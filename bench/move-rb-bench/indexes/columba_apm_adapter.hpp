/**
 * part of LukasNalbach/Move-r
 *
 * MIT License
 *
 * Copyright (c) Lukas Nalbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

// columba's index (flavor chosen by RUN_LENGTH_COMPRESSION, as in columba_adapter.hpp)
#ifdef RUN_LENGTH_COMPRESSION
#include <bmove/bmove.h>
#else
#include <fmindex/fmindex.h>
#endif
#include <definitions.h>
#include <indexhelpers.h>

#define umul128 gtl_umul128_shim
#include <misc/apm.hpp>
#include <algorithms/apm/hamming_distance.tpp>
#include <algorithms/apm/edit_distance.tpp>
#undef umul128

#include "index_adapter.hpp" // the interface this adapter implements (checked by the static_assert at the bottom)

/**
 * @brief Adapter that exposes columba's index (the run-length-compressed BMove or the FM-index) through the same
 *        query interface that move_rb provides, so that move_r's approximate-pattern-matching algorithms
 *        (apm_hamming and apm_edit) can drive columba's index unmodified -- i.e. the SAME algorithm run on
 *        columba's index instead of move_rb's, isolating the index from the algorithm.
 *
 * The adapter wraps a columba index and re-implements the bidirectional search context (extend / extend_phase /
 * extend_next) on top of columba's findRangesWithExtraChar routines (via the extendRange wrapper), and the locate
 * phase on top of columba's toehold: for the RLC flavor via BMove::locateCentered (a phi/phi^-1 walk), for the
 * FM-index via getTextPositionsFromSARange.
 */
class columba_apm_adapter {
  public:
#ifdef RUN_LENGTH_COMPRESSION
    using columba_index_t = BMove;
#else
    using columba_index_t = FMIndex;
#endif

    using pos_t = uint64_t;
    using pos_type = uint64_t;
    using sym_t = char;
    using sym_type = char;
    using i_sym_t = uint8_t;
    using inp_t = std::string;

    static constexpr bool supports_locate = true;
    static constexpr bool supports_multiple_locate = true;

    template <typename, cigar_mode_t> friend class apm_hamming;
    template <typename, cigar_mode_t> friend class apm_edit;

    struct locate_context_t;
    template <query_support_t query_support> struct extend_context_t;
    template <query_support_t query_support> struct search_context_t;

  protected:
    columba_index_t& _index; // the wrapped columba index (query routines are const)
    std::vector<uint8_t> _map_int; // maps a character to its effective symbol (1..sigma-1); 0 if absent/sentinel

  public:
    pos_t n; // length of the indexed text (including the sentinel)
    pos_t sigma; // size of the alphabet (including the sentinel)

    columba_apm_adapter(columba_index_t& index) : _index(index)
    {
        n = _index.getTextSize();
        sigma = _index.alphabetSize();

        _map_int.assign(256, 0);
        for (int32_t c = 0; c < 256; c++) {
            int32_t i = _index.charToIdx((char) c);
            _map_int[c] = i >= 1 ? (uint8_t) i : 0;
        }
    }

    inline columba_index_t& index() const { return _index; }

    inline i_sym_t map_symbol(sym_t sym) const { return _map_int[(uint8_t) sym]; }
    inline sym_t unmap_symbol(i_sym_t sym) const { return (sym_t) _index.idxToChar((int32_t) sym); }
    inline const std::vector<uint8_t>& map_int() const { return _map_int; }
    inline const columba_apm_adapter& forward_index() const { return *this; }
    inline pos_t input_size() const { return n - 1; }

    template <query_support_t query_support>
    inline search_context_t<query_support> empty_context() const;

    inline pos_t count(const inp_t& P) const;
    template <typename report_fnc_t>
    inline void locate(const inp_t& P, report_fnc_t report) const;

    inline pos_t count_hamming_dist(const inp_t& P, const search_scheme_t& scheme) const;
    template <distance_metric_t dist_metr, typename report_fnc_t>
    inline void locate(const inp_t& P, const search_scheme_t& scheme, report_fnc_t report) const;
    template <distance_metric_t dist_metr>
    inline std::vector<aprx_occ_t<pos_t>> locate(const inp_t& P, const search_scheme_t& scheme) const;
};

// ############################# search_context_t #############################

template <query_support_t query_support>
struct columba_apm_adapter::search_context_t {
    template <typename, cigar_mode_t> friend class apm_hamming;
    template <typename, cigar_mode_t> friend class apm_edit;
    friend struct columba_apm_adapter::locate_context_t;
    template <query_support_t _qs> friend struct columba_apm_adapter::extend_context_t;

    // locate-only bookkeeping (mirrors move_rb's)
    struct locate_sample_t {
        direction_t d = LEFT;
        pos_t dpth = 0;
        pos_t shft = 0;
        bool rprtd = false;
    };

    const columba_apm_adapter* idx = nullptr;

    SARangePair ranges; // columba's forward + reverse SA-ranges (and toehold) of the matched pattern
    pos_t m = 0; // length of the currently matched pattern
    pos_t err = 0; // number of errors (approximate pattern matching only)
    sym_t sym_lst = 0; // last-added symbol
    direction_t dir_lst = NO_DIR; // last performed extension direction

    [[no_unique_address]] std::conditional_t<query_support == LOCATE, locate_sample_t, std::monostate> s;

    using extend_res_t = std::tuple<search_context_t, bool>;

    // LEFT (prepend) maps to columba's BACKWARD extension, RIGHT (append) to FORWARD
    static constexpr Direction col_dir(direction_t dir) { return dir == LEFT ? BACKWARD : FORWARD; }

    search_context_t() {}
    search_context_t(sym_t sym) { sym_lst = sym; }
    search_context_t(const columba_apm_adapter& idx) : idx(&idx) { reset(); }

    inline void reset()
    {
        ranges = idx->index().getCompleteRange();
        m = 0;
        err = 0;
        sym_lst = sym_t(0);
        dir_lst = NO_DIR;
        if constexpr (query_support == LOCATE) s = locate_sample_t {};
    }

    inline sym_t last_added_symbol() const { return sym_lst; }
    inline direction_t last_direction() const { return dir_lst; }
    inline pos_t length() const { return m; }
    inline pos_t num_occ() const { return ranges.getRangeSA().width(); }
    inline pos_t errors() const { return err; }
    inline void set_errors(pos_t err) { this->err = err; }

    inline pos_t depth() const requires(query_support == LOCATE) { return s.dpth; }
    inline void set_depth(pos_t depth) requires(query_support == LOCATE) { s.dpth = depth; }
    inline pos_t shift() const requires(query_support == LOCATE) { return s.shft; }
    inline void set_shift(pos_t shft) requires(query_support == LOCATE) { s.shft = shft; }
    inline bool reported() const requires(query_support == LOCATE) { return s.rprtd; }
    inline void set_reported(bool reported) requires(query_support == LOCATE) { s.rprtd = reported; }

    inline locate_context_t locate_phase() const requires(query_support == LOCATE);
    inline extend_context_t<query_support> extend_phase();

    inline std::tuple<pos_t, pos_t> forward_sa_interval() const
    {
        const Range& r = ranges.getRangeSA();
        return {(pos_t) r.getBegin(), (pos_t) r.getEnd() - 1};
    }
    inline std::tuple<pos_t, pos_t> backward_sa_interval() const
    {
        const Range& r = ranges.getRangeSARev();
        return {(pos_t) r.getBegin(), (pos_t) r.getEnd() - 1};
    }

    bool operator==(const search_context_t& other) const
    {
        auto [b, e] = forward_sa_interval();
        auto [ob, oe] = other.forward_sa_interval();
        bool res = b == ob && e == oe;
        if constexpr (query_support == LOCATE) res = res && s.shft == other.s.shft && s.dpth == other.s.dpth;
        return res;
    }

    bool operator<(const search_context_t& other) const
    {
        auto [b, e] = forward_sa_interval();
        auto [ob, oe] = other.forward_sa_interval();
        if (b != ob) return b < ob;
        if (e != oe) return e > oe;
        if constexpr (query_support == LOCATE) {
            if (s.shft != other.s.shft) return s.shft < other.s.shft;
            if (s.dpth != other.s.dpth) return s.dpth < other.s.dpth;
        }
        return m < other.m;
    }

    struct hash {
        inline static pos_t operator()(const search_context_t& ctx)
        {
            auto [b, e] = ctx.forward_sa_interval();
            pos_t h = pos_hash<pos_t>(b);
            hash_combine<pos_t>(h, pos_hash<pos_t>(e));
            if constexpr (query_support == LOCATE) {
                hash_combine<pos_t>(h, pos_hash<pos_t>(ctx.s.shft));
                hash_combine<pos_t>(h, pos_hash<pos_t>(ctx.s.dpth));
            }
            return h;
        }
    };

    inline extend_res_t extend(sym_t sym, direction_t dir) const
    {
        int32_t alph = idx->index().charToIdx(sym);
        if (alph < 1) [[unlikely]] return {search_context_t {}, false}; // sentinel (0) or absent (-1)

        SARangePair child;
        if (!idx->index().extendRange((length_t) alph, col_dir(dir), ranges, child))
            return {search_context_t {}, false};

        search_context_t ctx = *this;
        ctx.ranges = child;
        ctx.m = m + 1;
        ctx.err = 0;
        ctx.sym_lst = sym;
        ctx.dir_lst = dir;
        if constexpr (query_support == LOCATE) {
            ctx.s.d = dir;
            ctx.s.dpth = s.dpth + 1;
            ctx.s.shft = s.shft;
            ctx.s.rprtd = false;
        }
        return {ctx, true};
    }
};

// ############################# extend_context_t #############################

template <query_support_t query_support>
struct columba_apm_adapter::extend_context_t {
    friend struct columba_apm_adapter::search_context_t<query_support>;

    const columba_apm_adapter* idx = nullptr;
    search_context_t<query_support>* ctx = nullptr;
    direction_t dir = NO_DIR;

    int32_t count = 0;
    int32_t cursor = 0;
    std::array<SARangePair, 256> children;
    std::array<sym_t, 256> chars;

    extend_context_t() {}
    extend_context_t(const columba_apm_adapter& idx, search_context_t<query_support>& ctx) : idx(&idx), ctx(&ctx) {}

    sym_t current_symbol() const { return cursor < count ? chars[cursor] : sym_t(0); }
    bool can_extend() const { return cursor < count; }

    extend_context_t& prepare_extend_all(direction_t dir)
    {
        this->dir = dir;
        count = 0;
        cursor = 0;

        // columba has no batched extension, so try every real alphabet symbol (1..sigma-1)
        const Direction cd = search_context_t<query_support>::col_dir(dir);
        for (length_t a = 1; a < idx->sigma; a++) {
            SARangePair child;
            if (idx->index().extendRange(a, cd, ctx->ranges, child)) {
                children[count] = child;
                chars[count] = idx->index().idxToChar((int32_t) a);
                count++;
            }
        }
        return *this;
    }

    search_context_t<query_support> extend_next()
    {
        int32_t j = cursor++;
        search_context_t<query_support> nxt = *ctx;
        nxt.ranges = children[j];
        nxt.m = ctx->m + 1;
        nxt.err = 0;
        nxt.sym_lst = chars[j];
        nxt.dir_lst = dir;
        if constexpr (query_support == LOCATE) {
            nxt.s.d = dir;
            nxt.s.dpth = ctx->s.dpth + 1;
            nxt.s.shft = ctx->s.shft;
            nxt.s.rprtd = false;
        }
        return nxt;
    }
};

// ############################# locate_context_t #############################

struct columba_apm_adapter::locate_context_t {
  protected:
    const columba_apm_adapter* idx = nullptr;
    const search_context_t<LOCATE>* ctx = nullptr;

    pos_t occ_rem = 0;
    bool centered = false;
    pos_t center_sa = 0; // SA-index of the centered occurrence
    length_t first_pos = 0; // text position of the centered occurrence
    std::vector<length_t> left; // text positions left of the center (descending SA index)
    std::vector<length_t> right; // text positions right of the center (ascending SA index)

    inline void enumerate()
    {
        if (centered) return;
        centered = true;
#ifdef RUN_LENGTH_COMPRESSION
        center_sa = (pos_t) idx->index().locateCentered(ctx->ranges, first_pos, left, right);
#else
        // FM-index: enumerate all positions of the range; use the range begin as the (arbitrary) center
        const Range& r = ctx->ranges.getRangeSA();
        center_sa = (pos_t) r.getBegin();
        std::vector<length_t> positions;
        idx->index().getTextPositionsFromSARange(ctx->ranges, positions);
        if (!positions.empty()) {
            first_pos = positions[0];
            right.assign(positions.begin() + 1, positions.end()); // remaining positions to the right of the center
        }
#endif
    }

  public:
    locate_context_t() {}
    locate_context_t(const search_context_t<LOCATE>& ctx) : idx(ctx.idx), ctx(&ctx) { occ_rem = ctx.num_occ(); }

    inline pos_t num_occ_rem() const { return occ_rem; }

    inline pos_t compute_center()
    {
        enumerate();
        return center_sa;
    }

    template <typename report_fnc_t>
    inline void locate(report_fnc_t report)
    {
        static constexpr bool report_pos = function_traits<report_fnc_t>::arity > 1;
        enumerate();
        pos_t shft = ctx->s.shft;

        if constexpr (report_pos) report(center_sa, (pos_t) first_pos + shft);
        else report((pos_t) first_pos + shft);

        for (pos_t j = 0; j < left.size(); j++) {
            if constexpr (report_pos) report(center_sa - 1 - j, (pos_t) left[j] + shft);
            else report((pos_t) left[j] + shft);
        }
        for (pos_t j = 0; j < right.size(); j++) {
            if constexpr (report_pos) report(center_sa + 1 + j, (pos_t) right[j] + shft);
            else report((pos_t) right[j] + shft);
        }
        occ_rem = 0;
    }

    std::vector<pos_t> locate()
    {
        std::vector<pos_t> Occ;
        Occ.reserve(occ_rem);
        locate([&](pos_t occ) { Occ.emplace_back(occ); });
        return Occ;
    }
};

// ############################# out-of-line methods #############################

template <query_support_t query_support>
inline columba_apm_adapter::locate_context_t
columba_apm_adapter::search_context_t<query_support>::locate_phase() const requires(query_support == LOCATE)
{
    return locate_context_t(*this);
}

template <query_support_t query_support>
inline columba_apm_adapter::extend_context_t<query_support>
columba_apm_adapter::search_context_t<query_support>::extend_phase() { return extend_context_t<query_support>(*idx, *this); }

template <query_support_t query_support>
inline columba_apm_adapter::search_context_t<query_support> columba_apm_adapter::empty_context() const
{
    return search_context_t<query_support>(*this);
}

inline columba_apm_adapter::pos_t columba_apm_adapter::count(const inp_t& P) const
{
    auto ctx = empty_context<COUNT>();
    for (pos_t i = 0; i < P.size(); i++) {
        auto [ctx_nxt, ok] = ctx.extend(P[i], RIGHT);
        if (!ok) return 0;
        ctx = ctx_nxt;
    }
    return ctx.num_occ();
}

template <typename report_fnc_t>
inline void columba_apm_adapter::locate(const inp_t& P, report_fnc_t report) const
{
    auto ctx = empty_context<LOCATE>();
    for (pos_t i = 0; i < P.size(); i++) {
        auto [ctx_nxt, ok] = ctx.extend(P[i], RIGHT);
        if (!ok) return;
        ctx = ctx_nxt;
    }
    ctx.locate_phase().locate([&](pos_t occ) { report(occ); });
}

inline columba_apm_adapter::pos_t columba_apm_adapter::count_hamming_dist(const inp_t& P, const search_scheme_t& scheme) const
{
    return apm_hamming{*this}.count(P, scheme);
}

template <distance_metric_t dist_metr, typename report_fnc_t>
inline void columba_apm_adapter::locate(const inp_t& P, const search_scheme_t& scheme, report_fnc_t report) const
{
    if constexpr (dist_metr == HAMMING_DISTANCE) apm_hamming{*this}.locate(P, scheme, report);
    if constexpr (dist_metr == EDIT_DISTANCE) apm_edit{*this}.locate(P, scheme, report);
}

template <distance_metric_t dist_metr>
inline std::vector<aprx_occ_t<columba_apm_adapter::pos_t>>
columba_apm_adapter::locate(const inp_t& P, const search_scheme_t& scheme) const
{
    std::vector<aprx_occ_t<pos_t>> Occ;
    locate<dist_metr>(P, scheme, [&](aprx_occ_t<pos_t> occ) { Occ.emplace_back(occ); });
    return Occ;
}

static_assert(index_adapter<columba_apm_adapter>, "columba_apm_adapter must implement the index_adapter interface");
