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

#include <bmove.h>

#define umul128 gtl_umul128_shim
#include <misc/apm.hpp>
#include <algorithms/apm/hamming_distance.tpp>
#include <algorithms/apm/edit_distance.tpp>
#undef umul128

// the interface this adapter implements (checked by the static_assert at the bottom)
#include "index_adapter.hpp"

/**
 * @brief Adapter that exposes b-move's bidirectional move structure (BMove)
 *        through the same query interface that move_rb provides, so that move_r's
 *        approximate-pattern-matching algorithms (apm_hamming and
 *        apm_edit) can drive it unmodified.
 *
 * The adapter wraps a BMove index and re-implements the bidirectional search
 * context (extend / extend_phase / extend_next) on top of b-move's
 * findRangesWithExtraChar{Forward,Backward} routines (via BMove::extendRange),
 * and the locate phase on top of b-move's toehold (via BMove::locateCentered).
 */
class b_move_adapter {
  public:
    using pos_t = uint64_t; // index integer type (= b-move's length_t)
    using pos_type = uint64_t; // exposed for the generic APM helpers
    using sym_t = char; // value type
    using sym_type = char; // exposed for the generic APM helpers
    using i_sym_t = uint8_t; // internal (unsigned) symbol type
    using inp_t = std::string; // input container type

    static constexpr bool supports_locate = true; // the index supports locate
    static constexpr bool supports_multiple_locate = true; // it can locate multiple occurrences

    template <typename> friend class apm_hamming;
    template <typename> friend class apm_edit;

    // ############################# NESTED QUERY TYPES (declarations) #############################

    struct locate_context_t;
    template <query_support_t query_support> struct extend_context_t;
    template <query_support_t query_support> struct search_context_t;

    // ############################# INDEX STATE #############################

  protected:
    const BMove& _index; // the wrapped b-move index
    std::vector<uint8_t> _map_int; // maps a character to its (1-based) internal symbol; 0 if absent

  public:
    pos_t n; // length of the indexed text (including the sentinel)
    pos_t sigma; // size of the alphabet (including the sentinel)

    /**
     * @brief constructs an approximate-pattern-matching adapter for a b-move index
     * @param index the b-move index to query
     */
    b_move_adapter(const BMove& index) : _index(index)
    {
        n = _index.getTextSize();
        sigma = _index.alphabetSize();

        // b-move maps the sentinel to alphabet index 0 and the other characters
        // to 1..sigma-1, which is exactly the move_r map_int convention.
        _map_int.assign(256, 0);

        for (int c = 0; c < 256; c++) {
            int i = _index.charToIdx((char) c);
            _map_int[c] = (i < 0) ? 0 : (uint8_t) i;
        }
    }

    // ############################# ACCESS METHODS #############################

    /**
     * @brief returns the wrapped b-move index
     */
    inline const BMove& bmove() const { return _index; }

    /**
     * @brief maps a symbol to its corresponding internal effective-alphabet symbol
     * @return the internal symbol (0 if sym does not occur in the alphabet)
     */
    inline i_sym_t map_symbol(sym_t sym) const
    {
        int i = _index.charToIdx(sym);
        return i < 0 ? 0 : (i_sym_t) i;
    }

    /**
     * @brief maps an internal effective-alphabet symbol back to its input symbol
     */
    inline sym_t unmap_symbol(i_sym_t sym) const { return _index.idxToChar(sym); }

    /**
     * @brief returns the map from characters to internal symbols (used by the edit-distance matrix)
     */
    inline const std::vector<uint8_t>& map_int() const { return _map_int; }

    /**
     * @brief returns a reference to the forward index; the adapter answers exact
     *        count/locate queries itself, so it returns *this
     */
    inline const b_move_adapter& forward_index() const { return *this; }

    /**
     * @brief returns the size of the indexed text (without the sentinel), matching
     *        move_r's input_size() convention (n holds the length including the sentinel)
     */
    inline pos_t input_size() const { return n - 1; }

    // ############################# QUERY METHODS (declarations) #############################

    /**
     * @brief returns an empty query context bound to this index
     */
    template <query_support_t query_support>
    inline search_context_t<query_support> empty_context() const;

    /**
     * @brief counts the exact occurrences of P in the indexed text
     */
    inline pos_t count(const inp_t& P) const;

    /**
     * @brief locates the exact occurrences of P in the indexed text
     */
    template <typename report_fnc_t>
    inline void locate(const inp_t& P, report_fnc_t report) const;

    // ############################# APPROXIMATE PATTERN MATCHING (declarations) #############################

    /**
     * @brief counts a pattern with at most k mismatches (w.r.t. hamming distance)
     */
    inline pos_t count_hamming_dist(const inp_t& P, const search_scheme_t& scheme) const;

    /**
     * @brief locates a pattern with at most k errors (w.r.t. dist_metr)
     */
    template <distance_metric_t dist_metr, typename report_fnc_t>
    inline void locate(const inp_t& P, const search_scheme_t& scheme, report_fnc_t report) const;

    /**
     * @brief locates a pattern with at most k errors (w.r.t. dist_metr)
     */
    template <distance_metric_t dist_metr>
    inline std::vector<aprx_occ_t<pos_t>> locate(const inp_t& P, const search_scheme_t& scheme) const;
};

// ############################# search_context_t #############################

/**
 * @brief stores the variables needed to perform bidirectional pattern search and
 *        locate queries on the wrapped b-move index
 */
template <query_support_t query_support>
struct b_move_adapter::search_context_t {
    template <typename> friend class apm_hamming;
    template <typename> friend class apm_edit;
    friend struct b_move_adapter::locate_context_t;
    template <query_support_t _query_support> friend struct b_move_adapter::extend_context_t;

    // sample/locate information maintained during the search (LOCATE only)
    struct locate_sample_t {
        direction_t d = LEFT; // direction of the last extension (used to deduplicate)
        pos_t dpth = 0; // depth (additional length) of this context in the banded matrix
        pos_t shft = 0; // right shift (in the text) of the occurrences
        bool rprtd = false; // true <=> this context has already been reported
    };

    const b_move_adapter* idx = nullptr; // the index this context queries

    SARangePair ranges; // b-move ranges (forward + reverse SA range and toehold)
    pos_t m = 0; // length of the currently matched pattern
    pos_t err = 0; // number of errors (approximate pattern matching only)
    sym_t sym_lst = 0; // last-added symbol
    direction_t dir_lst = NO_DIR; // last performed extension direction

    std::conditional_t<query_support == LOCATE, locate_sample_t, std::monostate> s;

    using extend_res_t = std::tuple<search_context_t, bool>;

    search_context_t() {}

    /**
     * @brief constructs a new search context with last-added symbol sym
     */
    search_context_t(sym_t sym) { sym_lst = sym; }

    /**
     * @brief constructs a new query context for the index idx
     */
    search_context_t(const b_move_adapter& idx) : idx(&idx) { reset(); }

    /**
     * @brief resets the query context to an empty P
     */
    inline void reset()
    {
        ranges = idx->bmove().getCompleteRange();
        m = 0;
        err = 0;
        sym_lst = sym_t(0);
        dir_lst = NO_DIR;

        if constexpr (query_support == LOCATE) {
            s = locate_sample_t {};
        }
    }

    // ############################# accessors #############################

    inline sym_t last_added_symbol() const { return sym_lst; }
    inline direction_t last_direction() const { return dir_lst; }
    inline pos_t length() const { return m; }
    inline pos_t num_occ() const { return ranges.width(); }
    inline pos_t errors() const { return err; }
    inline void set_errors(pos_t err) { this->err = err; }

    inline pos_t depth() const requires(query_support == LOCATE) { return s.dpth; }
    inline void set_depth(pos_t depth) requires(query_support == LOCATE) { s.dpth = depth; }
    inline pos_t shift() const requires(query_support == LOCATE) { return s.shft; }
    inline void set_shift(pos_t shft) requires(query_support == LOCATE) { s.shft = shft; }
    inline bool reported() const requires(query_support == LOCATE) { return s.rprtd; }
    inline void set_reported(bool reported) requires(query_support == LOCATE) { s.rprtd = reported; }

    /**
     * @brief returns a locate context for this search context (defined out of line)
     */
    inline locate_context_t locate_phase() const requires(query_support == LOCATE);

    /**
     * @brief returns an extend context for this search context (defined out of line)
     */
    inline extend_context_t<query_support> extend_phase();

    /**
     * @brief returns the SA-interval in the (forward) text of the currently matched P
     */
    inline std::tuple<pos_t, pos_t> forward_sa_interval() const
    {
        const auto& r = ranges.getRangeSA();
        return {r.getBegin(), r.getEnd() - 1};
    }

    /**
     * @brief returns the SA-interval in the reversed text of the currently matched P^R
     */
    inline std::tuple<pos_t, pos_t> backward_sa_interval() const
    {
        const auto& r = ranges.getRangeSARev();
        return {r.getBegin(), r.getEnd() - 1};
    }

    template <direction_t dir>
    inline std::tuple<pos_t, pos_t> sa_interval() const
    {
        if constexpr (dir == LEFT) return forward_sa_interval();
        else return backward_sa_interval();
    }

    // ############################# comparison & hashing #############################

    bool operator==(const search_context_t& other) const
    {
        bool res = ranges.getRangeSA().getBegin() == other.ranges.getRangeSA().getBegin()
                && ranges.getRangeSA().getEnd() == other.ranges.getRangeSA().getEnd();

        if constexpr (query_support == LOCATE) {
            res = res && s.shft == other.s.shft && s.dpth == other.s.dpth;
        }

        return res;
    }

    bool operator<(const search_context_t& other) const
    {
        pos_t b = ranges.getRangeSA().getBegin();
        pos_t e = ranges.getRangeSA().getEnd();
        pos_t ob = other.ranges.getRangeSA().getBegin();
        pos_t oe = other.ranges.getRangeSA().getEnd();

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
            pos_t hash = pos_hash<pos_t>(b);
            hash_combine<pos_t>(hash, pos_hash<pos_t>(e));

            if constexpr (query_support == LOCATE) {
                hash_combine<pos_t>(hash, pos_hash<pos_t>(ctx.s.shft));
                hash_combine<pos_t>(hash, pos_hash<pos_t>(ctx.s.dpth));
            }

            return hash;
        }
    };

    // ############################# extension #############################

    /**
     * @brief extends the currently matched P with sym in direction dir
     * @return {extended context, whether the extended P occurs in the text}
     */
    inline extend_res_t extend(sym_t sym, direction_t dir) const
    {
        return dir == LEFT ? extend<LEFT>(sym) : extend<RIGHT>(sym);
    }

    /**
     * @brief extends the currently matched P with sym in direction dir
     * @return {extended context, whether the extended P occurs in the text}
     */
    template <direction_t dir>
    inline extend_res_t extend(sym_t sym) const
    {
        int pa = idx->bmove().charToIdx(sym);

        // sym is the sentinel (0) or does not occur in the alphabet (-1)
        if (pa <= 0) [[unlikely]] return {search_context_t {}, false};

        SARangePair child;
        Direction d = (dir == LEFT) ? BACKWARD : FORWARD;

        if (!idx->bmove().extendRange((length_t) pa, d, ranges, child)) {
            return {search_context_t {}, false};
        }

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

/**
 * @brief stores the precomputed children of a search context so that it can be
 *        extended with every applicable symbol in O(sigma) time via extend_next()
 */
template <query_support_t query_support>
struct b_move_adapter::extend_context_t {
    friend struct b_move_adapter::search_context_t<query_support>;

    const b_move_adapter* idx = nullptr; // the index this context queries
    search_context_t<query_support>* ctx = nullptr; // the search context being extended
    direction_t dir = NO_DIR; // direction the context is being extended in

    int count = 0; // number of applicable extensions
    int cursor = 0; // index of the next extension to apply
    std::array<SARangePair, 8> children; // the ranges of each applicable extension
    std::array<sym_t, 8> chars; // the symbol of each applicable extension

    extend_context_t() {}

    extend_context_t(const b_move_adapter& idx, search_context_t<query_support>& ctx)
        : idx(&idx), ctx(&ctx) {}

    /**
     * @brief returns the next symbol the context can be extended with
     */
    sym_t current_symbol() const { return cursor < count ? chars[cursor] : sym_t(0); }

    /**
     * @brief returns whether there is a symbol left to extend the context with
     */
    bool can_extend() const { return cursor < count; }

    /**
     * @brief prepares the context to be extended in direction dir; afterwards the
     *        search context can be extended in O(sigma) time via extend_next()
     * @return a reference to this extend context
     */
    extend_context_t& prepare_extend_all(direction_t dir)
    {
        this->dir = dir;
        count = 0;
        cursor = 0;

        Direction d = (dir == LEFT) ? BACKWARD : FORWARD;

        // skip the sentinel (alphabet index 0)
        for (int c = 1; c < (int) idx->sigma; c++) {
            SARangePair child;

            if (idx->bmove().extendRange((length_t) c, d, ctx->ranges, child)) {
                children[count] = child;
                chars[count] = idx->bmove().idxToChar(c);
                count++;
            }
        }

        return *this;
    }

    /**
     * @brief extends the context with the next applicable symbol
     */
    search_context_t<query_support> extend_next()
    {
        int j = cursor++;

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

/**
 * @brief structure storing the information for locating all occurrences of a search context
 */
struct b_move_adapter::locate_context_t {
  protected:
    const b_move_adapter* idx = nullptr; // the index this context queries
    const search_context_t<LOCATE>* ctx = nullptr; // the search context this locate context belongs to

    pos_t occ_rem = 0; // number of remaining occurrences to locate
    bool centered = false; // true <=> the occurrences have been enumerated already
    pos_t center_sa = 0; // SA-index of the centered (toehold) occurrence
    length_t first_pos = 0; // text position (SA value) of the centered occurrence
    std::vector<length_t> left; // text positions left of the center (descending SA index)
    std::vector<length_t> right; // text positions right of the center (ascending SA index)

    /**
     * @brief enumerates the occurrences of the context centered on the toehold occurrence
     */
    inline void enumerate()
    {
        if (!centered) {
            center_sa = idx->bmove().locateCentered(ctx->ranges, first_pos, left, right);
            centered = true;
        }
    }

  public:
    locate_context_t() {}

    /**
     * @brief constructs a new locate context for a given search context
     */
    locate_context_t(const search_context_t<LOCATE>& ctx) : idx(ctx.idx), ctx(&ctx)
    {
        occ_rem = ctx.num_occ();
    }

    /**
     * @brief returns the number of remaining (not yet reported) occurrences
     */
    inline pos_t num_occ_rem() const { return occ_rem; }

    /**
     * @brief computes the SA-index of the centered (toehold) occurrence of P
     */
    inline pos_t compute_center()
    {
        enumerate();
        return center_sa;
    }

    /**
     * @brief locates the remaining occurrences of the currently matched P; report is
     *        called with (occ) or, if it has arity > 1, with (sa_index, occ) per occurrence
     */
    template <typename report_fnc_t>
    inline void locate(report_fnc_t report)
    {
        static constexpr bool report_pos = function_traits<report_fnc_t>::arity > 1;
        enumerate();
        pos_t shft = ctx->s.shft;

        // centered (toehold) occurrence
        if constexpr (report_pos) report(center_sa, first_pos + shft);
        else report(first_pos + shft);

        // occurrences left of the center: SA[c-1], SA[c-2], ..., SA[begin]
        for (pos_t j = 0; j < left.size(); j++) {
            if constexpr (report_pos) report(center_sa - 1 - j, left[j] + shft);
            else report(left[j] + shft);
        }

        // occurrences right of the center: SA[c+1], SA[c+2], ..., SA[end-1]
        for (pos_t j = 0; j < right.size(); j++) {
            if constexpr (report_pos) report(center_sa + 1 + j, right[j] + shft);
            else report(right[j] + shft);
        }

        occ_rem = 0;
    }

    /**
     * @brief locates the remaining occurrences of the currently matched P
     * @return vector containing the occurrences
     */
    std::vector<pos_t> locate()
    {
        std::vector<pos_t> Occ;
        Occ.reserve(occ_rem);
        locate([&](pos_t occ) { Occ.emplace_back(occ); });
        return Occ;
    }
};

// ############################# out-of-line search_context_t methods #############################

template <query_support_t query_support>
inline b_move_adapter::locate_context_t
b_move_adapter::search_context_t<query_support>::locate_phase() const requires(query_support == LOCATE)
{
    return locate_context_t(*this);
}

template <query_support_t query_support>
inline b_move_adapter::extend_context_t<query_support>
b_move_adapter::search_context_t<query_support>::extend_phase()
{
    return extend_context_t<query_support>(*idx, *this);
}

// ############################# out-of-line b_move_adapter methods #############################

template <query_support_t query_support>
inline b_move_adapter::search_context_t<query_support> b_move_adapter::empty_context() const
{
    return search_context_t<query_support>(*this);
}

inline b_move_adapter::pos_t b_move_adapter::count(const inp_t& P) const
{
    auto ctx = empty_context<COUNT>();

    for (pos_t i = 0; i < P.size(); i++) {
        auto [ctx_nxt, result] = ctx.extend(P[i], RIGHT);
        if (!result) return 0;
        ctx = ctx_nxt;
    }

    return ctx.num_occ();
}

template <typename report_fnc_t>
inline void b_move_adapter::locate(const inp_t& P, report_fnc_t report) const
{
    auto ctx = empty_context<LOCATE>();

    for (pos_t i = 0; i < P.size(); i++) {
        auto [ctx_nxt, result] = ctx.extend(P[i], RIGHT);
        if (!result) return;
        ctx = ctx_nxt;
    }

    ctx.locate_phase().locate([&](pos_t occ) { report(occ); });
}

inline b_move_adapter::pos_t b_move_adapter::count_hamming_dist(const inp_t& P, const search_scheme_t& scheme) const
{
    return apm_hamming{*this}.count(P, scheme);
}

template <distance_metric_t dist_metr, typename report_fnc_t>
inline void b_move_adapter::locate(const inp_t& P, const search_scheme_t& scheme, report_fnc_t report) const
{
    if constexpr (dist_metr == HAMMING_DISTANCE) apm_hamming{*this}.locate(P, scheme, report);
    if constexpr (dist_metr == EDIT_DISTANCE) apm_edit{*this}.locate(P, scheme, report);
}

template <distance_metric_t dist_metr>
inline std::vector<aprx_occ_t<b_move_adapter::pos_t>> b_move_adapter::locate(const inp_t& P, const search_scheme_t& scheme) const
{
    std::vector<aprx_occ_t<pos_t>> Occ;
    locate<dist_metr>(P, scheme, [&](aprx_occ_t<pos_t> occ) { Occ.emplace_back(occ); });
    return Occ;
}

static_assert(index_adapter<b_move_adapter>,
    "b_move_adapter must implement the index_adapter interface");
