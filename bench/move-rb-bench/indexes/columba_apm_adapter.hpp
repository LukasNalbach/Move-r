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

#ifndef MOVE_R_BENCH_MOVE_RB_BENCH_INDEXES_COLUMBA_APM_ADAPTER_HPP
#define MOVE_R_BENCH_MOVE_RB_BENCH_INDEXES_COLUMBA_APM_ADAPTER_HPP

#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

// ALPHABET is set to 5 by the columba targets (see CMakeLists.txt); the fallback is for tools that parse this header
// on its own
#ifndef ALPHABET
#define ALPHABET 5
#endif

// columba's index (selected by RUN_LENGTH_COMPRESSION, see columba_adapter.hpp)
#ifdef RUN_LENGTH_COMPRESSION
#include <bmove/bmove.h>
#else
#include <fmindex/fmindex.h>
#endif
#include <definitions.h>
#include <indexhelpers.h>

// renames umul128 in the headers of Move-r (the name clashes with columba)
#define umul128 gtl_umul128_shim
#include <misc/apm.hpp>
#include <misc/portable_layout.hpp>
#include <algorithms/move_rb/apm/hamming.hpp>
#include <algorithms/move_rb/apm/edit.hpp>
#undef umul128

#include "index_adapter.hpp"

/**
 * @brief columba's index (the run-length compressed b-move or the FM-index) with the query interface of move_rb
 *        (see index_adapter), so that the approximate-pattern-matching algorithms of move_rb (apm_hamming and
 *        apm_edit) run on it
 *
 * A search context stores columba's SARangePair (the suffix array intervals of P in T and of P^R in T^R); it is
 * extended with columba's extendRange. The locate phase reports the occurrences around the toehold computed by
 * BMove::locateCentered (b-move), or all occurrences of the interval computed by getTextPositionsFromSARange,
 * starting with the first one (FM-index).
 */
class columba_apm_adapter {
  public:
#ifdef RUN_LENGTH_COMPRESSION
    using columba_index_t = BMove;
#else
    using columba_index_t = FMIndex;
#endif

    using inp_t = std::string; // type of a pattern

    static constexpr bool supports_locate = true;
    static constexpr bool supports_multiple_locate = true;
    // the number of byte values
    static constexpr uint64_t num_byte_values = uint64_t { std::numeric_limits<uint8_t>::max() } + 1;

    template <typename> friend class apm_hamming;
    template <typename> friend class apm_edit;

    struct locate_context_t;
    template <query_support_t query_support> struct extend_context_t;
    template <query_support_t query_support> struct search_context_t;

  protected:
    columba_index_t& _index; // the columba index
    // maps the symbols to the effective alphabet: [0..255] _map_int[c] = columba's index of c in [1..sigma-1], or 0 if
    // c is the terminator or does not occur in T
    std::vector<uint8_t> _map_int;

  public:
    uint64_t n; // length of the text, including the terminator
    uint64_t sigma; // size of the alphabet, including the terminator

    /**
     * @brief creates an adapter for a columba index
     * @param index the columba index
     */
    columba_apm_adapter(columba_index_t& index) : _index(index)
    {
        n = _index.getTextSize();
        sigma = _index.alphabetSize();
        _map_int.assign(num_byte_values, 0);

        for (int32_t c = 0; c < int32_t(num_byte_values); c++) {
            int32_t i = _index.charToIdx((char) c);
            _map_int[c] = i >= 1 ? (uint8_t) i : 0;
        }
    }

    // returns the columba index
    inline columba_index_t& index() const { return _index; }

    // returns the symbol of the effective alphabet of sym (0 if sym does not occur in T)
    inline uint8_t map_symbol(uint8_t sym) const { return _map_int[(uint8_t) sym]; }

    // returns the symbol of T of a symbol of the effective alphabet
    inline uint8_t unmap_symbol(uint8_t sym) const { return (uint8_t) _index.idxToChar((int32_t) sym); }

    // returns the map from the symbols to the effective alphabet (see _map_int)
    inline const std::vector<uint8_t>& map_int() const { return _map_int; }

    // returns the index that answers exact count and locate queries (the adapter itself)
    inline const columba_apm_adapter& forward_index() const { return *this; }

    // returns the length of the text without the terminator
    inline uint64_t input_size() const { return n - 1; }

    // the following methods are defined after the nested types

    // returns an empty search context
    template <query_support_t query_support>
    inline search_context_t<query_support> empty_context() const;

    // returns the number of occurrences of P in T
    inline uint64_t count(const inp_t& P) const;

    // reports the occurrences of P in T as report(i, occurrence), where i is the suffix array position
    template <typename report_fnc_t>
    inline void locate(const inp_t& P, report_fnc_t report) const;

    // returns the number of occurrences of P in T with at most k mismatches (k is given by the search scheme)
    inline uint64_t count_hamming_dist(const inp_t& P, const search_scheme_t& scheme, apm_params_t params = {}) const;

    // appends the occurrences of P in T with at most k errors w.r.t. dist_metr to occ
    inline void locate(const inp_t& P, const search_scheme_t& scheme, distance_metric_t dist_metr,
        std::vector<aprx_occ_t<>>& occ, apm_params_t params = {}) const;

    // returns the occurrences of P in T with at most k errors w.r.t. dist_metr
    inline std::vector<aprx_occ_t<>> locate(const inp_t& P, const search_scheme_t& scheme,
        distance_metric_t dist_metr, apm_params_t params = {}) const;
};

// ############################# search_context_t #############################

// a pattern P matched in the columba index: its suffix array intervals, its length and (LOCATE) the members of apm_edit
template <query_support_t query_support>
struct columba_apm_adapter::search_context_t {
    template <typename> friend class apm_hamming;
    template <typename> friend class apm_edit;
    friend struct columba_apm_adapter::locate_context_t;
    template <query_support_t> friend struct columba_apm_adapter::extend_context_t;

    using locate_sample_t = adapter_locate_sample_t;
    using extend_res_t = std::tuple<search_context_t, bool>;

    const columba_apm_adapter* idx = nullptr; // the adapter
    SARangePair ranges; // the suffix array intervals of P in T and of P^R in T^R (end exclusive)
    uint64_t m = 0; // length of P
    uint64_t err = 0; // number of errors (approximate pattern matching only)
    uint8_t sym_lst = 0; // the symbol added last
    STRING_UTILS_NO_UNIQUE_ADDRESS std::conditional_t<query_support == LOCATE, locate_sample_t, std::monostate> s;

    // returns columba's direction of an extension: BACKWARD for LEFT (prepending), FORWARD for RIGHT (appending)
    static constexpr Direction col_dir(direction_t dir) { return dir == LEFT ? BACKWARD : FORWARD; }

    search_context_t() {}

    // creates an empty search context for the adapter idx
    search_context_t(const columba_apm_adapter& idx) : idx(&idx) { reset(); }

    // sets P to the empty string
    inline void reset()
    {
        ranges = idx->index().getCompleteRange();
        m = 0;
        err = 0;
        sym_lst = 0;
        if constexpr (query_support == LOCATE) s = locate_sample_t {};
    }

    inline uint8_t last_added_symbol() const { return sym_lst; }
    inline uint64_t pattern_length() const { return m; }
    inline uint64_t num_occ() const { return ranges.getRangeSA().width(); }
    inline uint64_t errors() const { return err; }
    inline void set_errors(uint64_t err) { this->err = err; }

    inline uint64_t depth() const requires(query_support == LOCATE) { return s.dpth; }
    inline void set_depth(uint64_t depth) requires(query_support == LOCATE) { s.dpth = depth; }
    inline uint64_t shift() const requires(query_support == LOCATE) { return s.shft; }
    inline void set_shift(uint64_t shft) requires(query_support == LOCATE) { s.shft = shft; }
    inline bool reported() const requires(query_support == LOCATE) { return s.rprtd; }
    inline void set_reported(bool reported) requires(query_support == LOCATE) { s.rprtd = reported; }

    // returns a locate context for the occurrences of P (defined after locate_context_t)
    inline locate_context_t locate_phase() const requires(query_support == LOCATE);

    // returns an extend context for the extensions of P (defined after extend_context_t)
    inline extend_context_t<query_support> extend_phase();

    // returns the suffix array interval of P in T
    inline std::tuple<uint64_t, uint64_t> forward_sa_interval() const
    {
        const Range& r = ranges.getRangeSA();
        return { (uint64_t) r.getBegin(), (uint64_t) r.getEnd() - 1 };
    }

    // returns the suffix array interval of P^R in T^R
    inline std::tuple<uint64_t, uint64_t> backward_sa_interval() const
    {
        const Range& r = ranges.getRangeSARev();
        return { (uint64_t) r.getBegin(), (uint64_t) r.getEnd() - 1 };
    }

    bool operator==(const search_context_t& other) const { return adapter_contexts_equal<query_support>(*this, other); }
    bool operator<(const search_context_t& other) const { return adapter_context_less<query_support>(*this, other); }

    struct hash {
        inline static uint64_t operator()(const search_context_t& ctx) { return adapter_context_hash<query_support>(ctx); }
    };

    /**
     * @brief extends P by sym in direction dir
     * @return the extended context and whether the extended pattern occurs in T
     */
    inline extend_res_t extend(uint8_t sym, direction_t dir) const
    {
        // reject the terminator (index 0) and the symbols that do not occur in T (index -1)
        int32_t alph = idx->index().charToIdx(sym);
        if (alph < 1) [[unlikely]] return { search_context_t {}, false };

        SARangePair extended_ranges;
        if (!idx->index().extendRange((length_t) alph, col_dir(dir), ranges, extended_ranges))
            return { search_context_t {}, false };

        search_context_t extended = extended_context<query_support>(*this, sym, dir);
        extended.ranges = extended_ranges;
        return { extended, true };
    }
};

// ############################# extend_context_t #############################

// the extensions of a search context in one direction (computed with extendRange for each symbol)
template <query_support_t query_support>
struct columba_apm_adapter::extend_context_t {
    friend struct columba_apm_adapter::search_context_t<query_support>;

    const columba_apm_adapter* idx = nullptr; // the adapter
    search_context_t<query_support>* ctx = nullptr; // the extended search context
    direction_t dir = NO_DIR; // the direction of the extensions

    int32_t count = 0; // number of extensions
    int32_t cursor = 0; // index of the next extension returned by extend_next()
    // the extensions: [0..count-1] suffix array intervals of the i-th extension
    std::array<SARangePair, num_byte_values> children;
    // the symbols of the extensions: [0..count-1] symbol of the i-th extension
    std::array<uint8_t, num_byte_values> chars;

    extend_context_t() {}

    extend_context_t(const columba_apm_adapter& idx, search_context_t<query_support>& ctx) : idx(&idx), ctx(&ctx) {}

    // returns whether there is an extension left
    bool can_extend() const { return cursor < count; }

    /**
     * @brief computes the extensions of the search context in direction dir (with extendRange for each symbol
     *        1..sigma-1)
     * @tparam scope must be NO_TERM (extensions by the terminator are not supported)
     */
    template <extension_scope_t scope = TERM>
    extend_context_t& prepare_extend_all(direction_t dir)
    {
        static_assert(scope == NO_TERM);

        this->dir = dir;
        count = 0;
        cursor = 0;
        const Direction col_dir = search_context_t<query_support>::col_dir(dir);

        for (length_t a = 1; a < idx->sigma; a++) {
            SARangePair child;

            if (idx->index().extendRange(a, col_dir, ctx->ranges, child)) {
                children[count] = child;
                chars[count] = idx->index().idxToChar((int32_t) a);
                count++;
            }
        }

        return *this;
    }

    // returns the search context of the next extension
    search_context_t<query_support> extend_next()
    {
        int32_t j = cursor++;
        search_context_t<query_support> extended = extended_context<query_support>(*ctx, chars[j], dir);
        extended.ranges = children[j];
        return extended;
    }
};

// ############################# locate_context_t #############################

// the occurrences of a matched pattern (computed with BMove::locateCentered or getTextPositionsFromSARange)
struct columba_apm_adapter::locate_context_t {
  protected:
    const columba_apm_adapter* idx = nullptr; // the adapter
    const search_context_t<LOCATE>* ctx = nullptr; // the search context of the pattern

    uint64_t occ_rem = 0; // number of occurrences not reported yet
    bool centered = false; // true <=> left and right have been computed
    uint64_t center_sa = 0; // c, the suffix array position of the toehold (FM-index: the beginning of the interval)
    length_t first_pos = 0; // SA[c]
    std::vector<length_t> left; // the occurrences before c in suffix array order: [0..x-1] left[j] = SA[c-1-j]
    // the occurrences after c in suffix array order: [0..y-1] right[j] = SA[c+1+j] (FM-index: the other occurrences,
    // in the order of getTextPositionsFromSARange)
    std::vector<length_t> right;

    // computes center_sa, first_pos, left and right (once)
    inline void enumerate()
    {
        if (centered) return;
        centered = true;
#ifdef RUN_LENGTH_COMPRESSION
        center_sa = (uint64_t) idx->index().locateCentered(ctx->ranges, first_pos, left, right);
#else
        const Range& r = ctx->ranges.getRangeSA();
        center_sa = (uint64_t) r.getBegin();
        std::vector<length_t> positions;
        idx->index().getTextPositionsFromSARange(ctx->ranges, positions);

        if (!positions.empty()) {
            first_pos = positions[0];
            right.assign(positions.begin() + 1, positions.end());
        }
#endif
    }

  public:
    locate_context_t() {}

    // creates a locate context for a search context
    locate_context_t(const search_context_t<LOCATE>& ctx) : idx(ctx.idx), ctx(&ctx) { occ_rem = ctx.num_occ(); }

    inline uint64_t num_occ_rem() const { return occ_rem; }

    // returns the suffix array position of the toehold
    inline uint64_t compute_center()
    {
        enumerate();
        return center_sa;
    }

    // reports the occurrences (see report_centered_occurrences)
    template <typename report_fnc_t>
    inline void locate(report_fnc_t report)
    {
        enumerate();
        report_centered_occurrences(report, center_sa, (uint64_t) first_pos, left, right, ctx->s.shft);
        occ_rem = 0;
    }

    // returns the occurrences
    std::vector<uint64_t> locate()
    {
        std::vector<uint64_t> occurrences;
        occurrences.reserve(occ_rem);
        locate([&](uint64_t occ) { occurrences.emplace_back(occ); });
        return occurrences;
    }
};

// ############################# methods defined after the nested types #############################

template <query_support_t query_support>
inline columba_apm_adapter::locate_context_t
columba_apm_adapter::search_context_t<query_support>::locate_phase() const requires(query_support == LOCATE)
{
    return locate_context_t(*this);
}

template <query_support_t query_support>
inline columba_apm_adapter::extend_context_t<query_support>
columba_apm_adapter::search_context_t<query_support>::extend_phase()
{
    return extend_context_t<query_support>(*idx, *this);
}

template <query_support_t query_support>
inline columba_apm_adapter::search_context_t<query_support> columba_apm_adapter::empty_context() const
{
    return search_context_t<query_support>(*this);
}

inline uint64_t columba_apm_adapter::count(const inp_t& P) const { return count_by_appending(*this, P); }

template <typename report_fnc_t>
inline void columba_apm_adapter::locate(const inp_t& P, report_fnc_t report) const
{
    locate_by_appending(*this, P, report);
}

inline uint64_t columba_apm_adapter::count_hamming_dist(const inp_t& P, const search_scheme_t& scheme,
    apm_params_t params) const
{
    return apm_hamming { *this, params }.count(P, scheme);
}

inline void columba_apm_adapter::locate(const inp_t& P, const search_scheme_t& scheme, distance_metric_t dist_metr,
    std::vector<aprx_occ_t<>>& occ, apm_params_t params) const
{
    auto collect = [&occ](aprx_occ_t<> o) { occ.emplace_back(o); };
    if (dist_metr == HAMMING_DISTANCE) apm_hamming { *this, params }.locate(P, scheme, collect);
    else apm_edit { *this, params }.locate(P, scheme, collect);
}

inline std::vector<aprx_occ_t<>> columba_apm_adapter::locate(const inp_t& P, const search_scheme_t& scheme,
    distance_metric_t dist_metr, apm_params_t params) const
{
    std::vector<aprx_occ_t<>> occ;
    locate(P, scheme, dist_metr, occ, params);
    return occ;
}

static_assert(index_adapter<columba_apm_adapter>, "columba_apm_adapter must implement the index_adapter interface");

#endif // MOVE_R_BENCH_MOVE_RB_BENCH_INDEXES_COLUMBA_APM_ADAPTER_HPP
