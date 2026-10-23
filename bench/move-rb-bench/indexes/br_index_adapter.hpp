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

#ifndef MOVE_R_BENCH_MOVE_RB_BENCH_INDEXES_BR_INDEX_ADAPTER_HPP
#define MOVE_R_BENCH_MOVE_RB_BENCH_INDEXES_BR_INDEX_ADAPTER_HPP

#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <tuple>
#include <stdexcept>
#include <variant>
#include <vector>

#include <misc/apm.hpp>
#include <misc/portable_layout.hpp>
#include <misc/index_cli.hpp>
#include <algorithms/move_rb/apm/hamming.hpp>
#include <algorithms/move_rb/apm/edit.hpp>

#include <br_index.hpp>

#include "index_adapter.hpp"

/**
 * @brief the bidirectional r-index (bri::br_index) with the query interface of move_rb (see index_adapter), so that
 *        the approximate-pattern-matching algorithms of move_rb (apm_hamming and apm_edit) run on it
 *
 * A search context stores the toehold sample of br_index (the suffix array intervals of P in T and of P^R in T^R
 * and one occurrence); it is extended with left_extension() and right_extension(). The locate phase reports the
 * occurrences around the toehold (see br_index::locate_sample_centered).
 */
class br_index_adapter {
  public:
    using br_index_t = bri::br_index<>;
    using inp_t = std::string; // type of a pattern

    static constexpr bool supports_locate = true;
    static constexpr bool supports_multiple_locate = true;
    // the number of byte values
    static constexpr uint64_t num_byte_values = uint64_t { std::numeric_limits<uint8_t>::max() } + 1;
    // the smallest internal symbol of br_index of a symbol of T (0: does not occur, 1: the terminator)
    static constexpr bri::uchar min_text_symbol = 2;

    template <typename> friend class apm_hamming;
    template <typename> friend class apm_edit;

    struct locate_context_t;
    template <query_support_t query_support> struct extend_context_t;
    template <query_support_t query_support> struct search_context_t;

  protected:
    br_index_t& _index; // the br-index (its query methods are not const)
    // maps the symbols to the effective alphabet: [0..255] _map_int[c] = (internal symbol of c in br_index) - 1, or 0
    // if c does not occur in T
    std::vector<uint8_t> _map_int;

  public:
    uint64_t n; // length of the text, including the terminator
    uint64_t sigma; // size of the alphabet, including the terminator

    /**
     * @brief creates an adapter for a br-index
     * @param index the br-index
     */
    br_index_adapter(br_index_t& index) : _index(index)
    {
        n = _index.bwt_size();
        sigma = _index.alphabet_size();

        // br_index's internal symbols are 1 (the terminator) and 2..sigma; the edit-distance matrix expects map_int[c]
        // in [1, sigma-1] for the symbols of T
        _map_int.assign(num_byte_values, 0);

        for (int32_t c = 0; c < int32_t(num_byte_values); c++) {
            bri::uchar i = _index.char_to_internal((bri::uchar) c);
            _map_int[c] = (i < min_text_symbol) ? 0 : (uint8_t) (i - 1);
        }
    }

    // returns the br-index
    inline br_index_t& brindex() const { return _index; }

    // returns the symbol of the effective alphabet of sym (0 if sym does not occur in T)
    inline uint8_t map_symbol(uint8_t sym) const { return _map_int[(uint8_t) sym]; }

    // returns the symbol of T of a symbol of the effective alphabet
    inline uint8_t unmap_symbol(uint8_t sym) const { return (uint8_t) _index.internal_to_char((bri::uchar) (sym + 1)); }

    // returns the map from the symbols to the effective alphabet (see _map_int)
    inline const std::vector<uint8_t>& map_int() const { return _map_int; }

    // returns the index that answers exact count and locate queries (the adapter itself)
    inline const br_index_adapter& forward_index() const { return *this; }

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

// a pattern P matched in the br-index: its toehold sample, its length and (LOCATE) the members of apm_edit
template <query_support_t query_support>
struct br_index_adapter::search_context_t {
    template <typename> friend class apm_hamming;
    template <typename> friend class apm_edit;
    friend struct br_index_adapter::locate_context_t;
    template <query_support_t> friend struct br_index_adapter::extend_context_t;

    using locate_sample_t = adapter_locate_sample_t;
    using extend_res_t = std::tuple<search_context_t, bool>;

    const br_index_adapter* idx = nullptr; // the adapter
    bri::br_sample sample; // the toehold sample of P (range: SA-interval of P, rangeR: SA-interval of P^R in T^R)
    uint64_t m = 0; // length of P
    uint64_t err = 0; // number of errors (approximate pattern matching only)
    uint8_t sym_lst = 0; // the symbol added last
    STRING_UTILS_NO_UNIQUE_ADDRESS std::conditional_t<query_support == LOCATE, locate_sample_t, std::monostate> s;

    search_context_t() {}

    // creates an empty search context for the adapter idx
    search_context_t(const br_index_adapter& idx) : idx(&idx) { reset(); }

    // sets P to the empty string
    inline void reset()
    {
        sample = idx->brindex().get_initial_sample(false);
        m = 0;
        err = 0;
        sym_lst = 0;
        if constexpr (query_support == LOCATE) s = locate_sample_t {};
    }

    inline uint8_t last_added_symbol() const { return sym_lst; }
    inline uint64_t pattern_length() const { return m; }
    inline uint64_t num_occ() const { return sample.range.second + 1 - sample.range.first; }
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
        return { sample.range.first, sample.range.second };
    }

    // returns the suffix array interval of P^R in T^R
    inline std::tuple<uint64_t, uint64_t> backward_sa_interval() const
    {
        return { sample.rangeR.first, sample.rangeR.second };
    }

    bool operator==(const search_context_t& other) const { return adapter_contexts_equal<query_support>(*this, other); }
    bool operator<(const search_context_t& other) const { return adapter_context_less<query_support>(*this, other); }

    struct hash {
        inline uint64_t operator()(const search_context_t& ctx) const { return adapter_context_hash<query_support>(ctx); }
    };

    /**
     * @brief extends P by sym in direction dir
     * @return the extended context and whether the extended pattern occurs in T
     */
    inline extend_res_t extend(uint8_t sym, direction_t dir) const
    {
        // reject the terminator (internal symbol 1) and the symbols that do not occur in T (internal symbol 0)
        if (idx->brindex().char_to_internal((bri::uchar) sym) < min_text_symbol) [[unlikely]]
            return { search_context_t {}, false };

        bri::br_sample extended_sample = dir == LEFT
            ? idx->brindex().left_extension((bri::uchar) sym, sample)
            : idx->brindex().right_extension((bri::uchar) sym, sample);

        if (extended_sample.is_invalid()) return { search_context_t {}, false };
        search_context_t extended = extended_context<query_support>(*this, sym, dir);
        extended.sample = extended_sample;
        return { extended, true };
    }
};

// ############################# extend_context_t #############################

// the extensions of a search context in one direction (computed with left_extension_all or right_extension_all)
template <query_support_t query_support>
struct br_index_adapter::extend_context_t {
    friend struct br_index_adapter::search_context_t<query_support>;

    const br_index_adapter* idx = nullptr; // the adapter
    search_context_t<query_support>* ctx = nullptr; // the extended search context
    direction_t dir = NO_DIR; // the direction of the extensions

    int32_t count = 0; // number of extensions
    int32_t cursor = 0; // index of the next extension returned by extend_next()
    // the extensions: [0..count-1] toehold sample of the i-th extension
    std::array<bri::br_sample, num_byte_values> children;
    // the symbols of the extensions: [0..count-1] symbol of the i-th extension
    std::array<uint8_t, num_byte_values> chars;

    extend_context_t() {}

    extend_context_t(const br_index_adapter& idx, search_context_t<query_support>& ctx)
        : idx(&idx), ctx(&ctx) {}

    // returns whether there is an extension left
    bool can_extend() const { return cursor < count; }

    /**
     * @brief computes the extensions of the search context in direction dir
     * @tparam scope must be NO_TERM (extensions by the terminator are not supported)
     */
    template <extension_scope_t scope = TERM>
    extend_context_t& prepare_extend_all(direction_t dir)
    {
        static_assert(scope == NO_TERM);

        this->dir = dir;
        count = 0;
        cursor = 0;

        auto report = [this](bri::uchar c, const bri::br_sample& child) {
            children[count] = child;
            chars[count] = (uint8_t) c;
            count++;
        };

        if (dir == LEFT) idx->brindex().left_extension_all(ctx->sample, report);
        else idx->brindex().right_extension_all(ctx->sample, report);

        return *this;
    }

    // returns the search context of the next extension
    search_context_t<query_support> extend_next()
    {
        int32_t j = cursor++;
        search_context_t<query_support> extended = extended_context<query_support>(*ctx, chars[j], dir);
        extended.sample = children[j];
        return extended;
    }
};

// ############################# locate_context_t #############################

// the occurrences of a matched pattern, computed with br_index::locate_sample_centered
struct br_index_adapter::locate_context_t {
  protected:
    const br_index_adapter* idx = nullptr; // the adapter
    const search_context_t<LOCATE>* ctx = nullptr; // the search context of the pattern

    uint64_t occ_rem = 0; // number of occurrences not reported yet
    bool centered = false; // true <=> left and right have been computed
    uint64_t center_sa = 0; // c, the suffix array position of the toehold
    bri::ulint first_pos = 0; // SA[c]
    std::vector<bri::ulint> left; // the occurrences before c in suffix array order: [0..x-1] left[j] = SA[c-1-j]
    std::vector<bri::ulint> right; // the occurrences after c in suffix array order: [0..y-1] right[j] = SA[c+1+j]

    // computes center_sa, first_pos, left and right (once)
    inline void enumerate()
    {
        if (!centered) {
            center_sa = idx->brindex().locate_sample_centered(ctx->sample, first_pos, left, right);
            centered = true;
        }
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
        report_centered_occurrences(report, center_sa, first_pos, left, right, ctx->s.shft);
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
inline br_index_adapter::locate_context_t
br_index_adapter::search_context_t<query_support>::locate_phase() const requires(query_support == LOCATE)
{
    return locate_context_t(*this);
}

template <query_support_t query_support>
inline br_index_adapter::extend_context_t<query_support> br_index_adapter::search_context_t<query_support>::extend_phase()
{
    return extend_context_t<query_support>(*idx, *this);
}

template <query_support_t query_support>
inline br_index_adapter::search_context_t<query_support> br_index_adapter::empty_context() const
{
    return search_context_t<query_support>(*this);
}

inline uint64_t br_index_adapter::count(const inp_t& P) const { return count_by_appending(*this, P); }

template <typename report_fnc_t>
inline void br_index_adapter::locate(const inp_t& P, report_fnc_t report) const
{
    locate_by_appending(*this, P, report);
}

inline uint64_t br_index_adapter::count_hamming_dist(const inp_t& P, const search_scheme_t& scheme,
    apm_params_t params) const
{
    return apm_hamming { *this, params }.count(P, scheme);
}

inline void br_index_adapter::locate(const inp_t& P, const search_scheme_t& scheme, distance_metric_t dist_metr,
    std::vector<aprx_occ_t<>>& occ, apm_params_t params) const
{
    auto collect = [&occ](aprx_occ_t<> o) { occ.emplace_back(o); };
    if (dist_metr == HAMMING_DISTANCE) apm_hamming { *this, params }.locate(P, scheme, collect);
    else apm_edit { *this, params }.locate(P, scheme, collect);
}

inline std::vector<aprx_occ_t<>> br_index_adapter::locate(const inp_t& P, const search_scheme_t& scheme,
    distance_metric_t dist_metr, apm_params_t params) const
{
    std::vector<aprx_occ_t<>> occ;
    locate(P, scheme, dist_metr, occ, params);
    return occ;
}

static_assert(index_adapter<br_index_adapter>, "br_index_adapter must implement the index_adapter interface");

// ############################# br_index_native #############################

/**
 * @brief the approximate pattern matching of the br-index itself (hamming distance only): the pattern is divided into
 *        k + 1 parts, and the mismatches are placed with bidirectional search (search_with_mismatch) before the
 *        resulting samples are counted or located
 */
class br_index_native {
  public:
    using br_index_t = bri::br_index<>;

  protected:
    br_index_t _index; // the br-index (its query methods are not const)

  public:
    /**
     * @brief loads a br-index
     * @param path path of the br-index file (.bri)
     */
    br_index_native(const std::string& path)
    {
        direct_ifstream in(path);
        if (!in.good()) throw std::runtime_error("br-index file " + path + " not found");
        load_index_with_progress(_index, in, "native br-index");
        in.close();
    }

    // returns the length of the text without the terminator
    inline uint64_t input_size() { return _index.text_size(); }

    // returns the length of the text including the terminator
    inline uint64_t bwt_size() { return _index.bwt_size(); }

    // returns the number of occurrences of P with at most k mismatches
    inline uint64_t count_hamming_dist(const std::string& P, uint64_t k)
    {
        auto samples = _index.search_with_mismatch(P, (bri::ulint) k);
        return _index.count_samples(samples);
    }

    // locates the occurrences of P with at most k mismatches and returns their number (they are not deduplicated)
    inline uint64_t locate_hamming_dist(const std::string& P, uint64_t k)
    {
        auto samples = _index.search_with_mismatch(P, (bri::ulint) k);
        return _index.locate_samples(samples).size();
    }
};

#endif // MOVE_R_BENCH_MOVE_RB_BENCH_INDEXES_BR_INDEX_ADAPTER_HPP
