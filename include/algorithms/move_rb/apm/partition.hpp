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

/**
 * @file partition.hpp
 * @brief a partition of a pattern into parts with the search context of each part
 */

#ifndef MOVE_R_ALGORITHMS_MOVE_RB_APM_PARTITION_HPP
#define MOVE_R_ALGORITHMS_MOVE_RB_APM_PARTITION_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <ranges>
#include <span>
#include <vector>

#include <misc/apm.hpp>
#include <misc/search_schemes.hpp>
#include <misc/utils.hpp>
#include <algorithms/move_rb/apm/model.hpp>
#include <algorithms/move_rb/apm/step_prediction.hpp>

/**
 * @brief a partition of a pattern P into p parts, with the search context of each part: part i is
 *        P[beg(i), end(i)). The search context of a part is computed if no position of the part has a free symbol
 *        other than its own symbol (see apm_costs_t::has_free; always without an error model). The parts have equal
 *        lengths, are chosen by the dynamic partitioning of Columba (Renders, Depuydt, Fostier), or are chosen
 *        together with the direction of the search scheme such that the predicted number of steps of the searches is
 *        small (see apm_step_prediction).
 * @tparam idx_t the bidirectional index type
 * @tparam query_support whether the contexts support locate (LOCATE) or only count (COUNT)
 */
template <typename idx_t, query_support_t query_support>
class apm_partition {
public:
    using inp_t = typename idx_t::inp_t;
    using sym_t = std::ranges::range_value_t<inp_t>;
    // the search contexts of approximate pattern matching (see move_rb::apm_context_t)
    using ctx_t = typename idx_t::template apm_context_t<query_support>;

    // the table of the numbers of occurrences of the substrings of P: [i][l - 1] the number of occurrences of P[i, i + l)
    using occurrence_table_t = std::vector<std::vector<uint64_t>>;

    // whether the index has methods for its alphabet size, the number r of runs in its BWT and its text length n
    static constexpr bool can_optimize = requires(const idx_t& idx) {
        idx.alphabet_size();
        idx.num_bwt_runs();
        idx.input_size();
    };

protected:
    const idx_t& _index;

    std::vector<uint64_t> _bounds; // [0..p] part i is P[_bounds[i], _bounds[i + 1]); _bounds[0] = 0, _bounds[p] = m
    std::vector<ctx_t> _ctx;       // [0..p-1] the search context of part i (valid if _has_context[i] and _occurs[i])
    std::vector<uint8_t> _occurs;  // [0..p-1] 1 <=> part i occurs in the text (valid if _has_context[i])
    // [0..p-1] 1 <=> no position of part i has a free symbol other than its own, so _ctx[i] is the context of part i
    std::vector<uint8_t> _has_context;
    const apm_costs_t* _costs;     // the error model, or nullptr
    bool _mirrored = false;        // whether the mirror of the search scheme is used for the searches
    mutable std::vector<double> _sym_p; // see symbol_probabilities (computed by its first call)
    mutable std::vector<double> _distinct_share; // see distinct_share (extended by its calls)
    mutable occurrence_table_t _occ_table; // see occurrences (computed by its first call)
    mutable uint64_t _occ_steps = 0;       // the number of context extensions of occurrences
    std::unique_ptr<search_scheme_t> _mirror_scheme; // the mirror of the search scheme (set by partition_optimized)
    // the predictions of partition_optimized for the searches of the scheme ([0]) and of its mirror ([1])
    std::optional<apm_step_prediction> _prediction[2];

    // returns whether P[i] has a free symbol other than itself
    bool is_free(uint64_t i) const { return _costs != nullptr && _costs->has_free(i); }

    // returns the share of the occurrences of a string of length L with distinct continuations (see apm_distinct_shares)
    double distinct_share(uint64_t L) const
    {
        if (L >= _distinct_share.size()) {
            _distinct_share = apm_distinct_shares(_index, std::max<uint64_t>(L, 2 * _distinct_share.size()));
        }

        return _distinct_share[L];
    }

    /* returns whether part i is rare: it has a search context, and at most one of its occurrences has a distinct
       continuation by one symbol (occ * d(l + 1) / n <= 1 for its length l, see apm_distinct_shares) */
    bool is_rare(uint64_t i) const
    {
        return _has_context[i] && double(num_occ(i)) * distinct_share(_bounds[i + 1] - _bounds[i] + 1) <= 1;
    }

    // returns whether no position of P[beg, end) has a free symbol other than itself
    bool has_no_free_symbol(uint64_t beg, uint64_t end) const
    {
        for (uint64_t i = beg; i < end; i++) if (is_free(i)) return false;
        return true;
    }

    // whether the search context of a pattern can be looked up in context tables of the index (see move_rb::apm_context_of)
    static constexpr bool has_context_tables = requires(const idx_t& idx, bool& occurs) {
        idx.template apm_context_of<query_support>(std::span<const sym_t>(), occurs);
    };

    /**
     * @brief computes the search context of P[beg, end), with the context tables of the index if it has them
     * @param P the pattern
     * @param beg starting position of the part
     * @param end end position of the part (exclusive)
     * @param occurs set to whether the part occurs in the text
     * @return the search context of the part, or the context of its longest suffix that occurs
     */
    ctx_t context_of(const inp_t& P, uint64_t beg, uint64_t end, bool& occurs) const;

    // returns the number of occurrences of part i, or 0 if it does not occur (only if _has_context[i])
    uint64_t num_occ_of_part(uint64_t i) const { return _occurs[i] ? _ctx[i].num_occ() : 0; }

    // computes the search context of each part of P from _bounds
    void compute_contexts(const inp_t& P);

    /**
     * @brief returns the length of the initial parts of the dynamic partitioning: the maximum pattern length of the
     *        context tables of the index if it has them, else 1
     * @param m the pattern length
     * @param p the number of parts
     * @return the length, in [1, max(1, m / (2p))]
     */
    uint64_t initial_part_length(uint64_t m, uint64_t p) const;

    /**
     * @brief partitions the pattern into p parts of equal length and computes the search context of each part
     * @param P the pattern
     * @param p the number of parts
     */
    void partition_uniform(const inp_t& P, uint64_t p);

    /**
     * @brief partitions the pattern into p parts with few occurrences each and computes the search context
     *        of each part
     *
     * 1. Each part starts as an initial part of length initial_part_length(m, p): the outer ones at the two ends of
     *    the pattern, the inner ones evenly in between.
     * 2. The remaining symbols are assigned one by one: the part with the largest number of occurrences (those of
     *    the two outer parts count twice) is extended by one symbol, towards the neighbour with fewer occurrences
     *    if it can be extended in both directions. A part is not extended over a position i at which a symbol other
     *    than P[i] is free, and a part without a search context is not extended.
     * 3. If no part that can be extended occurs, the remaining symbols are assigned to the part on their left.
     *
     * @param P the pattern
     * @param p the number of parts
     */
    void partition_dynamic(const inp_t& P, uint64_t p);

    // partitions P into p parts of equal length, or with partition_dynamic unless params.partitioning is uniform
    void make_parts(const inp_t& P, uint64_t p, const apm_params_t& params);

    /**
     * @brief returns the bounds of the p parts of P that minimize the sum of predictor.part_steps over the parts, by
     *        dynamic programming over the bounds
     * @param m the pattern length
     * @param p the number of parts
     * @param predictor the prediction of the number of steps of the searches
     * @param occ_of occ_of(b, e) returns the number of occurrences of P[b, e)
     */
    template <typename occ_fnc_t>
    static std::vector<uint64_t> minimal_bounds(uint64_t m, uint64_t p, const apm_step_prediction& predictor, occ_fnc_t occ_of);

    /**
     * @brief partitions P into the parts of a search scheme or of its mirror (then _mirrored is set) such that the
     *        predicted number of steps of the searches (see apm_step_prediction) is small, and computes the search context
     *        of each part; the optimization counts its own steps (the context extensions of step 2, the cells of
     *        minimal_bounds and the steps of the predictions, see apm_step_prediction::work) and stops once they exceed
     *        the smallest prediction found so far. A part is rare if it has at most one occurrence with a distinct
     *        continuation (see is_rare; in step 1, a part without a search context is not rare).
     *
     * 1. P is partitioned with partition_dynamic. These parts are kept if every part is rare. They are also kept if the
     *    smaller prediction of the scheme and of the mirror is at most the largest number of context extensions of step 2
     *    (see occurrences_steps); then the direction with the smaller prediction is chosen.
     * 2. The numbers of occurrences of the substrings of P are computed (see occurrences). The number of occurrences
     *    of a part is estimated by the smallest number of occurrences of its maximal substrings without a free symbol;
     *    a substring longer than its row of the table is counted as its longest prefix in the row.
     * 3. For the scheme and for the mirror, a local search starts from the parts of step 1, from the uniform partition
     *    and from minimal_bounds (if the steps so far plus the p * m^2 cells of minimal_bounds are at most the smallest
     *    prediction). In each round, the bound between two parts whose move by one position lowers the prediction most
     *    is moved; the local search stops when no move lowers the prediction or when the best move lowers it by fewer
     *    steps than the round took. A bound next to a rare part is not moved.
     *
     * @param P the pattern
     * @param scheme the search scheme
     * @param metric the distance metric of the searches
     * @param params the query parameters (see apm_step_prediction)
     */
    void partition_optimized(const inp_t& P, const search_scheme_t& scheme, distance_metric_t metric,
        const apm_params_t& params);

public:
    /**
     * @brief partitions a pattern into p parts and computes the search context of each part
     * @param index the index to search in
     * @param P the pattern
     * @param p the number of parts (0 < p <= P.size())
     * @param params the query parameters (the partitioning strategy and the error model)
     */
    apm_partition(const idx_t& index, const inp_t& P, uint64_t p, const apm_params_t& params)
        : _index(index), _bounds(p + 1), _occurs(p, 0), _has_context(p, 1), _costs(params.costs)
    {
        assert(0 < p && p <= P.size());
        _ctx.reserve(p);
        for (uint64_t i = 0; i < p; i++) _ctx.emplace_back(index);
        make_parts(P, p, params);
    }

    /**
     * @brief partitions a pattern into the parts of a search scheme and computes the search context of each part; with
     *        _partition_optimized (if can_optimize), the parts are chosen with partition_optimized and the mirror of the
     *        scheme is used for the searches if mirrored(); otherwise, the parts are chosen as with the constructor that
     *        takes p
     * @param index the index to search in
     * @param P the pattern
     * @param scheme the search scheme (0 < scheme.p <= P.size())
     * @param metric the distance metric of the searches
     * @param params the query parameters (the partitioning strategy and the error model)
     */
    apm_partition(const idx_t& index, const inp_t& P, const search_scheme_t& scheme, distance_metric_t metric,
        const apm_params_t& params);

    /**
     * @brief partitions a pattern into the parts given by their bounds and computes the search context of each part
     * @param index the index to search in
     * @param P the pattern
     * @param bounds [0..p] part i is P[bounds[i], bounds[i + 1]) (bounds[0] = 0, bounds[p] = P.size())
     * @param params the query parameters (the error model)
     */
    apm_partition(const idx_t& index, const inp_t& P, std::vector<uint64_t> bounds, const apm_params_t& params)
        : _index(index), _bounds(std::move(bounds)), _occurs(_bounds.size() - 1, 0),
          _has_context(_bounds.size() - 1, 1), _costs(params.costs)
    {
        for (uint64_t i = 0; i + 1 < _bounds.size(); i++) _ctx.emplace_back(index);
        compute_contexts(P);
    }

    /**
     * @brief returns [0..255] the probability of each byte in the text (how often it occurs, divided by n) for an index
     *        with a byte alphabet, else an empty vector (computed by the first call)
     */
    const std::vector<double>& symbol_probabilities() const;

    /**
     * @brief returns the numbers of occurrences of the substrings of P that start at each position (computed by the
     *        first call): [i][l - 1] is the number of occurrences of P[i, i + l), for l = 1, 2, ... until the substring
     *        occurs at most once, the next position has free symbols (see is_free) or P ends; if position i has free
     *        symbols, [i][0] is the number of occurrences of all of them and the row ends
     * @param P the pattern
     */
    const occurrence_table_t& occurrences(const inp_t& P) const;

    /* returns the largest number of context extensions of occurrences for a pattern of length m: the row of position
       i has at most m - i entries */
    static double occurrences_steps(uint64_t m) { return double(m) * double(m + 1) / 2; }

    // returns the number of occurrences of P[beg, end) in the text
    uint64_t num_occ_of(const inp_t& P, uint64_t beg, uint64_t end) const
    {
        bool occurs;
        const ctx_t ctx = context_of(P, beg, end, occurs);
        return occurs ? ctx.num_occ() : 0;
    }

    // returns the number of parts
    uint64_t size() const { return _bounds.size() - 1; }

    // returns whether the mirror of the search scheme is used for the searches (see mirrored_scheme)
    bool mirrored() const { return _mirrored; }

    /* returns the prediction of partition_optimized for the searches of the scheme or of its mirror that are used (see
       mirrored), or nullptr if partition_optimized did not run */
    apm_step_prediction* prediction() { return _prediction[_mirrored].has_value() ? &*_prediction[_mirrored] : nullptr; }

    // returns the bounds of the parts ([0..p], part i is P[bounds()[i], bounds()[i + 1]))
    const std::vector<uint64_t>& bounds() const { return _bounds; }

    // returns the number of occurrences of part i, or the largest uint64_t if its context is not computed
    uint64_t num_occ(uint64_t i) const { return _has_context[i] ? num_occ_of_part(i) : std::numeric_limits<uint64_t>::max(); }

    // returns the starting position of part i
    uint64_t beg(uint64_t i) const { return _bounds[i]; }

    // returns the end position of part i (exclusive)
    uint64_t end(uint64_t i) const { return _bounds[i + 1]; }

    // returns the length of part i
    uint64_t length(uint64_t i) const { return _bounds[i + 1] - _bounds[i]; }

    // returns whether part i occurs in the text (true if its context is not computed)
    bool occurs(uint64_t i) const { return !_has_context[i] || _occurs[i]; }

    // returns whether the search context of part i is computed (see _has_context)
    bool has_context(uint64_t i) const { return _has_context[i]; }

    // returns the search context of part i (valid if has_context(i) and occurs(i))
    const ctx_t& context(uint64_t i) const { return _ctx[i]; }
};

#include <algorithms/move_rb/apm/partition.tpp>

#endif // MOVE_R_ALGORITHMS_MOVE_RB_APM_PARTITION_HPP
