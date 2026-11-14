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
    using ctx_t = typename idx_t::template search_context_t<query_support>;

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
    mutable occurrence_table_t _occ_table; // see occurrences (computed by its first call)
    mutable uint64_t _occ_steps = 0;       // the number of context extensions of occurrences

    // returns whether P[i] has a free symbol other than itself
    bool is_free(uint64_t i) const { return _costs != nullptr && _costs->has_free(i); }

    // returns whether no position of P[beg, end) has a free symbol other than itself
    bool has_no_free_symbol(uint64_t beg, uint64_t end) const
    {
        for (uint64_t i = beg; i < end; i++) if (is_free(i)) return false;
        return true;
    }

    // whether the search context of a pattern can be looked up in context tables of the index
    static constexpr bool has_context_tables = requires(const idx_t& idx, bool& occurs) {
        idx.template search<query_support>(std::span<const sym_t>(), occurs);
    };

    /**
     * @brief computes the search context of P[beg, end), with the context tables of the index if it has them
     * @param P the pattern
     * @param beg starting position of the part
     * @param end end position of the part (exclusive)
     * @param occurs set to whether the part occurs in the text
     * @return the search context of the part, or the context of its longest suffix that occurs
     */
    ctx_t context_of(const inp_t& P, uint64_t beg, uint64_t end, bool& occurs) const
    {
        if constexpr (has_context_tables) {
            return _index.template search<query_support>(
                std::span<const sym_t>(P.data() + beg, end - beg), occurs);
        } else {
            ctx_t ctx = _index.template empty_context<query_support>();
            occurs = true;

            for (uint64_t i = end; i > beg; i--) {
                auto [next, ok] = ctx.extend(P[i - 1], LEFT);
                if (!ok) { occurs = false; return ctx; }
                ctx = std::move(next);
            }

            return ctx;
        }
    }

    // returns the number of occurrences of part i, or 0 if it does not occur (only if _has_context[i])
    uint64_t num_occ_of_part(uint64_t i) const { return _occurs[i] ? _ctx[i].num_occ() : 0; }

    // computes the search context of each part of P from _bounds
    void compute_contexts(const inp_t& P)
    {
        for (uint64_t i = 0; i + 1 < _bounds.size(); i++) {
            _has_context[i] = has_no_free_symbol(_bounds[i], _bounds[i + 1]);
            if (!_has_context[i]) continue;
            bool occurs;
            _ctx[i] = context_of(P, _bounds[i], _bounds[i + 1], occurs);
            _occurs[i] = occurs;
        }
    }

    /**
     * @brief returns the length of the initial parts of the dynamic partitioning: the maximum pattern length of the
     *        context tables of the index if it has them, else 1
     * @param m the pattern length
     * @param p the number of parts
     * @return the length, in [1, max(1, m / (2p))]
     */
    uint64_t initial_part_length(uint64_t m, uint64_t p) const
    {
        uint64_t len = 1;

        if constexpr (requires(const idx_t& idx) { idx.has_context_tables(); idx.context_tables_max_length(); }) {
            if (_index.has_context_tables()) len = std::max<uint64_t>(1, _index.context_tables_max_length());
        }

        // the p initial parts cover at most half of the pattern
        return std::clamp<uint64_t>(len, 1, std::max<uint64_t>(1, m / (2 * p)));
    }

    /**
     * @brief partitions the pattern into p parts of equal length and computes the search context of each part
     * @param P the pattern
     * @param p the number of parts
     */
    void partition_uniform(const inp_t& P, uint64_t p)
    {
        const uint64_t m = P.size();
        _bounds[0] = 0;
        _bounds[p] = m;

        for (uint64_t i = 1; i < p; i++) {
            // each part contains at least one symbol
            _bounds[i] = std::clamp<uint64_t>((i * m) / p, _bounds[i - 1] + 1, m - (p - i));
        }

        compute_contexts(P);
    }

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
    void partition_dynamic(const inp_t& P, uint64_t p)
    {
        const uint64_t m = P.size();
        const uint64_t init_len = initial_part_length(m, p);
        std::vector<uint64_t> beg(p), end(p); // [0..p-1] part i is P[beg[i], end[i])

        // place the initial parts
        beg[0] = 0;
        beg[p - 1] = m - init_len;

        for (uint64_t i = 1; i + 1 < p; i++) {
            uint64_t centre = (i * m) / (p - 1);
            uint64_t start = centre > init_len / 2 ? centre - init_len / 2 : 0; // the initial part is centred at centre
            beg[i] = std::clamp<uint64_t>(start, beg[i - 1] + init_len, m - init_len * (p - i));
        }

        for (uint64_t i = 0; i < p; i++) {
            end[i] = beg[i] + init_len;
            _has_context[i] = has_no_free_symbol(beg[i], end[i]);
            if (!_has_context[i]) continue;
            bool occurs;
            _ctx[i] = context_of(P, beg[i], end[i], occurs);
            _occurs[i] = occurs;
        }

        // the weight of the numbers of occurrences of the two outer parts in the dynamic partitioning of Columba
        constexpr uint64_t outer_part_weight = 2;
        // [0..p-1] the weight of the number of occurrences of part i: outer_part_weight for the two outer parts, else 1
        std::vector<uint64_t> weight(p, 1);
        weight.front() = outer_part_weight;
        weight.back() = outer_part_weight;

        for (uint64_t assigned = p * init_len; assigned < m; assigned++) {
            // the part to extend (the one with the largest weighted number of occurrences), that number and the direction
            uint64_t max_part = 0, max_weighted_occ = 0;
            direction_t dir = RIGHT;

            for (uint64_t i = 0; i < p; i++) {
                // whether part i cannot be extended to the left (no_left) and to the right (no_right)
                bool no_left = i == 0 || beg[i] == end[i - 1] || is_free(beg[i] - 1);
                bool no_right = i + 1 == p || end[i] == beg[i + 1] || is_free(end[i]);
                if ((no_left && no_right) || !_has_context[i]) continue;

                uint64_t weighted_occ = num_occ_of_part(i) * weight[i];
                if (weighted_occ <= max_weighted_occ) continue;

                max_weighted_occ = weighted_occ;
                max_part = i;

                if (no_left) dir = RIGHT;
                else if (no_right) dir = LEFT;
                else dir = num_occ_of_part(i - 1) < num_occ_of_part(i + 1) ? LEFT : RIGHT;
            }

            // no part that can be extended occurs: assign the remaining symbols to the parts on their left
            if (max_weighted_occ == 0) {
                for (uint64_t i = 0; i + 1 < p; i++) {
                    for (; end[i] < beg[i + 1]; end[i]++) {
                        if (!_has_context[i] || !_occurs[i]) continue;
                        if (is_free(end[i])) { _has_context[i] = false; continue; }
                        auto [next, ok] = _ctx[i].extend(P[end[i]], RIGHT);
                        if (ok) _ctx[i] = std::move(next);
                        else _occurs[i] = false;
                    }
                }

                break;
            }

            uint64_t pos = dir == RIGHT ? end[max_part]++ : --beg[max_part];
            auto [next, ok] = _ctx[max_part].extend(P[pos], dir);
            if (ok) _ctx[max_part] = std::move(next);
            else _occurs[max_part] = false;
        }

        for (uint64_t i = 0; i < p; i++) _bounds[i] = beg[i];
        _bounds[p] = m;
    }

    // partitions P into p parts of equal length, or with partition_dynamic unless params.partitioning is uniform
    void make_parts(const inp_t& P, uint64_t p, const apm_params_t& params)
    {
        // the dynamic partitioning requires p >= 2 and m >= 2p
        if (p == 1 || params.partitioning == _partition_uniform || P.size() < 2 * p) {
            partition_uniform(P, p);
        } else {
            partition_dynamic(P, p);
        }
    }

    /**
     * @brief returns the bounds of the p parts of P that minimize the sum of predictor.part_steps over the parts, by
     *        dynamic programming over the bounds
     * @param m the pattern length
     * @param p the number of parts
     * @param predictor the prediction of the number of steps of the searches
     * @param occ_of occ_of(b, e) returns the number of occurrences of P[b, e)
     */
    template <typename occ_fnc_t>
    static std::vector<uint64_t> minimal_bounds(uint64_t m, uint64_t p, const apm_step_prediction& predictor, occ_fnc_t occ_of)
    {
        /* [0..p][0..m] f[i][e] = the smallest sum of part_steps over parts 0, ..., i - 1 that cover P[0, e) (infinity if
           there are none), and from[i][e] = the start of part i - 1 in these parts */
        const double inf = std::numeric_limits<double>::infinity();
        std::vector<std::vector<double>> f(p + 1, std::vector<double>(m + 1, inf));
        std::vector<std::vector<uint64_t>> from(p + 1, std::vector<uint64_t>(m + 1, 0));
        f[0][0] = 0;

        for (uint64_t i = 0; i < p; i++) {
            for (uint64_t b = i; b + (p - i) <= m; b++) {
                if (f[i][b] == inf) continue;

                // part i is P[b, e)
                for (uint64_t e = b + 1; e + (p - 1 - i) <= m; e++) {
                    const double sum = f[i][b] + predictor.part_steps(i, occ_of(b, e));
                    if (sum < f[i + 1][e]) { f[i + 1][e] = sum; from[i + 1][e] = b; }
                }
            }
        }

        std::vector<uint64_t> bounds(p + 1);
        bounds[p] = m;
        for (uint64_t i = p; i > 0; i--) bounds[i - 1] = from[i][bounds[i]];
        return bounds;
    }

    /**
     * @brief partitions P into the parts of a search scheme or of its mirror (then _mirrored is set) such that the
     *        predicted number of steps of the searches (see apm_step_prediction) is small, and computes the search context
     *        of each part; the optimization counts its own steps (the context extensions of step 2, the cells of
     *        minimal_bounds and the steps of the predictions, see apm_step_prediction::work) and stops once they exceed
     *        the smallest prediction found so far. A part is rare if it has at most one occurrence with a distinct
     *        continuation (occ * r / n <= 1; in step 1, a part without a search context is not rare).
     *
     * 1. P is partitioned with partition_dynamic. These parts are kept if every part is rare. They are also kept if one
     *    part is rare or if the smaller prediction of the scheme and of the mirror is at most the number of context
     *    extensions of step 2 (see occurrences_steps); then the direction with the smaller prediction is chosen.
     * 2. The numbers of occurrences of the substrings of P are computed (see occurrences). The number of occurrences
     *    of a part is estimated by the smallest number of occurrences of its maximal substrings without a free symbol;
     *    a substring longer than its row of the table is counted as its longest prefix in the row. From here on, the
     *    predictions count the occurrences near the pattern (see apm_step_prediction::set_occurrences).
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
        const apm_params_t& params)
    {
        const uint64_t m = P.size();
        const uint64_t p = scheme.p;
        const uint64_t r = _index.num_bwt_runs();
        const uint64_t n = _index.input_size();
        const search_scheme_t mirror = mirrored_scheme(scheme);

        // how often the anchor of the error model occurs (see apm_costs_t::anchor)
        uint64_t anchor_occ = 0;

        if (_costs != nullptr && params.anchored_search) {
            auto [a_beg, a_end] = _costs->anchor();
            if (a_beg < a_end) anchor_occ = num_occ_of(P, a_beg, a_end);
        }

        const std::vector<double>& sym_p = symbol_probabilities();
        // [0] predicts the steps of the searches of the scheme, [1] those of the searches of its mirror
        apm_step_prediction predictor[2] {
            { scheme, _costs, metric, P, _index.alphabet_size(), sym_p, r, n, anchor_occ, params },
            { mirror, _costs, metric, P, _index.alphabet_size(), sym_p, r, n, anchor_occ, params }
        };

        // 1. the parts of partition_dynamic
        partition_dynamic(P, p);
        const std::vector<uint64_t> bounds_dynamic = _bounds;

        // the parts are kept if every part is rare
        bool all_rare = true;
        for (uint64_t i = 0; i < p; i++) all_rare &= _has_context[i] && double(num_occ(i)) * double(r) <= double(n);
        if (all_rare) return;

        std::vector<uint64_t> occ(p); // [0..p-1] the number of occurrences of part i
        bool has_rare = false;

        for (uint64_t i = 0; i < p; i++) {
            occ[i] = num_occ(i);
            has_rare |= _has_context[i] && double(occ[i]) * double(r) <= double(n);
        }

        // with a rare part, or if the smaller prediction is at most the context extensions of step 2, only the direction is chosen
        const double steps_scheme = predictor[0].steps(_bounds, occ), steps_mirror = predictor[1].steps(_bounds, occ);
        double best = std::min(steps_scheme, steps_mirror); // the smallest prediction so far
        _mirrored = steps_mirror < steps_scheme;
        if (has_rare || best <= occurrences_steps(m)) return;

        // 2. the numbers of occurrences of the substrings of P
        const occurrence_table_t& table = occurrences(P);
        for (apm_step_prediction& c : predictor) c.set_occurrences(&table);
        double table_steps = double(_occ_steps); // the steps of the optimization other than the predictions
        double steps = table_steps + double(predictor[0].work() + predictor[1].work()); // all steps of the optimization

        // [0..m] next_free[i] = the first position >= i with a free symbol, or m
        std::vector<uint64_t> next_free(m + 1, m);
        for (uint64_t i = m; i-- > 0;) next_free[i] = is_free(i) ? i : next_free[i + 1];

        /* returns the estimated number of occurrences of P[b, e): the smallest number of occurrences of its maximal
           substrings without a free symbol, since P[b, e) occurs at most as often as each of them, or n if every
           position of P[b, e) has a free symbol */
        auto occ_of = [&](uint64_t b, uint64_t e) {
            uint64_t smallest = n;

            // P[sub_beg, sub_end) is a maximal substring without a free symbol
            for (uint64_t sub_beg = b; sub_beg < e;) {
                if (is_free(sub_beg)) { sub_beg++; continue; }
                const uint64_t sub_end = std::min(next_free[sub_beg], e);
                const std::vector<uint64_t>& row = table[sub_beg];
                smallest = std::min(smallest, row[std::min<uint64_t>(sub_end - sub_beg, row.size()) - 1]);
                sub_beg = sub_end;
            }

            return smallest;
        };

        // returns whether P[b, e) is rare
        auto is_rare = [&](uint64_t b, uint64_t e) { return double(occ_of(b, e)) * double(r) <= double(n); };

        /* returns the prediction of predictor[mirrored] for the parts with the bounds, and sets steps to the steps of the
           optimization so far */
        auto predict = [&](uint64_t mirrored, const std::vector<uint64_t>& bounds) {
            for (uint64_t i = 0; i < p; i++) occ[i] = occ_of(bounds[i], bounds[i + 1]);
            const double prediction = predictor[mirrored].steps(bounds, occ);
            steps = table_steps + double(predictor[0].work() + predictor[1].work());
            return prediction;
        };

        best = predict(0, bounds_dynamic);
        const double best_mirror = predict(1, bounds_dynamic);
        _mirrored = best_mirror < best;
        best = std::min(best, best_mirror);

        // 3. the local searches
        std::vector<uint64_t> bounds_uniform(p + 1);
        for (uint64_t i = 0; i <= p; i++) bounds_uniform[i] = (i * m) / p;
        std::vector<uint64_t> bounds_best = bounds_dynamic;
        bool mirrored_best = _mirrored;

        for (uint64_t mirrored = 0; mirrored < 2; mirrored++) {
            std::vector<std::vector<uint64_t>> starts { bounds_dynamic, bounds_uniform };
            const double cells_dp = double(p) * double(m) * double(m); // the cells of minimal_bounds

            if (steps + cells_dp <= best) {
                starts.emplace_back(minimal_bounds(m, p, predictor[mirrored], occ_of));
                table_steps += cells_dp;
                steps = table_steps + double(predictor[0].work() + predictor[1].work());
            }

            for (std::vector<uint64_t> bounds : starts) {
                double prediction = predict(mirrored, bounds);

                while (steps <= best) {
                    const double steps_before = steps;
                    // the move with the smallest prediction: bound i_move by one position to the right if move_right
                    double prediction_move = prediction;
                    uint64_t i_move = 0;
                    bool move_right = false;

                    for (uint64_t i = 1; i < p; i++) {
                        // a bound next to a rare part is not moved
                        if (is_rare(bounds[i - 1], bounds[i]) || is_rare(bounds[i], bounds[i + 1])) continue;

                        for (bool right : { false, true }) {
                            if (right ? bounds[i] + 1 == bounds[i + 1] : bounds[i] - 1 == bounds[i - 1]) continue;
                            bounds[i] = right ? bounds[i] + 1 : bounds[i] - 1;
                            const double prediction_i = predict(mirrored, bounds);
                            bounds[i] = right ? bounds[i] - 1 : bounds[i] + 1;
                            if (prediction_i < prediction_move) { prediction_move = prediction_i; i_move = i; move_right = right; }
                        }
                    }

                    if (i_move == 0) break;
                    bounds[i_move] = move_right ? bounds[i_move] + 1 : bounds[i_move] - 1;
                    const double gain = prediction - prediction_move;
                    prediction = prediction_move;
                    if (gain < steps - steps_before) break;
                }

                if (prediction < best) { best = prediction; bounds_best = bounds; mirrored_best = mirrored == 1; }
            }
        }

        _mirrored = mirrored_best;
        if (bounds_best == bounds_dynamic) return;
        _bounds = bounds_best;
        compute_contexts(P);
    }

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
        for (uint64_t i = 0; i < p; i++) _ctx.emplace_back(index.template empty_context<query_support>());
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
        const apm_params_t& params)
        : _index(index), _bounds(scheme.p + 1), _occurs(scheme.p, 0), _has_context(scheme.p, 1), _costs(params.costs)
    {
        const uint64_t p = scheme.p;
        assert(0 < p && p <= P.size());
        _ctx.reserve(p);
        for (uint64_t i = 0; i < p; i++) _ctx.emplace_back(index.template empty_context<query_support>());

        if constexpr (can_optimize) {
            if (params.partitioning == _partition_optimized && p >= 2 && P.size() >= 2 * p) {
                partition_optimized(P, scheme, metric, params);
                return;
            }
        }

        make_parts(P, p, params);
    }

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
        for (uint64_t i = 0; i + 1 < _bounds.size(); i++) _ctx.emplace_back(index.template empty_context<query_support>());
        compute_contexts(P);
    }

    /**
     * @brief returns [0..255] the probability of each byte in the text (how often it occurs, divided by n) for an index
     *        with a byte alphabet, else an empty vector (computed by the first call)
     */
    const std::vector<double>& symbol_probabilities() const
    {
        if constexpr (requires(const idx_t& idx) { idx.forward_index().map_int(); idx.input_size(); }) {
            const std::vector<uint8_t>& map_int = _index.forward_index().map_int();
            if (!_sym_p.empty() || map_int.empty()) return _sym_p;
            _sym_p.resize(map_int.size(), 0);

            for (uint64_t c = 0; c < map_int.size(); c++) {
                if (map_int[c] == 0) continue;
                auto [ctx, occurs] = ctx_t(_index).extend(sym_t(c), RIGHT);
                if (occurs) _sym_p[c] = double(ctx.num_occ()) / double(_index.input_size());
            }
        }

        return _sym_p;
    }

    /**
     * @brief returns the numbers of occurrences of the substrings of P that start at each position (computed by the
     *        first call): [i][l - 1] is the number of occurrences of P[i, i + l), for l = 1, 2, ... until at most one
     *        occurrence has a distinct continuation (occ * r / n <= 1, see apm_step_prediction), the substring does not
     *        occur or P ends; at a position with free symbols, the occurrences of all of them are counted
     * @param P the pattern
     */
    const occurrence_table_t& occurrences(const inp_t& P) const
    {
        if (!_occ_table.empty()) return _occ_table;
        const uint64_t m = P.size();
        const double runs_per_symbol = double(_index.num_bwt_runs()) / double(_index.input_size()); // r / n
        _occ_table.resize(m);
        std::vector<ctx_t> ctxs, next; // the contexts of the strings that match P[i, e), and a buffer

        for (uint64_t i = 0; i < m; i++) {
            ctxs.assign(1, ctx_t(_index));

            for (uint64_t e = i; e < m && !ctxs.empty(); e++) {
                uint64_t occ = 0;
                next.clear();

                auto extend = [&](ctx_t& ctx, sym_t c) {
                    auto [ext, ok] = ctx.extend(c, RIGHT);
                    _occ_steps++;
                    if (!ok) return;
                    occ += ext.num_occ();
                    next.emplace_back(std::move(ext));
                };

                for (ctx_t& ctx : ctxs) {
                    if (is_free(e)) for (uint8_t c : _costs->free_symbols(e)) extend(ctx, sym_t(c));
                    else extend(ctx, P[e]);
                }

                _occ_table[i].emplace_back(occ);
                std::swap(ctxs, next);
                if (double(occ) * runs_per_symbol <= 1) break;
            }
        }

        return _occ_table;
    }

    // returns whether occurrences has been computed
    bool has_occurrences() const { return !_occ_table.empty(); }

    /* returns an estimate of the number of context extensions of occurrences: m * log_sigma(r) with the effective
       alphabet size sigma (see apm_step_prediction), since a row of the table ends once occ * r / n <= 1, after about
       log_sigma(r) symbols */
    double occurrences_steps(uint64_t m) const
    {
        double sigma = double(_index.alphabet_size());
        double collision = 0; // sum_c p_c^2, the probability that two random symbols of the text are equal
        for (double p : symbol_probabilities()) collision += p * p;
        if (collision > 0) sigma = 1 / collision;
        if (sigma <= 1) return std::numeric_limits<double>::infinity();
        return double(m) * std::log(double(_index.num_bwt_runs())) / std::log(sigma);
    }

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

#endif // MOVE_R_ALGORITHMS_MOVE_RB_APM_PARTITION_HPP
