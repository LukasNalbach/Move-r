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
#include <algorithms/move_rb/apm/search_cost.hpp>

/**
 * @brief a partition of a pattern P into p parts, with the search context of each part: part i is
 *        P[beg(i), end(i)). The search context of a part is computed if no position of the part has a free symbol
 *        other than its own symbol (see apm_costs_t::has_free; always without an error model). The parts are of
 *        equal length, chosen by the dynamic partitioning of Columba (Renders, Depuydt, Fostier), or chosen together
 *        with the direction of the search scheme such that the predicted number of steps of the searches is small
 *        (see apm_search_cost).
 * @tparam idx_t the bidirectional index type
 * @tparam query_support whether the contexts support locate (LOCATE) or only count (COUNT)
 */
template <typename idx_t, query_support_t query_support>
class apm_partition {
public:
    using inp_t = typename idx_t::inp_t;
    using sym_t = std::ranges::range_value_t<inp_t>;
    using ctx_t = typename idx_t::template search_context_t<query_support>;

protected:
    const idx_t& _index;

    std::vector<uint64_t> _bounds; // [0..p] part i is P[_bounds[i], _bounds[i + 1]); _bounds[0] = 0, _bounds[p] = m
    std::vector<ctx_t> _ctx;       // [0..p-1] the search context of part i (valid if _has_context[i] and _occurs[i])
    std::vector<uint8_t> _occurs;  // [0..p-1] 1 <=> part i occurs in the text (valid if _has_context[i])
    // [0..p-1] 1 <=> no position of part i has a free symbol other than its own, so _ctx[i] is the context of part i
    std::vector<uint8_t> _has_context;
    const apm_costs_t* _costs;     // the error model, or nullptr
    bool _mirrored = false;        // whether the searches run with the mirror of the search scheme

    // returns whether P[i] has a free symbol other than itself
    bool is_free(uint64_t i) const { return _costs != nullptr && _costs->has_free(i); }

    // returns whether no position of P[beg, end) has a free symbol other than itself
    bool has_no_free_symbol(uint64_t beg, uint64_t end) const
    {
        for (uint64_t i = beg; i < end; i++) if (is_free(i)) return false;
        return true;
    }

    // whether the index can look up the search context of a pattern in its context tables
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

    // computes the search context of each part from _bounds
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
     * 2. The remaining symbols are assigned one by one: the part with the largest number of occurrences (that of
     *    the two outer parts weighted by 2) is extended by one symbol, towards the neighbour with fewer occurrences
     *    if it can be extended in both directions. A part is not extended over a position with a free symbol other
     *    than itself, nor once its context is not computed.
     * 3. If no part that can be extended occurs, the remaining symbols are assigned to the part on their left.
     *
     * @param P the pattern
     * @param p the number of parts
     */
    void partition_dynamic(const inp_t& P, uint64_t p)
    {
        const uint64_t m = P.size();
        const uint64_t w = initial_part_length(m, p);
        std::vector<uint64_t> beg(p), end(p);

        // place the initial parts
        beg[0] = 0;
        beg[p - 1] = m - w;

        for (uint64_t i = 1; i + 1 < p; i++) {
            uint64_t centre = (i * m) / (p - 1);
            uint64_t start = centre > w / 2 ? centre - w / 2 : 0; // the initial part is centred at centre
            beg[i] = std::clamp<uint64_t>(start, beg[i - 1] + w, m - w * (p - i));
        }

        for (uint64_t i = 0; i < p; i++) {
            end[i] = beg[i] + w;
            _has_context[i] = has_no_free_symbol(beg[i], end[i]);
            if (!_has_context[i]) continue;
            bool occurs;
            _ctx[i] = context_of(P, beg[i], end[i], occurs);
            _occurs[i] = occurs;
        }

        std::vector<uint64_t> weight(p, 1);
        weight.front() = 2;
        weight.back() = 2;

        for (uint64_t assigned = p * w; assigned < m; assigned++) {
            // the part to extend (the one with the largest weighted number of occurrences) and that number
            uint64_t max_part = 0, max_val = 0;
            direction_t dir = RIGHT;

            for (uint64_t i = 0; i < p; i++) {
                bool no_left = i == 0 || beg[i] == end[i - 1] || is_free(beg[i] - 1);
                bool no_right = i + 1 == p || end[i] == beg[i + 1] || is_free(end[i]);
                if ((no_left && no_right) || !_has_context[i]) continue;

                uint64_t val = num_occ_of_part(i) * weight[i];
                if (val <= max_val) continue;

                max_val = val;
                max_part = i;

                if (no_left) dir = RIGHT;
                else if (no_right) dir = LEFT;
                else dir = num_occ_of_part(i - 1) < num_occ_of_part(i + 1) ? LEFT : RIGHT;
            }

            // no part that can be extended occurs: assign the remaining symbols to the parts on their left
            if (max_val == 0) {
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

    // whether the index provides its alphabet size, the number of runs in its BWT and the length of its text
    static constexpr bool can_optimize = requires(const idx_t& idx) {
        idx.alphabet_size();
        idx.num_bwt_runs();
        idx.input_size();
    };

    /**
     * @brief returns the numbers of occurrences of the substrings of P that start at each position: [i][l - 1] is the
     *        number of occurrences of P[i, i + l), for l = 1, 2, ... until at most one occurrence has a distinct
     *        continuation (occ * share <= 1, see apm_search_cost), the substring does not occur, it has a position with
     *        a free symbol or P ends; a substring with a free symbol has no search context (a search that starts with it
     *        starts at the root), it counts as n occurrences
     * @param P the pattern
     * @param share r / n
     * @param steps increased by the number of extensions
     */
    std::vector<std::vector<uint64_t>> occurrence_table(const inp_t& P, double share, double& steps) const
    {
        const uint64_t m = P.size();
        std::vector<std::vector<uint64_t>> occ(m);

        for (uint64_t i = 0; i < m; i++) {
            ctx_t ctx = _index.template empty_context<query_support>();

            for (uint64_t e = i; e < m; e++) {
                if (is_free(e)) { occ[i].emplace_back(_index.input_size()); break; }
                auto [next, ok] = ctx.extend(P[e], RIGHT);
                steps++;
                if (!ok) { occ[i].emplace_back(0); break; }
                ctx = std::move(next);
                occ[i].emplace_back(ctx.num_occ());
                if (double(ctx.num_occ()) * share <= 1) break;
            }
        }

        return occ;
    }

    /**
     * @brief returns the bounds of the p parts of P that minimize the sum of cost.part_cost (dynamic programming over
     *        the bounds: f[i][e] is the smallest sum for parts 0..i-1 covering P[0, e))
     * @param m the pattern length
     * @param p the number of parts
     * @param cost the prediction
     * @param occ_of occ_of(b, e) returns the number of occurrences of P[b, e)
     */
    template <typename occ_fnc_t>
    static std::vector<uint64_t> minimal_bounds(uint64_t m, uint64_t p, const apm_search_cost& cost, occ_fnc_t occ_of)
    {
        const double inf = std::numeric_limits<double>::infinity();
        std::vector<std::vector<double>> f(p + 1, std::vector<double>(m + 1, inf));
        std::vector<std::vector<uint64_t>> from(p + 1, std::vector<uint64_t>(m + 1, 0));
        f[0][0] = 0;

        for (uint64_t i = 0; i < p; i++) {
            for (uint64_t b = i; b + (p - i) <= m; b++) {
                if (f[i][b] == inf) continue;

                for (uint64_t e = b + 1; e + (p - 1 - i) <= m; e++) {
                    const double c = f[i][b] + cost.part_cost(i, occ_of(b, e));
                    if (c < f[i + 1][e]) { f[i + 1][e] = c; from[i + 1][e] = b; }
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
     *        predicted number of steps of the searches (see apm_search_cost) is small, and computes the search context
     *        of each part; the optimization counts its own steps (an extension of a context, a cell or a transition of
     *        a dynamic program) and stops once they exceed the prediction for the best partition found:
     *
     * 1. The partition of partition_dynamic; it stays if each of its parts has at most one occurrence with a distinct
     *    continuation (occ * r / n <= 1, a rare part) or if apm_model_plan runs the anchored search for every partition;
     *    it also stays, with the scheme or the mirror, whichever has the smaller prediction, if one of its parts is rare
     *    or if the prediction is below the m * log_sigma(r) extensions of step 2.
     * 2. The numbers of occurrences of the substrings of P (see occurrence_table); a longer substring is assumed to
     *    occur as often as its longest prefix in the table.
     * 3. For the scheme and for the mirror: a local search from the uniform partition, from the partition of step 1
     *    and from minimal_bounds (if its cells do not exceed the budget), which moves one bound by one symbol, the move
     *    that lowers the prediction most, until no move lowers it by more than the steps of trying the moves; a bound
     *    next to a rare part does not move.
     *
     * A position has one symbol without an error, sigma - 1 substitutions and, under the edit distance, an insertion and
     * sigma deletions; with an error model, its free symbols and the errors whose costs are within its bound.
     *
     * @param P the pattern
     * @param scheme the search scheme
     * @param metric the distance metric of the searches
     */
    void partition_optimized(const inp_t& P, const search_scheme_t& scheme, distance_metric_t metric)
    {
        const uint64_t m = P.size();
        const uint64_t p = scheme.p;
        const double sigma = double(_index.alphabet_size());
        const uint64_t r = _index.num_bwt_runs();
        const uint64_t n = _index.input_size();
        const search_scheme_t mirror = mirrored_scheme(scheme);
        const bool edit = metric == EDIT_DISTANCE;
        std::vector<double> match(m, 1), sub(m, sigma - 1), ins(m, edit ? 1 : 0), del(m, edit ? sigma : 0);

        if (_costs != nullptr) {
            for (uint64_t i = 0; i < m; i++) {
                match[i] = double(_costs->free_symbols(i).size());
                sub[i] = 0;

                for (uint64_t c = 0; c < _costs->num_ranks(); c++) {
                    if (!_costs->is_free(i, c) && _costs->sub(i, c) <= _costs->bound()) sub[i]++;
                }

                ins[i] = edit && _costs->ins(i) <= _costs->bound() ? 1 : 0;
                del[i] = edit && _costs->del(i) <= _costs->bound() ? sigma : 0;
            }
        }

        // the occurrences of the anchor of the error model (see apm_costs_t::anchor)
        uint64_t anchor_occ = 0;

        if (_costs != nullptr) {
            auto [a_beg, a_end] = _costs->anchor();

            if (a_beg < a_end) {
                bool occurs;
                const ctx_t ctx = context_of(P, a_beg, a_end, occurs);
                if (occurs) anchor_occ = ctx.num_occ();
            }
        }

        const apm_search_cost cost[2] {
            { scheme, _costs, sigma, match, sub, ins, del, r, n, anchor_occ },
            { mirror, _costs, sigma, match, sub, ins, del, r, n, anchor_occ }
        };

        partition_dynamic(P, p);
        const std::vector<uint64_t> bounds_dynamic = _bounds;

        // if every part has at most one occurrence with a distinct continuation, every search starts with at most one context
        bool all_rare = true;
        for (uint64_t i = 0; i < p; i++) all_rare &= _has_context[i] && double(num_occ(i)) * double(r) <= double(n);
        if (all_rare) return;

        /* if a position has a free symbol and every part is the first part of a search, every partition has a search that
           starts without a search context, so apm_model_plan runs the anchored search if there is an anchor, whatever the
           parts are; then the dynamic partition stays */
        if (_costs != nullptr && _costs->anchor().first < _costs->anchor().second) {
            bool has_free_symbol = false;
            for (uint64_t i = 0; i < m; i++) has_free_symbol |= is_free(i);
            std::vector<uint8_t> is_first(p, 0);
            for (const search_t& S : scheme.S) is_first[S[0].part] = 1;
            if (has_free_symbol && std::find(is_first.begin(), is_first.end(), 0) == is_first.end()) return;
        }

        std::vector<uint64_t> occ(p);
        bool has_rare = false;

        for (uint64_t i = 0; i < p; i++) {
            occ[i] = num_occ(i);
            has_rare |= _has_context[i] && double(occ[i]) * double(r) <= double(n);
        }

        // with a rare part or a prediction below the steps of the table, only the direction is chosen
        const double cost_scheme = cost[0].cost(_bounds, occ), cost_mirror = cost[1].cost(_bounds, occ);
        double best = std::min(cost_scheme, cost_mirror);
        _mirrored = cost_mirror < cost_scheme;
        if (has_rare || best <= double(m) * std::log(double(r)) / std::log(sigma)) return;

        double table_steps = 0;
        const std::vector<std::vector<uint64_t>> table = occurrence_table(P, double(r) / double(n), table_steps);
        double steps = table_steps + double(cost[0].work() + cost[1].work());

        // next_free[i] = the first position >= i with a free symbol, or m
        std::vector<uint64_t> next_free(m + 1, m);
        for (uint64_t i = m; i-- > 0;) next_free[i] = is_free(i) ? i : next_free[i + 1];

        /* returns the number of occurrences of P[b, e): at most that of each of its substrings without a free symbol, so
           the smallest of those of its maximal ones (a substring longer than its row of the table occurs as often as the
           longest prefix in the row), or n if every position of it has a free symbol */
        auto occ_of = [&](uint64_t b, uint64_t e) {
            uint64_t smallest = n;

            for (uint64_t x = b; x < e;) {
                if (is_free(x)) { x++; continue; }
                const uint64_t y = std::min(next_free[x], e);
                const std::vector<uint64_t>& row = table[x];
                smallest = std::min(smallest, row[std::min<uint64_t>(y - x, row.size()) - 1]);
                x = y;
            }

            return smallest;
        };

        // returns whether P[b, e) has at most one occurrence with a distinct continuation
        auto is_rare = [&](uint64_t b, uint64_t e) { return double(occ_of(b, e)) * double(r) <= double(n); };

        // returns the prediction of cost[o] for the parts with the bounds; steps counts the work of all predictions
        auto predict = [&](uint64_t o, const std::vector<uint64_t>& bounds) {
            for (uint64_t i = 0; i < p; i++) occ[i] = occ_of(bounds[i], bounds[i + 1]);
            const double c = cost[o].cost(bounds, occ);
            steps = table_steps + double(cost[0].work() + cost[1].work());
            return c;
        };

        best = predict(0, bounds_dynamic);
        const double best_mirror = predict(1, bounds_dynamic);
        _mirrored = best_mirror < best;
        best = std::min(best, best_mirror);

        std::vector<uint64_t> bounds_uniform(p + 1);
        for (uint64_t i = 0; i <= p; i++) bounds_uniform[i] = (i * m) / p;
        std::vector<uint64_t> bounds_best = bounds_dynamic;
        bool mirrored_best = _mirrored;

        for (uint64_t o = 0; o < 2; o++) {
            std::vector<std::vector<uint64_t>> starts { bounds_dynamic, bounds_uniform };
            const double cells_dp = double(p) * double(m) * double(m);

            if (steps + cells_dp <= best) {
                starts.emplace_back(minimal_bounds(m, p, cost[o], occ_of));
                table_steps += cells_dp;
                steps = table_steps + double(cost[0].work() + cost[1].work());
            }

            for (std::vector<uint64_t> bounds : starts) {
                double c = predict(o, bounds);

                while (steps <= best) {
                    const double steps_before = steps;
                    double c_move = c;
                    uint64_t i_move = 0;
                    bool move_right = false;

                    for (uint64_t i = 1; i < p; i++) {
                        // a border between parts with at most one occurrence with a distinct continuation stays
                        if (is_rare(bounds[i - 1], bounds[i]) || is_rare(bounds[i], bounds[i + 1])) continue;

                        for (bool right : { false, true }) {
                            if (right ? bounds[i] + 1 == bounds[i + 1] : bounds[i] - 1 == bounds[i - 1]) continue;
                            bounds[i] = right ? bounds[i] + 1 : bounds[i] - 1;
                            const double c_i = predict(o, bounds);
                            bounds[i] = right ? bounds[i] - 1 : bounds[i] + 1;
                            if (c_i < c_move) { c_move = c_i; i_move = i; move_right = right; }
                        }
                    }

                    if (i_move == 0) break;
                    bounds[i_move] = move_right ? bounds[i_move] + 1 : bounds[i_move] - 1;
                    const double gain = c - c_move;
                    c = c_move;
                    if (gain < steps - steps_before) break;
                }

                if (c < best) { best = c; bounds_best = bounds; mirrored_best = o == 1; }
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
     * @param params query parameters (partitioning strategy, error model)
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
     *        _partition_optimized (if the index provides what it needs), the searches run with the mirror of the scheme
     *        if mirrored(), else as with the constructor that takes p
     * @param index the index to search in
     * @param P the pattern
     * @param scheme the search scheme (0 < scheme.p <= P.size())
     * @param metric the distance metric of the searches
     * @param params query parameters (partitioning strategy, error model)
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
                partition_optimized(P, scheme, metric);
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
     * @param params query parameters (the error model)
     */
    apm_partition(const idx_t& index, const inp_t& P, std::vector<uint64_t> bounds, const apm_params_t& params)
        : _index(index), _bounds(std::move(bounds)), _occurs(_bounds.size() - 1, 0),
          _has_context(_bounds.size() - 1, 1), _costs(params.costs)
    {
        for (uint64_t i = 0; i + 1 < _bounds.size(); i++) _ctx.emplace_back(index.template empty_context<query_support>());
        compute_contexts(P);
    }

    // returns the number of parts
    uint64_t size() const { return _bounds.size() - 1; }

    // returns whether the searches run with the mirror of the search scheme (see mirrored_scheme)
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
