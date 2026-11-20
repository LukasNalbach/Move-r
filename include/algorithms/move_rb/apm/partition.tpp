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
 * @file partition.tpp
 * @brief the methods of apm_partition: the search contexts and the uniform, dynamic and optimized partitions
 */

#ifndef MOVE_R_ALGORITHMS_MOVE_RB_APM_PARTITION_TPP
#define MOVE_R_ALGORITHMS_MOVE_RB_APM_PARTITION_TPP

#include <algorithms/move_rb/apm/partition.hpp>

template <typename idx_t, query_support_t query_support>
inline auto apm_partition<idx_t, query_support>::context_of(const inp_t& P, uint64_t beg, uint64_t end, bool& occurs) const -> ctx_t
{
    if constexpr (has_context_tables) {
        return _index.template apm_context_of<query_support>(
            std::span<const sym_t>(P.data() + beg, end - beg), occurs);
    } else {
        ctx_t ctx(_index);
        occurs = true;

        for (uint64_t i = end; i > beg; i--) {
            auto [next, ok] = ctx.extend(P[i - 1], LEFT);
            if (!ok) { occurs = false; return ctx; }
            ctx = std::move(next);
        }

        return ctx;
    }
}

template <typename idx_t, query_support_t query_support>
inline void apm_partition<idx_t, query_support>::compute_contexts(const inp_t& P)
{
    for (uint64_t i = 0; i + 1 < _bounds.size(); i++) {
        _has_context[i] = has_no_free_symbol(_bounds[i], _bounds[i + 1]);
        if (!_has_context[i]) continue;
        bool occurs;
        _ctx[i] = context_of(P, _bounds[i], _bounds[i + 1], occurs);
        _occurs[i] = occurs;
    }
}

template <typename idx_t, query_support_t query_support>
inline uint64_t apm_partition<idx_t, query_support>::initial_part_length(uint64_t m, uint64_t p) const
{
    uint64_t len = 1;

    if constexpr (requires(const idx_t& idx) { idx.has_context_tables(); idx.context_tables_max_length(); }) {
        if (_index.has_context_tables()) len = std::max<uint64_t>(1, _index.context_tables_max_length());
    }

    // the p initial parts cover at most half of the pattern
    return std::clamp<uint64_t>(len, 1, std::max<uint64_t>(1, m / (2 * p)));
}

template <typename idx_t, query_support_t query_support>
inline void apm_partition<idx_t, query_support>::partition_uniform(const inp_t& P, uint64_t p)
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

template <typename idx_t, query_support_t query_support>
inline void apm_partition<idx_t, query_support>::partition_dynamic(const inp_t& P, uint64_t p)
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

template <typename idx_t, query_support_t query_support>
inline void apm_partition<idx_t, query_support>::make_parts(const inp_t& P, uint64_t p, const apm_params_t& params)
{
    // the dynamic partitioning requires p >= 2 and m >= 2p
    if (p == 1 || params.partitioning == _partition_uniform || P.size() < 2 * p) {
        partition_uniform(P, p);
    } else {
        partition_dynamic(P, p);
    }
}

template <typename idx_t, query_support_t query_support>
template <typename occ_fnc_t>
inline std::vector<uint64_t> apm_partition<idx_t, query_support>::minimal_bounds(uint64_t m, uint64_t p,
    const apm_step_prediction& predictor, occ_fnc_t occ_of)
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

template <typename idx_t, query_support_t query_support>
inline void apm_partition<idx_t, query_support>::partition_optimized(const inp_t& P, const search_scheme_t& scheme,
    distance_metric_t metric, const apm_params_t& params)
{
    const uint64_t m = P.size();
    const uint64_t p = scheme.p;
    const uint64_t n = _index.input_size();
    _mirror_scheme = std::make_unique<search_scheme_t>(mirrored_scheme(scheme));

    // how often the anchor of the error model occurs (see apm_costs_t::anchor)
    uint64_t anchor_occ = 0;

    if (_costs != nullptr && params.anchored_search) {
        auto [a_beg, a_end] = _costs->anchor();
        if (a_beg < a_end) anchor_occ = num_occ_of(P, a_beg, a_end);
    }

    const std::vector<double>& sym_p = symbol_probabilities();
    // the strings of the searches have at most m + k symbols
    const std::vector<double> shares = apm_distinct_shares(_index, m + scheme.k);
    // [0] predicts the steps of the searches of the scheme, [1] those of the searches of its mirror
    _prediction[0].emplace(scheme, _costs, metric, P, _index.alphabet_size(), sym_p, shares, n, anchor_occ, params);
    _prediction[1].emplace(*_mirror_scheme, _costs, metric, P, _index.alphabet_size(), sym_p, shares, n, anchor_occ, params);
    apm_step_prediction* predictor[2] { &*_prediction[0], &*_prediction[1] };

    // 1. the parts of partition_dynamic
    partition_dynamic(P, p);
    const std::vector<uint64_t> bounds_dynamic = _bounds;

    // the parts are kept if every part is rare
    bool all_rare = true;
    for (uint64_t i = 0; i < p; i++) all_rare &= is_rare(i);
    if (all_rare) return;

    std::vector<uint64_t> occ(p); // [0..p-1] the number of occurrences of part i
    for (uint64_t i = 0; i < p; i++) occ[i] = num_occ(i);

    /* if the smaller prediction is at most the largest number of context extensions of step 2, only the direction is
       chosen */
    const double steps_scheme = predictor[0]->steps(_bounds, occ), steps_mirror = predictor[1]->steps(_bounds, occ);
    double best = std::min(steps_scheme, steps_mirror); // the smallest prediction so far
    _mirrored = steps_mirror < steps_scheme;
    if (best <= occurrences_steps(m)) return;

    // 2. the numbers of occurrences of the substrings of P
    const occurrence_table_t& table = occurrences(P);
    double table_steps = double(_occ_steps); // the steps of the optimization other than the predictions
    double steps = table_steps + double(predictor[0]->work() + predictor[1]->work()); // all steps of the optimization

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

    // returns whether P[b, e) is rare (see is_rare)
    auto is_rare_range = [&](uint64_t b, uint64_t e) { return double(occ_of(b, e)) * distinct_share(e - b + 1) <= 1; };

    /* returns the prediction of predictor[mirrored] for the parts with the bounds, and sets steps to the steps of the
       optimization so far */
    auto predict = [&](uint64_t mirrored, const std::vector<uint64_t>& bounds) {
        for (uint64_t i = 0; i < p; i++) occ[i] = occ_of(bounds[i], bounds[i + 1]);
        const double prediction = predictor[mirrored]->steps(bounds, occ);
        steps = table_steps + double(predictor[0]->work() + predictor[1]->work());
        return prediction;
    };

    // [0..1] the predictions for the parts of step 1 with the searches of the scheme and of its mirror
    const double prediction_dynamic[2] { predict(0, bounds_dynamic), predict(1, bounds_dynamic) };
    _mirrored = prediction_dynamic[1] < prediction_dynamic[0];
    best = std::min(prediction_dynamic[0], prediction_dynamic[1]);

    // 3. the local searches
    std::vector<uint64_t> bounds_uniform(p + 1);
    for (uint64_t i = 0; i <= p; i++) bounds_uniform[i] = (i * m) / p;
    std::vector<uint64_t> bounds_best = bounds_dynamic;
    bool mirrored_best = _mirrored;

    for (uint64_t mirrored = 0; mirrored < 2; mirrored++) {
        std::vector<std::vector<uint64_t>> starts { bounds_dynamic, bounds_uniform };
        const double cells_dp = double(p) * double(m) * double(m); // the cells of minimal_bounds

        if (steps + cells_dp <= best) {
            starts.emplace_back(minimal_bounds(m, p, *predictor[mirrored], occ_of));
            table_steps += cells_dp;
            steps = table_steps + double(predictor[0]->work() + predictor[1]->work());
        }

        for (std::vector<uint64_t> bounds : starts) {
            double prediction = bounds == bounds_dynamic ? prediction_dynamic[mirrored] : predict(mirrored, bounds);

            while (steps <= best) {
                const double steps_before = steps;
                // the move with the smallest prediction: bound i_move by one position to the right if move_right
                double prediction_move = prediction;
                uint64_t i_move = 0;
                bool move_right = false;

                for (uint64_t i = 1; i < p; i++) {
                    // a bound next to a rare part is not moved
                    if (is_rare_range(bounds[i - 1], bounds[i]) || is_rare_range(bounds[i], bounds[i + 1])) continue;

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

template <typename idx_t, query_support_t query_support>
inline apm_partition<idx_t, query_support>::apm_partition(const idx_t& index, const inp_t& P, const search_scheme_t& scheme,
    distance_metric_t metric, const apm_params_t& params)
    : _index(index), _bounds(scheme.p + 1), _occurs(scheme.p, 0), _has_context(scheme.p, 1), _costs(params.costs)
{
    const uint64_t p = scheme.p;
    assert(0 < p && p <= P.size());
    _ctx.reserve(p);
    for (uint64_t i = 0; i < p; i++) _ctx.emplace_back(index);

    if constexpr (can_optimize) {
        if (params.partitioning == _partition_optimized && p >= 2 && P.size() >= 2 * p) {
            partition_optimized(P, scheme, metric, params);
            return;
        }
    }

    make_parts(P, p, params);
}

template <typename idx_t, query_support_t query_support>
inline const std::vector<double>& apm_partition<idx_t, query_support>::symbol_probabilities() const
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

template <typename idx_t, query_support_t query_support>
inline auto apm_partition<idx_t, query_support>::occurrences(const inp_t& P) const -> const occurrence_table_t&
{
    if (!_occ_table.empty()) return _occ_table;
    const uint64_t m = P.size();
    _occ_table.resize(m);

    for (uint64_t i = 0; i < m; i++) {
        if (is_free(i)) {
            uint64_t occ = 0; // the occurrences of the free symbols of position i

            for (uint8_t c : _costs->free_symbols(i)) {
                auto [ext, ok] = ctx_t(_index).extend(sym_t(c), RIGHT);
                _occ_steps++;
                if (ok) occ += ext.num_occ();
            }

            _occ_table[i].emplace_back(occ);
            continue;
        }

        ctx_t ctx(_index); // the context of P[i, e)

        for (uint64_t e = i; e < m && !is_free(e); e++) {
            auto [ext, ok] = ctx.extend(P[e], RIGHT);
            _occ_steps++;
            const uint64_t occ = ok ? ext.num_occ() : 0;
            _occ_table[i].emplace_back(occ);
            if (occ <= 1) break;
            ctx = std::move(ext);
        }
    }

    return _occ_table;
}

#endif // MOVE_R_ALGORITHMS_MOVE_RB_APM_PARTITION_TPP
