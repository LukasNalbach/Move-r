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

#ifndef MOVE_R_ALGORITHMS_MOVE_RB_APM_MODEL_PLAN_HPP
#define MOVE_R_ALGORITHMS_MOVE_RB_APM_MODEL_PLAN_HPP

#include <algorithm>
#include <tuple>
#include <optional>
#include <stdexcept>
#include <vector>

#include <algorithms/move_rb/apm/model.hpp>
#include <algorithms/move_rb/apm/partition.hpp>

/**
 * @brief the partition of a pattern and the search scheme of a query with an error model
 *
 * The searches are those of the given scheme, with the upper bound of each step lowered to the number of errors a
 * valid alignment can have in the parts searched until then (see apm_costs_t::max_errors). A search whose lower bound
 * then exceeds its upper bound is dropped. For edit distance, all lower bounds are 0 once an upper bound is lowered or
 * the search bounds the cost (see apm_costs_t::bounds_costs).
 *
 * Instead, a single search runs if the pattern has a range of positions without errors and free symbols (an anchor)
 * that has at most as many occurrences as the first parts without errors of the searches together: it starts at the
 * anchor, extends it to the end of the pattern and then to its start.
 *
 * @tparam idx_t the bidirectional index type
 * @tparam query_support COUNT or LOCATE
 */
template <typename idx_t, query_support_t query_support>
class apm_model_plan {
public:
    using inp_t = typename idx_t::inp_t;
    using partition_t = apm_partition<idx_t, query_support>;

protected:
    partition_t _partition;
    std::optional<partition_t> _anchored; // the partition of the anchored search, if it runs
    search_scheme_t _scheme;

    /**
     * @brief returns the longest range [beg, end) of P whose positions are literal (see apm_costs_t::is_literal)
     *        and admit no deletion between them
     */
    static std::pair<uint64_t, uint64_t> anchor(const apm_costs_t& costs)
    {
        const uint64_t m = costs.length();
        uint64_t best_beg = 0, best_end = 0, beg = 0;

        for (uint64_t i = 0; i < m; i++) {
            if (!costs.is_literal(i)) { beg = i + 1; continue; }
            if (i > beg && costs.del(i) <= costs.bound()) beg = i;
            if (i + 1 - beg > best_end - best_beg) { best_beg = beg; best_end = i + 1; }
        }

        return {best_beg, best_end};
    }

public:
    /**
     * @brief plans a query with an error model
     * @param index the index to search in
     * @param P the pattern
     * @param scheme the search scheme (its k must be at least params.costs->max_errors())
     * @param params query parameters (the error model must be set)
     */
    apm_model_plan(const idx_t& index, const inp_t& P, const search_scheme_t& scheme, const apm_params_t& params)
        : _partition(index, P, scheme.p, params)
    {
        const apm_costs_t& costs = *params.costs;
        const uint64_t m = P.size();
        const uint64_t p = scheme.p;
        if (scheme.k < costs.max_errors()) {
            throw std::invalid_argument("the search scheme allows fewer errors than the model");
        }

        // the computed values of costs.max_errors(beg, end), as (beg, end, value)
        std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> max_errors;

        // returns costs.max_errors(beg, end)
        auto max_errors_in = [&](uint64_t beg, uint64_t end) {
            for (auto [b, e, v] : max_errors) if (b == beg && e == end) return v;
            max_errors.emplace_back(beg, end, costs.max_errors(beg, end));
            return std::get<2>(max_errors.back());
        };

        const uint64_t no_limit = std::numeric_limits<uint64_t>::max();
        auto add_bounded = [&](uint64_t a, uint64_t b) { return a > no_limit - b ? no_limit : a + b; };

        // the searches with lowered upper bounds; num_first_part_occ is the number of occurrences of their first
        // parts without errors (no_limit if a first part can have errors)
        _scheme.k = scheme.k;
        _scheme.p = p;
        _scheme.S.clear();
        bool lowered = false;
        uint64_t num_first_part_occ = 0;

        for (const search_t& S : scheme.S) {
            search_t T = S;
            uint64_t beg = _partition.beg(S[0].part), end = _partition.end(S[0].part);
            bool feasible = true;

            for (uint64_t t = 0; t < p; t++) {
                beg = std::min(beg, _partition.beg(S[t].part));
                end = std::max(end, _partition.end(S[t].part));
                T[t].k_max = uint8_t(std::min<uint64_t>(S[t].k_max, max_errors_in(beg, end)));
                lowered |= T[t].k_max < S[t].k_max;
                feasible &= T[t].k_min <= T[t].k_max;
            }

            if (!feasible) continue;
            num_first_part_occ = add_bounded(num_first_part_occ, T[0].k_max == 0 ? _partition.num_occ(T[0].part) : no_limit);
            _scheme.S.emplace_back(std::move(T));
        }

        // for edit distance, a lower bound can exclude a valid alignment if an upper bound or the costs exclude the
        // alignment of the same text string with the fewest errors; then the lower bounds are dropped
        if (costs.metric() == EDIT_DISTANCE && (lowered || costs.bounds_costs())) {
            for (search_t& S : _scheme.S) for (search_step_t& step : S) step.k_min = 0;
        }

        // the anchored search
        auto [a_beg, a_end] = anchor(costs);
        if (a_beg == a_end) return;

        std::vector<uint64_t> bounds { 0, a_beg, a_end, m };
        for (uint64_t x : _partition.bounds()) if (x < a_beg || x > a_end) bounds.emplace_back(x);
        std::sort(bounds.begin(), bounds.end());
        bounds.erase(std::unique(bounds.begin(), bounds.end()), bounds.end());

        partition_t anchored(index, P, bounds, params);
        const uint64_t q = std::find(bounds.begin(), bounds.end(), a_beg) - bounds.begin();
        if (anchored.num_occ(q) > num_first_part_occ) return;

        // the anchor, then the parts to its right, then the parts to its left
        search_t T;
        uint64_t beg = a_beg, end = a_end;

        auto add_step = [&](uint64_t part) {
            beg = std::min(beg, anchored.beg(part));
            end = std::max(end, anchored.end(part));
            const uint64_t k_max = std::min<uint64_t>(scheme.k, max_errors_in(beg, end));
            T.emplace_back(search_step_t{.part = uint8_t(part), .k_min = 0, .k_max = uint8_t(k_max)});
        };

        for (uint64_t part = q; part < anchored.size(); part++) add_step(part);
        for (uint64_t part = q; part-- > 0;) add_step(part);

        _scheme.p = anchored.size();
        _scheme.S = { std::move(T) };
        _anchored.emplace(std::move(anchored));
    }

    // returns the partition of the pattern
    const partition_t& partition() const { return _anchored ? *_anchored : _partition; }

    // returns the search scheme over the parts of partition()
    const search_scheme_t& scheme() const { return _scheme; }
};

#endif // MOVE_R_ALGORITHMS_MOVE_RB_APM_MODEL_PLAN_HPP
