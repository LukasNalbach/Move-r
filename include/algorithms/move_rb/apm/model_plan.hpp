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
 * @file model_plan.hpp
 * @brief the partition of a pattern and the search scheme of a query with an error model
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
#include <algorithms/move_rb/apm/step_prediction.hpp>

/**
 * @brief the partition of a pattern and the search scheme of a query with an error model
 *
 * One of three alternatives is run:
 * 1. The searches of the given scheme with lowered upper bounds: the upper bound of each step is lowered to the
 *    largest number of errors that a valid alignment can have in the parts searched until then (see
 *    apm_costs_t::max_errors), and a search with a lower bound above its upper bound is dropped. For the edit
 *    distance, all lower bounds are 0 if an upper bound is lowered or if the costs are computed in the search (see
 *    apm_costs_t::bounds_costs).
 * 2. The searches of the given scheme with the bounds of the scheme. For the edit distance, all lower bounds are 0 if
 *    the costs are computed in the search.
 * 3. The anchored search, if the pattern has an anchor, i.e., a range of positions without possible errors and without
 *    free symbols (see apm_costs_t::anchor): a single search that starts with the anchor, extends it to the end of the
 *    pattern and then to its start. Its lower bounds are 0, and its upper bounds are lowered as in alternative 1.
 *
 * The alternative with the smallest predicted number of steps is run (see apm_step_prediction::choose). Alternative 2 is
 * considered only for the edit distance if lowering an upper bound sets lower bounds to 0 that are kept otherwise.
 * With apm_params_t::lower_upper_bounds = false, alternative 2 replaces alternative 1, and the upper bounds of the
 * anchored search are not lowered; with apm_params_t::anchored_search = false, alternative 3 is not considered. For an
 * index without the alphabet size, n and r (see apm_partition::can_optimize), nothing is predicted: alternative 1 is
 * run, or alternative 2 if apm_params_t::lower_upper_bounds is false.
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
    partition_t _partition;               // the partition of the pattern for the searches of the scheme
    std::optional<partition_t> _anchored; // the partition of the anchored search, if it is run
    search_scheme_t _scheme;              // the searches that are run, over the parts of partition()

public:
    /**
     * @brief plans a query with an error model
     * @param index the index to search in
     * @param P the pattern
     * @param scheme the search scheme (its k must be at least params.costs->max_errors())
     * @param params the query parameters (params.costs must be set)
     */
    apm_model_plan(const idx_t& index, const inp_t& P, const search_scheme_t& scheme, const apm_params_t& params)
        : _partition(index, P, scheme, params.costs->metric(), params)
    {
        const apm_costs_t& costs = *params.costs;
        const uint64_t m = P.size();
        const uint64_t p = scheme.p;
        if (scheme.k < costs.max_errors()) {
            throw std::invalid_argument("the search scheme has fewer errors than the largest number of errors of the error model");
        }

        // the computed values of costs.max_errors(beg, end), as (beg, end, value)
        std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> max_errors;

        // returns costs.max_errors(beg, end)
        auto max_errors_in = [&](uint64_t beg, uint64_t end) {
            for (auto [b, e, v] : max_errors) if (b == beg && e == end) return v;
            max_errors.emplace_back(beg, end, costs.max_errors(beg, end));
            return std::get<2>(max_errors.back());
        };

        // the searches of the scheme or of its mirror (see apm_partition::mirrored)
        const search_scheme_t mirror = _partition.mirrored() ? mirrored_scheme(scheme) : search_scheme_t();
        const search_scheme_t& oriented = _partition.mirrored() ? mirror : scheme;

        // whether the upper bounds of the searches of the scheme are lowered
        bool lower = params.lower_upper_bounds;

        if constexpr (partition_t::can_optimize) {
            bool anchored = false; // whether the anchored search is run
            auto [a_beg, a_end] = costs.anchor();
            const bool has_anchor = a_beg < a_end && params.anchored_search;

            // the alternative with the smallest prediction, if there are several (see apm_step_prediction::choose)
            if (has_anchor || (lower && costs.metric() == EDIT_DISTANCE && !costs.bounds_costs())) {
                std::vector<uint64_t> occ(p); // [0..p-1] the number of occurrences of part i
                for (uint64_t i = 0; i < p; i++) occ[i] = _partition.num_occ(i);
                // the prediction of the partition, or a new one if the partition did not predict
                std::optional<apm_step_prediction> own;
                apm_step_prediction* predictor = _partition.prediction();

                if (predictor == nullptr) {
                    predictor = &own.emplace(oriented, &costs, costs.metric(), P, index.alphabet_size(),
                        _partition.symbol_probabilities(), apm_distinct_shares(index, m + scheme.k), index.input_size(),
                        has_anchor ? _partition.num_occ_of(P, a_beg, a_end) : 0, params);
                }

                const apm_step_prediction::choice_t choice = predictor->choose(_partition.bounds(), occ);
                lower = choice.lowered;
                anchored = choice.anchored;
            }

            /* the anchored search: the anchor, then the parts to its right, then the parts to its left, with the bounds of
               the parts of _partition outside the anchor */
            if (anchored) {
                std::vector<uint64_t> bounds { 0, a_beg, a_end, m };
                for (uint64_t b : _partition.bounds()) if (b < a_beg || b > a_end) bounds.emplace_back(b);
                std::sort(bounds.begin(), bounds.end());
                bounds.erase(std::unique(bounds.begin(), bounds.end()), bounds.end());
                const uint64_t anchor_part = std::find(bounds.begin(), bounds.end(), a_beg) - bounds.begin();
                const partition_t& parts = _anchored.emplace(index, P, std::move(bounds), params);
                search_t search;
                uint64_t beg = a_beg, end = a_end; // the positions searched so far are P[beg, end)

                // appends the step that searches the part with index part
                auto add_step = [&](uint64_t part) {
                    beg = std::min(beg, parts.beg(part));
                    end = std::max(end, parts.end(part));
                    const uint64_t k_max = lower ? std::min<uint64_t>(scheme.k, max_errors_in(beg, end)) : scheme.k;
                    search.emplace_back(search_step_t{.part = uint8_t(part), .k_min = 0, .k_max = uint8_t(k_max)});
                };

                for (uint64_t part = anchor_part; part < parts.size(); part++) add_step(part);
                for (uint64_t part = anchor_part; part-- > 0;) add_step(part);

                _scheme.k = scheme.k;
                _scheme.p = parts.size();
                _scheme.S = { std::move(search) };
                return;
            }
        }

        // the searches of the scheme, with lowered upper bounds if lower
        _scheme.k = scheme.k;
        _scheme.p = p;
        _scheme.S.clear();
        bool lowered = false; // whether an upper bound is lowered

        for (const search_t& S : oriented.S) {
            search_t search = S;
            uint64_t beg = _partition.beg(S[0].part), end = _partition.end(S[0].part); // the positions searched so far
            bool feasible = true; // whether no step has a lower bound above its upper bound

            for (uint64_t j = 0; j < p; j++) {
                beg = std::min(beg, _partition.beg(S[j].part));
                end = std::max(end, _partition.end(S[j].part));
                if (lower) search[j].k_max = uint8_t(std::min<uint64_t>(S[j].k_max, max_errors_in(beg, end)));
                lowered |= search[j].k_max < S[j].k_max;
                feasible &= search[j].k_min <= search[j].k_max;
            }

            if (feasible) _scheme.S.emplace_back(std::move(search));
        }

        /* for edit distance, a valid alignment can be missed because of a lower bound if the alignment of the same text
           string with the fewest errors exceeds an upper bound or the cost bound; then the lower bounds are set to 0 */
        if (costs.metric() == EDIT_DISTANCE && (lowered || costs.bounds_costs())) {
            for (search_t& S : _scheme.S) for (search_step_t& step : S) step.k_min = 0;
        }
    }

    // returns the partition of the pattern
    const partition_t& partition() const { return _anchored ? *_anchored : _partition; }

    // returns the search scheme over the parts of partition()
    const search_scheme_t& scheme() const { return _scheme; }
};

#endif // MOVE_R_ALGORITHMS_MOVE_RB_APM_MODEL_PLAN_HPP
