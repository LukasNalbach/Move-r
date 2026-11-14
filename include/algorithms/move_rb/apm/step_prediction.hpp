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
 * @file step_prediction.hpp
 * @brief predicts how many steps the searches of a search scheme take
 */

#ifndef MOVE_R_ALGORITHMS_MOVE_RB_APM_STEP_PREDICTION_HPP
#define MOVE_R_ALGORITHMS_MOVE_RB_APM_STEP_PREDICTION_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>
#include <unordered_map>
#include <vector>

#include <misc/apm.hpp>
#include <misc/search_schemes.hpp>
#include <algorithms/move_rb/apm/model.hpp>
#include <algorithms/move_rb/apm/search_base.hpp>

/**
 * @brief predicts how many steps the searches of a search scheme take for a partition of a pattern P into parts
 *
 * The number of steps of a search is the number of nodes of its search tree, as in the analysis of search schemes in
 * a random text (Kucherov, Salikhov, Tsur; Kianfar et al.). A search starts with the occurrences of its first part and
 * extends them by the strings whose numbers of errors lie within the error bounds of its steps. A(t) denotes the number
 * of these strings at depth t, i.e., after t further positions of P. At position i of P, a string continues in
 * match(i) ways without an error and in error(i) ways with an error (a substitution, an insertion or a deletion), where
 * a text symbol c counts as sigma * p_c ways, its weight (p_c is the probability of c in the text, sigma the effective
 * alphabet size, see the constructor).
 *
 * In a text with n symbols and r runs in its BWT, x = occ * r / n of the occ occurrences of the first part are assumed
 * to have distinct continuations. Each of the A(t) strings occurs among these continuations with a probability of at
 * most min(1, x / sigma^t), so the predicted number of steps of the search is the sum of A(t) * min(1, x / sigma^t)
 * over t >= 1. A search that matches its first part symbol by symbol (if the part has a free symbol or can have
 * errors, or if not apm_params_t::start_at_part_contexts) starts with the n occurrences of the empty string (x = r),
 * and the positions of its first part are counted in A(t). With an error model, the prediction is made for the
 * searches that apm_model_plan runs (see choose), and a string whose counters exceed their bounds is not counted.
 *
 * Occurrences near the pattern (only with the table of the numbers of occurrences of the substrings of P, see
 * set_occurrences): let P[beg, end) be the positions of P searched so far. The surplus of P[beg, end) is its number of
 * occurrences with distinct continuations minus the number the random text predicts, i.e., x times the probability
 * that a random string matches the positions of P[beg, end) after the start of the search. The surplus occurrences
 * are the occurrences near the pattern (in a repetitive text, P can resemble a string with many similar copies).
 * At the next position i, the share of them that match P[i] is the surplus after i divided by the surplus before i, at
 * most 1 (if the table does not contain both substrings, it is the probability that a random symbol matches P[i]); the
 * others continue with an error. Per number of errors e, the occurrences near the pattern are counted together with the
 * number of distinct strings they form: with an operation of weight w (w = 1 for a match), a string continues with at
 * most w strings and with at most as many strings as it has continuing occurrences. The predicted steps at
 * depth t are those of the random text plus the strings near the pattern; per number of errors e, they are at most the
 * number of all strings at depth t with e errors.
 */
class apm_step_prediction {
public:
    // the searches that apm_model_plan runs for a partition (see choose)
    struct choice_t {
        double steps;  // the predicted number of steps of the searches
        bool lowered;  // whether the upper bounds of the searches are lowered (see apm_model_plan)
        bool anchored; // whether the anchored search is run instead of the searches of the scheme
    };

private:
    // the error bounds of the searches
    struct bounds_mode_t {
        bool lowered;         // whether the upper bounds are lowered (see upper_bound)
        bool no_lower_bounds; // whether all lower bounds are 0
    };

    const search_scheme_t& _scheme;
    const apm_costs_t* _costs; // the error model, or nullptr
    double _sigma;             // the effective alphabet size (see the constructor)
    double _runs_per_symbol;   // r / n: of the occ occurrences of a string, occ * r / n have distinct continuations
    uint64_t _n;               // n, the text length
    /* the anchor P[_anchor_beg, _anchor_end) of the error model (empty if the anchored search is not considered) and
       its number of occurrences */
    uint64_t _anchor_beg = 0, _anchor_end = 0, _anchor_occ = 0;
    bool _may_lower;     // see apm_params_t::lower_upper_bounds
    bool _from_contexts; // see apm_params_t::start_at_part_contexts

    // [0..m-1] 1 <=> a symbol other than P[i] is free at position i (see apm_costs_t::has_free)
    std::vector<uint8_t> _free;

    /* the numbers of ways in which a string continues at position i of P (see the class description): [0..m-1]
       _match[i] without an error, _sub[i] with a substitution, _ins[i] with an insertion, _del[i] with a deletion and
       _error[i] = _sub[i] + _ins[i] + _del[i] with an error */
    std::vector<double> _match, _sub, _ins, _del, _error;

    /* the table of the numbers of occurrences of the substrings of P: [i][l - 1] the number of occurrences of
       P[i, i + l) (see apm_partition::occurrences), or nullptr */
    const std::vector<std::vector<uint64_t>>* _substring_occ = nullptr;

    // the operations with which a string continues at a position: a match, a substitution, an insertion and a deletion
    static constexpr uint64_t num_operations = 4;

    // a group of the strings at depth t that have the same number of errors, the same counters and the same last operation
    struct counted_t {
        uint64_t counters;   // the packed counters (see apm_costs_t)
        uint8_t errors;      // the number of errors
        uint8_t last;        // the last operation: 0 (match or substitution), 1 (insertion) or 2 (deletion)
        double num;          // the number of these strings divided by sigma^t
        double near_occ;     // the occurrences near the pattern of these strings (see the class description)
        double near_strings; // the number of these strings with occurrences near the pattern

        // returns the values that identify the group
        auto key() const { return std::tuple(errors, last, counters); }
    };

    // the occurrences near the pattern of a search (see the class description)
    struct near_t {
        uint64_t beg, end;           // the positions of P searched so far are P[beg, end)
        double x;                    // the occurrences of the start of the search with distinct continuations
        double match_prob = 1;       // the probability that a random string matches the positions searched after the start
        double surplus = 0;          // the occurrences of P[beg, end) with distinct continuations minus x * match_prob
        // whether the table contains P[beg, end) (else its longest prefix in the table is used, see advance)
        bool counted = true;
        double match_share = 1;      // the share of the occurrences near the pattern that match the last searched position
        bool error_free_kept = true; // whether the lower bounds of the steps so far keep the strings without an error
        double power = 1;            // sigma^t at depth t (at most the largest double)

        // returns the occurrences near the pattern without an error
        double error_free() const { return error_free_kept ? surplus : 0; }
    };

    /**
     * @brief appends position i to the positions searched so far: updates the surplus and the share of the occurrences
     *        near the pattern that match P[i]
     * @param near the occurrences near the pattern of a search
     * @param i the position of P
     */
    void advance(near_t& near, uint64_t i) const
    {
        near.beg = std::min(near.beg, i);
        near.end = std::max(near.end, i + 1);
        near.match_prob *= _match[i] / _sigma;
        if (near.power < std::numeric_limits<double>::max() / _sigma) near.power *= _sigma;
        // the occurrences of P[beg, end) with distinct continuations (those of its longest prefix in the table if not counted)
        const std::vector<uint64_t>& row = (*_substring_occ)[near.beg];
        const bool counted = near.end - near.beg <= row.size();
        const double occ = double(row[std::min<uint64_t>(near.end - near.beg, row.size()) - 1]) * _runs_per_symbol;
        const double surplus = std::max(0.0, occ - near.x * near.match_prob);
        near.match_share = near.counted && counted && near.surplus > 0 ? std::min(1.0, surplus / near.surplus) : _match[i] / _sigma;
        near.surplus = surplus;
        near.counted = counted;
    }

    /* returns the share of the occurrences near the pattern that continue at position i with errors of weight w:
       (1 - match_share) * w / max(error(i), sigma - match(i)), i.e., the occurrences that do not match P[i] are split
       by weight among the other symbols, or among the errors if they have the larger weight */
    double error_share(const near_t& near, uint64_t i, double w) const
    {
        const double other = std::max(_error[i], _sigma - _match[i]);
        return other > 0 ? (1 - near.match_share) * w / other : 0;
    }

    /* returns the number of strings with which num strings near the pattern with occ occurrences continue if a share of
       these occurrences continues with an operation of weight w (w = 1 for a match): a string continues with at most w
       strings and with at most as many strings as it has continuing occurrences (occ / num * share) */
    static double continuing_strings(double num, double occ, double share, double w)
    {
        return num > 0 ? num * std::min(w, occ / num * share) : 0;
    }

    /* returns the steps that num strings near the pattern add to the x * b steps that the random text predicts for the
       b * sigma^t strings at depth t with one number of errors; the sum is at most b * sigma^t */
    static double near_steps(double b, double num, const near_t& near)
    {
        if (num <= 0 || b <= 0) return 0;
        const double all = b * near.power;
        return std::min(all, near.x * b + num) - std::min(all, near.x * b);
    }

    // the computed values of _costs->max_errors(beg, end), by beg * (m + 1) + end
    mutable std::unordered_map<uint64_t, uint64_t> _max_errors;

    // [0..p-1] _searches_of[i] = the indices of the searches that start with part i
    std::vector<std::vector<uint64_t>> _searches_of;

    /* the sums of part_steps (computed by its first call), for the bounds of the uniform partition: [0..|S|-1][t]
       _sum_a[s][t] = A(1) + ... + A(t) of search s for the t with sigma^t <= r, and _sum_b[s][t] = the sum of
       A(t') / sigma^t' over t' > t */
    mutable std::vector<std::vector<double>> _sum_a, _sum_b;

    /* the steps of the predictions so far: the cells (position, number of errors) and the transitions (group of strings,
       operation) computed for A(t) */
    mutable uint64_t _work = 0;

    /* returns the upper bound of a step of a search that has searched P[beg, end) and has the upper bound k_max, lowered
       to _costs->max_errors(beg, end) if lowered */
    uint64_t upper_bound(uint64_t beg, uint64_t end, uint64_t k_max, bool lowered) const
    {
        if (_costs == nullptr || !lowered) return k_max;
        auto [it, inserted] = _max_errors.try_emplace(beg * (_match.size() + 1) + end, 0);
        if (inserted) it->second = _costs->max_errors(beg, end);
        return std::min(k_max, it->second);
    }

    /**
     * @brief appends position i to the groups of strings with counters (see for_each_counted_step): each group continues
     *        with a match, a substitution, an insertion and a deletion, whose counter increments (and those of opening a
     *        gap) are added to its counters; the strings whose counters exceed their bounds or that have more than k_max
     *        errors are dropped, and the occurrences near the pattern are split among the operations as in for_each_step
     * @param groups the groups at depth t - 1, set to the groups at depth t
     * @param next a buffer
     * @param i the position of P
     * @param k_max the upper bound on the number of errors
     * @param near the occurrences near the pattern of the search
     */
    void append_counted(std::vector<counted_t>& groups, std::vector<counted_t>& next, uint64_t i, uint64_t k_max,
        const near_t& near) const
    {
        next.clear();
        _work += num_operations * groups.size();

        for (const counted_t& g : groups) {
            // appends the strings that continue g with an operation of weight w, with a share of its occurrences near the pattern
            auto add = [&](uint64_t counters, uint8_t errors, uint8_t last, double w, double share) {
                if (w == 0 || errors > k_max || _costs->exceeds(counters)) return;
                next.emplace_back(counted_t{counters, errors, last, g.num * w / _sigma, g.near_occ * share,
                    continuing_strings(g.near_strings, g.near_occ, share, last == 0 && errors == g.errors ? 1 : w)});
            };

            add(g.counters, g.errors, 0, _match[i], near.match_share);
            add(g.counters + _costs->sub_count(i), g.errors + 1, 0, _sub[i], error_share(near, i, _sub[i]));
            add(g.counters + _costs->ins_count(i) + (g.last == 1 ? 0 : _costs->ins_open_count(i)), g.errors + 1, 1,
                _ins[i], error_share(near, i, _ins[i]));
            add(g.counters + _costs->del_count(i) + (g.last == 2 ? 0 : _costs->del_open_count(i)), g.errors + 1, 2,
                _del[i], error_share(near, i, _del[i]));
        }

        // the groups with the same key are merged
        std::sort(next.begin(), next.end(), [](const counted_t& a, const counted_t& b) { return a.key() < b.key(); });
        groups.clear();

        for (const counted_t& g : next) {
            if (!groups.empty() && groups.back().key() == g.key()) {
                groups.back().num += g.num;
                groups.back().near_occ += g.near_occ;
                groups.back().near_strings += g.near_strings;
            } else {
                groups.emplace_back(g);
            }
        }
    }

    /* returns the occurrences near the pattern at the start of the search S, whose start has x occurrences with distinct
       continuations: the positions searched so far are the first part of S, or none if from_root */
    near_t near_start(const search_t& S, const std::vector<uint64_t>& bounds, bool from_root, double x) const
    {
        const uint64_t part_beg = bounds[S[0].part], part_end = bounds[S[0].part + 1];
        if (!from_root) return { .beg = part_beg, .end = part_end, .x = x };
        const uint64_t start = apm_step_directions(S)[0] == RIGHT ? part_beg : part_end;
        return { .beg = start, .end = start, .x = x };
    }

    /**
     * @brief computes A(t) of the search S depth by depth (see the class description) and calls step(b, steps_near) for
     *        each depth t, with b = A(t) / sigma^t and steps_near the steps that the strings near the pattern add at
     *        depth t (0 without set_occurrences or if x = 0); with counters in the error model, see
     *        for_each_counted_step
     * @param S the search
     * @param bounds [0..p] the bounds of the parts
     * @param mode the error bounds of the steps
     * @param from_root whether S matches its first part symbol by symbol (see starts_at_root)
     * @param x the occurrences of the start of S with distinct continuations
     * @param step the function called for each depth
     */
    template <typename step_fnc_t>
    void for_each_step(const search_t& S, const std::vector<uint64_t>& bounds, bounds_mode_t mode, bool from_root,
        double x, step_fnc_t step) const
    {
        if (_costs != nullptr && _costs->num_counters() > 0) {
            for_each_counted_step(S, bounds, mode, from_root, x, step);
            return;
        }

        const std::vector<direction_t> dirs = apm_step_directions(S);
        const bool has_near = _substring_occ != nullptr && x > 0;
        near_t near = near_start(S, bounds, from_root, x);
        // [0..k] b_by_errors[e] = A_e(t) / sigma^t, where A_e(t) is the number of strings at depth t with e errors
        std::vector<double> b_by_errors(_scheme.k + 1, 0.0);
        // [0..k] near_occ[e] = the occurrences near the pattern with e errors, near_strings[e] = their number of strings
        std::vector<double> near_occ(_scheme.k + 1, 0.0), near_strings(_scheme.k + 1, 0.0);
        b_by_errors[0] = 1;
        uint64_t beg = bounds[S[0].part], end = bounds[S[0].part + 1]; // the positions searched so far are P[beg, end)

        for (uint64_t j = from_root ? 0 : 1; j < S.size(); j++) {
            const uint64_t part_beg = bounds[S[j].part], part_end = bounds[S[j].part + 1];
            beg = std::min(beg, part_beg);
            end = std::max(end, part_end);
            const uint64_t k_max = upper_bound(beg, end, S[j].k_max, mode.lowered);
            // the strings with more than k_max errors are dropped
            for (uint64_t e = k_max + 1; e < b_by_errors.size(); e++) b_by_errors[e] = near_occ[e] = near_strings[e] = 0;

            _work += (part_end - part_beg) * (k_max + 1);

            for (uint64_t y = 0; y < part_end - part_beg; y++) {
                const uint64_t i = dirs[j] == RIGHT ? part_beg + y : part_end - 1 - y;
                // depth t from depth t - 1: a match keeps the number of errors, an error increases it by 1
                for (uint64_t e = k_max; e > 0; e--) b_by_errors[e] = (b_by_errors[e] * _match[i] + b_by_errors[e - 1] * _error[i]) / _sigma;
                b_by_errors[0] = b_by_errors[0] * _match[i] / _sigma;
                double b = 0;
                for (double b_e : b_by_errors) b += b_e;
                double steps_near = 0;

                // the occurrences near the pattern and their strings at depth t
                if (has_near) {
                    advance(near, i);
                    const double share_error = error_share(near, i, _error[i]);

                    for (uint64_t e = k_max; e > 0; e--) {
                        near_strings[e] = continuing_strings(near_strings[e], near_occ[e], near.match_share, 1) +
                            continuing_strings(near_strings[e - 1], near_occ[e - 1], share_error, _error[i]);
                        near_occ[e] = near_occ[e] * near.match_share + near_occ[e - 1] * share_error;
                    }

                    near_occ[0] = near.error_free();
                    near_strings[0] = std::min(1.0, near_occ[0]);
                    for (uint64_t e = 0; e <= k_max; e++) steps_near += near_steps(b_by_errors[e], near_strings[e], near);
                }

                step(b, steps_near);
            }

            // the strings with fewer than S[j].k_min errors are dropped
            if (!mode.no_lower_bounds) for (uint64_t e = 0; e < S[j].k_min; e++) b_by_errors[e] = near_occ[e] = near_strings[e] = 0;
            if (!mode.no_lower_bounds && S[j].k_min > 0) near.error_free_kept = false;
        }
    }

    /* for_each_step with counters in the error model: the strings are grouped by their numbers of errors, their counters
       and their last operations (see append_counted) */
    template <typename step_fnc_t>
    void for_each_counted_step(const search_t& S, const std::vector<uint64_t>& bounds, bounds_mode_t mode,
        bool from_root, double x, step_fnc_t step) const
    {
        const std::vector<direction_t> dirs = apm_step_directions(S);
        const bool has_near = _substring_occ != nullptr && x > 0;
        near_t near = near_start(S, bounds, from_root, x);
        std::vector<counted_t> groups { counted_t{0, 0, 0, 1, 0, 0} }, next; // the groups at depth t, and a buffer
        uint64_t beg = bounds[S[0].part], end = bounds[S[0].part + 1]; // the positions searched so far are P[beg, end)

        for (uint64_t j = from_root ? 0 : 1; j < S.size(); j++) {
            const uint64_t part_beg = bounds[S[j].part], part_end = bounds[S[j].part + 1];
            beg = std::min(beg, part_beg);
            end = std::max(end, part_end);
            const uint64_t k_max = upper_bound(beg, end, S[j].k_max, mode.lowered);

            for (uint64_t y = 0; y < part_end - part_beg; y++) {
                const uint64_t i = dirs[j] == RIGHT ? part_beg + y : part_end - 1 - y;
                if (has_near) advance(near, i);
                append_counted(groups, next, i, k_max, near);
                double b = 0, steps_near = 0;

                for (counted_t& g : groups) {
                    b += g.num;
                    if (!has_near) continue;
                    if (g.errors == 0) {
                        g.near_occ = near.error_free();
                        g.near_strings = std::min(1.0, g.near_occ);
                    }

                    steps_near += near_steps(g.num, g.near_strings, near);
                }

                step(b, steps_near);
            }

            // the strings with fewer than S[j].k_min errors are dropped
            if (!mode.no_lower_bounds) std::erase_if(groups, [&](const counted_t& g) { return g.errors < S[j].k_min; });
            if (!mode.no_lower_bounds && S[j].k_min > 0) near.error_free_kept = false;
        }
    }

    /**
     * @brief returns the predicted number of steps of the search S if its start has occ occurrences: the sum of
     *        min(A(t), x * A(t) / sigma^t) over the depths t, with x = occ * r / n, plus the steps of the strings near
     *        the pattern (see for_each_step)
     */
    double search_steps(const search_t& S, const std::vector<uint64_t>& bounds, uint64_t occ, bounds_mode_t mode,
        bool from_root) const
    {
        const double x = double(occ) * _runs_per_symbol;
        double power = 1; // sigma^t while it is at most x
        double steps = 0;

        for_each_step(S, bounds, mode, from_root, x, [&](double b, double steps_near) {
            if (power <= x) power *= _sigma;
            steps += (power <= x ? b * power : x * b) + steps_near;
        });

        return steps;
    }

    // returns whether a symbol other than P[i] is free at a position i of P[beg, end) (then the part has no search context)
    bool has_free(uint64_t beg, uint64_t end) const
    {
        for (uint64_t i = beg; i < end; i++) if (_free[i]) return true;
        return false;
    }

    // returns whether the search S matches its first part symbol by symbol, starting with the n occurrences of the empty string
    bool starts_at_root(const search_t& S, const std::vector<uint64_t>& bounds, bool lowered) const
    {
        const uint64_t beg = bounds[S[0].part], end = bounds[S[0].part + 1];
        return !_from_contexts || has_free(beg, end) || upper_bound(beg, end, S[0].k_max, lowered) > 0;
    }

    // returns whether a step of the search S has a lower bound above its upper bound (lowered if lowered)
    bool is_infeasible(const search_t& S, const std::vector<uint64_t>& bounds, bool lowered) const
    {
        uint64_t beg = bounds[S[0].part], end = bounds[S[0].part + 1];

        for (const search_step_t& step : S) {
            beg = std::min(beg, bounds[step.part]);
            end = std::max(end, bounds[step.part + 1]);
            if (step.k_min > upper_bound(beg, end, step.k_max, lowered)) return true;
        }

        return false;
    }

    // returns whether lowering the upper bounds changes the upper bound of a step of a search
    bool lowers_a_bound(const std::vector<uint64_t>& bounds) const
    {
        for (const search_t& S : _scheme.S) {
            uint64_t beg = bounds[S[0].part], end = bounds[S[0].part + 1];

            for (const search_step_t& step : S) {
                beg = std::min(beg, bounds[step.part]);
                end = std::max(end, bounds[step.part + 1]);
                if (upper_bound(beg, end, step.k_max, true) < step.k_max) return true;
            }
        }

        return false;
    }

    /**
     * @brief returns the predicted number of steps of the searches of the scheme with the error bounds of mode; a search
     *        with a lower bound above an upper bound is not run and not counted
     * @param bounds [0..p] the bounds of the parts
     * @param occ [0..p-1] the numbers of occurrences of the parts
     * @param mode the error bounds
     */
    double scheme_steps(const std::vector<uint64_t>& bounds, const std::vector<uint64_t>& occ, bounds_mode_t mode) const
    {
        double steps = 0;

        for (const search_t& S : _scheme.S) {
            if (is_infeasible(S, bounds, mode.lowered)) continue;
            const bool from_root = starts_at_root(S, bounds, mode.lowered);
            steps += search_steps(S, bounds, from_root ? _n : std::min(occ[S[0].part], _n), mode, from_root);
        }

        return steps;
    }

    /**
     * @brief returns the predicted number of steps of the anchored search of apm_model_plan: the anchor, then the parts
     *        to its right, then the parts to its left, with the bounds of the parts outside the anchor
     * @param bounds [0..p] the bounds of the parts
     */
    double anchored_steps(const std::vector<uint64_t>& bounds) const
    {
        // the bounds of the parts of the anchored search: the bounds of the anchor and those in bounds outside the anchor
        std::vector<uint64_t> anchored_bounds { 0, _anchor_beg, _anchor_end, _match.size() };
        for (uint64_t b : bounds) if (b < _anchor_beg || b > _anchor_end) anchored_bounds.emplace_back(b);
        std::sort(anchored_bounds.begin(), anchored_bounds.end());
        anchored_bounds.erase(std::unique(anchored_bounds.begin(), anchored_bounds.end()), anchored_bounds.end());
        // the index of the part that is the anchor
        const uint64_t anchor_part = std::find(anchored_bounds.begin(), anchored_bounds.end(), _anchor_beg) - anchored_bounds.begin();
        search_t search;
        for (uint64_t part = anchor_part; part + 1 < anchored_bounds.size(); part++) search.emplace_back(search_step_t{uint8_t(part), 0, _scheme.k});
        for (uint64_t part = anchor_part; part-- > 0;) search.emplace_back(search_step_t{uint8_t(part), 0, _scheme.k});
        const bool from_root = starts_at_root(search, anchored_bounds, _may_lower);
        return search_steps(search, anchored_bounds, from_root ? _n : _anchor_occ, {_may_lower, true}, from_root);
    }

public:
    /**
     * @brief prepares the prediction for a pattern P: the effective alphabet size sigma is 1 / (sum_c p_c^2), where p_c
     *        is the probability of symbol c in the text, and a text symbol c counts as sigma * p_c ways (without the
     *        probabilities, sigma is the alphabet size and each symbol counts as 1 way). At position i, a string
     *        continues without an error with P[i], and with an error by a substitution of P[i] by another symbol and,
     *        for the edit distance, by an insertion of P[i] or a deletion of any text symbol. With an error model, a
     *        string continues without an error with the free symbols of position i, and the errors are those with a
     *        cost of at most the cost bound.
     * @param scheme the search scheme
     * @param costs the error model compiled for P, or nullptr
     * @param metric the distance metric of the searches
     * @param P the pattern
     * @param sigma the alphabet size of the text
     * @param sym_p [0..255] the probability of each byte in the text (see apm_partition::symbol_probabilities), or empty
     * @param r the number of runs in the BWT of the text
     * @param n the text length
     * @param anchor_occ how often the anchor of the error model occurs in the text (see apm_costs_t::anchor)
     * @param params the query parameters (see apm_params_t::start_at_part_contexts, apm_params_t::lower_upper_bounds and
     *        apm_params_t::anchored_search)
     */
    template <typename pat_t>
    apm_step_prediction(const search_scheme_t& scheme, const apm_costs_t* costs, distance_metric_t metric, const pat_t& P,
        uint64_t sigma, const std::vector<double>& sym_p, uint64_t r, uint64_t n, uint64_t anchor_occ,
        const apm_params_t& params)
        : _scheme(scheme), _costs(costs), _sigma(double(sigma)), _runs_per_symbol(double(r) / double(n)), _n(n),
          _anchor_occ(anchor_occ), _may_lower(params.lower_upper_bounds), _from_contexts(params.start_at_part_contexts),
          _free(P.size(), 0), _match(P.size(), 0), _sub(P.size(), 0), _ins(P.size(), 0), _del(P.size(), 0),
          _error(P.size()), _searches_of(scheme.p)
    {
        using sym_t = std::remove_cvref_t<decltype(P[0])>;
        const uint64_t m = P.size();
        const bool edit = metric == EDIT_DISTANCE;
        double collision = 0; // sum_c p_c^2, the probability that two random symbols of the text are equal
        for (double p : sym_p) collision += p * p;
        if (collision > 0) _sigma = 1 / collision;

        // returns the number of ways that byte c counts as: sigma * p_c, or 1 without the probabilities
        auto weight = [&](uint64_t c) { return sym_p.empty() ? 1.0 : c < sym_p.size() ? sym_p[c] * _sigma : 0.0; };

        // [0..num_ranks-1] rank_weight[x] = the weight of the byte with rank x in the error model
        std::vector<double> rank_weight;

        if (costs != nullptr) {
            rank_weight.resize(costs->num_ranks(), 0);

            for (uint64_t c = 0; c <= std::numeric_limits<uint8_t>::max(); c++) {
                const uint64_t rank = costs->rank(uint8_t(c));
                if (rank < costs->num_ranks()) rank_weight[rank] = weight(c);
            }
        }

        for (uint64_t i = 0; i < m; i++) {
            if (costs == nullptr) {
                _match[i] = weight(uint64_t(std::make_unsigned_t<sym_t>(P[i])));
                _sub[i] = _sigma - _match[i];
                _ins[i] = edit ? 1 : 0;
                _del[i] = edit ? _sigma : 0;
            } else {
                _free[i] = costs->has_free(i);

                for (uint64_t rank = 0; rank < costs->num_ranks(); rank++) {
                    if (costs->is_free(i, rank)) _match[i] += rank_weight[rank];
                    else if (costs->sub(i, rank) <= costs->bound()) _sub[i] += rank_weight[rank];
                }

                _ins[i] = edit && costs->ins(i) <= costs->bound() ? 1 : 0;
                _del[i] = edit && costs->del(i) <= costs->bound() ? _sigma : 0;
            }

            _error[i] = _sub[i] + _ins[i] + _del[i];
        }

        if (costs != nullptr && params.anchored_search) std::tie(_anchor_beg, _anchor_end) = costs->anchor();
        for (uint64_t s = 0; s < scheme.S.size(); s++) _searches_of[scheme.S[s][0].part].emplace_back(s);
    }

    // returns the steps of the predictions so far (see _work)
    uint64_t work() const { return _work; }

    /**
     * @brief sets the numbers of occurrences of the substrings of P, with which the predictions count the occurrences
     *        near the pattern (see the class description)
     * @param occ [i][l - 1] the number of occurrences of P[i, i + l) (see apm_partition::occurrences), or nullptr
     */
    void set_occurrences(const std::vector<std::vector<uint64_t>>* occ) { _substring_occ = occ; }

    /**
     * @brief returns the predicted number of steps of the searches that start with part i if part i has occ
     *        occurrences; the other parts have the bounds of the uniform partition, and the occurrences near the pattern
     *        are not counted (see apm_partition::minimal_bounds)
     */
    double part_steps(uint64_t i, uint64_t occ) const
    {
        if (_sum_a.empty()) {
            const uint64_t m = _match.size();
            const double r = _runs_per_symbol * double(_n);
            const bounds_mode_t mode { _may_lower, _costs != nullptr && _costs->metric() == EDIT_DISTANCE };
            // [0..p] the bounds of the parts of the uniform partition
            std::vector<uint64_t> bounds(_scheme.p + 1);
            for (uint64_t j = 0; j <= _scheme.p; j++) bounds[j] = (j * m) / _scheme.p;
            _sum_a.resize(_scheme.S.size());
            _sum_b.resize(_scheme.S.size());

            for (uint64_t s = 0; s < _scheme.S.size(); s++) {
                std::vector<double> b { 0 }; // b[t] = A(t) / sigma^t of search s for the depths t >= 1
                for_each_step(_scheme.S[s], bounds, mode, false, 0, [&](double b_t, double) { b.emplace_back(b_t); });
                double power = 1; // sigma^t
                _sum_a[s] = { 0 };

                for (uint64_t t = 1; t < b.size() && power * _sigma <= r; t++) {
                    power *= _sigma;
                    _sum_a[s].emplace_back(_sum_a[s].back() + b[t] * power);
                }

                _sum_b[s].resize(b.size(), 0);
                for (uint64_t t = b.size() - 1; t > 0; t--) _sum_b[s][t - 1] = _sum_b[s][t] + b[t];
            }
        }

        const double x = double(occ) * _runs_per_symbol;
        double steps = 0;

        for (uint64_t s : _searches_of[i]) {
            // the largest t with sigma^t <= x
            const uint64_t t = x < _sigma ? 0 : std::min<uint64_t>(uint64_t(std::log(x) / std::log(_sigma)),
                std::min(_sum_a[s].size(), _sum_b[s].size()) - 1);
            steps += _sum_a[s][t] + x * _sum_b[s][t];
        }

        return steps;
    }

    /**
     * @brief returns the searches that apm_model_plan runs for the parts with the bounds: without an error model, the
     *        searches of the scheme; with an error model, the alternative with the smallest prediction among
     *        1. the searches of the scheme with lowered upper bounds (if apm_params_t::lower_upper_bounds); for the edit
     *           distance, their lower bounds are 0 if an upper bound is lowered or if the costs are computed in the
     *           search (see apm_costs_t::bounds_costs),
     *        2. the searches of the scheme with the bounds of the scheme (if not apm_params_t::lower_upper_bounds, or for the
     *           edit distance if lowering an upper bound sets lower bounds to 0 that are kept otherwise); for the edit
     *           distance, their lower bounds are 0 if the costs are computed in the search, and
     *        3. the anchored search (if the error model has an anchor and apm_params_t::anchored_search).
     * @param bounds [0..p] the bounds of the parts
     * @param occ [0..p-1] the numbers of occurrences of the parts (read only for the first parts of the searches that
     *        start with the search context of their first part)
     */
    choice_t choose(const std::vector<uint64_t>& bounds, const std::vector<uint64_t>& occ) const
    {
        if (_costs == nullptr) return { scheme_steps(bounds, occ, {false, false}), false, false };

        const bool edit = _costs->metric() == EDIT_DISTANCE;
        const bool bounds_costs = _costs->bounds_costs();
        const bool lowers = _may_lower && lowers_a_bound(bounds);
        choice_t best { std::numeric_limits<double>::infinity(), false, false };

        // keeps an alternative if it has the smallest prediction so far
        auto consider = [&](double steps, bool lowered, bool anchored) {
            if (steps < best.steps) best = { steps, lowered, anchored };
        };

        if (_may_lower) consider(scheme_steps(bounds, occ, {true, edit && (lowers || bounds_costs)}), true, false);
        if (!_may_lower || (edit && lowers && !bounds_costs)) {
            consider(scheme_steps(bounds, occ, {false, edit && bounds_costs}), false, false);
        }
        if (_anchor_beg < _anchor_end) consider(anchored_steps(bounds), _may_lower, true);
        return best;
    }

    /**
     * @brief returns the predicted number of steps of the searches that apm_model_plan runs (see choose)
     * @param bounds [0..p] the bounds of the parts
     * @param occ [0..p-1] the numbers of occurrences of the parts
     */
    double steps(const std::vector<uint64_t>& bounds, const std::vector<uint64_t>& occ) const
    {
        return choose(bounds, occ).steps;
    }
};

#endif // MOVE_R_ALGORITHMS_MOVE_RB_APM_STEP_PREDICTION_HPP
