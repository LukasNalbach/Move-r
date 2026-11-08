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
 * @file search_cost.hpp
 * @brief the predicted number of steps of the searches of a search scheme
 */

#ifndef MOVE_R_ALGORITHMS_MOVE_RB_APM_SEARCH_COST_HPP
#define MOVE_R_ALGORITHMS_MOVE_RB_APM_SEARCH_COST_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <tuple>
#include <utility>
#include <unordered_map>
#include <vector>

#include <misc/search_schemes.hpp>
#include <algorithms/move_rb/apm/model.hpp>
#include <algorithms/move_rb/apm/search_base.hpp>

/**
 * @brief the predicted number of steps of the searches of a search scheme for a partition of a pattern (the number of
 *        nodes of the search trees in a random text, Kucherov et al., Kianfar et al.): after its first part, a search
 *        extends the occurrences of the part by the strings that respect its error bounds; A(t) of these strings have
 *        t symbols, where position i of the pattern has match(i) symbols without an error and error(i) errors. Of the
 *        occ occurrences of the first part, occ * r / n have distinct continuations (r runs in the BWT of a text with
 *        n symbols), and at most occ * r / n / sigma^t of the A(t) strings occur, so a search takes
 *        sum_t A(t) * min(1, occ * r / n / sigma^t) steps. With an error model, the prediction is that of the searches
 *        that apm_model_plan runs (with the error bounds it lowers), and the strings respect the counters of the model.
 */
class apm_search_cost {
    const search_scheme_t& _scheme;
    const apm_costs_t* _costs; // the error model, or nullptr
    double _sigma;             // the alphabet size
    double _share;             // r / n
    uint64_t _n;               // the length of the text
    uint64_t _anchor_beg = 0, _anchor_end = 0, _anchor_occ = 0; // the anchor of the error model and its occurrences

    /* [0..m-1] at position i: the numbers of symbols without an error, of substitutions, of insertions and of deletions,
       and of all errors */
    std::vector<double> _match, _sub, _ins, _del, _error;

    // the operations that append a position: a match, a substitution, an insertion and a deletion
    static constexpr uint64_t num_operations = 4;

    // a string of t symbols after the first part: its number of errors, its counters and its last operation
    struct counted_t {
        uint64_t counters;
        uint8_t errors;
        uint8_t last; // the operation before the next one: 0 (match or substitution), 1 (insertion) or 2 (deletion)
        double num;   // the number of these strings divided by sigma^t

        auto key() const { return std::tuple(errors, last, counters); }
    };

    // the computed values of _costs->max_errors(beg, end), by beg * (m + 1) + end
    mutable std::unordered_map<uint64_t, uint64_t> _max_errors;

    // for each part, the searches that start with it
    std::vector<std::vector<uint64_t>> _searches_of;

    /* for the uniform partition, per search (computed by the first call of part_cost): _sum_a[s][t] = A(1) + ... + A(t)
       for sigma^t <= r, and _sum_b[s][t] = the sum of A(t') / sigma^t' over t' > t */
    mutable std::vector<std::vector<double>> _sum_a, _sum_b;

    mutable uint64_t _work = 0; // the cells and the transitions of the dynamic programs computed so far

    // returns the upper bound of a step of a search that has searched P[beg, end) and has the upper bound k_max
    uint64_t upper_bound(uint64_t beg, uint64_t end, uint64_t k_max) const
    {
        if (_costs == nullptr) return k_max;
        auto [it, inserted] = _max_errors.try_emplace(beg * (_match.size() + 1) + end, 0);
        if (inserted) it->second = _costs->max_errors(beg, end);
        return std::min(k_max, it->second);
    }

    /**
     * @brief appends position i to the strings with counters (see for_each_step): a match, a substitution, an insertion
     *        and a deletion each add their counter increments (and those of opening a gap), and strings whose counters
     *        exceed their bounds or that have more than k_max errors are dropped
     */
    void append_counted(std::vector<counted_t>& cnt, std::vector<counted_t>& next, uint64_t i, uint64_t k_max) const
    {
        next.clear();
        _work += num_operations * cnt.size();

        for (const counted_t& c : cnt) {
            auto add = [&](uint64_t counters, uint8_t errors, uint8_t last, double num) {
                if (num == 0 || errors > k_max || _costs->exceeds(counters)) return;
                next.emplace_back(counted_t{counters, errors, last, num / _sigma});
            };

            add(c.counters, c.errors, 0, c.num * _match[i]);
            add(c.counters + _costs->sub_count(i), c.errors + 1, 0, c.num * _sub[i]);
            add(c.counters + _costs->ins_count(i) + (c.last == 1 ? 0 : _costs->ins_open_count(i)), c.errors + 1, 1,
                c.num * _ins[i]);
            add(c.counters + _costs->del_count(i) + (c.last == 2 ? 0 : _costs->del_open_count(i)), c.errors + 1, 2,
                c.num * _del[i]);
        }

        std::sort(next.begin(), next.end(), [](const counted_t& a, const counted_t& b) { return a.key() < b.key(); });
        cnt.clear();

        for (const counted_t& c : next) {
            if (!cnt.empty() && cnt.back().key() == c.key()) cnt.back().num += c.num;
            else cnt.emplace_back(c);
        }
    }

    /**
     * @brief calls step(b) for t = 1, 2, ... with b = A(t) / sigma^t of the search S, for the parts after its first one
     *        in the partition with the bounds; with counters in the error model, the strings are counted per number of
     *        errors, counters and last operation (see append_counted)
     */
    template <typename step_fnc_t>
    void for_each_step(const search_t& S, const std::vector<uint64_t>& bounds, step_fnc_t step) const
    {
        if (_costs != nullptr && _costs->num_counters() > 0) { for_each_counted_step(S, bounds, step); return; }

        const std::vector<direction_t> dirs = apm_step_directions(S);
        const bool no_lower_bounds = _costs != nullptr && _costs->metric() == EDIT_DISTANCE;
        std::vector<double> cnt(_scheme.k + 1, 0.0); // [0..k] A(t) / sigma^t per number of errors
        cnt[0] = 1;
        uint64_t beg = bounds[S[0].part], end = bounds[S[0].part + 1];

        for (uint64_t j = 1; j < S.size(); j++) {
            const uint64_t part_beg = bounds[S[j].part], part_end = bounds[S[j].part + 1];
            beg = std::min(beg, part_beg);
            end = std::max(end, part_end);
            const uint64_t k_max = upper_bound(beg, end, S[j].k_max);
            for (uint64_t e = k_max + 1; e < cnt.size(); e++) cnt[e] = 0;

            _work += (part_end - part_beg) * (k_max + 1);

            for (uint64_t x = 0; x < part_end - part_beg; x++) {
                const uint64_t i = dirs[j] == RIGHT ? part_beg + x : part_end - 1 - x;
                for (uint64_t e = k_max; e > 0; e--) cnt[e] = (cnt[e] * _match[i] + cnt[e - 1] * _error[i]) / _sigma;
                cnt[0] = cnt[0] * _match[i] / _sigma;
                double b = 0;
                for (double c : cnt) b += c;
                step(b);
            }

            if (!no_lower_bounds) for (uint64_t e = 0; e < S[j].k_min; e++) cnt[e] = 0;
        }
    }

    // for_each_step with counters (see append_counted)
    template <typename step_fnc_t>
    void for_each_counted_step(const search_t& S, const std::vector<uint64_t>& bounds, step_fnc_t step) const
    {
        const std::vector<direction_t> dirs = apm_step_directions(S);
        const bool no_lower_bounds = _costs->metric() == EDIT_DISTANCE;
        std::vector<counted_t> cnt { counted_t{0, 0, 0, 1} }, next;
        uint64_t beg = bounds[S[0].part], end = bounds[S[0].part + 1];

        for (uint64_t j = 1; j < S.size(); j++) {
            const uint64_t part_beg = bounds[S[j].part], part_end = bounds[S[j].part + 1];
            beg = std::min(beg, part_beg);
            end = std::max(end, part_end);
            const uint64_t k_max = upper_bound(beg, end, S[j].k_max);

            for (uint64_t x = 0; x < part_end - part_beg; x++) {
                append_counted(cnt, next, dirs[j] == RIGHT ? part_beg + x : part_end - 1 - x, k_max);
                double b = 0;
                for (const counted_t& c : cnt) b += c.num;
                step(b);
            }

            if (!no_lower_bounds) std::erase_if(cnt, [&](const counted_t& c) { return c.errors < S[j].k_min; });
        }
    }

    // returns the predicted number of steps of the search S whose first part has occ occurrences
    double search_cost(const search_t& S, const std::vector<uint64_t>& bounds, uint64_t occ) const
    {
        const double x = double(occ) * _share;
        double power = 1; // sigma^t while it is at most x
        double c = 0;

        for_each_step(S, bounds, [&](double b) {
            if (power <= x) power *= _sigma;
            c += power <= x ? b * power : x * b;
        });

        return c;
    }

    // returns whether a position of P[beg, end) has a free symbol (then the part has no search context)
    bool has_free(uint64_t beg, uint64_t end) const
    {
        for (uint64_t i = beg; i < end; i++) if (_match[i] > 1) return true;
        return false;
    }

    /**
     * @brief returns the predicted number of steps of the anchored search of apm_model_plan: the anchor, then the parts
     *        to its right, then the parts to its left, with the bounds of the parts outside the anchor
     */
    double anchored_cost(const std::vector<uint64_t>& bounds) const
    {
        std::vector<uint64_t> b { 0, _anchor_beg, _anchor_end, _match.size() };
        for (uint64_t x : bounds) if (x < _anchor_beg || x > _anchor_end) b.emplace_back(x);
        std::sort(b.begin(), b.end());
        b.erase(std::unique(b.begin(), b.end()), b.end());
        const uint64_t q = std::find(b.begin(), b.end(), _anchor_beg) - b.begin();
        search_t T;
        for (uint64_t part = q; part + 1 < b.size(); part++) T.emplace_back(search_step_t{uint8_t(part), 0, _scheme.k});
        for (uint64_t part = q; part-- > 0;) T.emplace_back(search_step_t{uint8_t(part), 0, _scheme.k});
        return search_cost(T, b, _anchor_occ);
    }

public:
    /**
     * @brief prepares the prediction for a pattern of length m
     * @param scheme the search scheme
     * @param costs the error model of the pattern, or nullptr
     * @param sigma the alphabet size
     * @param match [0..m-1] the numbers of symbols at position i of the pattern without an error
     * @param sub [0..m-1] the numbers of substitutions at position i
     * @param ins [0..m-1] the numbers of insertions of position i
     * @param del [0..m-1] the numbers of deletions before position i
     * @param r the number of runs in the BWT of the text
     * @param n the length of the text
     * @param anchor_occ the number of occurrences of the anchor of the error model (see apm_costs_t::anchor)
     */
    apm_search_cost(const search_scheme_t& scheme, const apm_costs_t* costs, double sigma, std::vector<double> match,
        std::vector<double> sub, std::vector<double> ins, std::vector<double> del, uint64_t r, uint64_t n,
        uint64_t anchor_occ)
        : _scheme(scheme), _costs(costs), _sigma(sigma), _share(double(r) / double(n)), _n(n), _anchor_occ(anchor_occ),
          _match(std::move(match)),
          _sub(std::move(sub)), _ins(std::move(ins)), _del(std::move(del)), _searches_of(scheme.p)
    {
        const uint64_t m = _match.size();
        if (costs != nullptr) std::tie(_anchor_beg, _anchor_end) = costs->anchor();
        _error.resize(m);
        for (uint64_t i = 0; i < m; i++) _error[i] = _sub[i] + _ins[i] + _del[i];
        for (uint64_t s = 0; s < scheme.S.size(); s++) _searches_of[scheme.S[s][0].part].emplace_back(s);
    }

    // returns the cells and the transitions of the dynamic programs computed so far
    uint64_t work() const { return _work; }

    /**
     * @brief returns the predicted number of steps of the searches that start with part i, if part i has occ occurrences
     *        (with the bounds of the parts of the uniform partition)
     */
    double part_cost(uint64_t i, uint64_t occ) const
    {
        if (_sum_a.empty()) {
            const uint64_t m = _match.size();
            const double r = _share * double(_n);
            std::vector<uint64_t> bounds(_scheme.p + 1);
            for (uint64_t j = 0; j <= _scheme.p; j++) bounds[j] = (j * m) / _scheme.p;
            _sum_a.resize(_scheme.S.size());
            _sum_b.resize(_scheme.S.size());

            for (uint64_t s = 0; s < _scheme.S.size(); s++) {
                std::vector<double> b { 0 };
                for_each_step(_scheme.S[s], bounds, [&](double b_t) { b.emplace_back(b_t); });
                double power = 1;
                _sum_a[s] = { 0 };

                for (uint64_t t = 1; t < b.size() && power * _sigma <= r; t++) {
                    power *= _sigma;
                    _sum_a[s].emplace_back(_sum_a[s].back() + b[t] * power);
                }

                _sum_b[s].resize(b.size(), 0);
                for (uint64_t t = b.size() - 1; t > 0; t--) _sum_b[s][t - 1] = _sum_b[s][t] + b[t];
            }
        }

        const double x = double(occ) * _share;
        double c = 0;

        for (uint64_t s : _searches_of[i]) {
            // the largest t with sigma^t <= x
            const uint64_t t = x < _sigma ? 0 : std::min<uint64_t>(uint64_t(std::log(x) / std::log(_sigma)),
                std::min(_sum_a[s].size(), _sum_b[s].size()) - 1);
            c += _sum_a[s][t] + x * _sum_b[s][t];
        }

        return c;
    }

    /**
     * @brief returns the predicted number of steps of the searches that run: all searches of the scheme, or with an
     *        error model those that apm_model_plan keeps, or the anchored search if the plan runs it instead (if the
     *        anchor occurs at most as often as the first parts of the kept searches together, where a first part that
     *        admits errors or has a free symbol counts as occurring unboundedly often)
     * @param bounds [0..p] the bounds of the parts
     * @param occ [0..p-1] the numbers of occurrences of the parts
     */
    double cost(const std::vector<uint64_t>& bounds, const std::vector<uint64_t>& occ) const
    {
        if (_costs == nullptr) {
            double c = 0;
            for (const search_t& S : _scheme.S) c += search_cost(S, bounds, occ[S[0].part]);
            return c;
        }

        // the searches that the plan keeps, and whether it runs the anchored search instead
        std::vector<const search_t*> kept;
        double first_occ = 0;
        bool unbounded = false;

        for (const search_t& S : _scheme.S) {
            const uint64_t first = S[0].part;
            uint64_t beg = bounds[first], end = bounds[first + 1];
            bool feasible = true;

            for (const search_step_t& step : S) {
                beg = std::min(beg, bounds[step.part]);
                end = std::max(end, bounds[step.part + 1]);
                feasible &= step.k_min <= upper_bound(beg, end, step.k_max);
            }

            if (!feasible) continue;
            kept.emplace_back(&S);
            unbounded |= has_free(bounds[first], bounds[first + 1]) || upper_bound(bounds[first], bounds[first + 1], S[0].k_max) > 0;
            first_occ += double(occ[first]);
        }

        if (_anchor_beg < _anchor_end && (unbounded || double(_anchor_occ) <= first_occ)) return anchored_cost(bounds);
        double c = 0;

        for (const search_t* S : kept) {
            const uint64_t first = (*S)[0].part;
            c += search_cost(*S, bounds, std::min(occ[first], _n));
        }

        return c;
    }
};

#endif // MOVE_R_ALGORITHMS_MOVE_RB_APM_SEARCH_COST_HPP
