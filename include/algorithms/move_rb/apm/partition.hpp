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

#ifndef MOVE_R_ALGORITHMS_MOVE_RB_APM_PARTITION_HPP
#define MOVE_R_ALGORITHMS_MOVE_RB_APM_PARTITION_HPP

#include <cstdint>
#include <limits>
#include <ranges>
#include <span>
#include <vector>

#include <misc/apm.hpp>
#include <misc/utils.hpp>
#include <algorithms/move_rb/apm/model.hpp>

/**
 * @brief a partition of a pattern P into p parts, with the search context of each part: part i is
 *        P[beg(i), end(i)). The search context of a part is computed if no position of the part has a free symbol
 *        other than its own symbol (see apm_costs_t::has_free; always without an error model). The parts are of
 *        equal length, or chosen by the dynamic partitioning of Columba (Renders, Depuydt, Fostier).
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

        // the dynamic partitioning requires p >= 2 and m >= 2p
        if (p == 1 || params.partitioning == _partition_uniform || P.size() < 2 * p) {
            partition_uniform(P, p);
        } else {
            partition_dynamic(P, p);
        }
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
