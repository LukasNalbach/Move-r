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
 * @file step_prediction.tpp
 * @brief the methods of apm_step_prediction and apm_distinct_shares
 */

#ifndef MOVE_R_ALGORITHMS_MOVE_RB_APM_STEP_PREDICTION_TPP
#define MOVE_R_ALGORITHMS_MOVE_RB_APM_STEP_PREDICTION_TPP

#include <algorithms/move_rb/apm/step_prediction.hpp>

template <typename idx_t>
inline std::vector<double> apm_distinct_shares(const idx_t& index, uint64_t max_len)
{
    const double n = double(index.input_size());

    if constexpr (requires { index.has_distinct_substrings(); index.distinct_substrings(uint64_t(0)); }) {
        if (index.has_distinct_substrings()) {
            std::vector<double> shares(max_len + 1);
            for (uint64_t L = 0; L <= max_len; L++) shares[L] = double(index.distinct_substrings(L)) / n;
            return shares;
        }
    }

    return std::vector<double>(max_len + 1, double(index.num_bwt_runs()) / n);
}

inline auto apm_step_prediction::search_order(const search_t& S, const std::vector<uint64_t>& bounds, bounds_mode_t mode,
    bool from_root) const -> const search_order_t&
{
    std::vector<direction_t>& dirs = _dirs;
    apm_step_directions(S, dirs);
    const bool band = has_band(mode);
    search_order_t& order = _order;
    order.pos.clear();
    order.k_max.clear();
    order.k_min.clear();
    order.boundary.clear();
    order.common_prev.clear();
    order.trail.clear();
    order.k_min.emplace_back(0);
    order.trail.emplace_back(no_trail);
    uint64_t beg = bounds[S[0].part], end = bounds[S[0].part + 1]; // the positions searched so far are P[beg, end)
    direction_t prev_dir = dirs[0]; // the direction of the last position in order.pos

    for (uint64_t j = from_root ? 0 : 1; j < S.size(); j++) {
        const uint64_t part_beg = bounds[S[j].part], part_end = bounds[S[j].part + 1];
        beg = std::min(beg, part_beg);
        end = std::max(end, part_end);
        const uint64_t k_max = upper_bound(beg, end, S[j].k_max, mode.lowered);
        _work += part_end - part_beg;

        for (uint64_t y = 0; y < part_end - part_beg; y++) {
            const uint64_t i = dirs[j] == RIGHT ? part_beg + y : part_end - 1 - y;
            const bool continues = !order.pos.empty() && dirs[j] == prev_dir;
            order.common_prev.emplace_back(!continues ? -1.0 : _common_next[dirs[j] == RIGHT ? i - 1 : i]);
            order.pos.emplace_back(i);
            order.k_max.emplace_back(k_max);
            order.k_min.emplace_back(0);
            order.boundary.emplace_back(dirs[j] == RIGHT ? i : i + 1);
            order.trail.emplace_back(no_trail);
            prev_dir = dirs[j];
        }

        /* the deletions after the step, if the next step changes the direction or there is none, at a part of P at
           the first or the last position or with the costs computed in the search (see apm_edit::search_part) */
        const bool outer = S[j].part == 0 || uint64_t(S[j].part) + 2 == bounds.size();

        if ((j + 1 == S.size() || dirs[j + 1] != dirs[j]) && (outer || band)) {
            order.trail.back() = dirs[j] == RIGHT ? part_end : part_beg;
        }

        // the strings with fewer than S[j].k_min errors at the end of step j are not extended
        if (!mode.no_lower_bounds && S[j].k_min > 0) order.k_min.back() = S[j].k_min;
    }

    return order;
}

inline void apm_step_prediction::merge_groups(std::vector<counted_t>& next, std::vector<counted_t>& groups) const
{
    _work += next.size() * std::bit_width(next.size());
    std::sort(next.begin(), next.end(), [](const counted_t& a, const counted_t& b) { return a.key() < b.key(); });
    groups.clear();

    for (const counted_t& g : next) {
        if (!groups.empty() && groups.back().key() == g.key()) {
            groups.back().num += g.num;
            groups.back().counted += g.counted;
        } else {
            groups.emplace_back(g);
        }
    }
}

template <typename put_fnc_t>
inline void apm_step_prediction::extend_counted(const search_order_t& order, bool counts, bool band, put_fnc_t put) const
{
    const uint64_t U = order.pos.size();
    // whether the operations are those of the error model (else those of the edit distance, see the class description)
    const bool model_ops = !_edit || band;

    for (const counted_t& g : _groups) {
        const uint64_t u = g.u;

        /* appends the strings that continue g with an operation of weight w to u_next positions with at most k_max
           errors; count is the weight of them that are counted as strings at depth t */
        auto add = [&](uint64_t counters, uint64_t u_next, uint8_t errors, last_op_t last, double w, double count,
            uint64_t k_max) {
            _work++;
            if (w <= 0 || errors > k_max || (counts && !band && _costs->exceeds(counters))) return;
            put(counters, u_next, errors, last, g.num * w / _sigma, g.num * count / _sigma);
        };

        // a deletion after the step that ends after u positions, of any symbol
        const uint64_t tb = order.trail[u];

        if (tb != no_trail && g.last != INS && g.last != DEL && (!model_ops || _del_allowed[tb])) {
            add(counts ? g.counters + _costs->del_count(tb) + (g.last == TRAIL ? 0 : _costs->del_open_count(tb)) : 0,
                u, g.errors + 1, TRAIL, _sigma, _sigma, order.k_max[u - 1]);
        }

        if (g.last == TRAIL || u >= U) continue;
        const uint64_t i = order.pos[u];
        // whether pos[u] continues the direction of the position before it, and pos[u + 1] that of pos[u]
        const bool continues = order.common_prev[u] >= 0;
        const bool next_continues = u + 1 < U && order.common_prev[u + 1] >= 0;
        const bool after_ins = g.last == INS && continues;

        // a match; after an insertion, of the symbols that do not match the position before
        const double match = after_ins ? std::max(0.0, _match[i] - order.common_prev[u]) : _match[i];
        const double match_counted = after_ins ? 0.0 : g.last == DEL && next_continues ?
            std::max(0.0, _match[i] - order.common_prev[u + 1]) : match;
        const last_op_t match_last = after_ins ? INS_MATCH : g.last == DEL ? DEL_MATCH : DIAG;
        add(g.counters, u + 1, g.errors, match_last, match, match_counted, order.k_max[u]);

        // a substitution, not after an insertion and not after a deletion if pos[u + 1] continues the direction
        if (!after_ins && !(g.last == DEL && next_continues)) {
            const double sub = model_ops ? _sub[i] : _sigma - _match[i];
            add(counts ? g.counters + _costs->sub_count(i) : 0, u + 1, g.errors + 1, DIAG, sub, sub, order.k_max[u]);
        }

        // a deletion of a symbol that does not match pos[u], not after an insertion (followed by a match)
        const uint64_t b = order.boundary[u];

        if ((!model_ops || _del_allowed[b]) && !(continues && (g.last == INS || g.last == INS_MATCH))) {
            const double del = std::max(0.0, _sigma - _match[i]);
            add(counts ? g.counters + _costs->del_count(b) + (g.last == DEL ? 0 : _costs->del_open_count(b)) : 0,
                u, g.errors + 1, DEL, del, 0, order.k_max[u]);
        }
    }
}

inline void apm_step_prediction::start_edit_depth() const
{
    _stamp++;

    if (_slots.empty()) {
        const uint64_t cells = 2 * uint64_t(_scheme.k) + 1;
        _slots.resize(cells * (_scheme.k + 1) * num_last_ops, slot_t { 0, 0, 0, 0, no_slot });
        _cell_head.resize(cells);
        _cell_stamp.resize(cells, 0);
    }
}

inline void apm_step_prediction::put_slot(uint64_t cell, uint8_t errors, last_op_t last, uint64_t counters, double num,
    double counted) const
{
    const uint32_t idx = uint32_t((cell * (_scheme.k + 1) + errors) * num_last_ops + last);
    slot_t& slot = _slots[idx];

    if (slot.stamp == _stamp) {
        slot.counters = apm_costs_t::min_counters(slot.counters, counters);
        slot.num += num;
        slot.counted += counted;
        return;
    }

    slot = { counters, num, counted, _stamp, _cell_stamp[cell] == _stamp ? _cell_head[cell] : no_slot };
    _cell_head[cell] = idx;
    _cell_stamp[cell] = _stamp;
}

inline void apm_step_prediction::close_edit_depth(uint64_t t, const search_order_t& order, bool counts, bool band) const
{
    const uint64_t U = order.pos.size();
    const uint64_t k = _scheme.k;
    _groups.clear();

    for (uint64_t cell = 0; cell <= 2 * k; cell++) {
        if (_cell_stamp[cell] != _stamp) continue;
        const uint64_t u = cell + t - k;

        if (band && counts) {
            uint64_t min_cnt = _slots[_cell_head[cell]].counters;

            for (uint32_t idx = _cell_head[cell]; idx != no_slot; idx = _slots[idx].next) {
                min_cnt = apm_costs_t::min_counters(min_cnt, _slots[idx].counters);
            }

            if (_costs->exceeds(min_cnt)) continue;
            for (uint32_t idx = _cell_head[cell]; idx != no_slot; idx = _slots[idx].next) _slots[idx].counters = min_cnt;
        }

        for (uint32_t idx = _cell_head[cell]; idx != no_slot; idx = _slots[idx].next) {
            const slot_t& slot = _slots[idx];
            const last_op_t last = last_op_t(idx % num_last_ops);
            const uint8_t errors = uint8_t((idx / num_last_ops) % (k + 1));
            _groups.emplace_back(counted_t{slot.counters, uint32_t(u), errors, last, slot.num, slot.counted});
            _work++;
            if (u >= U || last == DEL || last == TRAIL || (last == DEL_MATCH && order.common_prev[u] >= 0)) continue;
            if (errors < order.k_min[u]) continue;
            const uint64_t i = order.pos[u];
            const uint64_t counters =
                counts ? slot.counters + _costs->ins_count(i) + (last == INS ? 0 : _costs->ins_open_count(i)) : 0;
            const double ins = band ? _ins[i] : 1.0;
            if (ins <= 0 || uint64_t(errors) + 1 > order.k_max[u]) continue;
            put_slot(cell + 1, errors + 1, INS, counters, slot.num * ins, 0);
        }
    }
}

template <typename step_fnc_t>
inline void apm_step_prediction::for_each_step(const search_t& S, const std::vector<uint64_t>& bounds, bounds_mode_t mode, bool from_root,
    step_fnc_t step) const
{
    if (_edit || (_costs != nullptr && _costs->num_counters() > 0)) {
        for_each_counted_step(S, bounds, mode, from_root, step);
        return;
    }

    std::vector<direction_t>& dirs = _dirs;
    apm_step_directions(S, dirs);
    // [0..k] b_by_errors[e] = A_e(t) / sigma^t, where A_e(t) is the number of strings at depth t with e errors
    std::vector<double>& b_by_errors = _b_by_errors;
    b_by_errors.assign(_scheme.k + 1, 0.0);
    b_by_errors[0] = 1;
    uint64_t beg = bounds[S[0].part], end = bounds[S[0].part + 1]; // the positions searched so far are P[beg, end)

    for (uint64_t j = from_root ? 0 : 1; j < S.size(); j++) {
        const uint64_t part_beg = bounds[S[j].part], part_end = bounds[S[j].part + 1];
        beg = std::min(beg, part_beg);
        end = std::max(end, part_end);
        const uint64_t k_max = upper_bound(beg, end, S[j].k_max, mode.lowered);
        // the strings with more than k_max errors are dropped
        for (uint64_t e = k_max + 1; e < b_by_errors.size(); e++) b_by_errors[e] = 0;

        _work += (part_end - part_beg) * (k_max + 1);

        for (uint64_t y = 0; y < part_end - part_beg; y++) {
            const uint64_t i = dirs[j] == RIGHT ? part_beg + y : part_end - 1 - y;
            // depth t from depth t - 1: a match keeps the number of errors, an error increases it by 1
            for (uint64_t e = k_max; e > 0; e--) b_by_errors[e] = (b_by_errors[e] * _match[i] + b_by_errors[e - 1] * _error[i]) / _sigma;
            b_by_errors[0] = b_by_errors[0] * _match[i] / _sigma;
            double b = 0;
            for (double b_e : b_by_errors) b += b_e;
            step(b);
        }

        // the strings with fewer than S[j].k_min errors are dropped
        if (!mode.no_lower_bounds) for (uint64_t e = 0; e < S[j].k_min; e++) b_by_errors[e] = 0;
    }
}

template <typename step_fnc_t>
inline void apm_step_prediction::for_each_counted_step(const search_t& S, const std::vector<uint64_t>& bounds, bounds_mode_t mode,
    bool from_root, step_fnc_t step) const
{
    const search_order_t& order = search_order(S, bounds, mode, from_root);
    const uint64_t U = order.pos.size();
    const uint64_t k = _scheme.k;
    const bool band = has_band(mode);
    // whether counters are checked: per string for the Hamming distance, per cell in the dynamic program of apm_edit
    const bool counts = _costs != nullptr && _costs->num_counters() > 0 && (!_edit || band);
    // the groups at depth 0
    _groups.clear();

    if (_edit) {
        start_edit_depth();
        put_slot(k, 0, DIAG, 0, 1, 0);
        close_edit_depth(0, order, counts, band);
    } else {
        _groups.emplace_back(counted_t{0, 0, 0, DIAG, 1, 0});
    }

    // [0..k] the strings at depth t with e errors (divided by sigma^t)
    std::vector<double>& b_by_errors = _b_by_errors;
    b_by_errors.assign(k + 1, 0.0);
    double extended_before = 1; // the strings at depth t - 1 that the search extends, divided by sigma^(t - 1)

    for (uint64_t t = 1; !_groups.empty(); t++) {
        // the groups at depth t: in the slots for the edit distance, else merged by sorting
        if (_edit) {
            start_edit_depth();
            auto put = [&](uint64_t counters, uint64_t u, uint8_t errors, last_op_t last, double num, double counted) {
                put_slot(u + k - t, errors, last, counters, num, counted);
            };
            extend_counted(order, counts, band, put);
            close_edit_depth(t, order, counts, band);
        } else {
            _next.clear();
            auto put = [&](uint64_t counters, uint64_t u, uint8_t errors, last_op_t last, double num, double counted) {
                _next.emplace_back(counted_t{counters, uint32_t(u), errors, last, num, counted});
            };
            extend_counted(order, counts, band, put);
            merge_groups(_next, _groups);
        }

        // the strings at depth t and those of them that the search extends
        double strings = 0, extended = 0;
        std::fill(b_by_errors.begin(), b_by_errors.end(), 0.0);

        for (const counted_t& g : _groups) {
            strings += g.counted;
            b_by_errors[g.errors] += g.counted;
            if (g.u < U || order.trail[g.u] != no_trail) extended += g.counted;
        }

        step(_edit ? extended_before : strings);
        extended_before = std::min(1.0, extended);
        // the strings with fewer errors than the lower bound at the end of a step are not extended
        std::erase_if(_groups, [&](const counted_t& g) { return g.last != TRAIL && g.errors < order.k_min[g.u]; });
        for (counted_t& g : _groups) g.counted = 0;
    }
}

inline double apm_step_prediction::search_steps(const search_t& S, const std::vector<uint64_t>& bounds, uint64_t occ, bounds_mode_t mode,
    bool from_root) const
{
    const uint64_t l = from_root ? 0 : bounds[S[0].part + 1] - bounds[S[0].part]; // the length of the start
    double power = 1; // sigma^t (at most the largest double)
    double steps = 0;
    uint64_t t = 0;

    for_each_step(S, bounds, mode, from_root, [&](double b) {
        t++;
        if (power < std::numeric_limits<double>::max() / _sigma) power *= _sigma;
        steps += occurring(b, t, power, double(occ), l);
    });

    return steps;
}

inline bool apm_step_prediction::is_infeasible(const search_t& S, const std::vector<uint64_t>& bounds, bool lowered) const
{
    uint64_t beg = bounds[S[0].part], end = bounds[S[0].part + 1];

    for (const search_step_t& step : S) {
        beg = std::min(beg, bounds[step.part]);
        end = std::max(end, bounds[step.part + 1]);
        if (step.k_min > upper_bound(beg, end, step.k_max, lowered)) return true;
    }

    return false;
}

inline bool apm_step_prediction::lowers_a_bound(const std::vector<uint64_t>& bounds) const
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

inline double apm_step_prediction::scheme_steps(const std::vector<uint64_t>& bounds, const std::vector<uint64_t>& occ,
    bounds_mode_t mode) const
{
    double steps = 0;

    for (const search_t& S : _scheme.S) {
        if (is_infeasible(S, bounds, mode.lowered)) continue;
        const bool from_root = starts_at_root(S, bounds, mode.lowered);
        steps += search_steps(S, bounds, from_root ? _n : std::min(occ[S[0].part], _n), mode, from_root);
    }

    return steps;
}

inline double apm_step_prediction::anchored_steps(const std::vector<uint64_t>& bounds) const
{
    // the bounds of the parts of the anchored search: the bounds of the anchor and those in bounds outside the anchor
    std::vector<uint64_t>& anchored_bounds = _anchored_bounds;
    anchored_bounds.assign({ 0, _anchor_beg, _anchor_end, _match.size() });
    for (uint64_t b : bounds) if (b < _anchor_beg || b > _anchor_end) anchored_bounds.emplace_back(b);
    std::sort(anchored_bounds.begin(), anchored_bounds.end());
    anchored_bounds.erase(std::unique(anchored_bounds.begin(), anchored_bounds.end()), anchored_bounds.end());
    // the index of the part that is the anchor
    const uint64_t anchor_part = std::find(anchored_bounds.begin(), anchored_bounds.end(), _anchor_beg) - anchored_bounds.begin();
    search_t& search = _anchored_search;
    search.clear();

    for (uint64_t part = anchor_part; part + 1 < anchored_bounds.size(); part++) {
        search.emplace_back(search_step_t{uint8_t(part), 0, _scheme.k});
    }

    for (uint64_t part = anchor_part; part-- > 0;) search.emplace_back(search_step_t{uint8_t(part), 0, _scheme.k});
    const bool from_root = starts_at_root(search, anchored_bounds, _may_lower);
    return search_steps(search, anchored_bounds, from_root ? _n : _anchor_occ, {_may_lower, true}, from_root);
}

template <typename pat_t>
inline apm_step_prediction::apm_step_prediction(const search_scheme_t& scheme, const apm_costs_t* costs,
    distance_metric_t metric, const pat_t& P, uint64_t sigma, const std::vector<double>& sym_p,
    std::vector<double> distinct_shares, uint64_t n, uint64_t anchor_occ, const apm_params_t& params)
    : _scheme(scheme), _costs(costs), _sigma(double(sigma)), _edit(metric == EDIT_DISTANCE),
      _distinct_share(std::move(distinct_shares)), _n(n),
      _anchor_occ(anchor_occ), _may_lower(params.lower_upper_bounds), _from_contexts(params.start_at_part_contexts),
      _free(P.size(), 0), _match(P.size(), 0), _sub(P.size(), 0), _ins(P.size(), 0), _del(P.size(), 0),
      _error(P.size()), _del_allowed(P.size() + 1, 0), _common_next(P.size() > 0 ? P.size() - 1 : 0, 0),
      _searches_of(scheme.p)
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

    for (uint64_t b = 0; b <= m; b++) _del_allowed[b] = edit && (costs == nullptr || costs->del(b) <= costs->bound());

    for (uint64_t i = 0; i + 1 < m; i++) {
        if (costs == nullptr) {
            _common_next[i] = P[i] == P[i + 1] ? _match[i] : 0;
        } else {
            for (uint64_t rank = 0; rank < costs->num_ranks(); rank++) {
                if (costs->is_free(i, rank) && costs->is_free(i + 1, rank)) _common_next[i] += rank_weight[rank];
            }
        }
    }

    if (costs != nullptr && params.anchored_search) std::tie(_anchor_beg, _anchor_end) = costs->anchor();
    for (uint64_t s = 0; s < scheme.S.size(); s++) _searches_of[scheme.S[s][0].part].emplace_back(s);
    for (const search_t& S : scheme.S) for (const search_step_t& step : S) _scheme_lower_bounds |= step.k_min > 0;
}

inline double apm_step_prediction::part_steps(uint64_t i, uint64_t occ) const
{
    const uint64_t m = _match.size();

    if (_part_strings.empty()) {
        const bounds_mode_t mode { _may_lower, _costs != nullptr && _costs->metric() == EDIT_DISTANCE };
        // [0..p] the bounds of the parts of the uniform partition
        std::vector<uint64_t> bounds(_scheme.p + 1);
        for (uint64_t j = 0; j <= _scheme.p; j++) bounds[j] = (j * m) / _scheme.p;
        _part_strings.resize(_scheme.S.size());

        for (uint64_t s = 0; s < _scheme.S.size(); s++) {
            for_each_step(_scheme.S[s], bounds, mode, false, [&](double b) {
                _part_strings[s].emplace_back(b);
            });
        }
    }

    // the length of part i in the uniform partition
    const uint64_t l = ((i + 1) * m) / _scheme.p - (i * m) / _scheme.p;
    double steps = 0;

    for (uint64_t s : _searches_of[i]) {
        double power = 1; // sigma^t (at most the largest double)

        for (uint64_t t = 1; t <= _part_strings[s].size(); t++) {
            if (power < std::numeric_limits<double>::max() / _sigma) power *= _sigma;
            steps += occurring(_part_strings[s][t - 1], t, power, double(occ), l);
        }
    }

    return steps;
}

inline auto apm_step_prediction::choose(const std::vector<uint64_t>& bounds, const std::vector<uint64_t>& occ) const -> choice_t
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

#endif // MOVE_R_ALGORITHMS_MOVE_RB_APM_STEP_PREDICTION_TPP
