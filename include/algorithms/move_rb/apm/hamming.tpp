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

#ifndef MOVE_R_ALGORITHMS_MOVE_RB_APM_HAMMING_TPP
#define MOVE_R_ALGORITHMS_MOVE_RB_APM_HAMMING_TPP

#include <algorithms/move_rb/apm/hamming.hpp>

template <typename idx_t>
cigar_t<uint8_t> apm_hamming<idx_t>::build_cigar(std::vector<std::pair<uint64_t, uint8_t>> mismatches, uint64_t m)
{
    std::sort(mismatches.begin(), mismatches.end());
    cigar_t<uint8_t> cigar;
    uint64_t prev = 0;

    for (auto [p, sym] : mismatches) {
        if (p > prev) cigar.emplace_back(MATCH<uint8_t>(p - prev));
        cigar.emplace_back(MISMATCH(sym));
        prev = p + 1;
    }

    if (m > prev) cigar.emplace_back(MATCH<uint8_t>(m - prev));
    return cigar;
}

template <typename idx_t>
template <typename partition_t>
auto apm_hamming<idx_t>::search_positions(const search_t& S, const partition_t& partition, uint64_t m)
    -> std::vector<search_position_t>
{
    const uint64_t p = S.size();
    const std::vector<direction_t> dirs = apm_step_directions(S);
    std::vector<search_position_t> positions;
    positions.reserve(m);

    for (uint64_t step = 0; step < p; step++) {
        const uint64_t beg = partition.beg(S[step].part);
        const uint64_t end = partition.end(S[step].part);
        const uint64_t len = end - beg;
        if (step > 0) positions.back().num_left = len;

        for (uint64_t i = 0; i < len; i++) {
            const uint64_t pos = dirs[step] == RIGHT ? beg + i : end - 1 - i;
            positions.emplace_back(search_position_t{step, dirs[step], pos, len - 1 - i});
        }
    }

    return positions;
}

template <typename idx_t>
template <query_support_t query_support>
void apm_hamming<idx_t>::add_result(search_context_set_t<query_support>& result, search_context_t<query_support>& ctx,
    uint64_t m) const
{
    if constexpr (query_support == LOCATE) {
        located_context_t node{ctx, no_cigar};
        auto it = result.find(node);

        if (it != result.end()) {
            if (!(it->ctx.s.d() == RIGHT && ctx.s.d() == LEFT)) return;
            it = result.erase(it);
        }

        if (_computes_cigars) node.cigar_idx = add_cigar(m);
        result.emplace_hint(it, node);
    } else {
        result.emplace(ctx);
    }
}

template <typename idx_t>
uint64_t apm_hamming<idx_t>::count(const inp_t& P, const search_scheme_t& scheme) const
{
    uint64_t m = P.size();
    uint64_t k = scheme.k;

    if (_params.costs != nullptr) {
        uint64_t count = 0;
        for (const auto& ctx : search_model<COUNT>(P, scheme)) count += ctx.num_occ();
        return count;
    }

    if (k >= m) return _index.n - 1;
    if (k == 0) return _index.forward_index().count(P);

    uint64_t count = 0;

    auto result = search<COUNT>(P, scheme);
    for (const auto& ctx : result) count += ctx.num_occ();

    return count;
}

template <typename idx_t>
std::vector<uint64_t> apm_hamming<idx_t>::count_histogram(const inp_t& P, const search_scheme_t& scheme) const
{
    uint64_t m = P.size();
    uint64_t k = scheme.k;

    if (_params.costs != nullptr) {
        std::vector<uint64_t> histogram(_params.costs->bound() + 1, 0);
        for (const auto& ctx : search_model<COUNT>(P, scheme)) histogram[ctx.errors()] += ctx.num_occ();
        return histogram;
    }

    std::vector<uint64_t> histogram(k + 1, 0);

    // k >= m: every position matches; the total is counted with 0 errors
    if (k >= m) { histogram[0] = _index.n - 1; return histogram; }
    if (k == 0) { histogram[0] = _index.forward_index().count(P); return histogram; }

    // the occurrences of a search context have the same number of mismatches
    auto result = search<COUNT>(P, scheme);
    for (const auto& ctx : result) histogram[ctx.errors()] += ctx.num_occ();

    return histogram;
}

template <typename idx_t>
template <cigar_mode_t mode, typename report_fnc_t>
std::vector<cigar_t<uint8_t>> apm_hamming<idx_t>::locate(const inp_t& P, const search_scheme_t& scheme,
    report_fnc_t report) const
{
    uint64_t m = P.size();
    assert(m < (uint64_t(1) << aprx_occ_len_bits));
    uint64_t k = scheme.k;
    _computes_cigars = mode == CIGAR;
    if (k == 0 && _params.costs == nullptr) return this->template locate_exact<mode>(P, report);

    // in CIGAR mode, the search also fills _cigars
    auto result = _params.costs != nullptr ? search_model<LOCATE>(P, scheme) : search<LOCATE>(P, scheme);

    for (const auto& node : result) {
        if constexpr (mode == CIGAR) {
            node.ctx.locate_phase().locate([&](uint64_t occ){report({occ, node.ctx.pattern_length(), node.ctx.errors(), node.cigar_idx});});
        } else {
            node.ctx.locate_phase().locate([&](uint64_t occ){report({occ, node.ctx.pattern_length(), node.ctx.errors()});});
        }
    }

    if constexpr (mode == CIGAR) return std::move(_cigars);
    else return {};
}

template <typename idx_t>
template <query_support_t query_support>
auto apm_hamming<idx_t>::search(const inp_t& P, const search_scheme_t& scheme) const
    -> search_context_set_t<query_support>
{
    uint64_t m = P.size();
    uint64_t p = scheme.p;
    if (m == 0) return {_index.template empty_context<query_support>()};
    if (p >= m) throw std::runtime_error("p >= m");

    // an entry of the stack of the depth-first search
    struct search_state_t {
        search_context_t<query_support> ctx;
        uint32_t position_idx; // index in positions of the next position to match
        uint32_t added_pos; // the position of the mismatch this state added to its path, or no_mismatch (CIGAR mode)
        uint8_t num_errors; // number of mismatches so far (at most k)
        uint8_t num_path_mismatches; // number of mismatches on the path of this state (CIGAR mode)
        uint8_t added_sym; // the text symbol of the added mismatch (CIGAR mode)
    };

    search_context_set_t<query_support> result;
    auto empty_ctx = _index.template empty_context<query_support>();
    if constexpr (query_support == LOCATE) if (_computes_cigars) { _path_mismatches.clear(); _cigars.clear(); }
    auto add_ctx = [&](search_context_t<query_support>& ctx) { add_result<query_support>(result, ctx, m); };

    // partition P into p parts and compute the search context of each part
    apm_partition<idx_t, query_support> partition(_index, P, p, _params);

    for (const search_t& S : scheme.S) {
        std::vector<search_position_t> positions = search_positions(S, partition, m);
        std::vector<search_state_t> stack;
        stack.reserve(m);
        uint64_t len_first_part = partition.length(S[0].part);

        // a first part without errors starts at its search context (see apm_params_t::start_at_part_contexts); the search is
        // skipped if the part does not occur
        if (_params.start_at_part_contexts && S[0].k_max == 0 && len_first_part < m) {
            if (!partition.occurs(S[0].part)) continue;
            stack.emplace_back(search_state_t{partition.context(S[0].part),
                uint32_t(len_first_part), uint32_t(no_mismatch), 0, 0, uint8_t{}});
        } else {
            stack.emplace_back(search_state_t{empty_ctx, 0, uint32_t(no_mismatch), 0, 0, uint8_t{}});
        }

        while (!stack.empty()) {
            auto [ctx, position_idx, added_pos, num_errors, num_path_mismatches, added_sym] = stack.back(); stack.pop_back();

            if constexpr (query_support == LOCATE) if (_computes_cigars) {
                _path_mismatches.resize(num_path_mismatches - (added_pos != no_mismatch ? 1 : 0));
                if (added_pos != no_mismatch) _path_mismatches.emplace_back(added_pos, added_sym);
            }

            const search_position_t& cur = positions[position_idx];
            uint64_t next_position_idx = position_idx + 1;
            uint64_t next_step = next_position_idx == m ? p : positions[next_position_idx].step;

            // P[cur.pos] can be a mismatch: extend with every symbol
            if (num_errors < S[cur.step].k_max && (next_step == p || S[next_step].k_min <= num_errors + 1 + cur.num_left)) {
                auto ext_ctx = ctx.extend_phase();
                ext_ctx.template prepare_extend_all<NO_TERM>(cur.dir);

                while (ext_ctx.can_extend()) {
                    if (next_is_separator(ext_ctx)) { ext_ctx.extend_next(); continue; }
                    auto next_ctx = ext_ctx.extend_next();
                    bool is_mismatch = next_ctx.last_added_symbol() != P[cur.pos];
                    uint8_t sym = next_ctx.last_added_symbol();
                    uint64_t next_num_errors = num_errors + is_mismatch;

                    if (next_step == p) {
                        next_ctx.set_errors(next_num_errors);
                        if constexpr (query_support == LOCATE)
                            if (_computes_cigars && is_mismatch) _path_mismatches.emplace_back(cur.pos, sym);
                        add_ctx(next_ctx);
                        if constexpr (query_support == LOCATE)
                            if (_computes_cigars && is_mismatch) _path_mismatches.pop_back();
                    } else if (is_mismatch || S[next_step].k_min <= next_num_errors + cur.num_left) {
                        stack.emplace_back(search_state_t{next_ctx, uint32_t(next_position_idx),
                            uint32_t(is_mismatch ? cur.pos : no_mismatch), uint8_t(next_num_errors),
                            uint8_t(num_path_mismatches + (is_mismatch ? 1 : 0)), is_mismatch ? sym : uint8_t(0)});
                    }
                }
            } else {
                auto [next_ctx, occurs] = ctx.extend(P[cur.pos], cur.dir);

                if (occurs) {
                    if (next_step == p) {
                        next_ctx.set_errors(num_errors);
                        add_ctx(next_ctx);
                    } else {
                        stack.emplace_back(search_state_t{next_ctx, uint32_t(next_position_idx),
                            uint32_t(no_mismatch), num_errors, num_path_mismatches, uint8_t{}});
                    }
                }
            }
        }
    }

    return result;
}

template <typename idx_t>
template <query_support_t query_support>
auto apm_hamming<idx_t>::search_model(const inp_t& P, const search_scheme_t& scheme) const
    -> search_context_set_t<query_support>
{
    const apm_costs_t& costs = *_params.costs;
    uint64_t m = P.size();
    if (m == 0) return {_index.template empty_context<query_support>()};
    if (scheme.p >= m) throw std::runtime_error("p >= m");

    // an entry of the stack of the depth-first search
    struct search_state_t {
        search_context_t<query_support> ctx;
        uint64_t counters; // the packed counters (see apm_costs_t)
        uint32_t position_idx; // index in positions of the next position to match
        uint32_t added_pos; // the position at which this state differs from P, or no_mismatch (CIGAR mode)
        uint16_t cost; // the cost so far
        uint8_t num_errors; // number of errors so far (see apm_costs_t::tau)
        uint8_t num_path_mismatches; // number of positions on the path of this state that differ from P (CIGAR mode)
        uint8_t added_sym; // the text symbol at added_pos (CIGAR mode)
    };

    search_context_set_t<query_support> result;
    auto empty_ctx = _index.template empty_context<query_support>();
    if constexpr (query_support == LOCATE) if (_computes_cigars) { _path_mismatches.clear(); _cigars.clear(); }

    apm_model_plan<idx_t, query_support> plan(_index, P, scheme, _params);
    const auto& partition = plan.partition();
    const uint64_t p = plan.scheme().p;

    for (const search_t& S : plan.scheme().S) {
        std::vector<search_position_t> positions = search_positions(S, partition, m);
        std::vector<search_state_t> stack;
        stack.reserve(m);
        uint64_t len_first_part = partition.length(S[0].part);

        // a first part without errors and free symbols starts at its search context
        if (_params.start_at_part_contexts && S[0].k_max == 0 && partition.has_context(S[0].part) && len_first_part < m) {
            if (!partition.occurs(S[0].part)) continue;
            stack.emplace_back(search_state_t{partition.context(S[0].part), 0,
                uint32_t(len_first_part), uint32_t(no_mismatch), 0, 0, 0, 0});
        } else {
            stack.emplace_back(search_state_t{empty_ctx, 0, 0, uint32_t(no_mismatch), 0, 0, 0, 0});
        }

        while (!stack.empty()) {
            search_state_t st = std::move(stack.back());
            stack.pop_back();

            if constexpr (query_support == LOCATE) if (_computes_cigars) {
                _path_mismatches.resize(st.num_path_mismatches - (st.added_pos != no_mismatch ? 1 : 0));
                if (st.added_pos != no_mismatch) _path_mismatches.emplace_back(st.added_pos, st.added_sym);
            }

            const search_position_t& cur = positions[st.position_idx];
            uint64_t next_position_idx = st.position_idx + 1;
            uint64_t next_step = next_position_idx == m ? p : positions[next_position_idx].step;
            uint64_t next_k_min = next_step == p ? 0 : S[next_step].k_min;
            bool part_ends = next_step != cur.step;

            // extends the state by the text symbol last added to next_ctx
            auto add_child = [&](search_context_t<query_support>& next_ctx) {
                uint8_t sym = next_ctx.last_added_symbol();
                uint64_t r = costs.rank(sym);
                uint16_t cost = st.cost;
                uint64_t counters = st.counters;
                bool counted = false;

                if (!costs.matches(cur.pos, r)) {
                    uint16_t c = costs.sub(cur.pos, r);
                    if (c > costs.bound() - cost) return;
                    cost += c;
                    counters += costs.sub_count(cur.pos);
                    if (costs.exceeds(counters)) return;
                    counted = c >= costs.tau();
                }

                uint64_t next_num_errors = st.num_errors + counted;
                if (next_num_errors > S[cur.step].k_max || (part_ends && next_num_errors < S[cur.step].k_min)) return;
                bool differs = sym != P[cur.pos];

                if (next_step == p) {
                    next_ctx.set_errors(cost);
                    if constexpr (query_support == LOCATE) if (_computes_cigars && differs) _path_mismatches.emplace_back(cur.pos, sym);
                    add_result<query_support>(result, next_ctx, m);
                    if constexpr (query_support == LOCATE) if (_computes_cigars && differs) _path_mismatches.pop_back();
                } else if (next_k_min <= next_num_errors + cur.num_left) {
                    stack.emplace_back(search_state_t{next_ctx, counters, uint32_t(next_position_idx),
                        uint32_t(differs ? cur.pos : no_mismatch), cost, uint8_t(next_num_errors),
                        uint8_t(st.num_path_mismatches + (differs ? 1 : 0)), differs ? sym : uint8_t(0)});
                }
            };

            // P[cur.pos] can be an error: extend with every symbol; else with its free symbols only
            if (st.num_errors < S[cur.step].k_max && next_k_min <= st.num_errors + 1 + cur.num_left) {
                auto ext_ctx = st.ctx.extend_phase();
                ext_ctx.template prepare_extend_all<NO_TERM>(cur.dir);

                while (ext_ctx.can_extend()) {
                    if (next_is_separator(ext_ctx)) { ext_ctx.extend_next(); continue; }
                    auto next_ctx = ext_ctx.extend_next();
                    add_child(next_ctx);
                }
            } else {
                for (uint8_t c : costs.free_symbols(cur.pos)) {
                    auto [next_ctx, ok] = st.ctx.extend(c, cur.dir);
                    if (ok) add_child(next_ctx);
                }
            }
        }
    }

    return result;
}

#endif // MOVE_R_ALGORITHMS_MOVE_RB_APM_HAMMING_TPP
