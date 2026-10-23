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

#ifndef MOVE_R_ALGORITHMS_MOVE_RB_APM_EDIT_MODEL_TPP
#define MOVE_R_ALGORITHMS_MOVE_RB_APM_EDIT_MODEL_TPP

#include <algorithms/move_rb/apm/edit.hpp>

template <typename idx_t>
uint32_t apm_edit<idx_t>::ext_trie_child(uint32_t parent, uint8_t sym) const
{
    auto [it, added] = _ext_trie_child.try_emplace(ext_trie_key(parent, sym), uint32_t(_ext_trie_parent.size()));

    if (added) {
        _ext_trie_parent.emplace_back(parent);
        _ext_trie_symbol.emplace_back(sym);
    }

    return it->second;
}

template <typename idx_t>
auto apm_edit<idx_t>::entry_of(string_map_t& map, const node_t& node, uint64_t err, uint16_t cost, bool& inserted,
    uint64_t counters) const -> model_entry_t&
{
    auto [b, e] = node.ctx.forward_sa_interval();
    string_key_t key{b, e, node.ctx.pattern_length()};
    auto it = map.find(key);
    inserted = it == map.end();
    if (!inserted) return it->second;
    node_t copy = node;
    commit_copy(copy);
    copy.ctx.set_errors(err);
    copy.set_match_len(0);
    model_entry_t& entry = map.emplace(key, model_entry_t{std::move(copy), cost}).first->second;
    entry.counters = counters;
    return entry;
}

template <typename idx_t>
void apm_edit<idx_t>::add_member(model_entry_t& entry, const model_member_t& mem)
{
    auto it = std::find_if(entry.members.begin(), entry.members.end(),
        [&](const model_member_t& x) { return x.ext == mem.ext; });

    if (it == entry.members.end()) {
        entry.members.emplace_back(mem);
    } else {
        it->err = std::min(it->err, mem.err);
        it->cost = std::min(it->cost, mem.cost);
    }

    if (mem.err < entry.node.ctx.errors()) entry.node.ctx.set_errors(mem.err);
    entry.cost = std::min(entry.cost, mem.cost);
}

template <typename idx_t>
void apm_edit<idx_t>::insert_or_lower(string_map_t& map, const node_t& node, uint64_t err, uint16_t cost,
    uint64_t counters) const
{
    bool inserted;
    model_entry_t& entry = entry_of(map, node, err, cost, inserted, counters);
    if (inserted) return;

    if (!entry.members.empty()) {
        add_member(entry, model_member_t{0, uint16_t(err), cost});
        return;
    }

    if (err < entry.node.ctx.errors()) entry.node.ctx.set_errors(err);

    if (cost < entry.cost) {
        entry.cost = cost;
        entry.counters = counters;
    }
}

template <typename idx_t>
void apm_edit<idx_t>::insert_members(string_map_t& map, const node_t& node,
    const std::vector<model_member_t>& members) const
{
    bool inserted;
    model_entry_t& entry = entry_of(map, node, members.front().err, members.front().cost, inserted);

    // an entry without members stands for its own string: it becomes the member with the empty extension
    if (!inserted && entry.members.empty()) {
        entry.members.emplace_back(model_member_t{0, uint16_t(entry.node.ctx.errors()), entry.cost});
    }

    for (const model_member_t& mem : members) add_member(entry, mem);
}

template <typename idx_t>
template <cigar_mode_t mode, typename word_t>
void apm_edit<idx_t>::locate_at_text_ends_model(const inp_t& P, inp_t& strings, std::vector<uint64_t>& beg,
    std::vector<aprx_occ_t<mode>>& occ) const
{
    const apm_costs_t& costs = *_params.costs;
    const uint64_t m = P.size(), k = costs.max_errors();
    const uint64_t n = _index.n - 1;
    if (m == 0 || n == 0) return;
    auto [left, right] = text_end_windows(m, k);

    // adds the cheapest alignment of P against a prefix of a window W at position pos
    auto add = [&](uint64_t pos, const inp_t& W) {
        auto [cost, len] = costs.align_prefix(W);
        if (cost > costs.bound()) return;
        if constexpr (mode == CIGAR) occ.emplace_back(aprx_occ_t<CIGAR>{pos, len, cost, uint32_t(beg.size() - 1)});
        else                         occ.emplace_back(aprx_occ_t<NO_CIGAR>{pos, len, cost});
        strings.insert(strings.end(), W.begin(), W.end());
        beg.emplace_back(strings.size());
    };

    add(0, left);

    // the suffix of the last window with the fewest errors, from the matrix of the reversed pattern (a free symbol of
    // P[i] counts as equal to it)
    inp_t P_rev(P.rbegin(), P.rend()), right_rev(right.rbegin(), right.rend());
    directional_substring<inp_t> P_rev_col(P_rev, 0, m - 1, RIGHT);
    static thread_local edit_distance_matrix<word_t> mat;

    mat.set_input(P_rev_col, _index.sigma, _index.forward_index().map_int(), [&](uint64_t pos, auto report_rank) {
        for (uint8_t c : costs.free_symbols(m - 1 - pos)) {
            uint8_t r = _index.forward_index().map_int()[c];
            if (r != 0) report_rank(r);
        }
    });

    auto [err, len] = edit_dist_prefix(mat, right_rev, m, k);
    if (err <= k) add(n - len, inp_t(right.end() - len, right.end()));
}

template <typename idx_t>
template <cigar_mode_t mode, typename report_fnc_t>
std::vector<cigar_t<uint8_t>> apm_edit<idx_t>::locate_model(const inp_t& P, const search_scheme_t& scheme,
    report_fnc_t& report) const
{
    const apm_costs_t& costs = *_params.costs;
    std::vector<cigar_t<uint8_t>> cigars;
    if (P.empty()) return cigars;
    const bool uniform = costs.uniform_cost() != 0 && costs.max_errors() > 0 && scheme.k == costs.max_errors();
    std::vector<node_t> nodes;

    if (uniform) {
        node_set_t found = search(P, scheme);
        nodes.assign(found.begin(), found.end());
    } else {
        nodes = costs.max_errors() <= edit_distance_matrix<uint64_t>::max_supported_errors
            ? search_model<uint64_t>(P, scheme) : search_model<uint128_t>(P, scheme);
    }

    // the text strings of the contexts: strings[beg[x], beg[x + 1]) is the string of nodes[x]
    inp_t strings, match;
    std::vector<uint64_t> beg { 0 };

    for (const node_t& node : nodes) {
        decode_matched_string(P, node.last_piece_idx(), match);
        strings.insert(strings.end(), match.begin(), match.end());
        beg.emplace_back(strings.size());
    }

    // the cheapest valid alignment of P against a prefix of each string; ctxts contains each valid context, in CIGAR
    // mode with its index in nodes in place of a CIGAR index
    std::vector<std::pair<uint64_t, uint64_t>> aligned;
    costs.align_prefixes(strings, beg, aligned);
    std::vector<context_ref_t<mode>> ctxts;

    for (uint64_t x = 0; x < nodes.size(); x++) {
        auto [cost, len] = aligned[x];
        if (cost > costs.bound()) continue;
        nodes[x].ctx.set_errors(cost);
        nodes[x].set_match_len(len);
        if constexpr (mode == CIGAR) ctxts.emplace_back(&nodes[x], uint32_t(x));
        else                         ctxts.emplace_back(&nodes[x]);
    }

    std::vector<aprx_occ_t<mode>> occ;
    auto collect = [&](const aprx_occ_t<mode>& o) { occ.emplace_back(o); };
    report_contexts<mode>(ctxts, collect);
    const uint64_t located = occ.size();

    if (uniform) {
        if constexpr (can_decode_text_ends) {
            if (costs.max_errors() <= edit_distance_matrix<uint64_t>::max_supported_errors) {
                locate_at_text_ends_model<mode, uint64_t>(P, strings, beg, occ);
            } else {
                locate_at_text_ends_model<mode, uint128_t>(P, strings, beg, occ);
            }
        }
    }

    ips2ra::sort(occ.begin(), occ.end(), key_of_pos {});

    // one occurrence per position, the cheapest (a position occurs twice only if it is also found at a text end)
    if (occ.size() > located) {
        uint64_t write_idx = 0;

        for (uint64_t x = 0; x < occ.size(); x++) {
            if (write_idx > 0 && occ[write_idx - 1].pos() == occ[x].pos()) {
                if (occ[x].err() < occ[write_idx - 1].err()) occ[write_idx - 1] = occ[x];
            } else {
                occ[write_idx++] = occ[x];
            }
        }

        occ.resize(write_idx);
    }

    // every valid alignment starting at p is at most as cheap as the occurrence located at p, so the occurrences
    // kept by select_covering_occurrences cover all of them, unless contexts were dropped
    if (uniform || costs.prunes()) filter_edit_distance_occurrences<mode>(occ, costs.max_errors());
    else select_covering_occurrences<mode>(occ, 2 * costs.max_errors() + 1);

    if constexpr (mode == CIGAR) {
        std::vector<uint32_t> cigar_idx(beg.size() - 1, no_cigar);

        for (const aprx_occ_t<CIGAR>& o : occ) {
            uint32_t& c = cigar_idx[o.cig_idx];

            if (c == no_cigar) {
                const uint64_t x = o.cig_idx;
                costs.align_prefix(std::span(strings.data() + beg[x], beg[x + 1] - beg[x]), &cigars.emplace_back());
                c = uint32_t(cigars.size() - 1);
            }

            report({o.pos(), o.len(), o.err(), c});
        }
    } else {
        for (const aprx_occ_t<NO_CIGAR>& o : occ) report({o.pos(), o.len(), o.err()});
    }

    return cigars;
}

template <typename idx_t>
template <typename matrix_word_t>
auto apm_edit<idx_t>::search_model(const inp_t& P, const search_scheme_t& scheme) const -> std::vector<node_t>
{
    const apm_costs_t& costs = *_params.costs;
    const uint64_t m = P.size();
    if (scheme.p >= m) throw std::runtime_error("p >= m");
    _piece_owner = this; _matched_pieces.clear(); _free_pieces.clear();
    _max_diagonal_shift = costs.max_errors();
    _ext_trie_parent.assign(1, 0);
    _ext_trie_symbol.assign(1, 0);
    _ext_trie_child.clear();

    apm_model_plan<idx_t, LOCATE> plan(_index, P, scheme, _params);
    const auto& partition = plan.partition();
    const search_scheme_t& plan_scheme = plan.scheme();
    const uint64_t p = plan_scheme.p;

    std::vector<directional_substring<inp_t>> parts;
    parts.reserve(p);
    for (uint64_t part = 0; part < p; part++) parts.emplace_back(P, partition.beg(part), partition.end(part) - 1);

    // [0..2p-1] the matrix of part i read to the right (2i) and to the left (2i + 1) (see set_matrix_input)
    std::vector<edit_distance_matrix<matrix_word_t>> matrices(2 * p);

    string_map_t result, step_contexts, next_step_contexts;
    std::vector<model_entry_t*> order;

    // the dynamic program over the cost, if the model asks for it (see apm_costs_t::bounds_costs); if the model prunes
    // (see apm_costs_t::prunes), the last step and the steps that collect members prune their cells (see process_part)
    apm_cost_band cost_band;
    apm_cost_band* band = costs.bounds_costs() ? &cost_band : nullptr;
    const bool prune = band != nullptr && costs.prunes();

    for (const search_t& S : plan_scheme.S) {
        const std::vector<direction_t> dirs = apm_step_directions(S);
        for (uint64_t p_idx = 0; p_idx < p; p_idx++) parts[S[p_idx].part].set_direction(dirs[p_idx]);

        node_t root(_index.template empty_context<LOCATE>());
        root.ctx.set_errors(0);
        uint64_t p_idx = 0;

        // a first part without errors and free symbols starts at its search context
        if (_params.start_at_part_contexts && S[0].k_max == 0 && partition.has_context(S[0].part)) {
            if (!partition.occurs(S[0].part)) continue;
            const auto& part = parts[S[0].part];
            root.ctx = partition.context(S[0].part);
            root.ctx.set_errors(0);
            for (uint64_t i = 0; i < part.size(); i++) append_pattern_symbol(root, part.pos(i), part.direction());
            p_idx = 1;
        }

        step_contexts.clear();
        insert_or_lower(step_contexts, root, 0, 0);
        direction_t member_dir = NO_DIR;

        for (; p_idx < p && !step_contexts.empty(); p_idx++) {
            const auto& part = parts[S[p_idx].part];
            auto& matrix = matrices[2 * S[p_idx].part + (dirs[p_idx] == LEFT ? 1 : 0)];
            if (!matrix.is_initialized()) set_matrix_input(matrix, part);
            next_step_contexts.clear();

            // no later step extends the text beyond a part that ends the pattern in its direction
            bool outer = dirs[p_idx] == LEFT ? partition.beg(S[p_idx].part) == 0 : partition.end(S[p_idx].part) == m;
            const bool collect_as_members = outer && p_idx + 1 < p;
            if (collect_as_members) member_dir = dirs[p_idx];

            // the contexts in increasing order of their lengths; a search from a context continues the longer
            // contexts it reaches at most depth symbols after it (see process_part)
            order.clear();
            for (auto& entry : step_contexts) order.emplace_back(&entry.second);

            std::sort(order.begin(), order.end(), [](const model_entry_t* x, const model_entry_t* y) {
                return x->node.ctx.pattern_length() < y->node.ctx.pattern_length();
            });

            const uint64_t depth = std::min<uint64_t>(
                2 * edit_distance_matrix<matrix_word_t>::max_supported_errors - S[p_idx].k_max,
                order.back()->node.ctx.pattern_length() - order.front()->node.ctx.pattern_length());

            for (model_entry_t* entry : order) {
                if (entry->covered) continue;
                process_part<matrix_word_t>(*entry, part, partition.beg(S[p_idx].part),
                    partition.end(S[p_idx].part), dirs[p_idx], outer, collect_as_members,
                    prune && (p_idx + 1 == p || (collect_as_members && costs.num_counters() == 0)),
                    S[p_idx].k_min, S[p_idx].k_max, matrix, band, step_contexts, depth, next_step_contexts);
            }

            step_contexts.swap(next_step_contexts);
        }

        if (p_idx < p) continue;

        // the contexts of the last step and of their members
        inp_t ext;

        for (auto& [key, entry] : step_contexts) {
            if (entry.members.empty()) {
                insert_or_lower(result, entry.node, entry.node.ctx.errors(), entry.cost);
                continue;
            }

            for (const model_member_t& mem : entry.members) {
                // the symbols of the extension of the member, from the last one to the first one
                ext.clear();
                for (uint32_t x = mem.ext; x != 0; x = _ext_trie_parent[x]) ext.push_back(_ext_trie_symbol[x]);
                node_t node = entry.node;
                bool found = true;

                for (uint64_t x = ext.size(); x-- > 0 && found;) {
                    auto [ctx, ok] = node.ctx.extend(ext[x], member_dir);
                    found = ok;
                    node.ctx = ctx;
                    node.set_last_piece(append_literal(node.last_piece_idx(), ext[x], member_dir));
                }

                if (found) insert_or_lower(result, node, mem.err, mem.cost);
            }
        }
    }

    std::vector<node_t> nodes;
    nodes.reserve(result.size());
    for (auto& entry : result) nodes.emplace_back(std::move(entry.second.node));
    return nodes;
}

template <typename idx_t>
template <typename matrix_word_t>
void apm_edit<idx_t>::process_part(const model_entry_t& entry, const directional_substring<inp_t>& part,
    uint64_t beg, uint64_t end, direction_t dir, bool outer, bool collect_as_members, bool prune, uint64_t k_min,
    uint64_t k_max, edit_distance_matrix<matrix_word_t>& matrix, apm_cost_band* band, string_map_t& step_contexts,
    uint64_t depth, string_map_t& target) const
{
    const apm_costs_t& costs = *_params.costs;
    const uint64_t err_in = entry.node.ctx.errors();
    std::vector<node_t> stack;
    node_t start = entry.node;
    start.set_match_len(0);

    /* the members along the paths: member_lists[x] is a list of members, stack_member_list[y] the index of the list of
       the path of stack[y], or no_list if the path stands for the string of its node only (with the error count and
       the cost of its cell) */
    static constexpr uint32_t no_list = std::numeric_limits<uint32_t>::max();
    std::vector<std::vector<path_member_t>> member_lists;
    std::vector<uint32_t> stack_member_list;
    std::vector<model_member_t> cell_members;

    // appends the members of a context to a list, relative to the cell in which their alignments start
    auto append_members = [&](const model_entry_t& ctx, std::vector<path_member_t>& list) {
        const uint16_t err = ctx.node.ctx.errors(), cost = ctx.cost;
        if (ctx.members.empty()) list.emplace_back(path_member_t{0, 0, 0});
        for (const model_member_t& mem : ctx.members) {
            list.emplace_back(path_member_t{mem.ext, uint16_t(mem.err - err), uint16_t(mem.cost - cost)});
        }
    };

    uint32_t root_list = no_list;

    if (!entry.members.empty()) {
        append_members(entry, member_lists.emplace_back());
        root_list = 0;
    }

    // inserts the alignments of the members of a path that end in a cell with error count err, cost cost and counters
    auto insert_cell = [&](const node_t& node, uint32_t list, uint64_t err, uint64_t cost, uint64_t counters) {
        auto valid = [&](uint64_t e, uint64_t c) {
            return k_min <= e && e <= k_max && c <= costs.bound();
        };

        if (list == no_list) { if (valid(err, band != nullptr ? cost : 0)) insert_or_lower(target, node, err, cost, counters); return; }
        cell_members.clear();

        for (const path_member_t& mem : member_lists[list]) {
            const uint64_t mem_err = err + mem.err_offset, mem_cost = cost + mem.cost_offset;
            if (!valid(mem_err, band != nullptr ? mem_cost : 0)) continue;
            cell_members.emplace_back(model_member_t{mem.ext, uint16_t(mem_err), uint16_t(mem_cost)});
        }

        if (!cell_members.empty()) insert_members(target, node, cell_members);
    };

    // [0..stack.size()-1] 1 <=> the path of stack[y] has restarted (see the restart below)
    std::vector<uint8_t> stack_restarted;

    auto push_children = [&](node_t& node, uint32_t list, bool restarted) {
        auto ext_ctx = node.ctx.extend_phase();
        ext_ctx.template prepare_extend_all<NO_TERM>(dir);

        while (ext_ctx.can_extend()) {
            if (next_is_separator(ext_ctx)) { ext_ctx.extend_next(); continue; }
            stack.emplace_back(ext_ctx.extend_next());
            stack.back().set_match_len(node.match_len() + 1);
            copy_matched_string(stack.back(), node);
            stack_member_list.emplace_back(list);
            stack_restarted.emplace_back(restarted);
        }
    };

    // no further error: the part is matched by free symbols only, which cost what their substitutions cost (the
    // cost bounds the search only with the dynamic program over it)
    if (err_in == k_max) {
        if (err_in < k_min && root_list == no_list) return;
        std::vector<uint16_t> path_costs { entry.cost };
        std::vector<uint64_t> path_counters { entry.counters };
        stack.emplace_back(std::move(start));

        while (!stack.empty()) {
            node_t node = std::move(stack.back());
            uint16_t cost = path_costs.back();
            uint64_t counters = path_counters.back();
            stack.pop_back();
            path_costs.pop_back();
            path_counters.pop_back();
            if (node.match_len() > 0) append_last_symbol(node, part, dir);
            if (node.match_len() == part.size()) { insert_cell(node, root_list, err_in, cost, counters); continue; }
            uint64_t pos = part.pos(node.match_len());

            for (uint8_t c : costs.free_symbols(pos)) {
                uint64_t r = costs.rank(c);
                const bool accepted = costs.matches(pos, r);
                uint64_t next_cost = cost + (accepted ? 0 : costs.sub(pos, r));
                if (band != nullptr && next_cost > costs.bound()) continue;
                auto [next_ctx, ok] = node.ctx.extend(c, dir);
                if (!ok) continue;
                node_t child(next_ctx);
                child.set_match_len(node.match_len() + 1);
                copy_matched_string(child, node);
                stack.emplace_back(std::move(child));
                path_costs.emplace_back(uint16_t(std::min<uint64_t>(next_cost, apm_forbidden)));
                path_counters.emplace_back(counters + (accepted ? 0 : costs.sub_count(pos)));
            }
        }

        return;
    }

    matrix.init(k_max, {uint16_t(err_in)}, depth > 0 ? depth + k_max : 0);
    if (band != nullptr) band->init(costs, beg, end, dir, entry.cost, entry.counters, matrix.num_rows(), depth);

    // the cost of the cell of a row in which the whole part is aligned (0 if the cost is not bounded)
    auto final_cost = [&](uint64_t row) -> uint16_t { return band != nullptr ? band->final_cost(row) : 0; };

    // the context of the step that equals a context reached in a row in [1..depth], or nullptr
    auto context_of_step = [&](const node_t& node, uint64_t row) -> model_entry_t* {
        if (row > depth) return nullptr;
        auto [b, e] = node.ctx.forward_sa_interval();
        auto it = step_contexts.find(string_key_t{b, e, node.ctx.pattern_length()});
        return it == step_contexts.end() ? nullptr : &it->second;
    };

    /* with collect_as_members: the first row of the final column, the context of the path in it, and [0..num_rows-1]
       path_ext[row] the extension (see _ext_trie_parent) of the path from the first row of the final column to row */
    const uint64_t first_final_row = std::max<int64_t>(0,
        int64_t(matrix.num_rows()) - int64_t(matrix.last_col_size()));
    node_t first_final_node = start;
    std::vector<uint32_t> path_ext(matrix.num_rows(), 0);

    /* with prune: [0..num_rows-1] best_cost[row] the cost of the cheapest inserted cell of the final column on the
       path up to row (apm_forbidden if there is none), best_row[row] its row */
    const uint64_t max_shift = 2 * costs.max_errors() + 1;
    std::vector<uint16_t> best_cost(matrix.num_rows(), apm_forbidden);
    std::vector<uint64_t> best_row(matrix.num_rows(), 0);
    const uint16_t slack = costs.prune_slack();

    // whether an inserted cell above row on its path covers the cell of row with cost cost (see prune)
    auto is_covered = [&](uint64_t row, uint16_t cost) {
        return best_cost[row] != apm_forbidden && best_cost[row] + slack <= cost
            && (dir == RIGHT || row - best_row[row] <= max_shift);
    };

    // whether the path of a node with a list of members and the restart flag restarted is pruned (see prune)
    auto path_prunes = [&](uint32_t list, bool restarted) {
        return prune && list == no_list && !(restarted && costs.is_affine());
    };

    // inserts the alignments that end in the cell of a row in which the whole part is aligned
    auto insert_final = [&](const node_t& node, uint32_t list, bool restarted, uint64_t row) {
        const uint64_t err = matrix(row, part.size());
        const uint16_t cost = final_cost(row);
        const uint64_t counters = band != nullptr ? band->final_counters(row) : 0;
        if (cost == apm_forbidden || (path_prunes(list, restarted) && is_covered(row, cost))) return;

        if (path_prunes(list, restarted) && k_min <= err && err <= k_max && !costs.exceeds(counters)) {
            best_cost[row] = cost;
            best_row[row] = row;
        }

        if (!collect_as_members) { insert_cell(node, list, err, cost, counters); return; }
        assert(list == no_list);
        if (k_min > err || err > k_max) return;
        cell_members.assign(1, model_member_t{path_ext[row], uint16_t(err), cost});
        insert_members(target, first_final_node, cell_members);
    };

    // whether no deeper cell of the final column on the path can be inserted (see prune)
    auto is_pruned = [&](uint64_t row, uint32_t list, bool restarted) {
        return path_prunes(list, restarted) && row >= depth && best_cost[row] != apm_forbidden
            && best_cost[row] + slack <= band->min_partial_cost(row)
            && (dir == RIGHT || matrix.num_rows() - 1 - best_row[row] <= max_shift);
    };

    if (matrix.is_in_final_column(0)) insert_final(start, root_list, false, 0);
    push_children(start, root_list, false);

    while (!stack.empty()) {
        node_t node = std::move(stack.back());
        uint32_t list = stack_member_list.back();
        bool restarted = stack_restarted.back();
        stack.pop_back();
        stack_member_list.pop_back();
        stack_restarted.pop_back();
        append_last_symbol(node, part, dir);
        uint64_t row = node.match_len();
        if (row >= matrix.num_rows()) continue;
        uint8_t sym = node.ctx.last_added_symbol();

        if (collect_as_members && row == first_final_row) first_final_node = node;
        if (collect_as_members && row > first_final_row) path_ext[row] = ext_trie_child(path_ext[row - 1], sym);
        best_cost[row] = best_cost[row - 1];
        best_row[row] = best_row[row - 1];

        // a context of the step whose error count leaves its cells of this row in the band restarts the matrices;
        // its members join the members of the path
        model_entry_t* step_ctx = context_of_step(node, row);
        const uint64_t step_err = step_ctx != nullptr ? step_ctx->node.ctx.errors() : k_max + 1;
        const bool restart = step_err <= k_max && step_err + row >= err_in;
        bool alive = matrix.compute_row(row, sym);
        if (restart) alive = matrix.add_alignment_start(row, step_err);

        if (band != nullptr && (alive || restart)) {
            bool alive_band = band->compute_row(row, sym);
            if (restart) alive_band = band->restart(row, step_ctx->cost, step_ctx->counters);
            alive &= alive_band;
        }

        if (restart) {
            step_ctx->covered = true;
            restarted = true;

            if (list != no_list || !step_ctx->members.empty()) {
                std::vector<path_member_t> joined;
                if (list == no_list) joined.emplace_back(path_member_t{0, 0, 0});
                else joined = member_lists[list];
                append_members(*step_ctx, joined);
                list = member_lists.size();
                member_lists.emplace_back(std::move(joined));
            }
        }

        bool terminal = !alive;

        if (!terminal && matrix.is_in_final_column(row)) {
            insert_final(node, list, restarted, row);

            // with a bounded cost, a longer text string can have a cheaper alignment than the matrix shows; a context
            // of the step may continue in a longer one
            terminal = (!outer && band == nullptr && step_ctx == nullptr && matrix.only_vertical_gaps_left(row))
                || is_pruned(row, list, restarted);
        }

        if (!terminal && row + 1 < matrix.num_rows()) push_children(node, list, restarted);
    }
}

#endif // MOVE_R_ALGORITHMS_MOVE_RB_APM_EDIT_MODEL_TPP
