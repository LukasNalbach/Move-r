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

#pragma once

#include <ips4o.hpp>
#include <move_rb/move_rb.hpp>

/**
 * @brief approximate pattern matching (APM) helper for a move_rb index; holds a reference to the
 *        index and implements the edit-distance locate algorithm outside the
 *        move_rb class
 * @tparam support type of locate support (_locate_move or _locate_rlzsa)
 * @tparam sym_t value type (default: char for strings)
 * @tparam pos_t index integer type (use uint32_t if input size < UINT_MAX, else uint64_t)
 */
template <move_r_support support, typename sym_t, typename pos_t>
class move_rb_apm_edit {
protected:
    using move_rb_t = move_rb<support, sym_t, pos_t>;
    const move_rb_t& index; // the index to perform approximate pattern matching on

public:
    /**
     * @brief constructs an approximate-pattern-matching helper for the index
     * @param index the index to perform approximate pattern matching on
     */
    move_rb_apm_edit(const move_rb_t& index) : index(index) {}

protected:
    using inp_t = typename move_rb_t::inp_t;
    template <move_rb_query_support_t query_support>
    using search_context_t = typename move_rb_t::template search_context_t<query_support>;
    template <move_rb_query_support_t query_support>
    using search_context_set_t = tsl::sparse_set<search_context_t<query_support>, typename search_context_t<query_support>::hash>;

    template <typename matrix_word_t>
    using data_struct_ref_t = std::tuple<
        const edit_dist_search&,
        const std::vector<sub_string<pos_t>>&,
        std::vector<std::vector<search_context_t<LOCATE>>>&,
        std::vector<edit_distance_matrix<matrix_word_t>>&,
        std::vector<edit_distance_matrix<matrix_word_t>>&,
        search_context_set_t<LOCATE>&
    >;

    struct interval_t {pos_t beg; pos_t end; pos_t err; pos_t len;};

    class cluster_t
    {
      protected:
        std::vector<uint16_t> dists;
        std::vector<search_context_t<LOCATE>> nodes;

        pos_t last_cell;
        uint16_t k_max;
        pos_t start_depth;
        pos_t shift;

      public:
        cluster_t(pos_t size, pos_t k_max, pos_t start_depth, pos_t shift)
            : dists(size, k_max + 1), nodes(size), last_cell(-1), k_max(k_max),
            start_depth(start_depth), shift(shift) {}

        void set_cell(pos_t idx, search_context_t<LOCATE> node, uint16_t dist)
        {
            nodes[idx] = std::move(node);
            nodes[idx].set_reported(false);
            dists[idx] = dist;
            last_cell = idx;
        }

        pos_t size() const
        {
            return dists.size();
        }

        std::vector<search_context_t<LOCATE>> report_centers_at_end()
        {
            std::vector<search_context_t<LOCATE>> centers;
            centers.reserve(last_cell + 1);

            for (pos_t i = 0; i <= last_cell; i++) {
                if (!nodes[i].reported() &&
                    dists[i] <= k_max &&
                    (i == 0 || dists[i] <= dists[i - 1]) &&
                    (i == last_cell || dists[i] <= dists[i + 1] || nodes[i].num_occ() > nodes[i + 1].num_occ())
                ) {
                    nodes[i].set_reported(true);
                    search_context_t<LOCATE> ctx = nodes[i];
                    ctx.set_errors(dists[i]);
                    ctx.set_depth(ctx.depth() + start_depth);
                    ctx.set_shift(shift);
                    centers.emplace_back(ctx);
                }
            }

            return centers;
        }

        std::vector<search_context_t<LOCATE>> plateau_minima(direction_t dir)
        {
            std::vector<search_context_t<LOCATE>> result;

            for (pos_t a = 0; a <= last_cell;) {
                pos_t b = a;
                while (b + 1 <= last_cell && dists[b + 1] == dists[a]) b++;

                if (dists[a] <= k_max && !nodes[b].reported()) {
                    nodes[b].set_reported(true);
                    search_context_t<LOCATE> ctx = nodes[b];
                    ctx.set_errors(dists[a]);
                    ctx.set_depth(ctx.depth() + start_depth - (b - a));
                    ctx.set_shift(((dir == LEFT) ? (b - a) : 0) + shift);
                    result.emplace_back(ctx);
                }

                a = b + 1;
            }

            return result;
        }

        std::tuple<search_context_t<LOCATE>, bool> cluster_center(uint16_t k_min, std::vector<search_context_t<LOCATE>>& desc, std::vector<uint16_t>& init_dists)
        {
            desc.reserve(dists.size());
            dists.reserve(dists.size());
            search_context_t<LOCATE> ctx;
            pos_t i = last_cell + 1;

            for (pos_t j = 0; j <= last_cell; j++) {
                if (dists[j] >= k_min && dists[j] <= k_max) {
                    i = j;
                    break;
                }
            }

            if (i <= last_cell) {
                bool res = false;

                if (!nodes[i].reported()) {
                    ctx = nodes[i];
                    ctx.set_errors(dists[i]);
                    ctx.set_depth(ctx.depth() + start_depth);
                    ctx.set_shift(shift);
                    res = true;
                }
                
                init_dists.emplace_back(dists[i]);

                for (pos_t j = i + 1; j <= last_cell; j++) {
                    desc.emplace_back(nodes[j]);
                    init_dists.emplace_back(dists[j]);
                }

                for (pos_t k = 1; k < init_dists.size(); k++) {
                    if (init_dists[k] < k_min && init_dists[k] <= init_dists[k - 1] &&
                        (k == init_dists.size() - 1 || init_dists[k] <= init_dists[k + 1])
                    ) {
                        pos_t high = 0;
                        pos_t low = init_dists.size() - 1;

                        for (pos_t l = k; l-- > 0;) {
                            if (init_dists[l] != init_dists[l + 1] + 1) {
                                high = l + 1;
                                break;
                            }
                        }
                        
                        for (pos_t l = k + 1; l < init_dists.size(); l++) {
                            if (init_dists[l] != init_dists[l - 1] + 1) {
                                low = l - 1;
                                break;
                            }
                        }

                        if (high != 0 && low != init_dists.size() - 1) {
                            pos_t lC = low;
                            pos_t hC = high;
                            bool highest = true;

                            while (lC > hC) {
                                if (highest) {
                                    init_dists[hC] = std::min(k_max + 1, init_dists[hC - 1] + 1);
                                    hC++;
                                } else {
                                    init_dists[lC] = std::min(k_max + 1, init_dists[lC + 1] + 1);
                                    lC--;
                                }

                                highest = !highest;
                            }

                            if (lC == hC) {
                                init_dists[lC] = std::min(init_dists[lC + 1] + 1, init_dists[lC - 1] + 1);
                            }
                        } else if (high == 0 && low != init_dists.size() - 1) {
                            for (pos_t l = low; l-- > 0;) {
                                init_dists[l] = init_dists[l + 1] + 1;
                            }
                        } else if (high != 0 && low == init_dists.size() - 1) {
                            for (pos_t l = high; l < init_dists.size(); l++) {
                                init_dists[l] = init_dists[l - 1] + 1;
                            }
                        }
                    }
                }

                return {ctx, res};
            }

            return {ctx, false};
        }
    };

public:
    /**
     * @brief computes the maximal sub-intervals of the union of the given (sorted) SA-intervals s.t. all
     *        occurrences in a sub-interval share the same (minimum) error count and length
     * @param ctxts_sorted the LOCATE search contexts sorted by their forward SA-interval starting position
     * @return the effective intervals, each storing its starting/ending position, error count and length
     */
    static auto effective_intervals(const std::vector<const search_context_t<LOCATE>*>& ctxts_sorted) -> std::vector<interval_t>
    {
        std::vector<interval_t> iv_stack;
        std::vector<interval_t> ivs;

        iv_stack.reserve(16);
        ivs.reserve(ctxts_sorted.size());

        auto append_segment = [&](interval_t iv) {
            auto [b, e, err, len] = iv;
            if (b > e) return;

            if (!ivs.empty()) {
                interval_t& seg_lst = ivs.back();
                const auto& [b_lst, e_lst, err_lst, len_lst] = seg_lst;

                if (err_lst == err && len_lst == len && e_lst + 1 == b) {
                    seg_lst.end = e;
                    return;
                }
            }

            ivs.emplace_back(interval_t{b, e, err, len});
        };

        for (const auto* ctx : ctxts_sorted) {
            auto [beg, end] = ctx->forward_sa_interval();
            pos_t err = ctx->errors();
            pos_t len = ctx->depth();

            while (!iv_stack.empty() && iv_stack.back().end < beg) {
                interval_t iv = iv_stack.back(); iv_stack.pop_back();
                append_segment(iv);

                if (!iv_stack.empty()) {
                    iv_stack.back().beg = iv.end + 1;
                }
            }

            if (iv_stack.empty()) {
                iv_stack.emplace_back(interval_t{beg, end, err, len});
            } else {
                interval_t& iv_top = iv_stack.back();
                append_segment({iv_top.beg, beg - 1, iv_top.err, iv_top.len});
                iv_top.beg = beg;

                if (iv_top.err < err || (iv_top.err == err && iv_top.len <= len)) {
                    iv_stack.emplace_back(interval_t{beg, end, iv_top.err, iv_top.len});
                } else {
                    iv_stack.emplace_back(interval_t{beg, end, err, len});
                }
            }
        }

        while (!iv_stack.empty()) {
            interval_t iv = iv_stack.back(); iv_stack.pop_back();
            append_segment(iv);

            if (!iv_stack.empty()) {
                iv_stack.back().beg = iv.end + 1;
            }
        }

        return std::move(ivs);
    }

    /**
     * @brief locates a pattern with at most k errors (w.r.t. edit distance)
     * @tparam report_fnc_t type of the function report
     * @param P the pattern to search
     * @param scheme the search scheme to use (provides k)
     * @param report function that is called with every occurrence (occ, len, err) of P in T with at most k errors (w.r.t. edit distance)
     */
    template <typename report_fnc_t>
    void locate(const inp_t& P, const search_scheme_t& scheme, report_fnc_t report) const requires(move_rb_t::supports_locate)
    {
        pos_t k = scheme.k;
        pos_t m = P.size();

        if (k == 0) {
            index.forward_index().locate(P, [&](pos_t occ){report({occ, m, 0});});
            return;
        }

        auto ctxts_set = search_entry(P, scheme);
        std::vector<const search_context_t<LOCATE>*> ctxts_sorted;
        ctxts_sorted.reserve(ctxts_set.size());
        for (const auto& ctx : ctxts_set) ctxts_sorted.emplace_back(&ctx);
        ips4o::sort(ctxts_sorted.begin(), ctxts_sorted.end(), [](auto x, auto y){return *x < *y;});
        auto ivs = effective_intervals(ctxts_sorted);
        pos_t iv_idx = 0;

        for (auto it = ctxts_sorted.begin(); it != ctxts_sorted.end();) {
            const auto& ctx = **it;
            auto [beg, end] = ctx.forward_sa_interval();
            auto loc_ctx = ctx.locate_phase();
            pos_t cntr = loc_ctx.compute_center(index, ctx);
            pos_t iv_idx_cntr = exp_search_max_leq<pos_t, RIGHT>(cntr, iv_idx, ivs.size() - 1, [&](pos_t x) {return ivs[x].beg;});
            iv_idx = iv_idx_cntr;

            loc_ctx.locate(index, ctx, [&](pos_t pos, pos_t occ){
                if (pos <= cntr) {
                    if (pos < ivs[iv_idx].beg) [[unlikely]] iv_idx--;
                } else {
                    if (iv_idx < iv_idx_cntr) [[unlikely]] iv_idx = iv_idx_cntr;
                    if (ivs[iv_idx].end < pos) [[unlikely]] iv_idx++;
                }

                report({occ, ivs[iv_idx].len, ivs[iv_idx].err});
            });

            iv_idx = std::max<pos_t>(iv_idx, iv_idx_cntr);
            do {it++;} while (it != ctxts_sorted.end() && std::get<0>((*it)->forward_sa_interval()) <= end);
        }
    }

    /**
     * @brief executes a given search scheme for a pattern (w.r.t. edit distance), choosing the edit-distance
     *        matrix word type based on the maximum number of errors k
     * @param P the pattern to search
     * @param scheme the search scheme to use
     * @return all search contexts of P in T with at most k errors (assuming the search scheme covers all error configurations)
     */
    auto search_entry(const inp_t& P, const search_scheme_t& scheme) const
        -> search_context_set_t<LOCATE> requires(move_rb_t::supports_locate)
    {
        static constexpr uint16_t k_limit_64 = edit_distance_matrix<uint64_t>::k_limit; // 10
        static constexpr uint16_t k_limit_128 = edit_distance_matrix<uint128_t>::k_limit; // 20

        if (scheme.k <= k_limit_64) {
            return search<uint64_t>(P, scheme);
        } else if (scheme.k <= k_limit_128) {
            return search<uint128_t>(P, scheme);
        } else {
            throw std::runtime_error("cannot search patterns with more than " + std::to_string(k_limit_128) + " errors.");
        }
    }

    /**
     * @brief executes a given search scheme for a pattern (w.r.t. edit distance)
     * @tparam matrix_word_t word type of the edit-distance matrix
     * @param P the pattern to search
     * @param scheme the search scheme to use
     * @return all search contexts of P in T with at most k errors (assuming the search scheme covers all error configurations)
     */
    template <typename matrix_word_t>
    auto search(const inp_t& P, const search_scheme_t& scheme) const
        -> search_context_set_t<LOCATE> requires(move_rb_t::supports_locate)
    {
        pos_t m = P.size();
        pos_t p = scheme.p;
        if (m == 0) return {index.template empty_context<LOCATE>()};
        if (p >= m) throw std::runtime_error("p >= m");

        search_context_set_t<LOCATE> ctxts_set;
        std::vector<edit_dist_search> searches;
        std::vector<sub_string<pos_t>> parts;
        searches.reserve(scheme.S.size());
        parts.reserve(p);

        for (pos_t p_idx = 0; p_idx < p; p_idx++) {
            pos_t part = scheme.S[0][p_idx].part;
            pos_t beg = (part * m) / p;
            pos_t end = (part == p - 1 ? m : (((part + 1) * m) / p)) - 1;
            parts.emplace_back(sub_string<pos_t>(P, beg, end));
        }
    
        for (pos_t s_idx = 0; s_idx < scheme.S.size(); s_idx++) {
            searches.emplace_back(edit_dist_search(scheme, s_idx));
        }

        std::vector<std::vector<search_context_t<LOCATE>>> stacks;
        stacks.resize(p);
    
        std::vector<edit_distance_matrix<matrix_word_t>> matrices_fwd, matrices_bwd;
        matrices_fwd.resize(p);
        matrices_bwd.resize(p);

        for (const edit_dist_search& search : searches) {
            data_struct_ref_t<matrix_word_t> ds_ref(search, parts, stacks, matrices_fwd, matrices_bwd, ctxts_set);

            for (pos_t p_idx = 0; p_idx < p; p_idx++) {
                parts[search.part(p_idx)].set_direction(search.part_dir(p_idx));
            }
        
            auto ctx = index.template empty_context<LOCATE>();
            pos_t p_idx = 0;

            if (search.upper_bound(0) == 0) {
                bool match = true;

                for (; match && p_idx < parts.size(); p_idx++) {
                    if (search.upper_bound(p_idx) > 0) break;
                    const auto& part = parts[search.part(p_idx)];

                    for (pos_t i = 0; match && i < part.size(); i++) {
                        auto [ctx_nxt, match_nxt] = ctx.extend(index, part[i], part.direction());
                        if (match_nxt) ctx = ctx_nxt;
                        match = match_nxt;
                    }
                }

                if (!match) continue;
            }

            rec_search<matrix_word_t>(ctx, p_idx, ds_ref, {}, {}, {}, {});
        }

        return ctxts_set;
    }

    /**
     * @brief recursively processes the p_idx-th part of the search by extending the current context in the
     *        part's direction and branching/bounding on the edit-distance matrix
     * @tparam matrix_word_t word type of the edit-distance matrix
     * @param ctx the current search context
     * @param p_idx index of the part to process
     * @param ds_ref references to the data structures shared during the search
     * @param desc_prev_dir descendant contexts from the previous part in the same direction
     * @param desc_diff_dir descendant contexts from the previous part in the opposite direction
     * @param dists_prev_dir edit distances corresponding to desc_prev_dir
     * @param dists_diff_dir edit distances corresponding to desc_diff_dir
     */
    template <typename matrix_word_t>
    void rec_search(
        search_context_t<LOCATE>& ctx, pos_t p_idx,
        data_struct_ref_t<matrix_word_t>& ds_ref,
        const std::vector<search_context_t<LOCATE>>& desc_prev_dir,
        const std::vector<search_context_t<LOCATE>>& desc_diff_dir,
        const std::vector<uint16_t>& dists_prev_dir,
        const std::vector<uint16_t>& dists_diff_dir
    ) const requires(move_rb_t::supports_locate)
    {
        auto& [search, parts, stacks, matrices_fwd, matrices_bwd, ctxts_set] = ds_ref;

        const auto& part = parts[search.part(p_idx)];
        pos_t k_max = search.upper_bound(p_idx);
        direction_t dir = search.part_dir(p_idx);
        bool dir_switch = search.does_part_switch_dir(p_idx);
        auto& stack = stacks[p_idx];
        auto& matrix = dir == RIGHT ? matrices_fwd[search.part(p_idx)] : matrices_bwd[search.part(p_idx)];

        const auto& dists_dir = dir_switch ? dists_diff_dir : dists_prev_dir;
        const auto& desc_dir = dir_switch ? desc_diff_dir : desc_prev_dir;
        const auto& dists_rev_dir = dir_switch ? dists_prev_dir : dists_diff_dir;
        const auto& desc_rev_dir = dir_switch ? desc_prev_dir : desc_diff_dir;

        std::vector<uint16_t> dists_new;

        if (dists_dir.empty()) {
            dists_new.emplace_back(ctx.errors());
        } else {
            uint16_t prev_dist = dir_switch ? *min_element(dists_dir.begin(), dists_dir.end()) : dists_dir[0];
            uint16_t add = ctx.errors() - prev_dist;
            dists_new.reserve(dists_dir.size());

            for (pos_t i = 0; i < dists_dir.size(); i++) {
                dists_new.emplace_back(dists_dir[i] + add);
            }
        }

        if (!matrix.is_initialized()) matrix.set_input(part, index.sigma, index.forward_index().map_int());
        matrix.init(k_max, dists_new);
        cluster_t cluster(matrix.last_col_size(), k_max, ctx.depth(), ctx.shift());

        if (matrix.is_in_final_column(0)) {
            auto ctx_copy = ctx;
            ctx_copy.set_depth(0);
            cluster.set_cell(0, ctx_copy, matrix(0, part.size()));
        }

        pos_t old_depth = 0;

        if (!desc_dir.empty()) {
            for (pos_t i = 0; i < desc_dir.size() && desc_dir[i].depth() < matrix.num_rows(); i++) {
                if (branch_and_bound<matrix_word_t>(matrix, cluster, desc_dir[i], p_idx, ds_ref,
                    dists_rev_dir, desc_rev_dir, {desc_dir.begin() + i + 1, desc_dir.end()})
                ) return;
            }

            if (desc_dir.back().depth() == matrix.num_rows() - 1) return;
            if (!dir_switch) ctx = desc_dir.back();
            old_depth = desc_dir.back().depth();
        }

        auto ext_ctx = ctx.prepare_extend_all(index, dir);
        bool extended = false;

        while (ext_ctx.can_extend(index)) {
            stack.emplace_back(ctx.extend_next(index, ext_ctx, dir));
            stack.back().set_depth(old_depth + 1);
            extended = true;
        }

        if (!extended && matrix.is_in_final_column(0)) {
            flush_cluster<matrix_word_t>(cluster, p_idx, ds_ref, dists_rev_dir, desc_rev_dir, {});
            return;
        }

        while (!stack.empty()) {
            auto ctx_cur = stack.back(); stack.pop_back();
            if (branch_and_bound<matrix_word_t>(matrix, cluster, ctx_cur, p_idx, ds_ref, dists_rev_dir, desc_rev_dir, {})) continue;
            auto ext_ctx_cur = ctx_cur.prepare_extend_all(index, dir);

            while (ext_ctx_cur.can_extend(index)) {
                stack.emplace_back(ctx_cur.extend_next(index, ext_ctx_cur, dir));
            }
        }
    }

    /**
     * @brief flushes a completed cluster (the final column of a part's alignment matrix): for the
     *        last part it reports the occurrences, for an inner edge/non-edge part it recurses into
     *        the next part with the appropriate cluster center(s)
     * @tparam matrix_word_t word type of the edit-distance matrix
     * @param cluster the cluster collecting the cells of the matrix' final column
     * @param p_idx index of the current part
     * @param ds_ref references to the data structures shared during the search
     * @param dists_rev_dir edit distances of the descendants in the reverse direction
     * @param desc_rev_dir descendant contexts in the reverse direction
     * @param desc_dir_rem remaining descendant contexts in the current direction
     */
    template <typename matrix_word_t>
    void flush_cluster(
        cluster_t& cluster, pos_t p_idx,
        data_struct_ref_t<matrix_word_t>& ds_ref,
        const std::vector<uint16_t>& dists_rev_dir,
        const std::vector<search_context_t<LOCATE>>& desc_rev_dir,
        const std::vector<search_context_t<LOCATE>>& desc_dir_rem
    ) const requires(move_rb_t::supports_locate)
    {
        auto& [search, parts, stacks, matrices_fwd, matrices_bwd, ctxts_set] = ds_ref;
        const auto& k_min = search.lower_bound(p_idx);

        if (search.is_edge(p_idx)) {
            if (p_idx + 1 == parts.size()) {
                auto centers = cluster.report_centers_at_end();

                for (auto& ctx_cur : centers) {
                    if (ctx_cur.errors() >= k_min) {
                        auto it = ctxts_set.find(ctx_cur);

                        if (it != ctxts_set.end()) {
                            const auto& ctx_old = *it;
                            uint16_t cur_err = ctx_cur.errors();
                            uint16_t old_err = ctx_old.errors();

                            if (cur_err < old_err || (cur_err == old_err && ctx_old.s.d == RIGHT && ctx_cur.s.d == LEFT)) {
                                it = ctxts_set.erase(it);
                                ctxts_set.emplace_hint(it, ctx_cur);
                            }
                        } else {
                            ctxts_set.emplace(ctx_cur);
                        }
                    }
                }
            } else {
                auto minima = cluster.plateau_minima(search.part_dir(p_idx));

                for (auto& ctx_min : minima)
                    if (ctx_min.errors() >= k_min)
                        rec_search<matrix_word_t>(ctx_min, p_idx + 1, ds_ref, {}, desc_rev_dir, {}, dists_rev_dir);
            }

            return;
        }

        std::vector<search_context_t<LOCATE>> desc_dir;
        std::vector<uint16_t> dists_dir;
        auto [ctx_cntr, res] = cluster.cluster_center(k_min, desc_dir, dists_dir);
        if (!res) return;

        desc_dir.insert(desc_dir.end(), desc_dir_rem.begin(), desc_dir_rem.end());

        for (pos_t i = 0; i < desc_dir.size(); i++) {
            desc_dir[i].set_depth(i + 1);
        }

        pos_t k_max = search.upper_bound(p_idx + 1);
        while (dists_dir.back() > k_max) dists_dir.pop_back();

        if (search.does_part_switch_dir(p_idx + 1) && !desc_dir.empty()) {
            pos_t depth_cntr = ctx_cntr.depth();
            ctx_cntr = desc_dir.back();
            ctx_cntr.set_depth(depth_cntr);
            ctx_cntr.set_errors(*std::min_element(dists_dir.begin(), dists_dir.end()));
        }

        rec_search<matrix_word_t>(ctx_cntr, p_idx + 1, ds_ref, desc_dir, desc_rev_dir, dists_dir, dists_rev_dir);
    }

    /**
     * @brief computes the edit-distance matrix row for ctx and, when the final column is reached, either reports
     *        the occurrences or recurses into the next part; bounds the search if no valid distance remains
     * @tparam matrix_word_t word type of the edit-distance matrix
     * @param matrix the edit-distance matrix of the current part
     * @param cluster the cluster collecting the cells of the matrix' final column
     * @param ctx the current search context
     * @param p_idx index of the current part
     * @param ds_ref references to the data structures shared during the search
     * @param dists_rev_dir edit distances of the descendants in the reverse direction
     * @param desc_rev_dir descendant contexts in the reverse direction
     * @param desc_dir_rem remaining descendant contexts in the current direction
     * @return true <=> the subtree rooted at ctx has been fully processed (it does not need to be extended further)
     */
    template <typename matrix_word_t>
    bool branch_and_bound(
        edit_distance_matrix<matrix_word_t>& matrix, cluster_t& cluster,
        const search_context_t<LOCATE>& ctx, pos_t p_idx,
        data_struct_ref_t<matrix_word_t>& ds_ref,
        const std::vector<uint16_t>& dists_rev_dir,
        const std::vector<search_context_t<LOCATE>>& desc_rev_dir,
        const std::vector<search_context_t<LOCATE>>& desc_dir_rem
    ) const requires(move_rb_t::supports_locate)
    {
        auto& [search, parts, stacks, matrices_fwd, matrices_bwd, ctxts_set] = ds_ref;
        pos_t row = ctx.depth();

        if (row >= matrix.num_rows()) return true;
        bool valid_dist = matrix.compute_row(row, ctx.last_added_symbol());

        if (matrix.is_in_final_column(row)) {
            pos_t cluster_idx = cluster.size() + row - matrix.num_rows();
            cluster.set_cell(cluster_idx, ctx, matrix(row, matrix.num_cols() - 1));
            bool terminal = !valid_dist || matrix.only_vertical_gaps_left(row);

            if (!terminal) {
                search_context_t<LOCATE> ctx_ext = ctx;
                terminal = !ctx_ext.prepare_extend_all(index, search.part_dir(p_idx)).can_extend(index);
            }

            if (terminal) {
                flush_cluster<matrix_word_t>(cluster, p_idx, ds_ref, dists_rev_dir, desc_rev_dir, desc_dir_rem);
                return true;
            }
        }

        return !valid_dist;
    }
};

// deduction guide so that move_rb_apm_edit{index} deduces the index template parameters
template <move_r_support support, typename sym_t, typename pos_t>
move_rb_apm_edit(const move_rb<support, sym_t, pos_t>&) -> move_rb_apm_edit<support, sym_t, pos_t>;