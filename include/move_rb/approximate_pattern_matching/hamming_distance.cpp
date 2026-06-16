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

#include <move_rb/move_rb.hpp>

/**
 * @brief approximate pattern matching (APM) helper for a move_rb index; holds a reference to the
 *        index and implements all hamming-distance count- and locate algorithms outside the
 *        move_rb class
 * @tparam support type of locate support (_locate_move or _locate_rlzsa)
 * @tparam sym_t value type (default: char for strings)
 * @tparam pos_t index integer type (use uint32_t if input size < UINT_MAX, else uint64_t)
 */
template <move_r_support support, typename sym_t, typename pos_t>
class move_rb_apm_hamming {
protected:
    using move_rb_t = move_rb<support, sym_t, pos_t>;
    using inp_t = typename move_rb_t::inp_t;
    template <move_rb_query_support_t query_support>
    using search_context_t = typename move_rb_t::template search_context_t<query_support>;
    template <move_rb_query_support_t query_support>
    using search_context_set_t = tsl::sparse_set<search_context_t<query_support>, typename search_context_t<query_support>::hash>;

    const move_rb_t& index; // the index to perform approximate pattern matching on

public:
    /**
     * @brief constructs an approximate-pattern-matching helper for the index
     * @param index the index to perform approximate pattern matching on
     */
    move_rb_apm_hamming(const move_rb_t& index) : index(index) {}

    /**
     * @brief counts a pattern with at most k mismatches (w.r.t. hamming distance)
     * @param P the pattern to search
     * @param scheme the search scheme to use (provides k)
     * @return number of occurrences of P in T with at most k mismatches (w.r.t. hamming distance)
     */
    pos_t count(const inp_t& P, const search_scheme_t& scheme) const
    {
        pos_t m = P.size();
        pos_t k = scheme.k;

        if (k >= m) return index.n - 1;
        if (k == 0) return index.forward_index().count(P);

        pos_t count = 0;

        auto ctxts_set = search<COUNT>(P, scheme);
        for (const auto& ctx : ctxts_set) count += ctx.num_occ();

        return count;
    }

    /**
     * @brief locates a pattern with at most k mismatches (w.r.t. hamming distance)
     * @tparam report_fnc_t type of the function report
     * @param P the pattern to search
     * @param scheme the search scheme to use (provides k)
     * @param report function that is called with every occurrence (occ, len, err) of P in T with at most k mismatches (w.r.t. hamming distance)
     */
    template <typename report_fnc_t>
    void locate(const inp_t& P, const search_scheme_t& scheme, report_fnc_t report) const
    {
        pos_t m = P.size();
        pos_t k = scheme.k;

        if (k == 0) {
            index.forward_index().locate(P, [&](pos_t occ){report({occ, m, 0});});
        } else {
            auto ctxts_set = search<LOCATE>(P, scheme);

            for (const auto& ctx : ctxts_set) {
                ctx.locate_phase().locate(index, ctx, [&](pos_t occ){report({occ, ctx.length(), ctx.errors()});});
            }
        }
    }

    /**
     * @brief executes a given search scheme for a pattern (w.r.t. hamming distance)
     * @tparam query_support COUNT or LOCATE
     * @param P the pattern to search
     * @param scheme the search scheme to use
     * @return all search contexts of P in T with at most k mismatches (assuming the search scheme covers all error configurations)
     */
    template <move_rb_query_support_t query_support>
    auto search(const inp_t& P, const search_scheme_t& scheme) const -> search_context_set_t<query_support>
    {
        pos_t m = P.size();
        pos_t p = scheme.p;
        if (m == 0) return {index.template empty_context<query_support>()};
        if (p >= m) throw std::runtime_error("p >= m");

        using search_state_t = std::tuple<
            search_context_t<query_support>, // ctx
            pos_t, // k_cur
            pos_t // match_pos_idx
        >;

        using match_pos_t = std::tuple<
            pos_t, // p_idx
            direction_t, // dir
            pos_t, // part
            pos_t, // beg
            pos_t, // end
            pos_t, // pos
            pos_t // ext_rem
        >;

        search_context_set_t<query_support> ctxts_set;
        auto empty_ctx = index.template empty_context<query_support>();

        auto add_ctx = [&](search_context_t<query_support>& ctx){
            if constexpr (query_support == LOCATE) {
                auto it = ctxts_set.find(ctx);

                if (it != ctxts_set.end()) {
                    const auto& ctx_old = *it;

                    if (ctx_old.s.d == RIGHT && ctx.s.d == LEFT) {
                        it = ctxts_set.erase(it);
                    }
                }
            
                ctxts_set.emplace_hint(it, ctx);
            } else {
                ctxts_set.emplace(ctx);
            }
        };

        for (const search_t& S : scheme.S) {
            std::vector<match_pos_t> match_pos;
            match_pos.reserve(m);

            for (pos_t p_idx = 0; p_idx < p; p_idx++) {
                pos_t part = S[p_idx].part;
                pos_t beg = (part * m) / p;
                pos_t end = part == p - 1 ? m : (((part + 1) * m) / p);
                direction_t dir;

                if (p_idx == 0) {
                    dir = p == 1 || S[p_idx + 1].part > part ? RIGHT : LEFT;
                } else {
                    dir = part < S[p_idx - 1].part ? LEFT : RIGHT;
                }

                if (dir == LEFT) {
                    for (pos_t pos = end; pos > beg; pos--) {
                        pos_t cur_pos = pos - 1;
                        match_pos.emplace_back(match_pos_t{p_idx, dir, part, beg, end, cur_pos, cur_pos - beg});
                    }
                } else {
                    for (pos_t pos = beg; pos < end; pos++) {
                        match_pos.emplace_back(match_pos_t{p_idx, dir, part, beg, end, pos, end - pos - 1});
                    }
                }
            }

            for (pos_t i = 0; i < m - 1; i++) {
                auto& [p_idx_cur, dir_cur, part_cur, beg_cur, end_cur, pos_cur, ext_rem_cur] = match_pos[i];
                auto& [p_idx_nxt, dir_nxt, part_nxt, beg_nxt, end_nxt, pos_nxt, ext_rem_nxt] = match_pos[i + 1];

                if (p_idx_nxt != p_idx_cur) {
                    ext_rem_cur = end_nxt - beg_nxt;
                }
            }

            std::vector<search_state_t> states_stack;
            states_stack.reserve(m);
            states_stack.emplace_back(search_state_t{empty_ctx, 0, 0});

            while (!states_stack.empty()) {
                auto [ctx, k_cur, match_pos_idx] = states_stack.back(); states_stack.pop_back();
                const auto& [p_idx, dir, part, beg, end, pos, ext_rem] = match_pos[match_pos_idx];
                pos_t match_pos_idx_nxt = match_pos_idx + 1;
                pos_t p_idx_nxt = match_pos_idx_nxt == m ? p : std::get<0>(match_pos[match_pos_idx_nxt]);

                if (k_cur < S[p_idx].k_max && (p_idx_nxt == p || S[p_idx_nxt].k_min <= k_cur + 1 + ext_rem)) {
                    auto ext_ctx = ctx.prepare_extend_all(index, dir);

                    while (ext_ctx.can_extend(index)) {
                        auto ctx_nxt = ctx.extend_next(index, ext_ctx, dir);
                        bool is_mismatch = ctx_nxt.last_added_symbol() != P[pos];
                        pos_t k_nxt = k_cur + is_mismatch;

                        if (p_idx_nxt == p) {
                            ctx_nxt.set_errors(k_nxt);
                            add_ctx(ctx_nxt);
                        } else if (is_mismatch || S[p_idx_nxt].k_min <= k_nxt + ext_rem) {
                            states_stack.emplace_back(search_state_t{ctx_nxt, k_nxt, match_pos_idx_nxt});
                        }
                    }
                } else {
                    auto [ctx_nxt, result] = ctx.extend(index, P[pos], dir);

                    if (result) {
                        if (p_idx_nxt == p) {
                            ctx_nxt.set_errors(k_cur);
                            add_ctx(ctx_nxt);
                        } else {
                            states_stack.emplace_back(search_state_t{ctx_nxt, k_cur, match_pos_idx_nxt});
                        }
                    }
                }
            }
        }

        return ctxts_set;
    }
};

// deduction guide so that move_rb_apm_hamming{index} deduces the index template parameters
template <move_r_support support, typename sym_t, typename pos_t>
move_rb_apm_hamming(const move_rb<support, sym_t, pos_t>&) -> move_rb_apm_hamming<support, sym_t, pos_t>;