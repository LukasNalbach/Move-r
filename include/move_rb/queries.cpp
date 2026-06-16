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

// ############################# search_context_t METHODS #############################

template <move_r_support support, typename sym_t, typename pos_t>
template <move_rb_query_support_t query_support>
template <direction_t dir>
inline void move_rb<support, sym_t, pos_t>::search_context_t<query_support>::build_prev_next(
    const move_rb<support, sym_t, pos_t>& idx,
    std::array<int64_t, 256>& prev, std::array<int64_t, 256>& next,
    i_sym_t max_sym) const
{
    const auto& idx_dir = idx.index<dir>();
    pos_t blk_size = idx_dir.L_block_size();
    pos_t max = max_sym;

    std::fill_n(next.begin(), max_sym + 1, std::numeric_limits<int64_t>::max());
    std::fill_n(prev.begin(), max_sym + 1, std::numeric_limits<int64_t>::min());

    const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = const_vars<dir>();

    pos_t blk = div_ceil<pos_t>(b_, blk_size);
    int64_t beg = b_;
    int64_t end = std::min<pos_t>(blk * blk_size, e_);

    if (end != e_) [[likely]] {
        pos_t blk_beg = blk * idx.sigma;

        for (pos_t i = 0; i <= max; i++) {
            next[i] = idx_dir.L_next(blk_beg + i);
        }
    }

    for (int64_t i = end; i >= beg; i--) {
        next[idx_dir.L_(i)] = i;
    }

    blk = e_ / blk_size;
    beg = std::max<pos_t>(blk * blk_size, b_);
    end = e_;

    if (beg != b_) [[likely]] {
        pos_t blk_beg = blk * idx.sigma;

        for (pos_t i = 0; i <= max; i++) {
            prev[i] = idx_dir.L_prev(blk_beg + i);
        }
    }

    for (int64_t i = beg; i <= end; i++) {
        prev[idx_dir.L_(i)] = i;
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
template <move_rb_query_support_t query_support>
template <direction_t dir>
inline void move_rb<support, sym_t, pos_t>::search_context_t<query_support>::update_input_intervals(
    const move_rb<support, sym_t, pos_t>& idx)
{
    const auto& idx_dir = idx.index<dir>();
    auto&& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = vars<dir>();

    if (dir_lst != NO_DIR && dir != dir_lst) {
        if (!(idx_dir.M_LF().p(b_) <= b && b < idx_dir.M_LF().p(b_ + 1))) {
            b_ = input_interval(b, idx.S_MLF_p<dir>(), [&](pos_t x){return idx_dir.M_LF().p(x);});
        }

        if (!(idx_dir.M_LF().p(e_) <= e && e < idx_dir.M_LF().p(e_ + 1))) {
            e_ = input_interval(e, idx.S_MLF_p<dir>(), [&](pos_t x){return idx_dir.M_LF().p(x);});
        }
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
template <move_rb_query_support_t query_support>
template <direction_t dir>
inline void move_rb<support, sym_t, pos_t>::search_context_t<query_support>::update_sample(
    const move_rb<support, sym_t, pos_t>& idx, search_context_t& ctx) const
{
    const auto& idx_dir = idx.index<dir>();
    const auto& [b_old, e_old, b_R_old, e_R_old, b__old, e__old, b_R__old, e_R__old, s_old] = const_vars<dir>();
    auto&& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = ctx.vars<dir>();

    if (b_ != b__old) {
        s = { .d = dir, .t = BEG, .i = 1, .x = b_ };
    } else if (e_ != e__old) {
        s = { .d = dir, .t = END, .i = 1, .x = e_ };
    } else if (b_ != e_) {
        s = { .d = dir, .t = MID, .o = e - idx_dir.M_LF().p(e_), .i = 1, .x = e_ };
    } else {
        s = s_old;
        if (s_old.d == dir) s.i++;
    }

    s.shft = s_old.shft;
    s.dpth = s_old.dpth + 1;
    s.rprtd = false;
}

template <move_r_support support, typename sym_t, typename pos_t>
template <move_rb_query_support_t query_support>
inline auto move_rb<support, sym_t, pos_t>::search_context_t<query_support>::extend(
    const move_rb<support, sym_t, pos_t>& idx, sym_t sym, direction_t dir) -> extend_res_t
{
    if (dir == LEFT) {
        return extend<LEFT>(idx, sym);
    } else {
        return extend<RIGHT>(idx, sym);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
template <move_rb_query_support_t query_support>
template <direction_t dir>
inline auto move_rb<support, sym_t, pos_t>::search_context_t<query_support>::extend(
    const move_rb<support, sym_t, pos_t>& idx, sym_t sym) -> extend_res_t
{
    const auto& idx_dir = idx.index<dir>();
    i_sym_t i_sym = idx.idx_fwd.map_symbol(sym);

    // If i_sym does not occur in L', then P[i..m] does not occur in T
    if (i_sym == 0) [[unlikely]] {
        return {search_context_t{}, false};
    }

    update_input_intervals<dir>(idx);

    std::array<int64_t, 256> prev;
    std::array<int64_t, 256> next;
    build_prev_next<dir>(idx, prev, next, i_sym);

    if (next[i_sym] > prev[i_sym] || prev[i_sym] >= idx_dir.M_LF().num_intervals()) [[unlikely]] {
        return {search_context_t{}, false};
    }

    search_context_t ctx;
    const auto& [b_old, e_old, b_R_old, e_R_old, b__old, e__old, b_R__old, e_R__old, s_old] = const_vars<dir>();
    auto&& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = ctx.vars<dir>();

    b_ = next[i_sym];
    e_ = prev[i_sym];

    if (b_ != b__old) {
        b = idx_dir.M_LF().p(b_);
    } else {
        b = b_old;
    }

    if (e_ != e__old) {
        e = idx_dir.M_LF().p(e_ + 1) - 1;
    } else {
        e = e_old;
    }

    if constexpr (query_support == LOCATE) {
        update_sample<dir>(idx, ctx);
    }

    b_R = b_R_old;

    for (i_sym_t i = 0; i < i_sym; i++) {
        if (next[i] <= prev[i] && prev[i] < idx_dir.M_LF().num_intervals()) {
            pos_t x = next[i];
            pos_t y = prev[i];

            b_R += idx_dir.M_LF().p(y + 1) - idx_dir.M_LF().p(y);

            if (x < y) [[likely]] {
                b_R += idx_dir.M_LF().q(y) - idx_dir.M_LF().q(x);
            }
        }
    }

    if (idx_dir.L_(b__old) < i_sym) {
        b_R -= b_old - idx_dir.M_LF().p(b__old);
    }

    if (idx_dir.L_(e__old) < i_sym) {
        b_R -= idx_dir.M_LF().p(e__old + 1) - e_old - 1;
    }

    if (b_ == e_) {
        if (b == e) {
            idx_dir.M_LF().move(b, b_);
            e = b;
            e_ = b_;
        } else {
            pos_t diff_eb = e - b;
            idx_dir.M_LF().move(b, b_);
            e = b + diff_eb;
            e_ = b_;

            while (e >= idx_dir.M_LF().p(e_ + 1)) {
                e_++;
            }
        }
    } else {
        idx_dir.M_LF().move(b, b_);
        idx_dir.M_LF().move(e, e_);
    }

    e_R = b_R + (e - b);

    b_R_ = b_R__old;
    e_R_ = e_R__old;

    ctx.dir_lst = dir;
    ctx.sym_lst = sym;
    ctx.m = m + 1;
    ctx.err = 0;

    return {ctx, true};
}

template <move_r_support support, typename sym_t, typename pos_t>
template <move_rb_query_support_t query_support>
auto move_rb<support, sym_t, pos_t>::search_context_t<query_support>::prepare_extend_all(
    const move_rb<support, sym_t, pos_t>& idx, direction_t dir) -> extend_context_t
{
    extend_context_t ext_ctx;

    if (dir == LEFT) {
        prepare_extend_all<LEFT>(idx, ext_ctx);
    } else {
        prepare_extend_all<RIGHT>(idx, ext_ctx);
    }

    return ext_ctx;
}

template <move_r_support support, typename sym_t, typename pos_t>
template <move_rb_query_support_t query_support>
template <direction_t dir>
void move_rb<support, sym_t, pos_t>::search_context_t<query_support>::prepare_extend_all(
    const move_rb<support, sym_t, pos_t>& idx, extend_context_t& ext_ctx)
{
    const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = const_vars<dir>();
    auto& prev = ext_ctx.prev;
    auto& next = ext_ctx.next;

    update_input_intervals<dir>(idx);
    build_prev_next<dir>(idx, prev, next, idx.sigma);

    ext_ctx.sym_nxt = 0;
    ext_ctx.b_R_nxt = b_R + (next[0] <= prev[0] && prev[0] != idx.index<dir>().M_LF().num_intervals());
    ext_ctx.template next_symbol<dir>(idx);
}

template <move_r_support support, typename sym_t, typename pos_t>
template <move_rb_query_support_t query_support>
auto move_rb<support, sym_t, pos_t>::search_context_t<query_support>::extend_next(
    const move_rb<support, sym_t, pos_t>& idx, extend_context_t& ext_ctx, direction_t dir) -> search_context_t
{
    if (dir == LEFT) {
        return extend_next<LEFT>(idx, ext_ctx);
    } else {
        return extend_next<RIGHT>(idx, ext_ctx);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
template <move_rb_query_support_t query_support>
template <direction_t dir>
auto move_rb<support, sym_t, pos_t>::search_context_t<query_support>::extend_next(
    const move_rb<support, sym_t, pos_t>& idx, extend_context_t& ext_ctx) const -> search_context_t
{
    const auto& idx_dir = idx.index<dir>();
    i_sym_t i_sym = ext_ctx.sym_nxt;
    search_context_t ctx;

    const auto& [b_old, e_old, b_R_old, e_R_old, b__old, e__old, b_R__old, e_R__old, s_old] = const_vars<dir>();
    auto&& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = ctx.vars<dir>();

    b_ = ext_ctx.next[i_sym];
    e_ = ext_ctx.prev[i_sym];

    if (b_ != b__old) {
        b = idx_dir.M_LF().p(b_);
    } else {
        b = b_old;
    }

    if (e_ != e__old) {
        e = idx_dir.M_LF().p(e_ + 1) - 1;
    } else {
        e = e_old;
    }

    if constexpr (query_support == LOCATE) {
        update_sample<dir>(idx, ctx);
    }

    if (b_ == e_) {
        if (b == e) {
            idx_dir.M_LF().move(b, b_);
            e = b;
            e_ = b_;
        } else {
            pos_t diff_eb = e - b;
            idx_dir.M_LF().move(b, b_);
            e = b + diff_eb;
            e_ = b_;

            while (e >= idx_dir.M_LF().p(e_ + 1)) {
                e_++;
            }
        }
    } else {
        idx_dir.M_LF().move(b, b_);
        idx_dir.M_LF().move(e, e_);
    }

    b_R = ext_ctx.b_R_nxt;
    e_R = ext_ctx.b_R_nxt + (e - b);

    b_R_ = b_R__old;
    e_R_ = e_R__old;

    ctx.sym_lst = idx.unmap_symbol(i_sym);
    ctx.dir_lst = dir;
    ctx.m = m + 1;
    ctx.err = 0;

    ext_ctx.b_R_nxt = e_R + 1;
    ext_ctx.template next_symbol<dir>(idx);

    return ctx;
}

// ############################# locate_context_t METHODS #############################

template <move_r_support support, typename sym_t, typename pos_t>
inline pos_t move_rb<support, sym_t, pos_t>::locate_context_t::compute_center(
    const move_rb<support, sym_t, pos_t>& idx, const search_context_t<LOCATE>& ctx)
{
    const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = ctx.template const_vars<LEFT>();

    if (s.d == LEFT) {
        switch (s.t) {
            case BEG: c = b; break;
            case END: c = e; break;
            case MID: c = e - s.o; break;
            default: __builtin_unreachable();
        }
    } else /* if (s.d == RIGHT) */ {
        c = s.t != END ? idx._SA_sR_m1[s.x] : idx._SA_eR_m1[s.x];
        pos_t c_ = input_interval(c, idx._S_MLF_p_fwd, [&](pos_t x){return idx.idx_fwd.M_LF().p(x);});

        for (pos_t j = 0; j < ctx.m - s.i; j++) {
            idx.idx_fwd.M_LF().move(c, c_);
        }
    }

    assert(b <= c && c <= e);
    return c;
}

template <move_r_support support, typename sym_t, typename pos_t>
inline pos_t move_rb<support, sym_t, pos_t>::locate_context_t::first_occ(
    const move_rb<support, sym_t, pos_t>& idx, const search_context_t<LOCATE>& ctx)
{
    if (dir == NO_DIR) {
        compute_center(idx, ctx);
    }

    const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = ctx.template const_vars<LEFT>();
    i = c;

    if (s.d == LEFT) {
        SA_i = s.t != END ? idx.idx_fwd.SA_s(s.x) : idx.idx_fwd.SA_s_(s.x);
        SA_i -= s.i;
    } else /* if (s.d == RIGHT) */ {
        SA_i = s.t != END ? idx.idx_bwd.SA_s(s.x) : idx.idx_bwd.SA_s_(s.x);
        SA_i = idx.n - SA_i - (ctx.m - s.i + 1);
    }

    SA_c = SA_i;

    if (occ_rem > 1) [[likely]] {
        dir = i > b ? LEFT : RIGHT;

        if constexpr (support == _locate_move) {
            if (dir == LEFT) {
                s_ = input_interval(SA_i, idx._S_MPhi_p, [&](pos_t x){return idx.idx_fwd.M_Phi().p(x);});
            } else {
                s_ = input_interval(SA_i, idx._S_MPhi_m1_p, [&](pos_t x){return idx.idx_fwd.M_Phi_m1().p(x);});
            }
        } else if constexpr (support == _locate_rlzsa) {
            if (dir == LEFT) {
                idx.idx_fwd.init_rlzsa(i, rlz_l.x_p, rlz_l.x_lp, rlz_l.x_cp, rlz_l.x_r, rlz_l.s_p);

                if (c < e) [[likely]] {
                    rlz_r = rlz_l;
                    pos_t tmp_1 = c;
                    pos_t tmp_2;
                    idx.idx_fwd.next_rlzsa(tmp_1, tmp_2, rlz_r.x_p, rlz_r.x_lp, rlz_r.x_cp, rlz_r.x_r, rlz_r.s_p);
                }

                idx.idx_fwd.turn_rlzsa_left(rlz_l.x_p, rlz_l.x_lp, rlz_l.x_cp, rlz_l.x_r, rlz_l.s_p);
            } else {
                i++;
                idx.idx_fwd.init_rlzsa(i, rlz_r.x_p, rlz_r.x_lp, rlz_r.x_cp, rlz_r.x_r, rlz_r.s_p);
            }
        }
    }

    occ_rem--;
    return SA_i;
}

template <move_r_support support, typename sym_t, typename pos_t>
inline pos_t move_rb<support, sym_t, pos_t>::locate_context_t::next_occ(
    const move_rb<support, sym_t, pos_t>& idx, const search_context_t<LOCATE>& ctx)
{
    if (dir == NO_DIR) [[unlikely]] {
        return first_occ(idx, ctx);
    }

    const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = ctx.template const_vars<LEFT>();
    pos_t occ;

    if constexpr (support == _locate_move) {
        if (dir == LEFT) {
            idx.idx_fwd.M_Phi().move(SA_i, s_);
            occ = SA_i;
            i--;

            if (i == b) [[unlikely]] {
                i = c;
                SA_i = SA_c;
                dir = RIGHT;
                s_ = input_interval(SA_i, idx._S_MPhi_m1_p, [&](pos_t x){return idx.idx_fwd.M_Phi_m1().p(x);});
            }
        } else /* if (dir == RIGHT) */ {
            idx.idx_fwd.M_Phi_m1().move(SA_i, s_);
            occ = SA_i;
            i++;
        }
    } else if constexpr (support == _locate_rlzsa) {
        if (dir == LEFT) {
            idx.idx_fwd.prev_rlzsa(i, SA_i, rlz_l.x_p, rlz_l.x_lp, rlz_l.x_cp, rlz_l.x_r, rlz_l.s_p);
            occ = SA_i;

            if (i == b && c < e) [[unlikely]] {
                i = c + 1;
                dir = RIGHT;
                SA_i = SA_c;
            }
        } else /* if (dir == RIGHT) */ {
            idx.idx_fwd.next_rlzsa(i, SA_i, rlz_r.x_p, rlz_r.x_lp, rlz_r.x_cp, rlz_r.x_r, rlz_r.s_p);
            occ = SA_i;
        }
    }

    occ_rem--;
    return occ + s.shft;
}

template <move_r_support support, typename sym_t, typename pos_t>
template <typename report_fnc_t>
inline void move_rb<support, sym_t, pos_t>::locate_context_t::locate(
    const move_rb<support, sym_t, pos_t>& idx, const search_context_t<LOCATE>& ctx, report_fnc_t report)
{
    static constexpr bool report_pos = function_traits<report_fnc_t>::arity > 1;
    const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = ctx.template const_vars<LEFT>();

    if (occ_rem == ctx.num_occ()) {
        first_occ(idx, ctx);
        if constexpr (report_pos) report(c, SA_c + s.shft); else report(SA_c + s.shft);
    }

    if constexpr (support == _locate_move) {
        if (dir == LEFT) {
            while (i > b) {
                idx.idx_fwd.M_Phi().move(SA_i, s_);
                i--;
                if constexpr (report_pos) report(i, SA_i + s.shft); else report(SA_i + s.shft);
            }

            if (c < e) {
                i = c;
                SA_i = SA_c;
                s_ = input_interval(SA_i, idx._S_MPhi_m1_p, [&](pos_t x){return idx.idx_fwd.M_Phi_m1().p(x);});
            }
        }

        if (c < e) {
            while (i < e) {
                idx.idx_fwd.M_Phi_m1().move(SA_i, s_);
                i++;
                if constexpr (report_pos) report(i, SA_i + s.shft); else report(SA_i + s.shft);
            }
        }
    } else if constexpr (support == _locate_rlzsa) {
        SA_i += s.shft;

        if (dir == LEFT) {
            idx.idx_fwd.report_rlzsa_left(i, b, SA_i,
                rlz_l.x_p, rlz_l.x_lp, rlz_l.x_cp, rlz_l.x_r, rlz_l.s_p,
                report);

            if (c < e) [[likely]] {
                i = c + 1;
                dir = RIGHT;
                SA_i = SA_c + s.shft;
            }
        }

        if (c < e) [[likely]] {
            idx.idx_fwd.report_rlzsa_right(i, e, SA_i,
                rlz_r.x_p, rlz_r.x_lp, rlz_r.x_cp, rlz_r.x_r, rlz_r.s_p, report);
        }
    }

    occ_rem = 0;
}

// ############################# extend_context_t METHODS #############################

template <move_r_support support, typename sym_t, typename pos_t>
template <direction_t dir>
void move_rb<support, sym_t, pos_t>::extend_context_t::next_symbol(
    const move_rb<support, sym_t, pos_t>& idx)
{
    do {
        sym_nxt++;
    } while (sym_nxt < idx.sigma && (
             next[sym_nxt] > prev[sym_nxt] ||
             prev[sym_nxt] == idx.index<dir>().M_LF().num_intervals()));
}
