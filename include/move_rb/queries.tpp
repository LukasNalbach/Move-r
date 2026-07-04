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
template <query_support_t query_support>
template <direction_t dir>
inline void move_rb<support, sym_t, pos_t>::search_context_t<query_support>::build_prev_next(
    std::array<int64_t, 256>& prev, std::array<int64_t, 256>& next,
    i_sym_t max_sym) const
{
    const auto& idx_dir = idx->index<dir>();
    pos_t blk_size = idx_dir.L_block_size();
    pos_t max = max_sym;

    const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = const_vars<dir>();

    pos_t blk = div_ceil<pos_t>(b_, blk_size);
    int64_t beg = b_;
    int64_t end = std::min<pos_t>(blk * blk_size, e_);

    if (end != e_) [[likely]] {
        pos_t blk_beg = blk * idx->sigma;

        for (pos_t i = 0; i <= max; i++) {
            next[i] = idx_dir.L_next(blk_beg + i);
        }
    } else {
        std::fill_n(next.begin(), max_sym + 1, LONG_MAX);
    }

    for (int64_t i = end; i >= beg; i--) {
        next[idx_dir.L_(i)] = i;
    }

    blk = e_ / blk_size;
    beg = std::max<pos_t>(blk * blk_size, b_);
    end = e_;

    if (beg != b_) [[likely]] {
        pos_t blk_beg = blk * idx->sigma;

        for (pos_t i = 0; i <= max; i++) {
            prev[i] = idx_dir.L_prev(blk_beg + i);
        }
    } else {
        std::fill_n(prev.begin(), max_sym + 1, LONG_MIN);
    }

    for (int64_t i = beg; i <= end; i++) {
        prev[idx_dir.L_(i)] = i;
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
template <query_support_t query_support>
template <direction_t dir>
inline void move_rb<support, sym_t, pos_t>::search_context_t<query_support>::update_input_intervals()
{
    const auto& idx_dir = idx->index<dir>();
    auto&& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = vars<dir>();

    if (dir_lst != NO_DIR && dir != dir_lst) {
        if (!(idx_dir.M_LF().p(b_) <= b && b < idx_dir.M_LF().p(b_ + 1))) {
            b_ = input_interval(b, idx->S_MLF_p<dir>(), [&](pos_t x){return idx_dir.M_LF().p(x);});
        }

        if (!(idx_dir.M_LF().p(e_) <= e && e < idx_dir.M_LF().p(e_ + 1))) {
            e_ = input_interval(e, idx->S_MLF_p<dir>(), [&](pos_t x){return idx_dir.M_LF().p(x);});
        }
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
template <query_support_t query_support>
template <direction_t dir>
inline void move_rb<support, sym_t, pos_t>::search_context_t<query_support>::update_sample(
    search_context_t& ctx, const prime_tuple_t& prime
) const {
    const auto& idx_dir = idx->index<dir>();
    const auto& [b_old, e_old, b_R_old, e_R_old, b__old, e__old, b_R__old, e_R__old, s_old] = const_vars<dir>();
    auto&& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = ctx.vars<dir>();

    if (prime.b != b_old) {
        s = { .d = dir, .t = BEG, .o = 0, .i = 1, .x = prime.b_ };
    } else if (prime.e != e_old) {
        s = { .d = dir, .t = END, .o = e - b, .i = 1, .x = prime.e_ };
    } else if (e - b < e_old - b_old) {
        s = { .d = dir, .t = END, .o = idx_dir.M_LF().p(prime.b_ + 1) - 1 - prime.b, .i = 1, .x = prime.b_ };
    } else {
        s = s_old;
        if (s_old.d == dir) s.i++;
    }

}

template <move_r_support support, typename sym_t, typename pos_t>
template <query_support_t query_support>
inline auto move_rb<support, sym_t, pos_t>::search_context_t<query_support>::extend(
    sym_t sym, direction_t dir) -> extend_res_t
{
    if (dir == LEFT) {
        return extend<LEFT>(sym);
    } else {
        return extend<RIGHT>(sym);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
template <query_support_t query_support>
template <direction_t dir>
inline auto move_rb<support, sym_t, pos_t>::search_context_t<query_support>::extend(
    sym_t sym) -> extend_res_t
{
    const auto& idx_dir = idx->index<dir>();
    i_sym_t i_sym = idx->idx_fwd.map_symbol(sym);

    // If i_sym does not occur in L', then P[i..m] does not occur in T
    if (i_sym == 0) [[unlikely]] {
        return {search_context_t{}, false};
    }

    update_input_intervals<dir>();

    std::array<int64_t, 256> prev;
    std::array<int64_t, 256> next;
    build_prev_next<dir>(prev, next, i_sym);

    if (next[i_sym] > prev[i_sym] || prev[i_sym] >= idx_dir.M_LF().num_intervals()) [[unlikely]] {
        return {search_context_t{}, false};
    }

    search_context_t ctx;
    ctx.idx = idx;
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

    prime_tuple_t prime{b, e, b_, e_};
    bool is_term_in_iv = next[0] <= prev[0] && prev[0] < idx_dir.M_LF().num_intervals();
    b_R = b_R_old + is_term_in_iv;

    for (i_sym_t i = 1; i < i_sym; i++) {
        if (next[i] <= prev[i] && prev[i] < idx_dir.M_LF().num_intervals()) {
            pos_t x = next[i];
            pos_t y = prev[i];

            pos_t add_p = idx_dir.M_LF().p(y + 1) - idx_dir.M_LF().p(y);
            b_R += add_p;
            pos_t add_q = 0;
            if (x < y) [[likely]] {
                add_q = idx_dir.M_LF().q(y) - idx_dir.M_LF().q(x);
                b_R += add_q;
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

    if constexpr (query_support == LOCATE) {
        update_sample<dir>(ctx, prime);
    }

    ctx.dir_lst = dir;
    ctx.sym_lst = sym;
    ctx.m = m + 1;
    ctx.err = 0;

    return {ctx, true};
}

// ############################# extend_context_t METHODS #############################

template <move_r_support support, typename sym_t, typename pos_t>
template <query_support_t query_support>
auto move_rb<support, sym_t, pos_t>::extend_context_t<query_support>::prepare_extend_all(
    direction_t dir) -> extend_context_t&
{
    this->dir = dir;

    if (dir == LEFT) {
        prepare_extend_all<LEFT>();
    } else {
        prepare_extend_all<RIGHT>();
    }

    return *this;
}

template <move_r_support support, typename sym_t, typename pos_t>
template <query_support_t query_support>
template <direction_t dir>
void move_rb<support, sym_t, pos_t>::extend_context_t<query_support>::prepare_extend_all()
{
    const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = ctx->template const_vars<dir>();

    ctx->template update_input_intervals<dir>();
    ctx->template build_prev_next<dir>(prev, next, idx->sigma);

    bool is_term_in_iv = next[0] <= prev[0] && prev[0] < idx->template index<dir>().M_LF().num_intervals();
    b_R_nxt = b_R + is_term_in_iv;

    sym_nxt = 0;
    this->template advance_symbol<dir>();
}

template <move_r_support support, typename sym_t, typename pos_t>
template <query_support_t query_support>
auto move_rb<support, sym_t, pos_t>::extend_context_t<query_support>::extend_next() -> search_context_t<query_support>
{
    if (dir == LEFT) {
        return extend_next<LEFT>();
    } else {
        return extend_next<RIGHT>();
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
template <query_support_t query_support>
template <direction_t dir>
auto move_rb<support, sym_t, pos_t>::extend_context_t<query_support>::extend_next() -> search_context_t<query_support>
{
    const auto& idx_dir = idx->template index<dir>();
    i_sym_t i_sym = sym_nxt;
    search_context_t<query_support> new_ctx;
    new_ctx.idx = idx;

    const auto& [b_old, e_old, b_R_old, e_R_old, b__old, e__old, b_R__old, e_R__old, s_old] = ctx->template const_vars<dir>();
    auto&& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = new_ctx.template vars<dir>();

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
    
    prime_tuple_t prime{b, e, b_, e_};

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

    b_R = b_R_nxt;
    e_R = b_R_nxt + (e - b);

    b_R_ = b_R__old;
    e_R_ = e_R__old;

    if constexpr (query_support == LOCATE) {
        ctx->template update_sample<dir>(new_ctx, prime);
    }

    new_ctx.sym_lst = idx->unmap_symbol(i_sym);
    new_ctx.dir_lst = dir;
    new_ctx.m = ctx->m + 1;
    new_ctx.err = 0;

    b_R_nxt = e_R + 1;
    this->template advance_symbol<dir>();

    return new_ctx;
}

// ############################# locate_context_t METHODS #############################

template <move_r_support support, typename sym_t, typename pos_t>
inline pos_t move_rb<support, sym_t, pos_t>::locate_context_t::compute_center()
{
    const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = ctx->template const_vars<LEFT>();

    if (s.d == LEFT) {
        c = b + s.o;
    } else /* if (s.d == RIGHT) */ {
        c = s.t == BEG ? idx->_SA_s_pos[s.x] : idx->_SA_e_pos[s.x];
        pos_t c_ = input_interval(c, idx->_S_MLF_p_fwd, [&](pos_t x){return idx->idx_fwd.M_LF().p(x);});
        assert(idx->idx_fwd.M_LF().p(c_) <= c && c < idx->idx_fwd.M_LF().p(c_ + 1));

        for (pos_t j = 0; j < ctx->m - s.i; j++) {
            idx->idx_fwd.M_LF().move(c, c_);
        }
    }

    assert(b <= c && c <= e);
    return c;
}

template <move_r_support support, typename sym_t, typename pos_t>
inline pos_t move_rb<support, sym_t, pos_t>::locate_context_t::first_occ()
{
    if (dir == NO_DIR) {
        compute_center();
    }

    const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = ctx->template const_vars<LEFT>();
    i = c;

    if (s.d == LEFT) {
        SA_i = s.t == BEG ? idx->idx_fwd.SA_s(s.x) : idx->idx_fwd.SA_e(s.x);
        SA_i -= s.i;
    } else /* if (s.d == RIGHT) */ {
        SA_i = s.t == BEG ? idx->idx_bwd.SA_s(s.x) : idx->idx_bwd.SA_e(s.x);
        SA_i = idx->n - SA_i - (ctx->m - s.i + 1);
    }

    SA_c = SA_i;

    if (occ_rem > 1) [[likely]] {
        dir = i > b ? LEFT : RIGHT;

        if constexpr (support == _locate_move) {
            if (dir == LEFT) {
                s_ = input_interval(SA_i, idx->_S_MPhi_p, [&](pos_t x){return idx->idx_fwd.M_Phi().p(x);});
            } else {
                s_ = input_interval(SA_i, idx->_S_MPhi_m1_p, [&](pos_t x){return idx->idx_fwd.M_Phi_m1().p(x);});
            }
        } else if constexpr (support == _locate_rlzsa) {
            if (dir == LEFT) {
                rlz_l = idx->idx_fwd.rlzsa().decode();
                rlz_l.set_value(SA_i);
                rlz_l.init_left(c);
            }

            if (c < e) [[likely]] {
                rlz_r = idx->idx_fwd.rlzsa().decode();
                rlz_r.set_value(SA_i);
                rlz_r.init_right(c + 1);
            }
        }
    }

    occ_rem--;
    return SA_i;
}

template <move_r_support support, typename sym_t, typename pos_t>
inline pos_t move_rb<support, sym_t, pos_t>::locate_context_t::next_occ()
{
    if (dir == NO_DIR) [[unlikely]] {
        return first_occ();
    }

    const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = ctx->template const_vars<LEFT>();
    pos_t occ;

    if constexpr (support == _locate_move) {
        if (dir == LEFT) {
            idx->idx_fwd.M_Phi().move(SA_i, s_);
            occ = SA_i;
            i--;

            if (i == b) [[unlikely]] {
                i = c;
                SA_i = SA_c;
                dir = RIGHT;
                s_ = input_interval(SA_i, idx->_S_MPhi_m1_p, [&](pos_t x){return idx->idx_fwd.M_Phi_m1().p(x);});
            }
        } else /* if (dir == RIGHT) */ {
            idx->idx_fwd.M_Phi_m1().move(SA_i, s_);
            occ = SA_i;
            i++;
        }
    } else if constexpr (support == _locate_rlzsa) {
        if (dir == LEFT) {
            rlz_l.prev();
            occ = rlz_l.value();

            if (rlz_l.pos() == b && c < e) [[unlikely]] {
                dir = RIGHT;
                rlz_r.set_value(SA_c);
            }
        } else /* if (dir == RIGHT) */ {
            rlz_r.next();
            occ = rlz_r.value();
        }
    }

    occ_rem--;
    return occ;
}

template <move_r_support support, typename sym_t, typename pos_t>
template <typename report_fnc_t>
inline void move_rb<support, sym_t, pos_t>::locate_context_t::locate(report_fnc_t report)
{
    static constexpr bool report_pos = function_traits<report_fnc_t>::arity > 1;
    const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = ctx->template const_vars<LEFT>();

    if (occ_rem == ctx->num_occ()) {
        first_occ();
        if constexpr (report_pos) report(c, SA_c); else report(SA_c);
    }

    if constexpr (support == _locate_move) {
        if (dir == LEFT) {
            while (i > b) {
                idx->idx_fwd.M_Phi().move(SA_i, s_);
                i--;
                if constexpr (report_pos) report(i, SA_i); else report(SA_i);
            }

            if (c < e) {
                i = c;
                SA_i = SA_c;
                s_ = input_interval(SA_i, idx->_S_MPhi_m1_p, [&](pos_t x){return idx->idx_fwd.M_Phi_m1().p(x);});
            }
        }

        if (c < e) {
            while (i < e) {
                idx->idx_fwd.M_Phi_m1().move(SA_i, s_);
                i++;
                if constexpr (report_pos) report(i, SA_i); else report(SA_i);
            }
        }
    } else if constexpr (support == _locate_rlzsa) {
        if (dir == LEFT) {
            rlz_l.report_left(b, report);

            if (c < e) [[likely]] {
                dir = RIGHT;
                rlz_r.set_value(SA_c);
            }
        }

        if (c < e) [[likely]] {
            rlz_r.report_right(e, report);
        }
    }

    occ_rem = 0;
}

template <move_r_support support, typename sym_t, typename pos_t>
template <query_support_t query_support>
template <direction_t dir>
void move_rb<support, sym_t, pos_t>::extend_context_t<query_support>::advance_symbol()
{
    do {
        sym_nxt++;
    } while (sym_nxt < idx->sigma && (
             next[sym_nxt] > prev[sym_nxt] ||
             prev[sym_nxt] == idx->template index<dir>().M_LF().num_intervals()));
}
