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

#ifndef MOVE_R_MOVE_RB_SEARCH_CONTEXT_TPP
#define MOVE_R_MOVE_RB_SEARCH_CONTEXT_TPP

#include <move_rb/search_context.hpp>

template <move_r_support support>
template <query_support_t query_support>
auto move_rb<support>::search_context_t<query_support>::sample_at_run_boundary() const
    -> sample_info_pack_t requires(supports_locate)
{
    const auto& fwd = idx->idx_fwd;
    const auto& M_LF = fwd.M_LF();
    uint64_t x = idx->template input_interval_near<LEFT>(b, b_);

    // a BWT run starts at b
    if (b == M_LF.p(x) && fwd.run_starts_at(x)) return { LEFT, BEG, 0, 0, x };

    // else, use the SA-sample at the end of the run containing b
    uint64_t y = fwd.last_interval_of_run(x);
    uint64_t last = M_LF.p(y + 1) - 1;
    assert(last <= e);
    return { LEFT, END, last - b, 0, y };
}

template <move_r_support support>
template <query_support_t query_support>
inline uint64_t move_rb<support>::search_context_t<query_support>::hash::operator()(
    const search_context_t& ctx) const
{
    auto [b, e] = ctx.forward_sa_interval();
    uint64_t hash = hash_64(b);
    hash_combine(hash, hash_64(e));
    hash_combine(hash, hash_64(ctx.m));
    return hash;
}

template <move_r_support support>
template <query_support_t query_support>
inline void move_rb<support>::search_context_t<query_support>::reset()
{
    valid_interval_dirs = BOTH_DIRS;
    sym_lst = 0;
    m = 0;
    err = 0;
    // the empty pattern is the only entry of level 0 of the context tables
    ctx_entry = 0;
    ctx_entry_R = 0;

    b = 0;
    e = idx->idx_fwd.input_size();
    b_R = 0;
    e_R = e;
    b_ = 0;
    e_ = idx->idx_fwd.M_LF().num_intervals() - 1;
    b_R_ = 0;
    e_R_ = idx->idx_bwd.M_LF().num_intervals() - 1;

    if constexpr (query_support == LOCATE) {
        s = {
            LEFT,
            BEG,
            0,
            0,
            0
        };
    }
}

template <move_r_support support>
template <query_support_t query_support>
template <direction_t dir>
inline void move_rb<support>::search_context_t<query_support>::build_prev_next(
    int64_t* prev, int64_t* next, uint64_t max_sym) const
{
    const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = const_vars<dir>();
    build_prev_next<dir>(prev, next, max_sym, b_, e_);
}

template <move_r_support support>
template <query_support_t query_support>
template <direction_t dir>
inline void move_rb<support>::search_context_t<query_support>::build_prev_next(
    int64_t* prev, int64_t* next, uint64_t max_sym, uint64_t b_, uint64_t e_) const
{
    assert(max_sym < idx->sigma);

    const auto& idx_dir = idx->index<dir>();
    uint64_t blk_size = idx_dir.PSB_L_().block_size();
    uint64_t max = max_sym;

    uint64_t blk = div_ceil<uint64_t>(b_, blk_size);
    uint64_t end_pos = std::min<uint64_t>(blk * blk_size, e_);
    int64_t beg = b_;
    int64_t end = end_pos;

    if (end_pos != e_) [[likely]] {
        for (uint64_t i = 0; i <= max; i++) {
            next[i] = idx_dir.PSB_L_().succ_of_block(blk, i);
        }
    } else {
        std::fill_n(next, max_sym + 1, LONG_MAX);
    }

    for (int64_t i = end; i >= beg; i--) {
        next[idx_dir.L_(i)] = i;
    }

    blk = e_ / blk_size;
    uint64_t beg_pos = std::max<uint64_t>(blk * blk_size, b_);
    beg = beg_pos;
    end = e_;

    if (beg_pos != b_) [[likely]] {
        for (uint64_t i = 0; i <= max; i++) {
            prev[i] = idx_dir.PSB_L_().pred_of_block(blk, i);
        }
    } else {
        std::fill_n(prev, max_sym + 1, LONG_MIN);
    }

    for (int64_t i = beg; i <= end; i++) {
        prev[idx_dir.L_(i)] = i;
    }
}

template <move_r_support support>
template <query_support_t query_support>
template <direction_t dir>
inline uint64_t move_rb<support>::search_context_t<query_support>::num_smaller(
    const int64_t* prev, const int64_t* next,
    uint64_t sym, uint64_t b, uint64_t e, uint64_t b_, uint64_t e_) const
{
    const auto& idx_dir = idx->index<dir>();
    int64_t num_intervals = idx_dir.M_LF().num_intervals();

    // the terminator symbol occurs exactly once in L
    uint64_t num = next[0] <= prev[0] && prev[0] < num_intervals;

    /* for each symbol a < sym, LF maps the positions of a in the input intervals next[a] .. prev[a] onto one
       contiguous range, so their number is (q(prev[a]) + len(prev[a])) - q(next[a]) */
    for (uint64_t a = 1; a < sym; a++) {
        if (next[a] <= prev[a] && prev[a] < num_intervals) {
            uint64_t x = next[a];
            uint64_t y = prev[a];

            num += idx_dir.M_LF().p(y + 1) - idx_dir.M_LF().p(y);

            if (x < y) [[likely]] {
                num += idx_dir.M_LF().q(y) - idx_dir.M_LF().q(x);
            }
        }
    }

    // subtract the parts of the first and the last input interval that lie outside of [b, e]
    if (idx_dir.L_(b_) < sym) {
        num -= b - idx_dir.M_LF().p(b_);
    }

    if (idx_dir.L_(e_) < sym) {
        num -= idx_dir.M_LF().p(e_ + 1) - e - 1;
    }

    return num;
}

template <move_r_support support>
template <query_support_t query_support>
template <direction_t dir>
inline uint64_t move_rb<support>::search_context_t<query_support>::num_smaller_wavelet(
    uint64_t sym, uint64_t b, uint64_t e, uint64_t b_, uint64_t e_
) const {
    const auto& idx_dir = idx->template index<dir>();

    // a range of few input intervals is scanned
    if (e_ - b_ + 1 <= max_scanned_intervals_per_wavelet_level * idx_dir.WTPS_L_().num_levels()) {
        return num_smaller_scanned<dir>(sym, b, e, b_, e_);
    }

    uint64_t num = idx_dir.WTPS_L_().smaller(b_, e_ + 1, sym);

    // subtract the parts of the first and the last input interval that lie outside of [b, e]
    if (idx_dir.L_(b_) < sym) num -= b - idx_dir.M_LF().p(b_);
    if (idx_dir.L_(e_) < sym) num -= idx_dir.M_LF().p(e_ + 1) - e - 1;

    return num;
}

template <move_r_support support>
template <query_support_t query_support>
template <direction_t dir>
inline uint64_t move_rb<support>::search_context_t<query_support>::num_smaller_scanned(
    uint64_t sym, uint64_t b, uint64_t e, uint64_t b_, uint64_t e_
) const {
    const auto& idx_dir = idx->template index<dir>();
    const auto& M_LF = idx_dir.M_LF();
    uint64_t num = 0;
    uint64_t pos = M_LF.p(b_);

    // sum up the lengths of the input intervals in [b_, e_] with a symbol smaller than sym
    for (uint64_t x = b_; x <= e_; x++) {
        uint64_t next = M_LF.p(x + 1);
        if (idx_dir.L_(x) < sym) num += next - pos;
        pos = next;
    }

    // subtract the parts of the first and the last input interval that lie outside of [b, e]
    if (idx_dir.L_(b_) < sym) num -= b - M_LF.p(b_);
    if (idx_dir.L_(e_) < sym) num -= M_LF.p(e_ + 1) - e - 1;

    return num;
}

template <move_r_support support>
template <query_support_t query_support>
template <direction_t dir>
inline void move_rb<support>::search_context_t<query_support>::update_input_intervals()
{
    // b_ and e_ are up to date (see valid_interval_dirs)
    if (valid_interval_dirs == dir || valid_interval_dirs == BOTH_DIRS) [[likely]] return;

    auto&& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = vars<dir>();

    b_ = idx->template input_interval_near<dir>(b, b_);
    e_ = idx->template input_interval_near<dir>(e, e_);

    valid_interval_dirs = valid_interval_dirs == NO_DIR ? dir : BOTH_DIRS;
}

template <move_r_support support>
template <query_support_t query_support>
template <direction_t dir>
inline void move_rb<support>::search_context_t<query_support>::update_sample(
    search_context_t& ctx, const prime_tuple_t& prime
) const {
    const auto& idx_dir = idx->index<dir>();
    const auto& [b_old, e_old, b_R_old, e_R_old, b__old, e__old, b_R__old, e_R__old, s_kept] = const_vars<dir>();
    auto&& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = ctx.vars<dir>();
    // the SA-sample of this context (see sample)
    const sample_info_pack_t s_old = sample();

    /* a forward or an explicit SA-sample, whose position in the forward suffix array is known, stays the SA-sample if
       its occurrence is one of the extended pattern, i.e., if its position lies in the new forward interval */
    if constexpr (dir == RIGHT && stores_SA_pos) {
        if (s_old.d() == LEFT || s_old.d() == NO_DIR) {
            const uint64_t pos = b_R_old + s_old.o();

            if (b_R <= pos && pos <= e_R) {
                s = s_old;
                s.set_o(pos - b_R);
                return;
            }
        }
    }

    if (prime.b != b_old) {
        s = { dir, BEG, 0, 1, prime.b_ };
    } else if (prime.e != e_old) {
        s = { dir, END, e - b, 1, prime.e_ };
    } else if (e - b < e_old - b_old) {
        /* the symbols at both limits equal the extension symbol, but not all symbols between them, so the
           BWT run containing the left limit ends inside the old interval; use the SA-sample at its end */
        uint64_t y = idx_dir.last_interval_of_run(prime.b_);
        s = { dir, END, idx_dir.M_LF().p(y + 1) - 1 - prime.b, 1, y };
    } else {
        s = s_old;

        if (s_old.d() == dir) {
            s.set_i(s.i() + 1);
        }
    }
}

template <move_r_support support>
template <query_support_t query_support>
template <typename sym_t>
inline auto move_rb<support>::search_context_t<query_support>::extend(
    sym_t sym, direction_t dir) -> extend_res_t
{
    return dir == LEFT ? extend<LEFT>(sym) : extend<RIGHT>(sym);
}

template <move_r_support support>
template <query_support_t query_support>
template <direction_t dir, typename sym_t>
inline void move_rb<support>::search_context_t<query_support>::prefetch_extend(
    sym_t sym) const
{
    if (valid_interval_dirs != dir && valid_interval_dirs != BOTH_DIRS) return;
    if (ctx_entry != no_ctx_entry && m < idx->context_tables_max_length()) return;
    const auto& idx_dir = idx->index<dir>();
    if (idx_dir.PSO_L_().empty()) return;
    const uint64_t i_sym = idx->idx_fwd.map_symbol(sym);
    if (i_sym == 0) [[unlikely]] return;
    const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = const_vars<dir>();

    // the first and the last input interval with the symbol, as computed by extend
    const uint64_t first = idx_dir.L_(b_) == i_sym ? b_ : idx_dir.PSO_L_().succ(i_sym, b_, e_);
    if (first > e_) return;
    const uint64_t last = idx_dir.L_(e_) == i_sym ? e_ : idx_dir.PSO_L_().pred(i_sym, e_, b_);
    idx_dir.M_LF().prefetch_move(first);
    if (last != first) idx_dir.M_LF().prefetch_move(last);
}

template <move_r_support support>
template <query_support_t query_support>
template <direction_t dir>
inline void move_rb<support>::search_context_t<query_support>::extend_to_runs(
    search_context_t& ctx, uint64_t first, uint64_t last, uint64_t b_R_new) const
{
    static constexpr direction_t other_dir = dir == LEFT ? RIGHT : LEFT;
    const auto& M_LF = idx->template index<dir>().M_LF();
    const auto& [b_old, e_old, b_R_old, e_R_old, b__old, e__old, b_R__old, e_R__old, s_old] = const_vars<dir>();
    ctx.idx = idx;
    auto&& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = ctx.template vars<dir>();

    // [b, e] = the positions in [b_old, e_old] with the extension symbol, which lie in the input intervals first..last
    b_ = first;
    e_ = last;
    b = b_ != b__old ? M_LF.p(b_) : b_old;
    e = e_ != e__old ? M_LF.p(e_ + 1) - 1 : e_old;
    const prime_tuple_t prime { b, e, b_, e_ };

    // apply LF to [b, e]; for first = last, LF maps [b, e] to [LF(b), LF(b) + e - b]
    if (b_ == e_) {
        if (b == e) {
            M_LF.move(b, b_);
            e = b;
            e_ = b_;
        } else {
            uint64_t diff_eb = e - b;
            M_LF.move(b, b_);
            e = b + diff_eb;
            e_ = b_;

            while (e >= M_LF.p(e_ + 1)) {
                e_++;
            }
        }
    } else {
        M_LF.move(b, b_);
        M_LF.move(e, e_);
    }

    // both intervals have the same width; the input intervals of the other direction are not updated
    b_R = b_R_new;
    e_R = b_R + (e - b);
    b_R_ = b_R__old;
    e_R_ = e_R__old;

    if constexpr (query_support == LOCATE) {
        update_sample<dir>(ctx, prime);
    }

    ctx.valid_interval_dirs = dir;
    ctx.m = m + 1;
    ctx.err = 0;
}

template <move_r_support support>
template <query_support_t query_support>
template <direction_t dir, bool one_sided, typename sym_t>
inline auto move_rb<support>::search_context_t<query_support>::extend_impl(
    sym_t sym) -> extend_res_t
{
    const auto& idx_dir = idx->index<dir>();
    uint64_t i_sym = idx->idx_fwd.map_symbol(sym);

    // the extended pattern does not occur, if the symbol does not occur in the input
    if (i_sym == 0) [[unlikely]] {
        return {search_context_t{}, false};
    }

    // if the context tables store the extended pattern, load its context from them (one key lookup per order)
    constexpr context_extension_t key = ctx_key<dir>;

    if (ctx_entry != no_ctx_entry && m < idx->context_tables_max_length()
        && idx->context_tables().template has_keys<key, _forward_order>(m)) {
        const auto& tables = idx->context_tables();
        uint32_t entry = tables.template extend<key, _forward_order>(m, ctx_entry, i_sym);

        if (entry == no_ctx_entry) [[unlikely]] {
            return {search_context_t{}, false};
        }

        const uint32_t entry_R = tables.template extend<key, _reverse_order>(m, ctx_entry_R, i_sym);

        search_context_t ctx;
        idx->template load_context<query_support>(ctx, m + 1, entry, entry_R);
        ctx.sym_lst = sym;
        return {ctx, true};
    }

    update_input_intervals<dir>();

    const auto& [b_old, e_old, b_R_old, e_R_old, b__old, e__old, b_R__old, e_R__old, s_old] = const_vars<dir>();
    uint64_t first; // first occurrence of i_sym in L'[b__old..e__old]
    uint64_t last; // last occurrence of i_sym in L'[b__old..e__old]
    uint64_t num_smaller_syms = 0; // number of positions in [b_old, e_old] with a symbol smaller than i_sym

    // true <=> first, last and num_smaller_syms are computed from prev and next (see build_prev_next)
    bool use_prev_next = idx->has_PSB_L_() && idx_dir.WTPS_L_().empty()
        && e__old - b__old + 1 > idx->sigma;

    if (use_prev_next) {
        int64_t* prev = prev_next_buffer(idx->sigma);
        int64_t* next = prev + idx->sigma;
        build_prev_next<dir>(prev, next, i_sym);
        int64_t num_intervals = idx_dir.M_LF().num_intervals();

        if (next[i_sym] > prev[i_sym] || prev[i_sym] >= num_intervals) [[unlikely]] {
            return {search_context_t{}, false};
        }

        first = next[i_sym];
        last = prev[i_sym];
        if constexpr (!one_sided)
            num_smaller_syms = num_smaller<dir>(prev, next, i_sym, b_old, e_old, b__old, e__old);
    } else {
        // successor and predecessor queries on L'
        first = idx_dir.L_(b__old) == i_sym ? b__old
            : idx_dir.PSO_L_().succ(i_sym, b__old, e__old);

        if (first > e__old) [[unlikely]] {
            return {search_context_t{}, false};
        }

        last = idx_dir.L_(e__old) == i_sym ? e__old
            : idx_dir.PSO_L_().pred(i_sym, e__old, b__old);

        if constexpr (!one_sided) {
            if (!idx_dir.WTPS_L_().empty()) {
                num_smaller_syms = num_smaller_wavelet<dir>(i_sym, b_old, e_old, b__old, e__old);
            } else {
                num_smaller_syms = num_smaller_scanned<dir>(i_sym, b_old, e_old, b__old, e__old);
            }
        }
    }

    search_context_t ctx;
    extend_to_runs<dir>(ctx, first, last, b_R_old + num_smaller_syms);
    ctx.sym_lst = sym;

    return {ctx, true};
}

template <move_r_support support>
template <query_support_t query_support>
template <direction_t dir, extension_scope_t scope, typename fnc_t>
inline void move_rb<support>::search_context_t<query_support>::for_each_extension(
    fnc_t report)
{
    extend_context_t<query_support> ext_ctx(*idx, *this);
    ext_ctx.template prepare_extend_all<dir, scope>();

    while (ext_ctx.can_extend()) {
        uint64_t sym = ext_ctx.next_symbol();
        if (!report(sym, ext_ctx.template extend_next<dir>())) return;
    }
}

#endif // MOVE_R_MOVE_RB_SEARCH_CONTEXT_TPP
