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

#ifndef MOVE_R_MOVE_RB_EXTEND_CONTEXT_TPP
#define MOVE_R_MOVE_RB_EXTEND_CONTEXT_TPP

#include <move_rb/extend_context.hpp>

template <move_r_support support>
template <query_support_t query_support>
template <extension_scope_t scope>
auto move_rb<support>::extend_context_t<query_support>::prepare_extend_all(
    direction_t dir) -> extend_context_t&
{
    this->dir = dir;
    if (dir == LEFT) prepare_extend_all<LEFT, scope>(); else prepare_extend_all<RIGHT, scope>();
    return *this;
}

template <move_r_support support>
template <query_support_t query_support>
template <direction_t dir, extension_scope_t scope>
void move_rb<support>::extend_context_t<query_support>::prepare_extend_all()
{
    // true <=> the extension with the terminator is reported (only possible to the right, see extension_scope_t)
    constexpr bool terminator = scope == TERM && dir == RIGHT;

    constexpr context_extension_t key = ctx_key<dir>;
    constexpr context_order_t group_order = ctx_group_order<dir>;

    /* if the context tables store the extensions of the pattern, enumerate them from the tables (they do
       not store the extension with the terminator) */
    if constexpr (!terminator) {
        if (ctx->ctx_entry != no_ctx_entry && ctx->m < idx->context_tables_max_length()) {
            const auto& tables = idx->context_tables();

            // the extensions are consecutive entries in group_order
            if (tables.template has_extension_range<key, group_order>(ctx->m)) {
                const uint32_t in_group = group_order == _forward_order ? ctx->ctx_entry : ctx->ctx_entry_R;
                uint32_t start;
                uint32_t end;
                tables.template extension_range<key, group_order>(ctx->m, in_group, start, end);

                form = _table_group;
                tbl_len = ctx->m + 1;
                tbl_keys = &tables.template keys<key, group_order>(ctx->m);
                tbl_base = uint64_t(in_group) * tables.alphabet_size();
                cur = start;
                tbl_end = end;
                read_table_symbol();
                return;
            }
        }
    }

    ctx->template update_input_intervals<dir>();

    const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s] = ctx->template const_vars<dir>();
    const auto& idx_dir = idx->template index<dir>();
    cur = 0;

    // starts the enumeration of syms; skips the terminator, which occupies one position of the reverse interval
    auto start_symbol_list = [&]() {
        b_R_nxt = b_R;

        if constexpr (!terminator) {
            if (!syms.empty() && syms[0].sym == 0) {
                b_R_nxt++;
                cur = 1;
            }
        }
    };

    const uint64_t width = e_ - b_ + 1;
    const auto& L_prev_min = idx_dir.L_prev_min();

    /* collect the first occurrence in L'[b_..e_] of each symbol with L_prev_min (x is a first occurrence
       <=> the previous occurrence of L'[x] lies before b_); abort after more than cap symbols
       (cap = 0 <=> L_prev_min is not used) */
    const uint64_t cap = L_prev_min.empty() ? 0 : L_prev_min_max_symbols(idx->sigma, width);

    if (cap != 0) {
        take_buffer(syms);
        syms.clear();

        uint64_t x = L_prev_min.nsv(b_, b_ + 1);

        while (x <= e_ && syms.size() <= cap) {
            syms.emplace_back(symbol_runs_t { idx_dir.L_(x), x, x });

            if (x == e_) {
                x = L_prev_min.none;
                break;
            }

            x = L_prev_min.nsv(x + 1, b_ + 1);
        }

        const bool complete = x > e_ && syms.size() <= cap;

        if (complete) {
            form = _symbol_list;
            std::sort(syms.begin(), syms.end(),
                [](const symbol_runs_t& a, const symbol_runs_t& b) {
                    return a.sym < b.sym;
                });

            // compute the last occurrence in L'[b_..e_] of each symbol
            const uint64_t sym_at_end = idx_dir.L_(e_);

            for (symbol_runs_t& entry : syms) {
                entry.last = entry.sym == sym_at_end ? e_
                    : idx_dir.PSO_L_().pred(entry.sym, e_, b_);
            }

            start_symbol_list();
            return;
        }

        syms.clear();
    }

    // a range of more than scan_width input intervals is processed with prev and next, a smaller one is scanned
    const uint64_t scan_width = max_width_scan_extensions(idx->sigma);
    bool use_prev_next = idx->has_PSB_L_() && width > scan_width;

    if (use_prev_next) {
        form = _prev_next;
        resize_prev_next_buf();
        int64_t* prev = prev_arr();
        int64_t* next = next_arr();
        ctx->template build_prev_next<dir>(prev, next, idx->sigma - 1);

        int64_t num_intervals = idx_dir.M_LF().num_intervals();
        // true <=> the terminator occurs in L'[b_..e_]
        bool has_terminator = next[0] <= prev[0] && prev[0] < num_intervals;

        if constexpr (terminator) {
            b_R_nxt = b_R;
            if (has_terminator) return;
        } else {
            // the terminator occupies one position of the reverse interval
            b_R_nxt = b_R + has_terminator;
        }
    } else {
        // store one entry per input interval in [b_, e_], sort them by symbol and merge the entries of each symbol
        form = _symbol_list;
        take_buffer(syms);
        syms.clear();

        for (uint64_t x = b_; x <= e_; x++) {
            syms.emplace_back(symbol_runs_t { idx_dir.L_(x), x, x });
        }

        std::sort(syms.begin(), syms.end(),
            [](const symbol_runs_t& a, const symbol_runs_t& b) {
                return a.sym != b.sym ? a.sym < b.sym : a.first < b.first;
            });

        uint64_t num_syms = 0;

        for (uint64_t i = 0; i < syms.size(); i++) {
            if (num_syms > 0 && syms[num_syms - 1].sym == syms[i].sym) {
                syms[num_syms - 1].last = syms[i].first;
            } else {
                syms[num_syms++] = syms[i];
            }
        }

        syms.resize(num_syms);
        start_symbol_list();
        return;
    }

    this->template advance_symbol<dir>();
}

template <move_r_support support>
template <query_support_t query_support>
auto move_rb<support>::extend_context_t<query_support>::extend_next() -> search_context_t<query_support>
{
    return dir == LEFT ? extend_next<LEFT>() : extend_next<RIGHT>();
}

template <move_r_support support>
template <query_support_t query_support>
template <direction_t dir>
auto move_rb<support>::extend_context_t<query_support>::extend_next() -> search_context_t<query_support>
{
    // load the context of the extension from the context tables
    if (form == _table_group) {
        constexpr context_extension_t key = ctx_key<dir>;
        constexpr context_order_t group_order = ctx_group_order<dir>;
        const auto& tables = idx->context_tables();
        const uint64_t sym = next_internal_symbol();

        // cur is the entry in group_order, the entry in the other order is looked up with its keys
        uint32_t entry;
        uint32_t entry_R;

        if constexpr (group_order == _forward_order) {
            entry = uint32_t(cur);
            entry_R = tables.template extend_occurring<key, _reverse_order>(tbl_len - 1, ctx->ctx_entry_R, sym);
        } else {
            entry_R = uint32_t(cur);
            entry = tables.template extend_occurring<key, _forward_order>(tbl_len - 1, ctx->ctx_entry, sym);
        }

        search_context_t<query_support> new_ctx;
        idx->template load_context<query_support>(new_ctx, tbl_len, entry, entry_R);
        new_ctx.sym_lst = idx->unmap_symbol(sym);
        cur++;
        read_table_symbol();

        return new_ctx;
    }

    // extend with the symbol of the input intervals next_first_run()..next_last_run() (see extend_to_runs)
    const uint64_t i_sym = next_internal_symbol();
    search_context_t<query_support> new_ctx;
    ctx->template extend_to_runs<dir>(new_ctx, next_first_run(), next_last_run(), b_R_nxt);
    new_ctx.sym_lst = idx->unmap_symbol(i_sym);

    // the interval of the other direction of the next extension starts after the one of this extension
    b_R_nxt += new_ctx.num_occ();
    this->template advance_symbol<dir>();

    return new_ctx;
}

template <move_r_support support>
template <query_support_t query_support>
template <direction_t dir>
void move_rb<support>::extend_context_t<query_support>::advance_symbol()
{
    if (form == _symbol_list) {
        cur++;
        return;
    }

    int64_t num_intervals = idx->template index<dir>().M_LF().num_intervals();
    const int64_t* prev = prev_arr();
    const int64_t* next = next_arr();

    // advance cur to the next symbol that occurs in the range
    do {
        cur++;
    } while (cur < idx->sigma && (
             next[cur] > prev[cur] ||
             prev[cur] == num_intervals));
}

#endif // MOVE_R_MOVE_RB_EXTEND_CONTEXT_TPP
