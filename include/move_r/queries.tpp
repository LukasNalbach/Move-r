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

#ifndef MOVE_R_MOVE_R_QUERIES_TPP
#define MOVE_R_MOVE_R_QUERIES_TPP

#include <cmath>

#include <move_r/move_r.hpp>

// ############################# INTERNAL METHODS #############################

template <move_r_support support, move_pos_encoding_t mlf_enc>
inline void move_r<support, mlf_enc>::init_phi_inv_at_run_start(uint64_t x, uint64_t& s, uint64_t& s_) const
    requires(has_locate_move)
{
    if constexpr (support == _locate_move) {
        // the input interval of M_Phi^{-1} whose output interval starts with s = SA[M_LF.p[x]]
        uint64_t x_s_ = SA_Phi_inv(x);

        // s_ = the input interval of M_Phi^{-1} containing s
        s_ = M_Phi_inv().idx(x_s_);
        s = M_Phi_inv().p(s_) + M_Phi_inv().offs(x_s_);
    } else {
        // if the samples are stored in SA_Phi^{-1}, compute s and s_ as for _locate_move
        if (!_SA_Phi_inv.empty()) {
            uint64_t x_s_ = _SA_Phi_inv.template get<0, uint64_t>(x);
            s_ = M_Phi_inv().idx(x_s_);
            s = M_Phi_inv().p(s_) + M_Phi_inv().offs(x_s_);
            return;
        }

        // else, find s_ with a binary search over the starting positions of the input intervals of M_Phi^{-1}
        s = SA_run_start(x);
        s_ = bin_search_max_leq(s, 0, r__ - 1, [&](uint64_t y) { return M_Phi_inv().p(y); });
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
inline void move_r<support, mlf_enc>::init_phi_inv(
    uint64_t b, uint64_t e,
    uint64_t& s, uint64_t& s_,
    uint64_t hat_b_ap_y, int64_t y) const
    requires(has_locate_move)
{
    // SA[b] = SA[M_LF.p[\hat{b}'_y]] - (y + 1)
    init_phi_inv_at_run_start(hat_b_ap_y, s, s_);
    s -= y + 1;

    // if b < e and s < M_Phi^{-1}.p[s_], an input interval before the s_-th one contains s; find it with an
    // exponential search to the left
    if (b < e && s < M_Phi_inv().p(s_)) {
        s_ = exp_search_max_leq<LEFT>(s, 0, s_, [&](uint64_t x) { return M_Phi_inv().p(x); });
    }
}

// ############################# QUERY METHODS #############################

template <move_r_support support, move_pos_encoding_t mlf_enc>
inline uint64_t move_r<support, mlf_enc>::BWT(uint64_t i) const
{
    // L[i] is the symbol of the input interval of M_LF containing i
    return unmap_symbol(L_(_M_LF.interval_and_offset(i).first));
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
template <typename fnc_t>
void move_r<support, mlf_enc>::SA_range(uint64_t b, uint64_t e, fnc_t report) const
    requires(has_SA_access)
{
    // SA[b] by random access (see SA), then SA[b + 1..e] decoded from it
    const uint64_t first = SA(b);
    report(first);
    if (b == e) return;
    auto dec = _SA_rlz.decode_context();
    dec.init_right(b + 1);
    dec.set_value(first);
    dec.report_right(e, [&](uint64_t, uint64_t value) { report(value); });
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
inline uint64_t move_r<support, mlf_enc>::SA(uint64_t i) const
    requires(has_SA_access)
{
    // index of the first input interval of M_LF in the bwt run containing i (see SA_run_start)
    uint64_t x = first_interval_of_run(_M_LF.interval_and_offset(i).first);

    // position in the suffix array of the current suffix s
    uint64_t j = M_LF().p(x);

    // the current suffix (s = SA[j])
    uint64_t s = SA_run_start(x);

    if (j == i) return s;
    j++;

    // initialize the rlz decode context to position j (with s = SA[j-1])
    auto dec = _SA_rlz.decode_context();
    dec.init_right(j);
    dec.set_value(s);

    // compute SA[i]
    dec.skip_right(i + 1);

    return dec.value();
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
inline void move_r<support, mlf_enc>::query_context_t::set_pos(uint64_t i)
{
    if constexpr (has_SA_rlz) {
        dec = idx->_SA_rlz.decode_context();
        dec.set_pos(i);
    } else if constexpr (has_locate_move) {
        phi.i = i;
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
inline void move_r<support, mlf_enc>::query_context_t::ensure_be()
{
    if (be_valid) return;
    // for M_LF in EF, turn the offsets b,e within their input intervals b_,e_ into positions
    if constexpr (mlf_enc == EF) {
        b = idx->M_LF().pos(b_, b);
        e = idx->M_LF().pos(e_, e);
    }
    set_pos(b);
    be_valid = true;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
template <typename sym_t>
bool move_r<support, mlf_enc>::query_context_t::prepend(sym_t sym)
{
    // for M_LF in EF, turn the positions b,e back into offsets within their input intervals b_,e_
    if constexpr (mlf_enc == EF) {
        if (be_valid) {
            b -= idx->M_LF().p(b_);
            e -= idx->M_LF().p(e_);
            be_valid = false;
        }
    }

    query_context_t ctx_old = *this;

    if (idx->backward_search_step(sym, b, b_, e, e_, hat_b_ap_y, y, hat_e_ap_z, z)) {
        l++;
        be_valid = false; // b,e reconstructed lazily by ensure_be()
        return true;
    } else {
        *this = ctx_old;
        return false;
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
inline uint64_t move_r<support, mlf_enc>::query_context_t::next_occ()
    requires(supports_multiple_locate)
{
    ensure_be(); // reconstructs b,e and seeds the locate cursor on the first call after a prepend

    if constexpr (has_SA_rlz) {
        if (dec.pos() == b) [[unlikely]] {
            // compute the suffix array value at b
            dec.set_value(idx->SA_run_start(hat_b_ap_y) - (y + 1));
            dec.set_pos(b + 1);

            // check if there is more than one occurrence
            if (b < e) dec.init_right(dec.pos());
            return dec.value();
        } else {
            return dec.next();
        }
    } else {
        if (phi.i == b) [[unlikely]] {
            // compute the suffix array value at b
            idx->init_phi_inv(b, e, phi.s, phi.s_, hat_b_ap_y, y);
        } else {
            idx->M_Phi_inv().move(phi.s, phi.s_);
        }

        phi.i++;
        return phi.s;
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
inline uint64_t move_r<support, mlf_enc>::query_context_t::one_occ() const
    requires(supports_locate) { return idx->SA_run_start(hat_b_ap_y) - (y + 1); }

template <move_r_support support, move_pos_encoding_t mlf_enc>
std::vector<uint48_t> move_r<support, mlf_enc>::query_context_t::locate() requires(supports_multiple_locate)
{
    std::vector<uint48_t> Occ;
    Occ.reserve(num_occ_rem());
    locate([&](uint64_t, uint64_t occ){Occ.emplace_back(occ);});
    return Occ;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
template <typename report_fnc_t>
void move_r<support, mlf_enc>::query_context_t::locate(report_fnc_t report)
    requires(supports_multiple_locate)
{
    ensure_be(); // reconstructs b,e and seeds the locate cursor on the first call after a prepend

    if constexpr (has_SA_rlz) {
        if (dec.pos() == b) [[unlikely]] {
            // compute the suffix array value at b
            dec.set_value(idx->SA_run_start(hat_b_ap_y) - (y + 1));
            report(b, dec.value());
            dec.set_pos(b + 1);

            // check if there is more than one occurrence
            if (b < e) [[likely]] dec.init_right(b + 1);
        }

        // compute the remaining occurrences SA(b,e]
        if (dec.pos() <= e) [[likely]] {
            dec.report_right(e, report);
        }
    } else if constexpr (has_locate_move) {
        // compute the suffix array value at b
        if (phi.i == b) [[unlikely]] {
            idx->init_phi_inv(b, e, phi.s, phi.s_, hat_b_ap_y, y);
            report(phi.i, phi.s);
            phi.i++;
        }

        // compute the remaining occurrences SA(b,e]
        while (phi.i <= e) {
            idx->M_Phi_inv().move(phi.s, phi.s_);
            report(phi.i, phi.s);
            phi.i++;
        }
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
inline void move_r<support, mlf_enc>::init_backward_search(
    uint64_t& b, uint64_t& b_,
    uint64_t& e, uint64_t& e_,
    uint64_t& hat_b_ap_y, int64_t& y,
    uint64_t& hat_e_ap_z, int64_t& z) const
{
    // the full suffix array interval [0, n-1]
    b_ = 0;
    b = 0;
    e_ = r_ - 1;
    if constexpr (mlf_enc == POS) e = M_LF().p(r_) - 1;
    else e = M_LF().len(r_ - 1) - 1;
    hat_b_ap_y = 0;
    y = -1;
    hat_e_ap_z = r_ - 1;
    z = -1;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
inline bool move_r<support, mlf_enc>::backward_search_step_mapped(
    uint64_t i_sym,
    uint64_t& b, uint64_t& b_,
    uint64_t& e, uint64_t& e_,
    uint64_t& hat_b_ap_y, int64_t& y,
    uint64_t& hat_e_ap_z, int64_t& z) const
{
    // if sym does not occur in L', then P[i..m] does not occur in T
    if (i_sym == 0) [[unlikely]] return false;

    // b_ <- the first input interval x >= b_ of M_LF with L'[x] = P[i]; b <- its first position (offset 0)
    if (i_sym != L_(b_)) {
        b_ = PSO_L_().succ(i_sym, b_, e_);
        if (b_ > e_) [[unlikely]] return false;
        if constexpr (mlf_enc == POS) b = M_LF().p(b_);
        else b = 0;
        hat_b_ap_y = b_;
        y = 0;
    } else {
        y++;
    }

    // e_ <- the last input interval x <= e_ of M_LF with L'[x] = P[i]; e <- its last position (offset len(e_) - 1)
    if (i_sym != L_(e_)) {
        e_ = PSO_L_().pred(i_sym, e_, b_);
        if constexpr (mlf_enc == POS) e = M_LF().p(e_ + 1) - 1;
        else e = M_LF().len(e_) - 1;
        hat_e_ap_z = e_;
        z = 0;
    } else {
        z++;
    }

    // empty interval (b_ > e_, or b_ = e_ and b > e): P[i..m] does not occur in T
    if (b_ > e_ || (b_ == e_ && b > e)) [[unlikely]] return false;

    /* all suffixes in [b, e] start with P[i+1..m] and are preceded by P[i] in T, so the suffixes in [LF(b), LF(e)]
       are the suffixes starting with P[i..m]; b <- LF(b), e <- LF(e) */
    if (b_ == e_) {
        if (b == e) {
            LF_step(b, b_);
            e_ = b_; e = b;
        } else {
            // b and e are in the same input interval: move b, then advance e by e - b (O(a) time by a-balancedness)
            uint64_t diff_eb = e - b;
            LF_step(b, b_);
            e_ = b_; e = b + diff_eb;
            if constexpr (mlf_enc == POS) {
                while (e >= M_LF().p(e_ + 1)) e_++;
            } else {
                uint64_t l;
                while (e >= (l = M_LF().len(e_))) { e -= l; e_++; }
            }
        }
    } else {
        LF_step(b, b_);
        LF_step(e, e_);
    }

    return true;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
template <typename sym_t>
uint64_t move_r<support, mlf_enc>::count(move_r_pattern_view<sym_t> P) const
{
    uint64_t b_, b, e_, e, hat_b_ap_y, hat_e_ap_z;
    int64_t y, z;
    init_backward_search(b, b_, e, e_, hat_b_ap_y, y, hat_e_ap_z, z);
    int64_t i = int64_t(P.size()) - 1;

    for (; i >= 0; i--) {
        if (!backward_search_step(P[i], b, b_, e, e_, hat_b_ap_y, y, hat_e_ap_z, z)) {
            return 0;
        }
    }

    // e - b + 1 (M_LF in EF: computed from the input intervals and the offsets, see move_data_structure::distance)
    if constexpr (mlf_enc == POS) return e - b + 1;
    else return M_LF().distance(b_, b, e_, e) + 1;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
template <typename report_fnc_t, typename sym_t>
void move_r<support, mlf_enc>::locate(move_r_pattern_view<sym_t> P, report_fnc_t report) const
    requires(supports_multiple_locate)
{
    uint64_t b_, b, e_, e, hat_b_ap_y, hat_e_ap_z;
    int64_t y, z;
    init_backward_search(b, b_, e, e_, hat_b_ap_y, y, hat_e_ap_z, z);
    int64_t i = int64_t(P.size()) - 1;

    for (; i >= 0; i--) {
        if (!backward_search_step(P[i], b, b_, e, e_, hat_b_ap_y, y, hat_e_ap_z, z)) {
            return;
        }
    }

    // with M_LF in EF, compute the positions b and e from their input intervals and offsets
    if constexpr (mlf_enc == EF) { b = M_LF().pos(b_, b); e = M_LF().pos(e_, e); }

    if constexpr (has_SA_rlz) {
        uint64_t s = SA_run_start(hat_b_ap_y) - (y + 1);
        report(b, s);

        if (b < e) {
            auto dec = _SA_rlz.decode_context();
            dec.init_right(b + 1);
            dec.set_value(s);
            dec.report_right(e, report);
        }
    } else if constexpr (has_locate_move) {
        uint64_t s, s_;
        init_phi_inv(b, e, s, s_, hat_b_ap_y, y);
        report(b, s);

        if (b < e) {
            uint64_t i = b + 1;

            while (i <= e) {
                M_Phi_inv().move(s, s_);
                report(i, s);
                i++;
            }
        }
    }
}

// ############################# RETRIEVE-RANGE METHODS #############################

template <move_r_support support, move_pos_encoding_t mlf_enc>
inline void move_r<support, mlf_enc>::adjust_retrieve_params(retrieve_params& params, uint64_t range_max)
{
    if (params.l > params.r) {
        params.l = 0;
        params.r = range_max;
    }

    params.r = std::min(params.r, range_max);
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
template <typename fnc_t>
inline void move_r<support, mlf_enc>::run_parallel_threads(uint64_t p, fnc_t&& fnc)
{
    if (p > 1) {
        #pragma omp parallel num_threads(p)
        fnc(uint64_t(omp_get_thread_num()));
    } else {
        fnc(uint64_t(0));
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
inline uint64_t move_r<support, mlf_enc>::thread_range_start(uint64_t i_p, uint64_t p, uint64_t l, uint64_t r, uint64_t align)
{
    if (i_p == 0) return l;
    if (i_p >= p) return r + 1;
    uint64_t off = (uint64_t(i_p) * (r - l + 1)) / p; // offset from l before the alignment
    if (align > 1) off = div_ceil<uint64_t>(off, align) * align;
    return std::min<uint64_t>(l + off, r + 1);
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
template <typename report_fnc_t>
void move_r<support, mlf_enc>::revert(report_fnc_t report) const
    requires(std::invocable<report_fnc_t&, uint64_t, uint64_t>)
{
    if (n < 2) return;

    // start at position 0 in the bwt (L[0] = T[n-2]), which is the starting position of the 0-th input
    // interval of M_LF; each LF-move query leads to the previous character of T
    uint64_t x = 0; // index of the input interval in M_LF containing the current position in the bwt
    uint64_t d = 0; // the current position in the bwt (M_LF in POS) or its offset in the x-th input interval (EF)

    for (uint64_t j = n - 2;; j--) {
        report(j, unmap_symbol(L_(x)));
        if (j == 0) break;
        LF_step(d, x);
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
template <typename piece_fnc_t, typename report_fnc_t>
void move_r<support, mlf_enc>::revert_strands(uint64_t num_pieces, uint64_t num_threads, uint64_t num_strands,
    piece_fnc_t piece, report_fnc_t report) const
{
    if (num_strands == 1) {
        #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 1)
        for (uint64_t k = 0; k < num_pieces; k++) {
            revert_strand s;
            piece(k, s);
            const uint64_t w = omp_get_thread_num();
            uint64_t x = s.x;
            uint64_t d = s.d;

            for (uint64_t j = s.j;; j--) {
                report(w, j, unmap_symbol(L_(x)));
                if (j == s.first) break;
                LF_step(d, x);
            }
        }

        return;
    }

    std::atomic<uint64_t> next_piece = 0; // the next piece not taken by any thread

    #pragma omp parallel num_threads(num_threads)
    {
        // [0..num_strands-1] the pieces of this thread; [0..active-1] are not finished
        no_init_vector<revert_strand> strands;
        resize_for_overwrite(strands, num_strands);
        uint64_t active = 0;

        // sets s to the next piece not taken by any thread; returns false if there is none
        auto take_piece = [&](revert_strand& s) {
            const uint64_t k = next_piece.fetch_add(1, std::memory_order_relaxed);
            if (k >= num_pieces) return false;
            piece(k, s);
            return true;
        };

        while (active < num_strands && take_piece(strands[active])) {
            strands[active].w = omp_get_thread_num() * num_strands + active;
            active++;
        }

        while (active > 0) {
            for (uint64_t t = 0; t < active;) {
                revert_strand& s = strands[t];
                report(s.w, s.j, unmap_symbol(L_(s.x)));

                if (s.j > s.first) {
                    s.j--;
                    LF_step(s.d, s.x);
                } else if (!take_piece(s)) {
                    std::swap(s, strands[--active]);
                    continue;
                }

                t++;
            }
        }
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
template <typename report_fnc_t>
void move_r<support, mlf_enc>::revert_in_pieces(report_fnc_t report, uint64_t num_threads, uint64_t num_strands) const
{
    if (n < 2) return;
    const uint64_t m = _revert_starts.size();

    /* piece k in [0..m] is T[s_{k-1}..s_k - 1], where s_k is the text position of the k-th revert start, s_{-1} = 0
       and s_m = n - 1; it is reverted with LF-steps, starting at the BWT position of the k-th revert start (for k = m
       at position 0, L[0] = T[n-2]) */
    revert_strands(m + 1, revert_num_threads(num_threads), num_strands, [&](uint64_t k, revert_strand& s) {
        s.first = k == 0 ? 0 : _revert_starts.template get<1, uint64_t>(k - 1);
        s.j = (k == m ? n - 1 : _revert_starts.template get<1, uint64_t>(k)) - 1;
        const uint64_t i = k == m ? 0 : _revert_starts.template get<0, uint64_t>(k);
        const auto [x, offset] = _M_LF.interval_and_offset(i);
        s.x = x;
        s.d = mlf_enc == POS ? i : offset;
    }, report);
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
template <typename out_sym_t>
void move_r<support, mlf_enc>::revert(std::string file_name, uint64_t num_threads, uint64_t num_strands) const
{
    direct_positioned_file file(file_name, (n - 1) * sizeof(out_sym_t));

    {
        // [0..p * num_strands - 1] the writer of each strand slot
        std::vector<value_writer_backward<sizeof(out_sym_t)>> writers;
        const uint64_t p = revert_num_threads(num_threads);
        const uint64_t buffer_values = div_ceil<uint64_t>(huge_page_bytes / sizeof(out_sym_t), num_strands);
        writers.reserve(p * num_strands);
        for (uint64_t w = 0; w < p * num_strands; w++) writers.emplace_back(file, 0, buffer_values);
        revert_in_pieces([&](uint64_t w, uint64_t i, uint64_t c) { writers[w].put(i, c); }, p, num_strands);
    }

    if (!file.close_file()) throw std::runtime_error("cannot write " + file_name);
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
template <typename out_t>
out_t move_r<support, mlf_enc>::BWT_range(retrieve_params params) const
{
    adjust_retrieve_params(params, n - 1);
    out_t L;
    resize_for_overwrite(L, params.r - params.l + 1);
    BWT_range([&](uint64_t i, uint64_t c) { L[i - params.l] = terminator_as_zero(c); }, params);
    return L;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
template <typename report_fnc_t>
void move_r<support, mlf_enc>::BWT_range(report_fnc_t report, retrieve_params params) const
{
    adjust_retrieve_params(params, n - 1);

    uint64_t l = params.l;
    uint64_t r = params.r;

    // each thread reports at least one cache line of symbols
    uint64_t p = retrieve_num_threads(params.num_threads, (r - l + 1) / cache_line_bytes);

    run_parallel_threads(p, [&](uint64_t i_p) {
        // the range [b, e] of thread i_p (empty, if the alignment to the output buffer blocks moves b to e)
        uint64_t b = thread_range_start(i_p, p, l, r, params.buffer_align);
        uint64_t e = thread_range_start(i_p + 1, p, l, r, params.buffer_align);
        if (b >= e) return;
        e--;

        // the current position in the bwt
        uint64_t i = b;

        // the input interval x of M_LF containing i, and the offset d of i in it
        auto [x, d] = _M_LF.interval_and_offset(i);

        // the starting position M_LF.p[x + 1] of the next input interval
        uint64_t next_interval_beg = i + (M_LF().len(x) - d);

        // iterate until x is the input interval containing e
        while (next_interval_beg <= e) {
            // iterate over all positions in the x-th input interval
            while (i < next_interval_beg) {
                report(i, unmap_symbol(L_(x)));
                i++;
            }

            x++;
            next_interval_beg += M_LF().len(x);
        }

        // report the remaining characters
        while (i <= e) {
            report(i, unmap_symbol(L_(x)));
            i++;
        }
    });
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
no_init_vector<uint48_t> move_r<support, mlf_enc>::SA_range(retrieve_params params) const requires(supports_multiple_locate)
{
    adjust_retrieve_params(params, n - 1);
    no_init_vector<uint48_t> SA_rng;
    resize_no_init(SA_rng, params.r - params.l + 1);
    SA_range([&](uint64_t i, uint64_t v){SA_rng[i - params.l] = v;}, params);
    return SA_rng;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
template <typename report_fnc_t>
void move_r<support, mlf_enc>::SA_range_in_pieces(report_fnc_t report, retrieve_params params) const
    requires(supports_multiple_locate)
{
    adjust_retrieve_params(params, n - 1);
    const uint64_t l = params.l;
    const uint64_t r = params.r;
    const uint64_t p = SA_range_num_threads(params);
    const uint64_t num_strands = std::max<uint64_t>(1, params.num_strands);

    /* the pieces start at the positions starts[k] in [l, r], in increasing order; starts[0] = l, and every other
       start is the start of a BWT run; runs[k] is the first input interval of M_LF in the BWT run containing
       starts[k], so a piece is decoded from the SA-sample at the start of runs[k] */
    no_init_vector<uint64_t> starts;
    no_init_vector<uint64_t> runs;
    const uint64_t max_pieces = p * sa_pieces_per_thread;
    starts.reserve(max_pieces + 1);
    runs.reserve(max_pieces);

    for (uint64_t k = 0; k < max_pieces; k++) {
        const uint64_t target = l + (k * (r - l + 1)) / max_pieces;
        const uint64_t x = first_interval_of_run(_M_LF.interval_and_offset(target).first);
        const uint64_t start = k == 0 ? l : M_LF().p(x);

        if (k == 0 || start > starts.back()) {
            starts.push_back(start);
            runs.push_back(x);
        }
    }

    const uint64_t num_pieces = runs.size();
    starts.push_back(r + 1);

    if constexpr (has_locate_move) {
        // a piece decoded with Phi^{-1}-steps: s = SA[i], and s_ is the input interval of M_Phi^{-1} containing s
        struct strand_t {
            uint64_t i; // the position of the next value to report
            uint64_t e; // the last position of the piece
            uint64_t s;
            uint64_t s_;
            uint64_t w; // the strand slot
        };

        run_strands<strand_t>(num_pieces, p, num_strands, [&](uint64_t k, strand_t& st) {
            init_phi_inv_at_run_start(runs[k], st.s, st.s_);
            st.e = starts[k + 1] - 1;

            // only the first piece can start after the start of its BWT run
            for (st.i = M_LF().p(runs[k]); st.i < starts[k]; st.i++) M_Phi_inv().move(st.s, st.s_);
        }, [&](strand_t& st) {
            report(st.w, st.i, st.s);
            if (st.i == st.e) return false;
            M_Phi_inv().move(st.s, st.s_);
            st.i++;
            return true;
        });
    } else if constexpr (has_SA_rlz) {
        // a piece decoded with the rlz of the differential suffix array: s = SA[i]
        struct strand_t {
            typename rlz<_diff_access>::decode_context_t dec; // positioned after i (if i < e)
            uint64_t i; // the position of the next value to report
            uint64_t e; // the last position of the piece
            uint64_t s;
            uint64_t w; // the strand slot
        };

        // sets st to the beginning of the k-th piece
        auto start = [&](uint64_t k, strand_t& st) {
            // the start of the BWT run and its SA-sample
            const uint64_t j = M_LF().p(runs[k]);
            st.s = SA_run_start(runs[k]);
            st.i = starts[k];
            st.e = starts[k + 1] - 1;

            if (j < st.e) {
                st.dec = _SA_rlz.decode_context();
                st.dec.init_right(j + 1);
                st.dec.set_value(st.s);

                if (j < st.i) {
                    st.dec.skip_right(st.i);
                    st.s = st.dec.next();
                }
            }
        };

        // with one strand, a thread decodes the values of a piece in one go (see rlz::decode_context_t::report_right)
        if (num_strands == 1) {
            #pragma omp parallel for num_threads(p) schedule(dynamic, 1)
            for (uint64_t k = 0; k < num_pieces; k++) {
                strand_t st;
                start(k, st);
                const uint64_t w = omp_get_thread_num();
                report(w, st.i, st.s);
                if (st.i < st.e) st.dec.report_right(st.e, [&](uint64_t i, uint64_t s) { report(w, i, s); });
            }

            return;
        }

        run_strands<strand_t>(num_pieces, p, num_strands, start, [&](strand_t& st) {
            report(st.w, st.i, st.s);
            if (st.i == st.e) return false;
            st.s = st.dec.next();
            st.i++;
            return true;
        });
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::SA_range(std::string file_name, retrieve_params params) const
    requires(supports_multiple_locate)
{
    adjust_retrieve_params(params, n - 1);
    const uint64_t l = params.l;
    const uint64_t length = params.r - l + 1;
    const uint64_t num_slots = SA_range_strand_slots(params);
    // the output buffers take about 1 / retrieve_output_bytes_per_buffer_byte of the output, or max_bytes_alloc bytes
    const uint64_t buffer_bytes = std::max<uint64_t>(retrieve_min_buffer_bytes, params.max_bytes_alloc != -1
        ? params.max_bytes_alloc / num_slots
        : (length * sizeof(uint64_t)) / (num_slots * retrieve_output_bytes_per_buffer_byte));
    direct_positioned_file file(file_name, length * sizeof(uint64_t));

    {
        // [0..num_slots-1] the writer of each strand slot
        std::vector<value_writer<sizeof(uint64_t)>> writers;
        writers.reserve(num_slots);
        for (uint64_t w = 0; w < num_slots; w++) writers.emplace_back(file, 0, 0, buffer_bytes / sizeof(uint64_t));
        SA_range_in_pieces([&](uint64_t w, uint64_t i, uint64_t s) { writers[w].put(i - l, s); }, params);
    }

    if (!file.close_file()) throw std::runtime_error("cannot write " + file_name);
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
template <typename output_t, typename invoke_t>
void move_r<support, mlf_enc>::retrieve_range(
    invoke_t retrieve, std::string file_name, move_r<support, mlf_enc>::retrieve_params params) const
{
    uint64_t l = params.l;
    uint64_t r = params.r;
    uint64_t num_threads = params.num_threads;
    uint64_t buffer_bytes = std::max<uint64_t>(retrieve_min_buffer_bytes, params.max_bytes_alloc != -1
        ? params.max_bytes_alloc / num_threads
        : ((r - l + 1) * sizeof(output_t)) / (num_threads * retrieve_output_bytes_per_buffer_byte));

    // the ranges of the threads start at multiples of the buffer size, so that no two threads write to the same
    // block of the file
    uint64_t align = std::max<uint64_t>(1, buffer_bytes / sizeof(output_t));
    params.buffer_align = align;

    direct_positioned_file file(file_name, (r - l + 1) * sizeof(output_t));

    {
        // [0..num_threads-1] the writer of each thread
        std::vector<std::unique_ptr<value_writer<sizeof(output_t)>>> writers(num_threads);
        for (auto& writer : writers) writer = std::make_unique<value_writer<sizeof(output_t)>>(file, 0, 0, align);

        // the value at position pos is the (pos - l)-th value of the file
        retrieve([&](uint64_t pos, output_t val) {
            writers[omp_get_thread_num()]->put(pos - l, uint64_t(val));
        }, params);
    }

    if (!file.close_file()) throw std::runtime_error("cannot write " + file_name);
}

#endif // MOVE_R_MOVE_R_QUERIES_TPP
