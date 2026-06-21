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

#include <move_r/move_r.hpp>

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::setup_phi_m1_move_pair(pos_t& x, pos_t& s, pos_t& s_) const
    requires(has_locate_move)
{
    if constexpr (support == _locate_move) {
        // the index of the pair in M_Phi^{-1} creating the output interval with starting position s = SA[M_LF.p[x]]
        pos_t x_s_ = SA_Phi_m1(x);

        // set s_ to the index of the input interval in M_Phi^{-1} containing s
        s_ = M_Phi_m1().idx(x_s_);

        // compute s
        s = M_Phi_m1().p(s_) + M_Phi_m1().offs(x_s_);
    } else {
        s = SA_s(x);
        s_ = bin_search_max_leq<pos_t>(s, 0, r__ - 1, [&](pos_t x){ return M_Phi_m1().p(x); });
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
sym_t move_r<support, sym_t, pos_t>::BWT(pos_t i) const
{
    // find the index of the input interval in M_LF containing i with a binary search.
    return unmap_symbol(L_(bin_search_max_leq<pos_t>(i, 0, r_ - 1, [&](pos_t x) { return M_LF().p(x); })));
}

template <move_r_support support, typename sym_t, typename pos_t>
pos_t move_r<support, sym_t, pos_t>::SA(pos_t i) const
    requires(supports_multiple_locate)
{
    if constexpr (has_rlzsa) {
        // index of the input interval in M_LF containing i.
        pos_t x = bin_search_max_leq<pos_t>(i, 0, r_ - 1, [&](pos_t x_) { return M_LF().p(x_); });
        while (SA_s(x) == n) x--;

        // position in the suffix array of the current suffix s
        pos_t j = M_LF().p(x);

        // the current suffix (s = SA[j])
        pos_t s = SA_s(x);

        if (j == i) return s;
        j++;

        // initialize the rlzsa decode context to position j (with s = SA[j-1])
        auto dec = _rlzsa.decode();
        dec.init(j);
        dec.set_value(s);

        // compute SA[i]
        dec.skip_right(i + 1);

        return dec.value();
    } else if constexpr (has_locate_move) {
        // index of the input interval in M_LF containing i.
        pos_t x = bin_search_max_leq<pos_t>(i, 0, r_ - 1, [&](pos_t x_) { return M_LF().p(x_); });

        if constexpr (support == _locate_move) {
            /* if i is a bwt run end position (i = M_LF.p(x+1)-1) and SA_Phi^{-1}[x+1] != r'', then
                SA[i] = Phi(SA_s[(x+1) mod r'])
                    = Phi(M_Phi^{-1}.q(SA_Phi^{-1}[(x+1) mod r']))
                    =     M_Phi^{-1}.p(SA_Phi^{-1}[(x+1) mod r'])
            */
            if (i == M_LF().p(x + 1) - 1) [[unlikely]] {
                pos_t xp1 = (x + 1) == r_ ? 0 : (x + 1);

                if (SA_Phi_m1(xp1) != r__) [[unlikely]] {
                    return M_Phi_m1().p(SA_Phi_m1(xp1));
                }
            }

            // decrement x until the starting position of the x-th input interval of M_LF is a starting position of a bwt run
            while (SA_Phi_m1(x) == r__) x--;
        }

        // begin iterating at the start of the x-th run, because there is a
        // suffix array sample at the end position of the x-th input interval

        // position in the suffix array of the current suffix s
        pos_t j = M_LF().p(x);

        // index of the input interval in M_Phi^{-1} containing s
        pos_t s_;
        // the current suffix (s = SA[j])
        pos_t s;

        setup_phi_m1_move_pair(x, s, s_);

        // Perform Phi-move queries until s is the suffix at position
        // i; in each iteration, s = SA[j] = \Phi^{i-j}(SA[i]) holds.
        while (j < i) {
            // Set s = \Phi(s)
            M_Phi_m1().move(s, s_);
            j++;
        }

        // Since j = i, now s = SA[i] holds.
        return s;
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
bool move_r<support, sym_t, pos_t>::query_context_t::prepend(sym_t sym)
{
    query_context_t ctx_old = *this;

    if (idx->backward_search_step(sym, b, e, b_, e_, hat_b_ap_y, y, hat_e_ap_z, z)) {
        l++;
        set_pos(b);
        return true;
    } else {
        *this = ctx_old;
        return false;
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
pos_t move_r<support, sym_t, pos_t>::query_context_t::next_occ()
    requires(supports_multiple_locate)
{
    if constexpr (has_rlzsa) {
        if (dec.pos() == b) [[unlikely]] {
            // compute the suffix array value at b
            dec.set_value(idx->SA_s(hat_b_ap_y) - (y + 1));
            dec.set_pos(b + 1);

            // check if there is more than one occurrence
            if (b < e) dec.init(dec.pos());
            return dec.value();
        } else {
            return dec.next();
        }
    } else {
        if (phi.i == b) [[unlikely]] {
            // compute the suffix array value at b
            idx->init_phi_m1(b, e, phi.s, phi.s_, hat_b_ap_y, y);
        } else {
            idx->M_Phi_m1().move(phi.s, phi.s_);
        }

        phi.i++;
        return phi.s;
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
pos_t move_r<support, sym_t, pos_t>::query_context_t::one_occ() const
    requires(supports_locate)
{
    return idx->SA_s(hat_b_ap_y) - (y + 1);
}

template <move_r_support support, typename sym_t, typename pos_t>
template <typename report_fnc_t>
void move_r<support, sym_t, pos_t>::query_context_t::locate(report_fnc_t report)
    requires(supports_multiple_locate)
{
    if constexpr (has_rlzsa) {

        if (dec.pos() == b) [[unlikely]] {
            // compute the suffix array value at b
            dec.set_value(idx->SA_s(hat_b_ap_y) - (y + 1));
            report(dec.value());
            dec.set_pos(b + 1);

            // check if there is more than one occurrence
            if (b < e) [[likely]] dec.init(b + 1);
        }

        // compute the remaining occurrences SA(b,e]
        if (dec.pos() <= e) [[likely]] {
            dec.report_right(e, report);
        }
    } else if constexpr (has_locate_move) {
        // compute the suffix array value at b
        if (phi.i == b) [[unlikely]] {
            idx->init_phi_m1(b, e, phi.s, phi.s_, hat_b_ap_y, y);
            report(phi.s);
            phi.i++;
        }

        // compute the remaining occurrences SA(b,e]
        while (phi.i <= e) {
            idx->M_Phi_m1().move(phi.s, phi.s_);
            report(phi.s);
            phi.i++;
        }
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
bool move_r<support, sym_t, pos_t>::backward_search_step(
    sym_t sym,
    pos_t& b, pos_t& e,
    pos_t& b_, pos_t& e_,
    pos_t& hat_b_ap_y, int64_t& y,
    pos_t& hat_e_ap_z, int64_t& z) const
{
    // If the characters have been remapped internally, the pattern also has to be remapped.
    i_sym_t i_sym = map_symbol(sym);

    // If sym does not occur in L', then P[i..m] does not occur in T
    if (i_sym == 0) [[unlikely]] return false;

    // Find the lexicographically smallest suffix in the current suffix array interval that is prefixed by P[i]
    if (i_sym != L_(b_)) {
        /* To do so, we can at first find the first (sub-)run with character P[i] after the b_-th (sub-)run, save
        its index in b_ and set b to its start position M_LF.p(b_). */

        if constexpr (byte_alphabet) {
            pos_t blk = div_ceil<pos_t>(b_, L_block_size());
            pos_t max_b_ = std::min<pos_t>(blk * L_block_size(), e_);
            while (b_ <= max_b_ && i_sym != L_(b_)) b_++;

            if (b_ > max_b_ && i_sym != L_(b_)) [[likely]] {
                if (b_ > e_) [[unlikely]] return false;
                b_ = L_next(blk, i_sym);
            }
        } else {
            b_ = RS_L_().rank(i_sym, b_);
            if (b_ == RS_L_().frequency(i_sym)) [[unlikely]] return false;
            b_ = RS_L_().select(i_sym, b_ + 1);
        }

        if (b_ > e_) [[unlikely]] return false;
        b = M_LF().p(b_);
        hat_b_ap_y = b_;
        y = 0;
    } else {
        y++;
    }

    // Find the lexicographically largest suffix in the current suffix array interval that is prefixed by P[i]
    if (i_sym != L_(e_)) {
        /* To do so, we can at first find the (sub-)last run with character P[i] before the e_-th (sub-)run, save
        its index in e_ and set e to its end position M_LF.p(e_+1)-1. */

        if constexpr (byte_alphabet) {
            pos_t blk = e_ / L_block_size();
            pos_t min_e_ = std::max<pos_t>(blk * L_block_size(), b_);
            while (e_ >= min_e_ && i_sym != L_(e_)) e_--;

            if (e_ < min_e_ && i_sym != L_(e_)) [[likely]] {
                e_ = L_prev(blk, i_sym);
            }
        } else {
            e_ = RS_L_().select(i_sym, RS_L_().rank(i_sym, e_));
        }

        e = M_LF().p(e_ + 1) - 1;
        hat_e_ap_z = e_;
        z = 0;
    } else {
        z++;
    }

    // Else, because each suffix i in the previous suffix array interval starts with P[i+1..m] and the current
    // interval [b,e] contains all suffixes of it, before which there is a P[i] in T, all suffixes in the
    // interval SA[LF(b),LF(e)] start with P[i..m]

    /* If the suffix array interval [LF(b),LF(e)] of P[i..m] is empty, then b > e,
    because LF(i) is monotonic for a fixed L[i], hence it suffices to check, whether
    b <= e holds. */

    // If the suffix array interval is empty, P does not occur in T, so return false.
    if (b > e) [[unlikely]] return false;

    /* Else, set b <- LF(b) and e <- LF(e). The following two optimizations increase query throughput slightly
        if there are only few occurrences */
    if (b_ == e_) {
        if (b == e) {
            /* If \hat{b'}_i == \hat{e'}_i and b'_i = e'_i, then computing
            (e_i,\hat{e}_i) <- M_LF.move(e'_i,\hat{e'}_i) is redundant */
            M_LF().move(b, b_);
            e = b;
            e_ = b_;
        } else {
            /* If \hat{b'}_i == \hat{e'}_i, but b'_i != e'_i, then e_i = b_i + e'_i - b'_i and therefore
            \hat{b'}_i < \hat{e'}_i, hence we can compute \hat{e'}_i by setting e_ <- \hat{b'}_i = b_ and
            incrementing e_ until e < M_LF.p[e_+1] holds; This takes O(a) time because of the a-balancedness property */
            pos_t diff_eb = e - b;
            M_LF().move(b, b_);
            e = b + diff_eb;
            e_ = b_;
            while (e >= M_LF().p(e_ + 1)) e_++;
        }
    } else {
        M_LF().move(b, b_);
        M_LF().move(e, e_);
    }

    return true;
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::init_phi_m1(
    pos_t& b, pos_t& e,
    pos_t& s, pos_t& s_,
    pos_t& hat_b_ap_y, int64_t& y) const
    requires(has_locate_move)
{
    setup_phi_m1_move_pair(hat_b_ap_y, s, s_);
    s -= y + 1;

    // If there is more than one occurrence and s < M_Phi^{-1}.p[s_], now an input interval of M_Phi^{-1} before
    // the s_-th one contains s, so we have to decrease s_. To find the correct value for s_, we perform
    // an exponential search to the left over the input interval starting positions of M_Phi^{-1} starting at s_.
    if (b < e && s < M_Phi_m1().p(s_)) {
        s_ = exp_search_max_leq<pos_t, LEFT>(s, 0, s_, [&](pos_t x) { return M_Phi_m1().p(x); });
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
pos_t move_r<support, sym_t, pos_t>::count(const inp_t& P) const
{
    pos_t b, e, b_, e_, hat_b_ap_y, hat_e_ap_z;
    int64_t y, z;
    init_backward_search(b, e, b_, e_, hat_b_ap_y, y, hat_e_ap_z, z);

    for (int64_t i = P.size() - 1; i >= 0; i--) {
        if (!backward_search_step(P[i], b, e, b_, e_, hat_b_ap_y, y, hat_e_ap_z, z)) {
            return 0;
        }
    }

    return e - b + 1;
}

template <move_r_support support, typename sym_t, typename pos_t>
template <typename report_fnc_t>
void move_r<support, sym_t, pos_t>::locate(const inp_t& P, report_fnc_t report) const
    requires(supports_multiple_locate)
{
    pos_t b, e, b_, e_, hat_b_ap_y, hat_e_ap_z;
    int64_t y, z;
    init_backward_search(b, e, b_, e_, hat_b_ap_y, y, hat_e_ap_z, z);

    for (int64_t i = P.size() - 1; i >= 0; i--) {
        if (!backward_search_step(P[i], b, e, b_, e_, hat_b_ap_y, y, hat_e_ap_z, z)) {
            return;
        }
    }

    if constexpr (has_rlzsa) {
        pos_t s = SA_s(hat_b_ap_y) - (y + 1);
        report(s);

        if (b < e) {
            auto dec = _rlzsa.decode();
            dec.init(b + 1);
            dec.set_value(s);
            dec.report_right(e, report);
        }
    } else if constexpr (has_locate_move) {
        pos_t s, s_;
        init_phi_m1(b, e, s, s_, hat_b_ap_y, y);
        report(s);

        if (b < e) {
            pos_t i = b + 1;

            while (i <= e) {
                M_Phi_m1().move(s, s_);
                report(s);
                i++;
            }
        }
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
template <typename report_fnc_t>
void move_r<support, sym_t, pos_t>::revert_range(report_fnc_t report, retrieve_params params) const
{
    adjust_retrieve_params(params, n - 2);

    pos_t l = params.l;
    pos_t r = params.r;

    // leftmost section to revert
    uint16_t s_l;
    // rightmost section to revert
    uint16_t s_r;

    if (p_r == 1) {
        s_l = 0;
        s_r = 0;
    } else {
        s_l = bin_search_min_gt<pos_t>(l, 0, p_r - 1, [&](pos_t x) { return _D_e[x].second; });
        s_r = bin_search_min_geq<pos_t>(r, 0, p_r - 1, [&](pos_t x) { return _D_e[x].second; });
    }

    uint16_t p = std::max(
        (uint16_t)1, // use at least one thread
        std::min({
            (uint16_t)(s_r - s_l + 1), // use at most s_r-s_l+1 threads
            (uint16_t)omp_get_max_threads(), // use at most all threads
            params.num_threads // use at most the specified number of threads
        }));
    
    if (r == n - 1) {
        report(r, 0);
    }

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // leftmost section for thread i_p to revert
        uint16_t sl_ip = s_l + (i_p * (s_r - s_l + 1)) / p;
        // rightmost section for thread i_p to revert
        uint16_t sr_ip = i_p == p - 1 ? s_r : s_l + ((i_p + 1) * (s_r - s_l + 1)) / p - 1;

        // Iteration range start position of thread i_p.
        pos_t j_l = std::max(l, sl_ip == 0 ? 0 : (_D_e[sl_ip - 1].second + 1) % n);
        // Iteration range end position of thread i_p.
        pos_t j_r = sr_ip == p_r - 1 ? n - 2 : _D_e[sr_ip].second;

        // index of the input interval in M_LF containing i.
        pos_t x = sr_ip == p_r - 1 ? 0 : _D_e[sr_ip].first;
        // The position in the bwt of the current character in T.
        pos_t i = sr_ip == p_r - 1 ? 0 : M_LF().p(x);

        // start iterating at the right iteration range end position
        pos_t j = j_r;

        // iterate until j = r
        while (j > r) {
            // Set i <- LF(i) and j <- j-1.
            M_LF().move(i, x);
            j--;
        }

        if (j >= l) {
            // Report T[r] = T[j] = L[i] = L'[x]
            report(j, unmap_symbol(L_(x)));
        }

        // report T[l,r-1] from right to left
        while (j > j_l) {
            // Set i <- LF(i) and j <- j-1.
            M_LF().move(i, x);
            j--;
            // Report T[j] = L[i] = L'[x].
            report(j, unmap_symbol(L_(x)));
        }
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
template <typename report_fnc_t>
void move_r<support, sym_t, pos_t>::BWT_range(report_fnc_t report, retrieve_params params) const
{
    adjust_retrieve_params(params, n - 1);

    pos_t l = params.l;
    pos_t r = params.r;

    uint16_t p = std::max(
        (uint16_t)1, // use at least one thread
        std::min({
            (uint16_t)omp_get_max_threads(), // use at most all threads
            (uint16_t)((r - l + 1) / 10), // use at most (r-l+1)/100 threads
            params.num_threads // use at most the specified number of threads
        }));

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Iteration range start position of thread i_p.
        pos_t b = l + i_p * ((r - l + 1) / p);
        // Iteration range end position of thread i_p.
        pos_t e = i_p == p - 1 ? r : l + (i_p + 1) * ((r - l + 1) / p) - 1;

        // Current position in the bwt.
        pos_t i = b;

        // index of the input interval in M_LF containing i.
        pos_t x = bin_search_max_leq<pos_t>(i, 0, r_ - 1, [&](pos_t x_) { return M_LF().p(x_); });

        // start position of the next input interval in M_LF
        pos_t l_xp1;
 
        // iterate until x is the input interval containing e
        while ((l_xp1 = M_LF().p(x + 1)) <= e) {

            // iterate over all positions in the x-th input interval
            while (i < l_xp1) {
                report(i, unmap_symbol(L_(x)));
                i++;
            }

            x++;
        }

        // report the remaining characters
        while (i <= e) {
            report(i, unmap_symbol(L_(x)));
            i++;
        }
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
template <typename report_fnc_t>
void move_r<support, sym_t, pos_t>::SA_range(report_fnc_t report, retrieve_params params) const
    requires(supports_multiple_locate)
{
    adjust_retrieve_params(params, n - 1);

    pos_t l = params.l;
    pos_t r = params.r;

    uint16_t p = std::max(
        (uint16_t)1, // use at least one thread
        std::min({
            (uint16_t)omp_get_max_threads(), // use at most all threads
            params.num_threads, // use at most the specified number of threads
            (uint16_t)(((r - l + 1) * (double)r__) / (10.0 * (double)n)) // use at most (r-l+1)*(r/n)*(1/10) threads
        }));

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // iteration range start position
        pos_t b = l + i_p * ((r - l + 1) / p);
        // iteration range end position
        pos_t e = i_p == p - 1 ? r : l + (i_p + 1) * ((r - l + 1) / p) - 1;

        if constexpr (has_locate_move) {
            // the input interval of M_LF containing i
            pos_t x = bin_search_max_leq<pos_t>(b, 0, r_ - 1, [&](pos_t x_) { return M_LF().p(x_); });

            if constexpr (support == _locate_move) {
                // decrement x until the starting position of the x-th input interval of M_LF is a starting position of a bwt run
                while (SA_Phi_m1(x) == r__) {
                    x--;
                }
            }

            // current position in the suffix array, initially the starting position of the x-th interval of M_LF
            pos_t i = M_LF().p(x);

            // index of the input interval in M_Phi^{-1} containing s
            pos_t s_;
            /* the current suffix array value (SA[i]), initially the suffix array sample of the x-th run,
            initially the suffix array value at b */
            pos_t s;

            setup_phi_m1_move_pair(x, s, s_);

            // iterate up to the iteration range starting position
            while (i < b) {
                M_Phi_m1().move(s, s_);
                i++;
            }

            // report SA[b]
            report(i, s);

            // report the SA-values SA[b+1,e] from left to right
            while (i < e) {
                M_Phi_m1().move(s, s_);
                i++;
                report(i, s);
            }
        } else if constexpr (has_rlzsa) {
            // index of the input interval in M_LF containing i.
            pos_t x = bin_search_max_leq<pos_t>(b, 0, r_ - 1, [&](pos_t x_) { return M_LF().p(x_); });

            // position in the suffix array of the current suffix s
            pos_t j = M_LF().p(x);

            // the current suffix (s = SA[j])
            pos_t s = SA_s(x);
            if (j == b) report(b, s);

            if (j < e) [[likely]] {
                j++;

                // initialize the rlzsa decode context to position j (with s = SA[j-1])
                auto dec = _rlzsa.decode();
                dec.init(j);
                dec.set_value(s);

                // advance the rlzsa context to position b
                dec.skip_right(b);

                // decode and report SA[b..e]
                dec.report_right(e, report);
            }
        }
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
template <typename output_t, bool output_reversed>
void move_r<support, sym_t, pos_t>::retrieve_range(
    void (move_r<support, sym_t, pos_t>::*retrieve_method)(
        std::function<void(pos_t, output_t)>&&, move_r<support, sym_t, pos_t>::retrieve_params) const,
    std::string file_name, move_r<support, sym_t, pos_t>::retrieve_params params) const
{
    pos_t l = params.l;
    pos_t r = params.r;
    uint16_t num_threads = params.num_threads;
    uint64_t buffer_size_per_thread = std::max<uint64_t>(1024,
        params.max_bytes_alloc != -1 ? params.max_bytes_alloc / num_threads : ((r - l + 1) * sizeof(output_t)) / (num_threads * 500));

    std::filesystem::resize_file(std::filesystem::current_path() / file_name, (r - l + 1) * sizeof(output_t));
    std::vector<sdsl::int_vector_buffer<sizeof(output_t) * 8>> file_bufs;

    for (uint16_t i = 0; i < num_threads; i++) {
        file_bufs.emplace_back(sdsl::int_vector_buffer<sizeof(output_t) * 8>(file_name, std::ios::in, buffer_size_per_thread, sizeof(output_t) * 8, true));
    }
    
    (this->*retrieve_method)([&](pos_t pos, output_t val) {
        file_bufs[omp_get_thread_num()][pos] = *reinterpret_cast<uint64_t*>(&val);
    }, params);
}