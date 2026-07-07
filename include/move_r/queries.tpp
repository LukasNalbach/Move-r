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

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
void move_r<support, sym_t, pos_t, mlf_enc>::setup_phi_m1_move_pair(pos_t& x, pos_t& s, pos_t& s_) const
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

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
sym_t move_r<support, sym_t, pos_t, mlf_enc>::BWT(pos_t i) const
{
    // find the index of the input interval in M_LF containing i with a binary search.
    return unmap_symbol(L_(M_LF().interval_index(i)));
}

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
pos_t move_r<support, sym_t, pos_t, mlf_enc>::SA(pos_t i) const
    requires(supports_multiple_locate)
{
    if constexpr (has_rlzsa) {
        // index of the input interval in M_LF containing i.
        pos_t x = M_LF().interval_index(i);
        while (SA_s(x) == n) x--;

        // position in the suffix array of the current suffix s
        pos_t j = M_LF().p(x);

        // the current suffix (s = SA[j])
        pos_t s = SA_s(x);

        if (j == i) return s;
        j++;

        // initialize the rlzsa decode context to position j (with s = SA[j-1])
        auto dec = _rlzsa.decode();
        dec.init_right(j);
        dec.set_value(s);

        // compute SA[i]
        dec.skip_right(i + 1);

        return dec.value();
    } else if constexpr (has_locate_move) {
        // index of the input interval in M_LF containing i.
        pos_t x = M_LF().interval_index(i);

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

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
bool move_r<support, sym_t, pos_t, mlf_enc>::query_context_t::prepend(sym_t sym)
{
    if constexpr (mlf_enc != POS) {
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

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
pos_t move_r<support, sym_t, pos_t, mlf_enc>::query_context_t::next_occ()
    requires(supports_multiple_locate)
{
    ensure_be(); // reconstructs b,e and seeds the locate cursor on the first call after a prepend

    if constexpr (has_rlzsa) {
        if (dec.pos() == b) [[unlikely]] {
            // compute the suffix array value at b
            dec.set_value(idx->SA_s(hat_b_ap_y) - (y + 1));
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
            idx->init_phi_m1(b, e, phi.s, phi.s_, hat_b_ap_y, y);
        } else {
            idx->M_Phi_m1().move(phi.s, phi.s_);
        }

        phi.i++;
        return phi.s;
    }
}

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
pos_t move_r<support, sym_t, pos_t, mlf_enc>::query_context_t::one_occ() const
    requires(supports_locate) { return idx->SA_s(hat_b_ap_y) - (y + 1); }

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
template <typename report_fnc_t>
void move_r<support, sym_t, pos_t, mlf_enc>::query_context_t::locate(report_fnc_t report)
    requires(supports_multiple_locate)
{
    ensure_be(); // reconstructs b,e and seeds the locate cursor on the first call after a prepend

    if constexpr (has_rlzsa) {

        if (dec.pos() == b) [[unlikely]] {
            // compute the suffix array value at b
            dec.set_value(idx->SA_s(hat_b_ap_y) - (y + 1));
            report(dec.value());
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

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
bool move_r<support, sym_t, pos_t, mlf_enc>::backward_search_step(
    sym_t sym,
    pos_t& b, pos_t& b_,
    pos_t& e, pos_t& e_,
    pos_t& hat_b_ap_y, int64_t& y,
    pos_t& hat_e_ap_z, int64_t& z) const
{
    // [b,e] is kept as (position, run) endpoints (b,e absolute for positional M_LF, offset within run b_,e_ for
    // differential); no M_LF.p is reconstructed here (in the differential encoding that would cost O(delta)).

    i_sym_t i_sym = map_symbol(sym);

    // If sym does not occur in L', then P[i..m] does not occur in T
    if (i_sym == 0) [[unlikely]] return false;

    // first run with character P[i] at/after run b_ (character successor query)
    if (i_sym != L_(b_)) {
        if constexpr (byte_alphabet) b_ = PS_L_().succ(i_sym, b_, e_);
        else b_ = RS_L_().succ(i_sym, b_, e_);
        if (b_ > e_) [[unlikely]] return false;
        // b at the start of run b_ (absolute p(b_) / offset 0)
        if constexpr (mlf_enc == POS) b = M_LF().p(b_);
        else b = 0;
        hat_b_ap_y = b_;
        y = 0;
    } else {
        y++;
    }

    // last run with character P[i] at/before run e_ (character predecessor query)
    if (i_sym != L_(e_)) {
        if constexpr (byte_alphabet) e_ = PS_L_().pred(i_sym, e_, b_);
        else e_ = RS_L_().pred(i_sym, e_, b_);
        // e at the end of run e_ (absolute p(e_+1)-1 / offset len(e_)-1)
        if constexpr (mlf_enc == POS) e = M_LF().p(e_ + 1) - 1;
        else e = M_LF().len(e_) - 1;
        hat_e_ap_z = e_;
        z = 0;
    } else {
        z++;
    }

    // Else, because each suffix i in the previous suffix array interval starts with P[i+1..m] and the current
    // interval [b,e] contains all suffixes of it, before which there is a P[i] in T, all suffixes in the
    // interval SA[LF(b),LF(e)] start with P[i..m]

    // empty interval (different runs with b_ > e_, or the same run with b > e): P does not occur
    if (b_ > e_ || (b_ == e_ && b > e)) [[unlikely]] return false;

    // set b <- LF(b) and e <- LF(e) via one (position,run)-move query per endpoint; the two branches below
    // just avoid a redundant move when both endpoints lie in the same run
    if (b_ == e_) {
        if (b == e) {
            /* If \hat{b'}_i == \hat{e'}_i and b'_i = e'_i, then computing
            (e_i,\hat{e}_i) <- M_LF.move(e'_i,\hat{e'}_i) is redundant */
            M_LF().move(b, b_);
            e_ = b_; e = b;
        } else {
            // same run but b != e: move b, then advance e by the difference (O(a) by a-balancedness)
            pos_t diff_eb = e - b;
            M_LF().move(b, b_);
            e_ = b_; e = b + diff_eb;
            if constexpr (mlf_enc == POS) {
                while (e >= M_LF().p(e_ + 1)) e_++;
            } else {
                pos_t l;
                while (e >= (l = M_LF().len(e_))) { e -= l; e_++; }
            }
        }
    } else {
        M_LF().move(b, b_);
        M_LF().move(e, e_);
    }

    return true;
}

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
void move_r<support, sym_t, pos_t, mlf_enc>::init_phi_m1(
    pos_t& b_, pos_t& e_,
    pos_t& s, pos_t& s_,
    pos_t& hat_b_ap_y, int64_t& y) const
    requires(has_locate_move)
{
    setup_phi_m1_move_pair(hat_b_ap_y, s, s_);
    s -= y + 1;

    // If there is more than one occurrence and s < M_Phi^{-1}.p[s_], now an input interval of M_Phi^{-1} before
    // the s_-th one contains s, so we have to decrease s_. To find the correct value for s_, we perform
    // an exponential search to the left over the input interval starting positions of M_Phi^{-1} starting at s_.
    if (b_ < e_ && s < M_Phi_m1().p(s_)) {
        s_ = exp_search_max_leq<pos_t, LEFT>(s, 0, s_, [&](pos_t x) { return M_Phi_m1().p(x); });
    }
}

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
pos_t move_r<support, sym_t, pos_t, mlf_enc>::count(const inp_t& P) const
{
    pos_t b_, b, e_, e, hat_b_ap_y, hat_e_ap_z;
    int64_t y, z;
    init_backward_search(b, b_, e, e_, hat_b_ap_y, y, hat_e_ap_z, z);

    for (int64_t i = P.size() - 1; i >= 0; i--) {
        if (!backward_search_step(P[i], b, b_, e, e_, hat_b_ap_y, y, hat_e_ap_z, z)) {
            return 0;
        }
    }

    // interval size e - b + 1 (positional: direct subtraction of absolute components; differential: hybrid distance)
    if constexpr (mlf_enc == POS) return e - b + 1;
    else return M_LF().distance(b_, b, e_, e) + 1;
}

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
template <typename report_fnc_t>
void move_r<support, sym_t, pos_t, mlf_enc>::locate(const inp_t& P, report_fnc_t report) const
    requires(supports_multiple_locate)
{
    pos_t b_, b, e_, e, hat_b_ap_y, hat_e_ap_z;
    int64_t y, z;
    init_backward_search(b, b_, e, e_, hat_b_ap_y, y, hat_e_ap_z, z);

    for (int64_t i = P.size() - 1; i >= 0; i--) {
        if (!backward_search_step(P[i], b, b_, e, e_, hat_b_ap_y, y, hat_e_ap_z, z)) {
            return;
        }
    }

    // reconstruct the absolute interval [b,e] once (for differential M_LF; positional already has it), then locate
    if constexpr (mlf_enc != POS) { b = M_LF().pos(b_, b); e = M_LF().pos(e_, e); }

    if constexpr (has_rlzsa) {
        pos_t s = SA_s(hat_b_ap_y) - (y + 1);
        report(s);

        if (b < e) {
            auto dec = _rlzsa.decode();
            dec.init_right(b + 1);
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

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
template <typename report_fnc_t>
void move_r<support, sym_t, pos_t, mlf_enc>::revert_range(report_fnc_t report, retrieve_params params) const
{
    adjust_retrieve_params(params, n - 2);

    pos_t l = params.l;
    pos_t r = params.r;

    // the parallelism is capped at the number of sections overlapping [l, r]: each thread starts its backward walk
    // at a section boundary, so more threads than that would not add useful parallelism
    uint16_t s_l, s_r;

    if (p_r == 1) {
        s_l = 0;
        s_r = 0;
    } else {
        s_l = bin_search_min_gt<pos_t>(l, 0, p_r - 1, [&](pos_t x) { return _D_e[x].second; });
        s_r = bin_search_min_geq<pos_t>(r, 0, p_r - 1, [&](pos_t x) { return _D_e[x].second; });
    }

    uint16_t p = retrieve_num_threads(params.num_threads, s_r - s_l + 1);

    run_parallel_threads(p, [&](uint16_t i_p) {
        // the buffer-aligned position range [b, e] this thread reports (may be empty after alignment)
        pos_t b = thread_range_start(i_p, p, l, r, params.buffer_align);
        pos_t e = thread_range_start(i_p + 1, p, l, r, params.buffer_align);
        if (b >= e) return; // empty chunk (a buffer block wider than the chunk claimed it for another thread)
        e--; // make the right end inclusive

        // find the section boundary >= e to start the backward walk from (with a known BWT position)
        uint16_t sect = p_r == 1 ? 0 : bin_search_min_geq<pos_t>(e, 0, p_r - 1, [&](pos_t y) { return _D_e[y].second; });
        // walk-start (>= e) and its LF cursor (run x sits at a run start, i.e. offset 0 / absolute p(x))
        pos_t j = sect == p_r - 1 ? n - 2 : _D_e[sect].second;
        pos_t x = sect == p_r - 1 ? pos_t(0) : _D_e[sect].first;
        pos_t d; // position component of the LF cursor: absolute p(x) (positional) or offset 0 (differential)
        if constexpr (mlf_enc == POS) d = M_LF().p(x);
        else d = 0;

        // skip (without reporting) from the section boundary down to e (one LF step per iteration, move(component, run))
        while (j > e) {
            M_LF().move(d, x);
            j--;
        }

        // report T[e] = L'[x], then T[e-1,b] from right to left
        report(j, unmap_symbol(L_(x)));

        while (j > b) {
            M_LF().move(d, x);
            j--;
            report(j, unmap_symbol(L_(x)));
        }
    });
}

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
template <typename report_fnc_t>
void move_r<support, sym_t, pos_t, mlf_enc>::BWT_range(report_fnc_t report, retrieve_params params) const
{
    adjust_retrieve_params(params, n - 1);

    pos_t l = params.l;
    pos_t r = params.r;

    uint16_t p = retrieve_num_threads(params.num_threads, (r - l + 1) / 10);

    run_parallel_threads(p, [&](uint16_t i_p) {
        // the buffer-aligned position range [b, e] of thread i_p (may be empty after alignment)
        pos_t b = thread_range_start(i_p, p, l, r, params.buffer_align);
        pos_t e = thread_range_start(i_p + 1, p, l, r, params.buffer_align);
        if (b >= e) return; // empty chunk (a buffer block wider than the chunk claimed it for another thread)
        e--; // make the right end inclusive

        // Current position in the bwt.
        pos_t i = b;

        // the input interval of M_LF containing i, as (run,offset); the next run boundary p(x+1) is then
        // advanced incrementally via the interval lengths (no M_LF.p reconstruction per run)
        auto [x, d] = M_LF().run_and_offset(i);

        // start position of the next input interval in M_LF, i.e. p(x+1) = i + (len(x) - d)
        pos_t l_xp1 = i + (M_LF().len(x) - d);

        // iterate until x is the input interval containing e
        while (l_xp1 <= e) {

            // iterate over all positions in the x-th input interval
            while (i < l_xp1) {
                report(i, unmap_symbol(L_(x)));
                i++;
            }

            x++;
            l_xp1 += M_LF().len(x);
        }

        // report the remaining characters
        while (i <= e) {
            report(i, unmap_symbol(L_(x)));
            i++;
        }
    });
}

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
template <typename report_fnc_t>
void move_r<support, sym_t, pos_t, mlf_enc>::SA_range(report_fnc_t report, retrieve_params params) const
    requires(supports_multiple_locate)
{
    adjust_retrieve_params(params, n - 1);

    pos_t l = params.l;
    pos_t r = params.r;

    // use at most (r-l+1)*(r''/n)*(1/10) threads
    uint16_t p = retrieve_num_threads(params.num_threads,
        (uint64_t)(((r - l + 1) * (double) r__) / (10.0 * (double) n)));

    run_parallel_threads(p, [&](uint16_t i_p) {
        // the buffer-aligned position range [b, e] of thread i_p (may be empty after alignment)
        pos_t b = thread_range_start(i_p, p, l, r, params.buffer_align);
        pos_t e = thread_range_start(i_p + 1, p, l, r, params.buffer_align);
        if (b >= e) return; // empty chunk (a buffer block wider than the chunk claimed it for another thread)
        e--; // make the right end inclusive

        if constexpr (has_locate_move) {
            // the input interval of M_LF containing i
            pos_t x = M_LF().interval_index(b);

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
            pos_t x = M_LF().interval_index(b);

            // position in the suffix array of the current suffix s
            pos_t j = M_LF().p(x);

            // the current suffix (s = SA[j])
            pos_t s = SA_s(x);
            if (j == b) report(b, s);

            if (j < e) [[likely]] {
                j++;

                // initialize the rlzsa decode context to position j (with s = SA[j-1])
                auto dec = _rlzsa.decode();
                dec.init_right(j);
                dec.set_value(s);

                // advance the rlzsa context to position b
                dec.skip_right(b);

                // decode and report SA[b..e]
                dec.report_right(e, report);
            }
        }
    });
}

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
template <typename output_t, bool output_reversed>
void move_r<support, sym_t, pos_t, mlf_enc>::retrieve_range(
    void (move_r<support, sym_t, pos_t, mlf_enc>::*retrieve_method)(
        std::function<void(pos_t, output_t)>&&, move_r<support, sym_t, pos_t, mlf_enc>::retrieve_params) const,
    std::string file_name, move_r<support, sym_t, pos_t, mlf_enc>::retrieve_params params) const
{
    pos_t l = params.l;
    pos_t r = params.r;
    uint16_t num_threads = params.num_threads;
    uint64_t buffer_bytes = std::max<uint64_t>(1024,
        params.max_bytes_alloc != -1 ? params.max_bytes_alloc / num_threads : ((r - l + 1) * sizeof(output_t)) / (num_threads * 500));

    uint64_t align = std::max<uint64_t>(1, buffer_bytes / sizeof(output_t));
    params.buffer_align = align;
    uint64_t buffer_size_per_thread = align * sizeof(output_t);

    std::filesystem::resize_file(std::filesystem::current_path() / file_name, (r - l + 1) * sizeof(output_t));
    std::vector<sdsl::int_vector_buffer<sizeof(output_t) * 8>> file_bufs;

    for (uint16_t i = 0; i < num_threads; i++) {
        file_bufs.emplace_back(sdsl::int_vector_buffer<sizeof(output_t) * 8>(file_name, std::ios::in, buffer_size_per_thread, sizeof(output_t) * 8, true));
    }

    (this->*retrieve_method)([&](pos_t pos, output_t val) {
        file_bufs[omp_get_thread_num()][pos] = *reinterpret_cast<uint64_t*>(&val);
    }, params);
}