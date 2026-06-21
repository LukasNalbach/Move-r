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
void move_r<support, sym_t, pos_t>::construction::build_ilf()
{
    log_phase_start("building I_LF");

    no_init_resize(I_LF, r);

    #pragma omp parallel num_threads(p_)
    {
        // Index in [0..p'-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Iteration range start position of thread i_p.
        pos_t b_r = r_p[i_p];

        // Number of BWT runs in thread i_p's section.
        pos_t rp_diff = r_p[i_p + 1] - r_p[i_p];

        // i', Start position of the last-seen run.
        pos_t i_ = n_p[i_p];

        // Build I_LF
        for (pos_t i = 0; i < rp_diff; i++) {
            /* Write the pair (i',LF(i')) to the next position i in I_LF, where
            LF(i') = C[L[i']] + rank(L,L[i'],i'-1) = C[p'][L[i']] + C[i_p][L[i']]. */
            I_LF[b_r + i] = std::make_pair(i_, C[p_][run_sym(i_p, i)] + C[i_p][run_sym(i_p, i)]);

            /* Update the rank-function in C[i_p] to store C[i_p][c] = rank(L,c,i'-1),
            for each c in [0..255] */
            C[i_p][run_sym(i_p, i)] += run_len(i_p, i);

            // Update the position of the last-seen run.
            i_ += run_len(i_p, i);
        }
    }

    C.clear();
    C.shrink_to_fit();

    log_phase_end("time_build_ilf");
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::build_mlf()
{
    if (log && mf_mds != nullptr) {
        *mf_mds << "RESULT"
                << " algo=build_mlf"
                << " text=" << name_text_file
                << " num_threads=" << p
                << " a=" << idx.a;
    }
    log_phase_start("\nbuilding M_LF");

    uint8_t width_l_ = std::bit_width(idx.sigma);

    idx._M_LF = move_data_structure_l_<pos_t, i_sym_t>(
        std::move(I_LF), n, {
            .num_threads = p,
            .a = idx.a,
            .log = log,
            .mf = mf_mds,
        },
        width_l_);

    r_ = idx._M_LF.num_intervals();
    idx.r_ = r_;

    log_mds_phase_end("time_build_mlf", "r_", r_);
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::build_l__sas()
{
    bool build_sas = p == 1 && (support == _locate_move || support == _locate_rlzsa);

    log_phase_start(build_sas ? "building L' and SA_s" : "building L'");

    if (build_sas) {
        if constexpr (support == _locate_move) {
            no_init_resize(SA_s, r_);
        } else {
            idx._SA_s = interleaved_bit_aligned_vectors<pos_t>({ std::bit_width(uint64_t(n)) });
            idx._SA_s.resize_no_init(r_);
        }
    }

    // Simultaneously iterate over the input intervals of M_LF nad the bwt runs to build L'
    #pragma omp parallel num_threads(p_)
    {
        // Index in [0..p'-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Bwt range start position of thread i_p.
        pos_t b = n_p[i_p];

        // Iteration range start position of thread i_p.
        pos_t b_r = r_p[i_p];

        // Number of runs in thread i_p's section.
        pos_t rp_diff = r_p[i_p + 1] - r_p[i_p];

        // Index of the current input interval in M_LF, initially the index of the input interval of M_LF containing b
        pos_t j = bin_search_max_leq<pos_t>(b, 0, r_ - 1, [&](pos_t x) { return idx._M_LF.p(x); });
        // Starting position of the next bwt run.
        pos_t l_ = b;

        for (pos_t i = 0; i < rp_diff; i++) {
            idx._M_LF.template set_L_(j, run_sym(i_p, i));

            if (build_sas) {
                if constexpr (support == _locate_move) {
                    SA_s[j] = I_Phi_m1[b_r + i].second;
                } else {
                    idx._SA_s.template set_parallel<0, pos_t>(j, I_Phi_m1[b_r + i].second);
                }
            }

            j++;

            // update l_ to the next run start position
            l_ += run_len(i_p, i);

            // iterate over all input intervals in M_LF within the i-th bwt run in thread i_p's section that have been
            // created by the balancing algorithm
            while (idx._M_LF.p(j) < l_) {
                if (build_sas) {
                    if constexpr (support == _locate_move) {
                        SA_s[j] = n;
                    } else {
                        idx._SA_s.template set_parallel<0, pos_t>(j, n);
                    }
                }

                idx._M_LF.template set_L_(j, run_sym(i_p, i));
                j++;
            }
        }
    }

    n_p.clear();
    n_p.shrink_to_fit();

    r_p.clear();
    r_p.shrink_to_fit();

    RLBWT.clear();
    RLBWT.shrink_to_fit();

    log_phase_end("time_build_l__sas");
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::build_rsl_()
{
    log_phase_start("building RS_L'");

    idx._RS_L_ = rsl_t([&](pos_t i) { return idx.L_(i); }, idx.sigma, 0, r_ - 1);

    log_phase_end("time_build_rsl_");
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::build_l_prev_next()
{
    log_phase_start("building L'_prev & L'_next");

    uint64_t bits_prev_next = std::bit_width(uint64_t(r_)); // tight width; values are in [0,r_] (r_ = no occurrence)
    int64_t blk_size = std::min<pos_t>(_l_blk_size_factor * idx.sigma, r_);
    int64_t num_blks = div_ceil<int64_t>(r_, blk_size);
    idx._l_blk_size = blk_size;
    idx._num_blks_l_ = num_blks;

    #pragma omp parallel sections num_threads(p)
    {
        #pragma omp section
        {
            idx._L_next = interleaved_bit_aligned_vectors<pos_t>({ bits_prev_next });
            idx._L_next.resize_no_init((num_blks + 1) * idx.sigma);
            std::vector<pos_t> next(idx.sigma, r_);

            for (int64_t blk = num_blks; blk >= 0; blk--) {
                int64_t beg = blk * blk_size;
                int64_t end = std::min<int64_t>(beg + blk_size, r_);
                pos_t beg_next = blk * idx.sigma;

                for (int64_t i = end - 1; i >= beg; i--) {
                    next[idx.L_(i)] = i;
                }

                for (pos_t c = 0; c < idx.sigma; c++) {
                    idx._L_next.template set<0, pos_t>(beg_next + c, next[c]);
                }
            }
        }

        #pragma omp section
        {
            idx._L_prev = interleaved_bit_aligned_vectors<pos_t>({ bits_prev_next });
            idx._L_prev.resize_no_init(num_blks * idx.sigma);
            std::vector<pos_t> prev(idx.sigma, r_);

            for (int64_t blk = 0; blk < num_blks; blk++) {
                int64_t beg = blk * blk_size;
                int64_t end = std::min<int64_t>(beg + blk_size, r_);
                pos_t beg_prev = blk * idx.sigma;

                for (pos_t c = 0; c < idx.sigma; c++) {
                    idx._L_prev.template set<0, pos_t>(beg_prev + c, prev[c]);
                }

                for (int64_t i = beg; i < end; i++) {
                    prev[idx.L_(i)] = i;
                }
            }
        }
    }

    log_phase_end("time_build_l_prev_next");
}
