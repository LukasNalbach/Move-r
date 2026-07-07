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
#include <move_r/move_r.hpp>

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
template <rlbwt_build_mode mode, typename sa_sint_t>
void move_r<support, sym_t, pos_t, mlf_enc>::construction::build_rlbwt_c()
{
    log_phase_start("building RLBWT");

    std::vector<sa_sint_t>& SA = get_sa<sa_sint_t>(); // [0..n-1] The suffix array

    p_ = byte_alphabet ? p : 1;
    r_p.resize(p_ + 1, 0);
    uint8_t width_bwt = byte_alphabet ? 1 : byte_width(idx.sigma);
    RLBWT.resize(p_, interleaved_byte_aligned_vectors<uint32_t, uint32_t>({ width_bwt, 4 }));

    if constexpr (byte_alphabet) {
        C.resize(p_, std::vector<pos_t>(byte_alphabet ? 256 : idx.sigma, 0));
    }

    if constexpr (mode == _bwt_file) {
        for (uint16_t i = 0; i < p; i++) {
            BWT_file_bufs.emplace_back(sdsl::int_vector_buffer<8>(
                prefix_tmp_files + ".bwt", std::ios::in,
                128 * 1024, 8, true));
        }
    }

    for (uint16_t i = 0; i < p_; i++) {
        n_p.emplace_back(i * (n / p_));
    }

    n_p.emplace_back(n);

    #pragma omp parallel num_threads(p_)
    {
        // Index in [0..p_-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Iteration range start position of thread i_p.
        pos_t b = n_p[i_p];
        // Iteration range end position of thread i_p.
        pos_t e = n_p[i_p + 1];

        // Store in C[i_p][c] the number of occurrences of c in L[b..e], for each c in [0..255].

        i_sym_t prev_sym; // previous symbol in L.
        i_sym_t cur_sym; // current symbol in L.
        pos_t i_ = b; // start position of the last seen run in L.

        if constexpr (mode == _sa) {
            prev_sym = SA[b] == 0 ? 0 : T<i_sym_t>(SA[b] - 1);
        } else if constexpr (mode == _bwt) {
            prev_sym = char_to_uchar(L[b]);
        } else if constexpr (mode == _bwt_file) {
            prev_sym = BWT_file_bufs[i_p][b];

            if (prev_sym < 3) [[unlikely]] {
                prev_sym = 0;
            } else {
                prev_sym -= 2;
            }
        }

        // Iterate over the range L[b+1..e-1]
        for (pos_t i = b + 1; i < e; i++) {
            if constexpr (mode == _sa) {
                cur_sym = SA[i] == 0 ? 0 : T<i_sym_t>(SA[i] - 1);
            } else if constexpr (mode == _bwt) {
                cur_sym = char_to_uchar(L[i]);
            } else if constexpr (mode == _bwt_file) {
                cur_sym = BWT_file_bufs[i_p][i];

                if (cur_sym < 3) [[unlikely]] {
                    cur_sym = 0;
                } else {
                    cur_sym -= 2;
                }
            }

            // check if there is a run starting at L[i]
            if (cur_sym != prev_sym) {
                add_run(i_p, prev_sym, i - i_);
                if constexpr (byte_alphabet)
                    C[i_p][prev_sym] += i - i_;
                prev_sym = cur_sym;
                i_ = i;
            }
        }

        // add the run L[i'..e)
        add_run(i_p, prev_sym, e - i_);
        if constexpr (byte_alphabet)
            C[i_p][prev_sym] += e - i_;
        // Store in r_p[i_p] the number of runs starting in L[b..e).
        r_p[i_p] = RLBWT[i_p].size();
        RLBWT[i_p].shrink_to_fit();
    }

    // for i_p \in [1,p'-2], merge the last run in thread i_p's section with the first run in thread
    // i_p+1's section, if their characters are equal
    for (uint16_t i_p = 0; i_p < p_ - 1; i_p++) {
        i_sym_t c = run_sym(i_p, r_p[i_p] - 1);

        if (run_sym(i_p + 1, 0) == c) {
            pos_t l = run_len(i_p, r_p[i_p] - 1);
            set_run_len(i_p + 1, 0, run_len(i_p + 1, 0) + l);

            if constexpr (byte_alphabet) {
                C[i_p][c] -= l;
                C[i_p + 1][c] += l;
            }

            n_p[i_p + 1] -= l;
            r_p[i_p]--;
            RLBWT[i_p].resize(r_p[i_p]);
        }
    }

    if constexpr (mode == _bwt_file) {
        for (uint16_t i = 0; i < p; i++) {
            BWT_file_bufs[i].close();
        }

        BWT_file_bufs.clear();
        BWT_file_bufs.shrink_to_fit();
        std::filesystem::remove(prefix_tmp_files + ".bwt");
    }

    if (&L == &L_tmp) {
        L.clear();
        L.shrink_to_fit();
    }

    if (delete_T) {
        T_str.clear();
        T_str.shrink_to_fit();
        T_vec.clear();
        T_vec.shrink_to_fit();
    }

    if ((support == _count || support == _count_bi) && mode != _bwt) {
        SA.clear();
        SA.shrink_to_fit();
    }

    /* Now, r_p[i_p] stores the number of runs starting in the iteration range L[b..e] of thread
    i_p in [0..p'-1], and r_p[p'] = 0. We want r_p[i_p] to store the number of runs starting before
    the iteration range start position b of thread i_p in [0..p'-1]. Also, we want r_p[p'] to store
    the number of all runs. This way, we can build I_LF[r_p[i_p]..r_p[i_p+1]-1] with thread i_p in [0..p'-1]
    using the C-array in C[p'] while using and updating the rank-function in C[i_p] on-the-fly. */

    for (uint16_t i = p_; i > 0; i--) {
        r_p[i] = r_p[i - 1];
    }

    r_p[0] = 0;

    for (uint16_t i = 2; i <= p_; i++) {
        r_p[i] += r_p[i - 1];
    }

    /* Now, r_p[i_p] stores the number of runs starting before the iteration range start position b of
    thread i_p in [0..p'-1] and r_p[p'] stores the number of all runs, so we are done with r_p[0..p'] */

    r = r_p[p_];
    idx.r = r;

    if constexpr (int_alphabet) {
        C.emplace_back(std::vector<pos_t>(idx.sigma, 0));

        for (pos_t i = 0; i < r; i++) {
            C[0][run_sym(0, i)] += run_len(0, i);
        }
    }

    process_c();

    log_phase_end("time_build_rlbwt");
}

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
void move_r<support, sym_t, pos_t, mlf_enc>::construction::process_c()
{
    /* Now, C[i_p][c] is the number of occurrences of c in L[b..e], where [b..e] is the range of the
    thread i_p in [0..p'-1]. Also, we have C[p'][0..255] = 0. */

    /* We want to have C[i_p][c] = rank(L,c,b-1), where b is the iteration range start position of
    thread i_p in [0..p'-1]. Also, we want C[p'][0..255] to be the C-array, that is C[p'][c] stores
    the number of occurrences of all smaller characters c' < c in L[0..n-1], for c in [0..255]. */

    pos_t max_symbol = byte_alphabet ? 256 : idx.sigma;

    for (pos_t i = 1; i < p_; i++) {
        #pragma omp parallel for num_threads(p)
        for (uint64_t j = 0; j < max_symbol; j++) {
            C[i][j] += C[i - 1][j];
        }
    }

    /* Now, we have C[i_p][c] = rank(L,c,e), for each c in [0..255] and i_p in [0..p'-1],
    where e is the iteration range end position of thread i_p. */

    C.insert(C.begin(), std::vector<pos_t>());
    no_init_resize(C[0], max_symbol);

    #pragma omp parallel for num_threads(p)
    for (uint64_t j = 0; j < max_symbol; j++) {
        C[0][j] = 0;
    }

    /* Now, we have C[i_p][c] = rank(L,c,b-1), for each c in [0..255] and i_p in [0..p'-1],
    where b is the iteration range start position of thread i_p, so we are done with C[0..p'-1][0..255].
    Also, we have C[p'][c] = rank(L,c,n-1), for c in [0..255]. */

    pos_t cp_im1 = C[p_][0]; // stores at the start of the i-th iteration C[p'][i-1]
    C[p_][0] = 0;
    pos_t cp_im1_tmp;

    for (pos_t i = 1; i < max_symbol; i++) {
        cp_im1_tmp = C[p_][i];
        C[p_][i] = cp_im1 + C[p_][i - 1];
        cp_im1 = cp_im1_tmp;
    }

    // Now we are done with C, since C[p'] is the C-array.
}
