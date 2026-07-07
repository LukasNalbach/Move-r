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
template <bool bigbwt, typename sa_sint_t>
void move_r<support, sym_t, pos_t, mlf_enc>::construction::build_sas_from_sa()
{
    log_phase_start("building SA_s and SA_e");

    if constexpr (bigbwt) {
        for (uint16_t i = 0; i < p; i++) {
            SA_file_bufs.emplace_back(sdsl::int_vector_buffer<40>(
                prefix_tmp_files + ".sa", std::ios::in,
                128 * 1024, 40, true));
        }
    }

    idx._SA_s = interleaved_bit_aligned_vectors<pos_t>({ bit_width(n) });
    idx._SA_e = interleaved_bit_aligned_vectors<pos_t>({ bit_width(n) });

    idx._SA_s.resize_no_init(r_);
    idx._SA_e.resize_no_init(r_);

    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < r_; i++) {
        idx._SA_s.template set_parallel<0, pos_t>(i,
            SA<bigbwt, sa_sint_t>(omp_get_thread_num(), idx._M_LF.p(i)));

        idx._SA_e.template set_parallel<0, pos_t>(i,
            SA<bigbwt, sa_sint_t>(omp_get_thread_num(), idx._M_LF.p(i + 1) - 1));
    }

    log_phase_end("time_build_sa_s_sa_s_");
}

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
template <bool bigbwt, typename sa_sint_t>
void move_r<support, sym_t, pos_t, mlf_enc>::construction::build_iphim1_sas_from_sa()
{
    log_phase_start(is_bidirectional ? "building I_Phi^{-1} and SA_s and SA_e" : "building I_Phi^{-1} and SA_s");

    if constexpr (bigbwt) {
        for (uint16_t i = 0; i < p; i++) {
            SA_file_bufs.emplace_back(sdsl::int_vector_buffer<40>(
                prefix_tmp_files + ".sa", std::ios::in,
                128 * 1024, 40, true));
        }
    }

    no_init_resize(I_Phi_m1, r);

    if constexpr (support == _locate_move) {
        no_init_resize(SA_s, r_);
    } else {
        idx._SA_s = interleaved_bit_aligned_vectors<pos_t>({ bit_width(n) });
        idx._SA_s.resize_no_init(r_);

        if constexpr (is_bidirectional) {
            idx._SA_e = interleaved_bit_aligned_vectors<pos_t>({ bit_width(n) });
            idx._SA_e.resize_no_init(r_);
        }
    }

    I_Phi_m1[0] = std::make_pair(
        SA<bigbwt, sa_sint_t>(0, n - 1),
        SA<bigbwt, sa_sint_t>(0, 0));

    #pragma omp parallel num_threads(p_)
    {
        uint16_t i_p = omp_get_thread_num();

        // Iteration range start position of thread i_p.
        pos_t b_r = r_p[i_p];

        // Number of BWT runs within thread i_p's section.
        pos_t rp_diff = r_p[i_p + 1] - r_p[i_p];

        // start position of the current BWT run.
        pos_t j = n_p[i_p];

        // Index of the current input interval in M_LF, initially the index of the input interval of M_LF containing b
        pos_t x = bin_search_max_leq<pos_t>(j, 0, r_ - 1, [&](pos_t y) { return idx._M_LF.p(y); });

        for (pos_t i = 0; i < rp_diff; i++) {
            if (j != 0) [[likely]] {
                I_Phi_m1[b_r + i] = std::make_pair(
                    SA<bigbwt, sa_sint_t>(i_p, j - 1),
                    SA<bigbwt, sa_sint_t>(i_p, j));
            }

            pos_t k = j;
            j += run_len(i_p, i);
            bool first_interval = true;

            while (k < j) {
                if constexpr (support == _locate_move) {
                    if (first_interval) {
                        SA_s[x] = SA<bigbwt, sa_sint_t>(i_p, k);
                        first_interval = false;
                    } else {
                        SA_s[x] = n;
                    }
                } else {
                    idx._SA_s.template set_parallel<0, pos_t>(x,
                        SA<bigbwt, sa_sint_t>(i_p, k));
                }

                x++;
                k = idx._M_LF.p(x);

                if constexpr (is_bidirectional) {
                    idx._SA_e.template set_parallel<0, pos_t>(x - 1,
                        SA<bigbwt, sa_sint_t>(i_p, k - 1));
                }
            }
        }
    }

    if (!build_sa_and_l && !has_rlzsa && support != _locate_move_bi_fwd && support != _locate_rlzsa_bi_fwd) {
        std::vector<sa_sint_t>& SA = get_sa<sa_sint_t>(); // [0..n-1] The suffix array

        SA.clear();
        SA.shrink_to_fit();
    }

    log_phase_end("time_build_iphi");
}

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
void move_r<support, sym_t, pos_t, mlf_enc>::construction::build_iphi()
{
    log_phase_start("\nbuilding I_Phi");

    no_init_resize(I_Phi, r);

    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < r; i++) {
        I_Phi[i] = std::make_pair(I_Phi_m1[i].second, I_Phi_m1[i].first);
    }

    log_phase_end();
}

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
void move_r<support, sym_t, pos_t, mlf_enc>::construction::sort_iphim1()
{
    log_phase_start("sorting I_Phi^{-1}");

    // Sort I_Phi^{-1} by the starting positions of its input intervals.
    auto comp_I_Phi = [](std::pair<pos_t, pos_t> p1, std::pair<pos_t, pos_t> p2) { return p1.first < p2.first; };

    // Choose the correct sorting algorithm.
    if (p > 1) {
        ips4o::parallel::sort(I_Phi_m1.begin(), I_Phi_m1.end(), comp_I_Phi);
    } else {
        ips4o::sort(I_Phi_m1.begin(), I_Phi_m1.end(), comp_I_Phi);
    }

    if (log && mf_mds != nullptr) {
        *mf_mds << "RESULT"
                << " algo=build_mphi"
                << " text=" << name_text_file
                << " num_threads=" << p
                << " a=" << idx.a;
    }

    log_phase_end("time_sort_iphim1");
}

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
void move_r<support, sym_t, pos_t, mlf_enc>::construction::sort_iphi()
{
    log_phase_start("sorting I_Phi");

    // Sort I_Phi by the starting positions of its input intervals.
    auto comp_I_Phi = [](std::pair<pos_t, pos_t> p1, std::pair<pos_t, pos_t> p2) { return p1.first < p2.first; };

    // Choose the correct sorting algorithm.
    if (p > 1) {
        ips4o::parallel::sort(I_Phi.begin(), I_Phi.end(), comp_I_Phi);
    } else {
        ips4o::sort(I_Phi.begin(), I_Phi.end(), comp_I_Phi);
    }

    if (log && mf_mds != nullptr) {
        *mf_mds << "RESULT"
                << " algo=build_mphi"
                << " text=" << name_text_file
                << " num_threads=" << p
                << " a=" << idx.a;
    }

    log_phase_end("time_sort_iphi");
}

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
void move_r<support, sym_t, pos_t, mlf_enc>::construction::build_mphim1()
{
    log_phase_start("\nbuilding M_Phi^{-1}");

    idx._M_Phi_m1 = move_data_structure<pos_t, POS>(
        std::move(I_Phi_m1), n, {
            .num_threads = p,
            .a = idx.a,
            .log = log,
            .mf = mf_mds,
        },
        {}, &pi_mphi);

    r__ = idx._M_Phi_m1.num_intervals();
    idx.r__ = r__;

    log_mds_phase_end("time_build_mphim1", "r__", r__);
}

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
void move_r<support, sym_t, pos_t, mlf_enc>::construction::build_mphi()
{
    log_phase_start("\nbuilding M_Phi");

    idx._M_Phi = move_data_structure<pos_t, POS>(
        std::move(I_Phi), n, {
            .num_threads = p,
            .a = idx.a,
            .log = log,
            .mf = mf_mds,
        },
        {}, &pi_mphi);

    r___ = idx._M_Phi.num_intervals();
    idx.r___ = r___;

    log_mds_phase_end("time_build_mphi", "r___", r___);
}

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
void move_r<support, sym_t, pos_t, mlf_enc>::construction::build_saphim1()
{
    log_phase_start("building SA_Phi^{-1}");

    no_init_resize(pi_, r_);

    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < r_; i++) {
        pi_[i] = i;
    }

    auto comp_pi_ = [&](pos_t i, pos_t j) { return SA_s[i] < SA_s[j]; };
    if (p > 1) {
        ips4o::parallel::sort(pi_.begin(), pi_.end(), comp_pi_);
    } else {
        ips4o::sort(pi_.begin(), pi_.end(), comp_pi_);
    }

    idx.omega_idx = idx._M_Phi_m1.width_idx();
    idx._SA_Phi_m1 = interleaved_bit_aligned_vectors<pos_t>({ idx.omega_idx });
    idx._SA_Phi_m1.resize_no_init(r_);

    /* Now we will divide the range [0..n-1] up into p non-overlapping sub-ranges [s[i_p]..s[i_p+1]-1],
    for each i_p in [0..p-1], with 0 = s[0] < s[1] < ... < s[p] = n, where
    s[i_p] = min {s' in [0,n-1], s.t. x[i_p] + u[i_p] - 2 >= i_p * lfloor (r+r'')/p rfloor, where
                    x[i_p] = min {x' in [0,r''-1], s.t. M_Phi^{-1}.q(x') >= s'} and
                    u[i_p] = min {u' in [0,r-1], s.t. SA_s[u'] >= s'}
    }.
    By doing so, we ensure that the number of the output intervals of M_Phi^{-1} starting in the range
    [s[i_p]..s[i_p+1]-1] plus the number of suffix array samples in SA_s lying in the range
    [s[i_p]..s[i_p+1]-1] is lfloor (r+r'')/p rfloor +- 1. This property is useful, because it
    ensures that if with each thread i_p, we simultaneously iterate over those, then each thread
    iterates over almost exactly the same number lfloor (r+r'')/p rfloor +- 1 of entries in M_Phi^{-1}
    and SA_s combined. This way, we can acheive good load-balancing. Because we do not have to access
    s[0..p] later, we will not store those values in an array. */

    /* [0..p], x[i_p] = min {x' in [0,r''-1], s.t. M_Phi^{-1}.q(x') >= s'} stores the number of output
    intervals in M_Phi^{-1} starting before s[i_p]. */
    std::vector<pos_t> x(p + 1);
    x[0] = 0;
    x[p] = r__;

    /* [0..p], u[i_p] = min {u' in [0,r-1], s.t. SA_s[u'] >= s'} stores the number of suffix array
    samples in SA_s that are smaller than s[i_p]. */
    std::vector<pos_t> u(p + 1);
    u[0] = 0;
    u[p] = r;

    // Compute s[1..p-1], x[1..p-1] and u[1..p-1].
    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // The optimal value i_p * lfloor (r+r'')/p rfloor for s[i_p].
        pos_t o = i_p * ((r + r__) / p);

        // Left interval limit of the binary search for s[i_p].
        pos_t l_s;
        // Number of output intervals in M_Phi^{-1} starting before s[i_p].
        pos_t l_x;
        // Number of suffix array samples in SA_s that are smaller than s[i_p].
        pos_t l_u;
        // Candidate position in the binary search for s[i_p].
        pos_t m_s;
        // Right interval limit of the binary search for s[i_p].
        pos_t r_s;

        // Initialize the search range for s[i_p] to [0..n-1].
        l_s = 0;
        r_s = n - 1;

        // Perform a binary search over [0..n-1], to find s[i_p].
        while (true) {
            /* Set the Candidate position for s[i_p] to the position in the middle
            between l_s and r_s. */
            m_s = l_s + (r_s - l_s) / 2;

            // Find the minimum x' in [0,r''-1], s.t. M_Phi^{-1}.q(pi_mphi[x']) >= m_s.
            l_x = bin_search_min_geq<pos_t>(m_s, 0, r__ - 1,
                [&](pos_t y) { return idx._M_Phi_m1.q(pi_mphi[y]); });

            // Find the minimum u' in [0,r-1], s.t. SA_s[pi'[u']] >= m_s.
            l_u = bin_search_min_geq<pos_t>(m_s, 0, r - 1,
                [&](pos_t y) { return SA_s[pi_[y]]; });

            /* If l_s = r_s, then l_s is an optimal value for s[i_p] and l_x and
            l_u are valid values for x[i_p] and u[i_p], respectively. */
            if (l_s == r_s) {
                break;
            }

            // Else, adjust the range for the binary search over [0..n-1].
            if (l_x + l_u < o) {
                l_s = m_s + 1;
            } else {
                r_s = m_s;
            }
        }

        // Store l_x and l_u in x[i_p] and u[i_p], respectively.
        x[i_p] = l_x;
        u[i_p] = l_u;

        #pragma omp barrier

        // Iteration range start position in the output intervals of M_Phi^{-1}.
        pos_t i = x[i_p];
        // Iteration range start position in SA_s.
        pos_t j = u[i_p];
        // Iteration range end position + 1 in SA_s.
        pos_t j_ = u[i_p + 1];

        // Check if the range, over which the thread i_p has to iterate in SA_s, is empty
        if (j < j_) {
            while (idx._M_Phi_m1.q(pi_mphi[i]) != SA_s[pi_[j]]) {
                i++;
            }

            /* Iterate over SA_s[pi'[j]],SA_s[pi'[j+1]],...,SA_s[pi'[j'-1]] */
            while (j < j_) {
                // Skip the output intervals the balancing algorithm has added to I_Phi^{-1}
                while (idx._M_Phi_m1.q(pi_mphi[i]) != SA_s[pi_[j]]) {
                    i++;
                }

                idx.set_SA_Phi_m1(pi_[j], pi_mphi[i]);

                i++;
                j++;
            }
        }
    }

    /* Since we set SA_s[j] = n for each j-th input interval of M_LF, whiches starting position is not the starting position of a bwt run,
    * where j \in [0,r'), SA_s[pi'[r]],SA_s[pi'[r+1]],...,SA_s[pi'[r'-1]] = n holds, hence we set SA_Phi^{-1}[pi'[i]] = r'' for
    * i \in [r,r') to mark that we cannot recover SA_s[M_LF[x]] = M_Phi^{-1}.q[SA_Phi^{-1}[x]] for each x-th input interval of
    * M_LF, whiches starting position is not the starting position of a bwt run and where x \in [0,r'). */
    #pragma omp parallel for num_threads(p)
    for (uint64_t i = r; i < r_; i++) {
        idx.set_SA_Phi_m1(pi_[i], r__);
    }

    x.clear();
    x.shrink_to_fit();

    u.clear();
    u.shrink_to_fit();

    pi_mphi.clear();
    pi_mphi.shrink_to_fit();

    log_phase_end("time_build_saphim1");
}

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
void move_r<support, sym_t, pos_t, mlf_enc>::construction::build_de()
{
    if constexpr (supports_locate && !is_bidirectional) {
        idx.p_r = std::min<pos_t>(256, std::max<pos_t>(1, r / 100));
    }

    idx._D_e.resize(idx.p_r - 1);

    #pragma omp parallel for num_threads(p)
    for (uint16_t i = 0; i < idx.p_r - 1; i++) {
        pos_t x = bin_search_min_geq<pos_t>(
            (i + 1) * ((n - 1) / idx.p_r), 0, r - 1,
            [&](pos_t x) { return (((int64_t) SA_s[pi_[x]]) - 1) % n; });

        idx._D_e[i] = std::make_pair(pi_[x], (pos_t) ((((int64_t) SA_s[pi_[x]]) - 1) % n));
    }

    SA_s.clear();
    SA_s.shrink_to_fit();

    pi_.clear();
    pi_.shrink_to_fit();
}
