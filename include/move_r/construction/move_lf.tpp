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

#ifndef MOVE_R_MOVE_R_CONSTRUCTION_MOVE_LF_TPP
#define MOVE_R_MOVE_R_CONSTRUCTION_MOVE_LF_TPP

#include <algorithm>
#include <atomic>
#include <filesystem>

#include <move_r/move_r.hpp>

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::build_mlf()
{
    begin_mds_result("build_mlf");
    log_phase_start("\nbuilding M_LF");

    uint8_t width_l_ = std::bit_width(idx.sigma);

    mds_params mlf_params { .num_threads = p, .a = idx.a,
        .log = params.log, .mf = params.mf_mds };
    if (I_LF_sym.width_of(0) == 8) mlf_params.q_order_keys = I_LF_sym.data<uint8_t>();
    // L' is built from I_LF_sym by the construction of M_LF (see mds_pair_row)
    mds_pair_row l_row { .values = &I_LF_sym };
    mlf_params.pair_rows = &l_row;
    idx._M_LF = mlf_t(std::move(I_LF), n, mlf_params, { width_l_ });
    I_LF_sym = interleaved_bit_aligned_vectors();

    r_ = idx._M_LF.num_intervals();
    idx.r_ = r_;

    log_mds_phase_end("time_build_mlf", "r_", r_);
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::build_SA_run_from_iphi_inv()
{
    /* SA_run is built here from I_Phi^{-1} (its pairs contain the SA values at the run boundaries), if the RLBWT
       has been read from the output of the prefix-free parsing and the samples are stored in SA_run */
    const bool build = builds_SA_samples && rlbwt_from_pfp && samples_layout == _at_run_boundaries;
    const bool with_ends = samples_at_run_ends();

    if (build) {
        log_phase_start("building SA_run");
        prepare_sa_run(with_ends);

        // simultaneously iterate over the input intervals of M_LF and the BWT runs
        #pragma omp parallel num_threads(p_)
        {
            uint64_t i_p = omp_get_thread_num();
            uint64_t b = n_p[i_p]; // the first BWT position of thread i_p's section
            uint64_t first_run = r_p[i_p]; // the first run of thread i_p's section
            uint64_t num_runs = r_p[i_p + 1] - r_p[i_p]; // the number of runs in thread i_p's section

            // index of the current input interval of M_LF, initially the one containing b
            uint64_t j = bin_search_max_leq(b, 0, r_ - 1, [&](uint64_t x) { return idx.M_LF().p(x); });

            for (uint64_t i = 0; i < num_runs; i++) {
                const uint64_t first = j; // the first input interval of the run
                idx._SA_run.template set_parallel<0, uint64_t>(idx._SA_run_pos.start_sample_index(j), I_Phi_inv.get<1>(first_run + i));
                j++;

                // skip the input intervals in the run that have been created by the balancing algorithm
                while (j < r_ && !idx.run_starts_at(j)) j++;

                /* the SA value at the end of the run is the first component of the pair of the next run (of the
                   first run, for the last run); a run of length 1 has only a sample at its start (see
                   run_sample_positions) */
                if (with_ends && (j < r_ ? idx.M_LF().p(j) : n) - idx.M_LF().p(first) > 1) {
                    uint64_t next = first_run + i + 1 == r ? 0 : first_run + i + 1;
                    idx._SA_run.template set_parallel<0, uint64_t>(idx._SA_run_pos.end_sample_index(j - 1), I_Phi_inv.get<0>(next));
                }
            }
        }

        log_phase_end("time_build_sa_run");
    }

    n_p.clear();
    n_p.shrink_to_fit();
    r_p.clear();
    r_p.shrink_to_fit();
    release_rlbwt();
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::build_data_structures_over_l_()
{
    // true <=> the smaller symbols in a range of L' are counted with a wavelet tree
    const bool wavelet = is_bidirectional && idx.sigma > max_sigma_num_smaller_psb;

    log_phase_start("building PSO_L'");
    idx._PSO_L_ = pso_l_t(L_reader { &idx }, r_, idx.sigma, params.pso_occs_per_block, p);
    log_phase_end("time_build_pso_l_");

    if (enumerates_extensions) {
        log_phase_start("building PSB_L'");
        idx._PSB_L_ = psb_l_t(L_reader { &idx }, r_, idx.sigma, p, params.pred_succ_block_factor);
        log_phase_end("time_build_psb_l_");
    }

    if (wavelet) {
        log_phase_start("building WTPS_L'");
        // L' and the lengths of the input intervals of M_LF
        no_init_vector<uint32_t> wtps_sym;
        no_init_vector<uint64_t> wtps_weight;
        resize_no_init(wtps_sym, r_);
        resize_no_init(wtps_weight, r_);

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < r_; i++) {
            wtps_sym[i] = uint32_t(idx.L_(i));
            wtps_weight[i] = idx.M_LF().p(i + 1) - idx.M_LF().p(i);
        }

        idx._WTPS_L_ = wavelet_tree_prefix_sum(
            std::move(wtps_sym), std::move(wtps_weight), idx.sigma, p);
        log_phase_end("time_build_wtps_l_");
    }

    if (enumerates_extensions && idx.sigma >= params.min_sigma_L_prev_min) {
        log_phase_start("building L_prev_min");

        /* the tree stores for each i in [0..r'-1] the position of the previous occurrence of L'[i] in L' plus 1
           (0, if there is none, see move_r::_L_prev_min); [0..sigma-1] last[c] = the position of the last occurrence
           of c so far plus 1 */
        const uint8_t pos_width = std::max<uint8_t>(1, std::bit_width(r_));
        interleaved_bit_aligned_vectors last({ pos_width });
        last.resize(idx.sigma, p);

        idx._L_prev_min = L_prev_min_t::build_sequentially(r_, pos_width,
            [&](uint64_t i) {
                uint64_t c = idx.L_(i);
                uint64_t previous = last.template get<0, uint64_t>(c);
                last.template set<0, uint64_t>(c, i + 1);
                return previous;
            }, p);

        log_phase_end("time_build_l_prev_min");
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::build_ilf()
{
    log_phase_start("building I_LF");

    I_LF.reset(r, n);
    const bool byte_symbols = idx.sigma <= 256;
    I_LF_sym = interleaved_bit_aligned_vectors({ byte_symbols ? uint8_t(8) : uint8_t(std::bit_width(idx.sigma)) });
    I_LF_sym.resize_no_init(r, p);
    uint8_t* sym_bytes = I_LF_sym.data<uint8_t>();

    #pragma omp parallel num_threads(p_)
    {
        uint64_t i_p = omp_get_thread_num();
        uint64_t first_run = r_p[i_p]; // the first run of thread i_p's section
        uint64_t num_runs = r_p[i_p + 1] - r_p[i_p]; // the number of runs in it
        uint64_t run_beg = n_p[i_p]; // i', the starting position of the current run

        for (uint64_t i = 0; i < num_runs; i++) {
            /* the pair (i', LF(i')), where
               LF(i') = C[L[i']] + rank(L, L[i'], i' - 1) = C[p'][L[i']] + C[i_p][L[i']] */
            I_LF.set<0>(first_run + i, run_beg);
            I_LF.set<1>(first_run + i, C[p_][run_sym(i_p, i)] + C[i_p][run_sym(i_p, i)]);
            // bit-packed symbols written by different threads can share a word
            if (byte_symbols) sym_bytes[first_run + i] = run_sym(i_p, i);
            else I_LF_sym.template set_parallel<0, uint64_t>(first_run + i, run_sym(i_p, i));

            // C[i_p][c] <- rank(L, c, i' - 1) for the start i' of the next run
            C[i_p][run_sym(i_p, i)] = C[i_p][run_sym(i_p, i)] + run_len(i_p, i);
            run_beg += run_len(i_p, i);
        }
    }

    C.clear();
    C.shrink_to_fit();

    log_phase_end("time_build_ilf");
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::choose_revert_starts()
{
    if (!chooses_revert_starts()) return;
    log_phase_start("choosing the revert starts");

    // [0, n) is split into k ranges of equal length; range g contains the text positions s with s * k / n = g
    const uint64_t k = std::min<uint64_t>(r, params.revert_start_samples);
    auto range_of = [&](uint64_t s) { return uint64_t((uint128_t(s) * k) / n); };
    // returns whether a revert start can have the text position s, i.e., 1 <= s <= n-2 (see move_r::_revert_starts)
    auto is_candidate = [&](uint64_t s) { return s != 0 && s + 2 <= n; };
    const uint64_t threads = std::min<uint64_t>(p, r);

    /* calls visit(x, SA[i]) for each run x in the t-th of the sections [t * r / threads, (t + 1) * r / threads) of
       the runs, where i is the start of run x; SA[i] is read from the output of the prefix-free parsing or from sa */
    auto for_each_run = [&](uint64_t t, auto visit) {
        const uint64_t beg = (t * r) / threads;
        const uint64_t end = ((t + 1) * r) / threads;

        if (rlbwt_from_pfp) {
            pfp_rlbwt::run_stream runs(pfp_output, beg, end, true, false);
            pfp_rlbwt::run_record rec;

            for (uint64_t x = beg; x < end; x++) {
                runs.next(rec);
                visit(x, rec.sa_start);
            }
        } else {
            for (uint64_t x = beg; x < end; x++) visit(x, uint64_t(sa[I_LF.get<0>(x)]));
        }
    };

    // [0..k-1] smallest[g] = the smallest candidate in range g at a run start (n, if there is none)
    std::vector<uint64_t> smallest(k, n);

    #pragma omp parallel for num_threads(threads) schedule(static, 1)
    for (uint64_t t = 0; t < threads; t++) {
        for_each_run(t, [&](uint64_t, uint64_t s) {
            if (!is_candidate(s)) return;
            std::atomic_ref<uint64_t> best(smallest[range_of(s)]);
            uint64_t current = best.load(std::memory_order_relaxed);
            while (s < current && !best.compare_exchange_weak(current, s, std::memory_order_relaxed));
        });
    }

    // [0..k-1] start_of[g] = the run start i with SA[i] = smallest[g] (set if smallest[g] < n)
    no_init_vector<uint64_t> start_of;
    resize_no_init(start_of, k);

    #pragma omp parallel for num_threads(threads) schedule(static, 1)
    for (uint64_t t = 0; t < threads; t++) {
        for_each_run(t, [&](uint64_t x, uint64_t s) {
            if (is_candidate(s) && smallest[range_of(s)] == s) start_of[range_of(s)] = I_LF.get<0>(x);
        });
    }

    const uint8_t width = std::bit_width(n);
    idx._revert_starts = interleaved_bit_aligned_vectors({ width, width });
    idx._revert_starts.resize_no_init(k - std::count(smallest.begin(), smallest.end(), n));
    uint64_t j = 0;

    for (uint64_t g = 0; g < k; g++) {
        if (smallest[g] == n) continue;
        idx._revert_starts.template set<0, uint64_t>(j, start_of[g]);
        idx._revert_starts.template set<1, uint64_t>(j, smallest[g]);
        j++;
    }

    if (rlbwt_from_pfp) {
        release_pfp_output();
    } else if (!builds_SA_samples) {
        // the SA values at the run starts are not read again, if the index does not sample the suffix array
        sa.clear();
        sa.shrink_to_fit();
    }

    log_phase_end("time_choose_revert_starts");
}

#endif // MOVE_R_MOVE_R_CONSTRUCTION_MOVE_LF_TPP
