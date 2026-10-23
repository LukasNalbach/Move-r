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

#ifndef MOVE_R_MOVE_R_CONSTRUCTION_MOVE_PHI_TPP
#define MOVE_R_MOVE_R_CONSTRUCTION_MOVE_PHI_TPP

#include <misc/utils.hpp>

#include <misc/sort.hpp>
#include <misc/text_lce.hpp>
#include <move_r/move_r.hpp>

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::build_sa_run_from_sa()
{
    log_phase_start("building SA_run");

    // SA is sampled at the start of each run and at the end of each run longer than 1 (see run_sample_positions)
    prepare_sa_run(true);

    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < r_; i++) {
        const bool starts = idx.run_starts_at(i);
        if (starts) idx._SA_run.template set_parallel<0, uint64_t>(idx._SA_run_pos.start_sample_index(i), sa[idx.M_LF().p(i)]);

        if (idx.run_ends_at(i) && !(starts && idx.M_LF().len(i) == 1))
            idx._SA_run.template set_parallel<0, uint64_t>(idx._SA_run_pos.end_sample_index(i), sa[idx.M_LF().p(i + 1) - 1]);
    }

    log_phase_end("time_build_sa_run");
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::build_iphi_inv_from_sa()
{
    // true <=> SA_run is built (else the samples are stored as input intervals of M_Phi^{-1}, see sa_samples_layout_t)
    const bool at_run_boundaries = builds_SA_samples && samples_layout == _at_run_boundaries;
    log_phase_start(at_run_boundaries ? "building I_Phi^{-1} and SA_run" : "building I_Phi^{-1}");

    I_Phi_inv.reset(r, n);
    const bool with_ends = samples_at_run_ends();
    if (at_run_boundaries) prepare_sa_run(with_ends);

    I_Phi_inv.set<0>(0, sa[n - 1]);
    I_Phi_inv.set<1>(0, sa[0]);

    #pragma omp parallel num_threads(p_)
    {
        uint64_t i_p = omp_get_thread_num();
        uint64_t first_run = r_p[i_p]; // the first run of thread i_p's section
        uint64_t num_runs = r_p[i_p + 1] - r_p[i_p]; // the number of runs in it
        uint64_t j = n_p[i_p]; // the starting position of the current BWT run
        // the input interval of M_LF containing j
        uint64_t x = bin_search_max_leq(j, 0, r_ - 1, [&](uint64_t y) { return idx.M_LF().p(y); });

        for (uint64_t i = 0; i < num_runs; i++) {
            if (j != 0) [[likely]] {
                I_Phi_inv.set<0>(first_run + i, sa[j - 1]);
                I_Phi_inv.set<1>(first_run + i, sa[j]);
            }

            const uint64_t k = j; // the starting position of the run
            j += run_len(i_p, i);
            const uint64_t first = x; // the first input interval of the run
            while (x < r_ && idx.M_LF().p(x) < j) x++;

            if (at_run_boundaries) {
                idx._SA_run.template set_parallel<0, uint64_t>(idx._SA_run_pos.start_sample_index(first), sa[k]);
                // a run of length 1 has only a sample at its start (see run_sample_positions)
                if (with_ends && j - k > 1)
                    idx._SA_run.template set_parallel<0, uint64_t>(idx._SA_run_pos.end_sample_index(x - 1), sa[j - 1]);
            }
        }
    }

    // only the in-memory rlz construction reads the suffix array again (see build_SA_rlz)
    if constexpr (!has_SA_rlz) {
        sa.clear();
        sa.shrink_to_fit();
    }

    log_phase_end("time_build_iphi_inv");
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::build_mphi_inv()
{
    begin_mds_result("build_mphi_inv");
    log_phase_start("\nbuilding M_Phi^{-1}");

    /* the PLCP row stores the PLCP value at the start of the output interval of each input interval; plcp_pair is
       sorted by the starting positions of the output intervals of the pairs (see plcp_pair) */
    mds_pair_row plcp_row { .values = plcp_pair.empty() ? nullptr : &plcp_pair, .in_q_order = true, .subtract_offset = true };
    idx._M_Phi_inv = mphi_t(
        std::move(I_Phi_inv), n, {
            .num_threads = p,
            .a = idx.a,
            .pair_rows = &plcp_row,
            // marks the first input interval of each pair (see build_SA_Phi_inv)
            .first_pieces = samples_layout == _as_phi_inv_intervals ? &phi_inv_first_intervals : nullptr,
            .log = params.log,
            .mf = params.mf_mds,
        },
        idx.mphi_inv_row_widths(n), {}, nullptr);

    r__ = idx._M_Phi_inv.num_intervals();
    idx.r__ = r__;

    log_mds_phase_end("time_build_mphi_inv", "r__", r__);
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::store_run_ranks()
{
    store_to_tmp_file("the run ranks", ".run_ranks", "time_store_run_ranks", [&](direct_ofstream& file) {
        write_in_chunks(file, run_ranks.data(), run_ranks.size_in_bytes());
        run_ranks.clear();
    });

    run_ranks_on_disk = true;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::build_SA_Phi_inv()
{
    log_phase_start("building SA_Phi^{-1}");

    /* SA_Phi^{-1}[x] = index of the first input interval of M_Phi^{-1} of the pair of the run starting at the x-th
       input interval of M_LF, or r'', if no run starts there (see move_r::_SA_Phi_inv). The pair of a run with rank
       i (see run_ranks) starts at the i-th interval marked in phi_inv_first_intervals. */
    idx.omega_idx = idx._M_Phi_inv.width_idx();
    idx._SA_Phi_inv = interleaved_bit_aligned_vectors({ idx.omega_idx });
    idx._SA_Phi_inv.resize_no_init(r_, p);
    // number of ranges of input intervals of M_LF (at most one per 2^16 runs)
    const uint64_t num_ranges = std::clamp<uint64_t>(r >> 16, 1, p);
    // [0..num_ranges] first_run[t] = index of the first run starting in range t; first_run[num_ranges] = r
    std::vector<uint64_t> first_run(num_ranges + 1, 0);

    // count the run starts in each range, then compute the prefix sums
    #pragma omp parallel for num_threads(p) schedule(static, 1)
    for (uint64_t t = 0; t < num_ranges; t++) {
        uint64_t runs = 0;
        for (uint64_t x = (t * r_) / num_ranges; x < ((t + 1) * r_) / num_ranges; x++) runs += idx.run_starts_at(x);
        first_run[t + 1] = runs;
    }

    for (uint64_t t = 0; t < num_ranges; t++) first_run[t + 1] += first_run[t];
    const std::string file_name = prefix_tmp_files + ".run_ranks";
    // bytes per rank (also in the file, see store_run_ranks)
    const uint64_t rank_bytes = bytes_per_value(r);

    // iterate over the input intervals of each range and the ranks of its runs; ranks on disk are read in blocks
    #pragma omp parallel for num_threads(p) schedule(static, 1)
    for (uint64_t t = 0; t < num_ranges; t++) {
        static constexpr uint64_t block_length = 1 << 14; // number of ranks read from the file at once
        uint64_t i = first_run[t]; // index of the next run
        uint64_t num_left = first_run[t + 1] - i; // number of ranks of the range that have not been read from the file
        uint_bytes_array<1> block;
        std::optional<direct_ifstream> file;
        uint64_t block_pos = 0; // index of the next rank in block

        if (run_ranks_on_disk && num_left != 0) {
            file.emplace(file_name, std::ios::binary, num_left * rank_bytes);
            file->seekg(i * rank_bytes);
        }

        for (uint64_t x = (t * r_) / num_ranges; x < ((t + 1) * r_) / num_ranges; x++) {
            uint64_t entry = r__;

            if (idx.run_starts_at(x)) {
                if (run_ranks_on_disk) {
                    if (block_pos == block.size()) {
                        block.reset(std::min<uint64_t>(num_left, block_length), r);
                        read_in_chunks(*file, block.data(), block.size_in_bytes());
                        num_left -= block.size();
                        block_pos = 0;
                    }

                    entry = phi_inv_first_intervals.select(block.get(block_pos++));
                } else {
                    entry = phi_inv_first_intervals.select(run_ranks.get(i++));
                }
            }

            idx._SA_Phi_inv.template set_parallel<0, uint64_t>(x, entry);
        }
    }

    if (run_ranks_on_disk) std::filesystem::remove(file_name);
    run_ranks_on_disk = false;
    run_ranks.clear();
    phi_inv_first_intervals.clear();
    log_phase_end("time_build_saphi_inv");
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::build_mphi()
{
    begin_mds_result("build_mphi");
    log_phase_start("\nbuilding M_Phi");

    idx._M_Phi = mphi_t(
        std::move(I_Phi), n, {
            .num_threads = p,
            .a = idx.a,
            .log = params.log,
            .mf = params.mf_mds,
        },
        {}, {}, nullptr);

    r___ = idx._M_Phi.num_intervals();
    idx.r___ = r___;

    log_mds_phase_end("time_build_mphi", "r___", r___);
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::store_iphi_inv()
{
    store_to_tmp_file("I_Phi^{-1}", ".iphi_inv", "time_store_iphi_inv", [&](direct_ofstream& file) {
        const uint64_t iphi_inv_bytes = I_Phi_inv.size_in_bytes();
        phase_progress prog(params.log, iphi_inv_bytes);
        write_in_chunks(file, I_Phi_inv.data(), iphi_inv_bytes, &prog);
        prog.finish();
        I_Phi_inv.clear();
    });
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::load_iphi_inv(bool remove_file)
{
    load_from_tmp_file("I_Phi^{-1}", ".iphi_inv", "time_load_iphi_inv", remove_file, [&](direct_ifstream& file) {
        I_Phi_inv.reset(r, n);
        const uint64_t iphi_inv_bytes = I_Phi_inv.size_in_bytes();
        phase_progress prog(params.log, iphi_inv_bytes);
        read_in_chunks(file, I_Phi_inv.data(), iphi_inv_bytes, &prog);
        prog.finish();
    });
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::build_iphi()
{
    log_phase_start("\nbuilding I_Phi");

    I_Phi.reset(r, n);
    // true <=> the second component of each pair is the index of its run until the pairs are sorted (see sort_iphi)
    const bool with_runs = !LCP_run.empty();

    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < r; i++) {
        I_Phi.set<0>(i, I_Phi_inv.get<1>(i));
        I_Phi.set<1>(i, with_runs ? i : I_Phi_inv.get<0>(i));
    }

    log_phase_end("time_build_iphi");
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::build_plcp_run()
{
    log_phase_start("building the run-boundary LCP values");

    /* The pair of the x-th run is (SA[i-1], SA[i]), where i is the starting position of the run in L; LCP[i] is the
       length of the longest common prefix of the suffixes of T starting at SA[i-1] and SA[i]. The terminator occurs
       only once in T, so no comparison reads beyond the end of T. */
    if (T_inp.size() < n)
        throw std::runtime_error("move_r: T must be in memory to compute the LCP values at the run boundaries");
    const char* text = T_inp.data();
    const uint8_t sym_width = T_inp.symbol_width();
    LCP_run = interleaved_bit_aligned_vectors({ uint8_t(std::bit_width(n)) });
    LCP_run.resize_no_init(r, p);
    uint64_t max_val = 0;

    #pragma omp parallel for num_threads(p) schedule(dynamic, 1024) reduction(max : max_val)
    for (uint64_t x = 0; x < r; x++) {
        const uint64_t i = I_Phi_inv.get<0>(x);
        const uint64_t j = I_Phi_inv.get<1>(x);
        const uint64_t value = text_lce<RIGHT>(text, sym_width, i, j, (n - 1) - std::max(i, j));
        LCP_run.template set_parallel<0, uint64_t>(x, value);
        max_val = std::max(max_val, value);
    }

    // the maximum value in LCP_run is the maximum PLCP value, which determines the width of the PLCP rows
    idx._plcp_row_width = std::max<uint8_t>(1, std::bit_width(max_val));

    if (delete_T) {
        T_inp.clear();
        T_inp.shrink_to_fit();
    }

    log_phase_end("time_build_plcp_run");
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::sort_iphi_inv()
{
    log_phase_start("sorting I_Phi^{-1}");

    if (samples_layout == _as_phi_inv_intervals) {
        /* the pairs (p, x), where p is the starting position of the input interval of the x-th run, are sorted by p;
           then the rank of each pair is stored in the run ranks (see run_ranks), and x is replaced by the starting
           position q of the output interval, which is read from I_Phi^{-1} (still in the order of the runs) */
        mds_pairs sorted(r, n);

        #pragma omp parallel for num_threads(p)
        for (uint64_t x = 0; x < r; x++) {
            sorted.set<0>(x, I_Phi_inv.get<0>(x));
            sorted.set<1>(x, x);
        }

        sorted.sort_by<0>(p);

        run_ranks.reset(r, r);
        advise_huge_pages(run_ranks.data(), run_ranks.size_in_bytes());

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < r; i++) {
            const uint64_t x = sorted.get<1>(i);
            run_ranks.set(x, i);
            sorted.set<1>(i, I_Phi_inv.get<1>(x));
        }

        I_Phi_inv = std::move(sorted);
    } else {
        I_Phi_inv.sort_by<0>(p);
    }

    log_phase_end("time_sort_iphi_inv");
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::sort_iphi()
{
    log_phase_start("sorting I_Phi");

    I_Phi.sort_by<0>(p);

    /* the second component of each pair is the index x of its run (see build_iphi); it is replaced by the first
       component of the x-th pair in I_Phi^{-1}, and LCP_run[x] is written to plcp_pair */
    if (!LCP_run.empty()) {
        plcp_pair = interleaved_bit_aligned_vectors({ idx._plcp_row_width });
        plcp_pair.resize_no_init(r, p);

        #pragma omp parallel for num_threads(p)
        for (uint64_t k = 0; k < r; k++) {
            const uint64_t x = I_Phi.get<1>(k);
            I_Phi.set<1>(k, I_Phi_inv.get<0>(x));
            plcp_pair.template set_parallel<0, uint64_t>(k, LCP_run.template get<0, uint64_t>(x));
        }

        LCP_run = interleaved_bit_aligned_vectors();
    }

    log_phase_end("time_sort_iphi");
}

#endif // MOVE_R_MOVE_R_CONSTRUCTION_MOVE_PHI_TPP
