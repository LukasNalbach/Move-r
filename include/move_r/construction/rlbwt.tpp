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

#ifndef MOVE_R_MOVE_R_CONSTRUCTION_RLBWT_TPP
#define MOVE_R_MOVE_R_CONSTRUCTION_RLBWT_TPP

#include <filesystem>
#include <stdexcept>

#include <libsais40.hpp>
#include <misc/sort.hpp>
#include <move_r/construction/pfp_output.hpp>
#include <move_r/move_r.hpp>

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::build_sa()
{
    if (!byte_alphabet && idx.symbols_remapped) {
        log_phase_start("mapping T to its effective alphabet");
        int_alphabet->map_text(T_inp, n, p);
        log_phase_end("time_map_t_effective");
        release_int_alphabet();
    }

    log_phase_start("building SA");

    // fs, the extra space after the suffix array that libsais requires for a byte alphabet
    uint64_t fs = byte_alphabet ? 6 * 256 : 0;
    resize_no_init(sa, n + fs);
    move_r_build_sa(T_inp, sa.data(), n, byte_alphabet ? 0 : idx.sigma, fs, p);

    sa.resize(n);

    log_phase_end("time_build_sa");
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::build_rlbwt_c()
{
    log_phase_start(rlbwt_from_pfp ? "reading RLBWT" : "building RLBWT");

    p_ = byte_alphabet ? p : 1;
    // the runs written by the prefix-free parsing (see pfp_rlbwt::run_data_file)
    const pfp_rlbwt::run_data_fields& pfp_fields = pfp_output.fields();
    // each thread reads at least one of them (r >= 1, see pfp_rlbwt::run_data_file::open)
    if (rlbwt_from_pfp) p_ = std::min<uint64_t>(p_, pfp_fields.r);

    r_p.resize(p_ + 1, 0);
    uint8_t width_sym = rlbwt_from_pfp ? std::max<uint8_t>(1, std::bit_width(pfp_fields.max_symbol))
                                       : std::bit_width(idx.sigma - 1);
    uint8_t width_len = rlbwt_from_pfp ? std::bit_width(pfp_fields.max_run_length) : std::bit_width(n);
    RLBWT.resize(p_, interleaved_bit_aligned_vectors({ width_sym, width_len }));

    if (byte_alphabet) {
        C.resize(p_, std::vector<uint64_t>(num_symbol_values(), 0));
    }

    // read the runs written by the prefix-free parsing; thread i_p reads the runs [i_p * r / p', (i_p + 1) * r / p')
    if (rlbwt_from_pfp) {
        n_p.resize(p_ + 1);

        #pragma omp parallel num_threads(p_)
        {
            uint64_t i_p = omp_get_thread_num();
            const uint64_t beg = (i_p * pfp_fields.r) / p_;
            const uint64_t end = ((i_p + 1) * pfp_fields.r) / p_;
            RLBWT[i_p].reserve(end - beg);
            pfp_rlbwt::run_stream runs(pfp_output, beg, end, false, false);
            pfp_rlbwt::run_record rec;

            for (uint64_t x = beg; x < end; x++) {
                runs.next(rec);
                if (x == beg) n_p[i_p] = rec.start;
                uint64_t sym = rec.symbol;

                // the bytes have been mapped to [3..sigma+1] (see min_mapped_symbol); the terminator is 0
                if (sym < move_r_min_mapped_symbol_pfp) [[unlikely]] sym = 0;
                else sym -= move_r_min_mapped_symbol_pfp - move_r_min_mapped_symbol;

                add_run(i_p, sym, rec.length);
                if (byte_alphabet) C[i_p][sym] = C[i_p][sym] + rec.length;
            }

            r_p[i_p] = RLBWT[i_p].size();
            RLBWT[i_p].shrink_to_fit();
        }

        n_p[p_] = pfp_fields.n;
    } else {
        for (uint64_t i = 0; i < p_; i++) {
            n_p.emplace_back(i * (n / p_));
        }

        n_p.emplace_back(n);

        phase_progress prog(params.log, n);
        move_r_bwt_runs_from_sa(T_inp, sa, n_p, r_p, RLBWT, C, byte_alphabet, p_, prog);
        prog.finish();
    }

    // merge the last run in thread i_p's section with the first run in thread i_p+1's section, if their symbols are equal
    for (uint64_t i_p = 0; i_p < p_ - 1; i_p++) {
        uint64_t c = run_sym(i_p, r_p[i_p] - 1);

        if (run_sym(i_p + 1, 0) == c) {
            uint64_t l = run_len(i_p, r_p[i_p] - 1);
            set_run_len(i_p + 1, 0, run_len(i_p + 1, 0) + l);

            if (byte_alphabet) {
                C[i_p][c] = C[i_p][c] - l;
                C[i_p + 1][c] = C[i_p + 1][c] + l;
            }

            n_p[i_p + 1] -= l;
            r_p[i_p]--;
            RLBWT[i_p].resize(r_p[i_p]);
        }
    }

    if (rlbwt_from_pfp) release_pfp_output();

    // build_plcp_run reads T
    if (delete_T && !(has_mphi && !rlbwt_from_pfp)) {
        T_inp.clear();
        T_inp.shrink_to_fit();
    }

    // choose_revert_starts reads sa at the run starts
    if (!builds_SA_samples && !chooses_revert_starts()) {
        sa.clear();
        sa.shrink_to_fit();
    }

    // r_p[i_p] <- the number of runs before thread i_p's section (the index of its first run), r_p[p'] = r
    for (uint64_t i = p_; i > 0; i--) {
        r_p[i] = r_p[i - 1];
    }

    r_p[0] = 0;

    for (uint64_t i = 2; i <= p_; i++) {
        r_p[i] += r_p[i - 1];
    }

    r = r_p[p_];
    idx.r = r;
    // the move data structures have at least r intervals
    mds_check_num_intervals(r);

    if (!byte_alphabet) {
        C.emplace_back(std::vector<uint64_t>(idx.sigma, 0));

        for (uint64_t i = 0; i < r; i++) {
            C[0][run_sym(0, i)] = C[0][run_sym(0, i)] + run_len(0, i);
        }
    }

    compute_C();

    log_phase_end("time_build_rlbwt");
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::compute_C()
{
    uint64_t max_symbol = num_symbol_values();

    // C[i_p][c] <- the number of occurrences of c in the sections 0..i_p
    for (uint64_t i = 1; i < p_; i++) {
        #pragma omp parallel for num_threads(p)
        for (uint64_t j = 0; j < max_symbol; j++) {
            C[i][j] = C[i][j] + C[i - 1][j];
        }
    }

    // shift the rows by one; C[i_p][c] = the number of occurrences of c in L[0..n_p[i_p]-1], C[0][c] = 0
    C.insert(C.begin(), std::vector<uint64_t>(max_symbol, 0));

    // C[p'][c] is the number of occurrences of c in L; C[p'] <- its exclusive prefix sum, the C-array
    uint64_t count_prev = C[p_][0]; // the number of occurrences of the symbol i - 1 in L
    C[p_][0] = 0;

    for (uint64_t i = 1; i < max_symbol; i++) {
        const uint64_t count_cur = C[p_][i];
        C[p_][i] = count_prev + C[p_][i - 1];
        count_prev = count_cur;
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::store_rlbwt()
{
    store_to_tmp_file("RLBWT", ".rlbwt", "time_store_rlbwt", [&](direct_ofstream& file) {
        for (uint64_t i = 0; i < p_; i++) RLBWT[i].serialize(file);
        release_rlbwt();
    });
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::load_rlbwt()
{
    load_from_tmp_file("RLBWT", ".rlbwt", "time_load_rlbwt", true, [&](direct_ifstream& file) {
        RLBWT.resize(p_);
        for (uint64_t i = 0; i < p_; i++) RLBWT[i].load(file);
    });
}

#endif // MOVE_R_MOVE_R_CONSTRUCTION_RLBWT_TPP
