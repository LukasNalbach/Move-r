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

#ifndef MOVE_R_MOVE_R_CONSTRUCTION_SA_RLZ_TPP
#define MOVE_R_MOVE_R_CONSTRUCTION_SA_RLZ_TPP

#include <move_r/move_r.hpp>

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::SA_rlz_reader::step_right()
{
    auto& phi_inv = build->idx._M_Phi_inv;
    prev_value = value;
    has_prev_value = true;

    phi_inv.move(value, phi_interval);
    pos++;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::SA_rlz_reader::seek(uint64_t target)
{
    enter_run(first_interval_of_run_at(target));

    // start at the last sampled position before target, if it is in the run
    if (sampled != nullptr && (target - pos) >> sampled_shift != 0) {
        pos = (target >> sampled_shift) << sampled_shift;
        value = sampled[target >> sampled_shift];
    }

    locate_phi_interval();
    while (pos < target) step_right();

    has_pos = true;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
uint64_t move_r<support, mlf_enc>::construction::SA_rlz_reader::SAd(uint64_t i)
{
    if (i == 0) {
        if (!has_pos || pos != 0) seek(0);
        return value;
    }

    // move to position i with SA[i - 1] known
    if (!has_pos || pos != i || !has_prev_value) {
        if (has_pos && pos + 1 == i) {
            step_right();
        } else {
            seek(i - 1);
            step_right();
        }
    }

    return value + build->n - prev_value;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::build_SA_rlz() requires(has_SA_rlz)
{
    std::vector<SA_rlz_reader> readers; // [0..p-1] the reader of each thread
    elias_fano_bit_vector mlf_starts; // marks the starting positions of the input intervals of M_LF and position n
    elias_fano_bit_vector phi_inv_starts; // marks the starting positions of the input intervals of M_Phi^{-1} and position n
    no_init_vector<uint64_t> sampled; // the sampled SA values (see SA_rlz_reader::sampled)

    // true <=> SA^d is read from the suffix array (else it is decoded with a temporary M_Phi^{-1}, see SA_rlz_reader)
    const bool sa_in_memory = params.mode == _suffix_array && !sa.empty();

    if (!sa_in_memory) {
        // build_mphi_inv and the rlz construction both consume I_Phi^{-1}
        store_iphi_inv();
        load_iphi_inv(false);
        build_mphi_inv();
        load_iphi_inv();

        mlf_starts = elias_fano_bit_vector(
            [&](uint64_t x) { return x < idx.r_ ? idx._M_LF.p(x) : idx.n; },
            idx.r_ + 1, idx.n + 1, p);
        phi_inv_starts = elias_fano_bit_vector(
            [&](uint64_t x) { return x < idx.r__ ? idx._M_Phi_inv.p(x) : idx.n; },
            idx.r__ + 1, idx.n + 1, p);

        readers.reserve(p);

        for (uint64_t i_p = 0; i_p < p; i_p++)
            readers.emplace_back(*this, mlf_starts, phi_inv_starts, nullptr, 0);

        log_phase_start("sampling SA for the RLZ construction");

        /* SA is sampled at every 2^sampled_shift-th position inside the runs longer than 2^sampled_shift, where
           sampled_shift = max(10, floor(log_2(n/r))). Each of those runs is decoded once from its start to its last
           sampled position. The runs are distributed among the threads by their starting positions. */
        const uint8_t sampled_shift = std::max<uint8_t>(10, uint8_t(std::bit_width(std::max<uint64_t>(1, n / r)) - 1));
        const uint64_t spacing = uint64_t { 1 } << sampled_shift;
        resize_no_init(sampled, ((n - 1) >> sampled_shift) + 1);
        std::vector<uint64_t> first_run(p + 1, idx.r_);
        for (uint64_t i_p = 0; i_p < p; i_p++) first_run[i_p] = readers[0].first_interval_of_run_at(i_p * (n / p));

        #pragma omp parallel for num_threads(p) schedule(static, 1)
        for (uint64_t i_p = 0; i_p < p; i_p++) {
            SA_rlz_reader& reader = readers[i_p];

            for (uint64_t x = first_run[i_p]; x < first_run[i_p + 1];) {
                // next, the input interval of M_LF at which the next run starts
                uint64_t next = x + 1;
                while (next < idx.r_ && !idx.run_starts_at(next)) next++;
                const uint64_t end = next < idx.r_ ? idx._M_LF.p(next) : n;

                if (end - idx._M_LF.p(x) > spacing) {
                    const uint64_t last = ((end - 1) >> sampled_shift) << sampled_shift;
                    reader.enter_run(x);
                    reader.locate_phi_interval();

                    while (reader.pos < last) {
                        reader.step_right();
                        if ((reader.pos & (spacing - 1)) == 0) sampled[reader.pos >> sampled_shift] = reader.value;
                    }
                }

                x = next;
            }
        }

        for (SA_rlz_reader& reader : readers) {
            reader.sampled = sampled.data();
            reader.sampled_shift = sampled_shift;
        }

        log_phase_end("time_sample_sa_rlz");
    }

    // returns SA^d[i] = SA[i] - SA[i-1] + n (SA[0] = n - 1 for i = 0)
    auto SAd = [&, this](uint64_t i_p, uint64_t i) -> uint64_t {
        if (!sa_in_memory) return readers[i_p].SAd(i);
        return i == 0 ? n - 1 : (uint64_t(sa[i]) + n) - uint64_t(sa[i - 1]);
    };

    typename rlz<_diff_access>::template construction<decltype(SAd), mds_pairs> SA_rlz_construction(
        idx._SA_rlz, SAd, I_Phi_inv, n, r,
        {
            .num_threads = p,
            .space = params.mode == _suffix_array ? _in_memory : _on_disk,
            .log = params.log,
            .measurement_file = params.measurement_file,
            .measurement_key_suffix = params.measurement_key_suffix,
            .tmp_files_beside = tmp_files_beside(params),
            .reference_target_length = params.rlz_max_reference_len
        });
    SA_rlz_construction.build();

    if (!sa_in_memory) {
        readers.clear();
        sampled.clear();
        sampled.shrink_to_fit();
        mlf_starts = elias_fano_bit_vector();
        phi_inv_starts = elias_fano_bit_vector();

        // the temporary M_Phi^{-1} is not stored in the index
        idx._M_Phi_inv = mphi_t();
        idx.r__ = 0;
        r__ = 0;
    }
}

#endif // MOVE_R_MOVE_R_CONSTRUCTION_SA_RLZ_TPP
