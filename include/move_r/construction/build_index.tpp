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

#ifndef MOVE_R_MOVE_R_CONSTRUCTION_BUILD_INDEX_TPP
#define MOVE_R_MOVE_R_CONSTRUCTION_BUILD_INDEX_TPP

#include <libsais40_types.hpp>
#include <misc/serialized_size.hpp>
#include <move_r/move_r.hpp>

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::choose_sa_samples_layout()
{
    samples_layout = _at_run_boundaries;

    /* if the index stores M_Phi^{-1}, the samples are stored as its input intervals (in a bidirectional index
       only if SA is sampled at both ends of each run) */
    if constexpr (support == _locate_move) {
        samples_layout = _as_phi_inv_intervals;
    } else if constexpr (support == _locate_move_bi) {
        if (samples_at_run_ends()) samples_layout = _as_phi_inv_intervals;
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::start_construction()
{
    time = now();
    time_start = time;
    omp_set_num_threads(p);
    prefix_tmp_files = temporary_file_path("move-r_", tmp_files_beside(params));
    baseline_mem_usage = malloc_count_current();

    // the peak of an enclosing construction (e.g., of a move_rb) is restored when this one ends (see peak_scope)
    if (params.log) peak_measurement.emplace();
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::limit_num_threads()
{
    if (p > 1 && move_r_min_symbols_per_thread * p > n) {
        const uint64_t requested = p;
        p = std::max<uint64_t>(1, n / move_r_min_symbols_per_thread);
        p_ = std::min<uint64_t>(p_, p);
        log_note(params.log, "warning: using " + format_threads(p) + " instead of " + format_number(requested)
            + " (at most n/" + std::to_string(move_r_min_symbols_per_thread) + " threads for n = " + format_number(n)
            + ")");
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::log_phase_end(const std::string& mf_key)
{
    ::log_phase_end(params.log, time, params.measurement_file,
        mf_key.empty() ? mf_key : mf_key + params.measurement_key_suffix);

    if (params.log && params.measurement_file != nullptr && !mf_key.empty()) {
        *params.measurement_file << " live_" << mf_key << params.measurement_key_suffix
                       << "=" << idx.size_in_bytes();
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::log_phase_end_quiet(const std::string& mf_key)
{
    if (!params.log) return;
    const std::string key = mf_key + params.measurement_key_suffix;

    if (params.measurement_file != nullptr) {
        *params.measurement_file << " " << key << "=" << time_diff_ns(time, now()) << " peak_" << key << "=" << phase_peak()
                       << " held_" << key << "=" << phase_allocated() << " live_" << key << "=" << idx.size_in_bytes();
    }

    note_peak();
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::log_mds_phase_end(const std::string& time_key, const std::string& count_key, uint64_t count)
{
    if (!params.log) return;
    if (params.mf_mds != nullptr) *params.mf_mds << std::endl;
    if (params.measurement_file != nullptr)
        *params.measurement_file << " " << time_key << params.measurement_key_suffix << "=" << time_diff_ns(time, now())
                       << " peak_" << time_key << params.measurement_key_suffix << "=" << phase_peak()
                       << " " << count_key << params.measurement_key_suffix << "=" << count;
    std::cout << std::endl;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::log_peak_mem_usage()
{
    if (!params.log || !is_malloc_count_linked || params.peak_memory_usage != nullptr) return;
    log_peak_memory_usage("peak memory usage until now", overall_peak() - baseline_mem_usage, n);
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
move_r<support, mlf_enc>::construction::construction(move_r<support, mlf_enc>& index,
    move_r_input_view& T, bool delete_T, move_r_params params_)
    : inp_width(T.symbol_width())
    , byte_alphabet(T.byte_alphabet())
    , T_inp(T)
    , idx(index)
    , params(std::move(params_))
{
    this->delete_T = delete_T;
    p = params.num_threads;
    reversed = params.reversed;
    pfp_prefix = params.pfp_cache ? move_r<support, mlf_enc>::pfp_cache_path(
        params.input_file, params.mode, reversed) : "";
    pfp_cache = !pfp_prefix.empty();
    if (pfp_cache && std::filesystem::exists(pfp_prefix)) pfp_cache_descriptor.open(pfp_prefix);
    idx.byte_alphabet = byte_alphabet;
    idx.sigma = params.alphabet_size;
    idx.a = params.a;
    start_construction();

    const std::string& input_file = params.input_file;
    n = move_r<support, mlf_enc>::input_length(T, params);

    // the suffix array modes build a suffix array with 40-bit entries (see move_r_max_input_length_sa)
    if (params.mode != _pfp && n - 1 > move_r_max_input_length_sa) {
        throw std::runtime_error("the input is too long for the construction modes sa and sa-space (at most 2^"
            + std::to_string(std::bit_width(move_r_max_input_length_sa)) + " - 1 input symbols); use the "
            "construction mode pfp");
    }

    // the positions are stored with at most 48 bits (see move_r_max_input_length)
    if (n - 1 > move_r_max_input_length) {
        throw std::runtime_error("the input is too long (at most 2^"
            + std::to_string(std::bit_width(move_r_max_input_length)) + " - 2 input symbols)");
    }

    idx.n = n;

    if (params.mode != _pfp) {
        if (file_input()) read_t_from_file(input_file);
        else T.push_back(0);

        if (idx.sigma == 0) preprocess_t(input_file);
        build_index();

        if (!delete_T && !file_input()) {
            T.resize(n - 1);
            if (!is_bidirectional) unmap_t();
        }
    } else {
        // without a cache, the RLBWT and its samples are written to files with the prefix prefix_tmp_files
        if (!uses_pfp_cache()) pfp_prefix = prefix_tmp_files;
        if (idx.sigma == 0) preprocess_t(input_file);

        try {
            build_index();
        } catch (...) {
            if (!is_bidirectional && !file_input() && !delete_T) unmap_t();
            throw;
        }

        // a file input has not been remapped (it is read through the alphabet map, see preprocess_t)
        if (!is_bidirectional && !file_input() && !delete_T) unmap_t();
    }

    log_finished();

    if (params.peak_memory_usage != nullptr) {
        *params.peak_memory_usage = std::max(*params.peak_memory_usage, overall_peak() - baseline_mem_usage);
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::build_index()
{
    limit_num_threads();

    // sigma is stored with 32 bits (see serialize)
    if (idx.sigma > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("move_r: the input contains " + std::to_string(idx.sigma - 1) +
            " distinct symbols, but at most 2^32 - 2 are supported");
    }

    // true <=> data structures are stored on disk while they are not needed
    const bool spill = params.mode != _suffix_array;

    if (params.mode == _pfp) {
        rlbwt_from_pfp = true;
        build_rlbwt_pfp();
    }

    if (!rlbwt_from_pfp) build_sa();
    build_rlbwt_c();
    log_statistics();
    choose_sa_samples_layout();

    // after I_LF has been built, the RLBWT is only read to build I_Phi^{-1} from the suffix array
    const bool rlbwt_read_later = !rlbwt_from_pfp;
    const bool spill_rlbwt = spill && rlbwt_read_later;

    build_ilf();
    choose_revert_starts();
    if (!rlbwt_read_later) release_rlbwt();
    if (spill_rlbwt) store_rlbwt();
    build_mlf();
    if (spill_rlbwt) load_rlbwt();

    if (rlbwt_from_pfp && has_mphi) read_lcp_run_pfp();

    if constexpr (!builds_SA_samples) {
        build_SA_run_from_iphi_inv();
        build_data_structures_over_l_();
    } else {
        if (rlbwt_from_pfp) read_iphi_inv_pfp();
        else if constexpr (support == _locate_one_bi) build_sa_run_from_sa();
        else build_iphi_inv_from_sa();

        if constexpr (has_mphi) {
            if (!rlbwt_from_pfp) build_plcp_run();
        }

        build_SA_run_from_iphi_inv();

        if constexpr (has_locate_move) {
            // M_LF and the data structures over L' stay in memory while M_Phi and M_Phi^{-1} are built, if M_Phi is
            // built
            const bool spill_lf = spill && !has_mphi;

            build_data_structures_over_l_();
            if (spill_lf) {
                store_data_structures_over_l_();
                store_mlf();
            }

            if constexpr (has_mphi) {
                build_iphi();
                // sort_iphi reads I_Phi^{-1}, if LCP_run is not empty
                sort_iphi();
                build_mphi();
            }

            sort_iphi_inv();
            if (spill_lf && samples_layout == _as_phi_inv_intervals) store_run_ranks();
            build_mphi_inv();
            plcp_pair = interleaved_bit_aligned_vectors();

            if (spill_lf) load_mlf();
            if (spill_lf) load_data_structures_over_l_();

            if (samples_layout == _as_phi_inv_intervals) build_SA_Phi_inv();
        } else if constexpr (has_SA_rlz) {
            build_data_structures_over_l_();
            if (spill) store_data_structures_over_l_();

            // M_LF and SA_run stay in memory, the rlz construction reads them (see SA_rlz_reader)
            sort_iphi_inv();
            build_SA_rlz();

            if (spill) load_data_structures_over_l_();
        } else {
            // _locate_one and _locate_one_bi store only the samples
            I_Phi_inv.clear();
            build_data_structures_over_l_();
        }
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::log_finished()
{
    if (!params.log) return;
    uint64_t time_construction = time_diff_ns(time_start, now());
    uint64_t peak_mem_usage = overall_peak() - baseline_mem_usage;

    // a sub-index of a larger construction writes its summary only to the measurement file
    if (!is_bidirectional && params.peak_memory_usage == nullptr) {
        std::cout << std::endl;
        log_output_size("index", serialized_size_in_bytes(idx), n);
        idx.log_data_structure_sizes();
        std::cout << std::endl;
        // if the pfp output was read from a cache, the construction time does not include the BWT construction, so
        // no throughput is printed
        log_construction_summary(time_construction, n,
            is_malloc_count_linked ? peak_mem_usage : peak_not_measured, !idx._rlbwt_from_cache);
    }

    if (params.measurement_file != nullptr) {
        const std::string& sfx = params.measurement_key_suffix;
        *params.measurement_file << " time_construction" << sfx << "=" << time_construction;
        *params.measurement_file << " throughput_in_mb_per_s" << sfx << "=" << throughput_in_mb_per_s(n, time_construction);
        *params.measurement_file << " peak_mem_usage" << sfx << "=" << peak_mem_usage;
        idx.log_data_structure_sizes(*params.measurement_file, sfx);
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::log_statistics()
{
    if (!params.log) return;
    double n_r = std::round(100.0 * (n / (double)r)) / 100.0;

    if (params.measurement_file != nullptr) {
        const std::string& sfx = params.measurement_key_suffix;
        *params.measurement_file << " n" << sfx << "=" << n;
        *params.measurement_file << " sigma" << sfx << "=" << std::to_string(idx.sigma);
        *params.measurement_file << " r" << sfx << "=" << r;
    }

    std::cout << "n = " << format_number(n) << std::endl;
    std::cout << "r = " << format_number(r) << std::endl;
    std::cout << "n/r = " << format_ratio(n_r) << std::endl;
    std::cout << "sigma = " << format_number(idx.sigma) << std::endl;
}

#endif // MOVE_R_MOVE_R_CONSTRUCTION_BUILD_INDEX_TPP
