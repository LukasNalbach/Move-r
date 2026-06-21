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

#include <hash_table6.hpp>
#include <move_rb/move_rb.hpp>

template <move_r_support support, typename sym_t, typename pos_t>
template <bool bigbwt, typename sa_sint_t>
void move_rb<support, sym_t, pos_t>::build(inp_t& input, move_r_params params)
{
    auto time = now();
    auto time_start = time;
    bool log = params.log;
    uint64_t baseline_mem_usage = malloc_count_current();
    std::string prefix_tmp_files = "move-r_" + random_alphanumeric_string(10);
    std::string input_file_name;

    // file-/Big-BWT-based input handling only applies to the byte-string input (str_input); the integer-vector
    // input is always provided and processed in-memory
    if constexpr (str_input) {
        if (bigbwt && params.file_input) {
            input_file_name = input;
        } else if (!bigbwt && params.file_input) {
            log_phase_start(log, time, "reading input file from disk");

            input_file_name = input;
            no_init_resize(input, n - 1);
            std::ifstream input_ifile(input_file_name);
            read_from_file(input_ifile, input.data(), n - 1);
            input_ifile.close();

            log_phase_end(log, time);
        } else if (bigbwt && !params.file_input) {
            log_phase_start(log, time, "writing input file to disk");

            input_file_name = prefix_tmp_files + ".input";
            std::ofstream input_ofile(input_file_name);
            write_to_file(input_ofile, input.data(), n - 1);
            input_ofile.close();

            log_phase_end(log, time);
        }
    }

    log_phase_start(log, time, "preprocessing T");

    uint16_t p = params.num_threads;
    std::vector<std::vector<pos_t>> freq_thr(p, std::vector<pos_t>(256, 0));

    if constexpr (bigbwt) {
        std::vector<sdsl::int_vector_buffer<8>> T_buf;

        for (uint16_t i = 0; i < p; i++) {
            T_buf.emplace_back(sdsl::int_vector_buffer<8>(input_file_name, std::ios::in, 128 * 1024, 8, true));
        }

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < n - 1; i++) {
            freq_thr[omp_get_thread_num()][T_buf[omp_get_thread_num()][i]]++;
        }
    } else {
        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < n - 1; i++) {
            freq_thr[omp_get_thread_num()][*reinterpret_cast<i_sym_t*>(&input[i])]++;
        }
    }

    std::vector<pos_t> freq(256, 0);

    for (uint16_t c = 0; c < 256; c++) {
        for (uint16_t i_p = 0; i_p < p; i_p++) {
            freq[c] += freq_thr[i_p][c];
        }
    }

    freq_thr.clear();
    freq_thr.shrink_to_fit();
    std::vector<i_sym_t> uchars_sorted;

    for (uint16_t c = 0; c < 256; c++) {
        if (freq[c] != 0) {
            uchars_sorted.emplace_back(c);
        }
    }

    std::sort(uchars_sorted.begin(), uchars_sorted.end(),
        [&](i_sym_t c1, i_sym_t c2){
            return freq[c1] > freq[c2];
    });

    params.alphabet_size = uchars_sorted.size() + 1;
    sigma = params.alphabet_size;
    i_sym_t cur_uchar = params.mode == _bigbwt ? 3 : 1;
    std::vector<i_sym_t> map_int(256, 0);

    for (i_sym_t c : uchars_sorted) {
        map_int[c] = cur_uchar;
        cur_uchar++;
    }

    log_phase_end(log, time);
    log_message(log, "mapping T to its internal alphabet");

    if constexpr (bigbwt) {
        std::vector<sdsl::int_vector_buffer<8>> T_buf;

        for (uint16_t i = 0; i < p; i++) {
            T_buf.emplace_back(sdsl::int_vector_buffer<8>(input_file_name, std::ios::in, 128 * 1024, 8, true));
        }

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < n - 1; i++) {
            T_buf[omp_get_thread_num()][i] = map_int[T_buf[omp_get_thread_num()][i]];
        }
    } else {
        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < n - 1; i++) {
            input[i] = map_int[*reinterpret_cast<i_sym_t*>(&input[i])];
        }
    }

    std::vector<sym_t> map_ext(256, 0);

    for (uint16_t c = 0; c < 256; c++) {
        if (map_int[c] != 0) {
            if (params.mode == _bigbwt) {
                map_int[c] -= 2;
            }

            map_ext[map_int[c]] = c;
        }
    }

    log_phase_end(log, time);

    std::string idx_fwd_file_name;
    std::string sa_fwd_file_name;
    std::vector<sa_sint_t> SA_fwd;
    std::vector<sdsl::int_vector_buffer<40>> SA_fwd_file_bufs;

    uint64_t peak_memory_usage = 0;
    params.peak_memory_usage = &peak_memory_usage;

    if constexpr (bigbwt && str_input) {
        log_message(log, "building forward index:\n");

        move_r_params new_params = params;
        new_params.file_input = true;
        new_params.sa_file_name = &sa_fwd_file_name;

        idx_fwd = move_r_fwd_t(input_file_name, new_params);
        idx_fwd.set_alphabet_maps(map_int, map_ext);

        std::string old_sa_fwd_file_name = sa_fwd_file_name;
        sa_fwd_file_name = old_sa_fwd_file_name + "_fwd";
        std::filesystem::rename(old_sa_fwd_file_name, sa_fwd_file_name);

        log_phase_start(log, time, "\nstoring forward index to disk");

        idx_fwd_file_name = prefix_tmp_files + ".move_r_fwd";
        std::ofstream idx_fwd_ofile(idx_fwd_file_name);
        idx_fwd.serialize(idx_fwd_ofile);
        idx_fwd = move_r_fwd_t();
        idx_fwd_ofile.close();

        log_phase_end(log, time);

        reverse(input_file_name, p, false, log);

        log_message(log, "building backward index:\n");

        idx_bwd = move_r_bwd_t(input_file_name, new_params);
        idx_bwd.set_alphabet_maps(map_int, map_ext);

        log_phase_start(log, time, "\n");

        reverse(input_file_name, p, false, log);
    } else {
        move_r_params new_params = params;
        new_params.sa_vector = reinterpret_cast<void*>(&SA_fwd);
        new_params.file_input = false;
        
        idx_fwd = move_r_fwd_t(input, new_params);
        idx_fwd.set_alphabet_maps(map_int, map_ext);

        log_phase_start(log, time, "\n");

        reverse(input, p, true, log);

        idx_bwd = move_r_bwd_t(input, new_params);
        idx_bwd.set_alphabet_maps(map_int, map_ext);

        log_phase_start(log, time, "\n");

        reverse(input, p, true, log);
    }

    if constexpr (support != _count) {
        log_phase_start(log, time, "building H_SA");

        enum sample_type : uint8_t {
            _run_start,
            _run_end,
            _run_start_and_end
        };

        struct entry_t {
            pos_t x;
            sample_type type;
        };

        pos_t r_bwd = idx_bwd.M_LF().num_intervals();
        emhash6::HashMap<pos_t, entry_t, std::identity> H_SA;
        H_SA.reserve(2 * r_bwd);

        for (pos_t i = 0; i < r_bwd; i++) {
            if (idx_bwd.M_LF().p(i + 1) - idx_bwd.M_LF().p(i) == 1) {
                H_SA.emplace((n - idx_bwd.SA_s(i) - 1),  entry_t {.x = i, .type = _run_start_and_end });
            } else {
                H_SA.emplace((n - idx_bwd.SA_s(i) - 1),  entry_t {.x = i, .type = _run_start });
                H_SA.emplace((n - idx_bwd.SA_e(i) - 1), entry_t {.x = i, .type = _run_end   });
            }
        }
        
        log_phase_end(log, time);
        log_message(log, "building SA_s_pos and SA_e_pos");

        _SA_s_pos = interleaved_bit_aligned_vectors<pos_t>({ std::bit_width(uint64_t(n)) });
        _SA_e_pos = interleaved_bit_aligned_vectors<pos_t>({ std::bit_width(uint64_t(n)) });

        _SA_s_pos.resize_no_init(r_bwd);
        _SA_e_pos.resize_no_init(r_bwd);

        if constexpr (bigbwt) {
            for (uint16_t i = 0; i < p; i++) {
                SA_fwd_file_bufs.emplace_back(sdsl::int_vector_buffer<40>(
                    sa_fwd_file_name, std::ios::in,
                    128 * 1024, 40, true));
            }
        }

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < n; i++) {
            pos_t SA_fwd_i;

            if constexpr (bigbwt) {
                if (i == 0) [[unlikely]] {
                    SA_fwd_i = n - 1;
                } else {
                    SA_fwd_i = SA_fwd_file_bufs[omp_get_thread_num()][i - 1];
                }
            } else {
                SA_fwd_i = SA_fwd[i];
            }

            auto it = H_SA.find(SA_fwd_i);

            if (it != H_SA.end()) {
                if (it->second.type == _run_start || it->second.type == _run_start_and_end) [[likely]] {
                    _SA_s_pos.template set_parallel<0, pos_t>(it->second.x, i);
                }

                if (it->second.type == _run_end || it->second.type == _run_start_and_end) [[likely]] {
                    _SA_e_pos.template set_parallel<0, pos_t>(it->second.x, i);
                }
            }
        }

        H_SA.clear();
        H_SA.shrink_to_fit();

        if constexpr (bigbwt) {
            SA_fwd_file_bufs.clear();
            SA_fwd_file_bufs.shrink_to_fit();
            std::filesystem::remove(sa_fwd_file_name);
        }

        log_phase_end(log, time);

        if constexpr (bigbwt) {
            log_phase_start(log, time, "loading forward index from disk");

            std::ifstream idx_fwd_ifile(idx_fwd_file_name);
            idx_fwd.load(idx_fwd_ifile);
            idx_fwd_ifile.close();
            std::filesystem::remove(idx_fwd_file_name);
            
            log_phase_end(log, time);
        }

        if constexpr (support == _locate_move) {
            log_phase_start(log, time, "building S_MPhi_p");

            _S_MPhi_p = build_sampling(
                idx_fwd.M_Phi().num_intervals(), n, sample_rate_input_intervals,
                [&](pos_t i){return idx_fwd.M_Phi().p(i);});

            log_phase_end(log, time);
            log_message(log, "building S_MPhi^{-1}_p");

            _S_MPhi_m1_p = build_sampling(
                idx_fwd.M_Phi_m1().num_intervals(), n, sample_rate_input_intervals,
                [&](pos_t i){return idx_fwd.M_Phi_m1().p(i);});

            log_phase_end(log, time);
        }
    }

    log_phase_start(log, time, "building S_MLF_p_fwd");
    
    _S_MLF_p_fwd = build_sampling(
        idx_fwd.M_LF().num_intervals(), n, sample_rate_input_intervals,
        [&](pos_t i){return idx_fwd.M_LF().p(i);});

    log_phase_end(log, time);
    log_message(log, "building S_MLF_p_bwd");

    _S_MLF_p_bwd = build_sampling(
        idx_bwd.M_LF().num_intervals(), n, sample_rate_input_intervals,
        [&](pos_t i){return idx_bwd.M_LF().p(i);});

    if constexpr (str_input) {
        if (!bigbwt && params.file_input) {
            input = input_file_name;
        } else if (bigbwt && !params.file_input) {
            std::filesystem::remove(input_file_name);
        }
    }

    log_phase_end(log, time);
    log_message(log, "unmapping T from its internal alphabet");

    if (params.mode == _bigbwt) {
        for (i_sym_t c = 255; c >= 2; c--) {
            map_ext[c] = map_ext[c - 2];
        }

        map_ext[0] = 0;
        map_ext[1] = 0;
    }

    if (bigbwt && params.file_input) {
        std::vector<sdsl::int_vector_buffer<8>> T_buf;

        for (uint16_t i = 0; i < p; i++) {
            T_buf.emplace_back(sdsl::int_vector_buffer<8>(input_file_name, std::ios::in, 128 * 1024, 8, true));
        }

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < n - 1; i++) {
            T_buf[omp_get_thread_num()][i] = *reinterpret_cast<i_sym_t*>(&map_ext[T_buf[omp_get_thread_num()][i]]);
        }
    } else if (!bigbwt && !params.file_input) {
        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < n - 1; i++) {
            input[i] = map_ext[*reinterpret_cast<i_sym_t*>(&input[i])];
        }
    }

    if (log) {
        time = log_runtime(time);
        peak_memory_usage = std::max(peak_memory_usage, malloc_count_peak() - baseline_mem_usage);
        std::cout << std::endl;
        std::cout << "overall construction time: " << format_time(time_diff_ns(time_start, now())) << std::endl;
        std::cout << "overall peak memory usage: " << format_size(peak_memory_usage) << std::endl;
        log_data_structure_sizes();
        std::cout << std::endl;
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
template <typename seq_t>
void move_rb<support, sym_t, pos_t>::reverse(seq_t& T, uint16_t p, bool in_memory, bool log)
{
    auto time = now();
    log_phase_start(log, time, "reversing T");
    uint64_t n_2 = (n - 1) / 2;

    if (in_memory) {
        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < n_2; i++) {
            std::swap(T[i], T[(n - i) - 2]);
        }
    } else if constexpr (std::is_same_v<seq_t, std::string>) {
        // on-disk reversal (only used for the file-/Big-BWT-based construction, i.e. for str_input)
        std::vector<sdsl::int_vector_buffer<8>> T_buf_left;
        std::vector<sdsl::int_vector_buffer<8>> T_buf_right;

        for (uint16_t i = 0; i < p; i++) {
            T_buf_left.emplace_back(sdsl::int_vector_buffer<8>(T, std::ios::in, 128 * 1024, 8, true));
            T_buf_right.emplace_back(sdsl::int_vector_buffer<8>(T, std::ios::in, 128 * 1024, 8, true));
        }

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < n_2; i++) {
            uint64_t tmp = T_buf_left[omp_get_thread_num()][i];
            T_buf_left[omp_get_thread_num()][i] = T_buf_right[omp_get_thread_num()][(n - i) - 2];
            T_buf_right[omp_get_thread_num()][(n - i) - 2] = tmp;
        }
    }

    log_phase_end(log, time);
}

template <move_r_support support, typename sym_t, typename pos_t>
template <typename fnc_t>
sd_array<pos_t> move_rb<support, sym_t, pos_t>::build_sampling(pos_t num_values, pos_t max_value, pos_t sample_rate, fnc_t sample)
{
    pos_t num_samples = div_ceil(num_values, sample_rate);
    sdsl::sd_vector_builder builder(max_value + 1, num_samples + 1);

    for (pos_t i = 0; i < num_values; i += sample_rate) {
        builder.set(sample(i));
    }

    builder.set(max_value);
    return sd_array<pos_t>(sdsl::sd_vector<>(builder));
}
