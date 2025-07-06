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
#include <move_r/move_r.hpp>

/**
 * @brief bi-directional move-r index, size O(r*(a/(a-1)))
 * @tparam support type of locate support (_locate_move or _locate_rlzsa)
 * @tparam sym_t value type (default: char for strings)
 * @tparam pos_t index integer type (use uint32_t if input size < UINT_MAX, else uint64_t)
 */
template <move_r_support support = _locate_move, typename sym_t = char, typename pos_t = uint32_t>
class move_rb {
protected:
    static_assert(support != _locate_one);

    using move_r_fwd_t = constexpr_switch_t< // forward move_r index type
        constexpr_case<support == _locate_move,    move_r<_locate_move_bi_fwd,    sym_t, pos_t>>,
        constexpr_case<support == _locate_rlzsa,   move_r<_locate_rlzsa_bi_fwd,   sym_t, pos_t>>,
     /* constexpr_case<support == _count, */       move_r<_count_bi,              sym_t, pos_t>>;

    using move_r_bwd_t = constexpr_switch_t< // backward move_r index type
        constexpr_case<support == _locate_move,    move_r<_locate_bi_bwd, sym_t, pos_t>>,
        constexpr_case<support == _locate_rlzsa,   move_r<_locate_bi_bwd, sym_t, pos_t>>,
     /* constexpr_case<support == _count, */       move_r<_count_bi,      sym_t, pos_t>>;

    using i_sym_t = move_r_fwd_t::i_sym_t; // internal (unsigned) symbol type
    using inp_t = move_r_fwd_t::inp_t; // input container type

    static_assert(move_r_fwd_t::byte_alphabet);

    static constexpr bool str_input = move_r_fwd_t::str_input;
    static constexpr pos_t max_scan_l_ = move_r_fwd_t::max_scan_l_;
    static constexpr pos_t sample_rate_input_intervals = 4;

    // ############################# INDEX DATA STRUCTURES #############################

    pos_t n = 0;
    pos_t sigma = 0;

    move_r_fwd_t idx_fwd; // move_r index for the forward text
    move_r_bwd_t idx_bwd; // move_r index for the backward text

    sd_array<pos_t> _S_MLF_p_fwd;
    sd_array<pos_t> _S_MLF_p_bwd;

    sd_array<pos_t> _S_MPhi_p;
    sd_array<pos_t> _S_MPhi_m1_p;

    interleaved_byte_aligned_vectors<pos_t, pos_t> _SA_sR_m1;
    interleaved_byte_aligned_vectors<pos_t, pos_t> _SA_eR_m1;

    // ############################# INTERNAL METHODS #############################

    template <direction dir>
    const std::conditional_t<dir == LEFT, move_r_fwd_t, move_r_bwd_t>& index() const
    {
        if constexpr (dir == LEFT) {
            return idx_fwd;
        } else {
            return idx_bwd;
        }
    }

    template <direction dir>
    const sd_array<pos_t>& S_MLF_p() const
    {
        if constexpr (dir == LEFT) {
            return _S_MLF_p_fwd;
        } else {
            return _S_MLF_p_bwd;
        }
    }

    /**
     * @brief constructs a move_r index of the input
     * @tparam bigbwt true <=> use Big-BWT
     * @tparam sa_sint_t suffix array signed integer type
     * @param input the input
     * @param params construction parameters
     */
    template <bool bigbwt, typename sa_sint_t>
    void build(inp_t& input, move_r_params params = {})
    {
        auto time = now();
        auto time_start = time;
        bool log = params.log;
        uint64_t baseline_mem_usage = malloc_count_current();
        std::string prefix_tmp_files = std::filesystem::temp_directory_path().string() +
                "/" + random_alphanumeric_string(10);
        std::string input_file_name;

        if (bigbwt && params.file_input) {
            input_file_name = input;
        } else if (!bigbwt && params.file_input) {
            if (log) std::cout << "reading input file from disk" << std::flush;

            input_file_name = input;
            no_init_resize(input, n);
            std::ifstream input_ifile(input_file_name);
            read_from_file(input_ifile, input.data(), n);
            input_ifile.close();

            if (log) time = log_runtime(time);
        } else if (bigbwt && !params.file_input) {
            if (log) std::cout << "writing input file to disk" << std::flush;

            input_file_name = prefix_tmp_files + ".input";
            std::ofstream input_ofile(input_file_name);
            write_to_file(input_ofile, input.data(), n);
            input_ofile.close();

            if (log) time = log_runtime(time);
        }

        if (log) std::cout << "preprocessing T" << std::flush;

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

        if (log) {
            time = log_runtime(time);
            std::cout << "mapping T to its internal alphabet" << std::flush;
        }

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

        if (log) time = log_runtime(time);

        std::string idx_fwd_file_name;

        std::string sa_fwd_file_name;
        std::vector<sa_sint_t> SA_fwd;
        std::vector<sdsl::int_vector_buffer<40>> SA_fwd_file_bufs;

        uint64_t peak_memory_usage = 0;
        params.peak_memory_usage = &peak_memory_usage;

        if constexpr (bigbwt) {
            reverse(input_file_name, p, false, log);

            if (log) {
                std::cout << std::endl;
                std::cout << "building backward index:" << std::endl;
            }

            move_r_params new_params = params;
            new_params.file_input = true;
            idx_bwd = move_r_bwd_t(input_file_name, new_params);
            idx_bwd.set_alphabet_maps(map_int, map_ext);

            if (log) {
                std::cout << std::endl;
                std::cout << "storing backward index to disk" << std::flush;
                time = now();
            }

            std::string idx_bwd_file_name = prefix_tmp_files + ".move_r_bwd";
            std::ofstream idx_bwd_ofile(idx_bwd_file_name);
            idx_bwd.serialize(idx_bwd_ofile);
            idx_bwd = move_r_bwd_t();
            idx_bwd_ofile.close();
            
            if (log) time = log_runtime(time);

            reverse(input_file_name, p, false, log);
            
            if (log) {
                std::cout << std::endl;
                std::cout << "building forward index:" << std::endl;
            }

            new_params.sa_file_name = &sa_fwd_file_name;
            idx_fwd = move_r_fwd_t(input_file_name, new_params);
            idx_fwd.set_alphabet_maps(map_int, map_ext);

            if (log) {
                std::cout << std::endl;
                time = now();
            }

            if constexpr (support != _count) {
                if (log) std::cout << "storing forward index to disk" << std::flush;

                idx_fwd_file_name = prefix_tmp_files + ".move_r_fwd";
                std::ofstream idx_fwd_ofile(idx_fwd_file_name);
                idx_fwd.serialize(idx_fwd_ofile);
                idx_fwd = move_r_fwd_t();
                idx_fwd_ofile.close();

                if (log) time = log_runtime(time);
            }

            if (log) std::cout << "loading backward index from disk" << std::flush;

            std::ifstream idx_bwd_ifile(idx_bwd_file_name);
            idx_bwd.load(idx_bwd_ifile);
            idx_bwd_ifile.close();
            std::filesystem::remove(idx_bwd_file_name);
            
            if (log) time = log_runtime(time);
        } else {
            reverse(input, p, true, log);

            move_r_params new_params = params;
            new_params.file_input = false;
            idx_bwd = move_r_bwd_t(input, new_params);
            idx_bwd.set_alphabet_maps(map_int, map_ext);

            if (log) {
                std::cout << std::endl;
                time = now();
            }

            reverse(input, p, true, log);
            new_params.sa_vector = reinterpret_cast<void*>(&SA_fwd);
            idx_fwd = move_r_fwd_t(input, new_params);
            idx_fwd.set_alphabet_maps(map_int, map_ext);

            if (log) {
                std::cout << std::endl;
                time = now();
            }
        }


        if constexpr (support != _count) {
            if (log) std::cout << "building H_SA" << std::flush;

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
                    H_SA.emplace((n - idx_bwd.SA_s_(i) - 1), entry_t {.x = i, .type = _run_end   });
                }
            }
            
            if (log) {
                time = log_runtime(time);
                std::cout << "building SA^{-1}_sR and SA^{-1}_eR" << std::flush;
            }

            _SA_sR_m1 = interleaved_byte_aligned_vectors<pos_t, pos_t>({ byte_width(n) });
            _SA_eR_m1 = interleaved_byte_aligned_vectors<pos_t, pos_t>({ byte_width(n) });

            _SA_sR_m1.resize_no_init(r_bwd);
            _SA_eR_m1.resize_no_init(r_bwd);

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
                        _SA_sR_m1.template set_parallel<0, pos_t>(it->second.x, i);
                    }

                    if (it->second.type == _run_end || it->second.type == _run_start_and_end) [[likely]] {
                        _SA_eR_m1.template set_parallel<0, pos_t>(it->second.x, i);
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

            if (log) time = log_runtime(time);

            if constexpr (bigbwt) {
                if (log) std::cout << "loading forward index from disk" << std::flush;

                std::ifstream idx_fwd_ifile(idx_fwd_file_name);
                idx_fwd.load(idx_fwd_ifile);
                idx_fwd_ifile.close();
                std::filesystem::remove(idx_fwd_file_name);
                
                if (log) time = log_runtime(time);
            }

            if constexpr (support == _locate_move) {
                if (log) std::cout << "building S_MPhi_p" << std::flush;

                _S_MPhi_p = build_sampling(
                    idx_fwd.M_Phi().num_intervals(), n, sample_rate_input_intervals,
                    [&](pos_t i){return idx_fwd.M_Phi().p(i);});

                if (log) {
                    time = log_runtime(time);
                    std::cout << "building S_MPhi^{-1}_p" << std::flush;
                }

                _S_MPhi_m1_p = build_sampling(
                    idx_fwd.M_Phi_m1().num_intervals(), n, sample_rate_input_intervals,
                    [&](pos_t i){return idx_fwd.M_Phi_m1().p(i);});

                if (log) time = log_runtime(time);
            }
        }

        if (log) std::cout << "building S_MLF_p_fwd" << std::flush;
        
        _S_MLF_p_fwd = build_sampling(
            idx_fwd.M_LF().num_intervals(), n, sample_rate_input_intervals,
            [&](pos_t i){return idx_fwd.M_LF().p(i);});

        if (log) {
            time = log_runtime(time);
            std::cout << "building S_MLF_p_bwd" << std::flush;
        }

        _S_MLF_p_bwd = build_sampling(
            idx_bwd.M_LF().num_intervals(), n, sample_rate_input_intervals,
            [&](pos_t i){return idx_bwd.M_LF().p(i);});

        if (!bigbwt && params.file_input) {
            input = input_file_name;
        } else if (bigbwt && !params.file_input) {
            std::filesystem::remove(input_file_name);
        }

        if (log) {
            time = log_runtime(time);
            std::cout << "unmapping T from its internal alphabet" << std::flush;
        }

        if (params.mode == _bigbwt) {
            for (i_sym_t c = 255; c >= 2; c--) {
                map_ext[c] = map_ext[c - 2];
            }

            map_ext[0] = 0;
            map_ext[1] = 0;
        }

        if constexpr (bigbwt) {
            std::vector<sdsl::int_vector_buffer<8>> T_buf;

            for (uint16_t i = 0; i < p; i++) {
                T_buf.emplace_back(sdsl::int_vector_buffer<8>(input_file_name, std::ios::in, 128 * 1024, 8, true));
            }

            #pragma omp parallel for num_threads(p)
            for (uint64_t i = 0; i < n - 1; i++) {
                T_buf[omp_get_thread_num()][i] = *reinterpret_cast<i_sym_t*>(&map_ext[T_buf[omp_get_thread_num()][i]]);
            }
        } else {
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

    void reverse(std::string& T, uint16_t p, bool in_memory, bool log = false) {
        auto time = now();
        if (log) std::cout << "reversing T" << std::flush;
        uint64_t n_2 = (n - 1) / 2;

        if (in_memory) {
            #pragma omp parallel for num_threads(p)
            for (uint64_t i = 0; i < n_2; i++) {
                std::swap(T[i], T[(n - i) - 2]);
            }
        } else {
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

        if (log) log_runtime(time);
    }

    template <typename fnc_t>
    sd_array<pos_t> build_sampling(
        pos_t num_values,
        pos_t max_value,
        pos_t sample_rate,
        fnc_t sample
    ) {
        pos_t num_samples = div_ceil(num_values, sample_rate);
        sdsl::sd_vector_builder builder(max_value + 1, num_samples + 1);

        for (pos_t i = 0; i < num_values; i += sample_rate) {
            builder.set(sample(i));
        }

        builder.set(max_value);
        return sd_array<pos_t>(sdsl::sd_vector<>(builder));
    }

    template <typename fnc_t>
    inline static pos_t input_interval(pos_t i, const sd_array<pos_t>& sd_arr, fnc_t interval_start)
    {
        pos_t x = (std::max<pos_t>(1, sd_arr.rank_1(i)) - 1) * sample_rate_input_intervals;

        while (i >= interval_start(x + 1)) {
            x++;
        }

        return x;
    }

    // ############################# CONSTRUCTORS #############################

public:
    move_rb() = default;

    /**
     * @brief constructs a move_r index of the input
     * @param input the input
     * @param params construction parameters
     */
    move_rb(inp_t&& input, move_r_params params = {}) : move_rb(input, params) {}

    /**
     * @brief constructs a move_r index of the input
     * @param input the input
     * @param params construction parameters
     */
    move_rb(inp_t& input, move_r_params params = {})
    {
        if (params.file_input) {
            n = std::filesystem::file_size(input) + 1;
        } else {
            n = input.size() + 1;
        }

        if (params.mode == _suffix_array) {
            build<false, int32_t>(input, params);
        } else {
            if (n <= std::numeric_limits<int32_t>::max()) {
                build<true, int32_t>(input, params);
            } else {
                build<true, int64_t>(input, params);
            }
        }
    }

    // ############################# MISC PUBLIC METHODS #############################

    /**
     * @brief returns the size of the data structure in bytes
     * @return size of the data structure in bytes
     */
    uint64_t size_in_bytes() const
    {
        return
            idx_fwd.size_in_bytes() +
            idx_bwd.size_in_bytes() +
            _S_MLF_p_fwd.size_in_bytes() +
            _S_MLF_p_bwd.size_in_bytes() +
            _S_MPhi_p.size_in_bytes() +
            _S_MPhi_m1_p.size_in_bytes() +
            _SA_sR_m1.size_in_bytes() +
            _SA_eR_m1.size_in_bytes();
    }

    /**
     * @brief logs the index data structure sizes to cout
     */
    void log_data_structure_sizes(bool print_index_size = true) const
    {
        if (print_index_size) std::cout << "index size: " << format_size(size_in_bytes()) << std::endl << std::endl;

        std::cout << "Forward Index: " << std::endl;
        idx_fwd.log_data_structure_sizes(false);
        std::cout << "S_MLF_p: " << format_size(_S_MLF_p_fwd.size_in_bytes()) << std::endl;
        
        if constexpr (support == _locate_move) {
            std::cout << "S_MPhi_p: " << format_size(_S_MPhi_p.size_in_bytes()) << std::endl;
            std::cout << "S_MPhi^{-1}_p: " << format_size(_S_MPhi_m1_p.size_in_bytes()) << std::endl;
        }

        std::cout << std::endl << "Backward Index: " << std::endl;
        idx_bwd.log_data_structure_sizes(false);
        std::cout << "S_MLF_p: " << format_size(_S_MLF_p_bwd.size_in_bytes()) << std::endl;
        std::cout << "SA_^{-1}_sR: " << format_size(_SA_sR_m1.size_in_bytes()) << std::endl;
        std::cout << "SA_^{-1}_eR: " << format_size(_SA_eR_m1.size_in_bytes()) << std::endl;
    }

    /**
     * @brief logs all data structures
     */
    void log_data_structures() const
    {
        if (n > 50) {
            std::cout << "cannot log the contents for n > 50" << std::endl;
            return;
        }

        std::cout << "Forward Index: " << std::endl;
        idx_fwd.log_data_structures();

        std::cout << std::endl << "input interval sampling rate: "
            << sample_rate_input_intervals;

        std::cout << std::endl;
        std::cout << "S_MLF_p: ";
        log_contents(_S_MLF_p_fwd);
        
        if constexpr (support == _locate_move) {
            std::cout << "S_MPhi_p: ";
            log_contents(_S_MPhi_p);

            std::cout << "S_MPhi^{-1}_p: ";
            log_contents(_S_MPhi_m1_p);
        }

        std::cout << std::endl << "Backward Index: " << std::endl;
        idx_bwd.log_data_structures();

        std::cout << std::endl;
        std::cout << "S_MLF_p: ";
        log_contents(_S_MLF_p_bwd);
        
        std::cout << "SA_^{-1}_sR: ";
        log_contents(_SA_sR_m1);

        std::cout << "SA_^{-1}_eR: ";
        log_contents(_SA_eR_m1);
    }

    /**
     * @brief logs the index data structure sizes to the output stream out
     * @param out an output stream
     */
    void log_data_structure_sizes(std::ostream& out) const
    {
    }

    // ############################# PUBLIC ACCESS METHODS #############################

    /**
     * @brief returns a reference to the move_r index of the forward input
     * @return index the forward input
     */
    inline const move_r_fwd_t& forward_index() const
    {
        return idx_fwd;
    }

    /**
     * @brief returns a reference to the move_r index of the backward input
     * @return index the backward input
     */
    inline const move_r_bwd_t& backward_index() const
    {
        return idx_bwd;
    }

    // ############################# QUERY METHODS #############################

    /**
     * @brief stores the variables needed to perform bidirectional pattern search and locate queries
     */
    struct query_context {
    protected:
        // ############################# VARIABLES FOR THE SEARCH PHASE #############################

        direction last_dir; // last performed pattern extend direction
        pos_t m; // length of the currently matched pattern
        pos_t b, e; // [b, e] = SA-interval in T of the currently matched pattern
        pos_t b_R, e_R; // [b, e] = SA-interval in T^R of the currently matched pattern
        pos_t b_, e_; // indices of the input intervals in M_LF containing b and e
        pos_t b_R_, e_R_; // indices of the input intervals in M_LF^R containing b_R and e_R

        // ############################# VARIABLES FOR MAINTAINING SA-SAMPLE INFORMATION DURING SEARCH #############################

        enum sample_t : uint8_t {
            NO_SAMPLE = 0,
            RUN_START = 1,
            RUN_END = 2
        };

        sample_t s_1, s_2; // s_1/2 = RUN_START/RUN_END <=> the usable SA-sample of type 1/2 is at a run start/end in the forward index
        sample_t s_1_R, s_2_R; // s_1/2_R = RUN_START/RUN_END <=> the usable SA-sample of type 1/2 is at a run start/end in the backward index
        pos_t o_1, o_2; // o_1/2 = offset from the left/right of the SA-interval of the usable SA-sample of type 1/2 in the forward index
        pos_t o_1_R, o_2_R; // o_1/2_R = offset from the left/right of the SA-interval of the usable SA-sample of type 1/2 in the backward index
        pos_t i_1, i_2; // i_1/2 = number of iterations since s_1/2 has been set
        pos_t i_1_R, i_2_R; // i_1/2_R = number of iterations since s_1/2_R has been set
        pos_t x_1, x_2; // x_1/2 = index of the (sub-)run in L', whose start/end is the position of the usable SA-sample of type 1/2
        pos_t x_1_R, x_2_R; // x_1/2_R = index of the (sub-)run in L' of T^R, whose start/end is the position of the usable SA-sample of type 1/2

        // ############################# VARIABLES FOR THE LOCATE PHASE #############################
        
        direction locate_dir = NO_DIR; // current locate direction
        pos_t occ_rem; // number of remaining occurrences to locate
        pos_t c; // initial position in the suffix array
        pos_t SA_c; // initial suffix SA[c] in the suffix array interval
        pos_t i; // current position in the suffix array interval
        pos_t SA_i; // current suffix SA[i] in the suffix array interval

        struct empty {};

        std::conditional_t<support == _locate_move, pos_t, empty> s_; // index of the input inteval of M_Phi/M_Phi^{-1} containing SA_i

        struct rlzsa_ctx_t {
            pos_t x_p, x_lp, x_cp, x_r, s_p; // variables for decoding the rlzsa
        };

        std::conditional_t<support == _locate_rlzsa, rlzsa_ctx_t, empty> rlz_l; // rlzsa context for position i
        std::conditional_t<support == _locate_rlzsa, rlzsa_ctx_t, empty> rlz_r; // rlzsa context for position c

        const move_rb<support, sym_t, pos_t>* idx; // index to query

    public:
        /**
         * @brief constructs a new query context for the index idx
         * @param idx an index
         */
        query_context(const move_rb<support, sym_t, pos_t>& idx)
        {
            this->idx = &idx;
            reset();
        }

        /**
         * @brief resets the query context to an empty pattern
         */
        inline void reset()
        {
            last_dir = NO_DIR;
            m = 0;

            b = 0;
            e = idx->idx_fwd.input_size();
            b_R = 0;
            e_R = e;
            b_ = 0;
            e_ = idx->idx_fwd.M_LF().num_intervals() - 1;
            b_R_ = 0;
            e_R_ = idx->idx_bwd.M_LF().num_intervals() - 1;

            s_1 = RUN_START;
            s_2 = RUN_END;
            s_1_R = RUN_START;
            s_2_R = RUN_END;
            o_1 = 0;
            o_2 = 0;
            o_1_R = 0;
            o_2_R = 0;
            i_1 = 0;
            i_2 = 0;
            i_1_R = 0;
            i_2_R = 0;
            x_1 = 0;
            x_2 = e_;
            x_1_R = 0;
            x_2_R = e_R_;
        }

        /**
         * @brief resets the context to the beginning of the locate phase
         */
        void reset_locate()
        {
            occ_rem = e - b + 1;
            locate_dir = NO_DIR;
        }

        /**
         * @brief returns the length of the currently matched pattern
         * @return length of the currently matched pattern
         */
        inline pos_t length() const
        {
            return m;
        }

        /**
         * @brief returns the overall number of occurrences of the currently matched pattern
         * @return overall number of occurrences
         */
        inline pos_t num_occ() const
        {
            return e - b + 1;
        }

        /**
         * @brief returns the number of remaining (not yet reported) occurrences of the currently matched pattern
         * @return number of remaining occurrences
         */
        inline pos_t num_occ_rem() const
        {
            return occ_rem;
        }

        /**
         * @brief returns the suffix array interval of the currently matched pattern
         * @tparam dir text direction
         * @return suffix array interval in the forward text
         */
        template <direction dir>
        inline std::pair<pos_t, pos_t> sa_interval() const
        {
            if constexpr (dir == LEFT) {
                return forward_sa_interval();
            } else {
                return backward_sa_interval();
            }
        }

        /**
         * @brief returns the suffix array interval of the currently matched pattern in the forward text
         * @return suffix array interval in the forward text
         */
        inline std::pair<pos_t, pos_t> forward_sa_interval() const
        {
            return std::make_pair(b, e);
        }

        /**
         * @brief returns the suffix array interval of the currently matched pattern in the backward text
         * @return suffix array interval in the backward text
         */
        inline std::pair<pos_t, pos_t> backward_sa_interval() const
        {
            return std::make_pair(b_R, e_R);
        }

    protected:
        /**
         * @brief extends the currently matched pattern with P; if the extended pattern occurs in the input, true is
         * returned and the query context is adjusted to store the information for the extended pattern; else,
         * false is returned and the query context is not modified
         * @tparam dir extend direction
         * @param sym
         * @return whether the extended pattern occurs in the input
         */
        template <direction dir>
        bool extend(sym_t sym)
        {
            i_sym_t i_sym = idx->idx_fwd.map_symbol(sym);

            // If i_sym does not occur in L', then P[i..m] does not occur in T
            if (i_sym == 0) {
                return false;
            }

            query_context ctx_old = *this;
            bool result;

            if constexpr (dir == LEFT) {
                result = extend<LEFT>(
                    i_sym,
                    b, e, b_R, e_R,
                    b_, e_,
                    s_1, s_2, s_1_R, s_2_R,
                    o_1, o_2, o_1_R, o_2_R,
                    i_1, i_2,
                    x_1, x_2
                );
            } else {
                result = extend<RIGHT>(
                    i_sym,
                    b_R, e_R, b, e,
                    b_R_, e_R_,
                    s_1_R, s_2_R, s_1, s_2,
                    o_1_R, o_2_R, o_1, o_2,
                    i_1_R, i_2_R,
                    x_1_R, x_2_R
                );
            }

            if (result) {
                m++;
                i = 0;
                occ_rem = e - b + 1;
                locate_dir = NO_DIR;
                last_dir = dir;
            } else {
                *this = ctx_old;
            }

            return result;
        }

        /**
         * @brief extends the currently matched pattern with sym; if the extended pattern occurs in the input, true is
         * returned and the query context is adjusted to store the context for the extended pattern; else,
         * false is returned and the query context is not modified; let S be a string, then S(LEFT) = S and S(RIGHT) = S^R
         * and !LEFT = RIGHT = !!RIGHT
         * @tparam dir extend direction
         * @param sym next symbol to match
         * @param b Left interval limit of the suffix array interval of P(dir) in T(dir).
         * @param e Right interval limit of the suffix array interval of P(dir) in T(dir).
         * @param b_r Left interval limit of the suffix array interval of P(!dir) in T(!dir).
         * @param e_r Right interval limit of the suffix array interval of P(!dir) in T(!dir).
         * @param b_ index of the input interval in M_LF of T(dir) containing b.
         * @param e_ index of the input interval in M_LF of T(dir) containing e.
         * @param s_1
         * @param s_2
         * @param s_1_R
         * @param s_2_R
         * @param o_1
         * @param o_2
         * @param i_1
         * @param i_2
         * @param x_1
         * @param x_2
         * @return whether the extended pattern occurs in the input
         */
        template <direction dir>
        bool extend(
            i_sym_t i_sym,
            pos_t& b, pos_t& e, pos_t& b_R, pos_t& e_R,
            pos_t& b_, pos_t& e_,
            sample_t& s_1, sample_t& s_2, sample_t& s_1_R, sample_t& s_2_R,
            pos_t& o_1, pos_t& o_2, pos_t& o_1_R, pos_t& o_2_R,
            pos_t& i_1, pos_t& i_2,
            pos_t& x_1, pos_t& x_2
        ) const {
            const auto& idx_dir = idx->index<dir>();

            if (last_dir != NO_DIR && dir != last_dir) {
                if (!(idx_dir.M_LF().p(b_) <= b && b < idx_dir.M_LF().p(b_ + 1))) {
                    b_ = input_interval(b, idx->S_MLF_p<dir>(), [&](pos_t x){return idx_dir.M_LF().p(x);});
                }

                if (b == idx_dir.M_LF().p(b_)) {
                    s_1 = RUN_START;
                    o_1 = 0;
                    i_1 = 0;
                    x_1 = b_;
                } else if (b == idx_dir.M_LF().p(b_ + 1) - 1) {
                    s_1 = RUN_END;
                    o_1 = 0;
                    i_1 = 0;
                    x_1 = b_;
                }

                if (!(idx_dir.M_LF().p(e_) <= e && e < idx_dir.M_LF().p(e_ + 1))) {
                    e_ = input_interval(e, idx->S_MLF_p<dir>(), [&](pos_t x){return idx_dir.M_LF().p(x);});
                }

                if (e == idx_dir.M_LF().p(e_)) {
                    s_2 = RUN_START;
                    o_2 = 0;
                    i_2 = 0;
                    x_2 = e_;
                } else if (e == idx_dir.M_LF().p(e_ + 1) - 1) {
                    s_2 = RUN_END;
                    o_2 = 0;
                    i_2 = 0;
                    x_2 = e_;
                }
            }
            
            int64_t next[256];
            int64_t prev[256];

            std::fill_n(next, 256, std::numeric_limits<int64_t>::max());
            std::fill_n(prev, 256, std::numeric_limits<int64_t>::min());

            pos_t b_old = b;
            pos_t e_old = e;

            pos_t b__old = b_;
            pos_t e__old = e_;

            pos_t b_R_old = b_R;
            pos_t e_R_old = e_R;

            pos_t blk_size = idx_dir.L_block_size();

            {
                pos_t blk = div_ceil<pos_t>(b_, blk_size);
                int64_t beg = b_;
                int64_t end = std::min<pos_t>(blk * blk_size, e_);

                if (end != e_) {
                    pos_t blk_beg = blk * idx->sigma;
                    pos_t max_sym = i_sym;

                    for (pos_t i = 0; i <= max_sym; i++) {
                        next[i] = idx_dir.L_next(blk_beg + i);
                    }
                }

                for (int64_t i = end; i >= beg; i--) {
                    next[idx_dir.L_(i)] = i;
                }
            }

            {
                pos_t blk = e_ / blk_size;
                int64_t beg = std::max<pos_t>(blk * blk_size, b_);
                int64_t end = e_;

                if (beg != b_) {
                    pos_t blk_beg = blk * idx->sigma;
                    pos_t max_sym = i_sym;

                    for (pos_t i = 0; i <= max_sym; i++) {
                        prev[i] = idx_dir.L_prev(blk_beg + i);
                    }
                }

                for (int64_t i = beg; i <= end; i++) {
                    prev[idx_dir.L_(i)] = i;
                }
            }

            if (next[i_sym] > prev[i_sym]) [[unlikely]] {
                return false;
            }

            b_ = next[i_sym];
            e_ = prev[i_sym];

            if (b_ != b__old) {
                b = idx_dir.M_LF().p(b_);
                s_1 = RUN_START;
                o_1 = 0;
                x_1 = b_;
                i_1 = 1;
            } else if (s_1 == RUN_START) {
                i_1++;
            } else if (s_1 == RUN_END) {
                pos_t p_1 = b + o_1;
                pos_t y_1 = idx_dir.M_LF().p(b_ + 1);

                if (p_1 < y_1) {
                    i_1++;
                } else {
                    o_1 = y_1 - b - 1;
                    x_1 = b_;
                    i_1 = 1;
                }
            }

            if (e_ != e__old) {
                e = idx_dir.M_LF().p(e_ + 1) - 1;
                s_2 = RUN_END;
                o_2 = 0;
                x_2 = e_;
                i_2 = 1;
            } else if (s_2 == RUN_END) {
                i_2++;
            } else if (s_2 == RUN_START) {
                pos_t p_2 = e - o_2;
                pos_t y_2 = idx_dir.M_LF().p(e_);

                if (y_2 <= p_2) {
                    i_2++;
                } else {
                    o_2 = e - y_2;
                    x_2 = e_;
                    i_2 = 1;
                }
            }

            pos_t b_prime = b;
            pos_t e_prime = e;

            pos_t b__prime = b_;
            pos_t e__prime = e_;

            for (i_sym_t i = 0; i < i_sym; i++) {
                if (next[i] <= prev[i] && prev[i] != idx_dir.M_LF().num_intervals()) {
                    pos_t x = next[i];
                    pos_t y = prev[i];

                    b_R += idx_dir.M_LF().p(y + 1) - idx_dir.M_LF().p(y);

                    if (x < y) [[likely]] {
                        b_R += idx_dir.M_LF().q(y) - idx_dir.M_LF().q(x);
                    }
                }
            }

            if (idx_dir.L_(b__old) < i_sym) {
                b_R -= b_old - idx_dir.M_LF().p(b__old);
            }

            if (idx_dir.L_(e__old) < i_sym) {
                b_R -= idx_dir.M_LF().p(e__old + 1) - e_old - 1;
            }

            if (b_ == e_) {
                if (b == e) {
                    idx_dir.M_LF().move(b, b_);
                    e = b;
                    e_ = b_;
                } else {
                    pos_t diff_eb = e - b;
                    idx_dir.M_LF().move(b, b_);
                    e = b + diff_eb;
                    e_ = b_;

                    while (e >= idx_dir.M_LF().p(e_ + 1)) {
                        e_++;
                    }
                }
            } else {
                idx_dir.M_LF().move(b, b_);
                idx_dir.M_LF().move(e, e_);
            }

            e_R = b_R + (e - b);

            if (e - b != e_prime - b_prime) {
                if (b_prime == b_old) {
                    s_1 = RUN_END;
                    o_1 = idx_dir.M_LF().p(b__prime + 1) - 1 - b_prime;
                    x_1 = b__prime;
                    i_1 = 1;
                }

                if (e_prime == e_old) {
                    s_2 = RUN_START;
                    o_2 = e_prime - idx_dir.M_LF().p(e__prime);
                    x_2 = e__prime;
                    i_2 = 1;
                }
            }

            if (s_1_R) {
                pos_t p_1 = b_R_old + o_1_R;

                if (!(b_R <= p_1 && p_1 <= e_R)) {
                    s_1_R = NO_SAMPLE;
                } else {
                    o_1_R -= b_R - b_R_old;
                }
            }

            if (s_2_R) {
                pos_t p_2 = e_R_old - o_2_R;

                if (!(b_R <= p_2 && p_2 <= e_R)) {
                    s_2_R = NO_SAMPLE;
                } else {
                    o_2_R -= e_R_old - e_R;
                }
            }

            assert(s_1 || s_2 || s_1_R || s_2_R);

            return true;
        }

    public:
        /**
         * @brief prepends sym to the currently matched pattern P; if symP occurs in the input, true is
         * returned and the query context is adjusted to store the information for the pattern symP; else,
         * false is returned and the query context is not modified
         * @param sym
         * @return whether symP occurs in the input
         */
        bool prepend(sym_t sym)
        {
            return extend<LEFT>(sym);
        }

        /**
         * @brief appends sym to the currently matched pattern P; if Psym occurs in the input, true is
         * returned and the query context is adjusted to store the information for the pattern Psym; else,
         * false is returned and the query context is not modified
         * @param sym
         * @return whether Psym occurs in the input
         */
        bool append(sym_t sym)
        {
            return extend<RIGHT>(sym);
        }

    protected:
        inline pos_t first_occ()
        {
            if (s_1) {
                i = b + o_1;
                SA_i = (s_1 == RUN_START ? idx->idx_fwd.SA_s(x_1) : idx->idx_fwd.SA_s_(x_1)) - i_1;
            } else if (s_2) {
                i = e - o_2;
                SA_i = (s_2 == RUN_START ? idx->idx_fwd.SA_s(x_2) : idx->idx_fwd.SA_s_(x_2)) - i_2;
            } else if (s_1_R) {
                if (s_1_R == RUN_START) {
                    i = idx->_SA_sR_m1[x_1_R];
                    SA_i = idx->n - idx->idx_bwd.SA_s(x_1_R) - 1 + i_1_R - m;
                } else {
                    i = idx->_SA_eR_m1[x_1_R];
                    SA_i = idx->n - idx->idx_bwd.SA_s_(x_1_R) - 1 + i_1_R - m;
                }

                pos_t i_ = input_interval(i, idx->_S_MLF_p_fwd,
                    [&](pos_t x){return idx->idx_fwd.M_LF().p(x);});

                for (pos_t j = 0; j < m - i_1_R; j++) {
                    idx->idx_fwd.M_LF().move(i, i_);
                }
            } else /* if (s_2_R) */ {
                if (s_2_R == RUN_START) {
                    i = idx->_SA_sR_m1[x_2_R];
                    SA_i = idx->n - idx->idx_bwd.SA_s(x_2_R) - 1 + i_2_R - m;
                } else {
                    i = idx->_SA_eR_m1[x_2_R];
                    SA_i = idx->n - idx->idx_bwd.SA_s_(x_2_R) - 1 + i_2_R - m;
                }

                pos_t i_ = input_interval(i, idx->_S_MLF_p_fwd,
                    [&](pos_t x){return idx->idx_fwd.M_LF().p(x);});

                for (pos_t j = 0; j < m - i_2_R; j++) {
                    idx->idx_fwd.M_LF().move(i, i_);
                }
            }

            c = i;
            SA_c = SA_i;

            if (occ_rem > 1) [[likely]] {
                locate_dir = i > b ? LEFT : RIGHT;

                if constexpr (support == _locate_move) {
                    if (locate_dir == LEFT) {
                        s_ = input_interval(SA_i, idx->_S_MPhi_p,
                            [&](pos_t x){return idx->idx_fwd.M_Phi().p(x);});
                    } else {
                        s_ = input_interval(SA_i, idx->_S_MPhi_m1_p,
                            [&](pos_t x){return idx->idx_fwd.M_Phi_m1().p(x);});
                    }
                } else if constexpr (support == _locate_rlzsa) {
                    if (locate_dir == LEFT) {
                        idx->idx_fwd.init_rlzsa(i, rlz_l.x_p, rlz_l.x_lp, rlz_l.x_cp, rlz_l.x_r, rlz_l.s_p);

                        if (c < e) [[likely]] {
                            rlz_r = rlz_l;
                            pos_t tmp_1 = c;
                            pos_t tmp_2;
                            idx->idx_fwd.next_rlzsa(tmp_1, tmp_2, rlz_r.x_p, rlz_r.x_lp, rlz_r.x_cp, rlz_r.x_r, rlz_r.s_p);
                        }

                        idx->idx_fwd.turn_rlzsa_left(rlz_l.x_p, rlz_l.x_lp, rlz_l.x_cp, rlz_l.x_r, rlz_l.s_p);
                    } else {
                        i++;
                        idx->idx_fwd.init_rlzsa(i, rlz_r.x_p, rlz_r.x_lp, rlz_r.x_cp, rlz_r.x_r, rlz_r.s_p);
                    }
                }
            }

            occ_rem--;
            return SA_i;
        }

    public:
        /**
         * @brief reports the next occurrence of the currently matched pattern
         * @return next occurrence
         */
        inline pos_t next_occ()
        {
            if (locate_dir == NO_DIR) [[unlikely]] {
                return first_occ();
            }

            pos_t occ;

            if constexpr (support == _locate_move) {
                if (locate_dir == LEFT) {
                    idx->idx_fwd.M_Phi().move(SA_i, s_);
                    occ = SA_i;
                    i--;

                    if (i == b) [[unlikely]] {
                        i = c;
                        SA_i = SA_c;
                        locate_dir = RIGHT;

                        s_ = input_interval(SA_i, idx->_S_MPhi_m1_p,
                            [&](pos_t x){return idx->idx_fwd.M_Phi_m1().p(x);});
                    }
                } else /* if (locate_dir == RIGHT) */ {
                    idx->idx_fwd.M_Phi_m1().move(SA_i, s_);
                    occ = SA_i;
                    i++;
                }
            } else if constexpr (support == _locate_rlzsa) {
                if (locate_dir == LEFT) {
                    idx->idx_fwd.prev_rlzsa(i, SA_i, rlz_l.x_p, rlz_l.x_lp, rlz_l.x_cp, rlz_l.x_r, rlz_l.s_p);
                    occ = SA_i;

                    if (i == b && c < e) [[unlikely]] {
                        i = c + 1;
                        locate_dir = RIGHT;
                        SA_i = SA_c;
                    }
                } else /* if (locate_dir == RIGHT) */ {
                    idx->idx_fwd.next_rlzsa(i, SA_i, rlz_r.x_p, rlz_r.x_lp, rlz_r.x_cp, rlz_r.x_r, rlz_r.s_p);
                    occ = SA_i;
                }
            }

            occ_rem--;
            return occ;
        }

        /**
         * @brief locates the remaining (not yet reported) occurrences of the currently matched pattern
         * @param Occ vector to append the occurrences to
         */
        inline void locate(std::vector<pos_t>& Occ)
        {
            Occ.reserve(Occ.size() + occ_rem);

            if (occ_rem == num_occ()) {
                Occ.emplace_back(first_occ());
            }

            if constexpr (support == _locate_move) {
                if (locate_dir == LEFT) {
                    while (i > b) {
                        idx->idx_fwd.M_Phi().move(SA_i, s_);
                        Occ.emplace_back(SA_i);
                        i--;
                    }

                    if (c < e) {
                        i = c;
                        SA_i = SA_c;
                        s_ = input_interval(SA_i, idx->_S_MPhi_m1_p,
                            [&](pos_t x){return idx->idx_fwd.M_Phi_m1().p(x);});
                    }
                }

                if (c < e) {
                    while (i < e) {
                        idx->idx_fwd.M_Phi_m1().move(SA_i, s_);
                        Occ.emplace_back(SA_i);
                        i++;
                    }
                }
            } else if constexpr (support == _locate_rlzsa) {
                if (locate_dir == LEFT) {
                    idx->idx_fwd.report_rlzsa_left(i, b, SA_i,
                        rlz_l.x_p, rlz_l.x_lp, rlz_l.x_cp, rlz_l.x_r, rlz_l.s_p,
                        [&](pos_t, pos_t occ){Occ.emplace_back(occ);}
                    );

                    if (c < e) [[likely]] {
                        i = c + 1;
                        locate_dir = RIGHT;
                        SA_i = SA_c;
                    }
                }

                if (c < e) [[likely]] {
                    idx->idx_fwd.report_rlzsa_right(i, e, SA_i,
                        rlz_r.x_p, rlz_r.x_lp, rlz_r.x_cp, rlz_r.x_r, rlz_r.s_p,
                        [&](pos_t, pos_t occ){Occ.emplace_back(occ);}
                    );
                }
            }

            occ_rem = 0;
        }

        /**
         * @brief locates the remaining (not yet reported) occurrences of the currently matched pattern
         * @return vector containing the occurrences
         */
        std::vector<pos_t> locate()
        {
            std::vector<pos_t> Occ;
            locate(Occ);
            return Occ;
        }
    };

    /**
     * @brief returns a query context for the index
     * @return query_context
     */
    inline query_context query() const
    {
        return query_context(*this);
    }

    // ############################# SERIALIZATION METHODS #############################

    /**
     * @brief stores the index to an output stream
     * @param out output stream to store the index to
     */
    void serialize(std::ostream& out) const
    {
        idx_fwd.serialize(out);
        idx_bwd.serialize(out);

        _S_MLF_p_fwd.serialize(out);
        _S_MLF_p_bwd.serialize(out);

        _S_MPhi_p.serialize(out);
        _S_MPhi_m1_p.serialize(out);

        _SA_sR_m1.serialize(out);
        _SA_eR_m1.serialize(out);
    }

    /**
     * @brief reads a serialized index from an input stream
     * @param in an input stream storing a serialized index
     */
    void load(std::istream& in)
    {
        idx_fwd.load(in);
        idx_bwd.load(in);

        _S_MLF_p_fwd.load(in);
        _S_MLF_p_bwd.load(in);

        _S_MPhi_p.load(in);
        _S_MPhi_m1_p.load(in);

        _SA_sR_m1.load(in);
        _SA_eR_m1.load(in);
    }

    std::ostream& operator>>(std::ostream& os) const
    {
        serialize(os);
        return os;
    }

    std::istream& operator<<(std::istream& is)
    {
        load(is);
        return is;
    }
};