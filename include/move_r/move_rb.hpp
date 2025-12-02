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
#include <tsl/sparse_set.h>
#include <move_r/move_r.hpp>
#include <misc/search_schemes.hpp>

enum move_rb_query_support_t {
    COUNT = 0,
    LOCATE = 1
};

/**
 * @brief bi-directional move-r index, size O(r*(a/(a-1)))
 * @tparam support type of locate support (_locate_move or _locate_rlzsa)
 * @tparam sym_t value type (default: char for strings)
 * @tparam pos_t index integer type (use uint32_t if input size < UINT_MAX, else uint64_t)
 */
template <move_r_support support = _locate_move, typename sym_t = char, typename pos_t = uint32_t>
class move_rb {
    public:
    static_assert(support != _locate_one);

    using move_r_fwd_t = constexpr_switch_t< // forward move_r index type
        constexpr_case<support == _locate_move,    move_r<_locate_move_bi_fwd,    sym_t, pos_t>>,
        constexpr_case<support == _locate_rlzsa,   move_r<_locate_rlzsa_bi_fwd,   sym_t, pos_t>>,
     /* constexpr_case<support == _count, */       move_r<_count_bi,              sym_t, pos_t>>;

    using move_r_bwd_t = constexpr_switch_t< // backward move_r index type
        constexpr_case<support == _locate_move,    move_r<_locate_bi_bwd, sym_t, pos_t>>,
        constexpr_case<support == _locate_rlzsa,   move_r<_locate_bi_bwd, sym_t, pos_t>>,
     /* constexpr_case<support == _count, */       move_r<_count_bi,      sym_t, pos_t>>;

    static constexpr bool supports_locate = move_r_fwd_t::supports_locate; // true <=> the index supports locate
    static constexpr bool supports_multiple_locate = move_r_fwd_t::supports_multiple_locate; // true <=> the index supports locating multiple occurrences
    static constexpr bool supports_bwsearch = move_r_fwd_t::supports_bwsearch; // true <=> the index uses backward search for answering count queries
    static constexpr bool has_rlzsa = support == move_r_fwd_t::has_rlzsa; // true <=> the index has an rlzsa
    static constexpr bool has_lzendsa = support == move_r_fwd_t::has_lzendsa; // true <=> the index has an lzendsa
    static constexpr bool has_locate_move = move_r_fwd_t::has_locate_move; // true the index uses move data structures for Phi/Phi^{-1}
    static constexpr bool str_input = move_r_fwd_t::str_input; // true <=> the input is a string
    static constexpr bool int_input = move_r_fwd_t::int_input; // true <=> the input is an iteger vector
    static constexpr bool byte_alphabet = move_r_fwd_t::byte_alphabet; // true <=> the input uses a byte alphabet
    static constexpr bool int_alphabet = move_r_fwd_t::int_alphabet; // true <=> the input uses an integer alphabet
    static constexpr pos_t max_scan_l_ = move_r_fwd_t::max_scan_l_; // maximum distance to scan over L' to find the first and last occurrences of sym in L'[\hat{b},\hat{e}]
    static constexpr pos_t sample_rate_input_intervals = 4; // sample rate of the sd-arrays storing input interval starting positions

    static_assert(byte_alphabet);

    using map_int_t = move_r_fwd_t::map_int_t; // type of map_int
    using map_ext_t = move_r_fwd_t::map_ext_t; // type of map_ext
    using i_sym_t = move_r_fwd_t::i_sym_t; // internal (unsigned) symbol type
    using inp_t = move_r_fwd_t::inp_t; // input container type

    // ############################# INDEX DATA STRUCTURES #############################

protected:
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
            if (log) std::cout << "building forward index:" << std::endl;

            move_r_params new_params = params;
            new_params.file_input = true;
            new_params.sa_file_name = &sa_fwd_file_name;

            idx_fwd = move_r_fwd_t(input_file_name, new_params);
            idx_fwd.set_alphabet_maps(map_int, map_ext);

            std::string old_sa_fwd_file_name = sa_fwd_file_name;
            sa_fwd_file_name = prefix_tmp_files + ".sa_fwd";
            std::filesystem::rename(old_sa_fwd_file_name, sa_fwd_file_name);

            if (log) {
                std::cout << std::endl;
                time = now();
            }

            if (log) std::cout << "storing forward index to disk" << std::flush;

            idx_fwd_file_name = prefix_tmp_files + ".move_r_fwd";
            std::ofstream idx_fwd_ofile(idx_fwd_file_name);
            idx_fwd.serialize(idx_fwd_ofile);
            idx_fwd = move_r_fwd_t();
            idx_fwd_ofile.close();

            if (log) time = log_runtime(time);

            reverse(input_file_name, p, false, log);

            if (log) std::cout << "building backward index:" << std::endl;

            idx_bwd = move_r_bwd_t(input_file_name, new_params);
            idx_bwd.set_alphabet_maps(map_int, map_ext);

            if (log) {
                std::cout << std::endl;
                time = now();
            }

            reverse(input_file_name, p, false, log);
        } else {
            move_r_params new_params = params;
            new_params.sa_vector = reinterpret_cast<void*>(&SA_fwd);
            new_params.file_input = false;
            
            idx_fwd = move_r_fwd_t(input, new_params);
            idx_fwd.set_alphabet_maps(map_int, map_ext);

            if (log) {
                std::cout << std::endl;
                time = now();
            }

            reverse(input, p, true, log);

            idx_bwd = move_r_bwd_t(input, new_params);
            idx_bwd.set_alphabet_maps(map_int, map_ext);

            if (log) {
                std::cout << std::endl;
                time = now();
            }

            reverse(input, p, true, log);
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

    /**
     * @brief reverses T either in-memory or on disk
     * @param T string storing either T (if in_memory = true) or the name of the file containing T, else
     * @param p the number of threads to use
     * @param in_memory true <=> T is stored in-memory
     * @param log true <=> print log messages
     */
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

    /**
     * @brief computes an sd-array marking every sample_rate-th sample of all num_sample samples in the range [0, max_value)
     * @param num_values overall number of values
     * @param max_value maximum value
     * @param sample_rate sample rate
     * @param sample function, where sample(i) returns the ith sample
     */
    template <typename fnc_t>
    sd_array<pos_t> build_sampling(pos_t num_values, pos_t max_value, pos_t sample_rate, fnc_t sample)
    {
        pos_t num_samples = div_ceil(num_values, sample_rate);
        sdsl::sd_vector_builder builder(max_value + 1, num_samples + 1);

        for (pos_t i = 0; i < num_values; i += sample_rate) {
            builder.set(sample(i));
        }

        builder.set(max_value);
        return sd_array<pos_t>(sdsl::sd_vector<>(builder));
    }

    /**
     * @brief computes the index of the input interval containing position i
     * @param i position in [0, n)
     * @param sd_arr an sd-array marking exactly the starting positions of each x-th input interval of some move data structure
     * @param interval_start function, where interval_start(i) returns the starting position of the ith input interval of the same move data structure
     */
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
     * @brief maps a symbol to its corresponding symbol in the internal effective alphabet
     * @param sym symbol
     * @return its corresponding symbol in the internal effective alphabet
     */
    inline i_sym_t map_symbol(sym_t sym) const
    {
        return idx_fwd.map_symbol(sym);
    }

    /**
     * @brief maps a symbol that occurs in the internal effective alphabet to its corresponding
     *        symbol in the input
     * @param sym a symbol that occurs in the internal effective alphabet
     * @return its corresponding symbol in the input
     */
    inline sym_t unmap_symbol(i_sym_t sym) const
    {
        return idx_fwd.unmap_symbol(sym);
    }

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
     * @param print_index_size true <=> also prints the overall index size
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

    struct locate_context_t;
    struct extend_context_t;

    enum sample_t : uint8_t {
        NO_SAMPLE = 0,
        RUN_START = 1,
        RUN_END = 2
    };

    struct __attribute__((packed)) sample_info_pack_t {
        sample_t t; // RUN_START/RUN_END <=> the SA-sample is at a run start/end
        pos_t o; // offset from the left/right of the SA-interval of the SA-sample
        pos_t i; // number of iterations since s_idx has been set
        pos_t x; // index of the (sub-)run in L', whose start/end is the position of the SA-sample
    };

    /**
     * @brief stores the variables needed to perform bidirectional pattern search and locate queries
     */
    template <move_rb_query_support_t query_support>
    struct search_context_t {
        static_assert(!(support == _count && query_support == LOCATE));

        friend class locate_context_t;

    protected:
        // ############################# VARIABLES FOR THE SEARCH PHASE #############################
        
        direction dir_lst; // last performed P extend direction
        sym_t sym_lst; // last-added symbol
        pos_t m; // length of the currently matched P
        pos_t b, e; // [b, e] = SA-interval in T of the currently matched P
        pos_t b_R, e_R; // [b_R, e_R] = SA-interval in T^R of the currently matched P^R
        pos_t b_, e_; // indices of the input intervals in M_LF containing b and e
        pos_t b_R_, e_R_; // indices of the input intervals in M_LF^R containing b_R and e_R

        // ############################# VARIABLES FOR MAINTAINING SA-SAMPLE INFORMATION DURING SEARCH #############################

        using sample_info_t = std::conditional_t<query_support == LOCATE, sample_info_pack_t, empty_t>;
        
        sample_info_t s_b, s_e; // sample info for the beginning/end of the SA-interval in T of the currently matched P
        sample_info_t s_b_R, s_e_R; // sample info for the beginning/end of the SA-interval in T^R of the currently matched P^R

        template <direction dir> inline pos_t& beg() {if constexpr (dir == LEFT) {return b;} else {return b_R;}}
        template <direction dir> inline pos_t& end() {if constexpr (dir == LEFT) {return e;} else {return e_R;}}
        template <direction dir> inline pos_t& beg_run() {if constexpr (dir == LEFT) {return b_;} else {return b_R_;}}
        template <direction dir> inline pos_t& end_run() {if constexpr (dir == LEFT) {return e_;} else {return e_R_;}}
        template <direction dir> inline sample_info_t& smpl_beg() requires(query_support == LOCATE) {if constexpr (dir == LEFT) {return s_b;} else {return s_b_R;}}
        template <direction dir> inline sample_info_t& smpl_end() requires(query_support == LOCATE) {if constexpr (dir == LEFT) {return s_e;} else {return s_e_R;}}

    public:
        search_context_t() {}

        /**
         * @brief constructs a new search context with last-added symbol sym
         * @param sym the last-added symbol
         */
        search_context_t(sym_t sym)
        {
            sym_lst = sym;
        }
        
        /**
         * @brief constructs a new query context for the index idx
         * @param idx an index
         */
        search_context_t(const move_rb<support, sym_t, pos_t>& idx)
        {
            reset(idx);
        }

        /**
         * @brief equality operator on search contexts
         * @param other the other search context to test equality with
         * @return whether the other context other is equal to this search context; it is considered equal
         *         if it represents the same pattern (SA-interval); other context variables may vary
         */
        bool operator==(const search_context_t& other) const
        {
            return b == other.b && e == other.e;
        }

        /**
         * @brief hash function used to identify search contexts; the hash is computed only from the beginning of the forward-SA-interval,
         *        because the SA-intervals of each two patterns of the same length are either identical or disjoint
         */
        struct hash {
            /**
             * @brief function returning a hash for a given search context
             * @param ctx the given search context
             * @return the hash computed from just the beginning of the interval
             */
            inline std::size_t operator()(const search_context_t& ctx) const noexcept
            {
                return pos_hash<pos_t>(ctx.b);
            }
        };

        /**
         * @brief returns the last-added symbol
         * @return the last-added symbol
         */
        inline sym_t last_added_symbol() const
        {
            return sym_lst;
        }

        /**
         * @brief returns a locate context for this search context
         * @return a locate context for this search context
         */
        locate_context_t locate_phase() const requires (query_support == LOCATE)
        {
            return locate_context_t(*this);
        }

        /**
         * @brief resets the query context to an empty P
         * @param idx the index to query
         */
        inline void reset(const move_rb<support, sym_t, pos_t>& idx)
        {
            dir_lst = NO_DIR;
            m = 0;

            b = 0;
            e = idx.idx_fwd.input_size();
            b_R = 0;
            e_R = e;
            b_ = 0;
            e_ = idx.idx_fwd.M_LF().num_intervals() - 1;
            b_R_ = 0;
            e_R_ = idx.idx_bwd.M_LF().num_intervals() - 1;

            if constexpr (query_support == LOCATE) {
                s_b = { .t = RUN_START, .o = 0, .i = 0, .x = 0 };
                s_e = { .t = RUN_END, .o = 0, .i = 0, .x = e_ };
                s_b_R = { .t = RUN_START, .o = 0, .i = 0, .x = 0 };
                s_e_R = { .t = RUN_END, .o = 0, .i = 0, .x = e_R_ };
            }
        }

        /**
         * @brief returns the length of the currently matched P
         * @return length of the currently matched P
         */
        inline pos_t length() const
        {
            return m;
        }

        /**
         * @brief returns whether the context is valid (not empty)
         * @return whether the context is valid (not empty)
         */
        inline bool is_valid() const
        {
            return e >= b;
        }
        
        /**
         * @brief returns the overall number of occurrences of the currently matched P
         * @return overall number of occurrences
         */
        inline pos_t num_occ() const
        {
            return e - b + 1;
        }

        /**
         * @brief returns the SA-interval in T/T^R of the currently matched P/P^R
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
         * @brief returns the SA-interval in T of the currently matched P
         * @return suffix array interval in the forward text
         */
        inline std::pair<pos_t, pos_t> forward_sa_interval() const
        {
            return std::make_pair(b, e);
        }

        /**
         * @brief returns the SA-interval in T^R of the currently matched P^R
         * @return suffix array interval in the backward text
         */
        inline std::pair<pos_t, pos_t> backward_sa_interval() const
        {
            return std::make_pair(b_R, e_R);
        }

        /**
         * @brief prepends sym to the currently matched P; if symP occurs in the input, true is
         * returned and the query context is adjusted to store the information for symP; else,
         * false is returned and the query context is not modified
         * @param idx the index to query
         * @param sym the symbol to prepend to P
         * @return whether symP occurs in the input
         */
        inline bool prepend(const move_rb<support, sym_t, pos_t>& idx, sym_t sym)
        {
            return extend<LEFT>(idx, sym);
        }

        /**
         * @brief appends sym to the currently matched P; if Psym occurs in the input, true is
         * returned and the query context is adjusted to store the information for Psym; else,
         * false is returned and the query context is not modified
         * @param idx the index to query
         * @param sym the symbol to append to P
         * @return whether Psym occurs in the input
         */
        inline bool append(const move_rb<support, sym_t, pos_t>& idx, sym_t sym)
        {
            return extend<RIGHT>(idx, sym);
        }

        /**
         * @brief extends the currently matched P with sym; if the extended P occurs in the input, true is
         * returned and the query context is adjusted to store the information for the extended P; else,
         * false is returned and the query context is not modified
         * @tparam dir extend direction
         * @param idx the index to query
         * @param sym the symbol to extend the pattern with
         * @return whether the extended P occurs in the input
         */
        template <direction dir>
        inline bool extend(const move_rb<support, sym_t, pos_t>& idx, sym_t sym)
        {
            i_sym_t i_sym = idx.idx_fwd.map_symbol(sym);

            // If i_sym does not occur in L', then P[i..m] does not occur in T
            if (i_sym == 0) [[unlikely]] {
                return false;
            }

            search_context_t ctx_old = *this;
            bool result;

            if constexpr (dir == LEFT) {
                result = extend<LEFT>(
                    idx, i_sym,
                    b, e, b_R, e_R,
                    b_, e_,
                    s_b, s_e,
                    s_b_R, s_e_R
                );
            } else {
                result = extend<RIGHT>(
                    idx, i_sym,
                    b_R, e_R, b, e,
                    b_R_, e_R_,
                    s_b_R, s_e_R,
                    s_b, s_e
                );
            }

            if (result) [[likely]] {
                m++;
                dir_lst = dir;
            } else {
                *this = ctx_old;
            }

            return result;
        }

        /**
         * @brief prepares the context to be prepended
         * @param idx the index to query
         * @param ext_ctx 
         * @return extend context prepared for prepending the search context with every possible character in O(sigma) time
         */
        extend_context_t prepare_prepend(const move_rb<support, sym_t, pos_t>& idx) {
            return prepare_extend<LEFT>(idx);
        }

        /**
         * @brief prepares the context to be appended
         * @param idx the index to query
         * @param ext_ctx 
         * @return extend context prepared for appending the search context with every possible character in O(sigma) time
         */
        extend_context_t prepare_append(const move_rb<support, sym_t, pos_t>& idx) {
            return prepare_extend<RIGHT>(idx);
        }

        /**
         * @brief prepares the context to be extended
         * @tparam dir extend direction
         * @param idx the index to query
         * @param ext_ctx 
         * @return extend context prepared for extending the search context with every possible character in O(sigma) time
         */
        template <direction dir>
        extend_context_t prepare_extend(const move_rb<support, sym_t, pos_t>& idx) {
            extend_context_t ext_ctx;

            if constexpr (dir == LEFT) {
                prepare_extend<LEFT>(idx, ext_ctx, b, e, b_, e_, s_b, s_e);
            } else {
                prepare_extend<RIGHT>(idx, ext_ctx, b_R, e_R, b_R_, e_R_, s_b_R, s_e_R);
            }

            return ext_ctx;
        }

        /**
         * @brief prepends the context with the next character
         * @param idx the index to query
         * @param ext_ctx extend context to prepare
         */
        search_context_t prepend_next(const move_rb<support, sym_t, pos_t>& idx, extend_context_t& ext_ctx) {
            return extend_next<LEFT>(idx, ext_ctx);
        }
        
        /**
         * @brief appends the context with the next character
         * @param idx the index to query
         * @param ext_ctx extend context to prepare
         */
        search_context_t append_next(const move_rb<support, sym_t, pos_t>& idx, extend_context_t& ext_ctx) {
            return extend_next<RIGHT>(idx, ext_ctx);
        }

    protected:
        /**
         * @brief builds the array prev[0..max_sym] and next[0..max_sym],
         *        where prev[c] = select_c(L', rank_c(L', e_)) and
         *        next[c] = select_c(L', rank_c(L', b_ + 1) + 1)
         * @tparam dir extend direction
         * @param idx the index to query
         * @param prev output prev array
         * @param next output next array
         * @param max_sym maximum symbol to consider
         * @param b_ index of the input interval in M_LF of T(dir) containing b
         * @param e_ index of the input interval in M_LF of T(dir) containing e
         */
        template<direction dir>
        inline void build_prev_next(
            const move_rb<support, sym_t, pos_t>& idx,
            std::array<int64_t, 256>& prev, std::array<int64_t, 256>& next,
            i_sym_t max_sym, pos_t b_, pos_t e_
        ) const {
            const auto& idx_dir = idx.index<dir>();
            pos_t blk_size = idx_dir.L_block_size();
            pos_t max = max_sym;

            std::fill_n(next.begin(), max_sym + 1, std::numeric_limits<int64_t>::max());
            std::fill_n(prev.begin(), max_sym + 1, std::numeric_limits<int64_t>::min());

            // TODO: add optimized algorithm for small [b_, e_] intervals

            pos_t blk = div_ceil<pos_t>(b_, blk_size);
            int64_t beg = b_;
            int64_t end = std::min<pos_t>(blk * blk_size, e_);

            if (end != e_) [[likely]] {
                pos_t blk_beg = blk * idx.sigma;

                for (pos_t i = 0; i <= max; i++) {
                    next[i] = idx_dir.L_next(blk_beg + i);
                }
            }

            for (int64_t i = end; i >= beg; i--) {
                next[idx_dir.L_(i)] = i;
            }

            blk = e_ / blk_size;
            beg = std::max<pos_t>(blk * blk_size, b_);
            end = e_;

            if (beg != b_) [[likely]] {
                pos_t blk_beg = blk * idx.sigma;

                for (pos_t i = 0; i <= max; i++) {
                    prev[i] = idx_dir.L_prev(blk_beg + i);
                }
            }

            for (int64_t i = beg; i <= end; i++) {
                prev[idx_dir.L_(i)] = i;
            }
        }

        /**
         * @brief update the input interval indexes b_ and e_, and adjusts dir-SA-interval samples
         * @tparam dir extend direction
         * @param idx the index to query
         * @param b Left interval limit of the suffix array interval of P(dir) in T(dir)
         * @param e Right interval limit of the suffix array interval of P(dir) in T(dir)
         * @param b_r Left interval limit of the suffix array interval of P(!dir) in T(!dir)
         * @param e_r Right interval limit of the suffix array interval of P(!dir) in T(!dir)
         * @param b_ index of the input interval in M_LF of T(dir) containing b
         * @param e_ index of the input interval in M_LF of T(dir) containing e
         * @param s_b SA-sample information at the beginning of the dir-SA-interval
         * @param s_e SA-sample information at the end of the dir-SA-interval
         */
        template <direction dir>
        inline void update_input_intervals_and_samples(
            const move_rb<support, sym_t, pos_t>& idx,
            const pos_t& b, const pos_t& e,
            pos_t& b_, pos_t& e_,
            sample_info_t& s_b, sample_info_t& s_e
        ) const {
            const auto& idx_dir = idx.index<dir>();

            if (dir_lst != NO_DIR && dir != dir_lst) {
                if (b_ >= idx_dir.M_LF().num_intervals() || !(idx_dir.M_LF().p(b_) <= b && b < idx_dir.M_LF().p(b_ + 1))) {
                    b_ = input_interval(b, idx.S_MLF_p<dir>(), [&](pos_t x){return idx_dir.M_LF().p(x);});
                }

                if constexpr (query_support == LOCATE) {
                    if (b == idx_dir.M_LF().p(b_)) {
                        s_b = { .t = RUN_START, .o = 0, .i = 0, .x = b_ };
                    } else if (b == idx_dir.M_LF().p(b_ + 1) - 1) {
                        s_b = { .t = RUN_END, .o = 0, .i = 0, .x = b_ };
                    }
                }

                if (e_ >= idx_dir.M_LF().num_intervals() || !(idx_dir.M_LF().p(e_) <= e && e < idx_dir.M_LF().p(e_ + 1))) {
                    e_ = input_interval(e, idx.S_MLF_p<dir>(), [&](pos_t x){return idx_dir.M_LF().p(x);});
                }

                if constexpr (query_support == LOCATE) {
                    if (e == idx_dir.M_LF().p(e_)) {
                        s_e = { .t = RUN_START, .o = 0, .i = 0, .x = e_ };
                    } else if (e == idx_dir.M_LF().p(e_ + 1) - 1) {
                        s_e = { .t = RUN_END, .o = 0, .i = 0, .x = e_ };
                    }
                }
            }
        }

        /**
         * @brief extends the currently matched P with sym; if the extended P occurs in the input, true is
         * returned and the query context is adjusted to store the context for the extended P; else,
         * false is returned and the query context is not modified; let S be a string, then S(LEFT) = S and S(RIGHT) = S^R
         * and !LEFT = RIGHT = !!RIGHT
         * @tparam dir extend direction
         * @param idx the index to query
         * @param i_sym next (internal) symbol to match
         * @param b Left interval limit of the suffix array interval of P(dir) in T(dir)
         * @param e Right interval limit of the suffix array interval of P(dir) in T(dir)
         * @param b_r Left interval limit of the suffix array interval of P(!dir) in T(!dir)
         * @param e_r Right interval limit of the suffix array interval of P(!dir) in T(!dir)
         * @param b_ index of the input interval in M_LF of T(dir) containing b
         * @param e_ index of the input interval in M_LF of T(dir) containing e
         * @param s_b SA-sample information at the beginning of the dir-SA-interval
         * @param s_e SA-sample information at the end of the dir-SA-interval
         * @param s_b_R SA-sample information at the beginning of the !dir-SA-interval
         * @param s_e_R SA-sample information at the end of the !dir-SA-interval
         * @return whether the extended P occurs in the input
         */
        template <direction dir>
        bool extend(
            const move_rb<support, sym_t, pos_t>& idx,
            i_sym_t i_sym,
            pos_t& b, pos_t& e, pos_t& b_R, pos_t& e_R,
            pos_t& b_, pos_t& e_,
            sample_info_t& s_b, sample_info_t& s_e,
            sample_info_t& s_b_R, sample_info_t& s_e_R
        ) const {
            const auto& idx_dir = idx.index<dir>();

            update_input_intervals_and_samples<dir>(idx, b, e, b_, e_, s_b, s_e);

            std::array<int64_t, 256> prev;
            std::array<int64_t, 256> next;
            build_prev_next<dir>(idx, prev, next, i_sym, b_, e_);
            
            pos_t b_old = b;
            pos_t e_old = e;

            pos_t b__old = b_;
            pos_t e__old = e_;

            pos_t b_R_old = b_R;
            pos_t e_R_old = e_R;

            if (next[i_sym] > prev[i_sym] || prev[i_sym] >= idx_dir.M_LF().num_intervals()) [[unlikely]] {
                return false;
            }

            b_ = next[i_sym];
            e_ = prev[i_sym];

            if (b_ != b__old) {
                b = idx_dir.M_LF().p(b_);
            }

            if (e_ != e__old) {
                e = idx_dir.M_LF().p(e_ + 1) - 1;
            }

            if constexpr (query_support == LOCATE) {
                if (b_ != b__old) {
                    s_b = { .t = RUN_START, .o = 0, .i = 1, .x = b_ };
                } else if (s_b.t == RUN_START) {
                    s_b.i++;
                } else if (s_b.t == RUN_END) {
                    pos_t p_1 = b + s_b.o;
                    pos_t y_1 = idx_dir.M_LF().p(b_ + 1);

                    if (p_1 < y_1) {
                        s_b.i++;
                    } else {
                        s_b.o = y_1 - b - 1;
                        s_b.x = b_;
                        s_b.i = 1;
                    }
                }

                if (e_ != e__old) {
                    s_e = { .t = RUN_END, .o = 0, .i = 1, .x = e_ };
                } else if (s_e.t == RUN_END) {
                    s_e.i++;
                } else if (s_e.t == RUN_START) {
                    pos_t p_2 = e - s_e.o;
                    pos_t y_2 = idx_dir.M_LF().p(e_);

                    if (y_2 <= p_2) {
                        s_e.i++;
                    } else {
                        s_e.o = e - y_2;
                        s_e.x = e_;
                        s_e.i = 1;
                    }
                }
            }

            pos_t b_prime = b;
            pos_t e_prime = e;

            pos_t b__prime = b_;
            pos_t e__prime = e_;

            for (i_sym_t i = 0; i < i_sym; i++) {
                if (next[i] <= prev[i] && prev[i] < idx_dir.M_LF().num_intervals()) {
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

            if constexpr (query_support == LOCATE) {
                if (e - b != e_prime - b_prime) {
                    if (b_prime == b_old) {
                        s_b = {
                            .t = RUN_END,
                            .o = idx_dir.M_LF().p(b__prime + 1) - 1 - b_prime,
                            .i = 1,
                            .x = b__prime
                        };
                    }

                    if (e_prime == e_old) {
                        s_e = {
                            .t = RUN_START,
                            .o = e_prime - idx_dir.M_LF().p(e__prime),
                            .i = 1,
                            .x = e__prime
                        };
                    }
                }

                if (e - b != e_old - b_old) {
                    if (s_b_R.t != NO_SAMPLE) {
                        pos_t p_1 = b_R_old + s_b_R.o;

                        if (!(b_R <= p_1 && p_1 <= e_R)) {
                            s_b_R.t = NO_SAMPLE;
                        } else {
                            s_b_R.o -= b_R - b_R_old;
                        }
                    }

                    if (s_e_R.t != NO_SAMPLE) {
                        pos_t p_2 = e_R_old - s_e_R.o;

                        if (!(b_R <= p_2 && p_2 <= e_R)) {
                            s_e_R.t = NO_SAMPLE;
                        } else {
                            s_e_R.o -= e_R_old - e_R;
                        }
                    }
                }

                assert(s_b.t || s_e.t || s_b_R.t || s_e_R.t);
            }

            return true;
        }

        /**
         * @brief prepares the context to be extended
         * @tparam dir extend direction
         * @param idx the index to query
         * @param ext_ctx extend context to prepare
         * @param b Left interval limit of the suffix array interval of P(dir) in T(dir)
         * @param e Right interval limit of the suffix array interval of P(dir) in T(dir)
         * @param b_ index of the input interval in M_LF of T(dir) containing b
         * @param e_ index of the input interval in M_LF of T(dir) containing e
         * @param s_b SA-sample information at the beginning of the dir-SA-interval
         * @param s_e SA-sample information at the end of the dir-SA-interval
         */
        template <direction dir>
        void prepare_extend(
            const move_rb<support, sym_t, pos_t>& idx,
            extend_context_t& ext_ctx,
            const pos_t b, const pos_t e, pos_t& b_, pos_t& e_,
            sample_info_t& s_b, sample_info_t& s_e
        ) {
            static constexpr direction rev_dir = flip<dir>();

            update_input_intervals_and_samples<dir>(idx, b, e, b_, e_, s_b, s_e);
            build_prev_next<dir>(idx, ext_ctx.prev, ext_ctx.next, idx.sigma, b_, e_);

            ext_ctx.b_R_nxt = beg<rev_dir>() + (ext_ctx.next[0] <= ext_ctx.prev[0] &&
                                                ext_ctx.prev[0] != idx.index<dir>().M_LF().num_intervals());
            
            ext_ctx.sym_nxt = 0;
            ext_ctx.template next_symbol<dir>(idx);
        }

    public:
        /**
         * @brief extends the context with the next character
         * @tparam dir extend direction
         * @param idx the index to query
         * @param ext_ctx extend context to prepare
         */
        template <direction dir>
        search_context_t extend_next(const move_rb<support, sym_t, pos_t>& idx, extend_context_t& ext_ctx) {
            const auto& idx_dir = idx.index<dir>();
            static constexpr direction rev_dir = flip<dir>();

            pos_t b_old, e_old;

            if constexpr (dir == LEFT) {
                b_old = b;
                e_old = e;
            } else {
                b_old = b_R;
                e_old = e_R;
            }

            search_context_t ctx(idx.unmap_symbol(ext_ctx.sym_nxt));

            ctx.beg_run<dir>() = ext_ctx.next[ext_ctx.sym_nxt];
            ctx.end_run<dir>() = ext_ctx.prev[ext_ctx.sym_nxt];

            if (idx_dir.L_(beg_run<dir>()) == ext_ctx.sym_nxt) [[unlikely]] {
                ctx.beg<dir>() = beg<dir>();
                
                if constexpr (query_support == LOCATE) {
                    if (smpl_beg<dir>().t == RUN_START) {
                        ctx.smpl_beg<dir>() = smpl_beg<dir>();
                        ctx.smpl_beg<dir>().i++;
                    } else if (smpl_beg<dir>().t == RUN_END) {
                        pos_t p_1 = beg<dir>() + smpl_beg<dir>().o;
                        pos_t y_1 = idx_dir.M_LF().p(beg_run<dir>() + 1);

                        if (p_1 < y_1) {
                            ctx.smpl_beg<dir>() = smpl_beg<dir>();
                            ctx.smpl_beg<dir>().i++;
                        } else {
                            ctx.smpl_beg<dir>() = {
                                .t = RUN_END,
                                .o = y_1 - beg<dir>() - 1,
                                .i = 1,
                                .x = beg_run<dir>()
                            };
                        }
                    } else {
                        ctx.smpl_beg<dir>().t = NO_SAMPLE;
                    }
                }
            } else {
                ctx.beg<dir>() = idx_dir.M_LF().p(ctx.beg_run<dir>());
                if constexpr (query_support == LOCATE) ctx.smpl_beg<dir>() = { .t = RUN_START, .o = 0, .i = 1, .x = ctx.beg_run<dir>() };
            }

            if (idx_dir.L_(end_run<dir>()) == ext_ctx.sym_nxt) [[unlikely]] {
                ctx.end<dir>() = end<dir>();
                
                if constexpr (query_support == LOCATE) {
                    if (smpl_end<dir>().t == RUN_END) {
                        ctx.smpl_end<dir>() = smpl_end<dir>();
                        ctx.smpl_end<dir>().i++;
                    } else if (smpl_end<dir>().t == RUN_START) {
                        pos_t p_2 = end<dir>() - smpl_end<dir>().o;
                        pos_t y_2 = idx_dir.M_LF().p(end_run<dir>());

                        if (y_2 <= p_2) {
                            ctx.smpl_end<dir>() = smpl_end<dir>();
                            ctx.smpl_end<dir>().i++;
                        } else {
                            ctx.smpl_end<dir>() = {
                                .t = RUN_START,
                                .o = end<dir>() - y_2,
                                .i = 1,
                                .x = end_run<dir>()
                            };
                        }
                    } else {
                        ctx.smpl_end<dir>().t = NO_SAMPLE;
                    }
                }
            } else {
                ctx.end<dir>() = idx_dir.M_LF().p(ctx.end_run<dir>() + 1) - 1;
                if constexpr (query_support == LOCATE) ctx.smpl_end<dir>() = { .t = RUN_END, .o = 0, .i = 1, .x = ctx.end_run<dir>() };
            }
            
            pos_t b_prime = ctx.beg<dir>();
            pos_t e_prime = ctx.end<dir>();
            
            pos_t b__prime = ctx.beg_run<dir>();
            pos_t e__prime = ctx.end_run<dir>();

            if (ctx.beg_run<dir>() == ctx.end_run<dir>()) {
                if (ctx.beg<dir>() == ctx.end<dir>()) {
                    idx_dir.M_LF().move(ctx.beg<dir>(), ctx.beg_run<dir>());
                    ctx.end<dir>() = ctx.beg<dir>();
                    ctx.end_run<dir>() = ctx.beg_run<dir>();
                } else {
                    pos_t diff_eb = ctx.end<dir>() - ctx.beg<dir>();
                    idx_dir.M_LF().move(ctx.beg<dir>(), ctx.beg_run<dir>());
                    ctx.end<dir>() = ctx.beg<dir>() + diff_eb;
                    ctx.end_run<dir>() = ctx.beg_run<dir>();

                    while (ctx.end<dir>() >= idx_dir.M_LF().p(ctx.end_run<dir>() + 1)) {
                        ctx.end_run<dir>()++;
                    }
                }
            } else {
                idx_dir.M_LF().move(ctx.beg<dir>(), ctx.beg_run<dir>());
                idx_dir.M_LF().move(ctx.end<dir>(), ctx.end_run<dir>());
            }

            ctx.beg<rev_dir>() = ext_ctx.b_R_nxt;
            ctx.end<rev_dir>() = ext_ctx.b_R_nxt + (ctx.end<dir>() - ctx.beg<dir>());

            if constexpr (query_support == LOCATE) {
                if (ctx.end<dir>() - ctx.beg<dir>() != e_prime - b_prime) [[likely]] {
                    if (b_prime == b_old) {
                        ctx.smpl_beg<dir>() = {
                            .t = RUN_END,
                            .o = idx_dir.M_LF().p(b__prime + 1) - 1 - b_prime,
                            .i = 1,
                            .x = b__prime
                        };
                    }

                    if (e_prime == e_old) {
                        ctx.smpl_end<dir>() = {
                            .t = RUN_START,
                            .o = e_prime - idx_dir.M_LF().p(e__prime),
                            .i = 1,
                            .x = e__prime
                        };
                    }
                }

                if (ctx.end<dir>() - ctx.beg<dir>() != e_old - b_old) [[likely]] {
                    if (smpl_beg<rev_dir>().t != NO_SAMPLE) {
                        pos_t p_1 = beg<rev_dir>() + smpl_beg<rev_dir>().o;

                        if (!(ctx.beg<rev_dir>() <= p_1 && p_1 <= ctx.end<rev_dir>())) {
                            ctx.smpl_beg<rev_dir>().t = NO_SAMPLE;
                        } else {
                            ctx.smpl_beg<rev_dir>() = smpl_beg<rev_dir>();
                            ctx.smpl_beg<rev_dir>().o -= ctx.beg<rev_dir>() - beg<rev_dir>();
                        }
                    } else {
                        ctx.smpl_beg<rev_dir>().t = NO_SAMPLE;
                    }

                    if (smpl_end<rev_dir>().t != NO_SAMPLE) {
                        pos_t p_2 = end<rev_dir>() - smpl_end<rev_dir>().o;

                        if (!(ctx.beg<rev_dir>() <= p_2 && p_2 <= ctx.end<rev_dir>())) {
                            ctx.smpl_end<rev_dir>().t = NO_SAMPLE;
                        } else {
                            ctx.smpl_end<rev_dir>() = smpl_end<rev_dir>();
                            ctx.smpl_end<rev_dir>().o -= end<rev_dir>() - ctx.end<rev_dir>();
                        }
                    } else {
                        ctx.smpl_end<rev_dir>().t = NO_SAMPLE;
                    }
                } else {
                    ctx.smpl_beg<rev_dir>() = smpl_beg<rev_dir>();
                    ctx.smpl_end<rev_dir>() = smpl_end<rev_dir>();
                }

                assert(ctx.s_b.t || ctx.s_e.t || ctx.s_b_R.t || ctx.s_e_R.t);
            }
            
            ext_ctx.b_R_nxt = ctx.end<rev_dir>() + 1;
            ctx.dir_lst = dir;
            ctx.m = m + 1;
            ext_ctx.template next_symbol<dir>(idx);

            return ctx;
        }
    };

    /**
     * @brief structure storing the information for locating all occurrences of a search context
     */
    struct locate_context_t {
        protected:
        direction dir; // current locate direction
        pos_t occ_rem; // number of remaining occurrences to locate
        pos_t c; // initial position in the suffix array
        pos_t SA_c; // initial suffix SA[c] in the suffix array interval
        pos_t i; // current position in the suffix array interval
        pos_t SA_i; // current suffix SA[i] in the suffix array interval
        
        // index of the input inteval of M_Phi/M_Phi^{-1} containing SA_i
        std::conditional_t<support == _locate_move, pos_t, empty_t> s_; 

        struct rlzsa_ctx_t {
            pos_t x_p, x_lp, x_cp, x_r, s_p; // variables for decoding the rlzsa
        };

        std::conditional_t<support == _locate_rlzsa, rlzsa_ctx_t, empty_t> rlz_l; // rlzsa context for decoding to the left
        std::conditional_t<support == _locate_rlzsa, rlzsa_ctx_t, empty_t> rlz_r; // rlzsa context for decoding to the right

    public:
        /**
         * @brief constructs an empty locate context
         */
        locate_context_t() {}

        /**
         * @brief constructs a new locate context for a given search context
         * @param ctx a search context
         */
        locate_context_t(const search_context_t<LOCATE>& ctx)
        {
            occ_rem = ctx.num_occ();
            dir = NO_DIR;
        }

        /**
         * @brief returns the number of remaining (not yet reported) occurrences of the currently matched P
         * @return number of remaining occurrences
         */
        inline pos_t num_occ_rem() const
        {
            return occ_rem;
        }

    protected:
        /**
         * @brief computes a first occurrence (SA[i]) of P, and adjusts the locate context accordingly
         * @param idx the index to query
         * @param ctx the search context associated with this locate context
         * @return an occurrence (SA[i]) of P
         */
        inline pos_t first_occ(const move_rb<support, sym_t, pos_t>& idx, const search_context_t<LOCATE>& ctx)
        {
            if (ctx.s_b.t) {
                i = ctx.b + ctx.s_b.o;
                SA_i = (ctx.s_b.t == RUN_START ? idx.idx_fwd.SA_s(ctx.s_b.x) : idx.idx_fwd.SA_s_(ctx.s_b.x)) - ctx.s_b.i;
            } else if (ctx.s_e.t) {
                i = ctx.e - ctx.s_e.o;
                SA_i = (ctx.s_e.t == RUN_START ? idx.idx_fwd.SA_s(ctx.s_e.x) : idx.idx_fwd.SA_s_(ctx.s_e.x)) - ctx.s_e.i;
            } else if (ctx.s_b_R.t) {
                if (ctx.s_b_R.t == RUN_START) {
                    i = idx._SA_sR_m1[ctx.s_b_R.x];
                    SA_i = idx.n - idx.idx_bwd.SA_s(ctx.s_b_R.x) - (ctx.m - ctx.s_b_R.i + 1);
                } else {
                    i = idx._SA_eR_m1[ctx.s_b_R.x];
                    SA_i = idx.n - idx.idx_bwd.SA_s_(ctx.s_b_R.x) - (ctx.m - ctx.s_b_R.i + 1);
                }

                pos_t i_ = input_interval(i, idx._S_MLF_p_fwd,
                    [&](pos_t x){return idx.idx_fwd.M_LF().p(x);});

                for (pos_t j = 0; j < ctx.m - ctx.s_b_R.i; j++) {
                    idx.idx_fwd.M_LF().move(i, i_);
                }
            } else /* if (ctx.s_e_R.t) */ {
                if (ctx.s_e_R.t == RUN_START) {
                    i = idx._SA_sR_m1[ctx.s_e_R.x];
                    SA_i = idx.n - idx.idx_bwd.SA_s(ctx.s_e_R.x) - (ctx.m - ctx.s_e_R.i + 1);
                } else {
                    i = idx._SA_eR_m1[ctx.s_e_R.x];
                    SA_i = idx.n - idx.idx_bwd.SA_s_(ctx.s_e_R.x) - (ctx.m - ctx.s_e_R.i + 1);
                }

                pos_t i_ = input_interval(i, idx._S_MLF_p_fwd,
                    [&](pos_t x){return idx.idx_fwd.M_LF().p(x);});

                for (pos_t j = 0; j < ctx.m - ctx.s_e_R.i; j++) {
                    idx.idx_fwd.M_LF().move(i, i_);
                }
            }

            assert(ctx.b <= i && i <= ctx.e);

            c = i;
            SA_c = SA_i;

            if (occ_rem > 1) [[likely]] {
                dir = i > ctx.b ? LEFT : RIGHT;

                if constexpr (support == _locate_move) {
                    if (dir == LEFT) {
                        s_ = input_interval(SA_i, idx._S_MPhi_p,
                            [&](pos_t x){return idx.idx_fwd.M_Phi().p(x);});
                    } else {
                        s_ = input_interval(SA_i, idx._S_MPhi_m1_p,
                            [&](pos_t x){return idx.idx_fwd.M_Phi_m1().p(x);});
                    }
                } else if constexpr (support == _locate_rlzsa) {
                    if (dir == LEFT) {
                        idx.idx_fwd.init_rlzsa(i, rlz_l.x_p, rlz_l.x_lp, rlz_l.x_cp, rlz_l.x_r, rlz_l.s_p);

                        if (c < ctx.e) [[likely]] {
                            rlz_r = rlz_l;
                            pos_t tmp_1 = c;
                            pos_t tmp_2;
                            idx.idx_fwd.next_rlzsa(tmp_1, tmp_2, rlz_r.x_p, rlz_r.x_lp, rlz_r.x_cp, rlz_r.x_r, rlz_r.s_p);
                        }

                        idx.idx_fwd.turn_rlzsa_left(rlz_l.x_p, rlz_l.x_lp, rlz_l.x_cp, rlz_l.x_r, rlz_l.s_p);
                    } else {
                        i++;
                        idx.idx_fwd.init_rlzsa(i, rlz_r.x_p, rlz_r.x_lp, rlz_r.x_cp, rlz_r.x_r, rlz_r.s_p);
                    }
                }
            }

            occ_rem--;
            return SA_i;
        }

    public:
        /**
         * @brief reports the next occurrence of the currently matched P
         * @return next occurrence
         */
        inline pos_t next_occ(const move_rb<support, sym_t, pos_t>& idx, const search_context_t<LOCATE>& ctx)
        {
            if (dir == NO_DIR) [[unlikely]] {
                return first_occ(idx, ctx);
            }

            pos_t occ;

            if constexpr (support == _locate_move) {
                if (dir == LEFT) {
                    idx.idx_fwd.M_Phi().move(SA_i, s_);
                    occ = SA_i;
                    i--;

                    if (i == ctx.b) [[unlikely]] {
                        i = c;
                        SA_i = SA_c;
                        dir = RIGHT;

                        s_ = input_interval(SA_i, idx._S_MPhi_m1_p,
                            [&](pos_t x){return idx.idx_fwd.M_Phi_m1().p(x);});
                    }
                } else /* if (dir == RIGHT) */ {
                    idx.idx_fwd.M_Phi_m1().move(SA_i, s_);
                    occ = SA_i;
                    i++;
                }
            } else if constexpr (support == _locate_rlzsa) {
                if (dir == LEFT) {
                    idx.idx_fwd.prev_rlzsa(i, SA_i, rlz_l.x_p, rlz_l.x_lp, rlz_l.x_cp, rlz_l.x_r, rlz_l.s_p);
                    occ = SA_i;

                    if (i == ctx.b && c < ctx.e) [[unlikely]] {
                        i = c + 1;
                        dir = RIGHT;
                        SA_i = SA_c;
                    }
                } else /* if (dir == RIGHT) */ {
                    idx.idx_fwd.next_rlzsa(i, SA_i, rlz_r.x_p, rlz_r.x_lp, rlz_r.x_cp, rlz_r.x_r, rlz_r.s_p);
                    occ = SA_i;
                }
            }

            occ_rem--;
            return occ;
        }

        /**
         * @brief locates the remaining (not yet reported) occurrences of the currently matched P
         * @param Occ vector to append the occurrences to
         */
        inline void locate(
            const move_rb<support, sym_t, pos_t>& idx,
            const search_context_t<LOCATE>& ctx,
            std::vector<pos_t>& Occ
        ) {
            Occ.reserve(Occ.size() + occ_rem);

            if (occ_rem == ctx.num_occ()) {
                Occ.emplace_back(first_occ(idx, ctx));
            }

            if constexpr (support == _locate_move) {
                if (dir == LEFT) {
                    while (i > ctx.b) {
                        idx.idx_fwd.M_Phi().move(SA_i, s_);
                        Occ.emplace_back(SA_i);
                        i--;
                    }

                    if (c < ctx.e) {
                        i = c;
                        SA_i = SA_c;
                        s_ = input_interval(SA_i, idx._S_MPhi_m1_p,
                            [&](pos_t x){return idx.idx_fwd.M_Phi_m1().p(x);});
                    }
                }

                if (c < ctx.e) {
                    while (i < ctx.e) {
                        idx.idx_fwd.M_Phi_m1().move(SA_i, s_);
                        Occ.emplace_back(SA_i);
                        i++;
                    }
                }
            } else if constexpr (support == _locate_rlzsa) {
                if (dir == LEFT) {
                    idx.idx_fwd.report_rlzsa_left(i, ctx.b, SA_i,
                        rlz_l.x_p, rlz_l.x_lp, rlz_l.x_cp, rlz_l.x_r, rlz_l.s_p,
                        [&](pos_t, pos_t occ){Occ.emplace_back(occ);}
                    );

                    if (c < ctx.e) [[likely]] {
                        i = c + 1;
                        dir = RIGHT;
                        SA_i = SA_c;
                    }
                }

                if (c < ctx.e) [[likely]] {
                    idx.idx_fwd.report_rlzsa_right(i, ctx.e, SA_i,
                        rlz_r.x_p, rlz_r.x_lp, rlz_r.x_cp, rlz_r.x_r, rlz_r.s_p,
                        [&](pos_t, pos_t occ){Occ.emplace_back(occ);}
                    );
                }
            }

            occ_rem = 0;
        }

        /**
         * @brief locates the remaining (not yet reported) occurrences of the currently matched P
         * @return vector containing the occurrences
         */
        std::vector<pos_t> locate(const move_rb<support, sym_t, pos_t>& idx, const search_context_t<LOCATE>& ctx)
        {
            std::vector<pos_t> Occ;
            locate(idx, ctx, Occ);
            return Occ;
        }
    };
    
    /**
     * @brief stores the variables needed to extend a search_context with all passible characters in O(sigma) time
     */
    struct extend_context_t {
        friend class search_context_t<COUNT>;
        friend class search_context_t<LOCATE>;

        protected:
        std::array<int64_t, 256> prev; // the prev array for all extensions
        std::array<int64_t, 256> next; // the next array for all extensions

        i_sym_t sym_nxt; // next symbol to extend the context with
        pos_t b_R_nxt; // current left interval limit of the suffix array interval of P(!dir) extended with sym_nxt in T(!dir)

        public:
        /**
         * @brief returns the next symbol to extend the context with
         * @param idx the index to query
         * @return the next symbol to extend the context with
         */
        sym_t current_symbol(const move_rb<support, sym_t, pos_t>& idx) const {
            return idx.unmap_symbol(sym_nxt);
        }

        /**
         * @brief returns whether there is a symbol left to extend the context with
         * @param idx the index to query
         * @return whether there is a symbol left to extend the context with
         */
        bool can_extend(const move_rb<support, sym_t, pos_t>& idx) const {
            return sym_nxt < idx.sigma;
        }

        protected:
        /**
         * @brief advances sym_nxt to the next symbol the pattern can be extended with (or sigma, if it cannot be extended further)
         * @tparam dir direction the context should be extended with
         * @param idx the index to query
         */
        template <direction dir>
        void next_symbol(const move_rb<support, sym_t, pos_t>& idx) {
            do {
                sym_nxt++;
            } while (sym_nxt < idx.sigma && (
                     next[sym_nxt] > prev[sym_nxt] ||
                     prev[sym_nxt] == idx.index<dir>().M_LF().num_intervals()));
        }
    };

    /**
     * @brief returns a query context for the index
     * @return search_context_t
     */
    template <move_rb_query_support_t query_support>
    inline search_context_t<query_support> search() const
    {
        return search_context_t<query_support>(*this);
    }

    /**
     * @brief searches for a substring P[l,r) of a P
     * @param P the pattern to search
     * @param l left substring boundary
     * @param r right substring boundary
     * @return search context of P[l,r)
     */
    template <move_rb_query_support_t query_support>
    inline search_context_t<query_support> search_substring(const inp_t& P, pos_t l = 0, pos_t r = std::numeric_limits<pos_t>::max()) const
    {
        r = std::min<pos_t>(r, P.size() - 1);
        auto ctx = search<query_support>();

        for (int64_t i = r; i >= l; i--) {
            if (!ctx.prepend(*this, P[i])) {
                return search<query_support>();
            }
        }

        return ctx;
    }

    // ############################# APPROXIMATE PATTERN MATCHING #############################
    
    /**
     * @brief counts a pattern with at most k_max mismatches
     * @param P the pattern to search
     * @param scheme the search scheme to use (provides k_max)
     * @return number of occurrences of P in T with at most k_max mismatches
     */
    pos_t count_with_mismatches(const inp_t& P, const search_scheme_t& scheme) const
    {
        if (scheme.k_max == 0) {
            return forward_index().count(P);
        }

        auto ctxts = execute_search_scheme<COUNT>(P, scheme);
        pos_t count = 0;

        for (auto ctx: ctxts) {
            count += ctx.num_occ();
        }
        
        return count;
    }

    /**
     * @brief locaes a pattern with at most k_max mismatches
     * @param P the pattern to search
     * @param scheme the search scheme to use (provides k_max)
     * @return occurrences of P in T with at most k_max mismatches
     */
    std::vector<pos_t> locate_with_mismatches(const inp_t& P, const search_scheme_t& scheme) const
    {
        if (scheme.k_max == 0) {
            return forward_index().locate(P);
        }

        auto ctxts = execute_search_scheme<LOCATE>(P, scheme);
        pos_t count = 0;

        for (auto ctx: ctxts) {
            count += ctx.num_occ();
        }

        std::vector<pos_t> Occ;
        Occ.reserve(Occ.size() + count);

        for (auto ctx: ctxts) {
            ctx.locate_phase().locate(*this, ctx, Occ);
        }

        return Occ;
    }
    
    // type of hashset to store search contexts in
    template <move_rb_query_support_t query_support>
    using search_context_set_t = tsl::sparse_set<search_context_t<query_support>, typename search_context_t<query_support>::hash>;

    /**
     * @brief executes a given search scheme for a pattern
     * @param scheme the search scheme to use
     * @return all search contexts of P in T with at most k_max mismatches (assuming the schearch scheme covers all error configurations)
     */
    template <move_rb_query_support_t query_support>
    search_context_set_t<query_support> execute_search_scheme(const inp_t& P, const search_scheme_t& scheme) const
    {
        search_context_set_t<query_support> ctxts;
        pos_t m = P.size();
        if (m == 0) return ctxts;
        uint8_t k_max = scheme.k_max;
        const distance_metric_t& dist_metr = scheme.dist_metr;

        if (k_max == 0) [[unlikely]] {
            auto ctx = search_substring<query_support>(P, 0, m - 1);
            if (ctx.is_valid()) ctxts.emplace(ctx);
            return ctxts;
        }
        
        uint8_t parts = scheme.parts;
        assert(parts <= m);

        using search_state_t = std::tuple<
            search_context_t<query_support>, // ctx
            uint8_t, // k
            uint8_t, // p_idx
            pos_t // pos
        >;

        std::vector<search_state_t> nav_stack;
        auto init_ctx = search<query_support>();
        pos_t l = m / parts;

        for (uint8_t s_idx = 0; s_idx < scheme.searches.size(); s_idx++) {
            const search_t& search = scheme.searches[s_idx];

            pos_t init_beg = search[0].part * l;
            pos_t init_end = (search[0].part + 1) == parts ? m : ((search[0].part + 1) * l);
            direction init_dir = search[1].part < search[0].part ? LEFT : RIGHT;
            pos_t init_pos = init_dir == LEFT ? init_end - 1 : init_beg;

            nav_stack.emplace_back(search_state_t{init_ctx, 0, 0, init_pos});
            
            while (!nav_stack.empty()) {
                auto [ctx, k, p_idx, pos] = nav_stack.back();
                nav_stack.pop_back();
                
                pos_t part = search[p_idx].part;
                pos_t beg = part * l;
                pos_t end = (part + 1) == parts ? m : ((part + 1) * l);
                direction dir = p_idx == 0 ? init_dir : (part < search[p_idx - 1].part ? LEFT : RIGHT);

                uint8_t p_idx_nxt = p_idx;
                pos_t pos_nxt = pos + (dir == LEFT ? -1 : 1);
                pos_t beg_nxt = beg;
                pos_t end_nxt = end;
                direction dir_nxt = dir;
                pos_t ext_rem;

                if (!(beg <= pos_nxt && pos_nxt < end)) [[unlikely]] {
                    p_idx_nxt = p_idx + 1;

                    if (p_idx_nxt < parts) {
                        pos_t part_nxt = search[p_idx_nxt].part;
                        dir_nxt = part_nxt < part ? LEFT : RIGHT;
                        beg_nxt = part_nxt * l;
                        end_nxt = (part_nxt + 1) == parts ? m : ((part_nxt + 1) * l);
                        pos_nxt = dir_nxt == LEFT ? end_nxt - 1 : beg_nxt;
                        ext_rem = end_nxt - beg_nxt;
                    } else {
                        ext_rem = std::numeric_limits<pos_t>::max();
                    }
                } else {
                    ext_rem = dir == LEFT ? pos - beg : (end - pos - 1);
                }
                
                if (k < k_max && k < search[p_idx].k_max) {
                    auto ext_ctx = dir == LEFT ? ctx.prepare_prepend(*this) : ctx.prepare_append(*this);

                    while (ext_ctx.can_extend(*this)) {
                        auto ctx_nxt = dir == LEFT ? ctx.prepend_next(*this, ext_ctx) : ctx.append_next(*this, ext_ctx);
                        bool is_mismatch = ctx_nxt.last_added_symbol() != P[pos];
                        uint8_t k_nxt = k + is_mismatch;

                        if (p_idx_nxt == parts) {
                            ctxts.emplace(ctx_nxt);
                        } else if (search[p_idx_nxt].k_min <= k_nxt + ext_rem) {
                            nav_stack.emplace_back(search_state_t{ctx_nxt, k_nxt, p_idx_nxt, pos_nxt});
                        }
                    }
                } else {
                    uint8_t k_nxt = k;

                    if ((p_idx_nxt == parts || search[p_idx_nxt].k_min <= k_nxt + ext_rem) &&
                        (dir == LEFT ? ctx.prepend(*this, P[pos]) : ctx.append(*this, P[pos])))
                    {
                        if (p_idx_nxt == parts) {
                            ctxts.emplace(ctx);
                        } else {
                            nav_stack.emplace_back(search_state_t{ctx, k_nxt, p_idx_nxt, pos_nxt});
                        }
                    }
                }
            }
        }

        return ctxts;
    }

    // ############################# SERIALIZATION METHODS #############################

    /**
     * @brief stores the index to an output stream
     * @param out output stream to store the index to
     */
    void serialize(std::ostream& out) const
    {
        bool is_64_bit = std::is_same_v<pos_t, uint64_t>;
        out.write((char*) &is_64_bit, 1);
        move_r_support _support = support;
        out.write((char*) &_support, sizeof(move_r_support));

        out.write((char*) &n, sizeof(pos_t));
        out.write((char*) &sigma, sizeof(pos_t));

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
        bool is_64_bit;
        in.read((char*) &is_64_bit, 1);
        
        if (is_64_bit != std::is_same_v<pos_t, uint64_t>) {
            std::cout << "error: cannot load a" << (is_64_bit ? "64" : "32") << "-bit"
                      << " index into a " << (is_64_bit ? "32" : "64") << "-bit index-object" << std::flush;
            return;
        }

        move_r_support _support;
        in.read((char*) &_support, sizeof(move_r_support));

        in.read((char*) &n, sizeof(pos_t));
        in.read((char*) &sigma, sizeof(pos_t));

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