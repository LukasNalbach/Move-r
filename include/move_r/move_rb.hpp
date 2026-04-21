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
#include <hash_set2.hpp>
#include <hash_set3.hpp>
#include <hash_set4.hpp>
#include <hash_set8.hpp>
#include <tsl/sparse_set.h>

#include <misc/hash.hpp>
#include <misc/strings.hpp>
#include <move_r/move_r.hpp>
#include <misc/search_schemes.hpp>
#include <misc/sub_string.hpp>
#include <misc/edit_distance_matrix.hpp>

enum move_rb_query_support_t : uint8_t {
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

    template <direction_t dir>
    const std::conditional_t<dir == LEFT, move_r_fwd_t, move_r_bwd_t>& index() const
    {
        if constexpr (dir == LEFT) {
            return idx_fwd;
        } else {
            return idx_bwd;
        }
    }

    template <direction_t dir>
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
            if (std::is_same_v<pos_t, uint32_t> && n <= std::numeric_limits<int32_t>::max()) {
                build<false, int32_t>(input, params);
            } else {
                build<false, int64_t>(input, params);
            }
        } else {
            build<true, int32_t>(input, params);
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

    protected:

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

    public:
    /**
     * @brief stores the variables needed to perform bidirectional pattern search and locate queries
     */
    template <move_rb_query_support_t query_support>
    struct search_context_t {
        friend class move_rb;
        friend class locate_context_t;
        template <move_rb_query_support_t _query_support> friend class search_context_t;

    protected:
        // ############################# VARIABLES FOR THE SEARCH PHASE #############################
        
        direction_t dir_lst; // last performed pattern extend direction
        sym_t sym_lst; // last-added symbol
        pos_t err; // number of errors (only used for approximate pattern matching output)
        pos_t m; // length of the currently matched pattern P
        pos_t b, e; // [b, e] = SA-interval in T of the currently matched P
        pos_t b_R, e_R; // [b_R, e_R] = SA-interval in T^R of the reverse of the currently matched P
        pos_t b_, e_; // indices of the input intervals in M_LF containing b and e
        pos_t b_R_, e_R_; // indices of the input intervals in M_LF^R containing b_R and e_R

        // ############################# VARIABLES FOR MAINTAINING SA-SAMPLE INFORMATION DURING SEARCH #############################

        using sample_info_t = std::conditional_t<query_support == LOCATE, sample_info_pack_t, std::monostate>;
        using shift_t = std::conditional_t<query_support == LOCATE, pos_t, std::monostate>;
        using dpth_t = std::conditional_t<query_support == LOCATE, pos_t, std::monostate>;
        using rprtd_t = std::conditional_t<query_support == LOCATE, bool, std::monostate>;

        shift_t shft; // right shift (in the text) of the occurrences
        dpth_t dpth; // depth (additional length) of this context in the banded alignment matrix
        rprtd_t rprtd; // true <=> this context has already been reported
        
        sample_info_t s_b, s_e; // sample info for the beginning/end of the SA-interval in T of the currently matched P
        sample_info_t s_b_R, s_e_R; // sample info for the beginning/end of the SA-interval in T^R of the reverse of the currently matched P

        template <direction_t dir, typename this_t>
        static auto vars_impl(this_t&& t)
        {
            if constexpr (dir == LEFT) return std::forward_as_tuple(t.b, t.e, t.b_R, t.e_R, t.b_, t.e_, t.b_R_, t.e_R_, t.s_b, t.s_e, t.s_b_R, t.s_e_R);
            else                       return std::forward_as_tuple(t.b_R, t.e_R, t.b, t.e, t.b_R_, t.e_R_, t.b_, t.e_, t.s_b_R, t.s_e_R, t.s_b, t.s_e);
        }

        template <direction_t dir> auto vars() {return vars_impl<dir>(*this);}
        template <direction_t dir> auto const_vars() const {return vars_impl<dir>(*this);}

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
         * @return whether the other context other is equal to this search context; it is considered equal if it
         *         represents the same pattern (SA-interval and pattern length); other context variables may vary
         */
        bool operator==(const search_context_t& other) const
        {
            return b == other.b && e == other.e;
        }

        /**
         * @brief "less than" operator on search contexts
         * @param other the other search context to compare with
         * @return whether this context is considered "less than" the other
         */
        bool operator<(const search_context_t& other) const
        {
            return b != other.b ? b < other.b :
                  (e != other.e ? e > other.e :
                                  m < other.m);
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
            inline static pos_t operator()(const search_context_t& ctx)
            {
                auto [b, e] = ctx.forward_sa_interval();
                pos_t hash = pos_hash<pos_t>(b);
                hash_combine<pos_t>(hash, pos_hash<pos_t>(e));
                return hash;
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

        inline direction_t last_direction() const
        {
            return dir_lst;
        }

        /**
         * @brief returns a locate context for this search context
         * @return a locate context for this search context
         */
        locate_context_t locate_phase() const requires(query_support == LOCATE)
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
            sym_lst = 0;
            m = 0;
            err = 0;

            b = 0;
            e = idx.idx_fwd.input_size();
            b_R = 0;
            e_R = e;
            b_ = 0;
            e_ = idx.idx_fwd.M_LF().num_intervals() - 1;
            b_R_ = 0;
            e_R_ = idx.idx_bwd.M_LF().num_intervals() - 1;

            if constexpr (query_support == LOCATE) {
                shft = 0;
                dpth = 0;
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
         * @brief returns the depth (additional length) of this context in the banded alignment matrix
         * @return depth
         */
        inline pos_t depth() const requires(query_support == LOCATE)
        {
            return dpth;
        }

        /**
         * @brief sets the depth (additional length) of this context in the banded alignment matrix
         * @param depth depth
         */
        inline void set_depth(pos_t depth) requires(query_support == LOCATE)
        {
            this->dpth = depth;
        }

        /**
         * @brief returns true <=> this context has already been reported
         * @return whether this context has already been reported
         */
        inline bool reported() const requires(query_support == LOCATE)
        {
            return rprtd;
        }

        /**
         * @brief sets whether this context has already been reported
         * @param reported whether this context has already been reported
         */
        inline void set_reported(bool reported) requires(query_support == LOCATE)
        {
            this->rprtd = reported;
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
         * @brief returns the number of errors (only used for approximate pattern matching)
         * @return number of errors
         */
        inline pos_t errors() const
        {
            return err;
        }
        
        /**
         * @brief sets the number of errors to err (only used for approximate pattern matching)
         * @param err number of errors
         */
        inline void set_errors(pos_t err)
        {
            this->err = err;
        }
        
        /**
         * @brief sets the number of errors to err (only used for approximate pattern matching)
         * @param err number of errors
         */
        inline void set_shift(pos_t shft) requires(query_support == LOCATE)
        {
            this->shft = shft;
        }
        
        inline pos_t shift() const requires(query_support == LOCATE)
        {
            return shft;
        }

        /**
         * @brief returns the SA-interval in T/T^R of the currently matched P/P^R
         * @tparam dir text direction_t
         * @return suffix array interval in the forward text
         */
        template <direction_t dir>
        inline std::tuple<pos_t, pos_t> sa_interval() const
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
        inline std::tuple<pos_t, pos_t> forward_sa_interval() const
        {
            return std::tuple<pos_t, pos_t>{b, e};
        }

        /**
         * @brief returns the SA-interval in T^R of the currently matched P^R
         * @return suffix array interval in the backward text
         */
        inline std::tuple<pos_t, pos_t> backward_sa_interval() const
        {
            return std::tuple<pos_t, pos_t>{b_R, e_R};
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
        template<direction_t dir>
        inline void build_prev_next(
            const move_rb<support, sym_t, pos_t>& idx,
            std::array<int64_t, 256>& prev, std::array<int64_t, 256>& next,
            i_sym_t max_sym
        ) const {
            const auto& idx_dir = idx.index<dir>();
            pos_t blk_size = idx_dir.L_block_size();
            pos_t max = max_sym;

            std::fill_n(next.begin(), max_sym + 1, std::numeric_limits<int64_t>::max());
            std::fill_n(prev.begin(), max_sym + 1, std::numeric_limits<int64_t>::min());

            const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s_b, s_e, s_b_R, s_e_R] = const_vars<dir>();

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
        template <direction_t dir>
        inline void update_input_intervals_and_samples(const move_rb<support, sym_t, pos_t>& idx)
        {
            const auto& idx_dir = idx.index<dir>();
            auto&& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s_b, s_e, s_b_R, s_e_R] = vars<dir>();

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
        
        template <direction_t dir>
        inline void update_samples_stage_1(const move_rb<support, sym_t, pos_t>& idx, search_context_t& ctx) const
        {
            const auto& idx_dir = idx.index<dir>();
            const auto& [b_old, e_old, b_R_old, e_R_old,
                         b__old, e__old, b_R__old, e_R__old,
                         s_b_old, s_e_old, s_b_R_old, s_e_R_old] = const_vars<dir>();
            auto&& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s_b, s_e, s_b_R, s_e_R] = ctx.vars<dir>();

            if (b_ != b__old) [[likely]] {
                s_b = { .t = RUN_START, .o = 0, .i = 1, .x = b_ };
            } else if (s_b_old.t == RUN_START) {
                s_b = s_b_old;
                s_b.i++;
            } else if (s_b_old.t == RUN_END) {
                pos_t p_1 = b_old + s_b_old.o;
                pos_t y_1 = idx_dir.M_LF().p(b__old + 1);

                if (p_1 < y_1) {
                    s_b = s_b_old;
                    s_b.i++;
                } else {
                    s_b = {
                        .t = RUN_END,
                        .o = y_1 - b_old - 1,
                        .i = 1,
                        .x = b__old
                    };
                }
            } else {
                s_b.t = NO_SAMPLE;
            }

            if (e_ != e__old) [[likely]] {
                s_e = { .t = RUN_END, .o = 0, .i = 1, .x = e_ };
            } else if (s_e_old.t == RUN_END) {
                s_e = s_e_old;
                s_e.i++;
            } else if (s_e_old.t == RUN_START) {
                pos_t p_2 = e_old - s_e_old.o;
                pos_t y_2 = idx_dir.M_LF().p(e__old);

                if (y_2 <= p_2) {
                    s_e = s_e_old;
                    s_e.i++;
                } else {
                    s_e = {
                        .t = RUN_START,
                        .o = e_old - y_2,
                        .i = 1,
                        .x = e__old
                    };
                }
            } else {
                s_e.t = NO_SAMPLE;
            }
        }

        using prime_tuple_t = std::tuple<
            pos_t, // b_prime
            pos_t, // e_prime
            pos_t, // b__prime
            pos_t // e__prime
        >;

        template <direction_t dir>
        inline void update_samples_stage_2(const move_rb<support, sym_t, pos_t>& idx, search_context_t& ctx, prime_tuple_t& prime) const
        {
            const auto& idx_dir = idx.index<dir>();
            const auto& [b_old, e_old, b_R_old, e_R_old,
                         b__old, e__old, b_R__old, e_R__old,
                         s_b_old, s_e_old, s_b_R_old, s_e_R_old] = const_vars<dir>();
            auto&& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s_b, s_e, s_b_R, s_e_R] = ctx.vars<dir>();
            auto& [b_prime, e_prime, b__prime, e__prime] = prime;
            
            if (e - b != e_prime - b_prime) [[likely]] {
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

            if (e - b != e_old - b_old) [[likely]] {
                if (s_b_R_old.t) {
                    pos_t p_1 = b_R_old + s_b_R_old.o;

                    if (!(b_R <= p_1 && p_1 <= e_R)) {
                        s_b_R.t = NO_SAMPLE;
                    } else {
                        s_b_R = s_b_R_old;
                        s_b_R.o -= b_R - b_R_old;
                    }
                } else {
                    s_b_R.t = NO_SAMPLE;
                }

                if (s_e_R_old.t) {
                    pos_t p_2 = e_R_old - s_e_R_old.o;

                    if (!(b_R <= p_2 && p_2 <= e_R)) {
                        s_e_R.t = NO_SAMPLE;
                    } else {
                        s_e_R = s_e_R_old;
                        s_e_R.o -= e_R_old - e_R;
                    }
                } else {
                    s_e_R.t = NO_SAMPLE;
                }
            } else {
                s_b_R = s_b_R_old;
                s_e_R = s_e_R_old;
            }

            assert(s_b.t || s_e.t || s_b_R.t || s_e_R.t);
        }

    public:
        using extend_res_t = std::tuple<search_context_t, bool>;

        /**
         * @brief extends the currently matched P with sym; if the extended P occurs in the input, true is
         * returned and the query context is adjusted to store the information for the extended P; else,
         * false is returned and the query context is not modified
         * @param idx the index to query
         * @param sym the symbol to extend P with
         * @param dir extend direction
         * @return whether symP occurs in the input
         */
        inline extend_res_t extend(const move_rb<support, sym_t, pos_t>& idx, sym_t sym, direction_t dir)
        {
            if (dir == LEFT) {
                return extend<LEFT>(idx, sym);
            } else {
                return extend<RIGHT>(idx, sym);
            }
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
        template <direction_t dir>
        inline extend_res_t extend(const move_rb<support, sym_t, pos_t>& idx, sym_t sym)
        {
            const auto& idx_dir = idx.index<dir>();
            i_sym_t i_sym = idx.idx_fwd.map_symbol(sym);

            // If i_sym does not occur in L', then P[i..m] does not occur in T
            if (i_sym == 0) [[unlikely]] {
                return {search_context_t{}, false};
            }

            const auto& [b_old, e_old, b_R_old, e_R_old,
                         b__old, e__old, b_R__old, e_R__old,
                         s_b_old, s_e_old, s_b_R_old, s_e_R_old] = const_vars<dir>();
            
            update_input_intervals_and_samples<dir>(idx);

            std::array<int64_t, 256> prev;
            std::array<int64_t, 256> next;
            build_prev_next<dir>(idx, prev, next, i_sym);

            if (next[i_sym] > prev[i_sym] || prev[i_sym] >= idx_dir.M_LF().num_intervals()) [[unlikely]] {
                return {search_context_t{}, false};
            }

            search_context_t ctx;
            auto&& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s_b, s_e, s_b_R, s_e_R] = ctx.vars<dir>();

            b_ = next[i_sym];
            e_ = prev[i_sym];

            if (b_ != b__old) [[likely]] {
                b = idx_dir.M_LF().p(b_);
            } else {
                b = b_old;
            }
            
            if (e_ != e__old) [[likely]] {
                e = idx_dir.M_LF().p(e_ + 1) - 1;
            } else {
                e = e_old;
            }

            if constexpr (query_support == LOCATE) {
                update_samples_stage_1<dir>(idx, ctx);
            }

            prime_tuple_t prime{b, e, b_, e_};
            b_R = b_R_old;

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
                update_samples_stage_2<dir>(idx, ctx, prime);
                ctx.shft = shft;
                ctx.dpth = dpth + 1;
            }

            ctx.dir_lst = dir;
            ctx.sym_lst = sym;
            ctx.m = m + 1;
            ctx.err = 0;

            return {ctx, true};
        }

        /**
         * @brief prepares the context to be extended
         * @param idx the index to query
         * @param dir extend direction
         * @return extend context prepared for extending the search context with every possible symbol in O(sigma) time
         */
        extend_context_t prepare_extend_all(const move_rb<support, sym_t, pos_t>& idx, direction_t dir)
        {
            extend_context_t ext_ctx;

            if (dir == LEFT) {
                prepare_extend_all<LEFT>(idx, ext_ctx);
            } else {
                prepare_extend_all<RIGHT>(idx, ext_ctx);
            }

            return ext_ctx;
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
        template <direction_t dir>
        void prepare_extend_all(const move_rb<support, sym_t, pos_t>& idx, extend_context_t& ext_ctx)
        {
            const auto& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s_b, s_e, s_b_R, s_e_R] = const_vars<dir>();
            auto& prev = ext_ctx.prev;
            auto& next = ext_ctx.next;

            update_input_intervals_and_samples<dir>(idx);
            build_prev_next<dir>(idx, prev, next, idx.sigma);
            
            ext_ctx.sym_nxt = 0;
            ext_ctx.b_R_nxt = b_R + (next[0] <= prev[0] && prev[0] != idx.index<dir>().M_LF().num_intervals());
            ext_ctx.template next_symbol<dir>(idx);
        }
        
        /**
         * @brief extends the context with the next symbol
         * @param idx the index to query
         * @param ext_ctx extend context to prepare
         * @param dir extend direction
         */
        search_context_t extend_next(const move_rb<support, sym_t, pos_t>& idx, extend_context_t& ext_ctx, direction_t dir)
        {
            if (dir == LEFT) {
                return extend_next<LEFT>(idx, ext_ctx);
            } else {
                return extend_next<RIGHT>(idx, ext_ctx);
            }
        }

        /**
         * @brief extends the context with the next symbol
         * @tparam dir extend direction
         * @param idx the index to query
         * @param ext_ctx extend context to prepare
         */
        template <direction_t dir>
        search_context_t extend_next(const move_rb<support, sym_t, pos_t>& idx, extend_context_t& ext_ctx) const
        {
            const auto& idx_dir = idx.index<dir>();
            i_sym_t i_sym = ext_ctx.sym_nxt;
            search_context_t ctx;

            const auto& [b_old, e_old, b_R_old, e_R_old,
                         b__old, e__old, b_R__old, e_R__old,
                         s_b_old, s_e_old, s_b_R_old, s_e_R_old] = const_vars<dir>();

            auto&& [b, e, b_R, e_R, b_, e_, b_R_, e_R_, s_b, s_e, s_b_R, s_e_R] = ctx.vars<dir>();

            b_ = ext_ctx.next[i_sym];
            e_ = ext_ctx.prev[i_sym];

            if (b_ != b__old) [[likely]] {
                b = idx_dir.M_LF().p(b_);
            } else {
                b = b_old;
            }
            
            if (e_ != e__old) [[likely]] {
                e = idx_dir.M_LF().p(e_ + 1) - 1;
            } else {
                e = e_old;
            }

            if constexpr (query_support == LOCATE) {
                update_samples_stage_1<dir>(idx, ctx);
            }

            prime_tuple_t prime{b, e, b_, e_};

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

            b_R = ext_ctx.b_R_nxt;
            e_R = ext_ctx.b_R_nxt + (e - b);

            b_R_ = b_R__old;
            e_R_ = e_R__old;

            if constexpr (query_support == LOCATE) {
                update_samples_stage_2<dir>(idx, ctx, prime);
                ctx.shft = shft;
                ctx.dpth = dpth + 1;
            }
            
            ctx.sym_lst = idx.unmap_symbol(i_sym);
            ctx.dir_lst = dir;
            ctx.m = m + 1;
            ctx.err = 0;

            ext_ctx.b_R_nxt = e_R + 1;
            ext_ctx.template next_symbol<dir>(idx);

            return ctx;
        }
    };

    /**
     * @brief structure storing the information for locating all occurrences of a search context
     */
    struct locate_context_t {
        protected:
        direction_t dir; // current locate direction_t
        pos_t occ_rem; // number of remaining occurrences to locate
        pos_t c; // initial position in the suffix array
        pos_t SA_c; // initial suffix SA[c] in the suffix array interval
        pos_t i; // current position in the suffix array interval
        pos_t SA_i; // current suffix SA[i] in the suffix array interval
        
        // index of the input inteval of M_Phi/M_Phi^{-1} containing SA_i
        std::conditional_t<support == _locate_move, pos_t, std::monostate> s_; 

        struct rlzsa_ctx_t {
            pos_t x_p, x_lp, x_cp, x_r, s_p; // variables for decoding the rlzsa
        };

        std::conditional_t<support == _locate_rlzsa, rlzsa_ctx_t, std::monostate> rlz_l; // rlzsa context for decoding to the left
        std::conditional_t<support == _locate_rlzsa, rlzsa_ctx_t, std::monostate> rlz_r; // rlzsa context for decoding to the right

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

        inline pos_t compute_center(const move_rb<support, sym_t, pos_t>& idx, const search_context_t<LOCATE>& ctx)
        {
            if (ctx.s_b.t) {
                c = ctx.b + ctx.s_b.o;
            } else if (ctx.s_e.t) {
                c = ctx.e - ctx.s_e.o;
            } else if (ctx.s_b_R.t) {
                if (ctx.s_b_R.t == RUN_START) {
                    c = idx._SA_sR_m1[ctx.s_b_R.x];
                } else {
                    c = idx._SA_eR_m1[ctx.s_b_R.x];
                }

                pos_t c_ = input_interval(c, idx._S_MLF_p_fwd,
                    [&](pos_t x){return idx.idx_fwd.M_LF().p(x);});

                for (pos_t j = 0; j < ctx.m - ctx.s_b_R.i; j++) {
                    idx.idx_fwd.M_LF().move(c, c_);
                }
            } else /* if (ctx.s_e_R.t) */ {
                if (ctx.s_e_R.t == RUN_START) {
                    c = idx._SA_sR_m1[ctx.s_e_R.x];
                } else {
                    c = idx._SA_eR_m1[ctx.s_e_R.x];
                }

                pos_t c_ = input_interval(c, idx._S_MLF_p_fwd,
                    [&](pos_t x){return idx.idx_fwd.M_LF().p(x);});

                for (pos_t j = 0; j < ctx.m - ctx.s_e_R.i; j++) {
                    idx.idx_fwd.M_LF().move(c, c_);
                }
            }

            assert(ctx.b <= c && c <= ctx.e);

            return c;
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
            if (dir == NO_DIR) {
                compute_center(idx, ctx);
            }

            i = c;

            if (ctx.s_b.t) {
                SA_i = (ctx.s_b.t == RUN_START ? idx.idx_fwd.SA_s(ctx.s_b.x) : idx.idx_fwd.SA_s_(ctx.s_b.x)) - ctx.s_b.i;
            } else if (ctx.s_e.t) {
                SA_i = (ctx.s_e.t == RUN_START ? idx.idx_fwd.SA_s(ctx.s_e.x) : idx.idx_fwd.SA_s_(ctx.s_e.x)) - ctx.s_e.i;
            } else if (ctx.s_b_R.t) {
                if (ctx.s_b_R.t == RUN_START) {
                    SA_i = idx.n - idx.idx_bwd.SA_s(ctx.s_b_R.x) - (ctx.m - ctx.s_b_R.i + 1);
                } else {
                    SA_i = idx.n - idx.idx_bwd.SA_s_(ctx.s_b_R.x) - (ctx.m - ctx.s_b_R.i + 1);
                }
            } else /* if (ctx.s_e_R.t) */ {
                if (ctx.s_e_R.t == RUN_START) {
                    SA_i = idx.n - idx.idx_bwd.SA_s(ctx.s_e_R.x) - (ctx.m - ctx.s_e_R.i + 1);
                } else {
                    SA_i = idx.n - idx.idx_bwd.SA_s_(ctx.s_e_R.x) - (ctx.m - ctx.s_e_R.i + 1);
                }
            }

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
            return occ + ctx.shft;
        }

        /**
         * @brief locates the remaining (not yet reported) occurrences of the currently matched P
         * @param Occ vector to append the occurrences to
         */
        template <typename report_fnc_t>
        inline void locate(const move_rb<support, sym_t, pos_t>& idx, const search_context_t<LOCATE>& ctx, report_fnc_t report)
        {
            static constexpr bool report_pos = function_traits<report_fnc_t>::arity > 1;

            if (occ_rem == ctx.num_occ()) {
                first_occ(idx, ctx);
                if constexpr (report_pos) report(c, SA_c + ctx.shft); else report(SA_c + ctx.shft);
            }

            if constexpr (support == _locate_move) {
                if (dir == LEFT) {
                    while (i > ctx.b) {
                        idx.idx_fwd.M_Phi().move(SA_i, s_);
                        i--;
                        if constexpr (report_pos) report(i, SA_i + ctx.shft); else report(SA_i + ctx.shft);
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
                        i++;
                        if constexpr (report_pos) report(i, SA_i + ctx.shft); else report(SA_i + ctx.shft);
                    }
                }
            } else if constexpr (support == _locate_rlzsa) {
                SA_i += ctx.shft;

                if (dir == LEFT) {
                    idx.idx_fwd.report_rlzsa_left(i, ctx.b, SA_i,
                        rlz_l.x_p, rlz_l.x_lp, rlz_l.x_cp, rlz_l.x_r, rlz_l.s_p,
                        report);

                    if (c < ctx.e) [[likely]] {
                        i = c + 1;
                        dir = RIGHT;
                        SA_i = SA_c + ctx.shft;
                    }
                }

                if (c < ctx.e) [[likely]] {
                    idx.idx_fwd.report_rlzsa_right(i, ctx.e, SA_i,
                        rlz_r.x_p, rlz_r.x_lp, rlz_r.x_cp, rlz_r.x_r, rlz_r.s_p, report);
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
            Occ.reserve(occ_rem);
            locate(idx, ctx, [&](pos_t occ){Occ.emplace_back(occ);});
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
         * @tparam dir direction_t the context should be extended with
         * @param idx the index to query
         */
        template <direction_t dir>
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
    inline search_context_t<query_support> empty_context() const
    {
        return search_context_t<query_support>(*this);
    }

    // ############################# APPROXIMATE PATTERN MATCHING #############################
    
    // type of hashset to store search contexts for hamming distance search
    template <move_rb_query_support_t query_support>
    using search_context_set_t = tsl::sparse_set<search_context_t<query_support>, typename search_context_t<query_support>::hash>;

    /**
     * @brief counts a pattern with at most k mismatches (w.r.t. hamming distance)
     * @param P the pattern to search
     * @param scheme the search scheme to use (provides k)
     * @return number of occurrences of P in T with at most k mismatches (w.r.t. hamming distance)
     */
    pos_t count_hamming_dist(const inp_t& P, const search_scheme_t& scheme) const
    {
        pos_t m = P.size();
        pos_t k = scheme.k;

        if (k >= m) return n - 1;
        if (k == 0) return forward_index().count(P);

        pos_t count = 0;

        auto ctxts_set = search_hamming_dist<COUNT>(P, scheme);
        for (const auto& ctx : ctxts_set) count += ctx.num_occ();

        return count;
    }

    /**
     * @brief locaes a pattern with at most k mismatches (w.r.t. hamming distance)
     * @tparam report_fnc_t type of the function report
     * @param P the pattern to search
     * @param scheme the search scheme to use (provides k)
     * @param report function that is called with every tuple (occ, err) as parameters, where occ
     *               is an occurrence of P in T with at most k mismatches (w.r.t. hamming distance)
     */
    template <typename report_fnc_t>
    void locate_hamming_dist(const inp_t& P, const search_scheme_t& scheme, report_fnc_t report) const
    {
        pos_t m = P.size();
        pos_t k = scheme.k;

        if (k == 0) {
            forward_index().locate(P, [&](pos_t occ){report({occ, m, 0});});
        } else {
            auto ctxts_set = search_hamming_dist<LOCATE>(P, scheme);

            for (const auto& ctx : ctxts_set) {
                ctx.locate_phase().locate(*this, ctx, [&](pos_t occ){report({occ, ctx.length(), ctx.errors()});});
            }
        }
    }

    /**
     * @brief locaes a pattern with at most k mismatches (w.r.t. hamming distance)
     * @param P the pattern to search
     * @param scheme the search scheme to use (provides k)
     * @return occurrences of P in T with at most k mismatches (w.r.t. hamming distance)
     */
    std::vector<aprx_occ_t<pos_t>> locate_hamming_dist(const inp_t& P, const search_scheme_t& scheme) const
    {
        std::vector<aprx_occ_t<pos_t>> Occ;
        locate_hamming_dist(P, scheme, [&](aprx_occ_t<pos_t> occ){Occ.emplace_back(occ);});
        return Occ;
    }
    
    /**
     * @brief executes a given search scheme for a pattern (w.r.t. hamming distance)
     * @param scheme the search scheme to use
     * @return all search contexts of P in T with at most k mismatches (assuming the schearch scheme covers all error configurations)
     */
    template <move_rb_query_support_t query_support>
    search_context_set_t<query_support> search_hamming_dist(const inp_t& P, const search_scheme_t& scheme) const
    {
        pos_t m = P.size();
        pos_t p = scheme.p;
        if (m == 0) return {empty_context<query_support>()};
        if (p >= m) throw std::runtime_error("p >= m");

        using search_state_t = std::tuple<
            search_context_t<query_support>, // ctx
            pos_t // k_cur
        >;

        using match_pos_t = std::tuple<
            pos_t, // p_idx
            direction_t, // dir
            pos_t, // part
            pos_t, // beg
            pos_t, // end
            pos_t, // pos
            pos_t // ext_rem
        >;

        using states_arr_t = std::vector<search_state_t>;

        search_context_set_t<query_support> ctxts_set;
        states_arr_t states_cur;
        states_arr_t states_nxt;
        auto empty_ctx = empty_context<query_support>();
        pos_t part_len = m / p;

        auto add_ctx = [&](search_context_t<query_support>& ctx){
            if constexpr (query_support == LOCATE) {
                auto it = ctxts_set.find(ctx);

                if (it != ctxts_set.end()) {
                    const auto& ctx_old = *it;

                    if (!(ctx_old.s_b.t || ctx_old.s_e.t) && (ctx.s_b.t || ctx.s_e.t)) {
                        it = ctxts_set.erase(it);
                    }
                }
                
                ctxts_set.emplace_hint(it, ctx);
            } else {
                ctxts_set.emplace(ctx);
            }
        };

        for (const search_t& S : scheme.S) {
            match_pos_t match_pos_cur;
            auto& [p_idx, dir, part, beg, end, pos, ext_rem] = match_pos_cur;

            p_idx = 0;
            part = S[p_idx].part;
            beg = part * part_len;
            end = (part + 1) == p ? m : ((part + 1) * part_len);
            dir = S[p_idx + 1].part < part ? LEFT : RIGHT;
            pos = dir == LEFT ? end - 1 : beg;
            ext_rem = dir == LEFT ? pos - beg : (end - pos - 1);
            states_cur.emplace_back(search_state_t{empty_ctx, 0});

            match_pos_t match_pos_nxt {p_idx, dir, part, beg, end, pos, ext_rem};
            auto& [p_idx_nxt, dir_nxt, part_nxt, beg_nxt, end_nxt, pos_nxt, ext_rem_nxt] = match_pos_nxt;

            while (true) {
                pos_nxt = pos + (dir == LEFT ? -1 : 1);

                if (!(beg <= pos_nxt && pos_nxt < end)) [[unlikely]] {
                    p_idx_nxt = p_idx + 1;

                    if (p_idx_nxt < p) [[likely]] {
                        part_nxt = S[p_idx_nxt].part;
                        dir_nxt = part_nxt < part ? LEFT : RIGHT;
                        beg_nxt = part_nxt * part_len;
                        end_nxt = (part_nxt + 1) == p ? m : ((part_nxt + 1) * part_len);
                        pos_nxt = dir_nxt == LEFT ? end_nxt - 1 : beg_nxt;
                        ext_rem = end_nxt - beg_nxt;
                    }
                }
                
                while (!states_cur.empty()) {
                    auto [ctx, k_cur] = states_cur.back(); states_cur.pop_back();
                    
                    if (k_cur < S[p_idx].k_max && (p_idx_nxt == p || S[p_idx_nxt].k_min <= k_cur + 1 + ext_rem)) {
                        auto ext_ctx = ctx.prepare_extend_all(*this, dir);

                        while (ext_ctx.can_extend(*this)) {
                            auto ctx_nxt = ctx.extend_next(*this, ext_ctx, dir);
                            bool is_mismatch = ctx_nxt.last_added_symbol() != P[pos];
                            pos_t k_nxt = k_cur + is_mismatch;

                            if (p_idx_nxt == p) {
                                ctx_nxt.set_errors(k_nxt);
                                add_ctx(ctx_nxt);
                            } else if (is_mismatch || S[p_idx_nxt].k_min <= k_nxt + ext_rem) {
                                states_nxt.emplace_back(search_state_t{ctx_nxt, k_nxt});
                            }
                        }
                    } else {
                        auto [ctx_nxt, result] = ctx.extend(*this, P[pos], dir);

                        if (result) {
                            if (p_idx_nxt == p) {
                                ctx_nxt.set_errors(k_cur);
                                add_ctx(ctx_nxt);
                            } else {
                                states_nxt.emplace_back(search_state_t{ctx_nxt, k_cur});
                            }
                        }
                    }
                }

                if (p_idx == p) [[unlikely]] break;
                ext_rem_nxt = ext_rem - 1;
                match_pos_cur = match_pos_nxt;
                states_cur.swap(states_nxt);
            }
        }

        return ctxts_set;
    }

    struct interval_t {pos_t beg; pos_t end; pos_t err; pos_t len;};

    static std::vector<interval_t> effective_intervals(const std::vector<const search_context_t<LOCATE>*>& ctxts_sorted)
    {
        std::vector<interval_t> iv_stack;
        std::vector<interval_t> ivs;

        iv_stack.reserve(16);
        ivs.reserve(ctxts_sorted.size());
    
        auto append_segment = [&](interval_t iv) {
            auto [b, e, err, len] = iv;
            if (b > e) return;

            if (!ivs.empty()) {
                interval_t& seg_lst = ivs.back();
                const auto& [b_lst, e_lst, err_lst, len_lst] = seg_lst;

                if (err_lst == err && len_lst == len && e_lst + 1 == b) {
                    seg_lst.end = e;
                    return;
                }
            }

            ivs.emplace_back(interval_t{b, e, err, len});
        };

        for (const auto* ctx : ctxts_sorted) {
            auto [beg, end] = ctx->forward_sa_interval();
            pos_t err = ctx->errors();
            pos_t len = ctx->depth();

            while (!iv_stack.empty() && iv_stack.back().end < beg) {
                interval_t iv = iv_stack.back(); iv_stack.pop_back();
                append_segment(iv);

                if (!iv_stack.empty()) {
                    iv_stack.back().beg = iv.end + 1;
                }
            }

            if (iv_stack.empty()) {
                iv_stack.emplace_back(interval_t{beg, end, err, len});
            } else {
                interval_t& iv_top = iv_stack.back();
                append_segment({iv_top.beg, beg - 1, iv_top.err, iv_top.len});
                iv_top.beg = beg;

                if (iv_top.err < err || (iv_top.err == err && iv_top.len <= len)) {
                    iv_stack.emplace_back(interval_t{beg, end, iv_top.err, iv_top.len});
                } else {
                    iv_stack.emplace_back(interval_t{beg, end, err, len});
                }
            }
        }

        while (!iv_stack.empty()) {
            interval_t iv = iv_stack.back(); iv_stack.pop_back();
            append_segment(iv);

            if (!iv_stack.empty()) {
                iv_stack.back().beg = iv.end + 1;
            }
        }

        return std::move(ivs);
    }

    /**
     * @brief locaes a pattern with at most k mismatches (w.r.t. edit distance)
     * @tparam report_fnc_t type of the function report
     * @param P the pattern to search
     * @param scheme the search scheme to use (provides k)
     * @param report function that is called with every tuple (occ, err) as parameters, where occ
     *               is an occurrence of P in T with at most k mismatches (w.r.t. edit distance)
     */
    template <typename report_fnc_t>
    void locate_edit_dist(const inp_t& P, const search_scheme_t& scheme, report_fnc_t report) const requires(supports_locate)
    {
        pos_t k = scheme.k;
        pos_t m = P.size();

        if (k == 0) {
            forward_index().locate(P, [&](pos_t occ){report({occ, m, 0});});
            return;
        }

        auto ctxts_set = search_edit_dist_entry(P, scheme);
        std::vector<const search_context_t<LOCATE>*> ctxts_sorted;
        ctxts_sorted.reserve(ctxts_set.size());
        for (const auto& ctx : ctxts_set) ctxts_sorted.emplace_back(&ctx);
        ips4o::sort(ctxts_sorted.begin(), ctxts_sorted.end(), [](auto x, auto y){return *x < *y;});
        auto ivs = effective_intervals(ctxts_sorted);
        pos_t iv_idx = 0;

        for (auto it = ctxts_sorted.begin(); it != ctxts_sorted.end();) {
            const auto& ctx = **it;
            auto [beg, end] = ctx.forward_sa_interval();
            auto loc_ctx = ctx.locate_phase();
            pos_t cntr = loc_ctx.compute_center(*this, ctx);
            pos_t iv_idx_cntr = exp_search_max_leq<pos_t, RIGHT>(cntr, iv_idx, ivs.size() - 1, [&](pos_t x) {return ivs[x].beg;});
            iv_idx = iv_idx_cntr;

            loc_ctx.locate(*this, ctx, [&](pos_t pos, pos_t occ){
                if (pos <= cntr) {
                    if (pos < ivs[iv_idx].beg) [[unlikely]] iv_idx--;
                } else {
                    if (iv_idx < iv_idx_cntr) [[unlikely]] iv_idx = iv_idx_cntr;
                    if (ivs[iv_idx].end < pos) [[unlikely]] iv_idx++;
                }

                report({occ, ivs[iv_idx].len, ivs[iv_idx].err});
            });

            iv_idx = std::max<pos_t>(iv_idx, iv_idx_cntr);
            do {it++;} while (it != ctxts_sorted.end() && std::get<0>((*it)->forward_sa_interval()) <= end);
        }
    }

    /**
     * @brief locaes a pattern with at most k mismatches (w.r.t. edit distance)
     * @param P the pattern to search
     * @param scheme the search scheme to use (provides k)
     * @return occurrences of P in T with at most k mismatches (w.r.t. edit distance)
     */
    std::vector<aprx_occ_t<pos_t>> locate_edit_dist(const inp_t& P, const search_scheme_t& scheme) const
    {
        std::vector<aprx_occ_t<pos_t>> Occ;
        locate_edit_dist(P, scheme, [&](aprx_occ_t<pos_t> occ){Occ.emplace_back(occ);});
        return Occ;
    }

    class cluster_t
    {
      private:
        std::vector<uint16_t> dists;
        std::vector<search_context_t<LOCATE>> nodes;

        pos_t last_cell;
        uint16_t k_max;
        pos_t start_depth;
        pos_t shift;

      public:
        cluster_t(pos_t size, pos_t k_max, pos_t start_depth, pos_t shift)
            : dists(size, k_max + 1), nodes(size), last_cell(-1), k_max(k_max),
            start_depth(start_depth), shift(shift) {}

        void set_cell(pos_t idx, search_context_t<LOCATE> node, uint16_t dist)
        {
            nodes[idx] = std::move(node);
            dists[idx] = dist;
            last_cell = idx;
        }

        pos_t size() const
        {
            return dists.size();
        }

        std::vector<search_context_t<LOCATE>> report_centers_at_end()
        {
            std::vector<search_context_t<LOCATE>> centers;
            centers.reserve(last_cell + 1);

            for (pos_t i = 0; i <= last_cell; i++) {
                if (!nodes[i].reported() &&
                    dists[i] <= k_max &&
                    (i == 0 || dists[i] <= dists[i - 1]) &&
                    (i == last_cell || dists[i] <= dists[i + 1])
                ) {
                    nodes[i].set_reported(true);
                    search_context_t<LOCATE> ctx = nodes[i];
                    ctx.set_errors(dists[i]);
                    ctx.set_depth(ctx.depth() + start_depth);
                    ctx.set_shift(shift);
                    centers.emplace_back(ctx);
                }
            }

            return centers;
        }

        std::tuple<search_context_t<LOCATE>, bool> deepest_minimum(direction_t dir)
        {
            uint16_t k_min = k_max + 1;
            pos_t max_best = -1;
            pos_t min_best = -1;

            for (pos_t i = 0; i <= last_cell; i++) {
                if (dists[i] < k_min) {
                    k_min = dists[i];
                    max_best = i;
                    min_best = i;
                }
                
                if (dists[i] == k_min) {
                    min_best = i;
                }
            }

            search_context_t<LOCATE> ctx;

            if (k_min <= k_max && !nodes[min_best].reported()) {
                nodes[min_best].set_reported(true);
                ctx = nodes[min_best];
                ctx.set_errors(k_min);
                ctx.set_depth(ctx.depth() + start_depth - (min_best - max_best));
                ctx.set_shift(((dir == LEFT) ? (min_best - max_best) : 0) + shift);
                return {ctx, true};
            }
            
            return {ctx, false};
        }

        std::tuple<search_context_t<LOCATE>, bool> cluster_center(uint16_t k_min, std::vector<search_context_t<LOCATE>>& desc, std::vector<uint16_t>& init_dists)
        {
            desc.reserve(dists.size());
            dists.reserve(dists.size());
            search_context_t<LOCATE> ctx;
            
            for (pos_t i = 0; i <= last_cell; i++) {
                if (dists[i] > k_max || dists[i] < k_min) continue;

                if (((i == 0) || dists[i] <= dists[i - 1]) &&
                    ((i == last_cell) || dists[i] <= dists[i + 1])
                ) {
                    if (!nodes[i].reported()) {
                        ctx = nodes[i];
                        ctx.set_errors(dists[i]);
                        ctx.set_depth(ctx.depth() + start_depth);
                        ctx.set_shift(shift);
                    }
                    
                    init_dists.emplace_back(dists[i]);

                    for (pos_t j = i + 1; j <= last_cell; j++) {
                        desc.emplace_back(nodes[j]);
                        init_dists.emplace_back(dists[j]);
                    }

                    for (pos_t k = 1; k < init_dists.size(); k++) {
                        if (init_dists[k] < k_min && init_dists[k] <= init_dists[k - 1] &&
                           (k == init_dists.size() - 1 || init_dists[k] <= init_dists[k + 1])
                        ) {
                            pos_t high = 0;
                            pos_t low = init_dists.size() - 1;

                            for (pos_t l = k; l-- > 0;) {
                                if (init_dists[l] != init_dists[l + 1] + 1) {
                                    high = l + 1;
                                    break;
                                }
                            }
                            
                            for (pos_t l = k + 1; l < init_dists.size(); l++) {
                                if (init_dists[l] != init_dists[l - 1] + 1) {
                                    low = l - 1;
                                    break;
                                }
                            }

                            if (high != 0 && low != init_dists.size() - 1) {
                                pos_t lC = low;
                                pos_t hC = high;
                                bool highest = true;

                                while (lC > hC) {
                                    if (highest) {
                                        init_dists[hC] = std::min(k_max + 1, init_dists[hC - 1] + 1);
                                        hC++;
                                    } else {
                                        init_dists[lC] = std::min(k_max + 1, init_dists[lC + 1] + 1);
                                        lC--;
                                    }

                                    highest = !highest;
                                }

                                if (lC == hC) {
                                    init_dists[lC] = std::min(init_dists[lC + 1] + 1, init_dists[lC - 1] + 1);
                                }
                            } else if (high == 0 && low != init_dists.size() - 1) {
                                for (pos_t l = low; l-- > 0;) {
                                    init_dists[l] = init_dists[l + 1] + 1;
                                }
                            } else if (high != 0 && low == init_dists.size() - 1) {
                                for (pos_t l = high; l < init_dists.size(); l++) {
                                    init_dists[l] = init_dists[l - 1] + 1;
                                }
                            }
                        }
                    }

                    return {ctx, true};
                }
            }

            return {ctx, false};
        }
    };

    template <typename matrix_word_t>
    using data_struct_ref_t = std::tuple<
        const edit_dist_search&,
        const std::vector<sub_string<pos_t>>&,
        std::vector<std::vector<search_context_t<LOCATE>>>&,
        std::vector<edit_distance_matrix<matrix_word_t>>&,
        std::vector<edit_distance_matrix<matrix_word_t>>&,
        search_context_set_t<LOCATE>&
    >;

    /**
     * @brief executes a given search scheme for a pattern (w.r.t. edit distance)
     * @param scheme the search scheme to use
     * @return all search contexts of P in T with at most k mismatches (assuming the schearch scheme covers all error configurations)
     */
    search_context_set_t<LOCATE> search_edit_dist_entry(const inp_t& P, const search_scheme_t& scheme) const requires(supports_locate)
    {
        static constexpr uint16_t k_limit_64 = edit_distance_matrix<uint64_t>::k_limit; // 10
        static constexpr uint16_t k_limit_128 = edit_distance_matrix<uint128_t>::k_limit; // 20

        if (scheme.k <= k_limit_64) {
            return search_edit_dist<uint64_t>(P, scheme);
        } else if (scheme.k <= k_limit_128) {
            return search_edit_dist<uint128_t>(P, scheme);
        } else {
            throw std::runtime_error("cannot search patterns with more than " + std::to_string(k_limit_128) + " errors.");
        }
    }

    /**
     * @brief executes a given search scheme for a pattern (w.r.t. edit distance)
     * @param scheme the search scheme to use
     * @return all search contexts of P in T with at most k mismatches (assuming the schearch scheme covers all error configurations)
     */
    template <typename matrix_word_t>
    search_context_set_t<LOCATE> search_edit_dist(const inp_t& P, const search_scheme_t& scheme) const requires(supports_locate)
    {
        pos_t m = P.size();
        pos_t p = scheme.p;
        if (m == 0) return {empty_context<LOCATE>()};
        if (p >= m) throw std::runtime_error("p >= m");

        search_context_set_t<LOCATE> ctxts_set;
        std::vector<edit_dist_search> searches;
        std::vector<sub_string<pos_t>> parts;
        searches.reserve(scheme.S.size());
        parts.reserve(p);

        double part_len = (double) m / p;
        
        for (pos_t p_idx = 0; p_idx < p; p_idx++) {
            pos_t part = scheme.S[0][p_idx].part;
            pos_t beg = part * part_len;
            pos_t end = (p_idx == p - 1 ? m : (part + 1) * part_len) - 1;
            parts.emplace_back(sub_string<pos_t>(P, beg, end));
        }
        
        for (pos_t s_idx = 0; s_idx < scheme.S.size(); s_idx++) {
            searches.emplace_back(edit_dist_search(scheme, s_idx));
        }

        std::vector<std::vector<search_context_t<LOCATE>>> stacks;
        stacks.resize(p);
        
        std::vector<edit_distance_matrix<matrix_word_t>> matrices_fwd, matrices_bwd;
        matrices_fwd.resize(p);
        matrices_bwd.resize(p);

        for (const edit_dist_search& search : searches) {
            data_struct_ref_t<matrix_word_t> ds_ref(search, parts, stacks, matrices_fwd, matrices_bwd, ctxts_set);

            for (pos_t p_idx = 0; p_idx < p; p_idx++) {
                parts[search.part(p_idx)].set_direction(search.part_dir(p_idx));
            }
            
            auto ctx = empty_context<LOCATE>();
            pos_t p_idx = 0;

            if (search.upper_bound(0) == 0) {
                bool match = true;

                for (; match && p_idx < parts.size(); p_idx++) {
                    if (search.upper_bound(p_idx) > 0) break;
                    const auto& part = parts[search.part(p_idx)];

                    for (pos_t i = 0; match && i < part.size(); i++) {
                        auto [ctx_nxt, match_nxt] = ctx.extend(*this, part[i], part.direction());
                        if (match_nxt) ctx = ctx_nxt;
                        match = match_nxt;
                    }
                }

                if (!match) continue;
            }

            rec_edit_dist_search<matrix_word_t>(ctx, p_idx, ds_ref, {}, {}, {}, {});
        }

        return ctxts_set;
    }

    template <typename matrix_word_t>
    void rec_edit_dist_search(
        search_context_t<LOCATE>& ctx, pos_t p_idx,
        data_struct_ref_t<matrix_word_t>& ds_ref,
        const std::vector<search_context_t<LOCATE>>& desc_prev_dir,
        const std::vector<search_context_t<LOCATE>>& desc_diff_dir,
        const std::vector<uint16_t>& dists_prev_dir,
        const std::vector<uint16_t>& dists_diff_dir
    ) const requires(supports_locate) {
        auto& [search, parts, stacks, matrices_fwd, matrices_bwd, ctxts_set] = ds_ref;

        const auto& part = parts[search.part(p_idx)];
        pos_t k_max = search.upper_bound(p_idx);
        direction_t dir = search.part_dir(p_idx);
        bool dir_switch = search.does_part_switch_dir(p_idx);
        auto& stack = stacks[p_idx];
        auto& matrix = dir == RIGHT ? matrices_fwd[search.part(p_idx)] : matrices_bwd[search.part(p_idx)];

        const auto& dists_dir = dir_switch ? dists_diff_dir : dists_prev_dir;
        const auto& desc_dir = dir_switch ? desc_diff_dir : desc_prev_dir;
        const auto& dists_rev_dir = dir_switch ? dists_prev_dir : dists_diff_dir;
        const auto& desc_rev_dir = dir_switch ? desc_prev_dir : desc_diff_dir;

        std::vector<uint16_t> dists_new;

        if (dists_dir.empty()) {
            dists_new.emplace_back(ctx.errors());
        } else {
            uint16_t prev_dist = dir_switch ? *min_element(dists_dir.begin(), dists_dir.end()) : dists_dir[0];
            uint16_t add = ctx.errors() - prev_dist;
            dists_new.reserve(dists_dir.size());

            for (pos_t i = 0; i < dists_dir.size(); i++) {
                dists_new.emplace_back(dists_dir[i] + add);
            }
        }

        if (!matrix.is_initialized()) matrix.set_input(part, sigma, idx_fwd._map_int);
        matrix.init(k_max, dists_new);
        cluster_t cluster(matrix.last_col_size(), k_max, ctx.depth(), ctx.shift());

        if (matrix.is_in_final_column(0)) {
            auto ctx_copy = ctx;
            ctx_copy.set_depth(0);
            cluster.set_cell(0, ctx_copy, matrix(0, part.size()));
        }

        pos_t old_depth = 0;

        if (!desc_dir.empty()) {
            for (pos_t i = 0; i < desc_dir.size() && desc_dir[i].depth() < matrix.num_rows(); i++) {
                if (branch_and_bound<matrix_word_t>(matrix, cluster, desc_dir[i], p_idx, ds_ref,
                    dists_rev_dir, desc_rev_dir, {desc_dir.begin() + i + 1, desc_dir.end()})
                ) return;
            }

            if (desc_dir.back().depth() == matrix.num_rows() - 1) return;
            if (!dir_switch) ctx = desc_dir.back();
            old_depth = desc_dir.back().depth();
        }

        auto ext_ctx = ctx.prepare_extend_all(*this, dir);

        while (ext_ctx.can_extend(*this)) {
            stack.emplace_back(ctx.extend_next(*this, ext_ctx, dir));
            stack.back().set_depth(old_depth + 1);
        }

        while (!stack.empty()) {
            auto ctx_cur = stack.back(); stack.pop_back();
            if (branch_and_bound<matrix_word_t>(matrix, cluster, ctx_cur, p_idx, ds_ref, dists_rev_dir, desc_rev_dir, {})) continue;
            auto ext_ctx_cur = ctx_cur.prepare_extend_all(*this, dir);

            while (ext_ctx_cur.can_extend(*this)) {
                stack.emplace_back(ctx_cur.extend_next(*this, ext_ctx_cur, dir));
            }
        }
    }
    
    template <typename matrix_word_t>
    bool branch_and_bound(
        edit_distance_matrix<matrix_word_t>& matrix, cluster_t& cluster,
        const search_context_t<LOCATE>& ctx, pos_t p_idx,
        data_struct_ref_t<matrix_word_t>& ds_ref,
        const std::vector<uint16_t>& dists_rev_dir,
        const std::vector<search_context_t<LOCATE>>& desc_rev_dir,
        const std::vector<search_context_t<LOCATE>>& desc_dir_rem
    ) const requires(supports_locate) {
        auto& [search, parts, stacks, matrices_fwd, matrices_bwd, ctxts_set] = ds_ref;

        pos_t row = ctx.depth();
        bool valid_dist = matrix.compute_row(row, ctx.last_added_symbol());

        if (matrix.is_in_final_column(row)) {
            pos_t cluster_idx = cluster.size() + row - matrix.num_rows();
            cluster.set_cell(cluster_idx, ctx, matrix(row, matrix.num_cols() - 1));

            if (!valid_dist || matrix.only_vertical_gaps_left(row)) {
                const auto& k_min = search.lower_bound(p_idx);

                if (search.is_edge(p_idx)) {
                    if (p_idx + 1 == parts.size()) {
                        auto centers = cluster.report_centers_at_end();

                        for (auto& ctx_cur : centers) {
                            if (ctx_cur.errors() >= k_min) {
                                auto it = ctxts_set.find(ctx_cur);

                                if (it != ctxts_set.end()) {
                                    const auto& ctx_old = *it;

                                    if (ctx_cur.errors() < ctx_old.errors() ||
                                        (ctx_cur.errors() == ctx_old.errors() && ctx_cur.depth() > ctx_old.depth()) ||
                                        (ctx_cur.errors() == ctx_old.errors() && ctx_cur.depth() == ctx_old.depth() && !(ctx_old.s_b.t || ctx_old.s_e.t) && (ctx_cur.s_b.t || ctx_cur.s_e.t))
                                    ) {
                                        it = ctxts_set.erase(it);
                                        ctxts_set.emplace_hint(it, ctx_cur);
                                    }
                                } else {
                                    ctxts_set.emplace(ctx_cur);
                                }
                            }
                        }
                    } else {
                        auto [ctx_min, res] = cluster.deepest_minimum(search.part_dir(p_idx));

                        if (res && ctx_min.errors() >= k_min)
                            rec_edit_dist_search<matrix_word_t>(ctx_min, p_idx + 1, ds_ref, {}, desc_rev_dir, {}, dists_rev_dir);
                    }

                    return true;
                }

                std::vector<search_context_t<LOCATE>> desc_dir;
                std::vector<uint16_t> dists_dir;
                auto [ctx_cntr, res] = cluster.cluster_center(k_min, desc_dir, dists_dir);
                if (!res) return true;

                desc_dir.insert(desc_dir.end(), desc_dir_rem.begin(), desc_dir_rem.end());

                for (pos_t i = 0; i < desc_dir.size(); i++) {
                    desc_dir[i].set_depth(i + 1);
                }

                pos_t k_max = search.upper_bound(p_idx + 1);
                while (dists_dir.back() > k_max) dists_dir.pop_back();

                if (search.does_part_switch_dir(p_idx + 1) && !desc_dir.empty()) {
                    pos_t depth_cntr = ctx_cntr.depth();
                    ctx_cntr = desc_dir.back();
                    ctx_cntr.set_depth(depth_cntr);
                    ctx_cntr.set_errors(*std::min_element(dists_dir.begin(), dists_dir.end()));
                }

                rec_edit_dist_search<matrix_word_t>(ctx_cntr, p_idx + 1, ds_ref, desc_dir, desc_rev_dir, dists_dir, dists_rev_dir);
                return true;
            }
        }

        return !valid_dist;
    }

    /**
     * @brief locaes a pattern with at most k mismatches (w.r.t. dist_metr)
     * @tparam report_fnc_t type of the function report
     * @param P the pattern to search
     * @param scheme the search scheme to use (provides k)
     * @param report function that is called with every tuple (occ, err) as parameters, where occ
     *               is an occurrence of P in T with at most k mismatches (w.r.t. dist_metr)
     */
    template <distance_metric_t dist_metr, typename report_fnc_t>
    void locate(const inp_t& P, const search_scheme_t& scheme, report_fnc_t report) const requires(supports_locate)
    {
        if constexpr (dist_metr == HAMMING_DISTANCE) locate_hamming_dist(P, scheme, report);
        if constexpr (dist_metr == EDIT_DISTANCE) locate_edit_dist(P, scheme, report);
    }

    /**
     * @brief locaes a pattern with at most k mismatches (w.r.t. dist_metr)
     * @param P the pattern to search
     * @param scheme the search scheme to use (provides k)
     * @return occurrences of P in T with at most k mismatches (w.r.t. dist_metr)
     */
    template <distance_metric_t dist_metr>
    std::vector<aprx_occ_t<pos_t>> locate(const inp_t& P, const search_scheme_t& scheme) const requires(supports_locate)
    {
        std::vector<aprx_occ_t<pos_t>> Occ;
        locate<dist_metr>(P, scheme, [&](aprx_occ_t<pos_t> occ){Occ.emplace_back(occ);});
        return Occ;
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