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
        constexpr_case<support == _locate_lzendsa, move_r<_locate_lzendsa_bi_fwd, sym_t, pos_t>>,
     /* constexpr_case<support == _count, */       move_r<_count_bi,              sym_t, pos_t>>;

    using move_r_bwd_t = constexpr_switch_t< // backward move_r index type
        constexpr_case<support == _locate_move,    move_r<_locate_bi_bwd, sym_t, pos_t>>,
        constexpr_case<support == _locate_rlzsa,   move_r<_locate_bi_bwd, sym_t, pos_t>>,
        constexpr_case<support == _locate_lzendsa, move_r<_locate_bi_bwd, sym_t, pos_t>>,
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

    enum direction { _LEFT, _RIGHT, _NONE };

    template <direction dir>
    const std::conditional_t<dir == _LEFT, move_r_fwd_t, move_r_bwd_t>& index() const
    {
        if constexpr (dir == _LEFT) {
            return idx_fwd;
        } else {
            return idx_bwd;
        }
    }

    template <direction dir>
    const sd_array<pos_t>& S_MLF_p() const
    {
        if constexpr (dir == _LEFT) {
            return _S_MLF_p_fwd;
        } else {
            return _S_MLF_p_bwd;
        }
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
        auto time = now();
        auto time_start = time;
        std::cout << "preprocessing T" << std::flush;

        uint16_t p = params.num_threads;
        std::vector<std::vector<pos_t>> freq_thr(p, std::vector<pos_t>(256, 0));

        if (params.file_input) {
            n = std::filesystem::file_size(input) + 1;
            std::vector<sdsl::int_vector_buffer<8>> T_buf;

            for (uint16_t i = 0; i < p; i++) {
                T_buf.emplace_back(sdsl::int_vector_buffer<8>(input, std::ios::in, 128 * 1024, 8, true));
            }

            #pragma omp parallel for num_threads(p)
            for (uint64_t i = 0; i < n - 1; i++) {
                freq_thr[omp_get_thread_num()][T_buf[omp_get_thread_num()][i]]++;
            }
        } else {
            n = input.size() + 1;
            
            #pragma omp parallel for num_threads(p)
            for (uint64_t i = 0; i < n - 1; i++) {
                freq_thr[omp_get_thread_num()][char_to_uchar(input[i])]++;
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
        std::vector<uint8_t> uchars_sorted;

        for (uint16_t c = 0; c < 256; c++) {
            if (freq[c] != 0) {
                uchars_sorted.emplace_back(c);
            }
        }

        std::sort(uchars_sorted.begin(), uchars_sorted.end(),
            [&](uint8_t c1, uint8_t c2){
                return freq[c1] > freq[c2];
        });

        params.alphabet_size = uchars_sorted.size();
        sigma = params.alphabet_size + 1;
        uint8_t cur_uchar = params.mode == _bigbwt ? 3 : 1;
        std::vector<uint8_t> map_int(256, 0);

        for (uint8_t c : uchars_sorted) {
            map_int[c] = cur_uchar;
            cur_uchar++;
        }

        time = log_runtime(time);
        std::cout << "mapping T to its internal alphabet" << std::flush;

        if (params.file_input) {
            std::vector<sdsl::int_vector_buffer<8>> T_buf;

            for (uint16_t i = 0; i < p; i++) {
                T_buf.emplace_back(sdsl::int_vector_buffer<8>(input, std::ios::in, 128 * 1024, 8, true));
            }

            #pragma omp parallel for num_threads(p)
            for (uint64_t i = 0; i < n - 1; i++) {
                T_buf[omp_get_thread_num()][i] = map_int[T_buf[omp_get_thread_num()][i]];
            }
        } else {
            #pragma omp parallel for num_threads(p)
            for (uint64_t i = 0; i < n - 1; i++) {
                input[i] = uchar_to_char(map_int[char_to_uchar(input[i])]);
            }
        }

        std::vector<char> map_ext(256, 0);

        for (uint16_t c = 0; c < 256; c++) {
            if (map_int[c] != 0) {
                map_int[c] -= 2;
                map_ext[map_int[c]] = c;
            }
        }

        time = log_runtime(time);

        if (params.file_input) {
            reverse(input, p, false);

            std::cout << std::endl;
            std::cout << "building backward index:" << std::endl;

            idx_bwd = move_r_bwd_t(input, params);

            std::cout << std::endl;
            std::cout << "storing backward index to disk" << std::flush;
            time = now();

            std::string idx_bwd_file_name = random_alphanumeric_string(10);
            std::ofstream idx_bwd_ofile(idx_bwd_file_name);
            idx_bwd.set_alphabet_maps(map_int, map_ext);
            idx_bwd.serialize(idx_bwd_ofile);
            _S_MLF_p_bwd.serialize(idx_bwd_ofile);
            _SA_sR_m1.serialize(idx_bwd_ofile);
            _SA_eR_m1.serialize(idx_bwd_ofile);
            idx_bwd = move_r_bwd_t();
            _S_MLF_p_bwd = sd_array<pos_t>();
            _SA_sR_m1 = interleaved_byte_aligned_vectors<pos_t, pos_t>();
            _SA_eR_m1 = interleaved_byte_aligned_vectors<pos_t, pos_t>();
            idx_bwd_ofile.close();
            
            time = log_runtime(time);

            reverse(input, p, false);
            
            std::cout << std::endl;
            std::cout << "building forward index:" << std::endl;

            idx_fwd = move_r_fwd_t(input, params);
            idx_fwd.set_alphabet_maps(map_int, map_ext);

            std::cout << std::endl;
            std::cout << "loading backward index from disk" << std::flush;
            time = now();

            std::ifstream idx_bwd_ifile(idx_bwd_file_name);
            idx_bwd.load(idx_bwd_ifile);
            _S_MLF_p_bwd.load(idx_bwd_ifile);
            _SA_sR_m1.load(idx_bwd_ifile);
            _SA_eR_m1.load(idx_bwd_ifile);
            idx_bwd_ifile.close();
            
            time = log_runtime(time);
        } else {
            reverse(input, p, true);
            idx_bwd = move_r_bwd_t(input, params);
            idx_bwd.set_alphabet_maps(map_int, map_ext);
            reverse(input, p, true);
            idx_fwd = move_r_fwd_t(input, params);
            idx_fwd.set_alphabet_maps(map_int, map_ext);
        }
        
        time = now();
        std::cout << std::endl;
        std::cout << "unmapping T from its internal alphabet" << std::flush;

        for (uint8_t c = 255; c >= 2; c--) {
            map_ext[c] = map_ext[c - 2];
        }

        map_ext[0] = 0;
        map_ext[1] = 0;

        if (params.file_input) {
            std::vector<sdsl::int_vector_buffer<8>> T_buf;

            for (uint16_t i = 0; i < p; i++) {
                T_buf.emplace_back(sdsl::int_vector_buffer<8>(input, std::ios::in, 128 * 1024, 8, true));
            }

            #pragma omp parallel for num_threads(p)
            for (uint64_t i = 0; i < n - 1; i++) {
                T_buf[omp_get_thread_num()][i] = char_to_uchar(map_ext[T_buf[omp_get_thread_num()][i]]);
            }
        } else {
            #pragma omp parallel for num_threads(p)
            for (uint64_t i = 0; i < n - 1; i++) {
                input[i] = map_ext[char_to_uchar(input[i])];
            }
        }

        time = log_runtime(time);
        std::cout << "building S_MLF_p_fwd" << std::flush;
        
        pos_t n = idx_fwd.input_size() + 1;
        
        _S_MLF_p_fwd = build_sampling(
            idx_fwd.M_LF().num_intervals(), n, sample_rate_input_intervals,
            [&](pos_t i){return idx_fwd.M_LF().p(i);});

        time = log_runtime(time);
        std::cout << "building S_MLF_p_bwd" << std::flush;

        _S_MLF_p_bwd = build_sampling(
            idx_bwd.M_LF().num_intervals(), n, sample_rate_input_intervals,
            [&](pos_t i){return idx_bwd.M_LF().p(i);});

        time = log_runtime(time);

        if constexpr (support != _count) {
            if constexpr (support == _locate_move) {
                std::cout << "building S_MPhi_p" << std::flush;

                _S_MPhi_p = build_sampling(
                    idx_fwd.M_Phi().num_intervals(), n, sample_rate_input_intervals,
                    [&](pos_t i){return idx_fwd.M_Phi().p(i);});

                time = log_runtime(time);
                std::cout << "building S_MPhi^{-1}_p" << std::flush;

                _S_MPhi_m1_p = build_sampling(
                    idx_fwd.M_Phi_m1().num_intervals(), n, sample_rate_input_intervals,
                    [&](pos_t i){return idx_fwd.M_Phi_m1().p(i);});

                time = log_runtime(time);
            }
            
            std::vector<pos_t> SA_fwd = idx_fwd.SA_range();

            time = now();
            std::cout << "building H" << std::flush;

            enum sample_type {
                _run_start,
                _run_end,
                _run_start_and_end
            };

            struct entry_t {
                pos_t x;
                sample_type type;
            };

            pos_t r_bwd = idx_bwd.M_LF().num_intervals();
            emhash6::HashMap<pos_t, entry_t, std::identity> H;
            H.reserve(2 * r_bwd);

            for (pos_t i = 0; i < r_bwd; i++) {
                if (idx_bwd.M_LF().p(i + 1) - idx_bwd.M_LF().p(i) == 1) {
                    H.emplace((n - idx_bwd.SA_s(i)) - 1,  entry_t {.x = i, .type = _run_start_and_end });
                } else {
                    H.emplace((n - idx_bwd.SA_s(i)) - 1,  entry_t {.x = i, .type = _run_start });
                    H.emplace((n - idx_bwd.SA_s_(i)) - 1, entry_t {.x = i, .type = _run_end   });
                }
            }
                
            time = log_runtime(time);
            std::cout << "building SA^{-1}_sR and SA^{-1}_eR" << std::flush;

            _SA_sR_m1 = interleaved_byte_aligned_vectors<pos_t, pos_t>({ byte_width(n) });
            _SA_eR_m1 = interleaved_byte_aligned_vectors<pos_t, pos_t>({ byte_width(n) });

            _SA_sR_m1.resize_no_init(r_bwd);
            _SA_eR_m1.resize_no_init(r_bwd);

            #pragma omp parallel for num_threads(p)
            for (pos_t i = 0; i < n; i++) {
                auto it = H.find(SA_fwd[i]);

                if (it != H.end()) {
                    if (it->second.type == _run_start || it->second.type == _run_start_and_end) [[likely]] {
                        _SA_sR_m1.template set_parallel<0, pos_t>(it->second.x, i);
                    }

                    if (it->second.type == _run_end || it->second.type == _run_start_and_end) [[likely]] {
                        _SA_eR_m1.template set_parallel<0, pos_t>(it->second.x, i);
                    }
                }
            }

            time = log_runtime(time);
        }

        std::cout << std::endl;
        std::cout << "overall construction time: " << format_time(time_diff_ns(time_start, now())) << std::endl;
        log_data_structure_sizes();
    }

    void reverse(std::string& T, uint16_t p, bool in_memory) {
        auto time = now();
        std::cout << "reversing T" << std::flush;
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

        log_runtime(time);
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
    static pos_t input_interval(pos_t i, const sd_array<pos_t>& sd_arr, fnc_t interval_start)
    {
        pos_t x = (sd_arr.rank_1(i) - 1) * sample_rate_input_intervals;

        while (i >= interval_start(x + 1)) {
            x++;
        }

        return x;
    }

    /**
     * @brief constructs a move_r index from an input file
     * @param input_file input file
     * @param params construction parameters
     */
    move_rb(std::ifstream& input_file, move_r_params params = {})
        requires(str_input)
    {
        
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
            std::cout << "S_MPhi_m1_p: " << format_size(_S_MPhi_m1_p.size_in_bytes()) << std::endl;
            std::cout << "S_MPhi_p: " << format_size(_S_MPhi_p.size_in_bytes()) << std::endl;
        }

        std::cout << std::endl << "Backward Index: " << std::endl;
        idx_bwd.log_data_structure_sizes(false);
        std::cout << "S_MLF_p: " << format_size(_S_MLF_p_bwd.size_in_bytes()) << std::endl;
        std::cout << "SA_^{-1}_sR: " << format_size(_SA_sR_m1.size_in_bytes()) << std::endl;
        std::cout << "SA_^{-1}_eR: " << format_size(_SA_eR_m1.size_in_bytes()) << std::endl;
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
    inline const move_r<support, sym_t, pos_t>& forward_index() const
    {
        return idx_fwd;
    }

    /**
     * @brief returns a reference to the move_r index of the backward input
     * @return index the backward input
     */
    inline const move_r<support, sym_t, pos_t>& backward_index() const
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

        bool v_b, v_e; // true <=> SA[b/e] can be computed with SA-Samples
        bool v_b_R, v_e_R; // true <=> n - SA^R[b_R/e_R] + 1 can be computed with SA^R-Samples
        pos_t i_b, i_e; // numbers of iterations since b/e has been altered by a rank-select query
        pos_t i_b_R, i_e_R; // numbers of iterations since b_R/e_R has been altered by a rank-select query
        pos_t x_b, x_e; // indices of the (sub-)runs in L', to whose beginning/end b/e has been set i_b/i_e iterations ago
        pos_t x_b_R, x_e_R; // indices of the (sub-)runs in L', to whose beginning/end b_R/e_R has been set i_b_R/i_e_R iterations ago

        // ############################# VARIABLES FOR THE LOCATE PHASE #############################
        
        direction locate_dir = _NONE; // current locate direction
        pos_t occ_rem; // number of remaining occurrences to locate
        pos_t c; // initial position in the suffix array
        pos_t SA_c; // initial suffix SA[c] in the suffix array interval
        pos_t i; // current position in the suffix array interval
        pos_t SA_i; // current suffix SA[i] in the suffix array interval
        pos_t s_; // index of the input inteval of M_Phi/M_Phi^{-1} containing SA_i
        pos_t x_p, x_lp, x_cp, x_r, s_np; // variables for decoding the rlzsa

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
            last_dir = _NONE;
            m = 0;

            b = 0;
            e = idx->idx_fwd.input_size();
            b_R = 0;
            e_R = e;
            b_ = 0;
            e_ = idx->idx_fwd.M_LF().num_intervals() - 1;
            b_R_ = 0;
            e_R_ = idx->idx_bwd.M_LF().num_intervals() - 1;

            v_b = true;
            v_e = true;
            v_b_R = true;
            v_e_R = true;
            i_b = 0;
            i_e = 0;
            i_b_R = 0;
            i_e_R = 0;
            x_b = 0;
            x_e = e_ - 1;
            x_b_R = 0;
            x_e_R = e_R_ - 1;
            
            i = 0;
            c = 0;
            SA_i = 0;
            SA_c = 0;
            s_ = 0;
            x_p = 0;
            x_lp = 0;
            x_cp = 0;
            x_r = 0;
            s_np = 0;
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
            return e >= b ? e - b + 1 : 0;
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
            if constexpr (dir == _LEFT) {
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
            if (i_sym == 0) return false;

            pos_t b_tmp = b;
            pos_t e_tmp = e;
            pos_t b_R_tmp = b_R;
            pos_t e_R_tmp = e_R;
            pos_t b__tmp = b_;
            pos_t e__tmp = e_;
            pos_t b_R__tmp = b_R_;
            pos_t e_R__tmp = e_R_;

            bool v_b_tmp = v_b;
            bool v_e_tmp = v_e;
            bool v_b_R_tmp = v_b_R;
            bool v_e_R_tmp = v_e_R;
            pos_t i_l_tmp = i_b;
            pos_t i_r_tmp = i_e;
            pos_t i_l_R_tmp = i_b_R;
            pos_t i_r_R_tmp = i_e_R;
            pos_t x_b_tmp = x_b;
            pos_t x_e_tmp = x_e;
            pos_t x_b_R_tmp = x_b_R;
            pos_t x_e_R_tmp = x_e_R;

            bool result;

            if (dir == _LEFT) {
                result = extend<_LEFT>(
                    i_sym, last_dir,
                    v_b, v_e, v_b_R, v_e_R,
                    b, e, b_R, e_R,
                    b_, e_, b_R_, e_R_,
                    i_b, i_e, i_b_R, i_e_R,
                    x_b, x_e, x_b_R, x_e_R);
            } else {
                result = extend<_RIGHT>(
                    i_sym, last_dir,
                    v_b_R, v_e_R, v_b, v_e,
                    b_R, e_R, b, e,
                    b_R_, e_R_, b_, e_,
                    i_b_R, i_e_R, i_b, i_e,
                    x_b_R, x_e_R, x_b, x_e);
            }

            if (result) {
                m++;
                i = b;
                occ_rem = e - b + 1;
                locate_dir = _NONE;
                last_dir = dir;
            } else {
                b = b_tmp;
                e = e_tmp;
                b_R = b_R_tmp;
                e_R = e_R_tmp;
                b_ = b__tmp;
                e_ = e__tmp;
                b_R_ = b_R__tmp;
                e_R_ = e_R__tmp;

                v_b = v_b_tmp;
                v_e = v_e_tmp;
                v_b_R = v_b_R_tmp;
                v_e_R = v_e_R_tmp;
                i_b = i_l_tmp;
                i_e = i_r_tmp;
                i_b_R = i_l_R_tmp;
                i_e_R = i_r_R_tmp;
                x_b = x_b_tmp;
                x_e = x_e_tmp;
                x_b_R = x_b_R_tmp;
                x_e_R = x_e_R_tmp;
            }

            return result;
        }

        /**
         * @brief extends the currently matched pattern with sym; if the extended pattern occurs in the input, true is
         * returned and the query context is adjusted to store the information for the extended pattern; else,
         * false is returned and the query context is not modified; let S be a string, then S(_LEFT) = S; S(_RIGHT) = S^R
         * and !LEFT = RIGHT = !!LEFT
         * @tparam dir extend direction
         * @param sym next symbol to match
         * @param b Left interval limit of the suffix array interval of P(dir) in T(dir).
         * @param e Right interval limit of the suffix array interval of P(dir) in T(dir).
         * @param b_ index of the input interval in M_LF of T(dir) containing b.
         * @param e_ index of the input interval in M_LF of T(dir) containing e.
         * @param b_r Left interval limit of the suffix array interval of P(!dir) in T(!dir).
         * @param e_r Right interval limit of the suffix array interval of P(!dir) in T(!dir).
         * @param b_r_ index of the input interval in M_LF of T(!dir) containing b_r.
         * @param e_r_ index of the input interval in M_LF of T(!dir) containing e_r.
         * @param hat_b_ap_y \hat{b}'_y
         * @param y y
         * @param hat_e_ap_z \hat{e}'_z
         * @param z z
         * @param hat_b_r_ap_y \hat{b_r}'_y
         * @param y_r y_r
         * @param hat_e_r_ap_z \hat{e_r}'_y
         * @param z_r z_r
         * @return whether the extended pattern occurs in the input
         */
        template <direction dir>
        bool extend(
            i_sym_t i_sym,
            direction last_dir,
            bool& v_b, bool& v_e,
            bool& v_b_R, bool& v_e_R,
            pos_t& b, pos_t& e,
            pos_t& b_R, pos_t& e_R,
            pos_t& b_, pos_t& e_,
            pos_t& b_R_, pos_t& e_R_,
            pos_t& i_b, pos_t& i_e,
            pos_t& i_b_R, pos_t& i_e_R,
            pos_t& x_b, pos_t& x_e,
            pos_t& x_b_R, pos_t& x_e_R ) const
        {
            if (last_dir != _NONE && dir != last_dir) {
                b_ = input_interval(b, idx->S_MLF_p<dir>(), [&](pos_t x){return idx->index<dir>().M_LF().p(x);});

                if (b == idx->index<dir>().M_LF().p(b_)) {
                    v_b = true;
                    i_b = 0;
                    x_b = b_;
                } else {
                    v_b = false;
                }

                e_ = input_interval(e, idx->S_MLF_p<dir>(), [&](pos_t x){return idx->index<dir>().M_LF().p(x);});

                if (e == idx->index<dir>().M_LF().p(e_)) {
                    v_e = true;
                    i_e = 0;
                    x_e = e_;
                } else {
                    v_e = false;
                }
            }
            
            int64_t next[256];
            int64_t prev[256];

            std::fill_n(next, 256, LONG_MAX);
            std::fill_n(prev, 256, LONG_MIN);

            pos_t b__old = b_;
            pos_t e__old = e_;

            pos_t b_R_old = b_R;
            pos_t e_R_old = e_R;

            pos_t blk_size = idx->index<dir>().L_block_size();

            {
                pos_t blk = div_ceil<pos_t>(b_, blk_size);
                int64_t beg = b_;
                int64_t end = std::min<pos_t>(blk * blk_size, e_);

                for (int64_t i = end; i >= beg; i--) {
                    next[idx->index<dir>().L_(i)] = i;
                }

                if (end != e_) {
                    for (i_sym_t i = 0; i <= i_sym; i++) {
                        if (next[i] == LONG_MAX) {
                            next[i] = idx->index<dir>().L_next(blk, i);
                        }
                    }
                }
            }

            {
                pos_t blk = e_ / blk_size;
                int64_t beg = std::max<pos_t>(blk * blk_size, b_);
                int64_t end = e_;

                if (beg != b_) {
                    for (i_sym_t i = 0; i <= i_sym; i++) {
                        if (prev[i] == LONG_MIN) {
                            prev[i] = idx->index<dir>().L_prev(blk, i);
                        }
                    }
                }

                for (int64_t i = beg; i <= end; i++) {
                    prev[idx->index<dir>().L_(i)] = i;
                }
            }

            b_ = next[i_sym];
            e_ = prev[i_sym];

            if (b_ > e_) [[unlikely]] {
                return false;
            }

            if (b_ != b__old) {
                b = idx->index<dir>().M_LF().p(b_);
                v_b = true;
                x_b = b_;
                i_b = 0;
            } else if (v_b) {
                i_b++;
            }

            if (e_ != e__old) {
                e = idx->index<dir>().M_LF().p(e_ + 1) - 1;
                v_e = true;
                x_e = e_;
                i_e = 0;
            } else if (v_e) {
                i_e++;
            }

            if (b_ == e_) {
                if (b == e) {
                    idx->index<dir>().M_LF().move(b, b_);
                    e = b;
                    e_ = b_;
                } else {
                    pos_t diff_eb = e - b;
                    idx->index<dir>().M_LF().move(b, b_);
                    e = b + diff_eb;
                    e_ = b_;

                    while (e >= idx->index<dir>().M_LF().p(e_ + 1)) {
                        e_++;
                    }
                }
            } else {
                idx->index<dir>().M_LF().move(b, b_);
                idx->index<dir>().M_LF().move(e, e_);
            }

            for (i_sym_t i = 0; i < i_sym; i++) {
                pos_t x = next[i];
                pos_t y = prev[i];

                if (x <= y) {
                    b_R += (idx->index<dir>().M_LF().q(y) -     idx->index<dir>().M_LF().q(x)) +
                           (idx->index<dir>().M_LF().p(y + 1) - idx->index<dir>().M_LF().p(y));
                }
            }

            e_R = b_R + (e - b);

            if (b_R != b_R_old) {
                v_b_R = false;
                i_b_R = 0;
                x_b_R = 0;
                b_R_ = 0;
            } else if (v_b_R) {
                i_b_R++;
            }

            if (e_R != e_R_old) {
                v_e_R = false;
                i_e_R = 0;
                x_e_R = 0;
                e_R_ = 0;
            } else if (v_e_R) {
                i_e_R++;
            }

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
            return extend<_LEFT>(sym);
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
            return extend<_RIGHT>(sym);
        }

        /**
         * @brief reports the next occurrence of the currently matched pattern
         * @return next occurrence
         */
        inline pos_t next_occ()
        {
            if (locate_dir == _NONE) {
                if (v_b) {
                    c = b;
                    i = b;
                    SA_i = idx->idx_fwd.SA_s(x_b) - i_b;
                } else if (v_e) {
                    c = e;
                    i = e;
                    SA_i = idx->idx_fwd.SA_s_(x_e) - i_e;
                } else if (v_b_R) {
                    SA_i = (n - idx->idx_bwd.SA_s(x_b_R) + 1) + i_b_R - m;
                    i = idx->_SA_sR_m1(x_b_R);
                    pos_t i_ = input_interval(i, idx->_S_MLF_p_fwd, [&](pos_t x){return idx->idx_fwd.M_LF().p(x);});

                    for (pos_t j = 0; j < m - i_b_R; j++) {
                        idx->idx_fwd.M_LF().move(i, i_);
                    }

                    c = i;
                } else {
                    SA_i = (n - idx->idx_bwd.SA_s_(x_e_R) + 1) + i_e_R - m;
                    i = idx->_SA_sR_m1(x_e_R);
                    pos_t i_ = input_interval(i, idx->_S_MLF_p_fwd, [&](pos_t x){return idx->idx_fwd.M_LF().p(x);});

                    for (pos_t j = 0; j < m - i_e_R; j++) {
                        idx->idx_fwd.M_LF().move(i, i_);
                    }

                    c = i;
                }

                SA_c = SA_i;
                locate_dir = b == e ? _NONE : (i > b ? _LEFT : _RIGHT);

                if (locate_dir == _LEFT) {
                    s_ = input_interval(i, idx->_S_MPhi_p, [&](pos_t x){return idx->idx_fwd.M_Phi().p(x);});
                } else if (locate_dir == _RIGHT) {
                    s_ = input_interval(i, idx->_S_MPhi_m1_p, [&](pos_t x){return idx->idx_fwd.M_Phi_m1().p(x);});
                }
            } else if (locate_dir == _LEFT) {
                idx->idx_fwd.M_Phi().move(SA_i, s_);

                if (i > b) [[likely]] {
                    i--;
                } else {
                    i = c;
                    SA_i = SA_c;
                    s_ = input_interval(i, idx->_S_MPhi_m1_p, [&](pos_t x){return idx->idx_fwd.M_Phi_m1().p(x);});
                    locate_dir = _RIGHT;
                }
            } else if (locate_dir == _RIGHT) {
                idx->idx_fwd.M_Phi_m1().move(SA_i, s_);
                
                if (i < e) [[likely]] {
                    i++;
                } else {
                    locate_dir = _NONE;
                }
            }

            return SA_i;
        }

        /**
         * @brief locates the remaining (not yet reported) occurrences of the currently matched pattern
         * @param Occ vector to append the occurrences to
         */
        inline void locate(std::vector<pos_t>& Occ)
        {
        }

        /**
         * @brief locates the remaining (not yet reported) occurrences of the currently matched pattern
         * @return vector containing the occurrences
         */
        std::vector<pos_t> locate()
        {
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