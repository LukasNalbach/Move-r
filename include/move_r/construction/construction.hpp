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

#include <sdsl/int_vector_buffer.hpp>
#include <gtl/btree.hpp>
#include <hash_table5.hpp>
#include <move_r/move_r.hpp>
#include <misc/strings.hpp>
#include <lzendsa/lzendsa_construction.hpp>

enum rlbwt_build_mode {
    _sa, // the BWT is read by L[i] = T[(SA[i]-1) mod n]
    _bwt, // the BWT is read by accessing L[i]
    _bwt_file // the BWT is read from files output by Big-BWT
};

template <move_r_support support, typename sym_t, typename pos_t>
class move_r<support, sym_t, pos_t>::construction {
public:
    construction() = delete;
    construction(construction&&) = delete;
    construction(const construction&) = delete;
    construction& operator=(construction&&) = delete;
    construction& operator=(const construction&) = delete;
    ~construction() { }

    // ############################# MISC VARIABLES #############################

    std::string T_str_tmp;
    std::vector<sym_t> T_vec_tmp;
    std::string L_tmp;
    std::vector<int32_t> SA_32_tmp;
    std::vector<int64_t> SA_64_tmp;
    uint16_t p = 1; // the number of threads to use
    move_r_construction_mode mode = _suffix_array;
    /* the number of threads to use during the construction of the L,C and I_LF,I_Phi^{-1},L' and SA_s; when building move-r for an
     * integer alphabet, this needs (p+1)*sigma = O(p*n) words of space, so we limit the number of threads to p' to 1 to ensure
     * that we use O(sigma) = O(n) words of space; */
    uint16_t p_ = 1;
    bool build_sa_and_l = false; // controls whether the index should be built from the suffix array and the bwt
    bool delete_T = false; // controls whether T should be deleted when not needed anymore
    bool log = false; // controls, whether to print log messages
    std::ostream* mf_idx = nullptr; // file to write measurement data of the index construction to
    std::ostream* mf_mds = nullptr; // file to write measurement data of the move data structure construction to
    std::string name_text_file = ""; // name of the text file (only for measurement output)
    std::string prefix_tmp_files = ""; // prefix of temporary files
    std::chrono::steady_clock::time_point time; // time of the start of the last build phase
    std::chrono::steady_clock::time_point time_start; // time of the start of the whole build phase
    uint64_t baseline_mem_usage = 0; // memory allocation at the start of the construction
    uint64_t bigbwt_peak_mem_usage = 0; // peak memory usage during the execution of Big-BWT
    uint8_t min_valid_char = 0; // the minimum valid character that is allowed to occur in T
    void* sa_vector; // pointer to the suffix array (used only for bidirectional forward indexes)
    std::string* sa_file_name; // name of the suffix array file (used only for bidirectional forward indexes)

    // ############################# INDEX VARIABLES #############################

    pos_t n = 0; // the length of T
    pos_t r = 0; // r, the number of runs in L
    pos_t r_ = 0; // r', the number of input/output intervals in M_LF
    pos_t r__ = 0; // r'', the number of input/output intervals in M_Phi^{-1}
    pos_t r___ = 0; // r''', the number of input/output intervals in M_Phi

    // ############################# CONSTRUCTION DATA STRUCTURE VARIABLES #############################

    /** the string containing T */
    std::string& T_str;
    /** the vector containing T */
    std::vector<sym_t>& T_vec;
    /** The move-r index to construct */
    move_r<support, sym_t, pos_t>& idx;
    /** [0..n-1] The suffix array (32-bit) */
    std::vector<int32_t>& SA_32;
    /** [0..n-1] The suffix array (64-bit) */
    std::vector<int64_t>& SA_64;
    /** [0..n-1] The BWT */
    std::string& L;
    /** [0..p-1] file buffers (of each thread) for reading the BWT file output by bigbwt */
    std::vector<sdsl::int_vector_buffer<8>> BWT_file_bufs;
    /** [0..p-1] vectors that contain the RLBWT concatenated */
    std::vector<interleaved_byte_aligned_vectors<uint32_t, uint32_t>> RLBWT;
    /** [0..p] n_p[0] < n_p[1] < ... < n_p[p] = n; n_p[i] = start position of thread i's section in L and SA */
    std::vector<pos_t> n_p;
    /** [0..p] r_p[0] < r_p[1] < ... < r_p[p] = r; r_p[i] = index of the first run in L starting in
     * [n_p[i]..n_p[i+1]-1]; there is a run starting at n_p[i] */
    std::vector<pos_t> r_p;
    /** [0..p][0..255] see the code to see how this variable is used */
    std::vector<std::vector<pos_t>> C;
    /** The disjoint interval sequence for LF */
    std::vector<std::pair<pos_t, pos_t>> I_LF;
    /** The disjoint interval sequence for Phi */
    std::vector<std::pair<pos_t, pos_t>> I_Phi;
    /** The disjoint interval sequence for Phi^{-1} */
    std::vector<std::pair<pos_t, pos_t>> I_Phi_m1;
    /** [0..r'-1] SA_s[x] = SA[M_LF.p[x]]; if the starting position of the
     * x-th input interval of M_LF is not starting position of a BWT run, then SA_s[x] = n */
    std::vector<pos_t> SA_s;
    /** [0..r'-1] Permutation storing the order of the values in SA_s */
    std::vector<pos_t> pi_;
    /** [0..r''-1] Permutation storing the order of the output interval starting positions of M_Phi^{-1} */
    std::vector<pos_t> pi_mphi;
    /** [0..p-1] file buffers (of each thread) for reading the suffix array file output by bigbwt */
    std::vector<sdsl::int_vector_buffer<40>> SA_file_bufs;

    /**
     * @brief builds the rlzsa; the construction is embedded in the rlzsa data structure (see
     *        rlzsa_opt/construction.hpp) and does not depend on this construction object: it receives the data
     *        it needs through parameters (a function SAd(i_p, i) returning SA^d[i] for the in-memory case, a
     *        reference to I_Phi^{-1}, n, r, and the optional rlzsa_construction_params). For the bigbwt case
     *        it reads SA^d itself from the suffix-array file (sa_file_path). It builds idx._rlzsa.
     * @tparam bigbwt true <=> read suffix array values from the file output by bigbwt
     * @tparam sa_sint_t signed integer type to use for the suffix array entries
     */
    template <bool bigbwt, typename sa_sint_t>
    void construct_rlzsa()
    {
        // SA^d[i] = SA[i] - SA[i-1] + n (with SA[-1] := -1, i.e. SA^d[0] = n - 1) read from the in-memory
        // suffix array; only used for the in-memory case (for bigbwt the construction reads the SA file)
        auto SAd = [this](uint16_t /*i_p*/, pos_t i) -> uint64_t {
            if (i == 0) [[unlikely]] return uint64_t { n } - 1;
            return (get_sa<sa_sint_t>()[i] + n) - get_sa<sa_sint_t>()[i - 1];
        };

        typename rlzsa_opt<pos_t>::template construction<decltype(SAd)> rlzsa_construction(
            idx._rlzsa, SAd, I_Phi_m1, n, r,
            { .num_threads = p, .mode = mode, .log = log,
              .sa_file_path = bigbwt ? (prefix_tmp_files + ".sa") : std::string() });
        rlzsa_construction.build();
    }

    // ############################# COMMON MISC METHODS #############################

    /**
     * @brief returns T at index i interpreted as type
     * @tparam type type to reinterpret the symbol at index i as
     * @param i an index into T
     * @return T at index i
     */
    template <typename type>
    inline type& T(pos_t i)
    {
        if constexpr (str_input) {
            return *reinterpret_cast<type*>(&T_str[i]);
        } else {
            return *reinterpret_cast<type*>(&T_vec[i]);
        }
    }

    /**
     * @brief returns the suffix array
     * @tparam sa_sint_t suffix array signed integer type
     * @return the suffix array
     */
    template <typename sa_sint_t>
    constexpr std::vector<sa_sint_t>& get_sa()
    {
        if constexpr (std::is_same_v<sa_sint_t, int32_t>) {
            return SA_32;
        } else {
            return SA_64;
        }
    }

    /**
     * @brief sets the run length of the i-th BWT run in thread i_p's section to len
     * @param i_p [0..p-1] thread index
     * @param i [0..r-1] run index
     * @param len run length
     */
    inline void set_run_len(uint16_t i_p, pos_t i, pos_t len)
    {
        RLBWT[i_p].set_unsafe<1, uint32_t>(i, len);
    }

    /**
     * @brief returns the length of the i-th BWT run in thread i_p's section
     * @param i_p [0..p-1] thread index
     * @param i [0..r-1] run index
     * @return run length
     */
    inline pos_t run_len(uint16_t i_p, pos_t i)
    {
        return RLBWT[i_p].get_unsafe<1, uint32_t>(i);
    }

    /**
     * @brief sets the character of the i-th BWT run in thread i_p's section to c
     * @param i_p [0..p-1] thread index
     * @param i [0..r-1] run index
     * @param c character
     */
    inline void set_run_sym(uint16_t i_p, pos_t i, i_sym_t c)
    {
        if constexpr (byte_alphabet) {
            RLBWT[i_p].set_unsafe<0, i_sym_t>(i, c);
        } else {
            RLBWT[i_p].set_parallel<0, i_sym_t>(i, c);
        }
    }

    /**
     * @brief returns the character of the i-th BWT run in thread i_p's section
     * @param i_p [0..p-1] thread index
     * @param i [0..r-1] run index
     * @return character
     */
    inline i_sym_t run_sym(uint16_t i_p, pos_t i)
    {
        if constexpr (byte_alphabet) {
            return RLBWT[i_p].get_unsafe<0, i_sym_t>(i);
        } else {
            return RLBWT[i_p].get<0, i_sym_t>(i);
        }
    }

    /**
     * @brief adds a bwt run of length len with the symbol sym to the end of thread i_p's section of the RLBWT
     * @param i_p [0..p-1] thread index
     * @param sym run symbol
     * @param len run length
     */
    inline void add_run(uint16_t i_p, i_sym_t sym, pos_t len)
    {
        if constexpr (byte_alphabet) {
            RLBWT[i_p].emplace_back_unsafe<i_sym_t, uint32_t>({ sym, len });
        } else {
            RLBWT[i_p].emplace_back<i_sym_t>({ sym });
            set_run_len(i_p, RLBWT[i_p].size() - 1, len);
        }
    }

    // ############################# SUFFIX ARRAY ACCESS #############################

    /**
     * @brief returns SA[i]
     * @tparam bigbwt true <=> read SA[i] from SA_file
     * @tparam sa_sint_t suffix array signed integer type
     * @param i_p [0..p-1] thread index
     * @param i [0..n-1] suffix array index
     * @return SA[i]
     */
    template <bool bigbwt, typename sa_sint_t>
    inline pos_t SA(uint16_t i_p, pos_t i)
    {
        if constexpr (bigbwt) {
            if (i == 0) [[unlikely]] return n - 1;
            return (pos_t)(SA_file_bufs[i_p][i - 1]);
        } else {
            return get_sa<sa_sint_t>()[i];
        }
    }

    // ############################# COMMON MISC METHODS #############################

    /**
     * @brief sets some variables and logs
     */
    void prepare_phase_1()
    {
        time = now();
        time_start = time;
        omp_set_num_threads(p);
        prefix_tmp_files = "move-r_" + random_alphanumeric_string(10);
        baseline_mem_usage = malloc_count_current();
        if (log) malloc_count_reset_peak();
    }

    /**
     * @brief sets some variables
     */
    void prepare_phase_2()
    {
        if (p > 1 && 1000 * p > n) {
            p = std::max<pos_t>(1, n / 1000);
            p_ = std::min<uint16_t>(p_, p);
            log_message(log, "warning: p > n/1000, setting p to n/1000 ~ " + std::to_string(p) + "\n");
        }
    }

    /**
     * @brief starts a construction phase: resets the phase timer and prints msg (if logging)
     * @param msg message describing the phase
     */
    void log_phase_start(const std::string& msg)
    {
        ::log_phase_start(log, time, msg);
    }

    /**
     * @brief ends a construction phase: writes the elapsed time to mf_idx under the key mf_key
     *        (if mf_key is non-empty) and logs the runtime (if logging)
     * @param mf_key measurement-file key for the elapsed time (empty => not written)
     */
    void log_phase_end(const std::string& mf_key = "")
    {
        ::log_phase_end(log, time, mf_idx, mf_key);
    }

    /**
     * @brief ends a move-data-structure construction phase: closes the mf_mds RESULT line,
     *        writes the elapsed time and the resulting interval count to mf_idx, and prints a newline.
     *        The move data structure itself logs its own progress, so no runtime is logged here.
     * @param time_key measurement-file key for the elapsed time
     * @param count_key measurement-file key for the interval count
     * @param count number of intervals in the built move data structure
     */
    void log_mds_phase_end(const std::string& time_key, const std::string& count_key, pos_t count)
    {
        if (log) {
            if (mf_mds != nullptr) *mf_mds << std::endl;
            if (mf_idx != nullptr)
                *mf_idx << " " << time_key << "=" << time_diff_ns(time, now()) << " " << count_key << "=" << count;
            std::cout << std::endl;
        }
    }

    /**
     * @brief logs the current memory usage
     */
    void log_mem_usage()
    {
        std::cout << "current memory allocation: "
                  << format_size(malloc_count_current() - baseline_mem_usage)
                  << std::endl;
    }

    /**
     * @brief logs the peak memory usage until now
     */
    void log_peak_mem_usage()
    {
        std::cout << "peak memory allocation until now: "
                  << format_size(std::max(malloc_count_peak() - baseline_mem_usage, bigbwt_peak_mem_usage))
                  << std::endl;
    }

    /**
     * @brief logs a construction summary
     */
    void log_finished()
    {
        if (!log) return;
        uint64_t time_construction = time_diff_ns(time_start, now());
        uint64_t peak_mem_usage = std::max(malloc_count_peak() - baseline_mem_usage, bigbwt_peak_mem_usage);

        std::cout << std::endl;
        std::cout << "construction time: " << format_time(time_construction) << std::endl;
        std::cout << "peak memory usage: " << format_size(peak_mem_usage) << std::endl;
        if (!is_bidirectional) {
            std::cout << "construction throughput: " << format_construction_throughput(n, time_construction) << std::endl;
            idx.log_data_structure_sizes();
        }

        if (mf_idx != nullptr) {
            *mf_idx << " time_construction=" << time_construction;
            *mf_idx << " construction_throughput_mb_per_s=" << construction_throughput_mb_per_s(n, time_construction);
            *mf_idx << " peak_mem_usage=" << peak_mem_usage;
            idx.log_data_structure_sizes(*mf_idx);
        }
    }

    /**
     * @brief logs statistics of T
     */
    void log_statistics()
    {
        if (!log) return;
        double n_r = std::round(100.0 * (n / (double)r)) / 100.0;
        if (mf_idx != nullptr) {
            *mf_idx << " n=" << n;
            *mf_idx << " sigma=" << std::to_string(idx.sigma);
            *mf_idx << " r=" << r;
        }
        std::cout << "n = " << n << ", sigma = " << std::to_string(idx.sigma) << ", r = " << r << ", n/r = " << n_r << std::endl;
    }

    // ############################# CONSTRUCTORS #############################

    void read_parameters(move_r_params& params)
    {
        idx.sigma = params.alphabet_size;
        this->p = params.num_threads;
        this->mode = params.mode;
        idx.a = params.a;
        this->log = params.log;
        this->mf_idx = params.mf_idx;
        this->mf_mds = params.mf_mds;
        this->name_text_file = params.name_text_file;
        this->sa_vector = params.sa_vector;
        this->sa_file_name = params.sa_file_name;
    }

    /**
     * @brief constructs a move_r index of the string input
     * @param index The move-r index to construct
     * @param T the string containing T
     * @param delete_T controls whether T should be deleted once it is not needed anymore
     * @param params construction parameters
     */
    construction(move_r<support, sym_t, pos_t>& index, std::string& T, bool delete_T, move_r_params params)
        requires(str_input)
        : T_str(T)
        , T_vec(T_vec_tmp)
        , idx(index)
        , SA_32(SA_32_tmp)
        , SA_64(SA_64_tmp)
        , L(L_tmp)
    {
        this->delete_T = delete_T;
        read_parameters(params);
        prepare_phase_1();
        
        if (params.file_input) {
            n = std::filesystem::file_size(T) + 1;
        } else {
            n = T.size() + 1;
        }

        idx.n = n;

        if (mode == _suffix_array || mode == _suffix_array_space) {
            min_valid_char = 1;

            if (params.file_input) {
                read_t_from_file(T);
            } else {
                T.push_back(0);
            }
            
            if (idx.sigma == 0) preprocess_t(true, false, T);
            construct_from_sa();

            if (!delete_T && !params.file_input) {
                T.resize(n - 1);

                if (!is_bidirectional) {
                    unmap_t(true);
                }
            }
        } else {
            min_valid_char = 3;
            if (params.file_input) prefix_tmp_files = T;
            if (idx.sigma == 0) preprocess_t(!params.file_input, true, T);
            if (!params.file_input) store_t_in_file();
            construct_from_bigbwt(!params.file_input);
            if (!is_bidirectional && params.file_input) unmap_t(false);
        }

        log_finished();
        
        if (params.peak_memory_usage != nullptr) {
            *params.peak_memory_usage = std::max({
                *params.peak_memory_usage,
                bigbwt_peak_mem_usage,
                malloc_count_peak() - baseline_mem_usage});
        }
    }

    /**
     * @brief constructs a move_r index of the vector input
     * @param index The move-r index to construct
     * @param T the vector containing T
     * @param alphabet_size alphabet size of the input vector (= maximum value in the input vector)
     * @param delete_T controls whether T should be deleted once it is not needed anymore
     * @param params construction parameters
     */
    construction(move_r<support, sym_t, pos_t>& index, std::vector<sym_t>& T, bool delete_T, move_r_params params)
        requires(int_input)
        : T_str(T_str_tmp)
        , T_vec(T)
        , idx(index)
        , SA_32(SA_32_tmp)
        , SA_64(SA_64_tmp)
        , L(L_tmp)
    {
        this->delete_T = delete_T;
        read_parameters(params);
        prepare_phase_1();
        // reserve 0 for the suffix-array sentinel: for a byte alphabet preprocess_t remaps the symbols starting
        // at min_valid_char, so this keeps the smallest input symbol distinct from the sentinel (for an integer
        // alphabet the symbols are remapped to [1, sigma) regardless, so min_valid_char is not used)
        min_valid_char = 1;
        T.push_back(0);
        n = T.size();
        idx.n = n;
        if (idx.sigma == 0) preprocess_t(true, false, "");
        construct_from_sa();

        if (!delete_T) {
            T.resize(n - 1);
            if (idx.symbols_remapped && !delete_T) unmap_t(true);
        }

        log_finished();
        
        if (params.peak_memory_usage != nullptr) {
            *params.peak_memory_usage = std::max({
                *params.peak_memory_usage,
                bigbwt_peak_mem_usage,
                malloc_count_peak() - baseline_mem_usage});
        }
    }

    /**
     * @brief constructs a move_r index from an input file
     * @param index The move-r index to construct
     * @param suffix_array vector containing the suffix array of the input
     * @param bwt string containing the bwt of the input
     * @param params construction parameters
     */
    construction(move_r<support, sym_t, pos_t>& index, std::vector<int32_t>& suffix_array, std::string& bwt, move_r_params params)
        requires(str_input)
        : T_str(T_str_tmp)
        , T_vec(T_vec_tmp)
        , idx(index)
        , SA_32(suffix_array)
        , SA_64(SA_64_tmp)
        , L(bwt)
    {
        read_parameters(params);
        construct_from_sa_and_l<int32_t>();
    }

    /**
     * @brief constructs a move_r index from an input file
     * @param index The move-r index to construct
     * @param suffix_array vector containing the suffix array of the input
     * @param bwt string containing the bwt of the input
     * @param params construction parameters
     */
    construction(move_r<support, sym_t, pos_t>& index, std::vector<int64_t>& suffix_array, std::string& bwt, move_r_params params)
        requires(str_input)
        : T_str(T_str_tmp)
        , T_vec(T_vec_tmp)
        , idx(index)
        , SA_32(SA_32_tmp)
        , SA_64(suffix_array)
        , L(bwt)
    {
        read_parameters(params);
        construct_from_sa_and_l<int64_t>();
    }

    // ############################# CONSTRUCTION #############################

    /**
     * @brief constructs the index from a suffix array and a bwt
     * @tparam sa_sint_t signed integer type of the suffix array entries
     */
    template <typename sa_sint_t>
    void construct_from_sa_and_l()
    {
        build_sa_and_l = true;
        min_valid_char = 1;
        n = L.size();
        idx.n = n;

        prepare_phase_1();
        prepare_phase_2();
        build_rlbwt_c<_bwt, sa_sint_t>();
        log_statistics();

        build_ilf();
        build_mlf();

        if constexpr (supports_locate) {
            build_iphim1_sas_from_sa<false, sa_sint_t>();
            build_l__sas();
        } else {
            build_l__sas();
        }

        if constexpr (supports_multiple_locate) {
            if constexpr (has_locate_move) {
                if constexpr (support == _locate_move_bi_fwd) {
                    build_iphi();
                    sort_iphi();
                    build_mphi();
                }

                sort_iphim1();
                build_mphim1();
                if constexpr (support == _locate_move) build_saphim1();
                build_de();
            } else if constexpr (has_rlzsa) {
                construct_rlzsa<false, sa_sint_t>();
            }
        }

        build_l_prev_next();

        log_finished();
    }

    /**
     * @brief constructs the index in memory (uses libsais)
     */
    void construct_from_sa()
    {
        if constexpr (std::is_same_v<pos_t, uint64_t>) {
            construct_from_sa<int64_t>();
        } else if (n <= INT_MAX) {
            construct_from_sa<int32_t>();
        } else {
            construct_from_sa<int64_t>();
        }
    }

    /**
     * @brief constructs the index using the suffix array
     * @tparam sa_sint_t signed integer type to use for the suffix array entries
     */
    template <typename sa_sint_t>
    void construct_from_sa()
    {
        bool _space = mode == _suffix_array_space;

        prepare_phase_2();
        build_sa<sa_sint_t>();
        build_rlbwt_c<_sa, sa_sint_t>();
        log_statistics();

        build_ilf();
        if (_space) store_rlbwt();
        build_mlf();
        if (_space) load_rlbwt();

        if constexpr (supports_locate) {
            if (support == _locate_bi_bwd) {
                build_sas_from_sa<false, sa_sint_t>();
            } else {
                build_iphim1_sas_from_sa<false, sa_sint_t>();
            }

            build_l__sas();

            if constexpr (supports_multiple_locate) {
                if constexpr (has_locate_move) {
                    if (_space) store_sas();
                    if constexpr (byte_alphabet) build_l_prev_next();
                    else build_rsl_();
                    if (_space) {
                        if constexpr (byte_alphabet) store_l_prev_next();
                        else store_rsl_();
                    }
                    if (_space) store_mlf();
                    if constexpr (support == _locate_move_bi_fwd) {
                        build_iphi();
                        if (_space) store_iphim1();
                        sort_iphi();
                        build_mphi();
                        if (_space) store_mphi();
                        if (_space) load_iphim1();
                    }
                    sort_iphim1();
                    build_mphim1();
                    if (_space) load_sas();
                    if (!is_bidirectional) build_saphim1();
                    build_de();
                    if (_space) {
                        load_mlf();
                        if constexpr (support == _locate_move_bi_fwd) load_mphi();
                        if constexpr (byte_alphabet) load_l_prev_next();
                        else load_rsl_();
                    }
                } else if constexpr (has_rlzsa) {
                    if (_space && support != _locate_bi_bwd) {
                        store_sas_idx();
                        if (is_bidirectional) store_sas__idx();
                    }
                    if constexpr (byte_alphabet) build_l_prev_next();
                    else build_rsl_();
                    if (_space) {
                        if constexpr (byte_alphabet) store_l_prev_next();
                        else store_rsl_();
                    }
                    if (_space) store_mlf();

                    sort_iphim1();
                    construct_rlzsa<false, sa_sint_t>();

                    if (_space) {
                        load_sas_idx();
                        if (is_bidirectional) load_sas__idx();
                    }
                    if (_space) load_mlf();
                    if (_space) {
                        if constexpr (byte_alphabet) load_l_prev_next();
                        else load_rsl_();
                    }
                } else if constexpr (support == _locate_bi_bwd) {
                    build_l_prev_next();
                }

                if constexpr (support == _locate_move_bi_fwd || support == _locate_rlzsa_bi_fwd) {
                    *reinterpret_cast<std::vector<sa_sint_t>*>(sa_vector) = std::move(get_sa<sa_sint_t>());
                }
            } else {
                if constexpr (byte_alphabet) build_l_prev_next();
                else build_rsl_();
            }
        } else {
            build_l__sas();
            if constexpr (byte_alphabet) build_l_prev_next();
            else build_rsl_();
        }

        if constexpr (int_alphabet) {
            if (idx.symbols_remapped && mode == _suffix_array_space)
                load_mapintext();
        }
    }

    /**
     * @brief constructs the index from a file using Big-BWT
     * @param delete_T delete T file
     */
    void construct_from_bigbwt(bool delete_T)
    {
        prepare_phase_2();
        bigbwt(delete_T);
        build_rlbwt_c<_bwt_file, int32_t>();
        log_statistics();

        build_ilf();
        store_rlbwt();
        build_mlf();
        load_rlbwt();

        if constexpr (supports_locate) {
            if (is_bidirectional || has_rlzsa || p > 1) {
                if (support == _locate_bi_bwd) {
                    build_sas_from_sa<true, int32_t>();
                } else {
                    build_iphim1_sas_from_sa<true, int32_t>();
                }
            } else {
                read_iphim1_bigbwt();
            }

            build_l__sas();

            if constexpr (supports_multiple_locate) {
                if constexpr (has_locate_move) {
                    store_sas();
                    build_l_prev_next();
                    store_mlf();
                    store_l_prev_next();

                    if constexpr (support == _locate_move_bi_fwd) {
                        build_iphi();
                        store_iphim1();
                        sort_iphi();
                        build_mphi();
                        store_mphi();
                        load_iphim1();
                    }

                    sort_iphim1();
                    build_mphim1();
                    load_sas();
                    if constexpr (support == _locate_move) build_saphim1();
                    build_de();
                    load_mlf();
                    if constexpr (support == _locate_move_bi_fwd) load_mphi();
                    load_l_prev_next();
                } else if constexpr (has_rlzsa) {
                    store_sas_idx();
                    if (is_bidirectional) store_sas__idx();
                    build_l_prev_next();
                    store_l_prev_next();
                    store_mlf();

                    sort_iphim1();
                    construct_rlzsa<true, int32_t>();

                    load_mlf();
                    load_l_prev_next();
                    load_sas_idx();
                    if (is_bidirectional) load_sas__idx();
                } else if constexpr (support == _locate_bi_bwd) {
                    build_l_prev_next();
                }
            } else {
                build_l_prev_next();
            }

            if (has_rlzsa || p > 1 || is_bidirectional) {
                SA_file_bufs.clear();
                SA_file_bufs.shrink_to_fit();

                if constexpr (support == _locate_move_bi_fwd || support == _locate_rlzsa_bi_fwd) {
                    *sa_file_name = prefix_tmp_files + ".sa";
                } else {
                    std::filesystem::remove(prefix_tmp_files + ".sa");
                }
            }
        } else {
            build_l__sas();
            build_l_prev_next();
        }
    };

    // ############################# COMMON CONSTRUCTION METHODS #############################

    /**
     * @brief reads the input T and possibly remaps it to an internal alphabet, if it contains an invalid character
     * @param in_memory controls, whether the input should be processed in memory or read buffered from a file
     * @param bigbwt true <=> remap to a dense alphabet and make sure '\n' does not occur in the input
     * @param T_file_name name of the file containing T (for in_memory = false)
     */
    void preprocess_t(bool in_memory, bool bigbwt, const std::string& T_file_name);

    /**
     * @brief builds the RLBWT and C
     * @tparam mode determines how the BWT is read
     * @tparam sa_sint_t suffix array signed integer type
     */
    template <rlbwt_build_mode mode, typename sa_sint_t>
    void build_rlbwt_c();

    /**
     * @brief processes the C-array
     */
    void process_c();

    /**
     * @brief builds I_LF
     */
    void build_ilf();

    /**
     * @brief builds M_LF
     */
    void build_mlf();

    /**
     * @brief builds SA_s and SA_e from SA
     * @tparam bigbwt true <=> read I_Phi from the SA file output by Big-BWT
     * @tparam sa_sint_t suffix array signed integer type
     */
    template <bool bigbwt, typename sa_sint_t>
    void build_sas_from_sa();

    /**
     * @brief builds I_Phi^{m1} and SA_s and SA_e from SA
     * @tparam bigbwt true <=> read I_Phi from the SA file output by Big-BWT
     * @tparam sa_sint_t suffix array signed integer type
     */
    template <bool bigbwt, typename sa_sint_t>
    void build_iphim1_sas_from_sa();

    /**
     * @brief builds L' in M_LF (and SA_s)
     */
    void build_l__sas();

    /**
     * @brief builds I_Phi
     */
    void build_iphi();

    /**
     * @brief sorts I_Phi^{-1}
     */
    void sort_iphim1();

    /**
     * @brief sorts I_Phi
     */
    void sort_iphi();

    /**
     * @brief builds M_Phi^{-1}
     */
    void build_mphim1();

    /**
     * @brief builds M_Phi
     */
    void build_mphi();

    /**
     * @brief builds SA_Phi^{-1}
     */
    void build_saphim1();

    /**
     * @brief builds D_e
     */
    void build_de();

    /**
     * @brief builds RS_L'
     */
    void build_rsl_();

    /**
     * @brief builds L_prev and L_next
     */
    void build_l_prev_next();

    /**
     * @brief stores the RLBWT to disk
     */
    void store_rlbwt();

    /**
     * @brief loads the RLBWT from disk
     */
    void load_rlbwt();

    /**
     * @brief stores M_LF to disk
     */
    void store_mlf();

    /**
     * @brief loads M_LF from disk
     */
    void load_mlf();

    /**
     * @brief stores I_Phi^{-1} to disk
     */
    void store_iphim1();

    /**
     * @brief loads I_Phi^{-1} from disk
     */
    void load_iphim1();

    /**
     * @brief stores M_LF to disk
     */
    void store_mphi();

    /**
     * @brief loads M_LF from disk
     */
    void load_mphi();

    /**
     * @brief stores SA_s to disk
     */
    void store_sas();

    /**
     * @brief loads SA_s from disk
     */
    void load_sas();

    /**
     * @brief stores SA_s to disk
     */
    void store_sas_idx();

    /**
     * @brief loads SA_s from disk
     */
    void load_sas_idx();

    /**
     * @brief stores SA_e to disk
     */
    void store_sas__idx();

    /**
     * @brief loads SA_e from disk
     */
    void load_sas__idx();

    /**
     * @brief stores RS_L' to disk
     */
    void store_rsl_();

    /**
     * @brief stores L_prev and L_next to disk
     */
    void store_l_prev_next();

    /**
     * @brief loads RS_L' from disk
     */
    void load_rsl_();

    /**
     * @brief loads RS_L' from disk
     */
    void load_l_prev_next();
    
    /**
     * @brief reads I_Phi^{-1}
     */
    void read_iphim1_bigbwt();

    // ############################# IN-MEMORY CONSTRUCTION METHODS #############################

    /**
     * @brief reads T from t_file
     * @param T_file_name name of the file containing T
     */
    void read_t_from_file(std::string& T_file_name);

    /**
     * @brief builds the suffix array
     * @tparam sa_sint_t suffix array signed integer type
     */
    template <typename sa_sint_t>
    void build_sa();

    /**
     * @brief unmaps T from the internal alphabet
     * @param in_memory unmaps T in memory
     */
    void unmap_t(bool in_memory);

    /**
     * @brief stores map_int and map_ext to disk
     */
    void store_mapintext();

    /**
     * @brief loads map_int and map_ext from disk
     */
    void load_mapintext();

    // ############################# Big-BWT CONSTRUCTION METHODS #############################

    /**
     * @brief preprocesses T and stores it in a file
     */
    void store_t_in_file();

    /**
     * @brief computes the BWT (and suffix array samples) using Big-BWT
     * @param delete_T delete T file
     */
    void bigbwt(bool delete_T);
};

#include "common.tpp"
#include "bigbwt.tpp"
#include "rlbwt.tpp"
#include "move_lf.tpp"
#include "move_phi.tpp"
#include <rlzsa_opt/construction.hpp>
#include <rlzsa_opt/construction.tpp>
#include "suffix_array.tpp"
#include "load_store.tpp"