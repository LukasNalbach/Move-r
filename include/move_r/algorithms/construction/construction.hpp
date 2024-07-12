#pragma once

#include <move_r/move_r.hpp>
#include <gtl/phmap.hpp>
#include <gtl/btree.hpp>

enum rlbwt_build_mode {
    _sa, // the BWT is read by L[i] = T[(SA[i]-1) mod n]
    _bwt, // the BWT is read by accessing L[i]
    _bwt_file // the BWT is read from files output by Big-BWT
};

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
class move_r<locate_support,sym_t,pos_t>::construction {
    public:
    construction() = delete;
    construction(construction&&) = delete;
    construction(const construction&) = delete;
    construction& operator=(construction&&) = delete;
    construction& operator=(const construction&) = delete;
    ~construction() {}

    // ############################# MISC VARIABLES #############################

    std::string T_str_tmp;
    std::vector<sym_t> T_vec_tmp;
    std::string L_tmp;
    std::vector<int32_t> SA_32_tmp;
    std::vector<int64_t> SA_64_tmp;
    uint16_t p = 1; // the number of threads to use
    move_r_constr_mode mode = _suffix_array;
    /* the number of threads to use during the construction of the L,C and I_LF,I_Phi^{-1},L' and SA_s; when building move-r for an
     * integer alphabet, this needs (p+1)*sigma = O(p*n) words of space, so we limit the number of threads to p' to O(n/sigma)
     * in this case, so we choose p_ maximal, s.t. (p'+1)*sigma < n; this then ensures that we use O(n) words of space; */
    uint16_t p_ = 1;
    bool build_sa_and_l = false; // controls whether the index should be built from the suffix array and the bwt
    bool delete_T = false; // controls whether T should be deleted when not needed anymore
    bool log = false; // controls, whether to print log messages
    std::ostream* mf_idx = NULL; // file to write measurement data of the index construction to 
    std::ostream* mf_mds = NULL; // file to write measurement data of the move data structure construction to 
    std::string name_text_file = ""; // name of the text file (only for measurement output)
    std::string prefix_tmp_files = ""; // prefix of temporary files
    std::chrono::steady_clock::time_point time; // time of the start of the last build phase
    std::chrono::steady_clock::time_point time_start; // time of the start of the whole build phase
    uint64_t baseline_mem_usage = 0; // memory allocation at the start of the construction
    uint64_t bigbwt_peak_mem_usage = 0; // peak memory usage during the execution of Big-BWT
    bool build_count_support = false; // = true <=> build support for count (RS_L')
    // = true <=> build support for locate (either SA_Phi^{-1} and M_Phi^{-1}, or SA_s,R,PT,SP,SR and LP); and buildparallel revert support (D_e)
    bool build_locate_support = false;
    uint8_t min_valid_char = 0; // the minimum valid character that is allowed to occur in T
    uint8_t max_remapped_uchar = 0; // the maximum character in T that has been remapped
    uint8_t max_remapped_to_uchar = 0; // the maximum character in the effective alphabet of T that a character has been remapped to
    pos_t sigma_SAd = 0; // number of distinct values in SA^d (alphabet size of SA^d)
    pos_t size_R = 0; // size of the reference (R) in the rlzdsa
    pos_t size_R_target = 0; // target size for R
    pos_t seg_size = 0; // (maximum) size of each segment from SA^d to be included in R
    pos_t num_cand_segs = 0; // number of candidate segments that are considered in each iteration during the construction of R

    // ############################# INDEX VARIABLES #############################

    pos_t n = 0; // the length of T
    uint64_t n_u64 = 0; // the length of T
    pos_t r = 0; // r, the number of runs in L
    pos_t r_ = 0; // r', the number of input/output intervals in M_LF
    pos_t r__ = 0; // r'', the number of input/output intervals in M_Phi^{-1}

    // ############################# CONSTRUCTION DATA STRUCTURE VARIABLES #############################

    /** the string containing T */
    std::string& T_str;
    /** the vector containing T */
    std::vector<sym_t>& T_vec;
    /** The move-r index to construct */
    move_r<locate_support,sym_t,pos_t>& idx;
    /** [0..n-1] The suffix array (32-bit) */
    std::vector<int32_t>& SA_32;
    /** [0..n-1] The suffix array (64-bit) */
    std::vector<int64_t>& SA_64;
    /** [0..n-1] The BWT */
    std::string& L;
    /** [0..p-1] file buffers (of each thread) for reading the BWT file output by bigbwt */
    std::vector<sdsl::int_vector_buffer<>> BWT_file_bufs;
    /** [0..p-1] vectors that contain the RLBWT concatenated */
    std::vector<interleaved_vectors<uint32_t,uint32_t>> RLBWT;
    /** [0..p] n_p[0] < n_p[1] < ... < n_p[p] = n; n_p[i] = start position of thread i's section in L and SA */
    std::vector<pos_t> n_p;
    /** [0..p] r_p[0] < r_p[1] < ... < r_p[p] = r; r_p[i] = index of the first run in L starting in
     * [n_p[i]..n_p[i+1]-1]; there is a run starting at n_p[i] */
    std::vector<pos_t> r_p;
    /** [0..p][0..255] see the code to see how this variable is used */
    std::vector<std::vector<pos_t>> C;
    /** The disjoint interval sequence for LF */
    std::vector<std::pair<pos_t,pos_t>> I_LF;
    /** The disjoint interval sequence for Phi^{-1} */
    std::vector<std::pair<pos_t,pos_t>> I_Phi_m1;
    /** [0..r'-1] SA_s[x] = SA[M_LF.p[x]]; if locate_support = _mds, and the starting position of the
     * x-th input interval of M_LF is not starting position of a BWT run, then SA_s[x] = n */
    std::vector<pos_t> SA_s;
    /** [0..r'-1] Permutation storing the order of the values in SA_s */
    std::vector<pos_t> pi_;
    /** [0..r''-1] Permutation storing the order of the output interval starting positions of M_Phi^{-1} */
    std::vector<pos_t> pi_mphi;
    /** [0..p-1] file buffers (of each thread) for reading the suffix array file output by bigbwt */
    std::vector<sdsl::int_vector_buffer<>> SA_file_bufs;
    /** type of hash map for storing the frequencies of values in SA^d */
    template <typename sad_t> using sad_freq_t = gtl::flat_hash_map<sad_t,pos_t>;
    /** [0..p-1] hashmaps (for sad_t = uint32_t) mapping each value in SA^d (+n) to its frequency in SA^d; initially
     * stores at index i_p the hashmap that thread i_p uses to count the frequencies of all values in its section
     * SA^d[n_p[i_p]..n_p[i_p+1]-1]; then the hashmaps are merged into the hashmap at index 0 */
    std::vector<sad_freq_t<uint32_t>> SAd_freq_32;
    /** [0..p-1] hashmaps ... (for sad_t = uint64_t; see SAd_freq_32) */
    std::vector<sad_freq_t<uint64_t>> SAd_freq_64;
    using ts_t = gtl::btree_set<std::pair<pos_t,pos_t>>; // type of T_s
    using ts_it_t = ts_t::iterator; // type of iterator in T_s
    /** B-tree storing the selected segments from SA^d to be included in the reference (R) for the rlzdsa; stores pairs (pos,len)
     * to represent segments SA^d[pos..pos+len] in ascending order of their starting position (pos) in SA^d */
    ts_t T_s;
    /** [0..size_R-1] rev(R) (for sad_t = uint32_t) */
    std::vector<uint32_t> revR_32;
    /** [0..size_R-1] rev(R) (for sad_t = uint64_t) */
    std::vector<uint64_t> revR_64;
    /** type of move-r index for finding copy phrases in the rlzdsa factorization */
    template <typename sad_t, typename irr_pos_t> using idx_revr_t = move_r<_mds,sad_t,irr_pos_t>;
    /** index for finding maximum length copy phrases in the rlzdsa factorization (for sad_t = uint32_t and irr_pos_t = uint32_t) */
    idx_revr_t<uint32_t,uint32_t> idx_revR_32_32;
    /** index for finding maximum length copy phrases in the rlzdsa factorization (for sad_t = uint64_t and irr_pos_t = uint32_t) */
    idx_revr_t<uint64_t,uint32_t> idx_revR_64_32;
    /** index for finding maximum length copy phrases in the rlzdsa factorization (for sad_t = uint64_t and irr_pos_t = uint64_t) */
    idx_revr_t<uint64_t,uint64_t> idx_revR_64_64;
    /** [0..p-1] stores at position i_p in ascending order the starting positions (in SA^d) of rlzdsa
     * copy phrases in SA^d[n_p[i_p]..n_p[i_p+1]-1] */
    std::vector<interleaved_vectors<pos_t,pos_t>> SCP;
    /** [0..p-1] stores at position i_p in ascending order the starting positions (in R) of rlzdsa
     * copy phrases in SA^d[n_p[i_p]..n_p[i_p+1]-1] */
    std::vector<interleaved_vectors<pos_t,pos_t>> SR;
    /** [0..p-1] stores at position i_p in ascending order the literal phrases of the rlzdsa
     * in SA^d[n_p[i_p]..n_p[i_p+1]-1] */
    std::vector<interleaved_vectors<pos_t,pos_t>> LP;
    /** [0..p-1] stores at position i_p in ascending order the types of phrases of the rlzdsa in SA^d[n_p[i_p]..n_p[i_p+1]-1];
     * i.e, PT[i_p][i] = 0 <=> the i-th phrase in SA^d[n_p[i_p]..n_p[i_p+1]-1] is literal */
    std::vector<sdsl::bit_vector> PT;
    /** [0..p-1] stores at position i_p a file buffer with the same content as SCP[i_p]*/
    std::vector<sdsl::int_vector_buffer<>> SCP_file_bufs;
    /** [0..p-1] stores at position i_p a file buffer with the same content as SR[i_p]*/
    std::vector<sdsl::int_vector_buffer<>> SR_file_bufs;
    /** [0..p-1] stores at position i_p a file buffer with the same content as LP[i_p]*/
    std::vector<sdsl::int_vector_buffer<>> LP_file_bufs;
    /** [0..p-1] stores at position i_p a file buffer with the same content as PT[i_p]*/
    std::vector<sdsl::int_vector_buffer<>> PT_file_bufs;
    /** sd_vector builder to build idx._SCP from SCP[0..p-1] */
    sdsl::sd_vector_builder SCP_b;
    /** [0..p] z_p[0] < z_p[1] < ... < z_p[p] = z; z_p[i_p] = number of phrases in the rlzdsa starting before n_p[i_p] in SA^d */
    std::vector<pos_t> z_p;
    /** [0..p] zl_p[0] < zl_p[1] < ... < zl_p[p] = z; zl_p[i_p] = number of literal phrases in the rlzdsa starting before n_p[i_p] in SA^d */
    std::vector<pos_t> zl_p;
    /** [0..p] zc_p[0] < zc_p[1] < ... < zc_p[p] = z; zc_p[i_p] = number of copy phrases in the rlzdsa starting before n_p[i_p] in SA^d */
    std::vector<pos_t> zc_p;

    // ############################# COMMON MISC METHODS #############################

    /**
     * @brief reports the operations supported by this index
     */
    void report_supports() {
        std::cout << "building an index with the following support:" << std::endl;

        static auto report_support = [](move_r_supp sup){
            if (sup == _revert) {
                std::cout << "revert";
            } else if (sup == _count) {
                std::cout << "count";
            } else if (sup == _locate) {
                std::cout << "locate";
            }
        };

        for (uint8_t i=0; i<idx._support.size()-1; i++) {
            report_support(idx._support[i]);
            std::cout << ", ";
        }

        report_support(idx._support.back());
        std::cout << std::endl;
    }

    /**
     * @brief returns T at index i interpreted as type
     * @return T at index i
     */
    template <typename type>
    inline type& T(pos_t i) {
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
    constexpr std::vector<sa_sint_t>& get_sa() {
        if constexpr (std::is_same_v<sa_sint_t,int32_t>) {
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
    inline void set_run_len(uint16_t i_p, pos_t i, pos_t len) {
        RLBWT[i_p].set_unsafe<1,uint32_t>(i,len);
    }

    /**
     * @brief returns the length of the i-th BWT run in thread i_p's section
     * @param i_p [0..p-1] thread index
     * @param i [0..r-1] run index
     * @return run length
     */
    inline pos_t run_len(uint16_t i_p, pos_t i) {
        return RLBWT[i_p].get_unsafe<1,uint32_t>(i);
    }

    /**
     * @brief sets the character of the i-th BWT run in thread i_p's section to c
     * @param i_p [0..p-1] thread index
     * @param i [0..r-1] run index
     * @param c character
     */
    inline void set_run_sym(uint16_t i_p, pos_t i, i_sym_t c) {
        if constexpr (byte_alphabet) {
            RLBWT[i_p].set_unsafe<0,i_sym_t>(i,c);
        } else {
            RLBWT[i_p].set<0,i_sym_t>(i,c);
        }
    }

    /**
     * @brief returns the character of the i-th BWT run in thread i_p's section
     * @param i_p [0..p-1] thread index
     * @param i [0..r-1] run index
     * @return character
     */
    inline i_sym_t run_sym(uint16_t i_p, pos_t i) {
        if constexpr (byte_alphabet) {
            return RLBWT[i_p].get_unsafe<0,i_sym_t>(i);
        } else {
            return RLBWT[i_p].get<0,i_sym_t>(i);
        }
    }

    /**
     * @brief adds a bwt run of length len with the symbol sym to the end of thread i_p's section of the RLBWT
     * @param i_p [0..p-1] thread index
     * @param sym run symbol
     * @param len run length
     */
    inline void add_run(uint16_t i_p, i_sym_t sym, pos_t len) {
        if constexpr (byte_alphabet) {
            RLBWT[i_p].emplace_back_unsafe<i_sym_t,uint32_t>({sym,len});
        } else {
            RLBWT[i_p].emplace_back<i_sym_t>({sym});
            set_run_len(i_p,RLBWT[i_p].size()-1,len);
        }
    }
    
    // ############################# rlzdsa MISC METHODS #############################

    /**
     * @brief returns SA[i]
     * @tparam bigbwt true <=> read SA[i] from SA_file
     * @tparam sa_sint_t suffix array signed integer type
     * @return SA[i]
     */
    template <bool bigbwt, typename sa_sint_t>
    inline pos_t SA(uint16_t i_p, pos_t i){
        if constexpr (bigbwt) {
            if (i == 0) return n-1;
            return pos_t{SA_file_bufs[i_p][i-1]};
        } else {
            return get_sa<sa_sint_t>()[i];
        }
    }

    /**
     * @brief returns SA^d[i]
     * @tparam bigbwt true <=> read SA[i] and SA[i-1] from SA_file
     * @tparam sa_sint_t suffix array signed integer type
     * @return SA^d[i]
     */
    template <bool bigbwt, typename sa_sint_t>
    inline uint64_t SAd(uint16_t i_p, pos_t i){
        if constexpr (bigbwt) {
            if (i == 0) return n_u64-1;
            if (i == 1) return SA_file_bufs[i_p][0]+1;
            return (SA_file_bufs[i_p][i-1]+n_u64)-SA_file_bufs[i_p][i-2];
        } else {
            return i == 0 ? n_u64-1 : ((get_sa<sa_sint_t>()[i]+n)-get_sa<sa_sint_t>()[i-1]);
        }
    }

    /**
     * @brief returns the suffix array
     * @tparam sad_t type of the values in SA^d
     * @return SAd_freq
     */
    template <typename sad_t>
    constexpr std::vector<sad_freq_t<sad_t>>& get_SAd_freq() {
        if constexpr (std::is_same_v<sad_t,uint32_t>) {
            return SAd_freq_32;
        } else {
            return SAd_freq_64;
        }
    }

    /**
     * @brief returns idx_revR
     * @tparam sad_t type of the values in SA^d
     * @return idx_revR
     */
    template <typename sad_t, typename irr_pos_t>
    constexpr idx_revr_t<sad_t,irr_pos_t>& get_idx_revR() {
        if constexpr (std::is_same_v<sad_t,uint32_t>) {
            return idx_revR_32_32;
        } else {
            if constexpr (std::is_same_v<irr_pos_t,uint32_t>) {
                return idx_revR_64_32;
            } else {
                return idx_revR_64_64;
            }
        }
    }

    /**
     * @brief returns revR
     * @tparam sad_t type of the values in SA^d
     * @return revR
     */
    template <typename sad_t>
    constexpr std::vector<sad_t>& get_revR() {
        if constexpr (std::is_same_v<sad_t,uint32_t>) {
            return revR_32;
        } else {
            return revR_64;
        }
    }

    // ############################# COMMON MISC METHODS #############################

    /**
     * @brief sets some variables and logs
     */
    void prepare_phase_1() {
        time = now();
        time_start = time;
        omp_set_num_threads(p);

        baseline_mem_usage = malloc_count_current();
        if (log) malloc_count_reset_peak();

        adjust_supports(idx._support);

        if (!contains(idx._support,_revert)) {
            std::cout << "error: cannot build an index without revert support";
            return;
        }

        build_count_support = contains(idx._support,_count);
        build_locate_support = contains(idx._support,_locate);

        // print the operations to build support for
        if (log) {
            report_supports();
            std::cout << std::endl;
        }
    }

    /**
     * @brief sets some variables
     */
    void prepare_phase_2() {
        if (p > 1 && 1000*p > n) {
            p = std::max<pos_t>(1,n/1000);
            p_ = std::min<uint16_t>(p_,p);
            if (log) std::cout << "warning: p > n/1000, setting p to n/1000 ~ " << std::to_string(p) << std::endl;
        }
    }

    /**
     * @brief logs the current memory usage
     */
    void log_mem_usage() {
        std::cout << "current memory allocation: "
            << format_size(malloc_count_current()-baseline_mem_usage)
            << std::endl;
    }

    /**
     * @brief logs the peak memory usage until now
     */
    void log_peak_mem_usage() {
        std::cout << "peak memory allocation until now: "
            << format_size(std::max(malloc_count_peak()-baseline_mem_usage,bigbwt_peak_mem_usage))
            << std::endl;
    }

    /**
     * @brief logs a construction summary
     */
    void log_finished() {
        uint64_t time_construction = time_diff_ns(time_start,now());
        uint64_t peak_mem_usage = std::max(malloc_count_peak()-baseline_mem_usage,bigbwt_peak_mem_usage);
        
        std::cout << std::endl;
        std::cout << "construction time: " << format_time(time_construction) << std::endl;
        std::cout << "peak memory usage: " << format_size(peak_mem_usage) << std::endl;
        idx.log_data_structure_sizes();

        if (mf_idx != NULL) {
            *mf_idx << " time_construction=" << time_construction;
            *mf_idx << " peak_mem_usage=" << peak_mem_usage;
            idx.log_data_structure_sizes(*mf_idx);
            *mf_idx << std::endl;
        }
    }

    /**
     * @brief logs statistics of T
     */
    void log_statistics() {
        double n_r = std::round(100.0*(n/(double)r))/100.0;
        if (mf_idx != NULL) {
            *mf_idx << " n=" << n;
            *mf_idx << " sigma=" << std::to_string(idx.sigma);
            *mf_idx << " r=" << r;
        }
        std::cout << "n = " << n << ", sigma = " << std::to_string(idx.sigma) << ", r = " << r << ", n/r = " << n_r << std::endl;
    }

    // ############################# CONSTRUCTORS #############################

    void read_parameters(move_r_params& params) {
        idx._support = params.support;
        this->p = params.num_threads;
        this->mode = params.mode;
        idx.a = params.a;
        this->log = params.log;
        this->mf_idx = params.mf_idx;
        this->mf_mds = params.mf_mds;
        this->name_text_file = params.name_text_file;
    }

    /**
     * @brief constructs a move_r index of the string input
     * @param index The move-r index to construct
     * @param T the string containing T 
     * @param delete_T controls whether T should be deleted once it is not needed anymore
     * @param params construction parameters
     */
    construction(move_r<locate_support,sym_t,pos_t>& index, std::string& T, bool delete_T, move_r_params params)
    requires(str_input) : T_str(T), T_vec(T_vec_tmp), L(L_tmp), SA_32(SA_32_tmp), SA_64(SA_64_tmp), idx(index) {
        this->delete_T = delete_T;
        read_parameters(params);
        prepare_phase_1();

        if (mode == _suffix_array || mode == _suffix_array_space) {
            min_valid_char = 1;
            T.push_back(uchar_to_char((uint8_t)0));
            n = T.size();
            idx.n = n;
            preprocess_t(true);
            construct_from_sa();

            if (!delete_T) {
                T.resize(n-1);
                if (idx.symbols_remapped) unmap_t();
            }
        } else {
            min_valid_char = 3;
            n = T.size()+1;
            idx.n = n;
            preprocess_and_store_t_in_file();
            construct_from_bigbwt();
        }

        if (log) log_finished();
    }

    /**
     * @brief constructs a move_r index of the vector input
     * @param index The move-r index to construct
     * @param T the vector containing T
     * @param alphabet_size alphabet size of the input vector (= maximum value in the input vector)
     * @param delete_T controls whether T should be deleted once it is not needed anymore
     * @param params construction parameters
     */
    construction(move_r<locate_support,sym_t,pos_t>& index, std::vector<sym_t>& T, bool delete_T, move_r_params params)
    requires(int_input) : T_str(T_str_tmp), T_vec(T), L(L_tmp), SA_32(SA_32_tmp), SA_64(SA_64_tmp), idx(index) {
        this->delete_T = delete_T;
        read_parameters(params);
        prepare_phase_1();
        T.push_back(0);
        n = T.size();
        idx.n = n;
        preprocess_t(true);
        construct_from_sa();

        if (!delete_T) {
            T.resize(n-1);
            if (idx.symbols_remapped) unmap_t();
        }

        if (log) log_finished();
    }

    /**
     * @brief constructs a move_r index from an input file
     * @param index The move-r index to construct
     * @param T_ifile file containing T
     * @param params construction parameters
     */
    construction(move_r<locate_support,sym_t,pos_t>& index, std::ifstream& T_ifile, move_r_params params)
    requires(str_input) : T_str(T_str_tmp), T_vec(T_vec_tmp), idx(index), SA_32(SA_32_tmp), SA_64(SA_64_tmp), L(L_tmp) {
        read_parameters(params);
        prepare_phase_1();

        if (mode == _suffix_array || mode == _suffix_array_space) {
            min_valid_char = 1;
            read_t_from_file(T_ifile);
            preprocess_t(true);
            construct_from_sa();
        } else {
            min_valid_char = 3;
            preprocess_t(false,&T_ifile);
            construct_from_bigbwt();
        }

        if (log) log_finished();
    }

    /**
     * @brief constructs a move_r index from an input file
     * @param index The move-r index to construct
     * @param suffix_array vector containing the suffix array of the input
     * @param bwt string containing the bwt of the input
     * @param params construction parameters
     */
    construction(move_r<locate_support,sym_t,pos_t>& index, std::vector<int32_t>& suffix_array, std::string& bwt, move_r_params params)
    requires(str_input) : T_str(T_str_tmp), T_vec(T_vec_tmp), L(bwt), SA_32(suffix_array), SA_64(SA_64_tmp), idx(index) {
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
    construction(move_r<locate_support,sym_t,pos_t>& index, std::vector<int64_t>& suffix_array, std::string& bwt, move_r_params params)
    requires(str_input) : T_str(T_str_tmp), T_vec(T_vec_tmp), L(bwt), SA_32(SA_32_tmp), SA_64(suffix_array), idx(index) {
        read_parameters(params);
        construct_from_sa_and_l<int64_t>();
    }

    // ############################# CONSTRUCTION #############################

    /**
     * @brief constructs the index from a suffix array and a bwt
     */
    template <typename sa_sint_t>
    void construct_from_sa_and_l() {
        build_sa_and_l = true;
        min_valid_char = 1;
        n = L.size();
        idx.n = n;

        prepare_phase_1();
        prepare_phase_2();
        build_rlbwt_c<_bwt,sa_sint_t>();
        if (log) log_statistics();
        build_ilf();
        build_mlf();

        if (build_locate_support) {
            build_iphim1_sa<sa_sint_t>();
            build_l__sas<locate_support == _mds>();
        } else {
            build_l__sas<false>();
        }

        if (build_locate_support) {
            if constexpr (locate_support == _mds) {
                sort_iphim1();
                build_mphim1();
                build_saphim1();
                build_de();
            } else {
                construct_rlzdsa<false,sa_sint_t>();
                build_sas_rlzdsa<false,sa_sint_t>();
            }
        }

        if (build_count_support) build_rsl_();
        if (log) log_finished();
    }

    /**
     * @brief constructs the index in memory (uses libsais)
     */
    void construct_from_sa() {
        if (n <= INT_MAX && (sizeof(sym_t) < 8 || mode == _suffix_array_space)) {
            construct_from_sa<int32_t>();
        } else {
            construct_from_sa<int64_t>();
        }
    }

    /**
     * @brief constructs the index using the suffix array
     * @param sa_sint_t signed integer type to use for the suffix array entries
     */
    template <typename sa_sint_t>
    void construct_from_sa() {
        bool _space = mode == _suffix_array_space;

        prepare_phase_2();
        build_sa<sa_sint_t>();
        build_rlbwt_c<_sa,sa_sint_t>();
        if (log) log_statistics();
        build_ilf();

        if (build_locate_support) {
            build_iphim1_sa<sa_sint_t>();
            if (_space) store_iphim1();
        }

        build_mlf();

        if (build_locate_support) {
            if constexpr (locate_support == _mds) {
                if (_space) load_iphim1();
                build_l__sas<true>();
                if (_space) store_sas();
                build_rsl_();
                if (_space) store_rsl_();
                if (_space) store_mlf();
                sort_iphim1();
                build_mphim1();
                if (_space) load_sas();
                build_saphim1();
                build_de();
                if (_space) load_mlf();
                if (_space) load_rsl_();
            } else {
                build_l__sas<false>();
                build_rsl_();
                if (_space) store_rsl_();
                if (_space) store_mlf();
                construct_rlzdsa<false,sa_sint_t>();
                if (_space) load_mlf();
                build_sas_rlzdsa<false,sa_sint_t>();
                if (_space) load_rsl_();
            }
        } else {
            build_l__sas<false>();
            if (build_count_support) build_rsl_();
        }

        if constexpr (!byte_alphabet) if (mode == _suffix_array_space) load_mapintext();
    }

    /**
     * @brief constructs the index from a file using Big-BWT
     */
    void construct_from_bigbwt() {
        prepare_phase_2();
        bigbwt();
        build_rlbwt_c<_bwt_file,int32_t>();
        if (log) log_statistics();
        build_ilf();
        build_mlf();

        if (build_locate_support) {
            if constexpr (locate_support == _mds) {
                read_iphim1_bigbwt();
                build_l__sas<true>();
                store_sas();
                build_rsl_();
                store_mlf();
                store_rsl_();
                sort_iphim1();
                build_mphim1();
                load_sas();
                build_saphim1();
                build_de();
                load_mlf();
                load_rsl_();
            } else {
                build_l__sas<false>();
                build_rsl_();
                store_mlf();
                store_rsl_();
                construct_rlzdsa<true,int32_t>();
                load_mlf();
                build_sas_rlzdsa<true,int32_t>();
                load_rsl_();
            }
        } else {
            build_l__sas<false>();
            if (build_count_support) build_rsl_();
        }
    };

    /**
     * @brief builds the rlzdsa
     * @tparam bigbwt true <=> read suffix array values from the file output by bigbwt
     * @tparam sa_sint_t signed integer type to use for the suffix array entries
     */
    template <bool bigbwt, typename sa_sint_t>
    void construct_rlzdsa() {
        size_R_target = std::min(std::max<pos_t>(1,n/3),5*r_);
        seg_size = std::min<pos_t>(3072,size_R_target);

        // TODO: somehow detect if choosing seg_size >> 3072 increases the compression rate (optimum is between 2000 and 12000)
        // generally, the optimum is approx. 2000-4000, but sometimes it can be approx. 12000 (Horspool)

        if constexpr (std::is_same_v<pos_t,uint32_t>) {
            construct_rlzdsa<bigbwt,uint32_t,uint32_t,sa_sint_t>();
        } else if (2*n <= UINT_MAX) {
            construct_rlzdsa<bigbwt,uint32_t,uint32_t,sa_sint_t>();
        } else if (size_R_target+seg_size <= UINT_MAX) {
            construct_rlzdsa<bigbwt,uint64_t,uint32_t,sa_sint_t>();
        } else {
            construct_rlzdsa<bigbwt,uint64_t,uint64_t,sa_sint_t>();
        }
    }

    /**
     * @brief builds the rlzdsa
     * @tparam bigbwt true <=> read suffix array values from the file output by bigbwt
     * @tparam sad_t type of the values in SA^d
     * @tparam irr_pos_t position type (pos_t) for the index of rev(R)
     * @tparam sa_sint_t signed integer type to use for the suffix array entries
     */
    template <bool bigbwt, typename sad_t, typename irr_pos_t, typename sa_sint_t>
    void construct_rlzdsa() {
        bool _space = bigbwt || mode == _suffix_array_space;
        build_freq_sad<bigbwt,sad_t,sa_sint_t>();
        build_r_revR<bigbwt,sad_t,sa_sint_t>();
        if (_space) store_r();
        build_idx_rev_r<sad_t,irr_pos_t>();

        if (_space) {
            build_rlzdsa_factorization<bigbwt,true,sad_t,irr_pos_t,sa_sint_t>();
            load_r();
        } else {
            build_rlzdsa_factorization<bigbwt,false,sad_t,irr_pos_t,sa_sint_t>();
        }
    }

    // ############################# COMMON CONSTRUCTION METHODS #############################

    /**
     * @brief reads the input T and possibly remaps it to an internal alphabet, if it contains an invalid character
     * @param in_memory controls, whether the input should be processed in memory or read buffered from a file
     * @param t_file file containing T (for in_memory = false)
     */
    void preprocess_t(bool in_memory, std::ifstream* T_ifile = NULL);

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
     * @brief builds L' in M_LF (and SA_s)
     * @tparam build_sas_ controls, whether to build SA_s
     */
    template <bool build_sas_>
    void build_l__sas();

    /**
     * @brief sorts I_Phi^{-1}
     */
    void sort_iphim1();

    /**
     * @brief builds M_Phi^{-1}
     */
    void build_mphim1();

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
     * @brief stores M_LF to disk
     */
    void store_mlf();

    /**
     * @brief loads M_LF from disk
     */
    void load_mlf();

    /**
     * @brief stores SA_s to disk
     */
    void store_sas();

    /**
     * @brief loads SA_s from disk
     */
    void load_sas();

    /**
     * @brief stores RS_L' to disk
     */
    void store_rsl_();

    /**
     * @brief loads RS_L' from disk
     */
    void load_rsl_();

    // ############################# IN-MEMORY CONSTRUCTION METHODS #############################

    /**
     * @brief reads T from t_file
     * @param t_file file containing T
     */
    void read_t_from_file(std::ifstream& t_file);

    /**
     * @brief execute the correct libsais algorithm
     * @tparam sa_sint_t suffix array signed integer type
     * @tparam ls_sym_t libsais symbol type
     * @param T text
     * @param SA suffix array
     * @param fs free space at the end of the suffix array
     */
    template <typename ls_sym_t, typename sa_sint_t>
    void run_libsais(ls_sym_t* T, sa_sint_t* SA, pos_t fs);

    /**
     * @brief builds the suffix array
     * @tparam sa_sint_t suffix array signed integer type
     */
    template <typename sa_sint_t>
    void build_sa();

    /**
     * @brief builds I_Phi^{m1} from SA in memory
     * @tparam sa_sint_t suffix array signed integer type
     */
    template <typename sa_sint_t>
    void build_iphim1_sa();

    /**
     * @brief unmaps T from the internal alphabet
     */
    void unmap_t();

    /**
     * @brief stores I_Phi^{-1} to disk
     */
    void store_iphim1();

    /**
     * @brief loads I_Phi^{-1} from disk
     */
    void load_iphim1();

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
    void preprocess_and_store_t_in_file();

    /**
     * @brief computes the BWT (and suffix array samples) using Big-BWT
     */
    void bigbwt();

    /**
     * @brief reads I_Phi^{-1}
     */
    void read_iphim1_bigbwt();

    // ############################# rlzdsa CONSTRUCTION METHODS #############################

    /**
     * @brief builds freq_SAd
     * @tparam bigbwt true <=> read suffix array values from the file output by bigbwt
     * @tparam sad_t type of the values in SA^d
     * @tparam sa_sint_t signed integer type to use for the suffix array entries
     */
    template <bool bigbwt, typename sad_t, typename sa_sint_t>
    void build_freq_sad();

    /**
     * @brief builds R and rev(R)
     * @tparam bigbwt true <=> read suffix array values from the file output by bigbwt
     * @tparam sad_t type of the values in SA^d
     * @tparam sa_sint_t signed integer type to use for the suffix array entries
     */
    template <bool bigbwt, typename sad_t, typename sa_sint_t>
    void build_r_revR();

    /**
     * @brief builds the move-r index of rev(R)
     * @tparam sad_t type of the values in SA^d
     * @tparam irr_pos_t position type (pos_t) for the index of rev(R)
     */
    template <typename sad_t, typename irr_pos_t>
    void build_idx_rev_r();

    /**
     * @brief builds the rlzdsa factorization
     * @tparam bigbwt true <=> read suffix array values from the file output by bigbwt
     * @tparam space true <=> store rlzdsa data structures in files during the construction
     * @tparam sad_t type of the values in SA^d
     * @tparam irr_pos_t position type (pos_t) for the index of rev(R)
     * @tparam sa_sint_t signed integer type to use for the suffix array entries
     */
    template <bool bigbwt, bool space, typename sad_t, typename irr_pos_t, typename sa_sint_t>
    void build_rlzdsa_factorization();

    /**
     * @brief stores R to disk
     */
    void store_r();

    /**
     * @brief loads R from disk
     */
    void load_r();

    /**
     * @brief builds SA_s
     * @tparam bigbwt true <=> read suffix array values from the file output by bigbwt
     * @tparam sa_sint_t signed integer type to use for the suffix array entries
     */
    template <bool bigbwt, typename sa_sint_t>
    void build_sas_rlzdsa();
};

#include "modes/common.cpp"
#include "modes/sa.cpp"
#include "modes/bigbwt.cpp"
#include "modes/rlzdsa.cpp"