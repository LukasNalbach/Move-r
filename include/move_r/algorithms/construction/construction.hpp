#pragma once

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
    bool build_from_sa_and_l = false; // controls whether the index should be built from the suffix array and the bwt
    bool delete_T = false; // controls whether T should be deleted when not needed anymore
    bool log = false; // controls, whether to print log messages
    std::ostream* mf_idx = NULL; // file to write measurement data of the index construction to 
    std::ostream* mf_mds = NULL; // file to write measurement data of the move data structure construction to 
    std::string name_text_file = ""; // name of the text file (only for measurement output)
    std::string prefix_tmp_files = ""; // prefix of temporary files
    std::chrono::steady_clock::time_point time; // time of the start of the last build phase
    std::chrono::steady_clock::time_point time_start; // time of the start of the whole build phase
    uint64_t baseline_mem_usage = 0; // memory allocation at the start of the construction
    uint64_t bigbwt_peak_mem_usage = 0;
    bool build_count_support = false; // = true <=> build support for count (RS_L')
    // = true <=> build support for locate (either SA_Phi and M_Phi, or SA_s,R,PT,SP,SR and LP); and buildparallel revert support (D_e)
    bool build_locate_support = false;
    uint8_t min_valid_char = 0; // the minimum valid character that is allowed to occur in T
    uint8_t max_remapped_uchar = 0; // the maximum character in T that has been remappd
    uint8_t max_remapped_to_uchar = 0; // the maximum character in T that a character has been remappd to

    // ############################# INDEX VARIABLES #############################

    pos_t n = 0; // the length of T
    pos_t r = 0; // r, the number of runs in L
    pos_t r_ = 0; // r', the number of input/output intervals in M_LF
    pos_t r__ = 0; // r'', the number of input/output intervals in M_Phi

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
    /** The disjoint interval sequence for Phi */
    std::vector<std::pair<pos_t,pos_t>> I_Phi;
    /** [0..r'-1] if the end position of the i-th input interval of M_LF is the end position of a BWT run, then
     * SA_s'[i] is the suffix array sample at the end position of the i-th input interval of M_LF; else SA_s'[i] = n */
    std::vector<pos_t> SA_s_;
    /** [0..r'-1] Permutation storing the order of the values in SA_s' */
    std::vector<pos_t> pi_;
    /** [0..r''-1] Permutation storing the order of the output interval starting positions of M_Phi */
    std::vector<pos_t> pi_mphi;

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
     * @brief returns T at index i
     * @return T at index i
     */
    inline sym_t& T(pos_t i) {
        if constexpr (std::is_same<sym_t,char>::value) {
            return T_str[i];
        } else {
            return T_vec[i];
        }
    }

    /**
     * @brief returns the suffix array
     * @tparam sa_sint_t suffix array signed integer type
     * @return the suffix array
     */
    template <typename sa_sint_t>
    constexpr std::vector<sa_sint_t>& get_sa() {
        if constexpr (std::is_same<sa_sint_t,int32_t>::value) {
            return SA_32;
        } else {
            return SA_64;
        }
    }

    inline pos_t symbol_idx(sym_t sym) {
        if constexpr (std::is_same<sym_t,char>::value) {
            return char_to_uchar(sym);
        } else {
            return sym;
        }
    };

    /**
     * @brief sets the run length of the i-th BWT run in thread i_p's section to len
     * @param i_p [0..p-1] thread index
     * @param i [0..r-1] run index
     * @param len run length
     */
    inline void set_run_length(uint16_t i_p, pos_t i, uint32_t len) {
        RLBWT[i_p].set_unsafe<1,uint32_t>(i,len);
    }

    /**
     * @brief returns the length of the i-th BWT run in thread i_p's section
     * @param i_p [0..p-1] thread index
     * @param i [0..r-1] run index
     * @return run length
     */
    inline uint32_t run_length(uint16_t i_p, pos_t i) {
        return RLBWT[i_p].get_unsafe<1,uint32_t>(i);
    }

    /**
     * @brief sets the character of the i-th BWT run in thread i_p's section to c
     * @param i_p [0..p-1] thread index
     * @param i [0..r-1] run index
     * @param c character
     */
    inline void set_run_symbol(uint16_t i_p, pos_t i, sym_t c) {
        RLBWT[i_p].set<0,sym_t>(i,c);
    }

    /**
     * @brief returns the character of the i-th BWT run in thread i_p's section
     * @param i_p [0..p-1] thread index
     * @param i [0..r-1] run index
     * @return character
     */
    inline sym_t run_symbol(uint16_t i_p, pos_t i) {
        return RLBWT[i_p].get<0,sym_t>(i);
    }

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
            if (log) std::cout << "warning: p > n/1000, setting p to n/1000 ~ " << std::to_string(p) << std::endl;
        }
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
    requires(std::is_same<sym_t,char>::value)
    : T_str(T), T_vec(T_vec_tmp), L(L_tmp), SA_32(SA_32_tmp), SA_64(SA_64_tmp), idx(index) {
        this->delete_T = delete_T;
        read_parameters(params);
        prepare_phase_1();

        if (params.mode == _libsais) {
            min_valid_char = 1;
            T.push_back(uchar_to_char((uint8_t)0));
            n = T.size();
            idx.n = n;
            construct_libsais();

            if (!delete_T) {
                T.resize(n-1);
                if (idx.symbols_remapped) unmap_t();
            }
        } else {
            min_valid_char = 3;
            n = T.size()+1;
            idx.n = n;
            preprocess_and_store_t_in_file();
            construct_bigbwt();
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
    requires(!std::is_same<sym_t,char>::value)
    : T_str(T_str_tmp), T_vec(T), L(L_tmp), SA_32(SA_32_tmp), SA_64(SA_64_tmp), idx(index) {
        this->delete_T = delete_T;
        read_parameters(params);
        prepare_phase_1();
        T.push_back(0);
        n = T.size();
        idx.n = n;
        construct_libsais();

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
    requires(std::is_same<sym_t,char>::value)
    : T_str(T_str_tmp), T_vec(T_vec_tmp), idx(index), SA_32(SA_32_tmp), SA_64(SA_64_tmp), L(L_tmp) {
        read_parameters(params);
        prepare_phase_1();

        if (params.mode == _libsais) {
            min_valid_char = 1;
            read_t_from_file(T_ifile);
            construct_libsais();
        } else {
            min_valid_char = 3;
            preprocess_t(false,&T_ifile);
            construct_bigbwt();
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
    requires(std::is_same<sym_t,char>::value)
    : T_str(T_str_tmp), T_vec(T_vec_tmp), L(bwt), SA_32(suffix_array), SA_64(SA_64_tmp), idx(index) {
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
    requires(std::is_same<sym_t,char>::value)
    : T_str(T_str_tmp), T_vec(T_vec_tmp), L(bwt), SA_32(SA_32_tmp), SA_64(suffix_array), idx(index) {
        read_parameters(params);
        construct_from_sa_and_l<int64_t>();
    }

    // ############################# CONSTRUCTION #############################

    /**
     * @brief constructs the index from a suffix array and a bwt
     */
    template <typename sa_sint_t>
    void construct_from_sa_and_l() {
        build_from_sa_and_l = true;
        min_valid_char = 1;
        n = L.size();
        idx.n = n;

        prepare_phase_1();
        prepare_phase_2();
        build_rlbwt_c_libsais<true,sa_sint_t>();
        if (log) log_statistics();
        build_ilf();
        build_mlf();
        if (build_locate_support) build_iphi_from_sa<sa_sint_t>();
        build_l__sas_();

        if (build_locate_support) {
            if constexpr (locate_support == _phi) {
                sort_iphi();
                build_mphi();
                build_saphi();
                build_de();
            } else {
                build_rlzdsa<false,sa_sint_t>();
            }
        }

        if (build_count_support) build_rsl_();
        if (log) log_finished();
    }

    /**
     * @brief constructs the index from a string in memory (uses libsais)
     */
    void construct_libsais() {
        if constexpr (std::is_same<pos_t,uint32_t>::value) {
            if (n <= INT_MAX) {
                construct_libsais<int32_t>();
            } else {
                construct_libsais<int64_t>();
            }
        } else {
            construct_libsais<int64_t>();
        }
    }

    /**
     * @brief constructs the index from a string using libsais
     * @param sa_sint_t signed integer type to use for the suffix array entries
     */
    template <typename sa_sint_t = int32_t>
    void construct_libsais() {
        preprocess_t(true);
        prepare_phase_2();

        if constexpr (std::is_same<sym_t,char>::value) {
            build_sa<sa_sint_t,char>();
        } else {
            build_sa<sa_sint_t>();
        }
        
        build_rlbwt_c_libsais<false,sa_sint_t>();
        if (log) log_statistics();
        build_ilf();
        if (build_locate_support) build_iphi_from_sa<sa_sint_t>();
        build_mlf();
        build_l__sas_();

        if (build_locate_support) {
            if constexpr (locate_support == _phi) {
                sort_iphi();
                build_mphi();
                build_saphi();
                build_de();
            } else {
                build_rlzdsa<false,sa_sint_t>();
            }
        }

        if (build_count_support) build_rsl_();
    }

    /**
     * @brief constructs the index from a file using Big-BWT
     */
    void construct_bigbwt() {
        if constexpr (locate_support == _rlzdsa) {
            std::cout << "error: bigbwt is not supported when using rlzdsa";
            return;
        }

        prepare_phase_2();
        bigbwt();
        build_rlbwt_c_bigbwt();

        if (log) {
            std::cout
                << "peak memory allocation until now: "
                << format_size(std::max(malloc_count_peak()-baseline_mem_usage,bigbwt_peak_mem_usage))
                << std::endl;
            log_statistics();
        }

        build_ilf();
        build_mlf();
        if (build_locate_support) read_iphi_from_bigbwt();
        build_l__sas_();

        if (build_locate_support) {
            store_mlf();

            if constexpr (locate_support == _phi) {
                
                sort_iphi();
                build_mphi();
                build_saphi();
                build_de();
            } else {
                build_rlzdsa<false,int32_t>();
            }

            load_mlf();
        }

        if (build_count_support) build_rsl_();
    };

    // ############################# COMMON CONSTRUCTION METHODS #############################

    /**
     * @brief reads the input T and possibly remaps it to an internal alphabet, if it contains an invalid character
     * @param use_libsais controls, whether the input should be processed in memory or read buffered from a file
     * @param t_file file containing T (for use_libsais = false)
     */
    void preprocess_t(bool use_libsais, std::ifstream* T_ifile = NULL);

    /**
     * @brief builds the RLBWT and C
     * @tparam sa_sint_t suffix array signed integer type
     * @tparam read_l whether the RLBWT should be read from L
     */
    template <bool read_l, typename sa_sint_t = int32_t>
    void build_rlbwt_c_libsais();

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
     * @brief builds L' (and SA_s')
     */
    void build_l__sas_();

    /**
     * @brief sorts I_Phi
     */
    void sort_iphi();

    /**
     * @brief builds M_Phi
     */
    void build_mphi();

    /**
     * @brief builds SA_Phi
     */
    void build_saphi();

    /**
     * @brief builds D_e
     */
    void build_de();

    /**
     * @brief builds RS_L'
     */
    void build_rsl_();

    /**
     * @brief builds the rlzdsa
     * @tparam bigbwt controls, whether to use the suffix array file output by Big-BWT
     * @tparam sa_sint_t suffix array signed integer type
     */
    template <bool bigbwt, typename sa_sint_t>
    void build_rlzdsa();

    // ############################# IN-MEMORY CONSTRUCTION METHODS #############################

    /**
     * @brief reads T from t_file
     * @param t_file file containing T
     */
    void read_t_from_file(std::ifstream& t_file);

    /**
     * @brief execute the correct libsais algorithm
     * @tparam inp_t input type
     * @param T text
     * @param SA suffix array
     * @param fs free space at the end of the suffix array
     */
    template <typename inp_t>
    void execute_libsais(inp_t* T, int32_t* SA, pos_t fs);

    /**
     * @brief builds the suffix array
     * @tparam sa_sint_t suffix array signed integer type
     * @tparam inp_t input type
     */
    template <typename sa_sint_t, typename inp_t>
    void build_sa();

    /**
     * @brief builds the suffix array
     * @tparam sa_sint_t suffix array signed integer type
     */
    template <typename sa_sint_t>
    void build_sa() {
        uint8_t bytes = std::ceil(std::log2(idx.sigma+1)/(double)8);

        switch (bytes) {
            case 1:
                build_sa<sa_sint_t,uint8_t>();
                break;
            case 2:
                build_sa<sa_sint_t,uint16_t>();
                break;
            default:
                build_sa<sa_sint_t,int32_t>();
        }
    }

    /**
     * @brief builds SA_s' from SA in memory
     * @tparam sa_sint_t suffix array signed integer type
     */
    template <typename sa_sint_t>
    void build_iphi_from_sa();

    /**
     * @brief unmaps T from the internal alphabet
     */
    void unmap_t();

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
     * @brief builds the RLBWT and C
     */
    void build_rlbwt_c_bigbwt();

    /**
     * @brief reads I_Phi
     */
    void read_iphi_from_bigbwt();

    /**
     * @brief stores M_LF in a file
     */
    void store_mlf();

    /**
     * @brief loads M_LF from a file
     */
    void load_mlf();
};

#include "modes/common.cpp"
#include "modes/libsais.cpp"
#include "modes/bigbwt.cpp"
#include "modes/rlzdsa.cpp"