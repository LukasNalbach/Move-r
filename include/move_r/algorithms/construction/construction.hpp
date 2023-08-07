#pragma once

template <typename uint_t>
class move_r<uint_t>::construction {
    static_assert(std::is_same<uint_t,uint32_t>::value || std::is_same<uint_t,uint64_t>::value);

    public:
    construction() = delete;
    construction(construction&&) = delete;
    construction(const construction&) = delete;
    construction& operator=(construction&&) = delete;
    construction& operator=(const construction&) = delete;
    ~construction() {}

    // ############################# MISC VARIABLES #############################

    std::string T_tmp;
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
    std::string prefix_tempfiles = ""; // prefix of temporary files
    std::chrono::steady_clock::time_point time; // time of the start of the last build phase
    std::chrono::steady_clock::time_point time_start; // time of the start of the whole build phase
    uint64_t baseline_memory_allocation = 0; // memory allocation at the start of the construction
    bool build_count_support = false; // = true <=> build support for count (RS_L')
    bool build_locate_support = false; // = true <=> build support for locate (SA_idx and M_Phi) and build parallel revert support (D_e)
    bool read_rlbwt = false; // = true <=> read the run-length encoded BWT
    uint8_t terminator = 0; // terminator (dollar) symbol
    uint8_t min_valid_char = 0; // the minimum valid character that is allowed to occur in T
    uint8_t max_remapped_uchar = 0; // the maximum character in T, that has been remappd
    uint8_t max_remapped_to_uchar = 0; // the maximum character in T, that a character has been remappd to

    // ############################# INDEX VARIABLES #############################

    std::vector<move_r_support> support;
    uint_t n = 0; // the length of T
    uint8_t sigma = 0; // the number of distinct characters (including the terminator symbol 1) of T
    uint_t r = 0; // r, the number of runs in L
    uint_t r_ = 0; // r', the number of input/output intervals in M_LF
    uint_t r__ = 0; // r'', the number of input/output intervals in M_Phi
    uint16_t a = 0; // balancing parameter, restricts size to O(r*(a/(a-1))), 2 <= a
    uint16_t p_r = 0; // maximum possible number of threads to use while reverting the index
    uint8_t omega_idx = 0; // word width of SA_idx

    // ############################# CONSTRUCTION DATA STRUCTURE VARIABLES #############################

    /** the string containing T */
    std::string& T;
    /** The move-r index to construct */
    move_r<uint_t>& idx;
    /** [0..n-1] The suffix array (32-bit) */
    std::vector<int32_t>& SA_32;
    /** [0..n-1] The suffix array (64-bit) */
    std::vector<int64_t>& SA_64;
    /** [0..n-1] The BWT */
    std::string& L;
    /** [0..p-1] vectors that contain the RLBWT concatenated */
    std::vector<std::vector<std::pair<char,uint32_t>>> RLBWT_thr;
    /** [0..r-1] the RLBWT */
    interleaved_vectors<uint32_t> RLBWT;
    /** [0..256] marks at position c whether the character c occurs in T */
    std::vector<uint8_t> contains_uchar = {};
    /** [0..p] n_p[0] < n_p[1] < ... < n_p[p] = n; n_p[i] = start position of thread i's section in L and SA */
    std::vector<uint_t> n_p;
    /** [0..p] r_p[0] < r_p[1] < ... < r_p[p] = r; r_p[i] = index of the first run in L, that starts in
     * [n_p[i]..n_p[i+1]-1]; there is a run starting at n_p[i] */
    std::vector<uint_t> r_p;
    /** [0..p][0..255] see the code to see how this variable is used */
    std::vector<std::vector<uint_t>> C;
    /** The disjoint interval sequence for LF */
    std::vector<std::pair<uint_t,uint_t>> I_LF;
    /** The disjoint interval sequence for Phi */
    std::vector<std::pair<uint_t,uint_t>> I_Phi;
    /** [0..r'-1] if the end position of the i-th input interval of M_LF is the end position of a BWT run, then
     * SA_[i] is the suffix array sample at the end position of the i-th input interval of M_LF; else SA_s[i] = n */
    std::vector<uint_t> SA_s;
    /** [0..r'-1] Permutation, that stores the order of the values in SA_s */
    std::vector<uint_t> pi_;
    /** [0..r''-1] Permutation, that stores the order of the output interval starting positions of M_Phi */
    std::vector<uint_t> pi_mphi;

    // ############################# COMMON MISC METHODS #############################

    /**
     * @brief reports the operations, that are supported by this index
     */
    void report_supports() {
        std::cout << "building an index with the following support:" << std::endl;

        static auto report_support = [](move_r_support sup){
            if (sup == revert) {
                std::cout << "revert";
            } else if (sup == move_r_support::count) {
                std::cout << "count";
            } else if (sup == move_r_support::locate) {
                std::cout << "locate";
            }
        };

        for (uint8_t i=0; i<idx.support.size()-1; i++) {
            report_support(idx.support[i]);
            std::cout << ", ";
        }

        report_support(idx.support.back());
        std::cout << std::endl;
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

    /**
     * @brief sets the run length of the i-th BWT run to len
     * @param i [0..r-1] run index
     * @param len run length
     */
    void set_run_length(uint_t i, uint32_t len) {
        RLBWT.set_unsafe<1,uint32_t>(i,len);
    }

    /**
     * @brief returns the length of the i-th BWT run
     * @param i [0..r-1] run index
     * @return run length
     */
    uint32_t run_length(uint_t i) {
        return RLBWT.get_unsafe<1,uint32_t>(i);
    }

    /**
     * @brief sets the character of the i-th BWT run to c
     * @param i [0..r-1] run index
     * @param c character
     */
    void set_run_char(uint_t i, char c) {
        RLBWT.set_unsafe<0,char>(i,c);
    }

    /**
     * @brief returns the character of the i-th BWT run
     * @param i [0..r-1] run index
     * @return character
     */
    char run_char(uint_t i) {
        return RLBWT.get_unsafe<0,char>(i);
    }

    uint8_t run_uchar(uint_t i) {
        return RLBWT.get_unsafe<0,uint8_t>(i);
    }

    /**
     * @brief sets some variables and logs
     */
    void prepare_phase_1() {
        time = now();
        time_start = time;
        omp_set_num_threads(p);

        if (log) {
            baseline_memory_allocation = malloc_count_current();
            malloc_count_reset_peak();
        }

        adjust_supports(support);

        idx.support = support;
        idx.a = a;

        if (!contains(idx.support,revert)) {
            std::cout << "error: cannot build an index without revert support";
            return;
        }

        build_count_support = contains(support,move_r_support::count);
        build_locate_support = contains(support,move_r_support::locate);

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
            p = std::max((uint_t)1,n/1000);
            if (log) std::cout << "warning: p > n/1000, setting p to n/1000 ~ " << std::to_string(p) << std::endl;
        }

        if (build_locate_support) {
            p_r = std::min((uint_t)256,std::max((uint_t)1,n/1000));
        } else {
            p_r = 1;
        }
        
        idx.p_r = p_r;
    }

    /**
     * @brief logs a construction summary
     */
    void log_finished() {
        uint64_t time_construction = time_diff_ns(time_start,now());
        uint64_t peak_memory_allocation = malloc_count_peak() - baseline_memory_allocation;
        
        std::cout << std::endl;
        std::cout << "construction time: " << format_time(time_construction) << std::endl;
        std::cout << "peak memory allocation: " << format_size(peak_memory_allocation) << std::endl;
        idx.log_data_structure_sizes();

        if (mf_idx != NULL) {
            *mf_idx << " time_construction=" << time_construction;
            *mf_idx << " peak_memory_allocation=" << peak_memory_allocation;
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
            *mf_idx << " sigma=" << std::to_string(sigma);
            *mf_idx << " r=" << r;
        }
        std::cout << "n = " << n << ", sigma = " << std::to_string(sigma) << ", r = " << r << ", n/r = " << n_r << std::endl;
    }

    // ############################# CONSTRUCTORS #############################

    /**
     * @brief constructs a move_r index of the string input
     * @param index The move-r index to construct
     * @param T the string containing T 
     * @param support a vector containing move_r operations to build support for
     * @param construction_mode cosntruction mode to use (default: optimized for low runtime)
     * @param p the number of threads to use during the construction
     * @param a balancing parameter, O(r*(a/(a-1))), 2 <= a
     * @param log controls, whether to print log messages
     * @param mf_idx measurement file for the index construciton
     * @param mf_mds measurement file for the move data structure construction
     * @param name_text_file name of the input file (used only for measurement output)
     */
    construction(
        move_r<uint_t>& index,
        std::string& T,
        bool delete_T,
        std::vector<move_r_support> support,
        move_r_construction_mode construction_mode,
        uint16_t p,
        uint16_t a,
        bool log,
        std::ostream* mf_idx,
        std::ostream* mf_mds,
        std::string name_text_file
    ) : T(T), L(L_tmp), SA_32(SA_32_tmp), SA_64(SA_64_tmp), idx(index) {
        this->delete_T = delete_T;
        this->support = support;
        this->p = p;
        this->a = a;
        this->log = log;
        this->mf_idx = mf_idx;
        this->mf_mds = mf_mds;
        this->name_text_file = name_text_file;

        prepare_phase_1();

        if (construction_mode == move_r_construction_mode::runtime) {
            min_valid_char = 2;
            terminator = 1;

            T.push_back(uchar_to_char((uint8_t)1));
            n = T.size();
            idx.n = n;
            construct_in_memory();
            T.resize(n-1);
            if (idx.chars_remapped) unmap_t();
        } else {
            min_valid_char = 3;
            terminator = 0;
            n = T.size()+1;
            idx.n = n;

            std::ifstream t_file = preprocess_and_store_t_in_file();
            construct_space_saving(t_file,true);
        }

        if (log) log_finished();
    }

    /**
     * @brief constructs a move_r index from an input file
     * @param index The move-r index to construct
     * @param t_file file containing T
     * @param support a vector containing move_r operations to build support for
     * @param construction_mode cosntruction mode to use (default: optimized for low runtime)
     * @param p the number of threads to use during the construction
     * @param a balancing parameter, O(r*(a/(a-1))), 2 <= a
     * @param log controls, whether to print log messages
     * @param mf_idx measurement file for the index construciton
     * @param mf_mds measurement file for the move data structure construction
     * @param name_text_file name of the input file (used only for measurement output)
     */
    construction(
        move_r<uint_t>& index,
        std::ifstream& t_file,
        std::vector<move_r_support> support,
        move_r_construction_mode construction_mode,
        uint16_t p,
        uint16_t a,
        bool log,
        std::ostream* mf_idx,
        std::ostream* mf_mds,
        std::string name_text_file
    ) : T(T_tmp), L(L_tmp), SA_32(SA_32_tmp), SA_64(SA_64_tmp), idx(index) {
        this->support = support;
        this->p = p;
        this->a = a;
        this->log = log;
        this->mf_idx = mf_idx;
        this->mf_mds = mf_mds;
        this->name_text_file = name_text_file;

        prepare_phase_1();

        if (construction_mode == move_r_construction_mode::runtime) {
            min_valid_char = 2;
            terminator = 1;

            read_t_from_file_in_memory(t_file);
            construct_in_memory();
        } else {
            min_valid_char = 3;
            terminator = 0;

            construct_space_saving(t_file,false);
        }

        if (log) log_finished();
    }

    /**
     * @brief constructs a move_r index from an input file
     * @param index The move-r index to construct
     * @param suffix_array vector containing the suffix array of the input
     * @param bwt string containing the bwt of the input
     * @param support a vector containing move_r operations to build support for
     * @param p the number of threads to use during the construction
     * @param a balancing parameter, O(r*(a/(a-1))), 2 <= a
     * @param log controls, whether to print log messages
     * @param mf_idx measurement file for the index construciton
     * @param mf_mds measurement file for the move data structure construction
     * @param name_text_file name of the input file (used only for measurement output)
     */
    construction(
        move_r<uint_t>& index,
        std::vector<int32_t>& suffix_array,
        std::string& bwt,
        std::vector<move_r_support> support,
        uint16_t p,
        uint16_t a,
        bool log,
        std::ostream* mf_idx,
        std::ostream* mf_mds,
        std::string name_text_file
    ) : T(T_tmp), L(bwt), SA_32(suffix_array), SA_64(SA_64_tmp), idx(index) {
        this->support = support;
        this->p = p;
        this->a = a;
        this->log = log;
        this->mf_idx = mf_idx;
        this->mf_mds = mf_mds;
        this->name_text_file = name_text_file;
        
        construct_from_sa_and_l<int32_t>();
    }

    /**
     * @brief constructs a move_r index from an input file
     * @param index The move-r index to construct
     * @param suffix_array vector containing the suffix array of the input
     * @param bwt string containing the bwt of the input
     * @param support a vector containing move_r operations to build support for
     * @param p the number of threads to use during the construction
     * @param a balancing parameter, O(r*(a/(a-1))), 2 <= a
     * @param log controls, whether to print log messages
     * @param mf_idx measurement file for the index construciton
     * @param mf_mds measurement file for the move data structure construction
     * @param name_text_file name of the input file (used only for measurement output)
     */
    construction(
        move_r<uint_t>& index,
        std::vector<int64_t>& suffix_array,
        std::string& bwt,
        std::vector<move_r_support> support,
        uint16_t p,
        uint16_t a,
        bool log,
        std::ostream* mf_idx,
        std::ostream* mf_mds,
        std::string name_text_file
    ) : T(T_tmp), L(bwt), SA_32(SA_32_tmp), SA_64(suffix_array), idx(index) {
        this->support = support;
        this->p = p;
        this->a = a;
        this->log = log;
        this->mf_idx = mf_idx;
        this->mf_mds = mf_mds;
        this->name_text_file = name_text_file;
        
        construct_from_sa_and_l<int64_t>();
    }

    // ############################# CONSTRUCTION #############################

    /**
     * @brief constructs the index from a suffix array and a bwt
     */
    template <typename sa_sint_t>
    void construct_from_sa_and_l() {
        build_from_sa_and_l = true;
        min_valid_char = 2;
        terminator = 1;
        n = L.size();
        idx.n = n;

        prepare_phase_1();
        prepare_phase_2();
        build_rlbwt_c_in_memory<sa_sint_t,true>();

        if (log) log_statistics();

        build_ilf();
        build_mlf();

        if (build_locate_support) build_iphi_in_memory<sa_sint_t>();

        build_l__sas();

        if (build_locate_support) {
            sort_iphi();
            build_mphi();
            build_saidx();
            build_de();
        }

        if (build_count_support) build_rsl_();

        if (log) log_finished();
    }

    /**
     * @brief constructs the index from a string in memory (optimized for low runtime)
     */
    void construct_in_memory() {
        if constexpr (std::is_same<uint_t,uint32_t>::value) {
            if (n <= INT_MAX) {
                construct_in_memory<int32_t>();
            } else {
                construct_in_memory<int64_t>();
            }
        } else {
            construct_in_memory<int64_t>();
        }
    }

    /**
     * @brief constructs the index from a string in memory (optimized for low runtime)
     * @param sa_sint_t signed integer type to use for the suffix array entries
     */
    template <typename sa_sint_t = int32_t>
    void construct_in_memory() {
        preprocess_t(true,true);
        prepare_phase_2();
        build_sa_in_memory<sa_sint_t>();
        build_rlbwt_c_in_memory<sa_sint_t,false>();

        if (log) log_statistics();

        build_ilf();
        build_mlf();
        
        if (build_locate_support) build_iphi_in_memory<sa_sint_t>();

        build_l__sas();

        if (build_locate_support) {
            sort_iphi();
            build_mphi();
            build_saidx();
            build_de();
        }

        if (build_count_support) build_rsl_();
    }

    /**
     * @brief constructs the index from a file space saving (optimized for low peak memory usage)
     * @param t_file file containing T
     * @param delete_t_file controls, whether t_file will be deleted as soon as it is not needed anymore
     */
    void construct_space_saving(std::ifstream& t_file, bool delete_t_file) {
        if (!delete_t_file) preprocess_t_buffered_from_file(t_file);

        prepare_phase_2();
        pfp(t_file,delete_t_file);

        if (log) {
            log_peak_mem_usage(baseline_memory_allocation);
            log_statistics();
        }

        read_rlbwt_bwt();

        if (read_rlbwt) {
            preprocess_rlbwt_space_saving();
        } else {
            build_rlbwt_c_in_memory<int32_t,true>();
        }
        
        build_ilf();
        build_mlf();
        
        if (build_locate_support) read_iphi_space_saving();

        build_l__sas();

        if (build_locate_support) {
            store_mlf();
            sort_iphi();
            build_mphi();
            build_saidx();
            build_de();
            load_mlf();
        }

        if (build_count_support) build_rsl_();
    };

    // ############################# COMMON CONSTRUCTION METHODS #############################

    /**
     * @brief reads the input T and possibly remaps it to an internal alphabet, if it contains an invalid character
     * @param t_file file containing T (for in_memory = false)
     */
    void preprocess_t(bool in_memory, bool map_t, std::ifstream* t_file = NULL);

    /**
     * @brief builds the RLBWT and C in-memory
     * @tparam sa_sint_t suffix array signed integer type
     * @tparam whether the RLBWT should be read from L
     */
    template <typename sa_sint_t = int32_t, bool read_l>
    void build_rlbwt_c_in_memory();

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
     * @brief builds L' (and SA_s)
     */
    void build_l__sas();

    /**
     * @brief sorts I_Phi
     */
    void sort_iphi();

    /**
     * @brief builds M_Phi
     */
    void build_mphi();

    /**
     * @brief builds SA_idx
     */
    void build_saidx();

    /**
     * @brief builds D_e
     */
    void build_de();

    /**
     * @brief builds RS_L'
     */
    void build_rsl_();

    // ############################# IN-MEMORY CONSTRUCTION METHODS #############################

    /**
     * @brief reads T from t_file
     * @param t_file file containing T
     */
    void read_t_from_file_in_memory(std::ifstream& t_file);

    /**
     * @brief builds the suffix array in-memory
     * @tparam sa_sint_t suffix array signed integer type
     */
    template <typename sa_sint_t = int32_t>
    void build_sa_in_memory();

    /**
     * @brief builds SA_s from SA in memory
     * @tparam sa_sint_t suffix array signed integer type
     */
    template <typename sa_sint_t>
    void build_iphi_in_memory();

    /**
     * @brief unmaps T from the internal alphabet
     */
    void unmap_t();

    // ############################# SPACE-SAVING CONSTRUCTION METHODS #############################

    /**
     * @brief preprocesses T and stores it in a file
     * @return file containing T
     */
    std::ifstream preprocess_and_store_t_in_file();

    /**
     * @brief preprocesses T buffered while reading it from a file
     * @param t_file 
     */
    void preprocess_t_buffered_from_file(std::ifstream& t_file);

    /**
     * @brief computes the (RL-)BWT (and I_Phi) using prefix-free-parsing
     * @param t_file file containing T
     * @param delete_t_file controls, whether to delete t_file as soon as it is not needed anymore
     */
    void pfp(std::ifstream& t_file, bool delete_t_file);

    /**
     * @brief reads the RLBWT or BWT from a file
     */
    void read_rlbwt_bwt();

    /**
     * @brief builds the RLBWT and C
     */
    void preprocess_rlbwt_space_saving();

    /**
     * @brief reads I_Phi
     */
    void read_iphi_space_saving();

    /**
     * @brief stores M_LF in a file
     */
    void store_mlf();

    /**
     * @brief loads M_LF from a file
     */
    void load_mlf();

    /**
     * @brief computes the missing SA-samples
     */
    void compute_missing_sa_samples_space_saving();
};

#include "modes/common.cpp"
#include "modes/in_memory.cpp"
#include "modes/space_saving.cpp"