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
    bool build_from_sa_and_l = false; // 
    bool log = false; // controls, whether to print log messages
    std::ostream* measurement_file_index = NULL; // file to write measurement data of the index construction to 
    std::ostream* measurement_file_move_data_structures = NULL; // file to write measurement data of the move data structure construction to 
    std::string name_textfile = ""; // name of the text file (only for measurement output)
    std::string prefix_tempfiles = ""; // prefix of temporary files
    std::chrono::steady_clock::time_point time; // time of the start of the last build phase
    std::chrono::steady_clock::time_point time_start; // time of the start of the whole build phase
    uint64_t baseline_memory_allocation = 0; // memory allocation at the start of the construction
    bool build_count_support = false; // = true <=> build support for count (RS_L')
    bool build_locate_support = false; // = true <=> build support for locate (SA_idx, SA_offs and M_Phi) and build parallel revert support (D_e)
    bool read_rle_bwt = false; // = true <=> read the run-length encoded BWT
    bool compress_br = false; // = true <=> compress B_r using sd-arrays
    uint8_t terminator = 0; // terminator (dollar) symbol
    uint8_t min_valid_char = 0; // the minimum valid character that is allowed to occur in T

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
    uint8_t omega_offs = 0; // word width of SA_offs

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
    /** [0..r-1] characters of the bwt runs */
    std::string bwt_run_heads;
    /** [0..r-1] lengths of the bwt runs */
    std::vector<uint32_t> bwt_run_lengths;
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
    /** [0..r'-1] Suffix array samples at the end positions of the input intervals in M_LF; SA_s[x] = SA[M_LF.p(x+1)-1] */
    std::vector<uint_t> SA_s;
    /** [0..r'-1] Permutation, that stores the order of the values in SA_s */
    std::vector<uint_t> pi_;
    /** [...] vector that contains tuples of (i,j,l), where i marks the first input interval in a run of input intervals
     * of M_LF, whiches end position does not correspond to the end position of a run in the bwt */
    std::vector<std::tuple<uint_t,uint_t,uint_t>> SA_s_missing;
    /** [0..p-1] vector of vectors like SA_s_missing, where the i-th vector corresponds to threads i */
    std::vector<std::vector<std::tuple<uint_t,uint_t,uint_t>>> SA_s_missing_thr;
    /** [0..p] start positions of each threads section in SA_s_missing */
    std::vector<uint_t> SA_s_missing_sect;
    /** [0..p] stores at position i the number of Phi-queries, that is required to compute all missing SA-samples in
     *  the section of thread i */
    std::vector<uint_t> num_SA_s_missing_thr;
    /** [0..n-1] bitvector that marks the bwt run start positions */
    std::vector<sdsl::bit_vector> B_r;
    /** [0..n-1] bitvector that marks the bwt run start positions */
    std::vector<sd_array<uint_t>> B_r_sd;

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
            omp_set_num_threads(p);
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

        if (measurement_file_index != NULL) {
            *measurement_file_index << " time_construction=" << time_construction;
            *measurement_file_index << " peak_memory_allocation=" << peak_memory_allocation;
            idx.log_data_structure_sizes(*measurement_file_index);
            *measurement_file_index << std::endl;
        }
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
     * @param measurement_file_index measurement file for the index construciton
     * @param measurement_file_move_data_structures measurement file for the move data structure construction
     * @param name_textfile name of the input file (used only for measurement output)
     */
    construction(
        move_r<uint_t>& index,
        std::string& T,
        std::vector<move_r_support> support,
        move_r_construction_mode construction_mode,
        uint16_t p,
        uint16_t a,
        bool log,
        std::ostream* measurement_file_index,
        std::ostream* measurement_file_move_data_structures,
        std::string name_textfile
    ) : T(T), L(L_tmp), SA_32(SA_32_tmp), SA_64(SA_64_tmp), idx(index) {
        this->support = support;
        this->p = p;
        this->a = a;
        this->log = log;
        this->measurement_file_index = measurement_file_index;
        this->measurement_file_move_data_structures = measurement_file_move_data_structures;
        this->name_textfile = name_textfile;

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
     * @param measurement_file_index measurement file for the index construciton
     * @param measurement_file_move_data_structures measurement file for the move data structure construction
     * @param name_textfile name of the input file (used only for measurement output)
     */
    construction(
        move_r<uint_t>& index,
        std::ifstream& t_file,
        std::vector<move_r_support> support,
        move_r_construction_mode construction_mode,
        uint16_t p,
        uint16_t a,
        bool log,
        std::ostream* measurement_file_index,
        std::ostream* measurement_file_move_data_structures,
        std::string name_textfile
    ) : T(T_tmp), L(L_tmp), SA_32(SA_32_tmp), SA_64(SA_64_tmp), idx(index) {
        this->support = support;
        this->p = p;
        this->a = a;
        this->log = log;
        this->measurement_file_index = measurement_file_index;
        this->measurement_file_move_data_structures = measurement_file_move_data_structures;
        this->name_textfile = name_textfile;

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
     * @tparam sa_sint_t suffix array signed integer type
     * @param index The move-r index to construct
     * @param suffix_array vector containing the suffix array of the input
     * @param bwt string containing the bwt of the input
     * @param support a vector containing move_r operations to build support for
     * @param p the number of threads to use during the construction
     * @param a balancing parameter, O(r*(a/(a-1))), 2 <= a
     * @param log controls, whether to print log messages
     * @param measurement_file_index measurement file for the index construciton
     * @param measurement_file_move_data_structures measurement file for the move data structure construction
     * @param name_textfile name of the input file (used only for measurement output)
     */
    construction(
        move_r<uint_t>& index,
        std::vector<int32_t>& suffix_array,
        std::string& bwt,
        std::vector<move_r_support> support,
        uint16_t p,
        uint16_t a,
        bool log,
        std::ostream* measurement_file_index,
        std::ostream* measurement_file_move_data_structures,
        std::string name_textfile
    ) : T(T_tmp), L(bwt), SA_32(suffix_array), SA_64(SA_64_tmp), idx(index) {
        this->support = support;
        this->p = p;
        this->a = a;
        this->log = log;
        this->measurement_file_index = measurement_file_index;
        this->measurement_file_move_data_structures = measurement_file_move_data_structures;
        this->name_textfile = name_textfile;
        
        construct_from_sa_and_l<int32_t>();
    }

    /**
     * @brief constructs a move_r index from an input file
     * @tparam sa_sint_t suffix array signed integer type
     * @param index The move-r index to construct
     * @param suffix_array vector containing the suffix array of the input
     * @param bwt string containing the bwt of the input
     * @param support a vector containing move_r operations to build support for
     * @param p the number of threads to use during the construction
     * @param a balancing parameter, O(r*(a/(a-1))), 2 <= a
     * @param log controls, whether to print log messages
     * @param measurement_file_index measurement file for the index construciton
     * @param measurement_file_move_data_structures measurement file for the move data structure construction
     * @param name_textfile name of the input file (used only for measurement output)
     */
    construction(
        move_r<uint_t>& index,
        std::vector<int64_t>& suffix_array,
        std::string& bwt,
        std::vector<move_r_support> support,
        uint16_t p,
        uint16_t a,
        bool log,
        std::ostream* measurement_file_index,
        std::ostream* measurement_file_move_data_structures,
        std::string name_textfile
    ) : T(T_tmp), L(bwt), SA_32(SA_32_tmp), SA_64(suffix_array), idx(index) {
        this->support = support;
        this->p = p;
        this->a = a;
        this->log = log;
        this->measurement_file_index = measurement_file_index;
        this->measurement_file_move_data_structures = measurement_file_move_data_structures;
        this->name_textfile = name_textfile;
        
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
        build_l_and_c_in_memory<sa_sint_t>();
        process_c_array();
        process_rp_in_memory();
        build_ilf_iphi_and_bwt_run_heads_in_memory<sa_sint_t>();
        build_mlf();
        build_l__in_memory_from_l<sa_sint_t>();

        if (build_locate_support) {
            sort_iphi();
            build_mphi();
            build_sas_in_memory<sa_sint_t>();
            build_saidxoffs(r_);
            build_de(r_);
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
        build_l_and_c_in_memory<sa_sint_t>();
        process_c_array();
        process_rp_in_memory();
        build_ilf_iphi_and_bwt_run_heads_in_memory<sa_sint_t>();
        build_br_in_memory();
        build_mlf();
        build_l__and_iphi_in_memory<sa_sint_t>();

        if (build_locate_support) {
            build_sas_in_memory<sa_sint_t>();
            sort_iphi();
            build_mphi();
            build_saidxoffs(r_);
            build_de(r_);
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
        if (log) log_peak_mem_usage(baseline_memory_allocation);
        build_rlbwt_and_c_space_saving();
        process_c_array();
        build_ilf_space_saving();
        build_mlf();

        if (build_locate_support) read_iphi_space_saving();
        build_l__and_sas_space_saving();

        if (build_locate_support) {
            store_mlf();
            sort_iphi();
            build_mphi();
            build_saidxoffs(r);
            build_de(r);
            load_mlf();
            compute_missing_sa_samples_space_saving();
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
     * @brief processes the C-array
     */
    void process_c_array();

    /**
     * @brief builds M_LF
     */
    void build_mlf();

    /**
     * @brief sorts I_Phi
     */
    void sort_iphi();

    /**
     * @brief builds M_Phi
     */
    void build_mphi();

    /**
     * @brief builds SA_idxoffs
     * @param r_ r' allows to override r' := r during memory-saving construction
     */
    void build_saidxoffs(uint_t r_);

    /**
     * @brief builds D_e
     * @param r_ r' allows to override r' := r during memory-saving construction
     */
    void build_de(uint_t r_);

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
     * @brief builds L and C in-memory
     * @tparam sa_sint_t suffix array signed integer type
     */
    template <typename sa_sint_t = int32_t>
    void build_l_and_c_in_memory();

    /**
     * @brief processes r_p in-memory
     */
    void process_rp_in_memory();

    /**
     * @brief builds I_LF and the BWT run heads in-memory
     * @tparam sa_sint_t suffix array signed integer type
     */
    template <typename sa_sint_t>
    void build_ilf_iphi_and_bwt_run_heads_in_memory();

    /**
     * @brief builds B_r in-memory
     */
    void build_br_in_memory();

    /**
     * @brief builds L' and I_Phi in-memory
     * @tparam sa_sint_t suffix array signed integer type
     */
    template <typename sa_sint_t = int32_t>
    void build_l__and_iphi_in_memory();

    /**
     * @brief builds L' in-memory from L
     * @tparam sa_sint_t suffix array signed integer type 
     */
    template <typename sa_sint_t>
    void build_l__in_memory_from_l();

    /**
     * @brief builds SA_s in-memory
     * @tparam sa_sint_t suffix array signed integer type
     */
    template <typename sa_sint_t = int32_t>
    void build_sas_in_memory();

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
     * @brief builds the RLBWT and C
     */
    void build_rlbwt_and_c_space_saving();

    /**
     * @brief builds I_LF
     */
    void build_ilf_space_saving();

    /**
     * @brief reads I_Phi
     */
    void read_iphi_space_saving();

    /**
     * @brief builds L' (and SA_s)
     */
    void build_l__and_sas_space_saving();

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