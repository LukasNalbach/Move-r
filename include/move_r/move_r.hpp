#pragma once

#include <iostream>
#include <omp.h>
#include <move_r/misc/utils.hpp>
#include <move_r/data_structures/hybrid_bit_vector.hpp>
#include <move_r/data_structures/string_rank_select_support.hpp>
#include <move_r/data_structures/interleaved_vectors.hpp>
#include <move_r/data_structures/move_data_structure/move_data_structure.hpp>
#include <move_r/data_structures/move_data_structure/move_data_structure_str.hpp>

/**
 * @brief an operation that can be supported by a move_r object
 */
enum move_r_support {
    /* support for retrieving (a range of) the input string from the index (reverting
       the index); this also includes support for accessing or retrieving (a range in)
       the bwt (reverting in parallel requires the index to be built with locate support) */
    revert = 0,
    count = 1, // support for counting the occurrences of a pattern in the input string
    /* support for calculating the positions of occurrences of a pattern in the input
       string; this also adds support for accessing and retrieving (a range in)
       the suffix array (also in parallel) */
    locate = 2
};

/**
 * @brief a vector containing all operations supported by move_r
 */
static std::vector<move_r_support> full_support = {
    move_r_support::revert,
    move_r_support::count,
    move_r_support::locate
};

/**
 * @brief move-r construction mode
 */
enum move_r_construction_mode {
    space = 0, // optimized for low peak memory consumption (and low-runtime if the input is repetitive)
    runtime = 1 // optimized for low runtime if the input is non-repetitive
};

/**
 * @brief move-r index, size O(r*(a/(a-1)))
 * @tparam uint_t index integer type (use uint32_t if input size < UINT_MAX, else uint64_t)
 */
template <typename uint_t = uint32_t>
class move_r {
    static_assert(std::is_same<uint_t,uint32_t>::value || std::is_same<uint_t,uint64_t>::value);
    protected:

    // ############################# INDEX VARIABLES #############################

    uint_t n = 0; // the length of T
    uint8_t sigma = 0; // the number of distinct characters (including the terminator symbol 1) of T
    uint_t r = 0; // r, the number of runs in L
    uint_t r_ = 0; // r', the number of input/output intervals in M_LF
    uint_t r__ = 0; // r'', the number of input/output intervals in M_Phi
    uint16_t a = 0; // balancing parameter, restricts size to O(r*(a/(a-1))), 2 <= a
    uint16_t p_r = 0; // maximum possible number of threads to use while reverting the index
    uint8_t omega_idx = 0; // word width of SA_phi

    std::vector<move_r_support> support; // contains all supported operations
    /* true <=> the characters of the input string have been remapped internally, because the input
       contains 0 and/or 1*/
    bool chars_remapped = false;

    // ############################# INDEX DATA STRUCTURES #############################

    /* mapping function from the alphabet of T to the internal effective alphabet (which starts with 1) */
    std::vector<uint8_t> map_char;
    std::vector<uint8_t> unmap_char; // inverse function of map_char
    /* The Move Data Structure for LF. It also stores L', which can be accessed at
    position i with M_LF.L_(i). */
    move_data_structure_str<uint_t> M_LF;
    /* [0..p_r-1], where D_e[i] = <x,j>, x in [0,r'-1] and j is minimal, s.t. SA_s[x]=j > i* lfloor (n-1)/p rfloor;
    see the parallel revert algorithm to understand why this is useful. */
    std::vector<std::pair<uint_t,uint_t>> D_e;
    string_rank_select_support<uint_t> RS_L_; // rank-select data structure for L'
    move_data_structure<uint_t> M_Phi; // The Move Data Structure for Phi.
    interleaved_vectors<uint_t> SA_phi; // [0..r'-1] SA_phi

    // ############################# INTERNAL METHODS #############################

    /**
     * @brief returns SA_s[x]
     * @param x [0..r'-1] the end position of the x-th input interval in M_LF must be an end position of a bwt run
     * @return SA_s[x]
     */
    inline uint_t SA_s(uint_t x) {
        return M_Phi.q(SA_phi[x]);
    }

    /**
     * @brief sets SA_phi[x] to idx
     * @param x [0..r-1]
     * @param idx [0..r''-1]
     */
    inline void set_SA_phi(uint_t x, uint_t idx) {
        SA_phi.template set<0>(x,idx);
    }

    /**
     * @brief maps a character that occurs in the input string to the internal alphabet, if the characters
     *        have been remapped
     * @param c a character that occurs in the input string
     * @return the character in the internal effective alphabet, to which c has been mapped to
     */
    inline char map_to_internal(char c) {
        return uchar_to_char(map_char[char_to_uchar(c)]);
    }

    /**
     * @brief maps a character that occurs in the internal effective alphabet to its corresponding
     *          character in the input string
     * @param c a character that occurs in the internal effective alphabet
     * @return its corresponding character in the input string
     */
    inline char unmap_from_internal(char c) {
        return uchar_to_char(unmap_char[char_to_uchar(c)]);
    }

    /**
     * @brief returns L'[x]
     * @param x [0..r'-1]
     * @return L'[x]
     */
    inline char L_(uint_t x) {
        return M_LF.character(x);
    }

    /**
     * @brief Sets the up a Phi-move-pair for the suffix array sample at the end position of the x-th input interval in M_LF
     * @param x an input interval in M_LF (the end position of the x-th input interval in M_LF must be an end position of a BWT run)
     * @param s variable to store the suffix array sample at position l'_{x+1}-1
     * @param s_ variable to store the index of the input interval in M_Phi containing s
     */
    void setup_phi_move_pair(uint_t x, uint_t& s, uint_t& s_);

    /**
     * @brief adds implicitly supported move_r operations to the vector support and sorts it afterwards
     * @param support a vector containing move_r operations
     */
    static void adjust_supports(std::vector<move_r_support>& support) {
        if (contains(support,move_r_support::locate) && !(contains(support,move_r_support::count))) {
            support.emplace_back(move_r_support::count);
        }
        
        if (contains(support,move_r_support::count) && !(contains(support,move_r_support::revert))) {
            support.emplace_back(move_r_support::revert);
        }
        
        ips4o::sort(support.begin(),support.end());
        support.erase(std::unique(support.begin(),support.end()),support.end());
    }

    /**
     * @brief executes retrieve_method with the parameters l, r and num_threads, buffers the output in num_threads
     * buffers and num_threads temporary files and then writes the temporary files into the file out
     * @tparam output_t type of the output data
     * @tparam output_reversed controls, whether the output should be reversed
     * @param retrieve_method function, whiches output should be buffered
     * @param out file to write the output to
     * @param l left range limit
     * @param r right range limit
     * @param num_threads maximum number of threads to use
     * @param max_bytes_alloc maximum number of bytes to allocate
     */
    template <typename output_t, bool output_reversed>
    void retrieve_range(
        void(move_r<uint_t>::*retrieve_method)(const std::function<void(uint_t,output_t)>&,uint_t,uint_t,uint16_t),
        std::ofstream& out, uint_t l, uint_t r, uint16_t num_threads, int64_t max_bytes_alloc
    );

    class construction;

    // ############################# CONSTRUCTORS #############################

    public:
    move_r() = default;

    /**
     * @brief constructs a move_r index of the string input
     * @param input the input string
     * @param support a vector containing move_r operations to build support for
     * @param construction_mode cosntruction mode to use (default: optimized for low runtime)
     * @param num_threads maximum number of threads to use during the construction
     * @param a balancing parameter, 2 <= a
     * @param log controls, whether to print log messages
     * @param mf_idx measurement file for the index construciton
     * @param mf_mds measurement file for the move data structure construction
     * @param name_text_file name of the input file (used only for measurement output)
     */
    move_r(
        std::string& input,
        std::vector<move_r_support> support = full_support,
        move_r_construction_mode construction_mode = move_r_construction_mode::runtime,
        uint16_t num_threads = omp_get_max_threads(),
        uint16_t a = 8,
        bool log = false,
        std::ostream* mf_idx = NULL,
        std::ostream* mf_mds = NULL,
        std::string name_text_file = ""
    ) {
        construction mrc(*this,input,false,support,construction_mode,num_threads,a,log,mf_idx,mf_mds,name_text_file);
    }

    /**
     * @brief constructs a move_r index of the string input
     * @param input the input string
     * @param support a vector containing move_r operations to build support for
     * @param construction_mode cosntruction mode to use (default: optimized for low runtime)
     * @param num_threads maximum number of threads to use during the construction
     * @param a balancing parameter, 2 <= a
     * @param log controls, whether to print log messages
     * @param mf_idx measurement file for the index construciton
     * @param mf_mds measurement file for the move data structure construction
     * @param name_text_file name of the input file (used only for measurement output)
     */
    move_r(
        std::string&& input,
        std::vector<move_r_support> support = full_support,
        move_r_construction_mode construction_mode = move_r_construction_mode::runtime,
        uint16_t num_threads = omp_get_max_threads(),
        uint16_t a = 8,
        bool log = false,
        std::ostream* mf_idx = NULL,
        std::ostream* mf_mds = NULL,
        std::string name_text_file = ""
    ) {
        construction mrc(*this,input,true,support,construction_mode,num_threads,a,log,mf_idx,mf_mds,name_text_file);
    }

    /**
     * @brief constructs a move_r index from an input file
     * @param input_file input file
     * @param support a vector containing move_r operations to build support for
     * @param construction_mode cosntruction mode to use (default: optimized for low runtime)
     * @param num_threads maximum number of threads to use during the construction
     * @param a balancing parameter, 2 <= a
     * @param log controls, whether to print log messages
     * @param mf_idx measurement file for the index construciton
     * @param mf_mds measurement file for the move data structure construction
     * @param name_text_file name of the input file (used only for measurement output)
     */
    move_r(
        std::ifstream& input_file,
        std::vector<move_r_support> support = full_support,
        move_r_construction_mode construction_mode = move_r_construction_mode::runtime,
        uint16_t num_threads = omp_get_max_threads(),
        uint16_t a = 8,
        bool log = false,
        std::ostream* mf_idx = NULL,
        std::ostream* mf_mds = NULL,
        std::string name_text_file = ""
    ) {
        construction mrc(*this,input_file,support,construction_mode,num_threads,a,log,mf_idx,mf_mds,name_text_file);
    }

    /**
     * @brief constructs a move_r index from a suffix array and a bwt
     * @tparam sa_sint_t suffix array signed integer type
     * @param suffix_array vector containing the suffix array of the input
     * @param bwt string containing the bwt of the input
     * @param support a vector containing move_r operations to build support for
     * @param num_threads maximum number of threads to use during the construction
     * @param a balancing parameter, 2 <= a
     * @param log controls, whether to print log messages
     * @param mf_idx measurement file for the index construciton
     * @param mf_mds measurement file for the move data structure construction
     * @param name_text_file name of the input file (used only for measurement output)
     */
    template <typename sa_sint_t>
    move_r(
        std::vector<sa_sint_t>& suffix_array,
        std::string& bwt,
        std::vector<move_r_support> support = full_support,
        uint16_t num_threads = omp_get_max_threads(),
        uint16_t a = 8,
        bool log = false,
        std::ostream* mf_idx = NULL,
        std::ostream* mf_mds = NULL,
        std::string name_text_file = ""
    ) {
        construction mrc(*this,suffix_array,bwt,support,num_threads,a,log,mf_idx,mf_mds,name_text_file);
    }

    // ############################# MISC PUBLIC METHODS #############################

    /**
     * @brief returns the size of the input string
     * @return size of the input string
     */
    uint_t input_size() {
        return n-1;
    }

    /**
     * @brief returns the number of distinct characters in the input string (alphabet size)
     * @return alphabet_size 
     */
    uint8_t alphabet_size() {
        return sigma-1;
    }

    /**
     * @brief returns the number of runs in the bwt (+ up to p, where p is the number of threads,
     * that where used during the construction of the index)
     * @return number of runs in the bwt 
     */
    uint_t num_bwt_runs() {
        return r;
    }

    /**
     * @brief returns the number of intervals in M_LF
     * @return number of intervals in M_LF
     */
    uint_t num_intervals_m_lf() {
        return r_;
    }

    /**
     * @brief returns the number of intervals in M_Phi
     * @return number of intervals in M_Phi
     */
    uint_t num_intervals_m_phi() {
        return r__;
    }

    /**
     * @brief returns the balancing parameter the index has been built with
     * @return balancing parameter 
     */
    uint16_t balancing_parameter() {
        return a;
    }

    /**
     * @brief returns the number omega_idx of bits used by one entry in SA_phi (word width of SA_phi)
     * @return omega_idx
     */
    inline uint8_t width_saidx() {
        return omega_idx;
    }

    /**
     * @brief returns the maximum number of threads that can be used to revert the index
     * @return maximum number of threads that can be used to revert the index 
     */
    uint16_t max_revert_threads() {
        return p_r;
    }

    /**
     * @brief returns a vector containing the supported operations of this index
     * @return vector containing the operations
     */
    std::vector<move_r_support> supported_operations() {
        return support;
    }

    /**
     * @brief returns whether the provided operation is supported by this index
     * @param operation a move_r operation
     * @return whether operation is supported by this index
     */
    bool does_support(move_r_support operation) {
        return contains(support,operation);
    }

    /**
     * @brief returns the size of the data structure in bytes
     * @return size of the data structure in bytes
     */
    uint64_t size_in_bytes() {
        return
            support.size()*sizeof(move_r_support)+ // variables
            4*sizeof(uint_t)+3+2*sizeof(uint16_t)+ // ...
            chars_remapped*2*256+ // mapchar & unmap_char
            p_r*sizeof(uint_t)+ // D_e
            M_LF.size_in_bytes()+ // M_LF and L'
            RS_L_.size_in_bytes()+ // RS_L'
            M_Phi.size_in_bytes()+ // M_Phi
            SA_phi.size_in_bytes(); // SA_phi
    }

    /**
     * @brief logs the index data structure sizes to cout
     */
    void log_data_structure_sizes() {
        std::cout << "index size: " << format_size(size_in_bytes()) << std::endl;

        std::cout << "M_LF: " << format_size(M_LF.size_in_bytes()-(r_+1)) << std::endl;
        std::cout << "L': " << format_size(r_+1) << std::endl;
        
        if (does_support(move_r_support::count)) {
            std::cout << "RS_L': " << format_size(RS_L_.size_in_bytes()) << std::endl;

            if (does_support(move_r_support::locate)) {
                std::cout << "M_Phi: " << format_size(M_Phi.size_in_bytes()) << std::endl;
                std::cout << "SA_phi: " << format_size(SA_phi.size_in_bytes()) << std::endl;
            }
        }
    }

    /**
     * @brief logs the index data structure sizes to the output stream out
     * @param out an output stream
     */
    void log_data_structure_sizes(std::ostream& out) {
        out << " size_index=" << size_in_bytes();
        out << " size_m_lf=" << M_LF.size_in_bytes()-(r_+1);
        out << " size_l_=" << r_+1;
        
        if (does_support(move_r_support::count)) {
            out << " size_rs_l_=" << RS_L_.size_in_bytes();

            if (does_support(move_r_support::locate)) {
                out << " size_m_phi=" << M_Phi.size_in_bytes();
                out << " size_sa_idx=" << SA_phi.size_in_bytes();
            }
        }
    }

    // ############################# PUBLIC ACCESS METHODS #############################

    /**
     * @brief returns L'[x]
     * @param x [0..num_intervals_m_lf()-1]
     * @return L'[x]
     */
    inline char access_l_(uint_t x) {
        return chars_remapped ? unmap_from_internal(L_(x)) : L_(x);
    }

    /**
     * @brief returns L[i], where $ = 0, so if the input contained 0, the output is not equal to the real bwt
     * @param x [0..input size]
     * @return L[i]
     */
    char access_bwt(uint_t i);

    /**
     * @brief returns SA[i]
     * @param x [0..input size]
     * @return SA[i]
     */
    uint_t access_sa(uint_t i);

    /**
     * @brief returns D_e[i]
     * @param i [0..max_revert_threads()-2]
     * @return D_e[i]
     */
    uint_t access_d_e(uint16_t i) {
        return D_e[i];
    }

    /**
     * @brief returns a reference to M_LF
     * @return M_LF
     */
    const move_data_structure_str<uint_t>& m_lf() {
        return M_LF;
    }

    /**
     * @brief returns a reference to M_Phi
     * @return M_Phi
     */
    const move_data_structure<uint_t>& m_phi() {
        return M_Phi;
    }

    /**
     * @brief returns a reference to RS_L'
     * @return RS_L'
     */
    const string_rank_select_support<uint_t>& rs_l_() {
        return RS_L_;
    }

    // ############################# QUERY METHODS #############################

    /**
     * @brief returns the number of occurrences of P in the input string
     * @param P the pattern to count in the input string
     * @return the number of occurrences of P in the input string
     */
    uint_t count(const std::string& P) {
        return count(P.size(),[&P](uint_t i){return P[i];});
    }

    /**
     * @brief returns the number of occurrences of a pattern in the input string
     * @param pattern_length length of the query pattern
     * @param read read(i) must return the character of the query pattern at position i, for each i \in [0,pattern_length)
     * @return the number of occurrences of the pattern in the input string
     */
    uint_t count(uint_t pattern_length, const std::function<char(uint_t)>& read);

    /**
     * @brief locates the pattern P in the input string
     * @param P the pattern to locate in the input string
     * @return a vector containing the occurrences of P in the input string
     */
    std::vector<uint_t> locate(const std::string& P) {
        std::vector<uint_t> Occ;
        locate(P,Occ);
        return Occ;
    }

    /**
     * @brief locates the pattern P in the input string and appends the positions of the occurrences to Occ
     * @param P the pattern to locate in the input string
     */
    void locate(const std::string& P, std::vector<uint_t>& Occ) {        
        locate(P,[&Occ](uint_t o){Occ.emplace_back(o);});
    }

    /**
     * @brief locates the pattern P in the input string
     * @param P the pattern to locate in the input string
     * @param report function that is called with every occurrence of P in the input string as a parameter
     */
    void locate(const std::string& P, const std::function<void(uint_t)>& report) {
        locate(P.size(),[&P](uint_t i){return P[i];},report);
    }

    /**
     * @brief locates a pattern in the input string
     * @param pattern_length length of the query pattern
     * @param read read(i) must return the character of the query pattern at position i, for each i \in [0,pattern_length)
     * @param report function that is called with every occurrence of the pattern in the input string as a parameter
     */
    void locate(uint_t pattern_length, const std::function<char(uint_t)>& read, const std::function<void(uint_t)>& report);

    // ############################# RETRIEVE-RANGE METHODS #############################

    /**
     * @brief returns the bwt in the range [l,r] (0 <= l <= r <= input size), else
     * if l > r, then the whole bwt is returned (default); $ = 0, so if the input contained 0, the output is not
     * equal to the real bwt
     * @param report function that is called with every tuple (i,c) as a parameter, where i in [l,r] and c = L[i]
     * @param l left range limit
     * @param r right range limit
     * @param num_threads maximum number of threads to use
     */
    std::string retrieve_bwt_range(uint_t l = 1, uint_t r = 0, uint16_t num_threads = omp_get_max_threads()) {
        r = std::max(r,n-1);

        if (l > r) {
            l = 0;
            r = n-1;
        }

        std::string L;
        no_init_resize(L,r-l+1);
        retrieve_bwt_range([&L,&l](uint_t i, char c){L[i-l] = c;},l,r,num_threads);
        return L;
    }

    /**
     * @brief reports the characters in the bwt in the range [l,r] (0 <= l <= r <= input size), else if l > r, then all
     * characters of the bwt are reported (default); $ = 0, so if the input contained 0, the output is not equal to the real bwt
     * @param report function that is called with every tuple (i,c) as a parameter, where i in [l,r] and c = L[i]
     * @param l left range limit
     * @param r right range limit
     * @param num_threads maximum number of threads to use (if num_threads = 1, then the values are reported from left to right, if num_threads > 1,
     * the order may vary)
     */
    void retrieve_bwt_range(const std::function<void(uint_t,char)>& report, uint_t l = 1, uint_t r = 0, uint16_t num_threads = omp_get_max_threads());

    /**
     * @brief writes the characters in the bwt in the range [l,r] blockwise to the file out (0 <= l <= r <= input size), else if
     * l > r, then the whole bwt is written (default); $ = 0, so if the input contained 0, the output is not equal to the real bwt
     * @param out file to write the bwt to
     * @param l left range limit
     * @param r right range limit
     * @param num_threads maximum number of threads to use
     * @param max_bytes_alloc maximum number of bytes to allocate (default (if set to -1): ~ (r-l+1)/500)
     */
    void retrieve_bwt_range(std::ofstream& out, uint_t l = 1, uint_t r = 0, uint16_t num_threads = omp_get_max_threads(), int64_t max_bytes_alloc = -1) {
        retrieve_range<char,false>(&move_r<uint_t>::retrieve_bwt_range,out,l,r,num_threads,max_bytes_alloc);
    }

    /**
     * @brief returns the input in the range [l,r] (0 <= l <= r < input size), else
     * if l > r, then the whole input is returned (default)
     * @param l left range limit
     * @param r right range limit
     * @param num_threads maximum number of threads to use
     */
    std::string revert_range(uint_t l = 1, uint_t r = 0, uint16_t num_threads = omp_get_max_threads()) {
        r = std::max(r,n-2);

        if (l > r) {
            l = 0;
            r = n-2;
        }

        std::string T;
        no_init_resize(T,r-l+1);
        revert_range([&T,&l](uint_t i, char c){T[i-l] = c;},l,r,num_threads);
        return T;
    }

    /**
     * @brief writes the input in the range [l,r] to the range [b,e] in string;
     * (0 <= l <= r < input size), else if l > r, then l := 0 and r := input size-1 (default);
     * (0 <= b <= e), else b := 0 and e := r-l (default);
     * (e < string.size());
     * (r-l = e-b)
     * @param l left range limit of the range in the input to revert
     * @param r right range limit of the range in the input to revert
     * @param b left range limit of the output range in string 
     * @param e right range limit of the output range in string 
     * @param num_threads maximum number of threads to use
     */
    void revert_range(std::string& string, uint_t l = 1, uint_t r = 0, uint_t b = 1, uint_t e = 0, uint16_t num_threads = omp_get_max_threads()) {
        r = std::max(r,n-2);

        if (l > r) {
            l = 0;
            r = n-2;
        }

        if (b > e) {
            b = 0;
            e = r-l;
        }

        if (e-b != r-l) {
            std::cout << "error: e-b != r-l";
            return;
        }

        if (e >= string.size()) {
            std::cout << "error: e >= string.size()";
            return;
        }
        
        revert_range([&string,&b](uint_t i, char c){string[i-b] = c;},l,r,num_threads);
    }

    /**
     * @brief reports the characters in the input in the range [l,r] (0 <= l <= r < input size), else if l > r, then
     * all characters of the input are reported (default)
     * @param report function that is called with every tuple (i,c) as a parameter, where i in [l,r] and c = T[i]
     * @param l left range limit
     * @param r right range limit
     * @param num_threads maximum number of threads to use (if num_threads = 1, then the values are reported from right to left, if num_threads > 1,
     * the order may vary)
     */
    void revert_range(const std::function<void(uint_t,char)>& report, uint_t l = 1, uint_t r = 0, uint16_t num_threads = omp_get_max_threads());

    /**
     * @brief reverts the input in the range [l,r] blockwise and writes it to the file out (0 <= l <= r < input size),
     * else if l > r, then the whole input is reverted (default)
     * @param out file to write the reverted input to
     * @param l left range limit
     * @param r right range limit
     * @param num_threads maximum number of threads to use
     * @param max_bytes_alloc maximum number of bytes to allocate (default (if set to -1): ~ (r-l+1)/500)
     */
    void revert_range(std::ofstream& out, uint_t l = 1, uint_t r = 0, uint16_t num_threads = omp_get_max_threads(), int64_t max_bytes_alloc = -1) {
        retrieve_range<char,true>(&move_r<uint_t>::revert_range,out,l,r,num_threads,max_bytes_alloc);
    }
    
    /**
     * @brief rebuilds and returns the suffix array in the range [l,r] (0 <= l <= r <= input size),
     * else if l > r, then the whole suffix array is rebuilt (default)
     * @param l left range limit
     * @param r right range limit
     * @param num_threads maximum number of threads to use
     * @return the suffix array
     */
    std::vector<uint_t> retrieve_sa_range(uint_t l = 1, uint_t r = 0, uint16_t num_threads = omp_get_max_threads()) {
        r = std::max(r,n-1);

        if (l > r) {
            l = 0;
            r = n-1;
        }

        std::vector<uint_t> SA;
        no_init_resize(SA,r-l+1);
        retrieve_sa_range([&SA,&l](uint_t i, uint_t s){SA[i-l] = s;},l,r,num_threads);
        return SA;
    }

    /**
     * @brief reports the suffix array values in the range [l,r] (0 <= l <= r <= input size),
     * else if l > r, then the whole suffix array is reported (default)
     * @param report function that is called with every tuple (i,s) as a parameter, where i in [l,r] and s = SA[i]
     * @param l left range limit
     * @param r right range limit
     * @param num_threads maximum number of threads to use (if num_threads = 1, then the values are reported from right to left,
     * if num_threads > 1, the order may vary)
     */
    void retrieve_sa_range(const std::function<void(uint_t,uint_t)>& report, uint_t l = 1, uint_t r = 0, uint16_t num_threads = omp_get_max_threads());

    /**
     * @brief writes the values in the suffix array of the input in the range [l,r] blockwise to the file out (0 <= l <= r <= input size),
     * else if l > r, then the whole suffix array is written (default)
     * @param out file to write the suffix array to
     * @param l left range limit
     * @param r right range limit
     * @param num_threads maximum number of threads to use
     * @param max_bytes_alloc maximum number of bytes to allocate (default (if set to -1): ~ (r-l+1)/500)
     */
    void retrieve_sa_range(std::ofstream& out, uint_t l = 1, uint_t r = 0, uint16_t num_threads = omp_get_max_threads(), int64_t max_bytes_alloc = -1) {
        retrieve_range<uint_t,true>(&move_r<uint_t>::retrieve_sa_range,out,l,r,num_threads,max_bytes_alloc);
    }

    // ############################# SERIALIZATION METHODS #############################

    /**
     * @brief stores the index to an output stream
     * @param out output stream to store the index to
     * @param support supported operations to store data structures for
     */
    void serialize(std::ostream& out, std::vector<move_r_support> support = full_support) {
        adjust_supports(support);

        if (!is_subset_of(support,this->support)) {
            std::cout << "error: cannot store an index with support it has not been built with" << std::flush;
            return;
        }

        bool is_64_bit = std::is_same<uint_t,uint64_t>::value;
        out.write((char*)&is_64_bit,1);

        std::streampos pos_data_structure_offsets = out.tellp();
        out.seekp(pos_data_structure_offsets+(std::streamoff)sizeof(std::streamoff),std::ios::beg);

        uint8_t num_supports = support.size();
        out.write((char*)&num_supports,1);
        out.write((char*)&support[0],num_supports*sizeof(move_r_support));

        out.write((char*)&n,sizeof(uint_t));
        out.write((char*)&sigma,1);
        out.write((char*)&r,sizeof(uint_t));
        out.write((char*)&r_,sizeof(uint_t));
        out.write((char*)&a,sizeof(uint16_t));
        out.write((char*)&p_r,sizeof(uint16_t));

        if (p_r > 0) {
            out.write((char*)&D_e[0],(p_r-1)*2*sizeof(uint_t));
        }

        out.write((char*)&chars_remapped,1);
        if (chars_remapped) {
            out.write((char*)&map_char[0],256);
            out.write((char*)&unmap_char[0],256);
        }

        M_LF.serialize(out);

        if (contains(support,move_r_support::count)) {
            RS_L_.serialize(out);
        }

        if (contains(support,move_r_support::locate)) {
            out.write((char*)&r__,sizeof(uint_t));
            M_Phi.serialize(out);

            out.write((char*)&omega_idx,1);
            SA_phi.serialize(out);
        }

        std::streamoff offs_end = out.tellp()-pos_data_structure_offsets;
        out.seekp(pos_data_structure_offsets,std::ios::beg);
        out.write((char*)&offs_end,sizeof(std::streamoff));
        out.seekp(pos_data_structure_offsets+offs_end,std::ios::beg);
    }

    /**
     * @brief reads a serialized index from an input stream
     * @param in an input stream storing a serialized index
     * @param support supported operations to load data structures for
     */
    void load(std::istream& in, std::vector<move_r_support> support = full_support) {
        adjust_supports(support);

        bool is_64_bit;
        in.read((char*)&is_64_bit,1);

        if (is_64_bit != std::is_same<uint_t,uint64_t>::value) {
            std::cout << "error: cannot load a" << (is_64_bit ? "64" : "32") << "-bit"
            << " index into a " << (is_64_bit ? "32" : "64") << "-bit index-object" << std::flush;
            return;
        }

        std::streampos pos_data_structure_offsets = in.tellg();
        std::streamoff offs_end;
        in.read((char*)&offs_end,sizeof(std::streamoff));

        uint8_t num_supports;
        in.read((char*)&num_supports,1);
        this->support.resize(num_supports);
        in.read((char*)&this->support[0],num_supports*sizeof(move_r_support));

        if (!is_subset_of(support,this->support)) {
            std::cout << "error: cannot load an index with support it has not been built with" << std::flush;
            return;
        }

        this->support = support;

        in.read((char*)&n,sizeof(uint_t));
        in.read((char*)&sigma,1);
        in.read((char*)&r,sizeof(uint_t));
        in.read((char*)&r_,sizeof(uint_t));
        in.read((char*)&a,sizeof(uint16_t));
        in.read((char*)&p_r,sizeof(uint16_t));

        if (p_r > 0) {
            D_e.resize(p_r-1);
            in.read((char*)&D_e[0],(p_r-1)*2*sizeof(uint_t));
        }

        in.read((char*)&chars_remapped,1);
        if (chars_remapped) {
            map_char.resize(256);
            in.read((char*)&map_char[0],256);

            unmap_char.resize(256);
            in.read((char*)&unmap_char[0],256);
        }

        M_LF.load(in);

        if (contains(support,move_r_support::count)) {
            RS_L_.load(in);
        }

        if (contains(support,move_r_support::locate)) {
            in.read((char*)&r__,sizeof(uint_t));
            M_Phi.load(in);

            in.read((char*)&omega_idx,1);
            SA_phi.load(in);
        }

        in.seekg(pos_data_structure_offsets+offs_end,std::ios::beg);
    }

    std::ostream& operator>>(std::ostream& os) {
        serialize(os);
        return os;
    }

    std::istream& operator<<(std::istream& is) {
        load(is);
        return is;
    }
};

#include "algorithms/construction/construction.hpp"
#include "algorithms/misc.cpp"
#include "algorithms/queries.cpp"