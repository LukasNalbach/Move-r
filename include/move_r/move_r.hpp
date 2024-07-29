#pragma once

#include <iostream>
#include <type_traits>
#include <omp.h>
#include <move_r/misc/utils.hpp>
#include <move_r/data_structures/rank_select_support.hpp>
#include <move_r/data_structures/interleaved_vectors.hpp>
#include <move_r/data_structures/move_data_structure/move_data_structure.hpp>
#include <move_r/data_structures/move_data_structure/move_data_structure_l_.hpp>
#include <tsl/sparse_map.h>

/**
 * @brief an operation that can be supported by a move_r object
 */
enum move_r_supp {
    /* support for retrieving (a range of) the input from the index (reverting
    the index); this also includes support for accessing or retrieving (a range in)
    the bwt (reverting in parallel requires the index to be built with locate support) */
    _revert = 0,
    _count = 1, // support for counting the occurrences of a pattern in the input
    _locate = 2 // support for calculating the positions of occurrences of a pattern in the inputstring
};

/**
 * @brief a vector containing all operations supported by move_r
 */
static std::vector<move_r_supp> _full_support = {_revert, _count, _locate};

/**
 * @brief type of locate support
 */
enum move_r_locate_supp {
    _mds = 0, // locate support is implemented using a move data structure to answer Phi^{-1}-queries
    _rlzdsa = 1 // locate support is implemented by relative lepel-ziv encoding the differential suffix array
};

/**
 * @brief move-r construction mode
 */
enum move_r_constr_mode {
    _bigbwt = 0, // builds the bwt with Big-BWT and stores many data structures on disk to reduce peak memory usage
    _suffix_array = 1, // builds the suffix array in-memory and stores no data structures on disk
    _suffix_array_space = 2 // builds the suffix array and stores some data structures on disk
};

/**
 * @brief move-r construction parameters
 */
struct move_r_params {
    std::vector<move_r_supp> support = _full_support; // a vector containing move_r operations to build support for
    move_r_constr_mode mode = _suffix_array; // cosntruction mode to use (default: sa)
    uint16_t num_threads = omp_get_max_threads(); // maximum number of threads to use during the construction
    uint16_t a = 8; // balancing parameter, 2 <= a
    bool log = false; // controls, whether to print log messages
    std::ostream* mf_idx = NULL; // measurement file for the index construciton
    std::ostream* mf_mds = NULL; // measurement file for the move data structure construction
    std::string name_text_file = ""; // name of the input file (used only for measurement output)
};

/**
 * @brief move-r index, size O(r*(a/(a-1)))
 * @tparam locate_support type of locate support (_mds or _rlzdsa)
 * @tparam sym_t value type (default: char for strings)
 * @tparam pos_t index integer type (use uint32_t if input size < UINT_MAX, else uint64_t)
 */
template <move_r_locate_supp locate_support = _mds, typename sym_t = char, typename pos_t = uint32_t>
class move_r {
    protected:

    // check if the position type is supported
    static_assert(
        std::is_same_v<pos_t,uint32_t> ||
        std::is_same_v<pos_t,uint64_t>
    );

    // check if the type of input is supported
    static_assert(
        std::is_same_v<sym_t,char> ||
        std::is_same_v<sym_t,uint8_t> ||
        std::is_same_v<sym_t,uint16_t> ||
        std::is_same_v<sym_t,uint32_t> ||
        std::is_same_v<sym_t,uint64_t> ||
        std::is_same_v<sym_t,int8_t> ||
        std::is_same_v<sym_t,int16_t> ||
        std::is_same_v<sym_t,int32_t> ||
        std::is_same_v<sym_t,int64_t>
    );
    
    // internal (unsigned) symbol type
    using i_sym_t = constexpr_switch_t<
        constexpr_case<sizeof(sym_t) == 1,    uint8_t>,
        constexpr_case<sizeof(sym_t) == 2,    uint16_t>,
        constexpr_case<sizeof(sym_t) == 4,    uint32_t>,
     /* constexpr_case<sizeof(sym_t) == 8, */ uint64_t
    >;
    
    static constexpr bool str_input = std::is_same_v<sym_t,char>; // true <=> the input is a string
    static constexpr bool int_input = !str_input; // true <=> the input is an iteger vector
    static constexpr bool byte_alphabet = sizeof(sym_t) == 1; // true <=> the input uses a byte alphabet

    using map_int_t = std::conditional_t<byte_alphabet,std::vector<uint8_t>,tsl::sparse_map<sym_t,i_sym_t>>; // type of map_int
    using map_ext_t = std::vector<sym_t>; // type of map_ext
    using inp_t = std::conditional_t<str_input,std::string,std::vector<sym_t>>; // input container type
    using rsl_t = rank_select_support<i_sym_t,pos_t,true,true>; // type of RS_L'

    // sample rate of the copy phrases in the rlzdsa
    static constexpr pos_t sr_scp = 4;

    // ############################# INDEX VARIABLES #############################

    pos_t n = 0; // the length of the input
    uint64_t sigma = 0; // the number of distinct characters (including the terminator symbol 1) of T
    pos_t r = 0; // r, the number of runs in L
    pos_t r_ = 0; // r', the number of input/output intervals in M_LF
    pos_t r__ = 0; // r'', the number of input/output intervals in M_Phi^{-1}
    pos_t z = 0; // z, the number of phrases in the rlzdsa
    pos_t z_l = 0; // z_l, the number of literal phrases in the rlzdsa
    pos_t z_c = 0; // z_c, the number of copy-phrases in the rlzdsa
    uint16_t a = 0; // balancing parameter, restricts size to O(r*(a/(a-1))+z), 2 <= a
    uint16_t p_r = 1; // maximum possible number of threads to use while reverting the index
    uint8_t omega_idx = 0; // word width of SA_Phi^{-1}

    std::vector<move_r_supp> _support; // contains all supported operations
    /* true <=> the characters of the input have been remapped internally, because sym_t != char or
       the input invalid characters */
    bool symbols_remapped = false;
    uint64_t size_map_int = 0; // size of _map_int (for byte_alphabet = false)

    // ############################# INDEX DATA STRUCTURES #############################

    // mapping function from the alphabet of the input to the internal effective alphabet
    map_int_t _map_int;
    // mapping function from the internal effective alphabet to the alphabet of the input
    map_ext_t _map_ext;
    
    /* The Move Data Structure for LF. It also stores L', which can be accessed at
    position i with M_LF.L_(i). */
    move_data_structure_l_<pos_t,i_sym_t> _M_LF;
    // rank-select data structure for L'
    rsl_t _RS_L_;

    // The Move Data Structure for Phi^{-1}.
    move_data_structure<pos_t> _M_Phi_m1;
    // [0..r'-1] stores at position x the index of the output interval of M_Phi^{-1} that starts with SA_s[x] = SA[M_LF.p[x]]
    interleaved_vectors<pos_t,pos_t> _SA_Phi_m1;

    /* [0..p_r-1], where D_e[i] = <x,j>, where x in [0,r'-1] and j is minimal, s.t. SA_s[x]=j > i* lfloor (n-1)/p rfloor;
    see the parallel revert algorithm to understand why this is useful. */
    std::vector<std::pair<pos_t,pos_t>> _D_e;

    // stores the suffix array values at the starting positions of the input intervals of M_LF, i.e, SA_s[i] = SA[M_LF.p[i]]
    interleaved_vectors<pos_t,pos_t> _SA_s;
    // reference for SA^d (differential suffix array)
    interleaved_vectors<uint64_t,pos_t> _R;
    // bit vector storing the phrase types of the rlzdsa, i.e, PT[i] = 1 <=> phrase i is literal
    plain_bit_vector<pos_t,true,true,true> _PT;
    // compressed bit vector marking the sampled starting positions in SA^d of the copy phrases of the rlzdsa
    sd_array<pos_t> _SCP_S;
    // lengths of the copy phrases of the rlzdsa
    std::vector<uint16_t> _CPL;
    // starting positions in R of the copy phrases of the rlzdsa
    interleaved_vectors<pos_t,pos_t> _SR;
    // literal phrases of the rlzdsa
    interleaved_vectors<pos_t,pos_t> _LP;
    
    // ############################# INTERNAL METHODS #############################

    /**
     * @brief sets SA_Phi^{-1}[x] to idx
     * @param x [0..r-1]
     * @param idx [0..r''-1]
     */
    inline void set_SA_Phi_m1(pos_t x, pos_t idx) {
        _SA_Phi_m1.template set<0,pos_t>(x,idx);
    }

    /**
     * @brief adds implicitly supported move_r operations to the vector support and sorts it afterwards
     * @param support a vector containing move_r operations
     */
    static void adjust_supports(std::vector<move_r_supp>& support) {
        if (contains(support,_locate) && !(contains(support,_count))) {
            support.emplace_back(_count);
        }
        
        if (contains(support,_count) && !(contains(support,_revert))) {
            support.emplace_back(_revert);
        }
        
        ips4o::sort(support.begin(),support.end());
        support.erase(std::unique(support.begin(),support.end()),support.end());
    }

    class construction;

    // ############################# CONSTRUCTORS #############################

    public:
    move_r() = default;

    /**
     * @brief constructs a move_r index of the input
     * @param input the input
     * @param params construction parameters
     */
    move_r(inp_t& input, move_r_params params = {}) {
        construction(*this,input,false,params);
    }

    /**
     * @brief constructs a move_r index of the input
     * @param input the input
     * @param params construction parameters
     */
    move_r(inp_t&& input, move_r_params params = {}) {
        construction(*this,input,true,params);
    }

    /**
     * @brief constructs a move_r index from an input file
     * @param input_file input file
     * @param params construction parameters
     */
    move_r(std::ifstream& input_file, move_r_params params = {}) requires(str_input) {
        construction(*this,input_file,params);
    }

    /**
     * @brief constructs a move_r index from a suffix array and a bwt
     * @tparam sa_sint_t suffix array signed integer type
     * @param suffix_array vector containing the suffix array of the input
     * @param bwt string containing the bwt of the input, where $ = 1
     * @param params construction parameters
     */
    template <typename sa_sint_t>
    move_r(std::vector<sa_sint_t>& suffix_array, std::string& bwt, move_r_params params = {}) requires(str_input) {
        construction(*this,suffix_array,bwt,params);
    }

    // ############################# MISC PUBLIC METHODS #############################

    /**
     * @brief returns the size of the input
     * @return size of the input
     */
    inline pos_t input_size() const {
        return n-1;
    }

    /**
     * @brief returns the number of distinct characters in the input (alphabet size)
     * @return alphabet_size 
     */
    inline uint32_t alphabet_size() const {
        return sigma-1;
    }

    /**
     * @brief returns the number of runs in the bwt
     * @return number of runs in the bwt 
     */
    inline pos_t num_bwt_runs() const {
        return r;
    }

    /**
     * @brief returns the number of phrases in the rlzdsa
     * @return number of phrases in the rlzdsa
     */
    inline pos_t num_phrases_rlzdsa() const {
        return z;
    }

    /**
     * @brief returns the number of literal phrases in the rlzdsa
     * @return number of literal phrases in the rlzdsa
     */
    inline pos_t num_literal_phrases_rlzdsa() const {
        return z_l;
    }

    /**
     * @brief returns the number of copy phrases in the rlzdsa
     * @return number of copy phrases in the rlzdsa
     */
    inline pos_t num_copy_phrases_rlzdsa() const {
        return z_c;
    }

    /**
     * @brief returns the balancing parameter the index has been built with
     * @return balancing parameter 
     */
    inline uint16_t balancing_parameter() const {
        return a;
    }

    /**
     * @brief returns the number omega_idx of bits used by one entry in SA_Phi^{-1} (word width of SA_Phi^{-1})
     * @return omega_idx
     */
    inline uint8_t width_saphi() const {
        return omega_idx;
    }

    /**
     * @brief returns the maximum number of threads that can be used to revert the index
     * @return maximum number of threads that can be used to revert the index 
     */
    inline uint16_t max_revert_threads() const {
        return p_r;
    }

    /**
     * @brief returns a vector containing the supported operations of this index
     * @return vector containing the operations
     */
    inline std::vector<move_r_supp> supported_operations() const {
        return _support;
    }

    /**
     * @brief returns whether the provided operation is supported by this index
     * @param operation a move_r operation
     * @return whether operation is supported by this index
     */
    inline bool does_support(move_r_supp operation) const {
        return contains(_support,operation);
    }

    /**
     * @brief returns the size of the data structure in bytes
     * @return size of the data structure in bytes
     */
    uint64_t size_in_bytes() const {
        uint64_t size =
            _support.size()*sizeof(move_r_supp)+ // variables
            4*sizeof(pos_t)+3+2*sizeof(uint16_t)+ // ...
            p_r*sizeof(pos_t)+ // D_e
            _M_LF.size_in_bytes()+ // M_LF and L'
            size_map_int+ // map_int
            sizeof(sym_t)*sigma+ // map_ext
            _RS_L_.size_in_bytes(); // RS_L'

        if (contains(_support,_locate)) {
            if constexpr (locate_support == _mds) {
                size +=
                    _M_Phi_m1.size_in_bytes()+ // M_Phi^{-1}
                    _SA_Phi_m1.size_in_bytes(); // SA_Phi^{-1}
            } else {
                size +=
                    _SA_s.size_in_bytes()+ // SA_s
                    _R.size_in_bytes()+ // R
                    (z_c+2)*sizeof(uint16_t)+ // CPL
                    _SCP_S.size_in_bytes()+ // SCP_S
                    _SR.size_in_bytes()+ // SR
                    _LP.size_in_bytes()+ // LP
                    _PT.size_in_bytes(); // PT
            }
        }

        return size;
    }

    /**
     * @brief logs the index data structure sizes to cout
     */
    void log_data_structure_sizes() const {
        std::cout << "index size: " << format_size(size_in_bytes()) << std::endl;

        uint64_t size_l_ = (_M_LF.width_l_()/8)*(r_+1);
        std::cout << "M_LF: " << format_size(_M_LF.size_in_bytes()-size_l_) << std::endl;
        std::cout << "L': " << format_size(size_l_) << std::endl;

        if constexpr (int_input) {
            std::cout << "map_int: " << format_size(size_map_int) << std::endl;
            std::cout << "map_ext: " << format_size(sizeof(sym_t)*sigma) << std::endl;
        }

        if (does_support(_count)) {
            std::cout << "RS_L': " << format_size(_RS_L_.size_in_bytes()) << std::endl;
        }

        if (does_support(_locate)) {
            if constexpr (locate_support == _mds) {
                std::cout << "M_Phi^{-1}: " << format_size(_M_Phi_m1.size_in_bytes()) << std::endl;
                std::cout << "SA_Phi^{-1}: " << format_size(_SA_Phi_m1.size_in_bytes()) << std::endl;
            } else {
                std::cout << "SA_s: " << format_size(_SA_s.size_in_bytes()) << std::endl;
                std::cout << "R: " << format_size(_R.size_in_bytes()) << std::endl;
                std::cout << "CPL: " << format_size((z_c+2)*sizeof(uint16_t)) << std::endl;
                std::cout << "SCP_S: " << format_size(_SCP_S.size_in_bytes()) << std::endl;
                std::cout << "SR: " << format_size(_SR.size_in_bytes()) << std::endl;
                std::cout << "LP: " << format_size(_LP.size_in_bytes()) << std::endl;
                std::cout << "PT: " << format_size(_PT.size_in_bytes()) << std::endl;
            }
        }
    }

    /**
     * @brief logs the index data structure sizes to the output stream out
     * @param out an output stream
     */
    void log_data_structure_sizes(std::ostream& out) const {
        out << " size_index=" << size_in_bytes();
        uint64_t size_l_ = (_M_LF.width_l_()/8)*(r_+1);
        out << " size_m_lf=" << _M_LF.size_in_bytes()-size_l_;
        out << " size_l_=" << size_l_;

        if constexpr (int_input) {
            out << " size_map_int=" << size_map_int;
            out << " size_map_ext=" << sizeof(sym_t)*sigma;
        }
        
        if (does_support(_count)) {
            out << " size_rs_l_=" << _RS_L_.size_in_bytes();
        }

        if (does_support(_locate)) {
            if constexpr (locate_support == _mds) {
                out << " size_m_phim1=" << _M_Phi_m1.size_in_bytes();
                out << " size_sa_phim1=" << _SA_Phi_m1.size_in_bytes();
            } else {
                out << "size_sa_s: " << _SA_s.size_in_bytes();
                out << "size_r: " << _R.size_in_bytes();
                out << "size_cpl: " << (z_c+2)*sizeof(uint16_t);
                out << "size_scp: " << _SCP_S.size_in_bytes();
                out << "size_sr: " << _SR.size_in_bytes();
                out << "size_lp: " << _LP.size_in_bytes();
                out << "size_pt: " << _PT.size_in_bytes();
            }
        }
    }

    // ############################# PUBLIC ACCESS METHODS #############################

    /**
     * @brief returns a reference to M_LF
     * @return M_LF
     */
    inline const move_data_structure_l_<pos_t,i_sym_t>& M_LF() const {
        return _M_LF;
    }

    /**
     * @brief returns a reference to M_Phi^{-1}
     * @return M_Phi^{-1}
     */
    inline const move_data_structure<pos_t>& M_Phi_m1() const {
        return _M_Phi_m1;
    }

    /**
     * @brief returns a reference to RS_L'
     * @return RS_L'
     */
    inline const rsl_t& RS_L_() const {
        return _RS_L_;
    }

    /**
     * @brief returns a reference to R
     * @return R
     */
    inline const interleaved_vectors<uint64_t,pos_t>& R() const {
        return _R;
    }

    /**
     * @brief returns a reference to PT
     * @return PT
     */
    inline const plain_bit_vector<pos_t,true,true,true>& PT() const {
        return _PT;
    }

    /**
     * @brief returns a reference to CPL
     * @return CPL
     */
    inline const std::vector<uint16_t>& CPL() const {
        return _CPL;
    }

    /**
     * @brief returns a reference to SCP_S
     * @return SCP_S
     */
    inline const sd_array<pos_t>& SCP_S() const {
        return _SCP_S;
    }

    /**
     * @brief returns a reference to SR
     * @return SR
     */
    inline const interleaved_vectors<pos_t,pos_t>& SR() const {
        return _SR;
    }

    /**
     * @brief returns a reference to LP
     * @return LP
     */
    inline const interleaved_vectors<pos_t,pos_t>& LP() const {
        return _LP;
    }

    /**
     * @brief returns R[x]
     * @param x [0..|R|-1] index in R
     * @return R[x]
     */
    inline uint64_t R(pos_t x) const {
        return _R[x];
    }

    /**
     * @brief returns PT[x]
     * @param x [0..z-1] index in PT
     * @return PT[x]
     */
    inline bool PT(pos_t x) const {
        return _PT[x];
    }

    /**
     * @brief returns CPL[x]
     * @param x [0..z_c-1] index in CPL
     * @return CPL[x]
     */
    inline uint16_t CPL(pos_t x) const {
        return _CPL[x];
    }

    /**
     * @brief returns SCP_S[x]
     * @param x [0..z_c/sr_scp-1] index in SCP_S
     * @return SCP_S[x]
     */
    inline pos_t SCP_S(pos_t x) const {
        return _SCP_S.select_1(x+1);
    }

    /**
     * @brief returns SR[x]
     * @param x [0..z_r-1] index in SR
     * @return SR[x]
     */
    inline pos_t SR(pos_t x) const {
        return _SR[x];
    }

    /**
     * @brief returns LP[x]
     * @param x [0..z_l-1] index in LP
     * @return LP[x]
     */
    inline pos_t LP(pos_t x) const {
        return _LP[x];
    }

    /**
     * @brief returns SA_Phi^{-1}[x]
     * @param x [0..r''-1]
     * @return SA_Phi^{-1}[x]
     */
    inline pos_t SA_Phi_m1(pos_t x) const {
        return _SA_Phi_m1[x];
    }

    /**
     * @brief returns SA_s[x]
     * @param x [0..r'-1] the starting position of the x-th input
     * interval in M_LF must be a starting position of a bwt run
     * @return SA_s[x]
     */
    inline pos_t SA_s(pos_t x) const {
        if constexpr (locate_support == _mds) {
            return _M_Phi_m1.q(_SA_Phi_m1[x]);
        } else {
            return _SA_s[x];
        }
    }

    /**
     * @brief returns L'[x]
     * @param x [0..r'-1]
     * @return L'[x]
     */
    inline i_sym_t L_(pos_t x) const {
        return _M_LF.L_(x);
    }

    /**
     * @brief reinterprets an int8_t as a uint8_t
     * @param sym symbol
     * @return sym reinterpreted as uint8_t
     */
    uint8_t symbol_idx(uint8_t sym) {
        return sym;
    }

    /**
     * @brief reinterprets an int8_t as a uint8_t
     * @param sym symbol
     * @return sym reinterpreted as uint8_t
     */
    uint8_t symbol_idx(int8_t sym) {
        return *reinterpret_cast<uint8_t*>(&sym);
    }

    /**
     * @brief reinterprets an char as a uint8_t
     * @param sym symbol
     * @return sym reinterpreted as uint8_t
     */
    uint8_t symbol_idx(char sym) const {
        return *reinterpret_cast<uint8_t*>(&sym);
    }

    /**
     * @brief maps a symbol to its corresponding symbol in the internal effective alphabet
     * @param sym symbol
     * @return its corresponding symbol in the internal effective alphabet
     */
    inline i_sym_t map_symbol(sym_t sym) const {
        if constexpr (byte_alphabet) {
            return symbols_remapped ? _map_int[symbol_idx(sym)] : symbol_idx(sym);
        } else {
            auto res = _map_int.find(sym);

            if (res == _map_int.end()) {
                return 0;
            } else {
                return (*res).second;
            }
        }
    }

    /**
     * @brief maps a symbol that occurs in the internal effective alphabet to its corresponding
     *        symbol in the input
     * @param sym a symbol that occurs in the internal effective alphabet
     * @return its corresponding symbol in the input
     */
    inline sym_t unmap_symbol(i_sym_t sym) const {
        if constexpr (byte_alphabet) {
            return symbols_remapped ? _map_ext[sym] : sym;
        } else {
            return _map_ext[sym];
        }
    }

    /**
     * @brief returns D_e[i]
     * @param i [0..max_revert_threads()-2]
     * @return D_e[i]
     */
    inline std::pair<pos_t,pos_t> D_e(uint16_t i) const {
        return _D_e[i];
    }

    // ############################# QUERY METHODS #############################

    /**
     * @brief returns L[i], where $ = 0, so if the input contained 0, the output is not equal to the real bwt
     * @param x [0..input size]
     * @return L[i]
     */
    inline sym_t BWT(pos_t i) const;

    /**
     * @brief returns SA[i]
     * @param x [0..input size]
     * @return SA[i]
     */
    pos_t SA(pos_t i) const;

    /**
     * @brief stores the variables needed to perform count- and locate-queries
     */
    struct query_context {
        protected:

        pos_t l;  // length of the currently matched pattern
        pos_t b,e,b_,e_,hat_b_ap_y; // variables for backward search
        int64_t y; // variable for backward search
        pos_t i; // current position in the suffix array interval
        pos_t s; // current suffix s = SA[i] in the suffix array interval
        pos_t s_; // index of the input inteval of M_Phi^{-1} containing s
        pos_t x_p,x_lp,x_cp,x_r,s_np; // variables for decoding the rlzdsa

        const move_r<locate_support,sym_t,pos_t>* idx; // index to query

        public:
        /**
         * @brief constructs a new query context for the index idx
         * @param idx an index
         */
        query_context(const move_r<locate_support,sym_t,pos_t>& idx) {
            this->idx = &idx;
            reset();
        }

        /**
         * @brief resets the query context to an empty pattern
         */
        inline void reset() {
            idx->init_backward_search(b,e,b_,e_,hat_b_ap_y,y);
            l = 0;
            i = b;
        }

        /**
         * @brief returns the length of the currently matched pattern
         * @return length of the currently matched pattern
         */
        inline pos_t length() const {
            return l;
        }

        /**
         * @brief returns the overall number of occurrences of the currently matched pattern
         * @return overall number of occurrences
         */
        inline pos_t num_occ() const {
            return e >= b ? e-b+1 : 0;
        }

        /**
         * @brief returns the number of remaining (not yet reported) occurrences of the currently matched pattern
         * @return number of remaining occurrences
         */
        inline pos_t num_occ_rem() const {
            return e >= i ? e-i+1 : 0;
        }

        /**
         * @brief returns the suffix array interval of the currently matched pattern
         * @return suffix array interval
         */
        inline std::pair<pos_t,pos_t> sa_interval() const {
            return std::make_pair(b,e);
        }

        /**
         * @brief prepends sym to the currently matched pattern P; if symP occurs in the input, true is
         * returned and the query context is adjusted to store the information for the pattern symP; else,
         * false is returned and the query context is not modified
         * @param sym 
         * @return whether symP occurs in the input
         */
        bool prepend(sym_t sym);

        /**
         * @brief reports the next occurrence of the currently matched pattern
         * @return next occurrence
         */
        inline pos_t next_occ();

        /**
         * @brief locates the remaining (not yet reported) occurrences of the currently matched pattern
         * @param Occ vector to append the occurrences to
         */
        inline void locate(std::vector<pos_t>& Occ);

        /**
         * @brief locates the remaining (not yet reported) occurrences of the currently matched pattern
         * @return vector containing the occurrences
         */
        std::vector<pos_t> locate() {
            std::vector<pos_t> Occ;
            locate(Occ);
            return Occ;
        }
    };

    /**
     * @brief returns a query context for the index
     * @return query_context 
     */
    inline query_context query() const {
        return query_context(*this);
    }

    protected:
    /**
     * @brief initializes the variables to start a new backward search
     * @param b Left interval limit of the suffix array interval.
     * @param e Right interval limit of the suffix array interval.
     * @param b_ index of the input interval in M_LF containing b.
     * @param e_ index of the input interval in M_LF containing e.
     * @param hat_b_ap_y \hat{b}'_y
     * @param y y
     */
    inline void init_backward_search(
        pos_t& b, pos_t& e,
        pos_t& b_, pos_t& e_,
        pos_t& hat_b_ap_y, int64_t& y
    ) const {
        b = 0;
        e = n-1;
        b_ = 0;
        e_ = r_-1;
        hat_b_ap_y = 0;
        y = -1;
    }
    
    /**
     * @brief prepends sym to the currently matched pattern P, adjusts the variables to store
     * the query context for the pattern symP and returns whether symP occurs in the input
     * @param sym next symbol to match
     * @param b Left interval limit of the suffix array interval.
     * @param e Right interval limit of the suffix array interval.
     * @param b_ index of the input interval in M_LF containing b.
     * @param e_ index of the input interval in M_LF containing e.
     * @param hat_b_ap_y \hat{b}'_y
     * @param y y
     * @return whether symP occurs in the input
     */
    bool backward_search_step(
        sym_t sym,
        pos_t& b, pos_t& e,
        pos_t& b_, pos_t& e_,
        pos_t& hat_b_ap_y, int64_t& y
    ) const;

    /**
     * @brief Sets the up a Phi-move-pair for the suffix array sample at the end position of the x-th input interval in M_LF
     * @param x an input interval in M_LF (the end position of the x-th input interval in M_LF must be an end position of a BWT run)
     * @param s variable to store the suffix array sample at position l'_{x+1}-1
     * @param s_ variable to store the index of the input interval in M_Phi^{-1} containing s
     */
    inline void setup_phi_move_pair(pos_t& x, pos_t& s, pos_t& s_) const;

    /**
     * @brief prepares the variables to decode SA[b]
     * @param b left interval limit of the suffix array interval
     * @param e right interval limit of the suffix array interval
     * @param s variable to store SA[b] in
     * @param s_ index of the input interval in M_Phi^{-1} containing s
     * @param hat_b_ap_y \hat{b}'_y
     * @param y y
     */
    inline void init_phi(
        pos_t& b, pos_t& e,
        pos_t& s, pos_t& s_,
        pos_t& hat_b_ap_y, int64_t& y
    ) const;
    
    /**
     * @brief prepares the variables to decode SA[i]
     * @param i current position in the suffix array
     * @param x_p phrase-index of the phrase of the rlzdsa contianing i
     * @param x_lp literal-phrase index of the current or next literal phrase of the rlzdsa
     * @param x_cp copy-phrase index of the current or next copy-phrase of the rlzdsa
     * @param x_r position in R inside the current copy-phrase (or the starting position in R of the next copy phrase) of the rlzdsa
     * @param s_np starting position in the rlzdsa of the next phrase of the rlzdsa
     */
    inline void init_rlzdsa(
        pos_t& i,
        pos_t& x_p, pos_t& x_lp, pos_t& x_cp, pos_t& x_r, pos_t& s_np
    ) const;
    
    /**
     * @brief prepares the context to decode SA[i]; if there
     * is a literal phrase at position i, s stores SA[i], else SA[i-1]
     * @param i [0..n-1] position in the suffix array
     * @param s variable to (possibly) store SA[i] or SA[i-1] in
     * @param x_p phrase-index of the phrase of the rlzdsa contianing i
     * @param x_lp literal-phrase index of the current or next literal phrase of the rlzdsa
     * @param x_cp copy-phrase index of the current or next copy-phrase of the rlzdsa
     * @param x_r position in R inside the current copy-phrase (or the starting position in R of the next copy phrase) of the rlzdsa
     * @param s_np starting position in the rlzdsa of the next phrase of the rlzdsa
     */
    inline void init_rlzdsa(
        pos_t& i, pos_t& s,
        pos_t& x_p, pos_t& x_lp, pos_t& x_cp, pos_t& x_r, pos_t& s_np
    ) const;

    /**
     * @brief decodes and stores SA[i] in s and prepares the context to decode
     * SA[i+1]; the context must be prepared to decode SA[i]
     * @param i [0..n-1] position in the suffix array
     * @param s variable to store SA[i] in
     * @param x_p phrase-index of the phrase of the rlzdsa contianing i
     * @param x_lp literal-phrase index of the current or next literal phrase of the rlzdsa
     * @param x_cp copy-phrase index of the current or next copy-phrase of the rlzdsa
     * @param x_r position in R inside the current copy-phrase (or the starting position in R of the next copy phrase) of the rlzdsa
     * @param s_np starting position in the rlzdsa of the next phrase of the rlzdsa
     */
    inline void next_rlzdsa(
        pos_t& i, pos_t& s,
        pos_t& x_p, pos_t& x_lp, pos_t& x_cp, pos_t& x_r, pos_t& s_np
    ) const;

    /**
     * @brief locates the remaining (not yet reported) occurrences of the currently matched pattern
     * @param i current position in the suffix array
     * @param e right interval limit of the suffix array interval
     * @param s current suffix array value
     * @param x_p phrase-index of the phrase of the rlzdsa contianing i
     * @param x_lp literal-phrase index of the current or next literal phrase of the rlzdsa
     * @param x_cp copy-phrase index of the current or next copy-phrase of the rlzdsa
     * @param x_r position in R inside the current copy-phrase (or the starting position in R of the next copy phrase) of the rlzdsa
     * @param s_np starting position in the rlzdsa of the next phrase of the rlzdsa
     */
    inline void locate_rlzdsa(
        pos_t& i, pos_t& e, pos_t& s,
        pos_t& x_p, pos_t& x_lp, pos_t& x_cp, pos_t& x_r, pos_t& s_np,
        std::vector<pos_t>& Occ
    ) const;

    public:
    /**
     * @brief returns the number of occurrences of P in the input
     * @param P the pattern to count in the input
     * @return the number of occurrences of P in the input
     */
    inline pos_t count(const inp_t& P) const;

    /**
     * @brief locates the pattern P in the input
     * @param P the pattern to locate in the input
     * @return a vector containing the occurrences of P in the input
     */
    inline std::vector<pos_t> locate(const inp_t& P) const {
        std::vector<pos_t> Occ;
        locate(P,Occ);
        return Occ;
    }

    /**
     * @brief locates the pattern P in the input and appends the positions of the occurrences to Occ
     * @param P the pattern to locate in the input
     * @param Occ vector to append the occurrences of P in the input to
     */
    void locate(const inp_t& P, std::vector<pos_t>& Occ) const;

    // ############################# RETRIEVE-RANGE METHODS #############################

    struct retrieve_params {
        pos_t l = 1; // left range limit
        pos_t r = 0; // right range limit
        uint16_t num_threads = omp_get_max_threads(); // maximum number of threads to use
        // maximum number of bytes to allocate (only applicable if the method writes data to a file; default (if set to -1): ~ (r-l+1)/500)
        int64_t max_bytes_alloc = -1;
    };

    protected:
    /**
     * @brief adjusts retrieve parameters, ensures (0 <= l <= r <= range_max); if l > r, it sets [l,r] <- [0,range_max]
     * @param params retrieve parameters to adjust
     * @param range_max maximum value for r
     */
    inline static void adjust_retrieve_params(retrieve_params& params, pos_t range_max) {
        if (params.l > params.r) {
            params.l = 0;
            params.r = range_max;
        }

        params.r = std::max(params.r,range_max);
    }

    /**
     * @brief executes retrieve_method with the parameters l, r and num_threads, buffers the output in num_threads
     * buffers and num_threads temporary files and then writes the temporary files into the file out
     * @tparam output_t type of the output data
     * @tparam output_reversed controls, whether the output should be reversed
     * @param retrieve_method function, whiches output should be buffered
     * @param file_name name of the file to write the output to
     * @param params parameters
     */
    template <typename output_t, bool output_reversed>
    void retrieve_range(
        void(move_r<locate_support,sym_t,pos_t>::*retrieve_method)(const std::function<void(pos_t,output_t)>&,retrieve_params)const,
        std::string file_name, retrieve_params params
    ) const;

    public:
    /**
     * @brief returns the bwt in the range [l,r] (0 <= l <= r <= input size), else
     * if l > r, then the whole bwt is returned (default); $ = 0, so if the input contained 0, the output is not
     * equal to the real bwt
     * @param params parameters
     * @return the bwt range [l,r]
     */
    inp_t BWT(retrieve_params params = {}) const {
        adjust_retrieve_params(params,n-1);
        inp_t L;
        no_init_resize(L,params.r-params.l+1);
        BWT([&L,&params](pos_t i, sym_t c){L[i-params.l] = c;},params);
        return L;
    }

    /**
     * @brief reports the characters in the bwt in the range [l,r] (0 <= l <= r <= input size), else if l > r, then all
     * characters of the bwt are reported (default); $ = 0, so if the input contained 0, the output is not equal to the real bwt
     * @param report function that is called with every tuple (i,c) as a parameter, where i in [l,r] and c = L[i]; if num_threads = 1,
     * then the values are reported from left to right, if num_threads > 1, the order may vary
     * @param params parameters
     */
    void BWT(const std::function<void(pos_t,sym_t)>& report, retrieve_params params = {}) const;

    /**
     * @brief writes the characters in the bwt in the range [l,r] blockwise to the file out (0 <= l <= r <= input size), else if
     * l > r, then the whole bwt is written (default); $ = 0, so if the input contained 0, the output is not equal to the real bwt
     * @param file_name name of the file to write the bwt to
     * @param params parameters
     */
    void BWT(std::string file_name, retrieve_params params = {}) const {
        adjust_retrieve_params(params,n-1);
        retrieve_range<sym_t,false>(&move_r<locate_support,sym_t,pos_t>::BWT,file_name,params);
    }

    /**
     * @brief returns the input in the range [l,r] (0 <= l <= r < input size), else
     * if l > r, then the whole input is returned (default)
     * @param params parameters
     * @return the input range [l,r]
     */
    inp_t revert(retrieve_params params = {}) const {
        adjust_retrieve_params(params,n-2);
        inp_t input;
        no_init_resize(input,params.r-params.l+1);
        revert([&input,&params](pos_t i, sym_t c){input[i-params.l] = c;},params);
        return input;
    }

    /**
     * @brief reports the characters in the input in the range [l,r] (0 <= l <= r < input size), else if l > r, then
     * all characters of the input are reported (default); if num_threads = 1, then the values are reported from right
     * to left, if num_threads > 1, the order may vary
     * @param report function that is called with every tuple (i,c) as a parameter, where i in [l,r] and c = input[i]
     * @param params parameters
     */
    void revert(const std::function<void(pos_t,sym_t)>& report, retrieve_params params = {}) const;

    /**
     * @brief reverts the input in the range [l,r] blockwise and writes it to the file out (0 <= l <= r < input size),
     * else if l > r, then the whole input is reverted (default)
     * @param file_name name of the file to write the reverted input to
     * @param params parameters
     */
    void revert(std::string file_name, retrieve_params params = {}) const {
        adjust_retrieve_params(params,n-2);
        retrieve_range<sym_t,true>(&move_r<locate_support,sym_t,pos_t>::revert,file_name,params);
    }
    
    /**
     * @brief rebuilds and returns the suffix array in the range [l,r] (0 <= l <= r <= input size),
     * else if l > r, then the whole suffix array is rebuilt (default)
     * @param params parameters
     * @return the suffix array range [l,r]
     */
    std::vector<pos_t> SA(retrieve_params params = {}) const {
        adjust_retrieve_params(params,n-1);
        std::vector<pos_t> SA_range;
        no_init_resize(SA_range,params.r-params.l+1);
        SA([&SA_range,&params](pos_t i, pos_t s){SA_range[i-params.l] = s;},params);
        return SA_range;
    }

    /**
     * @brief reports the suffix array values in the range [l,r] (0 <= l <= r <= input size), else if l > r, then the
     * whole suffix array is reported (default); if num_threads = 1, then the values are reported from left to right,
     * if num_threads > 1, the order may vary
     * @param report function that is called with every tuple (i,s) as a parameter, where i in [l,r] and s = SA[i]
     * @param params parameters
     */
    void SA(const std::function<void(pos_t,pos_t)>& report, retrieve_params params = {}) const;

    /**
     * @brief writes the values in the suffix array of the input in the range [l,r] blockwise to the file out (0 <= l <= r <= input size),
     * else if l > r, then the whole suffix array is written (default)
     * @param file_name name of the file to write the suffix array to
     * @param params parameters
     */
    void SA(std::string file_name, retrieve_params params = {}) const {
        adjust_retrieve_params(params,n-1);
        retrieve_range<pos_t,false>(&move_r<locate_support,sym_t,pos_t>::SA,file_name,params);
    }

    // ############################# SERIALIZATION METHODS #############################

    /**
     * @brief stores the index to an output stream
     * @param out output stream to store the index to
     * @param support supported operations to store data structures for
     */
    void serialize(std::ostream& out, std::vector<move_r_supp> support = {}) const {
        adjust_supports(support);

        if (!is_subset_of(support,this->_support)) {
            std::cout << "error: cannot store an index with support it has not been built with" << std::flush;
            return;
        }

        if (support.empty()) {
            support = this->_support;
        }

        bool is_64_bit = std::is_same_v<pos_t,uint64_t>;
        out.write((char*)&is_64_bit,1);
        move_r_locate_supp _locate_support = locate_support;
        out.write((char*)&_locate_support,sizeof(move_r_locate_supp));

        std::streampos pos_data_structure_offsets = out.tellp();
        out.seekp(pos_data_structure_offsets+(std::streamoff)sizeof(std::streamoff),std::ios::beg);

        uint8_t num_supports = support.size();
        out.write((char*)&num_supports,1);
        out.write((char*)&support[0],num_supports*sizeof(move_r_supp));

        out.write((char*)&n,sizeof(pos_t));
        out.write((char*)&sigma,sizeof(uint32_t));
        out.write((char*)&r,sizeof(pos_t));
        out.write((char*)&r_,sizeof(pos_t));
        out.write((char*)&a,sizeof(uint16_t));
        out.write((char*)&p_r,sizeof(uint16_t));

        if (p_r > 0) {
            out.write((char*)&_D_e[0],(p_r-1)*2*sizeof(pos_t));
        }

        out.write((char*)&symbols_remapped,1);
        if (symbols_remapped) {
            if constexpr (byte_alphabet) {
                out.write((char*)&_map_int[0],256);
                out.write((char*)&_map_ext[0],256);
            } else {
                write_to_file(out,(char*)&_map_ext[0],sizeof(sym_t)*sigma);
                std::vector<std::pair<sym_t,i_sym_t>> map_int_vec(_map_int.begin(),_map_int.end());
                write_to_file(out,(char*)&map_int_vec[0],sizeof(std::pair<sym_t,i_sym_t>)*sigma);
            }
        }

        _M_LF.serialize(out);

        if (contains(support,_count)) {
            _RS_L_.serialize(out);
        }

        if (contains(support,_locate)) {
            if constexpr (locate_support == _mds) {
                out.write((char*)&r__,sizeof(pos_t));
                _M_Phi_m1.serialize(out);

                out.write((char*)&omega_idx,1);
                _SA_Phi_m1.serialize(out);
            } else {
                out.write((char*)&z,sizeof(pos_t));
                out.write((char*)&z_l,sizeof(pos_t));
                out.write((char*)&z_c,sizeof(pos_t));

                _SA_s.serialize(out);
                _R.serialize(out);
                _SCP_S.serialize(out);
                write_to_file(out,(char*)&_CPL[0],(z_c+2)*sizeof(uint16_t));
                _SR.serialize(out);
                _LP.serialize(out);
                _PT.serialize(out);
            }
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
    void load(std::istream& in, std::vector<move_r_supp> support = {}) {
        adjust_supports(support);

        bool is_64_bit;
        in.read((char*)&is_64_bit,1);

        if (is_64_bit != std::is_same_v<pos_t,uint64_t>) {
            std::cout << "error: cannot load a" << (is_64_bit ? "64" : "32") << "-bit"
            << " index into a " << (is_64_bit ? "32" : "64") << "-bit index-object" << std::flush;
            return;
        }

        move_r_locate_supp _locate_support;
        in.read((char*)&_locate_support,sizeof(move_r_locate_supp));

        if (_locate_support != locate_support) {
            std::cout << "error: cannot load an index with " << (_locate_support ? "phi" : "rlzdsa") << "-locate "
            <<  "support into an index object with " << (locate_support ? "phi" : "rlzdsa") << "-locate support" << std::flush;
            return;
        }

        std::streampos pos_data_structure_offsets = in.tellg();
        std::streamoff offs_end;
        in.read((char*)&offs_end,sizeof(std::streamoff));

        uint8_t num_supports;
        in.read((char*)&num_supports,1);
        this->_support.resize(num_supports);
        in.read((char*)&this->_support[0],num_supports*sizeof(move_r_supp));

        if (!is_subset_of(support,this->_support)) {
            std::cout << "error: cannot load an index with support it has not been built with" << std::flush;
            return;
        }

        if (support.empty()) {
            support = this->_support;
        }

        this->_support = support;

        in.read((char*)&n,sizeof(pos_t));
        in.read((char*)&sigma,sizeof(uint32_t));
        in.read((char*)&r,sizeof(pos_t));
        in.read((char*)&r_,sizeof(pos_t));
        in.read((char*)&a,sizeof(uint16_t));
        in.read((char*)&p_r,sizeof(uint16_t));
        
        if (p_r > 0) {
            _D_e.resize(p_r-1);
            in.read((char*)&_D_e[0],(p_r-1)*2*sizeof(pos_t));
        }

        in.read((char*)&symbols_remapped,1);
        if (symbols_remapped) {
            if constexpr (byte_alphabet) {
                _map_int.resize(256);
                in.read((char*)&_map_int[0],256);

                _map_ext.resize(256);
                in.read((char*)&_map_ext[0],256);
            } else {
                no_init_resize(_map_ext,sigma);
                read_from_file(in,(char*)&_map_ext[0],sizeof(sym_t)*sigma);

                std::vector<std::pair<sym_t,i_sym_t>> map_int_vec;
                no_init_resize(map_int_vec,sigma);
                read_from_file(in,(char*)&map_int_vec[0],sizeof(std::pair<sym_t,i_sym_t>)*sigma);
                uint64_t alloc_before = malloc_count_current();
                _map_int.insert(map_int_vec.begin(),map_int_vec.end());
                size_map_int = malloc_count_current()-alloc_before;
            }
        }

        _M_LF.load(in);

        if (contains(support,_count)) {
            _RS_L_.load(in);
        }

        if (contains(support,_locate)) {
            if constexpr (locate_support == _mds) {
                in.read((char*)&r__,sizeof(pos_t));
                _M_Phi_m1.load(in);

                in.read((char*)&omega_idx,1);
                _SA_Phi_m1.load(in);
            } else {
                in.read((char*)&z,sizeof(pos_t));
                in.read((char*)&z_l,sizeof(pos_t));
                in.read((char*)&z_c,sizeof(pos_t));

                _SA_s.load(in);
                _R.load(in);
                _SCP_S.load(in);
                no_init_resize(_CPL,z_c+2);
                read_from_file(in,(char*)&_CPL[0],(z_c+2)*sizeof(uint16_t));
                _SR.load(in);
                _LP.load(in);
                _PT.load(in);
            }
        }

        in.seekg(pos_data_structure_offsets+offs_end,std::ios::beg);
    }

    std::ostream& operator>>(std::ostream& os) const {
        serialize(os);
        return os;
    }

    std::istream& operator<<(std::istream& is) {
        load(is);
        return is;
    }
};

#include "algorithms/construction/construction.hpp"
#include "algorithms/queries.cpp"