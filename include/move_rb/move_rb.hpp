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
#include <misc/approximate_pattern_matching.cpp>
#include <misc/strings.hpp>
#include <move_r/move_r.hpp>
#include <move_rb/approximate_pattern_matching/search_schemes.hpp>
#include <move_rb/approximate_pattern_matching/sub_string.hpp>
#include <move_rb/approximate_pattern_matching/edit_distance_matrix.hpp>

enum move_rb_query_support_t : uint8_t {
    COUNT = 0,
    LOCATE = 1
};

// approximate-pattern-matching helper, defined completely outside the move_rb class
template <move_r_support support, typename sym_t, typename pos_t> class move_rb_apm_hamming;

// approximate-pattern-matching helper, defined completely outside the move_rb class
template <move_r_support support, typename sym_t, typename pos_t> class move_rb_apm_edit;

/**
 * @brief bi-directional move-r index, size O(r*(a/(a-1)))
 * @tparam support type of locate support (_locate_move or _locate_rlzsa)
 * @tparam sym_t value type (default: char for strings)
 * @tparam pos_t index integer type (use uint32_t if input size < UINT_MAX, else uint64_t)
 */
template <move_r_support support = _locate_move, typename sym_t = char, typename pos_t = uint32_t>
class move_rb {
    public:
    template <move_r_support, typename, typename> friend class move_rb_apm_hamming;
    template <move_r_support, typename, typename> friend class move_rb_apm_edit;
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
    static constexpr bool has_rlzsa = move_r_fwd_t::has_rlzsa; // true <=> the index has an rlzsa
    static constexpr bool has_lzendsa = move_r_fwd_t::has_lzendsa; // true <=> the index has an lzendsa
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
    void build(inp_t& input, move_r_params params = {});

    /**
     * @brief reverses T either in-memory or on disk
     * @param T string storing either T (if in_memory = true) or the name of the file containing T, else
     * @param p the number of threads to use
     * @param in_memory true <=> T is stored in-memory
     * @param log true <=> print log messages
     */
    void reverse(std::string& T, uint16_t p, bool in_memory, bool log = false);

    /**
     * @brief computes an sd-array marking every sample_rate-th sample of all num_values samples in the range [0, max_value)
     * @tparam fnc_t type of the sample function
     * @param num_values overall number of values
     * @param max_value maximum value
     * @param sample_rate sample rate
     * @param sample function, where sample(i) returns the ith sample
     * @return an sd-array marking every sample_rate-th sample (and max_value)
     */
    template <typename fnc_t>
    sd_array<pos_t> build_sampling(pos_t num_values, pos_t max_value, pos_t sample_rate, fnc_t sample);

    /**
     * @brief computes the index of the input interval containing position i
     * @tparam fnc_t type of the interval_start function
     * @param i position in [0, n)
     * @param sd_arr an sd-array marking exactly the starting positions of each x-th input interval of some move data structure
     * @param interval_start function, where interval_start(i) returns the starting position of the ith input interval of the same move data structure
     * @return the index of the input interval containing position i
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
        out << " sigma=" << sigma;
        out << " r=" << idx_fwd.num_bwt_runs();
        out << " r_rev=" << idx_bwd.num_bwt_runs();
        out << " r_=" << idx_fwd.M_LF().num_intervals();
        out << " r_rev_=" << idx_bwd.M_LF().num_intervals();

        if constexpr (supports_multiple_locate) {
            if constexpr (has_locate_move) {
                out << " r__=" << idx_fwd.M_Phi_m1().num_intervals();
                out << " r___=" << idx_fwd.M_Phi().num_intervals();
            } else if constexpr (has_rlzsa) {
                out << " z=" << idx_fwd.num_phrases_rlzsa();
                out << " z_l=" << idx_fwd.num_literal_phrases_rlzsa();
                out << " z_c=" << idx_fwd.num_copy_phrases_rlzsa();
            }
        }

        out << " size_index=" << size_in_bytes();
        idx_fwd.log_data_structure_sizes(out, "_fwd");
        out << " size_s_mlf_p_fwd=" << _S_MLF_p_fwd.size_in_bytes();
        
        if constexpr (support == _locate_move) {
            out << " size_s_mphi_p=" << _S_MPhi_p.size_in_bytes();
            out << " size_s_mphi_m1_p=" << _S_MPhi_m1_p.size_in_bytes();
        }

        idx_bwd.log_data_structure_sizes(out, "_bwd");
        out << " size_s_mlf_p_bwd=" << _S_MLF_p_bwd.size_in_bytes();
        out << " size_sa_m1_sr=" << _SA_sR_m1.size_in_bytes();
        out << " size_sa_m1_er=" << _SA_eR_m1.size_in_bytes();
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
        BEG = 0,
        MID = 1,
        END = 2
    };

    struct __attribute__((packed)) sample_info_pack_t {
        direction_t d; // LEFT/RIGHT <=> the sample can be computed from an SA-/Backward-SA-sample
        sample_t t; // BEG/MID/END <=> the SA-sample is at a run start/start/end
        pos_t o; // offset from the right of the SA-interval of the SA-sample (only it t = MID)
        pos_t i; // number of d-extensions since x has been set
        pos_t x; // index of the (sub-)run in L', whose start/end is the position of the SA-sample
        pos_t shft; // right shift (in the text) of the occurrences
        pos_t dpth; // depth (additional length) of this context in the banded alignment matrix
        bool rprtd; // true <=> this context has already been reported
    };

    public:
    /**
     * @brief stores the variables needed to perform bidirectional pattern search and locate queries
     */
    template <move_rb_query_support_t query_support>
    struct search_context_t {
        friend class move_rb;
        template <move_r_support, typename, typename> friend class move_rb_apm_hamming;
        template <move_r_support, typename, typename> friend class move_rb_apm_edit;
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

        std::conditional_t<query_support == LOCATE, sample_info_pack_t, std::monostate> s;

        template <direction_t dir, typename this_t>
        static auto vars_impl(this_t&& t)
        {
            if constexpr (dir == LEFT) return std::forward_as_tuple(t.b, t.e, t.b_R, t.e_R, t.b_, t.e_, t.b_R_, t.e_R_, t.s);
            else                       return std::forward_as_tuple(t.b_R, t.e_R, t.b, t.e, t.b_R_, t.e_R_, t.b_, t.e_, t.s);
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
            bool res = b == other.b && e == other.e;

            if constexpr (query_support == LOCATE) {
                res = res &&
                    s.shft == other.s.shft &&
                    s.dpth == other.s.dpth;
            }

            return res;
        }

        /**
         * @brief "less than" operator on search contexts
         * @param other the other search context to compare with
         * @return whether this context is considered "less than" the other
         */
        bool operator<(const search_context_t& other) const
        {
            if (b != other.b) return b < other.b;
            if (e != other.e) return e > other.e;

            if constexpr (query_support == LOCATE) {
                if (s.shft != other.s.shft) return s.shft < other.s.shft;
                if (s.dpth != other.s.dpth) return s.dpth < other.s.dpth;
            }

            return m < other.m;
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

                if constexpr (query_support == LOCATE) {
                    hash_combine<pos_t>(hash, pos_hash<pos_t>(ctx.s.shft));
                    hash_combine<pos_t>(hash, pos_hash<pos_t>(ctx.s.dpth));
                }
                
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
                s = {
                    .d = LEFT,
                    .t = BEG,
                    .o = 0,
                    .i = 0,
                    .x = 0,
                    .shft = 0,
                    .dpth = 0,
                    .rprtd = false
                };
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
            return s.dpth;
        }

        /**
         * @brief sets the depth (additional length) of this context in the banded alignment matrix
         * @param depth depth
         */
        inline void set_depth(pos_t depth) requires(query_support == LOCATE)
        {
            this->s.dpth = depth;
        }

        /**
         * @brief returns true <=> this context has already been reported
         * @return whether this context has already been reported
         */
        inline bool reported() const requires(query_support == LOCATE)
        {
            return s.rprtd;
        }

        /**
         * @brief sets whether this context has already been reported
         * @param reported whether this context has already been reported
         */
        inline void set_reported(bool reported) requires(query_support == LOCATE)
        {
            this->s.rprtd = reported;
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
         * @brief sets the right shift (in the text) of the occurrences to shft (only used for approximate pattern matching)
         * @param shft right shift (in the text) of the occurrences
         */
        inline void set_shift(pos_t shft) requires(query_support == LOCATE)
        {
            this->s.shft = shft;
        }
        
        inline pos_t shift() const requires(query_support == LOCATE)
        {
            return s.shft;
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
        ) const;

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
        inline void update_input_intervals(const move_rb<support, sym_t, pos_t>& idx);
        
        template <direction_t dir>
        inline void update_sample(const move_rb<support, sym_t, pos_t>& idx, search_context_t& ctx) const;

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
        inline extend_res_t extend(const move_rb<support, sym_t, pos_t>& idx, sym_t sym, direction_t dir);

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
        inline extend_res_t extend(const move_rb<support, sym_t, pos_t>& idx, sym_t sym);

        /**
         * @brief prepares the context to be extended
         * @param idx the index to query
         * @param dir extend direction
         * @return extend context prepared for extending the search context with every possible symbol in O(sigma) time
         */
        extend_context_t prepare_extend_all(const move_rb<support, sym_t, pos_t>& idx, direction_t dir);

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
        void prepare_extend_all(const move_rb<support, sym_t, pos_t>& idx, extend_context_t& ext_ctx);
        
        /**
         * @brief extends the context with the next symbol
         * @param idx the index to query
         * @param ext_ctx extend context to advance
         * @param dir extend direction
         * @return the search context extended by the next symbol
         */
        search_context_t extend_next(const move_rb<support, sym_t, pos_t>& idx, extend_context_t& ext_ctx, direction_t dir);

        /**
         * @brief extends the context with the next symbol
         * @tparam dir extend direction
         * @param idx the index to query
         * @param ext_ctx extend context to advance
         * @return the search context extended by the next symbol
         */
        template <direction_t dir>
        search_context_t extend_next(const move_rb<support, sym_t, pos_t>& idx, extend_context_t& ext_ctx) const;
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

        inline pos_t compute_center(const move_rb<support, sym_t, pos_t>& idx, const search_context_t<LOCATE>& ctx);

    protected:
        /**
         * @brief computes a first occurrence (SA[i]) of P, and adjusts the locate context accordingly
         * @param idx the index to query
         * @param ctx the search context associated with this locate context
         * @return an occurrence (SA[i]) of P
         */
        inline pos_t first_occ(const move_rb<support, sym_t, pos_t>& idx, const search_context_t<LOCATE>& ctx);

    public:
        /**
         * @brief reports the next occurrence of the currently matched P
         * @param idx the index to query
         * @param ctx the search context of the currently matched P
         * @return next occurrence
         */
        inline pos_t next_occ(const move_rb<support, sym_t, pos_t>& idx, const search_context_t<LOCATE>& ctx);

        /**
         * @brief locates the remaining (not yet reported) occurrences of the currently matched P
         * @tparam report_fnc_t type of the function report
         * @param idx the index to query
         * @param ctx the search context of the currently matched P
         * @param report function that is called with every (remaining) occurrence of the currently matched P
         */
        template <typename report_fnc_t>
        inline void locate(const move_rb<support, sym_t, pos_t>& idx, const search_context_t<LOCATE>& ctx, report_fnc_t report);

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
        void next_symbol(const move_rb<support, sym_t, pos_t>& idx);
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

    /**
     * @brief counts a pattern with at most k errors (w.r.t. hamming distance)
     * @param P the pattern to search
     * @param scheme the search scheme to use (provides k)
     * @return number of occurrences of P in T with at most k mismatches (w.r.t. hamming distance)
     */
    pos_t count_hamming_dist(const inp_t& P, const search_scheme_t& scheme) const
    {
        return move_rb_apm_hamming{*this}.count(P, scheme);
    }

    /**
     * @brief locates a pattern with at most k errors (w.r.t. dist_metr)
     * @tparam dist_metr distance metric (HAMMING_DISTANCE or EDIT_DISTANCE)
     * @tparam report_fnc_t type of the function report
     * @param P the pattern to search
     * @param scheme the search scheme to use (provides k)
     * @param report function that is called with every occurrence (occ, len, err) of P in T with at most k errors (w.r.t. dist_metr)
     */
    template <distance_metric_t dist_metr, typename report_fnc_t>
    void locate(const inp_t& P, const search_scheme_t& scheme, report_fnc_t report) const requires(supports_locate)
    {
        if constexpr (dist_metr == HAMMING_DISTANCE) move_rb_apm_hamming{*this}.locate(P, scheme, report);
        if constexpr (dist_metr == EDIT_DISTANCE) move_rb_apm_edit{*this}.locate(P, scheme, report);
    }

    /**
     * @brief locates a pattern with at most k errors (w.r.t. dist_metr)
     * @tparam dist_metr distance metric (HAMMING_DISTANCE or EDIT_DISTANCE)
     * @param P the pattern to search
     * @param scheme the search scheme to use (provides k)
     * @return occurrences of P in T with at most k errors (w.r.t. dist_metr)
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

#include "construction.cpp"
#include "queries.cpp"
#include "approximate_pattern_matching/hamming_distance.cpp"
#include "approximate_pattern_matching/edit_distance.cpp"
