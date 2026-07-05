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

#include <tsl/sparse_set.h>

#include <move_r/move_r.hpp>
// lightweight approximate-pattern-matching dependencies + query_support_t
#include <misc/apm.hpp>
#include <misc/fasta_sequence_data.hpp>

// approximate-pattern-matching helper, defined completely outside the move_rb class;
// generic over the index type idx_t (move_rb or a compatible bidirectional index)
template <typename idx_t, cigar_mode_t mode> class apm_hamming;

// approximate-pattern-matching helper, defined completely outside the move_rb class;
// generic over the index type idx_t (move_rb or a compatible bidirectional index)
template <typename idx_t, cigar_mode_t mode> class apm_edit;

/**
 * @brief bi-directional move-r index, size O(r*(a/(a-1)))
 * @tparam support type of locate support (_locate_move or _locate_rlzsa)
 * @tparam sym_t value type (default: char for strings)
 * @tparam pos_t index integer type (use uint32_t if input size < UINT_MAX, else uint64_t)
 */
template <move_r_support support = _locate_move, typename sym_t = char, typename pos_t = uint32_t>
class move_rb {
    public:
    template <typename, cigar_mode_t> friend class apm_hamming;
    template <typename, cigar_mode_t> friend class apm_edit;
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
    static constexpr bool has_rlzsa = move_r_fwd_t::has_rlzsa; // true <=> the index has an rlzsa
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

    using pos_type = pos_t; // index integer type (exposed for the generic APM helpers)
    using sym_type = sym_t; // value type (exposed for the generic APM helpers)

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

    interleaved_bit_aligned_vectors<pos_t> _SA_s_pos;
    interleaved_bit_aligned_vectors<pos_t> _SA_e_pos;

    // the multi-sequence (FASTA/DNA) metadata (start positions, names and separator symbol of the sequences), empty
    // unless the index was built from a multi-sequence input (see misc/fasta_sequence_data.hpp)
    fasta_sequence_data<pos_t, sym_t> _seq_data;

    // ############################# INTERNAL METHODS #############################

    /**
     * @brief returns the forward index (for dir == LEFT) or the backward index (for dir == RIGHT)
     * @tparam dir the search direction
     * @return the forward or backward move_r index
     */
    template <direction_t dir>
    const std::conditional_t<dir == LEFT, move_r_fwd_t, move_r_bwd_t>& index() const
    {
        if constexpr (dir == LEFT) {
            return idx_fwd;
        } else {
            return idx_bwd;
        }
    }

    /**
     * @brief returns the sd-array sampling the input-interval starting positions of M_LF in the given direction
     * @tparam dir the search direction
     * @return the sd-array of the forward (dir == LEFT) or backward (dir == RIGHT) M_LF
     */
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
     * @tparam seq_t the sequence type holding T in-memory (inp_t) or the file-name type (std::string) on disk
     * @param T the sequence storing T (if in_memory = true), or the name of the file containing T, else
     * @param p the number of threads to use
     * @param in_memory true <=> T is stored in-memory
     * @param log true <=> print log messages
     */
    template <typename seq_t>
    void reverse(seq_t& T, uint16_t p, bool in_memory, bool log = false);

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
        if constexpr (str_input) {
            if (params.file_input) {
                n = std::filesystem::file_size(input) + 1;
            } else {
                n = input.size() + 1;
            }
        } else {
            n = input.size() + 1;
        }

        if (params.mode == _suffix_array) {
            if (std::is_same_v<pos_t, uint32_t> && n <= INT_MAX) {
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
     * @brief returns the index's multi-sequence (FASTA/DNA) metadata (start positions, names and separator symbol of
     *        the sequences); empty (has_sequences() == false) unless the index was built from a multi-sequence input
     * @return the multi-sequence metadata
     */
    inline const fasta_sequence_data<pos_t, sym_t>& seq_data() const
    {
        return _seq_data;
    }

    /**
     * @brief attaches multi-sequence (FASTA/DNA) metadata to the index (as produced by process_fasta), enabling
     *        sequence-relative coordinates and names for SAM output; call after construction
     * @param seq_data the multi-sequence metadata
     */
    inline void set_fasta_sequence_data(fasta_sequence_data<pos_t, sym_t> seq_data)
    {
        _seq_data = std::move(seq_data);
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
            _SA_s_pos.size_in_bytes() +
            _SA_e_pos.size_in_bytes() +
            _seq_data.size_in_bytes();
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
        std::cout << "SA_s_pos: " << format_size(_SA_s_pos.size_in_bytes()) << std::endl;
        std::cout << "SA_e_pos: " << format_size(_SA_e_pos.size_in_bytes()) << std::endl;
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

        std::cout << std::endl << "Backward Index: " << std::endl;
        idx_bwd.log_data_structures();

        // log the forward index's SA_s / SA_e aligned in pairs with SA_s_pos / SA_e_pos
        if (!idx_fwd._SA_s.empty()) {
            std::cout << std::endl;
            aligned_log log;
            log.add_row("SA_s:", idx_fwd._SA_s.size(), [&](pos_t i) { return idx_fwd._SA_s[i]; });
            log.add_row("SA_s_pos:", _SA_s_pos.size(), [&](pos_t i) { return _SA_s_pos[i]; });
            log.print();
        }

        if (!idx_fwd._SA_e.empty()) {
            std::cout << std::endl;
            aligned_log log;
            log.add_row("SA_e:", idx_fwd._SA_e.size(), [&](pos_t i) { return idx_fwd._SA_e[i]; });
            log.add_row("SA_e_pos:", _SA_e_pos.size(), [&](pos_t i) { return _SA_e_pos[i]; });
            log.print();
        }
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
        out << " size_sa_s_pos=" << _SA_s_pos.size_in_bytes();
        out << " size_sa_e_pos=" << _SA_e_pos.size_in_bytes();
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
    template <query_support_t query_support> struct extend_context_t;

    protected:

    /** @brief whether an SA-sample sits at the start (BEG) or end (END) of a run */
    enum sample_t : uint8_t {BEG = 0, END = 1};

    /** @brief packed information needed to maintain a single SA-sample during a bidirectional search */
    struct __attribute__((packed)) sample_info_pack_t {
        direction_t d; // LEFT/RIGHT <=> the sample can be computed from an SA-/Backward-SA-sample
        sample_t t; // BEG/END <=> the SA-sample is at a run start/end
        pos_t o; // offset of the SA-sample from the left of the SA-interval
        pos_t i; // number of d-extensions since x has been set
        pos_t x; // index of the (sub-)run in L', whose start/end is the position of the SA-sample
    };

    /** @brief a tuple of a suffix array interval [b, e] and the indices b_, e_ of the input intervals containing b and e */
    struct prime_tuple_t {pos_t b, e, b_, e_;};

    public:

    /**
     * @brief stores the variables needed to perform bidirectional pattern search and locate queries
     * @tparam query_support whether the context supports COUNT or LOCATE queries
     */
    template <query_support_t query_support>
    struct search_context_t {
        friend class move_rb;
        template <typename, cigar_mode_t> friend class apm_hamming;
        template <typename, cigar_mode_t> friend class apm_edit;
        friend class locate_context_t;
        template <query_support_t _query_support> friend struct extend_context_t;
        template <query_support_t _query_support> friend class search_context_t;

    protected:
        // ############################# VARIABLES FOR THE SEARCH PHASE #############################

        const move_rb<support, sym_t, pos_t>* idx = nullptr; // the index this context queries

        direction_t dir_lst; // last performed pattern extend direction
        sym_t sym_lst; // last-added symbol
        pos_t err; // number of errors (only used for approximate pattern matching output)
        pos_t m; // length of the currently matched pattern P
        pos_t b, e; // [b, e] = SA-interval in T of the currently matched P
        pos_t b_R, e_R; // [b_R, e_R] = SA-interval in T^R of the reverse of the currently matched P
        pos_t b_, e_; // indices of the input intervals in M_LF containing b and e
        pos_t b_R_, e_R_; // indices of the input intervals in M_LF^R containing b_R and e_R

        // ############################# VARIABLES FOR MAINTAINING SA-SAMPLE INFORMATION DURING SEARCH #############################

        [[no_unique_address]] std::conditional_t<query_support == LOCATE, sample_info_pack_t, std::monostate> s;

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
            : idx(&idx)
        {
            reset();
        }

        /**
         * @brief equality operator on search contexts
         * @param other the other search context to test equality with
         * @return whether the other context other is equal to this search context; it is considered equal if it
         *         represents the same pattern (SA-interval and pattern length); other context variables may vary
         */
        bool operator==(const search_context_t& other) const
        {
            return b == other.b &&
                   e == other.e &&
                   m == other.m;
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
            inline pos_t operator()(const search_context_t& ctx) const
            {
                auto [b, e] = ctx.forward_sa_interval();
                pos_t hash = pos_hash<pos_t>(b);
                hash_combine<pos_t>(hash, pos_hash<pos_t>(e));
                hash_combine<pos_t>(hash, pos_hash<pos_t>(ctx.m));
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

        /**
         * @brief returns the direction of the last performed pattern extension
         * @return the last extension direction
         */
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
         * @brief returns an extend context for this search context, ready to be prepared for extending
         * @return an extend context for this search context
         */
        extend_context_t<query_support> extend_phase()
        {
            return extend_context_t<query_support>(*idx, *this);
        }

        /**
         * @brief resets the query context to an empty P
         */
        inline void reset()
        {
            dir_lst = NO_DIR;
            sym_lst = 0;
            m = 0;
            err = 0;

            b = 0;
            e = idx->idx_fwd.input_size();
            b_R = 0;
            e_R = e;
            b_ = 0;
            e_ = idx->idx_fwd.M_LF().num_intervals() - 1;
            b_R_ = 0;
            e_R_ = idx->idx_bwd.M_LF().num_intervals() - 1;

            if constexpr (query_support == LOCATE) {
                s = {
                    .d = LEFT,
                    .t = BEG,
                    .o = 0,
                    .i = 0,
                    .x = 0
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
         * @param prev output prev array
         * @param next output next array
         * @param max_sym maximum symbol to consider
         */
        template<direction_t dir>
        inline void build_prev_next(
            std::array<int64_t, 256>& prev, std::array<int64_t, 256>& next,
            i_sym_t max_sym
        ) const;

        /**
         * @brief update the input interval indexes b_ and e_, and adjusts dir-SA-interval samples
         * @tparam dir extend direction
         */
        template <direction_t dir>
        inline void update_input_intervals();

        template <direction_t dir>
        inline void update_sample(search_context_t& ctx, const prime_tuple_t& prime) const;

    public:
        using extend_res_t = std::tuple<search_context_t, bool>;

        /**
         * @brief extends the currently matched P with sym; if the extended P occurs in the input, true is
         * returned and the query context is adjusted to store the information for the extended P; else,
         * false is returned and the query context is not modified
         * @param sym the symbol to extend P with
         * @param dir extend direction
         * @return whether symP occurs in the input
         */
        inline extend_res_t extend(sym_t sym, direction_t dir);

        /**
         * @brief extends the currently matched P with sym; if the extended P occurs in the input, true is
         * returned and the query context is adjusted to store the information for the extended P; else,
         * false is returned and the query context is not modified
         * @tparam dir extend direction
         * @param sym the symbol to extend the pattern with
         * @return whether the extended P occurs in the input
         */
        template <direction_t dir>
        inline extend_res_t extend(sym_t sym);
    };

    /**
     * @brief structure storing the information for locating all occurrences of a search context
     */
    struct locate_context_t {
        protected:
        const move_rb<support, sym_t, pos_t>* idx = nullptr; // the index this context queries
        const search_context_t<LOCATE>* ctx = nullptr; // the search context this locate context belongs to

        direction_t dir; // current locate direction_t
        pos_t occ_rem; // number of remaining occurrences to locate
        pos_t c; // initial position in the suffix array
        pos_t SA_c; // initial suffix SA[c] in the suffix array interval
        pos_t i; // current position in the suffix array interval
        pos_t SA_i; // current suffix SA[i] in the suffix array interval

        // index of the input inteval of M_Phi/M_Phi^{-1} containing SA_i
        [[no_unique_address]] std::conditional_t<support == _locate_move, pos_t, std::monostate> s_;

        // rlzsa decode context (tracks its own position and suffix value while decoding)
        using rlzsa_ctx_t = typename rlzsa_opt<pos_t>::decode_context_t;

        [[no_unique_address]] std::conditional_t<support == _locate_rlzsa,
            rlzsa_ctx_t, std::monostate> rlz_l; // rlzsa context for decoding to the left
        [[no_unique_address]] std::conditional_t<support == _locate_rlzsa,
            rlzsa_ctx_t, std::monostate> rlz_r; // rlzsa context for decoding to the right

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
            : idx(ctx.idx), ctx(&ctx)
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

        /**
         * @brief computes the center position of the suffix array interval, from which locating starts
         * @return the center position of the suffix array interval
         */
        inline pos_t compute_center();

    protected:
        /**
         * @brief computes a first occurrence (SA[i]) of P, and adjusts the locate context accordingly
         * @return an occurrence (SA[i]) of P
         */
        inline pos_t first_occ();

    public:
        /**
         * @brief reports the next occurrence of the currently matched P
         * @return next occurrence
         */
        inline pos_t next_occ();

        /**
         * @brief locates the remaining (not yet reported) occurrences of the currently matched P
         * @tparam report_fnc_t type of the function report
         * @param report function that is called with every (remaining) occurrence of the currently matched P
         */
        template <typename report_fnc_t>
        inline void locate(report_fnc_t report);

        /**
         * @brief locates the remaining (not yet reported) occurrences of the currently matched P
         * @return vector containing the occurrences
         */
        std::vector<pos_t> locate()
        {
            std::vector<pos_t> Occ;
            Occ.reserve(occ_rem);
            locate([&](pos_t occ){Occ.emplace_back(occ);});
            return Occ;
        }
    };
    
    /**
     * @brief stores the variables needed to extend a search_context with all passible characters in O(sigma) time
     * @tparam query_support whether the underlying search context supports COUNT or LOCATE queries
     */
    template <query_support_t query_support>
    struct extend_context_t {
        friend class move_rb;
        friend struct search_context_t<query_support>;

        protected:
        const move_rb<support, sym_t, pos_t>* idx = nullptr; // the index this context queries
        search_context_t<query_support>* ctx = nullptr; // the search context this extend context belongs to
        direction_t dir; // direction the context is being extended in (set by prepare_extend_all)

        std::array<int64_t, 256> prev; // the prev array for all extensions
        std::array<int64_t, 256> next; // the next array for all extensions

        i_sym_t sym_nxt; // next symbol to extend the context with
        pos_t b_R_nxt; // current left interval limit of the suffix array interval of P(!dir) extended with sym_nxt in T(!dir)

        public:
        /**
         * @brief constructs an empty extend context
         */
        extend_context_t() {}

        /**
         * @brief constructs a new extend context for a given search context
         * @param idx the index to query
         * @param ctx the search context to extend
         */
        extend_context_t(const move_rb<support, sym_t, pos_t>& idx, search_context_t<query_support>& ctx)
            : idx(&idx), ctx(&ctx)
        { }

        /**
         * @brief returns the next symbol extend_next() will extend the context with (valid only while can_extend())
         * @return the next symbol to extend the context with
         */
        sym_t next_symbol() const {
            return idx->unmap_symbol(sym_nxt);
        }

        /**
         * @brief returns whether there is a symbol left to extend the context with
         * @return whether there is a symbol left to extend the context with
         */
        bool can_extend() const {
            return sym_nxt < idx->sigma;
        }

        /**
         * @brief prepares the context to be extended in direction dir; afterwards the search context can be
         *        extended with every possible symbol in O(sigma) time via extend_next()
         * @param dir extend direction
         * @return a reference to this extend context
         */
        extend_context_t& prepare_extend_all(direction_t dir);

        /**
         * @brief prepares the context to be extended
         * @tparam dir extend direction
         */
        template <direction_t dir>
        void prepare_extend_all();

        /**
         * @brief extends the context with the next symbol (in the direction set by prepare_extend_all)
         * @return the search context extended by the next symbol
         */
        search_context_t<query_support> extend_next();

        /**
         * @brief extends the context with the next symbol
         * @tparam dir extend direction
         * @return the search context extended by the next symbol
         */
        template <direction_t dir>
        search_context_t<query_support> extend_next();

        protected:
        /**
         * @brief advances sym_nxt to the next symbol the pattern can be extended with (or sigma, if it cannot be extended further)
         * @tparam dir direction_t the context should be extended with
         */
        template <direction_t dir>
        void advance_symbol();
    };

    /**
     * @brief returns a query context for the index
     * @tparam query_support whether the context supports COUNT or LOCATE queries
     * @return search_context_t
     */
    template <query_support_t query_support>
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
        return apm_hamming{*this}.count(P, scheme);
    }

    /**
     * @brief counts a pattern with at most k errors (w.r.t. hamming distance), broken down by number of mismatches
     * @param P the pattern to search
     * @param scheme the search scheme to use (provides k)
     * @return a histogram h of length k+1, where h[e] is the number of occurrences of P in T with exactly e mismatches
     */
    std::vector<pos_t> count_hamming_dist_histogram(const inp_t& P, const search_scheme_t& scheme) const
    {
        return apm_hamming{*this}.count_histogram(P, scheme);
    }

    /**
     * @brief locates a pattern with at most k errors (w.r.t. dist_metr)
     * @tparam dist_metr distance metric (HAMMING_DISTANCE or EDIT_DISTANCE)
     * @tparam report_fnc_t type of the function report
     * @param P the pattern to search
     * @param scheme the search scheme to use (provides k)
     * @param report function called with every occurrence (an aprx_occ_t<pos_t, mode>) of P in T with at most k
     *        errors (w.r.t. dist_metr). If mode == CIGAR the occurrence additionally carries a CIGAR alignment of P
     *        against its matched string (occ.cigar), and its error is the exact distance of P to that string.
     */
    template <distance_metric_t dist_metr, cigar_mode_t mode = NO_CIGAR, typename report_fnc_t>
    void locate(const inp_t& P, const search_scheme_t& scheme, report_fnc_t report) const requires(supports_locate)
    {
        if constexpr (dist_metr == HAMMING_DISTANCE) apm_hamming<move_rb, mode>{*this}.locate(P, scheme, report);
        if constexpr (dist_metr == EDIT_DISTANCE) apm_edit<move_rb, mode>{*this}.locate(P, scheme, report);
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

        _SA_s_pos.serialize(out);
        _SA_e_pos.serialize(out);

        _seq_data.serialize(out);
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

        _SA_s_pos.load(in);
        _SA_e_pos.load(in);

        _seq_data.load(in);
    }

    /**
     * @brief stores the index to an output stream
     * @param os output stream
     * @return the output stream
     */
    std::ostream& operator>>(std::ostream& os) const
    {
        serialize(os);
        return os;
    }

    /**
     * @brief reads a serialized index from an input stream
     * @param is input stream
     * @return the input stream
     */
    std::istream& operator<<(std::istream& is)
    {
        load(is);
        return is;
    }
};

#include "construction.tpp"
#include "queries.tpp"
#include <algorithms/apm/hamming_distance.tpp>
#include <algorithms/apm/edit_distance.tpp>
