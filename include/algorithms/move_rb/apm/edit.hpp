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

/**
 * @file edit.hpp
 * @brief approximate pattern matching with edit distance on a bidirectional index with search schemes
 */

#ifndef MOVE_R_ALGORITHMS_MOVE_RB_APM_EDIT_HPP
#define MOVE_R_ALGORITHMS_MOVE_RB_APM_EDIT_HPP

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <gtl/phmap.hpp>
#include <misc/apm.hpp>
#include <misc/sort.hpp>
#include <misc/utils.hpp>
#include <algorithms/move_rb/apm/model.hpp>
#include <algorithms/move_rb/apm/model_plan.hpp>
#include <algorithms/move_rb/apm/partition.hpp>
#include <algorithms/move_rb/apm/search_base.hpp>

// the steps of one search of a search scheme: their parts, error bounds and directions
class apm_search_steps {
protected:
    const search_scheme_t& _scheme; // the search scheme
    const std::vector<search_step_t>& _steps; // [0..p-1] the steps of the search
    std::vector<direction_t> _dirs; // [0..p-1] the direction of step i (see apm_step_directions)
    std::vector<bool> _switches_dir; // [0..p-1] true <=> step i has another direction than step i - 1 (false for i = 0)

public:
    /**
     * @brief constructs the steps of a search
     * @param scheme the search scheme
     * @param search_idx the index of the search in the scheme
     */
    apm_search_steps(const search_scheme_t& scheme, uint8_t search_idx)
        : _scheme(scheme), _steps(scheme.S[search_idx]), _dirs(apm_step_directions(_steps))
    {
        _switches_dir.reserve(scheme.p);
        _switches_dir.emplace_back(false);
        for (uint8_t i = 1; i < _dirs.size(); i++) _switches_dir.emplace_back(_dirs[i] != _dirs[i - 1]);
    }

    // returns the lower bound of step i (the minimum number of errors after step i)
    uint8_t k_min(uint8_t i) const { return _steps[i].k_min; }

    // returns the upper bound of step i (the maximum number of errors after step i)
    uint8_t k_max(uint8_t i) const { return _steps[i].k_max; }

    // returns the part of P that step i matches
    uint8_t part(uint8_t i) const { return _steps[i].part; }

    // returns the direction of step i
    direction_t dir(uint8_t i) const { return _dirs[i]; }

    // returns whether step i has another direction than step i - 1
    bool switches_dir(uint8_t i) const { return _switches_dir[i]; }

    // returns whether the part of step i is the first or the last part of P
    bool is_outer_part(uint8_t i) const { return _steps[i].part == 0 || _steps[i].part == _scheme.p - 1; }
};

/* sort key of a search context: its suffix array interval [beg, end] in T, its occurrence length (see recompute_errors)
   and its index in the unsorted vector of contexts; the keys are ordered by beg, then by decreasing end (an interval
   before the intervals nested in it), then by len and idx */
struct apm_context_order_t {
    uint64_t beg;
    uint64_t end;
    uint64_t len;
    uint32_t idx;

    inline bool operator<(const apm_context_order_t& other) const
    {
        if (beg != other.beg) return beg < other.beg;
        if (end != other.end) return end > other.end;
        if (len != other.len) return len < other.len;
        return idx < other.idx;
    }
};

/**
 * @brief approximate pattern matching with edit distance on a bidirectional index with search schemes (Lam et al.;
 *        Kucherov, Salikhov, Tsur; Kianfar et al.), with the banded bit-parallel matrices and the branch and bound of
 *        Columba (Renders, Depuydt, Fostier)
 *
 * 1. P is partitioned into the parts of the search scheme (see apm_partition).
 * 2. Each search matches its parts step by step: for each step, a depth-first search over the extensions of the
 *    search context in the direction of the step computes one row of the edit-distance matrix of the part per text
 *    symbol (see search_part and branch_and_bound). The cells of the last column of the matrix (the final column)
 *    of a path are collected in a cluster_t; the next step continues from the local minima of the final column (see
 *    flush_cluster).
 * 3. Each search context stores its matched text string as pieces (see matched_piece_t); after the search, the error
 *    count and the length of each context are recomputed by aligning P against its string (see recompute_errors). As
 *    with an error model, an alignment does not start or end with a deletion.
 * 4. Each suffix array position is reported once, with the context with the fewest errors that contains it (see
 *    report_contexts). The occurrences that overlap a text end are found separately (see locate_at_text_ends); of an
 *    occurrence found by the search and one found at a text end at the same position, only the one with fewer errors
 *    is reported. An occurrence that ends at a separator of the sequences is found by the search (see
 *    cluster_t::is_local_min).
 *
 * @tparam idx_t the bidirectional index type (move_rb, or a benchmark adapter with the same interface)
 */
template <typename idx_t>
class apm_edit : public apm_search_base<idx_t> {
protected:
    // ############################# TYPES AND CONSTANTS #############################

    using base_t = apm_search_base<idx_t>;
    using typename base_t::inp_t;
    using base_t::_index;
    using base_t::_params;
    using base_t::next_is_separator;

    template <query_support_t query_support>
    using search_context_t = typename base_t::template search_context_t<query_support>;

    // true <=> both text ends can be decoded with the index (see text_end); false for the benchmark adapters
    static constexpr bool can_decode_text_ends = requires(const idx_t& i) {
        i.backward_index().M_LF();
        i.forward_index().M_LF();
    };

    // number of bits of a piece index and of a position in P in a matched_piece_t
    static constexpr uint64_t piece_bits = 8 * uint40_t::width;
    // the largest piece index: the prev of the first piece of a chain, and the last piece of an empty chain
    static constexpr uint64_t no_piece = (uint64_t(1) << piece_bits) - 1;

    // the CIGAR index of an occurrence without a CIGAR
    static constexpr uint32_t no_cigar = std::numeric_limits<uint32_t>::max();

    /* a piece of a matched text string: a copy of P[pos, pos + len) or a single text symbol (a literal); the pieces
       of a string form a chain from its last appended piece back to its first one (see node_t) */
    struct __attribute__((packed)) matched_piece_t {
        using prev_field = next_bit_field_t<void, piece_bits>;
        using pos_field = next_bit_field_t<prev_field, piece_bits>;

        packed_bit_fields<bytes_of_bit_fields<pos_field>> bits; // prev and pos (in the layout of GCC bit-fields)
        uint32_t len;       // copy: its length (>= 1); literal: 0
        uint16_t refs;      // reference count: the number of pieces whose prev is this piece and of piece_ref objects that refer to it
        uint8_t sym;        // literal: the text symbol; copy: unused
        direction_t dir;    // LEFT: the piece was prepended to the string, RIGHT: appended

        matched_piece_t() = default;

        matched_piece_t(uint64_t prev, uint64_t pos, uint32_t len, uint16_t refs, uint8_t sym, direction_t dir)
            : len(len), refs(refs), sym(sym), dir(dir)
        { bits.template assign<prev_field, pos_field>(prev, pos); }

        // the index in _matched_pieces of the piece added to the string before this one, or no_piece
        uint64_t prev() const { return bits.template get<prev_field>(); }
        void set_prev(uint64_t value) { bits.template set<prev_field>(value); }

        // copy: the start of the copied substring of P; literal: unused
        uint64_t pos() const { return bits.template get<pos_field>(); }
        void set_pos(uint64_t value) { bits.template set<pos_field>(value); }
    };

    static_assert(sizeof(matched_piece_t) == 2 * uint40_t::width + sizeof(uint32_t) + sizeof(uint16_t)
        + sizeof(uint8_t) + sizeof(direction_t));

    /* a reference-counted index of a piece in _matched_pieces of _piece_owner (or no_piece): copying it increments and
       destroying it decrements the reference count of the piece; a piece whose count drops to 0 is freed, and the
       reference count of its prev is decremented in the same way (see dec_ref) */
    struct piece_ref {
        uint40_t idx = uint40_t(no_piece);

        piece_ref() = default;
        explicit piece_ref(uint64_t i) : idx(i) { if (idx != no_piece) _piece_owner->inc_ref(idx); }
        piece_ref(const piece_ref& o) : idx(o.idx) { if (idx != no_piece) _piece_owner->inc_ref(idx); }
        piece_ref(piece_ref&& o) noexcept : idx(o.idx) { o.idx = no_piece; }

        piece_ref& operator=(const piece_ref& o)
        {
            if (this != &o) {
                if (o.idx != no_piece) _piece_owner->inc_ref(o.idx);
                if (idx != no_piece) _piece_owner->dec_ref(idx);
                idx = o.idx;
            }
            return *this;
        }

        piece_ref& operator=(piece_ref&& o) noexcept
        {
            if (this != &o) {
                if (idx != no_piece) _piece_owner->dec_ref(idx);
                idx = o.idx;
                o.idx = no_piece;
            }
            return *this;
        }

        ~piece_ref() { if (idx != no_piece) _piece_owner->dec_ref(idx); }
    };

    // identifies the text string of a context: its suffix array interval [beg, end] in T and its length
    struct string_key_t {
        uint64_t beg;
        uint64_t end;
        uint64_t len;

        bool operator==(const string_key_t&) const = default;
    };

    // the hash function of a string_key_t
    struct string_key_hash {
        uint64_t operator()(const string_key_t& key) const
        {
            uint64_t h = hash_64(key.beg);
            hash_combine(h, hash_64(key.end));
            hash_combine(h, hash_64(key.len));
            return h;
        }
    };

    /* a search context with its matched text string: the chain of pieces ending at last_piece, together with the
       pending copy P[copy_pos, copy_pos + copy_len) that is not yet a piece (see commit_copy); the lengths and positions
       fit into 32 bits, since m + k <= aprx_occ_max_len */
    struct node_t {
        search_context_t<LOCATE> ctx;
        uint32_t len = 0;      // the matched length, in [0, m + k] (during a step, the row in the matrix of the part)
        uint32_t copy_pos = 0; // the start of the pending copy in P
        uint32_t copy_len = 0; // the length of the pending copy (0 <=> there is none)
        piece_ref last_piece;  // the last piece of the chain
        uint8_t copy_dir : 2;  // the direction in which the pending copy was matched (LEFT or RIGHT)
        uint8_t was_taken : 1; // whether the cell of the node in a cluster_t has been taken (see cluster_t::take_cell)

        node_t() : copy_dir(RIGHT), was_taken(false) {}
        node_t(const search_context_t<LOCATE>& c) : ctx(c), copy_dir(RIGHT), was_taken(false) {}

        uint64_t match_len() const { return len; }
        void set_match_len(uint64_t d) { len = d; }
        bool taken() const { return was_taken; }
        void set_taken(bool t) { was_taken = t; }
        uint64_t last_piece_idx() const { return last_piece.idx; }
        void set_last_piece(uint64_t i) { last_piece = piece_ref(i); }

        // returns the key of the node in a node_set_t: its suffix array interval in T and its matched length
        string_key_t key() const
        {
            auto [b, e] = ctx.forward_sa_interval();
            return { b, e, len };
        }
    };

    // the contexts found by a search, by their keys (see node_t::key)
    using node_set_t = apm_context_set<node_t, string_key_t, string_key_hash>;

    // a reference to a node, in CIGAR mode with the index of its CIGAR
    template <cigar_mode_t mode>
    using context_ref_t = std::conditional_t<mode == CIGAR, std::pair<const node_t*, uint32_t>, const node_t*>;

    // the data structures of one search of the scheme that search_part, flush_cluster and branch_and_bound share
    template <typename matrix_word_t>
    struct search_data_t {
        const apm_search_steps& steps; // the steps of the search
        const std::vector<directional_substring<inp_t>>& parts; // [0..p-1] the parts of P
        std::vector<std::vector<node_t>>& stacks; // [0..p-1] the stack of the depth-first search of step i
        std::vector<edit_distance_matrix<matrix_word_t>>& matrices_fwd; // [0..p-1] the matrix of part i read to the right
        std::vector<edit_distance_matrix<matrix_word_t>>& matrices_bwd; // [0..p-1] the matrix of part i read to the left
        node_set_t& result; // the contexts found by the searches
    };

    /* a range [beg, end] of suffix array positions in T that are reported with the same occurrence length, error
       count and CIGAR index (see effective_intervals) */
    struct __attribute__((packed)) occurrence_interval_t {
        using beg_field = next_bit_field_t<void, aprx_occ_pos_bits>;
        using end_field = next_bit_field_t<beg_field, aprx_occ_pos_bits>;
        using len_field = next_bit_field_t<end_field, aprx_occ_len_bits>;
        using err_field = next_bit_field_t<len_field, aprx_occ_err_bits>;

        // beg, end, len and err (in the layout of GCC bit-fields)
        packed_bit_fields<bytes_of_bit_fields<err_field>> bits;
        uint32_t cigar_idx = no_cigar;

        occurrence_interval_t(uint64_t beg, uint64_t end, uint64_t len, uint64_t err, uint32_t cigar_idx)
            : cigar_idx(cigar_idx)
        { bits.template assign<beg_field, end_field, len_field, err_field>(beg, end, len, err); }

        uint64_t beg() const { return bits.template get<beg_field>(); }
        uint64_t end() const { return bits.template get<end_field>(); }
        uint64_t len() const { return bits.template get<len_field>(); }
        uint64_t err() const { return bits.template get<err_field>(); }
    };

    static_assert(sizeof(occurrence_interval_t)
        == (2 * aprx_occ_pos_bits + aprx_occ_len_bits + aprx_occ_err_bits) / 8 + sizeof(uint32_t));

    /* the final column of the matrix of a part on one path: the cells of the rows num_rows() - size() to
       num_rows() - 1, each with the context of its row and its edit distance */
    class cluster_t {
    protected:
        uint64_t _start_depth; // the matched length before the part
        uint16_t _k_max; // the upper bound of the step
        uint64_t _last_cell; // the cell set last, or uint64_t(-1) if none is set (the cells after it belong to other paths)
        std::vector<node_t> _nodes; // [0..size()-1] the context of cell i
        std::vector<uint16_t> _dists; // [0..size()-1] the edit distance of cell i (_k_max + 1 if it is not set)

    public:
        // creates a final column with size cells, none of them set, for a step with the upper bound k_max
        cluster_t(uint64_t size, uint64_t k_max, uint64_t start_depth)
            : _start_depth(start_depth), _k_max(k_max), _last_cell(uint64_t(-1)),
              _nodes(size), _dists(size, k_max + 1) {}

        // sets cell idx to a context (not taken) and its edit distance
        void set_cell(uint64_t idx, node_t node, uint16_t dist)
        {
            _nodes[idx] = std::move(node);
            _nodes[idx].set_taken(false);
            _dists[idx] = dist;
            _last_cell = idx;
        }

        // returns the number of cells of the column
        uint64_t size() const { return _dists.size(); }

        /**
         * @brief returns whether the distance of cell i is at most the distances of its neighbours (up to _last_cell); a
         *        smaller distance in cell i + 1 is ignored if an occurrence of the string of cell i is followed by the
         *        separator in the direction of the step (the string of cell i + 1 does not contain this occurrence)
         * @tparam sep_fnc_t type of the function at_separator
         * @param i the cell
         * @param at_separator returns whether the string of a node can be extended with the separator in the
         *        direction of the step (see apm_search_base::extends_with_separator)
         */
        template <typename sep_fnc_t>
        bool is_local_min(uint64_t i, sep_fnc_t& at_separator)
        {
            if (i != 0 && _dists[i] > _dists[i - 1]) return false;
            if (i == _last_cell || _dists[i] <= _dists[i + 1]) return true;
            return ends_at_separator(i, at_separator);
        }

        /**
         * @brief returns whether an occurrence of the string of cell i (i < _last_cell) is followed by the separator in
         *        the direction of the step; the string of cell i + 1 then does not contain this occurrence
         * @tparam sep_fnc_t type of the function at_separator
         * @param i the cell
         * @param at_separator see is_local_min
         */
        template <typename sep_fnc_t>
        bool ends_at_separator(uint64_t i, sep_fnc_t& at_separator)
        {
            /* the string of cell i + 1 extends the string of cell i; if it has as many occurrences, no occurrence of the
               string of cell i is followed by the separator */
            return num_occ(_nodes[i]) > num_occ(_nodes[i + 1]) && at_separator(_nodes[i]);
        }

        // returns the number of occurrences of the string of a node
        static uint64_t num_occ(const node_t& node)
        {
            auto [b, e] = node.ctx.forward_sa_interval();
            return e - b + 1;
        }

        // marks cell i as taken and returns a copy of its node, with the distance of the cell as its error count and with
        // its matched length increased by _start_depth
        node_t take_cell(uint64_t i)
        {
            _nodes[i].set_taken(true);
            node_t node = _nodes[i];
            node.ctx.set_errors(_dists[i]);
            node.set_match_len(node.match_len() + _start_depth);
            return node;
        }

        // takes (see take_cell) and returns the cells up to _last_cell that are local minima (see is_local_min) with at
        // most _k_max errors and were not taken before
        template <typename sep_fnc_t>
        std::vector<node_t> local_minima(sep_fnc_t& at_separator);

        /**
         * @brief chooses the cell from which the next step continues in the same direction: the first local minimum
         *        (see is_local_min) with a distance in [k_min, _k_max]
         * @param k_min the lower bound of the step
         * @param descendants is set to the contexts of the cells after it
         * @param init_dists is set to the first column of the matrix of the next step: the distances of the cell and
         *        of the cells after it, with every local minimum below k_min raised (see the body)
         * @param at_separator see is_local_min
         * @return the cell (see take_cell) and true, or false if there is no such cell
         */
        template <typename sep_fnc_t>
        std::tuple<node_t, bool> choose_center(uint16_t k_min, std::vector<node_t>& descendants,
            std::vector<uint16_t>& init_dists, sep_fnc_t& at_separator);
    };

    /* a member of a context in search_model: an alignment against the string of the context extended by the extension
       ext (see _ext_trie_parent) in the direction of the members of the search (see search_model) */
    struct model_member_t {
        uint32_t ext; // the extension
        uint16_t err; // the error count of the alignment
        uint16_t cost; // the cost of the alignment
    };

    // a context reached after a step of a search in search_model
    struct model_entry_t {
        node_t node; // the context; node.ctx.errors() is its smallest error count
        uint16_t cost; // its smallest cost (see apm_cost_band)
        bool covered = false; // true <=> in the next step, the context is continued in the search from a shorter context
        // empty if the entry represents only the string of the context, else the alignments it represents (node.ctx.errors()
        // and cost are their smallest error count and cost)
        std::vector<model_member_t> members;
        uint64_t counters = 0; // without members: the packed counters of an alignment with the cost (see apm_cost_band)
        /* without members: the smallest value of each counter over the alignments of the string (see
           apm_cost_band::final_min_counters); 0 with members */
        uint64_t min_counters = 0;
    };

    // a member of the contexts on a path of process_part: its error count and its cost are those of the cell of the path
    // plus err_offset and cost_offset
    struct path_member_t {
        uint32_t ext; // the extension (see model_member_t)
        uint16_t err_offset;
        uint16_t cost_offset;
    };

    // the contexts reached after a step of a search, by the keys of their text strings (see string_key_t)
    using string_map_t = gtl::flat_hash_map<string_key_t, model_entry_t, string_key_hash>;

    // ############################# DATA MEMBERS #############################

    /* the apm_edit whose _matched_pieces the piece_ref objects of this thread refer to; set by locate, search and
       search_model before they create a piece_ref (a query is run on one thread) */
    static inline thread_local const apm_edit* _piece_owner = nullptr;

    // the pieces of the matched strings of all contexts of a query: [0..number of pieces - 1] (see matched_piece_t)
    mutable std::vector<matched_piece_t> _matched_pieces;

    // the indices of the freed pieces in _matched_pieces, which are reused first (a stack)
    mutable std::vector<uint64_t> _free_pieces;

    /* the largest distance between a position j of P and the position of P on the diagonal of a node at which
       append_last_symbol looks for a P[j] that equals the appended symbol; set by search (k) and search_model
       (max_errors()) */
    mutable uint64_t _max_diagonal_shift = 0;

    /* the extensions of the members of search_model (see model_member_t), a trie: [0..number of extensions - 1]
       extension 0 is empty, extension x > 0 is extension _ext_trie_parent[x] followed by the symbol
       _ext_trie_symbol[x]; _ext_trie_child maps ext_trie_key(parent, symbol) to the extension */
    mutable std::vector<uint32_t> _ext_trie_parent;
    mutable std::vector<uint8_t> _ext_trie_symbol;
    mutable gtl::flat_hash_map<uint64_t, uint32_t> _ext_trie_child;

    // ############################# INTERNAL METHODS: TEXT ENDS #############################

    /**
     * @brief returns the last symbols of the text of a sub-index, decoded with LF steps from BWT position 0
     * @tparam sub_t type of the sub-index
     * @param sub the sub-index
     * @param length number of symbols to decode
     * @return the last length symbols of the text, in text order
     */
    template <typename sub_t>
    static inp_t text_end(const sub_t& sub, uint64_t length);

    /**
     * @brief returns the windows of the text in which an occurrence with at most k errors can overlap a text end:
     *        the first and the last min(n, m + k) symbols (of the first and the last sequence), in text order; the
     *        first window ends before its first separator, the last one starts after its last separator
     * @param m the length of P
     * @param k the maximum number of errors
     * @return the position of the first window in the text, the first window and the last window (it ends at the text
     *         end)
     */
    std::tuple<uint64_t, inp_t, inp_t> text_end_windows(uint64_t m, uint64_t k) const;

    /**
     * @brief reports the occurrences that overlap a text end (the search can miss them) by aligning P against the
     *        windows of text_end_windows (as in recompute_errors)
     * @tparam mode whether a CIGAR is computed for each occurrence
     * @tparam report_fnc_t type of the report function
     * @param P the pattern
     * @param m length of P
     * @param k maximum number of errors
     * @param report function called with each occurrence
     * @param cigars the CIGARs are appended to it (in CIGAR mode)
     */
    template <cigar_mode_t mode, typename report_fnc_t>
    void locate_at_text_ends(const inp_t& P, uint64_t m, uint64_t k, report_fnc_t& report,
        std::vector<cigar_t<uint8_t>>& cigars) const;

    // locate_at_text_ends with a matrix of word type word_t
    template <cigar_mode_t mode, typename word_t, typename report_fnc_t>
    void locate_at_text_ends_impl(const inp_t& P, uint64_t m, uint64_t k, report_fnc_t& report,
        std::vector<cigar_t<uint8_t>>& cigars) const;

    // ############################# INTERNAL METHODS: MATCHED STRINGS #############################

    /**
     * @brief returns a locate context of a search context that reports each occurrence with its suffix array position
     *        (see move_rb::locate_output_t), if the index supports it
     * @tparam ctx_t type of the search context
     * @param ctx a search context
     */
    template <typename ctx_t>
    static auto locate_phase_with_positions(const ctx_t& ctx)
    {
        if constexpr (requires { ctx.template locate_phase<idx_t::_occ_and_sa_pos>(); }) {
            return ctx.template locate_phase<idx_t::_occ_and_sa_pos>();
        } else {
            return ctx.locate_phase();
        }
    }

    // returns the node of a context reference
    template <cigar_mode_t mode>
    static const node_t& node_of(const context_ref_t<mode>& e) { if constexpr (mode == CIGAR) return *e.first; else return *e; }

    // returns the CIGAR index of a context reference (no_cigar in NO_CIGAR mode)
    template <cigar_mode_t mode>
    static uint32_t cigar_idx_of(const context_ref_t<mode>& e) { if constexpr (mode == CIGAR) return e.second; else return no_cigar; }

    // returns the index of a new piece in _matched_pieces (a freed one if there is one)
    uint64_t allocate_piece() const
    {
        if (!_free_pieces.empty()) { uint64_t i = _free_pieces.back(); _free_pieces.pop_back(); return i; }
        _matched_pieces.emplace_back();
        return _matched_pieces.size() - 1;
    }

    // increments the reference count of a piece; a count of UINT16_MAX is not changed (the piece is never freed)
    void inc_ref(uint64_t idx) const { if (_matched_pieces[idx].refs != UINT16_MAX) _matched_pieces[idx].refs++; }

    /* decrements the reference count of a piece; a piece whose count drops to 0 is freed, and the reference count of
       its prev is decremented in the same way */
    void dec_ref(uint64_t idx) const;

    // appends a literal to the chain ending at prev (dir == LEFT: prepends it) and returns its index
    uint64_t append_literal(uint64_t prev, uint8_t sym, direction_t dir) const;

    // appends a copy of P[pos, pos + len) to the chain ending at prev (dir == LEFT: prepends it) and returns its index
    uint64_t append_copy(uint64_t prev, uint64_t pos, uint64_t len, direction_t dir) const;

    // copies the matched string (chain and pending copy) of a node to another node; the symbol of an extension is
    // appended when its node is taken from the stack (see append_last_symbol)
    void copy_matched_string(node_t& child, const node_t& parent) const
    {
        child.set_last_piece(parent.last_piece_idx());
        child.copy_pos = parent.copy_pos;
        child.copy_len = parent.copy_len;
        child.copy_dir = parent.copy_dir;
    }

    // returns the child of a node with the extended context ctx: its matched length is that of the node plus 1, its
    // matched string is that of the node (see copy_matched_string)
    node_t child_of(const node_t& node, const search_context_t<LOCATE>& ctx) const
    {
        node_t child(ctx);
        child.set_match_len(node.match_len() + 1);
        copy_matched_string(child, node);
        return child;
    }

    // makes a node whose context was extended in place its own child: its matched length is increased by 1, its
    // matched string is not changed (see child_of), and it is not taken (see cluster_t::take_cell)
    static void descend(node_t& node)
    {
        node.set_match_len(node.match_len() + 1);
        node.set_taken(false);
    }

    // appends the pending copy of a node (if there is one) to its chain as a piece
    void commit_copy(node_t& node) const
    {
        if (node.copy_len > 0) {
            node.set_last_piece(append_copy(node.last_piece_idx(), node.copy_pos, node.copy_len, direction_t(node.copy_dir)));
            node.copy_len = 0;
        }
    }

    // returns a copy of a node whose pending copy is appended to its chain (see commit_copy)
    node_t with_committed_copy(node_t node) const
    {
        commit_copy(node);
        return node;
    }

    // appends a matched symbol that equals P[pos] (dir == LEFT: prepends it): to the pending copy if pos continues it,
    // else as a new pending copy (after committing the old one)
    void append_pattern_symbol(node_t& node, uint64_t pos, direction_t dir) const;

    /**
     * @brief appends the symbol of the last extension of a node in direction dir: to the pending copy if it continues
     *        the copy, else as a new pending copy of P[j] for the position j with P[j] = the symbol that is closest to
     *        the position of P on the diagonal of the node (at most _max_diagonal_shift positions away), else as a
     *        literal
     * @param node the node
     * @param part the part of the current step
     * @param dir the direction of the step
     */
    void append_last_symbol(node_t& node, const directional_substring<inp_t>& part, direction_t dir) const;

    /**
     * @brief decodes the matched string of a chain of pieces
     * @param P the pattern
     * @param last the index of the last piece of the chain
     * @param out is set to the string, in text order
     */
    void decode_matched_string(const inp_t& P, uint64_t last, inp_t& out) const;

    /**
     * @brief computes the error count and the occurrence length of each context: P is aligned against the prefix of the
     *        matched string of the context with the fewest errors (the shortest one on ties, so the alignment does not
     *        end with a deletion), without a deletion before P[0]; the contexts with at most k errors are collected,
     *        except those whose matched string without its first symbol has fewer errors
     * @tparam mode whether a CIGAR is computed for each context
     * @tparam word_t word type of the edit-distance matrix
     * @param ctxts the contexts
     * @param P the pattern
     * @param m length of P
     * @param k maximum number of errors
     * @param ctxts_sorted the contexts with at most k errors are appended to it
     * @param cigars their CIGARs are appended to it (in CIGAR mode)
     */
    template <cigar_mode_t mode, typename word_t>
    void recompute_errors(node_set_t& ctxts, const inp_t& P, uint64_t m, uint64_t k,
        std::vector<context_ref_t<mode>>& ctxts_sorted, std::vector<cigar_t<uint8_t>>& cigars) const;

    // ############################# INTERNAL METHODS: REPORTING #############################

    /**
     * @brief reports the occurrences of the contexts: each suffix array position is reported once, with the context
     *        of the fewest errors (the longest on ties) that contains it
     * @tparam mode whether each occurrence is reported with a CIGAR index
     * @tparam report_fnc_t type of the report function
     * @param ctxts_sorted the contexts (their order is changed)
     * @param report the report function
     */
    template <cigar_mode_t mode, typename report_fnc_t>
    void report_contexts(std::vector<context_ref_t<mode>>& ctxts_sorted, report_fnc_t& report) const;

    /**
     * @brief splits the union of the suffix array intervals of contexts into the maximal intervals in which every
     *        position is reported with the same context (the one with the fewest errors, the longest on ties)
     * @param ctxts_sorted the contexts, sorted by the starting positions of their suffix array intervals
     * @return the intervals, sorted by their starting positions (strictly increasing)
     */
    template <cigar_mode_t mode>
    static std::vector<occurrence_interval_t> effective_intervals(const std::vector<context_ref_t<mode>>& ctxts_sorted);

    // ############################# INTERNAL METHODS: SEARCH #############################

    /**
     * @brief extends the context of a node in direction dir by each symbol except the separator: if it has one
     *        extension, the node becomes its child in place (see descend) and true is returned; else push(child) is
     *        called for each extension in the order of the symbols (see apm_dfs_children), where child is a child of
     *        the node (see child_of), and false is returned
     * @tparam push_fnc_t type of the function push
     * @param node the node
     * @param dir the direction
     * @param push the function called for each extension that is pushed
     * @return whether the node became its child
     */
    template <typename push_fnc_t>
    bool extend_node(node_t& node, direction_t dir, push_fnc_t& push) const;

    /**
     * @brief sets the input of a matrix to a part; with an error model, a free symbol of a position is treated as
     *        equal to its symbol in P
     * @tparam matrix_word_t word type of the matrix
     * @param matrix the matrix
     * @param part the part
     */
    template <typename matrix_word_t>
    void set_matrix_input(edit_distance_matrix<matrix_word_t>& matrix, const directional_substring<inp_t>& part) const;

    /**
     * @brief executes a search scheme
     * @tparam matrix_word_t word type of the edit-distance matrices
     * @param P the pattern
     * @param scheme the search scheme
     * @return the contexts found by the searches, without duplicates
     */
    template <typename matrix_word_t>
    node_set_t search(const inp_t& P, const search_scheme_t& scheme) const;

    /**
     * @brief matches the part of step p_idx from a context: a depth-first search over the extensions in the direction
     *        of the step computes the rows of the matrix of the part (see branch_and_bound); the carried descendants
     *        (see cluster_t::choose_center) are the first rows
     * @tparam matrix_word_t word type of the edit-distance matrix
     * @param node the context
     * @param p_idx the index of the step
     * @param data the data structures of the search
     * @param descendants_prev_dir the descendants carried in the direction of step p_idx - 1
     * @param descendants_other_dir the descendants carried in the other direction
     * @param dists_prev_dir the distances of descendants_prev_dir (the first column of the matrix)
     * @param dists_other_dir the distances of descendants_other_dir
     */
    template <typename matrix_word_t>
    void search_part(node_t& node, uint64_t p_idx, search_data_t<matrix_word_t>& data,
        const std::vector<node_t>& descendants_prev_dir, const std::vector<node_t>& descendants_other_dir,
        const std::vector<uint16_t>& dists_prev_dir, const std::vector<uint16_t>& dists_other_dir) const;

    /**
     * @brief adds a context at which the last step of a search ends to the result; of two contexts of the same string,
     *        the one with fewer errors is kept, on ties the one whose last extension is to the left
     * @param result the contexts found by the searches
     * @param node the context, without a pending copy (see commit_copy)
     */
    static void add_result(node_set_t& result, const node_t& node);

    /**
     * @brief continues the search from the final column of a part on one path
     *
     * After the last step, the local minima of the final column are added to the result. At an outer part or before a
     * direction switch, the next step starts at each local minimum, else at the cell chosen by
     * cluster_t::choose_center, with the cells after it as descendants.
     * @tparam matrix_word_t word type of the edit-distance matrix
     * @param cluster the final column
     * @param p_idx the index of the step
     * @param data the data structures of the search
     * @param dists_rev the distances of descendants_rev
     * @param descendants_rev the descendants carried in the direction other than that of step p_idx
     * @param remaining_descendants the descendants of step p_idx that are not yet rows of its matrix
     */
    template <typename matrix_word_t>
    void flush_cluster(cluster_t& cluster, uint64_t p_idx, search_data_t<matrix_word_t>& data,
        const std::vector<uint16_t>& dists_rev, const std::vector<node_t>& descendants_rev,
        const std::vector<node_t>& remaining_descendants) const;

    /**
     * @brief computes the row of a node in the matrix of the part (the row node.match_len()); if the row has a cell in
     *        the final column, the cell is set in the cluster, and the cluster is flushed (see flush_cluster) if the
     *        cell is terminal: no cell of the row has at most k_max errors, only deletions remain (see
     *        edit_distance_matrix::only_vertical_gaps_left), or the context cannot be extended
     * @tparam matrix_word_t word type of the edit-distance matrix
     * @param matrix the matrix of the part of step p_idx
     * @param cluster the final column
     * @param node the node
     * @param p_idx the index of the step
     * @param data the data structures of the search
     * @param dists_rev see flush_cluster
     * @param descendants_rev see flush_cluster
     * @param remaining_descendants see flush_cluster
     * @return true <=> the node is not extended (no cell of the row has at most k_max errors, or the cluster was
     *         flushed)
     */
    template <typename matrix_word_t>
    bool branch_and_bound(edit_distance_matrix<matrix_word_t>& matrix, cluster_t& cluster, const node_t& node,
        uint64_t p_idx, search_data_t<matrix_word_t>& data, const std::vector<uint16_t>& dists_rev,
        const std::vector<node_t>& descendants_rev, const std::vector<node_t>& remaining_descendants) const;

    // ############################# INTERNAL METHODS: ERROR MODELS #############################

    // returns the key of the extension parent followed by the symbol sym in _ext_trie_child
    static uint64_t ext_trie_key(uint32_t parent, uint8_t sym)
    {
        return (uint64_t(parent) << std::numeric_limits<uint8_t>::digits) | sym;
    }

    // returns the extension parent followed by the symbol sym (see _ext_trie_parent), adding it if it is new
    uint32_t ext_trie_child(uint32_t parent, uint8_t sym) const;

    /**
     * @brief returns the entry of the string of a context in a map, inserting the context with an error count and a
     *        cost if it is absent
     * @param map the map
     * @param node the context
     * @param err the error count of an inserted entry
     * @param cost the cost of an inserted entry
     * @param inserted is set to whether the entry was inserted
     * @param counters the packed counters of an inserted entry
     * @param min_counters the smallest value of each counter over the alignments of an inserted entry
     */
    model_entry_t& entry_of(string_map_t& map, const node_t& node, uint64_t err, uint16_t cost, bool& inserted,
        uint64_t counters = 0, uint64_t min_counters = 0) const;

    /**
     * @brief adds a member to an entry, or lowers the error count and the cost of its member with the same extension
     * @param entry the entry (with members)
     * @param mem the member
     */
    static void add_member(model_entry_t& entry, const model_member_t& mem);

    /**
     * @brief inserts a context with an error count and a cost into a map, or lowers the error count and the cost of
     *        the entry of its string
     * @param map the map
     * @param node the context
     * @param err the error count
     * @param cost the cost
     * @param counters the packed counters of an alignment with the cost (kept with the smallest cost)
     * @param min_counters the smallest value of each counter over the alignments (the entry keeps the smallest values)
     */
    void insert_or_lower(string_map_t& map, const node_t& node, uint64_t err, uint16_t cost, uint64_t counters = 0,
        uint64_t min_counters = 0) const;

    /**
     * @brief inserts a context with members into a map, or adds the members to those of the entry of its string
     * @param map the map
     * @param node the context
     * @param members the members (at least one)
     */
    void insert_members(string_map_t& map, const node_t& node, const std::vector<model_member_t>& members) const;

    /**
     * @brief extends each context with members in a map by the longest common prefix of the extensions of its members
     *        and removes this prefix from their extensions; a context whose only member has an empty extension then
     *        represents its own string
     * @param map the map
     * @param dir the direction of the extensions of the members
     */
    void narrow_to_members(string_map_t& map, direction_t dir) const;

    /**
     * @brief adds the occurrences at the text ends that the search with unit costs misses (see locate_at_text_ends)
     *        for a model in which all errors have the same cost: the cheapest alignment against a prefix of the first
     *        window, and the one with the fewest errors against a suffix of the last window
     * @tparam mode whether the occurrences have a CIGAR index
     * @tparam word_t word type of the edit-distance matrix
     * @param P the pattern
     * @param strings the strings of the occurrences (see locate_model); the windows are appended
     * @param beg [0..number of strings] beg[x] is the start of string x in strings, beg.back() is strings.size()
     * @param occ the occurrences; in CIGAR mode, the CIGAR index of an added one is the index of its window in strings
     */
    template <cigar_mode_t mode, typename word_t>
    void locate_at_text_ends_model(const inp_t& P, inp_t& strings, std::vector<uint64_t>& beg,
        std::vector<aprx_occ_t<mode>>& occ) const;

    /**
     * @brief locates P with the error model params.costs: for every reported occurrence (pos, len, err), P has a valid
     *        alignment against T[pos, pos + len) with cost err, and for every valid alignment that starts at a position
     *        p, an occurrence with at most its cost is reported at most 2 * max_errors() + 1 positions away; the
     *        occurrences are reported in increasing order of their positions
     *
     * If all errors have the same cost (see apm_costs_t::uniform_cost), the contexts are computed with search, which
     * drops the contexts covered by others; else with search_model, which discards contexts by their costs if
     * apm_costs_t::prunes() is true. The occurrences are then filtered as with unit costs (see
     * filter_edit_distance_occurrences); if no context was dropped, the occurrences that cover all others are kept
     * instead (see select_covering_occurrences).
     * @tparam mode whether each occurrence is reported with a CIGAR index
     * @tparam report_fnc_t type of the report function
     * @param P the pattern
     * @param scheme the search scheme (its k must be at least params.costs->max_errors())
     * @param report the report function
     * @return in CIGAR mode, the CIGARs (one per context with a reported occurrence); empty in NO_CIGAR mode
     */
    template <cigar_mode_t mode, typename report_fnc_t>
    std::vector<cigar_t<uint8_t>> locate_model(const inp_t& P, const search_scheme_t& scheme,
        report_fnc_t& report) const;

    /**
     * @brief executes the searches of an error model (see apm_model_plan); in every step, the search is continued from
     *        each cell of the final column of the part within the bounds of the step (a context reached by the search
     *        from a shorter context is continued in that search, see process_part), so the result contains every text
     *        string to which P has a valid alignment
     *
     * After a step whose part is at the end of P in its direction, the later steps only extend the other end: the cells
     * of the final column of a path become members of the context of its first cell in the final column (see
     * model_member_t), this context is extended by the longest common prefix of the extensions of its members (see
     * narrow_to_members), and the later steps continue it once for all of them; the members of the contexts of the last
     * step are then extended to their own contexts.
     * @tparam matrix_word_t word type of the edit-distance matrix
     * @param P the pattern
     * @param scheme the search scheme
     * @return the contexts of the text strings, with their matched strings
     */
    template <typename matrix_word_t>
    std::vector<node_t> search_model(const inp_t& P, const search_scheme_t& scheme) const;

    /**
     * @brief matches a part from a context with an error count and a cost: inserts every context at which the part
     *        ends with an error count in [k_min, k_max] and a cost of at most the bound of the model into target
     *
     * With prune, cells of the final column of a path are discarded by their costs:
     * 1. A cell is not inserted if an inserted cell above it on the path has a cost of at most its cost minus
     *    apm_costs_t::prune_slack and, if dir is LEFT, is at most 2 * max_errors() + 1 rows away. Every occurrence of
     *    the string of the cell then has an occurrence of the string of the inserted cell at the same position (dir is
     *    RIGHT) or at most that far to the right (dir is LEFT).
     * 2. A path is not extended further once every cell of its row outside the final column has a cost of at least
     *    the cost of such an inserted cell plus apm_costs_t::prune_slack.
     * Only an inserted cell whose counters (see apm_cost_band::final_counters) do not exceed their bounds covers other
     * cells, and with affine gaps (see apm_costs_t::is_affine), no cell of a path is discarded after a restart.
     *
     * @tparam matrix_word_t word type of the edit-distance matrix
     * @param entry the context, its error count (at most k_max) and its cost
     * @param part the part, read in the direction dir
     * @param beg the first position of the part in P
     * @param end the position after the last position of the part in P
     * @param dir the direction of the step
     * @param outer true <=> the part is at the end of P in direction dir; then a context is extended past the end of
     *        the part although only deletions remain (an alignment does not end with a deletion)
     * @param collect_as_members true <=> a cell of the final column is inserted as a member of the context of the
     *        first cell of the final column on its path (see search_model); entry and the contexts of step_contexts
     *        have no members then
     * @param prune true <=> cells are discarded by their costs in this step (see above); only if the step is the last
     *        one or collect_as_members is true (then the model has no counters), and apm_costs_t::prunes()
     * @param k_min the lower bound of the step
     * @param k_max the upper bound of the step
     * @param matrix the matrix of the part read in direction dir
     * @param band the dynamic program over the costs (see apm_cost_band), or nullptr if the costs are not computed in
     *        the search (see search_model)
     * @param step_contexts the contexts of the step, which contain entry; a context of it reached by this search at
     *        most depth symbols after entry is continued in this search: the matrices are restarted in its row with its
     *        error count and its cost, and it is marked as covered
     * @param depth see step_contexts
     * @param target the map the contexts are inserted into
     */
    template <typename matrix_word_t>
    void process_part(const model_entry_t& entry, const directional_substring<inp_t>& part, uint64_t beg, uint64_t end,
        direction_t dir, bool outer, bool collect_as_members, bool prune, uint64_t k_min, uint64_t k_max,
        edit_distance_matrix<matrix_word_t>& matrix, apm_cost_band* band, string_map_t& step_contexts, uint64_t depth,
        string_map_t& target) const;

public:
    // ############################# CONSTRUCTORS #############################

    /**
     * @brief constructs the queries on an index
     * @param index the index
     * @param params query parameters
     */
    apm_edit(const idx_t& index, apm_params_t params = {}) : base_t(index, params)
    {
        // in the edit-distance matrices, the symbols are ranked with the map of the bytes (see move_r::map_int)
        if (!index.forward_index().byte_alphabet)
            throw std::invalid_argument("apm_edit: the index does not have a byte alphabet (see move_r::byte_alphabet)");
    }

    // ############################# QUERY METHODS #############################

    /**
     * @brief locates P with at most k errors (with an error model, see locate_model); each position is reported at most
     *        once
     *
     * Every reported (pos, len, err) has ed(P, T[pos, pos + len)) = err <= k, with an alignment that does not start or
     * end with a deletion (as with an error model). For every p' and L with ed(P, T[p', p' + L)) = e <= k, an occurrence
     * with at most e errors is reported at most 2k + 1 positions away from p'.
     * @tparam mode whether each occurrence is reported with a CIGAR index
     * @tparam report_fnc_t type of the report function (called with {pos, len, err} or, in CIGAR mode,
     *         {pos, len, err, cigar index})
     * @param P the pattern; throws std::invalid_argument if P.size() + k > aprx_occ_max_len
     * @param scheme the search scheme (provides k)
     * @param report the report function
     * @return in CIGAR mode, the CIGARs (one per suffix array interval); empty in NO_CIGAR mode
     */
    template <cigar_mode_t mode = NO_CIGAR, typename report_fnc_t>
    std::vector<cigar_t<uint8_t>> locate(const inp_t& P, const search_scheme_t& scheme, report_fnc_t report) const;

    /**
     * @brief returns the contexts found by the search scheme (without the occurrences at the text ends), without
     *        duplicates; each stores its matched string
     * @param P the pattern
     * @param scheme the search scheme
     */
    node_set_t search(const inp_t& P, const search_scheme_t& scheme) const
    {
        return with_edit_distance_word(scheme.k, [&]<typename word_t>() { return search<word_t>(P, scheme); });
    }
};

#include <algorithms/move_rb/apm/edit.tpp>
#include <algorithms/move_rb/apm/edit_model.tpp>

#endif // MOVE_R_ALGORITHMS_MOVE_RB_APM_EDIT_HPP
