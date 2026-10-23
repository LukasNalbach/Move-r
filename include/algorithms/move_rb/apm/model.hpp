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

#ifndef MOVE_R_ALGORITHMS_MOVE_RB_APM_MODEL_HPP
#define MOVE_R_ALGORITHMS_MOVE_RB_APM_MODEL_HPP

#include <algorithm>
#include <cstdint>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <misc/apm.hpp>
#include <misc/cigar.hpp>
#include <misc/fasta.hpp>

/**
 * @brief the operations of an alignment of a pattern P against a text window, as bits of a mask (the names follow
 *        the CIGAR operations)
 */
enum apm_op_t : uint8_t {
    APM_SUB = 1,  // a symbol of P aligned to a different text symbol that it does not match (X)
    APM_INS = 2,  // a symbol of P that is absent from the text (I)
    APM_DEL = 4,  // a text symbol that is absent from P (D)
    APM_OPEN = 8  // the first operation of a run of insertions or of deletions (a gap)
};

// the cost of an operation that is not allowed
static constexpr uint16_t apm_forbidden = std::numeric_limits<uint16_t>::max();

// the counters of a model are packed into one word, one byte lane each
static constexpr uint64_t apm_max_counters = sizeof(uint64_t);

// the largest bound of a counter; a counter one above its bound sets the highest bit of its lane
static constexpr uint64_t apm_max_counter_bound = std::numeric_limits<int8_t>::max();

/**
 * @brief returns a + b, or apm_forbidden if it exceeds a bound below apm_forbidden (then a forbidden summand exceeds
 *        it as well)
 */
inline uint16_t apm_bounded_sum(uint32_t a, uint32_t b, uint64_t bound)
{ return a + b > bound ? apm_forbidden : uint16_t(a + b); }

/**
 * @brief a range [beg, end) of positions of P; a negative value counts from the end of P (-1 is the last position),
 *        as in a Python slice
 */
struct apm_range_t {
    int64_t beg = 0;
    int64_t end = std::numeric_limits<int64_t>::max();

    /**
     * @brief returns the absolute range [beg, end) for a pattern of length m, clamped to [0, m]
     * @param m the pattern length
     */
    std::pair<uint64_t, uint64_t> absolute(uint64_t m) const;
};

/**
 * @brief sets the cost of the operations ops in a range of P; ops is a mask of apm_op_t
 *
 * A substitution of P[i] by the text symbol c is affected if i lies in range, P[i] is in from (or from is empty)
 * and c is in to (or to is empty); an insertion of P[i] if i lies in range and P[i] is in from (or from is empty);
 * a deletion between P[b - 1] and P[b] if b lies strictly inside range (see apm_costs_t::in_range); a gap opened by
 * an insertion or a deletion like that insertion or deletion.
 */
struct apm_cost_rule_t {
    uint8_t ops = 0;
    uint16_t cost = 1; // apm_forbidden: the operations are not allowed
    apm_range_t range;
    std::string from; // symbols of P (substitutions and insertions)
    std::string to;   // text symbols (substitutions)
};

/** @brief bounds the number of operations ops (a mask of apm_op_t) in a range of P (see apm_cost_rule_t) */
struct apm_counter_t {
    uint8_t ops = 0;
    uint8_t max = 0; // at most apm_max_counter_bound
    apm_range_t range;
};

/** @brief a symbol of P in from that lies in range also matches the text symbols in to */
struct apm_class_t {
    std::string from;
    std::string to;
    apm_range_t range;
};

/** @brief whether the edit-distance search bounds the cost of the alignments (see apm_cost_band) */
enum apm_cost_dp_t : uint8_t {
    // PRUNE if the search can prune (see apm_costs_t::can_prune) and neither gap openings nor counters apply, else ON
    // if the costs allow at most one insertion or deletion but the search more errors, else OFF
    APM_COST_DP_AUTO = 0,
    APM_COST_DP_ON = 1,   // the search bounds the cost
    APM_COST_DP_OFF = 2,  // the search does not compute the cost
    APM_COST_DP_PRUNE = 3  // ON, and the search prunes its contexts by their costs if it can
};

/**
 * @brief an error model for approximate pattern matching: an alignment of P against a text window costs the sum of
 *        the costs of its operations and is valid if it costs at most the bound and no counter exceeds its bound
 *
 * By default, every substitution, insertion and deletion costs 1, opening a gap costs 0, there are no counters and a
 * symbol only matches itself; the bound is the maximum number of errors k of the query. An alignment does not start
 * or end with a deletion.
 */
struct apm_model_t {
    std::vector<apm_cost_rule_t> costs; // applied in this order (a later rule overrides an earlier one)
    std::vector<apm_counter_t> counters; // at most apm_max_counters
    std::vector<apm_class_t> classes;
    bool iupac = false;   // an IUPAC code in P matches each of its bases (A, C, G, T)
    bool quality = false; // a substitution of P[i] costs the Phred quality of P[i] (before the cost rules)
    uint16_t tau = 0;     // the smallest cost of a substitution that the search counts as an error (see
                          // apm_costs_t::tau); 0: chosen for each pattern (see apm_costs_t::compile)
    apm_cost_dp_t cost_dp = APM_COST_DP_AUTO; // whether the edit-distance search bounds the cost

    // returns whether the model is the default model
    bool is_default() const
    {
        return costs.empty() && counters.empty() && classes.empty() && !iupac && !quality && tau == 0 &&
            cost_dp == APM_COST_DP_AUTO;
    }
};

/**
 * @brief parses a non-negative integer or "inf" (apm_forbidden)
 * @param str the value
 * @param what the name of the value in the error message
 */
uint64_t parse_apm_value(std::string_view str, const std::string& what);

/**
 * @brief parses a range "<beg>:<end>", "<beg>:", ":<end>" or "<pos>" (see apm_range_t)
 * @param str the range
 */
apm_range_t parse_apm_range(std::string_view str);

/**
 * @brief parses a cost rule "<ops>=<cost>[@<range>][/<from>,<to>]", where ops is a subset of the letters s, i, d, o
 *        (substitution, insertion, deletion, gap opening) and cost is a non-negative integer or inf
 * @param str the rule
 */
apm_cost_rule_t parse_apm_cost_rule(std::string_view str);

/**
 * @brief parses a counter "<ops>=<max>[@<range>]" (ops as in parse_apm_cost_rule)
 * @param str the counter
 */
apm_counter_t parse_apm_counter(std::string_view str);

/**
 * @brief parses a symbol class "<from>,<to>[@<range>]"
 * @param str the class
 */
apm_class_t parse_apm_class(std::string_view str);

/**
 * @brief an error model compiled for one pattern P and the alphabet of an index: the cost of every operation at every
 *        position of P, the counters each operation increments, and the bounds the search derives from them
 *
 * Text symbols are addressed by their rank in the effective alphabet of the index (alph_map[c] - 1). A deletion lies
 * at a boundary b in [1, m - 1], between P[b - 1] and P[b]. The search counts a substitution as an error if it costs
 * at least tau (see tau()), and every insertion and deletion; the other substitutions and the matches by a class are
 * free. The counters of an alignment are packed into one word with one byte lane per counter (see
 * apm_max_counters).
 */
class apm_costs_t {
protected:
    // ############################# TYPES AND CONSTANTS #############################

    static constexpr uint64_t lane_bits = std::numeric_limits<uint8_t>::digits; // the width of a counter lane
    // the highest bit of each lane
    static constexpr uint64_t lane_high_bits = (~uint64_t(0) / std::numeric_limits<uint8_t>::max()) << (lane_bits - 1);

    // a label of the alignment DP with counters: the cost and the counters of an alignment of P[0, i) against M[0, j)
    struct label_t {
        uint16_t cost;
        uint8_t step;     // the operation into this label (see align_reusing)
        uint8_t state;    // the last operation: 0 (match or substitution), 1 (insertion) or 2 (deletion)
        uint32_t prev;    // index of the preceding label, or no_label
        uint64_t counters; // packed counters
    };

    static constexpr uint32_t no_label = std::numeric_limits<uint32_t>::max();

    // ############################# DATA MEMBERS #############################

    uint64_t _m = 0; // m, the length of P
    uint16_t _bound = 0; // the maximum cost of an alignment
    uint16_t _tau = 1; // see tau()
    distance_metric_t _metric = EDIT_DISTANCE;
    uint64_t _num_ranks = 0; // the size of the effective alphabet of the index
    std::vector<uint8_t> _pattern; // [0..m-1] P
    const std::vector<uint8_t>* _alph_map = nullptr; // the alphabet map of the last call of compile
    // the ranks of the bytes: [0..255] the zero-based rank of byte c in the effective alphabet, or _num_ranks
    std::vector<uint8_t> _rank;
    // the bytes of the ranks: [0.._num_ranks-1] the byte of rank r
    std::vector<uint8_t> _unrank;
    // the matches: [0..m * _num_ranks - 1] _matches[i * _num_ranks + r] = 1 <=> P[i] matches the text symbol of rank r
    std::vector<uint8_t> _matches;
    // the substitutions: [0..m * _num_ranks - 1] _sub[i * _num_ranks + r] = cost of aligning P[i] to rank r (0 if it
    // matches)
    std::vector<uint16_t> _sub;
    /* the costs of the insertions, the deletions and of opening their gaps: [0..m-1] _ins[i] and _ins_open[i] of
       P[i]; [0..m] _del[b] and _del_open[b] of boundary b (_del[0] = _del[m] = apm_forbidden) */
    std::vector<uint16_t> _ins;
    std::vector<uint16_t> _del;
    std::vector<uint16_t> _ins_open;
    std::vector<uint16_t> _del_open;
    /* the packed counter increments of the operations, indexed as their costs: [0..m-1] those of a substitution of
       P[i], of an insertion of P[i] and of opening a gap of insertions at P[i]; [0..m] those of a deletion and of
       opening a gap of deletions at boundary b */
    std::vector<uint64_t> _sub_counts;
    std::vector<uint64_t> _ins_counts;
    std::vector<uint64_t> _del_counts;
    std::vector<uint64_t> _ins_open_counts;
    std::vector<uint64_t> _del_open_counts;
    // the free bytes of all positions: [0.._free_beg[m]-1]; those of P[i] are _free_syms[_free_beg[i], _free_beg[i + 1]),
    // P[i] first
    std::vector<uint8_t> _free_syms;
    std::vector<uint32_t> _free_beg; // [0..m] see _free_syms
    // [0..m-1] the cost of the cheapest error at P[i] (a counted substitution or an insertion) within the bound, or
    // apm_forbidden
    std::vector<uint16_t> _min_error_cost;
    // [0..m-1] 1 <=> every error at P[i] within the bound increments a counter
    std::vector<uint8_t> _all_errors_counted;
    uint64_t _num_counters = 0;
    /* the counters: [0.._num_counters-1] counter j counts the operations _counter_ops[j] (a mask of apm_op_t) in the
       range [_counter_beg[j], _counter_end[j]) of P, at most _counter_max[j] times */
    uint8_t _counter_ops[apm_max_counters] = {};
    uint64_t _counter_beg[apm_max_counters] = {};
    uint64_t _counter_end[apm_max_counters] = {};
    uint64_t _counter_max[apm_max_counters] = {};
    // lane j: apm_max_counter_bound - _counter_max[j], so a counter above its bound sets the highest bit of its lane
    uint64_t _counter_bias = 0;
    bool _affine = false; // see is_affine
    bool _bounds_costs = false; // see bounds_costs
    bool _prunes = false; // see prunes
    uint16_t _prune_slack = 0; // see prune_slack
    uint16_t _uniform_cost = 0; // see uniform_cost
    uint64_t _max_errors = 0; // max_errors(0, m)
    uint64_t _max_indels = 0; // max_num_errors(0, m, true)

    // ############################# INTERNAL METHODS #############################

    /**
     * @brief returns an upper bound on the number of errors (or, if indels_only, of insertions and deletions) of a
     *        valid alignment in the range [beg, end) of P, counting the deletions at the boundaries strictly inside it;
     *        with the counters, unless counters is false
     */
    uint64_t max_num_errors(uint64_t beg, uint64_t end, bool indels_only, bool counters = true) const;

    // sets tau and computes the free symbols, the cheapest errors and the bounds of the search that depend on it
    void set_tau(uint16_t tau);

    // sets whether the edit-distance search bounds the cost of the alignments
    void choose_cost_dp(apm_cost_dp_t cost_dp);

    // returns the cost of every error if all substitutions, insertions and deletions cost the same and neither
    // counters nor gap openings apply, else 0
    uint16_t compute_uniform_cost() const;

    /**
     * @brief align_prefix, reusing the cells of the dynamic program in the columns [0..reuse] of the previous string
     *
     * The cells are kept in thread_local storage: with reuse > 0, the previous call of align_reusing on this thread
     * must have been for this object and a string whose first reuse symbols are those of M (as align_prefixes calls
     * it).
     * @tparam str_t contiguous range type of M (one byte per symbol)
     * @param M the text window
     * @param cigar if not nullptr, the CIGAR of the alignment is written to it
     * @param reuse the number of first symbols that M shares with the previous string (0: nothing is reused)
     */
    template <std::ranges::contiguous_range str_t>
    std::pair<uint64_t, uint64_t> align_reusing(const str_t& M, cigar_t<uint8_t>* cigar, uint64_t reuse) const;

    /**
     * @brief align_reusing without the counters: one cost per cell and state of the DP
     * @tparam str_t contiguous range type of M (one byte per symbol)
     * @param M the text window
     * @param cigar if not nullptr, the CIGAR of the alignment is written to it
     * @param counters if not nullptr, the packed counters of the alignment are written to it
     * @param reuse see align_reusing
     */
    template <std::ranges::contiguous_range str_t>
    std::pair<uint64_t, uint64_t> align_without_counters(const str_t& M, cigar_t<uint8_t>* cigar, uint64_t* counters,
        uint64_t reuse) const;

public:
    // ############################# CONSTRUCTORS #############################

    /**
     * @brief compiles a model for a pattern
     * @tparam pat_t contiguous range type of the pattern (one byte per symbol)
     * @param model the error model
     * @param P the pattern as it is searched
     * @param quality the Phred+33 qualities of the read P was taken from (if model.quality)
     * @param reverse_strand true <=> P is the reverse complement of the read the ranges, the qualities and the
     *        symbols of the rules refer to
     * @param bound the maximum cost of an alignment (at most the largest error count of an occurrence)
     * @param metric HAMMING_DISTANCE (only substitutions) or EDIT_DISTANCE
     * @param alph_map [0..255] maps a byte to its one-based rank in the effective alphabet of the index (0: absent)
     */
    template <std::ranges::contiguous_range pat_t>
    apm_costs_t(const apm_model_t& model, const pat_t& P, std::string_view quality, bool reverse_strand,
        uint16_t bound, distance_metric_t metric, const std::vector<uint8_t>& alph_map)
    { compile(model, P, quality, reverse_strand, bound, metric, alph_map); }

    apm_costs_t() = default;

    // ############################# CONSTRUCTION METHODS #############################

    /**
     * @brief compiles a model for a pattern, reusing the memory of this object (see the constructor); the ranks are
     *        recomputed only if alph_map is another vector than in the previous call. If model.tau is 0, tau is the
     *        candidate (1, or the cost of a substitution plus 1) that maximizes the number of positions without free
     *        symbols other than their own per part of a search with max_errors() + 1 parts.
     */
    template <std::ranges::contiguous_range pat_t>
    void compile(const apm_model_t& model, const pat_t& P, std::string_view quality, bool reverse_strand,
        uint16_t bound, distance_metric_t metric, const std::vector<uint8_t>& alph_map);

    // ############################# ACCESS METHODS #############################

    // returns the pattern length
    uint64_t length() const { return _m; }

    // returns P[i]
    uint8_t symbol(uint64_t i) const { return _pattern[i]; }

    // returns the maximum cost of an alignment
    uint16_t bound() const { return _bound; }

    // returns the distance metric the model was compiled for
    distance_metric_t metric() const { return _metric; }

    // returns tau, the smallest cost of a substitution that the search counts as an error
    uint16_t tau() const { return _tau; }

    // returns the zero-based rank of a byte in the effective alphabet of the index (or num_ranks() if it is absent)
    uint64_t rank(uint8_t c) const { return _rank[c]; }

    // returns the size of the effective alphabet of the index
    uint64_t num_ranks() const { return _num_ranks; }

    // returns whether P[i] matches the text symbol of rank r
    bool matches(uint64_t i, uint64_t r) const { return r < _num_ranks && _matches[i * _num_ranks + r]; }

    // returns the cost of substituting P[i] by the text symbol of rank r (not a match)
    uint16_t sub(uint64_t i, uint64_t r) const { return r < _num_ranks ? _sub[i * _num_ranks + r] : apm_forbidden; }

    // returns whether P[i] aligned to the text symbol of rank r is free in the search (a match or a substitution
    // cheaper than tau)
    bool is_free(uint64_t i, uint64_t r) const { return matches(i, r) || sub(i, r) < _tau; }

    // returns the bytes that are free at P[i] in the search, P[i] first
    std::span<const uint8_t> free_symbols(uint64_t i) const
    { return {_free_syms.data() + _free_beg[i], _free_syms.data() + _free_beg[i + 1]}; }

    // returns whether a byte other than P[i] is free at P[i]
    bool has_free(uint64_t i) const { return _free_beg[i + 1] - _free_beg[i] > 1; }

    // returns the cost of inserting P[i]
    uint16_t ins(uint64_t i) const { return _ins[i]; }

    // returns the cost of a deletion at boundary b
    uint16_t del(uint64_t b) const { return _del[b]; }

    // returns the cost of opening a gap of insertions at P[i]
    uint16_t ins_open(uint64_t i) const { return _ins_open[i]; }

    // returns the cost of opening a gap of deletions at boundary b
    uint16_t del_open(uint64_t b) const { return _del_open[b]; }

    // returns the packed counter increments of a substitution of P[i]
    uint64_t sub_count(uint64_t i) const { return _sub_counts[i]; }

    // returns the packed counter increments of an insertion of P[i]
    uint64_t ins_count(uint64_t i) const { return _ins_counts[i]; }

    // returns the packed counter increments of a deletion at boundary b
    uint64_t del_count(uint64_t b) const { return _del_counts[b]; }

    // returns the packed counter increments of opening a gap of insertions at P[i]
    uint64_t ins_open_count(uint64_t i) const { return _ins_open_counts[i]; }

    // returns the packed counter increments of opening a gap of deletions at boundary b
    uint64_t del_open_count(uint64_t b) const { return _del_open_counts[b]; }

    // returns the number of counters
    uint64_t num_counters() const { return _num_counters; }

    // returns the bound of counter j
    uint64_t counter_bound(uint64_t j) const { return _counter_max[j]; }

    // returns whether packed counters exceed the bound of a counter
    bool exceeds(uint64_t counters) const { return ((counters + _counter_bias) & lane_high_bits) != 0; }

    // returns whether opening a gap has a cost or is counted (then insertions and deletions depend on their predecessor)
    bool is_affine() const { return _affine; }

    // returns whether the edit-distance search bounds the cost of the alignments (see apm_cost_band and apm_cost_dp_t)
    bool bounds_costs() const { return _bounds_costs; }

    /**
     * @brief returns the cost of every error if every substitution, insertion and deletion costs the same (at least
     *        tau), no counter applies and opening a gap is free, else 0; then the cost of an alignment is this cost
     *        times its number of errors, where a match by a class is no error
     */
    uint16_t uniform_cost() const { return _uniform_cost; }

    /**
     * @brief returns whether the edit-distance search can prune its contexts by their costs: the errors do not all
     *        cost the same (else it prunes them by their errors), and gap openings and counters do not both apply
     */
    bool can_prune() const { return _uniform_cost == 0 && !(_affine && _num_counters > 0); }

    /**
     * @brief returns the largest cost by which apm_cost_band can underestimate the cost of an alignment of a part
     *        compared to the cost of the whole alignment: the openings of the gaps at the two ends of the part
     *        (2 * the largest cost of opening a gap of insertions + the largest cost of opening a gap of deletions)
     */
    uint16_t prune_slack() const { return _prune_slack; }

    // returns whether the edit-distance search prunes its contexts by their costs (see apm_cost_dp_t)
    bool prunes() const { return _prunes; }

    /**
     * @brief returns an upper bound on the number of errors in the search (see the class description) of a valid
     *        alignment in the range [beg, end) of P, counting the deletions at the boundaries strictly inside it
     * @param beg first position of the range
     * @param end end of the range (exclusive)
     */
    uint64_t max_errors(uint64_t beg, uint64_t end) const { return max_num_errors(beg, end, false); }

    // returns an upper bound on the number of errors in the search of a valid alignment of P
    uint64_t max_errors() const { return _max_errors; }

    // returns an upper bound on the number of insertions and deletions of a valid alignment of P
    uint64_t max_indels() const { return _max_indels; }

    // returns whether P[i] admits no error and only matches P[i] itself
    bool is_literal(uint64_t i) const { return _min_error_cost[i] == apm_forbidden && !has_free(i); }

    /**
     * @brief returns whether boundary b (between P[b - 1] and P[b]) lies in the range [beg, end) of positions: it
     *        does if beg < b < end
     */
    static bool in_range(uint64_t b, uint64_t beg, uint64_t end) { return beg < b && b < end; }

    // ############################# QUERY METHODS #############################

    /**
     * @brief computes the cheapest valid alignment of P against a prefix M[0, len) of M (len >= 1), with the shortest
     *        len among the cheapest ones
     * @tparam str_t contiguous range type of M (one byte per symbol)
     * @param M the text window
     * @param cigar if not nullptr, the CIGAR of the alignment is written to it
     * @return (cost, len), or (bound() + 1, 0) if there is no valid alignment
     */
    template <std::ranges::contiguous_range str_t>
    std::pair<uint64_t, uint64_t> align_prefix(const str_t& M, cigar_t<uint8_t>* cigar = nullptr) const
    { return align_reusing(M, cigar, 0); }

    /**
     * @brief computes align_prefix for several strings; they are aligned in lexicographic order, and each reuses the
     *        cells of the dynamic program of the preceding one that lie in their longest common prefix
     * @tparam str_t contiguous range type of the strings (one byte per symbol)
     * @param strings the concatenated strings: strings[beg[x], beg[x + 1]) is the x-th one
     * @param beg the start of each string in strings, followed by the size of strings
     * @param out (output) out[x] is align_prefix of the x-th string
     */
    template <std::ranges::contiguous_range str_t>
    void align_prefixes(const str_t& strings, const std::vector<uint64_t>& beg,
        std::vector<std::pair<uint64_t, uint64_t>>& out) const;
};

/**
 * @brief a banded dynamic program over the cost of the alignments of a part P[beg, end) of a pattern, read in a
 *        direction, against a text that grows by one symbol per row; it omits the counters and the costs of opening
 *        the gaps that touch an end of the part, so each cell is at most the cost of the model of every alignment
 *        that ends in it
 *
 * Opening a gap of insertions is paid at its first position: when it is entered (read to the right) or left (read to
 * the left). Opening a gap of deletions is paid at its boundary. A row is a text symbol, a column c a number of
 * symbols of the part (in apm_costs_t, the rows are the positions of P).
 */
class apm_cost_band {
protected:
    // ############################# DATA MEMBERS #############################

    const apm_costs_t* _costs = nullptr;
    uint64_t _bound = 0;
    uint64_t _beg = 0, _len = 0; // the part is P[_beg, _beg + _len)
    uint64_t _num_ranks = 0; // costs.num_ranks() + 1 (the last rank is a symbol that does not occur in the index)
    // row i stores the columns [i - _cols_before, i + _cols_after]; _band = _cols_before + _cols_after + 1
    uint64_t _cols_before = 0, _cols_after = 0, _band = 0;
    uint64_t _states = 1; // 1, or 3 with gaps (0: match or substitution, 1: insertion, 2: deletion)
    bool _reads_right = true; // whether the part is read to the right
    // the cells: [0..num_rows * _band * _states - 1] _rows[(i * _band + c + _cols_before - i) * _states + q] is the
    // cost of cell (i, c) in state q
    std::vector<uint16_t> _rows;
    bool _counts = false; // whether the counters of the cells are kept (see final_counters)
    // the counters of the cells: [0..num_rows * _band - 1] _counters[i * _band + c + _cols_before - i] of cell (i, c)
    std::vector<uint64_t> _counters;
    /* the counters of the operations on the part, by the number c of symbols of the part in reading direction:
       [0.._len * _num_ranks - 1] _sub_count_vals[c * _num_ranks + r] of aligning the c-th symbol to rank r;
       [0.._len-1] _ins_count_vals[c] of inserting the c-th symbol; [0.._len] _del_count_vals[c] of a deletion after c
       symbols */
    std::vector<uint64_t> _sub_count_vals;
    std::vector<uint64_t> _ins_count_vals;
    std::vector<uint64_t> _del_count_vals;
    // the costs of the operations on the part, indexed as the counters
    std::vector<uint16_t> _sub_vals;
    std::vector<uint16_t> _ins_vals;
    std::vector<uint16_t> _del_vals;
    std::vector<uint16_t> _ins_open_in;  // [0.._len-1] paid when a gap of insertions is entered at the c-th symbol
    std::vector<uint16_t> _ins_open_out; // [0.._len] paid when a gap of insertions is left after c symbols
    std::vector<uint16_t> _del_open;     // [0.._len] paid when a gap of deletions is opened after c symbols

    // ############################# INTERNAL METHODS #############################

    // returns the position of P of the c-th symbol of the part in reading direction
    uint64_t pos(uint64_t c) const { return _reads_right ? _beg + c : _beg + _len - 1 - c; }

    // returns the boundary of P after c symbols of the part in reading direction
    uint64_t boundary(uint64_t c) const { return _reads_right ? _beg + c : _beg + _len - c; }

    // returns a + b + c, or apm_forbidden if it exceeds the bound
    uint16_t add_bounded(uint32_t a, uint32_t b, uint32_t c = 0) const { return apm_bounded_sum(a, b + c, _bound); }

public:
    // ############################# CONSTRUCTION METHODS #############################

    /**
     * @brief starts the dynamic program
     * @param costs the compiled model
     * @param beg first position of the part
     * @param end end of the part (exclusive)
     * @param dir the direction the part is read in
     * @param cost_in the cost before the part
     * @param counters_in the packed counters of an alignment before the part with cost cost_in (see final_counters)
     * @param num_rows the number of rows (text symbols + 1)
     * @param restart_rows the rows [1..restart_rows] can be restarted (see restart)
     */
    void init(const apm_costs_t& costs, uint64_t beg, uint64_t end, direction_t dir, uint16_t cost_in,
        uint64_t counters_in, uint64_t num_rows, uint64_t restart_rows = 0);

    /**
     * @brief computes row i from row i - 1 for the next text symbol
     * @param i [1..num_rows-1] the row
     * @param sym the text symbol
     * @return whether some cell of row i is at most the bound of the model
     */
    bool compute_row(uint64_t i, uint8_t sym);

    /**
     * @brief adds an alignment that starts in cell (i, 0) with a cost, i.e., before the part, after i text symbols
     * @param i [1..restart_rows] the row (computed)
     * @param cost the cost of the alignment in cell (i, 0)
     * @param counters the packed counters of the alignment
     * @return whether some cell of row i is at most the bound of the model
     */
    bool restart(uint64_t i, uint16_t cost, uint64_t counters);

    // ############################# QUERY METHODS #############################

    // returns the cost of the cell of row i in which the whole part is aligned, or apm_forbidden
    uint16_t final_cost(uint64_t i) const;

    /**
     * @brief returns the packed counters of an alignment with the cost of final_cost(i), if the model has counters
     *        and no gap opening (else 0): the dynamic program keeps, for each cell, the counters of one of its
     *        cheapest alignments, from counters_in and the counters of a restart on
     */
    uint64_t final_counters(uint64_t i) const;

    // returns the smallest cost of the cells of row i in which the part is not aligned completely, or apm_forbidden
    uint16_t min_partial_cost(uint64_t i) const;
};

#include "model.tpp"

#endif // MOVE_R_ALGORITHMS_MOVE_RB_APM_MODEL_HPP
