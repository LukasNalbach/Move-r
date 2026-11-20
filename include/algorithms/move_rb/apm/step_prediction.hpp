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
 * @file step_prediction.hpp
 * @brief predicts how many steps the searches of a search scheme take
 */

#ifndef MOVE_R_ALGORITHMS_MOVE_RB_APM_STEP_PREDICTION_HPP
#define MOVE_R_ALGORITHMS_MOVE_RB_APM_STEP_PREDICTION_HPP

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <gtl/phmap.hpp>
#include <misc/apm.hpp>
#include <misc/search_schemes.hpp>
#include <algorithms/move_rb/apm/model.hpp>
#include <algorithms/move_rb/apm/search_base.hpp>

/**
 * @brief returns the share of the occurrences of a string of each length L in [0..max_len] with distinct continuations
 *        (see apm_step_prediction): d(L) / n with the number d(L) of the distinct substrings of length L of the text
 *        (see move_rb::distinct_substrings), or r / n if the index does not store them
 * @param index the index
 * @param max_len the largest length
 */
template <typename idx_t>
std::vector<double> apm_distinct_shares(const idx_t& index, uint64_t max_len);

/**
 * @brief predicts how many steps the searches of a search scheme take for a partition of a pattern P into parts
 *
 * The number of steps of a search is the number of nodes of its search tree, as in the analysis of search schemes in
 * a random text (Kucherov, Salikhov, Tsur; Kianfar et al.). A search starts with the occurrences of its first part and
 * extends them by the strings whose numbers of errors lie within the error bounds of its steps. A(t) denotes the number
 * of these strings at depth t, i.e., after t further positions of P. At position i of P, a string continues in
 * match(i) ways without an error and in error(i) ways with an error (a substitution, an insertion or a deletion), where
 * a text symbol c counts as sigma * p_c ways, its weight (p_c is the probability of c in the text, sigma the effective
 * alphabet size, see the constructor).
 *
 * Let d(L) be the number of distinct substrings of length L of the text (without a separator of the sequences, see
 * move_rb::distinct_substrings). Of the occ occurrences of the first part, with length l, x(L) = occ * d(L) / n are
 * assumed to have distinct continuations to the length L (with r / n instead of d(L) / n if the index does not store
 * d). A string at depth t occurs among these continuations with a probability of at most min(1, x(l + t) / sigma^t),
 * and the predicted number of steps of the search is the sum of these probabilities over the strings at the depths
 * t >= 1. A search that matches its first part symbol by symbol (if the part has a free symbol or can have errors, or
 * if not apm_params_t::start_at_part_contexts) starts with the n occurrences of the empty string (l = 0, x(L) = d(L)),
 * and the positions of its first part are counted in A(t). With an error model, the prediction is made for the
 * searches that apm_model_plan runs (see choose); for the Hamming distance, a string whose counters exceed their
 * bounds is not counted.
 *
 * For the edit distance (see for_each_counted_step), the depth t of a string is its number of text symbols after the
 * start, and the positions of P are taken in the order of the search. The search computes the row of the
 * edit-distance matrix of each extension of a string whose row has a cell within the bounds, so the steps at depth t
 * are the sigma * A(t - 1) children of the strings at depth t - 1 that the search extends. A string continues at the
 * next position i of P with a match or a substitution of P[i], with an insertion (P[i] without a text symbol) or with a
 * deletion (a text symbol without a position of P). Of two alignments of the same string that end at the same position
 * of P, the one whose insertions and deletions lie further in the direction of the search is counted, if it has at most
 * as many errors:
 * - a deletion before P[i] is of a symbol c that does not match P[i], and it is followed by a deletion or by a match of
 *   P[i] (a deletion of c followed by a substitution of P[i] by c' is replaced by a substitution of P[i] by c followed
 *   by a deletion of c'),
 * - an insertion of P[i] is followed by an insertion or by a match of P[i'] (the next position) with a symbol that does
 *   not match P[i],
 * - an insertion does not follow a deletion and a deletion does not follow an insertion (a substitution), and an
 *   insertion (deletion) does not follow a deletion (insertion) and one match.
 * A string is counted at its alignments that end with a match or a substitution, except those that end with an
 * insertion followed by a match (the string also ends at the position before) and those that end with a deletion
 * followed by a match of a symbol that also matches the next position (the string also ends at the next position). A
 * string with alignments that end at different positions of P is counted once per end, and A(t) is at most sigma^t.
 * After the last position of a step, if the next step has the other direction or there is none, the search also
 * extends the strings by deletions of any symbol at a part at the first or the last position of P and if it computes
 * the costs (see has_band). If it computes the costs, the operations are those of the error model, and as in its
 * dynamic program, the strings that end at the same position of P are dropped if the lane-wise minimum of their
 * counters exceeds a bound, else they continue with this minimum (the search checks the bounds of the steps on the
 * edit distance and the cost bound on the error model; the prediction applies the bounds of the steps to the
 * operations of the error model); otherwise the search uses the edit distance (any symbol can be substituted, inserted
 * and deleted) and does not check the counters.
 */
class apm_step_prediction {
public:
    // the searches that apm_model_plan runs for a partition (see choose)
    struct choice_t {
        double steps;  // the predicted number of steps of the searches
        bool lowered;  // whether the upper bounds of the searches are lowered (see apm_model_plan)
        bool anchored; // whether the anchored search is run instead of the searches of the scheme
    };

private:
    // the error bounds of the searches
    struct bounds_mode_t {
        bool lowered;         // whether the upper bounds are lowered (see upper_bound)
        bool no_lower_bounds; // whether all lower bounds are 0
    };

    const search_scheme_t& _scheme;
    const apm_costs_t* _costs; // the error model, or nullptr
    double _sigma;             // the effective alphabet size (see the constructor)
    bool _edit;                // whether the searches use the edit distance
    /* the share of the occurrences of a string of length L with distinct continuations (see the class description):
       [0..|_distinct_share| - 1] d(L) / n, or r / n (see apm_distinct_shares); the last value holds for longer strings */
    std::vector<double> _distinct_share;
    uint64_t _n;               // n, the text length
    /* the anchor P[_anchor_beg, _anchor_end) of the error model (empty if the anchored search is not considered) and
       its number of occurrences */
    uint64_t _anchor_beg = 0, _anchor_end = 0, _anchor_occ = 0;
    bool _scheme_lower_bounds = false; // whether a step of the scheme has a lower bound above 0
    bool _may_lower;     // see apm_params_t::lower_upper_bounds
    bool _from_contexts; // see apm_params_t::start_at_part_contexts

    // [0..m-1] 1 <=> a symbol other than P[i] is free at position i (see apm_costs_t::has_free)
    std::vector<uint8_t> _free;

    /* the numbers of ways in which a string continues at position i of P (see the class description): [0..m-1]
       _match[i] without an error, _sub[i] with a substitution, _ins[i] with an insertion, _del[i] with a deletion and
       _error[i] = _sub[i] + _ins[i] + _del[i] with an error */
    std::vector<double> _match, _sub, _ins, _del, _error;

    // [0..m] 1 <=> a deletion is allowed at boundary b of P (before P[b])
    std::vector<uint8_t> _del_allowed;

    // [0..m-2] _common_next[i] = the number of ways of the symbols that match both P[i] and P[i + 1]
    std::vector<double> _common_next;

    // the last operations of a group of strings (see the class description and for_each_counted_step)
    enum last_op_t : uint8_t {
        DIAG,       // a match or a substitution
        INS,        // an insertion
        DEL,        // a deletion
        DEL_MATCH,  // a deletion followed by a match
        INS_MATCH,  // an insertion followed by a match
        TRAIL       // a deletion after the last position of a step (see search_order_t::trail)
    };

    /* a group of the strings at depth t that have searched the same number of positions of P and have the same number
       of errors, the same counters and the same last operations */
    struct counted_t {
        uint64_t counters;   // the packed counters (see apm_costs_t)
        uint32_t u;          // the number of positions of P searched after the start (see search_order_t)
        uint8_t errors;      // the number of errors
        last_op_t last;      // the last operations
        double num;          // the number of these strings divided by sigma^t
        double counted;      // the part of num counted as distinct strings at depth t (see the class description)

        // returns the values that identify the group
        auto key() const { return std::tuple(u, errors, last, counters); }
    };

    // returns the share of the occurrences of a string of length L with distinct continuations (see _distinct_share)
    double distinct_share(uint64_t L) const { return _distinct_share[std::min<uint64_t>(L, _distinct_share.size() - 1)]; }

    /* returns the number of the steps at depth t of a search whose start has occ occurrences and length l, if the
       steps are b * sigma^t strings (see the class description); power = sigma^t */
    double occurring(double b, uint64_t t, double power, double occ, uint64_t l) const
    {
        return std::min(b, 1.0) * std::min(power, occ * distinct_share(l + t));
    }

    // the computed values of _costs->max_errors(beg, end), by beg * (m + 1) + end
    mutable gtl::flat_hash_map<uint64_t, uint64_t> _max_errors;

    // [0..p-1] _searches_of[i] = the indices of the searches that start with part i
    std::vector<std::vector<uint64_t>> _searches_of;

    /* the steps of the searches of part_steps (computed by its first call), for the bounds of the uniform partition:
       [0..|S|-1][t - 1] _part_strings[s][t - 1] * sigma^t are the steps of search s at depth t >= 1 (see for_each_step) */
    mutable std::vector<std::vector<double>> _part_strings;

    /* the steps of the predictions so far: the positions of the search orders, the cells (position, number of errors),
       the operations applied to the groups of strings, the groups at the end of each depth and n log n per sorting of
       n groups */
    mutable uint64_t _work = 0;

    // the number of values of last_op_t
    static constexpr uint64_t num_last_ops = TRAIL + 1;

    /* a slot of the groups of the edit distance at the depth t being computed (see for_each_counted_step): the groups
       with u positions of P, a number of errors and last operations are merged in slot
       ((u + k - t) * (k + 1) + errors) * num_last_ops + last, with the lane-wise minimum of their counters (see
       close_edit_depth); the slots with the same u form a list */
    struct slot_t {
        uint64_t counters; // the lane-wise minimum of the counters of the groups
        double num;        // the sum of counted_t::num of the groups
        double counted;    // the sum of counted_t::counted of the groups
        uint64_t stamp;    // the slot holds groups if stamp = _stamp
        uint32_t next;     // the next slot of the list of the slots with the same u, or no_slot
    };

    // the end of a list of slots
    static constexpr uint32_t no_slot = std::numeric_limits<uint32_t>::max();

    /* the buffers of the predictions (the predictions of one object do not run concurrently): the slots, the first slot
       of the list of each u + k - t and its stamp (see slot_t), and the current stamp */
    mutable std::vector<slot_t> _slots;
    mutable std::vector<uint32_t> _cell_head;
    mutable std::vector<uint64_t> _cell_stamp;
    mutable uint64_t _stamp = 0;
    // the groups at the current depth, and the groups at the next depth before merging (the Hamming distance)
    mutable std::vector<counted_t> _groups, _next;
    // see for_each_step
    mutable std::vector<double> _b_by_errors;
    mutable std::vector<direction_t> _dirs;
    // the bounds of the parts and the steps of the anchored search (see anchored_steps)
    mutable std::vector<uint64_t> _anchored_bounds;
    mutable search_t _anchored_search;

    /* returns the upper bound of a step of a search that has searched P[beg, end) and has the upper bound k_max, lowered
       to _costs->max_errors(beg, end) if lowered */
    uint64_t upper_bound(uint64_t beg, uint64_t end, uint64_t k_max, bool lowered) const
    {
        if (_costs == nullptr || !lowered) return k_max;
        auto [it, inserted] = _max_errors.try_emplace(beg * (_match.size() + 1) + end, 0);
        if (inserted) it->second = _costs->max_errors(beg, end);
        return std::min(k_max, it->second);
    }

    // the positions of P that a search matches symbol by symbol after its start, in the order of the search
    struct search_order_t {
        std::vector<uint64_t> pos;   // [0..U-1] the positions
        std::vector<uint64_t> k_max; // [0..U-1] the upper bound of the step of pos[u]
        // [0..U] the lower bound of the strings that have searched u positions (that of the step ending there, else 0)
        std::vector<uint64_t> k_min;
        // [0..U-1] the boundary of P before pos[u] in the direction of its step (where a deletion before pos[u] lies)
        std::vector<uint64_t> boundary;
        /* [0..U-1] the number of ways of the symbols that match both pos[u - 1] and pos[u], or -1 if pos[u] does not
           continue the direction of pos[u - 1] (at the first position and after a change of the direction) */
        std::vector<double> common_prev;
        /* [0..U] the boundary of P after pos[u - 1] in the direction of its step if the search extends the strings that
           end there by deletions after a step (see the class description), else no_trail */
        std::vector<uint64_t> trail;
    };

    // the value of search_order_t::trail at positions without deletions after a step
    static constexpr uint64_t no_trail = std::numeric_limits<uint64_t>::max();

    // the buffer of search_order
    mutable search_order_t _order;

    /* returns whether the edit-distance search computes the costs and the counters of the error model (the dynamic program
       of apm_edit with an error model, see apm_cost_band): if apm_costs_t::bounds_costs(), or with counters if the
       searches have no lower bounds above 0 (with mode.no_lower_bounds or if the scheme has none) */
    bool has_band(bounds_mode_t mode) const
    {
        return _edit && _costs != nullptr && (_costs->bounds_costs() ||
            (_costs->num_counters() > 0 && (mode.no_lower_bounds || !_scheme_lower_bounds)));
    }

    /**
     * @brief returns the positions of P that the search S matches symbol by symbol after its start (in _order)
     * @param S the search
     * @param bounds [0..p] the bounds of the parts
     * @param mode the error bounds of the steps
     * @param from_root whether S matches its first part symbol by symbol (see starts_at_root)
     */
    const search_order_t& search_order(const search_t& S, const std::vector<uint64_t>& bounds, bounds_mode_t mode,
        bool from_root) const;

    // merges the groups in next with the same key into groups (for the Hamming distance)
    void merge_groups(std::vector<counted_t>& next, std::vector<counted_t>& groups) const;

    /**
     * @brief computes the groups of strings at depth t that end with a text symbol from the groups at depth t - 1 in
     *        _groups (see the class description): each group continues at its next position pos[u] of the search order
     *        with a match, a substitution or a deletion before pos[u] by the rules of the class description, and after
     *        a step with deletions (see search_order_t::trail); the counter increments (and those of opening a gap) are
     *        added to the counters if counters are checked; with checks per string (if not band), the strings whose
     *        counters exceed their bounds are dropped; the strings with more errors than the upper bound are dropped
     * @param order the search order
     * @param counts whether counters are checked
     * @param band whether the counters are checked per cell (see close_edit_depth)
     * @param put called with the counters, the number of positions, the number of errors, the last operations, the
     *        number and the counted number of each group of continuations (see counted_t)
     */
    template <typename put_fnc_t>
    void extend_counted(const search_order_t& order, bool counts, bool band, put_fnc_t put) const;

    // starts the groups of a new depth of the edit distance: all slots are empty (see slot_t)
    void start_edit_depth() const;

    /* adds a group of strings at depth t with u positions of P, a number of errors, last operations, counters, a number
       num and a counted number to its slot (see slot_t); cell = u + k - t */
    void put_slot(uint64_t cell, uint8_t errors, last_op_t last, uint64_t counters, double num, double counted) const;

    /**
     * @brief completes the groups of strings at depth t of the edit distance (see the class description) in the slots
     *        by positions u in increasing order and sets _groups to them: if band, the counters are checked per cell as
     *        in the dynamic program of apm_edit (the groups at u are dropped if the lane-wise minimum of their counters
     *        exceeds a bound, else their counters are set to it); then the groups at u continue with an insertion of
     *        pos[u] (to position u + 1, no text symbol), not after a deletion or a deletion followed by a match of the
     *        position before pos[u], with at most the upper bound of pos[u] errors and at least the lower bound at the
     *        end of a step
     * @param t the depth
     * @param order the search order
     * @param counts whether counters are checked
     * @param band whether the counters are checked per cell
     */
    void close_edit_depth(uint64_t t, const search_order_t& order, bool counts, bool band) const;

    /**
     * @brief computes A(t) of the search S depth by depth (see the class description) and calls step(b) for each depth
     *        t, with b * sigma^t the steps at depth t (the A(t) strings); for the edit distance or with counters in the
     *        error model, see for_each_counted_step
     * @param S the search
     * @param bounds [0..p] the bounds of the parts
     * @param mode the error bounds of the steps
     * @param from_root whether S matches its first part symbol by symbol (see starts_at_root)
     * @param step the function called for each depth
     */
    template <typename step_fnc_t>
    void for_each_step(const search_t& S, const std::vector<uint64_t>& bounds, bounds_mode_t mode, bool from_root,
        step_fnc_t step) const;

    /* for_each_step for the edit distance or with counters in the error model: the depth t is the number of text symbols
       after the start, the strings are grouped by their numbers of positions of P, their numbers of errors, their
       counters and their last operations (see extend_counted and close_edit_depth), and for the edit distance, the steps at
       depth t are the children of the strings at depth t - 1 that the search extends (see the class description) */
    template <typename step_fnc_t>
    void for_each_counted_step(const search_t& S, const std::vector<uint64_t>& bounds, bounds_mode_t mode,
        bool from_root, step_fnc_t step) const;

    /**
     * @brief returns the predicted number of steps of the search S if its start has occ occurrences: the strings at the
     *        depths t that occur in the text (see occurring)
     */
    double search_steps(const search_t& S, const std::vector<uint64_t>& bounds, uint64_t occ, bounds_mode_t mode,
        bool from_root) const;

    // returns whether a symbol other than P[i] is free at a position i of P[beg, end) (then the part has no search context)
    bool has_free(uint64_t beg, uint64_t end) const
    {
        for (uint64_t i = beg; i < end; i++) if (_free[i]) return true;
        return false;
    }

    // returns whether the search S matches its first part symbol by symbol, starting with the n occurrences of the empty string
    bool starts_at_root(const search_t& S, const std::vector<uint64_t>& bounds, bool lowered) const
    {
        const uint64_t beg = bounds[S[0].part], end = bounds[S[0].part + 1];
        return !_from_contexts || has_free(beg, end) || upper_bound(beg, end, S[0].k_max, lowered) > 0;
    }

    // returns whether a step of the search S has a lower bound above its upper bound (lowered if lowered)
    bool is_infeasible(const search_t& S, const std::vector<uint64_t>& bounds, bool lowered) const;

    // returns whether lowering the upper bounds changes the upper bound of a step of a search
    bool lowers_a_bound(const std::vector<uint64_t>& bounds) const;

    /**
     * @brief returns the predicted number of steps of the searches of the scheme with the error bounds of mode; a search
     *        with a lower bound above an upper bound is not run and not counted
     * @param bounds [0..p] the bounds of the parts
     * @param occ [0..p-1] the numbers of occurrences of the parts
     * @param mode the error bounds
     */
    double scheme_steps(const std::vector<uint64_t>& bounds, const std::vector<uint64_t>& occ, bounds_mode_t mode) const;

    /**
     * @brief returns the predicted number of steps of the anchored search of apm_model_plan: the anchor, then the parts
     *        to its right, then the parts to its left, with the bounds of the parts outside the anchor
     * @param bounds [0..p] the bounds of the parts
     */
    double anchored_steps(const std::vector<uint64_t>& bounds) const;

public:
    /**
     * @brief prepares the prediction for a pattern P: the effective alphabet size sigma is 1 / (sum_c p_c^2), where p_c
     *        is the probability of symbol c in the text, and a text symbol c counts as sigma * p_c ways (without the
     *        probabilities, sigma is the alphabet size and each symbol counts as 1 way). At position i, a string
     *        continues without an error with P[i], and with an error by a substitution of P[i] by another symbol and,
     *        for the edit distance, by an insertion of P[i] or a deletion of any text symbol. With an error model, a
     *        string continues without an error with the free symbols of position i, and the errors are those with a
     *        cost of at most the cost bound.
     * @param scheme the search scheme
     * @param costs the error model compiled for P, or nullptr
     * @param metric the distance metric of the searches
     * @param P the pattern
     * @param sigma the alphabet size of the text
     * @param sym_p [0..255] the probability of each byte in the text (see apm_partition::symbol_probabilities), or empty
     * @param distinct_shares [0..L_max] the share of the occurrences of a string of each length with distinct
     *        continuations (see apm_distinct_shares)
     * @param n the text length
     * @param anchor_occ how often the anchor of the error model occurs in the text (see apm_costs_t::anchor)
     * @param params the query parameters (see apm_params_t::start_at_part_contexts, apm_params_t::lower_upper_bounds and
     *        apm_params_t::anchored_search)
     */
    template <typename pat_t>
    apm_step_prediction(const search_scheme_t& scheme, const apm_costs_t* costs, distance_metric_t metric, const pat_t& P,
        uint64_t sigma, const std::vector<double>& sym_p, std::vector<double> distinct_shares, uint64_t n,
        uint64_t anchor_occ, const apm_params_t& params);

    // returns the steps of the predictions so far (see _work)
    uint64_t work() const { return _work; }

    /**
     * @brief returns the predicted number of steps of the searches that start with part i if part i has occ
     *        occurrences; the other parts have the bounds of the uniform partition (see apm_partition::minimal_bounds)
     */
    double part_steps(uint64_t i, uint64_t occ) const;

    /**
     * @brief returns the searches that apm_model_plan runs for the parts with the bounds: without an error model, the
     *        searches of the scheme; with an error model, the alternative with the smallest prediction among
     *        1. the searches of the scheme with lowered upper bounds (if apm_params_t::lower_upper_bounds); for the edit
     *           distance, their lower bounds are 0 if an upper bound is lowered or if the costs are computed in the
     *           search (see apm_costs_t::bounds_costs),
     *        2. the searches of the scheme with the bounds of the scheme (if not apm_params_t::lower_upper_bounds, or for the
     *           edit distance if lowering an upper bound sets lower bounds to 0 that are kept otherwise); for the edit
     *           distance, their lower bounds are 0 if the costs are computed in the search, and
     *        3. the anchored search (if the error model has an anchor and apm_params_t::anchored_search).
     * @param bounds [0..p] the bounds of the parts
     * @param occ [0..p-1] the numbers of occurrences of the parts (read only for the first parts of the searches that
     *        start with the search context of their first part)
     */
    choice_t choose(const std::vector<uint64_t>& bounds, const std::vector<uint64_t>& occ) const;

    /**
     * @brief returns the predicted number of steps of the searches that apm_model_plan runs (see choose)
     * @param bounds [0..p] the bounds of the parts
     * @param occ [0..p-1] the numbers of occurrences of the parts
     */
    double steps(const std::vector<uint64_t>& bounds, const std::vector<uint64_t>& occ) const
    {
        return choose(bounds, occ).steps;
    }
};

#include <algorithms/move_rb/apm/step_prediction.tpp>

#endif // MOVE_R_ALGORITHMS_MOVE_RB_APM_STEP_PREDICTION_HPP
