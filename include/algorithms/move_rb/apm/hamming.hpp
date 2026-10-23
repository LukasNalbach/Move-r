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

#ifndef MOVE_R_ALGORITHMS_MOVE_RB_APM_HAMMING_HPP
#define MOVE_R_ALGORITHMS_MOVE_RB_APM_HAMMING_HPP

#include <algorithm>
#include <limits>

#include <misc/apm.hpp>
#include <algorithms/move_rb/apm/model.hpp>
#include <algorithms/move_rb/apm/model_plan.hpp>
#include <algorithms/move_rb/apm/partition.hpp>
#include <algorithms/move_rb/apm/search_base.hpp>

/**
 * @brief approximate pattern matching with Hamming distance on a bidirectional index with search schemes: each search
 *        of the scheme matches the positions of P in the order of its steps (see search_positions) with a depth-first
 *        search over the extensions of the search contexts, within the error bounds of its steps
 * @tparam idx_t the bidirectional index type (move_rb, or a benchmark adapter with the same interface)
 */
template <typename idx_t>
class apm_hamming : public apm_search_base<idx_t> {
protected:
    // ############################# TYPES AND CONSTANTS #############################

    using base_t = apm_search_base<idx_t>;
    using typename base_t::inp_t;
    using base_t::_index;
    using base_t::_params;
    using base_t::next_is_separator;

    template <query_support_t query_support>
    using search_context_t = typename idx_t::template search_context_t<query_support>;

    // marks a search state that added no mismatch; the largest value of search_state_t::added_pos
    static constexpr uint64_t no_mismatch = std::numeric_limits<uint32_t>::max();

    // marks a located context without a CIGAR
    static constexpr uint32_t no_cigar = std::numeric_limits<uint32_t>::max();

    // a search context of a LOCATE search and the index of its CIGAR in _cigars (or no_cigar)
    struct located_context_t {
        search_context_t<LOCATE> ctx;
        uint32_t cigar_idx = no_cigar;

        located_context_t() = default;
        located_context_t(const search_context_t<LOCATE>& c, uint32_t ci = no_cigar) : ctx(c), cigar_idx(ci) {}

        // two located contexts are equal if their suffix array intervals in T start at the same position
        bool operator==(const located_context_t& o) const {
            return std::get<0>(ctx.forward_sa_interval()) ==
                   std::get<0>(o.ctx.forward_sa_interval());
        }

        struct hash {
            uint64_t operator()(const located_context_t& node) const {
                return hash_64(std::get<0>(node.ctx.forward_sa_interval()));
            }
        };
    };

    // the result of a search: a set of located contexts (LOCATE) or of search contexts (COUNT)
    template <query_support_t query_support>
    using search_context_set_t = std::conditional_t<query_support == LOCATE,
        tsl::sparse_set<located_context_t, typename located_context_t::hash>,
        tsl::sparse_set<search_context_t<query_support>, typename search_context_t<query_support>::hash>>;

    // a position of P in the order in which a search matches P (see search_positions)
    struct search_position_t {
        uint64_t step; // the step of the search that matches pos
        direction_t dir; // the direction of the step
        uint64_t pos; // the position in P
        uint64_t num_left; // the number of positions after pos that the step of the next position matches
    };

    // ############################# DATA MEMBERS #############################

    /* set by locate before it calls search or search_model: true <=> they compute the CIGAR of each located context
       (CIGAR mode); they then clear _path_mismatches and _cigars first */
    mutable bool _computes_cigars = false;

    /* the mismatches on the current path of the depth-first search (CIGAR mode only): [0..number of mismatches on the
       path - 1] (position in P, text symbol) of the i-th mismatch on the path */
    mutable std::vector<std::pair<uint64_t, uint8_t>> _path_mismatches;

    /* the CIGARs of the located contexts (CIGAR mode only): [0..number of located contexts - 1] the CIGAR of the
       located context with cigar_idx = i; returned by locate */
    mutable std::vector<cigar_t<uint8_t>> _cigars;

    // ############################# INTERNAL METHODS #############################

    /**
     * @brief builds the CIGAR of an occurrence without insertions and deletions
     * @param mismatches the mismatches of the occurrence, as pairs (position in P, text symbol)
     * @param m length of P
     * @return the CIGAR
     */
    static cigar_t<uint8_t> build_cigar(std::vector<std::pair<uint64_t, uint8_t>> mismatches, uint64_t m);

    // appends the CIGAR of the mismatches in _path_mismatches to _cigars and returns its index
    uint32_t add_cigar(uint64_t m) const
    {
        _cigars.emplace_back(build_cigar(_path_mismatches, m));
        return uint32_t(_cigars.size() - 1);
    }

    /**
     * @brief returns the positions of P in the order a search matches them: step by step, each part from its end
     *        that is adjacent to the previous parts; num_left of the last position of a step is the length of the part
     *        of the next step, else the number of positions of its step after it
     * @tparam partition_t type of the partition
     * @param S the search
     * @param partition the partition of P into the parts of the search
     * @param m the length of P
     * @return [0..m-1] the i-th position the search matches
     */
    template <typename partition_t>
    static std::vector<search_position_t> search_positions(const search_t& S, const partition_t& partition, uint64_t m);

    /**
     * @brief adds a search context that matches all of P to the result; of two located contexts with the same suffix
     *        array interval, the one whose last extension is to the left is kept
     * @tparam query_support COUNT or LOCATE
     * @param result the result of a search
     * @param ctx the search context
     * @param m the length of P
     */
    template <query_support_t query_support>
    void add_result(search_context_set_t<query_support>& result, search_context_t<query_support>& ctx, uint64_t m) const;

public:
    // ############################# CONSTRUCTORS #############################

    /**
     * @brief constructs the queries on an index
     * @param index the index
     * @param params query parameters
     */
    apm_hamming(const idx_t& index, apm_params_t params = {}) : base_t(index, params) {}

    // ############################# QUERY METHODS #############################

    /**
     * @brief counts the occurrences of P with at most k mismatches (with an error model: of cost at most its bound)
     * @param P the pattern
     * @param scheme the search scheme (provides k)
     * @return the number of occurrences
     */
    uint64_t count(const inp_t& P, const search_scheme_t& scheme) const;

    /**
     * @brief counts the occurrences of P with at most k mismatches by their number of mismatches (with an error
     *        model: by their costs)
     * @param P the pattern
     * @param scheme the search scheme (provides k)
     * @return [0..k] (with an error model [0..bound]) h[e] is the number of occurrences with e mismatches (cost e);
     *         for k >= m, all n - 1 occurrences are counted in h[0]
     */
    std::vector<uint64_t> count_histogram(const inp_t& P, const search_scheme_t& scheme) const;

    /**
     * @brief locates the occurrences of P with at most k mismatches (with an error model: of cost at most its bound)
     * @tparam mode whether each occurrence is reported with a CIGAR index
     * @tparam report_fnc_t type of the report function
     * @param P the pattern (P.size() < 2^aprx_occ_len_bits)
     * @param scheme the search scheme (provides k)
     * @param report called as report({pos, len, err}) or, in CIGAR mode, report({pos, len, err, cigar index}) for
     *        each occurrence
     * @return in CIGAR mode, the CIGARs (one per suffix array interval); empty in NO_CIGAR mode
     */
    template <cigar_mode_t mode = NO_CIGAR, typename report_fnc_t>
    std::vector<cigar_t<uint8_t>> locate(const inp_t& P, const search_scheme_t& scheme, report_fnc_t report) const;

    /**
     * @brief executes a search scheme for P
     * @tparam query_support COUNT or LOCATE
     * @param P the pattern
     * @param scheme the search scheme
     * @return the search contexts of the strings with at most k mismatches to P
     */
    template <query_support_t query_support>
    auto search(const inp_t& P, const search_scheme_t& scheme) const -> search_context_set_t<query_support>;

    /**
     * @brief executes a search scheme for P with the error model params.costs (see apm_model_plan); the error count
     *        of a search context is the cost of its alignment
     * @tparam query_support COUNT or LOCATE
     * @param P the pattern
     * @param scheme the search scheme (its k must be at least params.costs->max_errors())
     * @return the search contexts of the strings with a valid alignment of P
     */
    template <query_support_t query_support>
    auto search_model(const inp_t& P, const search_scheme_t& scheme) const -> search_context_set_t<query_support>;
};

#include <algorithms/move_rb/apm/hamming.tpp>

#endif // MOVE_R_ALGORITHMS_MOVE_RB_APM_HAMMING_HPP
