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

#ifndef MOVE_R_MOVE_RB_SEARCH_CONTEXT_HPP
#define MOVE_R_MOVE_RB_SEARCH_CONTEXT_HPP

#include <move_rb/move_rb.hpp>

/**
 * @brief the state of a bidirectional search for a pattern P of length m: the suffix array interval [b, e] of P
 *        in T, the interval [b_R, e_R] of the reversed pattern in T^R (both have the same width), the input
 *        intervals of M_LF containing their limits and an SA-sample of the occurrences (for LOCATE); an
 *        extension returns a new context
 * @tparam query_support whether the context supports COUNT or LOCATE queries
 */
template <move_r_support support>
template <query_support_t query_support>
struct move_rb<support>::search_context_t {
    friend class move_rb;
    template <typename> friend class apm_hamming;
    template <typename> friend class apm_edit;
    template <locate_output_t> friend struct locate_context_t;
    template <query_support_t _query_support> friend struct extend_context_t;
    template <query_support_t _query_support> friend struct search_context_t;

    using extend_res_t = std::tuple<search_context_t, bool>;

    // ############################# DATA MEMBERS #############################

protected:
    const move_rb<support>* idx = nullptr; // the index this context queries

    uint64_t b, e; // [b, e] = suffix array interval of P in T
    uint64_t b_R, e_R; // [b_R, e_R] = suffix array interval of the reversed P in T^R
    uint64_t b_, e_; // indices of the input intervals of M_LF (of idx_fwd) containing b and e
    uint64_t b_R_, e_R_; // indices of the input intervals of M_LF of idx_bwd containing b_R and e_R

    // the SA-sample of the context (see sample_info_pack_t)
    STRING_UTILS_NO_UNIQUE_ADDRESS std::conditional_t<query_support == LOCATE, sample_info_pack_t, absent_t<10>> s;

    uint48_t m; // length of P
    uint64_t sym_lst; // the last-added symbol

    // the directions, for which the indices of the input intervals (b_, e_ resp. b_R_, e_R_) are valid:
    // NO_DIR (none), LEFT, RIGHT or BOTH_DIRS; an operation in direction dir updates only the indices
    // of direction dir, the others are recomputed when they are needed (see update_input_intervals)
    direction_t valid_interval_dirs = NO_DIR;

    uint8_t err; // number of errors (only used by approximate pattern matching)

    // entry of P in the forward order of the context tables (no_ctx_entry if it has none); it is set by an
    // extension of a context that has an entry (see context_tables)
    uint32_t ctx_entry = no_ctx_entry;
    // entry of P in the reverse order of the context tables
    uint32_t ctx_entry_R = no_ctx_entry;

    // ############################# INTERNAL METHODS #############################

protected:
    // returns the SA-sample of the context (computed with sample_at_run_boundary, if s.d() == BOTH_DIRS)
    sample_info_pack_t sample() const requires(query_support == LOCATE)
    { return s.d() == BOTH_DIRS ? sample_at_run_boundary() : s; }

    /**
     * @brief computes the SA-sample of a context whose interval contains an SA-sample of the index (see
     *        sample_inside): the SA-sample at the run start b or at the end of the BWT run containing b
     * @return the SA-sample
     */
    sample_info_pack_t sample_at_run_boundary() const requires(supports_locate);

    /**
     * @brief returns whether the suffix array interval [b, e] contains an SA-sample of the index, i.e., a
     *        BWT run starts at b or the BWT run containing b ends in [b, e] (see run_sample_positions)
     * @param index the index
     * @param b left limit of the interval
     * @param e right limit of the interval
     * @param b_ index of the input interval of M_LF containing b
     */
    static bool sample_inside(const move_rb& index, uint64_t b, uint64_t e, uint64_t b_)
    {
        const auto& fwd = index.idx_fwd;
        const auto& M_LF = fwd.M_LF();
        if (b == M_LF.p(b_) && fwd.run_starts_at(b_)) return true;
        return M_LF.p(fwd.last_interval_of_run(b_) + 1) - 1 <= e;
    }

    /**
     * @brief returns references to the interval variables of a search context; for dir == RIGHT, the
     *        variables of the forward and the backward interval are swapped
     * @tparam dir the direction
     * @tparam this_t type of the context (const or non-const)
     * @param t the context
     * @return (b, e, b_R, e_R, b_, e_, b_R_, e_R_, s) for dir == LEFT
     */
    template <direction_t dir, typename this_t>
    static auto vars_impl(this_t&& t)
    {
        if constexpr (dir == LEFT) return std::forward_as_tuple(t.b, t.e, t.b_R, t.e_R, t.b_, t.e_, t.b_R_, t.e_R_, t.s);
        else                       return std::forward_as_tuple(t.b_R, t.e_R, t.b, t.e, t.b_R_, t.e_R_, t.b_, t.e_, t.s);
    }

    template <direction_t dir> auto vars() { return vars_impl<dir>(*this); }
    template <direction_t dir> auto const_vars() const { return vars_impl<dir>(*this); }

    /**
     * @brief computes prev[0..max_sym] and next[0..max_sym] for the input intervals [b_, e_] of M_LF of
     *        direction dir: next[c] is the first and prev[c] the last input interval x in [b_, e_] with
     *        L'[x] = c. A symbol c occurs in L'[b_..e_] <=> next[c] <= prev[c] < r'; for any other symbol,
     *        the values come from the block samples of PSB_L' and lie outside of [b_, e_] (prev[c] = r' if
     *        there is none), or are LONG_MAX resp. LONG_MIN, if no block sample is used
     * @tparam dir extend direction
     * @param prev [0..max_sym] output
     * @param next [0..max_sym] output
     * @param max_sym largest symbol to compute the values for (max_sym < sigma)
     */
    template <direction_t dir>
    inline void build_prev_next(int64_t* prev, int64_t* next, uint64_t max_sym) const;

    // same as above for the range [b_, e_] of input intervals of M_LF
    template <direction_t dir>
    inline void build_prev_next(int64_t* prev, int64_t* next, uint64_t max_sym, uint64_t b_, uint64_t e_) const;

    // returns a buffer of the current thread for prev[0..sigma-1] followed by next[0..sigma-1] (see build_prev_next)
    static int64_t* prev_next_buffer(uint64_t sigma)
    {
        static thread_local std::vector<int64_t> buffer;
        buffer.resize(2 * sigma);
        return buffer.data();
    }

    /**
     * @brief computes sum_{a < sym} |{ i in [b, e] : L[i] = a }|, the number of positions in the
     *        dir-interval [b, e] whose symbol in L is smaller than sym; an extension with sym moves the left
     *        limit of the interval of the other direction by this number
     * @tparam dir extend direction
     * @param prev the array prev built for [b_, e_] (see build_prev_next)
     * @param next the array next built for [b_, e_]
     * @param sym the symbol
     * @param b left limit of the dir-interval
     * @param e right limit of the dir-interval
     * @param b_ index of the input interval of M_LF containing b
     * @param e_ index of the input interval of M_LF containing e
     */
    template <direction_t dir>
    inline uint64_t num_smaller(const int64_t* prev, const int64_t* next,
        uint64_t sym, uint64_t b, uint64_t e, uint64_t b_, uint64_t e_) const;

    // computes the same value as num_smaller with the wavelet tree over L' (see wavelet_tree_prefix_sum)
    template <direction_t dir>
    inline uint64_t num_smaller_wavelet(uint64_t sym, uint64_t b, uint64_t e, uint64_t b_, uint64_t e_) const;

    // computes the same value as num_smaller by scanning L'[b_..e_] (see max_scanned_intervals_per_wavelet_level)
    template <direction_t dir>
    inline uint64_t num_smaller_scanned(uint64_t sym, uint64_t b, uint64_t e, uint64_t b_, uint64_t e_) const;

    /**
     * @brief updates the indices of the input intervals of M_LF containing the limits of the dir-interval,
     *        if they are not valid (see valid_interval_dirs); the previous indices are the starting points of
     *        the search (see input_interval_near)
     * @tparam dir extend direction
     */
    template <direction_t dir>
    inline void update_input_intervals();

    /**
     * @brief computes the SA-sample of the context ctx, which was produced by a dir-extension of this context
     * @tparam dir extend direction
     * @param ctx the extended search context
     * @param prime the dir-interval of this context restricted to the runs of the extension symbol (see
     *        prime_tuple_t)
     */
    template <direction_t dir>
    inline void update_sample(search_context_t& ctx, const prime_tuple_t& prime) const;

    /**
     * @brief sets the context ctx to the extension of this context by a symbol c in direction dir: the
     *        dir-interval of ctx is LF applied to the positions with symbol c in the dir-interval of this
     *        context, which lie in the input intervals first..last of M_LF; also sets the SA-sample and the
     *        limits of the other direction. sym_lst is not set.
     * @tparam dir extend direction
     * @param ctx the context to set
     * @param first the first input interval of M_LF with the symbol c in the dir-interval of this context
     * @param last the last such input interval
     * @param b_R_new left limit of the interval of the other direction of ctx
     */
    template <direction_t dir>
    inline void extend_to_runs(search_context_t& ctx, uint64_t first, uint64_t last,
        uint64_t b_R_new) const;

    /**
     * @brief extends P with sym (see extend and extend_one_sided)
     * @tparam dir extend direction
     * @tparam one_sided true <=> only the interval of direction dir is maintained
     * @param sym the symbol to extend P with
     */
    template <direction_t dir, bool one_sided, typename sym_t>
    inline extend_res_t extend_impl(sym_t sym);

    // ############################# CONSTRUCTORS #############################

public:
    search_context_t() {}

    // constructs the context of the empty pattern for the index idx
    search_context_t(const move_rb<support>& idx) : idx(&idx) { reset(); }

    // ############################# MISC PUBLIC METHODS #############################

    // returns whether both contexts have the same interval [b, e] and pattern length (the other members may differ)
    bool operator==(const search_context_t& other) const
    {
        return b == other.b && e == other.e && m == other.m;
    }

    // orders the contexts by b, then by decreasing e, then by m
    bool operator<(const search_context_t& other) const
    {
        if (b != other.b) return b < other.b;
        if (e != other.e) return e > other.e;
        return m < other.m;
    }

    // hash of a search context, computed from [b, e] and the pattern length
    struct hash {
        // returns the hash of ctx
        inline uint64_t operator()(const search_context_t& ctx) const;
    };

    // resets the context to the empty pattern
    inline void reset();

    /**
     * @brief sets the number of errors (only used by approximate pattern matching)
     * @param err number of errors
     */
    inline void set_errors(uint64_t err) { this->err = err; }

    // ############################# ACCESS METHODS #############################

    // returns the last-added symbol
    inline uint8_t last_added_symbol() const { return sym_lst; }

    // returns the length of P
    inline uint64_t pattern_length() const { return m; }

    // returns the width of [b, e]
    inline uint64_t width() const { return e - b + 1; }

    // returns the number of occurrences of P
    inline uint64_t num_occ() const { return width(); }

    /**
     * @brief returns the number of occurrences of P computed from its dir-interval (the dir-interval is the valid one
     *        after a one-sided extension, see extend_one_sided)
     * @tparam dir the text direction
     */
    template <direction_t dir>
    inline uint64_t num_occ_in() const
    {
        const auto [beg, end] = sa_interval<dir>();
        return end - beg + 1;
    }

    // returns the number of errors (only used by approximate pattern matching)
    inline uint64_t errors() const { return err; }

    // returns (b_, e_, b_R_, e_R_), the input intervals containing the limits of the intervals; computes them
    // first, if they are not valid (see valid_interval_dirs)
    inline std::array<uint64_t, 4> input_intervals()
    {
        update_input_intervals<LEFT>();
        update_input_intervals<RIGHT>();
        return { b_, e_, b_R_, e_R_ };
    }

    /**
     * @brief returns the suffix array interval of P in T (dir == LEFT) or of the reversed P in T^R (dir == RIGHT)
     * @tparam dir the text direction
     */
    template <direction_t dir>
    inline std::tuple<uint64_t, uint64_t> sa_interval() const
    {
        if constexpr (dir == LEFT) return forward_sa_interval(); else return backward_sa_interval();
    }

    // returns [b, e], the suffix array interval of P in T
    inline std::tuple<uint64_t, uint64_t> forward_sa_interval() const { return std::tuple<uint64_t, uint64_t> { b, e }; }

    // returns [b_R, e_R], the suffix array interval of the reversed P in T^R
    inline std::tuple<uint64_t, uint64_t> backward_sa_interval() const { return std::tuple<uint64_t, uint64_t> { b_R, e_R }; }

    /**
     * @brief returns whether P is a prefix of the pattern of the other context, i.e., [b, e] contains the
     *        interval of the other context in T
     * @param other the other search context
     */
    inline bool is_prefix_of(const search_context_t& other) const { return b <= other.b && other.e <= e; }

    /**
     * @brief returns whether P is a suffix of the pattern of the other context, i.e., [b_R, e_R] contains the
     *        interval of the other context in T^R
     * @param other the other search context
     */
    inline bool is_suffix_of(const search_context_t& other) const
    { return b_R <= other.b_R && other.e_R <= e_R; }

    // ############################# QUERY METHODS #############################

    // returns a locate context for this search context (see locate_output_t)
    template <locate_output_t output = _occ>
    locate_context_t<output> locate_phase() const requires(query_support == LOCATE)
    { return locate_context_t<output>(*this); }

    // returns an extend context for this search context (see extend_context_t::prepare_extend_all)
    extend_context_t<query_support> extend_phase() { return extend_context_t<query_support>(*idx, *this); }

    /**
     * @brief extends P with sym
     * @param sym the symbol to extend P with
     * @param dir extend direction
     * @return (the search context of the extended P, whether the extended P occurs in the input); if it
     *         does not occur, the returned context is invalid
     */
    template <typename sym_t>
    inline extend_res_t extend(sym_t sym, direction_t dir);

    /**
     * @brief extends P with sym
     * @tparam dir extend direction
     * @param sym the symbol to extend P with
     * @return (the search context of the extended P, whether the extended P occurs in the input); if it
     *         does not occur, the returned context is invalid
     */
    template <direction_t dir, typename sym_t>
    inline extend_res_t extend(sym_t sym) { return extend_impl<dir, false>(sym); }

    /**
     * @brief extends P with sym and maintains only the interval of direction dir (T for LEFT, T^R for
     *        RIGHT); the interval of the other direction is not valid afterwards (counting contexts only)
     * @tparam dir extend direction
     * @param sym the symbol to extend P with
     * @return (the search context of the extended P, whether the extended P occurs in the input)
     */
    template <direction_t dir, typename sym_t>
    inline extend_res_t extend_one_sided(sym_t sym) requires(query_support == COUNT)
    { return extend_impl<dir, true>(sym); }

    /**
     * @brief fetches the data into the cache that extend<dir>(sym) reads in M_LF after the input intervals
     *        of the limits; fetches nothing if the extension reads the context tables or if the input
     *        intervals of direction dir are not valid (see update_input_intervals)
     * @tparam dir extend direction
     * @param sym the symbol to extend P with
     */
    template <direction_t dir, typename sym_t>
    inline void prefetch_extend(sym_t sym) const;

    /**
     * @brief reports all extensions of P by one symbol that occur in the input, in the order of the internal
     *        alphabet
     * @tparam dir extend direction
     * @tparam scope which symbols are reported (see extension_scope_t)
     * @tparam fnc_t type of the function report
     * @param report function called as report(symbol, extended context); if it returns false, the
     *        enumeration stops
     */
    template <direction_t dir, extension_scope_t scope = TERM, typename fnc_t>
    inline void for_each_extension(fnc_t report);
};

#include <move_rb/search_context.tpp>

#endif // MOVE_R_MOVE_RB_SEARCH_CONTEXT_HPP
