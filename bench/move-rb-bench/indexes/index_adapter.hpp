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
 * @file index_adapter.hpp
 * @brief the interface through which the approximate-pattern-matching algorithms of move_rb (apm_hamming and
 *        apm_edit) query a foreign index (br_index_adapter, columba_apm_adapter), and the parts both adapters share
 */

#ifndef MOVE_R_BENCH_MOVE_RB_BENCH_INDEXES_INDEX_ADAPTER_HPP
#define MOVE_R_BENCH_MOVE_RB_BENCH_INDEXES_INDEX_ADAPTER_HPP

#include <concepts>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

#include <misc/apm.hpp>

/**
 * @brief extracts the arity and argument types of a callable (e.g. a lambda)
 * @tparam T the callable type
 */
template <typename T>
struct function_traits : public function_traits<decltype(&T::operator())> {};

template <typename class_t, typename return_t, typename... args>
struct function_traits<return_t(class_t::*)(args...) const> {
    static constexpr uint64_t arity = sizeof...(args);

    template <uint64_t N>
    using argument_type = std::tuple_element_t<N, std::tuple<args...>>;
};

/**
 * @brief a locate context: enumerates the occurrences of a matched pattern, starting at the occurrence of its
 *        toehold (an occurrence whose suffix array position is known) and continuing to both sides in suffix array
 *        order
 */
template <typename L>
concept locate_context = requires(L lc) {
    { lc.num_occ_rem() } -> std::convertible_to<uint64_t>; // number of occurrences not reported yet
    { lc.compute_center() } -> std::convertible_to<uint64_t>; // suffix array position of the toehold
    { lc.locate() } -> std::convertible_to<std::vector<uint64_t>>; // the occurrences not reported yet
};

/**
 * @brief an extend context: returns, one by one, the search contexts obtained by extending a search context with
 *        each symbol in a given direction for which the extended pattern occurs
 */
template <typename E, typename search_context_t>
concept extend_context = requires(E ec, direction_t dir) {
    // prepare the extensions in direction dir, the terminator excluded (see extension_scope_t)
    { ec.template prepare_extend_all<NO_TERM>(dir) } -> std::same_as<E&>;
    { ec.can_extend() } -> std::convertible_to<bool>; // is there another applicable symbol?
    { ec.extend_next() } -> std::same_as<search_context_t>; // the context extended by the next symbol
};

/**
 * @brief a bidirectional search context: a matched pattern P with its suffix array interval in T and the suffix
 *        array interval of P^R in T^R
 */
template <typename C>
concept search_context = requires(C ctx, const C cctx, uint8_t sym, direction_t dir, uint64_t p) {
    typename C::hash; // hashes the context by its forward SA-interval (and shift/depth)
    { cctx.forward_sa_interval() } -> std::convertible_to<std::tuple<uint64_t, uint64_t>>;
    { cctx.backward_sa_interval() } -> std::convertible_to<std::tuple<uint64_t, uint64_t>>;
    { cctx.num_occ() } -> std::convertible_to<uint64_t>;
    { cctx.pattern_length() } -> std::convertible_to<uint64_t>;
    { cctx.errors() } -> std::convertible_to<uint64_t>;
    ctx.set_errors(p);
    { cctx.last_added_symbol() } -> std::convertible_to<uint8_t>;
    { cctx.extend(sym, dir) } -> std::convertible_to<std::tuple<C, bool>>;
    ctx.extend_phase();
    { cctx == cctx } -> std::convertible_to<bool>;
    { cctx < cctx } -> std::convertible_to<bool>;
};

/**
 * @brief a search context that also supports the locate phase and stores the depth, shift and reported flag used by
 *        the edit-distance algorithm (see adapter_locate_sample_t)
 */
template <typename C>
concept locate_search_context = search_context<C> &&
    requires(C ctx, const C cctx, uint64_t p, bool b) {
        { cctx.depth() } -> std::convertible_to<uint64_t>;
        ctx.set_depth(p);
        { cctx.shift() } -> std::convertible_to<uint64_t>;
        ctx.set_shift(p);
        { cctx.reported() } -> std::convertible_to<bool>;
        ctx.set_reported(b);
        cctx.locate_phase();
    };

/**
 * @brief the interface of an index adapter: the text length, the alphabet maps, the creation of empty query
 *        contexts, exact and approximate count and locate, and the nested search, extend and locate context types
 */
template <typename T>
concept index_adapter =
    // associated types and capability flags
    requires {
        typename T::inp_t;
        requires std::convertible_to<decltype(T::supports_locate), bool>;
        requires std::convertible_to<decltype(T::supports_multiple_locate), bool>;
    } &&
    // index-level query surface
    requires(const T idx, typename T::inp_t P, search_scheme_t scheme, uint8_t sym, uint8_t isym) {
        { idx.n } -> std::convertible_to<uint64_t>; // text length (incl. sentinel)
        { idx.sigma } -> std::convertible_to<uint64_t>; // alphabet size
        { idx.input_size() } -> std::convertible_to<uint64_t>; // text length (excl. sentinel)
        { idx.map_symbol(sym) } -> std::convertible_to<uint8_t>;
        { idx.unmap_symbol(isym) } -> std::convertible_to<uint8_t>;
        { idx.map_int() } -> std::convertible_to<const std::vector<uint8_t>&>;
        idx.forward_index(); // index answering the exact count/locate queries
        { idx.count(P) } -> std::convertible_to<uint64_t>; // exact count
        { idx.count_hamming_dist(P, scheme) } -> std::convertible_to<uint64_t>;
        idx.template empty_context<COUNT>();
        idx.template empty_context<LOCATE>();
        { idx.locate(P, scheme, HAMMING_DISTANCE) }
            -> std::convertible_to<std::vector<aprx_occ_t<>>>;
        { idx.locate(P, scheme, EDIT_DISTANCE) }
            -> std::convertible_to<std::vector<aprx_occ_t<>>>;
    } &&
    // nested query-context types
    search_context<typename T::template search_context_t<COUNT>> &&
    locate_search_context<typename T::template search_context_t<LOCATE>> &&
    extend_context<typename T::template extend_context_t<LOCATE>, typename T::template search_context_t<LOCATE>> &&
    locate_context<typename T::locate_context_t>;

// ############################# PARTS SHARED BY THE ADAPTERS #############################

// the members of a LOCATE search context of an adapter (the same as in move_rb's search_context_t)
struct adapter_locate_sample_t {
    direction_t dir = LEFT; // direction of the last extension
    uint64_t dpth = 0; // depth, the number of symbols added to the context (used by apm_edit)
    uint64_t shft = 0; // shift, added to each reported occurrence
    bool rprtd = false; // true <=> the occurrences of the context have been reported

    // the direction of the last extension (read and set as the field d of move_rb's sample_info_pack_t)
    direction_t d() const { return dir; }
    void set_d(direction_t value) { dir = value; }
};

/**
 * @brief returns a copy of a search context of an adapter whose pattern is extended by one symbol; the caller
 *        sets the suffix array intervals of the extended pattern
 * @tparam query_support COUNT or LOCATE
 * @param ctx a search context with the members m, err, sym_lst and (LOCATE) s
 * @param sym the added symbol
 * @param dir the direction it is added in
 */
template <query_support_t query_support, typename ctx_t>
inline ctx_t extended_context(const ctx_t& ctx, uint8_t sym, direction_t dir)
{
    ctx_t extended = ctx;
    extended.m = ctx.m + 1;
    extended.err = 0;
    extended.sym_lst = sym;

    if constexpr (query_support == LOCATE) {
        extended.s.set_d(dir);
        extended.s.dpth = ctx.s.dpth + 1;
        extended.s.rprtd = false;
    }

    return extended;
}

// returns whether two search contexts of an adapter have the same forward suffix array interval (and, LOCATE, the
// same shift and depth)
template <query_support_t query_support, typename ctx_t>
inline bool adapter_contexts_equal(const ctx_t& ctx, const ctx_t& other)
{
    bool equal = ctx.forward_sa_interval() == other.forward_sa_interval();
    if constexpr (query_support == LOCATE) equal = equal && ctx.s.shft == other.s.shft && ctx.s.dpth == other.s.dpth;
    return equal;
}

// orders search contexts of an adapter by the beginning of the forward suffix array interval (ascending), its end
// (descending), (LOCATE) the shift and the depth, and the pattern length
template <query_support_t query_support, typename ctx_t>
inline bool adapter_context_less(const ctx_t& ctx, const ctx_t& other)
{
    auto [b, e] = ctx.forward_sa_interval();
    auto [ob, oe] = other.forward_sa_interval();
    if (b != ob) return b < ob;
    if (e != oe) return e > oe;

    if constexpr (query_support == LOCATE) {
        if (ctx.s.shft != other.s.shft) return ctx.s.shft < other.s.shft;
        if (ctx.s.dpth != other.s.dpth) return ctx.s.dpth < other.s.dpth;
    }

    return ctx.m < other.m;
}

// returns a hash of the forward suffix array interval (and, LOCATE, the shift and the depth) of a search context
template <query_support_t query_support, typename ctx_t>
inline uint64_t adapter_context_hash(const ctx_t& ctx)
{
    auto [b, e] = ctx.forward_sa_interval();
    uint64_t hash = hash_64(b);
    hash_combine(hash, hash_64(e));

    if constexpr (query_support == LOCATE) {
        hash_combine(hash, hash_64(ctx.s.shft));
        hash_combine(hash, hash_64(ctx.s.dpth));
    }

    return hash;
}

/**
 * @brief reports the occurrences of a pattern in suffix array order around its toehold: SA[c], then SA[c-1],
 *        SA[c-2], ..., then SA[c+1], SA[c+2], ...
 * @param report called as report(occurrence) or, if it takes two arguments, as report(i, occurrence), where i is
 *        the suffix array position of the occurrence
 * @param center_sa c, the suffix array position of the toehold
 * @param center_occ SA[c]
 * @param left [0..x-1] left[j] = SA[c-1-j]
 * @param right [0..y-1] right[j] = SA[c+1+j]
 * @param shift added to each occurrence
 */
template <typename report_fnc_t, typename vector_t>
inline void report_centered_occurrences(report_fnc_t report, uint64_t center_sa, uint64_t center_occ,
    const vector_t& left, const vector_t& right, uint64_t shift)
{
    static constexpr bool report_pos = function_traits<report_fnc_t>::arity > 1;

    if constexpr (report_pos) report(center_sa, center_occ + shift);
    else report(center_occ + shift);

    for (uint64_t j = 0; j < left.size(); j++) {
        if constexpr (report_pos) report(center_sa - 1 - j, (uint64_t) left[j] + shift);
        else report((uint64_t) left[j] + shift);
    }

    for (uint64_t j = 0; j < right.size(); j++) {
        if constexpr (report_pos) report(center_sa + 1 + j, (uint64_t) right[j] + shift);
        else report((uint64_t) right[j] + shift);
    }
}

/**
 * @brief matches P from left to right with a COUNT search context of an adapter and returns its number of
 *        occurrences
 * @param index the adapter
 * @param P the pattern
 */
template <typename adapter_t>
inline uint64_t count_by_appending(const adapter_t& index, const std::string& P)
{
    auto ctx = index.template empty_context<COUNT>();

    for (uint64_t i = 0; i < P.size(); i++) {
        auto [extended, occurs] = ctx.extend(P[i], RIGHT);
        if (!occurs) return 0;
        ctx = extended;
    }

    return ctx.num_occ();
}

/**
 * @brief matches P from left to right with a LOCATE search context of an adapter and reports its occurrences
 * @param index the adapter
 * @param P the pattern
 * @param report called as report(i, occurrence), where i is the suffix array position of the occurrence
 */
template <typename adapter_t, typename report_fnc_t>
inline void locate_by_appending(const adapter_t& index, const std::string& P, report_fnc_t report)
{
    auto ctx = index.template empty_context<LOCATE>();

    for (uint64_t i = 0; i < P.size(); i++) {
        auto [extended, occurs] = ctx.extend(P[i], RIGHT);
        if (!occurs) return;
        ctx = extended;
    }

    ctx.locate_phase().locate([&](uint64_t i, uint64_t occ) { report(i, occ); });
}

#endif // MOVE_R_BENCH_MOVE_RB_BENCH_INDEXES_INDEX_ADAPTER_HPP
