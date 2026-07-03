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

#include <concepts>
#include <tuple>
#include <vector>

#include <misc/apm.hpp>

/**
 * @brief a locate context: enumerates the occurrences of a matched pattern, centered
 *        on a toehold occurrence, in suffix-array order
 */
template <typename L, typename pos_t>
concept locate_context = requires(L lc) {
    { lc.num_occ_rem() } -> std::convertible_to<pos_t>; // remaining (not yet reported) occurrences
    { lc.compute_center() } -> std::convertible_to<pos_t>; // SA-index of the centered occurrence
    { lc.locate() } -> std::convertible_to<std::vector<pos_t>>; // all (remaining) occurrences
};

/**
 * @brief an extend context: yields, one by one, the search contexts obtained by
 *        extending a search context with every applicable symbol in a given direction
 */
template <typename E, typename search_context_t>
concept extend_context = requires(E ec, direction_t dir) {
    { ec.prepare_extend_all(dir) } -> std::same_as<E&>; // prepare extensions in direction dir
    { ec.can_extend() } -> std::convertible_to<bool>; // is there another applicable symbol?
    { ec.extend_next() } -> std::same_as<search_context_t>; // the context extended by the next symbol
};

/**
 * @brief a bidirectional search context: a matched pattern together with its SA-intervals
 *        and the bookkeeping the search/locate phases need
 */
template <typename C, typename sym_t, typename pos_t>
concept search_context = requires(C ctx, const C cctx, sym_t sym, direction_t dir, pos_t p) {
    typename C::hash; // hashes the context by its forward SA-interval (and shift/depth)
    { cctx.forward_sa_interval() } -> std::convertible_to<std::tuple<pos_t, pos_t>>;
    { cctx.backward_sa_interval() } -> std::convertible_to<std::tuple<pos_t, pos_t>>;
    { cctx.num_occ() } -> std::convertible_to<pos_t>;
    { cctx.length() } -> std::convertible_to<pos_t>;
    { cctx.errors() } -> std::convertible_to<pos_t>;
    ctx.set_errors(p);
    { cctx.last_added_symbol() } -> std::convertible_to<sym_t>;
    { cctx.last_direction() } -> std::convertible_to<direction_t>;
    { cctx.extend(sym, dir) } -> std::convertible_to<std::tuple<C, bool>>;
    ctx.extend_phase();
    { cctx == cctx } -> std::convertible_to<bool>;
    { cctx < cctx } -> std::convertible_to<bool>;
};

/**
 * @brief a search context that additionally supports the locate phase (error/shift/depth
 *        bookkeeping used by the edit-distance algorithm)
 */
template <typename C, typename sym_t, typename pos_t>
concept locate_search_context = search_context<C, sym_t, pos_t> &&
    requires(C ctx, const C cctx, pos_t p, bool b) {
        { cctx.depth() } -> std::convertible_to<pos_t>;
        ctx.set_depth(p);
        { cctx.shift() } -> std::convertible_to<pos_t>;
        ctx.set_shift(p);
        { cctx.reported() } -> std::convertible_to<bool>;
        ctx.set_reported(b);
        cctx.locate_phase();
    };

/**
 * @brief the full interface an index adapter exposes to the approximate-pattern-matching
 *        algorithms: alphabet/text accessors, query-context factories and the nested
 *        search/extend/locate context types
 */
template <typename T>
concept index_adapter =
    // associated types and capability flags
    requires {
        typename T::pos_type;
        typename T::sym_type;
        typename T::i_sym_t;
        typename T::inp_t;
        requires std::convertible_to<decltype(T::supports_locate), bool>;
        requires std::convertible_to<decltype(T::supports_multiple_locate), bool>;
    } &&
    // index-level query surface
    requires(const T idx, typename T::inp_t P, search_scheme_t scheme,
             typename T::sym_type sym, typename T::i_sym_t isym) {
        { idx.n } -> std::convertible_to<typename T::pos_type>; // text length (incl. sentinel)
        { idx.sigma } -> std::convertible_to<typename T::pos_type>; // alphabet size
        { idx.input_size() } -> std::convertible_to<typename T::pos_type>; // text length (excl. sentinel)
        { idx.map_symbol(sym) } -> std::convertible_to<typename T::i_sym_t>;
        { idx.unmap_symbol(isym) } -> std::convertible_to<typename T::sym_type>;
        { idx.map_int() } -> std::convertible_to<const std::vector<uint8_t>&>;
        idx.forward_index(); // index answering the exact count/locate queries
        { idx.count(P) } -> std::convertible_to<typename T::pos_type>; // exact count
        { idx.count_hamming_dist(P, scheme) } -> std::convertible_to<typename T::pos_type>;
        idx.template empty_context<COUNT>();
        idx.template empty_context<LOCATE>();
        { idx.template locate<HAMMING_DISTANCE>(P, scheme) }
            -> std::convertible_to<std::vector<aprx_occ_t<typename T::pos_type>>>;
        { idx.template locate<EDIT_DISTANCE>(P, scheme) }
            -> std::convertible_to<std::vector<aprx_occ_t<typename T::pos_type>>>;
    } &&
    // nested query-context types
    search_context<typename T::template search_context_t<COUNT>, typename T::sym_type, typename T::pos_type> &&
    locate_search_context<typename T::template search_context_t<LOCATE>, typename T::sym_type, typename T::pos_type> &&
    extend_context<typename T::template extend_context_t<LOCATE>, typename T::template search_context_t<LOCATE>> &&
    locate_context<typename T::locate_context_t, typename T::pos_type>;
