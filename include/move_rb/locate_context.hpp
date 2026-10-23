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

#ifndef MOVE_R_MOVE_RB_LOCATE_CONTEXT_HPP
#define MOVE_R_MOVE_RB_LOCATE_CONTEXT_HPP

#include <move_rb/move_rb.hpp>

/**
 * @brief the state of locating the occurrences SA[b..e] of a search context. Locating starts at a position c in
 *        [b, e] with a known value SA[c] (see compute_center) and reports the occurrences to its right and to
 *        its left: with M_Phi^{-1} and M_Phi (_locate_move) or with the rlz of SA^d (_locate_rlzsa).
 * @tparam output whether each occurrence is reported with its position in the suffix array (see
 *         locate_output_t)
 */
template <move_r_support support>
template <locate_output_t output>
struct move_rb<support>::locate_context_t {
    // true <=> each occurrence is reported with its position in the suffix array (see locate_output_t)
    static constexpr bool reports_sa_positions = output == _occ_and_sa_pos;

    // decoder of the rlz of SA^d at a suffix array position
    using SA_rlz_ctx_t = rlz<_diff_access>::decode_context_t;

    // ############################# DATA MEMBERS #############################

protected:
    const move_rb<support>* idx = nullptr; // the index this context queries
    const search_context_t<LOCATE>* ctx = nullptr; // the search context whose occurrences are located

    // direction of the current scan from c: RIGHT (towards e) or LEFT (towards b); NO_DIR before the first
    // occurrence
    direction_t scan_dir;
    bool centered; // true <=> c and SA_c are computed (see compute_center)
    uint64_t occ_rem; // number of occurrences that are not reported yet
    uint64_t c; // [b..e] the position from which the occurrences are reported to both sides (see compute_center)
    uint64_t SA_c; // SA[c]
    uint64_t i; // [b..e] the current position
    uint64_t SA_i; // SA[i]

    // index of the input interval of M_Phi^{-1} (scan_dir = RIGHT) resp. M_Phi (scan_dir = LEFT) containing SA_i
    STRING_UTILS_NO_UNIQUE_ADDRESS std::conditional_t<!stores_SA_pos, uint64_t, absent_t<13>> x_phi;

    // decoders of the rlz of SA^d for the positions left resp. right of c
    STRING_UTILS_NO_UNIQUE_ADDRESS std::conditional_t<stores_SA_pos, SA_rlz_ctx_t, absent_t<15>> SA_rlz_l;
    STRING_UTILS_NO_UNIQUE_ADDRESS std::conditional_t<stores_SA_pos, SA_rlz_ctx_t, absent_t<16>> SA_rlz_r;

    // index in buf of the next occurrence to report
    STRING_UTILS_NO_UNIQUE_ADDRESS std::conditional_t<!stores_SA_pos && reports_sa_positions,
        uint64_t, absent_t<17>> buf_pos;
    /* the occurrences right of c computed by compute_center: [0..k-1] buf[j] = SA[c + 1 + j], where k is the number
       of Phi^{-1} steps compute_center made (empty, if it made none) */
    STRING_UTILS_NO_UNIQUE_ADDRESS std::conditional_t<!stores_SA_pos && reports_sa_positions,
        std::vector<uint64_t>, absent_t<18>> buf;

    // ############################# INTERNAL METHODS #############################

    // computes the first occurrence SA_i to report and sets the context to it
    inline uint64_t first_occ() requires(supports_multiple_locate);

    // returns whether i = e, i.e., LCP[i + 1] < m
    inline bool is_right_limit() const
        requires(supports_multiple_locate && !stores_SA_pos);

    // reports the remaining occurrences to the left and to the right of c (see locate)
    template <typename report_fnc_t>
    inline void locate_from_sample(report_fnc_t report) requires(supports_multiple_locate);

    // ############################# CONSTRUCTORS #############################

public:
    locate_context_t() {}

    // constructs the locate context of the search context ctx
    locate_context_t(const search_context_t<LOCATE>& ctx)
        : idx(ctx.idx), ctx(&ctx)
    {
        occ_rem = ctx.num_occ();
        scan_dir = NO_DIR;
        centered = false;
        if constexpr (!stores_SA_pos && reports_sa_positions) buf_pos = 0;
    }

    // ############################# ACCESS METHODS #############################

    // returns the number of occurrences that are not reported yet
    inline uint64_t num_occ_rem() const { return occ_rem; }

    // ############################# QUERY METHODS #############################

    /**
     * @brief returns the occurrence given by the SA-sample of the search context in O(1) time; the locate
     *        context is not modified
     */
    inline uint64_t sampled_occ() const;

    /**
     * @brief computes c and SA_c = SA[c]: c is the position of the SA-sample of the search context, or b, if
     *        the sample is a backward one and b - (start of its BWT run) < m - i (SA[b] is then decoded from the
     *        start of the run); if the position of the SA-sample is not known, c = e - k, where k is the number
     *        of Phi^{-1} steps from SA_c to SA[e]
     * @return c
     */
    inline uint64_t compute_center() requires(supports_multiple_locate);

    // returns the next occurrence
    inline uint64_t next_occ() requires(supports_multiple_locate);

    /**
     * @brief reports the remaining occurrences
     * @tparam report_fnc_t type of the function report
     * @param report function called as report(occ) for _occ and as report(i, occ) with occ = SA[i] for
     *        _occ_and_sa_pos (see locate_output_t)
     */
    template <typename report_fnc_t>
    inline void locate(report_fnc_t report) requires(supports_multiple_locate);

    // returns the remaining occurrences
    std::vector<uint48_t> locate() requires(supports_multiple_locate);
};

#include <move_rb/locate_context.tpp>

#endif // MOVE_R_MOVE_RB_LOCATE_CONTEXT_HPP
