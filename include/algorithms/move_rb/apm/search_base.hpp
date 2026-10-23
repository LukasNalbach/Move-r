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

#ifndef MOVE_R_ALGORITHMS_MOVE_RB_APM_SEARCH_BASE_HPP
#define MOVE_R_ALGORITHMS_MOVE_RB_APM_SEARCH_BASE_HPP

#include <cstdint>
#include <vector>

#include <misc/apm.hpp>
#include <misc/cigar.hpp>

/**
 * @brief returns the direction in which each step of a search matches its part: RIGHT if its part lies right of the
 *        part of the previous step, else LEFT; the first step has the direction of the second one (RIGHT if p = 1)
 * @param S the search
 * @return [0..p-1] the direction of step i
 */
inline std::vector<direction_t> apm_step_directions(const search_t& S)
{
    std::vector<direction_t> dirs(S.size());
    dirs[0] = S.size() == 1 || S[1].part > S[0].part ? RIGHT : LEFT;
    for (uint64_t i = 1; i < S.size(); i++) dirs[i] = S[i].part > S[i - 1].part ? RIGHT : LEFT;
    return dirs;
}

/**
 * @brief the members of apm_hamming and apm_edit: the index, the query parameters, the locate query without errors
 *        and the test for the separator of the sequences
 * @tparam idx_t the bidirectional index type (move_rb, or a benchmark adapter with the same interface)
 */
template <typename idx_t>
class apm_search_base {
protected:
    using inp_t = typename idx_t::inp_t;

    template <query_support_t query_support>
    using search_context_t = typename idx_t::template search_context_t<query_support>;

    const idx_t& _index; // the index
    apm_params_t _params; // the query parameters

    apm_search_base(const idx_t& index, apm_params_t params) : _index(index), _params(params) {}

    /**
     * @brief returns whether the next symbol of an extension context is the separator of the sequences of the index
     *        (always false for an index without sequences); a search does not extend with it
     * @tparam ext_ctx_t type of the extension context
     * @param ext_ctx an extension context
     */
    template <typename ext_ctx_t>
    bool next_is_separator(const ext_ctx_t& ext_ctx) const
    {
        if constexpr (requires (const idx_t& ix) { ix.seq_data(); }) {
            return _index.seq_data().has_sequences() && ext_ctx.can_extend()
                && ext_ctx.next_symbol() == _index.seq_data().separator_sym();
        } else {
            return false;
        }
    }

    /**
     * @brief locates the occurrences of P without errors
     * @tparam mode whether each occurrence is reported with a CIGAR index
     * @tparam report_fnc_t type of the report function
     * @param P the pattern
     * @param report called as report({pos, m, 0}) or, in CIGAR mode, report({pos, m, 0, 0}) for each occurrence
     * @return in CIGAR mode, the CIGAR of the occurrences (m matches); empty in NO_CIGAR mode
     */
    template <cigar_mode_t mode, typename report_fnc_t>
    std::vector<cigar_t<uint8_t>> locate_exact(const inp_t& P, report_fnc_t& report) const
    {
        const uint64_t m = P.size();
        std::vector<cigar_t<uint8_t>> cigars;

        if constexpr (mode == CIGAR) {
            cigars.push_back(m > 0 ? cigar_t<uint8_t>{MATCH<uint8_t>(m)} : cigar_t<uint8_t>{});
            _index.forward_index().locate(P, [&](uint64_t, uint64_t occ){report({occ, m, 0, uint32_t(0)});});
        } else {
            _index.forward_index().locate(P, [&](uint64_t, uint64_t occ){report({occ, m, 0});});
        }

        return cigars;
    }
};

#endif // MOVE_R_ALGORITHMS_MOVE_RB_APM_SEARCH_BASE_HPP
