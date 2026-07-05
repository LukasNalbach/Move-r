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

#include <functional>
#include <cstdint>
#include <vector>
#include <type_traits>

#include <misc/utils.hpp>

/**
 * @brief computes the length of the longest common prefix of the ranges [beg_1, end_1) and [beg_2, end_2)
 * @tparam val_t type of the values in the ranges
 * @param beg_1 pointer to the first element of the first range
 * @param end_1 pointer past the last element of the first range
 * @param beg_2 pointer to the first element of the second range
 * @param end_2 pointer past the last element of the second range
 * @return the length (number of elements) of the longest common prefix
 */
template <typename val_t>
inline static uint64_t lce(const val_t* beg_1, const val_t* end_1, const val_t* beg_2, const val_t* end_2)
{
    static constexpr uint64_t blk_size = sizeof(uint64_t) / sizeof(val_t);
    const uint64_t max_lce = std::min<uint64_t>(end_1 - beg_1, end_2 - beg_2);
    const uint64_t max_blks = max_lce / blk_size;
    uint64_t const* const blks_1 = reinterpret_cast<uint64_t const*>(beg_1);
    uint64_t const* const blks_2 = reinterpret_cast<uint64_t const*>(beg_2);
    uint64_t lce_val = 0;

    while (lce_val < max_blks && blks_1[lce_val] == blks_2[lce_val]) {
        lce_val++;
    }

    if (lce_val == max_blks) [[unlikely]] {
        lce_val *= blk_size;

        while (lce_val < max_lce && beg_1[lce_val] == beg_2[lce_val]) {
            lce_val++;
        }

        return lce_val;
    }

    return lce_val * blk_size + std::countr_zero(blks_1[lce_val] ^ blks_2[lce_val]) / (8 * sizeof(val_t));
}

/**
 * @brief computes the length of the longest common prefix of the suffix of input starting at input_pos and pattern
 * @tparam inp_t container type of input and pattern
 * @param input the input
 * @param pattern the pattern
 * @param input_pos position in the input at which the suffix starts
 * @param lce_min a known lower bound on the longest common prefix (default: 0)
 * @return the length of the longest common prefix
 */
template <typename inp_t>
inline uint64_t lce(
    const inp_t& input, const inp_t& pattern,
    uint64_t input_pos, uint64_t lce_min = 0
) {
    return lce_min + lce<typename inp_t::value_type>(
        input.data() + input_pos + lce_min, input.data() + input.size(),
        pattern.data() + lce_min, pattern.data() + pattern.size());
}

/**
 * @brief returns whether the suffix of input starting at input_pos is lexicographically smaller than pattern
 * @tparam inp_t container type of input and pattern
 * @param input the input
 * @param pattern the pattern
 * @param input_pos position in the input at which the suffix starts
 * @param lce the length of the longest common prefix of the suffix and the pattern
 * @return whether the suffix is lexicographically smaller than the pattern
 */
template <typename inp_t>
inline bool lex_less(
    const inp_t& input, const inp_t& pattern,
    uint64_t input_pos, uint64_t lce
) {
    if (lce == pattern.size()) [[unlikely]] return false;
    if (input_pos + lce == input.size()) [[unlikely]] return true;
    using cmp_t = std::make_unsigned_t<typename inp_t::value_type>;
    return *reinterpret_cast<const cmp_t*>(&input[input_pos + lce]) <
           *reinterpret_cast<const cmp_t*>(&pattern[lce]);
}

/**
 * @brief finds, by binary search over the sampled suffixes in [beg, end], the interval of suffixes that have pattern as a prefix
 * @tparam inp_t container type of input and pattern
 * @tparam fnc_t type of the get_sa_sample function
 * @param input the input
 * @param pattern the pattern
 * @param beg index of the first sampled suffix to consider
 * @param end index of the last sampled suffix to consider
 * @param get_sa_sample function mapping a sample index to the starting position in the input of the corresponding suffix
 * @return the interval [b, e] of sample indices whose suffixes have pattern as a prefix (b > e if the pattern does not occur)
 */
template <typename inp_t, typename fnc_t>
std::tuple<uint64_t, uint64_t> sa_interval(
    const inp_t& input, const inp_t& pattern,
    uint64_t beg, uint64_t end,
    fnc_t get_sa_sample
) {
    uint64_t len = pattern.size();

    uint64_t l = beg;
    uint64_t r = end;
    uint64_t m;

    uint64_t lce_l = lce(input, pattern, get_sa_sample(l));
    uint64_t lce_r = lce(input, pattern, get_sa_sample(r));

    uint64_t l_ = l;
    uint64_t r_ = r;
    uint64_t lce_l_ = lce_l;
    uint64_t lce_r_ = lce_r;

    uint64_t lce_m;
    uint64_t pos_m;

    while (r - l > 1) {
        m = l + (r - l) / 2;
        pos_m = get_sa_sample(m);
        lce_m = lce(input, pattern, pos_m, std::min(lce_l, lce_r));

        if (lce_m >= len) {
            r = m;
            lce_r = lce_m;

            if (m > l_) {
                l_ = m;
                lce_l_ = lce_m;
            }
        } else if (lex_less(input, pattern, pos_m, lce_m)) {
            l = m;
            lce_l = lce_m;

            if (m > l_) {
                l_ = m;
                lce_l_ = lce_m;
            }
        } else {
            r = m;
            lce_r = lce_m;

            if (m < r_) {
                r_ = m;
                lce_r_ = lce_m;
            }
        }
    }

    uint64_t b;

    if (lce_l < len) {
        if (lce_r < len) {
            return {l + 1, l};
        }

        b = r;
    } else {
        b = l;
    }

    l = l_;
    r = r_;
    lce_l = lce_l_;
    lce_r = lce_r_;

    while (r - l > 1) {
        m = l + (r - l) / 2;
        lce_m = lce(input, pattern, get_sa_sample(m), std::min(lce_l, lce_r));

        if (lce_m >= len) {
            l = m;
            lce_l = lce_m;
        } else {
            r = m;
            lce_r = lce_m;
        }
    }

    uint64_t e;

    if (lce_r >= len) {
        e = r;
    } else {
        e = l;
    }

    return {b, e};
}

/**
 * @brief which lexicographic boundary lex_boundary() should search for
 */
enum lex_boundary_t {
    _min_lex_geq, // search for the minimum lexicographically greater than or equal suffix
    _max_lex_leq // search for the maximum lexicographically less than or equal suffix
};

/**
 * @brief finds, by binary search over the sampled suffixes in [beg, end], the lexicographic boundary of pattern
 * @tparam boundary_type whether to search for the minimum lex. greater-or-equal or maximum lex. less-or-equal boundary
 * @tparam inp_t container type of input and pattern
 * @param input the input
 * @param pattern the pattern
 * @param beg index of the first sampled suffix to consider
 * @param end index of the last sampled suffix to consider
 * @param get_sa_sample function mapping a sample index to the starting position in the input of the corresponding suffix
 * @return the sample index of the requested lexicographic boundary
 */
template <lex_boundary_t boundary_type, typename inp_t>
uint64_t lex_boundary(
    const inp_t& input, const inp_t& pattern,
    uint64_t beg, uint64_t end,
    const std::function<uint64_t(uint64_t)>& get_sa_sample
) {
    uint64_t len = pattern.size();

    uint64_t l = beg;
    uint64_t r = end;
    uint64_t m;

    uint64_t lce_l = lce(input, pattern, get_sa_sample(l));
    uint64_t lce_r = lce(input, pattern, get_sa_sample(r));
    uint64_t lce_m;

    while (r - l > 1) {
        m = l + (r - l) / 2;
        lce_m = lce(input, pattern, get_sa_sample(m), std::min(lce_l, lce_r));

        if ((lce_m >= len) == (boundary_type == _min_lex_geq)) {
            r = m;
            lce_r = lce_m;
        } else {
            l = m;
            lce_l = lce_m;
        }
    }

    if constexpr (boundary_type == _min_lex_geq) {
        return lce_l >= len ? l : r;
    } else {
        return lce_r >= len ? r : l;
    }
}

/**
 * @brief locates the SA-interval of pattern using binary search over sampled SA-positions and, optionally,
 *        extracts and reports the SA-values in that interval
 * @tparam out_t type of the extracted SA-values
 * @tparam inp_t container type of input and pattern
 * @tparam smpl_fnc_t type of the get_sa_sample function
 * @tparam smpl_pos_fnc_t type of the get_sa_smpl_pos function
 * @tparam rng_report_fnc_t type of the report_sa_rng function
 * @tparam report_fnc_t type of the report function
 * @param input the input
 * @param pattern the pattern
 * @param num_samples the number of sampled SA-positions
 * @param get_sa_sample function mapping a sample index to the starting position in the input of the corresponding suffix
 * @param get_sa_smpl_pos function mapping a sample index to its position in the (full) suffix array
 * @param report_sa_rng function that reports the SA-values in a given range of SA-positions
 * @param report function that reports a single extracted SA-value
 * @param reverse_rng whether the reported SA-ranges have to be reversed
 * @param extract_sa_values whether to extract and report the SA-values in the located interval
 * @return the located SA-interval [beg, end]
 */
template <typename out_t, typename inp_t, typename smpl_fnc_t, typename smpl_pos_fnc_t, typename rng_report_fnc_t, typename report_fnc_t>
std::tuple<uint64_t, uint64_t> binary_sa_search_and_extract(
    const inp_t& input, const inp_t& pattern, uint64_t num_samples,
    smpl_fnc_t get_sa_sample, smpl_pos_fnc_t get_sa_smpl_pos,
    rng_report_fnc_t report_sa_rng, report_fnc_t report,
    bool reverse_rng, bool extract_sa_values
) {
    // binary search over sampled SA-positions
    auto [fst_smpl_idx, lst_smpl_idx] = sa_interval(input, pattern, 0, num_samples - 1, get_sa_sample);

    if (fst_smpl_idx > lst_smpl_idx || lst_smpl_idx - fst_smpl_idx <= 2) {
        // the unsampled regions surrounding the (lexicographical) boundary-samples are connected

        if (fst_smpl_idx > 0) fst_smpl_idx--;
        if (lst_smpl_idx < num_samples - 1) lst_smpl_idx++;

        uint64_t fst_smpl_pos = get_sa_smpl_pos(fst_smpl_idx);
        uint64_t lst_smpl_pos = get_sa_smpl_pos(lst_smpl_idx);

        std::vector<out_t> sa_rng;
        sa_rng.reserve(lst_smpl_pos - fst_smpl_pos + 1);
        report_sa_rng(fst_smpl_pos, lst_smpl_pos,
            get_sa_sample(fst_smpl_idx), get_sa_sample(lst_smpl_idx),
            [&](uint64_t val){sa_rng.emplace_back(val);});
        if (reverse_rng) std::reverse(sa_rng.begin(), sa_rng.end());

        auto [beg_offs, end_offs] = sa_interval(
            input, pattern, 0, lst_smpl_pos - fst_smpl_pos,
            [&](uint64_t i){return sa_rng[i];});

        uint64_t beg = fst_smpl_pos + beg_offs;
        uint64_t end = fst_smpl_pos + end_offs;

        if (extract_sa_values && beg <= end) {
            // report the SA-values of the located interval; sa_rng[i] holds the SA-value of
            // SA-position fst_smpl_pos + i, so the interval [beg, end] maps to [beg_offs, end_offs]
            for (uint64_t i = beg_offs; i <= end_offs; i++) {
                report(sa_rng[i]);
            }
        }

        return {beg, end};
    } else {
        // the unsampled regions surrounding the (lexicographical) boundary-samples are not connected

        uint64_t prev_fst_smpl_idx = std::max<uint64_t>(fst_smpl_idx - 1, 0);
        uint64_t next_fst_smpl_idx = std::min<uint64_t>(fst_smpl_idx + 1, num_samples - 1);
        uint64_t prev_lst_smpl_idx = std::max<uint64_t>(lst_smpl_idx - 1, 0);
        uint64_t next_lst_smpl_idx = std::min<uint64_t>(lst_smpl_idx + 1, num_samples - 1);

        uint64_t prev_fst_smpl_pos = get_sa_smpl_pos(prev_fst_smpl_idx);
        uint64_t next_fst_smpl_pos = get_sa_smpl_pos(next_fst_smpl_idx);
        uint64_t prev_lst_smpl_pos = get_sa_smpl_pos(prev_lst_smpl_idx);
        uint64_t next_lst_smpl_pos = get_sa_smpl_pos(next_lst_smpl_idx);

        // extract the left surrounding region
        std::vector<out_t> sa_rng_left;
        sa_rng_left.reserve(next_fst_smpl_pos - prev_fst_smpl_pos + 1);
        report_sa_rng(prev_fst_smpl_pos, next_fst_smpl_pos,
            get_sa_sample(prev_fst_smpl_idx), get_sa_sample(next_fst_smpl_idx),
            [&](uint64_t val){sa_rng_left.emplace_back(val);});
        if (reverse_rng) std::reverse(sa_rng_left.begin(), sa_rng_left.end());

        std::vector<out_t> sa_rng_right;
        sa_rng_right.reserve(next_lst_smpl_pos - prev_lst_smpl_pos + 1);
        report_sa_rng(prev_lst_smpl_pos, next_lst_smpl_pos,
            get_sa_sample(prev_lst_smpl_idx), get_sa_sample(next_lst_smpl_idx),
            [&](uint64_t val){sa_rng_right.emplace_back(val);});
        if (reverse_rng) std::reverse(sa_rng_right.begin(), sa_rng_right.end());

        uint64_t beg_offs = lex_boundary<_min_lex_geq>(
            input, pattern, 0, next_fst_smpl_pos - prev_fst_smpl_pos,
            [&](uint64_t i){return sa_rng_left[i];});

        uint64_t end_offs = lex_boundary<_max_lex_leq>(
            input, pattern, 0, next_lst_smpl_pos - prev_lst_smpl_pos,
            [&](uint64_t i){return sa_rng_right[i];});

        uint64_t beg = prev_fst_smpl_pos + beg_offs;
        uint64_t end = prev_lst_smpl_pos + end_offs;

        if (extract_sa_values && beg <= end) {
            uint64_t left_len = next_fst_smpl_pos - beg;
            uint64_t left_offs = beg - prev_fst_smpl_pos;
            uint64_t right_len = end - prev_lst_smpl_pos;

            for (uint64_t i = 0; i < left_len; i++) {
                report(sa_rng_left[left_offs + i]);
            }

            report_sa_rng(next_fst_smpl_pos, prev_lst_smpl_pos,
                get_sa_sample(next_fst_smpl_idx), get_sa_sample(prev_lst_smpl_idx),
                [&](uint64_t val){report(val);});

            for (uint64_t i = 0; i < right_len; i++) {
                report(sa_rng_right[1 + i]);
            }
        }

        return {beg, end};
    }
}