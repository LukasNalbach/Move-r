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

#include <misc/utils.hpp>

template <typename char_t>
inline static uint64_t lce(const char_t* beg_1, const char_t* end_1, const char_t* beg_2, const char_t* end_2)
{
    static constexpr uint64_t blk_size = sizeof(uint64_t) / sizeof(char_t);
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

    return lce_val * blk_size + std::countr_zero(blks_1[lce_val] ^ blks_2[lce_val]) / (8 * sizeof(char_t));
}

inline uint64_t lce(
    const std::string& input, const std::string& pattern,
    uint64_t input_pos, uint64_t lce_min = 0
) {
    return lce_min + lce<char>(
        input.data() + input_pos + lce_min, input.data() + input.size(),
        pattern.data() + lce_min, pattern.data() + pattern.size());
}

inline bool lex_less(
    const std::string& input, const std::string& pattern,
    uint64_t input_pos, uint64_t lce
) {
    if (input_pos + lce == input.size() || lce == pattern.size()) [[unlikely]] return false;
    return char_to_uchar(input[input_pos + lce]) < char_to_uchar(pattern[lce]);
}

std::tuple<uint64_t, uint64_t> sa_interval(
    const std::string& input, const std::string& pattern,
    uint64_t beg, uint64_t end,
    const std::function<uint64_t(uint64_t)>& get_sa_sample
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

enum lex_boundary_t {
    _min_lex_geq, // search for the minimum lexicographically greater than or equal suffix
    _max_lex_leq // search for the maximum lexicographically less than or equal suffix
};

template <lex_boundary_t boundary_type>
uint64_t lex_boundary(
    const std::string& input, const std::string& pattern,
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

template <typename out_t>
std::tuple<uint64_t, uint64_t, std::vector<out_t>> binary_sa_search_and_extract(
    const std::string& input, const std::string& pattern, uint64_t num_samples,
    std::function<uint64_t(uint64_t)>&& get_sa_sample,
    std::function<uint64_t(uint64_t)>&& get_sa_smpl_pos,
    std::function<void(uint64_t, uint64_t, uint64_t, uint64_t, std::function<void(uint64_t, uint64_t)>&&)>&& report_sa_rng,
    bool extract_sa_values = true
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
        no_init_resize(sa_rng, lst_smpl_pos - fst_smpl_pos + 1);
        report_sa_rng(fst_smpl_pos, lst_smpl_pos,
            get_sa_sample(fst_smpl_idx), get_sa_sample(lst_smpl_idx),
            [&](uint64_t pos, uint64_t val){sa_rng[pos - fst_smpl_pos] = val;});

        auto [beg_offs, end_offs] = sa_interval(
            input, pattern, 0, lst_smpl_pos - fst_smpl_pos,
            [&](uint64_t i){return sa_rng[i];});

        uint64_t beg = fst_smpl_pos + beg_offs;
        uint64_t end = fst_smpl_pos + end_offs;

        if (!extract_sa_values) {
            return {beg, end, {}};
        } else {
            uint64_t len = end - beg + 1;

            if (beg != fst_smpl_pos) {
                uint64_t shft_left = beg - fst_smpl_pos;

                for (uint64_t i = 0; i < len; i++) {
                    sa_rng[i] = sa_rng[i + shft_left];
                }
            }

            sa_rng.resize(len);
            return {beg, end, sa_rng};
        }
    } else {
        // the unsampled regions surrounding the (lexicographical) boundary-samples are not connected

        uint64_t prev_fst_smpl_idx = std::max<int64_t>(int64_t{fst_smpl_idx} - 1, 0);
        uint64_t next_fst_smpl_idx = std::min<int64_t>(int64_t{fst_smpl_idx} + 1, num_samples - 1);
        uint64_t prev_lst_smpl_idx = std::max<int64_t>(int64_t{lst_smpl_idx} - 1, 0);
        uint64_t next_lst_smpl_idx = std::min<int64_t>(int64_t{lst_smpl_idx} + 1, num_samples - 1);

        uint64_t prev_fst_smpl_pos = get_sa_smpl_pos(prev_fst_smpl_idx);
        uint64_t next_fst_smpl_pos = get_sa_smpl_pos(next_fst_smpl_idx);
        uint64_t prev_lst_smpl_pos = get_sa_smpl_pos(prev_lst_smpl_idx);
        uint64_t next_lst_smpl_pos = get_sa_smpl_pos(next_lst_smpl_idx);

        // extract the left surrounding region
        std::vector<out_t> sa_rng_left;
        no_init_resize(sa_rng_left, next_fst_smpl_pos - prev_fst_smpl_pos + 1);
        report_sa_rng(prev_fst_smpl_pos, next_fst_smpl_pos,
            get_sa_sample(prev_fst_smpl_idx), get_sa_sample(next_fst_smpl_idx),
            [&](uint64_t pos, uint64_t val){sa_rng_left[pos - prev_fst_smpl_pos] = val;});

        std::vector<out_t> sa_rng_right;
        no_init_resize(sa_rng_right, next_lst_smpl_pos - prev_lst_smpl_pos + 1);
        report_sa_rng(prev_lst_smpl_pos, next_lst_smpl_pos,
            get_sa_sample(prev_lst_smpl_idx), get_sa_sample(next_lst_smpl_idx),
            [&](uint64_t pos, uint64_t val){sa_rng_right[pos - prev_lst_smpl_pos] = val;});

        uint64_t beg_offs = lex_boundary<_min_lex_geq>(
            input, pattern, 0, next_fst_smpl_pos - prev_fst_smpl_pos,
            [&](uint64_t i){return sa_rng_left[i];});

        uint64_t end_offs = lex_boundary<_max_lex_leq>(
            input, pattern, 0, next_lst_smpl_pos - prev_lst_smpl_pos,
            [&](uint64_t i){return sa_rng_right[i];});

        uint64_t beg = prev_fst_smpl_pos + beg_offs;
        uint64_t end = prev_lst_smpl_pos + end_offs;

        if (!extract_sa_values) {
            return {beg, end, {}};
        } else {
            uint64_t len = end - beg + 1;
            std::vector<out_t> sa_rng;
            no_init_resize(sa_rng, len);

            uint64_t left_len = next_fst_smpl_pos - beg + 1;
            uint64_t left_offs = beg - prev_fst_smpl_pos;
            
            for (uint64_t i = 0; i < left_len; i++) {
                sa_rng[i] = sa_rng_left[left_offs + i];
            }

            uint64_t right_len = end - prev_lst_smpl_pos + 1;
            uint64_t right_offs = prev_lst_smpl_pos - beg;
            
            for (uint64_t i = 0; i < right_len; i++) {
                sa_rng[right_offs + i] = sa_rng_right[i];
            }

            report_sa_rng(next_fst_smpl_pos, prev_lst_smpl_pos,
                get_sa_sample(next_fst_smpl_idx), get_sa_sample(prev_lst_smpl_idx),
                [&](uint64_t pos, uint64_t val){sa_rng[pos - beg] = val;});

            return {beg, end, sa_rng};
        }
    }
}