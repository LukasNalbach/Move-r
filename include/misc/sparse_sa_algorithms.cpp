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
    uint64_t const* const blks_2 = reinterpret_cast<uint64_t const*>(end_1);
    uint64_t lce_val = 0;

    while (lce_val < max_blks && blks_1[lce_val] == blks_2[lce_val]) {
        lce_val++;
    }

    if (lce_val == max_blks) [[unlikely]] {
        lce_val *= blk_size;

        while (lce_val < max_lce && beg_1[lce_val] == end_1[lce_val]) {
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

    uint64_t l_ = l;
    uint64_t r_ = r;
    uint64_t lce_l_;
    uint64_t lce_r_;

    uint64_t lce_l = lce(input, pattern, get_sa_sample(l));
    uint64_t lce_r = lce(input, pattern, get_sa_sample(r));
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
            return {1, 0};
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
    _min_lex_geq,
    _max_lex_leq
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

        if (lce_m >= len) {
            if constexpr (boundary_type == _min_lex_geq) {
                l = m;
                lce_l = lce_m;
            } else {
                r = m;
                lce_r = lce_m;
            }
        } else {
            if constexpr (boundary_type == _min_lex_geq) {
                r = m;
                lce_r = lce_m;
            } else {
                l = m;
                lce_l = lce_m;
            }
        }
    }

    if constexpr (boundary_type == _min_lex_geq) {
        return lce_r >= len ? r : l;
    } else {
        return lce_r >= len ? l : r;
    }
}

template <typename out_t>
std::tuple<uint64_t, uint64_t, std::vector<out_t>> binary_sa_search_and_extract(
    const std::string& input, const std::string& pattern, uint64_t num_samples,
    std::function<uint64_t(uint64_t)>&& get_sa_sample,
    std::function<uint64_t(uint64_t)>&& get_sa_sample_pos,
    std::function<void(uint64_t, uint64_t, std::function<void(uint64_t, uint64_t)>&&)>&& report_sa_range,
    bool extract_sa_values = true
) {
    // binary search over sampled SA-positions
    auto [sparse_beg, sparse_end] = sa_interval(input, pattern, 0, num_samples - 1, get_sa_sample);

    if (sparse_end - sparse_beg <= 2) {
        // the unsampled regions surrounding the (lexicographical) boundary-samples are connected

        uint64_t new_beg = std::max<int64_t>(int64_t{sparse_beg} - 1, 0);
        uint64_t new_end = std::min<int64_t>(int64_t{sparse_end} + 1, num_samples - 1);

        uint64_t rough_beg = get_sa_sample_pos(new_beg);
        uint64_t rough_end = get_sa_sample_pos(new_end);

        std::vector<out_t> sa_range;
        no_init_resize(sa_range, rough_end - rough_beg + 1);
        report_sa_range(rough_beg, rough_end,
            [&](uint64_t pos, uint64_t val){sa_range[pos - rough_beg] = val;});

        auto [exact_beg_offset, exact_end_offset] = sa_interval(
            input, pattern, 0, rough_end - rough_beg - 1,
            [&](uint64_t i){return sa_range[i];});

        uint64_t exact_beg = rough_beg + exact_beg_offset;
        uint64_t exact_end = rough_beg + exact_end_offset;

        if (!extract_sa_values) {
            return {exact_beg, exact_end, {}};
        } else {
            uint64_t len_interval = exact_end - exact_beg + 1;

            if (exact_beg != rough_beg) {
                uint64_t shift_left = exact_beg - rough_beg;

                for (uint64_t i = 0; i < len_interval; i++) {
                    sa_range[i] = sa_range[i + shift_left];
                }
            }

            sa_range.resize(len_interval);
            return {exact_beg, exact_end, sa_range};
        }
    } else {
        // the unsampled regions surrounding the (lexicographical) boundary-samples are not connected

        uint64_t new_sparse_left_beg = std::max<int64_t>(int64_t{sparse_beg} - 1, 0);
        uint64_t new_sparse_left_end = std::min<int64_t>(int64_t{sparse_beg} + 1, num_samples - 1);
        uint64_t new_sparse_right_beg = std::max<int64_t>(int64_t{sparse_end} - 1, 0);
        uint64_t new_sparse_right_end = std::min<int64_t>(int64_t{sparse_end} + 1, num_samples - 1);

        uint64_t rough_left_beg = get_sa_sample_pos(new_sparse_left_beg);
        uint64_t rough_left_end = get_sa_sample_pos(new_sparse_left_end);
        uint64_t rough_right_beg = get_sa_sample_pos(new_sparse_right_beg);
        uint64_t rough_right_end = get_sa_sample_pos(new_sparse_right_end);

        // extract the left surrounding region
        std::vector<out_t> sa_range_left;
        no_init_resize(sa_range_left, rough_left_end - rough_left_beg + 1);
        report_sa_range(rough_left_beg, rough_left_end,
            [&](uint64_t pos, uint64_t val){sa_range_left[pos] = val;});

        std::vector<out_t> sa_range_right;
        no_init_resize(sa_range_right, rough_right_end - rough_right_beg + 1);
        report_sa_range(rough_right_beg, rough_right_end,
            [&](uint64_t pos, uint64_t val){sa_range_right[pos] = val;});

        uint64_t exact_beg_offset = lex_boundary<_min_lex_geq>(
            input, pattern, 0, rough_left_end - rough_left_beg - 1,
            [&](uint64_t i){return sa_range_left[i];});

        uint64_t exact_end_offset = lex_boundary<_max_lex_leq>(
            input, pattern, 0, rough_right_end - rough_right_beg - 1,
            [&](uint64_t i){return sa_range_right[i];});

        uint64_t exact_beg = rough_left_beg + exact_beg_offset;
        uint64_t exact_end = rough_right_beg + exact_end_offset;

        if (!extract_sa_values) {
            return {exact_beg, exact_end, {}};
        } else {
            uint64_t len_interval = exact_end - exact_beg + 1;
            std::vector<out_t> result;
            no_init_resize(result, len_interval);

            uint64_t left_range_length = rough_left_end - exact_beg + 1;
            uint64_t left_range_offset = exact_beg - rough_left_beg;
            
            for (uint64_t i = 0; i < left_range_length; i++) {
                result[i] = sa_range_left[left_range_offset + i];
            }

            uint64_t right_range_length = exact_end - rough_right_beg + 1;
            uint64_t right_range_offset = rough_right_beg - exact_beg;
            
            for (uint64_t i = 0; i < right_range_length; i++) {
                result[right_range_offset + i] = sa_range_right[i];
            }

            report_sa_range(rough_left_end + 1, rough_right_beg - 1,
                [&](uint64_t pos, uint64_t val){result[pos - exact_beg] = val;});

            return {exact_beg, exact_end, result};
        }
    }
}

template <typename out_t>
std::vector<out_t> extract_range_using_samples(
    uint64_t beg, uint64_t end, uint64_t num_samples,
    std::function<uint64_t(uint64_t)>&& get_sample,
    std::function<uint64_t(uint64_t)>&& get_sample_pos,
    std::function<std::vector<out_t>(int64_t, int64_t)>&& extract_differential_range
) {
    int64_t interval_length = end - beg + 1;
    int64_t sample_index = bin_search_max_leq<uint64_t>(end, 0, num_samples - 1, [&](uint64_t x) { return get_sample_pos(x); });
    int64_t sample_pos = get_sample_pos(sample_index);
    std::vector<out_t> result;
    
    if (sample_pos < beg) {
        int64_t next_sample_pos = get_sample_pos(sample_index + 1);

        if (next_sample_pos <= end || next_sample_pos - end < beg - sample_pos) {
            sample_index++;
            sample_pos = next_sample_pos;
        }
    }

    if (sample_pos <= beg) {
        result = extract_differential_range(sample_pos, end);
        int64_t val = get_sample(sample_index);
        int64_t dist = beg - sample_pos;

        for (int64_t i = 1; i <= dist; i++) {
            val += result[i];
        }

        int64_t last_val;

        for (int64_t i = 0; i < interval_length; i++) {
            last_val = val;
            val += result[i + dist + 1];
            result[i] = last_val;
        }

        result.resize(interval_length);
    } else if (sample_pos <= end) {
        result = extract_differential_range(beg, end);
        int64_t sample = get_sample(sample_index);
        int64_t val = sample;

        for (int64_t i = sample_pos - beg + 1; i < interval_length; i++) {
            val += result[i];
            result[i] = val;
        }

        val = sample;
        int64_t last_val;

        for (int64_t i = sample_pos - beg; i >= 0; i--) {
            last_val = val;
            val -= result[i];
            result[i] = last_val;
        }
    } else {
        result = extract_differential_range(beg, sample_pos);
        int64_t val = get_sample(sample_index);

        for (int64_t i = result.size() - 1; i >= interval_length; i--) {
            val -= result[i];
        }

        int64_t last_val;

        for (int64_t i = interval_length - 1; i >= 1; i--) {
            last_val = val;
            val -= result[i];
            result[i] = last_val;
        }

        result[0] = val;
        result.resize(interval_length);
    }

    return result;
}