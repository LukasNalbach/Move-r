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

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <tuple>
#include <vector>
#include <cmath>
#include <sstream>
#include <string>

#include <data_structures/interleaved_byte_aligned_vectors.hpp>

template <typename word_t>
class edit_distance_matrix
{
    static_assert(std::is_same_v<word_t, uint64_t> || std::is_same_v<word_t, __uint128_t>);

  public:
    static constexpr uint16_t word_size = sizeof(word_t) * 8;
    static constexpr uint16_t blk_size = word_size / 2;
    static constexpr uint16_t k_limit = (word_size - blk_size - 2) / 3;
    static constexpr uint16_t max_first_col_size = 2 * k_limit + 1;
    static constexpr uint16_t two_k_limit = 2 * k_limit;
    
  protected:
    uint16_t sigma;
    uint16_t k_max;
    uint16_t m;
    uint16_t n;
    uint16_t band_width;
    uint16_t band_height;

    const std::vector<uint8_t>* map_int;

    struct bit_vectors_t {
        word_t HP;
        word_t HN;
        word_t D0;
        word_t RAC;
        uint16_t score;
    };

    std::vector<bit_vectors_t> bit_vectors;
    std::vector<word_t> match_vectors;

    uint16_t last_col_idx(uint16_t i) const
    {
        return std::min(n - 1, i + band_height);
    }

    uint8_t map_sym(char sym)
    {
        return (*map_int)[char_to_uchar(sym)] - 1;
    }

  public:
    edit_distance_matrix() {}

    template <typename pos_t>
    void set_input(const sub_string<pos_t>& input, uint16_t sigma, const std::vector<uint8_t>& map_int)
    {
        this->sigma = sigma;
        this->map_int = &map_int;
        n = input.size() + 1;
        m = 2 * k_limit + n;
        uint16_t match_vectors_size = (m + blk_size - 1) / blk_size;
        match_vectors.resize(match_vectors_size * sigma);
        word_t init = (word_t(1) << max_first_col_size) - word_t(1);
        std::fill(match_vectors.begin(), match_vectors.begin() + sigma, init);
        word_t bitmask = word_t(1) << max_first_col_size;
        uint16_t j_end = std::min<uint16_t>(input.size(), word_size - max_first_col_size);

        for (uint16_t j = 0; j < j_end; j++) {
            match_vectors[map_sym(input[j])] |= bitmask;
            bitmask <<= 1;
        }

        for (uint16_t b = 1; b < match_vectors_size; b++) {
            uint64_t beg_b = b * sigma;
            uint64_t beg_bm1 = beg_b - sigma;

            for (uint64_t i = 0; i < sigma; i++) {
                match_vectors[beg_b + i] = match_vectors[beg_bm1 + i] >> blk_size;
            }

            bitmask = word_t(1) << (word_size - blk_size);
            uint16_t j_beg = word_size - max_first_col_size + (b - 1) * blk_size;
            uint16_t j_end = std::min<uint16_t>(input.size(), j_beg + blk_size);

            for (uint16_t j = j_beg; j < j_end; j++) {
                match_vectors[beg_b + map_sym(input[j])] |= bitmask;
                bitmask <<= 1;
            }
        }
    }

    void init(uint16_t k_max, const std::vector<uint16_t>& init_dists = {})
    {
        assert(k_max <= k_limit);
        assert(is_initialized());
        assert(init_dists.empty() || init_dists.front() <= k_max);
        assert(init_dists.empty() || init_dists.back() <= k_max);

        this->k_max = k_max;
        band_width = (init_dists.empty()) ? k_max : init_dists.size() - 1 + k_max - init_dists.back();
        assert(band_width <= 2 * k_limit);
        m = band_width + n;
        bit_vectors.resize(m);
        bit_vectors[0].score = init_dists.empty() ? 0 : init_dists[0];
        band_height = k_max - bit_vectors[0].score;
        
        if (band_width + band_height + 1 > m) {
            m = band_width + band_height + 1;
            bit_vectors.resize(m);
        }

        bit_vectors[0].HP = (~word_t(0)) << max_first_col_size;
        bit_vectors[0].HN = ~bit_vectors[0].HP;
        uint16_t i_max = std::min<uint16_t>(init_dists.size(), max_first_col_size + 1);
        word_t xor_bit = word_t(1) << max_first_col_size;

        for (uint16_t i = 1; i < i_max; i++) {
            xor_bit >>= 1;

            if (init_dists[i] < init_dists[i - 1]) {
                bit_vectors[0].HP ^= xor_bit;
                bit_vectors[0].HN ^= xor_bit;
            } else if (init_dists[i] == init_dists[i - 1]) {
                bit_vectors[0].HN ^= xor_bit;
            }
        }

        bit_vectors[0].RAC = word_t(1) << (two_k_limit + band_height);
    }

    bool compute_row(uint16_t i, char Y)
    {
        assert(i > 0);
        assert(i < m);
        assert(is_initialized());

        const uint16_t b = i / blk_size;
        const uint16_t l = i % blk_size;
        word_t& HP = bit_vectors[i].HP;
        word_t& HN = bit_vectors[i].HN;
        word_t& D0 = bit_vectors[i].D0;
        word_t& RAC = bit_vectors[i].RAC;
        word_t& M = match_vectors[b * sigma + map_sym(Y)];
        HP = bit_vectors[i - 1].HP;
        HN = bit_vectors[i - 1].HN;
        RAC = bit_vectors[i - 1].RAC << 1;

        if (l == 0) {
            HP >>= blk_size;
            HN >>= blk_size;
            RAC >>= blk_size;
        }

        D0 = (((M & HP) + HP) ^ HP) | M | HN;
        word_t VP = HN | ~(D0 | HP);
        word_t VN = D0 & HP;
        HP = (VN << 1) | ~(D0 | (VP << 1));
        HN = (D0 & (VP << 1));
        const uint16_t diag_bit = l + two_k_limit;
        bit_vectors[i].score = bit_vectors[i - 1].score + (D0 & (word_t(1) << diag_bit) ? 0 : 1);

        if (!(D0 & RAC)) {
            uint16_t val = 1;

            while (val > 0) {
                if (HP & RAC) val--;
                if (HN & RAC) val++;
                if (RAC == (word_t(1) << (diag_bit - band_width))) return false;
                RAC >>= 1;
            }
        }

        return true;
    }

    bool is_in_final_column(const uint16_t i) const
    {
        return i >= num_rows() - last_col_size();
    }

    uint16_t operator()(uint16_t i, uint16_t j) const
    {
        assert(i < m);
        assert(j < n);
        const uint16_t bit = (i % blk_size) + two_k_limit;
        uint16_t b = (i > j) ? bit - (i - j) + 1 : bit + 1;
        uint16_t e = (i > j) ? bit + 1 : bit + (j - i) + 1;
        word_t mask = ((word_t(1) << (e - b)) - word_t(1)) << b;
        int neg = popcount(bit_vectors[i].HN & mask);
        int pos = popcount(bit_vectors[i].HP & mask);
        uint16_t score = bit_vectors[i].score;
        score += (i > j) ? (neg - pos) : (pos - neg);
        return score;
    }

    bool only_vertical_gaps_left(uint16_t i) const
    {
        assert(i < m);
        if (i + max_first_col_size < n) return false;
        const uint16_t blk = i / blk_size;
        const uint16_t blk_offs = i % blk_size;
        uint16_t bb = two_k_limit - band_width + blk_offs + 1;
        uint16_t be = two_k_limit + n - blk * blk_size;
        return (((~bit_vectors[i].HN >> bb) << bb) << (word_size - be)) == word_t(0);
    }

    uint16_t num_cols() const
    {
        return n;
    }

    uint16_t num_rows() const
    {
        return m;
    }

    bool is_initialized() const
    {
        return !match_vectors.empty();
    }

    uint16_t last_col_size() const
    {
        return band_height + band_width + 1;
    }
};