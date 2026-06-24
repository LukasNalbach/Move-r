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
 * @file edit_distance_matrix.hpp
 *
 * Banded, word-packed edit-distance matrix.
 *
 * Algorithmic background (the techniques below are published and not the
 * contribution of this file; only their realization here is original):
 *   - the bit-parallel edit-distance recurrence (the D0 / VP / VN / HP / HN
 *     update) is G. Myers, "A fast bit-vector algorithm for approximate string
 *     matching based on dynamic programming", J. ACM 46(3), 1999, with the
 *     refinements of H. Hyyrö.
 *   - restricting the computation to a diagonal band of width 2k+1 (so that one
 *     machine word suffices for up to k errors) is the approach used by the
 *     search-scheme aligners around Columba (L. Renders, L. Depuydt,
 *     J. Fostier, Ghent University).
 * The recurrence itself is fixed by the algorithm; the band layout, the storage
 * scheme and the accessor/pruning helpers are an original implementation written
 * for Move-r.
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

#include <misc/utils.hpp>

/**
 * @brief banded edit-distance dynamic-programming matrix between a fixed input and a query string,
 *        evaluated one query symbol (row) at a time; every row keeps only the 2k+1-wide diagonal band
 *        of the column, bit-packed into a single word using Myers' bit-parallel recurrence
 * @tparam word_t unsigned word type that bit-encodes one banded matrix row (uint64_t or __uint128_t)
 */
template <typename word_t>
class edit_distance_matrix
{
    static_assert(std::is_same_v<word_t, uint64_t> || std::is_same_v<word_t, __uint128_t>);

  public:
    // a row word stores the 2k+1-bit first column followed by two further k-bit expansions of the
    // band as it is shifted; that budget of 3k+... bits per word caps the number of errors at k_limit
    static constexpr uint16_t bits_per_word = sizeof(word_t) * 8;
    static constexpr uint16_t half_word = bits_per_word / 2;
    static constexpr uint16_t k_limit = (bits_per_word - half_word - 2) / 3;
    static constexpr uint16_t first_col_bits = 2 * k_limit + 1; // width budget of the band's first column
    static constexpr uint16_t diag0 = 2 * k_limit;              // bit position of the main diagonal in a block

  protected:
    uint16_t sigma;       // size of the effective alphabet
    uint16_t k_max;       // error budget of the current query
    uint16_t m;           // number of rows (band length)
    uint16_t n;           // number of columns (input length + 1)
    uint16_t band_width;  // horizontal extent of the band
    uint16_t band_height; // vertical extent of the band

    const std::vector<uint8_t>* byte_to_rank; // byte -> one-based symbol rank in the effective alphabet

    /** @brief one banded row: the zero-diagonal mask, the two horizontal-delta masks, the frontier bit
     *         marking the rightmost band cell still within the error budget, and the diagonal's distance */
    struct row_state_t {
        word_t diag_zero; // D0: positions whose diagonal delta is 0
        word_t dh_pos;    // HP: positions whose horizontal delta is +1
        word_t dh_neg;    // HN: positions whose horizontal delta is -1
        word_t frontier;  // single set bit at the rightmost band cell with value <= k_max
        uint16_t score;   // edit distance stored on the diagonal of this row
    };

    std::vector<row_state_t> rows;   // one entry per matrix row
    std::vector<word_t> eq_masks;    // [block][symbol] -> bitmask of input positions equal to that symbol

    /** @brief resizes all per-row storage to @p num rows */
    void resize_rows(uint16_t num)
    {
        rows.resize(num);
    }

    /**
     * @brief index of the last column inside the band on row @p i
     * @param i row index
     * @return the band's rightmost column on row i
     */
    uint16_t last_col_idx(uint16_t i) const
    {
        return std::min(n - 1, i + band_height);
    }

    /**
     * @brief maps a query symbol to its zero-based index into the eq_masks of a block
     * @tparam sym_t the (one-byte) symbol type
     * @param sym a query symbol
     * @return the zero-based symbol rank
     */
    template <typename sym_t>
    uint8_t map_sym(sym_t sym) const
    {
        return (*byte_to_rank)[sym_to_uchar(sym)] - 1;
    }

  public:
    edit_distance_matrix() {}

    /**
     * @brief fixes the (column-)input that query rows are matched against and precomputes its eq-masks
     * @tparam pos_t unsigned position type of the substring view
     * @tparam inp_t underlying sequence type of the substring view
     * @param input the input sequence
     * @param sigma the alphabet size
     * @param alph_map byte -> one-based symbol rank in the effective alphabet
     */
    template <typename pos_t, typename inp_t>
    void set_input(const sub_string<pos_t, inp_t>& input, uint16_t sigma, const std::vector<uint8_t>& alph_map)
    {
        this->sigma = sigma;
        this->byte_to_rank = &alph_map;
        n = input.size() + 1;
        m = 2 * k_limit + n;
        uint16_t num_blocks = (m + half_word - 1) / half_word;
        eq_masks.resize(num_blocks * sigma);

        // block 0 spans the first-column budget plus the input symbols that fit above it; the low
        // first_col_bits bits are set for every symbol (the first column always "matches")
        word_t low_ones = (word_t(1) << first_col_bits) - word_t(1);
        std::fill(eq_masks.begin(), eq_masks.begin() + sigma, low_ones);
        word_t bit = word_t(1) << first_col_bits;
        uint16_t j_end = std::min<uint16_t>(input.size(), bits_per_word - first_col_bits);

        for (uint16_t j = 0; j < j_end; j++) {
            uint8_t rank = alph_map[sym_to_uchar(input[j])];
            if (rank != 0) eq_masks[rank - 1] |= bit;
            bit <<= 1;
        }

        // every later block is the previous block shifted down by one half-word, with the input
        // symbols that newly enter the window folded into the top half-word
        for (uint16_t blk = 1; blk < num_blocks; blk++) {
            uint64_t base = uint64_t(blk) * sigma;
            uint64_t prev = base - sigma;

            for (uint64_t c = 0; c < sigma; c++) {
                eq_masks[base + c] = eq_masks[prev + c] >> half_word;
            }

            bit = word_t(1) << (bits_per_word - half_word);
            uint16_t j_beg = bits_per_word - first_col_bits + (blk - 1) * half_word;
            uint16_t blk_end = std::min<uint16_t>(input.size(), j_beg + half_word);

            for (uint16_t j = j_beg; j < blk_end; j++) {
                uint8_t rank = alph_map[sym_to_uchar(input[j])];
                if (rank != 0) eq_masks[base + rank - 1] |= bit;
                bit <<= 1;
            }
        }
    }

    /**
     * @brief starts a fresh query allowing at most @p k_max errors
     * @param k_max the error budget
     * @param init_dists optional first-column distances when continuing a previous (adjacent) band
     */
    void init(uint16_t k_max, const std::vector<uint16_t>& init_dists = {})
    {
        assert(k_max <= k_limit);
        assert(is_initialized());
        assert(init_dists.empty() || init_dists.front() <= k_max);
        assert(init_dists.empty() || init_dists.back() <= k_max);

        this->k_max = k_max;
        band_width = init_dists.empty() ? k_max : init_dists.size() - 1 + k_max - init_dists.back();
        assert(band_width <= 2 * k_limit);
        m = band_width + n;
        resize_rows(m);
        rows[0].score = init_dists.empty() ? 0 : init_dists[0];
        band_height = k_max - rows[0].score;

        if (band_width + band_height + 1 > m) {
            m = band_width + band_height + 1;
            resize_rows(m);
        }

        // first column: by default the distances are 0, 1, 2, ..., i.e. every horizontal step is +1
        rows[0].dh_pos = (~word_t(0)) << first_col_bits;
        rows[0].dh_neg = ~rows[0].dh_pos;

        // when continuing a previous band, replay the given first-column distances as horizontal deltas
        uint16_t i_max = std::min<uint16_t>(init_dists.size(), first_col_bits + 1);
        word_t bit = word_t(1) << first_col_bits;

        for (uint16_t i = 1; i < i_max; i++) {
            bit >>= 1;

            if (init_dists[i] < init_dists[i - 1]) {
                rows[0].dh_pos ^= bit;
                rows[0].dh_neg ^= bit;
            } else if (init_dists[i] == init_dists[i - 1]) {
                rows[0].dh_neg ^= bit;
            }
        }

        rows[0].frontier = word_t(1) << (diag0 + band_height);
    }

    /**
     * @brief derives row @p i from row i-1 for the next query symbol @p Y
     * @tparam sym_t the (one-byte) symbol type
     * @param i the row to compute (i > 0)
     * @param Y the next query symbol
     * @return false if the whole band on row i exceeds k_max (the search may be pruned), true otherwise
     */
    template <typename sym_t>
    bool compute_row(uint16_t i, sym_t Y)
    {
        assert(i > 0);
        assert(i < m);
        assert(is_initialized());

        const uint16_t blk = i / half_word;
        const uint16_t off = i % half_word;

        // HP/HN/D0/VP/VN below follow the Myers/Hyyrö notation; bind the row's persistent words
        word_t& HP = rows[i].dh_pos;
        word_t& HN = rows[i].dh_neg;
        word_t& D0 = rows[i].diag_zero;
        word_t& frontier = rows[i].frontier;
        const word_t& Eq = eq_masks[blk * sigma + map_sym(Y)];

        // pull in the previous row; descending one row advances the band one bit to the right
        HP = rows[i - 1].dh_pos;
        HN = rows[i - 1].dh_neg;
        frontier = rows[i - 1].frontier << 1;

        // crossing into the next half-word block, realign by discarding the now-vacated lower block
        if (off == 0) {
            HP >>= half_word;
            HN >>= half_word;
            frontier >>= half_word;
        }

        // Myers/Hyyrö bit-parallel column update (fixed by the algorithm)
        D0 = (((Eq & HP) + HP) ^ HP) | Eq | HN;
        word_t VP = HN | ~(D0 | HP);
        word_t VN = D0 & HP;
        HP = (VN << 1) | ~(D0 | (VP << 1));
        HN = (D0 & (VP << 1));

        // the diagonal cell grows by one unless its diagonal delta is zero
        const uint16_t diag_bit = off + diag0;
        rows[i].score = rows[i - 1].score + (D0 & (word_t(1) << diag_bit) ? 0 : 1);

        // move the frontier inward while the band still holds a value <= k_max; if it would leave the
        // bottom of the band, every band cell on this row exceeds k_max and the subtree can be pruned
        if (!(D0 & frontier)) {
            uint16_t surplus = 1;

            while (surplus > 0) {
                if (HP & frontier) surplus--;
                if (HN & frontier) surplus++;
                if (frontier == (word_t(1) << (diag_bit - band_width))) return false;
                frontier >>= 1;
            }
        }

        return true;
    }

    /**
     * @brief whether row @p i reaches the matrix' final column (a full match of the input)
     * @param i row index
     * @return whether row i belongs to the final column
     */
    bool is_in_final_column(const uint16_t i) const
    {
        return i >= num_rows() - last_col_size();
    }

    /**
     * @brief edit distance stored in cell (@p i, @p j); reconstructed from the diagonal score by
     *        counting the horizontal deltas between the diagonal and column j  -- O(1)
     * @param i row index (i < num_rows())
     * @param j column index (j < num_cols())
     * @return the edit distance in cell (i, j)
     */
    uint16_t operator()(uint16_t i, uint16_t j) const
    {
        assert(i < m);
        assert(j < n);
        // bit range [lo, hi) of the horizontal deltas lying between the diagonal and column j
        const uint16_t diag = (i % half_word) + diag0;
        uint16_t lo = (i > j) ? diag - (i - j) + 1 : diag + 1;
        uint16_t hi = (i > j) ? diag + 1 : diag + (j - i) + 1;
        word_t window = ((word_t(1) << (hi - lo)) - word_t(1)) << lo;
        int neg = ::popcount(rows[i].dh_neg & window);
        int pos = ::popcount(rows[i].dh_pos & window);
        uint16_t s = rows[i].score;
        s += (i > j) ? (neg - pos) : (pos - neg);
        return s;
    }

    /**
     * @brief whether, from row @p i on, the band can only grow by vertical gaps (input insertions)
     * @param i row index (i < num_rows())
     * @return whether only vertical gaps remain
     */
    bool only_vertical_gaps_left(uint16_t i) const
    {
        assert(i < m);
        if (i + first_col_bits < n) return false; // the final column n is not reached on row i yet
        const uint16_t blk = i / half_word;
        const uint16_t off = i % half_word;
        // the band only descends if, from the diagonal to the final column, every horizontal delta is +1
        uint16_t lo = diag0 - band_width + off + 1;
        uint16_t hi = diag0 + n - blk * half_word;
        return (((~rows[i].dh_neg >> lo) << lo) << (bits_per_word - hi)) == word_t(0);
    }

    /**
     * @brief number of columns of the matrix (input length + 1)
     * @return the number of columns
     */
    uint16_t num_cols() const
    {
        return n;
    }

    /**
     * @brief number of rows of the matrix
     * @return the number of rows
     */
    uint16_t num_rows() const
    {
        return m;
    }

    /**
     * @brief whether an input has been set (and thus the eq-masks have been precomputed)
     * @return whether the matrix has an input
     */
    bool is_initialized() const
    {
        return !eq_masks.empty();
    }

    /**
     * @brief number of rows that fall into the final column of the matrix
     * @return the size of the final column
     */
    uint16_t last_col_size() const
    {
        return band_height + band_width + 1;
    }
};
