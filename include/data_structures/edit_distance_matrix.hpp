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

#ifndef MOVE_R_DATA_STRUCTURES_EDIT_DISTANCE_MATRIX_HPP
#define MOVE_R_DATA_STRUCTURES_EDIT_DISTANCE_MATRIX_HPP

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <ranges>
#include <tuple>
#include <vector>
#include <cmath>
#include <sstream>
#include <string>

#include <misc/directional_substring.hpp>
#include <misc/utils.hpp>

// the matches of an input of an edit-distance matrix without further equal symbols (see set_input)
struct no_extra_matches_t {
    template <typename report_fnc_t>
    void operator()(uint64_t, report_fnc_t&&) const {}
};

/**
 * @brief banded edit-distance matrix between a fixed input (columns) and a query (rows), computed row by
 *        row with the bit-parallel algorithm of Myers (J. ACM 46(3), 1999) and Hyyrö; a row stores only
 *        the band of width 2k+1 around the diagonal in one word (as in Columba)
 *
 * The rows are grouped into blocks of half_word rows. Bit x of the word of row i in block b = i / half_word
 * corresponds to column x + b * half_word - first_col_bits + 1 (i.e., to input position x + b * half_word -
 * first_col_bits); the cell (i, i) of the main diagonal is at bit (i mod half_word) + block_diag_bit.
 * @tparam word_t unsigned type of the word storing a row of the band (uint64_t or uint128_t)
 */
template <typename word_t>
class edit_distance_matrix {
    static_assert(std::is_same_v<word_t, uint64_t> || std::is_same_v<word_t, uint128_t>);

public:
    // ############################# TYPES AND PARAMETERS #############################

    static constexpr uint16_t bits_per_word = sizeof(word_t) * 8;
    static constexpr uint16_t half_word = bits_per_word / 2; // number of rows of a block
    static constexpr uint16_t max_supported_errors = (bits_per_word - half_word - 2) / 3; // maximum supported number of errors
    static constexpr uint16_t first_col_bits = 2 * max_supported_errors + 1; // number of bits of the first column of the band
    static constexpr uint16_t block_diag_bit = 2 * max_supported_errors; // bit of the main diagonal in the first row of a block

protected:
    // a row of the band (HP, HN and D0 of Myers/Hyyrö)
    struct row_state_t {
        word_t diag_zero; // D0: positions whose diagonal delta is 0
        word_t dh_pos; // HP: positions whose horizontal delta is +1
        word_t dh_neg; // HN: positions whose horizontal delta is -1
        word_t last_active_bit; // single set bit at the rightmost band cell with value <= k_max
        uint16_t diag_distance; // edit distance in the cell of this row on the main diagonal
    };

    // ############################# VARIABLES #############################

    uint16_t sigma = 0; // size of the effective alphabet
    uint16_t k_max = 0; // maximum number of errors of the current query
    uint16_t _num_rows = 0; // number of rows (band length)
    uint16_t _num_cols = 0; // number of columns (input length + 1)
    uint16_t band_width = 0; // number of columns left of the diagonal in the band
    uint16_t band_height = 0; // number of rows below the diagonal whose value can still be <= k_max

    // maps a byte to its one-based rank in the effective alphabet (0 if it does not occur)
    const std::vector<uint8_t>* symbol_rank = nullptr;

    // the rows of the band: [0.._num_rows-1] rows[i] is row i (after i symbols of the query)
    std::vector<row_state_t> rows;

    /* the positions of the input that match each symbol, per block: [0..(number of blocks) * sigma - 1]
       bit x of match_masks[b * sigma + c] is 1 <=> the input position of bit x in block b (see the class
       description) has rank c + 1, or matches it (see set_input), or lies before the input (first column) */
    std::vector<word_t> match_masks;

    // ############################# INTERNAL METHODS #############################

    void resize_rows(uint16_t num) { rows.resize(num); }

    // maps a query symbol to its zero-based index in the match_masks of a block
    template <typename sym_t>
    uint8_t map_sym(sym_t sym) const { return (*symbol_rank)[std::bit_cast<uint8_t>(sym)] - 1; }

public:
    // ############################# CONSTRUCTORS #############################

    edit_distance_matrix() = default;

    // ############################# MISC PUBLIC METHODS #############################

    // returns the number of columns of the matrix (input length + 1)
    uint16_t num_cols() const { return _num_cols; }

    // returns the number of rows of the matrix
    uint16_t num_rows() const { return _num_rows; }

    // returns whether an input has been set (see set_input)
    bool is_initialized() const { return !match_masks.empty(); }

    // returns the number of rows with a cell in the last column of the matrix
    uint16_t last_col_size() const { return band_height + band_width + 1; }

    // returns whether the band in row i contains a cell of the last column (i.e., the whole input is matched)
    bool is_in_final_column(const uint16_t i) const { return i >= num_rows() - last_col_size(); }

    // ############################# CONSTRUCTION METHODS #############################

    /**
     * @brief sets the input (the columns of the matrix) and precomputes its match_masks
     * @tparam inp_t type of the sequence underlying the input
     * @tparam match_fnc_t type of matches
     * @param input the input
     * @param sigma size of the effective alphabet
     * @param symbol_rank byte -> one-based symbol rank in the effective alphabet (0 if the byte does not occur)
     * @param matches called as matches(pos, report) for the position pos of an input symbol in the underlying
     *        sequence; it calls report(rank) for each further one-based rank that is equal to that symbol
     */
    template <typename inp_t, typename match_fnc_t = no_extra_matches_t>
    void set_input(const directional_substring<inp_t>& input, uint16_t sigma, const std::vector<uint8_t>& symbol_rank,
        match_fnc_t matches = {})
    {
        this->sigma = sigma;
        this->symbol_rank = &symbol_rank;
        _num_cols = input.size() + 1;
        _num_rows = 2 * max_supported_errors + _num_cols;
        uint16_t num_blocks = (_num_rows + half_word - 1) / half_word;
        match_masks.resize(num_blocks * sigma);

        // block 0: the lowest first_col_bits bits (the first column) are set for every symbol, the
        // remaining bits belong to the first input symbols
        word_t low_ones = (word_t(1) << first_col_bits) - word_t(1);
        std::fill(match_masks.begin(), match_masks.begin() + sigma, low_ones);
        word_t bit = word_t(1) << first_col_bits;
        uint16_t j_end = std::min<uint16_t>(input.size(), bits_per_word - first_col_bits);

        for (uint16_t j = 0; j < j_end; j++) {
            uint8_t rank = symbol_rank[std::bit_cast<uint8_t>(input[j])];
            if (rank != 0) match_masks[rank - 1] |= bit;
            matches(input.pos(j), [&](uint8_t r) { match_masks[r - 1] |= bit; });
            bit <<= 1;
        }

        // each further block is the previous block shifted right by half_word bits, with the next
        // half_word input symbols in its upper half
        for (uint16_t block = 1; block < num_blocks; block++) {
            uint64_t base = uint64_t(block) * sigma;
            uint64_t prev = base - sigma;

            for (uint64_t c = 0; c < sigma; c++) {
                match_masks[base + c] = match_masks[prev + c] >> half_word;
            }

            bit = word_t(1) << (bits_per_word - half_word);
            uint16_t j_beg = bits_per_word - first_col_bits + (block - 1) * half_word;
            uint16_t block_end = std::min<uint16_t>(input.size(), j_beg + half_word);

            for (uint16_t j = j_beg; j < block_end; j++) {
                uint8_t rank = symbol_rank[std::bit_cast<uint8_t>(input[j])];
                if (rank != 0) match_masks[base + rank - 1] |= bit;
                matches(input.pos(j), [&](uint8_t r) { match_masks[base + r - 1] |= bit; });
                bit <<= 1;
            }
        }
    }

    /**
     * @brief initializes the matrix for a query with at most k_max errors
     * @param k_max maximum number of errors (k_max <= max_supported_errors)
     * @param init_dists (optional) distances in the first column, to continue the band of an adjacent matrix
     * @param min_band_width (optional) the band extends at least this many columns left of the diagonal (see add_alignment_start)
     */
    void init(uint16_t k_max, const std::vector<uint16_t>& init_dists = {}, uint16_t min_band_width = 0)
    {
        assert(k_max <= max_supported_errors);
        assert(is_initialized());
        assert(init_dists.empty() || init_dists.front() <= k_max);
        assert(init_dists.empty() || init_dists.back() <= k_max);

        this->k_max = k_max;
        band_width = init_dists.empty() ? k_max : init_dists.size() - 1 + k_max - init_dists.back();
        band_width = std::max(band_width, min_band_width);
        assert(band_width <= 2 * max_supported_errors);
        _num_rows = band_width + _num_cols;
        resize_rows(_num_rows);
        rows[0].diag_distance = init_dists.empty() ? 0 : init_dists[0];
        band_height = k_max - rows[0].diag_distance;

        if (band_width + band_height + 1 > _num_rows) {
            _num_rows = band_width + band_height + 1;
            resize_rows(_num_rows);
        }

        // first column: the distances are 0, 1, 2, ... by default, i.e., every horizontal delta is +1
        rows[0].dh_pos = (~word_t(0)) << first_col_bits;
        rows[0].dh_neg = ~rows[0].dh_pos;

        // encode the given first-column distances as horizontal deltas
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

        rows[0].last_active_bit = word_t(1) << (block_diag_bit + band_height);
    }

    // ############################# COMPUTING ROWS #############################

    /**
     * @brief computes row i from row i-1 for the next query symbol
     * @tparam sym_t symbol type (one byte)
     * @param i [1..num_rows()-1] row to compute
     * @param sym next query symbol
     * @return false if all values of the band in row i exceed k_max, else true
     */
    template <typename sym_t>
    bool compute_row(uint16_t i, sym_t sym)
    {
        assert(i > 0);
        assert(i < _num_rows);
        assert(is_initialized());

        const uint16_t block = i / half_word;
        const uint16_t offset = i % half_word;

        // HP, HN, D0, VP and VN as in Myers/Hyyrö
        word_t& HP = rows[i].dh_pos;
        word_t& HN = rows[i].dh_neg;
        word_t& D0 = rows[i].diag_zero;
        word_t& last_active_bit = rows[i].last_active_bit;
        const uint8_t rank = map_sym(sym);
        const word_t Eq = (rank == uint8_t(-1)) ? word_t(0) : match_masks[block * sigma + rank];

        // copy row i-1; the band advances by one bit per row
        HP = rows[i - 1].dh_pos;
        HN = rows[i - 1].dh_neg;
        last_active_bit = rows[i - 1].last_active_bit << 1;

        // row i is the first row of a block: shift out the lower half of the word
        if (offset == 0) {
            HP >>= half_word;
            HN >>= half_word;
            last_active_bit >>= half_word;
        }

        // bit-parallel update of Myers/Hyyrö
        D0 = (((Eq & HP) + HP) ^ HP) | Eq | HN;
        word_t VP = HN | ~(D0 | HP);
        word_t VN = D0 & HP;
        HP = (VN << 1) | ~(D0 | (VP << 1));
        HN = (D0 & (VP << 1));

        // the value on the diagonal increases by one, unless its diagonal delta is 0
        const uint16_t diag_bit = offset + block_diag_bit;
        rows[i].diag_distance = rows[i - 1].diag_distance + (D0 & (word_t(1) << diag_bit) ? 0 : 1);

        // move last_active_bit to the rightmost band cell with a value <= k_max; return false if there is none
        if (!(D0 & last_active_bit)) {
            // the value of the current cell minus k_max
            uint16_t excess = 1;

            while (excess > 0) {
                if (HP & last_active_bit) excess--;
                if (HN & last_active_bit) excess++;
                if (last_active_bit == (word_t(1) << (diag_bit - band_width))) return false;
                last_active_bit >>= 1;
            }
        }

        return true;
    }

    /**
     * @brief lowers every cell (i, j) of row i to dist + |j|, i.e., adds an alignment that starts in cell (i, 0) with
     *        distance dist (row i must lie in the first block, i.e., i < half_word); the band must contain cell (i, 0)
     *        and the cells (i, j) with j <= k_max - dist
     * @param i [0..num_rows()-1] row (computed)
     * @param dist the distance of the new alignment in cell (i, 0)
     * @return false if all values of the band in row i exceed k_max, else true
     */
    bool add_alignment_start(uint16_t i, uint16_t dist)
    {
        assert(i < _num_rows && i < half_word);
        assert(i <= band_width && dist + i + band_height >= k_max);
        row_state_t& row = rows[i];
        const uint16_t diag_bit = i + block_diag_bit;
        auto delta = [&](uint16_t b) { return int32_t((row.dh_pos >> b) & 1) - int32_t((row.dh_neg >> b) & 1); };

        // the values of row i at all bits of the word; bit b holds the cell in column i + b - diag_bit
        std::array<int32_t, bits_per_word> vals;
        vals[diag_bit] = row.diag_distance;
        for (uint16_t b = diag_bit + 1; b < bits_per_word; b++) vals[b] = vals[b - 1] + delta(b);
        for (uint16_t b = diag_bit; b > 0; b--) vals[b - 1] = vals[b] - delta(b);
        const int32_t v0 = vals[0];

        for (uint16_t b = 0; b < bits_per_word; b++) {
            int32_t j = int32_t(i) + b - diag_bit;
            vals[b] = std::min<int32_t>(vals[b], dist + (j < 0 ? -j : j));
        }

        // re-encode the horizontal deltas; left of bit 0, a lowered value is followed by one that is one higher
        word_t dh_pos = row.dh_pos & word_t(1), dh_neg = row.dh_neg & word_t(1);
        if (vals[0] < v0) { dh_pos = 0; dh_neg = 1; }

        for (uint16_t b = 1; b < bits_per_word; b++) {
            if (vals[b] > vals[b - 1]) dh_pos |= word_t(1) << b;
            else if (vals[b] < vals[b - 1]) dh_neg |= word_t(1) << b;
        }

        row.dh_pos = dh_pos;
        row.dh_neg = dh_neg;
        row.diag_distance = vals[diag_bit];

        // the rightmost cell of the band with a value <= k_max
        for (uint16_t b = diag_bit + band_height + 1; b-- > diag_bit - band_width;) {
            if (vals[b] <= k_max) {
                row.last_active_bit = word_t(1) << b;
                return true;
            }
        }

        return false;
    }

    // ############################# QUERY METHODS #############################

    /**
     * @brief returns the edit distance in cell (i, j), computed in O(1) time from the distance on the diagonal
     *        and the horizontal deltas between the diagonal and column j
     * @param i [0..num_rows()-1] row
     * @param j [0..num_cols()-1] column
     */
    uint16_t operator()(uint16_t i, uint16_t j) const
    {
        assert(i < _num_rows);
        assert(j < _num_cols);
        // bit range [lo, hi) of the horizontal deltas lying between the diagonal and column j
        const uint16_t diag = (i % half_word) + block_diag_bit;
        uint16_t lo = (i > j) ? diag - (i - j) + 1 : diag + 1;
        uint16_t hi = (i > j) ? diag + 1 : diag + (j - i) + 1;
        word_t window = ((word_t(1) << (hi - lo)) - word_t(1)) << lo;
        int32_t neg = ::popcount(rows[i].dh_neg & window);
        int32_t pos = ::popcount(rows[i].dh_pos & window);
        uint16_t s = rows[i].diag_distance;
        s += (i > j) ? (neg - pos) : (pos - neg);
        return s;
    }

    /**
     * @brief returns whether from row i on, the matrix can only be extended by vertical gaps (insertions)
     * @param i [0..num_rows()-1] row
     */
    bool only_vertical_gaps_left(uint16_t i) const
    {
        assert(i < _num_rows);
        if (i + first_col_bits < _num_cols) return false; // the band in row i does not reach the last column
        const uint16_t block = i / half_word;
        const uint16_t offset = i % half_word;
        // check whether dh_neg is set at all positions [lo, hi) of the band up to the last column
        uint16_t lo = block_diag_bit - band_width + offset + 1;
        uint16_t hi = block_diag_bit + _num_cols - block * half_word;
        return (((~rows[i].dh_neg >> lo) << lo) << (bits_per_word - hi)) == word_t(0);
    }
};

static constexpr int64_t edit_distance_max_errors_64 = edit_distance_matrix<uint64_t>::max_supported_errors;
static constexpr int64_t edit_distance_max_errors_128 = edit_distance_matrix<uint128_t>::max_supported_errors;
// the largest length of the input and of the compared string of a matrix: the numbers of rows and columns are uint16_t,
// and a row index plus first_col_bits must fit into a uint16_t
static constexpr int64_t edit_distance_max_length =
    std::numeric_limits<uint16_t>::max() - edit_distance_matrix<uint128_t>::first_col_bits;

#endif // MOVE_R_DATA_STRUCTURES_EDIT_DISTANCE_MATRIX_HPP
