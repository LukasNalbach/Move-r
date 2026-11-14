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
 * @file cigar.hpp
 * @brief CIGAR strings of approximate occurrences: operations, runs and the traceback of an edit-distance matrix
 */

#ifndef MOVE_R_MISC_CIGAR_HPP
#define MOVE_R_MISC_CIGAR_HPP

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <ostream>
#include <ranges>
#include <utility>
#include <vector>

#include <data_structures/edit_distance_matrix.hpp>

/**
 * @brief a CIGAR operation: a match (=), a mismatch (X, substitution), an insertion (I, a pattern symbol absent from
 *        the text) or a deletion (D, a text symbol absent from the pattern)
 */
enum class cigar_op_t : uint8_t { MATCH, MISMATCH, INS, DEL };

/**
 * @brief a run of a CIGAR (packed): a MATCH run consists of len >= 1 matches; a MISMATCH, INS or DEL run consists of
 *        one operation (len = 1) and stores its symbol in sym (the text symbol of a MISMATCH or a DEL, the pattern
 *        symbol of an INS)
 * @tparam sym_t the symbol type of the aligned sequences
 */
#pragma pack(push, 1)
template <typename sym_t>
struct cigar_run_t {
    cigar_op_t op;  // the operation
    uint32_t len;   // the number of operations of the run (1 for MISMATCH, INS and DEL)
    sym_t sym;      // MISMATCH, DEL: the text symbol; INS: the pattern symbol; MATCH: unused ({})

    bool operator==(const cigar_run_t&) const = default;

    /** @brief returns the one-character code of the operation (=, X, I, D) */
    char code() const
    {
        switch (op) {
            case cigar_op_t::MATCH:    return '=';
            case cigar_op_t::MISMATCH: return 'X';
            case cigar_op_t::INS:      return 'I';
            case cigar_op_t::DEL:      return 'D';
            default:      __builtin_unreachable();
        }
    }
};
#pragma pack(pop)

static_assert(sizeof(cigar_run_t<uint8_t>) == sizeof(cigar_op_t) + sizeof(uint32_t) + sizeof(uint8_t));

/**
 * @brief a CIGAR alignment: a sequence of runs (see cigar_run_t)
 * @tparam sym_t the symbol type of the aligned sequences
 */
template <typename sym_t>
using cigar_t = std::vector<cigar_run_t<sym_t>>;

/**
 * @name CIGAR runs: MATCH<sym_t>(len) returns a MATCH run with len matches; MISMATCH(sym), INS(sym) and DEL(sym)
 *       return a run with one operation and the symbol sym
 * @{
 */
template <typename sym_t> inline cigar_run_t<sym_t> MATCH(uint32_t len = 1) { return {cigar_op_t::MATCH, len, sym_t{}}; }
template <typename sym_t> inline cigar_run_t<sym_t> MISMATCH(sym_t sym)     { return {cigar_op_t::MISMATCH, 1, sym}; }
template <typename sym_t> inline cigar_run_t<sym_t> INS(sym_t sym)          { return {cigar_op_t::INS, 1, sym}; }
template <typename sym_t> inline cigar_run_t<sym_t> DEL(sym_t sym)          { return {cigar_op_t::DEL, 1, sym}; }
/** @} */

/** @brief whether an approximate-pattern-matching locate also computes a CIGAR alignment per occurrence */
enum cigar_mode_t : uint8_t { NO_CIGAR = 0, CIGAR = 1 };

/**
 * @brief appends a run to a CIGAR; a MATCH run is merged into a preceding MATCH run (MISMATCH, INS and DEL runs keep
 *        one operation each)
 * @tparam sym_t the symbol type of the CIGAR
 * @param cigar the CIGAR
 * @param run the run
 */
template <typename sym_t>
inline void cigar_push(cigar_t<sym_t>& cigar, cigar_run_t<sym_t> run)
{
    if (run.op == cigar_op_t::MATCH && !cigar.empty() && cigar.back().op == cigar_op_t::MATCH)
        cigar.back().len += run.len;
    else cigar.emplace_back(run);
}

/**
 * @brief reconstructs the substring of the text that a pattern is aligned against from the pattern and its CIGAR: a
 *        MATCH copies symbols of the pattern, a MISMATCH and a DEL append the text symbol stored in the run, and an
 *        INS skips a symbol of the pattern (see cigar_run_t)
 * @tparam seq_t contiguous range type of the pattern and of the result
 * @tparam sym_t the symbol type of the CIGAR
 * @param P the pattern
 * @param cigar the CIGAR of P against the text
 * @return the substring of the text that P is aligned against
 */
template <typename seq_t, typename sym_t>
inline seq_t apply_cigar(const seq_t& P, const cigar_t<sym_t>& cigar)
{
    using out_sym_t = typename seq_t::value_type;
    seq_t text;
    uint64_t i = 0; // position in the pattern

    for (const cigar_run_t<sym_t>& run : cigar) {
        switch (run.op) {
            case cigar_op_t::MATCH:
                for (uint32_t r = 0; r < run.len && i < P.size(); r++) text.push_back(P[i++]);
                break;

            case cigar_op_t::MISMATCH:
                text.push_back(out_sym_t(run.sym));
                i++;
                break;

            case cigar_op_t::INS: // in the pattern, not in the text
                i++;
                break;

            case cigar_op_t::DEL: // in the text, not in the pattern
                text.push_back(out_sym_t(run.sym));
                break;
        }
    }

    return text;
}

/**
 * @brief assigns the ranks 1, 2, ..., sigma to the distinct symbols of R[0..nR-1] that occur in Q[0..nQ-1], in the
 *        order of their first occurrences in R, for an edit-distance matrix with the input R, the query Q and the
 *        default first column (see edit_distance_matrix::init); the other symbols have no rank: a symbol of R without a
 *        rank does not occur in Q, and a symbol of Q without a rank matches no symbol of R. A rank is a byte, so at
 *        most 255 distinct symbols occur in both R and Q (e.g., if R or Q is a string of a text with a byte alphabet,
 *        see move_r::byte_alphabet).
 * @tparam ref_t container type of R
 * @tparam qry_t container type of Q
 * @param R the input of the matrix
 * @param nR length of R
 * @param Q the query of the matrix
 * @param nQ length of Q
 * @param alph_map (output) [0..255] alph_map[c] is the rank of the byte c, or 0 if c has no rank
 * @return sigma, the number of ranks of the matrix: the number of symbols with a rank, at least 1 (a matrix has the
 *         match masks of at least one rank, see edit_distance_matrix::is_initialized)
 */
template <std::ranges::contiguous_range ref_t, std::ranges::contiguous_range qry_t>
inline uint16_t build_alphabet_map(const ref_t& R, int64_t nR, const qry_t& Q, int64_t nQ, std::vector<uint8_t>& alph_map)
{
    // [0..255] in_query[c] <=> the byte c occurs in Q
    std::array<bool, 256> in_query {};
    for (int64_t i = 0; i < nQ; i++) in_query[std::bit_cast<uint8_t>(Q[i])] = true;
    alph_map.assign(256, 0);
    uint16_t sigma = 0;

    for (int64_t j = 0; j < nR; j++) {
        uint8_t b = std::bit_cast<uint8_t>(R[j]);
        if (alph_map[b] == 0 && in_query[b]) alph_map[b] = ++sigma;
    }

    assert(sigma <= std::numeric_limits<uint8_t>::max());
    return std::max<uint16_t>(sigma, 1);
}

/**
 * @brief sets R as the input (columns) of an edit-distance matrix and computes its rows for the query Q[0..nQ-1];
 *        afterwards, the cells of the matrix can be read, but no further row can be computed (the alphabet map that
 *        the matrix refers to is a local variable)
 * @tparam word_t word type of the edit-distance matrix
 * @tparam ref_t container type of R
 * @tparam qry_t container type of Q
 * @param mat the edit-distance matrix
 * @param R the input of the matrix
 * @param nR length of R
 * @param Q the query of the matrix
 * @param nQ length of Q
 * @param k the maximum number of errors
 */
template <typename word_t, std::ranges::contiguous_range ref_t, std::ranges::contiguous_range qry_t>
inline void setup_edit_matrix(edit_distance_matrix<word_t>& mat, const ref_t& R, int64_t nR, const qry_t& Q, int64_t nQ, int64_t k)
{
    std::vector<uint8_t> alph_map;
    uint16_t sigma = build_alphabet_map(R, nR, Q, nQ, alph_map);
    directional_substring<ref_t> R_col(R, 0, nR - 1, RIGHT);
    mat.set_input(R_col, sigma, alph_map);
    mat.init(k);
    for (int64_t i = 1; i <= nQ; i++) mat.compute_row(i, Q[i - 1]);
}

/**
 * @brief returns the distances in the first column of an edit-distance matrix whose input is P[1..m) and whose query is
 *        a string M, such that its alignments of P against prefixes of M do not start with a deletion: in row 0, P[0]
 *        is inserted (distance 1); in row i >= 1, P[0] is aligned to M[0] and M[1..i) is deleted (the distances of
 *        rows 0 and 1 are returned, see edit_distance_matrix::init)
 * @param first_matches whether P[0] matches M[0]
 */
inline const std::vector<uint16_t>& first_column_without_deletion(bool first_matches)
{
    static const std::vector<uint16_t> match { 1, 0 }, mismatch { 1, 1 };
    return first_matches ? match : mismatch;
}

/**
 * @brief aligns P against the prefix M[0..len) of M with the fewest errors e, without a deletion at the start (the
 *        shortest prefix on ties, as edit_dist_anchored_prefix), and computes the CIGAR of the alignment (INS: a symbol
 *        of P that is not in M, DEL: a symbol of M that is not in P); the CIGAR neither starts nor ends with a deletion
 * @tparam seq1_t random-access sequence type of the pattern
 * @tparam seq2_t random-access sequence type of M
 * @param P the pattern
 * @param n1 the length of P
 * @param M the string whose prefixes are aligned against P (it can be longer than the occurrence); M[0..n1 + k) has at
 *        most 255 distinct symbols (e.g., a string of a text with a byte alphabet, see move_r::byte_alphabet)
 * @param n2 the length of M
 * @param k the maximum number of errors, at least 1
 * @return the tuple (e, len, CIGAR of P against M[0..len)); the CIGAR is empty if e > k
 */
template <typename seq1_t, typename seq2_t>
inline std::tuple<uint64_t, uint64_t, cigar_t<std::ranges::range_value_t<seq2_t>>> edit_cigar_prefix(
    const seq1_t& P, int64_t n1, const seq2_t& M, int64_t n2, int64_t k)
{
    using sym_t = std::ranges::range_value_t<seq2_t>;
    assert(1 <= k && k <= edit_distance_max_errors);
    if (n1 == 0) return {uint64_t(0), uint64_t(0), cigar_t<sym_t>{}}; // empty pattern: empty prefix, no edits

    /* P[1..n1) is the input of the matrix and M its query; the first column contains the alignments of P[0] (see
       first_column_without_deletion). An alignment against M[0..len) has at least |len - n1| errors, and at least
       len - n1 + 1 if P[0] does not match M[0], so len lies in [len_min, len_max] */
    const bool first_matches = n2 > 0 && P[0] == M[0];
    const int64_t len_min = n1 > k ? n1 - k : 0;
    const int64_t len_max = std::min<int64_t>(n2, n1 + k - (first_matches ? 0 : 1));
    assert(len_max <= edit_distance_max_length);

    auto run = [&]<typename word_t>() -> std::tuple<uint64_t, uint64_t, cigar_t<sym_t>> {
        /* the ranks of M[0] and of the symbols of P[1..n1) that occur in M[0..len_max); the other symbols have no rank
           (see build_alphabet_map). The query symbol of row 1 needs a rank, since the distance of row 1 in the first
           column does not exceed that of row 0, and a query symbol without a rank does not match the bits of the first
           column (see edit_distance_matrix::_match_masks); in the rows below, the first column grows by 1 per row */
        std::array<bool, 256> in_query {}; // [0..255] in_query[c] <=> the byte c occurs in M[0..len_max)
        for (int64_t x = 0; x < len_max; x++) in_query[std::bit_cast<uint8_t>(M[x])] = true;
        std::vector<uint8_t> alph_map(256, 0);
        uint16_t sigma = 0;

        auto add_rank = [&](auto c) {
            uint8_t& rank = alph_map[std::bit_cast<uint8_t>(c)];
            if (rank == 0) rank = ++sigma;
        };

        if (len_max > 0) add_rank(M[0]);
        for (int64_t x = 1; x < n1; x++) if (in_query[std::bit_cast<uint8_t>(P[x])]) add_rank(P[x]);
        assert(sigma <= std::numeric_limits<uint8_t>::max());

        edit_distance_matrix<word_t> mat;
        // a matrix has the match masks of at least one rank (see edit_distance_matrix::is_initialized)
        mat.set_input(directional_substring<seq1_t>(P, 1, n1 - 1, RIGHT), std::max<uint16_t>(sigma, 1), alph_map);
        mat.init(k, first_column_without_deletion(first_matches));

        // compute the rows 1..len_max; once compute_row returns false, no later row has a cell with at most k errors
        int64_t last_row = 0;

        while (last_row < len_max) {
            last_row++;
            if (!mat.compute_row(last_row, M[last_row - 1])) break;
        }

        // the best prefix M[0..len) has the fewest errors in the last column of row len; on ties, the shortest prefix
        uint64_t best_err = uint64_t(k) + 1;
        int64_t best_len = 0;

        for (int64_t len = len_min; len <= last_row; len++) {
            uint64_t e = mat(len, n1 - 1); // the fewest errors of an alignment against M[0..len)
            if (e < best_err) { best_err = e; best_len = len; }
        }

        if (best_err > uint64_t(k)) return {best_err, uint64_t(best_len), cigar_t<sym_t>{}};

        // the traceback from the cell (best_len, n1 - 1); cells outside the band are skipped
        const int64_t left = mat.band_left(), right = mat.band_right();
        auto in_band = [&](int64_t i, int64_t j) { return i >= 0 && j >= 0 && i - j <= left && j - i <= right; };
        cigar_t<sym_t> cigar;
        int64_t i = best_len, j = n1 - 1;

        // column j > 0 (the column of P[j]): P[j] is aligned to M[i - 1] or inserted, or M[i - 1] is deleted
        while (j > 0) {
            const int64_t cur = mat(i, j);

            if (i > 0 && in_band(i - 1, j - 1) && mat(i - 1, j - 1) + (P[j] != M[i - 1] ? 1 : 0) == cur) {
                cigar_push(cigar, P[j] == M[i - 1] ? MATCH<sym_t>() : MISMATCH(sym_t(M[i - 1])));
                i--; j--;
            } else if (in_band(i, j - 1) && mat(i, j - 1) + 1 == cur) {
                cigar_push(cigar, INS(sym_t(P[j])));
                j--;
            } else {
                assert(i > 0 && in_band(i - 1, j) && mat(i - 1, j) + 1 == cur);
                cigar_push(cigar, DEL(sym_t(M[i - 1])));
                i--;
            }
        }

        // column 0: P[0] is inserted (row 0), or aligned to M[0] and followed by the deletion of M[1..i)
        if (i == 0) {
            cigar_push(cigar, INS(sym_t(P[0])));
        } else {
            for (int64_t x = i - 1; x > 0; x--) cigar_push(cigar, DEL(sym_t(M[x])));
            cigar_push(cigar, P[0] == M[0] ? MATCH<sym_t>() : MISMATCH(sym_t(M[0])));
        }

        std::reverse(cigar.begin(), cigar.end());
        return {best_err, uint64_t(best_len), std::move(cigar)};
    };

    return with_edit_distance_word(k, run);
}

/**
 * @brief counts the edits of a CIGAR and checks that it is an alignment of P against M: every MATCH is on equal
 *        symbols, every MISMATCH on different symbols, and the runs together consume exactly P and M
 * @tparam sym_t symbol type of the CIGAR
 * @tparam seq_sym_t symbol type of P and M; can differ from sym_t in signedness
 * @param P the pattern
 * @param n1 the length of P
 * @param M the string P is aligned against
 * @param n2 the length of M
 * @param cigar the CIGAR to check
 * @return the number of edits (mismatches + insertions + deletions) if the CIGAR aligns P against M, else -1
 */
template <typename sym_t, typename seq_sym_t>
inline int64_t num_cigar_edits(const seq_sym_t* P, uint64_t n1, const seq_sym_t* M, uint64_t n2, const cigar_t<sym_t>& cigar)
{
    uint64_t p1 = 0, p2 = 0; // the next positions in P and M
    int64_t edits = 0;

    // compares a sequence symbol and a CIGAR symbol as unsigned values (one can be a char, the other a uint8_t)
    auto same = [](seq_sym_t a, sym_t b) {
        return uint64_t(std::make_unsigned_t<seq_sym_t>(a)) == uint64_t(std::make_unsigned_t<sym_t>(b));
    };

    for (const cigar_run_t<sym_t>& run : cigar) {
        for (uint32_t r = 0; r < run.len; r++) {
            switch (run.op) {
                case cigar_op_t::MATCH:    if (p1 >= n1 || p2 >= n2 || P[p1] != M[p2]) return -1; p1++; p2++; break;
                case cigar_op_t::MISMATCH: if (p1 >= n1 || p2 >= n2 || P[p1] == M[p2] || !same(M[p2], run.sym)) return -1; p1++; p2++; edits++; break;
                case cigar_op_t::INS:      if (p1 >= n1 || !same(P[p1], run.sym)) return -1; p1++; edits++; break;
                case cigar_op_t::DEL:      if (p2 >= n2 || !same(M[p2], run.sym)) return -1; p2++; edits++; break;
                default:                   __builtin_unreachable();
            }
        }
    }

    return (p1 == n1 && p2 == n2) ? edits : -1;
}

/**
 * @brief writes the SAM CIGAR string of a CIGAR with the operations =, X, I and D to out, e.g. "5=1X2I3=";
 *        consecutive runs of the same operation are merged (two adjacent mismatches are written as "2X")
 * @tparam out_t type of the output (an output stream or a string_appender)
 * @tparam sym_t the symbol type of the CIGAR
 * @param out the output the CIGAR is written to
 * @param cigar the CIGAR alignment
 */
template <typename out_t, typename sym_t>
inline void append_cigar(out_t& out, const cigar_t<sym_t>& cigar)
{
    for (uint64_t i = 0; i < cigar.size();) {
        cigar_op_t op = cigar[i].op;
        uint64_t len = 0;
        uint64_t j = i;
        while (j < cigar.size() && cigar[j].op == op) len += cigar[j++].len;
        out << len << cigar[i].code();
        i = j;
    }
}

/**
 * @brief writes the value of the SAM MD tag of a CIGAR to out: each run of matches as its length, each mismatch as its
 *        text symbol, and each maximal run of deletions as '^' followed by the deleted text symbols; insertions are not
 *        part of the MD tag
 *
 * The text symbols are the ones stored in the runs of the CIGAR. For example, "10=1X5=2D6=" with the text symbol G of
 * the mismatch and the deleted text symbols TA is written as "10G5^TA6".
 * @tparam out_t type of the output (an output stream or a string_appender)
 * @tparam sym_t the symbol type of the CIGAR
 * @param out the output the MD-tag value is written to
 * @param cigar the CIGAR alignment
 */
template <typename out_t, typename sym_t>
inline void append_md(out_t& out, const cigar_t<sym_t>& cigar)
{
    uint64_t matches = 0; // length of the current run of matches
    bool in_del = false;  // whether the previous run was a deletion

    for (const cigar_run_t<sym_t>& run : cigar) {
        switch (run.op) {
            case cigar_op_t::MATCH:    matches += run.len; in_del = false; break;
            case cigar_op_t::MISMATCH: out << matches << run.sym; matches = 0; in_del = false; break;
            case cigar_op_t::INS:      in_del = false; break; // insertions do not appear in the MD tag
            case cigar_op_t::DEL:
                if (!in_del) { out << matches << '^'; matches = 0; in_del = true; }
                out << run.sym;
                break;
            default: __builtin_unreachable();
        }
    }

    out << matches; // the MD tag ends with a (possibly zero) match count
}

#endif // MOVE_R_MISC_CIGAR_HPP
