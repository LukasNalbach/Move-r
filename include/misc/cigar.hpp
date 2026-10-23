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

#ifndef MOVE_R_MISC_CIGAR_HPP
#define MOVE_R_MISC_CIGAR_HPP

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <ostream>
#include <ranges>
#include <utility>
#include <vector>

#include <data_structures/edit_distance_matrix.hpp>

/** @brief a CIGAR operation: a match (=), a mismatch (X, substitution), an insertion (I, a pattern symbol absent
 *         from the text) or a deletion (D, a text symbol absent from the pattern) */
enum class cigar_op_t : uint8_t { MATCH, MISMATCH, INS, DEL };

/** @brief one run of a CIGAR (packed): a MATCH run spans len (>= 1) positions; a MISMATCH/INS/DEL run spans one
 *         position (len == 1) and stores its symbol in @ref sym: the reference symbol for a MISMATCH or DEL, the
 *         pattern symbol for an INS
 * @tparam sym_t the symbol type of the aligned sequences */
#pragma pack(push, 1)
template <typename sym_t>
struct cigar_run_t {
    cigar_op_t op;  // the operation
    uint32_t len;   // number of positions the run spans (always 1 for MISMATCH/INS/DEL)
    sym_t sym;      // MISMATCH/DEL: the reference symbol; INS: the pattern symbol; unused ({}) for MATCH

    bool operator==(const cigar_run_t&) const = default;

    /** @brief the one-character code of the operation (=, X, I, D) */
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

/** @brief a CIGAR alignment as a sequence of runs (MISMATCH/INS/DEL runs are one position each and carry a symbol)
 * @tparam sym_t the symbol type of the aligned sequences */
template <typename sym_t>
using cigar_t = std::vector<cigar_run_t<sym_t>>;

/** @name CIGAR run factories: a MATCH takes a run length (and the symbol type); a MISMATCH/INS/DEL takes its single
 *        symbol (its length is 1), e.g. MATCH<char>(3), MISMATCH('A'), DEL('C').
 *  @{ */
template <typename sym_t> inline cigar_run_t<sym_t> MATCH(uint32_t len = 1) { return {cigar_op_t::MATCH, len, sym_t{}}; }
template <typename sym_t> inline cigar_run_t<sym_t> MISMATCH(sym_t sym)     { return {cigar_op_t::MISMATCH, 1, sym}; }
template <typename sym_t> inline cigar_run_t<sym_t> INS(sym_t sym)          { return {cigar_op_t::INS, 1, sym}; }
template <typename sym_t> inline cigar_run_t<sym_t> DEL(sym_t sym)          { return {cigar_op_t::DEL, 1, sym}; }
/** @} */

/** @brief whether an approximate-pattern-matching locate also computes a CIGAR alignment per occurrence */
enum cigar_mode_t : uint8_t { NO_CIGAR = 0, CIGAR = 1 };

/** @brief appends a CIGAR run; a MATCH run extends a preceding MATCH run (MISMATCH/INS/DEL runs stay one position
 *         each) */
template <typename sym_t>
inline void cigar_push(cigar_t<sym_t>& cigar, cigar_run_t<sym_t> run)
{
    if (run.op == cigar_op_t::MATCH && !cigar.empty() && cigar.back().op == cigar_op_t::MATCH)
        cigar.back().len += run.len;
    else cigar.emplace_back(run);
}

/**
 * @brief reconstructs the text a pattern is aligned against from the pattern and its CIGAR: a MATCH copies the
 *        symbols of the pattern, a MISMATCH and a DEL append the reference symbol stored in the run and an INS
 *        skips a symbol of the pattern (see cigar_run_t)
 * @tparam seq_t contiguous range type of the pattern and of the result
 * @tparam sym_t the symbol type of the CIGAR
 * @param P the pattern
 * @param cigar the CIGAR of P against the text
 * @return the text P is aligned against
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
 * @brief maps the symbols occurring in R[0..nR-1] to [1..sigma]
 * @tparam ref_t container type of the reference
 * @param R the reference
 * @param nR length of R
 * @param alph_map [0..255] (output) alph_map[c] = index in [1..sigma] of the symbol c, or 0 if c does not occur in R
 * @return sigma, the number of distinct symbols in R
 */
template <std::ranges::contiguous_range ref_t>
static uint16_t build_alphabet_map(const ref_t& R, int64_t nR, std::vector<uint8_t>& alph_map)
{
    alph_map.assign(256, 0);
    uint16_t sigma = 0;

    for (int64_t j = 0; j < nR; j++) {
        uint8_t b = std::bit_cast<uint8_t>(R[j]);
        if (alph_map[b] == 0) alph_map[b] = ++sigma;
    }

    return sigma;
}

/**
 * @brief sets R as the input (columns) of an edit-distance matrix and computes its rows for Q[0..nQ-1]; afterwards,
 *        the cells of the matrix can be read, but no further row can be computed (the alphabet map that the matrix
 *        refers to is a local variable)
 * @tparam word_t word type of the edit-distance matrix
 * @tparam ref_t container type of the reference
 * @tparam qry_t container type of the query
 * @param mat the edit-distance matrix
 * @param R the reference
 * @param nR length of R
 * @param Q the query
 * @param nQ length of Q
 * @param k the maximum number of errors
 */
template <typename word_t, std::ranges::contiguous_range ref_t, std::ranges::contiguous_range qry_t>
static void setup_edit_matrix(edit_distance_matrix<word_t>& mat, const ref_t& R, int64_t nR, const qry_t& Q, int64_t nQ, int64_t k)
{
    std::vector<uint8_t> alph_map;
    uint16_t sigma = build_alphabet_map(R, nR, alph_map);
    directional_substring<ref_t> R_col(R, 0, nR - 1, RIGHT);
    mat.set_input(R_col, sigma, alph_map);
    mat.init(k);
    for (int64_t i = 1; i <= nQ; i++) mat.compute_row(i, Q[i - 1]);
}

/* computes the CIGAR of P against M by a traceback from cell (i, j) of a banded matrix with M in the columns and P in
   the rows; cells with |i - j| > k are skipped */
template <typename word_t, typename seq1_t, typename seq2_t>
static cigar_t<std::ranges::range_value_t<seq2_t>> cigar_traceback(
    const edit_distance_matrix<word_t>& mat, const seq1_t& P, const seq2_t& M, int64_t i, int64_t j, int64_t k)
{
    using sym_t = std::ranges::range_value_t<seq2_t>;
    auto in_band = [&](int64_t a, int64_t b) { return a >= 0 && b >= 0 && std::abs(a - b) <= k; };
    cigar_t<sym_t> cigar;

    while (i > 0 || j > 0) {
        int64_t cur = mat(i, j);

        if (i > 0 && j > 0 && in_band(i - 1, j - 1) && mat(i - 1, j - 1) + (P[i - 1] != M[j - 1] ? 1 : 0) == cur) {
            cigar_push(cigar, P[i - 1] == M[j - 1] ? MATCH<sym_t>() : MISMATCH(sym_t(M[j - 1]))); // reference symbol
            i--; j--;
        } else if (i > 0 && in_band(i - 1, j) && mat(i - 1, j) + 1 == cur) {
            cigar_push(cigar, INS(sym_t(P[i - 1]))); // a P-symbol absent from M
            i--;
        } else if (j > 0 && in_band(i, j - 1) && mat(i, j - 1) + 1 == cur) {
            cigar_push(cigar, DEL(sym_t(M[j - 1]))); // an M-symbol absent from P
            j--;
        } else if (i > 0) {
            cigar_push(cigar, INS(sym_t(P[i - 1]))); // j == 0: only insertions of P remain
            i--;
        } else {
            cigar_push(cigar, DEL(sym_t(M[j - 1]))); // i == 0: only deletions of M remain
            j--;
        }
    }

    std::reverse(cigar.begin(), cigar.end());
    return cigar;
}

/**
 * @brief computes ed(P, M) and a CIGAR of the pattern P against the matched string M with a banded bit-parallel
 *        dynamic program (Myers) and a traceback; the CIGAR has M as the reference (INS = a symbol of P missing in M,
 *        DEL = a symbol of M missing in P). Requires ed(P, M) <= k, k <= edit_distance_max_errors_128 and n2 <= edit_distance_max_length.
 * @tparam seq1_t random-access sequence type of the pattern
 * @tparam seq2_t random-access sequence type of the matched string
 * @param P the pattern
 * @param n1 the length of P
 * @param M the matched string
 * @param n2 the length of M
 * @param k the band half-width (and error bound)
 * @return the pair (ed(P, M), CIGAR of P against M)
 */
template <typename seq1_t, typename seq2_t>
static std::pair<uint64_t, cigar_t<std::ranges::range_value_t<seq2_t>>> edit_cigar(
    const seq1_t& P, int64_t n1, const seq2_t& M, int64_t n2, int64_t k)
{
    using sym_t = std::ranges::range_value_t<seq2_t>;
    assert(k <= edit_distance_max_errors_128);
    assert(n2 <= edit_distance_max_length);

    // an empty side aligns trivially: only insertions of P (M empty) or deletions of M (P empty) remain
    if (n1 == 0 || n2 == 0) {
        cigar_t<sym_t> cigar;
        for (int64_t i = 0; i < n1; i++) cigar.emplace_back(INS(sym_t(P[i]))); // M empty: all of P is inserted
        for (int64_t j = 0; j < n2; j++) cigar.emplace_back(DEL(sym_t(M[j]))); // P empty: all of M is deleted
        return {uint64_t(n1 + n2), std::move(cigar)};
    }

    // builds the banded Myers matrix with reference = M (columns), query = P (rows), then traces back from cell (n1, n2)
    auto run = [&]<typename word_t>() -> std::pair<uint64_t, cigar_t<sym_t>> {
        edit_distance_matrix<word_t> mat;
        setup_edit_matrix<word_t>(mat, M, n2, P, n1, k); // reference = M, query = P
        return {mat(n1, n2), cigar_traceback(mat, P, M, n1, n2, k)};
    };

    if (k <= edit_distance_max_errors_64)
         return run.template operator()<uint64_t>();
    else return run.template operator()<uint128_t>();
}

/**
 * @brief finds the prefix M[0..len) of the matched string M with the smallest e = ed(P, M[0..len)) (the shortest one on
 *        ties, as edit_dist_prefix) and computes a CIGAR of the pattern P against it (M as the reference, see
 *        edit_cigar)
 * @tparam seq1_t random-access sequence type of the pattern
 * @tparam seq2_t random-access sequence type of the matched string
 * @param P the pattern
 * @param n1 the length of P
 * @param M the matched string (may be longer than the occurrence; only its best prefix is aligned)
 * @param n2 the length of M
 * @param k the band half-width (and error bound)
 * @return the tuple (e, len, CIGAR of P against M[0..len)); the CIGAR is empty if e > k
 */
template <typename seq1_t, typename seq2_t>
static std::tuple<uint64_t, uint64_t, cigar_t<std::ranges::range_value_t<seq2_t>>> edit_cigar_prefix(
    const seq1_t& P, int64_t n1, const seq2_t& M, int64_t n2, int64_t k)
{
    using sym_t = std::ranges::range_value_t<seq2_t>;
    assert(k <= edit_distance_max_errors_128);

    // only a prefix of length in [n1 - k, n1 + k] can be within edit distance k of P (a length difference of d costs
    // at least d indels); prefixes past n1 + k leave the band, so the reference is capped there
    const int64_t j_min = n1 > k ? n1 - k : 0;
    const int64_t j_max = std::min<int64_t>(n2, n1 + k);
    assert(j_max <= edit_distance_max_length);

    if (n1 == 0) return {uint64_t(0), uint64_t(0), cigar_t<sym_t>{}}; // empty pattern: empty prefix, no edits

    auto run = [&]<typename word_t>() -> std::tuple<uint64_t, uint64_t, cigar_t<sym_t>> {
        edit_distance_matrix<word_t> mat;
        setup_edit_matrix<word_t>(mat, M, j_max, P, n1, k); // reference = M[0..j_max), query = P

        // the best prefix is the column j (in the band at row n1) with the smallest ed(P, M[0..j)); ties keep the
        // shortest length
        uint64_t best_err = uint64_t(k) + 1;
        int64_t best_len = 0;

        for (int64_t j = j_min; j <= j_max; j++) {
            uint64_t e = mat(n1, j); // = ed(P, M[0..j))
            if (e < best_err) { best_err = e; best_len = j; }
        }

        if (best_err > uint64_t(k)) return {best_err, uint64_t(best_len), cigar_t<sym_t>{}};
        return {best_err, uint64_t(best_len), cigar_traceback(mat, P, M, n1, best_len, k)};
    };

    if (k <= edit_distance_max_errors_64)
         return run.template operator()<uint64_t>();
    else return run.template operator()<uint128_t>();
}

/**
 * @brief counts the edits of a CIGAR and checks that it is an alignment of P against M: every MATCH is on equal
 *        symbols, every MISMATCH on different symbols, and the runs together consume exactly P and M
 * @tparam sym_t symbol type of the CIGAR
 * @tparam seq_sym_t symbol type of P and M; can differ from sym_t in signedness
 * @param P the pattern
 * @param n1 the length of P
 * @param M the matched string
 * @param n2 the length of M
 * @param cigar the CIGAR to check
 * @return the number of edits (mismatches + insertions + deletions) if the CIGAR aligns P against M, else -1
 */
template <typename sym_t, typename seq_sym_t>
static int64_t num_cigar_edits(const seq_sym_t* P, uint64_t n1, const seq_sym_t* M, uint64_t n2, const cigar_t<sym_t>& cigar)
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
 * @brief writes the SAM MD-tag value of a CIGAR (the alignment of a read against the reference) to out: runs of
 *        matching bases as counts, each mismatch as the reference base, and each maximal run of deletions as '^'
 *        followed by the deleted reference bases (insertions are not part of the MD tag). E.g. "10=1X5=2D6=" over
 *        reference bases G / TA yields "10G5^TA6"; the reference bases are the symbols stored in the CIGAR runs
 * @tparam out_t type of the output (an output stream or a string_appender)
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
