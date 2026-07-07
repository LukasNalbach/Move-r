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
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <ostream>
#include <utility>
#include <vector>

#include <misc/edit_distance_matrix.hpp>

/** @brief a CIGAR operation: a match (=), a mismatch (X, substitution), an insertion (I, a pattern symbol absent
 *         from the text) or a deletion (D, a text symbol absent from the pattern) */
enum class cigar_op_t : uint8_t { MATCH, MISMATCH, INS, DEL };

/** @brief one run of a CIGAR. A MATCH run spans len (>= 1) positions. A MISMATCH/INS/DEL run always spans one
 *         position (len == 1) and additionally carries the involved symbol in @ref sym: the reference symbol for a
 *         MISMATCH or DEL, the pattern symbol for an INS. Together with the pattern this lets the reference be
 *         reconstructed and the SAM MD tag be derived. The struct is packed (6 bytes) to keep CIGARs compact. */
#pragma pack(push, 1)
struct cigar_run_t {
    cigar_op_t op;  // the operation
    uint32_t len;   // number of positions the run spans (always 1 for MISMATCH/INS/DEL)
    char sym;       // MISMATCH/DEL: the reference symbol; INS: the pattern symbol; unused (0) for MATCH

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

/** @brief a CIGAR alignment as a sequence of runs (MISMATCH/INS/DEL runs are one position each and carry a symbol) */
using cigar_t = std::vector<cigar_run_t>;

/** @name CIGAR run factories: a MATCH takes a run length; a MISMATCH/INS/DEL takes its single symbol (its length is
 *        1), e.g. MATCH(3), MISMATCH('A'), DEL('C').
 *  @{ */
inline cigar_run_t MATCH(uint32_t len = 1) { return {cigar_op_t::MATCH, len, 0}; }
inline cigar_run_t MISMATCH(char sym)      { return {cigar_op_t::MISMATCH, 1, sym}; }
inline cigar_run_t INS(char sym)           { return {cigar_op_t::INS, 1, sym}; }
inline cigar_run_t DEL(char sym)           { return {cigar_op_t::DEL, 1, sym}; }
/** @} */

/** @brief whether an approximate-pattern-matching locate also computes a CIGAR alignment per occurrence */
enum cigar_mode_t : uint8_t { NO_CIGAR = 0, CIGAR = 1 };

/** @brief appends a CIGAR run, extending the last run only for MATCH (MISMATCH/INS/DEL runs stay one position each,
 *         since each carries its own symbol) */
inline void cigar_push(cigar_t& cigar, cigar_run_t run)
{
    if (run.op == cigar_op_t::MATCH && !cigar.empty() && cigar.back().op == cigar_op_t::MATCH)
        cigar.back().len += run.len;
    else cigar.emplace_back(run);
}

/**
 * @brief computes the edit distance and a CIGAR alignment of the pattern P against the matched string M with a banded
 *        bit-parallel dynamic program (Myers) plus traceback. The returned edit distance is the exact ed(P, M) (which
 *        corrects the relaxed search's possibly-overestimated error), and the CIGAR is expressed relative to M as the
 *        reference (INS = a P-symbol absent from M, DEL = an M-symbol absent from P). Requires ed(P, M) <= k (so the
 *        optimal alignment stays within the band), k <= bp_k_limit_128 and n2 <= bp_max_ref_len.
 * @tparam pos_t unsigned integer type for the distance
 * @tparam seq1_t random-access sequence type of the pattern
 * @tparam seq2_t random-access sequence type of the matched string
 * @param P the pattern
 * @param n1 the length of P
 * @param M the matched string
 * @param n2 the length of M
 * @param k the band half-width (and error bound)
 * @return the pair (ed(P, M), CIGAR of P against M)
 */
// traces back the CIGAR of P against M from the (already materialized) banded matrix, starting at cell (i, j)
// (reference = M in the columns, query = P in the rows). Only cells within the band |i - j| <= k are meaningful; the
// optimal (<= k) path stays inside it, so band-out neighbours are skipped. Reads cells via the matrix' O(1) accessor.
template <typename word_t, typename seq1_t, typename seq2_t>
static cigar_t cigar_traceback(const edit_distance_matrix<word_t>& mat, const seq1_t& P, const seq2_t& M,
                               int64_t i, int64_t j, int64_t k)
{
    auto in_band = [&](int64_t a, int64_t b) { return a >= 0 && b >= 0 && std::abs(a - b) <= k; };
    cigar_t cigar;

    while (i > 0 || j > 0) {
        int64_t cur = mat(i, j);

        if (i > 0 && j > 0 && in_band(i - 1, j - 1) && mat(i - 1, j - 1) + (P[i - 1] != M[j - 1] ? 1 : 0) == cur) {
            cigar_push(cigar, P[i - 1] == M[j - 1] ? MATCH() : MISMATCH((char) M[j - 1])); // reference symbol
            i--; j--;
        } else if (i > 0 && in_band(i - 1, j) && mat(i - 1, j) + 1 == cur) {
            cigar_push(cigar, INS((char) P[i - 1])); // a P-symbol absent from M
            i--;
        } else if (j > 0 && in_band(i, j - 1) && mat(i, j - 1) + 1 == cur) {
            cigar_push(cigar, DEL((char) M[j - 1])); // an M-symbol absent from P
            j--;
        } else if (i > 0) {
            cigar_push(cigar, INS((char) P[i - 1])); // j == 0: only insertions of P remain
            i--;
        } else {
            cigar_push(cigar, DEL((char) M[j - 1])); // i == 0: only deletions of M remain
            j--;
        }
    }

    std::reverse(cigar.begin(), cigar.end());
    return cigar;
}

template <typename pos_t, typename seq1_t, typename seq2_t>
static std::pair<pos_t, cigar_t> edit_cigar(
    const seq1_t& P, int64_t n1, const seq2_t& M, int64_t n2, int64_t k)
{
    assert(k <= bp_k_limit_128);
    assert(n2 <= bp_max_ref_len);

    // an empty side aligns trivially: only insertions of P (M empty) or deletions of M (P empty) remain
    if (n1 == 0 || n2 == 0) {
        cigar_t cigar;
        for (int64_t i = 0; i < n1; i++) cigar.emplace_back(INS((char) P[i])); // M empty: all of P is inserted
        for (int64_t j = 0; j < n2; j++) cigar.emplace_back(DEL((char) M[j])); // P empty: all of M is deleted
        return {pos_t(n1 + n2), std::move(cigar)};
    }

    // builds the banded Myers matrix with reference = M (columns), query = P (rows), then traces back from cell (n1, n2)
    auto run = [&]<typename word_t>() -> std::pair<pos_t, cigar_t> {
        edit_distance_matrix<word_t> mat;
        setup_edit_matrix<word_t>(mat, M, n2, P, n1, k); // reference = M, query = P
        return {mat(n1, n2), cigar_traceback(mat, P, M, n1, n2, k)};
    };

    if (k <= bp_k_limit_64)
         return run.template operator()<uint64_t>();
    else return run.template operator()<uint128_t>();
}

/**
 * @brief aligns the pattern P against its BEST PREFIX of the matched string M -- the prefix M[0..len) minimizing
 *        ed(P, M[0..len)) -- and returns that edit distance, the winning prefix length len and the CIGAR of P against
 *        M[0..len). Combines the best-prefix search (which trims a possibly-overshooting matched string to its true
 *        occurrence, cf. apm_edit::edit_dist_prefix) with the CIGAR traceback into a SINGLE banded DP over M, so the
 *        matched string is scanned only once. The returned distance is the exact ed (correcting the search's relaxed
 *        estimate); the CIGAR is relative to M (INS = a P-symbol absent from M, DEL = an M-symbol absent from P).
 * @tparam pos_t unsigned integer type for the distance / length
 * @tparam seq1_t random-access sequence type of the pattern
 * @tparam seq2_t random-access sequence type of the matched string
 * @param P the pattern
 * @param n1 the length of P
 * @param M the matched string (may be longer than the occurrence; only its best prefix is aligned)
 * @param n2 the length of M
 * @param k the band half-width (and error bound)
 * @return the tuple (ed(P, M[0..len)), len, CIGAR of P against M[0..len)); the CIGAR is empty if the best prefix
 *         exceeds k (the caller drops such a context)
 */
template <typename pos_t, typename seq1_t, typename seq2_t>
static std::tuple<pos_t, pos_t, cigar_t> edit_cigar_prefix(
    const seq1_t& P, int64_t n1, const seq2_t& M, int64_t n2, int64_t k)
{
    assert(k <= bp_k_limit_128);

    // only a prefix of length in [n1 - k, n1 + k] can be within edit distance k of P (a length difference of d costs
    // at least d indels); prefixes past n1 + k leave the band, so the reference is capped there
    const int64_t j_min = n1 > k ? n1 - k : 0;
    const int64_t j_max = std::min<int64_t>(n2, n1 + k);
    assert(j_max <= bp_max_ref_len);

    if (n1 == 0) return {pos_t(0), pos_t(0), cigar_t{}}; // empty pattern: empty prefix, no edits

    auto run = [&]<typename word_t>() -> std::tuple<pos_t, pos_t, cigar_t> {
        edit_distance_matrix<word_t> mat;
        setup_edit_matrix<word_t>(mat, M, j_max, P, n1, k); // reference = M[0..j_max), query = P

        // the best prefix is the column j (in the band at row n1) with the smallest ed(P, M[0..j)); ties keep the
        // shortest length, matching edit_dist_prefix
        pos_t best_err = pos_t(k) + 1;
        int64_t best_len = 0;
        for (int64_t j = j_min; j <= j_max; j++) {
            pos_t e = mat(n1, j); // = ed(P, M[0..j))
            if (e < best_err) { best_err = e; best_len = j; }
        }

        if (best_err > pos_t(k)) return {best_err, pos_t(best_len), cigar_t{}}; // > k window: caller drops it
        return {best_err, pos_t(best_len), cigar_traceback(mat, P, M, n1, best_len, k)};
    };

    if (k <= bp_k_limit_64)
         return run.template operator()<uint64_t>();
    else return run.template operator()<uint128_t>();
}

/**
 * @brief counts the edits of a CIGAR while checking it is a valid alignment of S1 against S2: every MATCH must sit
 *        on equal symbols, every MISMATCH on differing symbols, and the operations together must consume exactly
 *        all of S1 and all of S2
 * @tparam sym_t symbol type of the pattern
 * @tparam sym2_t symbol type of the matched string
 * @param S1 pointer to the pattern symbols
 * @param n1 the length of S1
 * @param S2 pointer to the matched-string symbols
 * @param n2 the length of S2
 * @param cigar the CIGAR to validate
 * @return the number of edits (mismatches + insertions + deletions) if the CIGAR validly aligns S1 against S2, else -1
 */
template <typename sym_t, typename sym2_t>
static int64_t num_cigar_edits(const sym_t* S1, uint64_t n1, const sym2_t* S2, uint64_t n2, const cigar_t& cigar)
{
    uint64_t p1 = 0, p2 = 0;
    int64_t edits = 0;

    for (const cigar_run_t& run : cigar) {
        for (uint32_t r = 0; r < run.len; r++) {
            switch (run.op) {
                case cigar_op_t::MATCH:    if (p1 >= n1 || p2 >= n2 || S1[p1] != S2[p2]) return -1; p1++; p2++; break;
                case cigar_op_t::MISMATCH: if (p1 >= n1 || p2 >= n2 || S1[p1] == S2[p2] || (char) S2[p2] != run.sym) return -1; p1++; p2++; edits++; break;
                case cigar_op_t::INS:      if (p1 >= n1 || (char) S1[p1] != run.sym) return -1; p1++; edits++; break;
                case cigar_op_t::DEL:      if (p2 >= n2 || (char) S2[p2] != run.sym) return -1; p2++; edits++; break;
                default:                   __builtin_unreachable();
            }
        }
    }

    return (p1 == n1 && p2 == n2) ? edits : -1;
}

/**
 * @brief writes the SAM CIGAR string of a CIGAR to out, e.g. "5=1X2I3="; consecutive runs of the same operation are
 *        merged (each MISMATCH/INS/DEL run is a single position, so two adjacent mismatches print as "2X"). The
 *        extended operations = and X (rather than M) are used, since the alignment distinguishes matches from
 *        mismatches.
 * @param out the stream the CIGAR is written to
 * @param cigar the CIGAR alignment
 */
inline void append_cigar(std::ostream& out, const cigar_t& cigar)
{
    for (size_t i = 0; i < cigar.size();) {
        cigar_op_t op = cigar[i].op;
        uint64_t len = 0;
        size_t j = i;
        while (j < cigar.size() && cigar[j].op == op) len += cigar[j++].len;
        out << len << cigar[i].code();
        i = j;
    }
}

/**
 * @brief writes the SAM MD-tag value of a CIGAR (the alignment of a read against the reference) to out: runs of
 *        matching bases as counts, each mismatch as the reference base, and each maximal run of deletions as '^'
 *        followed by the deleted reference bases (insertions are not part of the MD tag). E.g. "10=1X5=2D6=" over
 *        reference bases G / TA yields "10G5^TA6". Uses the reference symbols stored in the CIGAR (only present in a
 *        CIGAR built by cigar_traceback or the hamming mismatch capture).
 * @param out the stream the MD-tag value is written to
 * @param cigar the CIGAR alignment
 */
inline void append_md(std::ostream& out, const cigar_t& cigar)
{
    uint64_t matches = 0; // current run of matching (or, via '=', reference-equal) bases
    bool in_del = false;  // whether the previous emitted op was a deletion (to group a '^'-run)

    for (const cigar_run_t& run : cigar) {
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
