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
 * @file apm.hpp
 * @brief approximate pattern matching: the occurrences, the query parameters, the naive search with Hamming and edit
 *        distance and the filters of edit-distance occurrences
 */

#ifndef MOVE_R_MISC_APM_HPP
#define MOVE_R_MISC_APM_HPP

#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

#include <misc/utils.hpp>
#include <misc/memory.hpp>
#include <misc/portable_layout.hpp>
#include <misc/hash.hpp>
#include <misc/log.hpp>
#include <misc/files.hpp>
#include <misc/search.hpp>
#include <misc/search_schemes.hpp>
#include <data_structures/edit_distance_matrix.hpp>
#include <misc/cigar.hpp>

/** @brief whether an approximate-pattern-matching query only counts occurrences or also locates them */
enum query_support_t : uint8_t {
    COUNT = 0,
    LOCATE = 1
};

/**
 * @brief whether an enumeration of the right extensions of a pattern reports the terminator symbol; a left
 *        extension with the terminator is never reported
 */
enum extension_scope_t : uint8_t {
    NO_TERM = 0, // all symbols of the text except the terminator
    TERM = 1     // all symbols of the text and the terminator (as the smallest symbol, reported as TERMINATOR_SYMBOL)
};

// bit-field widths of a packed aprx_occ_t: position, length and error count (10 bytes together)
static constexpr uint32_t aprx_occ_pos_bits = 48; // positions up to 2^48 - 1 (see uint48_t)
static constexpr uint32_t aprx_occ_len_bits = 24; // lengths up to 2^24 - 1 (an occurrence has at most m + k symbols)
static constexpr uint32_t aprx_occ_err_bits = 8;  // error counts up to 2^8 - 1 (an occurrence has at most k errors)

// the largest length of an approximate occurrence (see aprx_occ_t)
static constexpr uint64_t aprx_occ_max_len = (uint64_t(1) << aprx_occ_len_bits) - 1;

/**
 * @brief throws std::invalid_argument if an occurrence can be longer than aprx_occ_max_len
 * @param max_len the largest length of an occurrence (m for the Hamming distance, m + k for the edit distance)
 */
inline void require_storable_occurrence_length(uint64_t max_len)
{
    if (max_len <= aprx_occ_max_len) return;
    throw std::invalid_argument("the occurrences of a pattern can be " + std::to_string(max_len) + " symbols long "
        "(its length, plus the number of errors with the edit distance); at most " + std::to_string(aprx_occ_max_len)
        + " are supported");
}

/**
 * @brief an approximate occurrence: its starting position, length and error count, in bit-fields with the layout of a
 *        packed GCC struct (see packed_bit_fields)
 *
 * With a CIGAR, the occurrence also stores the index cig_idx of its CIGAR in the vector of CIGARs returned by locate;
 * all occurrences of an SA interval have the same CIGAR. An occurrence has 10 bytes without and 14 bytes with a CIGAR
 * index.
 * @tparam mode whether the occurrence also stores the index of its CIGAR
 */
template <cigar_mode_t mode = NO_CIGAR>
struct __attribute__((packed)) aprx_occ_t {
    using pos_field = next_bit_field_t<void, aprx_occ_pos_bits>;
    using len_field = next_bit_field_t<pos_field, aprx_occ_len_bits>;
    using err_field = next_bit_field_t<len_field, aprx_occ_err_bits>;

    packed_bit_fields<bytes_of_bit_fields<err_field>> bits; // pos, len and err
    STRING_UTILS_NO_UNIQUE_ADDRESS std::conditional_t<mode == CIGAR, uint32_t, std::monostate> cig_idx = {};

    aprx_occ_t() = default;

    aprx_occ_t(uint64_t pos, uint64_t len, uint64_t err)
    { bits.template assign<pos_field, len_field, err_field>(pos, len, err); }

    aprx_occ_t(uint64_t pos, uint64_t len, uint64_t err, uint32_t cigar_idx) requires(mode == CIGAR)
        : aprx_occ_t(pos, len, err)
    { cig_idx = cigar_idx; }

    uint64_t pos() const { return bits.template get<pos_field>(); } // the starting position
    uint64_t len() const { return bits.template get<len_field>(); } // the length
    uint64_t err() const { return bits.template get<err_field>(); } // the error count

    void set_pos(uint64_t value) { bits.template set<pos_field>(value); }
    void set_len(uint64_t value) { bits.template set<len_field>(value); }
    void set_err(uint64_t value) { bits.template set<err_field>(value); }

    /** @brief orders occurrences by position, then by error count, then by length (the CIGAR is ignored) */
    bool operator<(const aprx_occ_t& other) const
    {
        return pos() != other.pos() ? pos() < other.pos() :
              (err() != other.err() ? err() < other.err() :
                                      len() < other.len());
    }

    bool operator==(const aprx_occ_t& other) const
    {
        return pos() == other.pos() &&
               len() == other.len() &&
               err() == other.err();
    }

    /** @brief writes an occurrence to an output stream as "(pos, len, err)" */
    friend std::ostream& operator<<(std::ostream& os, const aprx_occ_t& occ) {
        return os << "(" << occ.pos() << ", " << occ.len() << ", " << occ.err() << ")";
    }
};

// IntelliSense (EDG) lays out empty [[no_unique_address]] members of packed structs differently than GCC and clang, so
// the size check is only done by the compiler
#ifndef __INTELLISENSE__
static_assert(sizeof(aprx_occ_t<NO_CIGAR>) == (aprx_occ_pos_bits + aprx_occ_len_bits + aprx_occ_err_bits) / 8);
#endif
static_assert(sizeof(aprx_occ_t<CIGAR>) == sizeof(aprx_occ_t<NO_CIGAR>) + sizeof(uint32_t));

// key function for radix_sort: returns the starting position of an occurrence (see aprx_occ_t)
struct key_of_pos {
    template <typename occ_t>
    inline auto operator()(const occ_t& occurrence) const { return occurrence.pos(); }
};

/**
 * @brief computes the Hamming distance of two equally long sequences; the computation stops once the distance exceeds k
 * @tparam inp_t contiguous range type of the sequences (e.g. std::string, std::vector, std::span)
 * @param str_1 the first sequence
 * @param str_2 the second sequence (as long as str_1)
 * @param k the largest distance that is computed exactly
 * @return the Hamming distance of str_1 and str_2, or a value > k if it exceeds k
 */
template <std::ranges::contiguous_range inp_t>
inline uint64_t hamming_dist_bounded(const inp_t& str_1, const inp_t& str_2, uint64_t k)
{
    using sym_t = std::ranges::range_value_t<inp_t>;
    assert(str_1.size() == str_2.size());
    const uint64_t n = str_1.size();
    uint64_t dist = 0;

    /* for a byte alphabet, compare the symbols of one word per step: XOR two words, then count the mismatching (nonzero)
       bytes by OR-reducing the bits of each byte into its low bit and popcounting */
    if constexpr (sizeof(sym_t) == 1) {
        constexpr uint64_t syms_per_word = sizeof(uint64_t);
        const auto* b1 = reinterpret_cast<const uint8_t*>(std::ranges::data(str_1));
        const auto* b2 = reinterpret_cast<const uint8_t*>(std::ranges::data(str_2));
        uint64_t i = 0;

        for (; i + syms_per_word <= n && dist <= k; i += syms_per_word) {
            uint64_t w1, w2;
            std::memcpy(&w1, b1 + i, sizeof(uint64_t));
            std::memcpy(&w2, b2 + i, sizeof(uint64_t));
            uint64_t x = w1 ^ w2;
            x |= x >> 1;
            x |= x >> 2;
            x |= x >> 4;
            x &= 0x0101010101010101ULL; // low bit of each byte set iff that byte differed
            dist += ::popcount(x);
        }

        for (; i < n && dist <= k; i++) dist += b1[i] != b2[i];
    } else {
        const auto* p1 = std::ranges::data(str_1);
        const auto* p2 = std::ranges::data(str_2);

        for (uint64_t i = 0; i < n && dist <= k; i++) dist += p1[i] != p2[i];
    }

    return dist;
}

/**
 * @brief locates all occurrences of P in T with at most k mismatches (Hamming distance) by comparing P with every
 *        window of T
 * @tparam inp_t contiguous range type of the text and the pattern
 * @param T the text
 * @param P the pattern
 * @param k the maximum number of mismatches
 * @return the occurrences of P in T with at most k mismatches, in increasing order of their positions
 */
template <std::ranges::contiguous_range inp_t>
inline std::vector<aprx_occ_t<>> locate_hamming_dist(const inp_t& T, const inp_t& P, uint64_t k)
{
    using sym_t = std::ranges::range_value_t<inp_t>;
    std::span<const sym_t> Ts(T), Ps(P); // spans of T and P (a window of T is a subspan of Ts)
    uint64_t m = Ps.size();
    uint64_t n = Ts.size();

    std::vector<aprx_occ_t<>> Occ;

    for (uint64_t i = 0; i + m <= n; i++) {
        uint64_t dist = hamming_dist_bounded(Ts.subspan(i, m), Ps, k);
        if (dist <= k) Occ.emplace_back(aprx_occ_t<>(i, m, dist));
    }

    return Occ;
}

/**
 * @brief computes the edit distance of two sequences with a banded bit-parallel dynamic program (Myers), or k + 1 if it
 *        exceeds k; requires k <= edit_distance_max_errors and at most edit_distance_max_length symbols in the longer
 *        sequence
 * @tparam inp_t contiguous range type of the sequences (a byte alphabet, since the matrix indexes symbols by byte)
 * @param T the first sequence
 * @param P the second sequence
 * @param k the largest distance that is computed exactly
 * @return the edit distance of T and P, or k + 1 if it exceeds k
 */
template <std::ranges::contiguous_range inp_t>
inline int64_t edit_dist_bounded(const inp_t& T, const inp_t& P, int64_t k)
{
    const int64_t nT = T.size(), nP = P.size();
    if (std::abs(nT - nP) > k) return k + 1;

    // the longer sequence R is the input (columns) of the matrix, the shorter sequence Q its query (rows)
    const bool p_longer = nP > nT;
    const int64_t nR = p_longer ? nP : nT, nQ = p_longer ? nT : nP;
    if (nQ == 0) return nR;
    assert(k <= edit_distance_max_errors);
    assert(nR <= edit_distance_max_length);

    auto run = [&]<typename word_t>() -> int64_t {
        edit_distance_matrix<word_t> mat;
        if (p_longer) setup_edit_matrix<word_t>(mat, P, nP, T, nT, k);
        else          setup_edit_matrix<word_t>(mat, T, nT, P, nP, k);
        return std::min<int64_t>(mat(nQ, nR), k + 1);
    };

    return with_edit_distance_word(k, run);
}

/**
 * @brief finds the prefix M[0..len) of M with the fewest errors against P, where max(1, m - k) <= len and
 *        len <= min(|M|, m + k, l_max); on ties, the shortest prefix
 * @tparam word_t word type of the edit-distance matrix (chosen by the caller to fit k)
 * @tparam inp_t contiguous range type of M (a byte alphabet, since the matrix indexes symbols by byte)
 * @param mat the banded edit-distance matrix, with P already set as its input (columns)
 * @param M the string whose prefixes are aligned against P (the query of mat)
 * @param m the pattern length
 * @param k the maximum number of errors
 * @param first_col the distances in the first column of mat (see edit_distance_matrix::init; empty: 0, 1, 2, ...)
 * @param l_max the largest prefix length that is considered
 * @return the pair (e, len) of the fewest errors e and the shortest prefix length len with e errors, or (k + 1, 0) if
 *         every considered prefix has more than k errors
 */
template <typename word_t, std::ranges::contiguous_range inp_t>
std::pair<uint64_t, uint64_t> edit_dist_prefix(edit_distance_matrix<word_t>& mat, const inp_t& M, uint64_t m,
    uint64_t k, const std::vector<uint16_t>& first_col = {}, uint64_t l_max = std::numeric_limits<uint64_t>::max())
{
    mat.init(k, first_col);
    const uint64_t l_min = m > k ? m - k : 0;
    const uint64_t i_max = std::min<uint64_t>({uint64_t(M.size()), m + k, l_max});
    const uint64_t final_col = mat.num_cols() - 1;
    uint64_t best_err = k + 1, best_len = 0;

    for (uint64_t i = 1; i <= i_max; i++) {
        const bool alive = mat.compute_row(i, M[i - 1]);

        if (i >= l_min) {
            uint64_t err = mat(i, final_col);
            if (err < best_err) { best_err = err; best_len = i; }
        }

        if (!alive) break;
    }

    return {best_err, best_len};
}

/**
 * @brief computes edit_dist_prefix for the alignments that do not start with a deletion, i.e. that insert P[0] or align
 *        it to M[0] (see first_column_without_deletion); requires k >= 1
 * @param mat the banded edit-distance matrix, with P[1..m) already set as its input (columns)
 * @param first_matches whether P[0] matches M[0]
 * @see edit_dist_prefix for the other parameters and the result
 */
template <typename word_t, std::ranges::contiguous_range inp_t>
std::pair<uint64_t, uint64_t> edit_dist_anchored_prefix(edit_distance_matrix<word_t>& mat, const inp_t& M, uint64_t m,
    uint64_t k, bool first_matches)
{
    // an alignment against a prefix of length len > m has at least len - m errors, and len - m + 1 if P[0] does not
    // match M[0]
    return edit_dist_prefix(mat, M, m, k, first_column_without_deletion(first_matches), m + k - (first_matches ? 0 : 1));
}

/**
 * @brief locates all occurrences of P in T with at most k errors (edit distance) with a banded bit-parallel matrix
 *        (Myers) whose input (columns) is P: for every starting position i of T, the shortest window T[i..i+len) with
 *        the fewest errors against P is reported if it has at most k errors
 * @tparam inp_t contiguous range type of the text and the pattern (a byte alphabet, since the matrix indexes symbols by
 *         byte)
 * @param T the text
 * @param P the pattern (1 <= |P| <= edit_distance_max_length)
 * @param k the maximum number of errors (k <= edit_distance_max_errors)
 * @return the occurrences of P in T with at most k errors, in increasing order of their positions
 */
template <std::ranges::contiguous_range inp_t>
inline std::vector<aprx_occ_t<>> locate_edit_dist(const inp_t& T, const inp_t& P, uint64_t k)
{
    static_assert(sizeof(std::ranges::range_value_t<inp_t>) == 1);
    const int64_t n = T.size(), m = P.size();
    assert(m > 0 && m <= edit_distance_max_length && int64_t(k) <= edit_distance_max_errors);
    const int64_t l_min = std::max<int64_t>(1, m - int64_t(k)); // occurrences have positive length (no empty window)
    const int64_t l_max = m + k;

    auto run = [&]<typename word_t>() {
        edit_distance_matrix<word_t> mat;
        std::vector<uint8_t> alph_map;
        uint16_t sigma = build_alphabet_map(P, m, T, n, alph_map);
        directional_substring<inp_t> P_col(P, 0, m - 1, RIGHT);
        mat.set_input(P_col, sigma, alph_map);

        std::vector<aprx_occ_t<>> Occ;

        for (int64_t i = 0; i < n; i++) {
            mat.init(k);
            aprx_occ_t<> occ(uint64_t(i), 0, uint64_t(k + 1));
            const int64_t j_max = std::min(l_max, n - i);

            for (int64_t j = 1; j <= j_max; j++) {
                // once all cells of the band exceed k, every longer window has more than k errors
                if (!mat.compute_row(j, T[i + j - 1])) [[unlikely]] break;

                // a cell of the final column contains ed(T[i..i+j), P), which is exact if it is at most k
                if (j >= l_min && mat.is_in_final_column(j)) {
                    int64_t err = mat(j, m);

                    if (err <= int64_t(k) && (err < int64_t(occ.err()) ||
                       (err == int64_t(occ.err()) && j < int64_t(occ.len())))
                    ) {
                        occ.set_err(err);
                        occ.set_len(j);
                    }
                }
            }

            if (occ.err() <= k) Occ.emplace_back(occ);
        }

        return Occ;
    };

    return with_edit_distance_word(k, run);
}

/**
 * @brief removes redundant approximate occurrences in place: an occurrence is removed if the previous kept occurrence
 *        and the next occurrence both have at most as many errors as it and their positions differ by at most 2k + 3
 * @tparam mode whether the occurrences store a CIGAR index (moved along with the kept occurrences)
 * @param Occ the occurrences, sorted by position; filtered in place
 * @param k the maximum number of errors
 */
template <cigar_mode_t mode = NO_CIGAR>
inline void filter_edit_distance_occurrences(std::vector<aprx_occ_t<mode>>& Occ, uint64_t k)
{
    if (Occ.size() <= 2) return;
    uint64_t write_idx = 0;

    for (uint64_t read_idx = 0; read_idx < Occ.size(); read_idx++) {
        auto occ = std::move(Occ[read_idx]);

        // remove the last kept occurrence (mid) while it lies between two occurrences that cover it (prev and occ)
        while (write_idx >= 2) {
            const auto& mid = Occ[write_idx - 1];
            const auto& prev = Occ[write_idx - 2];

            if (occ.pos() - prev.pos() <= 2 * k + 3 &&
                prev.err() <= mid.err() &&
                occ.err() <= mid.err()
            ) {
                write_idx--;
            } else {
                break;
            }
        }

        Occ[write_idx++] = std::move(occ);
    }

    Occ.resize(write_idx);
}

/**
 * @brief returns whether every occurrence in occ_all is covered by an occurrence in occ_filtered, i.e. by one with at
 *        most as many errors at most 2k + 1 positions away (so the filtering removed no required occurrence)
 * @param occ_filtered the filtered occurrences, sorted by position
 * @param occ_all all occurrences, sorted by position
 * @param k the maximum number of errors
 */
inline bool verify_edit_distance_coverage(
    const std::vector<aprx_occ_t<>>& occ_filtered,
    const std::vector<aprx_occ_t<>>& occ_all,
    uint64_t k
) {
    uint64_t max_dist = 2 * k + 1;
    uint64_t left = 0, right = 0;

    for (const auto& o : occ_all) {
        while (right < occ_filtered.size() && occ_filtered[right].pos() <= o.pos() + max_dist) right++;
        while (left < right && occ_filtered[left].pos() + max_dist < o.pos()) left++;
        if (left == right) return false;
        bool valid = false;

        for (uint64_t i = left; i < right; i++) {
            if (occ_filtered[i].err() <= o.err()) {
                valid = true;
                break;
            }
        }

        if (!valid) return false;
    }

    return true;
}

/**
 * @brief keeps a subset of the occurrences such that every occurrence has a kept one with at most as many errors at
 *        most max_dist positions away: the error counts are processed in increasing order, and the leftmost occurrence
 *        of the current error count that is not covered yet is covered by keeping the rightmost occurrence at most
 *        max_dist positions after it with at most as many errors
 * @tparam mode whether the occurrences store a CIGAR index (moved along with the kept occurrences)
 * @param Occ the occurrences, sorted by position; filtered in place
 * @param max_dist the maximum distance between an occurrence and a kept occurrence covering it
 */
template <cigar_mode_t mode = NO_CIGAR>
inline void select_covering_occurrences(std::vector<aprx_occ_t<mode>>& Occ, uint64_t max_dist)
{
    enum : uint8_t { UNCOVERED, COVERED, KEPT };
    constexpr uint64_t num_errs = uint64_t(1) << aprx_occ_err_bits;
    const uint64_t n = Occ.size();

    // whether each occurrence is uncovered, covered or kept: [0..n-1] state[i] is the state of Occ[i]
    std::vector<uint8_t> state(n, UNCOVERED);

    /* the indices in Occ of the occurrences, sorted by error count, then by position: [0..n-1] by_err[x] is the index of
       the x-th occurrence; the occurrences with e errors are by_err[err_beg[e]..err_beg[e + 1] - 1] */
    std::array<uint64_t, num_errs + 1> err_beg {};
    for (const auto& occ : Occ) err_beg[occ.err() + 1]++;
    for (uint64_t err = 0; err < num_errs; err++) err_beg[err + 1] += err_beg[err];
    std::array<uint64_t, num_errs + 1> err_end = err_beg;
    no_init_vector<uint64_t> by_err;
    resize_no_init(by_err, n);
    for (uint64_t i = 0; i < n; i++) by_err[err_end[Occ[i].err()]++] = i;

    for (uint64_t err = 0; err < num_errs; err++) {
        // the index of the last occurrence at most max_dist positions after Occ[i] (it only increases, since i does)
        uint64_t hi = 0;

        for (uint64_t x = err_beg[err]; x < err_beg[err + 1]; x++) {
            const uint64_t i = by_err[x];
            if (state[i] != UNCOVERED) continue;
            hi = std::max(hi, i);
            while (hi + 1 < n && Occ[hi + 1].pos() <= Occ[i].pos() + max_dist) hi++;

            // keep the last occurrence Occ[j], j in [i..hi], with at most err errors, and mark the occurrences it covers
            uint64_t j = hi;
            while (Occ[j].err() > err) j--;
            state[j] = KEPT;

            for (uint64_t t = j; t-- > 0 && Occ[j].pos() - Occ[t].pos() <= max_dist;) {
                if (state[t] == UNCOVERED && Occ[t].err() >= Occ[j].err()) state[t] = COVERED;
            }

            for (uint64_t t = j + 1; t < n && Occ[t].pos() - Occ[j].pos() <= max_dist; t++) {
                if (state[t] == UNCOVERED && Occ[t].err() >= Occ[j].err()) state[t] = COVERED;
            }
        }
    }

    uint64_t write_idx = 0;
    for (uint64_t i = 0; i < n; i++) if (state[i] == KEPT) Occ[write_idx++] = std::move(Occ[i]);
    Occ.resize(write_idx);
}

/** @brief the distance metric used for approximate pattern matching */
enum distance_metric_t : int8_t {
    HAMMING_DISTANCE = 0,
    EDIT_DISTANCE = 1
};

/**
 * @brief locates all occurrences of P in T with at most k errors with the given distance metric
 * @tparam inp_t contiguous range type of the text and pattern
 * @param T the text
 * @param P the pattern
 * @param k the maximum number of errors
 * @param metric the distance metric
 * @return the occurrences of P in T with at most k errors
 */
template <std::ranges::contiguous_range inp_t>
inline std::vector<aprx_occ_t<>> locate_naive(const inp_t& T, const inp_t& P, uint64_t k,
    distance_metric_t metric)
{
    assert(metric == HAMMING_DISTANCE || metric == EDIT_DISTANCE);
    if (metric == HAMMING_DISTANCE) return locate_hamming_dist(T, P, k);
    return locate_edit_dist(T, P, k);
}

// ############################# QUERY PARAMETERS #############################

/** @brief strategy for partitioning a pattern into the parts of a search scheme */
enum partitioning_t : uint8_t {
    _partition_uniform = 0,  // p parts of equal length
    _partition_dynamic = 1,  // part boundaries chosen per pattern such that each part has few exact occurrences
    _partition_optimized = 2 // part boundaries and the direction of the search scheme chosen per pattern such that the
                             // searches have a small predicted number of steps (see apm_step_prediction)
};

class apm_costs_t;

/** @brief the options of an approximate-pattern-matching query that are not part of the search scheme */
struct apm_params_t {
    // how P is partitioned into the parts of the search scheme
    partitioning_t partitioning = _partition_optimized;

    /* true <=> a search whose first step has the upper bound 0 starts at the search context of its first part, which
       the partitioning computes (if the part has one, see apm_partition::has_context); false <=> the search matches
       its first part symbol by symbol */
    bool start_at_part_contexts = true;

    /* the error model compiled for the pattern (see apm_costs_t), or nullptr. Without a model, each error has cost 1,
       and an alignment is valid if it has at most k errors (the k of the search scheme) and does not start or end
       with a deletion. With a model, the k of the search scheme must be at least costs->max_errors(), and each
       occurrence is reported with its cost as its error count. With edit distance and a model, the occurrences are
       reported in increasing order of their positions, and for every valid alignment that starts at a position i, an
       occurrence with at most the same cost is reported that starts at most 2 * costs->max_errors() + 1 positions
       away from i */
    const apm_costs_t* costs = nullptr;

    /* with an error model: true <=> apm_model_plan may lower the upper bound of a step of a search to the largest
       number of errors that a valid alignment can have in the parts matched up to this step, if this lowers the
       predicted number of steps; false <=> the bounds of the search scheme are used */
    bool lower_upper_bounds = true;

    /* with an error model: true <=> apm_model_plan runs the anchored search instead of the searches of the scheme if
       the anchored search has a smaller predicted number of steps; false <=> the searches of the scheme are run */
    bool anchored_search = true;
};

// ############################# ALL-BEST MAPPING #############################

/* locate_best increases the error bound k by bestmap_small_step for k < bestmap_large_k, else by bestmap_large_step
   (plus the number of further error counts that are reported); the steps of Columba (Renders et al.) */
inline constexpr uint64_t bestmap_large_k = 5;
inline constexpr uint64_t bestmap_small_step = 2;
inline constexpr uint64_t bestmap_large_step = 4;

/**
 * @brief locates the occurrences of a pattern with at most best + x errors, where best is the minimum number
 *        of errors of an occurrence with at most max_k errors (all-best mapping)
 *
 * locate_k is called with increasing error bounds k, starting with k = 0, until it reports an occurrence;
 * k is increased by x + 2 for k < 5 and by x + 4 for k >= 5 (see bestmap_large_k). If then
 * best + x > k, locate_k is called once more with k = min(best + x, max_k).
 *
 * @tparam mode whether the occurrences have a CIGAR alignment
 * @tparam locate_fnc_t type of locate_k
 * @param max_k the maximum number of errors
 * @param x the number of error counts above best that are also reported
 * @param occ (output) the occurrences; the vector is cleared before each call of locate_k
 * @param locate_k function called as locate_k(k, occ), which appends all occurrences with at most k errors to occ
 * @return best, or max_k + 1 if there is no occurrence with at most max_k errors
 */
template <cigar_mode_t mode = NO_CIGAR, typename locate_fnc_t>
uint64_t locate_best(uint64_t max_k, uint64_t x, std::vector<aprx_occ_t<mode>>& occ, locate_fnc_t locate_k)
{
    for (uint64_t k = 0;;) {
        occ.clear();
        locate_k(k, occ);

        if (!occ.empty()) {
            uint64_t best = occ.front().err();
            for (const aprx_occ_t<mode>& o : occ) best = std::min<uint64_t>(best, o.err());

            // the occurrences with more than k and at most best + x errors have not been located yet
            if (best + x > k && k < max_k) { k = std::min(best + x, max_k); continue; }

            std::erase_if(occ, [&](const aprx_occ_t<mode>& o) { return o.err() > best + x; });
            return best;
        }

        if (k >= max_k) return max_k + 1;
        k = std::min(max_k, k + x + (k < bestmap_large_k ? bestmap_small_step : bestmap_large_step));
    }
}

#endif // MOVE_R_MISC_APM_HPP
