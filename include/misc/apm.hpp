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

#include <cstdint>
#include <cstring>
#include <ranges>
#include <span>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

#include <tsl/sparse_set.h>

#include <misc/utils.hpp>
#include <misc/hash.hpp>
#include <misc/log.hpp>
#include <misc/files.hpp>
#include <misc/search.hpp>
#include <misc/strings.hpp>
#include <misc/search_schemes.hpp>
#include <misc/directional_substring.hpp>
#include <misc/edit_distance_matrix.hpp>
#include <misc/cigar.hpp>

/** @brief whether an approximate-pattern-matching query only counts occurrences or also locates them */
enum query_support_t : uint8_t {
    COUNT = 0,
    LOCATE = 1
};

/**
 * @brief an approximate occurrence: its starting position, length and error count. When mode == CIGAR it also
 *        owns the occurrence's CIGAR alignment.
 * @tparam pos_t unsigned integer position type
 * @tparam mode whether the occurrence also carries a CIGAR alignment
 */
template <typename pos_t, cigar_mode_t mode = NO_CIGAR>
struct aprx_occ_t {
    pos_t pos;
    pos_t len;
    pos_t err;
    [[no_unique_address]] std::conditional_t<mode == CIGAR, cigar_t, std::monostate> cigar;

    /** @brief orders occurrences by position, then by error count, then by length (the CIGAR is ignored) */
    bool operator<(const aprx_occ_t& other) const
    {
        return pos != other.pos ? pos < other.pos :
              (err != other.err ? err < other.err :
                                  len < other.len);
    }

    bool operator==(const aprx_occ_t& other) const
    {
        return pos == other.pos &&
               len == other.len &&
               err == other.err;
    }

    /** @brief writes an occurrence to an output stream as "(pos, len, err)" */
    friend std::ostream& operator<<(std::ostream& os, const aprx_occ_t& occ) {
        return os << "(" << occ.pos << ", " << occ.len << ", " << occ.err << ")";
    }
};

/**
 * @brief computes the hamming distance of two equally long sequences, stopping early once it exceeds k
 * @tparam pos_t unsigned integer type for the distance
 * @tparam inp_t contiguous range type of the sequences (e.g. std::string, std::vector, std::span)
 * @param str_1 the first sequence
 * @param str_2 the second sequence (must have the same length as str_1)
 * @param k the distance threshold
 * @return the hamming distance of str_1 and str_2 (or a value > k if it exceeds k)
 */
template <typename pos_t, std::ranges::contiguous_range inp_t>
static pos_t hamming_dist_bounded(const inp_t& str_1, const inp_t& str_2, pos_t k)
{
    using sym_t = std::ranges::range_value_t<inp_t>;
    assert(str_1.size() == str_2.size());
    const uint64_t n = str_1.size();
    pos_t dist = 0;

    // for a byte alphabet, compare eight symbols per step: XOR two 8-byte words, then count the
    // mismatching (nonzero bytes by OR-reducing each byte's bits into its low bit and popcounting
    if constexpr (sizeof(sym_t) == 1) {
        const auto* b1 = reinterpret_cast<const uint8_t*>(std::ranges::data(str_1));
        const auto* b2 = reinterpret_cast<const uint8_t*>(std::ranges::data(str_2));
        uint64_t i = 0;

        for (; i + 8 <= n && dist <= k; i += 8) {
            uint64_t w1, w2;
            std::memcpy(&w1, b1 + i, 8);
            std::memcpy(&w2, b2 + i, 8);
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
 * @brief naively locates all occurrences of P in T with at most k mismatches (w.r.t. hamming distance)
 * @tparam pos_t unsigned integer position type
 * @tparam inp_t contiguous range type of the text and pattern
 * @param T the text
 * @param P the pattern
 * @param k the maximum number of mismatches
 * @return the occurrences of P in T with at most k mismatches
 */
template <typename pos_t, std::ranges::contiguous_range inp_t>
static std::vector<aprx_occ_t<pos_t>> locate_hamming_dist(const inp_t& T, const inp_t& P, pos_t k)
{
    using sym_t = std::ranges::range_value_t<inp_t>;
    std::span<const sym_t> Ts(T), Ps(P); // read-only views for the windowed distance calls
    pos_t m = Ps.size();
    pos_t n = Ts.size();

    std::vector<aprx_occ_t<pos_t>> Occ;

    for (pos_t i = 0; i + m <= n; i++) {
        pos_t dist = hamming_dist_bounded<pos_t>(Ts.subspan(i, m), Ps, k);
        if (dist <= k) Occ.emplace_back(aprx_occ_t<pos_t>{.pos = i, .len = m, .err = dist});
    }

    return Occ;
}

template <std::ranges::contiguous_range ref_t>
static uint16_t build_alph_map(const ref_t& R, int64_t nR, std::vector<uint8_t>& alph_map)
{
    alph_map.assign(256, 0);
    uint16_t sigma = 0;

    for (int64_t j = 0; j < nR; j++) {
        uint8_t b = sym_to_uchar(R[j]);
        if (alph_map[b] == 0) alph_map[b] = ++sigma;
    }

    return sigma;
}

template <typename word_t, std::ranges::contiguous_range ref_t, std::ranges::contiguous_range qry_t>
static void setup_edit_matrix(edit_distance_matrix<word_t>& mat, const ref_t& R, int64_t nR, const qry_t& Q, int64_t nQ, int64_t k)
{
    std::vector<uint8_t> alph_map;
    uint16_t sigma = build_alph_map(R, nR, alph_map);
    directional_substring<int64_t, ref_t> R_col(R, 0, nR - 1, RIGHT);
    mat.set_input(R_col, sigma, alph_map);
    mat.init(k);
    for (int64_t i = 1; i <= nQ; i++) mat.compute_row(i, Q[i - 1]);
}

/**
 * @brief computes the edit distance of two sequences with a banded bit-parallel dynamic program (Myers), capped at
 *        k + 1. Requires k <= bp_k_limit_128 and the longer sequence to be at most bp_max_ref_len long.
 * @tparam pos_t unsigned integer type for the distance
 * @tparam inp_t contiguous range type of the sequences (a byte alphabet, since the matrix indexes symbols by byte)
 * @param T the first sequence
 * @param P the second sequence
 * @param k the distance threshold
 * @return the edit distance of T and P, or k + 1 if it exceeds k
 */
template <typename pos_t, std::ranges::contiguous_range inp_t>
static int64_t edit_dist_bounded(const inp_t& T, const inp_t& P, int64_t k)
{
    const int64_t nT = T.size(), nP = P.size();
    if (std::abs(nT - nP) > k) return k + 1;
    const bool p_longer = nP > nT;
    const int64_t nR = p_longer ? nP : nT, nQ = p_longer ? nT : nP;
    if (nQ == 0) return nR;
    assert(k <= bp_k_limit_128);
    assert(nR <= bp_max_ref_len);

    auto run = [&]<typename word_t>() -> int64_t {
        edit_distance_matrix<word_t> mat;
        if (p_longer) setup_edit_matrix<word_t>(mat, P, nP, T, nT, k);
        else          setup_edit_matrix<word_t>(mat, T, nT, P, nP, k);
        return std::min<int64_t>(mat(nQ, nR), k + 1);
    };

    return k <= bp_k_limit_64 ? run.template operator()<uint64_t>() : run.template operator()<__uint128_t>();
}

/**
 * @brief over all prefixes M[0..len) of M (len in the band [m-k, m+k]), finds the one minimizing ed(P, M[0..len))
 * @tparam pos_t unsigned integer type for the distance / length
 * @tparam word_t word type of the edit-distance matrix (chosen by the caller to fit k)
 * @tparam inp_t contiguous range type of the matched string (a byte alphabet, since the matrix indexes symbols by byte)
 * @param mat the banded edit-distance matrix, with the pattern P already set as its input (columns)
 * @param M the matched string whose best pattern-length prefix is sought
 * @param m the pattern length (= mat's number of columns minus one)
 * @param k the error threshold
 * @return (minimum error, shortest length achieving it), or (k + 1, 0) if no prefix stays within k errors
 */
template <typename pos_t, typename word_t, std::ranges::contiguous_range inp_t>
std::pair<pos_t, pos_t> edit_dist_prefix(
    edit_distance_matrix<word_t>& mat, const inp_t& M, pos_t m, pos_t k)
{
    mat.init(k);
    const pos_t n2 = M.size();
    const pos_t l_min = m > k ? m - k : 0;
    const pos_t i_max = std::min<pos_t>(n2, m + k);
    const uint16_t final_col = mat.num_cols() - 1;
    pos_t best_err = k + 1, best_len = 0;

    for (pos_t i = 1; i <= i_max; i++) {
        bool alive = mat.compute_row(i, M[i - 1]);
        if (i >= l_min) {
            pos_t err = mat(i, final_col);
            if (err < best_err) { best_err = err; best_len = i; }
        }
        if (!alive) break;
    }

    return {best_err, best_len};
}

/**
 * @brief locates all occurrences of P in T with at most k errors (w.r.t. edit distance) by running the banded
 *        bit-parallel matrix (Myers; columns = P, set up once) over every starting position of T, keeping for
 *        each starting position the shortest minimum-error occurrence
 * @tparam pos_t unsigned integer position type
 * @tparam inp_t contiguous range type of the text and pattern (a byte alphabet, since the matrix indexes symbols by byte)
 * @param T the text
 * @param P the pattern (1 <= |P| <= bp_max_ref_len)
 * @param k the maximum number of errors (k <= bp_k_limit_128)
 * @return the occurrences of P in T with at most k errors
 */
template <typename pos_t, std::ranges::contiguous_range inp_t>
static std::vector<aprx_occ_t<pos_t>> locate_edit_dist(const inp_t& T, const inp_t& P, pos_t k)
{
    static_assert(sizeof(std::ranges::range_value_t<inp_t>) == 1);
    const int64_t n = T.size(), m = P.size();
    assert(m > 0 && m <= bp_max_ref_len && int64_t(k) <= bp_k_limit_128);
    const int64_t l_min = std::max<int64_t>(1, m - int64_t(k)); // occurrences have positive length (no empty window)
    const int64_t l_max = m + k;

    auto run = [&]<typename word_t>() {
        edit_distance_matrix<word_t> mat;
        std::vector<uint8_t> alph_map;
        uint16_t sigma = build_alph_map(P, m, alph_map);
        directional_substring<int64_t, inp_t> P_col(P, 0, m - 1, RIGHT);
        mat.set_input(P_col, sigma, alph_map);

        std::vector<aprx_occ_t<pos_t>> Occ;

        for (int64_t i = 0; i < n; i++) {
            mat.init(k);
            aprx_occ_t<pos_t> occ{.pos = pos_t(i), .len = 0, .err = pos_t(k + 1)};
            const int64_t j_max = std::min(l_max, n - i);

            for (int64_t j = 1; j <= j_max; j++) {
                // once the whole band exceeds k, no longer window can come back under the budget
                if (!mat.compute_row(j, T[i + j - 1])) [[unlikely]] break;

                // a final-column cell holds ed(T[i..i+j), P) and is exact whenever it is within the budget
                if (j >= l_min && mat.is_in_final_column(j)) {
                    int64_t err = mat(j, m);

                    if (err <= int64_t(k) && (err < int64_t(occ.err) ||
                       (err == int64_t(occ.err) && j < int64_t(occ.len)))
                    ) {
                        occ.err = err;
                        occ.len = j;
                    }
                }
            }

            if (occ.err <= k) Occ.emplace_back(occ);
        }

        return Occ;
    };

    return int64_t(k) <= bp_k_limit_64 ?
        run.template operator()<uint64_t>() :
        run.template operator()<__uint128_t>();
}

/**
 * @brief removes redundant approximate occurrences in place: an occurrence is dropped when a nearby occurrence
 *        (within 4k+3 positions) on each side has an equal or smaller error count
 * @tparam pos_t unsigned integer position type
 * @tparam mode whether the occurrences carry a CIGAR alignment (moved along with the surviving occurrences)
 * @param Occ the occurrences (sorted by position); filtered in place
 * @param k the maximum number of errors
 */
template <typename pos_t, cigar_mode_t mode = NO_CIGAR>
static void filter_edit_distance_occurrences(std::vector<aprx_occ_t<pos_t, mode>>& Occ, pos_t k)
{
    if (Occ.size() <= 2) return;
    pos_t write_idx = 0;

    for (pos_t read_idx = 0; read_idx < Occ.size(); read_idx++) {
        auto occ = std::move(Occ[read_idx]);

        while (write_idx >= 2) {
            const auto& left = Occ[write_idx - 1];
            const auto& right = Occ[write_idx - 2];

            if (occ.pos - right.pos <= 4 * k + 3 &&
                right.err <= left.err &&
                occ.err <= left.err
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
 * @brief verifies that every occurrence in occ_all is covered by an occurrence in occ_filtered with an equal or
 *        smaller error count within 2k+1 positions (i.e. filtering did not drop any required occurrence)
 * @tparam pos_t unsigned integer position type
 * @param occ_filtered the filtered occurrences (sorted by position)
 * @param occ_all all occurrences (sorted by position)
 * @param k the maximum number of errors
 * @return whether occ_filtered covers all occurrences in occ_all
 */
template <typename pos_t>
bool verify_edit_distance_coverage(
    const std::vector<aprx_occ_t<pos_t>>& occ_filtered,
    const std::vector<aprx_occ_t<pos_t>>& occ_all,
    pos_t k
) {
    pos_t max_dist = 2 * k + 1;
    pos_t left = 0, right = 0;

    for (const auto& o : occ_all) {
        while (right < occ_filtered.size() && occ_filtered[right].pos <= o.pos + max_dist) right++;
        while (left < right && occ_filtered[left].pos + max_dist < o.pos) left++;
        if (left == right) return false;
        bool valid = false;

        for (pos_t i = left; i < right; i++) {
            if (occ_filtered[i].err <= o.err) {
                valid = true;
                break;
            }
        }

        if (!valid) return false;
    }

    return true;
}

/** @brief the distance metric used for approximate pattern matching */
enum distance_metric_t : int8_t {
    NO_METRIC = -1,
    HAMMING_DISTANCE = 0,
    EDIT_DISTANCE = 1
};

/**
 * @brief locates all occurrences of P in T with at most k errors w.r.t. the given distance metric
 * @tparam pos_t unsigned integer position type
 * @tparam dist_metr the distance metric (HAMMING_DISTANCE or EDIT_DISTANCE)
 * @tparam inp_t contiguous range type of the text and pattern
 * @param T the text
 * @param P the pattern
 * @param k the maximum number of errors
 * @return the occurrences of P in T with at most k errors
 */
template <typename pos_t, distance_metric_t dist_metr, std::ranges::contiguous_range inp_t>
static std::vector<aprx_occ_t<pos_t>> locate(const inp_t& T, const inp_t& P, pos_t k)
{
    if constexpr (dist_metr == HAMMING_DISTANCE) return locate_hamming_dist<pos_t>(T, P, k);
    if constexpr (dist_metr == EDIT_DISTANCE) return locate_edit_dist<pos_t>(T, P, k);
}