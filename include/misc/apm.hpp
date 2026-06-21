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

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include <tsl/sparse_set.h>

#include <misc/utils.hpp>
#include <misc/hash.hpp>
#include <misc/log.hpp>
#include <misc/files.hpp>
#include <misc/search.hpp>
#include <misc/strings.hpp>
#include <algorithms/apm/search_schemes.hpp>
#include <algorithms/apm/sub_string.hpp>
#include <algorithms/apm/edit_distance_matrix.hpp>

/** @brief whether an approximate-pattern-matching query only counts occurrences or also locates them */
enum query_support_t : uint8_t {
    COUNT = 0,
    LOCATE = 1
};

/**
 * @brief an approximate occurrence: its starting position, length and error count
 * @tparam pos_t unsigned integer position type
 */
template <typename pos_t>
struct aprx_occ_t {
    pos_t pos;
    pos_t len;
    pos_t err;

    /**
     * @brief orders occurrences by position, then by error count, then by length
     * @param other another occurrence
     * @return whether this occurrence is ordered before other
     */
    bool operator<(const aprx_occ_t& other) const
    {
        return pos != other.pos ? pos < other.pos :
              (err != other.err ? err < other.err :
                                  len < other.len);
    }

    bool operator==(const aprx_occ_t& other) const = default;

    /**
     * @brief writes an occurrence to an output stream as "(pos, len, err)"
     * @param os the output stream
     * @param occ the occurrence to write
     * @return the output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const aprx_occ_t& occ) {
        return os << "(" << occ.pos << ", " << occ.len << ", " << occ.err << ")";
    }
};

/**
 * @brief computes the hamming distance of two equally long strings
 * @tparam pos_t unsigned integer type for the distance
 * @param str_1 the first string
 * @param str_2 the second string (must have the same length as str_1)
 * @return the hamming distance of str_1 and str_2
 */
template <typename pos_t>
static pos_t hamming_dist(const std::string_view str_1, const std::string_view str_2)
{
    assert(str_1.size() == str_2.size());
    pos_t dist = 0;

    for (pos_t i = 0; i < str_1.size(); i++) {
        dist += str_1[i] != str_2[i];
    }

    return dist;
}

/**
 * @brief computes the hamming distance of two equally long strings, stopping early once it exceeds k
 * @tparam pos_t unsigned integer type for the distance
 * @param str_1 the first string
 * @param str_2 the second string (must have the same length as str_1)
 * @param k the distance threshold
 * @return the hamming distance of str_1 and str_2 (or a value > k if it exceeds k)
 */
template <typename pos_t>
static pos_t hamming_dist_bounded(const std::string_view str_1, const std::string_view str_2, pos_t k)
{
    assert(str_1.size() == str_2.size());
    pos_t dist = 0;

    for (pos_t i = 0; i < str_1.size() && dist <= k; i++) {
        dist += str_1[i] != str_2[i];
    }

    return dist;
}

/**
 * @brief naively locates all occurrences of P in T with at most k mismatches (w.r.t. hamming distance)
 * @tparam pos_t unsigned integer position type
 * @param T the text
 * @param P the pattern
 * @param k the maximum number of mismatches
 * @return the occurrences of P in T with at most k mismatches
 */
template <typename pos_t>
static std::vector<aprx_occ_t<pos_t>> locate_hamming_dist(const std::string& T, const std::string& P, pos_t k)
{
    pos_t m = P.length();
    pos_t n = T.length();

    std::vector<aprx_occ_t<pos_t>> Occ;

    for (pos_t i = 0; i <= n - m; i++) {
        pos_t dist = hamming_dist_bounded<pos_t>(std::string_view(T.c_str() + i, m), P, k);
        if (dist <= k) Occ.emplace_back(aprx_occ_t<pos_t>{.pos = i, .len = m, .err = dist});
    }

    return Occ;
}

/**
 * @brief computes the edit (Levenshtein) distance of two strings using a full dynamic-programming matrix
 * @tparam pos_t unsigned integer type for the distance
 * @param T the first string
 * @param P the second string
 * @return the edit distance of T and P
 */
template <typename pos_t>
static pos_t edit_dist(const std::string_view T, const std::string_view P)
{
    pos_t n = T.size();
    pos_t m = P.size();

    std::vector<std::vector<pos_t>> dist(n + 1, std::vector<pos_t>(m + 1));

    for (pos_t i = 0; i <= n; i++) dist[i][0] = i;
    for (pos_t j = 0; j <= m; j++) dist[0][j] = j;

    for (pos_t i = 1; i <= n; i++) {
        char c1 = T[i - 1];

        for (pos_t j = 1; j <= m; j++) {
            char c2 = P[j - 1];
            pos_t cost = c1 != c2;

            dist[i][j] = std::min({
                dist[i - 1][j] + 1,        // delete
                dist[i][j - 1] + 1,        // insert
                dist[i - 1][j - 1] + cost  // match / mismatch
            });
        }
    }

    return dist[n][m];
}

/**
 * @brief computes the edit distance of two strings using a banded dynamic program, capped at k + 1
 * @tparam pos_t unsigned integer type for the distance
 * @param T the first string
 * @param P the second string
 * @param k the distance threshold
 * @return the edit distance of T and P, or k + 1 if it exceeds k
 */
template <typename pos_t>
static int64_t edit_dist_bounded(const std::string_view T, const std::string_view P, int64_t k)
{
    if (T.size() < P.size()) return edit_dist_bounded<pos_t>(P, T, k);

    int64_t n = T.size();
    int64_t m = P.size();
    const pos_t infty = k + 1;
    if (n - m > k) return infty;
    std::vector<pos_t> prev(m + 1), curr(m + 1);

    for (int64_t j = 0; j <= m; j++) {
        prev[j] = (j <= k) ? j : infty;
    }

    for (int64_t i = 1; i <= n; i++) {
        int64_t x = std::max<int64_t>(1, i - k);
        int64_t y = std::min<int64_t>(m, i + k);
        bool abort = true;

        curr[x - 1] = infty;
        if (i <= k) [[unlikely]] curr[0] = i;

        for (int64_t j = x; j <= y; j++) {
            curr[j] = std::min({
                infty,                               // limit to k + 1
                prev[j] + 1,                         // deletion
                curr[j - 1] + 1,                     // insertion
                prev[j - 1] + (T[i - 1] != P[j - 1]) // match / mismatch
            });

            if (curr[j] <= k) abort = false;
        }
        
        if (abort) [[unlikely]] return infty;
        if (y < m) [[likely]] curr[y + 1] = infty;
        prev.swap(curr);
    }

    return prev[m];
}

/**
 * @brief naively locates all occurrences of P in T with at most k errors (w.r.t. edit distance), keeping for
 *        each starting position the shortest minimum-error occurrence
 * @tparam pos_t unsigned integer position type
 * @param T the text
 * @param P the pattern
 * @param k the maximum number of errors
 * @return the occurrences of P in T with at most k errors
 */
template <typename pos_t>
static std::vector<aprx_occ_t<pos_t>> locate_edit_dist_naive(const std::string& T, const std::string& P, pos_t k)
{
    pos_t n = T.size();
    pos_t m = P.size();

    pos_t l_min = std::max<int64_t>(0, int64_t{m} - int64_t{k});
    pos_t l_max = m + k;

    std::vector<aprx_occ_t<pos_t>> Occ;

    for (pos_t i = 0; i < n; i++) {
        aprx_occ_t<pos_t> occ{.pos = i, .len = 0, .err = k + 1};

        for (pos_t l = l_min; l <= l_max; l++) {
            if (i + l > n) [[unlikely]] break;
            pos_t err = edit_dist_bounded<pos_t>(std::string_view(T.c_str() + i, l), P, k);
            if (err < occ.err || (err == occ.err && l < occ.len)) {occ.err = err; occ.len = l;}
        }

        if (occ.err <= k) Occ.emplace_back(occ);
    }

    return Occ;
}

/**
 * @brief locates all occurrences of P in T with at most k errors (w.r.t. edit distance) using Myers' bit-parallel
 *        banded algorithm (a band of width 2k+1 encoded per word), keeping for each starting position the shortest
 *        minimum-error occurrence
 * @tparam pos_t unsigned integer position type
 * @tparam word_t word type holding the band (its width must be at least 2k+1)
 * @param T the text
 * @param P the pattern
 * @param k the maximum number of errors
 * @param Peq for each byte, the bitmask of the positions in P at which it occurs
 * @return the occurrences of P in T with at most k errors
 */
template <typename pos_t, typename word_t>
static std::vector<aprx_occ_t<pos_t>> locate_edit_dist_bp(
    const std::string& T, const std::string& P, int64_t k,
    const std::array<__uint128_t, 256>& Peq
) {
    int64_t n = T.size();
    int64_t m = P.size();
    int64_t l_min = std::max<int64_t>(0, m - k);
    int64_t l_max = m + k;
    int64_t l_start = std::max<int64_t>(1, l_min);
    int64_t W = 2 * k + 1;
    const int64_t infty = k + 1;
    const word_t mask_W = W >= int64_t(8 * sizeof(word_t)) ?
        ~word_t(0) : ((word_t(1) << W) - 1);

    std::vector<aprx_occ_t<pos_t>> Occ;
    std::vector<int64_t> prev(W + 1), curr(W + 1);

    for (int64_t i = 0; i < n; i++) {
        aprx_occ_t<pos_t> occ{.pos = pos_t(i), .len = 0, .err = pos_t(k + 1)};
        int64_t j_max = std::min(l_max, n - i);
        int64_t j_0 = std::min(k, j_max);

        auto consider = [&](int64_t j, int64_t err) {
            if (j >= l_start && (err < int64_t(occ.err) ||
               (err == int64_t(occ.err) && j < int64_t(occ.len)))
            ) {
                occ.err = err;
                occ.len = j;
            }
        };

        for (int64_t b = 0; b <= W; b++) {
            int64_t r = b - k;
            prev[b] = (r >= 0 && r <= m) ? std::min(r, infty) : infty;
        }

        for (int64_t j = 1; j <= j_0; j++) {
            char c = T[i + j - 1];

            for (int64_t b = 0; b < W; b++) {
                int64_t r = j - k + b;

                if (r < 0 || r > m) {
                    curr[b] = infty;
                } else if (r == 0) {
                    curr[b] = std::min(j, infty);
                } else {
                    curr[b] = std::min({infty,
                        prev[b] + (P[r - 1] != c),           // match / mismatch
                        prev[b + 1] + 1,                     // deletion
                        (b > 0 ? curr[b - 1] : infty) + 1}); // insertion
                }
            }

            curr[W] = infty;
            int64_t b_m = m - (j - k);
            if (b_m >= 0 && b_m < W) consider(j, curr[b_m]);
            prev.swap(curr);
        }

        int64_t base = prev[0];
        word_t vp = 0, vn = 0;

        for (int64_t b = 1; b < W; b++) {
            int64_t y = prev[b] - prev[b - 1];
            if (y == 1) vp |= word_t(1) << b;
            else if (y == -1) vn |= word_t(1) << b;
        }

        for (int64_t j = j_0 + 1; j <= j_max; j++) {
            uint8_t c = char_to_uchar(T[i + j - 1]);
            int64_t sh = j - k - 1;
            word_t eq = word_t(Peq[c] >> sh) & mask_W;
            word_t z = eq | (vn >> 1);
            word_t prop = vp;

            for (int64_t s = 1; s < W; s <<= 1) {
                z |= prop & (z << s);
                prop &= (prop << s);
                z &= mask_W;
                prop &= mask_W;
            }

            word_t hp = ~z & mask_W;
            word_t hp_s = hp << 1;
            word_t d_pos = hp & ~hp_s & mask_W;
            word_t d_neg = ~hp & hp_s & mask_W;
            word_t d_zero = (~(d_pos | d_neg)) & mask_W;
            word_t v_zero = (~vp & ~vn) & mask_W;
            word_t vp_n = ((vp & d_zero) | (v_zero & d_pos)) & mask_W;
            word_t vn_n = ((vn & d_zero) | (v_zero & d_neg)) & mask_W;
            vp = vp_n & ~word_t(1);
            vn = vn_n & ~word_t(1);
            base += hp & 1;
            int64_t b_m = m - (j - k);

            if (b_m >= 0 && b_m < W) {
                word_t lo = b_m >= 1 ? (word_t(2) << b_m) - 2 : 0;
                int64_t err = base + popcount<word_t>(vp & lo) - popcount<word_t>(vn & lo);
                consider(j, std::min(err, infty));
            }
        }

        if (occ.err <= k) Occ.emplace_back(occ);
    }

    return Occ;
}

/**
 * @brief locates all occurrences of P in T with at most k errors (w.r.t. edit distance), choosing the bit-parallel
 *        band word type based on k (or falling back to the naive algorithm when the band does not fit)
 * @tparam pos_t unsigned integer position type
 * @param T the text
 * @param P the pattern
 * @param k the maximum number of errors
 * @return the occurrences of P in T with at most k errors
 */
template <typename pos_t>
static std::vector<aprx_occ_t<pos_t>> locate_edit_dist(const std::string& T, const std::string& P, pos_t k)
{
    int64_t m = P.size();
    if (int64_t(k) > 63 || m > 128) return locate_edit_dist_naive<pos_t>(T, P, k);

    std::array<__uint128_t, 256> Peq{};
    for (int64_t p = 0; p < m; p++) Peq[uint8_t(P[p])] |= (__uint128_t)(1) << p;

    if (k <= 15) return locate_edit_dist_bp<pos_t, uint32_t>(T, P, k, Peq);
    if (k <= 31) return locate_edit_dist_bp<pos_t, uint64_t>(T, P, k, Peq);
    return locate_edit_dist_bp<pos_t, __uint128_t>(T, P, k, Peq);
}

/**
 * @brief removes redundant approximate occurrences in place: an occurrence is dropped when a nearby occurrence
 *        (within 4k+3 positions) on each side has an equal or smaller error count
 * @tparam pos_t unsigned integer position type
 * @param Occ the occurrences (sorted by position); filtered in place
 * @param k the maximum number of errors
 */
template <typename pos_t>
static void filter_aprx_occurrences(std::vector<aprx_occ_t<pos_t>>& Occ, pos_t k)
{
    if (Occ.size() <= 2) return;
    pos_t write_idx = 0;

    for (pos_t read_idx = 0; read_idx < Occ.size(); read_idx++) {
        const auto occ = Occ[read_idx]; 

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

        Occ[write_idx++] = occ;
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
    int64_t max_dist = 2 * k + 1;
    uint64_t left = 0, right = 0;

    for (const auto& o : occ_all) {
        while (right < occ_filtered.size() && int64_t{occ_filtered[right].pos} <= int64_t{o.pos} + max_dist) right++;
        while (left < right && int64_t{occ_filtered[left].pos} < int64_t{o.pos} - max_dist) left++;
        if (left == right) return false;

        bool valid = false;
        for (uint64_t i = left; i < right; ++i) {
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
 * @param T the text
 * @param P the pattern
 * @param k the maximum number of errors
 * @return the occurrences of P in T with at most k errors
 */
template <typename pos_t, distance_metric_t dist_metr>
static std::vector<aprx_occ_t<pos_t>> locate(const std::string& T, const std::string& P, pos_t k)
{
    if constexpr (dist_metr == HAMMING_DISTANCE) return locate_hamming_dist<pos_t>(T, P, k);
    if constexpr (dist_metr == EDIT_DISTANCE) return locate_edit_dist<pos_t>(T, P, k);
}