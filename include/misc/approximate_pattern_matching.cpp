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
#include <string>
#include <vector>

template <typename pos_t>
struct aprx_occ_t {
    pos_t pos;
    pos_t len;
    pos_t err;

    bool operator<(const aprx_occ_t& other) const
    {
        return pos != other.pos ? pos < other.pos :
              (err != other.err ? err < other.err :
                                  len < other.len);
    }

    bool operator==(const aprx_occ_t& other) const = default;

    friend std::ostream& operator<<(std::ostream& os, const aprx_occ_t& occ) {
        return os << "(" << occ.pos << ", " << occ.len << ", " << occ.err << ")";
    }
};

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

template <typename pos_t>
static std::vector<aprx_occ_t<pos_t>> locate_edit_dist(const std::string& T, const std::string& P, pos_t k)
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

enum distance_metric_t : int8_t {
    NO_METRIC = -1,
    HAMMING_DISTANCE = 0,
    EDIT_DISTANCE = 1
};

template <typename pos_t, distance_metric_t dist_metr>
static std::vector<aprx_occ_t<pos_t>> locate(const std::string& T, const std::string& P, pos_t k)
{
    if constexpr (dist_metr == HAMMING_DISTANCE) return locate_hamming_dist<pos_t>(T, P, k);
    if constexpr (dist_metr == EDIT_DISTANCE) return locate_edit_dist<pos_t>(T, P, k);
}