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
#include <random>
#include <fstream>
#include <ips4o.hpp>

#include "utils.hpp"

static std::string random_alphanumeric_string(uint64_t length)
{
    static std::string possible_chars = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint8_t> char_idx_distrib(0, possible_chars.size() - 1);
    
    std::string str_rand;
    str_rand.reserve(length);

    for (uint64_t i = 0; i < length; i++) {
        str_rand.push_back(possible_chars[char_idx_distrib(gen)]);
    }

    return str_rand;
}

static uint64_t malloc_count_peak_memory_usage(std::ifstream& log_file)
{
    std::string log_file_content;
    log_file.seekg(0, std::ios::end);
    no_init_resize(log_file_content, log_file.tellg());
    log_file.seekg(0, std::ios::beg);
    log_file.read((char*) &log_file_content[0], log_file_content.size());
    int32_t pos = 0;
    uint64_t cur_peak = 0;
    std::string str_cur_peak;

    while ((pos = log_file_content.find(", peak", pos)) != -1) {
        while (!('0' <= log_file_content[pos] && log_file_content[pos] <= '9')) {
            pos++;
        }

        while (('0' <= log_file_content[pos] && log_file_content[pos] <= '9') || log_file_content[pos] == '.') {
            if (log_file_content[pos] != '.') {
                str_cur_peak.push_back(log_file_content[pos]);
            }

            pos++;
        }

        cur_peak = std::max(cur_peak, (uint64_t)stol(str_cur_peak));
        str_cur_peak.clear();
    }

    return cur_peak;
}

template <typename inp_t>
static inp_t random_repetitive_input(
    uint64_t min_size, uint64_t max_size,
    typename inp_t::value_type min_sym = std::numeric_limits<typename inp_t::value_type>::min(),
    typename inp_t::value_type max_sym = std::numeric_limits<typename inp_t::value_type>::max()
) {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);
    std::uniform_int_distribution<uint64_t> input_size_distrib(min_size, max_size);
    std::uniform_int_distribution<typename inp_t::value_type> sym_distrib(min_sym, max_sym);

    uint64_t target_input_size = input_size_distrib(mt);
    enum construction_operation { new_symbol = 0, repetition = 1, run = 2 };
    double repetition_repetitiveness = prob_distrib(mt);
    double run_repetitiveness = prob_distrib(mt);

    std::uniform_int_distribution<uint64_t> repetition_length_distrib(
        1, (repetition_repetitiveness * target_input_size) / 100);

    std::uniform_int_distribution<uint64_t> run_length_distrib(
        1, (run_repetitiveness * target_input_size) / 200);

    std::discrete_distribution<uint8_t> next_operation_distrib({
        2 - (repetition_repetitiveness + run_repetitiveness),
        repetition_repetitiveness,
        run_repetitiveness });

    inp_t input;
    no_init_resize(input, target_input_size);
    input.clear();
    input.push_back(sym_distrib(mt));

    while (input.size() < target_input_size) {
        switch (next_operation_distrib(mt)) {
            case new_symbol: {
                input.push_back(sym_distrib(mt));
                break;
            }
            case repetition: {
                uint64_t repetition_length = std::min<uint64_t>(
                    target_input_size - input.size(), repetition_length_distrib(mt));
                uint64_t repstition_source = std::uniform_int_distribution<uint64_t>(0, input.size() - 1)(mt);

                for (uint64_t i = 0; i < repetition_length; i++)
                    input.push_back(input[repstition_source + i]);

                break;
            }
            case run: {
                uint64_t run_length = std::min<uint64_t>(
                    target_input_size - input.size(), run_length_distrib(mt));
                typename inp_t::value_type run_sym = sym_distrib(mt);

                for (uint64_t i = 0; i < run_length; i++)
                    input.push_back(run_sym);

                break;
            }
        }
    }

    return input;
}

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
static pos_t edit_dist(const std::string_view T, const std::string_view P) {
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
static int64_t edit_dist_bounded(const std::string_view T, const std::string_view P, int64_t k) {
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