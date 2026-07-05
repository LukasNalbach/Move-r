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
#include <vector>
#include <string>

#include <misc/files.hpp>

/**
 * @brief prints an error message about a malformed search scheme file and exits
 */
static void print_search_scheme_error()
{
    std::cout << "Error: malformed header in search scheme file" << std::endl;
    exit(0);
}

/**
 * @brief a single step of a search: process pattern part `part` while allowing between k_min and k_max errors
 */
struct search_step_t {
    uint8_t part;
    uint8_t k_min;
    uint8_t k_max;
};

/** @brief a search, that is the ordered sequence of search steps processing the pattern parts */
using search_t = std::vector<search_step_t>;

/**
 * @brief a search scheme: a set S of searches over a pattern split into p parts, allowing up to k errors
 */
struct search_scheme_t {
    uint8_t k = 0;
    uint8_t p = 1;
    std::vector<search_t> S = {{{0, 0, 0}}};
};

/**
 * @brief builds the pigeonhole search scheme for k errors (k + 1 parts, k + 1 searches)
 * @param k the maximum number of errors
 * @return the pigeonhole search scheme for k errors
 */
inline const search_scheme_t pigeon_hole_scheme(uint8_t k)
{
    uint8_t p = k + 1;
    std::vector<search_t> S;
    S.reserve(k + 1);

    for (int16_t i = 0; i < k + 1; i++) {
        search_t s;
        s.reserve(p);
        s.emplace_back(search_step_t{.part = (uint8_t)i, .k_min = 0, .k_max = 0});

        for (int16_t j = i + 1; j < p; j++) {
            s.emplace_back(search_step_t{.part = (uint8_t)j, .k_min = 0, .k_max = k});
        }

        for (int16_t j = i - 1; j >= 0; j--) {
            s.emplace_back(search_step_t{.part = (uint8_t)j, .k_min = 0, .k_max = k});
        }

        S.emplace_back(std::move(s));
    }

    return search_scheme_t {
        .k = k,
        .p = p,
        .S = std::move(S)
    };
}

/**
 * @brief builds the suffix-filter search scheme for k errors (k + 1 parts, k + 1 searches)
 * @param k the maximum number of errors
 * @return the suffix-filter search scheme for k errors
 */
inline const search_scheme_t suffix_filter_scheme(uint8_t k)
{
    uint8_t p = k + 1;
    std::vector<search_t> S;
    S.reserve(k + 1);

    for (int16_t i = 0; i < k + 1; i++) {
        search_t s;
        s.reserve(p);

        for (int16_t j = i; j < p; j++) {
            s.emplace_back(search_step_t{.part = (uint8_t)j, .k_min = 0, .k_max = (uint8_t)(j - i)});
        }

        for (int16_t j = i - 1; j >= 0; j--) {
            s.emplace_back(search_step_t{.part = (uint8_t)j, .k_min = 0, .k_max = k});
        }

        S.emplace_back(std::move(s));
    }

    return search_scheme_t {
        .k = k,
        .p = p,
        .S = std::move(S)
    };
}

/**
 * @brief builds the 0/1-seeds search scheme for k errors (k + 2 parts, k + 1 searches)
 * @param k the maximum number of errors
 * @return the 0/1-seeds search scheme for k errors
 */
inline const search_scheme_t zero_one_scheme(uint8_t k)
{
    uint8_t p = k + 2;
    std::vector<search_t> S;
    S.reserve(k + 1);

    for (int16_t i = 0; i < k + 1; i++) {
        search_t s;
        s.reserve(p);
        s.emplace_back(search_step_t{.part = (uint8_t)i, .k_min = 0, .k_max = 0});
        s.emplace_back(search_step_t{.part = (uint8_t)(i + 1), .k_min = 0, .k_max = (uint8_t)(i == k ? 0 : 1)});

        for (int16_t j = i + 2; j < p; j++) {
            s.emplace_back(search_step_t{.part = (uint8_t)j, .k_min = 0, .k_max = k});
        }

        for (int16_t j = i - 1; j >= 0; j--) {
            s.emplace_back(search_step_t{.part = (uint8_t)j, .k_min = 0, .k_max = k});
        }

        S.emplace_back(std::move(s));
    }

    return search_scheme_t {
        .k = k,
        .p = p,
        .S = std::move(S)
    };
}

/**
 * @brief parses a comma-separated list of part indices from a single bracket of a search scheme file
 * @param content the content between a pair of braces
 * @param p the number of parts (each parsed value must be in [0, p))
 * @return the parsed part indices
 */
static std::vector<uint8_t> parse_bracket(std::string content, int64_t p)
{
    std::replace(content.begin(), content.end(), ',', ' ');
    
    std::stringstream ss(content);
    std::vector<uint8_t> numbers;
    int64_t num;
    
    while (ss >> num) {
        if (!(0 <= num && num < p)) {
            print_search_scheme_error();
        }

        numbers.push_back(num);
    }

    return numbers;
}

/**
 * @brief parses a search scheme from its textual representation
 * @param str the textual representation of the search scheme (a header line "p=.. k=.." followed by one line per search)
 * @return the parsed search scheme
 */
inline search_scheme_t parse_search_scheme(const std::string& str)
{
    std::stringstream ss_input(str);
    std::string line;
    std::getline(ss_input, line);

    int64_t p = value_from_key(line, "p=");
    int64_t k = value_from_key(line, "k=");

    if (p == -1 || k == -1) {
        print_search_scheme_error();
    }

    std::vector<search_t> S;

    // each non-empty line describes one search as three consecutive braced arrays: {pi} {L} {U}, each holding p
    // comma-separated values (the part order, the lower and the upper error bounds)
    while (std::getline(ss_input, line)) {
        if (line.find('{') == std::string::npos) continue; // skip blank lines

        std::vector<std::vector<uint8_t>> arrays; // the parsed pi, L and U arrays
        for (size_t pos = 0; arrays.size() < 3;) {
            size_t open = line.find('{', pos);
            size_t close = open == std::string::npos ? std::string::npos : line.find('}', open);
            if (close == std::string::npos) print_search_scheme_error();
            arrays.emplace_back(parse_bracket(line.substr(open + 1, close - open - 1), p));
            if ((int64_t) arrays.back().size() != p) print_search_scheme_error();
            pos = close + 1;
        }

        search_t s;
        s.reserve(p);
        for (int64_t i = 0; i < p; i++)
            s.emplace_back(search_step_t{.part = arrays[0][i], .k_min = arrays[1][i], .k_max = arrays[2][i]});
        S.emplace_back(std::move(s));
    }

    return search_scheme_t {
        .k = (uint8_t)k,
        .p = (uint8_t)p,
        .S = std::move(S)
    };
}

#include "min_u_schemes.tpp"