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

#include "files.hpp"

static void print_search_scheme_error()
{
    std::cout << "Error: malformed header in search scheme file" << std::endl;
    exit(0);
}

struct search_step_t {
    uint8_t part;
    uint8_t k_min;
    uint8_t k_max;
};

using search_t = std::vector<search_step_t>;

struct search_scheme_t {
    uint8_t k = 0;
    uint8_t p = 1;
    std::vector<search_t> S = {{{0, 0, 0}}};
};

static search_scheme_t pigeon_hole_scheme(uint8_t k)
{
    uint8_t p = k + 1;
    std::vector<search_t> S;
    S.reserve(k + 1);

    for (int16_t i = 0; i < k + 1; i++) {
        search_t s;
        s.reserve(p);
        s.emplace_back(search_step_t{.part = i, .k_min = 0, .k_max = 0});

        for (int16_t j = i + 1; j < p; j++) {
            s.emplace_back(search_step_t{.part = j, .k_min = 0, .k_max = k});
        }

        for (int16_t j = i - 1; j >= 0; j--) {
            s.emplace_back(search_step_t{.part = j, .k_min = 0, .k_max = k});
        }

        S.emplace_back(std::move(s));
    }

    return search_scheme_t {
        .k = k,
        .p = p,
        .S = std::move(S)
    };
}

static search_scheme_t suffix_filter_scheme(uint8_t k)
{
    uint8_t p = k + 1;
    std::vector<search_t> S;
    S.reserve(k + 1);

    for (int16_t i = 0; i < k + 1; i++) {
        search_t s;
        s.reserve(p);

        for (int16_t j = i; j < p; j++) {
            s.emplace_back(search_step_t{.part = j, .k_min = 0, .k_max = j - i});
        }

        for (int16_t j = i - 1; j >= 0; j--) {
            s.emplace_back(search_step_t{.part = j, .k_min = 0, .k_max = k});
        }

        S.emplace_back(std::move(s));
    }

    return search_scheme_t {
        .k = k,
        .p = p,
        .S = std::move(S)
    };
}

static search_scheme_t zero_one_scheme(uint8_t k)
{
    uint8_t p = k + 2;
    std::vector<search_t> S;
    S.reserve(k + 1);

    for (int16_t i = 0; i < k + 1; i++) {
        search_t s;
        s.reserve(p);
        s.emplace_back(search_step_t{.part = i, .k_min = 0, .k_max = 0});
        s.emplace_back(search_step_t{.part = i + 1, .k_min = 0, .k_max = (i == k ? 0 : 1)});

        for (int16_t j = i + 2; j < p; j++) {
            s.emplace_back(search_step_t{.part = j, .k_min = 0, .k_max = k});
        }

        for (int16_t j = i - 1; j >= 0; j--) {
            s.emplace_back(search_step_t{.part = j, .k_min = 0, .k_max = k});
        }

        S.emplace_back(std::move(s));
    }

    return search_scheme_t {
        .k = k,
        .p = p,
        .S = std::move(S)
    };
}

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

static search_scheme_t parse_search_scheme(const std::string& str)
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

    while (std::getline(ss_input, line)) {
        if (line.empty()) {
            print_search_scheme_error();
        }

        std::vector<std::vector<uint8_t>> bracket_contents;
        uint64_t pos = 0;

        while ((pos = line.find('{', pos)) != std::string::npos) {
            uint64_t end = line.find('}', pos);
            if (end == std::string::npos) break;

            std::string content = line.substr(pos + 1, end - pos - 1);
            bracket_contents.emplace_back(parse_bracket(content, p));
            
            pos = end + 1;
        }

        if (bracket_contents.size() != 3 ||
            bracket_contents[0].size() != p ||
            bracket_contents[1].size() != p ||
            bracket_contents[2].size() != p
        ) {
            print_search_scheme_error();
        }

        search_t s;
        s.reserve(p);

        for (int64_t i = 0; i < p; i++) {
            s.emplace_back(search_step_t{
                .part =  bracket_contents[0][i],
                .k_min = bracket_contents[1][i],
                .k_max = bracket_contents[2][i]
            });
        }

        S.emplace_back(s);
    }

    return search_scheme_t {
        .k = k,
        .p = p,
        .S = std::move(S)
    };
}