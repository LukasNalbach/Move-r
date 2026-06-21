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
static const search_scheme_t pigeon_hole_scheme(uint8_t k)
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

/**
 * @brief builds the suffix-filter search scheme for k errors (k + 1 parts, k + 1 searches)
 * @param k the maximum number of errors
 * @return the suffix-filter search scheme for k errors
 */
static const search_scheme_t suffix_filter_scheme(uint8_t k)
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

/**
 * @brief builds the 0/1-seeds search scheme for k errors (k + 2 parts, k + 1 searches)
 * @param k the maximum number of errors
 * @return the 0/1-seeds search scheme for k errors
 */
static const search_scheme_t zero_one_scheme(uint8_t k)
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

/**
 * @brief precomputes, for a single search of a search scheme, the per-step information needed to execute it
 *        (search directions, direction switches and the range of previously processed parts)
 */
class edit_dist_search
{
  protected:
    const search_scheme_t& scheme;
    const std::vector<search_step_t>& search_arr;
    uint8_t search_idx;
    std::vector<direction_t> dirs;
    std::vector<bool> is_dir_switch;
    std::vector<uint8_t> leftmost_prev_part;
    std::vector<uint8_t> rightmost_prev_part;

  public:
    /**
     * @brief precomputes the per-step information for the search with index search_idx in scheme
     * @param scheme the search scheme
     * @param search_idx the index of the search within the scheme
     */
    edit_dist_search(const search_scheme_t& scheme, uint8_t search_idx)
        : scheme(scheme), search_arr(scheme.S[search_idx]), search_idx(search_idx)
    {
        dirs.reserve(scheme.p);
        dirs.emplace_back((search_arr[1].part > search_arr[0].part) ? RIGHT : LEFT);

        for (uint8_t i = 1; i < scheme.p; i++) {
            dirs.emplace_back((search_arr[i].part > search_arr[i - 1].part) ? RIGHT : LEFT);
        }

        is_dir_switch.reserve(scheme.p);
        is_dir_switch.emplace_back(false);

        for (uint8_t i = 1; i < dirs.size(); i++) {
            is_dir_switch.emplace_back(dirs[i] != dirs[i - 1]);
        }

        leftmost_prev_part.reserve(scheme.p);
        rightmost_prev_part.reserve(scheme.p);
        leftmost_prev_part.emplace_back(search_arr[0].part);
        rightmost_prev_part.emplace_back(search_arr[0].part);        

        for (uint8_t i = 1; i < scheme.p; i++) {
            uint8_t cur_part = search_arr[i].part;

            if (cur_part < leftmost_prev_part[i - 1]) {
                leftmost_prev_part.emplace_back(cur_part);
                rightmost_prev_part.emplace_back(rightmost_prev_part[i - 1]);
            } else {
                leftmost_prev_part.emplace_back(leftmost_prev_part[i - 1]);
                rightmost_prev_part.emplace_back(cur_part);
            }
        }
    }
    
    /**
     * @brief returns the minimum number of errors allowed after the i-th search step
     * @param i index of the search step
     * @return the lower error bound of the i-th search step
     */
    uint8_t lower_bound(uint8_t i) const
    {
        return search_arr[i].k_min;
    }

    /**
     * @brief returns the maximum number of errors allowed after the i-th search step
     * @param i index of the search step
     * @return the upper error bound of the i-th search step
     */
    uint8_t upper_bound(uint8_t i) const
    {
        return search_arr[i].k_max;
    }

    /**
     * @brief returns the pattern part processed in the i-th search step
     * @param i index of the search step
     * @return the index of the pattern part processed in the i-th search step
     */
    uint8_t part(uint8_t i) const
    {
        return search_arr[i].part;
    }

    /**
     * @brief returns the leftmost pattern part processed before the idx-th search step
     * @param idx index of the search step
     * @return the index of the leftmost previously processed part
     */
    uint8_t leftmost_previous_part(uint8_t idx) const
    {
        return leftmost_prev_part[idx - 1];
    }

    /**
     * @brief returns the rightmost pattern part processed before the idx-th search step
     * @param idx index of the search step
     * @return the index of the rightmost previously processed part
     */
    uint8_t rightmost_previous_part(uint8_t idx) const
    {
        return rightmost_prev_part[idx - 1];
    }

    /**
     * @brief returns the direction in which the i-th search step extends the pattern
     * @param i index of the search step
     * @return the search direction of the i-th search step
     */
    direction_t part_dir(uint8_t i) const
    {
        return dirs[i];
    }

    /**
     * @brief returns whether the i-th search step switches the search direction
     * @param i index of the search step
     * @return whether the i-th search step switches the search direction
     */
    bool does_part_switch_dir(uint8_t i) const
    {
        return is_dir_switch[i];
    }

    /**
     * @brief returns the number of pattern parts
     * @return the number of pattern parts
     */
    uint8_t num_parts() const
    {
        return scheme.p;
    }

    /**
     * @brief returns the index of this search within its search scheme
     * @return the index of this search within its search scheme
     */
    uint8_t search_index() const
    {
        return search_idx;
    }

    /**
     * @brief returns whether the part processed in the i-th search step is an outer (edge) part of the pattern
     * @param i index of the search step
     * @return whether the i-th part is the first or the last pattern part
     */
    bool is_edge(uint8_t i) const
    {
        return search_arr[i].part == 0 || search_arr[i].part == scheme.p - 1;
    }
};

#include "min_u_schemes.tpp"