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
 * @file search_schemes.hpp
 * @brief search schemes: searches, steps, the built-in schemes and the parsing of scheme files
 */

#ifndef MOVE_R_MISC_SEARCH_SCHEMES_HPP
#define MOVE_R_MISC_SEARCH_SCHEMES_HPP

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>
#include <sstream>
#include <string>
#include <cstdlib>
#include <stdexcept>

#include <misc/pattern_reader.hpp>

// throws the error about a malformed search scheme file
[[noreturn]] inline void throw_search_scheme_error()
{
    throw std::runtime_error("malformed search scheme file");
}

/**
 * @brief one step of a search: matches the pattern part `part`; afterwards, the number of errors of the parts matched
 *        so far must be in [k_min, k_max]
 */
struct search_step_t {
    uint8_t part;  // [0..p-1] the index of the part
    uint8_t k_min; // the lower bound on the number of errors after this step
    uint8_t k_max; // the upper bound on the number of errors after this step
};

/** @brief a search: [0..p-1] the steps in the order in which the parts are matched */
using search_t = std::vector<search_step_t>;

/**
 * @brief a search scheme: a set of searches over a pattern split into p parts, which together find every occurrence
 *        with at most k errors
 */
struct search_scheme_t {
    uint8_t k = 0; // the maximum number of errors
    uint64_t p = 1; // the number of parts
    std::vector<search_t> S = {{{0, 0, 0}}}; // the searches (by default, the exact search of one part)
};

// returns whether a pattern with m symbols is long enough for the search scheme: the searches need more symbols than parts
inline bool fits_search_scheme(uint64_t m, const search_scheme_t& scheme) { return m > scheme.p; }

/**
 * @brief returns the error message for a pattern that is too short for the search scheme (see fits_search_scheme)
 * @param pattern how the pattern is named in the message, e.g. "the pattern read_1"
 * @param m the pattern length
 * @param scheme the search scheme
 */
inline std::string too_short_for_scheme(const std::string& pattern, uint64_t m, const search_scheme_t& scheme)
{
    return pattern + " has " + std::to_string(m) + " symbols; the search scheme with " + std::to_string(scheme.p)
        + " parts needs at least " + std::to_string(scheme.p + 1);
}

/**
 * @brief builds the pigeonhole search scheme for k errors (k + 1 parts, k + 1 searches)
 * @param k the maximum number of errors
 * @return the pigeonhole search scheme for k errors
 */
inline search_scheme_t pigeon_hole_scheme(uint8_t k)
{
    uint64_t p = k + 1;
    std::vector<search_t> S;
    S.reserve(k + 1);

    for (int16_t i = 0; i < k + 1; i++) {
        search_t s;
        s.reserve(p);
        s.emplace_back(search_step_t{.part = (uint8_t)i, .k_min = 0, .k_max = 0});

        for (int16_t j = i + 1; uint64_t(j) < p; j++) {
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
inline search_scheme_t suffix_filter_scheme(uint8_t k)
{
    uint64_t p = k + 1;
    std::vector<search_t> S;
    S.reserve(k + 1);

    for (int16_t i = 0; i < k + 1; i++) {
        search_t s;
        s.reserve(p);

        for (int16_t j = i; uint64_t(j) < p; j++) {
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
inline search_scheme_t zero_one_scheme(uint8_t k)
{
    uint64_t p = k + 2;
    std::vector<search_t> S;
    S.reserve(k + 1);

    for (int16_t i = 0; i < k + 1; i++) {
        search_t s;
        s.reserve(p);
        s.emplace_back(search_step_t{.part = (uint8_t)i, .k_min = 0, .k_max = 0});
        s.emplace_back(search_step_t{.part = (uint8_t)(i + 1), .k_min = 0, .k_max = (uint8_t)(i == k ? 0 : 1)});

        for (int16_t j = i + 2; uint64_t(j) < p; j++) {
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
 * @brief returns the mirror of a search scheme: part i of each search becomes part p - 1 - i (the same
 *        occurrences are found, with the parts in the reverse order)
 * @param scheme the search scheme
 * @return the mirrored search scheme
 */
inline search_scheme_t mirrored_scheme(const search_scheme_t& scheme)
{
    search_scheme_t mirror = scheme;
    for (search_t& S : mirror.S) for (search_step_t& step : S) step.part = uint8_t(scheme.p - 1 - step.part);
    return mirror;
}

/**
 * @brief parses a comma-separated list of numbers between a pair of braces of a search scheme file
 * @param content the content between the braces
 * @param max_value the largest allowed value (throws if a value is larger)
 * @return the parsed numbers
 */
inline std::vector<uint8_t> parse_bracket(std::string content, int64_t max_value)
{
    std::replace(content.begin(), content.end(), ',', ' ');

    std::stringstream ss(content);
    std::vector<uint8_t> numbers;
    int64_t num;

    while (ss >> num) {
        if (!(0 <= num && num <= max_value)) {
            throw_search_scheme_error();
        }

        numbers.push_back(num);
    }

    return numbers;
}

/**
 * @brief parses a search scheme from its textual representation
 *
 * A search scheme for k errors with p parts must be complete for the edit distance (this is not checked): for every
 * distribution of at most k errors over the p parts and the p - 1 boundaries of adjacent parts, some search accepts
 * the numbers of errors it counts per part, i.e., after each step of the search, the errors counted in the parts matched
 * so far sum to a value in [k_min, k_max] of the step. An error inside a part (a substitution, an insertion of a
 * symbol of the part, or a deletion between two symbols of the part) is counted in this part; a deletion of text symbols
 * between the last symbol of a part and the first symbol of the next part is counted in the one of the two parts that
 * the search matches later. The built-in schemes and the files in search_schemes/ are complete.
 * @param str the textual representation of the search scheme (a header line "p=.. k=.." followed by one line per search)
 * @return the parsed search scheme
 * @throws std::runtime_error if the search scheme is malformed, e.g., if a part of a search is not adjacent to the parts
 *         searched before it, if a lower bound of a step is larger than its upper bound, or if a lower or an upper bound
 *         of a search decreases from one step to the next
 */
inline search_scheme_t parse_search_scheme(const std::string& str)
{
    std::stringstream ss_input(str);
    std::string line;
    std::getline(ss_input, line);

    int64_t p = header_value(line, "p=");
    int64_t k = header_value(line, "k=");

    // a part index and an error count are stored in a byte
    constexpr int64_t max_byte = std::numeric_limits<uint8_t>::max();
    if (p < 1 || p > max_byte + 1 || k < 0 || k > max_byte) {
        throw_search_scheme_error();
    }

    std::vector<search_t> S;

    /* each line with a '{' describes one search as three braced arrays {pi} {L} {U} with p comma-separated values
       each: the order of the parts (values in [0, p-1]), the lower and the upper bounds on the errors (in [0, k]) */
    while (std::getline(ss_input, line)) {
        if (line.find('{') == std::string::npos) continue;

        std::vector<std::vector<uint8_t>> arrays; // [0..2] the parsed arrays pi, L and U

        for (uint64_t pos = 0; arrays.size() < 3;) {
            uint64_t open = line.find('{', pos);
            uint64_t close = open == std::string::npos ? std::string::npos : line.find('}', open);
            if (close == std::string::npos) throw_search_scheme_error();
            arrays.emplace_back(parse_bracket(line.substr(open + 1, close - open - 1), arrays.empty() ? p - 1 : k));
            if ((int64_t) arrays.back().size() != p) throw_search_scheme_error();
            pos = close + 1;
        }

        search_t s;
        s.reserve(p);
        int64_t lo = arrays[0][0], hi = arrays[0][0]; // the parts searched so far are lo, ..., hi

        for (int64_t i = 0; i < p; i++) {
            if (arrays[1][i] > arrays[2][i]) throw_search_scheme_error();

            /* each part after the first is adjacent to the parts searched so far (so pi is a permutation of 0, ..., p - 1),
               and the bounds of a step are at least those of the previous step */
            if (i > 0) {
                if (arrays[1][i] < arrays[1][i - 1] || arrays[2][i] < arrays[2][i - 1]) throw_search_scheme_error();
                if (arrays[0][i] == lo - 1) lo--;
                else if (arrays[0][i] == hi + 1) hi++;
                else throw_search_scheme_error();
            }

            s.emplace_back(search_step_t{.part = arrays[0][i], .k_min = arrays[1][i], .k_max = arrays[2][i]});
        }

        S.emplace_back(std::move(s));
    }

    return search_scheme_t {
        .k = (uint8_t)k,
        .p = uint64_t(p),
        .S = std::move(S)
    };
}

#include "min_u_schemes.hpp"

#endif // MOVE_R_MISC_SEARCH_SCHEMES_HPP
