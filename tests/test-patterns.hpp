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

// random patterns for the tests and their occurrences

#ifndef MOVE_R_TESTS_TEST_PATTERNS_HPP
#define MOVE_R_TESTS_TEST_PATTERNS_HPP

#include <algorithm>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

// apply_random_edits draws max_edits from [0, |pattern| / pattern_symbols_per_max_edit]
inline constexpr uint64_t pattern_symbols_per_max_edit = 2;

/**
 * @brief applies k random edits (substitution, insertion or deletion of one symbol) to a pattern, with k drawn from
 *        [0, max_edits] and max_edits drawn from [0, |pattern| / pattern_symbols_per_max_edit]; the inserted and
 *        substituted symbols are drawn from the positions of input
 * @tparam inp_t container of the pattern and of the input (e.g. std::string or std::vector<int32_t>)
 * @tparam gen_t random number generator type
 * @param pattern the pattern, edited in place (a non-empty pattern stays non-empty)
 * @param input the input the symbols are drawn from
 * @param gen random number generator
 */
template <typename inp_t, typename gen_t>
inline void apply_random_edits(inp_t& pattern, const inp_t& input, gen_t& gen)
{
    uint64_t input_size = input.size();
    uint64_t max_edits = std::uniform_int_distribution<uint64_t>(0, pattern.size() / pattern_symbols_per_max_edit)(gen);
    uint64_t num_edits = std::uniform_int_distribution<uint64_t>(0, max_edits)(gen);

    for (uint64_t k = 0; k < num_edits; k++) {
        uint64_t pos = std::uniform_int_distribution<uint64_t>(0, pattern.size() - 1)(gen);
        typename inp_t::value_type sym = input[std::uniform_int_distribution<uint64_t>(0, input_size - 1)(gen)];
        // 0: deletion, 1: insertion, 2: substitution (also if the deletion or insertion is not possible)
        uint64_t edit_type = gen() % 3;

        if (edit_type == 0 && pattern.size() > 1)               pattern.erase(pattern.begin() + pos);
        else if (edit_type == 1 && pattern.size() < input_size) pattern.insert(pattern.begin() + pos, sym);
        else                                                    pattern[pos] = sym;
    }
}

/**
 * @brief returns the starting positions of the occurrences of a pattern in an input, by comparing it at every position
 * @tparam inp_t container of the input and of the pattern (e.g. std::string or std::vector<int32_t>)
 * @param input the input
 * @param pattern the pattern
 * @return the starting positions of the occurrences of pattern in input, in increasing order
 */
template <typename inp_t>
inline std::vector<uint64_t> locate_naive(const inp_t& input, const inp_t& pattern)
{
    std::vector<uint64_t> occurrences;
    uint64_t input_size = input.size();
    uint64_t pattern_length = pattern.size();

    for (uint64_t i = 0; i + pattern_length <= input_size; i++) {
        bool match = true;

        for (uint64_t j = 0; j < pattern_length; j++) {
            if (input[i + j] != pattern[j]) {
                match = false;
                break;
            }
        }

        if (match) occurrences.emplace_back(i);
    }

    return occurrences;
}

/**
 * @brief draws a random substring of an input as a pattern, applies random edits to it (if mutate is true) and returns
 *        it together with the starting positions of its occurrences in the input
 * @tparam inp_t container of the input and of the pattern (e.g. std::string or std::vector<int32_t>)
 * @tparam gen_t random number generator type
 * @param input the input (|input| >= 1)
 * @param max_pattern_length the maximum length of the substring
 * @param mutate whether random edits are applied to the substring (see apply_random_edits)
 * @param gen random number generator
 * @return the pattern and the starting positions of its occurrences in input, in increasing order
 */
template <typename inp_t, typename gen_t>
inline std::pair<inp_t, std::vector<uint64_t>> random_pattern_and_occurrences(
    const inp_t& input, uint64_t max_pattern_length, bool mutate, gen_t& gen)
{
    uint64_t input_size = input.size();
    std::uniform_int_distribution<uint64_t> pattern_pos_distrib(0, input_size - 1);
    std::uniform_int_distribution<uint64_t> pattern_length_distrib(
        1, std::min<uint64_t>(max_pattern_length, input_size));

    uint64_t pattern_pos = pattern_pos_distrib(gen);
    uint64_t pattern_length = std::min<uint64_t>(input_size - pattern_pos, pattern_length_distrib(gen));
    inp_t pattern(input.begin() + pattern_pos, input.begin() + pattern_pos + pattern_length);
    if (mutate) apply_random_edits(pattern, input, gen);
    std::vector<uint64_t> occurrences = locate_naive(input, pattern);
    return { std::move(pattern), std::move(occurrences) };
}

#endif // MOVE_R_TESTS_TEST_PATTERNS_HPP
