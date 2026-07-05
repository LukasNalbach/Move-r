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

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>
#include <random>
#include <fstream>

#include "utils.hpp"

/**
 * @brief generates a random string of alphanumeric characters
 * @param length number of characters in the generated string
 * @return a random string of the given length over the alphabet [0-9A-Za-z]
 */
inline std::string random_alphanumeric_string(uint64_t length)
{
    static std::string possible_chars = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> char_idx_distrib(0, possible_chars.size() - 1);
    
    std::string str_rand;
    str_rand.reserve(length);

    for (uint64_t i = 0; i < length; i++) {
        str_rand.push_back(possible_chars[char_idx_distrib(gen)]);
    }

    return str_rand;
}

/**
 * @brief extracts the peak memory usage from a malloc_count log file
 * @param leftg_file input stream of the malloc_count log file
 * @return the maximum peak memory usage reported in the log file
 */
inline uint64_t malloc_count_peak_memory_usage(std::ifstream& leftg_file)
{
    std::string leftg_file_content;
    leftg_file.seekg(0, std::ios::end);
    no_init_resize(leftg_file_content, leftg_file.tellg());
    leftg_file.seekg(0, std::ios::beg);
    leftg_file.read((char*) &leftg_file_content[0], leftg_file_content.size());
    int32_t pos = 0;
    uint64_t cur_peak = 0;
    std::string str_cur_peak;

    while ((pos = leftg_file_content.find(", peak", pos)) != -1) {
        while (!('0' <= leftg_file_content[pos] && leftg_file_content[pos] <= '9')) {
            pos++;
        }

        while (('0' <= leftg_file_content[pos] && leftg_file_content[pos] <= '9') || leftg_file_content[pos] == '.') {
            if (leftg_file_content[pos] != '.') {
                str_cur_peak.push_back(leftg_file_content[pos]);
            }

            pos++;
        }

        cur_peak = std::max(cur_peak, (uint64_t)stol(str_cur_peak));
        str_cur_peak.clear();
    }

    return cur_peak;
}

/**
 * @brief draws an integer in [min_size, max_size] log-uniformly, so inputs of every magnitude
 *        (including very small ones) are generated about equally often — unlike a plain uniform draw,
 *        which almost never picks a small size when the range spans several orders of magnitude
 * @tparam gen_t random number generator type
 * @param min_size minimum size (inclusive)
 * @param max_size maximum size (inclusive)
 * @param gen random number generator
 * @return a log-uniformly distributed size in [min_size, max_size]
 */
template <typename gen_t>
inline uint64_t random_log_uniform_size(uint64_t min_size, uint64_t max_size, gen_t& gen)
{
    std::uniform_real_distribution<double> log_distrib(
        std::log((double) std::max<uint64_t>(1, min_size)),
        std::log((double) std::max<uint64_t>(1, max_size)));
    return std::clamp<uint64_t>((uint64_t) std::llround(std::exp(log_distrib(gen))), min_size, max_size);
}

/**
 * @brief generates a random, repetitive input sequence (composed of new symbols, repetitions and runs)
 * @tparam inp_t sequence type to generate (e.g. std::string or std::vector<int32_t>)
 * @param min_size minimum length of the generated sequence
 * @param max_size maximum length of the generated sequence
 * @param min_sym smallest symbol value that may occur in the sequence
 * @param max_sym largest symbol value that may occur in the sequence
 * @return a random sequence with a length in [min_size, max_size]
 */
template <typename inp_t>
static inp_t random_repetitive_input(
    uint64_t min_size, uint64_t max_size,
    typename inp_t::value_type min_sym = std::numeric_limits<typename inp_t::value_type>::min(),
    typename inp_t::value_type max_sym = std::numeric_limits<typename inp_t::value_type>::max()
) {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);
    using sym_t = typename inp_t::value_type;
    using sym_dist_t = std::conditional_t<sizeof(sym_t) == 1,
        std::conditional_t<std::is_signed_v<sym_t>, int, unsigned int>, sym_t>;
    std::uniform_int_distribution<sym_dist_t> sym_distrib(min_sym, max_sym);

    uint64_t target_input_size = random_log_uniform_size(min_size, max_size, mt);
    enum construction_operation { new_symbol = 0, repetition = 1, run = 2 };
    double repetition_repetitiveness = prob_distrib(mt);
    double run_repetitiveness = prob_distrib(mt);

    std::uniform_int_distribution<uint64_t> repetition_length_distrib(
        1, std::max<double>(1.0, (repetition_repetitiveness * target_input_size) / 100));

    std::uniform_int_distribution<uint64_t> run_length_distrib(
        1, std::max<double>(1.0, ((run_repetitiveness * target_input_size) / 200)));

    std::discrete_distribution<uint32_t> next_operation_distrib({
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
                input.push_back((sym_t)sym_distrib(mt));
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

/**
 * @brief applies a random number of edits (substitutions, insertions and deletions) to pattern, drawing the
 * replacement symbols from input; the number of edits is skewed towards few (at most half the pattern
 * length), so that the mutated pattern is usually no longer an exact substring of input, while exact and
 * lightly-modified patterns are still produced frequently
 * @tparam pos_t integer type used for positions and counts
 * @tparam inp_t sequence type of the pattern and input (e.g. std::string or std::vector<int32_t>)
 * @tparam gen_t random number generator type
 * @param pattern the pattern to mutate in place
 * @param input the sequence to draw the replacement symbols from
 * @param gen random number generator
 */
template <typename pos_t, typename inp_t, typename gen_t>
static void mutate_pattern(inp_t& pattern, const inp_t& input, gen_t& gen)
{
    pos_t input_size = input.size();
    pos_t max_mutations = std::uniform_int_distribution<pos_t>(0, pattern.size() / 2)(gen);
    pos_t num_mutations = std::uniform_int_distribution<pos_t>(0, max_mutations)(gen);

    for (pos_t m = 0; m < num_mutations; m++) {
        pos_t pos = std::uniform_int_distribution<pos_t>(0, pattern.size() - 1)(gen);
        typename inp_t::value_type sym = input[std::uniform_int_distribution<pos_t>(0, input_size - 1)(gen)];
        pos_t op = gen() % 3;

        if (op == 0 && pattern.size() > 1)               pattern.erase(pattern.begin() + pos);
        else if (op == 1 && pattern.size() < input_size) pattern.insert(pattern.begin() + pos, sym);
        else                                             pattern[pos] = sym;
    }
}

/**
 * @brief naively scans input for all exact occurrences of pattern
 * @tparam pos_t integer type used for positions
 * @tparam inp_t sequence type of the input and pattern (e.g. std::string or std::vector<int32_t>)
 * @param input the sequence to search in
 * @param pattern the pattern to search for
 * @return the ascending starting positions of all exact occurrences of pattern in input
 */
template <typename pos_t, typename inp_t>
static std::vector<pos_t> locate_naive(const inp_t& input, const inp_t& pattern)
{
    std::vector<pos_t> occurrences;
    pos_t input_size = input.size();
    pos_t pattern_length = pattern.size();

    for (pos_t i = 0; i + pattern_length <= input_size; i++) {
        bool match = true;

        for (pos_t j = 0; j < pattern_length; j++) {
            if (input[i + j] != pattern[j]) {
                match = false;
                break;
            }
        }

        if (match) occurrences.emplace_back(i);
    }

    return occurrences;
}