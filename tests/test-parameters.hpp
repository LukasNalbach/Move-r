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

// test parameters shared by several test files

#ifndef MOVE_R_TESTS_TEST_PARAMETERS_HPP
#define MOVE_R_TESTS_TEST_PARAMETERS_HPP

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <string_view>
#include <vector>

#include <misc/memory.hpp>
#include <move_data_structure/move_data_structure.hpp>
#include "test-inputs.hpp"

// the number of byte values
inline constexpr uint64_t byte_alphabet_size = uint64_t { std::numeric_limits<uint8_t>::max() } + 1;

// the smallest number of distinct symbols of a random integer input whose alphabet size is drawn
inline constexpr uint64_t min_random_alphabet_size = 2;

// a small random byte alphabet has at most as many bytes as A, C, G, T and N (see random_byte_alphabet)
inline constexpr uint64_t max_small_byte_alphabet_size = std::string_view("ACGTN").size();

// the alphabet of a random byte input (see random_byte_alphabet); it may contain the byte 0 and all 256 byte values
struct random_input_alphabet {
    std::vector<uint8_t> bytes; // the distinct bytes of the alphabet in random order

    template <typename gen_t>
    explicit random_input_alphabet(gen_t& gen) : bytes(random_byte_alphabet(max_small_byte_alphabet_size, gen)) { }

    // returns a random repetitive input over the alphabet with a length in [min_size, max_size] (see
    // random_repetitive_input), in which every byte of the alphabet occurs if the input is long enough
    template <typename inp_t>
    inp_t input(uint64_t min_size, uint64_t max_size, uint64_t seed) const
    {
        std::vector<uint8_t> ranks = random_repetitive_input<std::vector<uint8_t>>(
            min_size, max_size, uint8_t { 0 }, uint8_t(bytes.size() - 1), seed);

        if (ranks.size() >= bytes.size()) {
            for (uint64_t c = 0; c < bytes.size(); c++) ranks[c * ranks.size() / bytes.size()] = uint8_t(c);
        }

        inp_t input;
        resize_for_overwrite(input, ranks.size());
        for (uint64_t i = 0; i < ranks.size(); i++) input[i] = typename inp_t::value_type(bytes[ranks[i]]);
        return input;
    }

    // returns a byte that is not in the alphabet (one that is, if the alphabet has all byte values)
    uint8_t absent() const
    {
        std::array<bool, byte_alphabet_size> in_alphabet {};
        for (uint8_t c : bytes) in_alphabet[c] = true;
        const uint64_t c = std::find(in_alphabet.begin(), in_alphabet.end(), false) - in_alphabet.begin();
        return c < byte_alphabet_size ? uint8_t(c) : bytes[0];
    }
};

// the largest balancing parameter a random index is built with (the largest int16_t)
inline constexpr uint64_t max_balancing_parameter = std::numeric_limits<int16_t>::max();
// a random balancing parameter is min_balancing_parameter plus a log-normal value with these parameters (capped at
// max_balancing_parameter)
inline constexpr double balancing_parameter_log_mean = 2.0;
inline constexpr double balancing_parameter_log_deviation = 3.0;

// the smallest input that is built with the prefix-free parsing (_pfp)
inline constexpr uint64_t min_pfp_input_size = 60;
// an input of at least min_pfp_input_size symbols is built with _pfp with probability 1 / pfp_inverse_probability
inline constexpr uint64_t pfp_inverse_probability = 8;

// the values of move_r_params::context_table_share a random index is built with (0: no context tables)
inline constexpr std::array<double, 3> context_table_shares { 0.0, 0.1, 0.9 };

#endif
