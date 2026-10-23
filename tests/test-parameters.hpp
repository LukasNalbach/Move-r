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

#include <array>
#include <cstdint>
#include <limits>

#include <move_data_structure/move_data_structure.hpp>
#include <move_r/construction/text_alphabet.hpp>

// the number of byte values
inline constexpr uint64_t byte_alphabet_size = uint64_t { std::numeric_limits<uint8_t>::max() } + 1;

// the smallest internal symbol a byte is mapped to with _pfp and with the in-memory modes (see
// move_r::construction::min_mapped_symbol); the values below it are reserved (0 is the terminator)
inline constexpr uint64_t pfp_min_mapped_symbol = move_r_min_mapped_symbol_pfp;
inline constexpr uint64_t in_memory_min_mapped_symbol = move_r_min_mapped_symbol;

// the smallest number of distinct symbols of a random input whose alphabet size is drawn
inline constexpr uint64_t min_random_alphabet_size = 2;

// the largest number of distinct symbols of a byte input that every construction mode accepts
inline constexpr uint64_t max_byte_alphabet_size = byte_alphabet_size - pfp_min_mapped_symbol;

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
