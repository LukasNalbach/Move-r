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

#include <gtest/gtest.h>
#include <move_r/move_rb.hpp>

std::random_device rd;
std::mt19937 gen(rd());
uint16_t max_num_threads = omp_get_max_threads();
uint8_t mismatches_limit = 10;

std::lognormal_distribution<double> avg_input_rep_length_distrib(4.0, 2.0);
std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);
std::uniform_int_distribution<uint8_t> alphabet_size_distrib(1, 254);
std::uniform_int_distribution<uint8_t> uchar_distrib(0, 255);
std::uniform_int_distribution<uint32_t> input_size_distrib(1, 200000);
std::uniform_int_distribution<uint16_t> num_threads_distrib(1, max_num_threads);
std::lognormal_distribution<double> a_distrib(2.0, 3.0);

uint32_t input_size;
uint8_t alphabet_size;
std::string input;
uint32_t max_pattern_length;
uint32_t num_queries;

template <move_r_support support>
void test_move_rb()
{
    // choose a random input length
    input_size = input_size_distrib(gen);

    // choose a random alphabet size
    alphabet_size = alphabet_size_distrib(gen);

    // choose a random alphabet
    std::vector<uint8_t> alphabet;
    alphabet.reserve(alphabet_size);
    uint8_t uchar;

    for (uint32_t i = 0; i < alphabet_size; i++) {
        do {
            uchar = uchar_distrib(gen);
        } while (contains(alphabet, uchar));

        alphabet.push_back(uchar);
    }

    // choose a random input based on the alphabet
    std::uniform_int_distribution<uint8_t> char_idx_distrib(0, alphabet_size - 1);
    double avg_input_rep_length = 1.0 + avg_input_rep_length_distrib(gen);
    uint8_t cur_uchar = alphabet[char_idx_distrib(gen)];
    input.reserve(input_size);

    for (uint32_t i = 0; i < input_size; i++) {
        if (prob_distrib(gen) < 1 / avg_input_rep_length)
            cur_uchar = alphabet[char_idx_distrib(gen)];

        input.push_back(uchar_to_char(cur_uchar));
    }

    // build move-r and choose a random number of threads and balancing parameter, but always use libsais,
    // because there are bugs in Big-BWT that come through during fuzzing but not really in practice
    move_rb<support, char, uint32_t> index(input, {
        .mode = _suffix_array,
        .num_threads = num_threads_distrib(gen),
        .a = std::min<uint16_t>(2 + a_distrib(gen), 32767)
    });

    // generate patterns from the input and test count- and locate queries
    std::uniform_int_distribution<uint32_t> pattern_pos_distrib(0, input_size - 1);
    max_pattern_length = std::min<uint32_t>(10000, std::max<uint32_t>(100, input_size / 1000));
    std::uniform_int_distribution<uint32_t> pattern_length_distrib(1, max_pattern_length);
    num_queries = std::min<uint32_t>(1000, std::max<uint32_t>(100, input_size / 1000));
    num_queries = std::max<uint32_t>(1, num_queries / max_num_threads);

    #pragma omp parallel num_threads(max_num_threads)
    {
        std::random_device rd_thr;
        std::mt19937 gen_thr(rd_thr());
        uint32_t pattern_pos;
        uint32_t pattern_length;
        std::string pattern;
        std::vector<uint32_t> correct_occurrences;
        std::vector<uint32_t> occurrences;

        for (uint32_t cur_query = 0; cur_query < num_queries; cur_query++) {
            pattern_pos = pattern_pos_distrib(gen_thr);
            pattern_length = std::min<uint32_t>(input_size - pattern_pos, pattern_length_distrib(gen_thr));
            uint32_t max_mismatches = std::uniform_int_distribution<uint32_t>(0,
                std::min<uint32_t>(mismatches_limit, pattern_length - 1))(gen_thr);
            no_init_resize(pattern, pattern_length);

            for (uint32_t i = 0; i < pattern_length; i++)
                pattern[i] = input[pattern_pos + i];

            for (uint32_t i = 0; i <= input_size - pattern_length; i++) {
                uint32_t mismatches = 0;

                for (uint32_t j = 0; j < pattern_length; j++) {
                    if (input[i + j] != pattern[j]) mismatches++;
                    if (mismatches > max_mismatches) break;
                }

                if (mismatches <= max_mismatches) correct_occurrences.emplace_back(i);
            }

            occurrences = index.locate_with_mismatches(pattern, suffix_filter_scheme(max_mismatches, HAMMING_DISTANCE));
            EXPECT_EQ(occurrences.size(), correct_occurrences.size());
            ips4o::sort(occurrences.begin(), occurrences.end());
            EXPECT_EQ(occurrences, correct_occurrences);

            correct_occurrences.clear();
            occurrences.clear();
        }
    }

    input.clear();
}

TEST(test_move_rb, fuzzy_test)
{
    auto start_time = now();

    while (time_diff_min(start_time, now()) < 60) {
        if (prob_distrib(gen) < 0.5) {
          test_move_rb<_locate_move>();
        } else {
            test_move_rb<_locate_rlzsa>();
        }
    }
}