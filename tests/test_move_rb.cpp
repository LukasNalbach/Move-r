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
void test_move_r()
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
    num_queries = std::min<uint32_t>(10000, std::max<uint32_t>(1000, input_size / 100));
    num_queries = std::max<uint32_t>(1, num_queries / max_num_threads);

    #pragma omp parallel num_threads(max_num_threads)
    {
        std::random_device rd_thr;
        std::mt19937 gen_thr(rd());
        uint32_t pattern_pos;
        uint32_t pattern_length;
        std::string pattern;
        std::vector<uint32_t> correct_occurrences;
        std::vector<uint32_t> occurrences;
        bool match;

        for (uint32_t cur_query = 0; cur_query < num_queries; cur_query++) {
            pattern_pos = pattern_pos_distrib(gen_thr);
            pattern_length = std::min<uint32_t>(input_size - pattern_pos, pattern_length_distrib(gen));
            no_init_resize(pattern, pattern_length);

            for (uint32_t i = 0; i < pattern_length; i++)
                pattern[i] = input[pattern_pos + i];

            for (uint32_t i = 0; i <= input_size - pattern_length; i++) {
                match = true;

                for (uint32_t j = 0; j < pattern_length; j++) {
                    if (input[i + j] != pattern[j]) {
                        match = false;
                        break;
                    }
                }

                if (match) correct_occurrences.emplace_back(i);
            }

            uint32_t start_pos = std::uniform_int_distribution<uint32_t>(0, pattern_length - 1)(gen);
            std::vector<direction> dirs;
            dirs.reserve(pattern_length);
            for (uint32_t i = 0; i <= start_pos; i++) dirs.emplace_back(LEFT);
            for (uint32_t i = start_pos + 1; i < pattern_length; i++) dirs.emplace_back(RIGHT);
            std::shuffle(dirs.begin(), dirs.end(), gen);
            auto query = index.query();
            uint32_t pos_prev_sym = start_pos;
            uint32_t pos_next_sym = start_pos + 1;

            for (uint32_t i = 0; i < pattern_length; i++) {
                if (dirs[i] == LEFT) {
                    EXPECT_TRUE(query.prepend(pattern[pos_prev_sym]));
                    pos_prev_sym--;
                } else {
                    EXPECT_TRUE(query.append(pattern[pos_next_sym]));
                    pos_next_sym++;
                }
            }
            
            EXPECT_EQ(query.num_occ(), correct_occurrences.size());

            for (uint32_t i = 0; i < query.num_occ(); i++) {
                if (prob_distrib(gen) < 1 / (double) query.num_occ()) {
                    std::vector<uint32_t> remaining_occurrences = query.locate();
                    occurrences.insert(occurrences.end(),
                        remaining_occurrences.begin(), remaining_occurrences.end());
                    break;
                }

                occurrences.emplace_back(query.next_occ());
            }

            ips4o::sort(occurrences.begin(), occurrences.end());
            EXPECT_EQ(occurrences, correct_occurrences);
            correct_occurrences.clear();
            occurrences.clear();
        }
    }

    input.clear();
}

TEST(test_move_r, fuzzy_test)
{
    auto start_time = now();

    while (time_diff_min(start_time, now()) < 60) {
        if (prob_distrib(gen) < 0.5) {
          test_move_r<_locate_move>();
        } else {
            test_move_r<_locate_rlzsa>();
        }
    }
}