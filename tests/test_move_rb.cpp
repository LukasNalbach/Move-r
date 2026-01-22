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
char min_uchar = std::numeric_limits<char>::min();
char max_uchar = std::numeric_limits<char>::max() - 2;
uint64_t min_input_size = 1;
uint64_t max_input_size = 100000;
uint64_t max_pattern_length = 20;
uint8_t mismatches_limit = 10;

std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);
std::uniform_int_distribution<uint16_t> num_threads_distrib(1, max_num_threads);
std::lognormal_distribution<double> a_distrib(2.0, 3.0);

template <typename pos_t, move_r_support support>
void test_move_rb()
{
    // generate a random input
    std::string input = random_repetitive_input<std::string>(min_input_size, max_input_size, min_uchar, max_uchar);
    pos_t input_size = input.size();

    // build move-r and choose a random number of threads and balancing parameter, but always use libsais,
    // because there are bugs in Big-BWT that come through during fuzzing but not really in practice
    move_rb<support, char, pos_t> index(input, {
        .mode = _suffix_array,
        .num_threads = num_threads_distrib(gen),
        .a = std::min<uint16_t>(2 + a_distrib(gen), 32767)
    });

    // generate patterns from the input and test count- and locate queries
    std::uniform_int_distribution<pos_t> pattern_pos_distrib(0, input_size - 1);
    std::uniform_int_distribution<pos_t> pattern_length_distrib(1, std::min<pos_t>(max_pattern_length, input_size));
    pos_t num_queries = std::min<pos_t>(1000, std::max<pos_t>(100, input_size / 1000));
    num_queries = std::max<pos_t>(1, num_queries / max_num_threads);

    #pragma omp parallel num_threads(max_num_threads)
    {
        std::random_device rd_thr;
        std::mt19937 gen_thr(rd_thr());
        pos_t pattern_pos;
        pos_t pattern_length;
        std::string pattern;
        std::vector<aprx_occ_t<pos_t>> correct_occurrences;
        std::vector<aprx_occ_t<pos_t>> occurrences;

        for (pos_t cur_query = 0; cur_query < num_queries; cur_query++) {
            pattern_pos = pattern_pos_distrib(gen_thr);
            pattern_length = std::min<pos_t>(input_size - pattern_pos, pattern_length_distrib(gen_thr));
            pos_t max_mismatches = std::uniform_int_distribution<pos_t>(0,
                std::min<pos_t>(mismatches_limit, pattern_length - 1))(gen_thr);
            no_init_resize(pattern, pattern_length);
            for (pos_t i = 0; i < pattern_length; i++)
                pattern[i] = input[pattern_pos + i];

            for_each_constexpr<HAMMING_DISTANCE, EDIT_DISTANCE>([&](auto dist_metr){
                correct_occurrences = locate<pos_t, dist_metr>(input, pattern, max_mismatches);
                search_scheme_t scheme = suffix_filter_scheme(max_mismatches);
                EXPECT_EQ(index.template count<dist_metr>(pattern, scheme), correct_occurrences.size());
                occurrences = index.template locate<dist_metr>(pattern, scheme);
                EXPECT_EQ(occurrences.size(), correct_occurrences.size());
                ips4o::sort(occurrences.begin(), occurrences.end());
                EXPECT_EQ(occurrences, correct_occurrences);
                correct_occurrences.clear();
                occurrences.clear();
            });
        }
    }
}

TEST(test_move_rb, fuzzy_test)
{
    auto start_time = now();

    while (time_diff_min(start_time, now()) < 60) {
        if (prob_distrib(gen) < 0.5) {
            if (prob_distrib(gen) < 0.5) {
                test_move_rb<uint32_t, _locate_move>();
            } else {
                test_move_rb<uint64_t, _locate_move>();
            }
        } else {
            if (prob_distrib(gen) < 0.5) {
                test_move_rb<uint32_t, _locate_rlzsa>();
            } else {
                test_move_rb<uint64_t, _locate_rlzsa>();
            }
        }
    }
}