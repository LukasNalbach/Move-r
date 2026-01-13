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
#include <move_r/move_r.hpp>

std::random_device rd;
std::mt19937 gen(rd());
uint16_t max_num_threads = omp_get_max_threads();
char min_uchar = std::numeric_limits<char>::min();
char max_uchar = std::numeric_limits<char>::max() - 2;
uint64_t min_input_size = 1;
uint64_t max_input_size = 200000;

std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);
std::uniform_int_distribution<uint16_t> num_threads_distrib(1, max_num_threads);
std::lognormal_distribution<double> a_distrib(2.0, 3.0);

template <typename pos_t, move_r_support support>
void test_move_r()
{
    // generate a random input
    std::string input = random_repetitive_input<std::string>(min_input_size, max_input_size, min_uchar, max_uchar);
    pos_t input_size = input.size();

    // build move-r and choose a random number of threads and balancing parameter, but always use libsais,
    // because there are bugs in Big-BWT that come through during fuzzing but not really in practice
    move_r<support, char, pos_t> index(input, {
        .mode = _suffix_array,
        .num_threads = num_threads_distrib(gen),
        .a = std::min<uint16_t>(2 + a_distrib(gen), 32767)
    });

    // revert the index and compare the output with the input string
    std::string input_reverted = index.revert_range({ .num_threads = num_threads_distrib(gen) });

    #pragma omp parallel for num_threads(max_num_threads)
    for (pos_t i = 0; i < input_size; i++)
        EXPECT_EQ(input[i], input_reverted[i]);

    // retrieve the suffix array and compare it with the correct suffix array; if the input contains 0,
    // then temporarily remap the characters of the input string s.t. it does not contain 0
    std::vector<uint8_t> contains_uchar(256, 0);
    #pragma omp parallel for num_threads(max_num_threads)
    for (pos_t i = 0; i < input_size; i++)
        contains_uchar[char_to_uchar(input[i])] = 1;

    std::vector<uint8_t> map_uchar;
    std::vector<uint8_t> unmap_uchar;

    if (contains_uchar[0]) {
        map_uchar.resize(256, 0);
        unmap_uchar.resize(256, 0);
        uint8_t next_uchar = 1;

        for (uint16_t i = 0; i < 256; i++) {
            if (contains_uchar[i]) {
                map_uchar[i] = next_uchar;
                unmap_uchar[next_uchar] = i;
                next_uchar++;
            }
        }

        #pragma omp parallel for num_threads(max_num_threads)
        for (pos_t i = 0; i < input_size; i++)
            input[i] = uchar_to_char(map_uchar[char_to_uchar(input[i])]);
    }

    input.push_back(uchar_to_char((uint8_t)0));
    std::vector<int64_t> suffix_array;
    no_init_resize(suffix_array, input_size + 1);
    libsais64_omp((uint8_t*)&input[0], &suffix_array[0], input_size + 1, 0, NULL, max_num_threads);

    if (contains_uchar[0]) {
        #pragma omp parallel for num_threads(max_num_threads)
        for (pos_t i = 0; i < input_size; i++)
            input[i] = uchar_to_char(unmap_uchar[char_to_uchar(input[i])]);
    }

    input.resize(input_size);
    std::vector<pos_t> suffix_array_retrieved = index.SA_range({ .num_threads = num_threads_distrib(gen) });

    #pragma omp parallel for num_threads(max_num_threads)
    for (pos_t i = 0; i <= input_size; i++)
        EXPECT_EQ(suffix_array[i], suffix_array_retrieved[i]);

    // compute each suffix array value separately and check if it is correct
    #pragma omp parallel for num_threads(max_num_threads)
    for (pos_t i = 0; i <= input_size; i++)
        EXPECT_EQ(index.SA(i), suffix_array[i]);

    // retrieve the bwt and compare it with the correct bwt
    std::string bwt;
    no_init_resize(bwt, input_size + 1);
    
    #pragma omp parallel for num_threads(max_num_threads)
    for (pos_t i = 0; i <= input_size; i++)
        bwt[i] = suffix_array[i] == 0 ? 0 : input[suffix_array[i] - 1];

    std::string bwt_retrieved = index.BWT_range({ .num_threads = num_threads_distrib(gen) });

    #pragma omp parallel for num_threads(max_num_threads)
    for (pos_t i = 0; i <= input_size; i++)
        EXPECT_EQ(bwt[i], bwt_retrieved[i]);

    // compute each bwt character separately and check if it is correct
    #pragma omp parallel for num_threads(max_num_threads)
    for (pos_t i = 0; i <= input_size; i++)
        EXPECT_EQ(index.BWT(i), bwt[i]);

    // generate patterns from the input and test count- and locate queries
    std::uniform_int_distribution<pos_t> pattern_pos_distrib(0, input_size - 1);
    pos_t max_pattern_length = std::min<pos_t>(10000, std::max<pos_t>(100, input_size / 1000));
    std::uniform_int_distribution<pos_t> pattern_length_distrib(1, max_pattern_length);
    pos_t num_queries = std::min<pos_t>(10000, std::max<pos_t>(1000, input_size / 100));
    num_queries = std::max<pos_t>(1, num_queries / max_num_threads);

    #pragma omp parallel num_threads(max_num_threads)
    {
        std::random_device rd_thr;
        std::mt19937 gen_thr(rd_thr());
        pos_t pattern_pos;
        pos_t pattern_length;
        std::string pattern;
        std::vector<pos_t> correct_occurrences;
        std::vector<pos_t> occurrences;
        bool match;

        for (pos_t cur_query = 0; cur_query < num_queries; cur_query++) {
            pattern_pos = pattern_pos_distrib(gen_thr);
            pattern_length = std::min<pos_t>(input_size - pattern_pos, pattern_length_distrib(gen_thr));
            no_init_resize(pattern, pattern_length);

            for (pos_t i = 0; i < pattern_length; i++)
                pattern[i] = input[pattern_pos + i];

            for (pos_t i = 0; i <= input_size - pattern_length; i++) {
                match = true;

                for (pos_t j = 0; j < pattern_length; j++) {
                    if (input[i + j] != pattern[j]) {
                        match = false;
                        break;
                    }
                }

                if (match) correct_occurrences.emplace_back(i);
            }
            
            EXPECT_EQ(index.count(pattern), correct_occurrences.size());
            occurrences = index.locate(pattern);
            ips4o::sort(occurrences.begin(), occurrences.end());
            EXPECT_EQ(occurrences, correct_occurrences);
            occurrences.clear();
            auto query = index.query();
            for (int32_t i = pattern.size() - 1; i >= 0; i--)
                query.prepend(pattern[i]);
            EXPECT_EQ(query.num_occ(), correct_occurrences.size());

            if constexpr (support == _locate_lzendsa) {
                occurrences = query.locate();
            } else {
                for (int32_t i = 0; i < query.num_occ(); i++) {
                    if (prob_distrib(gen_thr) < 1 / (double) query.num_occ()) {
                        std::vector<pos_t> remaining_occurrences = query.locate();
                        occurrences.insert(occurrences.end(),
                            remaining_occurrences.begin(), remaining_occurrences.end());
                        break;
                    }
    
                    occurrences.emplace_back(query.next_occ());
                }
            }

            ips4o::sort(occurrences.begin(), occurrences.end());
            EXPECT_EQ(occurrences, correct_occurrences);
            correct_occurrences.clear();
            occurrences.clear();
        }
    }
}

TEST(test_move_r, fuzzy_test)
{
    auto start_time = now();

    while (time_diff_min(start_time, now()) < 60) {
        double val = prob_distrib(gen);

        if (val < 1.0 / 3.0) {
            if (prob_distrib(gen) < 0.5) {
                test_move_r<uint32_t, _locate_move>();
            } else {
                test_move_r<uint64_t, _locate_move>();
            }
        } else if (val < 2.0 / 3.0) {
            if (prob_distrib(gen) < 0.5) {
                test_move_r<uint32_t, _locate_rlzsa>();
            } else {
                test_move_r<uint64_t, _locate_rlzsa>();
            }
        } else {
            if (prob_distrib(gen) < 0.5) {
                test_move_r<uint32_t, _locate_lzendsa>();
            } else {
                test_move_r<uint64_t, _locate_lzendsa>();
            }
        }
    }
}