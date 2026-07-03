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
#include <gtl/phmap.hpp>
#include <move_rb/move_rb.hpp>
#include <misc/apm.hpp>

std::random_device rd;
std::mt19937 gen(rd());
uint16_t max_num_threads = omp_get_max_threads();
uint64_t min_input_size = 1;
uint64_t max_input_size = 10'000;
uint64_t max_pattern_length = 100;
uint64_t mismatches_limit = 20;

std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);
std::uniform_int_distribution<uint16_t> num_threads_distrib(1, max_num_threads);
std::lognormal_distribution<double> a_distrib(2.0, 3.0);

/**
 * @brief builds a move_rb index of a random input and verifies its approximate (hamming-/edit-distance)
 * count- and locate-queries against naive references
 * @tparam sym_t symbol type of the input (a one-byte type: char, uint8_t or int8_t)
 * @tparam inp_t sequence type of the input (std::string for char, else std::vector<sym_t>)
 * @tparam pos_t index integer type
 * @tparam support type of locate support
 */
template <typename sym_t, typename inp_t, typename pos_t, move_r_support support>
void test_move_rb()
{
    // a byte alphabet reserves the two largest values (for the sentinel and the BWT-escape)
    sym_t min_sym = std::numeric_limits<sym_t>::min();
    sym_t max_sym = std::numeric_limits<sym_t>::max() - 2;

    // generate a random input
    inp_t input = random_repetitive_input<inp_t>(min_input_size, max_input_size, min_sym, max_sym);
    pos_t input_size = input.size();

    // build move-r and choose a random number of threads and balancing parameter, but always use libsais,
    // because there are bugs in Big-BWT that come through during fuzzing but not really in practice
    move_rb<support, sym_t, pos_t> index(input, {
        .mode = _suffix_array,
        .num_threads = num_threads_distrib(gen),
        .a = std::min<uint16_t>(2 + a_distrib(gen), 32767)
    });

    // a read-only view for the naive reference algorithms; created AFTER the build, which reallocates input's
    // buffer (its contents are preserved) -- a view captured earlier would dangle
    std::span<const sym_t> input_span(input);

    // generate patterns from the input and test count- and locate queries
    std::uniform_int_distribution<pos_t> pattern_pos_distrib(0, input_size - 1);
    std::uniform_int_distribution<pos_t> pattern_length_distrib(1, std::min<pos_t>(max_pattern_length, input_size));
    pos_t num_queries = std::min<pos_t>(1'000, std::max<pos_t>(100, input_size / 10));
    num_queries = std::max<pos_t>(1, num_queries / max_num_threads);

    #pragma omp parallel for num_threads(max_num_threads), schedule(dynamic,1)
    for (pos_t cur_query = 0; cur_query < num_queries; cur_query++) {
        std::random_device rd_thr;
        std::mt19937 gen_thr(rd_thr());
        pos_t pattern_pos = pattern_pos_distrib(gen_thr);
        pos_t pattern_length = std::min<pos_t>(input_size - pattern_pos, pattern_length_distrib(gen_thr));
        inp_t pattern;
        no_init_resize(pattern, pattern_length);
        for (pos_t i = 0; i < pattern_length; i++)
            pattern[i] = input[pattern_pos + i];

        mutate_pattern<pos_t>(pattern, input, gen_thr);
        pattern_length = pattern.size();
        std::span<const sym_t> pattern_span(pattern);
        pos_t max_mismatches = std::uniform_int_distribution<pos_t>(0,
            std::max<int32_t>(std::min<int32_t>(pattern_length / 2 - 1, mismatches_limit), 0))(gen_thr);

        for_each_constexpr<HAMMING_DISTANCE, EDIT_DISTANCE>([&](auto dist_metr){
            auto correct_occurrences = locate<pos_t, dist_metr>(input_span, pattern_span, max_mismatches);
            const search_scheme_t scheme = min_u_scheme(max_mismatches);
            auto occurrences = index.template locate<dist_metr>(pattern, scheme);
            ips4o::sort(occurrences.begin(), occurrences.end());

            if constexpr (dist_metr == HAMMING_DISTANCE) {
                EXPECT_EQ(index.count_hamming_dist(pattern, scheme), correct_occurrences.size());
                EXPECT_EQ(occurrences.size(), correct_occurrences.size());
                EXPECT_EQ(occurrences, correct_occurrences);
            } else {
                gtl::flat_hash_map<std::span<const sym_t>, pos_t, span_hash, span_eq> dist_cache;
                for (const auto& occ : occurrences) {
                    if (occ.pos + occ.len > input_size) { ADD_FAILURE(); continue; }
                    std::span<const sym_t> match = input_span.subspan(occ.pos, occ.len);
                    auto [it, inserted] = dist_cache.try_emplace(match, 0);
                    if (inserted) it->second = edit_dist_bounded<pos_t>(pattern_span, match, max_mismatches);
                    if constexpr (dist_metr == HAMMING_DISTANCE) EXPECT_EQ(it->second, occ.err);
                    else EXPECT_TRUE(it->second == occ.err && occ.err <= max_mismatches);
                }
                EXPECT_TRUE(verify_edit_distance_coverage(occurrences, correct_occurrences, max_mismatches));
                filter_edit_distance_occurrences(occurrences, max_mismatches);
                EXPECT_TRUE(verify_edit_distance_coverage(occurrences, correct_occurrences, max_mismatches));
            }

            std::vector<aprx_occ_t<pos_t, CIGAR>> cigar_results;
            index.template locate<dist_metr, CIGAR>(pattern, scheme, [&](aprx_occ_t<pos_t, CIGAR> occ){
                cigar_results.emplace_back(std::move(occ));
            });

            std::vector<aprx_occ_t<pos_t>> cigar_occurrences;
            gtl::flat_hash_map<std::span<const sym_t>, pos_t, span_hash, span_eq> dist_cache;
            for (const auto& occ : cigar_results) {
                cigar_occurrences.push_back({occ.pos, occ.len, occ.err});
                if (occ.pos + occ.len > input_size) { ADD_FAILURE(); continue; }
                std::span<const sym_t> match = input_span.subspan(occ.pos, occ.len);
                EXPECT_EQ(num_cigar_edits(pattern_span.data(), pattern_span.size(), match.data(), match.size(), occ.cigar), occ.err);
                auto [it, inserted] = dist_cache.try_emplace(match, 0);

                if (inserted) {
                    if constexpr (dist_metr == HAMMING_DISTANCE) {
                        EXPECT_EQ(match.size(), pattern_span.size());
                        it->second = hamming_dist_bounded<pos_t>(pattern_span, match, max_mismatches);
                    } else {
                        it->second = edit_dist_bounded<pos_t>(pattern_span, match, max_mismatches);
                    }
                }

                EXPECT_TRUE(it->second == occ.err && occ.err <= max_mismatches);
            }

            ips4o::sort(cigar_occurrences.begin(), cigar_occurrences.end());
            if constexpr (dist_metr == HAMMING_DISTANCE) {
                EXPECT_EQ(cigar_occurrences, correct_occurrences);
            } else {
                EXPECT_TRUE(verify_edit_distance_coverage(cigar_occurrences, correct_occurrences, max_mismatches));
                filter_edit_distance_occurrences(cigar_occurrences, max_mismatches);
                EXPECT_TRUE(verify_edit_distance_coverage(cigar_occurrences, correct_occurrences, max_mismatches));
            }
        });
    }
}

/**
 * @brief builds and verifies one random move_rb index over the given input value type, picking a random
 * locate-support type and a random 32-/64-bit index integer type
 * @tparam sym_t symbol type of the input (a one-byte type: char, uint8_t or int8_t)
 * @tparam inp_t sequence type of the input (std::string for char, else std::vector<sym_t>)
 */
template <typename sym_t, typename inp_t>
void run_move_rb_once()
{
    bool use_64_bit = prob_distrib(gen) < 0.5;

    auto run = [&]<move_r_support support>() {
        if (use_64_bit) test_move_rb<sym_t, inp_t, uint64_t, support>();
        else            test_move_rb<sym_t, inp_t, uint32_t, support>();
    };

    if (prob_distrib(gen) < 0.5) run.template operator()<_locate_move>();
    else                         run.template operator()<_locate_rlzsa>();
}

// move_rb over every supported input value type (the byte alphabets char/uint8_t/int8_t), rotating over the
// locate-support types and the 32-/64-bit index integer type
TEST(test_move_rb, fuzzy_test)
{
    std::uniform_int_distribution<int> sym_type_distrib(0, 2);
    auto start_time = now();

    while (time_diff_min(start_time, now()) < 60) {
        switch (sym_type_distrib(gen)) {
            case 0: run_move_rb_once<char,    std::string>();          break;
            case 1: run_move_rb_once<uint8_t, std::vector<uint8_t>>(); break;
            case 2: run_move_rb_once<int8_t,  std::vector<int8_t>>();  break;
        }
    }
}