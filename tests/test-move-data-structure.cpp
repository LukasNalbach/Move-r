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
#include <misc/strings.hpp>
#include "test-progress.hpp"
#include <move_data_structure/move_data_structure.hpp>
#include <move_data_structure/move_data_structure_l_.hpp>

static std::mt19937 gen(std::random_device{}());
static uint16_t max_num_threads = omp_get_max_threads();

// a random disjoint interval sequence together with the move data structure built from it
struct random_mds {
    move_data_structure<uint32_t> mds;
    std::vector<std::pair<uint32_t, uint32_t>> interval_sequence;
    uint32_t input_size;
    uint32_t num_intervals;
    uint16_t a; // the balancing parameter used for the construction
};

// generates a random disjoint interval sequence (a random permutation of random-length input intervals into
// output intervals) and builds a move data structure from it
static random_mds build_random_mds()
{
    std::lognormal_distribution<double> avg_interval_length_distrib(4.0, 2.0);
    std::lognormal_distribution<double> a_distrib(2.0, 3.0);

    uint32_t input_size = random_log_uniform_size(550000, 1800000, gen);
    uint32_t num_intervals = 0;
    std::vector<std::pair<uint32_t, uint32_t>> interval_sequence;
    std::vector<uint32_t> interval_permutation;

    // choose input intervals of random lengths
    std::uniform_int_distribution<uint32_t> interval_length_distrib(1, 2 * avg_interval_length_distrib(gen));
    uint32_t interval_length_prefix_sum = 0;
    interval_sequence.emplace_back(std::make_pair(0, 0));

    while (interval_length_prefix_sum < input_size) {
        interval_length_prefix_sum = std::min<uint32_t>(
            interval_length_prefix_sum + interval_length_distrib(gen), input_size);
        interval_sequence.emplace_back(std::make_pair(interval_length_prefix_sum, 0));
        num_intervals++;
    }

    // permute the input intervals randomly into the output intervals
    no_init_resize(interval_permutation, num_intervals);

    #pragma omp parallel for num_threads(max_num_threads)
    for (uint32_t i = 0; i < num_intervals; i++)
        interval_permutation[i] = i;

    std::shuffle(interval_permutation.begin(), interval_permutation.end(), gen);
    interval_length_prefix_sum = 0;

    for (uint32_t i = 0; i < num_intervals; i++) {
        interval_sequence[interval_permutation[i]].second = interval_length_prefix_sum;
        interval_length_prefix_sum += interval_sequence[interval_permutation[i] + 1].first - interval_sequence[interval_permutation[i]].first;
    }

    interval_sequence.pop_back();

    uint16_t a = std::min<uint16_t>(2 + a_distrib(gen), 32767);
    move_data_structure<uint32_t> mds(interval_sequence, input_size, { .num_threads = max_num_threads, .a = a });
    return random_mds { std::move(mds), std::move(interval_sequence), input_size, num_intervals, a };
}

// the invariants of the resulting interval sequence that move_data_structure::construction used to check itself
// in debug builds (verify_correctness), moved here and verified through the public interface
static void verify_mds_structure(uint64_t)
{
    random_mds r = build_random_mds();
    move_data_structure<uint32_t>& mds = r.mds;

    #pragma omp parallel for num_threads(max_num_threads)
    for (uint32_t i = 0; i < mds.num_intervals(); i++) {
        // the input interval starting positions ascend
        EXPECT_LT(mds.p(i), mds.p(i + 1)) << "input interval starting positions must ascend (i = " << i << ")";

        // D_idx maps each output-interval start q(i) into the input interval that contains it
        EXPECT_TRUE(mds.p(mds.idx(i)) <= mds.q(i) && mds.q(i) < mds.p(mds.idx(i) + 1))
            << "D_idx does not map q(i) into its input interval (i = " << i << ")";

        // no output interval overlaps more than 2a input intervals (the a-balancing invariant)
        uint32_t j = bin_search_min_geq<uint32_t>(mds.q(i), 0, mds.num_intervals() - 1, [&mds](auto x) { return mds.p(x); });
        uint16_t num_intervals_in_output_interval = 0;

        while (mds.p(j) < mds.q(i) + (mds.p(i + 1) - mds.p(i))) {
            num_intervals_in_output_interval++;
            j++;
        }

        EXPECT_LE(num_intervals_in_output_interval, 2 * r.a) << "a-heavy output interval (i = " << i << ")";
    }
}

// the move queries: mds.move() must reproduce the interval sequence's mapping, and report the correct
// input-interval index of the output value
static void verify_mds_move(uint64_t)
{
    random_mds r = build_random_mds();
    move_data_structure<uint32_t>& mds = r.mds;
    uint32_t input_size = r.input_size;
    std::vector<std::pair<uint32_t, uint32_t>>& interval_sequence = r.interval_sequence;

    // roughly evenly distributed input values in [0, input_size)
    uint32_t avg_step_size = std::max<uint32_t>(2, input_size / 10000);
    std::uniform_int_distribution<uint32_t> step_size_distrib(avg_step_size / 1.5, 1.5 * avg_step_size);

    #pragma omp parallel for num_threads(max_num_threads)
    for (uint32_t i = 0; i < input_size; i += step_size_distrib(gen)) {
        std::pair<uint32_t, uint32_t> ix_mds { i,
            bin_search_max_leq<uint32_t>(i, 0, mds.num_intervals() - 1, [&mds](uint32_t x) { return mds.p(x); }) };

        std::pair<uint32_t, uint32_t> ix_is { i,
            bin_search_max_leq<uint32_t>(i, 0, r.num_intervals - 1, [&interval_sequence](uint32_t x) { return interval_sequence[x].first; }) };

        ix_mds = mds.move(ix_mds);
        ix_is.first = interval_sequence[ix_is.second].second + (ix_is.first - interval_sequence[ix_is.second].first);
        EXPECT_EQ(ix_mds.first, ix_is.first);

        // the index of the input interval containing the output value must be correct
        EXPECT_TRUE(mds.p(ix_mds.second) <= ix_mds.first && ix_mds.first < mds.p(ix_mds.second + 1));
    }
}

TEST(test_move_data_structure, structure) { run_fuzz("move-data-structure", { { "structure", verify_mds_structure } }, fuzz_iterations(1200)); }
TEST(test_move_data_structure, move)      { run_fuzz("move-data-structure", { { "move",      verify_mds_move } }, fuzz_iterations(1200)); }