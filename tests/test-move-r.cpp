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
#include "test-index-common.hpp"

/**
 * @brief index adapter for the generic index fuzz-test (test_index_instance) that tests the move_r index;
 * besides count/locate it verifies the full-text queries (revert, SA and BWT) and the incremental query()
 * interface
 * @tparam pos_t index integer type
 * @tparam support type of locate support
 * @tparam sym_t symbol type of the input (a one-byte type for a byte alphabet, or a wider integer alphabet)
 * @tparam inp_t sequence type of the input (std::string for char, else std::vector<sym_t>)
 */
template <typename pos_t_, move_r_support support, typename sym_t_, typename inp_t_>
struct move_r_adapter {
    using pos_t = pos_t_;
    using sym_t = sym_t_;
    using inp_t = inp_t_;
    using index_t = move_r<support, sym_t_, pos_t_>;

    // move_r is verified via its full-text queries instead of a serialize/load round-trip
    static constexpr bool serializable = false;

    static std::string name() { return "move_r"; }

    // a byte alphabet (a one-byte symbol type) reserves the two largest values (for the sentinel and the
    // BWT-escape), which bounds the number of distinct symbols to at most 254
    static sym_t min_sym() { return std::numeric_limits<sym_t>::min(); }
    static sym_t max_sym()
    {
        return sizeof(sym_t) == 1 ?
            sym_t(std::numeric_limits<sym_t>::max() - 2) : std::numeric_limits<sym_t>::max();
    }

    static uint64_t max_pattern_length(uint64_t n) { return std::min<uint64_t>(10000, std::max<uint64_t>(100, n / 1000)); }
    static uint64_t num_queries(uint64_t n) { return std::min<uint64_t>(10000, std::max<uint64_t>(1000, n / 100)); }

    template <typename gen_t>
    static index_t build(inp_t& input, gen_t& gen, uint16_t max_num_threads)
    {
        std::uniform_int_distribution<uint16_t> num_threads_distrib(1, max_num_threads);
        std::lognormal_distribution<double> a_distrib(2.0, 3.0);

        // always use libsais, because there are bugs in Big-BWT that come through during fuzzing but not really in practice
        return index_t(input, {
            .mode = _suffix_array,
            .num_threads = num_threads_distrib(gen),
            .a = std::min<uint16_t>(2 + a_distrib(gen), 32767)
        });
    }

    static void after_build(index_t&, const inp_t&) { }

    static uint64_t count(index_t& index, inp_t& pattern) { return index.count(pattern); }
    static std::vector<pos_t> locate(index_t& index, inp_t& pattern) { return index.locate(pattern); }

    /**
     * @brief verifies count(), the one-shot locate() and the incremental query() interface (which is
     * exercised either fully via locate() or by drawing next_occ() until a random cutoff)
     */
    template <typename gen_t>
    static void verify_pattern(
        index_t& index, inp_t& pattern, const std::vector<pos_t>& correct_occurrences, gen_t& gen)
    {
        std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);

        // count() and the one-shot locate()
        verify_count_locate_basic<move_r_adapter>(index, pattern, correct_occurrences);

        // incrementally backward-search the pattern; prepend() is transactional and returns false (leaving
        // the context unchanged) as soon as a symbol cannot extend the match, so the pattern occurs in the
        // input iff every prepend() succeeds
        auto query = index.query();
        bool fully_matched = true;
        for (int32_t i = pattern.size() - 1; i >= 0; i--) {
            if (!query.prepend(pattern[i])) {
                fully_matched = false;
                break;
            }
        }

        // if the full pattern does not occur, the context holds the longest matching suffix instead, so
        // only the absence of exact occurrences can be verified here
        if (!fully_matched) {
            EXPECT_TRUE(correct_occurrences.empty());
            return;
        }

        EXPECT_EQ(query.num_occ(), correct_occurrences.size());
        std::vector<pos_t> occurrences;

        for (pos_t i = 0; i < query.num_occ(); i++) {
            if (prob_distrib(gen) < 1 / (double) query.num_occ()) {
                std::vector<pos_t> remaining_occurrences = query.locate();
                occurrences.insert(occurrences.end(),
                    remaining_occurrences.begin(), remaining_occurrences.end());
                break;
            }

            occurrences.emplace_back(query.next_occ());
        }

        ips4o::sort(occurrences.begin(), occurrences.end());
        EXPECT_EQ(occurrences, correct_occurrences);
    }

    /**
     * @brief verifies the revert-, SA- and BWT-queries (both the range- and the single-position variants)
     * against naive references
     */
    template <typename gen_t>
    static void verify_build(index_t& index, const inp_t& input, gen_t& gen, uint16_t max_num_threads)
    {
        std::uniform_int_distribution<uint16_t> num_threads_distrib(1, max_num_threads);
        pos_t input_size = input.size();

        // revert the index and compare the output with the input
        inp_t input_reverted = index.revert_range({ .num_threads = num_threads_distrib(gen) });

        #pragma omp parallel for num_threads(max_num_threads)
        for (pos_t i = 0; i < input_size; i++)
            EXPECT_EQ(input[i], input_reverted[i]);

        // retrieve the suffix array and compare it with the correct suffix array
        std::vector<int64_t> suffix_array = compute_reference_sa<sym_t, inp_t>(input, input_size, max_num_threads);
        std::vector<pos_t> suffix_array_retrieved = index.SA_range({ .num_threads = num_threads_distrib(gen) });

        #pragma omp parallel for num_threads(max_num_threads)
        for (pos_t i = 0; i <= input_size; i++)
            EXPECT_EQ(suffix_array[i], suffix_array_retrieved[i]);

        // compute each suffix array value separately and check if it is correct
        #pragma omp parallel for num_threads(max_num_threads)
        for (pos_t i = 0; i <= input_size; i++)
            EXPECT_EQ(index.SA(i), suffix_array[i]);

        // retrieve the bwt and compare it with the correct bwt
        inp_t bwt;
        no_init_resize(bwt, input_size + 1);

        #pragma omp parallel for num_threads(max_num_threads)
        for (pos_t i = 0; i <= input_size; i++)
            bwt[i] = suffix_array[i] == 0 ? sym_t(0) : input[suffix_array[i] - 1];

        inp_t bwt_retrieved = index.BWT_range({ .num_threads = num_threads_distrib(gen) });

        #pragma omp parallel for num_threads(max_num_threads)
        for (pos_t i = 0; i <= input_size; i++)
            EXPECT_EQ(bwt[i], bwt_retrieved[i]);

        // compute each bwt character separately and check if it is correct
        #pragma omp parallel for num_threads(max_num_threads)
        for (pos_t i = 0; i <= input_size; i++)
            EXPECT_EQ(index.BWT(i), bwt[i]);
    }
};

/**
 * @brief builds and verifies one random move_r index over the given input value type via the generic
 * framework, picking a random locate-support type and a random 32-/64-bit index integer type
 * @tparam sym_t symbol type of the input (a one-byte type for a byte alphabet, or a wider integer alphabet)
 * @tparam inp_t sequence type of the input (std::string for char, else std::vector<sym_t>)
 * @tparam gen_t random number generator type
 * @param gen random number generator
 * @param max_num_threads the maximum number of threads to use
 * @param min_input_size minimum input length
 * @param max_input_size maximum input length
 */
template <typename sym_t, typename inp_t, typename gen_t>
static void run_move_r_once(
    gen_t& gen, uint16_t max_num_threads, uint64_t min_input_size, uint64_t max_input_size)
{
    std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);
    double val = prob_distrib(gen);
    bool use_64_bit = prob_distrib(gen) < 0.5;

    auto run = [&]<move_r_support support>() {
        if (use_64_bit) test_index_instance<move_r_adapter<uint64_t, support, sym_t, inp_t>>(gen, max_num_threads, min_input_size, max_input_size);
        else            test_index_instance<move_r_adapter<uint32_t, support, sym_t, inp_t>>(gen, max_num_threads, min_input_size, max_input_size);
    };

    if (val < 0.5) run.template operator()<_locate_move>();
    else           run.template operator()<_locate_rlzsa>();
}

std::random_device rd;
std::mt19937 gen(rd());
uint16_t max_num_threads = omp_get_max_threads();
uint64_t min_input_size = 1;
uint64_t max_input_size = 200000;

// move_r over every supported input value type (byte alphabets: char/uint8_t/int8_t; integer alphabets:
// (u)int16/32/64), rotating over the locate-support types and the 32-/64-bit index integer type
TEST(test_move_r, fuzzy_test)
{
    std::uniform_int_distribution<int> sym_type_distrib(0, 8);
    auto start_time = now();

    while (time_diff_min(start_time, now()) < 60) {
        switch (sym_type_distrib(gen)) {
            case 0: run_move_r_once<char,     std::string>          (gen, max_num_threads, min_input_size, max_input_size); break;
            case 1: run_move_r_once<uint8_t,  std::vector<uint8_t>> (gen, max_num_threads, min_input_size, max_input_size); break;
            case 2: run_move_r_once<int8_t,   std::vector<int8_t>>  (gen, max_num_threads, min_input_size, max_input_size); break;
            case 3: run_move_r_once<uint16_t, std::vector<uint16_t>>(gen, max_num_threads, min_input_size, max_input_size); break;
            case 4: run_move_r_once<int16_t,  std::vector<int16_t>> (gen, max_num_threads, min_input_size, max_input_size); break;
            case 5: run_move_r_once<uint32_t, std::vector<uint32_t>>(gen, max_num_threads, min_input_size, max_input_size); break;
            case 6: run_move_r_once<int32_t,  std::vector<int32_t>> (gen, max_num_threads, min_input_size, max_input_size); break;
            case 7: run_move_r_once<uint64_t, std::vector<uint64_t>>(gen, max_num_threads, min_input_size, max_input_size); break;
            case 8: run_move_r_once<int64_t,  std::vector<int64_t>> (gen, max_num_threads, min_input_size, max_input_size); break;
        }
    }
}
