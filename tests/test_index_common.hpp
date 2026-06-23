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

// shared framework for the index fuzz tests (move_r, move_r_int, lzendsa, rlzsa): reference suffix array
// construction, pattern generation, count-/locate-occurrence verification, and a single generic driver
// (test_index_instance) that builds a random index via an adapter and verifies its queries. Each index type
// provides a small adapter describing how to build it, how to run count/locate, and which extra checks apply
// (full-text SA/BWT/revert/query, suffix-array random access, or none).

#include <omp.h>

#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <ips4o.hpp>
#include <libsais64.h>
#include <gtl/phmap.hpp>

#include <misc/strings.hpp>
#include <misc/utils.hpp>

// the standalone SA-based indexes (lzendsa, rlzsa, ...) are designed for byte values >= 2: the value 0 is
// reserved as the suffix-array sentinel and 1 is used internally as a BWT-escape, so their inputs are drawn
// from the remaining alphabet (matching how the indexes are used on real text)
static constexpr char sa_index_min_uchar = 2;
static constexpr char sa_index_max_uchar = std::numeric_limits<char>::max() - 2;

/**
 * @brief computes the suffix array of input with an appended 0-sentinel (i.e. over the n + 1 suffixes), which
 * is the suffix array the indexes encode internally; the symbols are remapped to an effective alphabet so
 * that libsais sees no 0-byte (other than the sentinel) regardless of the symbols occurring in input
 *
 * one-byte symbols (char, uint8_t, int8_t) are ordered by their unsigned byte value, matching how move_r (and
 * the standalone indexes) order a byte alphabet; wider symbol types are ordered by their native value
 * @tparam sym_t symbol type of the input (a one-byte type for a byte alphabet, or a wider integer alphabet)
 * @tparam inp_t sequence type of the input (e.g. std::string or std::vector<int32_t>)
 * @param input the input sequence (taken by value, as it is modified internally)
 * @param input_size the length of the input
 * @param num_threads the number of threads to use for the (parallel) suffix array construction
 * @return the suffix array of input + the 0-sentinel (of size input_size + 1)
 */
template <typename sym_t, typename inp_t>
static std::vector<int64_t> compute_reference_sa(inp_t input, uint64_t input_size, uint16_t num_threads)
{
    std::vector<int64_t> suffix_array;
    no_init_resize(suffix_array, input_size + 1);

    if constexpr (sizeof(sym_t) == 1) {
        // if the input contains a 0-byte, remap the symbols (order-preserving over the unsigned byte values,
        // so the suffix array is unchanged) so that the input does not contain 0 before appending the sentinel
        std::vector<uint8_t> contains_uchar(256, 0);

        #pragma omp parallel for num_threads(num_threads)
        for (uint64_t i = 0; i < input_size; i++)
            contains_uchar[sym_to_uchar(input[i])] = 1;

        if (contains_uchar[0]) {
            std::vector<uint8_t> map_uchar(256, 0);
            uint8_t next_uchar = 1;

            for (uint16_t i = 0; i < 256; i++)
                if (contains_uchar[i]) map_uchar[i] = next_uchar++;

            #pragma omp parallel for num_threads(num_threads)
            for (uint64_t i = 0; i < input_size; i++)
                input[i] = static_cast<sym_t>(map_uchar[sym_to_uchar(input[i])]);
        }

        input.push_back(sym_t(0));
        libsais64_omp((uint8_t*) &input[0], &suffix_array[0], input_size + 1, 0, NULL, num_threads);
    } else {
        // remap the symbols to the effective alphabet [1, alphabet_size] and build the suffix array of the
        // remapped sequence (with an appended 0-sentinel)
        std::vector<sym_t> alphabet(input.begin(), input.end());
        ips4o::parallel::sort(alphabet.begin(), alphabet.end());
        alphabet.resize(std::unique(alphabet.begin(), alphabet.end()) - alphabet.begin());
        gtl::flat_hash_map<sym_t, int64_t> map_int;
        int64_t sym_cur = 1;

        for (uint64_t i = 0; i < alphabet.size(); i++)
            map_int[alphabet[i]] = sym_cur++;

        std::vector<int64_t> input_libsais;
        no_init_resize(input_libsais, input_size + 1);

        #pragma omp parallel for num_threads(num_threads)
        for (uint64_t i = 0; i < input_size; i++)
            input_libsais[i] = map_int[input[i]];

        input_libsais[input_size] = 0;
        suffix_array[input_size] = 0;
        libsais64_long_omp(&input_libsais[0], &suffix_array[0], input_size + 1, alphabet.size() + 1, 0, num_threads);
    }

    return suffix_array;
}

/**
 * @brief verifies the occurrences reported by a count- and a locate-query against the naive ground-truth
 * @tparam pos_t integer type of the reported occurrences
 * @param count_num_occ the number of occurrences reported by count() (end - beg + 1, or 0 if beg > end)
 * @param located the occurrences reported by locate() (reordered in place)
 * @param correct_occurrences the ascending starting positions of all exact occurrences of the pattern
 * @param label a short label identifying the index/query for failure messages
 */
template <typename pos_t>
static void verify_occurrences(
    uint64_t count_num_occ, std::vector<pos_t>& located,
    const std::vector<pos_t>& correct_occurrences, const std::string& label)
{
    EXPECT_EQ(count_num_occ, correct_occurrences.size()) << label << ": count";
    ips4o::sort(located.begin(), located.end());
    EXPECT_EQ(located.size(), correct_occurrences.size()) << label << ": locate size";
    EXPECT_EQ(located, correct_occurrences) << label << ": locate values";
}

/**
 * @brief draws a pattern from input, optionally mutates it, and returns it together with the ascending
 * starting positions of all of its exact occurrences in input
 * @tparam pos_t integer type used for positions
 * @tparam inp_t sequence type of the input and pattern (e.g. std::string or std::vector<int32_t>)
 * @tparam gen_t random number generator type
 * @param input the input sequence to draw the pattern from and search in
 * @param max_pattern_length the maximum length of the drawn (pre-mutation) pattern
 * @param mutate whether to apply random edits to the drawn pattern
 * @param gen random number generator
 * @return a pair of the pattern and the ascending positions of its exact occurrences in input
 */
template <typename pos_t, typename inp_t, typename gen_t>
static std::pair<inp_t, std::vector<pos_t>> random_pattern_and_occurrences(
    const inp_t& input, uint64_t max_pattern_length, bool mutate, gen_t& gen)
{
    uint64_t input_size = input.size();
    std::uniform_int_distribution<uint64_t> pattern_pos_distrib(0, input_size - 1);
    std::uniform_int_distribution<uint64_t> pattern_length_distrib(
        1, std::min<uint64_t>(max_pattern_length, input_size));

    uint64_t pattern_pos = pattern_pos_distrib(gen);
    uint64_t pattern_length = std::min<uint64_t>(input_size - pattern_pos, pattern_length_distrib(gen));
    inp_t pattern(input.begin() + pattern_pos, input.begin() + pattern_pos + pattern_length);

    if (mutate) mutate_pattern<uint64_t>(pattern, input, gen);
    if (pattern.empty()) pattern.push_back(input[pattern_pos]);

    std::vector<pos_t> correct_occurrences = locate_naive<pos_t>(input, pattern);
    return { std::move(pattern), std::move(correct_occurrences) };
}

/**
 * @brief verifies a single count- and locate-query of an index against the naive occurrences (the part of
 * the verification that is common to every index); used by the adapters' verify_pattern()
 * @tparam adapter_t the index adapter type (see test_index_instance)
 * @param index the index to query
 * @param pattern the pattern to query for
 * @param correct_occurrences the ascending starting positions of all exact occurrences of pattern
 */
template <typename adapter_t>
static void verify_count_locate_basic(
    typename adapter_t::index_t& index, typename adapter_t::inp_t& pattern,
    const std::vector<typename adapter_t::pos_t>& correct_occurrences)
{
    using pos_t = typename adapter_t::pos_t;
    std::vector<pos_t> located = adapter_t::locate(index, pattern);
    verify_occurrences<pos_t>(adapter_t::count(index, pattern), located, correct_occurrences, adapter_t::name());
}

/**
 * @brief the generic index fuzz-test: builds a random index via the adapter, runs the adapter's build-time
 * checks (verify_build) and a parallel fuzz of count-/locate-queries (verify_pattern), and, if the index is
 * serializable, repeats the checks after a serialize/load round-trip
 *
 * an adapter is a type providing:
 *   - types: index_t, inp_t (input/pattern sequence), pos_t (occurrence type)
 *   - static constexpr bool serializable
 *   - static std::string name()
 *   - static <sym> min_sym(); static <sym> max_sym()           // alphabet bounds for the random input
 *   - static uint64_t max_pattern_length(uint64_t n)
 *   - static uint64_t num_queries(uint64_t n)
 *   - static index_t build(inp_t& input, gen_t&, uint16_t max_num_threads)
 *   - static void after_build(index_t&, const inp_t& input)     // e.g. set_input (or a no-op)
 *   - static uint64_t count(index_t&, inp_t& pattern)
 *   - static std::vector<pos_t> locate(index_t&, inp_t& pattern)
 *   - static void verify_pattern(index_t&, inp_t& pattern, const std::vector<pos_t>& correct, gen_t&)
 *   - static void verify_build(index_t&, const inp_t& input, gen_t&, uint16_t max_num_threads)
 *
 * @tparam adapter_t the index adapter type
 * @tparam gen_t random number generator type
 * @param gen random number generator
 * @param max_num_threads the maximum number of threads to use
 * @param min_input_size minimum input length
 * @param max_input_size maximum input length
 */
template <typename adapter_t, typename gen_t>
static void test_index_instance(
    gen_t& gen, uint16_t max_num_threads, uint64_t min_input_size, uint64_t max_input_size)
{
    using inp_t = typename adapter_t::inp_t;
    using pos_t = typename adapter_t::pos_t;
    using index_t = typename adapter_t::index_t;

    inp_t input = random_repetitive_input<inp_t>(
        min_input_size, max_input_size, adapter_t::min_sym(), adapter_t::max_sym());
    uint64_t n = input.size();

    index_t index = adapter_t::build(input, gen, max_num_threads);
    adapter_t::after_build(index, input);

    uint64_t max_pattern_length = adapter_t::max_pattern_length(n);
    uint64_t total_queries = adapter_t::num_queries(n);

    auto run_queries = [&](index_t& idx) {
        uint64_t queries_per_thread = std::max<uint64_t>(1, total_queries / max_num_threads);

        #pragma omp parallel num_threads(max_num_threads)
        {
            std::random_device rd_thr;
            std::mt19937 gen_thr(rd_thr());

            for (uint64_t q = 0; q < queries_per_thread; q++) {
                auto [pattern, correct_occurrences] =
                    random_pattern_and_occurrences<pos_t, inp_t>(input, max_pattern_length, true, gen_thr);
                adapter_t::verify_pattern(idx, pattern, correct_occurrences, gen_thr);
            }
        }
    };

    adapter_t::verify_build(index, input, gen, max_num_threads);
    run_queries(index);

    if constexpr (adapter_t::serializable) {
        std::stringstream stream;
        index.serialize(stream);
        index_t index_reloaded;
        index_reloaded.load(stream);
        adapter_t::after_build(index_reloaded, input);
        adapter_t::verify_build(index_reloaded, input, gen, max_num_threads);
        run_queries(index_reloaded);
    }
}
