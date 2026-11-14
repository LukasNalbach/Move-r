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

/**
 * @file test-move-r-common.hpp
 * @brief the checks of an index built by the fuzz tests of move_r, shared by tests/test-move-r.cpp and
 *        tests/test-move-r-instantiations.cpp
 */

#ifndef MOVE_R_TESTS_TEST_MOVE_R_COMMON_HPP
#define MOVE_R_TESTS_TEST_MOVE_R_COMMON_HPP

#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>

#include <gtest/gtest.h>
#include <misc/build_cli.hpp>
#include <move_r/move_r.hpp>
#include "test-references.hpp"
#include "test-patterns.hpp"
#include "test-fuzz.hpp"
#include "test-inputs.hpp"
#include "test-parameters.hpp"
#include "test-revert.hpp"

// one generator per thread (run_fuzz runs parallel functionalities on all threads)
inline thread_local std::mt19937 gen(test_seed() + omp_get_thread_num());
inline uint64_t max_num_threads = omp_get_max_threads();
inline uint64_t min_input_size = 1;
// the largest length of a random input; run_fuzz reduces max_input_size, the inputs of full size are not reduced (see
// verify_index)
inline constexpr uint64_t full_input_size = 16000;
inline uint64_t max_input_size = full_input_size;

/* verify_index draws an input of input_length_per_thread symbols per thread of the construction (such that the
   construction uses all of its threads, see construction::limit_num_threads) with probability 1 /
   parallel_input_inverse_probability, and an input of a length between max_input_size and full_input_size with
   probability 1 / full_size_inverse_probability */
inline constexpr uint64_t input_length_per_thread = 2 * move_r_min_symbols_per_thread;
inline constexpr uint64_t parallel_input_inverse_probability = 8;
inline constexpr uint64_t full_size_inverse_probability = 8;

/* TEST(test_move_r, index) checks the same support and input type in iterations_per_combination consecutive iterations;
   the input type shifts by one after each round over the num_supports supports, so that every combination of a
   support and an input type occurs in iterations_per_combination * num_supports * num_input_types consecutive
   iterations */
inline constexpr uint64_t iterations_per_combination = 2;
inline constexpr uint64_t num_supports = 4;
inline constexpr uint64_t num_input_types = 5;

// max_pattern_length(n) of move_r_adapter is n / input_length_per_pattern_length, clamped to
// [min_max_pattern_length..max_max_pattern_length]
inline constexpr uint64_t input_length_per_pattern_length = 1000;
inline constexpr uint64_t min_max_pattern_length = 100;
inline constexpr uint64_t max_max_pattern_length = 10000;
// num_queries(n) of move_r_adapter is n / input_length_per_query, clamped to [min_num_queries..max_num_queries]
inline constexpr uint64_t input_length_per_query = 4;
inline constexpr uint64_t min_num_queries = 16;
inline constexpr uint64_t max_num_queries = 25000;

/**
 * @brief returns a random input (see random_text): over a random byte alphabet (see random_input_alphabet) for symbols
 *        of one byte, else over a random integer alphabet (see random_integer_alphabet)
 * @tparam inp_t the container of the input
 * @param min_size the smallest length
 * @param max_size the largest length
 * @param gen random number generator
 */
template <typename inp_t, typename gen_t>
inline inp_t random_input(uint64_t min_size, uint64_t max_size, gen_t& gen)
{
    using sym_t = typename inp_t::value_type;

    if constexpr (sizeof(sym_t) == 1) {
        return random_input_alphabet(gen).template input<inp_t>(min_size, max_size, gen);
    } else {
        const integer_alphabet<sym_t> alphabet = random_integer_alphabet<sym_t>(1, max_size, gen);
        return random_text<inp_t>(min_size, max_size, alphabet.sigma, alphabet, gen);
    }
}

/**
 * @brief returns the contents of a file as an array of values of type value_t
 * @param file the path of the file
 * @return [0..file size / sizeof(value_t) - 1] the values
 */
template <typename value_t>
std::vector<value_t> read_values(const std::string& file)
{
    std::ifstream in(file, std::ios::binary);
    std::vector<value_t> values(std::filesystem::file_size(file) / sizeof(value_t));
    in.read(reinterpret_cast<char*>(values.data()), values.size() * sizeof(value_t));
    return values;
}

/**
 * @brief writes an input to a new temporary file
 * @param input the input
 * @return the path of the file
 */
template <typename inp_t>
std::string write_temporary_input(const inp_t& input)
{
    std::string input_file = temporary_file_path("move-r-test-");
    direct_ofstream out(input_file);
    out.write(reinterpret_cast<const char*>(input.data()), input.size() * sizeof(typename inp_t::value_type));
    out.close();
    return input_file;
}

/**
 * @brief removes a temporary input file and the pfp cache the construction with _pfp writes next to it (see
 *        move_r::pfp_cache_path)
 * @param input_file the path of the input file
 */
inline void remove_input_file_and_pfp_cache(const std::string& input_file)
{
    std::filesystem::remove(input_file);
    std::filesystem::remove(move_r<>::pfp_cache_path(input_file, _pfp));
}

/**
 * @brief checks count() and locate() of an index for a pattern
 * @tparam adapter_t the index adapter type (see verify_index)
 * @param index the index
 * @param pattern the pattern
 * @param correct_occurrences the starting positions of the occurrences of pattern in the input, in increasing order
 */
template <typename adapter_t>
inline void verify_count_and_locate(
    typename adapter_t::index_t& index, typename adapter_t::inp_t& pattern,
    const std::vector<uint64_t>& correct_occurrences)
{
    const std::string& label = adapter_t::label;
    EXPECT_EQ(adapter_t::count(index, pattern), correct_occurrences.size()) << label << ": count";

    auto reported = adapter_t::locate(index, pattern);
    std::vector<uint64_t> located(reported.begin(), reported.end());
    std::sort(located.begin(), located.end());
    EXPECT_EQ(located.size(), correct_occurrences.size()) << label << ": locate size";
    EXPECT_EQ(located, correct_occurrences) << label << ": locate values";
}

/**
 * @brief builds one index of a random input (see random_input) with the adapter and checks its queries with random
 *        patterns (see random_pattern_and_occurrences), also after serializing and loading the index
 *
 * An adapter is a type with:
 *   - types index_t (the index) and inp_t (the container of the input and of the patterns)
 *   - static uint64_t max_pattern_length(uint64_t n), static uint64_t num_queries(uint64_t n)
 *   - static index_t build(inp_t& input, gen_t&, uint64_t num_threads)
 *   - static uint64_t count(index_t&, inp_t& pattern), static <occurrences> locate(index_t&, inp_t& pattern)
 *   - static void verify_pattern(index_t&, inp_t& pattern, const std::vector<uint64_t>& occurrences, gen_t&)
 *   - static void verify_retrieval(index_t&, const inp_t& input, gen_t&, uint64_t max_num_threads)
 *
 * The number of threads of the construction is 1 in a parallel functionality (run_fuzz runs it on all threads), else
 * drawn from [1..max_num_threads]; the input has input_length_per_thread symbols per thread of the construction with
 * probability 1 / parallel_input_inverse_probability, a length in [max_input_size, full_input_size] with probability
 * 1 / full_size_inverse_probability, else a length in [min_input_size, max_input_size].
 *
 * @tparam adapter_t the index adapter type
 * @tparam gen_t random number generator type
 * @param gen random number generator
 * @param max_num_threads the maximum number of threads to use
 * @param min_input_size minimum input length
 * @param max_input_size maximum input length
 */
template <typename adapter_t, typename gen_t>
inline void verify_index(
    gen_t& gen, uint64_t max_num_threads, uint64_t min_input_size, uint64_t max_input_size)
{
    using inp_t = typename adapter_t::inp_t;
    using index_t = typename adapter_t::index_t;

    const uint64_t num_threads = omp_in_parallel() ? 1 : std::uniform_int_distribution<uint64_t>(1, max_num_threads)(gen);
    uint64_t min_n = min_input_size;
    uint64_t max_n = max_input_size;

    if (num_threads > 1 && gen() % parallel_input_inverse_probability == 0) {
        min_n = max_n = input_length_per_thread * num_threads;
    } else if (gen() % full_size_inverse_probability == 0) {
        min_n = max_input_size;
        max_n = std::max(max_input_size, full_input_size);
    }

    inp_t input = random_input<inp_t>(min_n, max_n, gen);
    uint64_t n = input.size();

    index_t index = adapter_t::build(input, gen, num_threads);

    uint64_t max_pattern_length = adapter_t::max_pattern_length(n);
    uint64_t num_queries = adapter_t::num_queries(n);
    // the symbols the edits of the patterns draw from
    const std::vector<typename inp_t::value_type> symbols = edit_symbols(input);

    // seed of the per-thread generators, drawn from gen (the same in both passes of verify_everything)
    uint64_t query_seed = gen();

    // checks the retrieval queries and num_queries random patterns on idx
    auto verify_everything = [&](index_t& idx) {
        adapter_t::verify_retrieval(idx, input, gen, max_num_threads);
        uint64_t queries_per_thread = std::max<uint64_t>(1, num_queries / max_num_threads);
        const std::string& index_label = adapter_t::label;

        #pragma omp parallel num_threads(max_num_threads)
        {
            // the failure messages of every thread of the team contain the label of the index
            if (&adapter_t::label != &index_label) adapter_t::label = index_label;
            std::mt19937_64 gen_thr(query_seed + omp_get_thread_num());

            for (uint64_t q = 0; q < queries_per_thread; q++) {
                auto [pattern, correct_occurrences] =
                    random_pattern_and_occurrences<inp_t>(input, max_pattern_length, symbols, gen_thr);
                adapter_t::verify_pattern(idx, pattern, correct_occurrences, gen_thr);
            }
        }
    };

    verify_everything(index);

    // check the index again after serializing and loading it
    std::stringstream stream;
    index.serialize(stream);
    index_t index_reloaded;
    index_reloaded.load(stream);
    verify_everything(index_reloaded);
}

/**
 * @brief the adapter of move_r for verify_index
 * @tparam support the locate support of the index
 * @tparam sym_t_ the symbol type of the input
 * @tparam inp_t_ the container of the input
 */
template <move_r_support support, typename sym_t_, typename inp_t_>
struct move_r_adapter {
    // the index type does not depend on the symbol type of the input
    using index_t = move_r<support>;
    using sym_t = sym_t_;
    using inp_t = inp_t_;
    // container type of the symbols reported by the index (unsigned symbols of the width of sym_t)
    using out_t = std::vector<std::make_unsigned_t<sym_t_>>;

    // parameters of the last index built by the thread, appended to the failure messages
    inline static thread_local std::string label;

    static uint64_t max_pattern_length(uint64_t n)
    {
        return std::min<uint64_t>(max_max_pattern_length,
            std::max<uint64_t>(min_max_pattern_length, n / input_length_per_pattern_length));
    }

    static uint64_t num_queries(uint64_t n)
    {
        return std::min<uint64_t>(max_num_queries, std::max<uint64_t>(min_num_queries, n / input_length_per_query));
    }

    /**
     * @brief builds an index of an input with random construction parameters and sets label to them
     * @param input the input
     * @param gen random number generator
     * @param num_threads the number of threads of the construction
     * @return the index
     */
    template <typename gen_t>
    static index_t build(inp_t& input, gen_t& gen, uint64_t num_threads)
    {
        std::lognormal_distribution<double> a_distrib(balancing_parameter_log_mean, balancing_parameter_log_deviation);
        std::uniform_int_distribution<int32_t> share_distrib(0, int32_t(context_table_shares.size()) - 1);

        move_r_params params {
            .mode = _suffix_array,
            .num_threads = num_threads,
            .a = min_balancing_parameter
                + uint64_t(std::min(a_distrib(gen), double(max_balancing_parameter - min_balancing_parameter))),
            .reorder_alphabet = gen() % 2 == 0,
            .revert_start_samples = revert_start_sample_counts[gen() % revert_start_sample_counts.size()],
            .context_table_share = context_table_shares[share_distrib(gen)]
        };

        // the prefix-free parsing for large enough inputs, else one of the two modes with the suffix array
        bool use_pfp = input.size() >= min_pfp_input_size && gen() % pfp_inverse_probability == 0;
        params.mode = use_pfp ? _pfp
            : (gen() % 2 == 0 ? _suffix_array : _suffix_array_space);

        // half of the constructions read the input from a file; the symbol type is then passed to the constructor with
        // the empty container symbol_type
        std::string input_file;
        inp_t symbol_type;

        if (gen() % 2 == 0) {
            input_file = write_temporary_input(input);
            params.input_file = input_file;
        }

        label = std::string("support ") + std::to_string(int(support))
            + ", " + std::to_string(8 * sizeof(sym_t)) + "-bit symbols"
            + ", n = " + std::to_string(input.size())
            + ", a = " + std::to_string(params.a)
            + ", mode " + (params.mode == _pfp ? "pfp"
                : (params.mode == _suffix_array ? "suffix-array" : "suffix-array-space"))
            + (input_file.empty() ? ", from memory" : ", from a file")
            + ", " + std::to_string(params.num_threads) + " threads"
            + (params.reorder_alphabet ? ", reordered" : ", not reordered")
            + ", revert starts " + std::to_string(params.revert_start_samples)
            + ", context share " + std::to_string(params.context_table_share);

        index_t index(input_file.empty() ? input : symbol_type, params);
        if (!input_file.empty()) remove_input_file_and_pfp_cache(input_file);
        return index;
    }

    static uint64_t count(index_t& index, inp_t& pattern) { return index.count(pattern); }
    static std::vector<uint48_t> locate(index_t& index, inp_t& pattern) requires(index_t::supports_multiple_locate) { return index.locate(pattern); }

    /**
     * @brief checks count(), locate() and the query context of index.query() for a pattern; the occurrences of the
     *        context are reported with next_occ() up to a random number of them, the remaining ones with locate()
     * @param index the index
     * @param pattern the pattern
     * @param correct_occurrences the starting positions of the occurrences of pattern in the input, in increasing order
     * @param gen random number generator
     */
    template <typename gen_t>
    static void verify_pattern(
        index_t& index, inp_t& pattern, const std::vector<uint64_t>& correct_occurrences, gen_t& gen)
    {
        if constexpr (index_t::supports_multiple_locate) {
            verify_count_and_locate<move_r_adapter>(index, pattern, correct_occurrences);
        } else {
            EXPECT_EQ(count(index, pattern), correct_occurrences.size()) << label;
        }

        // prepend() returns false and leaves the context unchanged if the extended pattern does not occur
        auto query = index.query();
        bool fully_matched = true;

        for (int32_t i = pattern.size() - 1; i >= 0; i--) {
            if (!query.prepend(pattern[i])) {
                fully_matched = false;
                break;
            }
        }

        // if a prepend() failed, the context is the context of the longest suffix of the pattern that occurs
        if (!fully_matched) {
            EXPECT_TRUE(correct_occurrences.empty());
            return;
        }

        EXPECT_EQ(query.num_occ(), correct_occurrences.size());

        // one_occ() reports one of the occurrences (see _locate_one)
        if constexpr (index_t::supports_locate && !index_t::supports_multiple_locate) {
            EXPECT_NE(std::find(correct_occurrences.begin(), correct_occurrences.end(),
                query.one_occ()), correct_occurrences.end()) << label;
        }

        if constexpr (index_t::supports_multiple_locate) {
            std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);
            std::vector<uint64_t> occurrences;

            for (uint64_t i = 0; i < query.num_occ(); i++) {
                if (prob_distrib(gen) < 1 / (double) query.num_occ()) {
                    // locate() reports the remaining occurrences
                    auto remaining_occurrences = query.locate();
                    occurrences.insert(occurrences.end(),
                        remaining_occurrences.begin(), remaining_occurrences.end());
                    break;
                }

                occurrences.emplace_back(query.next_occ());
            }

            std::sort(occurrences.begin(), occurrences.end());
            EXPECT_EQ(occurrences, correct_occurrences);
        }
    }

    /**
     * @brief checks revert(), BWT(i), SA(i) and, for the whole range and a random range, BWT_range() and SA_range() in
     *        memory and to a file against the naive suffix array
     * @param index the index
     * @param input the input of the index
     * @param gen random number generator
     * @param max_num_threads the maximum number of threads
     */
    template <typename gen_t>
    static void verify_retrieval(index_t& index, const inp_t& input, [[maybe_unused]] gen_t& gen,
        uint64_t max_num_threads)
    {
        uint64_t input_size = input.size();
        using usym_t = std::make_unsigned_t<sym_t>;
        // the index reports unsigned symbols
        auto expected = [&](uint64_t i) { return index_t::symbol_idx(input[i]); };

        verify_parallel_revert<usym_t>(input_size, expected,
            [&](uint64_t p, uint64_t k) { return index.template revert<out_t>(p, k); },
            [&](const std::string& file, uint64_t p, uint64_t k) { index.template revert<usym_t>(file, p, k); }, label);

        // revert(report) reports the symbols from right to left
        out_t reported;
        resize_for_overwrite(reported, input_size);
        uint64_t next = input_size; // the position revert(report) reports next, plus 1

        index.revert([&](uint64_t i, uint64_t c) {
            EXPECT_EQ(i + 1, next) << "revert(report) order, " << label;
            next = i;
            reported[i] = c;
        });

        compare_reverted_text(reported, input_size, expected, "revert(report), " + label);

        /* the input mapped to the internal symbols of the index (in the index, the suffixes are ordered by the
           internal symbols): [0..input_size-1] mapped[i] = map_symbol(input[i]) */
        std::vector<uint64_t> mapped;
        resize_for_overwrite(mapped, input_size);

        #pragma omp parallel for num_threads(max_num_threads)
        for (uint64_t i = 0; i < input_size; i++) mapped[i] = index.map_symbol(input[i]);

        // [0..input_size] the suffix array of mapped with the sentinel
        std::vector<int64_t> suffix_array = compute_reference_sa<uint64_t>(mapped, input_size, max_num_threads);
        // the label of the calling thread (label is per thread)
        const std::string& index_label = label;
        // the BWT of the input as the index reports it: [0..input_size] bwt[i] = input[SA[i] - 1], 0 if SA[i] = 0
        std::vector<uint64_t> bwt(input_size + 1);

        #pragma omp parallel for num_threads(max_num_threads)
        for (uint64_t i = 0; i <= input_size; i++)
            bwt[i] = suffix_array[i] == 0 ? 0 : index_t::symbol_idx(input[suffix_array[i] - 1]);

        // BWT(i) reports the terminator as TERMINATOR_SYMBOL
        #pragma omp parallel for num_threads(max_num_threads)
        for (uint64_t i = 0; i <= input_size; i++)
            EXPECT_EQ(index.BWT(i), suffix_array[i] == 0 ? TERMINATOR_SYMBOL : bwt[i]) << "BWT(" << i << "), " << index_label;

        if constexpr (index_t::has_SA_access) {
            #pragma omp parallel for num_threads(max_num_threads)
            for (uint64_t i = 0; i <= input_size; i++)
                EXPECT_EQ(index.SA(i), suffix_array[i]) << "SA(" << i << "), " << index_label;
        }

        // the whole range (l > r) and a random range [l, r]
        const uint64_t l = std::uniform_int_distribution<uint64_t>(0, input_size)(gen);
        const uint64_t r = std::uniform_int_distribution<uint64_t>(l, input_size)(gen);
        const uint64_t range_threads = std::uniform_int_distribution<uint64_t>(1, max_num_threads)(gen);

        for (const typename index_t::retrieve_params range : { typename index_t::retrieve_params { .num_threads = max_num_threads },
                typename index_t::retrieve_params { .l = l, .r = r, .num_threads = range_threads } }) {
            const uint64_t first = range.l <= range.r ? range.l : 0;
            const uint64_t last = range.l <= range.r ? range.r : input_size;
            const std::string what = "[" + std::to_string(first) + ", " + std::to_string(last) + "], "
                + std::to_string(range.num_threads) + " threads, " + index_label;
            const std::vector<uint64_t> bwt_range(bwt.begin() + first, bwt.begin() + last + 1);
            const auto bwt_retrieved = index.template BWT_range<out_t>(range);
            EXPECT_EQ(std::vector<uint64_t>(bwt_retrieved.begin(), bwt_retrieved.end()), bwt_range) << "BWT_range " << what;
            const std::string file = temporary_file_path("move-r-test-");
            index.template BWT_range<usym_t>(file, range);
            const std::vector<usym_t> bwt_from_file = read_values<usym_t>(file);
            EXPECT_EQ(std::vector<uint64_t>(bwt_from_file.begin(), bwt_from_file.end()), bwt_range)
                << "BWT_range to a file " << what;

            if constexpr (index_t::supports_multiple_locate) {
                const std::vector<uint64_t> sa_range(suffix_array.begin() + first, suffix_array.begin() + last + 1);
                const no_init_vector<uint48_t> sa_retrieved = index.SA_range(range);
                EXPECT_EQ(std::vector<uint64_t>(sa_retrieved.begin(), sa_retrieved.end()), sa_range) << "SA_range " << what;
                index.SA_range(file, range);
                EXPECT_EQ(read_values<uint64_t>(file), sa_range) << "SA_range to a file " << what;
            }

            std::filesystem::remove(file);
        }
    }
};

/**
 * @brief builds an index of a random input and checks its queries (see verify_index), with the support of iteration i
 *        and its default encoding of M_LF (see move_r)
 * @tparam sym_t the symbol type of the input
 * @tparam inp_t the container of the input
 * @param i the iteration
 * @param min_n minimum input length
 * @param max_n maximum input length
 */
template <typename sym_t, typename inp_t>
void verify_support_of_iteration(uint64_t i, uint64_t min_n, uint64_t max_n)
{
    auto run = [&]<move_r_support support>() {
        verify_index<move_r_adapter<support, sym_t, inp_t>>(gen, max_num_threads, min_n, max_n);
    };

    switch ((i / iterations_per_combination) % num_supports) {
        case 0: run.template operator()<_locate_move>(); break;
        case 1: run.template operator()<_locate_rlzsa>(); break;
        case 2: run.template operator()<_count>(); break;
        case 3: run.template operator()<_locate_one>(); break;
    }
}

/* the instantiations of verify_support_of_iteration, declared (prefix extern) or defined (empty prefix);
   tests/test-move-r-instantiations.cpp is compiled once per input type */
#define TEST_MOVE_R_VERIFY_SUPPORT_OF_ITERATION(prefix, sym_t, inp_t) \
    prefix template void verify_support_of_iteration<sym_t, inp_t>(uint64_t, uint64_t, uint64_t);
TEST_MOVE_R_VERIFY_SUPPORT_OF_ITERATION(extern, char, std::string)
TEST_MOVE_R_VERIFY_SUPPORT_OF_ITERATION(extern, uint8_t, std::vector<uint8_t>)
TEST_MOVE_R_VERIFY_SUPPORT_OF_ITERATION(extern, int8_t, std::vector<int8_t>)
TEST_MOVE_R_VERIFY_SUPPORT_OF_ITERATION(extern, uint32_t, std::vector<uint32_t>)
TEST_MOVE_R_VERIFY_SUPPORT_OF_ITERATION(extern, int32_t, std::vector<int32_t>)

#endif // MOVE_R_TESTS_TEST_MOVE_R_COMMON_HPP
