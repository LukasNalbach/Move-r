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

// the checks of an index built by the fuzz tests of move_r, shared by tests/test-move-r.cpp and
// tests/test-move-r-instantiations.cpp

#ifndef MOVE_R_TESTS_TEST_MOVE_R_COMMON_HPP
#define MOVE_R_TESTS_TEST_MOVE_R_COMMON_HPP

#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>

#include <ips2ra.hpp>
#include <gtest/gtest.h>
#include <misc/build_cli.hpp>
#include <move_r/move_r.hpp>
#include "test-references.hpp"
#include "test-patterns.hpp"
#include "test-fuzz.hpp"
#include "test-inputs.hpp"
#include "test-parameters.hpp"
#include "test-revert.hpp"

inline std::mt19937 gen(test_seed());
inline uint64_t max_num_threads = omp_get_max_threads();
inline uint64_t min_input_size = 1;
inline uint64_t max_input_size = 16000;

/* number of largest values of the symbol type that do not occur in the random inputs, so that a byte input has at most
   max_byte_alphabet_size distinct symbols (as many as the construction with _pfp reserves) */
inline constexpr uint64_t num_excluded_symbol_values = pfp_min_mapped_symbol;

/* TEST(test_move_r, index) checks the same support and input type in iterations_per_combination consecutive iterations;
   the input type shifts by one after each round over the num_supports supports, so that every combination of a
   support and an input type occurs in iterations_per_combination * num_supports * num_input_types consecutive
   iterations */
inline constexpr uint64_t iterations_per_combination = 2;
inline constexpr uint64_t num_supports = 8;
inline constexpr uint64_t num_input_types = 5;

// max_pattern_length(n) of move_r_adapter is n / input_length_per_pattern_length, clamped to
// [min_max_pattern_length..max_max_pattern_length]
inline constexpr uint64_t input_length_per_pattern_length = 1000;
inline constexpr uint64_t min_max_pattern_length = 100;
inline constexpr uint64_t max_max_pattern_length = 10000;
// num_queries(n) of move_r_adapter is n / input_length_per_query, clamped to [min_num_queries..max_num_queries]
inline constexpr uint64_t input_length_per_query = 4;
inline constexpr uint64_t min_num_queries = 12000;
inline constexpr uint64_t max_num_queries = 25000;

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
    sort_sequential(located.begin(), located.end(), std::less<>());
    EXPECT_EQ(located.size(), correct_occurrences.size()) << label << ": locate size";
    EXPECT_EQ(located, correct_occurrences) << label << ": locate values";
}

/**
 * @brief builds one index of a random input with the adapter and checks its queries with random patterns, also after
 *        serializing and loading the index
 *
 * An adapter is a type with:
 *   - types index_t (the index) and inp_t (the container of the input and of the patterns)
 *   - static <sym> min_sym(), static <sym> max_sym(): the smallest and largest symbol of the random input
 *   - static uint64_t max_pattern_length(uint64_t n), static uint64_t num_queries(uint64_t n)
 *   - static index_t build(inp_t& input, gen_t&, uint64_t max_num_threads)
 *   - static uint64_t count(index_t&, inp_t& pattern), static <occurrences> locate(index_t&, inp_t& pattern)
 *   - static void verify_pattern(index_t&, inp_t& pattern, const std::vector<uint64_t>& occurrences, gen_t&)
 *   - static void verify_retrieval(index_t&, const inp_t& input, gen_t&, uint64_t max_num_threads)
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

    // seed of the input, drawn from gen (the input is reproduced with MOVE_R_TEST_SEED)
    uint64_t input_seed = std::uniform_int_distribution<uint64_t>(1, ~uint64_t(0))(gen);
    inp_t input = random_repetitive_input<inp_t>(
        min_input_size, max_input_size, adapter_t::min_sym(), adapter_t::max_sym(), input_seed);
    uint64_t n = input.size();

    index_t index = adapter_t::build(input, gen, max_num_threads);

    uint64_t max_pattern_length = adapter_t::max_pattern_length(n);
    uint64_t num_queries = adapter_t::num_queries(n);

    // seed of the per-thread generators, drawn from gen (the same in both passes of verify_everything)
    uint64_t query_seed = gen();

    // checks the retrieval queries and num_queries random patterns on idx
    auto verify_everything = [&](index_t& idx) {
        adapter_t::verify_retrieval(idx, input, gen, max_num_threads);
        uint64_t queries_per_thread = std::max<uint64_t>(1, num_queries / max_num_threads);

        #pragma omp parallel num_threads(max_num_threads)
        {
            std::mt19937_64 gen_thr(query_seed + omp_get_thread_num());

            for (uint64_t q = 0; q < queries_per_thread; q++) {
                auto [pattern, correct_occurrences] =
                    random_pattern_and_occurrences<inp_t>(input, max_pattern_length, true, gen_thr);
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
 * @tparam mlf_enc the encoding of M_LF
 */
template <move_r_support support, typename sym_t_, typename inp_t_, move_pos_encoding_t mlf_enc = EF>
struct move_r_adapter {
    // the index type does not depend on the symbol type of the input
    using index_t = move_r<support, mlf_enc>;
    using sym_t = sym_t_;
    using inp_t = inp_t_;
    // container type of the symbols reported by the index (unsigned symbols of the width of sym_t)
    using out_t = std::vector<std::make_unsigned_t<sym_t_>>;

    // parameters of the last built index, appended to the failure messages
    inline static std::string label;

    static sym_t min_sym() { return std::numeric_limits<sym_t>::min(); }
    static sym_t max_sym() { return std::numeric_limits<sym_t>::max() - num_excluded_symbol_values; }

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
     * @param max_num_threads the maximum number of threads of the construction
     * @return the index
     */
    template <typename gen_t>
    static index_t build(inp_t& input, gen_t& gen, uint64_t max_num_threads)
    {
        std::lognormal_distribution<double> a_distrib(balancing_parameter_log_mean, balancing_parameter_log_deviation);
        std::uniform_int_distribution<int32_t> share_distrib(0, int32_t(context_table_shares.size()) - 1);

        move_r_params params {
            .mode = _suffix_array,
            .num_threads = std::uniform_int_distribution<uint64_t>(1, max_num_threads)(gen),
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

        // half of the constructions read the input from a file; the empty container symbol_type then passes the
        // symbol type to the constructor
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
    static std::vector<uint48_t> locate(index_t& index, inp_t& pattern)
        requires(index_t::supports_multiple_locate) { return index.locate(pattern); }

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

        // if a prepend() failed, the context belongs to the longest suffix of the pattern that occurs
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

            ips2ra::sort(occurrences.begin(), occurrences.end(), key_of_value {});
            EXPECT_EQ(occurrences, correct_occurrences);
        }
    }

    /**
     * @brief checks revert(), SA_range(), SA(i), BWT_range() and BWT(i) against the naive suffix array
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

        if constexpr (index_t::supports_multiple_locate) {
            /* the input mapped to the internal symbols of the index (the index orders the suffixes by the internal
               symbols): [0..input_size-1] mapped[i] = map_symbol(input[i]) */
            inp_t mapped;
            resize_for_overwrite(mapped, input_size);

            #pragma omp parallel for num_threads(max_num_threads)
            for (uint64_t i = 0; i < input_size; i++) mapped[i] = sym_t(index.map_symbol(input[i]));

            // [0..input_size] the suffix array of mapped with the sentinel
            std::vector<int64_t> suffix_array = compute_reference_sa<sym_t, inp_t>(mapped, input_size, max_num_threads);
            const no_init_vector<uint48_t> suffix_array_retrieved = index.SA_range({ .num_threads = max_num_threads });

            #pragma omp parallel for num_threads(max_num_threads)
            for (uint64_t i = 0; i <= input_size; i++)
                EXPECT_EQ(suffix_array[i], suffix_array_retrieved[i]) << "SA[" << i << "], " << label;

            if constexpr (index_t::has_SA_access) {
                #pragma omp parallel for num_threads(max_num_threads)
                for (uint64_t i = 0; i <= input_size; i++)
                    EXPECT_EQ(index.SA(i), suffix_array[i]) << "SA(" << i << "), " << label;
            }

            // the BWT of the input: [0..input_size] bwt[i] = input[SA[i] - 1], 0 if SA[i] = 0
            inp_t bwt;
            resize_for_overwrite(bwt, input_size + 1);

            #pragma omp parallel for num_threads(max_num_threads)
            for (uint64_t i = 0; i <= input_size; i++)
                bwt[i] = suffix_array[i] == 0 ? sym_t(0) : input[suffix_array[i] - 1];

            auto bwt_retrieved = index.template BWT_range<out_t>({ .num_threads = max_num_threads });

            #pragma omp parallel for num_threads(max_num_threads)
            for (uint64_t i = 0; i <= input_size; i++)
                EXPECT_EQ(index_t::symbol_idx(bwt[i]), uint64_t(bwt_retrieved[i]));

            // BWT(i) reports the terminator as TERMINATOR_SYMBOL
            #pragma omp parallel for num_threads(max_num_threads)
            for (uint64_t i = 0; i <= input_size; i++)
                EXPECT_EQ(index.BWT(i), suffix_array[i] == 0 ? TERMINATOR_SYMBOL : index_t::symbol_idx(bwt[i]));
        }
    }
};

/**
 * @brief builds an index of a random input and checks its queries (see verify_index), with the support and the encoding
 *        of M_LF of iteration i
 * @tparam sym_t the symbol type of the input
 * @tparam inp_t the container of the input
 * @param i the iteration
 * @param min_n minimum input length
 * @param max_n maximum input length
 */
template <typename sym_t, typename inp_t>
void verify_support_of_iteration(uint64_t i, uint64_t min_n, uint64_t max_n)
{
    auto run = [&]<move_r_support support, move_pos_encoding_t mlf_enc>() {
        verify_index<move_r_adapter<support, sym_t, inp_t, mlf_enc>>(gen, max_num_threads, min_n, max_n);
    };

    switch ((i / iterations_per_combination) % num_supports) {
        case 0: run.template operator()<_locate_move,  EF>();   break;
        case 1: run.template operator()<_locate_rlzsa, EF>();   break;
        case 2: run.template operator()<_locate_move,  POS>();  break;
        case 3: run.template operator()<_locate_rlzsa, POS>();  break;
        case 4: run.template operator()<_count,        EF>();   break;
        case 5: run.template operator()<_locate_one,   EF>();   break;
        case 6: run.template operator()<_count,        POS>();  break;
        case 7: run.template operator()<_locate_one,   POS>();  break;
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
