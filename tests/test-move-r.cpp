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

#include "test-move-r-common.hpp"

// smallest value run_fuzz reduces max_input_size to in TEST(test_move_r, index)
constexpr uint64_t min_reduced_max_input_size = 512;

// TEST(test_move_r, parallel_construction): inputs that give every thread twice the minimum number of symbols per
// thread of the construction (see construction::limit_num_threads)
constexpr uint64_t input_length_per_thread = 2 * move_r_min_symbols_per_thread;

// TEST(test_move_r, retrieve_range_to_file): the retrieved range [range_beg, range_end]
constexpr uint64_t range_beg = 10;
constexpr uint64_t range_end = 20;

// the number of iterations of the fuzz tests other than TEST(test_move_r, index) (see fuzz_iterations)
constexpr uint64_t num_iterations_parallel_construction = 20;
constexpr uint64_t num_iterations_pfp_integer_alphabet = 30;
constexpr uint64_t num_iterations_input_widths = 10;
constexpr uint64_t num_iterations_retrieve_range_to_file = 20;

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
 * @brief builds an index of a random input and checks its queries (see verify_index), with the input type and the
 *        support of iteration i; the input types are one-byte symbols (byte alphabet) and four-byte symbols (integer
 *        alphabet), each signed and unsigned
 * @param i the iteration
 * @param min_n minimum input length
 * @param max_n maximum input length
 */
static void verify_index_of_iteration(uint64_t i, uint64_t min_n = min_input_size, uint64_t max_n = max_input_size)
{
    const uint64_t combination = i / iterations_per_combination;

    switch ((combination + combination / num_supports) % num_input_types) {
        case 0: verify_support_of_iteration<char,     std::string>          (i, min_n, max_n); break;
        case 1: verify_support_of_iteration<uint8_t,  std::vector<uint8_t>> (i, min_n, max_n); break;
        case 2: verify_support_of_iteration<int8_t,   std::vector<int8_t>>  (i, min_n, max_n); break;
        case 3: verify_support_of_iteration<uint32_t, std::vector<uint32_t>>(i, min_n, max_n); break;
        case 4: verify_support_of_iteration<int32_t,  std::vector<int32_t>> (i, min_n, max_n); break;
    }
}

// builds one index per iteration and checks every query on it
TEST(test_move_r, index)
{
    run_fuzz("move-r", { { "index", [](uint64_t i) { verify_index_of_iteration(i); } } },
        fuzz_iterations(), &max_input_size, min_reduced_max_input_size);
}

// inputs of length input_length_per_thread * max_num_threads, which the construction builds with all threads
TEST(test_move_r, parallel_construction)
{
    run_fuzz("move-r", { { "parallel-construction", [](uint64_t i) {
        uint64_t n = input_length_per_thread * max_num_threads;
        verify_index_of_iteration(i, n, n);
    } } }, fuzz_iterations(num_iterations_parallel_construction), nullptr, 0);
}

/* the index of a text with an integer alphabet is the same with the prefix-free parsing (of the text in memory and of
   a file) and with the suffix array */
TEST(test_move_r, pfp_integer_alphabet)
{
    run_fuzz("move-r", { { "pfp-integer-alphabet", [](uint64_t i) {
        auto run = [&]<typename sym_t>() {
            /* the symbols are drawn from span consecutive values (up to the maximum input length, so that the symbols
               of the parsing are often wider than a byte), starting at 1 or close to the largest value of sym_t */
            const uint64_t span = random_log_uniform_size(min_random_alphabet_size, max_input_size, gen);
            const uint64_t first = gen() % 2 == 0 ? 1
                : uint64_t(std::numeric_limits<sym_t>::max()) - num_excluded_symbol_values - span;
            std::vector<sym_t> input = random_repetitive_input<std::vector<sym_t>>(min_pfp_input_size,
                max_input_size, sym_t(first), sym_t(first + span - 1), gen() | 1);

            move_r_params params { .mode = _suffix_array, .num_threads = 1 + gen() % max_num_threads,
                .reorder_alphabet = gen() % 2 == 0 };
            const move_r<_locate_move> reference(input, params);
            std::stringstream expected;
            reference.serialize(expected);
            params.mode = _pfp;

            for (bool from_file : { false, true }) {
                const std::string label = std::to_string(8 * sizeof(sym_t)) + "-bit symbols, n = "
                    + std::to_string(input.size()) + ", sigma = " + std::to_string(reference.alphabet_size())
                    + (from_file ? ", from a file" : ", from memory") + ", "
                    + std::to_string(params.num_threads) + " threads"
                    + (params.reorder_alphabet ? ", reordered" : ", not reordered");
                // the empty container symbol_type passes the symbol type to the constructor
                std::vector<sym_t> symbol_type;
                params.input_file = from_file ? write_temporary_input(input) : "";
                const move_r<_locate_move> index(from_file ? symbol_type : input, params);
                std::stringstream actual;
                index.serialize(actual);
                if (from_file) remove_input_file_and_pfp_cache(params.input_file);

                ASSERT_TRUE(actual.str() == expected.str()) << "the index differs from the one built with the "
                    << "suffix array, " << label;
            }
        };

        switch (i % 3) {
            case 0: run.template operator()<uint16_t>(); break;
            case 1: run.template operator()<uint32_t>(); break;
            default: run.template operator()<uint64_t>(); break;
        }
    } } }, fuzz_iterations(num_iterations_pfp_integer_alphabet), nullptr, 0);
}

/* the index of a text with symbols of 2 to 7 bytes (built with the suffix array or the prefix-free parsing) is the same
   as the one of the same values in 64 bits (the symbols of 1 byte are a byte alphabet) */
TEST(test_move_r, input_widths)
{
    run_fuzz("move-r", { { "input-widths", [](uint64_t i) {
        auto run = [&]<typename sym_t>() {
            if constexpr (sizeof(sym_t) > sizeof(uint8_t) && sizeof(sym_t) < sizeof(uint64_t)) {
                // the symbols are drawn from span consecutive values, starting at 1 or close to the largest value
                const uint64_t span = random_log_uniform_size(min_random_alphabet_size, max_input_size, gen);
                const uint64_t first = gen() % 2 == 0 ? 1 : max_value_of<sym_t>() - num_excluded_symbol_values - span;
                std::vector<uint64_t> values = random_repetitive_input<std::vector<uint64_t>>(min_pfp_input_size,
                    max_input_size, first, first + span - 1, gen() | 1);
                std::vector<sym_t> input(values.begin(), values.end());

                move_r_params params { .mode = _suffix_array, .num_threads = 1 + gen() % max_num_threads,
                    .reorder_alphabet = gen() % 2 == 0 };
                const move_r<_locate_move> reference(values, params);
                std::stringstream expected;
                reference.serialize(expected);
                params.mode = i % 2 == 0 ? _suffix_array : _pfp;
                const move_r<_locate_move> index(input, params);
                std::stringstream actual;
                index.serialize(actual);

                ASSERT_TRUE(actual.str() == expected.str()) << "the index differs from the one of 64-bit symbols, "
                    << sizeof(sym_t) << "-byte symbols, n = " << input.size() << ", sigma = "
                    << reference.alphabet_size() << (params.mode == _pfp ? ", prefix-free parsing" : ", suffix array");
            }
        };

        #define MOVE_R_TEST_INPUT_WIDTH(sym_t) run.template operator()<sym_t>();
        MOVE_R_FOR_EACH_INPUT_WIDTH(MOVE_R_TEST_INPUT_WIDTH)
        #undef MOVE_R_TEST_INPUT_WIDTH
    } } }, fuzz_iterations(num_iterations_input_widths), nullptr, 0);
}

// BWT_range and SA_range with a file name write the values of the range [l, r] to the file
TEST(test_move_r, retrieve_range_to_file)
{
    run_fuzz("move-r", { { "retrieve-range-to-file", [](uint64_t) {
        using index_t = move_r<_locate_move>;

        // n - 1 = input size >= range_end
        std::string input = random_repetitive_input<std::string>(range_end, max_input_size, 'a', 'd', gen() | 1);
        const index_t index(input, move_r_params { .mode = _suffix_array, .num_threads = 1 + gen() % max_num_threads });
        const index_t::retrieve_params range { .l = range_beg, .r = range_end, .num_threads = 1 + gen() % max_num_threads };
        const std::string label = "n = " + std::to_string(input.size() + 1) + ", "
            + std::to_string(range.num_threads) + " threads";
        const std::string file = temporary_file_path("move-r-test-");

        index.BWT_range(file, range);
        const std::vector<uint8_t> bwt_from_file = read_values<uint8_t>(file);
        EXPECT_EQ(bwt_from_file.size(), range_end - range_beg + 1) << "BWT, " << label;
        EXPECT_EQ(bwt_from_file, index.BWT_range<std::vector<uint8_t>>(range)) << "BWT, " << label;

        index.SA_range(file, range);
        const std::vector<uint64_t> sa_from_file = read_values<uint64_t>(file);
        const no_init_vector<uint48_t> sa = index.SA_range(range);
        EXPECT_EQ(sa_from_file.size(), range_end - range_beg + 1) << "SA, " << label;
        EXPECT_EQ(sa_from_file, std::vector<uint64_t>(sa.begin(), sa.end())) << "SA, " << label;

        std::filesystem::remove(file);
    } } }, fuzz_iterations(num_iterations_retrieve_range_to_file), nullptr, 0);
}

// short texts (fewer runs than revert starts requested, a single run of a symbol), with a byte and an integer
// alphabet, reverted with every number of threads of revert_thread_counts
TEST(test_move_r, revert_short_texts)
{
    const std::vector<std::string> texts { "a", "aa", "ab", "aaaaaaaa", "abababab", "abracadabra" };

    auto check = [&]<move_r_support support, typename inp_t>(inp_t text, uint64_t starts) {
        using usym_t = std::make_unsigned_t<typename inp_t::value_type>;
        const move_r<support> index(text, move_r_params { .mode = _suffix_array, .revert_start_samples = starts });
        const std::string label = "support " + std::to_string(int(support)) + ", " + std::to_string(8 * sizeof(usym_t))
            + "-bit symbols, n = " + std::to_string(text.size()) + ", revert starts " + std::to_string(starts);
        verify_parallel_revert<usym_t>(text.size(), [&](uint64_t i) { return uint64_t(usym_t(text[i])); },
            [&](uint64_t p, uint64_t k) { return index.template revert<std::vector<usym_t>>(p, k); },
            [&](const std::string& file, uint64_t p, uint64_t k) { index.template revert<usym_t>(file, p, k); }, label);
    };

    for (const std::string& text : texts) {
        const std::vector<uint32_t> integers(text.begin(), text.end());

        for (uint64_t starts : revert_start_sample_counts) {
            check.template operator()<_count>(text, starts);
            check.template operator()<_locate_move>(text, starts);
            check.template operator()<_count>(integers, starts);
            check.template operator()<_locate_rlzsa>(integers, starts);
        }
    }
}

// a failed pfp construction is repeated with sa-space only without -c and if the input is short enough for sa-space
TEST(test_move_r, build_fallback_limit)
{
    move_r_construction_mode last_mode = _pfp; // the construction mode of the last call of build

    auto build = [&](move_r_construction_mode mode) {
        last_mode = mode;
        if (mode == _pfp) throw std::runtime_error("the pfp construction failed");
        return mode;
    };

    EXPECT_EQ(build_with_fallback(_pfp, false, move_r_max_input_length_sa, build), _suffix_array_space);
    EXPECT_THROW(build_with_fallback(_pfp, false, move_r_max_input_length_sa + 1, build), std::runtime_error);
    EXPECT_EQ(last_mode, _pfp);
    EXPECT_THROW(build_with_fallback(_pfp, true, 0, build), std::runtime_error);
    EXPECT_EQ(last_mode, _pfp);
}
