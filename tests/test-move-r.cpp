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
 * @file test-move-r.cpp
 * @brief fuzz tests of move_r: the index, the construction from inputs of every symbol width and the fallback of the
 *        construction mode
 */

#include "test-move-r-common.hpp"

// smallest value run_fuzz reduces max_input_size to in TEST(test_move_r, index)
constexpr uint64_t min_reduced_max_input_size = 512;

// the number of iterations of TEST(test_move_r, construction) and TEST(test_move_r, build_fallback) (see fuzz_iterations)
constexpr uint64_t num_iterations_construction = 40;
constexpr uint64_t num_iterations_build_fallback = 1000;

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

// builds one index per iteration and checks every query on it (serial: the checks are run on all threads, see verify_index)
TEST(test_move_r, index)
{
    run_fuzz("move-r", { { "index", [](uint64_t i) { verify_index_of_iteration(i); } } },
        fuzz_iterations(), &max_input_size, min_reduced_max_input_size);
}

/* the index of a random text (see random_text) over a random integer alphabet (see random_integer_alphabet) with
   symbols of 2 to 8 bytes (see MOVE_R_FOR_EACH_INPUT_WIDTH; the symbols of 1 byte are a byte alphabet), built with a
   random construction mode from memory or from a file, is the same as the one of the same values as 64-bit symbols
   built with the suffix array */
TEST(test_move_r, construction)
{
    run_fuzz("move-r", { { "construction", [](uint64_t) {
        auto run = [&]<typename sym_t>() {
            if constexpr (sizeof(sym_t) > sizeof(uint8_t)) {
                const integer_alphabet<uint64_t> alphabet = random_integer_alphabet<uint64_t>(1, max_input_size, gen,
                    max_value_of<sym_t>());
                std::vector<uint64_t> values = random_text<std::vector<uint64_t>>(min_input_size, max_input_size,
                    alphabet.sigma, alphabet, gen);
                std::vector<sym_t> input;
                for (uint64_t value : values) input.push_back(sym_t(value));

                move_r_params params { .mode = _suffix_array, .num_threads = 1, .reorder_alphabet = gen() % 2 == 0 };
                const move_r<_locate_move> reference(values, params);
                std::stringstream expected;
                reference.serialize(expected);

                // the prefix-free parsing (for inputs of at least min_pfp_input_size symbols), else one of the two
                // modes with the suffix array
                const bool use_pfp = input.size() >= min_pfp_input_size && gen() % 2 == 0;
                params.mode = use_pfp ? _pfp : (gen() % 2 == 0 ? _suffix_array : _suffix_array_space);
                const bool from_file = gen() % 2 == 0;
                // the symbol type is passed to the constructor with the empty container symbol_type
                std::vector<sym_t> symbol_type;
                params.input_file = from_file ? write_temporary_input(input) : "";
                const move_r<_locate_move> index(from_file ? symbol_type : input, params);
                std::stringstream actual;
                index.serialize(actual);
                if (from_file) remove_input_file_and_pfp_cache(params.input_file);

                ASSERT_TRUE(actual.str() == expected.str()) << "the index differs from the one of 64-bit symbols built "
                    << "with the suffix array, " << sizeof(sym_t) << "-byte symbols, n = " << input.size() << ", sigma = "
                    << reference.alphabet_size() << ", mode " << int(params.mode)
                    << (from_file ? ", from a file" : ", from memory")
                    << (params.reorder_alphabet ? ", reordered" : ", not reordered");
            }
        };

        // the symbol widths of MOVE_R_FOR_EACH_INPUT_WIDTH, one of them drawn
        std::vector<std::function<void()>> widths;
        #define MOVE_R_TEST_INPUT_WIDTH(sym_t) widths.emplace_back([&] { run.template operator()<sym_t>(); });
        MOVE_R_FOR_EACH_INPUT_WIDTH(MOVE_R_TEST_INPUT_WIDTH)
        #undef MOVE_R_TEST_INPUT_WIDTH
        widths[gen() % widths.size()]();
    }, true } }, fuzz_iterations(num_iterations_construction), nullptr, 0);
}

/* build_with_fallback with a random construction mode, with and without -c, with an input length around
   move_r_max_input_length_sa or below it, and a build that fails for a random set of modes: the result is the mode of
   the build with the given mode, or after a failed _pfp build without -c and with at most move_r_max_input_length_sa
   symbols the one of the build with _suffix_array_space; otherwise std::runtime_error is thrown */
TEST(test_move_r, build_fallback)
{
    run_fuzz("move-r", { { "build-fallback", [](uint64_t) {
        const move_r_construction_mode modes[] { _suffix_array, _suffix_array_space, _pfp };
        const move_r_construction_mode mode = modes[gen() % std::size(modes)];
        const bool mode_given = gen() % 2 == 0;
        const uint64_t input_length = gen() % 2 == 0 ? move_r_max_input_length_sa - 1 + gen() % 3
            : random_log_uniform_size(1, move_r_max_input_length_sa, gen);
        // fails[m] <=> the build with the mode m fails
        std::map<move_r_construction_mode, bool> fails;
        for (move_r_construction_mode m : modes) fails[m] = gen() % 2 == 0;

        auto build = [&](move_r_construction_mode m) {
            if (fails[m]) throw std::runtime_error("the build failed");
            return m;
        };

        const bool falls_back = fails[mode] && mode == _pfp && !mode_given && input_length <= move_r_max_input_length_sa;
        const move_r_construction_mode used = falls_back ? _suffix_array_space : mode;
        const std::string label = "mode " + std::to_string(int(mode)) + (mode_given ? ", given" : ", not given")
            + ", n = " + std::to_string(input_length);

        if (fails[used]) {
            EXPECT_THROW(build_with_fallback(mode, mode_given, input_length, build), std::runtime_error) << label;
        } else {
            EXPECT_EQ(build_with_fallback(mode, mode_given, input_length, build), used) << label;
        }
    } } }, fuzz_iterations(num_iterations_build_fallback));
}
