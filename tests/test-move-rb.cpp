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
 * @file test-move-rb.cpp
 * @brief fuzz tests of move_rb: the indexes, the parallel construction, the documents and approximate pattern matching
 */

#include "test-move-rb-common.hpp"

// the largest input of the approximate-matching checks with error models (their reference aligns the pattern at every
// starting position), and the smallest value run_fuzz reduces it to
inline uint64_t max_model_input_size = 400;
inline constexpr uint64_t min_model_fuzz_input_size = 64;

// TEST(test_move_rb, parallel_construction): the number of inputs
inline constexpr uint64_t num_parallel_construction_inputs = 10;

// the index operations, for every support (see cycle_profiles)
TEST(test_move_rb, index)
{
    run_fuzz("move-rb", { { "index", [](uint64_t i) {
        cycle_profiles(i, _test_index_operations, _all_supports);
    }, true } }, fuzz_iterations(), &max_input_size, min_fuzz_input_size);
}

// the index operations on inputs of input_symbols_per_thread * max_num_threads symbols, for which the construction uses
// all threads (it uses at most n / move_r_min_symbols_per_thread threads, see move_r::construction::limit_num_threads)
TEST(test_move_rb, parallel_construction)
{
    run_fuzz("move-rb", { { "parallel-construction", [](uint64_t i) {
        const uint64_t n = input_symbols_per_thread * max_num_threads;
        cycle_profiles(i, _test_index_operations, _all_supports, n, n);
    } } }, fuzz_iterations(num_parallel_construction_inputs), nullptr, 0);
}

// the documents and the document queries, for the supports that locate all occurrences (see move_r_params::documents)
TEST(test_move_rb, documents)
{
    run_fuzz("move-rb", { { "documents", [](uint64_t i) {
        cycle_profiles(i, _test_documents, _multiple_locate_supports);
    }, true } }, fuzz_iterations(), &max_input_size, min_fuzz_input_size);
}

/* approximate pattern matching with the Hamming and the edit distance (see verify_apm_query), with random error models
   and without them (serial with the edit distance: one input can take long; its queries are checked on all threads,
   see verify_algorithm) */
TEST(test_move_rb, apm)
{
    run_fuzz("move-rb", {
        { "apm-hamming-model", [&](uint64_t i) {
            cycle_profiles(i, _test_apm_hamming_model, _all_supports, min_input_size, max_model_input_size);
        }, true },
        { "apm-edit-model", [&](uint64_t i) {
            cycle_profiles(i, _test_apm_edit_model, _multiple_locate_supports, min_input_size, max_model_input_size);
        } }
    }, fuzz_iterations(), &max_model_input_size, min_model_fuzz_input_size);

    run_fuzz("move-rb", {
        { "apm-hamming", [&](uint64_t i) { cycle_profiles(i, _test_apm_hamming, _all_supports); }, true },
        { "apm-edit", [&](uint64_t i) { cycle_profiles(i, _test_apm_edit, _multiple_locate_supports); } }
    }, fuzz_iterations(), &max_input_size, min_fuzz_input_size);
}
