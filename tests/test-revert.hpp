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

// checks of the parallel reverts of the indexes

#ifndef MOVE_R_TESTS_TEST_REVERT_HPP
#define MOVE_R_TESTS_TEST_REVERT_HPP

#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <misc/temp_files.hpp>
#include <move_r/move_r.hpp>

// the numbers of threads the reverts are checked with
inline constexpr std::array<uint64_t, 4> revert_thread_counts { 1, 2, 3, 8 };

// the numbers of pieces each thread reverts alternately the reverts are checked with
inline constexpr std::array<uint64_t, 3> revert_strand_counts { 1, 2, 3 };

// the values of move_r_params::revert_start_samples a random index is built with (0: none)
inline constexpr std::array<uint64_t, 4> revert_start_sample_counts { 0, 1, 7, default_revert_start_samples };

/**
 * @brief compares a reverted text with the expected one; reports the first difference only
 * @param reverted the reverted text
 * @param length the length of the text
 * @param expected function returning the expected symbol at a position
 * @param what the revert that is checked, appended to the failure message
 */
template <typename text_t, typename expected_t>
void compare_reverted_text(const text_t& reverted, uint64_t length, expected_t expected, const std::string& what)
{
    ASSERT_EQ(uint64_t(reverted.size()), length) << what;

    for (uint64_t i = 0; i < length; i++) {
        if (uint64_t(reverted[i]) != expected(i)) {
            ADD_FAILURE() << what << ": " << uint64_t(reverted[i]) << " instead of " << expected(i) << " at " << i;
            return;
        }
    }
}

/**
 * @brief checks the parallel reverts of an index with each number of threads in revert_thread_counts
 *        and each number of strands in revert_strand_counts, into a container and into a file (sizeof(usym_t) bytes
 *        per symbol)
 * @tparam usym_t the unsigned symbol type of the text
 * @param length the length of the text
 * @param expected function returning the expected symbol at a position
 * @param revert function called as revert(p, k) that returns the text reverted with p threads and k strands
 * @param revert_to_file function called as revert_to_file(file, p, k) that writes the text reverted with p threads
 *        and k strands
 * @param label the parameters of the index, appended to the failure messages
 */
template <typename usym_t, typename expected_t, typename revert_t, typename revert_to_file_t>
void verify_parallel_revert(uint64_t length, expected_t expected, revert_t revert, revert_to_file_t revert_to_file,
    const std::string& label)
{
    for (uint64_t p : revert_thread_counts) {
        for (uint64_t k : revert_strand_counts) {
            const std::string what = "revert with " + std::to_string(p) + " threads and " + std::to_string(k)
                + " strands, " + label;
            compare_reverted_text(revert(p, k), length, expected, what);

            const std::string file = temporary_file_path("move-r-test-revert-");
            revert_to_file(file, p, k);
            std::vector<usym_t> from_file(std::filesystem::file_size(file) / sizeof(usym_t));
            std::ifstream(file, std::ios::binary).read(reinterpret_cast<char*>(from_file.data()),
                from_file.size() * sizeof(usym_t));
            std::filesystem::remove(file);
            compare_reverted_text(from_file, length, expected, what + ", to a file");
        }
    }
}

#endif
