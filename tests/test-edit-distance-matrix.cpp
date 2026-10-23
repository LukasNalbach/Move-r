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

#include <random>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <data_structures/edit_distance_matrix.hpp>
#include <misc/apm.hpp>
#include "test-inputs.hpp"
#include "test-patterns.hpp"
#include <misc/utils.hpp>
#include <misc/log.hpp>
#include "test-fuzz.hpp"
#include "test-parameters.hpp"

// one generator per thread (run_fuzz runs parallel functionalities on all threads)
thread_local std::mt19937 gen(test_seed() + omp_get_thread_num());
thread_local std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);

// the random inputs have a length in [min_input_size, max_input_size]; max_input_size is calibrated by run_fuzz to its
// time budget, down to min_fuzz_input_size
uint64_t min_input_size = 1;
uint64_t max_input_size = 80000;
constexpr uint64_t min_fuzz_input_size = 2000;

// number of random inputs of each fuzz test (see fuzz_iterations)
constexpr uint64_t num_iterations_matrix = 320000;
constexpr uint64_t num_iterations_distance = 320000;
constexpr uint64_t num_iterations_locate = 430;

// the two largest symbol values do not occur in the random inputs: the ranks 1..sigma of the alphabet map stay below
// uint8_t(-1), which edit_distance_matrix reserves for a symbol that does not occur in the pattern, next to rank 0
constexpr uint64_t num_excluded_symbols = 2;

// verify_edit_distance_matrix: the pattern and the text have a length of at most matrix_length_per_error times the
// largest number of errors the matrix supports; the text is a mutated copy of the pattern with probability
// prob_text_from_pattern
constexpr uint64_t matrix_length_per_error = 2;
constexpr double prob_text_from_pattern = 0.75;
// verify_edit_dist_bounded_and_cigar: P and T are empty with probability prob_empty each; otherwise P is a mutated
// substring of length at most bounded_max_length, and T is a mutated copy of P with probability prob_copy, else a
// mutated substring
constexpr uint64_t bounded_max_length = 40;
constexpr double prob_empty = 0.05;
constexpr double prob_copy = 0.5;
// verify_locate_edit_dist: k is drawn from [0, locate_max_errors] and is at most |P| / pattern_length_per_error for a
// pattern of length at most locate_max_pattern_length
constexpr uint64_t locate_max_errors = 6;
constexpr uint64_t pattern_length_per_error = 3;
constexpr uint64_t locate_max_pattern_length = 20;
// TEST(test_edit_distance_matrix, matrix) checks the 64-bit matrix with probability prob_64_bit_matrix, else the
// 128-bit one
constexpr double prob_64_bit_matrix = 0.5;

/**
 * @brief returns the edit distance of T and P (full DP matrix)
 * @tparam inp_t contiguous range type of the sequences
 * @param T the first sequence
 * @param P the second sequence
 */
template <std::ranges::contiguous_range inp_t>
static uint64_t reference_edit_distance(const inp_t& T, const inp_t& P)
{
    uint64_t n = T.size();
    uint64_t m = P.size();

    // [0..n][0..m] dist[i][j] = edit distance of T[0..i-1] and P[0..j-1]
    std::vector<std::vector<uint64_t>> dist(n + 1, std::vector<uint64_t>(m + 1));

    for (uint64_t i = 0; i <= n; i++) dist[i][0] = i;
    for (uint64_t j = 0; j <= m; j++) dist[0][j] = j;

    for (uint64_t i = 1; i <= n; i++) {
        auto c1 = T[i - 1];

        for (uint64_t j = 1; j <= m; j++) {
            auto c2 = P[j - 1];
            uint64_t cost = c1 != c2;

            dist[i][j] = std::min({
                dist[i - 1][j] + 1,        // delete
                dist[i][j - 1] + 1,        // insert
                dist[i - 1][j - 1] + cost  // match / mismatch
            });
        }
    }

    return dist[n][m];
}

/**
 * @brief returns sigma and alph_map, where alph_map[c] is the rank of byte c among the distinct bytes of input
 *        (1-based, in increasing order of the unsigned byte values), or 0 if c does not occur
 * @tparam inp_t sequence type of the input (std::string for char, else std::vector<sym_t>)
 * @param input the input sequence
 */
template <typename inp_t>
static std::pair<uint16_t, std::vector<uint8_t>> compute_alphabet_map(const inp_t& input)
{
    std::vector<uint8_t> alph_map(byte_alphabet_size, 0);

    for (uint64_t i = 0; i < input.size(); i++)
        alph_map[std::bit_cast<uint8_t>(input[i])] = 1;

    uint16_t sigma = 0;

    for (uint16_t c = 0; c < byte_alphabet_size; c++)
        if (alph_map[c] != 0) alph_map[c] = ++sigma;

    return { sigma, std::move(alph_map) };
}

/**
 * @brief returns a random substring of input of length at most max_length, with random edits (see apply_random_edits)
 *        if mutate is true
 * @tparam inp_t sequence type of the input (std::string for char, else std::vector<sym_t>)
 * @tparam gen_t random number generator type
 * @param input the input sequence to draw the substring from
 * @param max_length the maximum length of the substring before the edits
 * @param mutate whether to apply random edits to the substring
 * @param gen random number generator
 */
template <typename inp_t, typename gen_t>
static inp_t random_substring(const inp_t& input, uint64_t max_length, bool mutate, gen_t& gen)
{
    std::uniform_int_distribution<uint64_t> pos_distrib(0, input.size() - 1);
    std::uniform_int_distribution<uint64_t> length_distrib(1, std::min<uint64_t>(max_length, input.size()));

    uint64_t pos = pos_distrib(gen);
    uint64_t length = std::min<uint64_t>(input.size() - pos, length_distrib(gen));
    inp_t substring(input.begin() + pos, input.begin() + pos + length);

    if (mutate) apply_random_edits(substring, input, gen);
    return substring;
}

/**
 * @brief checks every cell of the last column of an edit_distance_matrix with value <= k against
 *        reference_edit_distance; the pattern is a mutated substring of input, the text a mutated copy of the pattern
 *        or of another substring
 * @tparam word_t the matrix word type (uint64_t or uint128_t)
 * @tparam inp_t sequence type of the input (std::string for char, else std::vector<sym_t>)
 * @tparam gen_t random number generator type
 * @param input the input sequence to draw the pattern and the text from
 * @param gen random number generator
 */
template <typename word_t, typename inp_t, typename gen_t>
static void verify_edit_distance_matrix(const inp_t& input, gen_t& gen)
{
    using matrix_t = edit_distance_matrix<word_t>;
    constexpr uint16_t max_supported_errors = matrix_t::max_supported_errors;

    auto [sigma, alph_map] = compute_alphabet_map(input);
    std::uniform_int_distribution<int32_t> k_distrib(0, max_supported_errors);
    uint16_t k = k_distrib(gen);

    inp_t P = random_substring(input, matrix_length_per_error * max_supported_errors, true, gen);
    int32_t m = P.size();

    inp_t T = prob_distrib(gen) < prob_text_from_pattern ? P
        : random_substring(input, matrix_length_per_error * max_supported_errors, false, gen);
    apply_random_edits(T, input, gen);
    int32_t n = T.size();

    matrix_t M;
    directional_substring<inp_t> pattern_view(P, 0, m - 1, RIGHT);
    M.set_input(pattern_view, sigma, alph_map);
    M.init(k);
    for (int32_t i = 1; i <= n && i < M.num_rows(); i++)
        M.compute_row(i, T[i - 1]);

    for (uint16_t i = 0; i < M.num_rows() && i <= n; i++) {
        if (!M.is_in_final_column(i)) continue;
        uint16_t value = M(i, M.num_cols() - 1);

        if (value <= k) {
            inp_t prefix(T.begin(), T.begin() + i);
            EXPECT_EQ(value, reference_edit_distance(prefix, P)) << "matrix value at final column, row " << i
                << " (plen=" << m << " tlen=" << n << " sigma=" << sigma << " k=" << k << ")";
        }
    }
}

/**
 * @brief checks edit_dist_bounded against reference_edit_distance (capped at k + 1) and, if the distance is at most
 *        k, that edit_cigar returns the distance and a CIGAR with that many edits, for two random strings P and T and
 *        a random k
 * @tparam inp_t sequence type of the input (std::string for char, else std::vector<sym_t>)
 * @tparam gen_t random number generator type
 * @param input the input sequence to draw the two strings from
 * @param gen random number generator
 */
template <typename inp_t, typename gen_t>
static void verify_edit_dist_bounded_and_cigar(const inp_t& input, gen_t& gen)
{
    constexpr uint64_t max_length = bounded_max_length;

    // k is drawn up to the maximum of the larger matrix word type
    std::uniform_int_distribution<int32_t> k_distrib(0, edit_distance_max_errors_128);
    int32_t k = k_distrib(gen);

    inp_t P = prob_distrib(gen) < prob_empty ? inp_t() : random_substring(input, max_length, true, gen);
    inp_t T = prob_distrib(gen) < prob_empty ? inp_t() : (prob_distrib(gen) < prob_copy ? P : random_substring(input, max_length, false, gen));
    apply_random_edits(T, input, gen);
    int32_t len_P = P.size(), len_T = T.size();

    int64_t dist = edit_dist_bounded(P, T, k);
    int64_t dist_ref = std::min<int64_t>(reference_edit_distance(P, T), k + 1);
    EXPECT_EQ(dist, dist_ref) << "len_P=" << len_P << " len_T=" << len_T << " k=" << k;

    if (dist <= k && len_P > 0 && len_T > 0) {
        auto [err, cigar] = edit_cigar(P, len_P, T, len_T, k);
        EXPECT_EQ(int64_t(err), dist) << "len_P=" << len_P << " len_T=" << len_T << " k=" << k;
        EXPECT_EQ(num_cigar_edits(P.data(), len_P, T.data(), len_T, cigar), int64_t(err))
            << "len_P=" << len_P << " len_T=" << len_T << " k=" << k;
    }
}

/**
 * @brief checks locate_edit_dist against reference_edit_distance for a mutated pattern drawn from input
 * @tparam inp_t sequence type of the input (std::string for char, else std::vector<sym_t>)
 * @tparam gen_t random number generator type
 * @param input the input sequence to draw the pattern from and search in
 * @param gen random number generator
 */
template <typename inp_t, typename gen_t>
static void verify_locate_edit_dist(const inp_t& input, gen_t& gen)
{
    std::uniform_int_distribution<uint64_t> k_distrib(0, locate_max_errors);
    inp_t P = random_substring(input, locate_max_pattern_length, true, gen);
    uint64_t k = std::min<uint64_t>(k_distrib(gen), P.size() / pattern_length_per_error);
    int64_t n = input.size(), m = P.size();

    // the occurrences, sorted by pos: for each pos, the shortest window with the minimum number of errors, if it is <= k
    std::vector<aprx_occ_t<>> ref;

    for (int64_t i = 0; i < n; i++) {
        aprx_occ_t<> occ(uint64_t(i), 0, k + 1);
        int64_t l_max = std::min<int64_t>(m + k, n - i);

        for (int64_t l = std::max<int64_t>(1, m - int64_t(k)); l <= l_max; l++) {
            inp_t window(input.begin() + i, input.begin() + i + l);
            uint64_t err = reference_edit_distance(window, P);
            if (err < occ.err()) { occ.set_err(err); occ.set_len(l); }
        }

        if (occ.err() <= k) ref.emplace_back(occ);
    }

    EXPECT_EQ(locate_edit_dist(input, P, k), ref) << "n=" << n << " m=" << m << " k=" << k;
}

/**
 * @brief generates a random repetitive input over the symbol type char, uint8_t or int8_t (chosen by i % 3) and calls
 *        verify(input)
 * @tparam verify_t type of verify
 * @param i the iteration
 * @param verify the check to run on the input
 */
template <typename verify_t>
static void for_random_input_type(uint64_t i, verify_t verify)
{
    auto run = [&]<typename sym_t, typename inp_t>() {
        sym_t min_sym = std::numeric_limits<sym_t>::min();
        sym_t max_sym = std::numeric_limits<sym_t>::max() - num_excluded_symbols;
        inp_t input = random_repetitive_input<inp_t>(
            min_input_size, max_input_size, min_sym, max_sym, gen());
        verify(input);
    };

    switch (i % 3) {
        case 0: run.template operator()<char,    std::string>();          break;
        case 1: run.template operator()<uint8_t, std::vector<uint8_t>>(); break;
        case 2: run.template operator()<int8_t,  std::vector<int8_t>>();  break;
    }
}

TEST(test_edit_distance_matrix, matrix)
{
    run_fuzz("edit-distance-matrix", { { "matrix", [](uint64_t i) { for_random_input_type(i, [](auto& input) {
        if (prob_distrib(gen) < prob_64_bit_matrix) verify_edit_distance_matrix<uint64_t>(input, gen);
        else                         verify_edit_distance_matrix<uint128_t>(input, gen);
    }); }, true } }, fuzz_iterations(num_iterations_matrix), &max_input_size, min_fuzz_input_size);
}

TEST(test_edit_distance_matrix, distance)
{
    run_fuzz("edit-distance-matrix", { { "distance", [](uint64_t i) {
        for_random_input_type(i, [](auto& input) { verify_edit_dist_bounded_and_cigar(input, gen); }); }, true } },
        fuzz_iterations(num_iterations_distance), &max_input_size, min_fuzz_input_size);
}

TEST(test_edit_distance_matrix, locate)
{
    run_fuzz("edit-distance-matrix", { { "locate", [](uint64_t i) {
        for_random_input_type(i, [](auto& input) { verify_locate_edit_dist(input, gen); }); }, true } },
        fuzz_iterations(num_iterations_locate), &max_input_size, min_fuzz_input_size);
}
