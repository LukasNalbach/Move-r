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

#include <misc/directional_substring.hpp>
#include <misc/edit_distance_matrix.hpp>
#include <misc/apm.hpp>
#include <misc/strings.hpp>
#include <misc/utils.hpp>
#include <misc/log.hpp>

std::random_device rd;
std::mt19937 gen(rd());
uint64_t min_input_size = 1;
uint64_t max_input_size = 200;

std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);

/**
 * @brief computes the edit (Levenshtein) distance of two sequences using a full (unbanded) dynamic-programming
 *        matrix -- the independent scalar reference all bit-parallel routines are verified against
 * @tparam pos_t unsigned integer type for the distance
 * @tparam inp_t contiguous range type of the sequences
 * @param T the first sequence
 * @param P the second sequence
 * @return the edit distance of T and P
 */
template <typename pos_t, std::ranges::contiguous_range inp_t>
static pos_t edit_dist(const inp_t& T, const inp_t& P)
{
    pos_t n = T.size();
    pos_t m = P.size();

    std::vector<std::vector<pos_t>> dist(n + 1, std::vector<pos_t>(m + 1));

    for (pos_t i = 0; i <= n; i++) dist[i][0] = i;
    for (pos_t j = 0; j <= m; j++) dist[0][j] = j;

    for (pos_t i = 1; i <= n; i++) {
        auto c1 = T[i - 1];

        for (pos_t j = 1; j <= m; j++) {
            auto c2 = P[j - 1];
            pos_t cost = c1 != c2;

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
 * @brief computes the effective byte alphabet of input: its size and the byte -> one-based rank map (with the
 * ranks ordered by unsigned byte value, matching how move_r orders a byte alphabet)
 * @tparam inp_t sequence type of the input (std::string for char, else std::vector<sym_t>)
 * @param input the input sequence
 * @return the pair (alphabet size, byte -> one-based rank map)
 */
template <typename inp_t>
static std::pair<uint16_t, std::vector<uint8_t>> effective_alphabet(const inp_t& input)
{
    std::vector<uint8_t> alph_map(256, 0);

    for (uint64_t i = 0; i < input.size(); i++)
        alph_map[sym_to_uchar(input[i])] = 1;

    uint16_t sigma = 0;

    for (uint16_t c = 0; c < 256; c++)
        if (alph_map[c] != 0) alph_map[c] = ++sigma;

    return { sigma, std::move(alph_map) };
}

/**
 * @brief draws a random substring of input (of length at most max_length before mutation) and optionally
 * applies random edits to it, drawing the replacement symbols from input (so the result stays within
 * input's effective alphabet)
 * @tparam inp_t sequence type of the input (std::string for char, else std::vector<sym_t>)
 * @tparam gen_t random number generator type
 * @param input the input sequence to draw the substring from
 * @param max_length the maximum length of the drawn (pre-mutation) substring
 * @param mutate whether to apply random edits to the drawn substring
 * @param gen random number generator
 * @return the (possibly mutated) substring
 */
template <typename inp_t, typename gen_t>
static inp_t random_substring(const inp_t& input, uint64_t max_length, bool mutate, gen_t& gen)
{
    std::uniform_int_distribution<uint64_t> pos_distrib(0, input.size() - 1);
    std::uniform_int_distribution<uint64_t> length_distrib(1, std::min<uint64_t>(max_length, input.size()));

    uint64_t pos = pos_distrib(gen);
    uint64_t length = std::min<uint64_t>(input.size() - pos, length_distrib(gen));
    inp_t substring(input.begin() + pos, input.begin() + pos + length);

    if (mutate) mutate_pattern<uint64_t>(substring, input, gen);
    return substring;
}

/**
 * @brief builds an edit_distance_matrix over a pattern drawn from input, feeds it a mutated copy of the
 *        pattern (or an unrelated substring) as the query, and verifies that every final-column cell's edit
 *        distance (when within the error budget) matches the unbanded reference; the alphabet map covers all of
 *        input (as an index' map_int does), so every pattern and query symbol has a rank
 * @tparam word_t the matrix word type (uint64_t or __uint128_t)
 * @tparam inp_t sequence type of the input (std::string for char, else std::vector<sym_t>)
 * @tparam gen_t random number generator type
 * @param input the input sequence to draw the pattern and query from
 * @param gen random number generator
 */
template <typename word_t, typename inp_t, typename gen_t>
static void verify_edit_distance_matrix(const inp_t& input, gen_t& gen)
{
    using matrix_t = edit_distance_matrix<word_t>;
    constexpr uint16_t k_limit = matrix_t::k_limit;

    auto [sigma, alph_map] = effective_alphabet(input);
    std::uniform_int_distribution<int32_t> k_distrib(0, k_limit);
    uint16_t k = k_distrib(gen);

    inp_t pat = random_substring(input, 2 * k_limit, true, gen);
    int32_t plen = pat.size();

    // the query text is a further-mutated copy of the pattern most of the time (so that, even over the
    // exhaustive alphabet, many band cells stay within the error budget) and an unrelated substring otherwise
    inp_t text = prob_distrib(gen) < 0.75 ? pat : random_substring(input, 2 * k_limit, false, gen);
    mutate_pattern<uint64_t>(text, input, gen);
    int32_t tlen = text.size();

    matrix_t M;
    directional_substring<uint32_t, inp_t> pat_col(pat, 0, plen - 1, RIGHT);
    M.set_input(pat_col, sigma, alph_map);
    M.init(k);
    for (int32_t i = 1; i <= tlen && i < M.num_rows(); i++)
        M.compute_row(i, text[i - 1]);

    for (uint16_t i = 0; i < M.num_rows() && i <= tlen; i++) {
        if (!M.is_in_final_column(i)) continue;
        uint16_t mv = M(i, M.num_cols() - 1);

        // the matrix value must match the unbanded reference whenever it is within the error budget
        if (mv <= k) {
            inp_t prefix(text.begin(), text.begin() + i);
            EXPECT_EQ(mv, edit_dist<uint16_t>(prefix, pat)) << "matrix value at final column, row " << i
                << " (plen=" << plen << " tlen=" << tlen << " sigma=" << sigma << " k=" << k << ")";
        }
    }
}

/**
 * @brief verifies the bit-parallel edit-distance routines against the unbanded reference over a pair of
 *        strings drawn from input and a random error budget: edit_dist_bounded must equal the clamped unbanded
 *        reference, and (when the distance is within the budget) edit_cigar must report the exact distance and
 *        a CIGAR that validly aligns the two strings with exactly that many edits
 * @tparam inp_t sequence type of the input (std::string for char, else std::vector<sym_t>)
 * @tparam gen_t random number generator type
 * @param input the input sequence to draw the two strings from
 * @param gen random number generator
 */
template <typename inp_t, typename gen_t>
static void verify_edit_dist_routines(const inp_t& input, gen_t& gen)
{
    using pos_t = uint32_t;
    std::uniform_int_distribution<int32_t> k_distrib(0, 20); // spans both matrix word types (k_limit 10 / 20)
    int32_t k = k_distrib(gen);

    // M starts as a copy of P half the time (so that pairs within the error budget occur frequently) and as
    // an independent substring otherwise, and is always further mutated; either may be empty (the routines
    // must handle empty sequences)
    inp_t P = prob_distrib(gen) < 0.05 ? inp_t() : random_substring(input, 40, true, gen);
    inp_t M = prob_distrib(gen) < 0.05 ? inp_t() : (prob_distrib(gen) < 0.5 ? P : random_substring(input, 40, false, gen));
    mutate_pattern<uint64_t>(M, input, gen);
    int32_t n1 = P.size(), n2 = M.size();

    // the bit-parallel bounded distance must match the clamped unbanded reference
    int64_t bp = edit_dist_bounded<pos_t>(P, M, k);
    int64_t full = std::min<int64_t>(edit_dist<pos_t>(P, M), k + 1);
    EXPECT_EQ(bp, full) << "n1=" << n1 << " n2=" << n2 << " k=" << k;

    // within the budget, edit_cigar reports the exact distance and a CIGAR that validly aligns P against M
    if (bp <= k && n1 > 0 && n2 > 0) {
        auto [err, cigar] = edit_cigar<pos_t>(P, n1, M, n2, k);
        EXPECT_EQ(int64_t(err), bp) << "n1=" << n1 << " n2=" << n2 << " k=" << k;
        EXPECT_EQ(num_cigar_edits(P.data(), n1, M.data(), n2, cigar), int64_t(err))
            << "n1=" << n1 << " n2=" << n2 << " k=" << k;
    }
}

/**
 * @brief verifies the bit-parallel edit-distance locate against an inline unbanded reference, searching input
 *        for a mutated pattern drawn from it
 * @tparam inp_t sequence type of the input (std::string for char, else std::vector<sym_t>)
 * @tparam gen_t random number generator type
 * @param input the input sequence to draw the pattern from and search in
 * @param gen random number generator
 */
template <typename inp_t, typename gen_t>
static void verify_locate_edit_dist(const inp_t& input, gen_t& gen)
{
    using pos_t = uint32_t;
    std::uniform_int_distribution<pos_t> k_distrib(0, 6);
    inp_t P = random_substring(input, 20, true, gen);
    pos_t k = k_distrib(gen);
    int64_t n = input.size(), m = P.size();

    // reference: for each starting position, the shortest window with the minimum (within-budget) edit distance
    std::vector<aprx_occ_t<pos_t>> ref;

    for (int64_t i = 0; i < n; i++) {
        aprx_occ_t<pos_t> occ{.pos = pos_t(i), .len = 0, .err = k + 1};
        int64_t l_max = std::min<int64_t>(m + k, n - i);

        for (int64_t l = std::max<int64_t>(1, m - int64_t(k)); l <= l_max; l++) {
            inp_t window(input.begin() + i, input.begin() + i + l);
            pos_t err = edit_dist<pos_t>(window, P);
            if (err < occ.err) { occ.err = err; occ.len = l; } // ascending l: ties keep the shortest window
        }

        if (occ.err <= k) ref.emplace_back(occ);
    }

    EXPECT_EQ(locate_edit_dist<pos_t>(input, P, k), ref) << "n=" << n << " m=" << m << " k=" << k;
}

/**
 * @brief runs one round of every verification over the given input value type: generates a random repetitive
 * input over the exhaustive alphabet and verifies the matrix (alternating the two word types) and the
 * edit-distance routines on strings drawn from it
 * @tparam sym_t symbol type of the input (a one-byte type: char, uint8_t or int8_t)
 * @tparam inp_t sequence type of the input (std::string for char, else std::vector<sym_t>)
 */
template <typename sym_t, typename inp_t>
static void run_edit_distance_once()
{
    // a byte alphabet reserves the two largest values (for the sentinel and the BWT-escape), matching the
    // alphabet the apm machinery sees in the move-r indexes
    sym_t min_sym = std::numeric_limits<sym_t>::min();
    sym_t max_sym = std::numeric_limits<sym_t>::max() - 2;

    inp_t input = random_repetitive_input<inp_t>(min_input_size, max_input_size, min_sym, max_sym);

    // alternate the two matrix word types (k up to 10 / up to 20)
    if (prob_distrib(gen) < 0.5) verify_edit_distance_matrix<uint64_t>(input, gen);
    else                         verify_edit_distance_matrix<__uint128_t>(input, gen);

    verify_edit_dist_routines(input, gen);
    verify_locate_edit_dist(input, gen);
}

// the edit-distance machinery over every supported input value type (the byte alphabets char/uint8_t/int8_t),
// rotating over the two matrix word types
TEST(test_edit_distance_matrix, fuzzy_test)
{
    std::uniform_int_distribution<int> sym_type_distrib(0, 2);
    auto start_time = now();

    while (time_diff_min(start_time, now()) < 60) {
        switch (sym_type_distrib(gen)) {
            case 0: run_edit_distance_once<char,    std::string>();          break;
            case 1: run_edit_distance_once<uint8_t, std::vector<uint8_t>>(); break;
            case 2: run_edit_distance_once<int8_t,  std::vector<int8_t>>();  break;
        }
    }
}
