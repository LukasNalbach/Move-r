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

#include <map>
#include <random>
#include <sstream>
#include <vector>

#include <gtest/gtest.h>

#include <data_structures/plain_bit_vector.hpp>
#include <data_structures/sd_array.hpp>
#include <data_structures/hybrid_bit_vector.hpp>
#include <data_structures/rank_select_support.hpp>
#include <data_structures/interleaved_bit_aligned_vectors.hpp>
#include <data_structures/interleaved_byte_aligned_vectors.hpp>
#include <misc/utils.hpp>
#include <misc/log.hpp>

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);

using pos_t = uint32_t;

/**
 * @brief draws a random bit vector whose density is biased to occasionally produce sparse, dense, all-zero
 * and all-one vectors (so that the hybrid bit vector exercises both its compressed and its plain branch)
 * @param gen_t random number generator type
 * @param gen random number generator
 * @return a random bit vector
 */
template <typename gen_t>
static sdsl::bit_vector random_bit_vector(gen_t& gen)
{
    std::uniform_int_distribution<uint64_t> size_distrib(1, 20000);
    uint64_t n = size_distrib(gen);
    sdsl::bit_vector bits(n, 0);

    double r = prob_distrib(gen);
    double density = r < 0.1 ? 0.0 : (r < 0.2 ? 1.0 : prob_distrib(gen));

    for (uint64_t i = 0; i < n; i++)
        bits[i] = prob_distrib(gen) < density ? 1 : 0;

    return bits;
}

/**
 * @brief builds a bit-vector data structure from bits and verifies its rank-, select- and access-queries
 * (incl. a serialize/load round-trip) against naive references
 * @tparam bit_vector_t the bit-vector data structure type (plain_bit_vector, sd_array or hybrid_bit_vector)
 * @param bits the bit vector to build the data structure from
 */
template <typename bit_vector_t>
static void verify_bit_vector(const sdsl::bit_vector& bits)
{
    uint64_t n = bits.size();

    // naive references: prefix rank_1 and the positions of all ones/zeros
    std::vector<uint64_t> rank_1(n + 1, 0);
    std::vector<uint64_t> ones_positions;
    std::vector<uint64_t> zeros_positions;

    for (uint64_t i = 0; i < n; i++) {
        rank_1[i + 1] = rank_1[i] + (bits[i] ? 1 : 0);
        if (bits[i]) ones_positions.emplace_back(i);
        else         zeros_positions.emplace_back(i);
    }

    uint64_t num_ones = ones_positions.size();
    uint64_t num_zeros = zeros_positions.size();

    auto verify = [&](bit_vector_t& bv, const std::string& label) {
        ASSERT_EQ(bv.size(), n) << label << ": size";
        EXPECT_EQ(bv.num_ones(), num_ones) << label << ": num_ones";
        EXPECT_EQ(bv.num_zeros(), num_zeros) << label << ": num_zeros";

        for (uint64_t i = 0; i < n; i++)
            EXPECT_EQ((bool) bv[i], (bool) bits[i]) << label << ": operator[" << i << "]";

        for (uint64_t i = 0; i <= n; i++) {
            EXPECT_EQ(bv.rank_1(i), rank_1[i]) << label << ": rank_1(" << i << ")";
            EXPECT_EQ(bv.rank_0(i), i - rank_1[i]) << label << ": rank_0(" << i << ")";
        }

        // select_1/select_0 are 1-based (the i-th one/zero, i in [1, num_ones]/[1, num_zeros])
        for (uint64_t i = 1; i <= num_ones; i++)
            EXPECT_EQ(bv.select_1(i), ones_positions[i - 1]) << label << ": select_1(" << i << ")";

        for (uint64_t i = 1; i <= num_zeros; i++)
            EXPECT_EQ(bv.select_0(i), zeros_positions[i - 1]) << label << ": select_0(" << i << ")";
    };

    bit_vector_t bv(bits);
    verify(bv, "built");

    // verify that the data structure survives a serialize/load round-trip
    std::stringstream stream;
    bv.serialize(stream);
    bit_vector_t bv_reloaded;
    bv_reloaded.load(stream);
    verify(bv_reloaded, "reloaded");
}

/**
 * @brief builds an interleaved_bit_aligned_vectors from random values and verifies its get-, set- and
 * emplace_back-operations (incl. a serialize/load round-trip) against naive references
 * @tparam gen_t random number generator type
 * @param gen random number generator
 */
template <typename gen_t>
static void verify_interleaved_bit_aligned_vectors(gen_t& gen)
{
    constexpr uint8_t num_vectors = 4;
    std::array<uint64_t, num_vectors> widths = { 1, 9, 21, 33 }; // widths in bits
    interleaved_bit_aligned_vectors<uint64_t, num_vectors> vectors(widths);

    std::array<uint64_t, num_vectors> masks;
    for (uint8_t v = 0; v < num_vectors; v++)
        masks[v] = (uint64_t{1} << widths[v]) - 1;

    std::uniform_int_distribution<uint64_t> size_distrib(1, 5000);
    uint64_t n = size_distrib(gen);
    std::uniform_int_distribution<uint64_t> val_distrib(0, ULONG_MAX);

    std::array<std::vector<uint64_t>, num_vectors> reference;

    for (uint64_t i = 0; i < n; i++) {
        std::array<uint64_t, num_vectors> vals;
        for (uint8_t v = 0; v < num_vectors; v++) {
            vals[v] = val_distrib(gen) & masks[v];
            reference[v].emplace_back(vals[v]);
        }
        vectors.emplace_back(vals);
    }

    ASSERT_EQ(vectors.size(), n);

    auto verify = [&](interleaved_bit_aligned_vectors<uint64_t, num_vectors>& vecs, const std::string& label) {
        ASSERT_EQ(vecs.size(), n) << label << ": size";
        for (uint64_t i = 0; i < n; i++) {
            EXPECT_EQ(vecs.template get<0>(i), reference[0][i]) << label << ": get<0>(" << i << ")";
            EXPECT_EQ(vecs.template get<1>(i), reference[1][i]) << label << ": get<1>(" << i << ")";
            EXPECT_EQ(vecs.template get<2>(i), reference[2][i]) << label << ": get<2>(" << i << ")";
            EXPECT_EQ(vecs.template get<3>(i), reference[3][i]) << label << ": get<3>(" << i << ")";
            EXPECT_EQ(vecs[i], reference[0][i]) << label << ": operator[" << i << "]"; // operator[] == get<0>
        }
    };

    verify(vectors, "emplaced");

    // overwrite some entries with set<>() and re-verify
    std::uniform_int_distribution<uint64_t> pos_distrib(0, n - 1);
    for (uint64_t q = 0; q < n; q++) {
        uint64_t i = pos_distrib(gen);
        uint64_t v0 = val_distrib(gen) & masks[0];
        uint64_t v2 = val_distrib(gen) & masks[2];
        vectors.template set<0>(i, v0);
        vectors.template set<2>(i, v2);
        reference[0][i] = v0;
        reference[2][i] = v2;
    }

    verify(vectors, "after set");

    // overwrite every entry concurrently with set_parallel<>() and re-verify; neighbouring fields
    // share bytes, so this exercises the race-free aligned-word writes
    std::array<std::vector<uint64_t>, num_vectors> par_vals;
    for (uint8_t v = 0; v < num_vectors; v++) {
        par_vals[v].resize(n);
        for (uint64_t i = 0; i < n; i++) {
            par_vals[v][i] = val_distrib(gen) & masks[v];
            reference[v][i] = par_vals[v][i];
        }
    }

    #pragma omp parallel for
    for (uint64_t i = 0; i < n; i++) {
        vectors.template set_parallel<0>(i, par_vals[0][i]);
        vectors.template set_parallel<1>(i, par_vals[1][i]);
        vectors.template set_parallel<2>(i, par_vals[2][i]);
        vectors.template set_parallel<3>(i, par_vals[3][i]);
    }

    verify(vectors, "after set_parallel");

    // verify that the data structure survives a serialize/load round-trip
    std::stringstream stream;
    vectors.serialize(stream);
    interleaved_bit_aligned_vectors<uint64_t, num_vectors> reloaded(widths);
    reloaded.load(stream);
    verify(reloaded, "reloaded");
}

/**
 * @brief builds an interleaved_byte_aligned_vectors from random values and verifies its get-, set- and
 * emplace_back-operations (incl. a serialize/load round-trip) against naive references
 * @tparam gen_t random number generator type
 * @param gen random number generator
 */
template <typename gen_t>
static void verify_interleaved_byte_aligned_vectors(gen_t& gen)
{
    constexpr uint8_t num_vectors = 4;
    std::array<uint8_t, num_vectors> widths = { 1, 2, 3, 5 }; // widths in bytes
    interleaved_byte_aligned_vectors<uint64_t, uint64_t, num_vectors> vectors(widths);

    std::array<uint64_t, num_vectors> masks;
    for (uint8_t v = 0; v < num_vectors; v++)
        masks[v] = (uint64_t{1} << (8 * widths[v])) - 1;

    std::uniform_int_distribution<uint64_t> size_distrib(1, 5000);
    uint64_t n = size_distrib(gen);
    std::uniform_int_distribution<uint64_t> val_distrib(0, ULONG_MAX);

    std::array<std::vector<uint64_t>, num_vectors> reference;

    for (uint64_t i = 0; i < n; i++) {
        uint64_t v0 = val_distrib(gen) & masks[0];
        uint64_t v1 = val_distrib(gen) & masks[1];
        uint64_t v2 = val_distrib(gen) & masks[2];
        uint64_t v3 = val_distrib(gen) & masks[3];
        reference[0].emplace_back(v0);
        reference[1].emplace_back(v1);
        reference[2].emplace_back(v2);
        reference[3].emplace_back(v3);
        vectors.emplace_back(std::make_tuple(v0, v1, v2, v3));
    }

    ASSERT_EQ(vectors.size(), n);

    auto verify = [&](interleaved_byte_aligned_vectors<uint64_t, uint64_t, num_vectors>& vecs, const std::string& label) {
        ASSERT_EQ(vecs.size(), n) << label << ": size";
        for (uint64_t i = 0; i < n; i++) {
            EXPECT_EQ(vecs.template get<0>(i), reference[0][i]) << label << ": get<0>(" << i << ")";
            EXPECT_EQ(vecs.template get<1>(i), reference[1][i]) << label << ": get<1>(" << i << ")";
            EXPECT_EQ(vecs.template get<2>(i), reference[2][i]) << label << ": get<2>(" << i << ")";
            EXPECT_EQ(vecs.template get<3>(i), reference[3][i]) << label << ": get<3>(" << i << ")";
            EXPECT_EQ(vecs[i], reference[0][i]) << label << ": operator[" << i << "]";
        }
    };

    verify(vectors, "emplaced");

    // overwrite some entries with set<>() and re-verify
    std::uniform_int_distribution<uint64_t> pos_distrib(0, n - 1);
    for (uint64_t q = 0; q < n; q++) {
        uint64_t i = pos_distrib(gen);
        uint64_t v1 = val_distrib(gen) & masks[1];
        uint64_t v3 = val_distrib(gen) & masks[3];
        vectors.template set<1>(i, v1);
        vectors.template set<3>(i, v3);
        reference[1][i] = v1;
        reference[3][i] = v3;
    }

    verify(vectors, "after set");

    // verify that the data structure survives a serialize/load round-trip
    std::stringstream stream;
    vectors.serialize(stream);
    interleaved_byte_aligned_vectors<uint64_t, uint64_t, num_vectors> reloaded(widths);
    reloaded.load(stream);
    verify(reloaded, "reloaded");
}

/**
 * @brief verifies the rank-, select-, frequency- and contains-queries (incl. a serialize/load round-trip) of
 * a rank_select_support against naive references, for each symbol that occurs in the input
 * @tparam sym_t symbol type of the input
 * @tparam rank_select_t the rank_select_support type
 * @param input the input the rank_select_support was built from
 * @param index the rank_select_support to verify
 * @param label a short label identifying the instance for failure messages
 */
template <typename sym_t, typename rank_select_t, typename input_t>
static void verify_rank_select_queries(const input_t& input, rank_select_t& index, const std::string& label)
{
    pos_t n = input.size();
    ASSERT_EQ(index.size(), n) << label << ": size";

    // naive references: for each symbol, the ascending positions of its occurrences
    std::map<sym_t, std::vector<pos_t>> occurrences;
    for (pos_t i = 0; i < n; i++)
        occurrences[input[i]].emplace_back(i);

    for (const auto& [sym, positions] : occurrences) {
        EXPECT_TRUE(index.contains(sym)) << label << ": contains";
        EXPECT_EQ(index.frequency(sym), positions.size()) << label << ": frequency";

        // rank(sym, i) = number of occurrences of sym before index i (i in [1, n])
        pos_t next = 0;
        for (pos_t i = 1; i <= n; i++) {
            while (next < positions.size() && positions[next] < i) next++;
            EXPECT_EQ(index.rank(sym, i), next) << label << ": rank(" << (int64_t) sym << ", " << i << ")";
        }

        // select(sym, i) = index of the i-th occurrence of sym (i in [1, frequency])
        for (pos_t i = 1; i <= positions.size(); i++)
            EXPECT_EQ(index.select(sym, i), positions[i - 1]) << label << ": select(" << (int64_t) sym << ", " << i << ")";
    }
}

/**
 * @brief builds a rank_select_support over a random byte-alphabet input and verifies it
 * @tparam gen_t random number generator type
 * @param gen random number generator
 */
template <typename gen_t>
static void verify_rank_select_support_byte(gen_t& gen)
{
    std::uniform_int_distribution<uint64_t> size_distrib(1, 20000);
    std::uniform_int_distribution<int> sigma_distrib(1, 60);
    uint64_t n = size_distrib(gen);
    int sigma = sigma_distrib(gen);
    std::uniform_int_distribution<int> sym_distrib(1, sigma);

    std::string input;
    no_init_resize(input, n);
    for (uint64_t i = 0; i < n; i++)
        input[i] = (char) sym_distrib(gen);

    rank_select_support<char, pos_t> index(input);
    verify_rank_select_queries<char>(input, index, "rank_select (byte)");

    std::stringstream stream;
    index.serialize(stream);
    rank_select_support<char, pos_t> index_reloaded;
    index_reloaded.load(stream);
    verify_rank_select_queries<char>(input, index_reloaded, "rank_select (byte, reloaded)");
}

/**
 * @brief builds a rank_select_support over a random integer-alphabet input and verifies it
 * @tparam gen_t random number generator type
 * @param gen random number generator
 */
template <typename gen_t>
static void verify_rank_select_support_int(gen_t& gen)
{
    // rank_select_support requires an unsigned symbol type for its integer alphabet
    std::uniform_int_distribution<uint64_t> size_distrib(1, 20000);
    std::uniform_int_distribution<uint32_t> sigma_distrib(1, 200);
    uint64_t n = size_distrib(gen);
    uint32_t sigma = sigma_distrib(gen);
    std::uniform_int_distribution<uint32_t> sym_distrib(0, sigma - 1);

    std::vector<uint32_t> input;
    no_init_resize(input, n);
    for (uint64_t i = 0; i < n; i++)
        input[i] = sym_distrib(gen);

    rank_select_support<uint32_t, pos_t> index(input, sigma);
    verify_rank_select_queries<uint32_t>(input, index, "rank_select (int)");

    std::stringstream stream;
    index.serialize(stream);
    rank_select_support<uint32_t, pos_t> index_reloaded;
    index_reloaded.load(stream);
    verify_rank_select_queries<uint32_t>(input, index_reloaded, "rank_select (int, reloaded)");
}

TEST(test_data_structures, fuzzy_test)
{
    auto start_time = now();

    while (time_diff_min(start_time, now()) < 60) {
        double val = prob_distrib(gen);

        if (val < 1.0 / 7.0) {
            verify_bit_vector<plain_bit_vector<pos_t, true, true, true>>(random_bit_vector(gen));
        } else if (val < 2.0 / 7.0) {
            verify_bit_vector<sd_array<pos_t>>(random_bit_vector(gen));
        } else if (val < 3.0 / 7.0) {
            verify_bit_vector<hybrid_bit_vector<pos_t, true, true, true>>(random_bit_vector(gen));
        } else if (val < 4.0 / 7.0) {
            verify_interleaved_bit_aligned_vectors(gen);
        } else if (val < 5.0 / 7.0) {
            verify_interleaved_byte_aligned_vectors(gen);
        } else if (val < 6.0 / 7.0) {
            verify_rank_select_support_byte(gen);
        } else {
            verify_rank_select_support_int(gen);
        }
    }
}
