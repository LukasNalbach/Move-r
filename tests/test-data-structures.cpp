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

#include <algorithm>
#include <array>
#include <clocale>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <vector>

#include <gtest/gtest.h>

#include <data_structures/plain_bit_vector.hpp>
#include <data_structures/elias_fano_bit_vector.hpp>
#include <data_structures/rank_select_set.hpp>
#include <data_structures/blocked_positions.hpp>
#include <data_structures/wavelet_tree_prefix_sum.hpp>
#include <data_structures/pred_succ_occ.hpp>
#include <data_structures/symbol_groups.hpp>
#include <data_structures/pred_succ_block.hpp>
#include <data_structures/interleaved_bit_aligned_vectors.hpp>
#include <data_structures/uint_bytes_array.hpp>
#include <misc/utils.hpp>
#include <misc/parse.hpp>
#include <misc/portable_layout.hpp>
#include "test-inputs.hpp"
#include <misc/log.hpp>
#include "test-fuzz.hpp"

// one generator per thread (run_fuzz runs parallel functionalities on all threads)
thread_local std::mt19937 gen(test_seed() + omp_get_thread_num());
thread_local std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);

// maximum input size, calibrated by run_fuzz to its time budget
uint64_t max_input_size = 250000;
uint64_t max_num_threads = omp_get_max_threads();
// smallest value run_fuzz reduces max_input_size to
constexpr uint64_t min_fuzz_input_size = 2000;

// number of random inputs of each fuzz test (see fuzz_iterations; the tests without an entry use the default)
constexpr uint64_t num_iterations_plain_bit_vector = 138000;
constexpr uint64_t num_iterations_elias_fano_bit_vector = 30000;
constexpr uint64_t num_iterations_interleaved_bit_aligned_vectors = 96000;
constexpr uint64_t num_iterations_wavelet_tree_prefix_sum = 500;
constexpr uint64_t num_iterations_pred_succ_occ = 600;
constexpr uint64_t num_iterations_symbol_groups = 600;
constexpr uint64_t num_iterations_pred_succ_block = 6000;
constexpr uint64_t num_iterations_packed_bit_fields = 10000;

// maximum input size of a fuzz test relative to max_input_size (the tests without an entry use max_input_size)
constexpr double size_factor_interleaved_bit_aligned_vectors = 0.36;
constexpr double size_factor_wavelet_tree_prefix_sum = 0.08;
constexpr double size_factor_pred_succ_occ = 0.08;
constexpr double size_factor_symbol_groups = 0.08;
constexpr double size_factor_pred_succ_block = 0.80;

// random_bit_vector: all zeros with probability prob_all_zeros, all ones with probability prob_all_ones
constexpr double prob_all_zeros = 0.1;
constexpr double prob_all_ones = 0.1;
// verify_uint_bytes_array: the largest value is less than max_value_slack below the largest value of its width; the
// array is sorted with at most max_sort_threads threads
constexpr uint64_t max_value_slack = 1024;
constexpr uint64_t max_sort_threads = 8;
// verify_interleaved_bit_aligned_vectors: the widths in bits of the tested vectors (the remaining vectors are unused)
constexpr std::array<uint8_t, 4> tested_vector_widths { 1, 9, 21, 33 };
// verify_wavelet_tree_prefix_sum: the largest alphabet size and weight, and the number of queries
constexpr uint64_t prefix_sum_max_sigma = 300;
constexpr uint64_t prefix_sum_max_weight = 1'000'000;
constexpr uint64_t prefix_sum_num_queries = 200;
// verify_pred_succ_occ: the input has at least pso_min_input_size symbols; sigma is drawn from [pso_min_sigma,
// pso_max_sigma] (segmented input: [pso_min_sigma, pso_max_segmented_sigma]); every pso_frequent_symbol_period-th
// symbol of an unsegmented input is drawn from the max(pso_min_sigma, sigma / pso_sigma_per_frequent_symbol) smallest
// symbols; at most pso_max_occs_per_block occurrences per block and pso_max_num_queries queries
constexpr uint64_t pso_min_input_size = 2;
constexpr uint64_t pso_min_sigma = 2;
constexpr uint64_t pso_max_sigma = 3000;
constexpr uint64_t pso_max_segmented_sigma = 8;
constexpr uint64_t pso_frequent_symbol_period = 3;
constexpr uint64_t pso_sigma_per_frequent_symbol = 100;
constexpr uint64_t pso_max_occs_per_block = 4;
constexpr uint64_t pso_max_num_queries = 2000;
// verify_symbol_groups: the largest alphabet size; the group size is 2^e for a random e in [0, max_log_group_size]
constexpr uint64_t symbol_groups_max_sigma = 20;
constexpr uint64_t max_log_group_size = 8;
// verify_pred_succ_block: the largest alphabet size and number of queries
constexpr int32_t psb_max_sigma = 60;
constexpr uint64_t psb_max_num_queries = 4000;

// serializes ds and loads it into loaded
template <typename ds_t>
static void serialize_and_load(const ds_t& ds, ds_t& loaded)
{
    std::stringstream stream;
    ds.serialize(stream);
    loaded.load(stream);
}

/**
 * @brief returns a random bit vector of random length in [1, max_size]: all zeros with probability prob_all_zeros, all
 *        ones with probability prob_all_ones, else each bit is 1 with a random probability
 * @tparam gen_t random number generator type
 * @param gen random number generator
 * @param max_size maximum size of the bit vector
 */
template <typename gen_t>
static sdsl::bit_vector random_bit_vector(gen_t& gen, uint64_t max_size)
{
    uint64_t n = random_log_uniform_size(1, max_size, gen);
    sdsl::bit_vector bits(n, 0);

    double r = prob_distrib(gen);
    double density = r < prob_all_zeros ? 0.0 : (r < prob_all_zeros + prob_all_ones ? 1.0 : prob_distrib(gen));

    for (uint64_t i = 0; i < n; i++)
        bits[i] = prob_distrib(gen) < density ? 1 : 0;

    return bits;
}

/**
 * @brief verifies rank_select_set built from the ones of bits (dense or sparse, depending on the density) against
 *        naive references, before and after serialize/load
 * @param bits the bit vector whose ones are the elements
 */
static void verify_rank_select_set(const sdsl::bit_vector& bits)
{
    const uint64_t n = bits.size();
    // [0..num_ones-1] the positions of the ones of bits, increasing
    std::vector<uint64_t> values;
    for (uint64_t i = 0; i < n; i++) if (bits[i]) values.push_back(i);
    const uint64_t threads = 1 + gen() % max_num_threads;
    rank_select_set built(values, n, threads);
    rank_select_set loaded;
    serialize_and_load(built, loaded);

    for (const rank_select_set* set : { &built, &loaded }) {
        ASSERT_EQ(set->num_ones(), values.size());
        ASSERT_EQ(set->size(), n);
        uint64_t rank = 0;

        for (uint64_t i = 0; i < n; i++) {
            ASSERT_EQ(set->rank_1(i), rank) << "rank_1(" << i << "), dense " << set->is_dense();
            ASSERT_EQ(set->contains(i), bool(bits[i])) << "contains(" << i << ")";
            ASSERT_EQ(set->index_of(i), bits[i] ? rank : values.size()) << "index_of(" << i << ")";
            rank += bits[i];
        }

        ASSERT_EQ(set->rank_1(n), values.size()) << "rank_1(n)";
        for (uint64_t r = 1; r <= values.size(); r++) ASSERT_EQ(set->select_1(r), values[r - 1]) << "select_1(" << r << ")";
    }
}

/**
 * @brief sets random groups of values (of 4, 5 or 6 bytes each) in a uint_bytes_array and verifies them before and
 *        after sorting by the last field
 * @tparam fields number of values per group
 * @param max_size maximum number of groups
 */
template <uint8_t fields>
static void verify_uint_bytes_array(uint64_t max_size)
{
    // the values have min_value_bytes to max_value_bytes bytes (the widths of uint_bytes_array, see bytes_per_value)
    constexpr uint64_t min_value_bytes = sizeof(uint32_t);
    constexpr uint64_t max_value_bytes = sizeof(uint48_t);

    const uint64_t size = random_log_uniform_size(1, max_size, gen);
    const uint8_t bytes = uint8_t(min_value_bytes + gen() % (max_value_bytes - min_value_bytes + 1));
    const uint64_t max_value = (uint64_t { 1 } << (8 * bytes)) - 1 - gen() % max_value_slack;
    // [0..size-1] groups[i][f] = value f of group i
    std::vector<std::array<uint64_t, fields>> groups(size);
    uint_bytes_array<fields> array(size, max_value);
    if (max_value > max_narrow_value) ASSERT_EQ(array.bytes_per_value(), max_value_bytes);

    for (uint64_t i = 0; i < size; i++) {
        for (uint8_t f = 0; f < fields; f++) groups[i][f] = gen() % (max_value + 1);
        array.template set<0>(i, groups[i][0]);
        if constexpr (fields > 1) array.template set<1>(i, groups[i][1]);
    }

    auto check = [&](const char* label) {
        for (uint64_t i = 0; i < size; i++) {
            ASSERT_EQ(array.template get<0>(i), groups[i][0]) << label << " [" << i << "]";
            if constexpr (fields > 1) ASSERT_EQ(array.template get<1>(i), groups[i][1]) << label << " [" << i << "]";
        }
    };

    ASSERT_NO_FATAL_FAILURE(check("set"));
    array.template sort_by<fields - 1>(1 + gen() % max_sort_threads);
    std::stable_sort(groups.begin(), groups.end(),
        [](const auto& a, const auto& b) { return a[fields - 1] < b[fields - 1]; });

    // the keys are sorted; groups with equal keys may be in any order, so the groups are compared as a multiset
    for (uint64_t i = 0; i < size; i++) {
        ASSERT_EQ(array.template get<fields - 1>(i), groups[i][fields - 1]) << "sorted [" << i << "]";
    }

    // [0..size-1] groups_read[i][f] = value f of group i of the sorted array
    std::vector<std::array<uint64_t, fields>> groups_read(size);

    for (uint64_t i = 0; i < size; i++) {
        groups_read[i][0] = array.template get<0>(i);
        if constexpr (fields > 1) groups_read[i][1] = array.template get<1>(i);
    }

    std::sort(groups_read.begin(), groups_read.end());
    std::sort(groups.begin(), groups.end());
    ASSERT_EQ(groups_read, groups) << "sorted groups";
}

/**
 * @brief verifies blocked_positions built from random increasing positions (with gaps that fit into 16 or 32 bits or
 *        not) against naive references, before and after serialize/load
 * @param max_size maximum number of positions
 */
static void verify_blocked_positions(uint64_t max_size)
{
    const uint64_t num = random_log_uniform_size(1, max_size, gen);
    const uint64_t max_position = (uint64_t { 1 } << (8 * sizeof(uint48_t))) - 1;
    const uint64_t max_gap = random_log_uniform_size(1, max_position / num, gen);
    // [0..num-1] the positions, strictly increasing, below 2^48
    std::vector<uint64_t> values(num);
    uint64_t position = gen() % max_gap;

    for (uint64_t i = 0; i < num; i++) {
        values[i] = position;
        position += 1 + gen() % max_gap;
    }

    blocked_positions built(values);
    blocked_positions loaded;
    serialize_and_load(built, loaded);

    for (const blocked_positions* seq : { &built, &loaded }) {
        ASSERT_EQ(seq->size(), num);

        for (uint64_t i = 0; i < num; i++) {
            ASSERT_EQ((*seq)[i], values[i]) << "[" << i << "], maximum gap " << max_gap;
            if (i + 1 < num) ASSERT_EQ(seq->pair_at(i), std::make_pair(values[i], values[i + 1])) << "pair_at(" << i << ")";
            ASSERT_EQ(seq->find(values[i]), i) << "find(" << values[i] << ")";
            if (i + 1 < num && values[i] + 1 < values[i + 1]) ASSERT_EQ(seq->find(values[i] + 1), num) << "find absent";
        }

        if (values[0] > 0) ASSERT_EQ(seq->find(values[0] - 1), num) << "find before the first position";
    }
}

/**
 * @brief builds a bit vector data structure from bits and verifies its rank, select and access queries against
 *        naive references, before and after serialize/load
 * @tparam bit_vector_t the bit vector data structure (plain_bit_vector, elias_fano_bit_vector)
 * @param bits the bit vector to build the data structure from
 */
template <typename bit_vector_t>
static void verify_bit_vector(const sdsl::bit_vector& bits)
{
    uint64_t n = bits.size();

    // [0..n] rank_1[i] = number of ones in bits[0..i-1]
    std::vector<uint64_t> rank_1(n + 1, 0);
    // [0..num_ones-1] the positions of the ones, increasing
    std::vector<uint64_t> ones_positions;
    // [0..num_zeros-1] the positions of the zeros, increasing
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

        // select_1(i) / select_0(i) return the position of the i-th one / zero, i in [1, num_ones] / [1, num_zeros]
        for (uint64_t i = 1; i <= num_ones; i++)
            EXPECT_EQ(bv.select_1(i), ones_positions[i - 1]) << label << ": select_1(" << i << ")";

        for (uint64_t i = 1; i <= num_zeros; i++)
            EXPECT_EQ(bv.select_0(i), zeros_positions[i - 1]) << label << ": select_0(" << i << ")";
    };

    // builds the bit vector with all rank and select supports (see bit_vector_support_t)
    auto build = [&]() {
        if constexpr (requires { bit_vector_t(bits, RANK_1 | SELECT_0 | SELECT_1); })
             return bit_vector_t(bits, RANK_1 | SELECT_0 | SELECT_1);
        else return bit_vector_t(bits);
    };

    bit_vector_t bv = build();
    ASSERT_NO_FATAL_FAILURE(verify(bv, "built"));
    bit_vector_t bv_reloaded;
    serialize_and_load(bv, bv_reloaded);
    verify(bv_reloaded, "reloaded");
}

/**
 * @brief builds an interleaved_bit_aligned_vectors from random values and verifies get, set, set_parallel and
 *        emplace_back against naive references, before and after serialize/load
 * @tparam gen_t random number generator type
 * @param gen random number generator
 * @param max_size maximum number of entries
 */
template <typename gen_t>
static void verify_interleaved_bit_aligned_vectors(gen_t& gen, uint64_t max_size)
{
    constexpr uint8_t num_vectors = tested_vector_widths.size();
    // widths in bits of the tested vectors; the remaining vectors are unused (width 0)
    std::array<uint8_t, interleaved_bit_aligned_vectors::num_vectors> widths {};
    std::copy(tested_vector_widths.begin(), tested_vector_widths.end(), widths.begin());
    interleaved_bit_aligned_vectors vectors(widths);

    // [0..3] masks[v] = 2^widths[v] - 1
    std::array<uint64_t, num_vectors> masks;
    for (uint8_t v = 0; v < num_vectors; v++)
        masks[v] = (uint64_t{1} << widths[v]) - 1;

    uint64_t n = random_log_uniform_size(1, max_size, gen);
    std::uniform_int_distribution<uint64_t> val_distrib(0, UINT64_MAX);

    // [0..3][0..n-1] reference[v][i] = expected value of vector v at index i
    std::array<std::vector<uint64_t>, num_vectors> reference;

    for (uint64_t i = 0; i < n; i++) {
        std::array<uint64_t, interleaved_bit_aligned_vectors::num_vectors> vals {};
        for (uint8_t v = 0; v < num_vectors; v++) {
            vals[v] = val_distrib(gen) & masks[v];
            reference[v].emplace_back(vals[v]);
        }
        vectors.emplace_back(vals);
    }

    ASSERT_EQ(vectors.size(), n);

    auto verify = [&](interleaved_bit_aligned_vectors& vecs, const std::string& label) {
        ASSERT_EQ(vecs.size(), n) << label << ": size";
        for (uint64_t i = 0; i < n; i++) {
            EXPECT_EQ(vecs.template get<0>(i), reference[0][i]) << label << ": get<0>(" << i << ")";
            EXPECT_EQ(vecs.template get<1>(i), reference[1][i]) << label << ": get<1>(" << i << ")";
            EXPECT_EQ(vecs.template get<2>(i), reference[2][i]) << label << ": get<2>(" << i << ")";
            EXPECT_EQ(vecs.template get<3>(i), reference[3][i]) << label << ": get<3>(" << i << ")";
            EXPECT_EQ(vecs[i], reference[0][i]) << label << ": operator[" << i << "]"; // operator[] == get<0>
        }
    };

    ASSERT_NO_FATAL_FAILURE(verify(vectors, "emplaced"));

    // overwrite n random entries of the vectors 0 and 2 with set
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

    ASSERT_NO_FATAL_FAILURE(verify(vectors, "after set"));

    // [0..3][0..n-1] new_values[v][i] = value written to vector v at index i by set_parallel
    std::array<std::vector<uint64_t>, num_vectors> new_values;
    for (uint8_t v = 0; v < num_vectors; v++) {
        new_values[v].resize(n);
        for (uint64_t i = 0; i < n; i++) {
            new_values[v][i] = val_distrib(gen) & masks[v];
            reference[v][i] = new_values[v][i];
        }
    }

    // overwrite every entry with set_parallel from several threads
    #pragma omp parallel for
    for (uint64_t i = 0; i < n; i++) {
        vectors.template set_parallel<0>(i, new_values[0][i]);
        vectors.template set_parallel<1>(i, new_values[1][i]);
        vectors.template set_parallel<2>(i, new_values[2][i]);
        vectors.template set_parallel<3>(i, new_values[3][i]);
    }

    ASSERT_NO_FATAL_FAILURE(verify(vectors, "after set_parallel"));

    std::stringstream stream;
    vectors.serialize(stream);
    interleaved_bit_aligned_vectors reloaded(widths);
    reloaded.load(stream);
    verify(reloaded, "reloaded");
}

/**
 * @brief builds a wavelet_tree_prefix_sum over a random weighted sequence and verifies smaller(i, j, c) against
 *        sums computed by scanning the range, before and after serialize/load
 * @tparam gen_t random number generator type
 * @param gen random number generator
 * @param max_size maximum length of the generated sequence
 */
template <typename gen_t>
static void verify_wavelet_tree_prefix_sum(gen_t& gen, uint64_t max_size)
{
    constexpr uint64_t max_sigma = prefix_sum_max_sigma;
    constexpr uint64_t max_weight = prefix_sum_max_weight;
    constexpr uint64_t num_queries = prefix_sum_num_queries;

    std::uniform_int_distribution<uint64_t> sigma_distrib(1, max_sigma);
    uint64_t n = random_log_uniform_size(1, max_size, gen);
    uint64_t sigma = sigma_distrib(gen);
    std::uniform_int_distribution<uint64_t> sym_distrib(0, sigma - 1);
    std::uniform_int_distribution<uint64_t> weight_distrib(1, max_weight);

    // [0..n-1] the symbol in [0..sigma-1] of position i
    std::vector<uint64_t> syms(n);
    // [0..n-1] the weight in [1..max_weight] of position i
    std::vector<uint64_t> weights(n);

    for (uint64_t i = 0; i < n; i++) {
        syms[i] = sym_distrib(gen);
        weights[i] = weight_distrib(gen);
    }

    // checks smaller(i, j, c), the sum of the weights of the positions in [i, j) with a symbol smaller than c
    auto check_queries = [&](const wavelet_tree_prefix_sum& index, const std::string& name) {
        for (uint64_t q = 0; q < num_queries; q++) {
            uint64_t i = std::uniform_int_distribution<uint64_t>(0, n)(gen);
            uint64_t j = std::uniform_int_distribution<uint64_t>(0, n)(gen);
            if (i > j) std::swap(i, j);
            uint64_t c = std::uniform_int_distribution<uint64_t>(0, sigma - 1)(gen);

            uint64_t expected = 0;
            for (uint64_t k = i; k < j; k++) if (syms[k] < c) expected += weights[k];

            EXPECT_EQ(index.smaller(i, j, c), expected)
                << name << ": n " << n << ", sigma " << sigma << ", range [" << i << "," << j
                << "), symbol " << c;
        }
    };

    no_init_vector<uint32_t> sym_arg(n);
    for (uint64_t i = 0; i < n; i++) sym_arg[i] = uint32_t(syms[i]);
    wavelet_tree_prefix_sum index(std::move(sym_arg), no_init_vector<uint64_t>(weights.begin(), weights.end()), sigma);
    check_queries(index, "wavelet_tree_prefix_sum");

    wavelet_tree_prefix_sum index_reloaded;
    serialize_and_load(index, index_reloaded);
    check_queries(index_reloaded, "wavelet_tree_prefix_sum (reloaded)");
}

/**
 * @brief compares succ and pred of idx with a scan of input for random queries
 * @tparam sym_t type of the query symbols
 * @param idx the data structure built over input (pred_succ_occ or pred_succ_block)
 * @param input [0..n-1] the input
 * @param sym_distrib distribution of the symbols of the succ queries
 * @param num_queries number of queries
 * @param name name of idx in the failure messages
 * @param gen random number generator
 */
template <typename sym_t, typename idx_t, typename input_t, typename sym_distrib_t, typename gen_t>
static void verify_pred_succ_queries(const idx_t& idx, const input_t& input, sym_distrib_t& sym_distrib,
    uint64_t num_queries, const std::string& name, gen_t& gen)
{
    const uint64_t n = input.size();
    std::uniform_int_distribution<uint64_t> pos_distrib(0, n - 1);

    for (uint64_t q = 0; q < num_queries; q++) {
        const uint64_t x = pos_distrib(gen);

        // succ(sym, x, max_pos): the smallest j in [x, max_pos] with input[j] = sym, else a value > max_pos
        const sym_t sym = sym_t(sym_distrib(gen));
        const uint64_t max_pos = x + pos_distrib(gen) % (n - x);
        uint64_t ref_succ = max_pos + 1;
        for (uint64_t j = x; j <= max_pos; j++) if (input[j] == sym) { ref_succ = j; break; }
        const uint64_t succ = idx.succ(sym, x, max_pos);
        if (ref_succ <= max_pos) EXPECT_EQ(succ, ref_succ) << name << " succ(" << uint64_t(sym) << "," << x << "," << max_pos << ")";
        else EXPECT_GT(succ, max_pos) << name << " succ(" << uint64_t(sym) << "," << x << "," << max_pos << ") should report none";

        // pred(sym, x, min_pos): the largest j in [min_pos, x] with input[j] = sym; sym occurs in input[min_pos..x]
        const uint64_t min_pos = pos_distrib(gen) % (x + 1);
        const sym_t sym_pred = input[min_pos + pos_distrib(gen) % (x - min_pos + 1)];
        uint64_t ref_pred = x;
        for (uint64_t j = x; ; j--) { if (input[j] == sym_pred) { ref_pred = j; break; } if (j == min_pos) break; }
        EXPECT_EQ(idx.pred(sym_pred, x, min_pos), ref_pred)
            << name << " pred(" << uint64_t(sym_pred) << "," << x << "," << min_pos << ")";
    }
}

/**
 * @brief builds a pred_succ_occ over a random input and verifies its succ and pred queries against a scan, before
 *        and after serialize/load
 * @tparam gen_t random number generator type
 * @param gen random number generator
 * @param max_size maximum length of the generated input
 */
template <typename gen_t>
static void verify_pred_succ_occ(gen_t& gen, uint64_t max_size)
{
    constexpr uint64_t min_sigma = pso_min_sigma;
    constexpr uint64_t max_sigma = pso_max_sigma;
    constexpr uint64_t max_segmented_sigma = pso_max_segmented_sigma;
    constexpr uint64_t max_occs_per_block = pso_max_occs_per_block;
    constexpr uint64_t max_num_queries = pso_max_num_queries;

    std::uniform_int_distribution<uint64_t> sigma_distrib(min_sigma, max_sigma);
    uint64_t n = random_log_uniform_size(pso_min_input_size, max_size, gen);
    // segmented <=> the input alternates between segments over the lower and over the upper half of the alphabet
    const bool segmented = gen() % 2 == 0;
    uint64_t sigma = segmented ? min_sigma + gen() % (max_segmented_sigma - min_sigma + 1) : sigma_distrib(gen);
    std::uniform_int_distribution<uint64_t> sym_distrib(0, sigma - 1);
    std::uniform_int_distribution<uint64_t> occs_per_block_distrib(1, max_occs_per_block);
    // [0..n-1] the input, symbols in [0..sigma-1]
    no_init_vector<uint32_t> input;
    resize_no_init(input, n);

    if (segmented) {
        const uint64_t half = sigma / 2;
        bool upper = false;

        for (uint64_t i = 0; i < n; upper = !upper) {
            const uint64_t end = std::min(n, i + random_log_uniform_size(1, n, gen));
            for (; i < end; i++) input[i] = uint32_t(upper ? half + gen() % (sigma - half) : gen() % half);
        }
    } else {
        // every pso_frequent_symbol_period-th symbol is drawn from the frequent symbols (see
        // pso_sigma_per_frequent_symbol), the others from the whole alphabet
        for (uint64_t i = 0; i < n; i++) {
            input[i] = uint32_t(i % pso_frequent_symbol_period == 0
                ? sym_distrib(gen) % std::max<uint64_t>(min_sigma, sigma / pso_sigma_per_frequent_symbol)
                : sym_distrib(gen));
        }
    }

    const uint64_t num_queries = std::min<uint64_t>(n, max_num_queries);
    auto read = [&input](uint64_t i) -> uint64_t { return input[i]; };
    uint64_t occs_per_block = occs_per_block_distrib(gen);

    // build with a random number of threads
    uint64_t threads = std::uniform_int_distribution<uint64_t>(1, max_num_threads)(gen);
    pred_succ_occ<decltype(read)> index(read, n, sigma, occs_per_block, threads);
    verify_pred_succ_queries<uint64_t>(index, input, sym_distrib, num_queries,
        "pred_succ_occ (" + std::to_string(threads) + " threads)", gen);

    pred_succ_occ<decltype(read)> index_reloaded;
    serialize_and_load(index, index_reloaded);
    index_reloaded.set_read(read);
    verify_pred_succ_queries<uint64_t>(index_reloaded, input, sym_distrib, num_queries, "pred_succ_occ (reloaded)", gen);
}

/**
 * @brief builds symbol_groups over a random input, with a random group size and a random set of stored symbols, and
 *        verifies next and previous against a scan, before and after serialize/load
 * @tparam gen_t random number generator type
 * @param gen random number generator
 * @param max_size maximum length of the generated input
 */
template <typename gen_t>
static void verify_symbol_groups(gen_t& gen, uint64_t max_size)
{
    constexpr uint64_t max_sigma = symbol_groups_max_sigma;

    const uint64_t n = random_log_uniform_size(1, max_size, gen);
    const uint64_t sigma = 1 + gen() % max_sigma;
    const uint64_t group_size = uint64_t { 1 } << (gen() % (max_log_group_size + 1));
    // [0..n-1] the input, symbols in [0..sigma-1]
    std::vector<uint64_t> input(n);
    for (uint64_t& c : input) c = gen() % sigma;
    // [0..sigma-1] stored[c] <=> c is among the symbols passed to symbol_groups to store
    std::vector<bool> stored(sigma);
    for (uint64_t c = 0; c < sigma; c++) stored[c] = gen() % 2 == 0;
    const uint64_t threads = std::uniform_int_distribution<uint64_t>(1, max_num_threads)(gen);
    auto read = [&input](uint64_t i) { return input[i]; };
    symbol_groups groups(read, n, group_size, stored, threads);
    const uint64_t num_groups = div_ceil<uint64_t>(n, group_size);
    // [0..sigma-1][0..num_groups-1] contains[c][h] <=> the group h contains the symbol c
    std::vector<std::vector<bool>> contains(sigma, std::vector<bool>(num_groups, false));
    for (uint64_t i = 0; i < n; i++) contains[input[i]][i / group_size] = true;

    auto check = [&](const symbol_groups& sg, const std::string& name) {
        for (uint64_t c = 0; c < sigma; c++) {
            const bool occurs = std::find(contains[c].begin(), contains[c].end(), true) != contains[c].end();
            EXPECT_EQ(sg.stores(c), stored[c] && occurs) << name << " stores(" << c << ")";
            if (!sg.stores(c)) continue;

            for (uint64_t h = 0; h < num_groups; h++) {
                uint64_t next = symbol_groups::none;
                for (uint64_t g = h; g < num_groups; g++) if (contains[c][g]) { next = g; break; }
                uint64_t previous = symbol_groups::none;
                for (uint64_t g = h + 1; g > 0; g--) if (contains[c][g - 1]) { previous = g - 1; break; }
                EXPECT_EQ(sg.next(c, h), next) << name << " next(" << c << "," << h << ")";
                EXPECT_EQ(sg.previous(c, h), previous) << name << " previous(" << c << "," << h << ")";
            }
        }
    };

    check(groups, "symbol_groups (" + std::to_string(threads) + " threads)");
    symbol_groups reloaded;
    serialize_and_load(groups, reloaded);
    check(reloaded, "symbol_groups (reloaded)");
}

/**
 * @brief builds a pred_succ_block over a random byte input and verifies its succ and pred queries against a scan,
 *        before and after serialize/load
 * @tparam gen_t random number generator type
 * @param gen random number generator
 * @param max_size maximum length of the generated input
 */
template <typename gen_t>
static void verify_pred_succ_block(gen_t& gen, uint64_t max_size)
{
    constexpr int32_t max_sigma = psb_max_sigma;
    constexpr uint64_t max_num_queries = psb_max_num_queries;

    std::uniform_int_distribution<int32_t> sigma_distrib(1, max_sigma);
    uint64_t n = random_log_uniform_size(1, max_size, gen);
    uint64_t sigma = sigma_distrib(gen);
    std::uniform_int_distribution<int32_t> sym_distrib(0, sigma - 1);

    // [0..n-1] the input, symbols in [0..sigma-1]
    no_init_vector<uint8_t> input;
    resize_no_init(input, n);
    for (uint64_t i = 0; i < n; i++)
        input[i] = (uint8_t) sym_distrib(gen);

    const uint64_t num_queries = std::min<uint64_t>(n, max_num_queries);
    auto read = [&input](uint64_t i) -> uint8_t { return input[i]; };
    pred_succ_block<decltype(read)> index(read, n, sigma);
    verify_pred_succ_queries<uint8_t>(index, input, sym_distrib, num_queries, "pred_succ_block", gen);

    pred_succ_block<decltype(read)> index_reloaded;
    serialize_and_load(index, index_reloaded);
    index_reloaded.set_read(read);
    verify_pred_succ_queries<uint8_t>(index_reloaded, input, sym_distrib, num_queries, "pred_succ_block (reloaded)", gen);
}

TEST(test_data_structures, parse_double_ignores_locale)
{
    const std::string before = std::setlocale(LC_NUMERIC, nullptr);
    const bool comma_locale = std::setlocale(LC_NUMERIC, "de_DE.UTF-8") != nullptr;
    double value = 0;
    EXPECT_TRUE(parse_double("0.1", value));
    EXPECT_EQ(value, 0.1);
    EXPECT_TRUE(parse_double("-2.5e3", value));
    EXPECT_EQ(value, -2500.0);
    EXPECT_FALSE(parse_double("0,1", value));
    EXPECT_FALSE(parse_double("1.5x", value));
    EXPECT_FALSE(parse_double("", value));
    std::setlocale(LC_NUMERIC, before.c_str());
    if (!comma_locale) GTEST_SKIP() << "the locale de_DE.UTF-8 is not installed; only the C locale was tested";
}

TEST(test_data_structures, plain_bit_vector)
{
    run_fuzz("plain-bit-vector", {
        { "verify", [](uint64_t) { verify_bit_vector<plain_bit_vector>(random_bit_vector(gen, max_input_size)); }, true },
    }, fuzz_iterations(num_iterations_plain_bit_vector), &max_input_size, min_fuzz_input_size);
}

TEST(test_data_structures, elias_fano_bit_vector)
{
    run_fuzz("elias-fano-bit-vector", {
        { "verify", [](uint64_t) { verify_bit_vector<elias_fano_bit_vector>(random_bit_vector(gen, max_input_size)); }, true },
    }, fuzz_iterations(num_iterations_elias_fano_bit_vector), &max_input_size, min_fuzz_input_size);
}

TEST(test_data_structures, rank_select_set)
{
    run_fuzz("rank-select-set", {
        { "verify", [](uint64_t) { verify_rank_select_set(random_bit_vector(gen, max_input_size)); }, true },
    }, fuzz_iterations(), &max_input_size, min_fuzz_input_size);
}

TEST(test_data_structures, uint_bytes_array)
{
    run_fuzz("uint-bytes-array", {
        { "single values", [](uint64_t) { verify_uint_bytes_array<1>(max_input_size); }, true },
        { "pairs", [](uint64_t) { verify_uint_bytes_array<2>(max_input_size); }, true },
    }, fuzz_iterations(), &max_input_size, min_fuzz_input_size);
}

TEST(test_data_structures, blocked_positions)
{
    run_fuzz("blocked-positions", {
        { "verify", [](uint64_t) { verify_blocked_positions(max_input_size); }, true },
    }, fuzz_iterations(), &max_input_size, min_fuzz_input_size);
}

TEST(test_data_structures, interleaved_bit_aligned_vectors)
{
    run_fuzz("interleaved-bit-aligned-vectors", {
        { "verify", [](uint64_t) { verify_interleaved_bit_aligned_vectors(gen,
            uint64_t(size_factor_interleaved_bit_aligned_vectors * double(max_input_size))); }, true },
    }, fuzz_iterations(num_iterations_interleaved_bit_aligned_vectors), &max_input_size, min_fuzz_input_size);
}

TEST(test_data_structures, wavelet_tree_prefix_sum)
{
    run_fuzz("wavelet-tree-prefix-sum", {
        { "verify", [](uint64_t) { verify_wavelet_tree_prefix_sum(gen,
            uint64_t(size_factor_wavelet_tree_prefix_sum * double(max_input_size))); }, true },
    }, fuzz_iterations(num_iterations_wavelet_tree_prefix_sum), &max_input_size, min_fuzz_input_size);
}

TEST(test_data_structures, pred_succ_occ)
{
    run_fuzz("pred-succ-occ", {
        { "verify", [](uint64_t) { verify_pred_succ_occ(gen,
            uint64_t(size_factor_pred_succ_occ * double(max_input_size))); }, true },
    }, fuzz_iterations(num_iterations_pred_succ_occ), &max_input_size, min_fuzz_input_size);
}

TEST(test_data_structures, symbol_groups)
{
    run_fuzz("symbol-groups", {
        { "verify", [](uint64_t) { verify_symbol_groups(gen,
            uint64_t(size_factor_symbol_groups * double(max_input_size))); }, true },
    }, fuzz_iterations(num_iterations_symbol_groups), &max_input_size, min_fuzz_input_size);
}

TEST(test_data_structures, pred_succ_block)
{
    run_fuzz("pred-succ-block", {
        { "verify", [](uint64_t) { verify_pred_succ_block(gen,
            uint64_t(size_factor_pred_succ_block * double(max_input_size))); }, true },
    }, fuzz_iterations(num_iterations_pred_succ_block), &max_input_size, min_fuzz_input_size);
}

// packed_bit_fields stores bit-fields at the bit positions of the bit-fields of a packed struct of GCC and clang on an
// Itanium-ABI target (the Microsoft ABI lays out bit-fields differently, so there is nothing to compare with there)
#ifndef _MSC_VER
namespace {

// widths of the fields of the struct below: those of move_rb's sample_info_pack_t (two fields sharing a byte, three
// 48-bit values)
constexpr uint64_t sample_dir_bits = 2, sample_kind_bits = 8 - sample_dir_bits, sample_value_bits = 8 * uint48_t::width;

// a packed struct with fields that share a byte and fields that cross 64-bit words
struct __attribute__((packed)) gcc_sample_layout_t {
    uint64_t d : sample_dir_bits, t : sample_kind_bits, o : sample_value_bits, i : sample_value_bits,
        x : sample_value_bits;
};

} // namespace

TEST(test_data_structures, packed_bit_fields_layout)
{
    using d_field = next_bit_field_t<void, sample_dir_bits>;
    using t_field = next_bit_field_t<d_field, sample_kind_bits>;
    using o_field = next_bit_field_t<t_field, sample_value_bits>;
    using i_field = next_bit_field_t<o_field, sample_value_bits>;
    using x_field = next_bit_field_t<i_field, sample_value_bits>;

    static_assert(sizeof(packed_bit_fields<bytes_of_bit_fields<x_field>>) == sizeof(gcc_sample_layout_t));

    // a random value of a field of width bits
    auto random_value = [](uint64_t width) {
        return std::uniform_int_distribution<uint64_t>(0, low_bits_mask(width))(gen);
    };

    for (uint64_t iteration = 0; iteration < num_iterations_packed_bit_fields; iteration++) {
        gcc_sample_layout_t gcc_sample;
        std::memset(&gcc_sample, 0, sizeof(gcc_sample));
        packed_bit_fields<bytes_of_bit_fields<x_field>> sample {};
        const uint64_t d = random_value(d_field::width), t = random_value(t_field::width);
        const uint64_t o = random_value(o_field::width), i = random_value(i_field::width);
        const uint64_t x = random_value(x_field::width);
        gcc_sample.d = d; gcc_sample.t = t; gcc_sample.o = o; gcc_sample.i = i; gcc_sample.x = x;
        sample.set<d_field>(d); sample.set<t_field>(t); sample.set<o_field>(o); sample.set<i_field>(i);
        sample.set<x_field>(x);
        ASSERT_EQ(std::memcmp(&gcc_sample, sample.bytes, sizeof(gcc_sample)), 0);
        EXPECT_EQ(sample.get<d_field>(), d);
        EXPECT_EQ(sample.get<t_field>(), t);
        EXPECT_EQ(sample.get<o_field>(), o);
        EXPECT_EQ(sample.get<i_field>(), i);
        EXPECT_EQ(sample.get<x_field>(), x);

        // setting a field again changes only its bits
        const uint64_t o_new = random_value(o_field::width);
        gcc_sample.o = o_new;
        sample.set<o_field>(o_new);
        ASSERT_EQ(std::memcmp(&gcc_sample, sample.bytes, sizeof(gcc_sample)), 0);
    }
}
#endif
