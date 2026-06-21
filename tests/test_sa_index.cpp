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

#include <gtest/gtest.h>
#include <lzendsa/lzendsa.hpp>
#include <lzendsa/r_index_lzendsa.hpp>
#include <rlzsa/rlzsa.hpp>
#include <rlzsa/r_index_rlzsa.hpp>
#include "test_index_common.hpp"

// adapters for the generic index fuzz-test (test_index_instance) that test the standalone SA-based indexes.
// each concrete index supplies a small config type:
//   - self SA indexes (lzendsa, rlzsa):   { using index_t; static index_t make(std::string&); static std::string name(); }
//   - r-index variants (r_index_lzendsa, r_index_rlzsa): additionally { static std::vector<int_t> do_locate(index_t&, std::string&); }

/**
 * @brief returns the number of pattern (count/locate) queries to run for an input of the given size
 * @param input_size the size of the input
 * @return the number of queries
 */
static uint64_t num_sa_index_queries(uint64_t input_size)
{
    return std::min<uint64_t>(1000, std::max<uint64_t>(100, input_size / 10));
}

/**
 * @brief index adapter for a standalone SA self-index (lzendsa, rlzsa); besides count/locate it verifies the
 * random-access (sa_values) queries against a reference suffix array
 * @tparam int_t signed integer type of the suffix array entries
 * @tparam config_t the index config (provides index_t, make() and name())
 */
template <typename int_t, typename config_t>
struct sa_self_index_adapter {
    using pos_t = int_t;
    using inp_t = std::string;
    using index_t = typename config_t::index_t;

    static constexpr bool serializable = true;

    static std::string name() { return config_t::name(); }
    static char min_sym() { return sa_index_min_uchar; }
    static char max_sym() { return sa_index_max_uchar; }
    static uint64_t max_pattern_length(uint64_t) { return 100; }
    static uint64_t num_queries(uint64_t n) { return num_sa_index_queries(n); }

    template <typename gen_t>
    static index_t build(std::string& input, gen_t&, uint16_t) { return config_t::make(input); }

    static void after_build(index_t& index, const std::string& input) { index.set_input(input); }

    static uint64_t count(index_t& index, std::string& pattern)
    {
        auto [beg, end] = index.count(pattern);
        return end >= beg ? end - beg + 1 : 0;
    }

    static std::vector<pos_t> locate(index_t& index, std::string& pattern)
    {
        return index.template locate<int_t>(pattern);
    }

    template <typename gen_t>
    static void verify_pattern(index_t& index, std::string& pattern, const std::vector<pos_t>& correct, gen_t&)
    {
        verify_count_locate_basic<sa_self_index_adapter>(index, pattern, correct);
    }

    template <typename gen_t>
    static void verify_build(index_t& index, const std::string& input, gen_t& gen, uint16_t num_threads)
    {
        std::vector<int64_t> sa = compute_reference_sa<char, std::string>(input, input.size(), num_threads);
        uint64_t n = sa.size();
        std::uniform_int_distribution<uint64_t> sa_pos_distrib(0, n - 1);

        for (uint64_t q = 0; q < num_sa_index_queries(n); q++) {
            uint64_t beg = sa_pos_distrib(gen);
            uint64_t end = sa_pos_distrib(gen);
            if (beg > end) std::swap(beg, end);
            std::vector<int_t> values = index.template sa_values<int_t>(beg, end);
            ASSERT_EQ(values.size(), end - beg + 1) << name() << ": sa_values size";

            for (uint64_t i = 0; i < values.size(); i++)
                EXPECT_EQ((int64_t) values[i], sa[beg + i]) << name() << ": sa_values[" << (beg + i) << "]";
        }

        std::vector<int_t> all_values = index.template sa_values<int_t>(0, n - 1);
        ASSERT_EQ(all_values.size(), n) << name() << ": sa_values(full) size";
        for (uint64_t i = 0; i < n; i++)
            EXPECT_EQ((int64_t) all_values[i], sa[i]) << name() << ": sa_values(full)[" << i << "]";
    }
};

/**
 * @brief index adapter for an r-index variant (r_index_lzendsa, r_index_rlzsa); it verifies count/locate
 * (the r-indexes have no random-access or full-text queries to check at build time)
 * @tparam int_t signed integer type of the suffix array entries
 * @tparam config_t the index config (provides index_t, make(), do_locate() and name())
 */
template <typename int_t, typename config_t>
struct sa_r_index_adapter {
    using pos_t = int_t;
    using inp_t = std::string;
    using index_t = typename config_t::index_t;

    static constexpr bool serializable = true;

    static std::string name() { return config_t::name(); }
    static char min_sym() { return sa_index_min_uchar; }
    static char max_sym() { return sa_index_max_uchar; }
    static uint64_t max_pattern_length(uint64_t) { return 100; }
    static uint64_t num_queries(uint64_t n) { return num_sa_index_queries(n); }

    template <typename gen_t>
    static index_t build(std::string& input, gen_t&, uint16_t) { return config_t::make(input); }

    static void after_build(index_t&, const std::string&) { }

    static uint64_t count(index_t& index, std::string& pattern)
    {
        auto [beg, end] = index.count(pattern);
        return end >= beg ? end - beg + 1 : 0;
    }

    static std::vector<pos_t> locate(index_t& index, std::string& pattern)
    {
        return config_t::do_locate(index, pattern);
    }

    template <typename gen_t>
    static void verify_pattern(index_t& index, std::string& pattern, const std::vector<pos_t>& correct, gen_t&)
    {
        verify_count_locate_basic<sa_r_index_adapter>(index, pattern, correct);
    }

    template <typename gen_t>
    static void verify_build(index_t&, const std::string&, gen_t&, uint16_t) { }
};

// ############################# lzendsa #############################

// config for the standalone lzendsa self-index
template <typename int_t>
struct lzendsa_config {
    using index_t = lzendsa<int_t>;
    static index_t make(std::string& input) { return lzendsa<int_t>(input, -1, 8192, false, false); }
    static std::string name() { return "lzendsa"; }
};

// config for the r-index-lzendsa
template <typename int_t>
struct r_index_lzendsa_config {
    using index_t = r_index_lzendsa<int_t>;
    static index_t make(std::string& input) { return r_index_lzendsa<int_t>(input, 8192, false, false); }
    static std::vector<int_t> do_locate(index_t& index, std::string& pattern) { return index.locate(pattern); }
    static std::string name() { return "r_index_lzendsa"; }
};

// ############################# rlzsa #############################

// config for the standalone rlzsa self-index
template <typename int_t>
struct rlzsa_config {
    using index_t = rlzsa<int_t>;
    static index_t make(std::string& input) { return rlzsa<int_t>(input, -1, false, false); }
    static std::string name() { return "rlzsa"; }
};

// config for the r-index-rlzsa
template <typename int_t>
struct r_index_rlzsa_config {
    using index_t = r_index_rlzsa<int_t>;
    static index_t make(std::string& input) { return r_index_rlzsa<int_t>(input, false, false, false); }
    static std::vector<int_t> do_locate(index_t& index, std::string& pattern) { return index.template locate<int_t>(pattern); }
    static std::string name() { return "r_index_rlzsa"; }
};

std::random_device rd;
std::mt19937 gen(rd());
uint16_t max_num_threads = omp_get_max_threads();
uint64_t min_input_size = 1;
uint64_t max_input_size = 100000;

std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);
std::uniform_int_distribution<int> sa_index_distrib(0, 3);

// all standalone SA-based indexes (lzendsa & r-index-lzendsa & rlzsa & r-index-rlzsa), rotating over the
// index type and the 32-/64-bit suffix-array integer type
TEST(test_sa_index, fuzzy_test)
{
    auto start_time = now();

    while (time_diff_min(start_time, now()) < 60) {
        bool use_64_bit = prob_distrib(gen) < 0.5;

        auto run = [&]<template <typename> typename config_t, template <typename, typename> typename adapter_t>() {
            if (use_64_bit) test_index_instance<adapter_t<int64_t, config_t<int64_t>>>(gen, max_num_threads, min_input_size, max_input_size);
            else            test_index_instance<adapter_t<int32_t, config_t<int32_t>>>(gen, max_num_threads, min_input_size, max_input_size);
        };

        switch (sa_index_distrib(gen)) {
            case 0: run.template operator()<lzendsa_config,         sa_self_index_adapter>(); break;
            case 1: run.template operator()<r_index_lzendsa_config, sa_r_index_adapter>();    break;
            case 2: run.template operator()<rlzsa_config,           sa_self_index_adapter>(); break;
            case 3: run.template operator()<r_index_rlzsa_config,   sa_r_index_adapter>();    break;
        }
    }
}
