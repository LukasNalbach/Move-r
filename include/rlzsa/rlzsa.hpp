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

#pragma once

#include <bit>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <ostream>
#include <type_traits>
#include <vector>

#include "rlzsa_encoding.hpp"
#include <algorithms/build_sa_and_bwt.hpp>
#include <algorithms/sparse_sa_bin_search.tpp>
#include <data_structures/interleaved_bit_aligned_vectors.hpp>

/**
 * @brief a standalone self-index based on a relative-Lempel-Ziv encoding of the differential suffix array, with
 *        evenly-spaced SA-samples (unoptimized; the move-r optimized variant is rlzsa_opt in include/rlzsa_opt/)
 * @tparam int_t signed integer type of the suffix array entries
 */
template <typename int_t = int32_t>
class rlzsa {

protected:
    // default overall size of the SA-samples relative to the index size
    static constexpr double default_relative_sampling_size = 0.1;

    rlzsa_encoding<int_t> rlzsa_enc; // rlzsa encoding
    int64_t d; // sampling parameter
    int_t last_sa; // last SA-value (SA[n - 1])

    // evenly-spaced SA-samples (every d-th position is sampled)
    interleaved_bit_aligned_vectors<uint64_t, 1> sa_samples;

    const std::string* input;

    /**
     * @brief builds the evenly-spaced SA-sampling
     * @param sa the suffix array
     * @param d the sampling parameter (every d-th position is sampled); if d <= 0, it is chosen so that the
     *          samples take about default_relative_sampling_size of the index size
     * @param log whether to print log messages
     */
    void build_sampling(const std::vector<int_t>& sa, int64_t d, bool log)
    {
        auto time = now();
        assert(d != 0);
        uint64_t n = rlzsa_enc.input_size();

        if (d <= 0) {
            uint64_t target_sampling_size_in_bits = (rlzsa_enc.size_in_bytes() * 8) * default_relative_sampling_size;
            uint64_t divisor = std::max<uint64_t>(1, target_sampling_size_in_bits / std::bit_width(uint64_t{n}));
            d = std::max<uint64_t>(1, n / divisor);
        }

        this->d = d;

        // construct the SA sampling
        log_phase_start(log, time, "building SA-Samples, d = " + std::to_string(d));
        uint64_t num_samples = n / d;
        sa_samples = interleaved_bit_aligned_vectors<uint64_t, 1>({ bit_width(n) });
        sa_samples.resize_no_init(num_samples);

        for (uint64_t i = 0; i < num_samples; i++) {
            sa_samples.set<0>(i, sa[i * d]);
        }

        log_phase_end(log, time);
    }
    
public:
    rlzsa() = default;

    /**
     * @brief builds the index from the input text (computing the suffix array internally)
     * @param input the input text
     * @param d the SA-sampling parameter (if <= 0, it is chosen automatically)
     * @param use_bigbwt whether to use Big-BWT to build the suffix array and BWT
     * @param log whether to print log messages
     */
    rlzsa(std::string& input, int64_t d = -1, bool use_bigbwt = false, bool log = false)
    {
        auto time = now();
        assert(d != 0);

        // compute SA and BWT
        auto [sa, bwt] = build_sa_and_bwt<int_t>(input, true, use_bigbwt, log);
        uint64_t n = sa.size();
        last_sa = sa[n - 1];
        uint64_t r = 1;

        for (uint64_t i = 1; i < n; i++) {
            if (bwt[i] != bwt[i - 1]) {
                r++;
            }
        }

        // construct differential suffix array (DSA)
        log_phase_start(log, time, "building Differential Suffix Array (DSA)");
        std::vector<int_t> dsa;
        no_init_resize(dsa, n);
        dsa[0] = sa[0];

        for (uint64_t i = 1; i < n; i++) {
            dsa[i] = sa[i] - sa[i - 1];
        }

        log_phase_end(log, time);

        // construct rlzsa
        log_message(log, "building rlzsa:\n");
        uint64_t reference_size = std::min<uint64_t>(n / 3, 5.2 * r);
        rlzsa_enc = rlzsa_encoding<int_t>(sa, std::move(dsa), reference_size, false, log);
        log_phase_end(log, time);

        build_sampling(sa, d, log);
    }

    /**
     * @brief builds the index from a precomputed rlzsa encoding and suffix array
     * @param rlzsa_enc the rlzsa encoding of the differential suffix array
     * @param sa the suffix array
     * @param d the SA-sampling parameter (if <= 0, it is chosen automatically)
     * @param log whether to print log messages
     */
    rlzsa(const rlzsa_encoding<int_t>& rlzsa_enc, const std::vector<int_t>& sa, int64_t d = -1, bool log = false)
    {
        this->rlzsa_enc = rlzsa_enc;
        build_sampling(sa, d, log);
    }

    /**
     * @brief stores a pointer to the input text (needed for count/locate queries)
     * @param str the input text
     */
    void set_input(const std::string& str) { input = &str; }

    // returns a reference to the rlzsa encoding
    const rlzsa_encoding<int_t>& sa_encoding() const { return rlzsa_enc; }

    /**
     * @brief returns the index-th suffix array sample
     * @param index a sample index
     * @return the index-th suffix array sample
     */
    uint64_t sample(uint64_t index) const { return index == sa_samples.size() ? last_sa : sa_samples[index]; }

    /**
     * @brief returns the position in the suffix array of the index-th sample
     * @param index a sample index
     * @return the position in SA of the index-th sample
     */
    uint64_t sample_pos(uint64_t index) const { return index == sa_samples.size() ? input_size() - 1 : (index * d); }

    // returns the size of the input text
    uint64_t input_size() const { return rlzsa_enc.input_size(); }

    // returns the number of phrases in the rlzsa encoding
    uint64_t num_phrases() const { return rlzsa_enc.num_phrases(); }

    // returns the SA-sampling parameter
    int64_t delta() const { return d; }

    /**
     * @brief counts the occurrences of pattern in the input
     * @param pattern the pattern to count
     * @return the suffix array interval [beg, end] of the pattern (beg > end if it does not occur)
     */
    std::tuple<uint64_t, uint64_t> count(const std::string& pattern) const
    {
        auto [beg, end] = binary_sa_search_and_extract<int_t>(*input, pattern, num_samples(),
            [&](uint64_t i){return sample(i);}, [&](uint64_t i){return sample_pos(i);},
            [&](uint64_t b, uint64_t e, uint64_t sa_b, uint64_t /*sa_e*/, auto report){
                rlzsa_enc.template extract<int_t>(b, e, report, sa_b);}, [](int_t){}, false, false);
        
        return {beg, end};
    }

    /**
     * @brief locates the occurrences of pattern in the input
     * @tparam out_t output value type
     * @tparam report_fnc_t type of the report function
     * @param pattern the pattern to locate
     * @param report function that is called with every occurrence of pattern
     */
    template <typename out_t, typename report_fnc_t>
    void locate(const std::string& pattern, report_fnc_t report) const
    {
        binary_sa_search_and_extract<out_t>(*input, pattern, num_samples(),
            [&](uint64_t i){return sample(i);}, [&](uint64_t i){return sample_pos(i);},
            [&](uint64_t b, uint64_t e, uint64_t sa_b, uint64_t /*sa_e*/, auto report){
                rlzsa_enc.template extract<out_t>(b, e, report, sa_b);}, report, false, true);
    }

    /**
     * @brief locates the occurrences of pattern in the input
     * @tparam out_t output value type
     * @param pattern the pattern to locate
     * @return the occurrences of pattern in the input
     */
    template <typename out_t>
    std::vector<out_t> locate(const std::string& pattern) const
    {
        std::vector<out_t> result;
        locate<out_t>(pattern, [&](uint64_t val){result.emplace_back(val);});
        return result;
    }

    /**
     * @brief returns the suffix array values in the range [beg, end]
     * @tparam out_t output value type
     * @param beg left range limit
     * @param end right range limit
     * @return the suffix array values in the range [beg, end]
     */
    template <typename out_t>
    std::vector<out_t> sa_values(uint64_t beg, uint64_t end) const
    {
        uint64_t smpl_idx = std::min<uint64_t>(beg / d, sa_samples.size() - 1);
        uint64_t smpl_pos = smpl_idx * d;
        uint64_t smpl = sa_samples[smpl_idx];

        std::vector<out_t> result;
        result.reserve(end - beg + 1);
        uint64_t pos = smpl_pos;

        rlzsa_enc.template extract<out_t>(smpl_pos, end,
            [&](uint64_t val){if (pos++ >= beg) result.emplace_back(val);}, smpl);

        return result;
    }

    // returns the number of suffix array samples (including the last-SA sample)
    uint64_t num_samples() const { return sa_samples.size() + 1; }

    // returns the size of the whole data structure in bytes
    size_t size_in_bytes() const { return sizeof(this) + rlzsa_enc.size_in_bytes() + sa_samples.size_in_bytes(); }

    /**
     * @brief loads the index from an input stream
     * @param in input stream
     */
    void load(std::istream& in)
    {
        in.read((char*) &d, sizeof(d));
        in.read((char*) &last_sa, sizeof(last_sa));
        rlzsa_enc.load(in);
        sa_samples.load(in);
    }

    /**
     * @brief serializes the index to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out)
    {
        out.write((char*) &d, sizeof(d));
        out.write((char*) &last_sa, sizeof(last_sa));
        rlzsa_enc.serialize(out);
        sa_samples.serialize(out);
    }
};