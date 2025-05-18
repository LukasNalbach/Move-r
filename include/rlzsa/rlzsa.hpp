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
#include <algorithms/sparse_sa_bin_search.cpp>
#include <data_structures/interleaved_bit_aligned_vectors.hpp>

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

public:
    rlzsa() = default;

    // call when the suffix array of the input text is not yet calculated
    rlzsa(std::string& input, int64_t d = -1, bool use_bigbwt = false, bool log = false)
    {
        auto time_start = now();
        auto time = time_start;
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

        bwt.clear();
        bwt.shrink_to_fit();

        // construct differential suffix array (DSA)
        if (log) time = now();
        if (log) std::cout << "building Differential Suffix Array (DSA)" << std::flush;
        std::vector<int_t> dsa;
        no_init_resize(dsa, n);
        dsa[0] = sa[0];

        for (uint64_t i = 1; i < n; i++) {
            dsa[i] = sa[i] - sa[i - 1];
        }

        if (log) time = log_runtime(time);

        // construct rlzsa
        if (log) std::cout << "building rlzsa:" << std::endl;
        uint64_t reference_size = std::min<uint64_t>(n / 3, 5.2 * r);
        rlzsa_enc = rlzsa_encoding<int_t>(sa, std::move(dsa), reference_size, log);

        if (d <= 0) {
            uint64_t target_sampling_size_in_bits = (rlzsa_enc.size_in_bytes() * 8) * default_relative_sampling_size;
            d = n / (target_sampling_size_in_bits / std::bit_width(uint64_t{n}));
            this->d = d;
        }

        // construct the SA sampling
        if (log) std::cout << "building SA-Samples, d = " << d << std::flush;
        
        time = now();
        uint64_t num_samples = n / d;
        sa_samples = interleaved_bit_aligned_vectors<uint64_t, 1>({ std::bit_width(uint64_t{n}) });
        sa_samples.resize_no_init(num_samples);

        for (uint64_t i = 0; i < num_samples; i++) {
            sa_samples.set<0>(i, sa[i * d]);
        }

        if (log) time = log_runtime(time);
    }

    void set_input(const std::string& str)
    {
        input = &str;
    }

    // return a referrlzsa_ence to the rlzsa_encoding
    const rlzsa_encoding<int_t>& sa_encoding() const
    {
        return rlzsa_enc;
    }

    // return the index-th suffix array sample
    uint64_t sample(uint64_t index) const
    {
        return index == sa_samples.size() ? last_sa : sa_samples[index];
    }

    // return the position in SA of the index-th suffix array sample
    uint64_t sample_pos(uint64_t index) const
    {
        return index == sa_samples.size() ? input_size() - 1 : (index * d);
    }

    uint64_t input_size() const
    {
        return rlzsa_enc.input_size();
    }

    uint64_t num_phrases() const
    {
        return rlzsa_enc.num_phrases();
    }

    int64_t delta() const
    {
        return d;
    }

    std::tuple<uint64_t, uint64_t> count(const std::string& pattern) const
    {
        auto [beg, end, result] = binary_sa_search_and_extract<int_t>(*input, pattern, num_samples(),
            [&](uint64_t i){return sample(i);}, [&](uint64_t i){return sample_pos(i);},
            [&](uint64_t b, uint64_t e, uint64_t sa_b, uint64_t sa_e, auto report){
                rlzsa_enc.template extract<int_t>(b, e, report, sa_b);}, false);
        
        return {beg, end};
    }

    template <typename out_t>
    std::vector<out_t> locate(const std::string& pattern) const
    {
        auto [beg, end, result] = binary_sa_search_and_extract<out_t>(*input, pattern, num_samples(),
            [&](uint64_t i){return sample(i);}, [&](uint64_t i){return sample_pos(i);},
            [&](uint64_t b, uint64_t e, uint64_t sa_b, uint64_t sa_e, auto report){
                rlzsa_enc.template extract<out_t>(b, e, report, sa_b);}, true);
        
        return result;
    }

    // returns multiple consecutive suffix array values
    template <typename out_t>
    std::vector<out_t> sa_values(uint64_t beg, uint64_t end) const
    {
        uint64_t smpl_idx = beg / d;
        uint64_t smpl_pos = smpl_idx * d;
        uint64_t smpl = sa_samples[smpl_idx];

        std::vector<out_t> result;
        no_init_resize(result, end - beg + 1);

        rlzsa_enc.template extract<out_t>(smpl_pos, end,
            [&](uint64_t pos, uint64_t val){if (pos >= beg) result[pos - beg] = val;}, smpl);

        return result;
    }

    uint64_t num_samples() const
    {
        return sa_samples.size() + 1;
    }

    // return the size of the whole data structure in bytes
    size_t size_in_bytes() const
    {
        return sizeof(this) + rlzsa_enc.size_in_bytes() + sa_samples.size_in_bytes();
    }

    // load the rlzsa construction from a stream
    void load(std::istream& in)
    {
        in.read((char*) &d, sizeof(d));
        in.read((char*) &last_sa, sizeof(last_sa));
        rlzsa_enc.load(in);
        sa_samples.load(in);
    }

    // serialize the rlzsa construction to the output stream
    void serialize(std::ostream& out)
    {
        out.write((char*) &d, sizeof(d));
        out.write((char*) &last_sa, sizeof(last_sa));
        rlzsa_enc.serialize(out);
        sa_samples.serialize(out);
    }
};