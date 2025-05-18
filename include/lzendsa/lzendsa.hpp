/**
 * part of LukasNalbach/Move-r
 *
 * MIT License
 *
 * Copyright (c) Jan Zumbrink, Lukas Nalbach
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

#include "lzendsa_construction.hpp"
#include "lzendsa_encoding.hpp"

#include <algorithms/build_sa_and_bwt.hpp>
#include <algorithms/sparse_sa_bin_search.cpp>

template <typename int_t = int32_t>
class lzendsa {

protected:
    // default overall size of the SA-samples relative to the index size
    static constexpr double default_relative_sampling_size = 0.1;

    lzendsa_encoding lzendsa_enc; // LZ-End lzendsa_encoding
    int_t d; // sampling parameter
    int_t last_sa; // last SA-value (SA[n - 1])

    // evenly-spaced SA-samples (every d-th position is sampled, if d > 0)
    interleaved_bit_aligned_vectors<uint64_t, 1> sa_samples;

    const std::string* input;

public:
    lzendsa() = default;

    // call when the suffix array of the input text is not yet calculated
    lzendsa(std::string& input, int_t d = -1, int_t h = 8192, bool use_bigbwt = false, bool log = false) : d(d)
    {
        auto time_start = now();
        auto time = time_start;

        // compute SA and BWT
        auto [sa, bwt] = build_sa_and_bwt<int_t>(input, false, use_bigbwt, log);
        uint64_t n = sa.size();
        last_sa = sa[n - 1];

        // build DSA from SA
        if (log) time = now();
        if (log) std::cout << "building DSA from SA" << std::flush;
        std::vector<int_t> dsa;
        no_init_resize(dsa, n);
        dsa[0] = sa[0];

        for (uint64_t i = 1; i < n; i++) {
            dsa[i] = sa[i] - sa[i - 1];
        }

        if (log) time = log_runtime(time);

        // construct lzend parsing of DSA
        auto lzend_phrases = construct_lzend_of_reverse<int_t>(dsa, h, log);

        // discard DSA
        dsa.clear();
        dsa.shrink_to_fit();

        // lzendsa_encode the lzend parsing of DSA
        lzendsa_enc = lzendsa_encoding(lzend_phrases, n, h);

        if (d <= 0) {
            uint64_t target_sampling_size_in_bits = (lzendsa_enc.size_in_bytes() * 8) * default_relative_sampling_size;
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

    // return a reference to the lzendsa_encoding
    const lzendsa_encoding& sa_encoding() const
    {
        return lzendsa_enc;
    }

    uint64_t num_samples() const
    {
        return sa_samples.size() + 1;
    }

    // return the index-th suffix array sample
    int_t sample(uint64_t index) const
    {
        return index == sa_samples.size() ? last_sa : sa_samples[index];
    }

    // return the position in SA of the index-th suffix array sample
    int_t sample_pos(uint64_t index) const
    {
        return index == sa_samples.size() ? input_size() - 1 : (index * d);
    }

    uint64_t input_size() const
    {
        return lzendsa_enc.input_size();
    }

    uint64_t num_phrases() const
    {
        return lzendsa_enc.num_phrases();
    }

    int_t delta() const
    {
        return d;
    }

    std::tuple<uint64_t, uint64_t> count(const std::string& pattern) const
    {
        auto [beg, end, result] = binary_sa_search_and_extract<int_t>(*input, pattern, num_samples(),
            [&](uint64_t i){return sample(i);}, [&](uint64_t i){return sample_pos(i);},
            [&](uint64_t b, uint64_t e, uint64_t sa_b, uint64_t sa_e, auto report){
                lzendsa_enc.template extract(b, e, sa_e, report);}, false);
        
        return {beg, end};
    }

    template <typename out_t>
    std::vector<out_t> locate(const std::string& pattern) const
    {
        auto [beg, end, result] = binary_sa_search_and_extract<out_t>(*input, pattern, num_samples(),
            [&](uint64_t i){return sample(i);}, [&](uint64_t i){return sample_pos(i);},
            [&](uint64_t b, uint64_t e, uint64_t sa_b, uint64_t sa_e, auto report){
                lzendsa_enc.template extract(b, e, sa_e, report);}, true);
        
        return result;
    }

    // returns multiple consecutive suffix array values
    template <typename out_t>
    std::vector<out_t> sa_values(int64_t beg, int64_t end) const
    {
        int64_t len = end - beg + 1;
        int64_t smpl_idx = beg / d;
        int64_t smpl_pos = smpl_idx * d;
        std::vector<out_t> rng;
        
        if (smpl_pos < beg) {
            int64_t next_smpl_pos = sample_pos(smpl_idx + 1);
    
            if (next_smpl_pos <= end || next_smpl_pos - end < beg - smpl_pos) {
                smpl_idx++;
                smpl_pos = next_smpl_pos;
            }
        }
    
        if (smpl_pos <= beg) {
            rng = lzendsa_enc.template extract_deltas<out_t>(smpl_pos, end);
            int64_t val = sample(smpl_idx);
            int64_t dist = beg - smpl_pos;
    
            for (int64_t i = 1; i <= dist; i++) {
                val += rng[i];
            }
    
            int64_t last_val;
    
            for (int64_t i = 0; i < len; i++) {
                last_val = val;
                val += rng[i + dist + 1];
                rng[i] = last_val;
            }
    
            rng.resize(len);
        } else if (smpl_pos <= end) {
            rng = lzendsa_enc.template extract_deltas<out_t>(beg, end);
            int64_t smpl = sample(smpl_idx);
            int64_t val = smpl;
    
            for (int64_t i = smpl_pos - beg + 1; i < len; i++) {
                val += rng[i];
                rng[i] = val;
            }
    
            val = smpl;
            int64_t last_val;
    
            for (int64_t i = smpl_pos - beg; i >= 0; i--) {
                last_val = val;
                val -= rng[i];
                rng[i] = last_val;
            }
        } else {
            rng = lzendsa_enc.template extract_deltas<out_t>(beg, smpl_pos);
            int64_t val = sample(smpl_idx);
    
            for (int64_t i = rng.size() - 1; i >= len; i--) {
                val -= rng[i];
            }
    
            int64_t last_val;
    
            for (int64_t i = len - 1; i >= 1; i--) {
                last_val = val;
                val -= rng[i];
                rng[i] = last_val;
            }
    
            rng[0] = val;
            rng.resize(len);
        }
    
        return rng;
    }

    // return the size of the whole data structure in bytes
    uint64_t size_in_bytes() const
    {
        return sizeof(this) + lzendsa_enc.size_in_bytes() + sa_samples.size_in_bytes();
    }

    // load the lzendsa construction from a stream
    void load(std::istream& in)
    {
        in.read((char*) &d, sizeof(d));
        in.read((char*) &last_sa, sizeof(last_sa));
        lzendsa_enc.load(in);
        sa_samples.load(in);
    }

    // serialize the lzendsa construction to the output stream
    void serialize(std::ostream& out)
    {
        out.write((char*) &d, sizeof(d));
        out.write((char*) &last_sa, sizeof(last_sa));
        lzendsa_enc.serialize(out);
        sa_samples.serialize(out);
    }
};