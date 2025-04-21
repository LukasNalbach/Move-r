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
#include <misc/build_sa_and_bwt.hpp>

template <typename int_t = int32_t>
class lzendsa {

protected:
    lzendsa_encoding enc; // LZ-End encoding
    int_t d = 0; // sampling parameter

    // evenly-spaced SA-samples (every d-th position is sampled, if d != -1)
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

        // encode the lzend parsing of DSA
        enc = lzendsa_encoding(lzend_phrases, sa, n, d == -1);

        if (d != -1) {
            // construct the SA sampling
            if (log) std::cout << "building SA-Samples" << std::flush;

            uint64_t num_samples = n / d;
            sa_samples = interleaved_bit_aligned_vectors<uint64_t, 1>({ std::bit_width(uint64_t{n}) });
            sa_samples.resize_no_init(num_samples);

            for (uint64_t i = 0; i < num_samples; i++) {
                sa_samples.set<0>(i, sa[i * d]);
            }

            if (log) time = log_runtime(time);
        }
    }

    void set_input(const std::string& str) const
    {
        input = &str;
    }

    // return a reference to the encoding
    const lzendsa_encoding& encoding() const
    {
        return enc;
    }

    uint64_t num_samples() const
    {
        if (d == -1) return enc.num_phrases();
        return sa_samples.size();
    }

    // return the index-th suffix array sample
    int_t sample(uint64_t index) const
    {
        if (d == -1) return enc.sample(index);
        return sa_samples.get<0>(index);
    }

    // return the suffix array value at position index
    int_t operator[](uint64_t index) const
    {
        return enc[index];
    }

    uint64_t input_size() const
    {
        return enc.input_size();
    }

    uint64_t num_phrases() const
    {
        return enc.num_phrases();
    }

    int_t delta() const
    {
        return d;
    }

    // returns multiple consecutive suffix array values
    std::vector<int_t> sa_values(int64_t beg, int64_t end) const
    {
        // TODO

        /*
        
        #ifndef NDEBUG
        for (auto occ : result) {
            assert(0 <= occ && occ < input_size());
        }
        #endif

        return result;
        */
    }

    // return the size of the whole data structure in bytes
    size_t size_in_bytes() const
    {
        return sizeof(this) + enc.size_in_bytes() + sa_samples.size_in_bytes();
    }

    // load the lzendsa construction from a stream
    void load(std::istream& in)
    {
        enc.load(in);
        sa_samples.load(in);
    }

    // serialize the lzendsa construction to the output stream
    void serialize(std::ostream& out)
    {
        enc.serialize(out);
        sa_samples.serialize(out);
    }
};