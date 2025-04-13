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

#include <memory>
#include <utility>
#include <ostream>
#include <fstream>
#include <filesystem>

#include <misc/build_sa_and_bwt.hpp>
#include <custom_r_index/custom_r_index.hpp>
#include "rlzsa_encoding.hpp"

template<typename int_t = int32_t>
class r_index_rlzsa {

protected:
    rlzsa_encoding<int_t> enc;
    custom_r_index::index<_run_heads> r_index;

public:
    r_index_rlzsa() = default;

    r_index_rlzsa(std::string& input, bool use_bigbwt = false, bool use_rindex_samples = true, bool log = false)
    {
        auto time_start = now();
        auto time = time_start;

        // compute SA and BWT
        auto [sa, bwt] = build_sa_and_bwt<int_t>(input, true, use_bigbwt, log);
        uint64_t n = sa.size();

        // build r-index
        if (log) time = now();
        if (log) std::cout << "building r-index" << std::flush;
        r_index = custom_r_index::index<_run_heads>(bwt, sa, use_rindex_samples);
        bwt.clear();
        bwt.shrink_to_fit();
        if (log) time = log_runtime(time);

        // construct differential suffix array (DSA)
        if (log) std::cout << "building Differential Suffix Array (DSA)" << std::flush;
        std::vector<int_t> dsa;
        no_init_resize(dsa, n);
        dsa[0] = sa[0];

        for (uint64_t i = 1; i < n; i++) {
            dsa[i] = sa[i] - sa[i - 1];
        }

        if (log) time = log_runtime(time);
        
        // compute rlzsa
        if (log) std::cout << "building rlzsa:" << std::endl;
        uint64_t r = r_index.number_of_runs();
        uint64_t reference_size = std::min<uint64_t>(n / 3, 5.2 * r);
        enc = rlzsa_encoding<int_t>(sa, std::move(dsa), reference_size, !use_rindex_samples, log);
        if (log) std::cout << "r-index-rlzsa built in " << format_time(time_diff_ns(time_start, now())) << std::endl;
    }

    // return a reference to the encoding
    const rlzsa_encoding<int_t>& encoding() const
    {
        return enc;
    }

    template <typename out_t>
    std::vector<out_t> locate(const std::string &pattern) const
    {
        std::vector<out_t> Occ;

        if (r_index.has_sa_samples()) {
            auto [beg, end, first_value] = r_index.count_and_get_occ(pattern);
            if (end < beg) return {};
            Occ = enc.template extract<out_t>(beg, end, first_value);
        } else {
            auto [beg, end] = r_index.count(pattern);
            if (end < beg) return {};
            Occ = enc.template extract<out_t>(beg, end);
        }

        #ifndef NDEBUG
        for (auto occ : Occ) {
            assert(0 <= occ && occ < input_size());
        }
        #endif

        return Occ;
        
    }

    std::pair<uint64_t, uint64_t> count(std::string &pattern) const
    {
        return r_index.count(pattern);
    }

    uint64_t input_size() const
    {
        return enc.input_size();
    }

    uint64_t num_phrases() const
    {
        return enc.num_phrases();
    }

    uint64_t size_in_bytes() const
    {
        return sizeof(this) + enc.size_in_bytes() + r_index.size_in_bytes();
    }

    void serialize(std::ostream &out) const
    {
        enc.serialize(out);
        r_index.serialize(out);
    }

    void load(std::istream &in)
    {
        enc.load(in);
        r_index.load(in);
    }
};