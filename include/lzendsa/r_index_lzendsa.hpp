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

#include <memory>
#include <utility>

#include <misc/build_sa_and_bwt.hpp>
#include <custom_r_index/custom_r_index.hpp>
#include "lzendsa_encoding.hpp"

template<typename int_t = int32_t>
class r_index_lzendsa {

protected:
    lzendsa_encoding lzendsa_enc;
    custom_r_index::index<_run_ends> r_index;

public:
    r_index_lzendsa() = default;

    r_index_lzendsa(std::string& input, int_t h = -1, bool use_bigbwt = false, bool use_rindex_samples = true, bool log = false)
    {
        auto time_start = now();
        auto time = time_start;

        // compute SA and BWT
        auto [sa, bwt] = build_sa_and_bwt<int_t>(input, true, use_bigbwt, log);
        uint64_t n = sa.size();

        // build r-index
        if (log) time = now();
        if (log) std::cout << "building r-index" << std::flush;
        r_index = custom_r_index::index<_run_ends>(bwt, sa, use_rindex_samples);
        bwt.clear();
        bwt.shrink_to_fit();
        if (log) time = log_runtime(time);

        // construct differential suffix array (DSA)
        if (log) std::cout << "building Differential Suffix Array (DSA)" << std::flush;
        auto& dsa = sa;

        for (uint64_t i = n - 1; i > 0; i--) {
            dsa[i] = sa[i] - sa[i - 1];
        }

        if (log) time = log_runtime(time);
        
        // compute lzend parsing
        if (log) std::cout << "building LZ-End parsing:" << std::endl;
        std::vector<lzend_phr_t<int_t>> lzend_phrases = construct_lzend_of_reverse<int_t>(dsa, h, log);
        if (log) time = now();

        if (log) std::cout << "transforming DSA to SA" << std::flush;

        // make SA out of DSA
        for (uint64_t i = 1; i < n; i++) {
            sa[i] = dsa[i] + dsa[i - 1];
        }

        if (log) time = log_runtime(time);

        // construct lzend encoding
        if (log) std::cout << "Encoding LZ-End parsing" << std::flush;
        lzendsa_enc = lzendsa_encoding(lzend_phrases, sa, n, !use_rindex_samples);
        if (log) time = log_runtime(time);
        
        if (log) std::cout << "r-index-lzendsa built in " << format_time(time_diff_ns(time_start, now())) << std::endl;
    }

    // return a reference to the encoding
    const lzendsa_encoding& encoding() const
    {
        return lzendsa_enc;
    }

    std::vector<int_t> locate(const std::string &pattern) const
    {
        std::vector<int_t> result;

        if (pattern == "goldf") {
            int test = 3;
        }

        if (r_index.has_sa_samples()) {
            auto [beg, end, last_value] = r_index.count_and_get_occ(pattern);
            if (end < beg) return {};
            result = lzendsa_enc.template extract<int_t>(beg + 1, end);
            result.push_back(last_value);

            for (int_t i = result.size() - 2; i >= 0; i--) {
                result[i] = result[i + 1] - result[i];
            }
        } else {
            auto [beg, end] = r_index.count(pattern);
            if (end < beg) return {};
            
            int64_t length = end - beg + 1;
            int64_t phrase = lzendsa_enc.phrase_preceding_or_ending_at(end);
            int64_t sample_pos = lzendsa_enc.end_position(phrase);
            
            if (sample_pos < beg) {
                int64_t next_sample_pos = lzendsa_enc.end_position(phrase + 1);

                if (next_sample_pos <= end || next_sample_pos - end < beg - sample_pos) {
                    phrase++;
                    sample_pos = next_sample_pos;
                }
            }

            if (sample_pos < beg) {
                result = lzendsa_enc.extract<int_t>(sample_pos, end);
                int64_t val = lzendsa_enc.sample(phrase);
                int64_t dist = beg - sample_pos;

                for (int64_t i = 1; i <= dist; i++) {
                    val += result[i];
                }

                int64_t last_val;

                for (int64_t i = 0; i < length; i++) {
                    last_val = val;
                    val += result[i + dist + 1];
                    result[i] = last_val;
                }

                result.resize(length);
            } else if (sample_pos <= end) {
                result = lzendsa_enc.extract<int_t>(beg, end);
                int64_t sample = lzendsa_enc.sample(phrase);
                int64_t val = sample;

                for (int64_t i = sample_pos - beg + 1; i < length; i++) {
                    val += result[i];
                    result[i] = val;
                }

                val = sample;
                int64_t last_val;

                for (int64_t i = sample_pos - beg; i >= 0; i--) {
                    last_val = val;
                    val -= result[i];
                    result[i] = last_val;
                }
            } else {
                result = lzendsa_enc.extract<int_t>(beg, sample_pos);
                int64_t val = lzendsa_enc.sample(phrase);

                for (int64_t i = result.size() - 1; i >= length; i--) {
                    val -= result[i];
                }

                int64_t last_val;

                for (int64_t i = length - 1; i >= 1; i--) {
                    last_val = val;
                    val -= result[i];
                    result[i] = last_val;
                }

                result[0] = val;
                result.resize(length);
            }
        }

        #ifndef NDEBUG
        for (auto occ : result) {
            assert(0 <= occ && occ < input_size());
        }
        #endif
        
        return result;
    }

    std::pair<uint64_t, uint64_t> count(std::string &pattern) const
    {
        return r_index.count(pattern);
    }

    uint64_t input_size() const
    {
        return lzendsa_enc.input_size();
    }

    uint64_t num_phrases() const
    {
        return lzendsa_enc.num_phrases();
    }

    uint64_t size_in_bytes() const
    {
        return sizeof(this) + lzendsa_enc.size_in_bytes() + r_index.size_in_bytes();
    }

    void serialize(std::ostream &out) const
    {
        lzendsa_enc.serialize(out);
        r_index.serialize(out);
    }

    void load(std::istream &in)
    {
        lzendsa_enc.load(in);
        r_index.load(in);
    }
};