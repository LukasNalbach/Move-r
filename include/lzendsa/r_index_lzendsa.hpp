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

#include "lzendsa_encoding.hpp"

#include <algorithms/build_sa_and_bwt.hpp>
#include <custom_r_index/custom_r_index.hpp>
#include <algorithms/sparse_sa_bin_search.tpp>

/**
 * @brief an r-index combined with an LZ-End encoding of the differential suffix array for locating
 * @tparam int_t signed integer type of the suffix array entries
 */
template<typename int_t = int32_t>
class r_index_lzendsa {

protected:
    lzendsa_encoding lzendsa_enc;
    custom_r_index::index<_run_ends> r_idx;

public:
    r_index_lzendsa() = default;

    /**
     * @brief builds the index from the input text
     * @param input the input text
     * @param h the maximum phrase length (default: 8192)
     * @param use_bigbwt whether to use Big-BWT to build the suffix array and BWT
     * @param log whether to print log messages
     */
    r_index_lzendsa(std::string& input, int_t h = 8192, bool use_bigbwt = false, bool log = false)
    {
        auto time_start = now();
        auto time = time_start;

        // compute SA and BWT
        auto [sa, bwt] = build_sa_and_bwt<int_t>(input, true, use_bigbwt, log);
        uint64_t n = sa.size();

        // build r-index
        log_phase_start(log, time, "building r-index");
        r_idx = custom_r_index::index<_run_ends>(bwt, sa, true);
        bwt.clear();
        bwt.shrink_to_fit();
        log_phase_end(log, time);

        // construct differential suffix array (DSA)
        log_phase_start(log, time, "building Differential Suffix Array (DSA)");
        auto& dsa = sa;

        for (uint64_t i = n - 1; i > 0; i--) {
            dsa[i] = sa[i] - sa[i - 1];
        }

        log_phase_end(log, time);
        
        // compute lzend parsing
        log_message(log, "building LZ-End parsing:\n");
        std::vector<lzend_phr_t<int_t>> lzend_phrases = construct_lzend_of_reverse<int_t>(dsa, h, log);

        log_phase_start(log, time, "transforming DSA to SA");

        // make SA out of DSA
        for (uint64_t i = 1; i < n; i++) {
            sa[i] = dsa[i] + dsa[i - 1];
        }

        log_phase_end(log, time);

        // construct lzend encoding
        log_phase_start(log, time, "Encoding LZ-End parsing");
        lzendsa_enc = lzendsa_encoding(lzend_phrases, n, h);
        log_phase_end(log, time);
        
        log_message(log, "r-index-lzendsa built in " + format_time(time_diff_ns(time_start, now())) + "\n");
    }

    /**
     * @brief locates the occurrences of pattern in the input
     * @param pattern the pattern to locate
     * @return the starting positions of the occurrences of pattern in the input
     */
    std::vector<int_t> locate(const std::string &pattern) const
    {
        auto [beg, end, last_value] = r_idx.count_and_get_occ(pattern);
        if (end < beg) return {};
        std::vector<int_t> result = lzendsa_enc.template extract_deltas<int_t>(beg + 1, end);
        result.push_back(last_value);

        for (int_t i = result.size() - 2; i >= 0; i--) {
            result[i] = result[i + 1] - result[i];
        }
        
        return result;
    }

    // returns a reference to the r-index
    const custom_r_index::index<_run_ends>& r_index() const { return r_idx; }

    // returns a reference to the LZ-End encoding
    const lzendsa_encoding& sa_encoding() const { return lzendsa_enc; }

    /**
     * @brief counts the occurrences of pattern in the input
     * @param pattern the pattern to count
     * @return the suffix array interval [beg, end] of the pattern
     */
    std::pair<uint64_t, uint64_t> count(std::string &pattern) const { return r_idx.count(pattern); }

    // returns the size of the input text
    uint64_t input_size() const { return lzendsa_enc.input_size(); }

    // returns the number of phrases in the LZ-End encoding
    uint64_t num_phrases() const { return lzendsa_enc.num_phrases(); }

    // returns the number of suffix array samples (one per BWT run)
    uint64_t num_samples() const { return r_idx.num_bwt_runs(); }

    /**
     * @brief returns the i-th suffix array sample
     * @param i a sample index
     * @return the i-th suffix array sample
     */
    uint64_t sample(uint64_t i) const { return r_idx.sample(i); }

    /**
     * @brief returns the position in the suffix array of the i-th sample
     * @param i a sample index
     * @return the position in SA of the i-th sample
     */
    uint64_t sample_pos(uint64_t i) const { return r_idx.sample_pos(i); }

    // returns the size of the data structure in bytes
    uint64_t size_in_bytes() const { return sizeof(this) + lzendsa_enc.size_in_bytes() + r_idx.size_in_bytes(); }

    /**
     * @brief serializes the index to an output stream
     * @param out output stream
     */
    void serialize(std::ostream &out) const
    {
        lzendsa_enc.serialize(out);
        r_idx.serialize(out);
    }

    /**
     * @brief loads the index from an input stream
     * @param in input stream
     */
    void load(std::istream &in)
    {
        lzendsa_enc.load(in);
        r_idx.load(in);
    }
};