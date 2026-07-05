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

#include <algorithms/build_sa_and_bwt.hpp>
#include <custom_r_index/custom_r_index.hpp>
#include "rlzsa_encoding.hpp"

/**
 * @brief an r-index combined with a standalone rlzsa encoding of the differential suffix array for locating
 * @tparam int_t signed integer type of the suffix array entries
 */
template<typename int_t = int32_t>
class r_index_rlzsa {

protected:
    rlzsa_encoding<int_t> rlzsa_enc;
    custom_r_index::index<_run_heads> r_idx;

public:
    r_index_rlzsa() = default;

    /**
     * @brief builds the index from the input text
     * @param input the input text
     * @param use_bigbwt whether to use Big-BWT to build the suffix array and BWT
     * @param use_r_index_samples whether to use the r-index's suffix array samples for locating
     * @param log whether to print log messages
     */
    r_index_rlzsa(std::string& input, bool use_bigbwt = false, bool use_r_index_samples = false, bool log = false)
    {
        auto time_start = now();
        auto time = time_start;

        // compute SA and BWT
        auto [sa, bwt] = build_sa_and_bwt<int_t>(input, true, use_bigbwt, log);
        uint64_t n = sa.size();

        // build r-index
        log_phase_start(log, time, "building r-index");
        r_idx = custom_r_index::index<_run_heads>(bwt, sa, use_r_index_samples);
        bwt.clear();
        bwt.shrink_to_fit();
        log_phase_end(log, time);

        // construct differential suffix array (DSA)
        log_phase_start(log, time, "building Differential Suffix Array (DSA)");
        std::vector<int_t> dsa;
        no_init_resize(dsa, n);
        dsa[0] = sa[0];

        for (uint64_t i = 1; i < n; i++) {
            dsa[i] = sa[i] - sa[i - 1];
        }

        log_phase_end(log, time);
        
        // compute rlzsa
        log_message(log, "building rlzsa:\n");
        uint64_t r = r_idx.number_of_runs();
        uint64_t reference_size = std::min<uint64_t>(n / 3, 5.2 * r);
        rlzsa_enc = rlzsa_encoding<int_t>(sa, std::move(dsa), reference_size, !use_r_index_samples, log);

        if (log) {
            std::cout << "r-index-rlzsa built in " << format_time(time_diff_ns(time_start, now())) << std::endl;
            std::cout << "relative reference size: " << reference_size / (double) n << std::endl;
        }
    }

    /**
     * @brief returns a reference to the r-index
     * @return the r-index
     */
    const custom_r_index::index<_run_heads>& r_index() const
    {
        return r_idx;
    }

    /**
     * @brief returns a reference to the rlzsa encoding
     * @return the rlzsa encoding
     */
    const rlzsa_encoding<int_t>& sa_encoding() const
    {
        return rlzsa_enc;
    }

    /**
     * @brief locates the occurrences of pattern in the input
     * @tparam out_t output value type
     * @param pattern the pattern to locate
     * @return the starting positions of the occurrences of pattern in the input
     */
    template <typename out_t>
    std::vector<out_t> locate(const std::string &pattern) const
    {
        std::vector<out_t> Occ;

        if (r_idx.has_sa_samples()) {
            auto [beg, end, sa_beg] = r_idx.count_and_get_occ(pattern);
            if (end < beg) return {};
            Occ = rlzsa_enc.template extract<out_t>(beg, end, (int64_t)sa_beg);
        } else {
            auto [beg, end] = r_idx.count(pattern);
            if (end < beg) return {};
            Occ = rlzsa_enc.template extract<out_t>(beg, end);
        }

        return Occ;
    }

    /**
     * @brief counts the occurrences of pattern in the input
     * @param pattern the pattern to count
     * @return the suffix array interval [beg, end] of the pattern
     */
    std::pair<uint64_t, uint64_t> count(std::string &pattern) const
    {
        return r_idx.count(pattern);
    }

    /**
     * @brief returns the size of the input text
     * @return the size of the input text
     */
    uint64_t input_size() const
    {
        return rlzsa_enc.input_size();
    }

    /**
     * @brief returns the number of phrases in the rlzsa encoding
     * @return the number of phrases
     */
    uint64_t num_phrases() const
    {
        return rlzsa_enc.num_phrases();
    }

    /**
     * @brief returns the number of suffix array samples (one per BWT run)
     * @return the number of suffix array samples
     */
    uint64_t num_samples() const
    {
        return r_idx.num_bwt_runs();
    }

    /**
     * @brief returns the i-th suffix array sample
     * @param i a sample index
     * @return the i-th suffix array sample
     */
    uint64_t sample(uint64_t i) const
    {
        return r_idx.sample(i);
    }

    /**
     * @brief returns the size of the data structure in bytes
     * @return the size of the data structure in bytes
     */
    uint64_t size_in_bytes() const
    {
        return sizeof(this) + rlzsa_enc.size_in_bytes() + r_idx.size_in_bytes();
    }

    /**
     * @brief serializes the index to an output stream
     * @param out output stream
     */
    void serialize(std::ostream &out) const
    {
        rlzsa_enc.serialize(out);
        r_idx.serialize(out);
    }

    /**
     * @brief loads the index from an input stream
     * @param in input stream
     */
    void load(std::istream &in)
    {
        rlzsa_enc.load(in);
        r_idx.load(in);
    }
};