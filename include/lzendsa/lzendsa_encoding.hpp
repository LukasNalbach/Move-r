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

#include <algorithm>
#include <bit>
#include <omp.h>
#include <vector>

#include "lzendsa_construction.hpp"
#include <data_structures/interleaved_bit_aligned_vectors.hpp>
#include <data_structures/interleaved_byte_aligned_vectors.hpp>

/**
 * This lzendsa_encoding is designed to be rather small (3 words per phrase) while being simple and enabling efficient extraction.
 */
class lzendsa_encoding {

protected:
    uint64_t n; // the number of symbols in the uncompressed text
    uint64_t z; // the number of phrases of the lzend factorization
    int64_t min_ext; // the minimum over all extensions in the lzend parsing

    /* Stores three bit-packed vectors sources and extensions interleaved with each other.
     * The vector sources stores at position i the index of the phrase the ith phrase extends.
     * The vector extensions stores at position i the extension of the ith phrase - min_ext. */
    interleaved_bit_aligned_vectors<uint64_t, 3> sources_extensions_samples;

    /* Stores at position i the end position of the ith phrase.
     * It essentially serves as a constant-time select operation (index of the j-th 1) and a
     * log(z)-time rank operation over the end positions (number of 1s in a given interval).
     * The first phrase ends at position 0, so end_positions[0] = 0 holds. */
    interleaved_byte_aligned_vectors<uint64_t, uint64_t, 1> end_positions;

public:
    lzendsa_encoding() = default;

    template <typename int_t>
    lzendsa_encoding(
        const std::vector<lzend_phr_t<int_t>>& lzend_phrases,
        const std::vector<int_t>& sa,
        int64_t n, bool sample_phrase_ends = false
    ) : n(n), z(lzend_phrases.size())
    {
        min_ext = lzend_phrases[0].ext;
        int64_t max_ext = lzend_phrases[0].ext;

        for (uint64_t i = 1; i < z; i++) {
            if (min_ext > lzend_phrases[i].ext) min_ext = lzend_phrases[i].ext;
            if (max_ext < lzend_phrases[i].ext) max_ext = lzend_phrases[i].ext;
        }

        sources_extensions_samples = interleaved_bit_aligned_vectors<uint64_t, 3>({
            std::bit_width(uint64_t{z}), // bit-width of the sources vector
            std::bit_width(uint64_t{max_ext - min_ext}), // bit-width of the extensions vector
            sample_phrase_ends ? std::bit_width(uint64_t{n}) : 0
        });

        end_positions = interleaved_byte_aligned_vectors<uint64_t, uint64_t, 1>({
            byte_width(uint64_t{n}) // bit-width of the end_positions vector
        });

        sources_extensions_samples.resize_no_init(z);
        end_positions.resize_no_init(z);

        // current phrase end position
        int64_t cur_end_pos = -1;

        for (uint64_t i = 0; i < z; i++) {
            cur_end_pos += lzend_phrases[i].len; // update current end position

            sources_extensions_samples.set<0>(i, lzend_phrases[i].lnk); // store the sources of phrases
            sources_extensions_samples.set<1>(i, lzend_phrases[i].ext - min_ext); // store the extensions of phrases

            if (sample_phrase_ends) {
                sources_extensions_samples.set<2>(i, sa[cur_end_pos]); // store samples at phrase ends
            }

            end_positions.set<0>(i, cur_end_pos); // store the end positions of phrases
        }
    }

    inline bool has_sa_samples() const
    {
        return sources_extensions_samples.width<2>() != 0;
    }

    inline int64_t source(int64_t i) const
    {
        return sources_extensions_samples.get<0>(i);
    }

    inline int64_t extension(int64_t i) const
    {
        return min_ext + int64_t{sources_extensions_samples.get<1>(i)};
    }

    inline int64_t sample(int64_t i) const
    {
        return sources_extensions_samples.get<2>(i);
    }

    inline int64_t end_position(int64_t i) const
    {
        return end_positions.get<0>(i);
    }
    
    inline uint64_t num_phrases() const
    {
        return z;
    }

    inline uint64_t input_size() const
    {
        return n;
    }

    uint64_t size_in_bytes() const
    {
        return sizeof(this) + sources_extensions_samples.size_in_bytes() + end_positions.size_in_bytes();
    }

    // returns the phrase that contains index
    int64_t phrase_containing(int64_t index) const
    {
        return bin_search_min_geq<uint64_t>(index, 0, z - 1, [&](uint64_t x) { return end_position(x); });
    }

    // returns the phrase that ends at or precedes index
    int64_t phrase_preceding_or_ending_at(int64_t index) const
    {
        return bin_search_max_leq<uint64_t>(index, 0, z - 1, [&](uint64_t x) { return end_position(x); });
    }

    // extraction method
    void extract(int64_t beg, int64_t end, const std::function<void(int64_t, int64_t)>& report) const
    {
        if (end < beg) return;
        int64_t phrase_id = phrase_containing(end);
        int64_t pos_in_sa = end;

        struct SearchInterval {
            int64_t beg;
            int64_t end;
            int64_t phrase_id;
        };

        std::vector<SearchInterval> intervals;
        intervals.reserve(32);
        intervals.push_back({ beg, end, phrase_id });

        // to avoid recusion, we simply store the starting- and endpositions and the corresponding phrase id in a stack
        int64_t pos, src, phrase_pos_shift;

        while (!intervals.empty()) {
            SearchInterval current_interval = intervals.back();
            intervals.pop_back();

            beg = current_interval.beg;
            end = current_interval.end;
            phrase_id = current_interval.phrase_id;

            if (phrase_id > 0) {
                pos = end_position(phrase_id - 1) + 1; // first position of the current phrase
            } else {
                pos = 0;
            }

            if (beg < pos) {
                intervals.push_back({ beg, pos - 1, phrase_id - 1 });
                beg = pos;
            }

            while (end >= beg) {
                if (end_position(phrase_id) == end) {
                    report(pos_in_sa, extension(phrase_id));
                    --pos_in_sa;

                    if (phrase_id != 0 && end_position(phrase_id - 1) == end - 1) {
                        --phrase_id;
                    }

                    --end;
                } else {
                    src = source(phrase_id);
                    phrase_pos_shift = - end_position(phrase_id) + end_position(src) + 1;
                    beg += phrase_pos_shift;
                    end += phrase_pos_shift;

                    while (src != 0 && end_position(src - 1) >= end) {
                        src--;
                    }

                    phrase_id = src;
                    int64_t prev_end_pos = end_position(phrase_id - 1);

                    // make pos
                    if (phrase_id > 0 && beg <= prev_end_pos) {
                        pos = prev_end_pos + 1; // first position of the current phrase
                        intervals.push_back({ beg, pos - 1, phrase_id - 1 });
                        beg = pos;
                    }
                }
            }
        }
    }

    // extraction method
    template <typename out_t>
    std::vector<out_t> extract(int64_t beg, int64_t end) const
    {
        std::vector<out_t> result;
        no_init_resize(result, end - beg + 1);
        extract(beg, end, [&](int64_t pos, int64_t val){result[pos - beg] = val;});
        return result;
    }

    int64_t operator[](int64_t i) const
    {
        int64_t val;
        extract(i, i, [&](int64_t, int64_t v){val = v;});
        return val;
    }

    void load(std::istream& in)
    {
        in.read((char*) &n, sizeof(n));
        in.read((char*) &z, sizeof(z));
        in.read((char*) &min_ext, sizeof(min_ext));
        sources_extensions_samples.load(in);
        end_positions.load(in);
    }

    void serialize(std::ostream& out) const
    {
        out.write((char*) &n, sizeof(n));
        out.write((char*) &z, sizeof(z));
        out.write((char*) &min_ext, sizeof(min_ext));
        sources_extensions_samples.serialize(out);
        end_positions.serialize(out);
    }

    void log_data_structure_sizes() const {
        std::cout << "sources: " << format_size((sources_extensions_samples.width<0>() * z) / 8) << std::endl;
        std::cout << "extensions: " << format_size((sources_extensions_samples.width<1>() * z) / 8) << std::endl;
        std::cout << "end_positions: " << format_size(end_positions.size_in_bytes()) << std::endl;

        if (sources_extensions_samples.width<2>() != 0) {
            std::cout << "samples: " << format_size((sources_extensions_samples.width<2>() * z) / 8) << std::endl;
        }
        
    }
};