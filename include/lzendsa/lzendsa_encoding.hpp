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
    int64_t h; // maximum phrase length
    int64_t min_ext; // the minimum over all extensions in the lzend parsing

    /* Stores two bit-packed vectors sources and extensions interleaved with each other.
     * The vector sources stores at position i the index of the phrase the ith phrase extends.
     * The vector extensions stores at position i the extension of the ith phrase - min_ext. */
    interleaved_bit_aligned_vectors<uint64_t, 2> sources_extensions;

    /* Stores at position i the end position of the ith phrase.
     * It essentially serves as a constant-time select operation (index of the j-th 1) and a
     * log(z)-time rank operation over the end positions (number of 1s in a given interval).
     * The first phrase ends at position 0, so end_positions[0] = 0 holds. */
    interleaved_byte_aligned_vectors<uint64_t, uint64_t, 1> end_positions;

public:
    lzendsa_encoding() = default;

    template <typename int_t>
    lzendsa_encoding(const std::vector<lzend_phr_t<int_t>>& lzend_phrases,int64_t n, int64_t h = 8192)
        : n(n), h(h), z(lzend_phrases.size())
    {
        min_ext = lzend_phrases[0].ext;
        int64_t max_ext = lzend_phrases[0].ext;

        for (uint64_t i = 1; i < z; i++) {
            if (min_ext > lzend_phrases[i].ext) min_ext = lzend_phrases[i].ext;
            if (max_ext < lzend_phrases[i].ext) max_ext = lzend_phrases[i].ext;
        }

        sources_extensions = interleaved_bit_aligned_vectors<uint64_t, 2>({
            std::bit_width(uint64_t{z}), // bit-width of the sources vector
            std::bit_width(uint64_t{max_ext - min_ext}) // bit-width of the extensions vector
        });

        end_positions = interleaved_byte_aligned_vectors<uint64_t, uint64_t, 1>({
            byte_width(uint64_t{n})
        });

        sources_extensions.resize_no_init(z);
        end_positions.resize_no_init(z);

        // current phrase end position
        int64_t cur_end_pos = -1;

        for (uint64_t i = 0; i < z; i++) {
            cur_end_pos += lzend_phrases[i].len; // update current end position
            sources_extensions.set<0>(i, lzend_phrases[i].lnk); // store the sources of phrases
            sources_extensions.set<1>(i, lzend_phrases[i].ext - min_ext); // store the extensions of phrases
            end_positions.set<0>(i, cur_end_pos); // store the end positions of phrases
        }
    }

    inline int64_t maximum_phrase_length() const
    {
        return h;
    }

    inline int64_t source(int64_t i) const
    {
        return sources_extensions.get<0>(i);
    }

    inline int64_t extension(int64_t i) const
    {
        return min_ext + int64_t{sources_extensions.get<1>(i)};
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
        return sizeof(this) + sources_extensions.size_in_bytes() + end_positions.size_in_bytes();
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
    void extract_deltas(int64_t beg, int64_t end, const std::function<void(int64_t, int64_t)>& report) const
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
                int64_t end_pos = end_position(phrase_id);

                if (end_pos == end) {
                    report(pos_in_sa, extension(phrase_id));
                    --pos_in_sa;

                    if (phrase_id != 0 && end_position(phrase_id - 1) == end - 1) {
                        --phrase_id;
                    }

                    --end;
                } else {
                    src = source(phrase_id);
                    phrase_pos_shift = - end_pos + end_position(src) + 1;
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
    std::vector<out_t> extract_deltas(int64_t beg, int64_t end) const
    {
        std::vector<out_t> result;
        no_init_resize(result, end - beg + 1);
        extract_deltas(beg, end, [&](int64_t pos, int64_t val){result[pos - beg] = val;});
        return result;
    }

    void extract(int64_t beg, int64_t end, int64_t sa_end, const std::function<void(int64_t, int64_t)>& report) const
    {
        report(end, sa_end);
        int64_t cur_val = sa_end;

        extract_deltas(beg + 1, end, [&](int64_t pos, int64_t val){
            cur_val -= val;
            report(pos - 1, cur_val);
        });
    }

    // extraction method
    template <typename out_t>
    std::vector<out_t> extract(int64_t beg, int64_t end, int64_t sa_end) const
    {
        std::vector<out_t> result;
        no_init_resize(result, end - beg + 1);
        extract(beg, end, [&](int64_t pos, int64_t val){result[pos - beg] = val;});
        return result;
    }

    void load(std::istream& in)
    {
        in.read((char*) &n, sizeof(n));
        in.read((char*) &z, sizeof(z));
        in.read((char*) &h, sizeof(h));
        in.read((char*) &min_ext, sizeof(min_ext));
        sources_extensions.load(in);
        end_positions.load(in);
    }

    void serialize(std::ostream& out) const
    {
        out.write((char*) &n, sizeof(n));
        out.write((char*) &z, sizeof(z));
        out.write((char*) &h, sizeof(h));
        out.write((char*) &min_ext, sizeof(min_ext));
        sources_extensions.serialize(out);
        end_positions.serialize(out);
    }

    void log_data_structure_sizes() const {
        std::cout << "sources: " << format_size((sources_extensions.width<0>() * z) / 8) << std::endl;
        std::cout << "extensions: " << format_size((sources_extensions.width<1>() * z) / 8) << std::endl;
        std::cout << "end_positions: " << format_size(end_positions.size_in_bytes()) << std::endl;
    }
};