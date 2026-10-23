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

#ifndef MOVE_R_DATA_STRUCTURES_RANK_SELECT_SET_HPP
#define MOVE_R_DATA_STRUCTURES_RANK_SELECT_SET_HPP

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <data_structures/elias_fano_bit_vector.hpp>
#include <istream>
#include <limits>
#include <misc/utils.hpp>
#include <ostream>
#include <vector>

/**
 * @brief a set of integers in [0, universe) with rank and select. A dense set is stored as a bit vector in blocks of
 *        one cache line (the layout of cache_line_rank_bit_vector): the first word of a block is the number of
 *        elements before the block, the other words store the bits; rank and membership read one block. A sparse set
 *        is stored with Elias-Fano coding (see elias_fano_bit_vector); the smaller of the two is used.
 */
class rank_select_set {
protected:
    // ############################# TYPES AND PARAMETERS #############################

    // the number of bits of a word
    static constexpr uint64_t word_bits = 8 * sizeof(uint64_t);
    // the number of words of a block, which fills one cache line
    static constexpr uint64_t words_per_block = cache_aligned_no_init_allocator<uint64_t>::alignment / sizeof(uint64_t);
    // the number of bits of a block (its first word stores the number of elements before it)
    static constexpr uint64_t bits_per_block = (words_per_block - 1) * word_bits;
    // the number of elements between two select samples (as many as a block can hold)
    static constexpr uint64_t select_sample_rate = bits_per_block;

    // ############################# VARIABLES #############################

    bool dense = false; // true <=> the set is stored in words, else in sparse
    uint64_t universe = 0; // all elements are smaller than universe
    uint64_t num_elements = 0; // number of elements

    /* the blocks of a dense set: [0..words_per_block * (universe / bits_per_block + 1) - 1] block b is
       words[words_per_block * b..words_per_block * (b + 1) - 1]; words[words_per_block * b] is the number of elements
       smaller than b * bits_per_block, and element i is bit (i mod bits_per_block) mod word_bits of
       words[words_per_block * (i / bits_per_block) + 1 + (i mod bits_per_block) / word_bits] */
    std::vector<uint64_t, cache_aligned_no_init_allocator<uint64_t>> words;

    /* the select samples of a dense set: [0..(num_elements - 1) / select_sample_rate] select_samples[j] is the block
       containing the element with rank 1 + j * select_sample_rate */
    no_init_vector<uint32_t> select_samples;

    elias_fano_bit_vector sparse; // a sparse set

    // ############################# INTERNAL METHODS #############################

    // returns the r-th smallest element of a dense set
    uint64_t select_dense(uint64_t r) const
    {
        // the block of the element: the last block with fewer than r elements before it
        const uint64_t num_blocks = words.size() / words_per_block;
        uint64_t b = select_samples[(r - 1) / select_sample_rate];
        while (b + 1 < num_blocks && words[words_per_block * (b + 1)] < r) b++;
        uint64_t k = r - words[words_per_block * b]; // the rank of the element in the block (from 1)

        for (uint64_t w = 1; w < words_per_block; w++) {
            uint64_t word = words[words_per_block * b + w];
            const uint64_t count = std::popcount(word);

            if (k > count) {
                k -= count;
                continue;
            }

            // the k-th set bit of the word, by halving the searched range of bits
            uint64_t pos = 0;

            for (uint64_t width = word_bits / 2; width != 0; width /= 2) {
                const uint64_t low = word & ((uint64_t { 1 } << width) - 1);
                const uint64_t c = std::popcount(low);

                if (k > c) {
                    k -= c;
                    word >>= width;
                    pos += width;
                } else {
                    word = low;
                }
            }

            return b * bits_per_block + (w - 1) * word_bits + pos;
        }

        return universe;
    }

public:
    // ############################# CONSTRUCTORS #############################

    rank_select_set() = default;

    /**
     * @brief builds the set
     * @param values the elements in increasing order
     * @param universe all elements are smaller than universe
     * @param num_threads the number of threads
     */
    rank_select_set(const std::vector<uint64_t>& values, uint64_t universe, uint64_t num_threads = 1)
        : universe(universe), num_elements(values.size())
    {
        if (num_elements == 0) return;

        /* the bit vector is used if it is at most as large as the Elias-Fano coding and its block indices fit into
           the select samples */
        const double ef_bits = double(num_elements) * (2.0 + std::ceil(std::log2(std::max(1.0, double(universe) / double(num_elements)))));
        const uint64_t num_blocks = universe / bits_per_block + 1;
        dense = num_blocks - 1 <= std::numeric_limits<uint32_t>::max()
            && double(num_blocks * words_per_block * word_bits + (num_elements / select_sample_rate + 1) * 8 * sizeof(uint32_t))
                <= ef_bits;

        if (!dense) {
            sparse = elias_fano_bit_vector(values, universe, num_threads);
            return;
        }

        // set the bits, then the numbers of elements before the blocks
        words.assign(words_per_block * num_blocks, 0);

        for (uint64_t v : values) {
            const uint64_t offset = v % bits_per_block;
            words[words_per_block * (v / bits_per_block) + 1 + offset / word_bits] |= uint64_t { 1 } << (offset % word_bits);
        }

        uint64_t count = 0;

        for (uint64_t b = 0; b < num_blocks; b++) {
            words[words_per_block * b] = count;
            for (uint64_t w = 1; w < words_per_block; w++) count += std::popcount(words[words_per_block * b + w]);
        }

        select_samples.reserve(num_elements / select_sample_rate + 1);
        for (uint64_t r = 0; r < num_elements; r += select_sample_rate) select_samples.push_back(uint32_t(values[r] / bits_per_block));
    }

    // ############################# MISC PUBLIC METHODS #############################

    // returns whether the set is empty
    inline bool empty() const { return num_elements == 0; }

    // returns the number of elements
    inline uint64_t num_ones() const { return num_elements; }

    // returns the universe
    inline uint64_t size() const { return universe; }

    // returns whether the set is stored as a bit vector
    inline bool is_dense() const { return dense; }

    // returns the size in bytes
    uint64_t size_in_bytes() const
    {
        return words.size() * sizeof(uint64_t) + select_samples.size() * sizeof(uint32_t)
            + (dense ? 0 : sparse.size_in_bytes());
    }

    // fetches the data that rank_1(i) and contains(i) read into the cache
    inline void prefetch(uint64_t i) const
    {
        if (dense) __builtin_prefetch(words.data() + words_per_block * (i / bits_per_block));
    }

    // ############################# QUERY METHODS #############################

    /**
     * @brief returns the number of elements smaller than i
     * @param i [0..universe]
     */
    inline uint64_t rank_1(uint64_t i) const
    {
        if (!dense) return sparse.rank_1(i);
        if (i >= universe) [[unlikely]] return num_elements;
        const uint64_t* B = words.data() + words_per_block * (i / bits_per_block);
        const uint64_t offset = i % bits_per_block;
        const uint64_t full = offset / word_bits;
        uint64_t rank = B[0];
        for (uint64_t w = 1; w <= full; w++) rank += std::popcount(B[w]);
        return rank + std::popcount(B[1 + full] & ((uint64_t { 1 } << (offset % word_bits)) - 1));
    }

    /**
     * @brief returns whether i is an element
     * @param i [0..universe - 1]
     */
    inline bool contains(uint64_t i) const
    {
        if (!dense) return sparse.rank_1(i + 1) != sparse.rank_1(i);
        const uint64_t offset = i % bits_per_block;
        return (words[words_per_block * (i / bits_per_block) + 1 + offset / word_bits] >> (offset % word_bits)) & 1;
    }

    /**
     * @brief returns the number of elements smaller than i, if i is an element, else num_ones()
     * @param i [0..universe - 1]
     */
    inline uint64_t index_of(uint64_t i) const
    {
        if (!dense) {
            const uint64_t rank = sparse.rank_1(i);
            return rank < num_elements && sparse.select_1(rank + 1) == i ? rank : num_elements;
        }

        const uint64_t offset = i % bits_per_block;
        const uint64_t* B = words.data() + words_per_block * (i / bits_per_block);
        if (!((B[1 + offset / word_bits] >> (offset % word_bits)) & 1)) return num_elements;
        return rank_1(i);
    }

    /**
     * @brief returns the r-th smallest element
     * @param r [1..num_ones()]
     */
    inline uint64_t select_1(uint64_t r) const
    {
        if (!dense) return sparse.select_1(r);
        return select_dense(r);
    }

    // ############################# SERIALIZATION METHODS #############################

    /**
     * @brief writes the set to an output stream
     * @param out the output stream
     */
    void serialize(std::ostream& out) const
    {
        out.write((const char*) &dense, sizeof(bool));
        out.write((const char*) &universe, sizeof(uint64_t));
        out.write((const char*) &num_elements, sizeof(uint64_t));

        if (!dense) {
            if (num_elements != 0) sparse.serialize(out);
            return;
        }

        const uint64_t num_words = words.size();
        const uint64_t num_samples = select_samples.size();
        out.write((const char*) &num_words, sizeof(uint64_t));
        out.write((const char*) words.data(), num_words * sizeof(uint64_t));
        out.write((const char*) &num_samples, sizeof(uint64_t));
        out.write((const char*) select_samples.data(), num_samples * sizeof(uint32_t));
    }

    /**
     * @brief reads a set written by serialize
     * @param in the input stream
     */
    void load(std::istream& in)
    {
        in.read((char*) &dense, sizeof(bool));
        in.read((char*) &universe, sizeof(uint64_t));
        in.read((char*) &num_elements, sizeof(uint64_t));
        words.clear();
        select_samples.clear();
        sparse = elias_fano_bit_vector();

        if (!dense) {
            if (num_elements != 0) sparse.load(in);
            return;
        }

        uint64_t num_words;
        uint64_t num_samples;
        in.read((char*) &num_words, sizeof(uint64_t));
        words.resize(num_words);
        in.read((char*) words.data(), num_words * sizeof(uint64_t));
        in.read((char*) &num_samples, sizeof(uint64_t));
        resize_no_init(select_samples, num_samples);
        in.read((char*) select_samples.data(), num_samples * sizeof(uint32_t));
    }
};

#endif // MOVE_R_DATA_STRUCTURES_RANK_SELECT_SET_HPP
