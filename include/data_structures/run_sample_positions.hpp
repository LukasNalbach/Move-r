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

#ifndef MOVE_R_DATA_STRUCTURES_RUN_SAMPLE_POSITIONS_HPP
#define MOVE_R_DATA_STRUCTURES_RUN_SAMPLE_POSITIONS_HPP

#include <bit>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <new>
#include <utility>
#include <vector>

#include <omp.h>

#include <misc/files.hpp>
#include <misc/utils.hpp>

/**
 * @brief maps the input intervals of M_LF to the positions of their suffix array samples in an array that
 *        stores the samples at the BWT run boundaries in run order
 *
 * For each BWT run, the array stores the sample at its start and, if the run is longer than 1, the sample at
 * its end. Two bits per input interval x mark its samples: b[2x] = 1 <=> x starts a run, and b[2x + 1] = 1
 * <=> x ends a run of length at least 2. With rank(i) the number of ones before position i, the sample at the
 * start of a run starting at x is the rank(2x)-th, and the sample at the end of a run ending at x is the
 * (rank(2x + 2) - 1)-th (for a run of length 1, this is the sample at its start).
 *
 * The bits are stored in blocks of one cache line (the layout of cache_line_rank_bit_vector): the first word of a
 * block stores the number of ones before the block, the other seven words store 448 bits.
 */
class run_sample_positions {
protected:
    // ############################# TYPES AND PARAMETERS #############################

    static constexpr uint64_t word_bits = 8 * sizeof(uint64_t); // number of bits of a word
    static constexpr uint64_t words_per_block = cache_aligned_no_init_allocator<uint64_t>::alignment / sizeof(uint64_t); // one cache line
    static constexpr uint64_t bits_per_block = word_bits * (words_per_block - 1); // number of bits of b per block
    static constexpr uint64_t intervals_per_block = bits_per_block / 2; // number of input intervals per block

    // ############################# VARIABLES #############################

    /* the blocks: [0..words_per_block * ((2 r') / bits_per_block + 1) - 1] block k is blocks[words_per_block * k..
       words_per_block * (k + 1) - 1]; blocks[words_per_block * k] is the number of ones in b before the block, and
       b[i] is bit (i mod bits_per_block) mod word_bits of blocks[words_per_block * (i / bits_per_block) + 1 +
       (i mod bits_per_block) / word_bits] */
    std::vector<uint64_t, cache_aligned_no_init_allocator<uint64_t>> blocks;

    uint64_t num_intervals = 0; // r', the number of input intervals of M_LF
    uint64_t _num_samples = 0; // number of samples in the array

    // ############################# INTERNAL METHODS #############################

    // returns the number of ones before position i
    inline uint64_t rank(uint64_t i) const
    {
        const uint64_t b = i / bits_per_block;
        const uint64_t offset = i - b * bits_per_block;
        const uint64_t* B = blocks.data() + words_per_block * b;
        uint64_t ones = B[0];
        const uint64_t full = offset / word_bits;
        for (uint64_t w = 1; w <= full; w++) ones += std::popcount(B[w]);
        if (offset % word_bits != 0) ones += std::popcount(B[1 + full] & ((uint64_t { 1 } << (offset % word_bits)) - 1));
        return ones;
    }

public:
    // ############################# CONSTRUCTORS #############################

    run_sample_positions() = default;

    /**
     * @brief builds the data structure
     * @param num_intervals r', the number of input intervals of M_LF
     * @param starts_run starts_run(x) = true <=> the x-th input interval starts a BWT run (x = 0 does)
     * @param len len(x) = the length of the x-th input interval
     * @param num_threads number of threads to use
     * @param with_ends false <=> the array stores only the samples at the starts of the runs
     *        (see move_r::has_run_end_samples); end_sample_index must not be called then
     */
    template <typename starts_run_fnc_t, typename len_fnc_t>
    run_sample_positions(uint64_t num_intervals, starts_run_fnc_t starts_run, len_fnc_t len, uint64_t num_threads,
        bool with_ends = true)
        : num_intervals(num_intervals)
    {
        // one additional block, s.t. rank(2 * num_intervals) is in bounds
        const uint64_t num_blocks = (2 * num_intervals) / bits_per_block + 1;
        blocks.resize(num_blocks * words_per_block); // not zeroed here, but below by all threads

        // set the bits of each block and count its ones
        #pragma omp parallel for num_threads(num_threads) schedule(static)
        for (uint64_t b = 0; b < num_blocks; b++) {
            uint64_t* B = blocks.data() + words_per_block * b;
            for (uint64_t w = 0; w < words_per_block; w++) B[w] = 0;
            const uint64_t x_end = std::min(num_intervals, (b + 1) * intervals_per_block);

            for (uint64_t x = b * intervals_per_block; x < x_end; x++) {
                const bool starts = starts_run(x);
                // x ends a run of length >= 2 <=> x + 1 starts a run and x is not a run of length 1
                const bool ends_long = with_ends && (x + 1 == num_intervals || starts_run(x + 1)) &&
                    !(starts && len(x) == 1);
                const uint64_t bit = 2 * (x - b * intervals_per_block);
                B[1 + bit / word_bits] |= (uint64_t(starts) << (bit % word_bits)) | (uint64_t(ends_long) << (bit % word_bits + 1));
            }

            for (uint64_t w = 1; w < words_per_block; w++) B[0] += std::popcount(B[w]);
        }

        // replace the number of ones in each block by the number of ones before it
        for (uint64_t b = 0; b < num_blocks; b++) {
            uint64_t* B = blocks.data() + words_per_block * b;
            const uint64_t ones = B[0];
            B[0] = _num_samples;
            _num_samples += ones;
        }
    }

    // ############################# MISC PUBLIC METHODS #############################

    // returns the number of samples in the array (the number of BWT runs plus the number of BWT runs longer than 1)
    inline uint64_t num_samples() const { return _num_samples; }

    // returns whether the data structure is empty (i.e., has not been built)
    inline bool empty() const { return blocks.empty(); }

    // returns the size in bytes
    uint64_t size_in_bytes() const { return blocks.size() * sizeof(uint64_t); }

    // ############################# QUERY METHODS #############################

    /**
     * @brief returns the position in the array of the sample at the start of the x-th input interval
     * @param x [0..r'-1] an input interval that starts a BWT run
     */
    inline uint64_t start_sample_index(uint64_t x) const { return rank(2 * x); }

    /**
     * @brief returns the position in the array of the sample at the end of the x-th input interval
     * @param x [0..r'-1] an input interval that ends a BWT run
     */
    inline uint64_t end_sample_index(uint64_t x) const { return rank(2 * x + 2) - 1; }

    // ############################# SERIALIZATION METHODS #############################

    // writes the data structure to an output stream
    void serialize(std::ostream& out) const
    {
        out.write((char*) &num_intervals, sizeof(uint64_t));
        out.write((char*) &_num_samples, sizeof(uint64_t));
        uint64_t num_words = blocks.size();
        out.write((char*) &num_words, sizeof(uint64_t));
        write_in_chunks(out, (const char*) blocks.data(), num_words * sizeof(uint64_t));
    }

    // reads the data structure written by serialize
    void load(std::istream& in)
    {
        in.read((char*) &num_intervals, sizeof(uint64_t));
        in.read((char*) &_num_samples, sizeof(uint64_t));
        uint64_t num_words;
        in.read((char*) &num_words, sizeof(uint64_t));
        blocks.resize(num_words);
        read_in_chunks(in, (char*) blocks.data(), num_words * sizeof(uint64_t));
    }
};

#endif // MOVE_R_DATA_STRUCTURES_RUN_SAMPLE_POSITIONS_HPP
