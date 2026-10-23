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

#ifndef MOVE_R_DATA_STRUCTURES_BLOCKED_POSITIONS_HPP
#define MOVE_R_DATA_STRUCTURES_BLOCKED_POSITIONS_HPP

#include <algorithm>
#include <cstdint>
#include <istream>
#include <limits>
#include <misc/utils.hpp>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>

/**
 * @brief an increasing sequence of positions smaller than 2^63, stored in blocks of 32 bytes: a block stores its first
 *        position and the offsets of the following positions from it, either 12 offsets of 16 bits (13 positions per
 *        block) or 6 offsets of 32 bits (7 positions per block), whichever needs less space. The offsets of a block
 *        that do not fit are stored in a separate array of 64-bit offsets. A position, and the one after it, is read
 *        from one block (the blocks are aligned to 32 bytes).
 */
class blocked_positions {
protected:
    // ############################# TYPES AND PARAMETERS #############################

    // number of bits of block_t::first_and_flag that store the first position of the block
    static constexpr uint64_t first_pos_bits = 8 * sizeof(uint64_t) - 1;
    static constexpr uint64_t mask_first_pos = (uint64_t { 1 } << first_pos_bits) - 1;

    /* a block of positions_per_block consecutive positions; the low first_pos_bits bits of first_and_flag store the
       first position of the block, its highest bit is 1 <=> the offsets of the block are stored in wide_offsets; then
       offsets[0] is the index of the block among the blocks with offsets in wide_offsets, else offsets stores the
       positions_per_block - 1 offsets of the other positions from the first one, as 16- or 32-bit integers */
    struct alignas(4 * sizeof(uint64_t)) block_t {
        uint64_t first_and_flag;
        uint64_t offsets[3];
    };

    static_assert(sizeof(block_t) == alignof(block_t));

    // number of positions of a block with 16-bit and with 32-bit offsets
    static constexpr uint64_t positions_per_block_16 = 1 + (sizeof(block_t::offsets)) / sizeof(uint16_t);
    static constexpr uint64_t positions_per_block_32 = 1 + (sizeof(block_t::offsets)) / sizeof(uint32_t);

    // ############################# VARIABLES #############################

    uint64_t num_positions = 0; // number of positions
    uint64_t positions_per_block = positions_per_block_16; // positions_per_block_16 or positions_per_block_32

    /* the blocks: [0..ceil(num_positions / positions_per_block) - 1] blocks[b] stores the positions with indices
       b * positions_per_block .. (b + 1) * positions_per_block - 1 */
    no_init_vector<block_t> blocks;

    /* the offsets of the blocks whose largest offset does not fit into the offsets of a block: [0..(number of such
       blocks) * (positions_per_block - 1) - 1] the offsets of the j-th such block (in block order) are at
       [(positions_per_block - 1) j..(positions_per_block - 1) (j + 1) - 1]; the unused offsets of the last block are 0 */
    no_init_vector<uint64_t> wide_offsets;

    // ############################# INTERNAL METHODS #############################

    // returns the offset of the j-th position of a block from its first position
    inline uint64_t offset(const block_t& B, uint64_t j) const
    {
        if (j == 0) return 0;

        if ((B.first_and_flag >> first_pos_bits) != 0) [[unlikely]]
            return wide_offsets[B.offsets[0] * (positions_per_block - 1) + j - 1];

        if (positions_per_block == positions_per_block_16) return reinterpret_cast<const uint16_t*>(B.offsets)[j - 1];
        return reinterpret_cast<const uint32_t*>(B.offsets)[j - 1];
    }

public:
    // ############################# CONSTRUCTORS #############################

    blocked_positions() = default;

    /**
     * @brief builds the sequence
     * @param values the positions in increasing order (each smaller than 2^63)
     */
    blocked_positions(const std::vector<uint64_t>& values)
        : num_positions(values.size())
    {
        if (num_positions == 0) return;
        if (values.back() > mask_first_pos) throw std::invalid_argument("blocked_positions: a position is at least 2^63");

        // returns the size in bytes of the layout with per_block positions per block and offsets of at most max_offset
        auto size_of = [&](uint64_t per_block, uint64_t max_offset) {
            uint64_t num_blocks = 0;
            uint64_t num_wide = 0;

            for (uint64_t first = 0; first < num_positions; first += per_block) {
                num_blocks++;
                num_wide += values[std::min<uint64_t>(num_positions, first + per_block) - 1] - values[first] > max_offset;
            }

            return num_blocks * sizeof(block_t) + num_wide * (per_block - 1) * sizeof(uint64_t);
        };

        // choose the layout that needs less space
        constexpr uint64_t max_offset_16 = std::numeric_limits<uint16_t>::max();
        constexpr uint64_t max_offset_32 = std::numeric_limits<uint32_t>::max();
        positions_per_block = size_of(positions_per_block_32, max_offset_32) < size_of(positions_per_block_16, max_offset_16) ?
            positions_per_block_32 : positions_per_block_16;
        const uint64_t max_offset = positions_per_block == positions_per_block_16 ? max_offset_16 : max_offset_32;
        blocks.resize((num_positions + positions_per_block - 1) / positions_per_block, block_t { 0, { 0, 0, 0 } });

        for (uint64_t b = 0; b < blocks.size(); b++) {
            const uint64_t first = b * positions_per_block;
            const uint64_t last = std::min<uint64_t>(num_positions, first + positions_per_block) - 1;
            block_t& B = blocks[b];
            B.first_and_flag = values[first];

            // the offsets do not fit: store them in wide_offsets
            if (values[last] - values[first] > max_offset) {
                B.first_and_flag |= uint64_t { 1 } << first_pos_bits;
                B.offsets[0] = wide_offsets.size() / (positions_per_block - 1);
                for (uint64_t j = 1; j < positions_per_block; j++)
                    wide_offsets.push_back(first + j <= last ? values[first + j] - values[first] : 0);
                continue;
            }

            for (uint64_t j = 1; first + j <= last; j++) {
                if (positions_per_block == positions_per_block_16) reinterpret_cast<uint16_t*>(B.offsets)[j - 1] = uint16_t(values[first + j] - values[first]);
                else reinterpret_cast<uint32_t*>(B.offsets)[j - 1] = uint32_t(values[first + j] - values[first]);
            }
        }
    }

    // ############################# MISC PUBLIC METHODS #############################

    // returns the number of positions
    inline uint64_t size() const { return num_positions; }

    // returns whether the sequence is empty
    inline bool empty() const { return num_positions == 0; }

    // returns the size in bytes
    uint64_t size_in_bytes() const
    { return blocks.size() * sizeof(block_t) + wide_offsets.size() * sizeof(uint64_t); }

    // fetches the block of the i-th position into the cache
    inline void prefetch(uint64_t i) const { __builtin_prefetch(blocks.data() + i / positions_per_block); }

    // ############################# ACCESS METHODS #############################

    /**
     * @brief returns the i-th position
     * @param i [0..size() - 1]
     */
    inline uint64_t operator[](uint64_t i) const
    {
        const block_t& B = blocks[i / positions_per_block];
        return (B.first_and_flag & mask_first_pos) + offset(B, i % positions_per_block);
    }

    /**
     * @brief returns the i-th and the (i + 1)-th position
     * @param i [0..size() - 2]
     */
    inline std::pair<uint64_t, uint64_t> pair_at(uint64_t i) const
    {
        const uint64_t j = i % positions_per_block;
        const block_t& B = blocks[i / positions_per_block];
        const uint64_t first = (B.first_and_flag & mask_first_pos) + offset(B, j);
        if (j + 1 < positions_per_block) [[likely]] return { first, (B.first_and_flag & mask_first_pos) + offset(B, j + 1) };
        return { first, blocks[i / positions_per_block + 1].first_and_flag & mask_first_pos };
    }

    // ############################# QUERY METHODS #############################

    /**
     * @brief returns the index of a position in the sequence
     * @param value a position
     * @return its index, or size() if it is not in the sequence
     */
    uint64_t find(uint64_t value) const
    {
        if (num_positions == 0) return 0;

        // binary search for the last block whose first position is at most value
        uint64_t lo = 0;
        uint64_t hi = blocks.size();

        while (hi - lo > 1) {
            const uint64_t mid = (lo + hi) / 2;
            if ((blocks[mid].first_and_flag & mask_first_pos) <= value) lo = mid;
            else hi = mid;
        }

        const block_t& B = blocks[lo];
        const uint64_t first = B.first_and_flag & mask_first_pos;
        if (value < first) return num_positions;

        // scan the positions of the block
        for (uint64_t j = 0; j < positions_per_block && lo * positions_per_block + j < num_positions; j++) {
            const uint64_t pos = first + offset(B, j);
            if (pos == value) return lo * positions_per_block + j;
            if (pos > value) break;
        }

        return num_positions;
    }

    // ############################# SERIALIZATION METHODS #############################

    /**
     * @brief writes the sequence to an output stream
     * @param out the output stream
     */
    void serialize(std::ostream& out) const
    {
        out.write((const char*) &num_positions, sizeof(uint64_t));
        out.write((const char*) &positions_per_block, sizeof(uint64_t));
        const uint64_t num_blocks = blocks.size();
        const uint64_t num_wide = wide_offsets.size();
        out.write((const char*) &num_blocks, sizeof(uint64_t));
        out.write((const char*) blocks.data(), num_blocks * sizeof(block_t));
        out.write((const char*) &num_wide, sizeof(uint64_t));
        out.write((const char*) wide_offsets.data(), num_wide * sizeof(uint64_t));
    }

    /**
     * @brief reads a sequence written by serialize
     * @param in the input stream
     */
    void load(std::istream& in)
    {
        in.read((char*) &num_positions, sizeof(uint64_t));
        in.read((char*) &positions_per_block, sizeof(uint64_t));
        uint64_t num_blocks;
        uint64_t num_wide;
        in.read((char*) &num_blocks, sizeof(uint64_t));
        resize_no_init(blocks, num_blocks);
        in.read((char*) blocks.data(), num_blocks * sizeof(block_t));
        in.read((char*) &num_wide, sizeof(uint64_t));
        resize_no_init(wide_offsets, num_wide);
        in.read((char*) wide_offsets.data(), num_wide * sizeof(uint64_t));
    }
};

#endif // MOVE_R_DATA_STRUCTURES_BLOCKED_POSITIONS_HPP
