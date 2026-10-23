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

#ifndef MOVE_R_DATA_STRUCTURES_PRED_SUCC_BLOCK_HPP
#define MOVE_R_DATA_STRUCTURES_PRED_SUCC_BLOCK_HPP

#include <functional>
#include <limits>
#include <optional>

#include <data_structures/interleaved_bit_aligned_vectors.hpp>
#include <misc/utils.hpp>

/**
 * @brief character predecessor/successor queries over an input read through a stored function read(i),
 *        using per-block lookup tables plus a scan of the block containing the query
 *
 * A block spans block_size_factor * sigma positions, so the tables store 2 * input_size / block_size_factor
 * entries, and a query scans at most block_size_factor * sigma positions of the input.
 *
 * @tparam read_fnc_t type of the read function (read(i) returns the symbol at index i)
 */
template <typename read_fnc_t>
class pred_succ_block {
protected:
    // ############################# TYPES AND PARAMETERS #############################

    // default upper bound of succ (no bound)
    static constexpr uint64_t unbounded = std::numeric_limits<uint64_t>::max();

    // ############################# VARIABLES #############################

    mutable std::optional<read_fnc_t> read; // the stored read function (rebindable with set_read)
    uint64_t block_size_factor = 8; // a block spans block_size_factor * sigma input positions
    uint64_t input_size = 0; // number of input positions
    uint64_t sigma = 0; // alphabet size
    uint64_t block_length = 0; // number of input positions per block
    uint64_t num_blocks = 0; // number of blocks

    /* the next occurrence of each symbol at each block start: [0..(num_blocks + 1) * sigma - 1] _succ[b * sigma + c] is
       the first position >= b * block_length with symbol c (input_size if none); bit_width(input_size) bits each */
    interleaved_bit_aligned_vectors _succ;

    /* the previous occurrence of each symbol before each block start: [0..num_blocks * sigma - 1] _pred[b * sigma + c]
       is the last position < b * block_length with symbol c (input_size if none); bit_width(input_size) bits each */
    interleaved_bit_aligned_vectors _pred;

public:
    // ############################# CONSTRUCTORS #############################

    pred_succ_block() = default;

    /**
     * @brief builds the data structure over input[0..input_size-1] read through the function read
     * @param read function returning the symbol at index i
     * @param input_size number of input positions
     * @param alphabet_size alphabet size
     * @param num_threads number of threads to use during the construction
     * @param block_size_factor a block spans block_size_factor * sigma input positions
     */
    pred_succ_block(read_fnc_t read, uint64_t input_size, uint64_t alphabet_size,
        uint64_t num_threads = 1, uint64_t block_size_factor = 8)
    {
        this->input_size = input_size;
        this->sigma = alphabet_size;
        this->block_size_factor = std::max<uint64_t>(1, block_size_factor);
        this->read.emplace(read);
        const read_fnc_t& rd = *this->read;
        uint8_t bits = std::bit_width(input_size); // values in [0, input_size], input_size = no occurrence
        block_length = std::min<uint64_t>(this->block_size_factor * sigma, input_size);
        num_blocks = div_ceil<uint64_t>(input_size, block_length);

        #pragma omp parallel sections num_threads(num_threads)
        {
            // _succ: one pass from right to left
            #pragma omp section
            {
                _succ = interleaved_bit_aligned_vectors({ bits });
                _succ.resize_no_init((num_blocks + 1) * sigma);
                std::vector<uint64_t> next(sigma, input_size);

                for (int64_t block = num_blocks; block >= 0; block--) {
                    int64_t beg = block * (int64_t) block_length;
                    int64_t end = std::min<int64_t>(beg + block_length, input_size);
                    uint64_t first_entry = block * sigma;

                    for (int64_t i = end - 1; i >= beg; i--) next[symbol_bit_pattern(rd(i))] = i;
                    for (uint64_t c = 0; c < sigma; c++) _succ.template set<0, uint64_t>(first_entry + c, next[c]);
                }
            }

            // _pred: one pass from left to right
            #pragma omp section
            {
                _pred = interleaved_bit_aligned_vectors({ bits });
                _pred.resize_no_init(num_blocks * sigma);
                std::vector<uint64_t> prev(sigma, input_size);

                for (int64_t block = 0; block < (int64_t) num_blocks; block++) {
                    int64_t beg = block * (int64_t) block_length;
                    int64_t end = std::min<int64_t>(beg + block_length, input_size);
                    uint64_t first_entry = block * sigma;

                    for (uint64_t c = 0; c < sigma; c++) _pred.template set<0, uint64_t>(first_entry + c, prev[c]);
                    for (int64_t i = beg; i < end; i++) prev[symbol_bit_pattern(rd(i))] = i;
                }
            }
        }
    }

    // ############################# MISC PUBLIC METHODS #############################

    // (re-)binds the read function; call after loading, or after the object read references has moved
    inline void set_read(read_fnc_t read) const { this->read.emplace(read); }

    // returns the number of input positions per block
    inline uint64_t block_size() const { return block_length; }

    // returns whether the data structure is empty (i.e., has not been built)
    inline bool empty() const { return _succ.empty() && _pred.empty(); }

    // returns the size of the data structure in bytes
    uint64_t size_in_bytes() const { return _succ.size_in_bytes() + _pred.size_in_bytes(); }

    // ############################# ACCESS METHODS #############################

    // returns the first index >= block * block_size() whose symbol is sym (input_size if none)
    inline uint64_t succ_of_block(uint64_t block, uint64_t sym) const { return _succ.template get<0, uint64_t>(block * sigma + sym); }

    // returns the last index < block * block_size() whose symbol is sym (input_size if none)
    inline uint64_t pred_of_block(uint64_t block, uint64_t sym) const { return _pred.template get<0, uint64_t>(block * sigma + sym); }

    // ############################# QUERY METHODS #############################

    /**
     * @brief returns the smallest index in [x, max] whose symbol is sym (a value > max if none)
     * @param sym a symbol
     * @param x lower bound
     * @param max upper bound
     */
    inline uint64_t succ(uint64_t sym, uint64_t x, uint64_t max = unbounded) const
    {
        const read_fnc_t& rd = *read;
        // scan up to the next block start, then look the result up in _succ
        uint64_t block = div_ceil<uint64_t>(x, block_length);
        uint64_t max_x = std::min<uint64_t>({block * block_length, max, input_size - 1});
        while (x <= max_x && sym != symbol_bit_pattern(rd(x))) x++;

        if (x > max_x) [[likely]] {
            if (x > max) [[unlikely]] return x;
            x = succ_of_block(block, sym);
        }

        return x;
    }

    /**
     * @brief returns the largest index in [min, x] whose symbol is sym (sym must occur in [min, x]; if it does not occur
     *        in [0, x], input_size is returned)
     * @param sym a symbol
     * @param x upper bound
     * @param min lower bound
     */
    inline uint64_t pred(uint64_t sym, uint64_t x, uint64_t min = 0) const
    {
        const read_fnc_t& rd = *read;
        // scan down to the start of the block of x, then look the result up in _pred
        uint64_t block = x / block_length;
        uint64_t min_x = std::max<uint64_t>(block * block_length, min);
        while (x > min_x && sym != symbol_bit_pattern(rd(x))) x--;
        if (sym != symbol_bit_pattern(rd(x))) [[likely]] x = pred_of_block(block, sym);
        return x;
    }

    // ############################# SERIALIZATION METHODS #############################

    /**
     * @brief writes the data structure to an output stream (the read function is not serialized)
     * @param out output stream
     */
    void serialize(std::ostream& out) const
    {
        out.write((char*) &input_size, sizeof(uint64_t));
        out.write((char*) &sigma, sizeof(uint64_t));
        out.write((char*) &block_size_factor, sizeof(uint64_t));
        out.write((char*) &block_length, sizeof(uint64_t));
        out.write((char*) &num_blocks, sizeof(uint64_t));
        _succ.serialize(out);
        _pred.serialize(out);
    }

    // reads the data structure from an input stream; set_read has to be called afterwards
    void load(std::istream& in)
    {
        in.read((char*) &input_size, sizeof(uint64_t));
        in.read((char*) &sigma, sizeof(uint64_t));
        in.read((char*) &block_size_factor, sizeof(uint64_t));
        in.read((char*) &block_length, sizeof(uint64_t));
        in.read((char*) &num_blocks, sizeof(uint64_t));
        _succ.load(in);
        _pred.load(in);
    }
};

#endif // MOVE_R_DATA_STRUCTURES_PRED_SUCC_BLOCK_HPP
