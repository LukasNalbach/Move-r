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

#include <array>
#include <cstring>
#include <iostream>
#include <vector>

#include <misc/utils.hpp>

/**
 * @brief variable-width word-packed interleaved vectors (the maximum supported width is 64 bits)
 * @tparam val_t type of the values stored in the vectors
 * @tparam num_vectors maximum number of vectors that can be stored
 */
template <typename val_t, uint8_t num_vectors = 8>
class interleaved_bit_aligned_vectors {
protected:
    static_assert(num_vectors > 0);
    static_assert(
        std::is_same_v<val_t, uint8_t> ||
        std::is_same_v<val_t, uint16_t> ||
        std::is_same_v<val_t, uint32_t> ||
        std::is_same_v<val_t, uint64_t>);

    // size (in bits) of the blocks that are read/written at once
    static constexpr val_t data_block_bit_width = 8 * sizeof(val_t);

    // mask that is use to quickly compute the left padding within a block when accessing/writing to a block
    static constexpr uint64_t left_padding_mask = data_block_bit_width - 1;

    uint64_t size_vectors = 0; // size of each stored vector
    uint64_t capacity_vectors = 0; // capacity of each stored vector
    uint64_t width_entry = 0; // sum of the widths (in bits) of all vectors

    // [0 .. byte_size(capacity_vectors)] vector of bytes storing the interleaved vectors
    std::vector<char> data_vectors;

    // [0 .. num_vectors - 1] widths of the stored vectors; widths[i] = width (in bits) of vector i
    std::array<uint64_t, num_vectors> widths;

    // [0 .. num_vectors - 1] offsets[i] = offset (in bits) of the ith vectors entry within the region storing all ith entries
    std::array<uint64_t, num_vectors> offsets;
    
    // offset info
    struct offset_info_t {
        val_t left_excess; // the number of bits to the right of the data in the left block
        val_t right_excess; // the number of bits of the data in the right block
        val_t left_mask; // mask to mask of bits left and right of the data in the left block
    };

    /**
     * @brief [0 .. num_vectors - 1][0 .. data_block_bit_width - 1] stores at position [vec_idx][o] the offset info for 
     * accessing/writing to the vector with index vec_idx within a block with the offset o
     */
    std::array<std::array<offset_info_t, data_block_bit_width>, num_vectors> offset_info;

    /**
     * @brief initializes the interleaved_bit_aligned_vectors with the vector-widths stored in widths
     * @param widths vector containing the widths (in bits) of the interleaved arrays
     */
    void initialize(std::array<uint64_t, num_vectors> widths)
    {
        size_vectors = 0;
        capacity_vectors = 0;
        width_entry = 0;

        data_vectors.clear();
        data_vectors.shrink_to_fit();

        for (uint8_t vec_idx = 0; vec_idx < num_vectors; vec_idx++) {
            offsets[vec_idx] = width_entry;

            if (widths[vec_idx] == 0) {
                this->widths[vec_idx] = 0;
            } else {
                this->widths[vec_idx] = widths[vec_idx];
                width_entry += widths[vec_idx];

                for (uint8_t left_padding = 0; left_padding < data_block_bit_width; left_padding++) {
                    offset_info_t& info = offset_info[vec_idx][left_padding];

                    int64_t diff = int64_t{left_padding + widths[vec_idx]} - int64_t{data_block_bit_width};
                    info.left_excess = std::max<int64_t>(-diff, 0);
                    info.right_excess = std::max<int64_t>(diff, 0);

                    info.left_mask = std::numeric_limits<val_t>::max() >> left_padding;
                    info.left_mask &= std::numeric_limits<val_t>::max() << info.left_excess;
                    info.left_mask = ~info.left_mask;
                }
            }
        }
    }

    /**
     * @brief returns the combined size (in bytes) of all vectors for a given vector size
     */
    uint64_t byte_size(uint64_t size) const
    {
        return sizeof(val_t) * div_ceil(size * width_entry, 8 * sizeof(val_t));
    }

public:
    interleaved_bit_aligned_vectors() = default;

    /**
     * @brief Construct a new interleaved_bit_aligned_vectors
     * @param widths vector containing the widths (in bits) of the interleaved arrays
     */
    interleaved_bit_aligned_vectors(std::array<uint64_t, num_vectors> widths)
    {
        initialize(widths);
    }

    /**
     * @brief returns the size of each stored vector
     * @return the size of each stored vector
     */
    inline uint64_t size() const
    {
        return size_vectors;
    }

    /**
     * @brief returns whether the interleaved vectors are empty
     * @return whether the interleaved vectors are empty
     */
    inline bool empty() const
    {
        return size_vectors == 0;
    }

    /**
     * @brief returns the size of the data structure in bits
     * @return size of the data structure in bits
     */
    uint64_t size_in_bytes() const
    {
        return sizeof(this) + byte_size(size_vectors);
    }

    /**
     * @brief returns total width (number of bits) per entry, that is the (sum of all widths)
     * @return number of bits per entry
     */
    inline uint64_t bits_per_entry() const
    {
        return width_entry;
    }

    /**
     * @brief returns the width in bits of the vector with index vec_idx
     * @param vec_idx vector index
     * @return its witdth in bits
     */
    template <uint8_t vec_idx>
    inline uint64_t width() const
    {
        return widths[vec_idx];
    }

    inline val_t read_block(uint64_t blk_idx) const
    {
        return reinterpret_cast<const val_t*>(data_vectors.data())[blk_idx];
    }

    inline val_t& block(uint64_t blk_idx)
    {
        return reinterpret_cast<val_t*>(data_vectors.data())[blk_idx];
    }

    /**
     * @brief returns a pointer of type T to the ith byte in the data of the interleved vectors
     * @tparam T type
     * @return pointer of type T to the ith byte in the data of the interleved vectors
     */
    template <typename T = char>
    inline T* data(uint64_t i = 0)
    {
        return reinterpret_cast<T*>(data_vectors.data() + i);
    }

    /**
     * @brief returns the i-th entry of the vector with index 0
     * @param i entry index (0 <= i < size_vectors)
     * @return i-th entry of the vector with index 0
     */
    inline val_t operator[](uint64_t i) const
    {
        return get<0, val_t>(i);
    }

    /**
     * @brief reserves capacity entries in all stored vectors; if capacity is smaller than the
     *        current capacity, nothing happens (use shrink_to_fit() to lower the capacity)
     * @param capacity capacity
     * @param num_threads number of threads to use when copying entries (default: 1)
     */
    void reserve(uint64_t capacity, uint16_t num_threads = 1)
    {
        if (capacity_vectors < capacity) {
            uint64_t num_old_bytes = byte_size(capacity_vectors);
            uint64_t num_new_bytes = byte_size(capacity);

            std::vector<char> new_data_vectors;
            no_init_resize(new_data_vectors, num_new_bytes);
            uint64_t last_block_idx = num_new_bytes / sizeof(val_t) - 1;
            reinterpret_cast<val_t*>(new_data_vectors.data())[last_block_idx] = 0;

            #pragma omp parallel for num_threads(num_threads)
            for (uint64_t i = 0; i < num_old_bytes; i++) {
                new_data_vectors[i] = data_vectors[i];
            }

            std::swap(data_vectors, new_data_vectors);
            capacity_vectors = capacity;
        }
    }

    /**
     * @brief resizes all stored vectors to size and initializes new entries to 0
     * @param size size
     * @param num_threads number of threads to use when copying entries and initializing new entries
     *                    to 0 (default: 1)
     */
    void resize(uint64_t size, uint16_t num_threads = 1)
    {
        if (capacity_vectors < size) {
            reserve(size, num_threads);
        }

        uint64_t num_old_bytes = byte_size(size_vectors);
        uint64_t num_new_bytes = byte_size(size);

        #pragma omp parallel for num_threads(num_threads)
        for (uint64_t i = num_old_bytes; i < num_new_bytes; i++) {
            data_vectors[i] = 0;
        }
        
        size_vectors = size;
    }

    /**
     * @brief resizes all stored vectors to size without initializing new entries to 0
     * @param size size
     * @param num_threads number of threads to use when copying entries (default: 1)
     */
    void resize_no_init(uint64_t size, uint16_t num_threads = 1)
    {
        if (capacity_vectors < size) {
            reserve(size, num_threads);
        }

        size_vectors = size;
    }

    /**
     * @brief resizes all stored vectors to size 0
     */
    inline void clear()
    {
        resize(0);
    }

    /**
     * @brief shrinks all stored vectors to their size
     */
    void shrink_to_fit()
    {
        if (size_vectors < capacity_vectors) {
            capacity_vectors = std::max<uint64_t>(2, size_vectors);
            data_vectors.resize(byte_size(capacity_vectors));
            data_vectors.shrink_to_fit();
        }
    }

    /**
     * @brief sets the i-th entry in the vector with index vec_idx to v
     * @tparam vec_idx vector index (0 <= vec_idx < num_vectors)
     * @param i entry index (0 <= i < size_vectors)
     * @param v value to store
     */
    template <uint8_t vec_idx>
    inline void set(uint64_t i, val_t v)
    {
        static_assert(vec_idx < num_vectors);

        uint64_t bit_pos = i * width_entry + offsets[vec_idx];
        uint64_t left_block = bit_pos / data_block_bit_width;
        uint64_t right_block = left_block + 1;
        uint64_t left_padding = bit_pos & left_padding_mask;
        const offset_info_t& info = offset_info[vec_idx][left_padding];

        if (info.right_excess == 0) {
            block(left_block) = (read_block(left_block) & info.left_mask) | (v << info.left_excess);
        } else {
            val_t right_mask = std::numeric_limits<val_t>::max() >> info.right_excess;
            val_t right_shift = data_block_bit_width - info.right_excess;
            block(left_block) = (read_block(left_block) & info.left_mask) | (v >> info.right_excess);
            block(right_block) = (read_block(right_block) & right_mask) | (v << right_shift);
        }
    }

    /**
     * @brief returns the i-th entry in the vector with index vec_idx
     * @tparam vec_idx vector index (0 <= vec_idx < num_vectors)
     * @param i entry index (0 <= i < size_vectors)
     * @return the i-th entry in the vector with index vec_idx
     */
    template <uint8_t vec_idx>
    inline val_t get(uint64_t i) const
    {
        static_assert(vec_idx < num_vectors);

        uint64_t bit_pos = i * width_entry + offsets[vec_idx];
        uint64_t left_block = bit_pos / data_block_bit_width;
        uint64_t right_block = left_block + 1;
        uint64_t left_padding = bit_pos & left_padding_mask;
        const offset_info_t& info = offset_info[vec_idx][left_padding];
        val_t left_mask = std::numeric_limits<val_t>::max() >> left_padding;

        if (info.right_excess == 0) {
            return (read_block(left_block) & left_mask) >> info.left_excess;
        } else {
            return ((read_block(left_block) & left_mask) << info.right_excess) |
                    (read_block(right_block) >> (data_block_bit_width - info.right_excess));
        }
    }

    /**
     * @brief appends a tuple of values to the end of the interleaved vectors
     * @param vals tuple of values
     */
    inline void emplace_back(std::array<val_t, num_vectors> vals)
    {
        if (size_vectors == capacity_vectors) {
            reserve(1.5 * size_vectors);
        }

        for_constexpr<0, num_vectors, 1>([&](auto vec_idx) {
            set<vec_idx>(size_vectors, vals[vec_idx]);
        });

        size_vectors++;
    }

    /**
     * @brief serializes the interleaved vectors to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) const
    {
        out.write((char*) &size_vectors, sizeof(uint64_t));
        out.write((char*) widths.data(), num_vectors * sizeof(uint64_t));

        if (size_vectors > 0) {
            write_to_file(out, data_vectors.data(), data_vectors.size());
        }
    }

    /**
     * @brief loads the interleaved vectors from an input stream
     * @param in input stream
     */
    void load(std::istream& in)
    {
        uint64_t old_size;
        in.read((char*) &old_size, sizeof(uint64_t));
        in.read((char*) widths.data(), num_vectors * sizeof(uint64_t));
        initialize(widths);

        if (old_size > 0) {
            resize_no_init(old_size);
            read_from_file(in, data_vectors.data(), data_vectors.size());
        }
    }

    std::ostream& operator>>(std::ostream& os) const
    {
        serialize(os);
        return os;
    }

    std::istream& operator<<(std::istream& is)
    {
        load(is);
        return is;
    }
};