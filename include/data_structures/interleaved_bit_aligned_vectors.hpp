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

#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <vector>

#include <misc/utils.hpp>
#include <misc/files.hpp>

/**
 * @brief atomically overwrites the bits selected by submask of the 64-bit word with the
 *        corresponding bits of subval (all other bits of the word are preserved); this makes
 *        concurrent writes to disjoint bit-fields that share a 64-bit word race-free
 * @param word reference to the aligned 64-bit word to update
 * @param submask mask selecting the bits to overwrite
 * @param subval new values for the selected bits (must be zero outside submask)
 */
static inline void atomic_store_bits(uint64_t& word, uint64_t submask, uint64_t subval)
{
    uint64_t expected = __atomic_load_n(&word, __ATOMIC_RELAXED);
    uint64_t desired;

    do {
        desired = (expected & ~submask) | subval;
    } while (!__atomic_compare_exchange_n(&word, &expected, desired,
        true, __ATOMIC_RELAXED, __ATOMIC_RELAXED));
}

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

    using wide_t = constexpr_switch_t<
        constexpr_case<sizeof(val_t) == 1,    uint16_t>,
        constexpr_case<sizeof(val_t) == 2,    uint32_t>,
        constexpr_case<sizeof(val_t) == 4,    uint64_t>,
     /* constexpr_case<sizeof(val_t) == 8, */ __uint128_t>;

    uint64_t size_vectors = 0; // size of each stored vector
    uint64_t capacity_vectors = 0; // capacity of each stored vector
    uint64_t width_entry = 0; // sum of the widths (in bits) of all vectors

    // [0 .. byte_size(capacity_vectors)] vector of bytes storing the interleaved vectors
    std::vector<char> data_vectors;

    // [0 .. num_vectors - 1] widths of the stored vectors; widths[i] = width (in bits) of vector i
    std::array<uint64_t, num_vectors> widths {};

    // [0 .. num_vectors - 1] offsets[i] = offset (in bits) of the ith vectors entry within the region storing all ith entries
    std::array<uint64_t, num_vectors> offsets {};

    // [0 .. num_vectors - 1] masks that are used to mask off data of other vector entries when accessing a vector
    std::array<wide_t, num_vectors> masks {};

    // padding (in bytes) appended to data_vectors so that the wide (block) read/write and the aligned
    // 64-bit atomic-word writes (set_parallel) of the last entry stay in bounds (independent of size)
    static constexpr uint64_t padding = std::max<uint64_t>(sizeof(wide_t), sizeof(uint64_t));

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
            if (widths[vec_idx] == 0) {
                this->widths[vec_idx] = 0;
                offsets[vec_idx] = 0;
                masks[vec_idx] = 0;
            } else {
                offsets[vec_idx] = width_entry;
                this->widths[vec_idx] = widths[vec_idx];
                width_entry += widths[vec_idx];
                masks[vec_idx] = (wide_t{1} << widths[vec_idx]) - 1;
            }
        }
    }

    /**
     * @brief returns the combined size (in bytes) of all vectors for a given vector size
     * @param size number of entries per vector
     * @return the combined size (in bytes) of all vectors
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
     * @tparam vec_idx vector index
     * @return its witdth in bits
     */
    template <uint8_t vec_idx>
    inline uint64_t width() const
    {
        return widths[vec_idx];
    }

    /**
     * @brief reads the wide_t-sized block of data starting at byte byte_idx
     * @param byte_idx byte index into the interleaved data
     * @return the wide_t-sized block of data starting at byte byte_idx
     */
    inline wide_t read_block(uint64_t byte_idx) const
    {
        return *reinterpret_cast<const wide_t*>(&data_vectors[byte_idx]);
    }

    /**
     * @brief returns a writable reference to the wide_t-sized block of data starting at byte byte_idx
     * @param byte_idx byte index into the interleaved data
     * @return reference to the wide_t-sized block of data starting at byte byte_idx
     */
    inline wide_t& block(uint64_t byte_idx)
    {
        return *reinterpret_cast<wide_t*>(&data_vectors[byte_idx]);
    }

    /**
     * @brief returns a pointer of type T to the ith byte in the data of the interleved vectors
     * @tparam T type
     * @param i byte index into the interleaved data
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
        return get<0>(i);
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

            // pad by padding bytes so the wide (block) and aligned 64-bit (set_parallel) read/write of
            // the last entry stays in bounds
            std::vector<char> new_data_vectors;
            no_init_resize(new_data_vectors, num_new_bytes + padding);

            for (uint64_t i = num_old_bytes; i < num_new_bytes + padding; i++) {
                new_data_vectors[i] = 0;
            }

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
            data_vectors.resize(byte_size(capacity_vectors) + padding);
            data_vectors.shrink_to_fit();
        }
    }

    using block_info_t = std::tuple<
        uint64_t, // byte_idx
        uint64_t // bit_offs
    >;

    protected:

    /**
     * @brief returns the block index and bit offset of the i-th entry of the vector with index vec_idx within its block
     * @tparam vec_idx vector index (0 <= vec_idx < num_vectors)
     * @param i entry index (0 <= i < size_vectors)
     * @return the block index and bit offset of the i-th entry of the vector with index vec_idx within its block
     */
    template <uint8_t vec_idx>
    inline block_info_t block_info(uint64_t i) const
    {
        static_assert(vec_idx < num_vectors);
        uint64_t bit_pos = i * width_entry + offsets[vec_idx];
        uint64_t byte_idx = bit_pos / 8;
        uint64_t bit_offs = bit_pos % 8;
        return {byte_idx, bit_offs};
    }

    public:
    /**
     * @brief sets the i-th entry in the vector with index vec_idx to v
     * @tparam vec_idx vector index (0 <= vec_idx < num_vectors)
     * @tparam T type of the value to store
     * @param i entry index (0 <= i < size_vectors)
     * @param v value to store
     */
    template <uint8_t vec_idx, typename T = val_t>
    inline void set(uint64_t i, T v)
    {
        auto [byte_idx, bit_offs] = block_info<vec_idx>(i);
        block(byte_idx) = (read_block(byte_idx) & ~(masks[vec_idx] << bit_offs)) | ((wide_t(v) & masks[vec_idx]) << bit_offs);
    }

    /**
     * @brief alias for set(); provided for API-compatibility with interleaved_byte_aligned_vectors
     *        (bit-packed entries are always written safely, so there is no separate unsafe path)
     * @tparam vec_idx vector index (0 <= vec_idx < num_vectors)
     * @tparam T type of the value to store
     * @param i entry index (0 <= i < size_vectors)
     * @param v value to store
     */
    template <uint8_t vec_idx, typename T = val_t>
    inline void set_unsafe(uint64_t i, T v)
    {
        set<vec_idx, T>(i, v);
    }

    /**
     * @brief sets the i-th entry in the vector with index vec_idx to v; unlike set(), this may be
     *        called concurrently for distinct entries (or distinct vectors of the same entry), even
     *        when their bit-fields share a byte. Concurrent writes to the same (i, vec_idx) field are 
     *        not supported.
     * @tparam vec_idx vector index (0 <= vec_idx < num_vectors)
     * @tparam T type of the value to store
     * @param i entry index (0 <= i < size_vectors)
     * @param v value to store
     */
    template <uint8_t vec_idx, typename T = val_t>
    inline void set_parallel(uint64_t i, T v)
    {
        static_assert(vec_idx < num_vectors);

        // a field of <= 64 bits spans at most two 64-bit words; write each word's portion
        // atomically so concurrent writes to neighbouring fields sharing a word do not race
        uint64_t bit_pos = i * width_entry + offsets[vec_idx];
        uint64_t width = widths[vec_idx];
        uint64_t val = uint64_t(val_t(v)) & uint64_t(masks[vec_idx]);

        uint64_t* words = reinterpret_cast<uint64_t*>(data_vectors.data());
        uint64_t word_idx = bit_pos >> 6; // bit_pos / 64
        uint64_t bit_offs = bit_pos & 63; // bit_pos % 64
        uint64_t bits_lo = std::min<uint64_t>(64 - bit_offs, width); // bits stored in the first word

        uint64_t mask_lo = (bits_lo == 64 ? ~uint64_t{0} : ((uint64_t{1} << bits_lo) - 1)) << bit_offs;
        atomic_store_bits(words[word_idx], mask_lo, (val << bit_offs) & mask_lo);

        if (bits_lo < width) { // the field continues in the next 64-bit word
            uint64_t mask_hi = (uint64_t{1} << (width - bits_lo)) - 1;
            atomic_store_bits(words[word_idx + 1], mask_hi, (val >> bits_lo) & mask_hi);
        }
    }

    /**
     * @brief returns the i-th entry in the vector with index vec_idx
     * @tparam vec_idx vector index (0 <= vec_idx < num_vectors)
     * @tparam T type to return the entry as
     * @param i entry index (0 <= i < size_vectors)
     * @return the i-th entry in the vector with index vec_idx
     */
    template <uint8_t vec_idx, typename T = val_t>
    inline T get(uint64_t i) const
    {
        auto [byte_idx, bit_offs] = block_info<vec_idx>(i);
        return T((read_block(byte_idx) >> bit_offs) & masks[vec_idx]);
    }

    /**
     * @brief alias for get(); provided for API-compatibility with interleaved_byte_aligned_vectors
     * @tparam vec_idx vector index (0 <= vec_idx < num_vectors)
     * @tparam T type to return the entry as
     * @param i entry index (0 <= i < size_vectors)
     * @return the i-th entry in the vector with index vec_idx
     */
    template <uint8_t vec_idx, typename T = val_t>
    inline T get_unsafe(uint64_t i) const
    {
        return get<vec_idx, T>(i);
    }

    /**
     * @brief appends a tuple of values to the end of the interleaved vectors
     * @param vals tuple of values
     */
    inline void emplace_back(std::array<val_t, num_vectors> vals)
    {
        if (size_vectors == capacity_vectors) {
            reserve(size_vectors + std::max<uint64_t>(1, size_vectors / 2));
        }

        for_constexpr<0, num_vectors, 1>([&](auto vec_idx) {
            set<vec_idx>(size_vectors, vals[vec_idx]);
        });

        size_vectors++;
    }

    /**
     * @brief appends a single value to the end of the interleaved vectors, setting only the vector
     *        with index vec_idx (the other vectors of the new entry are left zero-initialized)
     * @tparam vec_idx vector index (0 <= vec_idx < num_vectors)
     * @tparam T type of the value to append
     * @param val value to append
     */
    template <uint8_t vec_idx = 0, typename T = val_t>
    inline void emplace_back(T val)
    {
        if (size_vectors == capacity_vectors) {
            reserve(size_vectors + std::max<uint64_t>(1, size_vectors / 2));
        }

        set<vec_idx, T>(size_vectors, val);
        size_vectors++;
    }

    /**
     * @brief logs the contents of all vectos
     */
    void log_contents() const
    {
        uint8_t last_used_vec = 0;
        uint8_t num_used_vecs = 0;

        for_constexpr<0, num_vectors, 1>([&](auto vec_idx){
            if (widths[vec_idx] != 0) {
                last_used_vec = vec_idx;
                num_used_vecs++;
            }
        });

        std::string opening_bracket_str = num_used_vecs <= 1 ? "" : "(";
        std::string closing_bracket_str = num_used_vecs <= 1 ? "" : ")";

        for (uint64_t i = 0; i < size_vectors; i++) {
            std::cout << opening_bracket_str;

            for_constexpr<0, num_vectors, 1>([&](auto vec_idx){
                if (widths[vec_idx] != 0) {
                    std::cout << get<vec_idx>(i) << (vec_idx == last_used_vec ? "" : " ");
                }
            });

            std::cout << closing_bracket_str;
            if (i != size_vectors - 1) std::cout << " ";
        }
        
        std::cout << std::endl;
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
            // write only the bytes actually holding entries (capacity and padding are scratch)
            write_to_file(out, data_vectors.data(), byte_size(size_vectors));
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
            read_from_file(in, data_vectors.data(), byte_size(old_size));
        }
    }

    /**
     * @brief serializes the interleaved vectors to an output stream
     * @param os output stream
     * @return the output stream
     */
    std::ostream& operator>>(std::ostream& os) const
    {
        serialize(os);
        return os;
    }

    /**
     * @brief loads the interleaved vectors from an input stream
     * @param is input stream
     * @return the input stream
     */
    std::istream& operator<<(std::istream& is)
    {
        load(is);
        return is;
    }
};