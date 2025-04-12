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
 * @brief variable-width interleaved vectors (widths are fixed to whole bytes)
 * @tparam val_t default type of the values stored in the vectors
 * @tparam pos_t unsigned integer position (vector index) type
 * @tparam num_vectors maximum number of vectors that can be stored
 */
template <typename val_t, typename pos_t = uint32_t, uint8_t num_vectors = 8>
class interleaved_byte_aligned_vectors {

    static_assert(std::is_same_v<pos_t, uint32_t> || std::is_same_v<pos_t, uint64_t>);

    static_assert(
        std::is_same_v<val_t, char> ||
        std::is_same_v<val_t, uint8_t> ||
        std::is_same_v<val_t, uint16_t> ||
        std::is_same_v<val_t, uint32_t> ||
        std::is_same_v<val_t, uint64_t>
    );

    static_assert(num_vectors > 0);

protected:
    static constexpr uint64_t block_size = sizeof(uint128_t); // size (in bytes) of the blocks that are copied at once

    uint64_t size_vectors = 0; // size of each stored vector
    uint64_t capacity_vectors = 0; // capacity of each stored vector
    uint64_t width_entry = 0; // sum of the widths of all vectors

    // [0 .. (capacity_vectors + 1) * width_entry - 1] vector of bytes storing the interleaved vectors
    std::vector<char> data_vectors;

    // [0 .. num_vectors - 1] or shorter; widths of the stored vectors; widths[i] = width of vector i
    std::array<uint64_t, num_vectors> widths;

    /** @brief [0 .. num_vectors - 1] pointers to the first entries of each vector; bases[vec_idx] = base
     *        of the first entry of the vector with index vec_idx */
    std::array<char*, num_vectors> bases;

    /** @brief [0 .. num_vectors - 1] masks that are used to mask off data of other vector entries when
     *        reading from a vector */
    std::array<val_t, num_vectors> masks_get;

    /** @brief [0 .. num_vectors - 1] masks that are used to mask off data of other vector entries when
     *        writing to a vector */
    std::array<val_t, num_vectors> masks_set;

    /**
     * @brief initializes the interleaved_byte_aligned_vectors with the vector-widths stored in widths
     * @param widths vector containing the widths (in bytes) of the interleaved arrays
     */
    void initialize(std::array<uint8_t, num_vectors> widths = { sizeof(val_t) })
    {
        size_vectors = 0;
        capacity_vectors = 0;
        width_entry = 0;

        data_vectors.clear();
        data_vectors.shrink_to_fit();

        for (uint8_t i = 0; i < num_vectors; i++) {
            if (widths[i] == 0) {
                masks_get[i] = 0;
                masks_set[i] = 0;
                this->widths[i] = 0;
            } else {
                this->widths[i] = widths[i];
                width_entry += widths[i];
                masks_get[i] = std::numeric_limits<val_t>::max() >> (8 * (sizeof(val_t) - widths[i]));
                masks_set[i] = ~masks_get[i];
            }
        }

        reserve(2);
    }

    /**
     * @brief builds the bases array
     * @param data pointer to the data storing the interleaved vectors
     */
    void set_bases(char* data)
    {
        bases[0] = data;

        for (uint8_t i = 1; i < num_vectors; i++) {
            if (widths[i] == 0) {
                bases[i] = nullptr;
            } else {
                bases[i] = bases[i - 1] + widths[i - 1];
            }
        }
    }

    /**
     * @brief copies another interleaved_byte_aligned_vectors object into this object
     * @param other another interleaved_byte_aligned_vectors object
     */
    void copy_from_other(const interleaved_byte_aligned_vectors& other)
    {
        size_vectors = other.size_vectors;
        capacity_vectors = other.capacity_vectors;
        width_entry = other.width_entry;
        data_vectors = other.data_vectors;
        widths = other.widths;
        masks_get = other.masks_get;
        masks_set = other.masks_set;

        set_bases(&data_vectors[0]);
    }

    /**
     * @brief moves another interleaved_byte_aligned_vectors object into this object
     * @param other another interleaved_byte_aligned_vectors object
     */
    void move_from_other(interleaved_byte_aligned_vectors&& other)
    {
        size_vectors = other.size_vectors;
        capacity_vectors = other.capacity_vectors;
        width_entry = other.width_entry;

        data_vectors = std::move(other.data_vectors);
        widths = other.widths;
        bases = other.bases;
        masks_get = other.masks_get;
        masks_set = other.masks_set;

        other.initialize();
    }

public:
    interleaved_byte_aligned_vectors() { initialize(); }
    interleaved_byte_aligned_vectors(interleaved_byte_aligned_vectors&& other) { move_from_other(std::move(other)); }
    interleaved_byte_aligned_vectors(const interleaved_byte_aligned_vectors& other) { copy_from_other(other); }
    interleaved_byte_aligned_vectors& operator=(interleaved_byte_aligned_vectors&& other) { move_from_other(std::move(other)); return *this; }
    interleaved_byte_aligned_vectors& operator=(const interleaved_byte_aligned_vectors& other) { copy_from_other(other); return *this; }

    ~interleaved_byte_aligned_vectors()
    {
        for (uint8_t i = 0; i < num_vectors; i++) {
            bases[i] = nullptr;
        }
    }

    /**
     * @brief Construct a new interleaved_byte_aligned_vectors
     * @param widths vector containing the widths (in bytes) of the interleaved arrays
     */
    interleaved_byte_aligned_vectors(std::array<uint8_t, num_vectors> widths)
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
     * @brief returns the size of the data structure in bytes
     * @return size of the data structure in bytes
     */
    uint64_t size_in_bytes() const
    {
        return sizeof(this) + size_vectors * width_entry;
    }

    /**
     * @brief returns total width (number of bytes) per entry, that is the (sum of all widths)
     * @return number of bytes per entry
     */
    inline uint8_t bytes_per_entry() const
    {
        return this->width_entry;
    }

    /**
     * @brief returns the width in bytes of the vector with index vec_idx
     * @param vec_idx vector index
     * @return its witdth in bytes
     */
    template <uint8_t vec_idx>
    inline uint8_t width() const
    {
        return widths[vec_idx];
    }

    /**
     * @brief returns a pointer to the data of the interleved vectors
     * @return pointer to the data of the interleved vectors
     */
    inline char* data() const
    {
        return bases[0];
    }

    /**
     * @brief returns the i-th entry of the vector with index 0
     * @param i entry index (0 <= i < size_vectors)
     * @return i-th entry of the vector with index 0
     */
    inline val_t operator[](pos_t i) const
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
            std::vector<char> new_data_vectors;
            no_init_resize(new_data_vectors, capacity * width_entry + block_size);

            #pragma omp parallel for num_threads(num_threads)
            for (uint64_t i = 0; i < size_vectors * width_entry; i++) {
                new_data_vectors[i] = data_vectors[i];
            }

            std::memset(&new_data_vectors[size_vectors * width_entry], 0, block_size);
            std::memset(&new_data_vectors[capacity * width_entry], 0, block_size);
            std::swap(data_vectors, new_data_vectors);
            set_bases(&data_vectors[0]);
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
        uint64_t size_vectors_old = size_vectors;

        if (capacity_vectors < size) {
            reserve(size, num_threads);
        }

        #pragma omp parallel for num_threads(num_threads)
        for (uint64_t i = size_vectors_old * width_entry; i < size * width_entry; i++) {
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
            data_vectors.resize(capacity_vectors * width_entry + block_size);
            data_vectors.shrink_to_fit();
            set_bases(&data_vectors[0]);
        }
    }

    /**
     * @brief sets the i-th entry in the vector with index vec_idx to v (cannot be executed in parallel for multiple
     * different positions in the vectors)
     * @tparam vec_idx vector index (0 <= vec_idx < num_vectors)
     * @param i entry index (0 <= i < size_vectors)
     * @param v value to store
     */
    template <uint8_t vec_idx, typename T = val_t>
    inline void set(pos_t i, T v)
    {
        static_assert(vec_idx < num_vectors);

        T* loc = reinterpret_cast<T*>(bases[vec_idx] + i * width_entry);
        *loc = v | (*loc & masks_set[vec_idx]);
    }

    /**
     * @brief sets the i-th entry in the vector with index vec_idx to v (can be executed in parallel for multiple
     * different positions in the vectors)
     * @tparam vec_idx vector index (0 <= vec_idx < num_vectors)
     * @param i entry index (0 <= i < size_vectors)
     * @param v value to store
     */
    template <uint8_t vec_idx, typename T = val_t>
    inline void set_parallel(pos_t i, T v)
    {
        static_assert(vec_idx < num_vectors);

        for (uint64_t byte = 0; byte < widths[vec_idx]; byte++) {
            *reinterpret_cast<char*>(bases[vec_idx] + i * width_entry + byte) = *(reinterpret_cast<char*>(&v) + byte);
        }
    }

    /**
     * @brief plainly stores the value v (of type T) at the memory location of the i-th entry in the vector
     * with index vec_idx (if sizeof(T) > widths[vec_idx], this is unsafe, since subsequnet entries are overwritten)
     * @tparam vec_idx vector index (0 <= vec_idx < num_vectors)
     * @tparam T type to treat the i-th entry in the vector with index vec_idx as
     * @param i entry index (0 <= i < size_vectors)
     * @param v value to store
     */
    template <uint8_t vec_idx, typename T = val_t>
    inline void set_unsafe(pos_t i, T v)
    {
        static_assert(vec_idx < num_vectors);
        *reinterpret_cast<T*>(bases[vec_idx] + i * width_entry) = v;
    }

    /**
     * @brief returns the i-th entry in the vector with index vec_idx
     * @tparam vec_idx vector index (0 <= vec_idx < num_vectors)
     * @param i entry index (0 <= i < size_vectors)
     * @return value
     */
    template <uint8_t vec_idx, typename T = val_t>
    inline T get(pos_t i) const
    {
        static_assert(vec_idx < num_vectors);
        return *reinterpret_cast<T*>(bases[vec_idx] + i * width_entry) & masks_get[vec_idx];
    }

    /**
     * @brief interprets the memory location of the i-th entry in the vector with index vec_idx as an object of type
     * T and returns it (if sizeof(T) > widths[vec_idx], the returned value contains data from subsequent entries)
     * @tparam vec_idx vector index (0 <= vec_idx < num_vectors)
     * @tparam pos_t type to treat the i-th entry in the vector with index vec_idx as
     * @param i entry index (0 <= i < size_vectors)
     * @return value
     */
    template <uint8_t vec_idx, typename T = val_t>
    inline T get_unsafe(pos_t i) const
    {
        static_assert(vec_idx < num_vectors);
        return *reinterpret_cast<T*>(bases[vec_idx] + i * width_entry);
    }

    /**
     * @brief appends val to the end of the interleaved vectors
     * @tparam vec_idx vector index (0 <= vec_idx < num_vectors)
     * @tparam T type of value to append
     * @param val value to append
     */
    template <uint8_t vec_idx = 0, typename T = val_t>
    inline void emplace_back(T val)
    {
        if (size_vectors == capacity_vectors) {
            reserve(1.5 * size_vectors);
        }

        set<vec_idx, T>(size_vectors, val);
        size_vectors++;
    }

    /**
     * @brief unsafely appends val to the end of the interleaved vectors
     * (if sizeof(T) > widths[vec_idx], this is unsafe, because
     * subsequent entries are overwritten)
     * @tparam vec_idx vector index (0 <= vec_idx < num_vectors)
     * @tparam T type of value to append
     * @param val value to append
     */
    template <uint8_t vec_idx = 0, typename T = val_t>
    inline void emplace_back_unsafe(T val)
    {
        if (size_vectors == capacity_vectors) {
            reserve(1.5 * size_vectors);
        }

        set_unsafe<vec_idx, T>(size_vectors, val);
        size_vectors++;
    }

    /**
     * @brief appends a tuple of values to the end of the interleaved vectors
     * @param vals tuple of values
     */
    template <typename... Ts>
    inline void emplace_back(std::tuple<Ts...> vals)
    {
        static_assert(std::tuple_size<std::tuple<Ts...>>::value <= num_vectors);

        if (size_vectors == capacity_vectors) {
            reserve(1.5 * size_vectors);
        }

        for_constexpr<0, sizeof...(Ts), 1>([&](auto vec_idx) {
            set<vec_idx, std::tuple_element_t<vec_idx, std::tuple<Ts...>>>(size_vectors, std::get<vec_idx>(vals));
        });

        size_vectors++;
    }

    /**
     * @brief unsafely appends a tuple of values to the end of the interleaved vectors (if
     * there exists a vec_idx \in [0,sizeof...(Ts...)) s.t. sizeof(Ts...[vec_idx]) > widths[vec_idx],
     * this is unsafe, because subsequent entries are overwritten)
     * @param vals tuple of values
     */
    template <typename... Ts>
    inline void emplace_back_unsafe(std::tuple<Ts...> vals)
    {
        static_assert(std::tuple_size<std::tuple<Ts...>>::value <= num_vectors);

        if (size_vectors == capacity_vectors) {
            reserve(1.5 * size_vectors);
        }

        for_constexpr<0, sizeof...(Ts), 1>([&](auto vec_idx) {
            set_unsafe<vec_idx, std::tuple_element_t<vec_idx, std::tuple<Ts...>>>(size_vectors, std::get<vec_idx>(vals));
        });

        size_vectors++;
    }

    /**
     * @brief reinterpret the memory at data as interleaved vectors of size size; do not perform any operations that
     * may change the size or the capacity of the interleaved vectors after using this method
     * @param data memory location storing interleaved vectors with the same byte-widths as stored in widths
     * @param size size of the interleaved vectors stored at data
     */
    void set_data(char* data, uint64_t size)
    {
        data_vectors.clear();
        data_vectors.shrink_to_fit();

        size_vectors = size;
        capacity_vectors = size;

        set_bases(data);
    }

    /**
     * @brief serializes the interleaved vectors to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) const
    {
        out.write((char*) &size_vectors, sizeof(uint64_t));
        out.write((char*) &width_entry, sizeof(uint64_t));

        if (num_vectors > 0) {
            out.write((char*) widths.data(), num_vectors * sizeof(uint64_t));
            out.write((char*) masks_get.data(), num_vectors * sizeof(val_t));
            out.write((char*) masks_set.data(), num_vectors * sizeof(val_t));
        }

        if (size_vectors > 0) {
            write_to_file(out, data(), size_vectors * width_entry);
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
        in.read((char*) &width_entry, sizeof(uint64_t));

        if (num_vectors > 0) {
            in.read((char*) widths.data(), num_vectors * sizeof(uint64_t));
            in.read((char*) masks_get.data(), num_vectors * sizeof(val_t));
            in.read((char*) masks_set.data(), num_vectors * sizeof(val_t));
        }

        if (old_size > 0) {
            resize_no_init(old_size);
            set_bases(&data_vectors[0]);
            read_from_file(in, data(), size_vectors * width_entry);
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