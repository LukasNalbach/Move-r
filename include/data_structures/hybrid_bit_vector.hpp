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

#include "plain_bit_vector.hpp"
#include "sd_array.hpp"
#include <iostream>
#include <optional>
#include <vector>

/**
 * @brief hybrid bit vector (either an sd_array or a plain bit vector)
 * @tparam pos_t unsigned integer type
 * @tparam build_rank_support whether to build rank support
 * @tparam build_select_0_support whether to build select-0 support
 * @tparam build_select_1_support whether to build select-1 support
 */
template <typename pos_t = uint32_t, bool build_rank_support = false, bool build_select_0_support = false, bool build_select_1_support = false>
class hybrid_bit_vector {
    static_assert(std::is_same_v<pos_t, uint32_t> || std::is_same_v<pos_t, uint64_t>);

public:
    using plain_bv_t = plain_bit_vector<pos_t, build_rank_support, build_select_0_support, build_select_1_support>;

    /* if there are less than (compression_threshold * size of the bitvector)
     * ones in the bit vector, then it should be compressed */
    static constexpr double compression_threshold = 0.1;

protected:
    std::optional<sd_array<pos_t>> sd_arr; // the sd_array
    std::optional<plain_bv_t> plain_bit_vec; // the plain bit vector

    /**
     * @brief returns, whether plain_bit_vec should be compressed
     * @param plain_bit_vec a bit vector
     * @return whether plain_bit_vec should be compressed
     */
    bool should_be_compressed(const sdsl::bit_vector& plain_bit_vec)
    {
        return sdsl::bit_vector::rank_1_type(&plain_bit_vec).rank(plain_bit_vec.size()) <
            plain_bit_vec.size() * compression_threshold;
    }

public:
    hybrid_bit_vector() = default;

    /**
     * @brief constructs a new hybrid_bit_vector from a bit vector
     * @param plain_bit_vec a bit vector
     */
    hybrid_bit_vector(const sdsl::bit_vector& plain_bit_vec)
    {
        if (should_be_compressed(plain_bit_vec))
             this->sd_arr = sd_array<pos_t>(plain_bit_vec);
        else this->plain_bit_vec = plain_bv_t(plain_bit_vec);
    }

    /**
     * @brief constructs a new hybrid_bit_vector from a bit vector
     * @param plain_bit_vec a bit vector
     */
    hybrid_bit_vector(sdsl::bit_vector&& plain_bit_vec)
    {
        if (should_be_compressed(plain_bit_vec))
             this->sd_arr = sd_array<pos_t>(plain_bit_vec);
        else this->plain_bit_vec = plain_bv_t(std::move(plain_bit_vec));
    }

    /**
     * @brief constructs a new hybrid_bit_vector from an sd_vector
     * @param sd_vec an sd_vector
     */
    hybrid_bit_vector(const sdsl::sd_vector<>& sd_vec) { sd_arr = sd_array<pos_t>(sd_vec); }

    /**
     * @brief constructs a new hybrid_bit_vector from an sd_vector
     * @param sd_vec an sd_vector
     */
    hybrid_bit_vector(sdsl::sd_vector<>&& sd_vec) { sd_arr = sd_array<pos_t>(std::move(sd_vec)); }

    // whether the bit vector has been initialized
    inline bool is_initialized() const { return sd_arr.has_value() || plain_bit_vec.has_value(); }

    // whether the bit vector is compressed
    inline bool is_compressed() const { return sd_arr.has_value(); }

    // the size of the bit vector
    inline pos_t size() const { return is_compressed() ? sd_arr.value().size() : plain_bit_vec.value().size(); }

    // whether the bit vector is empty
    inline bool empty() const { return size() == 0; }

    // the number of ones in the bit vector
    inline pos_t num_ones() const { return is_compressed() ? sd_arr.value().num_ones() : plain_bit_vec.value().num_ones(); }

    // the number of zeros in the bit vector
    inline pos_t num_zeros() const { return is_compressed() ? sd_arr.value().num_zeros() : plain_bit_vec.value().num_zeros(); }

    // the number of ones before index i in [0..size]
    inline pos_t rank_1(pos_t i) const
    {
        static_assert(build_rank_support);
        return is_compressed() ? sd_arr.value().rank_1(i) : plain_bit_vec.value().rank_1(i);
    }

    // the index of the i-th one, i in [1..number of ones]
    inline pos_t select_1(pos_t i) const
    {
        static_assert(build_select_1_support);
        return is_compressed() ? sd_arr.value().select_1(i) : plain_bit_vec.value().select_1(i);
    }

    // the number of zeros before index i in [0..size]
    inline pos_t rank_0(pos_t i) const
    {
        static_assert(build_rank_support);
        return is_compressed() ? sd_arr.value().rank_0(i) : plain_bit_vec.value().rank_0(i);
    }

    // the index of the i-th zero, i in [1..number of zeros]
    inline pos_t select_0(pos_t i) const
    {
        static_assert(build_select_0_support);
        return is_compressed() ? sd_arr.value().select_0(i) : plain_bit_vec.value().select_0(i);
    }

    // the index of the next one after index i in [1..size-1]
    inline pos_t next_1(pos_t i) const { return is_compressed() ? sd_arr.value().next_1(i) : plain_bit_vec.value().next_1(i); }

    // the index of the previous one before index i in [1..size-1]
    inline pos_t previous_1(pos_t i) const { return is_compressed() ? sd_arr.value().previous_1(i) : plain_bit_vec.value().previous_1(i); }

    // the index of the next zero after index i in [1..size-1]
    inline pos_t next_0(pos_t i) const { return is_compressed() ? sd_arr.value().next_0(i) : plain_bit_vec.value().next_0(i); }

    // the index of the previous zero before index i in [1..size-1]
    inline pos_t previous_0(pos_t i) const { return is_compressed() ? sd_arr.value().previous_0(i) : plain_bit_vec.value().previous_0(i); }

    // whether there is a one at index i in [0..size-1]
    inline bool operator[](pos_t i) const { return is_compressed() ? sd_arr.value()[i] : plain_bit_vec.value()[i]; }

    // the size of the data structure in bytes
    inline uint64_t size_in_bytes() const
    {
        if (!is_initialized()) return 0;
        return is_compressed() ? sd_arr.value().size_in_bytes() : plain_bit_vec.value().size_in_bytes();
    }

    // logs the contents of the bit vector
    void log_contents() const
    {
        if (empty()) return;
        for (uint64_t i = 0; i < size() - 1; i++)
            std::cout << operator[](i) << " ";
        std::cout << operator[](size() - 1) << std::endl;
    }

    /**
     * @brief serializes the hybrid_bit_vector to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) const
    {
        bool is_init = is_initialized();
        out.write((char*) &is_init, 1);
        if (!is_init) return;
        bool compressed = is_compressed();
        out.write((char*) &compressed, 1);
        if (compressed) sd_arr.value().serialize(out);
        else plain_bit_vec.value().serialize(out);
    }

    /**
     * @brief loads the hybrid_bit_vector from an input stream
     * @param in input stream
     */
    void load(std::istream& in)
    {
        bool is_init;
        in.read((char*) &is_init, 1);
        if (!is_init) return;
        bool compressed;
        in.read((char*) &compressed, 1);

        if (compressed) {
            sd_arr = sd_array<pos_t>();
            sd_arr.value().load(in);
        } else {
            plain_bit_vec = plain_bv_t();
            plain_bit_vec.value().load(in);
        }
    }

    /**
     * @brief serializes the hybrid_bit_vector to an output stream
     * @param os output stream
     * @return the output stream
     */
    std::ostream& operator>>(std::ostream& os) const
    {
        serialize(os);
        return os;
    }

    /**
     * @brief loads the hybrid_bit_vector from an input stream
     * @param is input stream
     * @return the input stream
     */
    std::istream& operator<<(std::istream& is)
    {
        load(is);
        return is;
    }
};