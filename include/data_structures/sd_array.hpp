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

#include <iostream>
#include <sdsl/bit_vectors.hpp>
#include <vector>

/**
 * @brief wrapper class for the sd_vector from sdsl
 * @tparam pos_t unsigned integer type
 */
template <typename pos_t = uint32_t>
class sd_array {
    static_assert(std::is_same_v<pos_t, uint32_t> || std::is_same_v<pos_t, uint64_t>);

protected:
    sdsl::sd_vector<> sd_vector;
    pos_t ones = 0;
    pos_t zeros = 0;

    // recomputes the number of ones/zeros (the sd_vector does not store them)
    void count() { if (size() > 0) { ones = rank_1(size()); zeros = size() - ones; } }

public:
    sd_array() = default;

    // constructs a new sd_array from a bit vector
    sd_array(const sdsl::bit_vector& bit_vector) { sd_vector = sdsl::sd_vector<>(bit_vector); count(); }

    // constructs a new sd_array from an sd_vector
    sd_array(const sdsl::sd_vector<>& sd_vector) { this->sd_vector = sd_vector; count(); }

    // constructs a new sd_array, moving the sd_vector into it
    sd_array(sdsl::sd_vector<>&& sd_vector) { this->sd_vector = std::move(sd_vector); count(); }

    inline pos_t size() const { return sd_vector.size(); }
    inline bool empty() const { return size() == 0; }
    inline pos_t num_ones() const { return ones; }
    inline pos_t num_zeros() const { return zeros; }
    inline pos_t rank_1(pos_t i) const { return sdsl::sd_vector<>::rank_1_type(&sd_vector).rank(i); } // ones before index i
    inline pos_t select_1(pos_t i) const { return sdsl::sd_vector<>::select_1_type(&sd_vector).select(i); } // index of the i-th one
    inline pos_t rank_0(pos_t i) const { return sdsl::sd_vector<>::rank_0_type(&sd_vector).rank(i); } // zeros before index i
    inline pos_t select_0(pos_t i) const { return sdsl::sd_vector<>::select_0_type(&sd_vector).select(i); } // index of the i-th zero
    inline pos_t next_1(pos_t i) const { return select_1(rank_1(i + 1) + 1); } // index of the next one after i
    inline pos_t previous_1(pos_t i) const { return select_1(rank_1(i)); } // index of the previous one before i
    inline pos_t next_0(pos_t i) const { return select_0(rank_0(i + 1) + 1); } // index of the next zero after i
    inline pos_t previous_0(pos_t i) const { return select_0(rank_0(i)); } // index of the previous zero before i
    inline bool operator[](pos_t i) const { return sd_vector[i]; } // whether there is a one at index i
    uint64_t size_in_bytes() const { return sdsl::size_in_bytes(sd_vector); }

    void serialize(std::ostream& out) const { sd_vector.serialize(out); }
    void load(std::istream& in) { sd_vector.load(in); count(); }
    std::ostream& operator>>(std::ostream& os) const { serialize(os); return os; }
    std::istream& operator<<(std::istream& is) { load(is); return is; }

    // logs the contents of the bit vector
    void log_contents() const
    {
        if (empty()) return;
        for (uint64_t i = 0; i < size() - 1; i++) std::cout << operator[](i) << " ";
        std::cout << operator[](size() - 1) << std::endl;
    }
};