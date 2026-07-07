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
#include <string>
#include <type_traits>
#include <vector>
#include <omp.h>

#include <data_structures/hybrid_bit_vector.hpp>
#include <misc/utils.hpp>

/**
 * @brief rank-select data structure for a byte alphabet: one hybrid bit vector (sd_array or plain bit vector,
 *        chosen by density) per used character, marking that character's occurrences in the input
 * @tparam sym_t symbol type (must be a byte type)
 * @tparam pos_t position type
 * @tparam build_rank_support whether to build rank support
 * @tparam build_select_support whether to build select support
 */
template <typename sym_t, typename pos_t = uint32_t, bool build_rank_support = true, bool build_select_support = true>
class rank_select_byte {
protected:
    static_assert(std::is_same_v<pos_t, uint32_t> || std::is_same_v<pos_t, uint64_t>);
    static_assert(build_rank_support || build_select_support);
    static_assert(sizeof(sym_t) == 1); // byte alphabet only

    static constexpr bool str_input = std::is_same_v<sym_t, char>; // true <=> the input is a string
    using inp_t = std::conditional_t<str_input, std::string, std::vector<sym_t>>; // input container type
    using hybrid_bv_t = hybrid_bit_vector<pos_t, build_rank_support, false, build_select_support>;
    // maximum relative number of occurrences to use an sd_array instead of a plain bit vector
    static constexpr double thrsh_sd_array = hybrid_bit_vector<pos_t>::compression_threshold;

    pos_t input_size = 0; // the size of the input
    pos_t sigma = 0; // the number of distinct symbols in the input

    std::vector<pos_t> freq; // [0..255] frequency of each symbol
    // [0..255] hyb_bit_vecs[c] marks (with ones) the occurrences of the character c in the input
    std::vector<hybrid_bv_t> hyb_bit_vecs;

    // maps a symbol to its (unsigned) internal index in [0..255]
    uint8_t symbol_idx(sym_t sym) const { return static_cast<uint8_t>(sym); }

public:
    rank_select_byte() = default;

    //builds the data structure
    rank_select_byte(const inp_t& input) { build([&](pos_t i){ return input[i]; }, input.size()); }

    template <typename read_fnc_t>
    rank_select_byte(read_fnc_t read, pos_t input_size) { build(read, input_size); }

protected:
    // builds the data structure by reading the input using the function read (input_size must be set beforehand)
    template <typename read_fnc_t>
    void build(read_fnc_t read, pos_t input_size)
    {
        this->input_size = input_size;
        freq.resize(256, 0);
        for (pos_t i = 0; i < input_size; i++) freq[symbol_idx(read(i))]++;
        for (pos_t v = 0; v < 256; v++) if (freq[v] != 0) sigma++;

        // for each character, mark its occurrences in a plain bit vector or (if sparse enough) an sd_array
        pos_t max_occ_sd_array = input_size * thrsh_sd_array;
        std::vector<sdsl::sd_vector_builder> sdv_builders(256);
        std::vector<sdsl::bit_vector> plain_bvs(256);

        for (pos_t v = 0; v < 256; v++) {
            if (freq[v] != 0) {
                if (freq[v] <= max_occ_sd_array) sdv_builders[v] = sdsl::sd_vector_builder(input_size, freq[v]);
                else plain_bvs[v] = sdsl::bit_vector(input_size);
            }
        }

        for (pos_t i = 0; i < input_size; i++) {
            pos_t v = symbol_idx(read(i));
            if (freq[v] <= max_occ_sd_array) sdv_builders[v].set(i);
            else plain_bvs[v][i] = 1;
        }

        hyb_bit_vecs.resize(256);
        for (pos_t v = 0; v < 256; v++) {
            if (plain_bvs[v].size() != 0) hyb_bit_vecs[v] = hybrid_bv_t(std::move(plain_bvs[v]));
            else if (sdv_builders[v].size() != 0) hyb_bit_vecs[v] = hybrid_bv_t(sdsl::sd_vector<>(sdv_builders[v]));
        }
    }

public:
    inline pos_t size() const { return input_size; }
    inline pos_t alphabet_size() const { return sigma; }
    inline bool empty() const { return size() == 0; }
    inline bool contains(sym_t v) const { return frequency(v) != 0; } // whether v occurs in the input
    inline pos_t frequency(sym_t v) const { return freq[symbol_idx(v)]; } // number of occurrences of v

    // number of occurrences of v before index i
    inline pos_t rank(sym_t v, pos_t i) const
    {
        static_assert(build_rank_support);
        return hyb_bit_vecs[symbol_idx(v)].rank_1(i);
    }

    // index of the i-th occurrence of v
    inline pos_t select(sym_t v, pos_t i) const
    {
        static_assert(build_select_support);
        return hyb_bit_vecs[symbol_idx(v)].select_1(i);
    }

    uint64_t size_in_bytes() const
    {
        uint64_t size = sizeof(*this) + freq.size() * sizeof(pos_t);
        for (const hybrid_bv_t& v : hyb_bit_vecs) size += v.size_in_bytes();
        return size;
    }

    void serialize(std::ostream& out) const
    {
        out.write((char*) &input_size, sizeof(pos_t));
        if (input_size != 0) {
            out.write((char*) &sigma, sizeof(pos_t));
            out.write((char*) &freq[0], 256 * sizeof(pos_t));
            for (pos_t v = 0; v < 256; v++) hyb_bit_vecs[v].serialize(out);
        }
    }

    void load(std::istream& in)
    {
        in.read((char*) &input_size, sizeof(pos_t));
        if (input_size != 0) {
            in.read((char*) &sigma, sizeof(pos_t));
            no_init_resize(freq, 256);
            in.read((char*) &freq[0], 256 * sizeof(pos_t));
            hyb_bit_vecs.resize(256);
            for (pos_t v = 0; v < 256; v++) hyb_bit_vecs[v].load(in);
        }
    }
};
