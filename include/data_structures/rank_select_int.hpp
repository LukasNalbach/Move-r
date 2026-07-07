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

#include <filesystem>
#include <iostream>
#include <omp.h>
#include <vector>

#include <data_structures/hybrid_bit_vector.hpp>
#include <data_structures/interleaved_bit_aligned_vectors.hpp>
#include <misc/utils.hpp>
#include <misc/search.hpp>

/**
 * @brief a rank-select data structure using hybrid bit vectors (either sd_array or plain bit vector), or the rs data structure
 * @tparam sym_t symbol type
 * @tparam pos_t position type
 * @tparam build_rank_support whether to build rank support
 * @tparam build_select_support whether to build select support
 */
template <typename sym_t, typename pos_t = uint32_t, bool build_rank_support = true, bool build_select_support = true>
class rank_select_int {
protected:
    static_assert(std::is_same_v<pos_t, uint32_t> || std::is_same_v<pos_t, uint64_t>);
    static_assert(build_rank_support || build_select_support);

    // rank_select_int is for integer alphabets only (byte alphabets use rank_select_byte / pred_succ_byte)
    static_assert(
        std::is_same_v<sym_t, uint16_t> ||
        std::is_same_v<sym_t, uint32_t> ||
        std::is_same_v<sym_t, uint64_t>
    );

    using i_sym_t = sym_t; // internal (unsigned) symbol type
    using inp_t = std::vector<sym_t>; // input container type
    using hybrid_bv_t = hybrid_bit_vector<pos_t, build_rank_support, false, build_select_support>; // hybrid bit vector type

    /* maximum number of occurrences to use scanning instead of binary search over the occurrences for answering rank */
    static constexpr pos_t max_occ_scan_rank = 16;
    /* minimum number of occurrences to use a bit vector for answering rank */
    static constexpr pos_t min_occ_vec_rank = 512;
    /* maximum relative number of occurrences to use an sd_array instead of a plain bitvector for answering rank (and select) */
    static constexpr double thrsh_sd_array = hybrid_bit_vector<pos_t>::compression_threshold;
    /* maximum relative number of occurrences to store all occurrences plainly to answer select */
    static constexpr double thrsh_plain_select = 0.01;

    /* minimum number of occurrences to build a bit vector that marks all occurrences */
    static constexpr pos_t min_occ_vec = min_occ_vec_rank;

    pos_t input_size = 0; // the size of the input
    pos_t sigma = 0; // the number of distinct symbols in the input
    pos_t num_vectors = 0; // the number of initialized vectors in hyb_bit_vecs

    // ############################# DATA STRUCTURES #############################

    // hyb_bit_vecs[vec_idx[v]] marks (with ones) the occurrences of v in the input, for each frequent value v
    std::vector<hybrid_bv_t> hyb_bit_vecs;
    // [0..sigma-1] vec_idx[v] = the position in hyb_bit_vecs of v's bit vector if freq(v) > min_occ_vec, else sigma
    interleaved_bit_aligned_vectors<pos_t> vec_idx;
    // [0..sigma-1] c_arr[v] = the sum of the frequencies of all v' < v with freq(v') <= max_occ_plain
    interleaved_bit_aligned_vectors<pos_t> c_arr;
    // occs[c_arr[v] + j] = the position of the (j+1)-th occurrence of v in the input
    interleaved_bit_aligned_vectors<pos_t> occs;

public:
    rank_select_int() = default;

    // builds the data structure over input
    rank_select_int(const inp_t& input, pos_t alphabet_size) { build([&](pos_t i){ return input[i]; }, input.size(), alphabet_size); }

    // builds the data structure over input, reading it with the function read
    template <typename read_fnc_t>
    rank_select_int(read_fnc_t read, pos_t input_size, pos_t alphabet_size) { build(read, input_size, alphabet_size); }

protected:
    // builds the data structure over input, reading it with the function read
    template <typename read_fnc_t>
    void build(read_fnc_t read, pos_t input_size, pos_t alphabet_size)
    {
        this->input_size = input_size;
        this->sigma = alphabet_size;
        pos_t max_occ_plain = build_select_support ? input_size * thrsh_plain_select : min_occ_vec;
        pos_t max_occ_sd_array = input_size * thrsh_sd_array;
        pos_t min_occ_vec_tmp = std::min<pos_t>(min_occ_vec, max_occ_plain);
        uint8_t bits_per_entry = bit_width(input_size);

        // first pass: c_arr[v] = frequency of v (histogram); used as the frequency array until it is converted
        // in place to the prefix sum of plain frequencies at the end
        c_arr = interleaved_bit_aligned_vectors<pos_t>({ bits_per_entry });
        c_arr.resize(sigma + 1); // zero-initialized
        for (pos_t i = 0; i < input_size; i++)
            c_arr.template set<0, pos_t>(read(i), c_arr[read(i)] + 1);

        interleaved_bit_aligned_vectors<pos_t> occ_idx({ bits_per_entry }); // occ_idx[v] = write offset into occs
        std::vector<sdsl::bit_vector> plain_bvs;
        std::vector<sdsl::sd_vector_builder> sdv_builders;
        num_vectors = 0;

        vec_idx = interleaved_bit_aligned_vectors<pos_t>({ bit_width(sigma) });
        vec_idx.resize_no_init(sigma);
        occ_idx.resize_no_init(sigma);

        // decide per symbol whether its occurrences are stored plainly (in occs) and/or in a bit vector, and
        // accumulate the number of plainly-stored occurrences in num_occs (= the size of occs)
        pos_t num_occs = 0;

        for (pos_t v = 0; v < sigma; v++) {
            pos_t freq_v = c_arr[v]; // c_arr still holds the histogram
            occ_idx.template set<0, pos_t>(v, num_occs);

            if (freq_v > min_occ_vec_tmp) {
                if (freq_v <= max_occ_sd_array) {
                    sdv_builders.emplace_back(sdsl::sd_vector_builder(input_size, freq_v));
                    plain_bvs.emplace_back(sdsl::bit_vector());
                } else {
                    sdv_builders.emplace_back(sdsl::sd_vector_builder());
                    plain_bvs.emplace_back(sdsl::bit_vector(input_size));
                }

                vec_idx.template set<0, pos_t>(v, num_vectors);
                num_vectors++;
            } else {
                vec_idx.template set<0, pos_t>(v, sigma);
            }

            if (freq_v <= max_occ_plain) num_occs += freq_v;
        }

        occs = interleaved_bit_aligned_vectors<pos_t>({ bits_per_entry });
        occs.resize_no_init(num_occs);

        // second pass: fill occs (the plainly-stored occurrences) and the bit vectors
        for (pos_t i = 0; i < input_size; i++) {
            pos_t v = read(i);
            pos_t freq_v = c_arr[v];

            if (freq_v <= max_occ_plain) {
                pos_t oi = occ_idx[v];
                occs.template set<0, pos_t>(oi, i);
                occ_idx.template set<0, pos_t>(v, oi + 1);
            }

            if (freq_v > min_occ_vec_tmp) {
                if (freq_v <= max_occ_sd_array) {
                    sdv_builders[vec_idx[v]].set(i);
                } else {
                    plain_bvs[vec_idx[v]][i] = 1;
                }
            }
        }

        occ_idx.clear();
        occ_idx.shrink_to_fit();
        pos_t prefix = 0;

        for (pos_t v = 0; v < sigma; v++) {
            pos_t freq_v = c_arr[v];
            c_arr.template set<0, pos_t>(v, prefix);
            if (freq_v <= max_occ_plain) prefix += freq_v;
        }

        c_arr.template set<0, pos_t>(sigma, prefix);
        hyb_bit_vecs.resize(num_vectors);

        for (pos_t i = 0; i < num_vectors; i++) {
            if (plain_bvs[i].size() != 0) {
                hyb_bit_vecs[i] = hybrid_bv_t(std::move(plain_bvs[i]));
            } else if (sdv_builders[i].size() != 0) {
                hyb_bit_vecs[i] = hybrid_bv_t(sdsl::sd_vector<>(sdv_builders[i]));
            }
        }
    }

public:
    // returns the size of the input
    inline pos_t size() const { return input_size; }

    // returns the number of distinct values in the input
    inline pos_t alphabet_size() const { return sigma; }

    // returns whether the input vector is empty
    inline bool empty() const { return size() == 0; }

    /**
     * @brief returns the size of the data structure in bytes
     * @return size of the data structure in bytes
     */
    uint64_t size_in_bytes() const
    {
        uint64_t size = 0;
        size += vec_idx.size_in_bytes();
        size += c_arr.size_in_bytes();
        size += occs.size_in_bytes();
        for (pos_t i = 0; i < num_vectors; i++)
            size += hyb_bit_vecs[i].size_in_bytes();
        return size;
    }

    /**
     * @brief returns whether v occurs in the input
     * @param v value
     * @return whether v occurs in the input
     */
    inline bool contains(sym_t v) const { return v < sigma; }

    /**
     * @brief returns the number of occurrences of the symbol v in the input (v must occur in the input)
     * @param v symbol that occurs in the input
     * @return the number of occurrences of v in the input
     */
    inline pos_t frequency(sym_t v) const
    {
        pos_t diff = c_arr[v + 1] - c_arr[v];
        return diff != 0 ? diff : hyb_bit_vecs[vec_idx[v]].num_ones();
    }


    /**
     * @brief returns the number of occurrences of v before index i
     * @param v [0..alphabet_size-1] value that occurs in the input
     * @param i [1..input size] position in the input
     * @return number of occurrences of v before index i
     */
    inline pos_t rank(sym_t v, pos_t i) const
    {
        static_assert(build_rank_support);

        pos_t v_s = c_arr[v];
        pos_t v_e = c_arr[v + 1];

        if (v_e == v_s || v_e - v_s > min_occ_vec_rank) {
            return hyb_bit_vecs[vec_idx[v]].rank_1(i);
        } else {
            pos_t pos;

            if (v_e - v_s <= max_occ_scan_rank) {
                if (occs[v_s] >= i) return 0;
                pos = v_s;
                v_e--;
                while (pos < v_e && occs[pos + 1] < i) pos++;
                return pos - v_s + 1;
            } else {
                pos = bin_search_max_lt<pos_t>(i, v_s, v_e - 1, [&](pos_t x) { return occs[x]; });
                return occs[pos] >= i ? 0 : pos - v_s + 1;
            }
        }
    }


    /**
     * @brief returns the index of the i-th occurrence of v
     * @param v a value that occurs in the input
     * @param i [1..number of occurrences of v in the input]
     * @return index of the i-th occurrence of v
     */
    inline pos_t select(sym_t v, pos_t i) const
    {
        static_assert(build_select_support);

        pos_t v_s = c_arr[v];
        pos_t v_e = c_arr[v + 1];
        return v_e != v_s ? occs[v_s + i - 1] : hyb_bit_vecs[vec_idx[v]].select_1(i);
    }

    // character successor: smallest index in [x, max] with symbol sym (a value > max if none)
    inline pos_t succ(sym_t sym, pos_t x, pos_t max) const
    {
        pos_t rnk = rank(sym, x);
        if (rnk == frequency(sym)) return max + 1;
        return select(sym, rnk + 1);
    }

    // character predecessor: largest index in [min, x] with symbol sym (sym must occur in [min, x])
    inline pos_t pred(sym_t sym, pos_t x, pos_t) const { return select(sym, rank(sym, x)); }

    /**
     * @brief serializes the rank_select_int to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) const
    {
        out.write((char*) &input_size, sizeof(pos_t));

        if (input_size != 0) {
            out.write((char*) &sigma, sizeof(pos_t));
            out.write((char*) &num_vectors, sizeof(pos_t));
            vec_idx.serialize(out);
            c_arr.serialize(out);
            occs.serialize(out);

            for (pos_t i = 0; i < num_vectors; i++) {
                hyb_bit_vecs[i].serialize(out);
            }
        }
    }

    /**
     * @brief loads the rank_select_int from an input stream
     * @param in input stream
     */
    void load(std::istream& in)
    {
        in.read((char*) &input_size, sizeof(pos_t));

        if (input_size != 0) {
            in.read((char*) &sigma, sizeof(pos_t));
            in.read((char*) &num_vectors, sizeof(pos_t));
            vec_idx.load(in);
            c_arr.load(in);
            occs.load(in);
            hyb_bit_vecs.resize(num_vectors);

            for (pos_t i = 0; i < num_vectors; i++) {
                hyb_bit_vecs[i].load(in);
            }
        }
    }

    /**
     * @brief serializes the rank_select_int to an output stream
     * @param os output stream
     * @return the output stream
     */
    std::ostream& operator>>(std::ostream& os) const
    {
        serialize(os);
        return os;
    }

    /**
     * @brief loads the rank_select_int from an input stream
     * @param is input stream
     * @return the input stream
     */
    std::istream& operator<<(std::istream& is)
    {
        load(is);
        return is;
    }
};