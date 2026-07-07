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

#include <functional>
#include <limits>
#include <optional>

#include <data_structures/interleaved_bit_aligned_vectors.hpp>
#include <misc/utils.hpp>

/**
 * @brief character predecessor/successor queries over a byte-alphabet input read through a stored function read(i),
 *        using per-block lookpup tables plus a scan of the block containing the query
 * @tparam sym_t symbol type (must be a byte type)
 * @tparam pos_t unsigned position type
 * @tparam read_fnc_t type of the read function (read(i) returns the symbol at index i)
 */
template <typename sym_t, typename pos_t = uint32_t, typename read_fnc_t = std::function<sym_t(pos_t)>>
class pred_succ_byte {
protected:
    static_assert(sizeof(sym_t) == 1); // byte alphabet only
    static_assert(std::is_same_v<pos_t, uint32_t> || std::is_same_v<pos_t, uint64_t>);

    static constexpr pos_t infty = std::numeric_limits<pos_t>::max();
    static constexpr pos_t block_size_factor = 4; // a block spans block_size_factor * sigma input positions
    mutable std::optional<read_fnc_t> read; // the stored read function (rebindable with set_read)
    pos_t input_size = 0; // number of input positions
    pos_t sigma = 0; // alphabet size
    pos_t blk_size = 0; // block size in input positions
    pos_t num_blks = 0; // number of blocks

    // per-block lookpup tables. _succ[blk*sigma + c] = first index >= blk*blk_size whose symbol is c (input_size if
    // none); _pred[blk*sigma + c] = last index < blk*blk_size whose symbol is c (input_size if none)
    interleaved_bit_aligned_vectors<pos_t> _succ;
    interleaved_bit_aligned_vectors<pos_t> _pred;

public:
    pred_succ_byte() = default;

    /**
     * @brief builds the structure over input[0..input_size-1] read through the function read
     * @param read function returning the symbol at index i
     * @param input_size number of input positions
     * @param alphabet_size alphabet size
     * @param num_threads number of threads to use during the construction
     */
    pred_succ_byte(read_fnc_t read, pos_t input_size, pos_t alphabet_size, uint16_t num_threads = 1)
    {
        this->read.emplace(read);
        this->input_size = input_size;
        this->sigma = alphabet_size;
        const read_fnc_t& rd = *this->read;
        uint8_t bits = bit_width(input_size); // values in [0, input_size], input_size = no occurrence
        blk_size = std::min<pos_t>(block_size_factor * sigma, input_size);
        num_blks = div_ceil<pos_t>(input_size, blk_size);

        #pragma omp parallel sections num_threads(num_threads)
        {
            #pragma omp section
            {
                _succ = interleaved_bit_aligned_vectors<pos_t>({ bits });
                _succ.resize_no_init((num_blks + 1) * sigma);
                std::vector<pos_t> next(sigma, input_size);

                for (int64_t blk = num_blks; blk >= 0; blk--) {
                    int64_t beg = blk * (int64_t) blk_size;
                    int64_t end = std::min<int64_t>(beg + blk_size, input_size);
                    pos_t beg_next = blk * sigma;

                    for (int64_t i = end - 1; i >= beg; i--) next[rd(i)] = i;
                    for (pos_t c = 0; c < sigma; c++) _succ.template set<0, pos_t>(beg_next + c, next[c]);
                }
            }

            #pragma omp section
            {
                _pred = interleaved_bit_aligned_vectors<pos_t>({ bits });
                _pred.resize_no_init(num_blks * sigma);
                std::vector<pos_t> prev(sigma, input_size);

                for (int64_t blk = 0; blk < (int64_t) num_blks; blk++) {
                    int64_t beg = blk * (int64_t) blk_size;
                    int64_t end = std::min<int64_t>(beg + blk_size, input_size);
                    pos_t beg_prev = blk * sigma;

                    for (pos_t c = 0; c < sigma; c++) _pred.template set<0, pos_t>(beg_prev + c, prev[c]);
                    for (int64_t i = beg; i < end; i++) prev[rd(i)] = i;
                }
            }
        }
    }

    // (re-)binds the read function; call after loading, or after the object read references has moved
    inline void set_read(read_fnc_t read) const { this->read.emplace(read); }

    inline pos_t block_size() const { return blk_size; }
    inline pos_t num_blocks() const { return num_blks; }

    // block-granularity jump tables: first/last index of symbol sym at/before block blk
    inline pos_t succ_blk(pos_t blk, pos_t sym) const { return _succ.template get<0, pos_t>(blk * sigma + sym); }
    inline pos_t pred_blk(pos_t blk, pos_t sym) const { return _pred.template get<0, pos_t>(blk * sigma + sym); }

    /**
     * @brief returns the smallest index in [x, max] whose symbol is sym (a value > max if none)
     * @param sym a symbol
     * @param x lower bound
     * @param max upper bound
     */
    inline pos_t succ(sym_t sym, pos_t x, pos_t max = infty) const
    {
        const read_fnc_t& rd = *read;
        pos_t blk = div_ceil<pos_t>(x, blk_size);
        pos_t max_x = std::min<pos_t>({blk * blk_size, max, input_size});
        while (x <= max_x && sym != rd(x)) x++;

        if (x > max_x && sym != rd(x)) [[likely]] {
            if (x > max) [[unlikely]] return x;
            x = succ_blk(blk, sym);
        }

        return x;
    }

    /**
     * @brief returns the largest index in [min, x] whose symbol is sym (sym must occur in [min, x])
     * @param sym a symbol
     * @param x upper bound
     * @param min lower bound
     */
    inline pos_t pred(sym_t sym, pos_t x, pos_t min = 0) const
    {
        const read_fnc_t& rd = *read;
        pos_t blk = x / blk_size;
        pos_t min_x = std::max<pos_t>(blk * blk_size, min);
        while (x >= min_x && sym != rd(x)) x--;
        if (x < min_x && sym != rd(x)) [[likely]] x = pred_blk(blk, sym);
        return x;
    }

    uint64_t size_in_bytes() const { return sizeof(*this) + _succ.size_in_bytes() + _pred.size_in_bytes(); }

    void serialize(std::ostream& out) const
    {
        out.write((char*) &input_size, sizeof(pos_t));
        out.write((char*) &sigma, sizeof(pos_t));
        out.write((char*) &blk_size, sizeof(pos_t));
        out.write((char*) &num_blks, sizeof(pos_t));
        _succ.serialize(out);
        _pred.serialize(out);
    }

    // note: the read function is not serialized; call set_read after loading
    void load(std::istream& in)
    {
        in.read((char*) &input_size, sizeof(pos_t));
        in.read((char*) &sigma, sizeof(pos_t));
        in.read((char*) &blk_size, sizeof(pos_t));
        in.read((char*) &num_blks, sizeof(pos_t));
        _succ.load(in);
        _pred.load(in);
    }
};
