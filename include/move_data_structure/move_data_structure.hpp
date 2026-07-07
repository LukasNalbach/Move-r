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

#include <cmath>
#include <tuple>
#include <utility>
#include <data_structures/interleaved_bit_aligned_vectors.hpp>
#include <data_structures/interleaved_byte_aligned_vectors.hpp>
#include <misc/utils.hpp>
#include <misc/search.hpp>
#include <omp.h>

/**
 * @brief encoding of the interval starting positions of a move data structure
 */
enum move_pos_encoding_t : uint8_t {
    // stores D_p, i.e. the absolute interval starting positions p_x (width omega_p = log n).
    // p(x) is a direct O(1) read; this is the classical move data structure.
    POS = 0,
    // stores D_len, i.e. the interval lengths len_x = p_{x+1} - p_x (width omega_offs = log l_max),
    // together with a sampled p-sample array S[j] = p_{j*d}. p(x) is reconstructed in O(d) time,
    // but the move query is performed on (run,offset)-pairs and never needs p on its hot path. This
    // trades the wide D_p slot (log n bits) for the narrow D_len slot (log l_max bits).
    DIFF = 1
};

struct mds_params {
    uint16_t num_threads = omp_get_max_threads(); // the number of threads to use during the construction
    uint16_t a = 8; // balancing parameter, restricts the number of intervals in the resulting move data structure to k*(a/(a-1))
    uint64_t d = 0; // p-sample sampling parameter (a p-sample S[j] = p_{j*d} is stored every d intervals); only used for the DIFF encoding
    bool log = false; // controls whether to print log messages during the construction
    std::ostream* mf = nullptr; // measurement file to write runtime data to
};

/**
 * @brief move data structure. Beyond the 3 core interleaved fields (D_p/D_len, D_idx, D_offs) it can interleave
 *        an arbitrary number of extra typed rows (e.g. L' for M_LF), exposed generically via row<J> / set_row<J>.
 * @tparam pos_t unsigned integer type of the interval starting positions
 * @tparam enc encoding of the interval starting positions (D_p vs. D_len + p-samples)
 * @tparam row_ts types of the extra interleaved rows (columns), stored at fields 3, 4, ...
 */
template <typename pos_t = uint32_t, move_pos_encoding_t enc = POS, typename... row_ts>
class move_data_structure {
    static_assert(std::is_same_v<pos_t, uint32_t> || std::is_same_v<pos_t, uint64_t>);

public:
    // the encoding of the interval starting positions (D_p vs. D_len + p-samples)
    static constexpr move_pos_encoding_t encoding = enc;
    // the number of extra interleaved rows, and the type of the J-th one
    static constexpr uint8_t num_rows = sizeof...(row_ts);
    template <uint8_t J>
    using row_t = std::tuple_element_t<J, std::tuple<row_ts...>>;

protected:
    class construction;

    using pair_t = std::pair<pos_t, pos_t>; // pair type
    using pair_arr_t = std::vector<pair_t>; // pair array type

    // the interleaved storage holds 3 core fields (D_p/D_len, D_idx, D_offs) followed by the extra rows
    static constexpr uint8_t num_data_fields = 8; // capacity of the interleaved storage
    static_assert(3 + num_rows <= num_data_fields);
    using row_widths_t = std::array<uint8_t, num_rows>; // widths (in bits) of the extra rows

    pos_t n = 0; // n = p_{k_'-1} + d_{k_'-1}, k_' <= n
    pos_t k = 0; // k, number of intervals in the original disjoint inteval sequence I
    pos_t k_ = 0; // k', number of intervals in the balanced disjoint inteval sequence B_a(I), k <= k_'
    uint16_t a = 0; // balancing parameter, restricts the number of intervals in the resulting move data structure to k*(a/(a-1))
    uint8_t omega_p = 0; // word width of D_p (and of the p-samples in the differential encoding)
    uint8_t omega_idx = 0; // word width of D_idx
    uint8_t omega_offs = 0; // word width of D_offs (and of D_len in the differential encoding)
    interleaved_bit_aligned_vectors<pos_t> data; // interleaved: (D_p | D_len), D_idx, D_offs, then the extra rows

    // builds the interleaved data layout from omega_p = log n, omega_idx = log k', omega_offs (already set) and the extra-row widths
    inline void init_data(row_widths_t row_widths)
    {
        omega_p = std::bit_width(uint64_t(n));
        omega_idx = std::bit_width(uint64_t(k_));
        std::array<uint8_t, num_data_fields> w {};
        w[0] = (enc == POS) ? omega_p : omega_offs; // D_p (positional) or D_len (differential)
        w[1] = omega_idx;
        w[2] = omega_offs;
        for (uint8_t j = 0; j < num_rows; j++) w[3 + j] = row_widths[j];
        data = interleaved_bit_aligned_vectors<pos_t>(w);
    }

    pos_t d = 0; // sampling parameter of the p-samples (only used in the differential encoding); S[j] = p_{j*d}
    interleaved_bit_aligned_vectors<pos_t> p_smpl; // [0..k_'/d] p-samples S[j] = p_{j*d} (only used in the differential encoding)

    // the j-th p-sample S[j] = p_{j*d}
    inline pos_t p_sample(pos_t j) const requires(enc == DIFF) { return p_smpl.template get<0, pos_t>(j); }

public:
    move_data_structure() = default;

    /**
     * @brief constructs a new move data structure from a disjoint interval sequence I (consumed if passed as an rvalue)
     * @param I a disjoint interval sequence
     * @param n n = p_k + d_j
     * @param params construction parameters (params.d = p-sample sampling parameter for the DIFF encoding)
     * @param row_widths widths (in bits) of the extra rows
     * @param pi_mphi vector to move pi into after the construction
     */
    move_data_structure(pair_arr_t& I, pos_t n, mds_params params = {}, row_widths_t row_widths = {}, std::vector<pos_t>* pi_mphi = nullptr)
    {
        if constexpr (enc == DIFF) this->d = params.d;
        construction(*this, I, n, false, row_widths, params, pi_mphi);
    }

    move_data_structure(pair_arr_t&& I, pos_t n, mds_params params = {}, row_widths_t row_widths = {}, std::vector<pos_t>* pi_mphi = nullptr)
    {
        if constexpr (enc == DIFF) this->d = params.d;
        construction(*this, I, n, true, row_widths, params, pi_mphi);
    }

    /**
     * @brief returns the size of the data structure in bytes
     * @return size of the data structure in bytes
     */
    uint64_t size_in_bytes() const
    {
        uint64_t size = sizeof(*this) + data.size_in_bytes();
        if constexpr (enc == DIFF) size += p_smpl.size_in_bytes();
        return size;
    }

    inline pos_t max_value() const { return n; } // n = p_k + d_k
    inline pos_t num_intervals() const { return k_; } // k', the number of intervals
    inline uint16_t balancing_parameter() const { return a; }
    inline bool empty() const { return num_intervals() == 0; }
    inline uint8_t width_p() const { return omega_p; } // bits per D_p entry (and per p-sample)
    inline uint8_t width_idx() const { return omega_idx; } // bits per D_idx entry
    inline uint8_t width_offs() const { return omega_offs; } // bits per D_offs entry
    inline pos_t sampling() const requires(enc == DIFF) { return d; } // p-sample sampling parameter

    // generic access to the extra interleaved rows: value / setter / bit width of the J-th row at position x
    template <uint8_t J> inline row_t<J> row(pos_t x) const { return data.template get<3 + J, row_t<J>>(x); }
    template <uint8_t J> inline void set_row(pos_t x, row_t<J> v) { data.template set_parallel<3 + J, row_t<J>>(x, v); }
    template <uint8_t J> inline uint8_t width_row() const { return data.template width<3 + J>(); }

protected:
    inline pos_t num_p_samples() const requires(enc == DIFF) { return p_smpl.size(); }

    // sets the sentinel entry (index k') of every extra row to a default-constructed value
    inline void init_row_sentinels()
    {
        for_constexpr<0, num_rows, 1>([&](auto i) { set_row<i>(k_, row_t<i> {}); });
    }

    /**
     * @brief resizes to k_ intervals and writes the sentinel entry (omega_offs and, differentially, d must
     *        already be set); the construction fills the interval data (differentially via finalize_differential)
     * @param n maximum value
     * @param k_ size
     * @param row_widths bit widths of the extra rows
     */
    void resize(pos_t n, pos_t k_, row_widths_t row_widths)
    {
        this->n = n;
        this->k_ = k_;
        init_data(row_widths);
        data.resize_no_init(k_ + 1);

        if constexpr (enc == POS) set_p(k_, n);
        else set_len(k_, 0);
        set_idx(k_, k_);
        set_offs(k_, 0);
        init_row_sentinels();

        if constexpr (enc == DIFF) {
            p_smpl = interleaved_bit_aligned_vectors<pos_t>({ omega_p });
            p_smpl.resize_no_init(k_ / d + 1);
        }
    }

    inline void set_p(pos_t j, pos_t p) requires(enc == POS) { data.template set_parallel<0, pos_t>(j, p); }
    inline void set_len(pos_t j, pos_t len) requires(enc == DIFF) { data.template set_parallel<0, pos_t>(j, len); }
    inline void set_idx(pos_t j, pos_t idx) { data.template set_parallel<1, pos_t>(j, idx); }
    inline void set_offs(pos_t j, pos_t offs) { data.template set_parallel<2, pos_t>(j, offs); }

public:
    /**
     * @brief returns D_p[x], the start of the x-th input interval; differentially reconstructed in O(d)
     * @param x [0..k_']
     * @return D_p[x]
     */
    inline pos_t p(pos_t x) const
    {
        if constexpr (enc == POS) {
            return data.template get<0, pos_t>(x);
        } else {
            pos_t j = x / d;
            pos_t p = p_sample(j);
            for (pos_t y = j * d; y < x; y++) p += len(y);
            return p;
        }
    }

    // len_x = p_{x+1} - p_x
    inline pos_t len(pos_t x) const
    {
        if constexpr (enc == POS) return data.template get<0, pos_t>(x + 1) - data.template get<0, pos_t>(x);
        else return data.template get<0, pos_t>(x);
    }

    inline pos_t q(pos_t x) const { return p(idx(x)) + offs(x); } // q_x
    inline pos_t idx(pos_t x) const { return data.template get<1, pos_t>(x); } // D_idx[x]
    inline pos_t offs(pos_t x) const { return data.template get<2, pos_t>(x); } // D_offs[x]
    inline pos_t pos(pos_t x, pos_t o) const { return p(x) + o; } // absolute position of the (run,offset)-pair (x,o)

    /**
     * @brief returns (p(x2) + o2) - (p(x1) + o1) for two (run,offset)-pairs with the second position >= the
     *        first. In the differential encoding, if the runs are close (<= 2*d apart) the intervening
     *        interval lengths are accumulated (cache-friendly, no p-sample search), otherwise p(x1) and
     *        p(x2) are reconstructed separately; in the positional encoding it is a direct subtraction.
     * @param x1 first input interval index
     * @param o1 offset within the x1-th input interval
     * @param x2 second input interval index (with p(x2) + o2 >= p(x1) + o1)
     * @param o2 offset within the x2-th input interval
     * @return (p(x2) + o2) - (p(x1) + o1)
     */
    inline pos_t distance(pos_t x1, pos_t o1, pos_t x2, pos_t o2) const
    {
        if constexpr (enc == POS) {
            return (p(x2) + o2) - (p(x1) + o1);
        } else {
            if (x2 - x1 <= 2 * d) {
                pos_t s = o2;
                for (pos_t r = x1; r < x2; r++) s += len(r);
                return s - o1;
            } else {
                return (p(x2) + o2) - (p(x1) + o1);
            }
        }
    }

    /**
     * @brief returns the (run, offset)-pair (x, i - p(x)) with p(x) <= i < p(x+1) of an absolute position i
     *        (binary search over D_p in the positional encoding; over the checkpoints + a D_len scan differentially)
     * @param i [0..n-1] an absolute position
     * @return the (run,offset)-pair (x, i - p(x))
     */
    inline pair_t run_and_offset(pos_t i) const
    {
        if constexpr (enc == POS) {
            pos_t x = interval_index(i);
            return { x, i - p(x) };
        } else {
            pos_t j = bin_search_max_leq<pos_t>(i, pos_t(0), num_p_samples() - 1, [&](pos_t s) { return p_sample(s); });
            pos_t x = d * j;
            pos_t p = p_sample(j);
            pos_t l;

            while (p + (l = len(x)) <= i) {
                p += l;
                x++;
            }

            return { x, i - p };
        }
    }

    /**
     * @brief returns the index of the input interval containing the absolute position i
     * @param i [0..n-1] an absolute position
     * @return x with p(x) <= i < p(x+1)
     */
    inline pos_t interval_index(pos_t i) const
    {
        if constexpr (enc == POS) {
            return bin_search_max_leq<pos_t>(i, pos_t(0), k_ - 1, [&](pos_t x) { return p(x); });
        } else {
            return run_and_offset(i).first;
        }
    }

    /**
     * @brief (positional encoding) move query on an absolute position: changes (i, x), i in [p(x), p(x+1)), to
     *        (i', x') with i' = f_I(i) and i' in [p(x'), p(x'+1))
     * @param i [0..n-1] absolute position
     * @param x [0..k_'-1] the input interval index containing i
     */
    inline void move(pos_t& i, pos_t& x) const
        requires(enc == POS)
    {
        i = q(x) + (i - p(x));
        x = idx(x);
        while (i >= p(x + 1)) {
            x++;
        }
    }

    /**
     * @brief (differential encoding) move query on an (offset, run)-pair: changes (o, x) to (o', x') with
     *        p(x') + o' = f_I(p(x) + o) and p(x') <= f_I(p(x) + o) < p(x'+1); never touches p -- it uses only
     *        offs, idx and len. The component-first order matches the positional move(i, x).
     * @param o offset within the x-th input interval
     * @param x [0..k_'-1] input interval index
     */
    inline void move(pos_t& o, pos_t& x) const
        requires(enc == DIFF)
    {
        o += offs(x);
        x = idx(x);
        pos_t l;

        while (o >= (l = len(x))) {
            o -= l;
            x++;
        }
    }

    /**
     * @brief move query on a (component, run)-pair -- (i, x) in the positional encoding, (o, x) in the differential
     *        encoding; forwards to the two-argument move and returns the updated pair
     * @param pair the (component, run)-pair to move
     * @return the moved pair
     */
    inline pair_t move(pair_t pair) const
    {
        move(pair.first, pair.second);
        return pair;
    }

    /**
     * @brief logs all data structures
     */
    void log_data_structures() const
    {
        aligned_log log;
        log.add_row("p:", k_, [&](pos_t i) { return p(i); });
        log.add_row("q:", k_, [&](pos_t i) { return q(i); });
        log.add_row("idx:", k_, [&](pos_t i) { return idx(i); });
        log.print();
    }

    /**
     * @brief serializes the move data structure to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) const
    {
        out.write((char*) &n, sizeof(pos_t));
        out.write((char*) &k, sizeof(pos_t));
        out.write((char*) &k_, sizeof(pos_t));
        out.write((char*) &a, sizeof(uint16_t));
        data.serialize(out);

        if constexpr (enc == DIFF) {
            out.write((char*) &d, sizeof(pos_t));
            p_smpl.serialize(out);
        }
    }

    /**
     * @brief loads the move data structure from an input stream
     * @param in input stream
     */
    void load(std::istream& in)
    {
        in.read((char*) &n, sizeof(pos_t));
        in.read((char*) &k, sizeof(pos_t));
        in.read((char*) &k_, sizeof(pos_t));
        in.read((char*) &a, sizeof(uint16_t));
        data.load(in);

        omega_p = std::bit_width(uint64_t(n));
        omega_idx = std::bit_width(uint64_t(k_));
        omega_offs = data.template width<2>();

        if constexpr (enc == DIFF) {
            in.read((char*) &d, sizeof(pos_t));
            p_smpl.load(in);
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

#include "construction.hpp"
