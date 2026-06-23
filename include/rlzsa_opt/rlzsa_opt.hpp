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

#include <cstdint>
#include <iostream>
#include <vector>

#include <data_structures/interleaved_bit_aligned_vectors.hpp>
#include <data_structures/plain_bit_vector.hpp>
#include <data_structures/sd_array.hpp>
#include <misc/utils.hpp>
#include <misc/files.hpp>
#include <misc/log.hpp>

/**
 * @brief relative-Lempel-Ziv encoding of the differential suffix array (SA^d), optimized for move-r.
 * @tparam pos_t index integer type (use uint32_t if input size < UINT_MAX, else uint64_t)
 */
template <typename pos_t = uint32_t>
class rlzsa_opt {
public:
    static_assert(
        std::is_same_v<pos_t, uint32_t> ||
        std::is_same_v<pos_t, uint64_t>);

    // sample rate of the copy phrases in the rlzsa
    static constexpr pos_t sample_rate_scp = 4;

protected:
    pos_t n = 0; // n, the length of the input (including the terminator symbol)
    pos_t z = 0; // z, the number of phrases in the rlzsa
    pos_t z_l = 0; // z_l, the number of literal phrases in the rlzsa
    pos_t z_c = 0; // z_c, the number of copy-phrases in the rlzsa

    // reference for SA^d (differential suffix array)
    interleaved_bit_aligned_vectors<uint64_t> _R;
    // bit vector storing the phrase types of the rlzsa, i.e, PT[i] = 1 <=> phrase i is literal
    plain_bit_vector<pos_t, true, true, true> _PT;
    // compressed bit vector marking the sampled starting positions in SA^d of the copy phrases of the rlzsa
    sd_array<pos_t> _SCP_S;
    // lengths of the copy phrases of the rlzsa
    std::vector<uint16_t> _CPL;
    // starting positions in R of the copy phrases of the rlzsa
    interleaved_bit_aligned_vectors<pos_t> _SR;
    // literal phrases of the rlzsa
    interleaved_bit_aligned_vectors<pos_t> _LP;

public:
    rlzsa_opt() = default;

    /**
     * @brief constructs the rlzsa from its (already built) components
     * @param n length of the input (including the terminator symbol)
     * @param z number of phrases
     * @param z_l number of literal phrases
     * @param z_c number of copy phrases
     * @param R reference for SA^d
     * @param PT phrase types
     * @param SCP_S sampled starting positions of the copy phrases
     * @param CPL copy phrase lengths
     * @param SR starting positions in R of the copy phrases
     * @param LP literal phrases
     */
    rlzsa_opt(
        pos_t n, pos_t z, pos_t z_l, pos_t z_c,
        interleaved_bit_aligned_vectors<uint64_t>&& R,
        plain_bit_vector<pos_t, true, true, true>&& PT,
        sd_array<pos_t>&& SCP_S,
        std::vector<uint16_t>&& CPL,
        interleaved_bit_aligned_vectors<pos_t>&& SR,
        interleaved_bit_aligned_vectors<pos_t>&& LP)
        : n(n), z(z), z_l(z_l), z_c(z_c)
        , _R(std::move(R))
        , _PT(std::move(PT))
        , _SCP_S(std::move(SCP_S))
        , _CPL(std::move(CPL))
        , _SR(std::move(SR))
        , _LP(std::move(LP))
    { }

    /**
     * @brief constructs the rlzsa
     * @tparam sad_func_t type of the caller-provided function returning SA^d[i] (callable as SAd(i_p, i))
     */
    template <typename sad_func_t> class construction;

    // ############################# ACCESS METHODS #############################

    /**
     * @brief returns whether the rlzsa is empty (i.e, has not been built)
     */
    inline bool empty() const
    {
        return z == 0;
    }

    /**
     * @brief returns the length of the input (including the terminator symbol)
     */
    inline pos_t input_size() const
    {
        return n;
    }

    /**
     * @brief returns the number of phrases in the rlzsa
     */
    inline pos_t num_phrases() const
    {
        return z;
    }

    /**
     * @brief returns the number of literal phrases in the rlzsa
     */
    inline pos_t num_literal_phrases() const
    {
        return z_l;
    }

    /**
     * @brief returns the number of copy phrases in the rlzsa
     */
    inline pos_t num_copy_phrases() const
    {
        return z_c;
    }

    /**
     * @brief returns a reference to R
     */
    inline const interleaved_bit_aligned_vectors<uint64_t>& R() const
    {
        return _R;
    }

    /**
     * @brief returns a reference to PT
     */
    inline const plain_bit_vector<pos_t, true, true, true>& PT() const
    {
        return _PT;
    }

    /**
     * @brief returns a reference to CPL
     */
    inline const std::vector<uint16_t>& CPL() const
    {
        return _CPL;
    }

    /**
     * @brief returns a reference to SCP_S
     */
    inline const sd_array<pos_t>& SCP_S() const
    {
        return _SCP_S;
    }

    /**
     * @brief returns a reference to SR
     */
    inline const interleaved_bit_aligned_vectors<pos_t>& SR() const
    {
        return _SR;
    }

    /**
     * @brief returns a reference to LP
     */
    inline const interleaved_bit_aligned_vectors<pos_t>& LP() const
    {
        return _LP;
    }

    /**
     * @brief returns R[x]
     * @param x [0..|R|-1] index in R
     */
    inline uint64_t R(pos_t x) const
    {
        return _R[x];
    }

    /**
     * @brief returns PT[x]
     * @param x [0..z-1] index in PT
     */
    inline bool PT(pos_t x) const
    {
        return _PT[x];
    }

    /**
     * @brief returns CPL[x]
     * @param x [0..z_c-1] index in CPL
     */
    inline uint16_t CPL(pos_t x) const
    {
        return _CPL[x];
    }

    /**
     * @brief returns SCP_S[x]
     * @param x [0..z_c/sample_rate_scp-1] index in SCP_S
     */
    inline pos_t SCP_S(pos_t x) const
    {
        return _SCP_S.select_1(x + 1);
    }

    /**
     * @brief returns SR[x]
     * @param x [0..z_c-1] index in SR
     */
    inline pos_t SR(pos_t x) const
    {
        return _SR[x];
    }

    /**
     * @brief returns LP[x]
     * @param x [0..z_l-1] index in LP
     */
    inline pos_t LP(pos_t x) const
    {
        return _LP[x];
    }

    // ############################# DECODING #############################

    /**
     * @brief stores the variables needed to decode consecutive suffix-array values from the rlzsa.
     *
     * A decode context tracks a current position i in the suffix array, the current value s = SA[i] and
     * the internal phrase-navigation state of the rlzsa. After initializing it to a position with
     * init_right(), the suffix array can be decoded to the right (next/skip_right/report_right), or, after
     * init_left(), to the left (prev/skip_left/report_left). The position i and the value s are readable from outside
     * through pos() and value() and writable through set_pos() and set_value() (e.g, to seed s before
     * decoding); the remaining navigation state is internal.
     */
    struct decode_context_t {
      protected:
        const rlzsa_opt<pos_t>* rlz = nullptr; // the rlzsa to decode from

        pos_t i; // current position in the suffix array
        pos_t s; // current suffix s = SA[i]

        pos_t x_p;  // phrase-index of the phrase of the rlzsa containing i
        pos_t x_lp; // literal-phrase index of the current or next literal phrase of the rlzsa
        pos_t x_cp; // copy-phrase index of the current or next copy-phrase of the rlzsa
        pos_t x_r;  // position in R inside the current copy-phrase (or the starting position in R of the current or next copy phrase) of the rlzsa
        pos_t s_np; // starting position in the rlzsa of the next phrase (or the current phrase when decoding to the left) of the rlzsa

      public:
        /**
         * @brief constructs an empty decode context (not bound to any rlzsa)
         */
        decode_context_t() = default;

        /**
         * @brief constructs a decode context for the rlzsa rlz
         * @param rlz an rlzsa
         */
        decode_context_t(const rlzsa_opt<pos_t>& rlz)
            : rlz(&rlz)
        { }

        /**
         * @brief returns the current position i in the suffix array
         * @return i (the position the context is currently at)
         */
        pos_t pos() const
        {
            return i;
        }

        /**
         * @brief sets the current position i in the suffix array
         * @param p the position to set i to
         */
        inline void set_pos(pos_t p)
        {
            i = p;
        }

        /**
         * @brief returns the current suffix-array value s = SA[i]
         * @return s (the current suffix array value)
         */
        pos_t value() const
        {
            return s;
        }

        /**
         * @brief sets the current suffix-array value s = SA[i] (e.g, to seed the context before decoding)
         * @param v the value to set s to
         */
        inline void set_value(pos_t v)
        {
            s = v;
        }

        /**
         * @brief prepares the context to decode SA[i] and advance to the right (i.e, with next()/
         * skip_right()/report_right()); seed s with set_value(SA[i-1]) before decoding
         * @param i position in the suffix array to initialize the context to
         */
        inline void init_right(pos_t i);

        /**
         * @brief prepares the context to decode SA[i] and advance to the left (i.e, with prev()/
         * skip_left()/report_left()); seed s with set_value(SA[i]) before decoding. Equivalent to
         * init_right(i) followed by turn_left().
         * @param i position in the suffix array to initialize the context to
         */
        inline void init_left(pos_t i);

        /**
         * @brief reverses the decode direction at the current position from rightward to leftward in O(1):
         * a right-ready context (init_right(), or a rightward step) becomes left-ready (prev()/skip_left()/
         * report_left()). pos() is left unchanged; seed s for the new direction with set_value(SA[pos()]).
         * The inverse of turn_right().
         */
        inline void turn_left();

        /**
         * @brief reverses the decode direction at the current position from leftward to rightward in O(1):
         * a left-ready context (init_left(), or a leftward step) becomes right-ready (next()/skip_right()/
         * report_right()). pos() is left unchanged; seed s for the new direction with set_value(SA[pos()-1]).
         * The inverse of turn_left().
         */
        inline void turn_right();

        /**
         * @brief advances the context to the right up to position e, without reporting the values in
         * between; the context must be prepared to advance to the right (i.e, init_right() with s seeded to
         * SA[i-1], or a previous right-decoding operation). Afterwards i == e and the context is prepared
         * to decode SA[e] (i.e, s == SA[e-1], and a following next()/report_right() yields SA[e]). This is
         * the rightward counterpart of skip_left().
         * @param e position in the suffix array to advance to (e >= i)
         */
        inline void skip_right(pos_t e);

        /**
         * @brief advances the context to the left down to position b, without reporting the values in
         * between; the context must be prepared to advance to the left (i.e, init_left(), or a previous
         * left-decoding operation). Afterwards i == b, s == SA[b] and the
         * context is prepared to decode SA[b-1] (i.e, a following prev()/report_left() yields SA[b-1]).
         * This is the leftward counterpart of skip_right().
         * @param b position in the suffix array to advance to (b <= i)
         */
        inline void skip_left(pos_t b);

        /**
         * @brief decodes and stores SA[i] in s and prepares the context to decode SA[i-1];
         * the context must be prepared to decode SA[i]
         */
        inline pos_t prev();

        /**
         * @brief decodes and stores SA[i] in s and prepares the context to decode SA[i+1];
         * the context must be prepared to decode SA[i]
         */
        inline pos_t next();

        /**
         * @brief decodes the suffix array values SA[b..i] from right to left, reporting each value
         * @tparam report_fnc_t type of the report function
         * @param b left position in the suffix array to decode to
         * @param report function called with each value SA[j] (or (j, SA[j])) for j in [b,i], from right to left
         */
        template <typename report_fnc_t>
        inline void report_left(pos_t b, report_fnc_t report);

        /**
         * @brief decodes the suffix array values SA[i..e] from left to right, reporting each value
         * @tparam report_fnc_t type of the report function
         * @param e right position in the suffix array to decode to
         * @param report function called with each value SA[j] (or (j, SA[j])) for j in [i,e], from left to right
         */
        template <typename report_fnc_t>
        inline void report_right(pos_t e, report_fnc_t report);
    };

    /**
     * @brief returns a decode context for the rlzsa
     * @return decode_context_t
     */
    inline decode_context_t decode() const
    {
        return decode_context_t(*this);
    }

    // ############################# MISC METHODS #############################

    /**
     * @brief returns the size of the data structure in bytes
     */
    uint64_t size_in_bytes() const
    {
        return sizeof(this) +
            _R.size_in_bytes() +
            (z_c + 2) * sizeof(uint16_t) + // CPL
            _SCP_S.size_in_bytes() +
            _SR.size_in_bytes() +
            _LP.size_in_bytes() +
            _PT.size_in_bytes();
    }

    /**
     * @brief logs the contents of the rlzsa (only useful for small inputs)
     */
    void log_data_structures() const
    {
        std::cout << "z = " << z << ", ";
        std::cout << "z_l = " << z_l << ", ";
        std::cout << "z_c = " << z_c << std::endl;
        std::cout << "SCP_S sampling rate = " << sample_rate_scp << std::endl;
        std::cout << std::endl;

        aligned_log log_r;
        log_r.add_row("R:", _R.size(), [&](pos_t i) { return int32_t{_R[i]} - int32_t{n}; });
        log_r.print();

        log_indexed("PT:", _PT);
        log_indexed("SCP_S:", _SCP_S);
        log_indexed("CPL:", _CPL);
        log_indexed("SR:", _SR);

        aligned_log log_lp;
        log_lp.add_row("LP:", _LP.size(), [&](pos_t i) { return int32_t{_LP[i]} - int32_t{n}; });
        log_lp.print();
    }

    /**
     * @brief logs the data structure sizes to cout
     */
    void log_data_structure_sizes() const
    {
        std::cout << "R: " << format_size(_R.size_in_bytes()) << std::endl;
        std::cout << "CPL: " << format_size((z_c + 2) * sizeof(uint16_t)) << std::endl;
        std::cout << "SCP_S: " << format_size(_SCP_S.size_in_bytes()) << std::endl;
        std::cout << "SR: " << format_size(_SR.size_in_bytes()) << std::endl;
        std::cout << "LP: " << format_size(_LP.size_in_bytes()) << std::endl;
        std::cout << "PT: " << format_size(_PT.size_in_bytes()) << std::endl;
    }

    /**
     * @brief logs the data structure sizes to the output stream out
     * @param out an output stream
     * @param suffix suffix appended to each measurement key
     */
    void log_data_structure_sizes(std::ostream& out, std::string suffix = "") const
    {
        out << " size_r" << suffix << "=" << _R.size_in_bytes();
        out << " size_cpl" << suffix << "=" << (z_c + 2) * sizeof(uint16_t);
        out << " size_scp" << suffix << "=" << _SCP_S.size_in_bytes();
        out << " size_sr" << suffix << "=" << _SR.size_in_bytes();
        out << " size_lp" << suffix << "=" << _LP.size_in_bytes();
        out << " size_pt" << suffix << "=" << _PT.size_in_bytes();
    }

    // ############################# SERIALIZATION METHODS #############################

    /**
     * @brief stores the rlzsa to an output stream
     * @param out output stream to store the rlzsa to
     */
    void serialize(std::ostream& out) const
    {
        out.write((char*) &n, sizeof(pos_t));
        out.write((char*) &z, sizeof(pos_t));
        out.write((char*) &z_l, sizeof(pos_t));
        out.write((char*) &z_c, sizeof(pos_t));

        _R.serialize(out);
        _SCP_S.serialize(out);
        write_to_file(out, (char*) &_CPL[0], (z_c + 2) * sizeof(uint16_t));
        _SR.serialize(out);
        _LP.serialize(out);
        _PT.serialize(out);
    }

    /**
     * @brief reads a serialized rlzsa from an input stream
     * @param in an input stream storing a serialized rlzsa
     */
    void load(std::istream& in)
    {
        in.read((char*) &n, sizeof(pos_t));
        in.read((char*) &z, sizeof(pos_t));
        in.read((char*) &z_l, sizeof(pos_t));
        in.read((char*) &z_c, sizeof(pos_t));

        _R.load(in);
        _SCP_S.load(in);
        no_init_resize(_CPL, z_c + 2);
        read_from_file(in, (char*) &_CPL[0], (z_c + 2) * sizeof(uint16_t));
        _SR.load(in);
        _LP.load(in);
        _PT.load(in);
    }
};

#include "queries.tpp"
