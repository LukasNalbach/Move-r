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

#ifndef MOVE_R_DATA_STRUCTURES_FASTA_SEQUENCE_DATA_HPP
#define MOVE_R_DATA_STRUCTURES_FASTA_SEQUENCE_DATA_HPP

#include <bit>
#include <cstdint>
#include <istream>
#include <ostream>
#include <string>
#include <vector>

#include <data_structures/interleaved_bit_aligned_vectors.hpp>
#include <misc/files.hpp>
#include <misc/search.hpp>

// the byte that replaces each header line and each masked run of a FASTA file (see pfp_rlbwt::fasta_transform)
inline constexpr char fasta_separator_symbol = (char) 0x7F;

/**
 * @brief stores the sequences of a multi-sequence (FASTA) input: the starting position of each sequence in the
 *        text, its name and the symbol separating the sequences; empty (has_sequences() is false) if the
 *        index was not built from a multi-sequence input
 */
class fasta_sequence_data {
protected:
    // ############################# VARIABLES #############################

    uint64_t _n = 0; // n, the length of the text

    /* the starting positions of the sequences: [0..num_sequences()-1] _seq_starts[i] is the position of the first
       symbol of sequence i in the text; strictly increasing, bit_width(n) bits per entry */
    interleaved_bit_aligned_vectors _seq_starts;

    // the names of the sequences: [0..num_sequences()-1] _seq_names[i] is the name of sequence i (RNAME in SAM output)
    std::vector<std::string> _seq_names;

    uint8_t _separator_sym {}; // symbol separating the sequences; a search never extends with it

public:
    // ############################# MISC PUBLIC METHODS #############################

    // returns whether the index was built from a multi-sequence input
    inline bool has_sequences() const { return !_seq_names.empty(); }

    // returns the number of sequences (0 if !has_sequences())
    inline uint64_t num_sequences() const { return _seq_names.size(); }

    // returns the separator symbol (undefined if !has_sequences())
    inline uint8_t separator_sym() const { return _separator_sym; }

    // returns the size in bytes
    uint64_t size_in_bytes() const
    {
        uint64_t bytes = _seq_starts.size_in_bytes();
        for (const std::string& name : _seq_names) bytes += name.size();
        return bytes;
    }

    // ############################# ACCESS METHODS #############################

    /**
     * @brief returns the starting position of sequence i in the text
     * @param i [0..num_sequences()-1] sequence index
     */
    inline uint64_t sequence_start(uint64_t i) const { return _seq_starts[i]; }

    /**
     * @brief returns the length of sequence i (without the separator symbol)
     * @param i [0..num_sequences()-1] sequence index
     */
    inline uint64_t sequence_length(uint64_t i) const
    {
        // position after sequence i (of the separator before sequence i+1, or of the terminator)
        uint64_t end = (i + 1 < num_sequences()) ? sequence_start(i + 1) - 1 : _n - 1;
        return end - sequence_start(i);
    }

    /**
     * @brief returns the name of sequence i (RNAME in SAM output)
     * @param i [0..num_sequences()-1] sequence index
     */
    inline const std::string& sequence_name(uint64_t i) const { return _seq_names[i]; }

    // ############################# QUERY METHODS #############################

    /**
     * @brief returns the index of the sequence containing a text position, using an exponential search that
     *        starts at the sequence of the previous query
     * @param pos [0..n-1] position in the text
     * @param hint index of the sequence returned by the previous query (0 before the first query); is set to the
     *        returned index
     */
    inline uint64_t sequence_index(uint64_t pos, uint64_t& hint) const
    {
        auto start = [&](uint64_t i){ return _seq_starts[i]; };
        hint = _seq_starts[hint] <= pos
            ? exp_search_max_leq<RIGHT>(pos, hint, num_sequences() - 1, start)  // pos is in sequence hint or later
            : exp_search_max_leq<LEFT>(pos, 0, hint, start);                    // pos is in an earlier sequence
        return hint;
    }

    // ############################# CONSTRUCTION METHODS #############################

    /**
     * @brief sets the sequences (see process_fasta)
     * @param n length of the text (input_size() + 1)
     * @param seq_starts starting positions of the sequences in the text (strictly increasing, each < n)
     * @param seq_names names of the sequences (same size as seq_starts)
     * @param separator symbol separating the sequences
     */
    void set_sequences(uint64_t n, const std::vector<uint64_t>& seq_starts, std::vector<std::string> seq_names, uint8_t separator)
    {
        *this = fasta_sequence_data();
        _n = n;
        _seq_names = std::move(seq_names);

        if (_seq_names.empty()) {
            return;
        }

        _separator_sym = separator;
        _seq_starts = interleaved_bit_aligned_vectors({ uint8_t(std::bit_width(n)) });
        _seq_starts.resize_no_init(seq_starts.size());
        for (uint64_t i = 0; i < seq_starts.size(); i++) _seq_starts.template set<0>(i, seq_starts[i]);
    }

    // ############################# SERIALIZATION METHODS #############################

    /**
     * @brief writes the sequence data to an output stream (only the number of sequences, if it is 0)
     * @param out output stream
     */
    void serialize(std::ostream& out) const
    {
        uint64_t num_seq = _seq_names.size();
        out.write((char*) &num_seq, sizeof(uint64_t));

        if (num_seq > 0) {
            out.write((char*) &_n, sizeof(uint64_t));
            out.write((char*) &_separator_sym, sizeof(uint8_t));
            _seq_starts.serialize(out);

            for (const std::string& name : _seq_names) write_string(out, name);
        }
    }

    /**
     * @brief reads the sequence data written by serialize from an input stream
     * @param in input stream
     */
    void load(std::istream& in)
    {
        *this = fasta_sequence_data();
        uint64_t num_seq;
        in.read((char*) &num_seq, sizeof(uint64_t));

        if (num_seq > 0) {
            in.read((char*) &_n, sizeof(uint64_t));
            in.read((char*) &_separator_sym, sizeof(uint8_t));
            _seq_starts.load(in);
            _seq_names.resize(num_seq);

            for (std::string& name : _seq_names) read_string(in, name);
        }
    }
};

#endif // MOVE_R_DATA_STRUCTURES_FASTA_SEQUENCE_DATA_HPP
