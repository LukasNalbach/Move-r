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

#include <bit>
#include <cstdint>
#include <istream>
#include <ostream>
#include <string>
#include <vector>

#include <data_structures/interleaved_bit_aligned_vectors.hpp>
#include <misc/search.hpp>

/**
 * @brief the multi-sequence (FASTA/DNA) metadata shared by move_r and move_rb, which each hold one as a member: for
 *        each sequence of a multi-FASTA input, the start position of its body in the text and its name (RNAME for SAM
 *        output), plus the separator symbol delimiting the sequences (which a search never
 *        extends into, so no match ever crosses a sequence boundary). The start positions are held in one bit-packed,
 *        strictly-increasing array; mapping a reported occurrence to its sequence is an exponential search over it,
 *        which -- resumed from the previous lookup via the cursor overload -- amortizes to ~O(1) per occurrence when
 *        consecutive occurrences fall in the same or nearby sequences.
 *
 *        The metadata is empty unless the index was built from a multi-sequence input; then the whole text is treated
 *        as one anonymous sequence, has_sequences() is false and no sequence-relative coordinates or names exist.
 * @tparam pos_t the index integer type (matching the owning index)
 * @tparam sym_t the external symbol type (matching the owning index)
 */
template <typename pos_t, typename sym_t>
class fasta_sequence_data {
protected:
    pos_t _seq_n = 0;                                   // the text length n (needed for the last sequence's end)
    interleaved_bit_aligned_vectors<pos_t> _seq_starts; // bit-packed, sorted body start positions of the sequences
    std::vector<std::string> _seq_names;                // _seq_names[i] = the name (SAM RNAME) of sequence i
    sym_t _separator_sym {};                            // the symbol delimiting the sequences; a search never extends
                                                        // into it, so no search ever crosses a sequence boundary

public:
    /**
     * @brief returns whether the index carries multi-sequence (FASTA/DNA) boundary metadata; if not, the input is
     *        treated as a single anonymous sequence and no sequence-relative coordinates or names are available
     * @return whether the index was built from a multi-sequence input
     */
    inline bool has_sequences() const
    {
        return !_seq_names.empty();
    }

    /**
     * @brief returns the number of sequences the index was built from (0 unless the index has multi-sequence metadata)
     * @return the number of sequences
     */
    inline pos_t num_sequences() const
    {
        return _seq_names.size();
    }

    /**
     * @brief returns the index of the sequence containing the global text position pos
     * @param pos a position in [0, n)
     * @return the index of the sequence containing pos
     */
    inline pos_t sequence_index(pos_t pos) const
    {
        // the sequence containing pos is the last one whose start is <= pos (galloping from the start)
        return exp_search_max_leq<pos_t, RIGHT>(pos, 0, num_sequences() - 1, [&](pos_t i){ return _seq_starts[i]; });
    }

    /**
     * @brief returns the index of the sequence containing pos, resuming the search from @p cursor (the sequence index
     *        of the previously looked-up position) with an incremental exponential search and updating it. Amortizes to
     *        ~O(1) per occurrence when consecutive occurrences fall in the same or nearby sequences, and is never worse
     *        than a binary search -- intended for the per-occurrence lookups of a reporting loop.
     * @param pos a position in [0, n)
     * @param cursor in/out: the previous sequence index; updated to the sequence containing pos
     * @return the index of the sequence containing pos
     */
    inline pos_t sequence_index(pos_t pos, pos_t& cursor) const
    {
        auto start = [&](pos_t i){ return _seq_starts[i]; };
        cursor = _seq_starts[cursor] <= pos
            ? exp_search_max_leq<pos_t, RIGHT>(pos, cursor, num_sequences() - 1, start)  // pos is in cursor or later
            : exp_search_max_leq<pos_t, LEFT>(pos, 0, cursor, start);                    // pos is in an earlier sequence
        return cursor;
    }

    /**
     * @brief returns the start position of sequence i's body in the text
     * @param i a sequence index in [0, num_sequences())
     * @return the start position of sequence i in the text
     */
    inline pos_t sequence_start(pos_t i) const
    {
        return _seq_starts[i];
    }

    /**
     * @brief returns the length of sequence i's body (excluding the separator symbol)
     * @param i a sequence index in [0, num_sequences())
     * @return the length of sequence i
     */
    inline pos_t sequence_length(pos_t i) const
    {
        // sequence i ends just before sequence i+1's separator symbol, or at the sentinel for the last sequence
        pos_t end = (i + 1 < num_sequences()) ? sequence_start(i + 1) - 1 : _seq_n - 1;
        return end - sequence_start(i);
    }

    /**
     * @brief returns the name (SAM RNAME) of sequence i
     * @param i a sequence index in [0, num_sequences())
     * @return the name of sequence i
     */
    inline const std::string& sequence_name(pos_t i) const
    {
        return _seq_names[i];
    }

    /**
     * @brief returns the symbol a search must not extend into (the sequence separator), so that no search crosses a
     *        sequence boundary; unspecified (and unused) if the index has no multi-sequence metadata
     * @return the separator symbol
     */
    inline sym_t separator_sym() const
    {
        return _separator_sym;
    }

    /**
     * @brief populates this metadata (typically from process_fasta): stores the sequences' names and bit-packs their
     *        body start positions, and records the separator symbol a search must not cross
     * @param n the length of the text the metadata is built for (n = input_size() + 1)
     * @param seq_starts the body start position of each sequence in the text (strictly increasing, each < n)
     * @param seq_names the name of each sequence (same length as seq_starts)
     * @param separator the separator symbol that delimits the sequences (a search never extends into it)
     */
    void set_sequences(pos_t n, const std::vector<pos_t>& seq_starts, std::vector<std::string> seq_names, sym_t separator)
    {
        _seq_n = n;
        _seq_names = std::move(seq_names);

        if (_seq_names.empty()) {
            return;
        }

        _separator_sym = separator;
        _seq_starts = interleaved_bit_aligned_vectors<pos_t>({ bit_width(n) });
        _seq_starts.resize_no_init(seq_starts.size());
        for (pos_t i = 0; i < seq_starts.size(); i++) _seq_starts.template set<0>(i, seq_starts[i]);
    }

    /**
     * @brief writes the multi-sequence metadata to an output stream (writes a leading count of 0 if there is none)
     * @param out the output stream
     */
    void serialize(std::ostream& out) const
    {
        pos_t num_seq = _seq_names.size();
        out.write((char*) &num_seq, sizeof(pos_t));

        if (num_seq > 0) {
            out.write((char*) &_seq_n, sizeof(pos_t));
            out.write((char*) &_separator_sym, sizeof(sym_t));
            _seq_starts.serialize(out);

            for (const std::string& name : _seq_names) {
                pos_t len = name.size();
                out.write((char*) &len, sizeof(pos_t));
                out.write(name.data(), len);
            }
        }
    }

    /**
     * @brief reads the multi-sequence metadata from an input stream (as written by serialize, which always writes at
     *        least the leading count, so this is called unconditionally after the index's data structures)
     * @param in the input stream
     */
    void load(std::istream& in)
    {
        pos_t num_seq;
        in.read((char*) &num_seq, sizeof(pos_t));

        if (num_seq > 0) {
            in.read((char*) &_seq_n, sizeof(pos_t));
            in.read((char*) &_separator_sym, sizeof(sym_t));
            _seq_starts.load(in);
            _seq_names.resize(num_seq);

            for (pos_t i = 0; i < num_seq; i++) {
                pos_t len;
                in.read((char*) &len, sizeof(pos_t));
                _seq_names[i].resize(len);
                in.read(_seq_names[i].data(), len);
            }
        }
    }

    /**
     * @brief returns the size of the multi-sequence metadata in bytes
     * @return size of the multi-sequence metadata in bytes
     */
    uint64_t size_in_bytes() const
    {
        return _seq_starts.size_in_bytes();
    }
};
