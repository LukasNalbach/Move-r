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
#include <cstdint>
#include <string>
#include <vector>

#include <sdsl/int_vector_buffer.hpp>

/** @brief the default reserved byte masking a non-allowed (e.g. non-ACGT / ambiguous) base */
inline constexpr char fasta_mask_symbol = 'N';
/** @brief the default reserved byte replacing each sequence's header line (the sequence separator). It must be the
 *         largest byte in the produced text so that it maps to the largest internal symbol, which the approximate
 *         search recognizes as the forbidden boundary and never extends into. 0x7F exceeds every printable base. */
inline constexpr char fasta_separator_symbol = (char) 0x7F;

/**
 * @brief the per-sequence metadata of a FASTA file (small -- one entry per sequence, not per base): for each
 *        sequence, the start position of its body in the produced text and its name (RNAME for SAM output)
 */
struct fasta_index_t {
    std::vector<uint64_t> seq_starts;    // seq_starts[i] = start position of sequence i's body in the text
    std::vector<std::string> seq_names;  // seq_names[i] = name of sequence i (header up to first whitespace)
};

/**
 * @brief preprocesses a FASTA file on disk, writing its preprocessed text to @p out_path and returning only the
 *        (small) per-sequence metadata. Both files are streamed through sdsl::int_vector_buffer, so the text is
 *        never held in memory and arbitrarily large genomes can be preprocessed within a compressed index's memory
 *        budget. The produced text is a plain byte file, read for construction like any other input.
 *
 *        Each record's header line (beginning with '>') becomes one reserved @p separator byte and gives the
 *        sequence name (the header up to the first whitespace). Within a body, every byte not in @p allowed is
 *        replaced with the reserved @p mask byte (after upper-casing if @p to_upper). Both reserved bytes lie
 *        outside @p allowed, so a query over the allowed alphabet cannot match a masked base or span a boundary.
 * @param in_path the path to the FASTA file
 * @param out_path the path the preprocessed text (a plain byte file) is written to
 * @param allowed the alphabet kept verbatim (default "ACGT"); any other byte in a body becomes @p mask
 * @param mask the reserved byte masking a non-allowed base (must not be in @p allowed)
 * @param separator the reserved byte replacing each header line / separating sequences (must not be in @p allowed)
 * @param to_upper whether to upper-case body bytes before testing membership in @p allowed
 * @param buffer_size the per-file buffer size in bytes
 * @return the per-sequence body start positions and names
 */
inline fasta_index_t process_fasta(const std::string& in_path, const std::string& out_path,
    const std::string& allowed = "ACGT", char mask = fasta_mask_symbol, char separator = fasta_separator_symbol,
    bool to_upper = true, uint64_t buffer_size = 8 * 1024 * 1024)
{
    // the separator must be the largest byte in the produced text: it then maps to the largest internal symbol, which
    // the search treats as the forbidden sequence boundary. So it must lie above the mask and every allowed byte.
    for (unsigned char c : allowed)
        if ((unsigned char) separator <= c)
            throw std::runtime_error("FASTA separator must be greater than every allowed byte");
    if ((unsigned char) separator <= (unsigned char) mask)
        throw std::runtime_error("FASTA separator must be greater than the mask byte");

    // both files are accessed through int_vector_buffer as plain (header-less) byte files, read/written sequentially
    sdsl::int_vector_buffer<8> in(in_path, std::ios::in, buffer_size, 8, true);
    sdsl::int_vector_buffer<8> out(out_path, std::ios::out, buffer_size, 8, true);

    std::array<uint8_t, 256> keep{};
    for (unsigned char c : allowed) keep[c] = 1;

    fasta_index_t index;
    uint64_t text_size = 0;     // number of bytes written to out so far
    bool in_sequence = false;   // whether any header has been seen yet
    bool at_line_start = true;  // whether the current byte begins a line
    bool header_line = false;   // whether the current line is a header ('>') line
    bool name_done = false;     // whether the current header's name (up to first whitespace) is complete
    std::string name;

    uint64_t size = in.size();

    for (uint64_t i = 0; i < size; i++) {
        char c = (char) in[i];

        if (c == '\r') continue; // tolerate CRLF

        if (c == '\n') {
            if (header_line) index.seq_names.emplace_back(std::move(name)), name.clear();
            at_line_start = true;
            header_line = false;
            name_done = false;
            continue;
        }

        if (at_line_start) {
            at_line_start = false;
            header_line = (c == '>');

            if (header_line) {
                // a header (and the newlines around it) becomes one separator symbol that precedes the body
                out.push_back((uint8_t) separator);
                text_size++;
                index.seq_starts.emplace_back(text_size);
                in_sequence = true;
                continue; // the '>' itself is not part of the name
            }
        }

        if (header_line) {
            // the name is the header up to the first whitespace
            if (!name_done) {
                if (c == ' ' || c == '\t') name_done = true;
                else name.push_back(c);
            }
            continue;
        }

        if (!in_sequence) continue; // sequence data before any header: skip

        char u = (to_upper && c >= 'a' && c <= 'z') ? char(c - 'a' + 'A') : c;
        out.push_back((uint8_t)(keep[(unsigned char) u] ? u : mask));
        text_size++;
    }

    if (header_line) index.seq_names.emplace_back(std::move(name)); // a final header line with no trailing newline

    return index;
}
