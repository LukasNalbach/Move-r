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

#ifndef MOVE_R_MISC_FASTA_HPP
#define MOVE_R_MISC_FASTA_HPP

#include <omp.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

#include <sdsl/int_vector_buffer.hpp>

#include <data_structures/fasta_sequence_data.hpp>
#include <misc/files.hpp>
#include <misc/log.hpp>
#include <misc/memory.hpp>
#include <pfp_rlbwt/fasta.hpp>

// ############################# DNA BASES #############################

// the DNA bases in the order of their 2-bit codes (A = 0, C = 1, G = 2, T = 3)
inline constexpr std::array<char, 4> dna_bases { 'A', 'C', 'G', 'T' };

// [0..255] the 2-bit code of the byte c, if it is a DNA base in either case (see dna_bases), else dna_bases.size()
inline constexpr std::array<uint8_t, 256> dna_base_codes = [] {
    std::array<uint8_t, 256> codes {};
    codes.fill(uint8_t(dna_bases.size()));

    for (uint8_t code = 0; code < dna_bases.size(); code++) {
        codes[uint8_t(dna_bases[code])] = code;
        codes[uint8_t(dna_bases[code] - 'A' + 'a')] = code;
    }

    return codes;
}();

// returns the bases (A, C, G, T) an IUPAC code stands for, or an empty view if c is no IUPAC code
inline std::string_view iupac_bases(uint8_t c)
{
    switch (c) {
        case 'R': return "AG";
        case 'Y': return "CT";
        case 'S': return "CG";
        case 'W': return "AT";
        case 'K': return "GT";
        case 'M': return "AC";
        case 'B': return "CGT";
        case 'D': return "AGT";
        case 'H': return "ACT";
        case 'V': return "ACG";
        case 'N': return "ACGT";
        default:  return {};
    }
}

/**
 * @brief returns the reverse complement of a DNA string: the string reversed, each base complemented (A<->T, C<->G)
 *        and each IUPAC code complemented (R<->Y, K<->M, B<->V, D<->H), case-preserving; other bytes are unchanged
 * @param s the DNA string
 * @return the reverse complement of s
 */
inline std::string reverse_complement(const std::string& s)
{
    static constexpr std::string_view sym = "ACGTRYKMBVDHacgtrykmbvdh";
    static constexpr std::string_view cmp = "TGCAYRMKVBHDtgcayrmkvbhd";
    std::string r(s.size(), 'N');

    for (uint64_t i = 0; i < s.size(); i++) {
        char c = s[s.size() - 1 - i];
        uint64_t j = sym.find(c);
        r[i] = j == std::string_view::npos ? c : cmp[j];
    }

    return r;
}

/** @brief the result of scan_fasta */
struct fasta_scan_t : pfp_rlbwt::fasta_scan {
    fasta_sequence_data data; // the names and the start positions of the sequences
};

/**
 * @brief the start positions and the names of the sequences whose header lines start in a chunk of a FASTA file (see
 *        pfp_rlbwt::scan_fasta_lines); the name of a sequence is its header line up to the first space or tab, without
 *        the '>'. A chunk can start and end inside a header line.
 */
struct fasta_chunk_t : pfp_rlbwt::fasta_chunk {
    /* the start positions of the sequences whose header line starts in the chunk: seq_starts[i] = position in the text
       of the chunk directly after the separator of the i-th such header */
    std::vector<uint64_t> seq_starts;
    // the names of the sequences whose header line starts in the chunk, in the order of the headers
    std::vector<std::string> seq_names;
    bool last_name_open = false; // the last name of seq_names does not end in the chunk
    // the bytes of the name of the header line that the chunk starts inside, as far as they are in the chunk
    std::string continued_name;
    bool continued_name_ends = false; // the name of the header line that the chunk starts inside ends in the chunk

    // the '>' of a header line; it yields a separator (see pfp_rlbwt::fasta_transform::separator_before_first_record)
    void header_start()
    {
        seq_starts.emplace_back(text_size + 1);
        seq_names.emplace_back();
        last_name_open = true;
    }

    // bytes of a header line
    void header_bytes(const uint8_t* bytes, uint64_t len)
    {
        const bool continued = seq_names.empty(); // the bytes belong to the header line the chunk starts inside
        if (continued ? continued_name_ends : !last_name_open) return;
        const uint8_t* end = std::find_if(bytes, bytes + len, [](uint8_t c) { return c == ' ' || c == '\t'; });
        (continued ? continued_name : seq_names.back()).append(reinterpret_cast<const char*>(bytes), end - bytes);
        if (end != bytes + len) end_name(continued);
    }

    // a line feed ends a header line (header = true) or another line
    void line_end(bool header, bool)
    {
        if (header) end_name(seq_names.empty());
    }

    // ends the name of the header line the chunk starts inside (continued = true) or of the last one of seq_names
    void end_name(bool continued)
    {
        if (continued) continued_name_ends = true;
        else last_name_open = false;
    }
};

/**
 * @brief scans a FASTA file in parallel (see pfp_rlbwt::scan_fasta_lines) and appends the start positions and the names
 *        of its sequences (see fasta_chunk_t)
 * @param path the FASTA file
 * @param transform the transformation; each header line has to yield a separator
 * @param prm the parameters of the scan
 * @param text_start the position in the text of the first symbol of the text of the file
 * @param seq_starts (appended) the positions in the text directly after the separators of the headers
 * @param seq_names (appended) the names of the sequences
 * @return the text of the file
 */
inline pfp_rlbwt::fasta_scan scan_fasta_file(const std::string& path, const pfp_rlbwt::fasta_transform& transform,
    const pfp_rlbwt::fasta_scan_params& prm, uint64_t text_start, std::vector<uint64_t>& seq_starts,
    std::vector<std::string>& seq_names)
{
    bool name_open = false; // the last name of seq_names continues in the next chunk

    return pfp_rlbwt::scan_fasta_lines<fasta_chunk_t>([&] { return direct_reader(path); }, transform, prm,
        [&](fasta_chunk_t& chunk, uint64_t shift) {
            if (name_open) {
                seq_names.back() += chunk.continued_name;
                name_open = !chunk.continued_name_ends;
            }

            for (uint64_t start : chunk.seq_starts) seq_starts.emplace_back(text_start + shift + start);
            for (std::string& name : chunk.seq_names) seq_names.emplace_back(std::move(name));
            if (!chunk.seq_names.empty()) name_open = chunk.last_name_open;
        });
}

/**
 * @brief writes the text of a FASTA file to a buffer in parallel, one range of the text per thread (see
 *        pfp_rlbwt::fasta_source)
 * @param path the FASTA file
 * @param transform the transformation
 * @param scan the result of a scan of the file with samples (see pfp_rlbwt::scan_fasta_lines)
 * @param dst [0..scan.text_size-1] the buffer
 * @param num_threads the maximum number of threads to use
 */
inline void read_fasta_text(const std::string& path, const pfp_rlbwt::fasta_transform& transform,
    const pfp_rlbwt::fasta_scan& scan, uint8_t* dst, uint64_t num_threads)
{
    const uint64_t n = scan.text_size;
    if (n == 0) return;
    const pfp_rlbwt::fasta_source<uint8_t> text(path, n, transform, scan.samples.data(), scan.samples.size());
    // each range has at least pfp_rlbwt::fasta_sample_period symbols
    const uint64_t p = std::max<uint64_t>(1, std::min<uint64_t>(num_threads, n / pfp_rlbwt::fasta_sample_period));

    #pragma omp parallel for num_threads(p) schedule(static, 1)
    for (uint64_t t = 0; t < p; t++) {
        const uint64_t beg = t * n / p;
        const uint64_t end = (t + 1) * n / p;
        text.new_reader()->read(beg, end - beg, dst + beg);
    }
}

/**
 * @brief computes the length, the symbol frequencies, the sequences and the samples of the text that FASTA files
 *        describe (see pfp_rlbwt::fasta_transform), without writing the text
 * @param in_paths the FASTA files; they are concatenated into one text as in process_fasta
 * @param alphabet the bytes that are kept verbatim; a run of other bytes becomes one @p separator
 * @param separator the byte a header line and a masked run each yield
 * @param to_upper whether to upper-case a body byte before looking it up in @p alphabet
 * @param num_threads the maximum number of threads to use
 * @param progress if not nullptr, counts the bytes scanned
 * @return the result of the scan (see fasta_scan_t); the samples are only computed for a single file
 */
inline fasta_scan_t scan_fasta(const std::vector<std::string>& in_paths,
    const std::string& alphabet = "ACGT", char separator = fasta_separator_symbol,
    bool to_upper = true, uint64_t num_threads = omp_get_max_threads(), phase_progress* progress = nullptr)
{
    const pfp_rlbwt::fasta_transform transform(alphabet, (uint8_t) separator, to_upper);
    fasta_scan_t result;
    // [0..num_seq-1] seq_starts[i] = position in the text directly after the separator of the header of sequence i
    std::vector<uint64_t> seq_starts;
    std::vector<std::string> seq_names; // [0..num_seq-1] seq_names[i] = name of sequence i

    // the samples store offsets in a file, so they are only computed for a single file
    const pfp_rlbwt::fasta_scan_params prm {
        .num_threads = num_threads,
        .sample_period = in_paths.size() == 1 ? pfp_rlbwt::fasta_sample_period : 0,
        .progress = progress };

    for (const std::string& in_path : in_paths) {
        pfp_rlbwt::fasta_scan file = scan_fasta_file(in_path, transform, prm, result.text_size, seq_starts, seq_names);
        for (uint64_t c = 0; c < result.freq.size(); c++) result.freq[c] += file.freq[c];
        result.text_size += file.text_size;
        if (prm.sample_period != 0) result.samples = std::move(file.samples);
    }

    // the index appends one sentinel to the text, so its length is text_size + 1
    result.data.set_sequences(result.text_size + 1, seq_starts, std::move(seq_names), separator);
    return result;
}

/**
 * @brief writes the text that FASTA files describe (see pfp_rlbwt::fasta_transform) to a buffer; the result is the
 *        text process_fasta writes to a file. Each file is scanned in parallel (see pfp_rlbwt::scan_fasta_lines), then
 *        its text is written in parallel (see read_fasta_text).
 * @param in_paths the FASTA files; they are concatenated into one text
 * @param dst the buffer; its size has to be at least the text size computed by scan_fasta for the same files
 * @param alphabet the bytes that are kept verbatim; a run of other bytes becomes one @p separator
 * @param separator the byte a header line and a masked run each yield
 * @param to_upper whether to upper-case a body byte before looking it up in @p alphabet
 * @param num_threads the maximum number of threads to use
 * @return the number of symbols written
 */
inline uint64_t read_fasta_into(const std::vector<std::string>& in_paths, uint8_t* dst,
    const std::string& alphabet = "ACGT", char separator = fasta_separator_symbol,
    bool to_upper = true, uint64_t num_threads = omp_get_max_threads())
{
    const pfp_rlbwt::fasta_transform transform(alphabet, (uint8_t) separator, to_upper);
    uint64_t written = 0;

    for (const std::string& in_path : in_paths) {
        const pfp_rlbwt::fasta_scan file = pfp_rlbwt::scan_fasta_lines<pfp_rlbwt::fasta_chunk>(
            [&] { return direct_reader(in_path); }, transform, { .num_threads = num_threads }, [](auto&&...) { });
        read_fasta_text(in_path, transform, file, dst + written, num_threads);
        written += file.text_size;
    }

    return written;
}

/**
 * @brief writes the text that FASTA files describe (see pfp_rlbwt::fasta_transform) to a file and returns its
 *        sequences; each file is scanned in parallel (see pfp_rlbwt::scan_fasta_lines), then its text is read in blocks
 *        of direct_reader::block symbols (see pfp_rlbwt::fasta_source)
 * @param in_paths the FASTA files; they are concatenated into one text, the sequences of a later file follow those
 *        of the earlier ones
 * @param out_path the file the text is written to
 * @param alphabet the bytes that are kept verbatim; a run of other bytes becomes one @p separator
 * @param separator the byte a header line and a masked run each yield; it must not be in @p alphabet
 * @param to_upper whether to upper-case a body byte before looking it up in @p alphabet
 * @param num_threads the maximum number of threads of the scans
 * @return the names and the start positions of the sequences
 */
inline fasta_sequence_data process_fasta(const std::vector<std::string>& in_paths,
    const std::string& out_path, const std::string& alphabet = "ACGT",
    char separator = fasta_separator_symbol, bool to_upper = true, uint64_t num_threads = omp_get_max_threads())
{
    if (alphabet.find(separator) != std::string::npos)
        throw std::runtime_error("the FASTA separator must not be in the alphabet");

    const pfp_rlbwt::fasta_transform transform(alphabet, (uint8_t) separator, to_upper);
    direct_ofstream out(out_path);
    // [0..num_seq-1] seq_starts[i] = position in the text directly after the separator of the header of sequence i
    std::vector<uint64_t> seq_starts;
    std::vector<std::string> seq_names; // [0..num_seq-1] seq_names[i] = name of sequence i
    uint64_t text_size = 0; // number of bytes written to out_path
    no_init_vector<uint8_t> block; // a block of the text of a file

    for (const std::string& in_path : in_paths) {
        const pfp_rlbwt::fasta_scan file = scan_fasta_file(in_path, transform, { .num_threads = num_threads },
            text_size, seq_starts, seq_names);

        if (file.text_size != 0) {
            const pfp_rlbwt::fasta_source<uint8_t> text(in_path, file.text_size, transform, file.samples.data(),
                file.samples.size());
            const auto reader = text.new_reader();
            resize_no_init(block, std::min<uint64_t>(file.text_size, direct_reader::block));

            for (uint64_t pos = 0; pos < file.text_size;) {
                const uint64_t len = std::min<uint64_t>(block.size(), file.text_size - pos);
                reader->read(pos, len, block.data());
                out.write(reinterpret_cast<const char*>(block.data()), len);
                pos += len;
            }
        }

        text_size += file.text_size;
    }

    out.close();

    // the index appends one sentinel to the text, so its length is text_size + 1
    fasta_sequence_data seq_data;
    seq_data.set_sequences(text_size + 1, seq_starts, std::move(seq_names), separator);
    return seq_data;
}

// writes the text that a FASTA file describes to a file and returns its sequences (see process_fasta)
inline fasta_sequence_data process_fasta(const std::string& in_path, const std::string& out_path,
    const std::string& alphabet = "ACGT", char separator = fasta_separator_symbol, bool to_upper = true,
    uint64_t num_threads = omp_get_max_threads())
{
    return process_fasta(std::vector<std::string>{in_path}, out_path, alphabet, separator, to_upper, num_threads);
}

#endif // MOVE_R_MISC_FASTA_HPP
