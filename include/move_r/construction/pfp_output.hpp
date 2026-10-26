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

#ifndef MOVE_R_MOVE_R_CONSTRUCTION_PFP_OUTPUT_HPP
#define MOVE_R_MOVE_R_CONSTRUCTION_PFP_OUTPUT_HPP

#include <algorithm>
#include <array>
#include <cstring>
#include <filesystem>
#include <limits>
#include <string>

#include <data_structures/fasta_sequence_data.hpp>
#include <misc/files.hpp>
#include <rlbwt/readers.hpp>
#include <pfp_rlbwt/pfp_rlbwt.hpp>

/* a pfp cache is the runs file of pfp_rlbwt::build (see pfp_rlbwt::ext_runs) with two further sections: the
   alphabet map (see pfp_cache_alphabet) and, if the text has been read from a FASTA file, the length of the text and
   its sequences (see write_fasta_to_pfp_cache) */
inline constexpr uint64_t pfp_cache_alphabet_section = pfp_rlbwt::first_extra_section;
inline constexpr uint64_t pfp_cache_fasta_section = pfp_rlbwt::first_extra_section + 1;

/**
 * @brief the alphabet map of a cache (the section pfp_cache_alphabet_section)
 */
struct pfp_cache_alphabet {
    // [0..255] the symbol (>= 3, see construction::min_mapped_symbol) of each byte in the cache, 0 if the byte does
    // not occur
    std::array<uint8_t, 256> map {};
};

static_assert(sizeof(pfp_cache_alphabet) == sizeof(pfp_cache_alphabet::map));

/**
 * @brief reads the table of contents and the alphabet map of a cache
 * @param cache_file name of the cache file
 * @param cache stores the table of contents
 * @param alphabet stores the alphabet map
 * @return whether the file is a cache, i.e., it is a runs file (see pfp_rlbwt::run_data_file::open) and stores the
 *         alphabet map
 */
inline bool open_pfp_cache(const std::string& cache_file, pfp_rlbwt::run_data_file& cache, pfp_cache_alphabet& alphabet)
{
    std::string bytes;
    if (!cache.open(cache_file) || !cache.read_section(pfp_cache_alphabet_section, bytes)
        || bytes.size() != sizeof(pfp_cache_alphabet))
        return false;
    std::memcpy(&alphabet, bytes.data(), sizeof(pfp_cache_alphabet));
    return true;
}

/**
 * @brief reads the length of the text and its sequences from a pfp cache for a FASTA file
 * @param cache_file name of the cache file
 * @param text_size stores the length of the text
 * @param sequences stores the sequences
 * @return whether the file is a cache that stores them (see pfp_cache_fasta_section)
 */
inline bool read_fasta_from_pfp_cache(const std::string& cache_file, uint64_t& text_size,
    fasta_sequence_data& sequences)
{
    pfp_rlbwt::run_data_file cache;
    pfp_cache_alphabet alphabet {};
    if (!open_pfp_cache(cache_file, cache, alphabet) || !cache.contains(pfp_cache_fasta_section)) return false;
    direct_ifstream in(cache_file, std::ios::binary);
    if (!in.good()) return false;
    in.seekg(std::streamoff(cache.section(pfp_cache_fasta_section).offset));
    in.read((char*) &text_size, sizeof(uint64_t));
    sequences.load(in);
    return in.good();
}

/**
 * @brief writes the length of the text and its sequences in the layout of the section pfp_cache_fasta_section
 * @param out output stream
 * @param text_size length of the text
 * @param sequences the sequences
 */
inline void write_fasta_to_pfp_cache(std::ostream& out, uint64_t text_size, const fasta_sequence_data& sequences)
{
    out.write((const char*) &text_size, sizeof(uint64_t));
    sequences.serialize(out);
}

#endif // MOVE_R_MOVE_R_CONSTRUCTION_PFP_OUTPUT_HPP
