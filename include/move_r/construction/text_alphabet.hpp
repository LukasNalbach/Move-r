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

#ifndef MOVE_R_MOVE_R_CONSTRUCTION_TEXT_ALPHABET_HPP
#define MOVE_R_MOVE_R_CONSTRUCTION_TEXT_ALPHABET_HPP

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>

#include <misc/files.hpp>
#include <misc/log.hpp>
#include <misc/memory.hpp>
#include <misc/utils.hpp>
#include <move_r/construction/pfp_output.hpp>

/**
 * @brief counts the occurrences of each byte in the first length bytes of a file; each thread reads its section
 *        of the file in blocks (see direct_reader)
 * @param file_name name of the file
 * @param length number of bytes to count (from the start of the file)
 * @param p the number of threads to use
 * @param freq_thr [0..p-1][0..255] the counts of each thread, to which the occurrences are added
 * @param prog progress of the phase, advanced by the number of bytes read
 */
inline void count_bytes_of_file(const std::string& file_name, uint64_t length, uint64_t p,
    std::vector<std::vector<uint64_t>>& freq_thr, phase_progress& prog)
{
    static constexpr uint64_t block = huge_page_bytes; // number of bytes read at once (one O_DIRECT read)

    #pragma omp parallel num_threads(p)
    {
        const uint64_t i_p = uint64_t(omp_get_thread_num());
        const uint64_t threads = uint64_t(omp_get_num_threads());
        std::vector<uint64_t>& freq = freq_thr[i_p];
        direct_reader file(file_name);

        // the section [from..to-1] of thread i_p
        const uint64_t from = (length * i_p) / threads;
        const uint64_t to = (length * (i_p + 1)) / threads;

        for (uint64_t pos = from; pos < to; pos += block) {
            const uint64_t num_bytes = std::min(block, to - pos);
            const uint8_t* data = file.read(pos, pos + num_bytes);
            for (uint64_t k = 0; k < num_bytes; k++) freq[data[k]]++;
            prog.advance(num_bytes);
        }
    }
}

/**
 * @brief computes the number of occurrences of each symbol in the BWT of a pfp cache (the sum of
 *        the lengths of its runs), in O(r) time
 * @param cache_file name of the cache file (see move_r::pfp_cache_path)
 * @param num_values exclusive upper bound of the symbols
 * @param p the number of threads to use
 * @param counts [0..num_values-1] stores the number of occurrences of each symbol
 * @param cache_map [0..255] stores the alphabet map of the cache (see pfp_cache_alphabet)
 * @return whether the file is a cache whose symbols are smaller than num_values; if not, counts and cache_map are
 *         unspecified
 */
inline bool read_symbol_counts_from_pfp_cache(const std::string& cache_file, uint64_t num_values, uint64_t p,
    std::vector<uint64_t>& counts, std::array<uint8_t, 256>& cache_map)
{
    pfp_rlbwt::run_data_file cache;
    pfp_cache_alphabet alphabet {};
    if (!open_pfp_cache(cache_file, cache, alphabet) || cache.fields().max_symbol >= num_values) return false;

    const uint64_t r = cache.r();
    const uint64_t threads = std::max<uint64_t>(1, std::min(p, r));
    // each thread adds up its own counts if they take at most r words, else all threads add to counts
    const bool own_counts = threads * num_values <= r;
    std::vector<std::vector<uint64_t>> counts_thr(own_counts ? threads : 0, std::vector<uint64_t>(num_values, 0));
    counts.assign(num_values, 0);
    bool in_range = true;

    #pragma omp parallel for num_threads(threads) schedule(static, 1)
    for (uint64_t t = 0; t < threads; t++) {
        // thread t reads the runs [t * r / threads, (t + 1) * r / threads)
        const uint64_t beg = (t * r) / threads;
        const uint64_t end = ((t + 1) * r) / threads;
        pfp_rlbwt::run_stream runs(cache, beg, end, false, false);
        pfp_rlbwt::run_record rec;

        for (uint64_t x = beg; x < end; x++) {
            runs.next(rec);
            const uint64_t sym = rec.symbol, len = rec.length;

            if (sym >= num_values) [[unlikely]] {
                in_range = false;
                break;
            }

            if (own_counts) counts_thr[t][sym] += len;
            else std::atomic_ref<uint64_t>(counts[sym]).fetch_add(len, std::memory_order_relaxed);
        }
    }

    if (!in_range) return false;

    for (const std::vector<uint64_t>& counts_t : counts_thr) {
        for (uint64_t sym = 0; sym < num_values; sym++) counts[sym] += counts_t[sym];
    }

    cache_map = alphabet.map;
    return true;
}

/**
 * @brief computes the byte frequencies of a text from a pfp cache for the text (see
 *        read_symbol_counts_from_pfp_cache). The caller must compare cache_map with the alphabet map it computes
 *        from freq; if they differ, the cache must not be used.
 * @param cache_file name of the cache file (see move_r::pfp_cache_path)
 * @param p the number of threads to use
 * @param freq [0..255] stores the number of occurrences of each byte
 * @param cache_map [0..255] stores the symbol (>= 3, see pfp_cache_alphabet) of each byte in the cache, 0 for a byte
 *        that does not occur
 * @return whether the file is a cache and the symbols in its BWT match its alphabet map; if not, freq and
 *         cache_map are unspecified
 */
inline bool read_alphabet_from_pfp_cache(const std::string& cache_file, uint64_t p, std::vector<uint64_t>& freq,
    std::array<uint8_t, 256>& cache_map)
{
    constexpr uint64_t num_values = 256;
    std::vector<uint64_t> counts;
    if (!read_symbol_counts_from_pfp_cache(cache_file, num_values, p, counts, cache_map)) return false;

    // byte_of[sym] = the byte that is mapped to the symbol sym (256, if there is none)
    std::array<uint16_t, 256> byte_of;
    byte_of.fill(256);

    for (uint16_t b = 0; b < 256; b++) {
        if (cache_map[b] == 0) continue;
        // the map is invalid, if it maps a byte to a reserved symbol or two bytes to the same symbol
        if (cache_map[b] < 3 || byte_of[cache_map[b]] != 256) return false;
        byte_of[cache_map[b]] = b;
    }

    // the number of occurrences of the symbols 0, 1 and 2, which no byte is mapped to (only the terminator 0 occurs, once)
    uint64_t terminators = 0;
    freq.assign(256, 0);

    for (uint64_t sym = 0; sym < num_values; sym++) {
        const uint64_t len = counts[sym];

        if (sym < 3) {
            terminators += len;
            continue;
        }

        // a symbol must occur in the BWT <=> a byte is mapped to it
        if ((len != 0) != (sym < 256 && byte_of[sym] < 256)) return false;
        if (len != 0) freq[byte_of[sym]] = len;
    }

    return terminators == 1;
}

/* the smallest internal symbol a byte is mapped to before the BWT is built (see move_r::construction::
   min_mapped_symbol): with the in-memory construction modes, the symbol after the terminator 0; with the construction
   mode pfp, the values 0, 1 and 2 are reserved (pfp-rlbwt uses the values below pfp_rlbwt::min_symbol internally), and the
   output of the prefix-free parsing and its cache store the bytes mapped this way */
inline constexpr uint8_t move_r_min_mapped_symbol = 1;
inline constexpr uint8_t move_r_min_mapped_symbol_pfp = 3;
static_assert(move_r_min_mapped_symbol_pfp >= pfp_rlbwt::min_symbol);

/**
 * @brief returns whether the distinct bytes of an input can be mapped to the internal symbols [min_symbol..255] (see
 *        compute_byte_alphabet_map); else the input is mapped like an integer input
 * @param freq_thr [0..p-1][0..255] freq_thr[i_p][c] = number of occurrences of the byte c counted by thread i_p
 * @param min_symbol the smallest symbol a byte is mapped to
 */
inline bool byte_alphabet_fits(const std::vector<std::vector<uint64_t>>& freq_thr, uint64_t min_symbol)
{
    uint64_t num_distinct = 0;

    for (uint16_t c = 0; c < 256; c++) {
        num_distinct += std::any_of(freq_thr.begin(), freq_thr.end(),
            [c](const std::vector<uint64_t>& freq_of_thread) { return freq_of_thread[c] != 0; });
    }

    return num_distinct <= 256 - min_symbol;
}

/**
 * @brief computes the alphabet map of a byte input from the numbers of occurrences of its bytes: the bytes that occur
 *        are mapped to min_symbol, min_symbol + 1, ... in increasing order, or by decreasing frequency (bytes with
 *        equal frequency in increasing order) if reorder is set; throws std::runtime_error if more than
 *        256 - min_symbol distinct bytes occur (see byte_alphabet_fits)
 * @param freq_thr [0..p-1][0..255] freq_thr[i_p][c] = number of occurrences of the byte c counted by thread i_p;
 *        freed
 * @param min_symbol the smallest symbol a byte is mapped to
 * @param reorder whether the bytes are mapped in the order of decreasing frequency (see
 *        move_r_params::reorder_alphabet)
 * @param map_int [0..255] stores the symbol of each byte, or 0 if the byte does not occur
 * @return the number of distinct bytes of the input
 */
inline uint64_t compute_byte_alphabet_map(std::vector<std::vector<uint64_t>>& freq_thr, uint64_t min_symbol,
    bool reorder, std::vector<uint8_t>& map_int)
{
    std::vector<uint64_t> freq(256, 0); // [0..255] the number of occurrences of each byte

    for (uint16_t c = 0; c < 256; c++) {
        for (const std::vector<uint64_t>& freq_of_thread : freq_thr) freq[c] += freq_of_thread[c];
    }

    freq_thr.clear();
    freq_thr.shrink_to_fit();

    // the distinct bytes of the input, in the order in which they are mapped
    std::vector<uint8_t> distinct_bytes;

    for (uint16_t c = 0; c < 256; c++) {
        if (freq[c] != 0) distinct_bytes.emplace_back(uint8_t(c));
    }

    if (distinct_bytes.size() > 256 - min_symbol) {
        throw std::runtime_error("the input contains " + std::to_string(distinct_bytes.size())
            + " distinct symbols, but this construction mode supports at most "
            + std::to_string(256 - min_symbol));
    }

    if (reorder) {
        std::sort(distinct_bytes.begin(), distinct_bytes.end(), [&](uint8_t c1, uint8_t c2) {
            return freq[c1] > freq[c2] || (freq[c1] == freq[c2] && c1 < c2);
        });
    }

    map_int.assign(256, 0);
    uint64_t next_symbol = min_symbol;
    for (uint8_t c : distinct_bytes) map_int[c] = uint8_t(next_symbol++);
    return distinct_bytes.size();
}

/**
 * @brief computes the inverse of an alphabet map of a byte input
 * @param map_int [0..255] the symbol of each byte, or 0 if the byte does not occur
 * @param map_ext [0..255] stores the byte of each symbol, or 0 if no byte is mapped to it
 */
inline void invert_byte_alphabet_map(const std::vector<uint8_t>& map_int, std::vector<uint8_t>& map_ext)
{
    map_ext.assign(256, 0);

    for (uint16_t c = 0; c < 256; c++) {
        if (map_int[c] != 0) map_ext[map_int[c]] = uint8_t(c);
    }
}

#endif // MOVE_R_MOVE_R_CONSTRUCTION_TEXT_ALPHABET_HPP
