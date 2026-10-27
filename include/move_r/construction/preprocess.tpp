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

#ifndef MOVE_R_MOVE_R_CONSTRUCTION_PREPROCESS_TPP
#define MOVE_R_MOVE_R_CONSTRUCTION_PREPROCESS_TPP

#include <misc/sort.hpp>
#include <misc/files.hpp>
#include <misc/utils.hpp>
#include <move_r/construction/text_alphabet.hpp>
#include <move_r/construction/construction.hpp>

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::map_byte_alphabet(std::vector<std::vector<uint64_t>>& freq_thr)
{
    // sigma includes the terminator
    idx.sigma = 1 + compute_byte_alphabet_map(freq_thr, min_mapped_symbol(), params.reorder_alphabet, idx._map_int);
    idx.symbols_remapped = true;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::finish_byte_alphabet()
{
    if (move_r_builds_runs_file(params.mode)) {
        // the symbols in the RLBWT are the mapped bytes decreased by 2 (see build_rlbwt_c)
        for (uint16_t i = 0; i < 256; i++) {
            if (idx._map_int[i] != 0) idx._map_int[i] -= 2;
        }
    }

    invert_byte_alphabet_map(idx._map_int, idx._map_ext);
    p_ = p;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::widen_t()
{
    move_r_widen_bytes(T_inp, T_wide, p);

    if (delete_T) {
        T_inp.clear();
        T_inp.shrink_to_fit();
    }

    T_inp = move_r_input_view(T_wide);
    inp_width = sizeof(uint16_t);
    T_widened = true;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::preprocess_t(const std::string& T_file_name)
{
    log_phase_start("preprocessing T");

    /* a byte alphabet is mapped to [min_mapped_symbol()..min_mapped_symbol()+sigma-2] (by decreasing frequency, if
       reorder_alphabet is set); an integer alphabet, and a byte input with more distinct bytes than that range holds,
       is mapped to [1..sigma-1] */
    if (byte_alphabet) {
        // [0..p-1][0..255] freq_thr[i_p][c] = number of occurrences of c that thread i_p has counted in T[0..n-2]
        std::vector<std::vector<uint64_t>> freq_thr(p, std::vector<uint64_t>(256, 0));
        bool alphabet_from_cache = false; // true <=> the frequencies have been read from the pfp cache
        std::vector<uint64_t> cache_freq; // [0..255] the frequencies read from the cache
        std::array<uint8_t, 256> cache_map {}; // [0..255] the alphabet map of the cache (see pfp_cache_alphabet)

        if (params.fasta.byte_freqs != nullptr) {
            // the frequencies of the symbols of the text of a FASTA file have been counted by scan_fasta
            for (uint16_t c = 0; c < 256; c++) freq_thr[0][c] = params.fasta.byte_freqs[c];
        } else if (input_in_memory()) {
            #pragma omp parallel for num_threads(p)
            for (uint64_t i = 0; i < n - 1; i++) {
                freq_thr[omp_get_thread_num()][T<uint8_t>(i)]++;
            }
        } else if (move_r_builds_runs_file(params.mode) && uses_pfp_cache() && !params.ignore_pfp_cache &&
                   std::filesystem::exists(pfp_cache_file()) &&
                   read_alphabet_from_pfp_cache(pfp_cache_file(), p, cache_freq, cache_map)) {
            // the frequencies have been computed from the run lengths of the cached BWT (see
            // read_alphabet_from_pfp_cache)
            freq_thr[0] = cache_freq;
            alphabet_from_cache = true;
        } else if (params.fasta.text_len != 0) {
            throw std::invalid_argument("move_r: the symbol frequencies of the FASTA input are missing");
        } else {
            phase_progress prog(params.log, n - 1);
            count_bytes_of_file(T_file_name, n - 1, p, freq_thr, prog);
            prog.finish();
        }

        // the text of a FASTA file is read only through the map of the bytes (see build_rlbwt_pfp)
        if (byte_alphabet_fits(freq_thr, min_mapped_symbol()) || params.fasta.text_len != 0) {
            map_byte_alphabet(freq_thr);

            // a cache that has been written with a different alphabet map is not used
            if (alphabet_from_cache && !std::equal(cache_map.begin(), cache_map.end(), idx._map_int.begin())) {
                params.ignore_pfp_cache = true;

                log_note(params.log, "the cache " + pfp_prefix + " was written with a different alphabet map, "
                    "rebuilding it");
            }

            if (input_in_memory()) {
                #pragma omp parallel for num_threads(p)
                for (uint64_t i = 0; i < n - 1; i++) {
                    T<uint8_t>(i) = idx._map_int[T<uint8_t>(i)];
                }
            } else {
                // the input file is not modified; the alphabet map is applied while the prefix-free parsing reads it
                params.pfp_alphabet_map_file = prefix_tmp_files + ".map";
                direct_ofstream map_ofile(params.pfp_alphabet_map_file, std::ios::binary);
                write_in_chunks(map_ofile, (char*) idx._map_int.data(), 256);
                map_ofile.close();
            }

            finish_byte_alphabet();
        } else {
            // a byte input with more distinct bytes than internal symbols of one byte is mapped like an integer input
            byte_alphabet = false;
            idx.byte_alphabet = false;
            // the internal symbols [1..256] of an input in memory with all 256 bytes do not fit into its symbols
            if (input_in_memory() && !byte_alphabet_fits(freq_thr, move_r_min_mapped_symbol)) widen_t();
        }
    }

    if (!byte_alphabet) {
        idx.symbols_remapped = true;
        uint64_t alloc_before = malloc_count_current();

        int_alphabet = move_r_int_alphabet::make(inp_width);

        // the prefix-free parsing reads a file input itself, so it is not read into memory (see build_index)
        if (move_r_builds_runs_file(params.mode) && file_input()) int_alphabet->collect_file(T_file_name, n, p);
        else int_alphabet->collect(T_inp, n);

        idx.sigma = int_alphabet->num_symbols() + 1;
        int_alphabet_bytes = malloc_count_current() - alloc_before;

        /* the distinct symbols in increasing order are stored in _alphabet (as a bit vector); the internal symbols
           are ordered by decreasing frequency, the permutation is stored in _alphabet_order (see
           move_r_order_alphabet) */
        std::vector<uint64_t> symbols;
        int_alphabet->sorted_symbols(symbols, p);

        if (std::binary_search(symbols.begin(), symbols.end(), TERMINATOR_SYMBOL))
            throw std::invalid_argument("move_r: the input contains the symbol TERMINATOR_SYMBOL (bit pattern 2^63)");

        std::vector<uint64_t> freq;
        int_alphabet->frequencies(symbols, freq, p);
        std::vector<uint64_t> effective;
        move_r_order_alphabet(freq, params.reorder_alphabet, effective, p);

        // the cache is used only if its symbol counts equal the counts of the effective symbols (offset min_mapped_symbol() - 1)
        if (move_r_builds_runs_file(params.mode) && uses_pfp_cache() && !params.ignore_pfp_cache &&
            std::filesystem::exists(pfp_cache_file())) {
            const uint64_t offset = min_mapped_symbol() - 1;
            std::vector<uint64_t> expected(idx.sigma + offset, 0);
            expected[0] = 1;
            for (uint64_t i = 0; i < freq.size(); i++) expected[(effective.empty() ? i + 1 : effective[i]) + offset] = freq[i];
            std::vector<uint64_t> counts;
            std::array<uint8_t, 256> cache_map;

            if (!read_symbol_counts_from_pfp_cache(pfp_cache_file(), expected.size(), p, counts, cache_map) ||
                counts != expected) {
                params.ignore_pfp_cache = true;

                log_note(params.log, "the cache " + pfp_prefix + " was written with a different alphabet order, "
                    "rebuilding it");
            }
        }

        freq.clear();
        freq.shrink_to_fit();

        int_alphabet->assign_effective_symbols(symbols, effective, p);

        idx.assign_alphabet(symbols, effective);
    }

    log_phase_end("time_preprocess_t");

    if (!byte_alphabet) {
        if (params.log) std::cout << "alphabet map: " << format_size(int_alphabet_bytes) << std::endl;
    }
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::read_t_from_file(
    const std::string& T_file_name)
{
    log_phase_start("reading T");

    T_inp.resize(n);
    T_inp.set(n - 1, 0);
    phase_progress prog(params.log, (n - 1) * inp_width);

    // the pages of the file are dropped from the page cache while it is read (see read_file_in_blocks)
    if (!read_file_in_blocks(T_file_name, (n - 1) * inp_width, T_inp.data(),
            ignore_block, 64 * 1024 * 1024, &prog)) {
        throw std::runtime_error("move_r: cannot read the input file " + T_file_name);
    }

    prog.finish();
    log_phase_end("time_read_t");
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::unmap_t()
{
    // a widened copy of T has been mapped instead of the input
    if (T_widened) return;
    log_phase_start("unmapping T");

    if (byte_alphabet) {
        const uint8_t offset = min_mapped_symbol() - 1;

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < n - 1; i++) {
            T<uint8_t>(i) = idx._map_ext[T<uint8_t>(i) - offset];
        }
    } else {
        // maps each internal symbol to its input symbol with _alphabet_order_inv and a select query on _alphabet
        move_r_unmap_text(T_inp, n, p, idx._alphabet, idx._alphabet_order_inv);
    }

    log_phase_end("time_unmap_t");
}

#endif // MOVE_R_MOVE_R_CONSTRUCTION_PREPROCESS_TPP
