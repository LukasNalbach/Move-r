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

#include <misc/files.hpp>
#include <move_r/move_r.hpp>

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
void move_r<support, sym_t, pos_t, mlf_enc>::construction::preprocess_t(bool in_memory, bool bigbwt, const std::string& T_file_name)
{
    log_message(log, "preprocessing T");

    if constexpr (byte_alphabet) {
        // contains_uchar_thr[i_p][c] = true <=> thread i_p found the character c in its section of T[0..n-2].
        std::vector<std::vector<uint8_t>> contains_uchar_thr(p, std::vector<uint8_t>(256, 0));

        if (in_memory) {
            // Iterate over T[0..n-2] and report the occurrence of each found character in T[0..n-2] in contains_uchar_thr.
            #pragma omp parallel for num_threads(p)
            for (uint64_t i = 0; i < n - 1; i++) {
                contains_uchar_thr[omp_get_thread_num()][T<i_sym_t>(i)] = 1;
            }
        } else {
            std::vector<sdsl::int_vector_buffer<8>> T_buf;

            for (uint16_t i = 0; i < p; i++) {
                T_buf.emplace_back(sdsl::int_vector_buffer<8>(T_file_name, std::ios::in, 128 * 1024, 8, true));
            }

            #pragma omp parallel for num_threads(p)
            for (uint64_t i = 0; i < n - 1; i++) {
                contains_uchar_thr[omp_get_thread_num()][T_buf[omp_get_thread_num()][i]] = 1;
            }
        }

        // contains_uchar[c] = 1 <=> c in T[0..n-2].
        std::vector<uint8_t> contains_uchar(256, 0);

        /* Combine the results of each thread's sections in contains_uchar_thr[0..i_p-1][0..255]
        into contains_uchar[0..255]. */
        for (uint16_t cur_uchar = 0; cur_uchar < 256; cur_uchar++) {
            for (uint16_t i_p = 0; i_p < p; i_p++) {
                if (contains_uchar_thr[i_p][cur_uchar] == 1) {
                    contains_uchar[cur_uchar] = 1;
                    break;
                }
            }
        }

        contains_uchar_thr.clear();
        contains_uchar_thr.shrink_to_fit();

        // The number of distinct characters in T[0..n-1].
        idx.sigma = 1;

        // Count the number of ones in contains_uchar[0..255] in idx._map_ext.
        for (uint16_t cur_uchar = 0; cur_uchar < 256; cur_uchar++) {
            if (contains_uchar[cur_uchar] == 1) {
                idx.sigma++;
            }
        }

        /* If T[0..n-2] contains more than 256 - min_valid_char distinct characters, we cannot remap them into the
            range [0..255] without using a character less than min_valid_char, hence we cannot build an index for T. */
        pos_t max_num_chars = 256 - min_valid_char;

        if (idx.sigma - 1 > max_num_chars) {
            std::cout << "Error: the input contains more than " << std::to_string(max_num_chars) << " distinct characters" << std::endl;
            return;
        }

        idx.symbols_remapped = true;

        // build the mapping function map_symbol that remaps the characters of T, s.t. it does
        // not contain 0 or 1; also build its inverse function unmap_symbol

        idx._map_int.resize(256, 0);
        idx._map_ext.resize(256, 0);

        /* To preserve the order among characters in T[0..n-2], we start by mapping smallest
        character in T[0..n-2] to 2, the second smallest to 3, ... . */

        // The character, to map the currently next largest character in T[0..n-2] to.
        uint16_t next_uchar_to_remap_to = min_valid_char;

        for (uint16_t cur_uchar = 0; cur_uchar < 256; cur_uchar++) {
            if (contains_uchar[cur_uchar] == 1) {
                idx._map_int[cur_uchar] = next_uchar_to_remap_to;
                idx._map_ext[next_uchar_to_remap_to] = cur_uchar;
                next_uchar_to_remap_to++;
            }
        }

        // Apply map_symbol to T.
        if (in_memory) {
            #pragma omp parallel for num_threads(p)
            for (uint64_t i = 0; i < n - 1; i++) {
                T<i_sym_t>(i) = idx._map_int[T<i_sym_t>(i)];
            }
        } else {
            transform_byte_file_parallel(T_file_name, n - 1, p,
                [&](uint8_t c) { return idx._map_int[c]; });
        }

        if (bigbwt) {
            for (uint16_t i = 0; i < 256; i++) {
                if (idx._map_int[i] != 0) {
                    idx._map_int[i] -= 2;
                }
            }

            for (uint16_t i = 0; i < 256; i++) {
                idx._map_ext[i] = 0;
            }

            for (uint16_t i = 0; i < 256; i++) {
                if (idx._map_int[i] != 0) {
                    idx._map_ext[idx._map_int[i]] = i;
                }
            }
        }

        p_ = p;
    } else {
        idx.symbols_remapped = true;
        uint64_t alloc_before = malloc_count_current();

        for (pos_t i = 0; i < n - 1; i++) {
            idx._map_int.try_emplace(T<sym_t>(i), 0);
        }

        idx.sigma = idx._map_int.size() + 1;
        idx.size_map_int = malloc_count_current() - alloc_before;
        no_init_resize(idx._map_ext, idx.sigma);
        idx._map_ext[0] = 0;
        pos_t sym_cur = 1;

        for (const auto& kv : idx._map_int) {
            idx._map_ext[sym_cur] = kv.first;
            sym_cur++;
        }

        if (p == 1) {
            ips4o::sort(idx._map_ext.begin() + 1, idx._map_ext.end());
        } else {
            ips4o::parallel::sort(idx._map_ext.begin() + 1, idx._map_ext.end());
        }

        if (idx.sigma > p) {
            #pragma omp parallel num_threads(p)
            {
                uint16_t i_p = omp_get_thread_num();

                i_sym_t b = 1 + i_p * ((idx.sigma - 1) / p);
                i_sym_t e = i_p == p - 1 ? idx.sigma - 1 : ((i_p + 1) * ((idx.sigma - 1) / p));

                for (i_sym_t i = b; i <= e; i++) {
                    idx._map_int[idx._map_ext[i]] = i;
                }
            }
        } else {
            for (i_sym_t i = 1; i < idx.sigma; i++) {
                idx._map_int[idx._map_ext[i]] = i;
            }
        }
    }

    log_phase_end("time_preprocess_t");
}

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
void move_r<support, sym_t, pos_t, mlf_enc>::construction::read_t_from_file(std::string& T_file_name)
{
    log_phase_start("reading T");

    no_init_resize(T_str, n);
    T<i_sym_t>(n - 1) = 0;
    std::ifstream T_file(T_file_name);
    read_from_file(T_file, T_str.data(), n - 1);

    log_phase_end();
}

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
void move_r<support, sym_t, pos_t, mlf_enc>::construction::unmap_t(bool in_memory)
{
    if constexpr (str_input) {
        if (in_memory) {
            #pragma omp parallel for num_threads(p)
            for (uint64_t i = 0; i < n - 1; i++) {
                T<i_sym_t>(i) = idx._map_ext[T<i_sym_t>(i)];
            }
        } else {
            std::vector<uint64_t> map_ext_adj(256, 0);

            for (uint16_t c = 2; c < 256; c++) {
                map_ext_adj[c] = char_to_uchar(idx._map_ext[c - 2]);
            }

            transform_byte_file_parallel(prefix_tmp_files, n - 1, p,
                [&](uint8_t c) { return map_ext_adj[c]; });
        }
    } else {
        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < n - 1; i++) {
            T<sym_t>(i) = idx._map_ext[T<i_sym_t>(i)];
        }
    }
}