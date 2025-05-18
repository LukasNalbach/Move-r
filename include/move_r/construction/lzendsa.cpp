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

#include <move_r/move_r.hpp>

template <move_r_support support, typename sym_t, typename pos_t>
template <bool bigbwt, typename sa_sint_t>
void move_r<support, sym_t, pos_t>::construction::construct_lzendsa()
{
    if (log) std::cout << "building DSA" << std::flush;

    std::vector<sa_sint_t> SA_vec = get_sa<sa_sint_t>();

    if constexpr (bigbwt) {
        no_init_resize(SA_vec, n);

        for (uint16_t i = 0; i < p; i++) {
            SA_file_bufs.emplace_back(sdsl::int_vector_buffer<40>(
                prefix_tmp_files + ".sa", std::ios::in,
                128 * 1024, 40, true));
        }

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < n; i++) {
            SA_vec[i] = SA<true, int32_t>(omp_get_thread_num(), i);
        }
    }
    
    SA_file_bufs.clear();
    SA_file_bufs.shrink_to_fit();

    // make DSA out of SA
    for (uint64_t i = n - 1; i > 0; i--) {
        SA_vec[i] = SA_vec[i] - SA_vec[i - 1];
    }

    if (log) time = log_runtime(time);
    if (log) std::cout << "building LZ-End parsing:" << std::endl;
    std::vector<lzend_phr_t<sa_sint_t>> lzend_phrases = construct_lzend_of_reverse<sa_sint_t>(SA_vec, -1, log);
    idx.z_end = lzend_phrases.size();
    SA_vec.clear();
    SA_vec.shrink_to_fit();
    if (log) time = now();

    if (log) std::cout << "Encoding LZ-End parsing" << std::flush;
    idx._lzendsa = lzendsa_encoding(lzend_phrases, n, 8192);
    if (log) time = log_runtime(time);
}