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

#include <gtl/btree.hpp>
#include <libsais.h>
#include <libsais16.h>
#include <libsais16x64.h>
#include <libsais64.h>
#include <move_r/move_r.hpp>

template <move_r_support support, typename sym_t, typename pos_t>
template <typename sa_sint_t>
void move_r<support, sym_t, pos_t>::construction::build_sa()
{
    std::vector<sa_sint_t>& SA = get_sa<sa_sint_t>(); // [0..n-1] The suffix array

    // Choose the correct suffix array construction algorithm.
    if constexpr (byte_alphabet) {
        log_message(log, "building SA");

        pos_t fs = 6 * 256;
        no_init_resize(SA, n + fs);

        if constexpr (std::is_same_v<sa_sint_t, int32_t>) {
            libsais_omp(&T<uint8_t>(0), SA.data(), n, fs, nullptr, p);
        } else {
            libsais64_omp(&T<uint8_t>(0), SA.data(), n, fs, nullptr, p);
        }
    } else {
        if (idx.symbols_remapped) {
            log_phase_start("mapping T to its effective alphabet");

            #pragma omp parallel for num_threads(p)
            for (uint64_t i = 0; i < n - 1; i++) {
                T<i_sym_t>(i) = (*idx._map_int.find(T<sym_t>(i))).second;
            }

            log_phase_end();
            if (mode == _suffix_array_space)
                store_mapintext();
        }

        log_message(log, "building SA");

        if constexpr (sizeof(i_sym_t) == 2) {
            no_init_resize(SA, n);

            // libsais16's 4th/5th arguments are the free space at the end of SA (0 is enough) and an optional
            // frequency table (NULL) - it has no separate alphabet-size argument, unlike libsais_int below
            if constexpr (sizeof(sa_sint_t) == 4) {
                libsais16_omp(&T<uint16_t>(0), SA.data(), n, 0, nullptr, p);
            } else {
                libsais16x64_omp(&T<uint16_t>(0), SA.data(), n, 0, nullptr, p);
            }
        } else if constexpr (sizeof(i_sym_t) == 4) {
            no_init_resize(SA, n);

            if constexpr (sizeof(sa_sint_t) == 4) {
                libsais_int_omp(&T<int32_t>(0), SA.data(), n, idx.sigma, 0, p);
            } else {
                std::vector<int64_t> T_tmp;
                no_init_resize(T_tmp, n);
            
                #pragma omp parallel for num_threads(p)
                for (pos_t i = 0; i < n; i++) {
                    T_tmp[i] = T<int32_t>(i);
                }

                libsais64_long_omp(T_tmp.data(), SA.data(), n, idx.sigma, 0, p);
            }
        } else if constexpr (sizeof(i_sym_t) == 8) {
            if constexpr (sizeof(sa_sint_t) == 4) {
                no_init_resize(SA, 2 * n);
                libsais64_long_omp(&T<int64_t>(0), (int64_t*) SA.data(), n, idx.sigma, 0, p);

                for (uint64_t i = 0; i < n; i++) {
                    SA[i] = SA[2 * i];
                }
            } else {
                no_init_resize(SA, n);
                libsais64_long_omp(&T<int64_t>(0), SA.data(), n, idx.sigma, 0, p);
            }
        }
    }

    SA.resize(n);

    log_phase_end("time_build_sa");
}