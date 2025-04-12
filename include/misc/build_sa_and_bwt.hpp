/**
 * part of LukasNalbach/Move-r
 * 
 * MIT License
 * 
 * Copyright (c) Jan Zumbrink, Lukas Nalbach
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

#include <fstream>
#include <filesystem>
#include <libsais.h>
#include <libsais64.h>
#include <vector>
#include <misc/utils.hpp>

template <typename int_t>
std::tuple<std::vector<int_t>, std::string> build_sa_and_bwt(std::string& input, bool build_bwt = true, bool use_bigbwt = false, bool log = false) {
    auto time = now();
    uint64_t n;

    std::vector<int_t> sa;
    std::string bwt;

    if (use_bigbwt) {
        n = std::filesystem::file_size(input) + 1;

        if (log) std::cout << "building Suffix Array and BWT using Big-BWT" << std::flush;
        system(("bigbwt -S " + input + (std::string)(log ? "" : " >log_1 >log_2")).c_str());
        
        if (!log) {
            std::filesystem::remove("log_1");
            std::filesystem::remove("log_2");
        }
        
        std::string sa_file_name = input + ".sa";
        std::string bwt_file_name = input + ".bwt";
        std::ifstream sa_file(sa_file_name);
        bool test = n == std::filesystem::file_size(bwt_file_name);
        sa.resize(n);
        sa[0] = n - 1;

        for (uint64_t i = 1; i < n; i++) {
            sa_file.read((char*) &sa[i], 5);
        }

        sa_file.close();
        std::filesystem::remove(sa_file_name);

        if (build_bwt) {
            no_init_resize(bwt, n);
            std::ifstream bwt_file(bwt_file_name);
            read_from_file(bwt_file, bwt.data(), n);
            bwt_file.close();
            
            for (uint64_t i = 0; i < n; i++) {
                if (bwt[i] == 0) [[unlikely]] {
                    bwt[i] = 1;
                }
            }
        }

        std::filesystem::remove(input + ".log");
        std::filesystem::remove(bwt_file_name);
        if (log) time = log_runtime(time);
    } else {
        input.push_back(uchar_to_char(0));
        n = input.size();

        if (log) std::cout << "building Suffix Array using libsais" << std::flush;
        no_init_resize(sa, n);

        if constexpr (std::is_same_v<int_t, int32_t>) {
            libsais_omp((uint8_t*) input.data(), sa.data(), n, 0, nullptr);
        } else {
            libsais64_omp((uint8_t*) input.data(), sa.data(), n, 0, nullptr);
        }

        if (log) time = log_runtime(time);

        if (build_bwt) {
            if (log) std::cout << "building BWT" << std::flush;
            no_init_resize(bwt, n);

            for (uint64_t i = 0; i < n; i++) {
                if (sa[i] == 0) [[unlikely]] {
                    bwt[i] = input[n - 1];
                } else {
                    bwt[i] = input[sa[i] - 1];
                }

                if (bwt[i] == 0) [[unlikely]] {
                    bwt[i] = 1;
                }
            }

            if (log) time = log_runtime(time);
        }
    }

    return { sa, bwt};
}