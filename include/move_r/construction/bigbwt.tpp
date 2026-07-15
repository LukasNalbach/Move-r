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

#include <filesystem>
#include <stdexcept>
#include <system_error>
#include <misc/bigbwt.hpp>
#include <move_r/move_r.hpp>

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
void move_r<support, sym_t, pos_t, mlf_enc>::construction::store_t_in_file()
{
    std::ofstream ofile(prefix_tmp_files);
    write_to_file(ofile, T_str.c_str(), n - 1);
    ofile.close();
}

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
void move_r<support, sym_t, pos_t, mlf_enc>::construction::bigbwt(bool delete_T)
{
    log_phase_start("running Big-BWT\n\n");

    ensure_bigbwt_on_path();
    system(("taskset -c 0-$(($(nproc)-1)) bigbwt " +
        (std::string)(supports_locate ? (bigbwt_full_sa ? "-S " : "-s -e ") : "") +
        (std::string)(p > 1 ? ("-t " + std::to_string(p) + " ") : "") +
        prefix_tmp_files + (std::string)(log ? "" : " >log_1 >log_2")).c_str());

    if (!std::filesystem::exists(prefix_tmp_files + ".bwt")) {
        for (const char* ext : { ".log", ".sa", ".ssa", ".esa", ".parse", ".parse_old", ".last",
                 ".bwlast", ".dict", ".ilist", ".occ", ".sai", ".bwsai", ".dicz", ".dicz.len" }) {
            std::error_code ec;
            std::filesystem::remove(prefix_tmp_files + ext, ec);
        }
        if (delete_T) std::filesystem::remove(prefix_tmp_files);
        if (!log) { std::filesystem::remove("log_1"); std::filesystem::remove("log_2"); }
        throw std::runtime_error("Big-BWT failed to build the BWT (the input likely parses into a "
            "single prefix-free phrase, which happens for very short or short-period inputs)");
    }

    std::ifstream log_ifile(prefix_tmp_files + ".log");
    bigbwt_peak_mem_usage = (malloc_count_current() - baseline_mem_usage) + malloc_count_peak_memory_usage(log_ifile);

#ifdef MOVE_R_BENCH
    external_peak_memory_usage = bigbwt_peak_mem_usage;
#endif

    log_ifile.close();
    std::filesystem::remove(prefix_tmp_files + ".log");
    if (delete_T) std::filesystem::remove(prefix_tmp_files);

    if (!log) {
        std::filesystem::remove("log_1");
        std::filesystem::remove("log_2");
    }

    if (log) {
        if (mf_idx != nullptr)
            *mf_idx << " time_pfp" << mf_suffix << "=" << time_diff_ns(time, now());
        std::cout << std::endl;
        time = log_runtime(time);
        log_peak_mem_usage();
    }
}

template <move_r_support support, typename sym_t, typename pos_t, move_pos_encoding_t mlf_enc>
void move_r<support, sym_t, pos_t, mlf_enc>::construction::read_iphim1_bigbwt()
{
    log_phase_start("reading I_Phi^{-1}");

    no_init_resize(I_Phi_m1, r);
    std::ifstream ssa_file(prefix_tmp_files + ".ssa");
    std::ifstream esa_file(prefix_tmp_files + ".esa");

    interleaved_byte_aligned_vectors<uint64_t, pos_t> ssa({ 5, 5 });
    interleaved_byte_aligned_vectors<uint64_t, pos_t> esa({ 5, 5 });
    ssa.resize_no_init(r);
    esa.resize_no_init(r);
    read_from_file(ssa_file, ssa.data(), 10 * r);
    read_from_file(esa_file, esa.data(), 10 * r);

    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < r - 1; i++) {
        I_Phi_m1[i].second = ssa.template get<1>(i);
        I_Phi_m1[i + 1].first = esa.template get<1>(i);
    }

    I_Phi_m1[r - 1].second = ssa.template get<1>(r - 1);
    I_Phi_m1[0].first = esa.template get<1>(r - 1);

    ssa.clear();
    esa.clear();
    ssa.shrink_to_fit();
    esa.shrink_to_fit();

    ssa_file.close();
    esa_file.close();
    std::filesystem::remove(prefix_tmp_files + ".ssa");
    std::filesystem::remove(prefix_tmp_files + ".esa");

    log_phase_end("time_read_iphi");
}