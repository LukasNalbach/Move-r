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
void move_r<support, sym_t, pos_t>::construction::store_rlbwt()
{
    if (log) {
        time = now();
        std::cout << "storing RLBWT to disk" << std::flush;
    }

    std::ofstream file_rlbwt(prefix_tmp_files + ".rlbwt");

    for (uint16_t i = 0; i < p_; i++) {
        RLBWT[i].serialize(file_rlbwt);
    }

    RLBWT.clear();
    RLBWT.shrink_to_fit();
    file_rlbwt.close();

    if (log) {
        if (mf_idx != nullptr)
            *mf_idx << " time_store_rlbwt=" << time_diff_ns(time, now());
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::load_rlbwt()
{
    if (log) {
        time = now();
        std::cout << "loading RLBWT from disk" << std::flush;
    }

    std::ifstream file_rlbwt(prefix_tmp_files + ".rlbwt");
    RLBWT.resize(p_);

    for (uint16_t i = 0; i < p_; i++) {
        RLBWT[i].load(file_rlbwt);
    }

    file_rlbwt.close();
    std::filesystem::remove(prefix_tmp_files + ".rlbwt");

    if (log) {
        if (mf_idx != nullptr)
            *mf_idx << " time_load_rlbwt=" << time_diff_ns(time, now());
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::store_mlf()
{
    if (log) {
        time = now();
        std::cout << "storing M_LF and L' to disk" << std::flush;
    }

    std::ofstream file_mlf(prefix_tmp_files + ".mlf");
    idx._M_LF.serialize(file_mlf);
    idx._M_LF = move_data_structure_l_<pos_t, i_sym_t>();
    file_mlf.close();

    if (log) {
        if (mf_idx != nullptr)
            *mf_idx << " time_store_mlf=" << time_diff_ns(time, now());
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::load_mlf()
{
    if (log) {
        time = now();
        std::cout << "loading M_LF and L' from disk" << std::flush;
    }

    std::ifstream file_mlf(prefix_tmp_files + ".mlf");
    idx._M_LF.load(file_mlf);
    file_mlf.close();
    std::filesystem::remove(prefix_tmp_files + ".mlf");

    if (log) {
        if (mf_idx != nullptr)
            *mf_idx << " time_load_mlf=" << time_diff_ns(time, now());
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::store_iphim1()
{
    if (log) {
        time = now();
        std::cout << "storing I_Phi^{-1} to disk" << std::flush;
    }

    std::ofstream file_iphim1(prefix_tmp_files + ".iphim1");
    write_to_file(file_iphim1, (char*) I_Phi_m1.data(), r * sizeof(std::pair<pos_t, pos_t>));
    I_Phi_m1.clear();
    I_Phi_m1.shrink_to_fit();
    file_iphim1.close();

    if (log) {
        if (mf_idx != nullptr)
            *mf_idx << " time_store_iphim1=" << time_diff_ns(time, now());
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::load_iphim1()
{
    if (log) {
        time = now();
        std::cout << "loading I_Phi^{-1} from disk" << std::flush;
    }

    std::ifstream file_iphim1(prefix_tmp_files + ".iphim1");
    no_init_resize(I_Phi_m1, r);
    read_from_file(file_iphim1, (char*) I_Phi_m1.data(), r * sizeof(std::pair<pos_t, pos_t>));
    file_iphim1.close();
    std::filesystem::remove(prefix_tmp_files + ".iphim1");

    if (log) {
        if (mf_idx != nullptr)
            *mf_idx << " time_load_iphim1=" << time_diff_ns(time, now());
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::store_mphi()
{
    if (log) {
        time = now();
        std::cout << "storing M_Phi to disk" << std::flush;
    }

    std::ofstream file_mphi(prefix_tmp_files + ".mphi");
    idx._M_Phi.serialize(file_mphi);
    idx._M_Phi = move_data_structure<pos_t>();
    file_mphi.close();

    if (log) {
        if (mf_idx != nullptr)
            *mf_idx << " time_store_mphi=" << time_diff_ns(time, now());
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::load_mphi()
{
    if (log) {
        time = now();
        std::cout << "loading M_Phi from disk" << std::flush;
    }

    std::ifstream file_mphi(prefix_tmp_files + ".mphi");
    idx._M_Phi.load(file_mphi);
    file_mphi.close();
    std::filesystem::remove(prefix_tmp_files + ".mphi");

    if (log) {
        if (mf_idx != nullptr)
            *mf_idx << " time_load_mphi=" << time_diff_ns(time, now());
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::store_sas()
{
    if (log) {
        time = now();
        std::cout << "storing SA_s" <<
            (is_bidirectional ? " and SA_s'" : "") <<
            " to disk" << std::flush;
    }

    std::ofstream file_sas(prefix_tmp_files + ".sas");

    if constexpr (support == _locate_move) {
        write_to_file(file_sas, (char*) &SA_s[0], r_ * sizeof(pos_t));
        SA_s.clear();
        SA_s.shrink_to_fit();
    } else {
        idx._SA_s.serialize(file_sas);
        idx._SA_s.clear();
        idx._SA_s.shrink_to_fit();

        idx._SA_s_.serialize(file_sas);
        idx._SA_s_.clear();
        idx._SA_s_.shrink_to_fit();
    }

    file_sas.close();

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::load_sas()
{
    if (log) {
        time = now();
        std::cout << "loading SA_s" <<
            (is_bidirectional ? " and SA_s'" : "") <<
            " from disk" << std::flush;
    }

    std::ifstream file_sas(prefix_tmp_files + ".sas");

    if constexpr (support == _locate_move) {
        no_init_resize(SA_s, r_);
        read_from_file(file_sas, (char*) &SA_s[0], r_ * sizeof(pos_t));
    } else {
        idx._SA_s.load(file_sas);
        idx._SA_s_.load(file_sas);
    }

    file_sas.close();
    std::filesystem::remove(prefix_tmp_files + ".sas");

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::store_sas_idx()
{
    if (log) {
        time = now();
        std::cout << "storing SA_s to disk" << std::flush;
    }

    std::ofstream file_sas(prefix_tmp_files + ".sas");
    idx._SA_s.serialize(file_sas);
    idx._SA_s.clear();
    idx._SA_s.shrink_to_fit();
    file_sas.close();

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::load_sas_idx()
{
    if (log) {
        time = now();
        std::cout << "loading SA_s from disk" << std::flush;
    }

    std::ifstream file_sas(prefix_tmp_files + ".sas");
    idx._SA_s.load(file_sas);
    file_sas.close();
    std::filesystem::remove(prefix_tmp_files + ".sas");

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::store_sas__idx()
{
    if (log) {
        time = now();
        std::cout << "storing SA_s' to disk" << std::flush;
    }

    std::ofstream file_sas_(prefix_tmp_files + ".sas_");
    idx._SA_s_.serialize(file_sas_);
    idx._SA_s_.clear();
    idx._SA_s_.shrink_to_fit();
    file_sas_.close();

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::load_sas__idx()
{
    if (log) {
        time = now();
        std::cout << "loading SA_s' from disk" << std::flush;
    }

    std::ifstream file_sas_(prefix_tmp_files + ".sas_");
    idx._SA_s_.load(file_sas_);
    file_sas_.close();
    std::filesystem::remove(prefix_tmp_files + ".sas_");

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::store_rsl_()
{
    if (log) {
        time = now();
        std::cout << "storing RS_L' to disk" << std::flush;
    }

    std::ofstream file_rsl_(prefix_tmp_files + ".rsl_");
    idx._RS_L_.serialize(file_rsl_);
    file_rsl_.close();
    idx._RS_L_ = rsl_t();

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::load_rsl_()
{
    if (log) {
        time = now();
        std::cout << "loading RS_L' from disk" << std::flush;
    }

    std::ifstream file_rsl_(prefix_tmp_files + ".rsl_");
    idx._RS_L_.load(file_rsl_);
    file_rsl_.close();
    std::filesystem::remove(prefix_tmp_files + ".rsl_");

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::store_l_prev_next()
{
    if (log) {
        time = now();
        std::cout << "storing L'_prev and L'_next to disk" << std::flush;
    }

    std::ofstream file_lprevnext(prefix_tmp_files + ".lprevnext");
    idx._L_prev.serialize(file_lprevnext);
    idx._L_next.serialize(file_lprevnext);
    file_lprevnext.close();
    idx._L_prev = interleaved_byte_aligned_vectors<pos_t, pos_t>();
    idx._L_next = interleaved_byte_aligned_vectors<pos_t, pos_t>();

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::load_l_prev_next()
{
    if (log) {
        time = now();
        std::cout << "loading L'_prev and L'_next from disk" << std::flush;
    }

    std::ifstream file_lprevnext(prefix_tmp_files + ".lprevnext");
    idx._L_prev.load(file_lprevnext);
    idx._L_next.load(file_lprevnext);
    file_lprevnext.close();
    std::filesystem::remove(prefix_tmp_files + ".lprevnext");

    if (log) {
        time = log_runtime(time);
    }
}