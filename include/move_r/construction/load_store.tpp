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

#ifndef MOVE_R_MOVE_R_CONSTRUCTION_LOAD_STORE_TPP
#define MOVE_R_MOVE_R_CONSTRUCTION_LOAD_STORE_TPP

#include <move_r/move_r.hpp>

template <move_r_support support, move_pos_encoding_t mlf_enc>
template <typename stream_t, typename fnc_t>
void move_r<support, mlf_enc>::construction::log_store_or_load(const char* verb, const char* name, stream_t& stream, fnc_t fnc)
{
    log_phase_start("");
    auto pos = [&] {
        if constexpr (std::is_base_of_v<std::ostream, stream_t>) return stream.num_bytes_written();
        else return stream.read_pos();
    };
    const uint64_t before = pos();
    fnc();
    const bool storing = std::string(verb) == "storing";
    end_spill_phase(std::string(storing ? "storing " : "reading ") + name, pos() - before, storing,
        std::string("time_") + (storing ? "store_" : "load_") + name);
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::store_mlf()
{
    store_to_tmp_file("M_LF and L'", ".mlf", "time_store_mlf", [&](direct_ofstream& file) {
        idx._M_LF.serialize(file);
        idx._M_LF = mlf_t();
    });
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::load_mlf()
{
    load_from_tmp_file("M_LF and L'", ".mlf", "time_load_mlf", true, [&](direct_ifstream& file) {
        idx._M_LF.load(file, idx.mlf_num_rows);
    });
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::store_data_structures_over_l_()
{
    direct_ofstream file(prefix_tmp_files + ".l_");
    log_store_or_load("storing", "PSB_L'", file, [&] {
        idx._PSB_L_.serialize(file); idx._PSB_L_ = psb_l_t(); });
    log_store_or_load("storing", "PSO_L'", file, [&] {
        idx._PSO_L_.serialize(file); idx._PSO_L_ = pso_l_t(); });
    log_store_or_load("storing", "L_prev_min", file, [&] {
        idx._L_prev_min.serialize(file); idx._L_prev_min = L_prev_min_t(); });
    file.close();
    drop_from_page_cache(prefix_tmp_files + ".l_");
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::load_data_structures_over_l_()
{
    direct_ifstream file(prefix_tmp_files + ".l_");
    log_store_or_load("loading", "PSB_L'", file, [&] {
        idx._PSB_L_.load(file); idx._PSB_L_.set_read(L_reader { &idx }); });
    log_store_or_load("loading", "PSO_L'", file, [&] {
        idx._PSO_L_.load(file); idx._PSO_L_.set_read(L_reader { &idx }); });
    log_store_or_load("loading", "L_prev_min", file, [&] { idx._L_prev_min.load(file); });
    file.close();
    std::filesystem::remove(prefix_tmp_files + ".l_");
}

#endif // MOVE_R_MOVE_R_CONSTRUCTION_LOAD_STORE_TPP
