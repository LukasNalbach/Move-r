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

#include <gtl/phmap.hpp>
#include <move_r/move_r.hpp>

template <move_r_support support, typename sym_t, typename pos_t>
template <bool bigbwt, typename sad_t, typename sa_sint_t>
void move_r<support, sym_t, pos_t>::construction::build_freq_sad()
{
    if (log) {
        time = now();
        std::cout << std::endl;
        std::cout << "building rlzsa:" << std::endl;
        std::cout << "computing frequencies of values in SA^d" << std::flush;
    }

    sad_freq_t<sad_t>& SAd_freq = get_SAd_freq<sad_t>();
    n_u64 = n;

    for (uint16_t i = 0; i < p; i++) {
        n_p.emplace_back(i * (n / p));
    }

    n_p.emplace_back(n);
    sad_t n_sad = n;
    SAd_freq.reserve(r);

    for (pos_t i = 0; i < r; i++) {
        sad_t sad = (sad_t { I_Phi_m1[i].second } + n_sad) - sad_t { I_Phi_m1[i].first };
        pos_t count = (i == r - 1 ? n : I_Phi_m1[i + 1].first) - I_Phi_m1[i].first;
        auto it = SAd_freq.find(sad);

        if (it == SAd_freq.end()) {
            SAd_freq.insert_unique(sad, count);
        } else {
            (*it).second += count;
        }
    }

    I_Phi_m1.clear();
    I_Phi_m1.shrink_to_fit();

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
template <bool bigbwt, typename sad_t, typename sa_sint_t>
void move_r<support, sym_t, pos_t>::construction::build_r()
{
    if (log) {
        std::cout << "choosing segments for R" << std::flush;
    }

    sad_freq_t<sad_t>& SAd_freq = get_SAd_freq<sad_t>();
    idx._R = interleaved_byte_aligned_vectors<uint64_t, pos_t>({ byte_width(2 * n + 1) });
    std::mt19937 mt(std::random_device {}());
    num_cand_segs = 5 * std::pow(n / (float) r, 0.45);
    pos_t num_cand_segs_thr = num_cand_segs / p + 1; // number of considered candidate segments per thread
    std::vector<gtl::flat_hash_set<sad_t, std::identity>> PV_thr(p);
    std::vector<std::uniform_int_distribution<pos_t>*> pos_distrib;

    for (uint16_t i_p = 0; i_p < p; i_p++) {
        pos_distrib.emplace_back(new std::uniform_int_distribution<pos_t>(
            n_p[i_p],
            std::max<int64_t>(n_p[i_p], int64_t { n_p[i_p + 1] } - seg_size)));
        if constexpr (bigbwt)
            SA_file_bufs[i_p].buffersize(8 * seg_size);
    }

    pos_t size_R_pre_target = 0.95 * size_R_target;

    while (size_R < size_R_pre_target) {
        pos_t beg_best = 0;
        pos_t end_best = 0;
        float score_best = 0;
        ts_it_t it_best = T_s.end();
        uint16_t ip_best = 0;

        #pragma omp parallel num_threads(p)
        {
            uint16_t i_p = omp_get_thread_num();

            pos_t beg_best_thr = 0;
            pos_t end_best_thr = 0;
            float score_best_thr = 0;
            ts_it_t it_best_thr = T_s.end();

            pos_t beg, end;
            float score;
            ts_it_t it = T_s.end();

            for (pos_t seg = 0; seg < num_cand_segs_thr; seg++) {
                uint8_t tries = 0;

                do {
                    tries++;
                    beg = (*pos_distrib[i_p])(mt);
                    end = beg + seg_size;
                    score = 0;
                    it = T_s.lower_bound(segment { beg, 0 });

                    if (it != T_s.end()) {
                        pos_t start_next = (*it).beg;

                        if (start_next < end) {
                            end = start_next;
                        }
                    }

                    if (it != T_s.begin()) {
                        auto it_prev = it;
                        --it_prev;
                        pos_t end_last = (*it_prev).end;

                        if (beg < end_last) {
                            if (end < end_last) {
                                end = beg;
                            } else {
                                beg = end_last;
                            }
                        }
                    }
                } while (tries < 8 && end == beg);

                for (pos_t i = beg; i < end; i++) {
                    sad_t val = SAd<bigbwt, sa_sint_t>(i_p, i);
                    pos_t freq = (*SAd_freq.find(val)).second;

                    if (freq != 0 && PV_thr[i_p].emplace(val).second) {
                        score += std::sqrt(freq);
                    }
                }

                score /= end - beg;
                PV_thr[i_p].clear();

                if (score > score_best_thr) {
                    beg_best_thr = beg;
                    end_best_thr = end;
                    score_best_thr = score;
                    it_best_thr = it;
                }
            }

            #pragma omp critical
            {
                if (score_best_thr > score_best) {
                    beg_best = beg_best_thr;
                    end_best = end_best_thr;
                    score_best = score_best_thr;
                    it_best = it_best_thr;
                    ip_best = i_p;
                }
            }
        }

        if (score_best == 0) {
            break;
        }

        bool merged = false;

        if (it_best != T_s.end() && end_best == (*it_best).beg) {
            merged = true;

            if (it_best != T_s.begin()) {
                auto it_prev = it_best;
                --it_prev;

                if ((*it_prev).end == beg_best) {
                    (*it_prev).end = (*it_best).end;
                    T_s.erase(it_best);
                } else {
                    (*it_best).beg = beg_best;
                }
            } else {
                (*it_best).beg = beg_best;
            }
        } else if (!T_s.empty() && it_best != T_s.begin()) {
            auto it_prev = it_best;
            --it_prev;

            if ((*it_prev).end == beg_best) {
                (*it_prev).end = end_best;
                merged = true;
            }
        }

        if (!merged) {
            T_s.emplace_hint(it_best, segment { beg_best, end_best });
        }

        for (pos_t i = beg_best; i < end_best; i++) {
            SAd_freq.find(SAd<bigbwt, sa_sint_t>(ip_best, i))->second = 0;
        }

        size_R += end_best - beg_best;
    }

    SAd_freq.clear();
    SAd_freq.shrink_to_fit();

    if (log) {
        time = log_runtime(time);
        std::cout << "num. of segments: " << T_s.size() << std::endl;
        std::cout << "closing gaps between segments" << std::flush;
    }

    {
        auto it = T_s.begin();

        while (it != T_s.end()) {
            pos_t beg_last = (*it).beg;
            pos_t end_last = (*it).end;
            ++it;

            if (it != T_s.end()) {
                T_g.emplace(gap {
                    .beg_prev = beg_last,
                    .score = ((*it).end - beg_last) / (float) ((*it).beg - end_last) });
            }
        }
    }

    while (!T_g.empty()) {
        auto tg_min = T_g.begin();
        auto ts_left = T_s.find(segment { (*tg_min).beg_prev, 0 });
        T_g.erase(tg_min);

        auto ts_right = ts_left;
        ++ts_right;
        pos_t len_gap = (*ts_right).beg - (*ts_left).end;

        if (size_R + len_gap <= size_R_target) {
            size_R += len_gap;

            if (ts_left != T_s.begin()) {
                auto ts_prev = ts_left;
                --ts_prev;

                pos_t len_prev_gap = (*ts_left).beg - (*ts_prev).end;
                float old_score = ((*ts_left).end - (*ts_prev).beg) / (float) len_prev_gap;

                if (T_g.erase(gap { (*ts_prev).beg, old_score }) && size_R + len_prev_gap <= size_R_target) {
                    float new_score = ((*ts_right).end - (*ts_prev).beg) / (float) len_prev_gap;
                    T_g.emplace(gap { (*ts_prev).beg, new_score });
                }
            }

            auto ts_next = ts_right;
            ++ts_next;

            if (ts_next != T_s.end()) {
                pos_t len_next_gap = (*ts_next).beg - (*ts_right).end;
                float old_score = ((*ts_next).end - (*ts_right).beg) / (float) len_next_gap;

                if (T_g.erase(gap { (*ts_right).beg, old_score }) && size_R + len_next_gap <= size_R_target) {
                    float new_score = ((*ts_next).end - (*ts_left).beg) / (float) len_next_gap;
                    T_g.emplace(gap { (*ts_left).beg, new_score });
                }
            }

            (*ts_left).end = (*ts_right).end;
            T_s.erase(ts_right);
        }
    }

    if (log) {
        time = log_runtime(time);
        std::cout << "num. of segments: " << T_s.size() << std::endl;
        std::cout << "building R" << std::flush;
    }

    idx._R.resize_no_init(size_R);
    uint64_t i = 0;

    for (auto seg : T_s) {
        #pragma omp parallel for num_threads(p)
        for (uint64_t j = seg.beg; j < seg.end; j++) {
            idx._R.template set_parallel<0, uint64_t>(i + (j - seg.beg), SAd<bigbwt, sa_sint_t>(omp_get_thread_num(), j));
        }

        i += seg.end - seg.beg;
    }

    T_s.clear();

    if (log) {
        time = log_runtime(time);
        std::cout << "building rev(R)" << std::flush;
    }

    std::vector<sad_t>& revR = get_revR<sad_t>();
    no_init_resize(revR, size_R);

    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < size_R; i++) {
        revR[i] = idx._R[size_R - i - 1];
    }

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::store_r()
{
    if (log) {
        std::cout << "storing R to disk" << std::flush;
    }

    std::ofstream tmp_file(prefix_tmp_files + "_R");
    idx._R.serialize(tmp_file);
    idx._R.clear();
    idx._R.shrink_to_fit();
    tmp_file.close();

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
template <typename sad_t, typename irr_pos_t>
void move_r<support, sym_t, pos_t>::construction::build_idx_rev_r()
{
    if (log) {
        std::cout << "building move-r of rev(R):" << std::endl
                  << std::endl;
    }

    std::vector<sad_t>& revR = get_revR<sad_t>();
    idx_revr_t<sad_t, irr_pos_t>& idx_revR = get_idx_revR<sad_t, irr_pos_t>();

    idx_revR = idx_revr_t<sad_t, irr_pos_t>(std::move(revR), {
        .mode = mode == _suffix_array ? _suffix_array : _suffix_array_space,
        .num_threads = p, .log = log
    });

    if (log) {
        std::cout << std::endl;
        time = log_runtime(time);
        log_peak_mem_usage();
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
template <bool bigbwt, bool space, typename sad_t, typename irr_pos_t, typename sa_sint_t>
void move_r<support, sym_t, pos_t>::construction::build_rlzsa_factorization()
{
    if (log) {
        std::cout << "computing rlzsa factorization" << std::flush;
    }

    idx_revr_t<sad_t, irr_pos_t>& idx_revR = get_idx_revR<sad_t, irr_pos_t>();

    uint8_t width_sr = byte_width(size_R) * 8;
    uint8_t width_lp = byte_width(n) * 8;

    if constexpr (space) {
        for (uint16_t i = 0; i < p; i++) {
            CPL_file_bufs.emplace_back(sdsl::int_vector_buffer<>(
                prefix_tmp_files + ".cpl_" + std::to_string(i),
                std::ios::out, 128 * 1024, 16, false));

            SR_file_bufs.emplace_back(sdsl::int_vector_buffer<>(
                prefix_tmp_files + ".sr_" + std::to_string(i),
                std::ios::out, 128 * 1024, width_sr, false));

            LP_file_bufs.emplace_back(sdsl::int_vector_buffer<>(
                prefix_tmp_files + ".lp_" + std::to_string(i),
                std::ios::out, 128 * 1024, width_lp, false));

            PT_file_bufs.emplace_back(sdsl::int_vector_buffer<>(
                prefix_tmp_files + ".pt_" + std::to_string(i),
                std::ios::out, 128 * 1024, 1, false));

            SA_file_bufs[i].buffersize(128 * 1024);
        }
    } else {
        _CPL.resize(p);
        _SR.resize(p, interleaved_byte_aligned_vectors<pos_t, pos_t>({ width_sr / 8 }));
        _LP.resize(p, interleaved_byte_aligned_vectors<pos_t, pos_t>({ width_lp / 8 }));
        _PT.resize(p);
    }

    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        pos_t b = n_p[i_p];
        pos_t e = n_p[i_p + 1];
        pos_t j = 1;
        pos_t pt_size = std::max<pos_t>(2, r / p);

        if constexpr (space) {
            PT_file_bufs[i_p].push_back(1);
            LP_file_bufs[i_p].push_back(SA<bigbwt, sa_sint_t>(i_p, b));
        } else {
            _PT[i_p].resize(pt_size);
            _PT[i_p][0] = 1;
            _LP[i_p].emplace_back(SA<bigbwt, sa_sint_t>(i_p, b));
        }

        auto query = idx_revR.query();

        for (pos_t i = b + 1; i < e;) {
            query.reset();
            pos_t max_query_len = std::min<pos_t>(e - i, 65535);

            while (query.prepend(SAd<bigbwt, sa_sint_t>(i_p, i + query.length())) &&
                   query.length() < max_query_len &&
                   (query.length() <= 1 || query.num_occ() > 1));

            if (query.length() <= 1) {
                if constexpr (space) {
                    PT_file_bufs[i_p].push_back(1);
                    LP_file_bufs[i_p].push_back(SA<bigbwt, sa_sint_t>(i_p, i));
                } else {
                    _PT[i_p][j] = 1;
                    _LP[i_p].emplace_back(SA<bigbwt, sa_sint_t>(i_p, i));
                }

                i++;
            } else {
                pos_t occ = size_R - query.one_occ() - query.length();
                pos_t len = query.length();

                if (query.num_occ() == 1) {
                    max_query_len = std::min<pos_t>(max_query_len, size_R - occ);

                    while (len < max_query_len &&
                           idx._R[occ + len] == SAd<bigbwt, sa_sint_t>(i_p, i + len))
                    {
                        len++;
                    }
                }

                if constexpr (space) {
                    PT_file_bufs[i_p].push_back(0);
                    CPL_file_bufs[i_p].push_back(len);
                    SR_file_bufs[i_p].push_back(occ);
                } else {
                    _PT[i_p][j] = 0;
                    _CPL[i_p].emplace_back(len);
                    _SR[i_p].emplace_back(occ);
                }

                i += len;
            }

            if constexpr (!space) {
                j++;

                if (j == pt_size) {
                    pt_size *= 2;
                    _PT[i_p].resize(pt_size);
                }
            }
        }

        if constexpr (!space)
            _PT[i_p].resize(j);
    }

    idx_revR = idx_revr_t<sad_t, irr_pos_t>();

    z_p.emplace_back(0);
    zl_p.emplace_back(0);
    zc_p.emplace_back(0);

    for (uint16_t i = 0; i < p; i++) {
        if constexpr (space) {
            idx.z += PT_file_bufs[i].size();
            idx.z_l += LP_file_bufs[i].size();
        } else {
            idx.z += _PT[i].size();
            idx.z_l += _LP[i].size();
        }

        z_p.emplace_back(idx.z);
        zl_p.emplace_back(idx.z_l);
        zc_p.emplace_back(idx.z - idx.z_l);
    }

    idx.z_c = zc_p[p];
    idx._CPL.reserve(idx.z_c + 2);
    no_init_resize(idx._CPL, idx.z_c);

    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        pos_t b_zc = zc_p[i_p];
        pos_t zc_diff = zc_p[i_p + 1] - b_zc;

        for (pos_t i = 0; i < zc_diff; i++) {
            idx._CPL[b_zc + i] = CPL<space>(i_p, i);
        }
    }

    idx._CPL.emplace_back(1);
    idx._CPL.emplace_back(1);
    _CPL.clear();
    _CPL.shrink_to_fit();

    idx._SR = interleaved_byte_aligned_vectors<pos_t, pos_t>({ width_sr / 8 });
    idx._SR.resize_no_init(idx.z_c + 1);

    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        pos_t b_zc = zc_p[i_p];
        pos_t zc_diff = zc_p[i_p + 1] - b_zc;

        for (pos_t i = 0; i < zc_diff; i++) {
            idx._SR.template set_parallel<0, pos_t>(b_zc + i, SR<space>(i_p, i));
        }
    }

    idx._SR.template set_parallel<0, pos_t>(idx.z_c, size_R);

    idx._LP = interleaved_byte_aligned_vectors<pos_t, pos_t>({ width_lp / 8 });
    idx._LP.resize_no_init(idx.z_l, p);

    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        pos_t b_zl = zl_p[i_p];
        pos_t zl_diff = zl_p[i_p + 1] - b_zl;

        for (pos_t i = 0; i < zl_diff; i++) {
            idx._LP.template set_parallel<0, pos_t>(b_zl + i, LP<space>(i_p, i));
        }
    }

    sdsl::bit_vector PT_tmp;
    PT_tmp.resize(idx.z + 2);

    for (uint16_t i_p = 0; i_p < p; i_p++) {
        pos_t b_z = z_p[i_p];
        pos_t z_diff = z_p[i_p + 1] - b_z;

        for (pos_t i = 0; i < z_diff; i++) {
            PT_tmp[b_z + i] = PT<space>(i_p, i);
        }
    }

    PT_tmp[idx.z] = 0;
    PT_tmp[idx.z + 1] = 0;
    idx._PT = plain_bit_vector<pos_t, true, true, true>(std::move(PT_tmp));

    if constexpr (space) {
        for (uint16_t i = 0; i < p; i++) {
            CPL_file_bufs[i].close(true);
            SR_file_bufs[i].close(true);
            LP_file_bufs[i].close(true);
            PT_file_bufs[i].close(true);
        }
    }

    if (idx.z_c == 0) {
        sdsl::sd_vector_builder SCP_S_b(n + 1, 1);
        SCP_S_b.set(n);
        idx._SCP_S = sd_array<pos_t>(sdsl::sd_vector<>(SCP_S_b));
    } else {
        pos_t n_s = std::max<pos_t>(1, idx.z_c / sample_rate_scp);
        sdsl::sd_vector_builder SCP_S_b(n + 1, n_s + 1);
        pos_t i_cp = 0;
        pos_t i_p = idx._PT.select_0(1);
        pos_t s_cp = i_p;

        for (pos_t i = 0; i < n_s - 1; i++) {
            SCP_S_b.set(s_cp);

            for (pos_t j = 0; j < sample_rate_scp; j++) {
                pos_t i_np = idx._PT.select_0(i_cp + 2);
                s_cp += idx._CPL[i_cp] + (i_np - i_p - 1);
                i_p = i_np;
                i_cp++;
            }
        }

        SCP_S_b.set(s_cp);
        SCP_S_b.set(n);
        idx._SCP_S = sd_array<pos_t>(sdsl::sd_vector<>(SCP_S_b));
    }

    if (log) {
        time = log_runtime(time);
        std::cout << "z: " << idx.z << std::endl;
        std::cout << "z_l/z: " << idx.z_l / (double)idx.z << std::endl;
        std::cout << "z/r: " << idx.z / (double)r << std::endl;
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
template <bool bigbwt, typename sa_sint_t>
void move_r<support, sym_t, pos_t>::construction::build_sa_delta()
{
    if (idx.delta == 0) {
        uint64_t target_sampling_size = idx.size_in_bytes() * default_relative_sampling_size;
        idx.delta = n / (target_sampling_size / byte_width(n));
    }

    if (log) std::cout << "building SA_delta, delta = " << idx.delta << std::flush;
    
    time = now();
    uint64_t num_samples = n / idx.delta;
    idx._SA_delta = interleaved_byte_aligned_vectors<pos_t, pos_t>({ byte_width(n) });
    idx._SA_delta.resize_no_init(num_samples);
    idx.last_sa = SA<bigbwt, sa_sint_t>(0, n - 1);

    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < num_samples; i++) {
        idx._SA_delta.template set_parallel<0, pos_t>(i, SA<bigbwt, sa_sint_t>(omp_get_thread_num(), i * idx.delta));
    }

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::construction::load_r()
{
    if (log) {
        std::cout << "reading R from disk" << std::flush;
    }

    std::ifstream tmp_file(prefix_tmp_files + "_R");
    idx._R.load(tmp_file);
    tmp_file.close();
    std::filesystem::remove(prefix_tmp_files + "_R");

    if (log) {
        time = log_runtime(time);
    }
}