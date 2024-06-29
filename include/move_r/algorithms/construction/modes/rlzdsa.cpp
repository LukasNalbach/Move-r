#pragma once

#include <move_r/move_r.hpp>

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
template <bool bigbwt, typename sad_t, typename sa_sint_t>
void move_r<locate_support,sym_t,pos_t>::construction::build_freq_sad() {
    if (log) {
        time = now();
        std::cout << std::endl;
        std::cout << "building rlzdsa:" << std::endl;
        std::cout << "computing frequencies of values in SA^d" << std::flush;
    }

    if constexpr (bigbwt) {
        for (uint16_t i=0; i<p; i++) {
            SA_file_bufs.emplace_back(sdsl::int_vector_buffer<>(
                prefix_tmp_files + ".sa", std::ios::in,
                128*1024, 40, true
            ));
        }
    }

    std::vector<sad_freq_t<sad_t>>& SAd_freq = get_SAd_freq<sad_t>();
    n_u64 = n;
    
    for (uint16_t i=0; i<p; i++) {
        n_p.emplace_back(i*(n/p));
    }

    n_p.emplace_back(n);

    SAd_freq.resize(p);
    sad_t n_sad = n;

    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        pos_t b = n_p[i_p];
        pos_t e = n_p[i_p+1];
        pos_t eb_diff = e-b;
        sad_t cur_sa,last_sa,cur_sad;

        if (i_p == 0) {
            SAd_freq[0].try_emplace(n-1,1);
            last_sa = n-1;
        } else {
            last_sa = SA<bigbwt,sa_sint_t>(i_p,b);
        }

        for (pos_t i=1; i<eb_diff; i++) {
            cur_sa = SA<bigbwt,sa_sint_t>(i_p,b+i);
            cur_sad = (cur_sa+n_sad)-last_sa;
            last_sa = cur_sa;

            auto it = SAd_freq[i_p].find(cur_sad);

            if (it == SAd_freq[i_p].end()) {
                SAd_freq[i_p].try_emplace(cur_sad,1);
            } else {
                (*it).second++;
            }
        }
    }

    for (uint16_t j=1; j<=std::ceil(std::log2(p)); j++) {
        #pragma omp parallel for num_threads(p)
        for (uint16_t i_p=0; i_p<p/(uint16_t)std::pow(2,j); i_p++) {
            std::vector<std::pair<uint16_t,uint16_t>> merges;
            merges.emplace_back(std::make_pair(std::pow(2,j)*i_p,std::pow(2,j)*i_p+std::pow(2,j-1)));

            if (p%(uint16_t)std::pow(2,j) != 0 && i_p == p/(uint16_t)std::pow(2,j)-1) {
                merges.emplace_back(std::make_pair(std::pow(2,j)*i_p,std::pow(2,j)*(i_p+1)));
            }

            for (auto m : merges) {
                for (auto p : SAd_freq[m.second]) {
                    auto it = SAd_freq[m.first].find(p.first);

                    if (it == SAd_freq[m.first].end()) {
                        SAd_freq[m.first].emplace(p);
                    } else {
                        (*it).second += p.second;
                    }
                }

                SAd_freq[m.second] = sad_freq_t<sad_t>();
            }
        }
    }

    SAd_freq.resize(1);
    SAd_freq.shrink_to_fit();
    sigma_SAd = SAd_freq[0].size();

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
template <bool bigbwt, typename sad_t, typename sa_sint_t>
void move_r<locate_support,sym_t,pos_t>::construction::build_r_revR() {
    if (log) {
        std::cout << "building R" << std::flush;
    }

    std::vector<sad_freq_t<sad_t>>& SAd_freq = get_SAd_freq<sad_t>();

    idx._R = std::move(interleaved_vectors<uint64_t,pos_t>({(uint8_t)std::ceil(std::log2(2*n+1)/(double)8)}));
    std::mt19937 mt(std::random_device{}());
    num_cand_segs = 4.5*powf(n/(float)r,0.6);
    pos_t num_cand_segs_thr = num_cand_segs/p+1; // number of considered candidate segments per thread
    std::vector<gtl::flat_hash_set<sad_t>> PV_thr(p);
    std::vector<std::uniform_int_distribution<pos_t>*> pos_distrib;

    for (uint16_t i_p=0; i_p<p; i_p++) {
        pos_distrib.emplace_back(new std::uniform_int_distribution<pos_t>(n_p[i_p],n_p[i_p+1]));
        if constexpr (bigbwt) SA_file_bufs[i_p].buffersize(8*seg_size); 
    }

    while (size_R < size_R_target) {
        pos_t pos_best = 0;
        pos_t len_best = 0;
        float score_best = 0;
        ts_it_t it_best = T_s.end();
        uint16_t ip_best = 0;

        #pragma omp parallel num_threads(p)
        {
            uint16_t i_p = omp_get_thread_num();

            pos_t pos_best_thr = 0;
            pos_t len_best_thr = 0;
            float score_best_thr = 0;
            ts_it_t it_best_thr = T_s.end();

            pos_t pos,len;
            float score;
            ts_it_t it = T_s.end();

            for (pos_t seg=0; seg<num_cand_segs_thr; seg++) {
                pos = (*pos_distrib[i_p])(mt);
                len = std::min<pos_t>(seg_size,n-pos);
                score = 0;
                it = T_s.end();

                if (!T_s.empty()) {
                    it = T_s.upper_bound(std::make_pair(pos,0));
                    
                    if (it != T_s.end()) {
                        pos_t start_next = (*it).first;

                        if (start_next < pos+len) {
                            len = start_next-pos;
                        }
                    }

                    if (it != T_s.begin()) {
                        auto it_before = it;
                        it_before--;
                        pos_t end_last = (*it_before).first+(*it_before).second;

                        if (pos < end_last) {
                            if (pos+len < end_last) {
                                len = 0;
                            } else {
                                len -= end_last-pos;
                                pos = end_last;
                            }
                        }
                    }
                }

                for (pos_t i=0; i<len; i++) {
                    sad_t val = SAd<bigbwt,sa_sint_t>(i_p,pos+i);
                    pos_t freq = SAd_freq[0][val];

                    if (freq != 0 && PV_thr[i_p].emplace(val).second) {
                        score += std::sqrt(freq);
                    }
                }

                score /= len;
                PV_thr[i_p].clear();

                if (score > score_best_thr) {
                    pos_best_thr = pos;
                    len_best_thr = len;
                    score_best_thr = score;
                    it_best_thr = it;
                }
            }

            #pragma omp critical
            {
                if (score_best_thr > score_best) {
                    pos_best = pos_best_thr;
                    len_best = len_best_thr;
                    score_best = score_best_thr;
                    it_best = it_best_thr;
                    ip_best = i_p;
                }
            }
        }

        if (score_best == 0) {
            break;
        }
        
        T_s.emplace_hint(it_best,std::make_pair(pos_best,len_best));

        for (pos_t i=0; i<len_best; i++) {
            SAd_freq[0][SAd<bigbwt,sa_sint_t>(ip_best,pos_best+i)] = 0;
        }

        size_R += len_best;
    }

    SAd_freq.clear();
    SAd_freq.shrink_to_fit();
    idx._R.resize_no_init(size_R);
    uint64_t i = 0;

    for (auto seg : T_s) {
        #pragma omp parallel for num_threads(p)
        for (uint64_t j=0; j<seg.second; j++) {
            idx._R.template set<0>(i+j,SAd<bigbwt,sa_sint_t>(omp_get_thread_num(),seg.first+j));
        }

        i += seg.second;
    }

    T_s.clear();

    if (log) {
        time = log_runtime(time);
        std::cout << "building rev(R)" << std::flush;
    }

    std::vector<sad_t>& revR = get_revR<sad_t>();
    no_init_resize(revR,size_R);
    
    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<size_R; i++) {
        revR[i] = idx._R[size_R-i-1];
    }

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
void move_r<locate_support,sym_t,pos_t>::construction::store_r() {
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

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
template <typename sad_t, typename irr_pos_t>
void move_r<locate_support,sym_t,pos_t>::construction::build_idx_rev_r() {
    if (log) {
        std::cout << "building move-r of rev(R)" << std::flush;
    }

    std::vector<sad_t>& revR = get_revR<sad_t>();
    idx_revr_t<sad_t,irr_pos_t>& idx_revR = get_idx_revR<sad_t,irr_pos_t>();

    idx_revR = idx_revr_t<sad_t,irr_pos_t>(std::move(revR),{
        .mode = (mode == _bigbwt || mode == _suffix_array_space) ? _suffix_array_space : _suffix_array,
        .num_threads = p
    });

    if (log) {
        std::cout << " (size: " << format_size(idx_revR.size_in_bytes()) << ")";
        time = log_runtime(time);
        log_peak_mem_usage();
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
template <bool bigbwt, bool space, typename sad_t, typename irr_pos_t, typename sa_sint_t>
void move_r<locate_support,sym_t,pos_t>::construction::build_rlzdsa_factorization() {
    if (log) {
        std::cout << "computing rlzdsa factorization" << std::flush;
    }

    idx_revr_t<sad_t,irr_pos_t>& idx_revR = get_idx_revR<sad_t,irr_pos_t>();

    uint8_t width_scp = std::ceil(std::log2(n+1)/(double)8)*8;
    uint8_t width_sr = std::ceil(std::log2(size_R+1)/(double)8)*8;
    uint8_t width_lp = std::ceil(std::log2(n+1)/(double)8)*8;

    if constexpr (space) {
        for (uint16_t i=0; i<p; i++) {
            SCP_file_bufs.emplace_back(sdsl::int_vector_buffer<>(
                prefix_tmp_files + ".scp_" + std::to_string(i),
                std::ios::out, 128*1024, width_scp, false
            ));

            SR_file_bufs.emplace_back(sdsl::int_vector_buffer<>(
                prefix_tmp_files + ".sr_" + std::to_string(i),
                std::ios::out, 128*1024, width_sr, false
            ));

            LP_file_bufs.emplace_back(sdsl::int_vector_buffer<>(
                prefix_tmp_files + ".lp_" + std::to_string(i),
                std::ios::out, 128*1024, width_lp, false
            ));

            PT_file_bufs.emplace_back(sdsl::int_vector_buffer<>(
                prefix_tmp_files + ".pt_" + std::to_string(i),
                std::ios::out, 128*1024, 1, false
            ));

            SA_file_bufs[i].buffersize(128*1024);
        }
    } else {
        SCP.resize(p,interleaved_vectors<pos_t,pos_t>({width_scp/8}));
        SR.resize(p,interleaved_vectors<pos_t,pos_t>({width_sr/8}));
        LP.resize(p,interleaved_vectors<pos_t,pos_t>({width_lp/8}));
        PT.resize(p);
    }

    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        pos_t b = n_p[i_p];
        pos_t e = n_p[i_p+1];
        pos_t j = 1;
        pos_t pt_size = std::max<pos_t>(2,r/p);

        if constexpr (space) {
            PT_file_bufs[i_p].push_back(1);
            LP_file_bufs[i_p].push_back(SA<bigbwt,sa_sint_t>(i_p,b));
        } else {
            PT[i_p].resize(pt_size);
            PT[i_p][0] = 1;
            LP[i_p].emplace_back(SA<bigbwt,sa_sint_t>(i_p,b));
        }
        
        auto query = idx_revR.query();
        
        for (pos_t i=b+1; i<e;) {
            query.reset();
            while (query.prepend(SAd<bigbwt,sa_sint_t>(i_p,i+query.length())) && i+query.length() < e);

            if (query.length() <= 1) {
                if constexpr (space) {
                    PT_file_bufs[i_p].push_back(1);
                    LP_file_bufs[i_p].push_back(SA<bigbwt,sa_sint_t>(i_p,i));
                } else {
                    PT[i_p][j] = 1;
                    LP[i_p].emplace_back(SA<bigbwt,sa_sint_t>(i_p,i));
                }

                i++;
            } else {
                if constexpr (space) {
                    PT_file_bufs[i_p].push_back(0);
                    SCP_file_bufs[i_p].push_back(i);
                    SR_file_bufs[i_p].push_back(size_R-query.next_occ()-query.length());
                } else {
                    PT[i_p][j] = 0;
                    SCP[i_p].emplace_back(i);
                    SR[i_p].emplace_back(size_R-query.next_occ()-query.length());
                }

                i += query.length();
            }
            
            if constexpr (!space) {
                j++;

                if (j == pt_size) {
                    pt_size *= 2;
                    PT[i_p].resize(pt_size);
                }
            }
        }

        if constexpr (!space) PT[i_p].resize(j);
    }

    idx_revR = idx_revr_t<sad_t,irr_pos_t>();

    z_p.emplace_back(0);
    zl_p.emplace_back(0);
    zc_p.emplace_back(0);

    for (uint16_t i=0; i<p; i++) {
        if constexpr (space) {
            idx.z += PT_file_bufs[i].size();
            idx.z_l += LP_file_bufs[i].size();
        } else {
            idx.z += PT[i].size();
            idx.z_l += LP[i].size();
        }

        z_p.emplace_back(idx.z);
        zl_p.emplace_back(idx.z_l);
        zc_p.emplace_back(idx.z-idx.z_l);
    }

    idx.z_c = zc_p[p];
    
    if constexpr (space) {
        SCP_b = sdsl::sd_vector_builder(n+2,idx.z_c+2);

        for (uint16_t i_p=0; i_p<p; i_p++) {
            pos_t zc_diff = zc_p[i_p+1]-zc_p[i_p];

            for (pos_t i=0; i<zc_diff; i++) {
                SCP_b.set(SCP_file_bufs[i_p][i]);
            }
        }

        SCP_b.set(n);
        SCP_b.set(n+1);
        idx._SCP = std::move(sd_array<pos_t>(std::move(sdsl::sd_vector<>(SCP_b))));

        idx._SR = interleaved_vectors<pos_t,pos_t>({width_sr/8});
        idx._SR.resize_no_init(idx.z_c+1);

        for (uint16_t i_p=0; i_p<p; i_p++) {
            pos_t b_zc = zc_p[i_p];
            pos_t zc_diff = zc_p[i_p+1]-b_zc;

            for (pos_t i=0; i<zc_diff; i++) {
                idx._SR.template set<0,pos_t>(b_zc+i,SR_file_bufs[i_p][i]);
            }
        }
        
        idx._SR.template set<0,pos_t>(idx.z_c,size_R);

        idx._LP = interleaved_vectors<pos_t,pos_t>({width_lp/8});
        idx._LP.resize_no_init(idx.z_l,p);

        for (uint16_t i_p=0; i_p<p; i_p++) {
            pos_t b_zl = zl_p[i_p];
            pos_t zl_diff = zl_p[i_p+1]-b_zl;

            for (pos_t i=0; i<zl_diff; i++) {
                idx._LP.template set<0,pos_t>(b_zl+i,LP_file_bufs[i_p][i]);
            }
        }

        sdsl::bit_vector PT_tmp;
        PT_tmp.resize(idx.z+2);

        for (uint16_t i_p=0; i_p<p; i_p++) {
            pos_t b_z = z_p[i_p];
            pos_t z_diff = z_p[i_p+1]-b_z;
            
            for (pos_t i=0; i<z_diff; i++) {
                PT_tmp[b_z+i] = PT_file_bufs[i_p][i];
            }
        }

        PT_tmp[idx.z] = 0;
        PT_tmp[idx.z+1] = 0;
        idx._PT = std::move(plain_bit_vector<pos_t,true,true,true>(std::move(PT_tmp)));

        for (uint16_t i=0; i<p; i++) {
            SCP_file_bufs[i].close(true);
            SR_file_bufs[i].close(true);
            LP_file_bufs[i].close(true);
            PT_file_bufs[i].close(true);
        }
    } else {
        idx._SR = std::move(SR[0]);
        idx._SR.resize_no_init(idx.z_c+1,p);

        #pragma omp parallel num_threads(p)
        {
            uint16_t i_p = omp_get_thread_num();

            if (i_p > 0) {
                pos_t b_zc = zc_p[i_p];
                pos_t zc_diff = zc_p[i_p+1]-b_zc;

                for (pos_t i=0; i<zc_diff; i++) {
                    idx._SR.template set<0>(b_zc+i,SR[i_p][i]);
                }
                
                SR[i_p].clear();
                SR[i_p].shrink_to_fit();
            }
        }
        
        idx._SR.template set<0>(idx.z_c,size_R);
        SR.clear();
        SR.shrink_to_fit();

        PT[0].resize(idx.z+2);

        for (uint16_t i_p=1; i_p<p; i_p++) {
            pos_t b_z = z_p[i_p];
            pos_t z_diff = z_p[i_p+1]-b_z;
            
            for (pos_t i=0; i<z_diff; i++) {
                PT[0][b_z+i] = PT[i_p][i];
            }

            PT[i_p] = sdsl::bit_vector();
        }

        PT[0][idx.z] = 0;
        PT[0][idx.z+1] = 0;
        idx._PT = std::move(plain_bit_vector<pos_t,true,true,true>(std::move(PT[0])));
        PT.clear();
        PT.shrink_to_fit();

        idx._LP = std::move(LP[0]);
        idx._LP.resize_no_init(idx.z_l,p);

        #pragma omp parallel num_threads(p)
        {
            uint16_t i_p = omp_get_thread_num();

            if (i_p > 0) {
                pos_t b_zl = zl_p[i_p];
                pos_t zl_diff = zl_p[i_p+1]-b_zl;

                for (pos_t i=0; i<zl_diff; i++) {
                    idx._LP.template set<0>(b_zl+i,LP[i_p][i]);
                }
                
                LP[i_p].clear();
                LP[i_p].shrink_to_fit();
            }
        }

        LP.clear();
        LP.shrink_to_fit();
        
        SCP_b = sdsl::sd_vector_builder(n+2,idx.z_c+2);

        for (uint16_t i_p=0; i_p<p; i_p++) {
            for (pos_t i=0; i<SCP[i_p].size(); i++) {
                SCP_b.set(SCP[i_p][i]);
            }

            SCP[i_p].clear();
            SCP[i_p].shrink_to_fit();
        }

        SCP_b.set(n);
        SCP_b.set(n+1);
        SCP.clear();
        SCP.shrink_to_fit();
        idx._SCP = std::move(sd_array<pos_t>(std::move(sdsl::sd_vector<>(SCP_b))));
    }

    if (log) {
        time = log_runtime(time);
        std::cout << "z: " << idx.z << std::endl;
        std::cout << "z_l/z: " << idx.z_l/(double)idx.z << std::endl;
        std::cout << "z/r': " << idx.z/(double)r_ << std::endl;
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
void move_r<locate_support,sym_t,pos_t>::construction::load_r() {
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

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
template <bool bigbwt, typename sa_sint_t>
void move_r<locate_support,sym_t,pos_t>::construction::build_sas_rlzdsa() {
    if (log) {
        std::cout << "building SA_s" << std::flush;
    }

    idx._SA_s = std::move(interleaved_vectors<pos_t,pos_t>({(uint8_t)std::ceil(std::log2(n+1)/(double)8)}));
    idx._SA_s.resize_no_init(idx.r_);

    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<r_; i++) {
        idx._SA_s.template set<0>(i,SA<bigbwt,sa_sint_t>(omp_get_thread_num(),idx.M_LF().p(i)));
    }

    if constexpr (bigbwt) {
        SA_file_bufs.clear();
        SA_file_bufs.shrink_to_fit();
        std::filesystem::remove(prefix_tmp_files + ".sa");
    }

    if (log) {
        time = log_runtime(time);
    }
}