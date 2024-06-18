template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
template <bool bigbwt, typename sa_sint_t>
void move_r<locate_support,sym_t,pos_t>::construction::build_rlzdsa() {
    if (log) {
        time = now();
        std::cout << "building SA_s" << std::flush;
    }
    
    std::vector<sa_sint_t>& SA = get_sa<sa_sint_t>(); // [0..n-1] The suffix array
    file_vector<uint64_t,pos_t> SA_file;
    std::ifstream SA_ifile;

    if constexpr (bigbwt) {
        p = 1;
        //SA_file = std::move(sdsl::int_vector_buffer<40>(prefix_tmp_files + ".sa",std::ios::in,1024*1024,40,false));
        SA_ifile.open(prefix_tmp_files + ".sa");
        SA_file = std::move(file_vector<uint64_t,pos_t>(SA_ifile,5));

        /*
        interleaved_vectors<uint64_t,pos_t> SA_40({5});
        SA_40.resize_no_init(n);
        std::ifstream SA_ifile;
        SA_ifile.open(prefix_tmp_files + ".sa");
        read_from_file(SA_ifile,SA_40.data(),5*n);
        no_init_resize(SA,n);

        for (uint64_t i=0; i<n; i++) {
            SA[i] = SA_40[i];
        }
        */
    }

    idx._SA_s = std::move(interleaved_vectors<pos_t,pos_t>({(uint8_t)std::ceil(std::log2(n+1)/(double)8)}));
    idx._SA_s.resize_no_init(idx.r_);

    if constexpr (bigbwt) {
        for (uint64_t i=0; i<r_; i++) {
            idx._SA_s.template set<0>(i,SA_file[idx.M_LF().p(i)]);
        }
    } else {
        #pragma omp parallel for num_threads(p)
        for (uint64_t i=0; i<r_; i++) {
            idx._SA_s.template set<0>(i,SA[idx.M_LF().p(i)]);
        }
    }

    if (log) {
        time = log_runtime(time);
        std::cout << "preprocessing SA_d" << std::flush;
    }

    std::function<int64_t(pos_t)> SA_d = [&](pos_t i){
        if constexpr (bigbwt) {
            return i == 0 ? SA_file[0] : (int64_t{SA_file[i]}-int64_t{SA_file[i-1]});
        } else {
            return i == 0 ? SA[0] : (int64_t{SA[i]}-int64_t{SA[i-1]});
        }
    };

    pos_t size_R_target = std::min(std::max<pos_t>(1,n/3),5*r_);
    using sad_map_t = ankerl::unordered_dense::map<int64_t,std::pair<uint32_t,pos_t>>;
    std::vector<sad_map_t> map_SAd_thr(p);

    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<n; i++) {
        int64_t val = SA_d(i);
        auto res = map_SAd_thr[omp_get_thread_num()].find(val);

        if (res == map_SAd_thr[omp_get_thread_num()].end()) {
            map_SAd_thr[omp_get_thread_num()].try_emplace(val,std::make_pair(0,1));
        } else {
            (*res).second.second++;
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
                for (auto p : map_SAd_thr[m.second]) {
                    auto res = map_SAd_thr[m.first].find(p.first);

                    if (res == map_SAd_thr[m.first].end()) {
                        map_SAd_thr[m.first].emplace(p);
                    } else {
                        (*res).second.second += p.second.second;
                    }
                }

                map_SAd_thr[m.second] = sad_map_t();
            }
        }
    }

    sad_map_t map_SAd = std::move(map_SAd_thr[0]);
    map_SAd_thr.clear();
    map_SAd_thr.shrink_to_fit();
    pos_t sigma_SAd = map_SAd.size();
    std::vector<int64_t> unmap_SAd;
    std::vector<pos_t> freq_SAd;
    no_init_resize(unmap_SAd,sigma_SAd);
    no_init_resize(freq_SAd,sigma_SAd);
    pos_t sym_cur = 0;

    for (auto p : map_SAd) {
        unmap_SAd[sym_cur] = p.first;
        freq_SAd[sym_cur] = p.second.second;
        sym_cur++;
    }

    if (sigma_SAd > p) {
        #pragma omp parallel num_threads(p)
        {
            uint16_t i_p = omp_get_thread_num();

            pos_t b = i_p*(sigma_SAd/p);
            pos_t e = i_p == p-1 ? sigma_SAd-1 : ((i_p+1)*(sigma_SAd/p)-1);

            for (pos_t i=b; i<=e; i++) {
                map_SAd[unmap_SAd[i]].first = i;
            }
        }
    } else {
        for (pos_t i=0; i<sigma_SAd; i++) {
            map_SAd[unmap_SAd[i]].first = i;
        }
    }

    unmap_SAd.clear();
    unmap_SAd.shrink_to_fit();

    std::function<int64_t(pos_t)> SA_d_mapped = [&](pos_t i){
        if constexpr (bigbwt) {
            return map_SAd[i == 0 ? SA_file[0] : (int64_t{SA_file[i]}-int64_t{SA_file[i-1]})].first;
        } else {
            return map_SAd[i == 0 ? SA[0] : (int64_t{SA[i]}-int64_t{SA[i-1]})].first;
        }
    };

    if (log) {
        time = log_runtime(time);
        std::cout << "sigma_SA^d: " << sigma_SAd << std::endl;
        std::cout << "building R" << std::flush;
    }

    idx._R = std::move(interleaved_vectors<uint64_t,pos_t>({(uint8_t)std::ceil(std::log2(2*n+1)/(double)8)}));
    sdsl::bit_vector occurs_in_R;
    occurs_in_R.resize(sigma_SAd);
    pos_t s = std::min<pos_t>(2048,size_R_target);
    float score_cutoff = 50.0/(float)s;
    uint64_t num_cand_str = 17*powf(n/(float)r,0.4);
    std::uniform_int_distribution<uint32_t> pos_distrib(0,n-1-s);
    absl::btree_set<std::pair<pos_t,pos_t>> segments;
    pos_t size_R_current = 0;
    //std::geometric_distribution<pos_t> geo_distrib(6.0/s);

    while (size_R_current < size_R_target) {
        pos_t best_pos = 0;
        pos_t best_len = 0;
        float best_score = 0;
        typename absl::btree_set<std::pair<pos_t,pos_t>>::iterator best_it = segments.end();

        #pragma omp parallel for num_threads(p)
        for (uint64_t str=0; str<num_cand_str; str++) {
            std::random_device rd;
            std::mt19937 gen(rd());
            pos_t pos = pos_distrib(gen);
            pos_t len = s;//-std::min<pos_t>(s-1,geo_distrib(gen));
            float score = 0;
            ankerl::unordered_dense::set<uint32_t> processed_values;
            typename absl::btree_set<std::pair<pos_t,pos_t>>::iterator it_thr = segments.end();

            if (!segments.empty()) {
                it_thr = segments.upper_bound(std::make_pair(pos,0));
                
                if (it_thr != segments.end()) {
                    pos_t start_next = (*it_thr).first;

                    if (start_next < pos+len) {
                        len = start_next-pos;
                    }
                }

                if (it_thr != segments.begin()) {
                    auto it_thr_before = it_thr;
                    it_thr_before--;
                    pos_t end_last = (*it_thr_before).first+(*it_thr_before).second;

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
                uint32_t val = SA_d_mapped(pos+i);

                if (occurs_in_R[val] == 0 && !processed_values.contains(val)) {
                    score += std::sqrt(freq_SAd[val]);
                    processed_values.emplace(val);
                }
            }

            score /= len;

            #pragma omp critical
            {
                if (score > best_score) {
                    best_pos = pos;
                    best_len = len;
                    best_score = score;
                    best_it = it_thr;
                }
            }
        }

        segments.emplace_hint(best_it,std::make_pair(best_pos,best_len));

        for (pos_t i=0; i<best_len; i++) {
            occurs_in_R[SA_d_mapped(best_pos+i)] = 1;
        }

        size_R_current += best_len;

        if (best_score < score_cutoff) {
            break;
        }
    }

    for (auto segment : segments) {
        for (pos_t i=0; i<segment.second; i++) {
            idx._R.emplace_back(std::make_tuple(SA_d(segment.first+i)+n));
        }
    }
    
    pos_t size_R = idx._R.size();

    /*
    // blockwise random reference generation

    pos_t num_blocks = 128;
    pos_t size_R = size_R_target;
    idx._R = std::move(interleaved_vectors<uint64_t,pos_t>({(uint8_t)std::ceil(std::log2(2*n+1)/(double)8)}));
    idx._R.resize_no_init(size_R);
    pos_t block_size = std::max(pos_t{1},size_R/num_blocks);
    num_blocks = std::ceil(size_R/(double)block_size);

    for (pos_t cur_block=0; cur_block<num_blocks; cur_block++) {
        pos_t cur_block_size = cur_block == num_blocks - 1 ? size_R - (num_blocks - 1) * block_size : block_size;
        pos_t b_R = cur_block * block_size;
        pos_t b_SAd = cur_block * (n/num_blocks);

        for (pos_t i=0; i<cur_block_size; i++) {
            idx._R.template set<0,pos_t>(b_R+i,SA_d(b_SAd+i)+n);
        }
    }

    */

    if (log) {
        time = log_runtime(time);
        std::cout << "building move-r of R^R" << std::flush;
    }

    std::vector<uint64_t> R_rev;
    R_rev.reserve(size_R);

    for (int64_t i=size_R-1; i>=0; i--) {
        R_rev.emplace_back(idx._R[i]);
    }
    
    move_r<_phi,uint64_t,uint32_t> idx_R_rev(std::move(R_rev),{.num_threads=p});

    if (log) {
        time = log_runtime(time);
        std::cout << "bulding rlzdsa" << std::flush;
    }

    interleaved_vectors<pos_t,pos_t> SP({(uint8_t)std::ceil(std::log2(n+1)/(double)8)});
    idx._SR = std::move(interleaved_vectors<pos_t,pos_t>({(uint8_t)std::ceil(std::log2(size_R+1)/(double)8)}));
    idx._LP = std::move(interleaved_vectors<pos_t,pos_t>({(uint8_t)std::ceil(std::log2(n+1)/(double)8)}));
    sdsl::bit_vector PT;
    PT.resize(1);
    PT[0] = 1;
    idx.z = 1;

    if constexpr (bigbwt) {
        idx._LP.emplace_back(std::make_tuple(SA_file[0]));
    } else {
        idx._LP.emplace_back(std::make_tuple(SA[0]));
    }

    for (pos_t i=1; i<n;) {
        auto qc = idx_R_rev.query();
        auto last_qc = qc;
        idx.z++;
        pos_t l = 0;

        do {
            last_qc = qc;
            qc.prepend(SA_d(i+l)+n);

            if (qc.num_occ() == 0) {
                break;
            }

            l++;

            if (i+l == n) {
                last_qc = qc;
                break;
            }
        } while (true);

        if (l <= 1) {
            PT.resize(idx.z);
            PT[idx.z-1] = 1;

            if constexpr (bigbwt) {
                idx._LP.emplace_back(std::make_tuple(SA_file[i]));
            } else {
                idx._LP.emplace_back(std::make_tuple(SA[i]));
            }

            i++;
        } else {
            SP.emplace_back(std::make_tuple(i));
            PT.resize(idx.z);
            PT[idx.z-1] = 0;
            idx._SR.emplace_back(std::make_tuple(size_R-last_qc.next_occ()-l));
            i += l;
        }
    }

    idx._SR.emplace_back(std::make_tuple(size_R));
    PT.resize(idx.z+1);
    PT[idx.z] = 0;
    idx._PT = std::move(plain_bit_vector<pos_t,true,true,true>(std::move(PT)));
    idx.z_c = idx._PT.num_zeros()-1;
    idx.z_l = idx._PT.num_ones();
    sdsl::sd_vector_builder SCP_builder(n+1,idx.z_c+1);

    for (pos_t i=0; i<idx.z_c; i++) {
        SCP_builder.set(SP[i]);
    }

    SCP_builder.set(n);
    idx._SCP = std::move(sd_array<pos_t>(std::move(sdsl::sd_vector<>(SCP_builder))));

    if (log) {
        time = log_runtime(time);
        std::cout << "z: " << idx.z << std::endl;
        std::cout << "z_l/z: " << idx.z_l/(double)idx.z << std::endl;
        std::cout << "z/r': " << idx.z/(double)r_ << std::endl;
    }

    /*
    std::cout << "SA: ";
    for (pos_t i=0; i<n; i++) {
        std::cout << SA[i] << ", ";
    }

    std::cout << std::endl << "SA^d: ";
    for (pos_t i=0; i<n; i++) {
        std::cout << SA_d(i) << ", ";
    }

    std::cout << std::endl << "R: ";
    for (pos_t i=0; i<size_R; i++) {
        std::cout << int64_t{idx._R[i]}-int64_t{n} << ", ";
    }

    std::cout << std::endl << "SR: ";
    for (pos_t i=0; i<idx._SR.size(); i++) {
        std::cout << idx._SR[i] << ", ";
    }

    std::cout << std::endl << "SP: ";
    for (pos_t i=0; i<SP.size(); i++) {
        std::cout << SP[i] << ", ";
    }
    std::cout << n << ", ";

    std::cout << std::endl << "PT: ";
    for (pos_t i=0; i<idx._PT.size(); i++) {
        std::cout << std::to_string(idx._PT[i]) << ", ";
    }

    std::cout << std::endl << "LP: ";
    for (pos_t i=0; i<idx._LP.size(); i++) {
        std::cout << std::to_string(idx._LP[i]) << ", ";
    }

    std::cout << std::endl;
    */
}