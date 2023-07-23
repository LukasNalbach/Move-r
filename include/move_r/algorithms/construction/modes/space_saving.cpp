#include <filesystem>
#include <move_r/algorithms/construction/Big-BWT/newscan.cpp>
#include <move_r/algorithms/construction/Big-BWT/pfbwt.cpp>

extern "C" {
    #include <move_r/algorithms/construction/Big-BWT/bwtparse.c>
}

template <typename uint_t>
std::ifstream move_r<uint_t>::construction::preprocess_and_store_t_in_file() {
    preprocess_t(true,false);
    
    prefix_tempfiles = "move-r_" + random_alphanumeric_string(10);
    std::ofstream T_output_file(prefix_tempfiles);
    write_to_file(T_output_file,T.c_str(),n-1);
    T_output_file.close();
    std::ifstream T_input_file(prefix_tempfiles);

    return T_input_file;
}

template <typename uint_t>
void move_r<uint_t>::construction::preprocess_t_buffered_from_file(std::ifstream& t_file) {
    t_file.seekg(0,std::ios::end);
    n = t_file.tellg()+(std::streamsize)+1;
    idx.n = n;
    t_file.seekg(0,std::ios::beg);
    preprocess_t(false,false,&t_file);
    t_file.seekg(0,std::ios::beg);
}

template <typename uint_t>
void move_r<uint_t>::construction::pfp(std::ifstream& t_file, bool delete_t_file) {
    if (log) {
        time = now();
        std::cout << "building PFP of T" << std::flush;
    }

    if (!delete_t_file) {
        prefix_tempfiles = "move-r_" + random_alphanumeric_string(10);
    }

    uint64_t dictionary_size;

    bigbwt_newscan({
        .input=t_file,
        .chars_remapped=idx.chars_remapped,
        .map_char=idx.map_char,
        .name_inputfile=prefix_tempfiles,
        .compute_sa_info=true
    },dictionary_size,false);

    if (delete_t_file) {
        t_file.close();
        std::filesystem::remove(prefix_tempfiles);
    }

    if (log) {
        if (measurement_file_index != NULL) *measurement_file_index << " time_pfp_phase_1=" << time_diff_ns(time,now());
        time = log_runtime(time);
        std::cout << "building BWT of PFP of T" << std::flush;
    }

    bigbwt_bwtparse({
        .name_inputfile=prefix_tempfiles.c_str(),
        .compute_sa_info=true
    },false);

    if (log) {
        if (measurement_file_index != NULL) *measurement_file_index << " time_pfp_phase_2=" << time_diff_ns(time,now());
        time = log_runtime(time);
        std::string msg = "building RLBWT";
        if (build_locate_support) {
            msg.append(" and I_Phi");
        }
        std::cout << msg << " of T" << std::flush;
    }

    uint64_t r = 0;

    bigbwt_pfbwt({
        .name_inputfile=prefix_tempfiles.c_str(),
        .n=n,
        .build_sa_samples=build_locate_support?3:0
    },r,false);

    this->r = r;
    idx.r = r;

    std::vector<std::string> suffixes = {".bwlast",".bwsai",".dict",".ilist",".last",".occ",".parse",".parse_old",".sai"};

    for (std::string& suffix : suffixes) {
        std::filesystem::remove(prefix_tempfiles + suffix);
    }

    if (log) {
        double n_r = std::round(100.0*(n/(double)r))/100.0;
        if (measurement_file_index != NULL) {
            *measurement_file_index << " time_pfp_phase_3=" << time_diff_ns(time,now());
            *measurement_file_index << " n=" << n;
            *measurement_file_index << " sigma=" << std::to_string(sigma);
            *measurement_file_index << " r=" << r;
        }
        time = log_runtime(time);
        std::cout << "n = " << n << ", sigma = " << std::to_string(sigma) << ", r = " << r << ", n/r = " << n_r << std::endl;
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::build_rlbwt_and_c_space_saving() {
    read_rle_bwt = std::filesystem::exists(prefix_tempfiles + ".bwtr");

    if (log) {
        time = now();
        std::string msg;
        if (read_rle_bwt) {
            msg = "reading RLBWT";
        } else {
            msg = "reading BWT";
        }
        std::cout << msg << std::flush;
    }
    
    no_init_resize(bwt_run_heads,r);

    r_p.resize(p+1);
    r_p[0] = 0;
    r_p[p] = r;

    for (uint16_t i=1; i<p; i++) {
        r_p[i] = i*(r/p);
    }
    
    C.resize(p+1,std::vector<uint_t>(256,0));

    if (read_rle_bwt) {
        std::ifstream bwt_run_heads_file(prefix_tempfiles + ".bwtr");
        read_from_file(bwt_run_heads_file,bwt_run_heads.c_str(),r);
        bwt_run_heads_file.close();
        std::filesystem::remove(prefix_tempfiles + ".bwtr");

        std::ifstream bwt_run_lengths_file(prefix_tempfiles + ".bwtrls");
        (*reinterpret_cast<std::vector<no_init<uint32_t>>*>(&bwt_run_lengths)).resize(r);
        read_from_file(bwt_run_lengths_file,(char*)&bwt_run_lengths[0],r*sizeof(uint32_t));
        bwt_run_lengths_file.close();
        std::filesystem::remove(prefix_tempfiles + ".bwtrls");

        n_p.resize(p+1,0);

        #pragma omp parallel num_threads(p)
        {
            uint16_t i_p = omp_get_thread_num();

            uint_t b = r_p[i_p];
            uint_t e = r_p[i_p+1];

            for (uint_t i=b; i<e; i++) {
                n_p[i_p] += bwt_run_lengths[i];
                C[i_p][char_to_uchar(bwt_run_heads[i])] += bwt_run_lengths[i];
            }
        }

        for (uint16_t i=1; i<p; i++) {
            n_p[i] += n_p[i-1];
        }

        for (int16_t i=p; i>0; i--) {
            n_p[i] = n_p[i-1];
        }

        n_p[0] = 0;

    } else {
        std::ifstream bwt_file(prefix_tempfiles + ".bwt");
        (*reinterpret_cast<std::vector<no_init<uint32_t>>*>(&bwt_run_lengths)).resize(r);
        
        char last_char = bwt_file.get();
        char cur_char = bwt_file.get();
        uint32_t cur_run_index;
        uint32_t cur_run_length;
        uint_t cur_pos_in_L;
        uint16_t cur_threads_section = 0;
        uint16_t next_threads_section = 1;

        if (cur_char == last_char) {
            cur_run_index = 0;
            cur_run_length = 2;
            cur_pos_in_L = 0;
        } else {
            bwt_run_heads[0] = last_char;
            bwt_run_lengths[0] = 1;
            C[0][char_to_uchar(last_char)] = 1;
            last_char = cur_char;
            cur_run_index = 1;
            cur_run_length = 1;
            cur_pos_in_L = 1;
        }

        n_p.emplace_back(0);

        while (cur_run_index < r) {
            while ((cur_char = bwt_file.get()) == last_char) {
                cur_run_length++;
            }

            bwt_run_heads[cur_run_index] = last_char;
            bwt_run_lengths[cur_run_index] = cur_run_length;
            C[cur_threads_section][char_to_uchar(last_char)] += cur_run_length;

            cur_run_index++;
            cur_pos_in_L += cur_run_length;
            last_char = cur_char;
            cur_run_length = 1;

            if (cur_run_index == r_p[next_threads_section]) {
                n_p.emplace_back(cur_pos_in_L);
                cur_threads_section++;
                next_threads_section++;
            }
        }

        bwt_file.close();
        std::filesystem::remove(prefix_tempfiles + ".bwt");
    }

    if (log) {
        if (measurement_file_index != NULL) *measurement_file_index << " time_read_rlbwt=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::build_ilf_space_saving() {
    if (log) {
        time = now();
        std::cout << "building I_LF" << std::flush;
    }

    (*reinterpret_cast<std::vector<std::pair<no_init<uint_t>,no_init<uint_t>>>*>(&I_LF)).resize(r);

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Build I_LF[r_p[i_p]..r_p[i_p+1]-1]

        // Iteration range start position of thread i_p.
        uint_t b = r_p[i_p];
        // Iteration range end position of thread i_p.
        uint_t e = r_p[i_p+1];

        // Start position of the current run.
        uint_t i_ = n_p[i_p];

        for (uint_t i=b; i<e; i++) {
            /* Write the pair (i_,LF(i_)) to the current position i in I_LF, where
            LF(i_) = C[L[i_]] + rank(L,L[i_],i_-1) = C[p][L[i_]] + C[i_p][L[i_]]. */
            I_LF[i] = std::make_pair(i_,C[p][char_to_uchar(bwt_run_heads[i])]+C[i_p][char_to_uchar(bwt_run_heads[i])]);

            // Update the start position i_ of the current run.
            i_ += bwt_run_lengths[i];

            /* Update the rank-function in C[i_p] to store C[i_p][c] = rank(L,c,i_-1),
            for each c in [0..255]. */
            C[i_p][char_to_uchar(bwt_run_heads[i])] += bwt_run_lengths[i];
        }
    }

    C.clear();
    C.shrink_to_fit();
    
    if (log) {
        if (measurement_file_index != NULL) *measurement_file_index << " build_ilf=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::read_iphi_space_saving() {
    if (log) {
        time = now();
        std::cout << "reading I_Phi" << std::flush;
    }

    (*reinterpret_cast<std::vector<std::pair<no_init<uint_t>,no_init<uint_t>>>*>(&I_Phi)).resize(r);
    std::ifstream I_Phi_file(prefix_tempfiles + ".iphi");
    
    if constexpr (std::is_same<uint_t,uint32_t>::value) {
        I_Phi_file.read((char*)&I_Phi[1],(r-1)*2*sizeof(uint32_t));
        I_Phi_file.read((char*)&I_Phi[0],2*sizeof(uint32_t));
    } else {
        #pragma omp parallel for num_threads(p)
        for (uint64_t i=0; i<r; i++) {
            std::memset(reinterpret_cast<uint8_t*>(&I_Phi[i].first)+5,0,3);
            std::memset(reinterpret_cast<uint8_t*>(&I_Phi[i].second)+5,0,3);
        }

        for (uint64_t i=1; i<r; i++) {
            I_Phi_file.read((char*)&I_Phi[i].first,5);
            I_Phi_file.read((char*)&I_Phi[i].second,5);
        }
    
        I_Phi_file.read((char*)&I_Phi[0].first,5);
        I_Phi_file.read((char*)&I_Phi[0].second,5);
    }

    I_Phi_file.close();
    std::filesystem::remove(prefix_tempfiles + ".iphi");
    
    if (log) {
        if (measurement_file_index != NULL) *measurement_file_index << " time_read_phi=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::build_l__and_sas_space_saving() {
    if (log) {
        time = now();
        std::string msg;
        if (build_locate_support) {
            msg = "building L' and SA_s";
        } else {
            msg = "building L'";
        }
        std::cout << msg << std::flush;
    }

    if (build_locate_support) {
        (*reinterpret_cast<std::vector<no_init<uint_t>>*>(&SA_s)).resize(r_);
        SA_s[r_-1] = I_Phi[0].second;

        SA_s_missing_thr.resize(p);
    }

    // Simultaneously iterate over the input intervals of M_LF nad the bwt runs to build SA_s and L'
    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Bwt range start position of thread i_p.
        uint_t b = n_p[i_p];

        // Bwt runs iteration range start position of thread i_p.
        uint_t b_r = r_p[i_p];
        // Bwt runs iteration range end position of thread i_p.
        uint_t e_r = r_p[i_p+1];

        // Index of the current input interval in M_LF, initially the index of the input interval of M_LF, in which b lies.
        uint_t j;
        {
            uint_t x,y,z;
            x = 0;
            z = r_-1;

            while (x != z) {
                y = (x+z)/2+1;
                if (idx.M_LF.p(y) <= b) {
                    x = y;
                } else {
                    z = y-1;
                }
            }

            j = x;
        }

        uint_t j_; // Index of the first input interval in M_LF starting in the current bwt run.
        uint_t l_ = b; // Starting position of the next bwt run.

        for (uint_t i=b_r; i<e_r; i++) {
            j_ = j;

            // update l_ to the next run start position
            l_ += bwt_run_lengths[i];

            // iterate over all input intervals in M_LF within the i-th bwt run, that have been created by the balancing algorithm
            do {
                idx.M_LF.template set_character<char>(j,bwt_run_heads[i]);
                j++;
            } while (idx.M_LF.p(j) < l_);

            if (build_locate_support) {
                if (j != r_) {
                    SA_s[j-1] = I_Phi[i+1].second;
                }

                // check if the i-th bwt run has been split during the balancing algorithm
                if (j-j_ > 1) {
                    // mark, that SA_s[j_..j-2] has yet to be computed with l_-M_LF.p(j_+1) Phi-move queries
                    SA_s_missing_thr[i_p].emplace_back(std::make_tuple(j_,j-1,l_-idx.M_LF.p(j_+1)));

                    for (uint_t x=j_; x<=j-2; x++) {
                        SA_s[x] = n;
                    }
                }
            }
        }
    }

    if (!build_locate_support && log) {
        if (measurement_file_index != NULL) *measurement_file_index << " time_build_l_=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
        
    if (build_locate_support) {
        SA_s_missing_sect.resize(p+1,0);
        num_SA_s_missing_thr.resize(p+1,0);

        for (uint16_t i=0; i<p; i++) {
            num_SA_s_missing_thr[i] = SA_s_missing_thr[i].size();
        }

        for (uint16_t i=1; i<p; i++) {
            num_SA_s_missing_thr[i] += num_SA_s_missing_thr[i-1];
        }

        for (uint16_t i=p; i>0; i--) {
            num_SA_s_missing_thr[i] = num_SA_s_missing_thr[i-1];
        }

        num_SA_s_missing_thr[0] = 0;
        SA_s_missing_sect[p] = num_SA_s_missing_thr[p];

        (*reinterpret_cast<std::vector<std::tuple<no_init<uint_t>,no_init<uint_t>,no_init<uint_t>>>*>(&SA_s_missing)).resize(num_SA_s_missing_thr[p]);
        
        #pragma omp parallel num_threads(p)
        {
            // Index in [0..p-1] of the current thread.
            uint16_t i_p = omp_get_thread_num();

            uint_t b = num_SA_s_missing_thr[i_p];
            uint_t e = num_SA_s_missing_thr[i_p+1];

            for (uint_t i=b; i<e; i++) {
                SA_s_missing[i] = SA_s_missing_thr[i_p][i-b];
            }
        }

        SA_s_missing_thr.clear();
        SA_s_missing_thr.shrink_to_fit();

        for (uint_t i=1; i<num_SA_s_missing_thr[p]; i++) {
            std::get<2>(SA_s_missing[i]) += std::get<2>(SA_s_missing[i-1]);
        }

        if (num_SA_s_missing_thr[p] != 0) {
            #pragma omp parallel num_threads(p)
            {
                // Index in [0..p-1] of the current thread.
                uint16_t i_p = omp_get_thread_num();

                uint_t x,y,z;
                x = 0;
                z = num_SA_s_missing_thr[p]-1;
                uint_t opt = i_p*(std::get<2>(SA_s_missing[num_SA_s_missing_thr[p]-1])/p);

                while (x != z) {
                    y = (x+z)/2+1;
                    if (std::get<2>(SA_s_missing[y]) <= opt) {
                        x = y;
                    } else {
                        z = y-1;
                    }
                }

                SA_s_missing_sect[i_p] = x;
            }
        }

        n_p.clear();
        n_p.shrink_to_fit();

        r_p.clear();
        r_p.shrink_to_fit();

        bwt_run_heads.clear();
        bwt_run_heads.shrink_to_fit();

        bwt_run_lengths.clear();
        bwt_run_lengths.shrink_to_fit();

        num_SA_s_missing_thr.clear();
        num_SA_s_missing_thr.shrink_to_fit();

        if (log) {
            if (measurement_file_index != NULL) *measurement_file_index << " time_build_l_sas=" << time_diff_ns(time,now());
            time = log_runtime(time);
        }
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::store_mlf() {
    if (log) {
        time = now();
        std::cout << "storing M_LF in a file" << std::flush;
    }

    std::ofstream file_mlf(prefix_tempfiles + ".mlf");
    idx.M_LF.serialize(file_mlf);
    idx.M_LF = std::move(move_data_structure_lf<uint_t>());
    file_mlf.close();

    if (log) {
        if (measurement_file_index != NULL) *measurement_file_index << " time_store_mlf=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::load_mlf() {
    if (log) {
        time = now();
        std::cout << "loading M_LF from a file" << std::flush;
    }

    std::ifstream file_mlf(prefix_tempfiles + ".mlf");
    idx.M_LF.load(file_mlf);
    file_mlf.close();
    std::filesystem::remove(prefix_tempfiles + ".mlf");

    if (log) {
        if (measurement_file_index != NULL) *measurement_file_index << " time_load_mlf=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::compute_missing_sa_samples_space_saving() {
    if (log) {
        time = now();
        std::cout << "computing the missing SA-samples" << std::flush;
    }

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        if (SA_s_missing_sect[i_p] < SA_s_missing_sect[i_p+1]) {
            int64_t b = SA_s_missing_sect[i_p];
            int64_t e = SA_s_missing_sect[i_p+1]-1;

            uint_t j; // current position in SA
            uint_t l_x; // start position of the current input interval in M_LF, that contains j

            uint_t i_s; // current SA-value at position j
            uint_t x_s; // Index in M_phi of i_s

            uint_t x_l; // index of the interval in M_LF to iterate to
            uint_t x; // index of the interval in M_LF to start iterating from

            for (int64_t i=e; i>=b; i--) {
                x_l = std::get<0>(SA_s_missing[i]);
                x = std::get<1>(SA_s_missing[i]);

                j = idx.M_LF.p(x+1)-1;

                x_s = idx.SA_idx(x);
                i_s = idx.M_Phi.p(x_s)+idx.SA_offs(x);

                while (x > x_l) {
                    l_x = idx.M_LF.p(x);

                    while (j >= l_x) {
                        idx.M_Phi.move(i_s,x_s);
                        j--;
                    }

                    x--;
                    idx.set_SA_idx(x,x_s);
                    idx.set_SA_offs(x,i_s-idx.M_Phi.p(x_s));
                }
            }
        }
    }

    if (log) {
        if (measurement_file_index != NULL) *measurement_file_index << " time_compute_missing_sa_samples=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}