#include <filesystem>

template <typename uint_t>
void move_r<uint_t>::construction::preprocess_and_store_t_in_file() {
    preprocess_t(true);
    
    prefix_tmp_files = "move-r_" + random_alphanumeric_string(10);
    std::ofstream T_ofile(prefix_tmp_files);
    write_to_file(T_ofile,T.c_str(),n-1);
    T_ofile.close();
}

template <typename uint_t>
void move_r<uint_t>::construction::pfp() {
    if (log) {
        time = now();
        std::cout << "executing Big-BWT" << std::endl << std::endl;
    }

    system((
        "external/Big-BWT/bigbwt " +
        (std::string)(build_locate_support ? "-s -e " : "") +
        (std::string)(p > 1 ? ("-t " + std::to_string(p) + " ") : "") +
        prefix_tmp_files +
        (std::string)(log ? "" : " >log_1 >log_2")
    ).c_str());

    std::ifstream log_ifile(prefix_tmp_files + ".log");
    bigbwt_peak_mem_usage = (malloc_count_current()-baseline_mem_usage)+malloc_count_peak_memory_usage(log_ifile);

#ifdef MOVE_R_BENCH
    external_peak_memory_usage = bigbwt_peak_mem_usage;
#endif

    log_ifile.close();
    std::filesystem::remove(prefix_tmp_files + ".log");
    std::filesystem::remove(prefix_tmp_files);

    if (log) {
        std::filesystem::remove("log_1");
        std::filesystem::remove("log_2");
    }

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_pfp=" << time_diff_ns(time,now());
        std::cout << std::endl;
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::build_rlbwt_c_space_saving() {
    if (log) {
        time = now();
        std::cout << "building RLBWT" << std::flush;
    }

    std::ifstream bwt_file(prefix_tmp_files + ".bwt");
    RLBWT.resize(p,std::move(interleaved_vectors<uint32_t>({1,4})));

    for (uint16_t i=0; i<p; i++) {
        n_p.emplace_back(i*(n/p));
    }

    n_p.emplace_back(n);
    r_p.resize(p+1,0);
    r_p[0] = 0;
    C.resize(p+1,std::vector<uint_t>(256,0));

    char cur_char;
    char prev_char;
    bwt_file.read(&prev_char,1);

    if (prev_char == uchar_to_char(2)) {
        prev_char = uchar_to_char(0);
    }

    uint_t i_ = 0;
    uint16_t i_p = 0;
    uint16_t ip_nxt = 1;
    uint_t np_nxt = n_p[ip_nxt];
    uint_t cur_bwt_buf_size;
    uint_t i = 1;
    std::string bwt_buf;
    uint_t max_bwt_buf_size = std::max((uint_t)1,n/500);
    no_init_resize(bwt_buf,max_bwt_buf_size);
    uint_t i_buf;

    while (i < n) {
        cur_bwt_buf_size = std::min(n-i,max_bwt_buf_size);
        read_from_file(bwt_file,bwt_buf.c_str(),cur_bwt_buf_size);
        i_buf = 0;
        
        while (i_buf < cur_bwt_buf_size) {
            cur_char = bwt_buf[i_buf] == uchar_to_char(2) ? uchar_to_char(0) : bwt_buf[i_buf];

            if (cur_char != prev_char) {
                C[i_p][char_to_uchar(prev_char)] += i-i_;
                RLBWT[0].emplace_back_unsafe<char,uint32_t>({prev_char,i-i_});
                prev_char = cur_char;
                i_ = i;

                if (i >= np_nxt) {
                    n_p[ip_nxt] = i;
                    r_p[ip_nxt] = RLBWT[0].size();
                    i_p++;
                    ip_nxt++;
                    np_nxt = n_p[ip_nxt];
                }
            }

            i++;
            i_buf++;
        }
    }

    RLBWT[0].emplace_back_unsafe<char,uint32_t>({prev_char,n-i_});
    C[i_p][char_to_uchar(prev_char)] += n-i_;
    RLBWT[0].shrink_to_fit();
    bwt_buf.clear();
    bwt_buf.shrink_to_fit();
    r = RLBWT[0].size();
    idx.r = r;
    r_p[i_p+1] = r;
    n_p[i_p+1] = n;
    i_p += 2;

    while (i_p <= p) {
        r_p[i_p] = r;
        n_p[i_p] = n;
        i_p++;
    }

    bwt_file.close();
    std::filesystem::remove(prefix_tmp_files + ".bwt");

    for (uint16_t i_p=1; i_p<p; i_p++) {
        RLBWT[i_p].set_data(RLBWT[0].data()+5*r_p[i_p],r_p[i_p+1]-r_p[i_p]);
    }

    process_c();

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_build_rlbwt=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::read_iphi_space_saving() {
    if (log) {
        time = now();
        std::cout << "reading I_Phi" << std::flush;
    }

    no_init_resize(I_Phi,r);
    std::ifstream ssa_file(prefix_tmp_files + ".ssa");
    std::ifstream esa_file(prefix_tmp_files + ".esa");
    
    interleaved_vectors<uint64_t> ssa({5});
    interleaved_vectors<uint64_t> esa({5});
    ssa.resize_no_init(r);
    esa.resize_no_init(r);
    read_from_file(ssa_file,ssa.data(),5*r);
    read_from_file(esa_file,esa.data(),5*r);

    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<r-1; i++) {
        if constexpr (std::is_same<uint_t,uint32_t>::value) {
            I_Phi[i].first = ssa.get_unsafe<0,uint32_t>(i);
            I_Phi[i+1].second = esa.get_unsafe<0,uint32_t>(i);
        } else {
            I_Phi[i].first = ssa[i];
            I_Phi[i+1].second = esa[i];
        }
    }

    I_Phi[r-1].first = ssa[r-1];
    I_Phi[0].second = esa[r-1];

    ssa.clear();
    esa.clear();
    ssa.shrink_to_fit();
    esa.shrink_to_fit();

    ssa_file.close();
    esa_file.close();
    std::filesystem::remove(prefix_tmp_files + ".ssa");
    std::filesystem::remove(prefix_tmp_files + ".esa");
    
    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_read_iphi=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::store_mlf() {
    if (log) {
        time = now();
        std::cout << "storing M_LF in a file" << std::flush;
    }

    std::ofstream file_mlf(prefix_tmp_files + ".mlf");
    idx.M_LF.serialize(file_mlf);
    idx.M_LF = std::move(move_data_structure_str<uint_t>());
    file_mlf.close();

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_store_mlf=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::load_mlf() {
    if (log) {
        time = now();
        std::cout << "loading M_LF from a file" << std::flush;
    }

    std::ifstream file_mlf(prefix_tmp_files + ".mlf");
    idx.M_LF.load(file_mlf);
    file_mlf.close();
    std::filesystem::remove(prefix_tmp_files + ".mlf");

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_load_mlf=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}