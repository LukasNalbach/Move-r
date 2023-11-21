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

    uint64_t num_dictionary_words;

    bigbwt_newscan({
        .input=t_file,
        .chars_remapped=idx.chars_remapped,
        .map_char=idx.map_char,
        .name_input_file=prefix_tempfiles,
        .compute_sa_info=true
    },num_dictionary_words,false);

    read_rlbwt = n > 100 * num_dictionary_words;

    if (delete_t_file) {
        t_file.close();
        std::filesystem::remove(prefix_tempfiles);
    }

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_pfp_phase_1=" << time_diff_ns(time,now());
        time = log_runtime(time);
        std::cout << "building BWT of PFP of T" << std::flush;
    }

    bigbwt_bwtparse({
        .name_input_file=prefix_tempfiles.c_str(),
        .compute_sa_info=true
    },false);

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_pfp_phase_2=" << time_diff_ns(time,now());
        time = log_runtime(time);
        std::string msg = "building " + (std::string)(read_rlbwt ? "RLBWT" : "BWT");
        if (build_locate_support) {
            msg.append(" and I_Phi");
        }
        std::cout << msg << " of T" << std::flush;
    }

    uint64_t r = 0;

    bigbwt_pfbwt({
        .name_input_file=prefix_tempfiles.c_str(),
        .n=n,
        .build_sa_samples=build_locate_support?3:0
    },r,false);

    this->r = r;
    idx.r = r;

    if (p*1000 > r) {
        p = std::max((uint64_t)1,r/1000);
    }

    std::vector<std::string> suffixes = {".bwlast",".bwsai",".dict",".ilist",".last",".occ",".parse",".parse_old",".sai"};

    for (std::string& suffix : suffixes) {
        std::filesystem::remove(prefix_tempfiles + suffix);
    }

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_pfp_phase_3=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::read_rlbwt_bwt() {
    read_rlbwt = std::filesystem::exists(prefix_tempfiles + ".rlbwt");

    if (log) {
        time = now();
        std::string msg;
        if (read_rlbwt) {
            msg = "reading RLBWT";
        } else {
            msg = "reading BWT";
        }
        std::cout << msg << std::flush;
    }

    std::string bwt_file_name = prefix_tempfiles + (read_rlbwt ? ".rlbwt" : ".bwt");
    std::ifstream bwt_file(bwt_file_name);

    if (read_rlbwt) {
        RLBWT = std::move(interleaved_vectors<uint32_t>({1,4},r,false));
        read_from_file(bwt_file,RLBWT.get_data(),5*r);
    } else {
        no_init_resize(L,n);
        read_from_file(bwt_file,&L[0],n);
    }

    bwt_file.close();
    std::filesystem::remove(bwt_file_name);

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_read_rlbwt_bwt=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::preprocess_rlbwt_space_saving() {
    if (log) {
        time = now();
        std::cout << "preprocessing RLBWT" << std::flush;
    }
    
    C.resize(p+1,std::vector<uint_t>(256));
    n_p.resize(p+1);

    for (uint16_t i=0; i<p; i++) {
        r_p.emplace_back(i*(r/p));
    }

    r_p.emplace_back(r);

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Iteration range start position of thread i_p.
        uint_t b_r = r_p[i_p];
        // Iteration range end position of thread i_p.
        uint_t e_r = r_p[i_p+1];

        for (uint_t i=b_r; i<e_r; i++) {
            n_p[i_p] += run_length(i);
            C[i_p][run_uchar(i)] += run_length(i);
        }
    }

    for (uint16_t i=1; i<p; i++) {
        n_p[i] += n_p[i-1];
    }

    for (uint16_t i=p; i>0; i--) {
        n_p[i] = n_p[i-1];
    }

    n_p[0] = 0;
    
    process_c();

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_preprocess_rlbwt=" << time_diff_ns(time,now());
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
    std::ifstream I_Phi_file(prefix_tempfiles + ".iphi");
    
    if constexpr (std::is_same<uint_t,uint32_t>::value) {
        read_from_file(I_Phi_file,(char*)&I_Phi[1],8*(r-1));
        read_from_file(I_Phi_file,(char*)&I_Phi[0],8);
    } else {
        interleaved_vectors<uint64_t> I_Phi_tmp({5,5},r,false);
        read_from_file(I_Phi_file,I_Phi_tmp.get_data()+10,10*(r-1));
        read_from_file(I_Phi_file,I_Phi_tmp.get_data(),10);

        #pragma omp parallel for num_threads(p)
        for (uint64_t i=0; i<r; i++) {
            I_Phi[i].first = I_Phi_tmp.get<0>(i);
            I_Phi[i].second = I_Phi_tmp.get<1>(i);
        }
    }

    I_Phi_file.close();
    std::filesystem::remove(prefix_tempfiles + ".iphi");
    
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

    std::ofstream file_mlf(prefix_tempfiles + ".mlf");
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

    std::ifstream file_mlf(prefix_tempfiles + ".mlf");
    idx.M_LF.load(file_mlf);
    file_mlf.close();
    std::filesystem::remove(prefix_tempfiles + ".mlf");

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_load_mlf=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}