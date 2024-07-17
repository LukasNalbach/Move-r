#pragma once

#include <libsais.h>
#include <libsais16.h>
#include <libsais16x64.h>
#include <libsais64.h>
#include <sais.hxx>
#include <move_r/move_r.hpp>
#include <gtl/btree.hpp>

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
void move_r<locate_support,sym_t,pos_t>::construction::read_t_from_file(std::ifstream& T_ifile) {
    time = now();
    if (log) std::cout << "reading T" << std::flush;

    T_ifile.seekg(0,std::ios::end);
    n = T_ifile.tellg()+(std::streamsize)+1;
    idx.n = n;
    T_ifile.seekg(0,std::ios::beg);
    no_init_resize(T_str,n);
    read_from_file(T_ifile,T_str.c_str(),n-1);
    T<i_sym_t>(n-1) = 0;

    if (log) time = log_runtime(time);
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
template <typename ls_sym_t, typename sa_sint_t>
void move_r<locate_support,sym_t,pos_t>::construction::run_libsais(ls_sym_t* T, sa_sint_t* SA, pos_t fs) {
    if (log) {
        time = now();
        std::cout << "building SA" << std::flush;
    }

    if constexpr (std::is_same_v<ls_sym_t,uint8_t>) {
        if constexpr (std::is_same_v<sa_sint_t,int32_t>) {
            libsais_omp(T,SA,n,fs,NULL,p);
        } else {
            libsais64_omp(T,SA,n,fs,NULL,p);
        }
    } else if constexpr (std::is_same_v<ls_sym_t,uint16_t>) {
        if constexpr (std::is_same_v<sa_sint_t,int32_t>) {
            libsais16_omp(T,SA,n,fs,NULL,p);
        } else {
            libsais16x64_omp(T,SA,n,fs,NULL,p);
        }
    } else if constexpr (std::is_same_v<ls_sym_t,int32_t>) {
        libsais_int_omp(T,SA,n,idx.sigma,fs,p);
    } else if constexpr (std::is_same_v<ls_sym_t,int64_t>) {
        saisxx(T,SA,(sa_sint_t)n,(sa_sint_t)idx.sigma);
        //libsais64_long_omp(T,SA,n,idx.sigma,fs,p); // buggy
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
template <typename sa_sint_t>
void move_r<locate_support,sym_t,pos_t>::construction::build_sa() {
    if (!byte_alphabet && log) {
        time = now();
        std::cout << "mapping T to its effective alphabet" << std::flush;
    }

    std::vector<sa_sint_t>& SA = get_sa<sa_sint_t>(); // [0..n-1] The suffix array
    pos_t fs = 6*idx.sigma;

    // Choose the correct suffix array construction algorithm.
    if constexpr (byte_alphabet) {
        no_init_resize(SA,n+6*256);
        run_libsais<uint8_t,sa_sint_t>(&T<i_sym_t>(0),&SA[0],6*256);
    } else if constexpr (sizeof(sym_t) == 2) {
        #pragma omp parallel for num_threads(p)
        for (uint64_t i=0; i<n-1; i++) {
            T<uint16_t>(i) = (*idx._map_int.find(T<sym_t>(i))).second;
        }

        if (mode == _suffix_array_space) store_mapintext();
        no_init_resize(SA,n+fs);
        run_libsais<uint16_t,sa_sint_t>((uint16_t*)&T_vec[0],&SA[0],fs);
    } else if (mode == _suffix_array_space) {
        #pragma omp parallel for num_threads(p)
        for (uint64_t i=0; i<n-1; i++) {
            T<i_sym_t>(i) = (*idx._map_int.find(T<sym_t>(i))).second;
        }

        if (log) time = log_runtime(time);
        store_mapintext();
        if (log) std::cout << "building SA" << std::flush;
        no_init_resize(SA,n);
        saisxx((i_sym_t*)&T_vec[0],&SA[0],(sa_sint_t)n,(sa_sint_t)idx.sigma);
    } else if constexpr (sizeof(sym_t) == sizeof(sa_sint_t)) {
        #pragma omp parallel for num_threads(p)
        for (uint64_t i=0; i<n-1; i++) {
            T<sa_sint_t>(i) = (sa_sint_t)(*idx._map_int.find(T<sym_t>(i))).second;
        }
        
        if (log) time = log_runtime(time);
        no_init_resize(SA,n+fs);
        run_libsais<sa_sint_t,sa_sint_t>((sa_sint_t*)&T_vec[0],&SA[0],fs);
    } else {
        std::vector<int64_t> T_vec_ls;
        no_init_resize(T_vec_ls,n);
        T_vec_ls[n-1] = 0;

        #pragma omp parallel for num_threads(p)
        for (uint64_t i=0; i<n-1; i++) {
            T<i_sym_t>(i) = (*idx._map_int.find(T<sym_t>(i))).second;
            T_vec_ls[i] = (int64_t)T<i_sym_t>(i);
        }

        if (log) time = log_runtime(time);
        no_init_resize(SA,n+fs);
        run_libsais<int64_t,sa_sint_t>((int64_t*)&T_vec_ls[0],&SA[0],fs);
    }

    no_init_resize(SA,n);

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_build_sa=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
template <typename sa_sint_t>
void move_r<locate_support,sym_t,pos_t>::construction::build_iphim1_sa() {
    if (log) {
        time = now();
        std::cout << "building I_Phi^{-1}" << std::flush;
    }

    std::vector<sa_sint_t>& SA = get_sa<sa_sint_t>(); // [0..n-1] The suffix array
    
    no_init_resize(I_Phi_m1,r);
    I_Phi_m1[0] = std::make_pair(SA[n-1],SA[0]);

    #pragma omp parallel num_threads(p_)
    {
        uint16_t i_p = omp_get_thread_num();

        // Iteration range start position of thread i_p.
        pos_t b_r = r_p[i_p];

        // Number of BWT runs within thread i_p's section.
        pos_t rp_diff = r_p[i_p+1]-r_p[i_p];

        // start position of the current BWT run.
        pos_t j = n_p[i_p];

        if (rp_diff > 0) {
            if (r_p[i_p] != 0) I_Phi_m1[b_r] = std::make_pair(SA[j-1],SA[j]);
            j += run_len(i_p,0);

            for (pos_t i=1; i<rp_diff; i++) {
                I_Phi_m1[b_r+i] = std::make_pair(SA[j-1],SA[j]);
                j += run_len(i_p,i);
            }
        }
    }

    if (!build_sa_and_l && locate_support != _rlzdsa) {
        SA.clear();
        SA.shrink_to_fit();
    }

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_build_iphi=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
void move_r<locate_support,sym_t,pos_t>::construction::store_iphim1() {
    if (log) {
        time = now();
        std::cout << "storing I_Phi^{-1} to disk" << std::flush;
    }

    std::ofstream file_iphim1(prefix_tmp_files + ".iphim1");
    write_to_file(file_iphim1,(char*)&I_Phi_m1[0].first,r*2*sizeof(pos_t));
    I_Phi_m1.clear();
    I_Phi_m1.shrink_to_fit();
    file_iphim1.close();

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
void move_r<locate_support,sym_t,pos_t>::construction::load_iphim1() {
    if (log) {
        time = now();
        std::cout << "loading I_Phi^{-1} from disk" << std::flush;
    }

    std::ifstream file_iphim1(prefix_tmp_files + ".iphim1");
    no_init_resize(I_Phi_m1,r);
    read_from_file(file_iphim1,(char*)&I_Phi_m1[0].first,r*2*sizeof(pos_t));
    file_iphim1.close();
    std::filesystem::remove(prefix_tmp_files + ".iphim1");

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
void move_r<locate_support,sym_t,pos_t>::construction::store_mapintext() {
    if (log) {
        time = now();
        std::cout << "storing map_int and map_ext to disk" << std::flush;
    }

    std::ofstream file_mapintext(prefix_tmp_files + ".mapintext");

    for (std::pair<sym_t,i_sym_t> p : idx._map_int) {
        file_mapintext.write((char*)&p,sizeof(std::pair<sym_t,i_sym_t>));
    }
    
    idx._map_int.clear();
    write_to_file(file_mapintext,(char*)&idx._map_ext[0],idx.sigma*sizeof(sym_t));
    idx._map_ext.clear();
    idx._map_ext.shrink_to_fit();
    file_mapintext.close();

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
void move_r<locate_support,sym_t,pos_t>::construction::load_mapintext() {
    if (log) {
        time = now();
        std::cout << "loading map_int and map_ext from disk" << std::flush;
    }

    std::ifstream file_mapintext(prefix_tmp_files + ".mapintext");
    std::pair<sym_t,i_sym_t> p;

    for (pos_t i=0; i<idx.sigma; i++) {
        file_mapintext.read((char*)&p,sizeof(std::pair<sym_t,i_sym_t>));
        idx._map_int.emplace(p);
    }

    no_init_resize(idx._map_ext,idx.sigma);
    read_from_file(file_mapintext,(char*)&idx._map_ext[0],idx.sigma*sizeof(sym_t));
    file_mapintext.close();
    std::filesystem::remove(prefix_tmp_files + ".mapintext");

    if (log) {
        time = log_runtime(time);
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
void move_r<locate_support,sym_t,pos_t>::construction::unmap_t() {
    if constexpr (str_input) {
        #pragma omp parallel for num_threads(p)
        for (uint64_t i=0; i<n-1; i++) {
            if (T<i_sym_t>(i) <= max_remapped_to_uchar) {
                T<i_sym_t>(i) = idx._map_ext[T<i_sym_t>(i)];
            }
        }
    } else {
        #pragma omp parallel for num_threads(p)
        for (uint64_t i=0; i<n-1; i++) {
            T<sym_t>(i) = idx._map_ext[T<i_sym_t>(i)];
        }
    }
}