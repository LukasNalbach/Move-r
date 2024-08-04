#pragma once

#include <libsais.h>
#include <libsais64.h>
#include <sais.hxx>
#include <move_r/move_r.hpp>
#include <gtl/btree.hpp>

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support,sym_t,pos_t>::construction::read_t_from_file(std::ifstream& T_ifile) {
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

template <move_r_support support, typename sym_t, typename pos_t>
template <typename sa_sint_t>
void move_r<support,sym_t,pos_t>::construction::build_sa() {
    std::vector<sa_sint_t>& SA = get_sa<sa_sint_t>(); // [0..n-1] The suffix array

    // Choose the correct suffix array construction algorithm.
    if constexpr (byte_alphabet) {
        pos_t fs = std::min<pos_t>(std::numeric_limits<sa_sint_t>::max()-n,6*256);
        no_init_resize(SA,n+fs);
        if (log) std::cout << "building SA" << std::flush;

        if constexpr (std::is_same_v<sa_sint_t,int32_t>) {
            libsais_omp(&T<uint8_t>(0),&SA[0],n,fs,NULL,p);
        } else {
            libsais64_omp(&T<uint8_t>(0),&SA[0],n,fs,NULL,p);
        }

        no_init_resize(SA,n);
    } else {
        if (idx.symbols_remapped) {
            if (log) {
                time = now();
                std::cout << "mapping T to its effective alphabet" << std::flush;
            }

            #pragma omp parallel for num_threads(p)
            for (uint64_t i=0; i<n-1; i++) {
                T<i_sym_t>(i) = (*idx._map_int.find(T<sym_t>(i))).second;
            }

            if (log) time = log_runtime(time);
            if (mode == _suffix_array_space) store_mapintext();
        }
        
        if (log) std::cout << "building SA" << std::flush;

        no_init_resize(SA,n);
        saisxx((i_sym_t*)&T_vec[0],&SA[0],(sa_sint_t)n,(sa_sint_t)idx.sigma);
    }

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_build_sa=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support,sym_t,pos_t>::construction::store_mapintext() {
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

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support,sym_t,pos_t>::construction::load_mapintext() {
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

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support,sym_t,pos_t>::construction::unmap_t() {
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