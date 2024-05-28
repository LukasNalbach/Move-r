#include <libsais.h>
#include <libsais16.h>
#include <libsais64.h>
#include <sais.hxx>
#include <absl/container/btree_map.h>

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
    T(n-1) = uchar_to_char(0);

    if (log) time = log_runtime(time);
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
template <typename inp_t>
void move_r<locate_support,sym_t,pos_t>::construction::execute_libsais(inp_t* T, int32_t* SA, pos_t fs) {
    if constexpr (std::is_same<inp_t,uint8_t>::value) {
        if (p == 1) {
            libsais(T,SA,n,fs,NULL);
        } else {
            libsais_omp(T,SA,n,fs,NULL,p);
        }
    } else if constexpr (std::is_same<inp_t,uint16_t>::value) {
        if (p == 1) {
            libsais16(T,SA,n,fs,NULL);
        } else {
            libsais16_omp(T,SA,n,fs,NULL,p);
        }
    } else {
        if (p == 1) {
            libsais_int(T,SA,n,idx.sigma,fs);
        } else {
            libsais_int_omp(T,SA,n,idx.sigma,fs,p);
        }
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
template <typename sa_sint_t, typename inp_t>
void move_r<locate_support,sym_t,pos_t>::construction::build_sa() {
    if (log) {
        time = now();
        std::cout << "building SA" << std::flush;
    }

    std::vector<sa_sint_t>& SA = get_sa<sa_sint_t>(); // [0..n-1] The suffix array

    // Choose the correct suffix array construction algorithm.
    if constexpr (std::is_same<sym_t,char>::value) {
        no_init_resize(SA,n+6*256);

        if constexpr (std::is_same<sa_sint_t,int32_t>::value) {
            if (p == 1) {
                libsais((uint8_t*)&T_str[0],&SA[0],n,6*256,NULL);
            } else {
                libsais_omp((uint8_t*)&T_str[0],&SA[0],n,6*256,NULL,p);
            }
        } else {
            if (p == 1) {
                libsais64((uint8_t*)&T_str[0],&SA[0],n,6*256,NULL);
            } else {
                libsais64_omp((uint8_t*)&T_str[0],&SA[0],n,6*256,NULL,p);
            }
        }

        no_init_resize(SA,n);
    } else {
        if constexpr (std::is_same<sa_sint_t,int32_t>::value) {
            pos_t fs = 6*idx.sigma;
            no_init_resize(SA,n+fs);

            if constexpr (std::is_same<sym_t,inp_t>::value) {
                execute_libsais<inp_t>(&T_vec[0],&SA[0],fs);
            } else {
                std::vector<inp_t> T_tmp;
                no_init_resize(T_tmp,n);

                #pragma omp parallel for num_threads(p)
                for (uint64_t i=0; i<n; i++) {
                    T_tmp[i] = T(i);
                }

                execute_libsais<inp_t>(&T_tmp[0],&SA[0],fs);
            }

            no_init_resize(SA,n);
        } else {
            no_init_resize(SA,n);
            saisxx(T_vec.begin(),SA.begin(),(sa_sint_t)n,(sa_sint_t)idx.sigma);
        }
    }

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_build_sa=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
template <bool read_l, typename sa_sint_t>
void move_r<locate_support,sym_t,pos_t>::construction::build_rlbwt_c_libsais() {
    if (log) {
        time = now();
        std::cout << "building RLBWT" << std::flush;
    }

    std::vector<sa_sint_t>& SA = get_sa<sa_sint_t>(); // [0..n-1] The suffix array

    r_p.resize(p+1,0);
    RLBWT.resize(p,std::move(interleaved_vectors<uint32_t,uint32_t>({(uint8_t)std::ceil(std::log2(idx.sigma+1)/(double)8),4})));

    if constexpr (std::is_same<sym_t,char>::value) {
        C.resize(p+1,std::vector<pos_t>(256,0));
    } else {
        C.resize(p+1,std::vector<pos_t>(idx.sigma,0));
    }

    for (uint16_t i=0; i<p; i++) {
        n_p.emplace_back(i*(n/p));
    }

    n_p.emplace_back(n);

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Iteration range start position of thread i_p.
        pos_t b = n_p[i_p];
        // Iteration range end position of thread i_p.
        pos_t e = n_p[i_p+1];

        // Store in C[i_p][c] the number of occurrences of c in L[b..e], for each c in [0..255].

        sym_t prev_sym; // previous symbol in L.
        sym_t cur_sym; // current symbol in L.
        pos_t i_ = b; // start position of the last seen run in L.

        if constexpr (read_l) {
            prev_sym = L[b];
        } else {
            prev_sym = SA[b] == 0 ? 0 : T(SA[b]-1);
        }

        // Iterate over the range L[b+1..e-1]
        for (pos_t i=b+1; i<e; i++) {
            if constexpr (read_l) {
                cur_sym = L[i];
            } else {
                cur_sym = SA[i] == 0 ? 0 : T(SA[i]-1);
            }

            // check if there is a run starting at L[i]
            if (cur_sym != prev_sym) {
                RLBWT[i_p].emplace_back_unsafe<sym_t,uint32_t>({prev_sym,i-i_});
                C[i_p][symbol_idx(prev_sym)] += i-i_;
                prev_sym = cur_sym;
                i_ = i;
            }
        }

        // add the run L[i'..e)
        RLBWT[i_p].emplace_back_unsafe<sym_t,uint32_t>({prev_sym,e-i_});
        C[i_p][symbol_idx(prev_sym)] += e-i_;
        // Store in r_p[i_p] the number of runs starting in L[b..e).
        r_p[i_p] = RLBWT[i_p].size();
        RLBWT[i_p].shrink_to_fit();
    }

    // for i_p \in [1,p-2], merge the last run in thread i_p's section with the first run in thread
    // i_p+1's section, if their characters are equal
    for (uint16_t i_p=0; i_p<p-1; i_p++) {
        sym_t c = run_symbol(i_p,r_p[i_p]-1);

        if (run_symbol(i_p+1,0) == c) {
            pos_t l = run_length(i_p,r_p[i_p]-1);
            set_run_length(i_p+1,0,run_length(i_p+1,0)+l);
            C[i_p][symbol_idx(c)] -= l;
            C[i_p+1][symbol_idx(c)] += l;
            n_p[i_p+1] -= l;
            r_p[i_p]--;
            RLBWT[i_p].resize(r_p[i_p]);
        }
    }

    if (&L == &L_tmp) {
        L.clear();
        L.shrink_to_fit();
    }

    if (delete_T) {
        T_str.clear();
        T_str.shrink_to_fit();
        T_vec.clear();
        T_vec.shrink_to_fit();
    }

    if (!build_locate_support && !read_l) {
        SA.clear();
        SA.shrink_to_fit();
    }

    /* Now, r_p[i_p] stores the number of runs starting in the iteration range L[b..e] of thread
    i_p in [0..p-1], and r_p[p] = 0. We want r_p[i_p] to store the number of runs starting before
    the iteration range start position b of thread i_p in [0..p-1]. Also, we want r_p[p] to store
    the number of all runs. This way, we can build I_LF[r_p[i_p]..r_p[i_p+1]-1] with thread i_p in [0..p-1]
    using the C-array in C[p] while using and updating the rank-function in C[i_p] on-the-fly. */

    for (uint16_t i=p; i>0; i--) {
        r_p[i] = r_p[i-1];
    }

    r_p[0] = 0;

    for (uint16_t i=2; i<=p; i++) {
        r_p[i] += r_p[i-1];
    }

    /* Now, r_p[i_p] stores the number of runs starting before the iteration range start position b of
    thread i_p in [0..p-1] and r_p[p] stores the number of all runs, so we are done with r_p[0..p] */
    
    r = r_p[p];
    idx.r = r;

    process_c();

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_build_rlbwt=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
template <typename sa_sint_t>
void move_r<locate_support,sym_t,pos_t>::construction::build_iphi_from_sa() {
    if (log) {
        time = now();
        std::cout << "building I_Phi" << std::flush;
    }

    std::vector<sa_sint_t>& SA = get_sa<sa_sint_t>(); // [0..n-1] The suffix array
    
    no_init_resize(I_Phi,r);
    I_Phi[0] = std::make_pair(SA[0],SA[n-1]);

    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        // Iteration range start position of thread i_p.
        pos_t b_r = r_p[i_p];

        // Number of BWT runs within thread i_p's section.
        pos_t rp_diff = r_p[i_p+1]-r_p[i_p];

        // start position of the current BWT run.
        pos_t j = n_p[i_p];

        if (rp_diff > 0) {
            if (r_p[i_p] != 0) I_Phi[b_r] = std::make_pair(SA[j],SA[j-1]);
            j += run_length(i_p,0);

            for (pos_t i=1; i<rp_diff; i++) {
                I_Phi[b_r+i] = std::make_pair(SA[j],SA[j-1]);
                j += run_length(i_p,i);
            }
        }
    }

    if (!build_from_sa_and_l && locate_support != _rlzdsa) {
        SA.clear();
        SA.shrink_to_fit();
    }

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_build_iphi=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
void move_r<locate_support,sym_t,pos_t>::construction::unmap_t() {
    if constexpr (std::is_same<sym_t,char>::value) {
        #pragma omp parallel for num_threads(p)
        for (uint64_t i=0; i<n-1; i++) {
            if (char_to_uchar(T(i)) <= max_remapped_to_uchar) {
                T(i) = idx.unmap_symbol(T(i));
            }
        }
    } else {
        #pragma omp parallel for num_threads(p)
        for (uint64_t i=0; i<n-1; i++) {
            T(i) = idx._unmap_int[T(i)];
        }
    }
}