#include <libsais.h>
#include <libsais64.h>

template <typename uint_t>
void move_r<uint_t>::construction::read_t_from_file_in_memory(std::ifstream& T_ifile) {
    time = now();
    if (log) std::cout << "reading T" << std::flush;

    T_ifile.seekg(0,std::ios::end);
    n = T_ifile.tellg()+(std::streamsize)+1;
    idx.n = n;
    T_ifile.seekg(0,std::ios::beg);
    no_init_resize(T,n);
    read_from_file(T_ifile,T.c_str(),n-1);
    T[n-1] = uchar_to_char(1);

    if (log) time = log_runtime(time);
}

template <typename uint_t>
template <typename sa_sint_t>
void move_r<uint_t>::construction::build_sa_in_memory() {
    if (log) {
        time = now();
        std::cout << "building SA" << std::flush;
    }

    std::vector<sa_sint_t>& SA = get_sa<sa_sint_t>(); // [0..n-1] The suffix array
    no_init_resize(SA,n);

    // Choose the correct suffix array construction algorithm.
    if constexpr (std::is_same<sa_sint_t,int32_t>::value) {
        if (p == 1) {
            libsais((uint8_t*)&T[0],&SA[0],n,0,NULL);
        } else {
            libsais_omp((uint8_t*)&T[0],&SA[0],n,0,NULL,p);
        }
    } else {
        if (p == 1) {
            libsais64((uint8_t*)&T[0],&SA[0],n,0,NULL);
        } else {
            libsais64_omp((uint8_t*)&T[0],&SA[0],n,0,NULL,p);
        }
    }
    
    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_build_sa=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
template <bool read_l, typename sa_sint_t>
void move_r<uint_t>::construction::build_rlbwt_c_in_memory() {
    if (log) {
        time = now();
        std::cout << "building RLBWT" << std::flush;
    }

    std::vector<sa_sint_t>& SA = get_sa<sa_sint_t>(); // [0..n-1] The suffix array

    r_p.resize(p+1,0);
    RLBWT.resize(p,std::move(interleaved_vectors<uint32_t>({1,4})));
    C.resize(p+1,std::vector<uint_t>(256,0));

    for (uint16_t i=0; i<p; i++) {
        n_p.emplace_back(i*(n/p));
    }

    n_p.emplace_back(n);

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Iteration range start position of thread i_p.
        uint_t b = n_p[i_p];
        // Iteration range end position of thread i_p.
        uint_t e = n_p[i_p+1];

        // Store in C[i_p][c] the number of occurrences of c in L[b..e], for each c in [0..255].

        char prev_char; // previous character in L.
        char cur_char; // current character in L.
        uint_t i_ = b; // start position of the last seen run in L.

        if constexpr (read_l) {
            prev_char = L[b] == uchar_to_char(1) ? uchar_to_char(0) : L[b];
        } else {
            prev_char = SA[b] == 0 ? uchar_to_char(0) : T[SA[b]-1];
        }

        // Iterate over the range L[b+1..e-1]
        for (uint_t i=b+1; i<e; i++) {
            if constexpr (read_l) {
                cur_char = L[i] == uchar_to_char(1) ? uchar_to_char(0) : L[i];
            } else {
                cur_char = SA[i] == 0 ? uchar_to_char(0) : T[SA[i]-1];
            }

            // check if there is a run starting at L[i]
            if (cur_char != prev_char) {
                RLBWT[i_p].emplace_back_unsafe<char,uint32_t>({prev_char,i-i_});
                C[i_p][char_to_uchar(prev_char)] += i-i_;
                prev_char = cur_char;
                i_ = i;
            }
        }

        // add the run L[i'..e)
        RLBWT[i_p].emplace_back_unsafe<char,uint32_t>({prev_char,e-i_});
        C[i_p][char_to_uchar(prev_char)] += e-i_;

        // Store in r_p[i_p] the number of runs starting in L[b..e).
        r_p[i_p] = RLBWT[i_p].size();
        RLBWT[i_p].shrink_to_fit();
    }

    // for i_p \in [1,p-2], merge the last run in thread i_p's section with the first run in thread
    // i_p+1's section, if their characters are equal
    for (uint16_t i_p=0; i_p<p-1; i_p++) {
        uint8_t c = run_uchar(i_p,r_p[i_p]-1);

        if (run_uchar(i_p+1,0) == c) {
            uint_t l = run_length(i_p,r_p[i_p]-1);
            set_run_length(i_p+1,0,run_length(i_p+1,0)+l);
            C[i_p][c] -= l;
            C[i_p+1][c] += l;
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
        T.clear();
        T.shrink_to_fit();
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

template <typename uint_t>
template <typename sa_sint_t>
void move_r<uint_t>::construction::build_iphi_in_memory() {
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
        uint_t b_r = r_p[i_p];

        // Number of BWT runs within thread i_p's section.
        uint_t rp_diff = r_p[i_p+1]-r_p[i_p];

        // start position of the current BWT run.
        uint_t j = n_p[i_p];

        if (rp_diff > 0) {
            if (r_p[i_p] != 0) I_Phi[b_r] = std::make_pair(SA[j],SA[j-1]);
            j += run_length(i_p,0);

            for (uint_t i=1; i<rp_diff; i++) {
                I_Phi[b_r+i] = std::make_pair(SA[j],SA[j-1]);
                j += run_length(i_p,i);
            }
        }
    }

    if (!build_from_sa_and_l) {
        SA.clear();
        SA.shrink_to_fit();
    }

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_build_iphi=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::unmap_t() {
    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<n-1; i++) {
        if (char_to_uchar(T[i]) <= max_remapped_to_uchar) {
            T[i] = idx.unmap_from_internal(T[i]);
        }
    }
}