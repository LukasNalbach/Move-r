#include <libsais.h>
#include <libsais64.h>

template <typename uint_t>
void move_r<uint_t>::construction::read_t_from_file_in_memory(std::ifstream& t_file) {
    time = now();
    if (log) std::cout << "reading T" << std::flush;

    t_file.seekg(0,std::ios::end);
    n = t_file.tellg()+(std::streamsize)+1;
    idx.n = n;
    t_file.seekg(0,std::ios::beg);

    no_init_resize(T,n);
    read_from_file(t_file,T.c_str(),n-1);
    T[n-1] = uchar_to_char(terminator);

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
template <typename sa_sint_t, bool read_l>
void move_r<uint_t>::construction::build_rlbwt_c_in_memory() {
    if (log) {
        time = now();
        std::cout << "building RLBWT" << std::flush;
    }

    std::vector<sa_sint_t>& SA = get_sa<sa_sint_t>(); // [0..n-1] The suffix array

    r_p.resize(p+1);
    RLBWT_thr.resize(p);

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
            prev_char = L[b];
        } else {
            prev_char = T[SA[b] == 0 ? n-1 : SA[b]-1];
        }

        // Iterate backwards in L until the next run start position.
        if (p > 0) {
            while (i_ > n_p[i_p-1] && (read_l ? L[i_-1] : T[SA[i_-1] == 0 ? n-1 : SA[i_-1]-1]) == prev_char) {
                i_--;
            }
        }

        #pragma omp barrier

        // communicate the new iteration range start positions in L.
        n_p[i_p] = i_;

        #pragma omp barrier

        // adjust the current threads iteration range end position.
        e = n_p[i_p+1];

        // Iterate over the range L[b+1..e-1]
        for (uint_t i=b+1; i<e; i++) {
            if constexpr (read_l) {
                cur_char = L[i];
            } else {
                cur_char = T[SA[i] == 0 ? n-1 : SA[i]-1];
            }

            if (cur_char != prev_char) {
                RLBWT_thr[i_p].emplace_back(std::make_pair(prev_char,i-i_));
                prev_char = cur_char;
                i_ = i;
            }
        }

        RLBWT_thr[i_p].emplace_back(std::make_pair(prev_char,e-i_));

        // Store in r_p[i_p] the number of runs starting in L[b..e].
        r_p[i_p] = RLBWT_thr[i_p].size();
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

    RLBWT = std::move(interleaved_vectors<uint32_t>({1,4},r,false));

    std::vector<uint_t> r_p_new;
    std::vector<uint_t> n_p_new;
    n_p_new.resize(p+1);
    C.resize(p+1);
    n_p_new[p] = n;

    for (uint16_t i=0; i<p; i++) {
        r_p_new.emplace_back(i*(r/p));
        C[i].resize(256);
    }

    r_p_new.emplace_back(r);

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Iteration range start position of thread i_p.
        uint_t b_r = r_p[i_p];
        // Iteration range end position of thread i_p.
        uint_t e_r = r_p[i_p+1];

        // Iterator pointing to the current RLBWT pair in RLBWT_thr[i_p].
        auto rlbwt_it = RLBWT_thr[i_p].begin();
        // Current position in L.
        uint_t j = n_p[i_p];
        // Index of the current section in r_p_new.
        uint_t cur_rp = bin_search_max_leq<uint_t>(b_r,0,p-1,[&r_p_new](uint_t x){return r_p_new[x];});

        if (b_r == r_p_new[cur_rp]) {
            n_p_new[cur_rp] = j;
        }
        
        // Index of the next section in r_p_new.
        uint_t next_rp = cur_rp+1;

        for (uint_t i=b_r; i<e_r; i++) {
            if (i == r_p_new[next_rp]) {
                n_p_new[next_rp] = j;
                cur_rp++;
                next_rp++;
            }

            // Count the occurrences of L[i-1] within the last run in C[i_p][L[i-1]].
            #pragma omp atomic
            C[cur_rp][char_to_uchar((*rlbwt_it).first)] += (*rlbwt_it).second;

            set_run_char(i,(*rlbwt_it).first);
            set_run_length(i,(*rlbwt_it).second);

            j += (*rlbwt_it).second;
            rlbwt_it++;
        }
    }

    RLBWT_thr.clear();
    RLBWT_thr.shrink_to_fit();

    r_p = std::move(r_p_new);
    n_p = std::move(n_p_new);

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

        // Bwt runs iteration range start position of thread i_p.
        uint_t b_r = r_p[i_p];
        // Bwt runs iteration range end position of thread i_p.
        uint_t e_r = r_p[i_p+1];

        // start position of the current BWT run.
        uint_t j = n_p[i_p];

        if (b_r != 0) I_Phi[b_r] = std::make_pair(SA[j],SA[j-1]);
        j += run_length(b_r);

        for (uint_t i=b_r+1; i<e_r; i++) {
            I_Phi[i] = std::make_pair(SA[j],SA[j-1]);
            j += run_length(i);
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