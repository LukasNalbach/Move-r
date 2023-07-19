#include <libsais.h>
#include <libsais64.h>

template <typename uint_t>
void move_r<uint_t>::construction::read_t_from_file_in_memory(std::ifstream& t_file) {
    time = now();
    if (log) std::cout << "reading T";

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
    (*reinterpret_cast<std::vector<no_init<sa_sint_t>>*>(&SA)).resize(n);

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
        if (measurement_file_index != NULL) *measurement_file_index << " time_build_sa=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
template <typename sa_sint_t>
void move_r<uint_t>::construction::build_l_and_c_in_memory() {
    if (log) {
        time = now();
        std::cout << "building L and C" << std::flush;
    }

    std::vector<sa_sint_t>& SA = get_sa<sa_sint_t>(); // [0..n-1] The suffix array
    no_init_resize(L,n);

    C.resize(p+1);
    r_p.resize(p+1);

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        C[i_p].resize(256,0);

        // Iteration range start position of thread i_p.
        uint_t b = i_p*(n/p);
        // Iteration range end position of thread i_p.
        uint_t e = i_p == p-1 ? n-1 : (i_p+1)*(n/p)-1;

        // Build L in the range [b..e].
        for (uint_t i=b; i<=e; i++) {
            L[i] = SA[i] == 0 ? T[n-1] : T[SA[i]-1];
        }

        // Store in C[i_p][c] the number of occurences of c in L[b..e], for each c in [0..255].
        // Store in r_p[i_p] the number of runs starting in L[b..e].

        // We add an artificial run starting at b.
        C[i_p][char_to_uchar(L[b])]++;
        r_p[i_p] = 1;
        
        // Iterate over the range L[b..e]
        for (uint_t i=b+1; i<=e; i++) {

            // Count the occurence of L[i] in C[i_p][L[i]].
            C[i_p][char_to_uchar(L[i])]++;

            // If there is a run starting at position i, count it in r_p[i_p].
            if (L[i] != L[i-1]) {
                r_p[i_p]++;
            }
        }
    }

    if (&T == &T_tmp) {
        T.clear();
        T.shrink_to_fit();
    }

    if (!build_locate_support) {
        SA.clear();
        SA.shrink_to_fit();
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::process_rp_in_memory() {
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

    if (log) {
        double n_r = std::round(100.0*(n/(double)r))/100.0;
        if (measurement_file_index != NULL) {
            *measurement_file_index << " time_build_l_c=" << time_diff_ns(time,now());
            *measurement_file_index << " n=" << n;
            *measurement_file_index << " sigma=" << std::to_string(sigma);
            *measurement_file_index << " r=" << r;
        }
        time = log_runtime(time);
        std::cout << "n = " << n << ", sigma = " << std::to_string(sigma) << ", r = " << r << ", n/r = " << n_r << std::endl;
    }

    if (p > 1 && 1000*p > r) {
        p = std::max((uint_t)1,r/1000);
        if (log) std::cout << "warning: p > r/1000, setting p to r/1000 ~ " << std::to_string(p) << std::endl;
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::build_ilf_and_bwt_run_heads_in_memory() {
    if (log) {
        time = now();
        std::cout << "building I_LF and bwt run heads" << std::flush;
    }

    (*reinterpret_cast<std::vector<std::pair<no_init<uint_t>,no_init<uint_t>>>*>(&I_LF)).resize(r);
    no_init_resize(bwt_run_heads,r);

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Iteration range start position of thread i_p.
        uint_t b = i_p*(n/p);
        // Iteration range end position of thread i_p.
        uint_t e = i_p == p-1 ? n-1 : (i_p+1)*(n/p)-1;

        // Start position of the last-seen run.
        uint_t i_ = b;
        // Position in I_LF and M_Phi to write the next pairs to.
        uint_t j = r_p[i_p];

        // Build I_LF[r_p[i_p]..r_p[i_p+1]-1]

        /* Each thread creates one input interval starting at b,
        whiches pair has to be placed at LF[r_p[i_p]]. */
        I_LF[j] = std::make_pair(i_,C[p][char_to_uchar(L[b])]+C[i_p][char_to_uchar(L[b])]);
        bwt_run_heads[j] = L[i_];

        j++;

        // Iterate over the range L[b+1..e]
        for (uint_t i=b+1; i<=e; i++) {

            // Check, if a run starts at a position i
            if (L[i] != L[i-1]) {
                bwt_run_heads[j] = L[i];

                /* Update the rank-function in C[i_p] to store C[i_p][c] = rank(L,c,i-1),
                for each c in [0..255]. */
                C[i_p][char_to_uchar(L[i-1])] += i-i_;

                // Update the position i_ of the last-seen run to i.
                i_ = i;

                /* Write the pair (i,LF(i)) to the next position j in I_LF, where
                LF(i) = C[L[i]] + rank(L,L[i],i-1) = C[p][L[i]] + C[i_p][L[i]]. */
                I_LF[j] = std::make_pair(i,C[p][char_to_uchar(L[i])]+C[i_p][char_to_uchar(L[i])]);

                j++;
            }
        }
    }

    C.clear();
    C.shrink_to_fit();

    L.clear();
    L.shrink_to_fit();

    if (log) {
        if (measurement_file_index != NULL) *measurement_file_index << " time_build_ilf_bwt_run_heads=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::build_br_in_memory() {
    if (log) {
        time = now();
        std::cout << "building B_r" << std::flush;
    }
    
    // build bitvectors, that mark the run head positions in each threads section
    // use an sd_vector, if the bitvectors are sparse (r << n)
    compress_br = (n/r) > 7;

    if (compress_br) {
        B_r_sdsl.resize(p);
    } else {
        B_r_pasta.resize(p,NULL);
    }

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Bwt range start position of thread i_p.
        uint_t b = i_p*(n/p);
        // Bwt range start position of thread i_p.
        uint_t e = i_p == p-1 ? n-1 : (i_p+1)*(n/p)-1;

        // Iteration range start position of thread i_p in I_LF.
        uint_t b_r = r_p[i_p];
        // Iteration range end position of thread i_p in I_LF.
        uint_t e_r = r_p[i_p+1]-1;

        if (compress_br) {
            // add one bit and mark the next section's first run start position
            B_r_sdsl[i_p] = std::move(sdsl::bit_vector(e-b+2));
            B_r_sdsl[i_p][e-b+1] = 1;

            for (uint_t i=b_r; i<=e_r; i++) {
                B_r_sdsl[i_p][I_LF[i].first-b] = 1;
            }
        } else {
            // add one bit and mark the next section's first run start position
            B_r_pasta[i_p] = new pasta::BitVector(e-b+2,0);
            (*B_r_pasta[i_p])[e-b+1] = 1;

            for (uint_t i=b_r; i<=e_r; i++) {
                (*B_r_pasta[i_p])[I_LF[i].first-b] = 1;
            }
        }
    }

    // compress B_r using sd_arrays
    if (compress_br) {
        B_r_compr.resize(p);

        #pragma omp parallel num_threads(p)
        {
            // Index in [0..p-1] of the current thread.
            uint16_t i_p = omp_get_thread_num();
            B_r_compr[i_p] = std::move(sd_array<uint_t>(B_r_sdsl[i_p]));
        }

        B_r_sdsl.clear();
        B_r_sdsl.shrink_to_fit();
    }

    if (log) {
        time = log_runtime(time);
    }
}

template <typename uint_t>
template <typename sa_sint_t>
void move_r<uint_t>::construction::build_l__and_iphi_in_memory() {
    if (log) {
        time = now();
        std::string msg;
        if (build_locate_support) {
            msg = "building L' and I_Phi";
        } else {
            msg = "building L'";
        }
        std::cout << msg << std::flush;
    }

    std::vector<sa_sint_t>& SA = get_sa<sa_sint_t>(); // [0..n-1] The suffix array

    if (build_locate_support) {
        (*reinterpret_cast<std::vector<std::pair<no_init<uint_t>,no_init<uint_t>>>*>(&I_Phi)).resize(r);
        I_Phi[0] = std::make_pair(SA[0],SA[n-1]);
    }

    // Simultaneously iterate over the input intervals of M_LF, the bwt run start positions and SA to build I_Phi and L'
    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        std::function<uint_t(uint_t)> Br_ip_select_1; // function that comutes select_1 on B_r[i_p]
        pasta::FlatRankSelect B_r_pasta_select_1;

        if (compress_br) {
            Br_ip_select_1 = [this,&i_p](uint_t idx){return B_r_compr[i_p].select_1(idx);};
        } else {
            B_r_pasta_select_1 = std::move(pasta::FlatRankSelect(*B_r_pasta[i_p]));
            Br_ip_select_1 = [&B_r_pasta_select_1](uint_t idx){return B_r_pasta_select_1.select1(idx);};
        }

        // Bwt range start position of thread i_p.
        uint_t b = i_p*(n/p);

        // Bwt runs iteration range start position of thread i_p.
        uint_t b_r = r_p[i_p];
        // Bwt runs iteration range end position of thread i_p.
        uint_t e_r = r_p[i_p+1]-1;

        /**
         * @brief Index of the current input interval in M_LF, initially
         *        the index of the input interval of M_LF, in which b lies.
         */
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

        uint_t l_ = b; // Starting position of the next bwt run.

        for (uint_t i=b_r; i<=e_r; i++) {
            if (build_locate_support && i != 0) {
                I_Phi[i] = std::make_pair(SA[l_],SA[l_-1]);
            }

            // update l_ to the next run start position
            l_ = b + Br_ip_select_1(i-b_r+2);

            // iterate over all input intervals in M_LF within the x-th bwt run, that have been created by the balancing algorithm
            do {
                idx.M_LF.template set_character<char>(j,bwt_run_heads[i]);
                j++;
            } while (idx.M_LF.p(j) < l_);
        }
    }

    r_p.clear();
    r_p.shrink_to_fit();

    bwt_run_heads.clear();
    bwt_run_heads.shrink_to_fit();

    if (compress_br) {
        B_r_compr.clear();
        B_r_compr.shrink_to_fit();
    } else {
        for (uint16_t i=0; i<p; i++) {
            delete B_r_pasta[i];
            B_r_pasta[i] = NULL;
        }

        B_r_pasta.clear();
        B_r_pasta.shrink_to_fit();
    }

    if (!build_locate_support) {
        if (compress_br) {
            B_r_sdsl.clear();
            B_r_sdsl.shrink_to_fit();
        } else {
            B_r_compr.clear();
            B_r_compr.shrink_to_fit();
        }
    }

    if (log) {
        if (measurement_file_index != NULL) *measurement_file_index << " time_build_l_=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
template <typename sa_sint_t>
void move_r<uint_t>::construction::build_sas_in_memory() {
    if (log) {
        time = now();
        std::cout << "building SA_s" << std::flush;
    }

    std::vector<sa_sint_t>& SA = get_sa<sa_sint_t>(); // [0..n-1] The suffix array

    (*reinterpret_cast<std::vector<no_init<uint_t>>*>(&SA_s)).resize(r_);

    // build SA_s
    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<r_; i++) {
        SA_s[i] = SA[idx.M_LF.p(i+1)-1];
    }

    // Now we do not need the suffix array anymore
    SA.clear();
    SA.shrink_to_fit();

    if (log) {
        if (measurement_file_index != NULL) *measurement_file_index << " time_build_sas=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::unmap_t() {
    #pragma omp parallel for num_threads(p)
    for (uint_t i=0; i<n-1; i++) {
        T[i] = idx.unmap_from_internal(T[i]);
    }
}