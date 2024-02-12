#include <ips4o.hpp>

template <typename uint_t>
void move_r<uint_t>::construction::preprocess_t(bool use_libsais, std::ifstream* T_ifile) {
    if (log) std::cout << "preprocessing T" << std::flush;

    if (!use_libsais) {
        T_ifile->seekg(0,std::ios::end);
        n = T_ifile->tellg()+(std::streamsize)+1;
        idx.n = n;
        T_ifile->seekg(0,std::ios::beg);
    }

    // contains_uchar_thr[i_p][c] = true <=> thread i_p found the character c in its section of T[0..n-2].
    std::vector<std::vector<uint8_t>> contains_uchar_thr(p,std::vector<uint8_t>(256,0));

    if (use_libsais) {
        // Iterate over T[0..n-2] and report the occurrence of each found character in T[0..n-2] in contains_uchar_thr.
        #pragma omp parallel for num_threads(p)
        for (uint64_t i=0; i<n-1; i++) {
            contains_uchar_thr[omp_get_thread_num()][char_to_uchar(T[i])] = 1;
        }
    } else {
        uint_t max_t_buf_size = std::max((uint_t)1,n/500);
        std::string T_buf;
        no_init_resize(T_buf,max_t_buf_size);
        uint_t cur_t_buf_size;
        uint_t n_ = n-1;

        while (n_ > 0) {
            cur_t_buf_size = std::min(n_,max_t_buf_size);
            read_from_file(*T_ifile,T_buf.c_str(),cur_t_buf_size);
            
            // Iterate over T[0..n-2] and report the occurrence of each found character in T[0..n-2] in contains_uchar_thr.
            #pragma omp parallel for num_threads(p)
            for (uint64_t i=0; i<cur_t_buf_size; i++) {
                contains_uchar_thr[omp_get_thread_num()][char_to_uchar(T_buf[i])] = 1;
            }

            n_ -= cur_t_buf_size;
        }

        T_ifile->seekg(0,std::ios::beg);
    }
    
    // contains_uchar[c] = 1 <=> c in T[0..n-2].
    std::vector<uint8_t> contains_uchar(256,0);

    /* Combine the results of each thread's sections in contains_uchar_thr[0..i_p-1][0..255]
    into contains_uchar[0..255]. */
    for (uint16_t cur_uchar=0; cur_uchar<256; cur_uchar++) {
        for (uint16_t i_p=0; i_p<p; i_p++) {
            if (contains_uchar_thr[i_p][cur_uchar] == 1) {
                contains_uchar[cur_uchar] = 1;
                break;
            }
        }
    }

    contains_uchar_thr.clear();
    contains_uchar_thr.shrink_to_fit();

    // The number of distinct characters in T[0..n-1].
    idx.sigma = 1;

    // Count the number of ones in contains_uchar[0..255] in sigma.
    for (uint16_t cur_uchar=0; cur_uchar<256; cur_uchar++) {
        if (contains_uchar[cur_uchar] == 1) {
            idx.sigma++;
        }
    }
    
    bool contains_invalid_char = false;

    for (uint8_t cur_uchar=0; cur_uchar<min_valid_char; cur_uchar++) {
        if (contains_uchar[cur_uchar] == 1) {
            contains_invalid_char = true;
            break;
        }
    }

    // If the input contains too many distinct characters, we have to remap the characters in T[0..n-2].
    if (contains_invalid_char) {
        if (idx.sigma > 256-min_valid_char) {
            /* If T[0..n-2] contains more than 256 - min_valid_char distinct characters, we cannot remap them into the 
            range [0..255] without using a character less than min_valid_char, hence we cannot build an index for T. */
            std::cout << "Error: the input contains more than " << std::to_string(256-min_valid_char-1) << " distinct characters" << std::endl;
            return;
        }

        idx.chars_remapped = true;
        
        // build the mapping function map_char that remaps the characters of T, s.t. it does
        // not contain 0 or 1; also build its inverse function unmap_char

        idx._map_char.resize(256,0);
        idx._unmap_char.resize(256,0);

        /* To preserve the order among characters in T[0..n-2], we start by mapping smallest
           character in T[0..n-2] to 2, the second smallest to 3, ... . */

        // The character, to map the currently next largest character in T[0..n-2] to.
        uint16_t next_uchar_to_remap_to = min_valid_char;
        max_remapped_uchar = 0;

        for (uint16_t cur_uchar=0; cur_uchar<next_uchar_to_remap_to; cur_uchar++) {
            if (contains_uchar[cur_uchar] == 1) {
                idx._map_char[cur_uchar] = next_uchar_to_remap_to;
                idx._unmap_char[next_uchar_to_remap_to] = cur_uchar;
                max_remapped_uchar = cur_uchar;
                next_uchar_to_remap_to++;
            }
        }

        max_remapped_to_uchar = next_uchar_to_remap_to - 1;

        for (uint16_t cur_uchar=max_remapped_to_uchar+1; cur_uchar<256; cur_uchar++) {
            if (contains_uchar[cur_uchar] == 1) {
                idx._map_char[cur_uchar] = cur_uchar;
                idx._unmap_char[cur_uchar] = cur_uchar;
            }
        }

        // Apply map_char to T.
        if (use_libsais) {
            #pragma omp parallel for num_threads(p)
            for (uint64_t i=0; i<n-1; i++) {
                if (char_to_uchar(T[i]) <= max_remapped_uchar) {
                    T[i] = idx.map_char(T[i]);
                }
            }
        }
    }

    if (!use_libsais) {
        prefix_tmp_files = "move-r_" + random_alphanumeric_string(10);
        std::ofstream T_ofile(prefix_tmp_files);
        uint_t max_t_buf_size = std::max((uint_t)1,n/500);
        std::string T_buf;
        no_init_resize(T_buf,max_t_buf_size);
        uint_t cur_t_buf_size;
        uint_t n_ = n-1;

        while (n_ > 0) {
            cur_t_buf_size = std::min(n_,max_t_buf_size);
            read_from_file(*T_ifile,T_buf.c_str(),cur_t_buf_size);

            if (idx.chars_remapped) {
                #pragma omp parallel for num_threads(p)
                for (uint64_t i=0; i<cur_t_buf_size; i++) {
                    if (char_to_uchar(T_buf[i]) <= max_remapped_uchar) {
                        T_buf[i] = idx.map_char(T_buf[i]);
                    }
                }
            }

            write_to_file(T_ofile,T_buf.c_str(),cur_t_buf_size);
            n_ -= cur_t_buf_size;
        }
        
        T_ifile->close();
        T_ofile.close();
    }

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_preprocess_t=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::process_c() {
    /* Now, C[i_p][c] is the number of occurrences of c in L[b..e], where [b..e] is the range of the
    thread i_p in [0..p-1]. Also, we have C[p][0..255] = 0. */

    /* We want to have C[i_p][c] = rank(L,c,b-1), where b is the iteration range start position of
    thread i_p in [0..p-1]. Also, we want C[p][0..255] to be the C-array, that is C[p][c] stores
    the number of occurrences of all smaller characters c' < c in L[0..n-1], for c in [0..255]. */

    for (uint16_t i=1; i<p; i++) {
        for (uint16_t j=0; j<256; j++) {
            C[i][j] += C[i-1][j];
        }
    }

    /* Now, we have C[i_p][c] = rank(L,c,e), for each c in [0..255] and i_p in [0..p-1],
    where e is the iteration range end position of thread i_p. */

    for (uint16_t i=p; i>0; i--) {
        for (uint16_t j=0; j<256; j++) {
            C[i][j] = C[i-1][j];
        }
    }
    
    for (uint16_t j=0; j<256; j++) {
        C[0][j] = 0;
    }

    /* Now, we have C[i_p][c] = rank(L,c,b-1), for each c in [0..255] and i_p in [0..p-1],
    where b is the iteration range start position of thread i_p, so we are done with C[0..p-1][0..255].
    Also, we have C[p][c] = rank(L,c,n-1), for c in [0..255]. */

    for (uint16_t i=255; i>0; i--) {
        C[p][i] = C[p][i-1];
    }

    C[p][0] = 0;

    // Now, we have C[p][c] = rank(L,c-1,n-1), for c in [1..255], and C[p][0] = 0.

    for (uint16_t i=2; i<256; i++) {
        C[p][i] += C[p][i-1];
    }

    // Now we are done with C, since C[p] is the C-array.
}

template <typename uint_t>
void move_r<uint_t>::construction::build_ilf() {
    if (log) {
        time = now();
        std::cout << "building I_LF" << std::flush;
    }

    no_init_resize(I_LF,r);

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Iteration range start position of thread i_p.
        uint_t b_r = r_p[i_p];

        // Number of BWT runs in thread i_p's section.
        uint_t rp_diff = r_p[i_p+1]-r_p[i_p];

        // i', Start position of the last-seen run.
        uint_t i_ = n_p[i_p];

        // Build I_LF
        for (uint_t i=0; i<rp_diff; i++) {
            /* Write the pair (i',LF(i')) to the next position i in I_LF, where
            LF(i') = C[L[i']] + rank(L,L[i'],i'-1) = C[p][L[i']] + C[i_p][L[i']]. */
            I_LF[b_r+i] = std::make_pair(i_,C[p][run_uchar(i_p,i)]+C[i_p][run_uchar(i_p,i)]);

            /* Update the rank-function in C[i_p] to store C[i_p][c] = rank(L,c,i'-1),
            for each c in [0..255] */
            C[i_p][run_uchar(i_p,i)] += run_length(i_p,i);

            // Update the position of the last-seen run.
            i_ += run_length(i_p,i);
        }
    }

    C.clear();
    C.shrink_to_fit();

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_build_ilf=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::build_mlf() {
    if (log) {
        if (mf_mds != NULL) {
            *mf_mds << "RESULT"
                << " type=build_mlf"
                << " text=" << name_text_file
                << " num_threads=" << p
                << " a=" << idx.a;
        }
        time = now();
        std::cout << std::endl << "building M_LF" << std::flush;
    }

    idx._M_LF = std::move(move_data_structure_str<uint_t>(std::move(I_LF),n,{
        .num_threads=p,
        .a=idx.a,
        .log=log,
        .mf=mf_mds,
    }));

    r_ = idx._M_LF.num_intervals();
    idx.r_ = r_;

    if (log) {
        if (mf_mds != NULL) *mf_mds << std::endl;
        if (mf_idx != NULL) {
            *mf_idx << " time_build_mlf=" << time_diff_ns(time,now())
                << " r_=" << r_;
        }
        std::cout << std::endl;
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::build_l__sas() {
    if (log) {
        time = now();
        std::cout << "building L'" << (std::string)(build_locate_support ? " and SA_s" : "") << std::flush;
    }

    if (build_locate_support) {
        no_init_resize(SA_s,r_);
        SA_s[r_-1] = I_Phi[0].second;
    }

    // Simultaneously iterate over the input intervals of M_LF nad the bwt runs to build L'
    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Bwt range start position of thread i_p.
        uint_t b = n_p[i_p];

        // Iteration range start position of thread i_p.
        uint_t b_r = r_p[i_p];

        // Number of runs in thread i_p's section.
        uint_t rp_diff = r_p[i_p+1]-r_p[i_p];

        // Index of the current input interval in M_LF, initially the index of the input interval of M_LF containing b
        uint_t j = bin_search_max_leq<uint_t>(b,0,r_-1,[this](uint_t x){return idx._M_LF.p(x);});
        // Starting position of the next bwt run.
        uint_t l_ = b;

        for (uint_t i=0; i<rp_diff; i++) {
            idx._M_LF.template set_character(j,run_char(i_p,i));
            j++;

            // update l_ to the next run start position
            l_ += run_length(i_p,i);

            // iterate over all input intervals in M_LF within the i-th bwt within thread i_p's section run that have been
            // created by the balancing algorithm
            while (idx._M_LF.p(j) < l_) {
                if (build_locate_support) SA_s[j-1] = n;
                idx._M_LF.template set_character(j,run_char(i_p,i));
                j++;
            }
            
            if (build_locate_support && b_r+i+1 != r) SA_s[j-1] = I_Phi[b_r+i+1].second;
        }
    }

    n_p.clear();
    n_p.shrink_to_fit();

    r_p.clear();
    r_p.shrink_to_fit();

    RLBWT.clear();
    RLBWT.shrink_to_fit();

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_build_l__sas=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::sort_iphi() {
    if (log) {
        time = now();
        std::cout << "sorting I_Phi" << std::flush;
    }

    // Sort I_Phi by the starting positions of its output intervals.
    auto comp_I_Phi = [](std::pair<uint_t,uint_t> p1, std::pair<uint_t,uint_t> p2) {return p1.first < p2.first;};

    // Choose the correct sorting algorithm.
    if (p > 1) {
        ips4o::parallel::sort(I_Phi.begin(),I_Phi.end(),comp_I_Phi);
    } else {
        ips4o::sort(I_Phi.begin(),I_Phi.end(),comp_I_Phi);
    }

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_sort_iphi=" << time_diff_ns(time,now());
        if (mf_mds != NULL) {
            *mf_mds << "RESULT"
                    << " type=build_mphi"
                    << " text=" << name_text_file
                    << " num_threads=" << p
                    << " a=" << idx.a;
        }
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::build_mphi() {
    if (log) {
        time = now();
        std::cout << std::endl << "building M_Phi" << std::flush;
    }

    idx._M_Phi = std::move(move_data_structure<uint_t>(std::move(I_Phi),n,{
        .num_threads=p,
        .a=idx.a,
        .log=log,
        .mf=mf_mds,
    },&pi_mphi));

    r__ = idx._M_Phi.num_intervals();
    idx.r__ = r__;

    if (log) {
        if (mf_mds != NULL) *mf_mds << std::endl;
        if (mf_idx != NULL) {
            *mf_idx << " time_build_mphi=" << time_diff_ns(time,now());
            *mf_idx << " r__=" << r__;
        }
        std::cout << std::endl;
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::build_saphi() {
    time = now();
    if (log) std::cout << "building SA_Phi" << std::flush;

    no_init_resize(pi_,r_);

    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<r_; i++) {
        pi_[i] = i;
    }

    auto comp_pi_ = [this](uint_t i, uint_t j){return SA_s[i] < SA_s[j];};
    if (p > 1) {
        ips4o::parallel::sort(pi_.begin(),pi_.end(),comp_pi_);
    } else {
        ips4o::sort(pi_.begin(),pi_.end(),comp_pi_);
    }

    idx.omega_idx = idx._M_Phi.width_idx();
    idx._SA_Phi = std::move(interleaved_vectors<uint_t>({(uint8_t)(idx.omega_idx/8)}));
    idx._SA_Phi.resize_no_init(r_);

    /* Now we will divide the range [0..n-1] up into p non-overlapping sub-ranges [s[i_p]..s[i_p+1]-1],
    for each i_p in [0..p-1], with 0 = s[0] < s[1] < ... < s[p] = n, where
    s[i_p] = min {s' in [0,n-1], s.t. x[i_p] + u[i_p] - 2 >= i_p * lfloor (r+r'')/p rfloor, where 
                    x[i_p] = min {x' in [0,r''-1], s.t. M_Phi.q(x') >= s'} and
                    u[i_p] = min {u' in [0,r-1], s.t. SA_s[u'] >= s'}
    }.
    By doing so, we ensure that the number of the output intervals of M_Phi starting in the range
    [s[i_p]..s[i_p+1]-1] plus the number of suffix array samples in SA_s lying in the range
    [s[i_p]..s[i_p+1]-1] is lfloor (r+r'')/p rfloor +- 1. This property is useful, because it
    ensures that if with each thread i_p, we simultaneously iterate over those, then each thread
    iterates over almost exactly the same number lfloor (r+r'')/p rfloor +- 1 of entries in M_Phi
    and SA_s combined. This way, we can acheive good load-balancing. Because we do not have to access
    s[0..p] later, we will not store those values in an array. */

    /* [0..p], x[i_p] = min {x' in [0,r''-1], s.t. M_Phi.q(x') >= s'} stores the number of output
    intervals in M_Phi starting before s[i_p]. */
    std::vector<uint_t> x(p+1);
    x[0] = 0;
    x[p] = r__;

    /* [0..p], u[i_p] = min {u' in [0,r-1], s.t. SA_s[u'] >= s'} stores the number of suffix array
    samples in SA_s that are smaller than s[i_p]. */
    std::vector<uint_t> u(p+1);
    u[0] = 0;
    u[p] = r;

    // Compute s[1..p-1], x[1..p-1] and u[1..p-1].
    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // The optimal value i_p * lfloor (r+r'')/p rfloor for s[i_p].
        uint_t o = i_p*((r+r__)/p);

        // Left interval limit of the binary search for s[i_p].
        uint_t l_s;
        // Left interval limit of the binary search for x[i_p].
        uint_t l_x;
        // Left interval limit of the binary search for u[i_p].
        uint_t l_u;
        // Candidate position in the binary search for s[i_p].
        uint_t m_s;
        // Candidate position in the binary search for x[i_p].
        uint_t m_x;
        // Candidate position in the binary search for u[i_p].
        uint_t m_u;
        // Right interval limit of the binary search for s[i_p].
        uint_t r_s;
        // Right interval limit of the binary search for x[i_p].
        uint_t r_x;
        // Right interval limit of the binary search for u[i_p].
        uint_t r_u;

        // Initialize the search range for s[i_p] to [0..n-1].
        l_s = 0;
        r_s = n-1;

        // Perform a binary search over [0..n-1], to find s[i_p].
        while (true) {
            /* Set the Candidate position for s[i_p] to the position in the middle
            between l_s and r_s. */
            m_s = l_s+(r_s-l_s)/2;

            // Find the minimum x' in [0,r''-1], s.t. M_Phi.q(x') >= m_s.
            l_x = 0;
            r_x = r__-1;
            while (l_x != r_x) {
                m_x = l_x+(r_x-l_x)/2;
                if (idx._M_Phi.q(pi_mphi[m_x]) < m_s) {
                    l_x = m_x+1;
                } else {
                    r_x = m_x;
                }
            }

            // Find the minimum u' in [0,r-1], s.t. SA_s[pi'[u']] >= m_s.
            l_u = 0;
            r_u = r-1;
            while (l_u != r_u) {
                m_u = l_u+(r_u-l_u)/2;
                if (SA_s[pi_[m_u]] < m_s) {
                    l_u = m_u+1;
                } else {
                    r_u = m_u;
                }
            }

            /* If l_s = r_s, then l_s is an optimal value for s[i_p] and l_x and
            l_u are valid values for x[i_p] and u[i_p], respectively. */
            if (l_s == r_s) {
                break;
            }

            // Else, adjust the range for the binary search over [0..n-1].
            if (l_x+l_u < o) {
                l_s = m_s + 1;
            } else {
                r_s = m_s;
            }
        }

        // Store l_x and l_u in x[i_p] and u[i_p], respectively.
        x[i_p] = l_x;
        u[i_p] = l_u;
        
        #pragma omp barrier

        // Iteration range start position in the output intervals of M_Phi.
        uint_t i = x[i_p];
        // Iteration range start position in SA_s.
        uint_t j = u[i_p];
        // Iteration range end position + 1 in SA_s.
        uint_t j_ = u[i_p+1];

        // Check if the range, over which the thread i_p has to iterate in SA_s, is empty
        if (j < j_) {
            while (idx._M_Phi.q(pi_mphi[i]) != SA_s[pi_[j]]) {
                i++;
            }

            /* Iterate over SA_s[pi'[j]],SA_s[pi'[j+1]],...,SA_s[pi'[j'-1]] */
            while (j < j_) {
                // Skip the output intervals the balancing algorithm has added to I_Phi
                while (idx._M_Phi.q(pi_mphi[i]) != SA_s[pi_[j]]) {
                    i++;
                }
                
                idx.set_SA_Phi(pi_[j],pi_mphi[i]);

                i++;
                j++;
            }
        }
    }

    /* Since we set SA_s[j] = n for each j-th input interval of M_LF, whiches end position is not the end position of a bwt run,
     * where j \in [0,r'), SA_s[pi'[r]],SA_s[pi'[r+1]],...,SA_s[pi'[r'-1]] = n holds, hence we set SA_Phi[pi'[i]] = r'' for
     * i \in [r,r') to mark that we cannot recover SA_s[M_LF[x+1]-1] = M_Phi.q[SA_Phi[x]] for each x-th input interval of
     * M_LF, whiches end position is not the end position of a bwt run and where x \in [0,r'). */
    #pragma omp parallel for num_threads(p)
    for (uint64_t i=r; i<r_; i++) {
        idx.set_SA_Phi(pi_[i],r__);
    }

    x.clear();
    x.shrink_to_fit();

    u.clear();
    u.shrink_to_fit();

    pi_mphi.clear();
    pi_mphi.shrink_to_fit();

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_build_saphi=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::build_de() {
    if (build_locate_support) {
        idx.p_r = std::min<uint_t>(256,std::max<uint_t>(1,r/100));
    } else {
        idx.p_r = 1;
    }
    
    idx._D_e.resize(idx.p_r-1);
    
    #pragma omp parallel for num_threads(p)
    for (uint16_t i=0; i<idx.p_r-1; i++) {
        uint_t x = bin_search_min_geq<uint_t>((i+1)*((n-1)/idx.p_r),0,r-1,[this](uint_t x){return (((int64_t)SA_s[pi_[x]])-1)%n;});
        idx._D_e[i] = std::make_pair(pi_[x],(uint_t)((((int64_t)SA_s[pi_[x]])-1)%n));
    }

    SA_s.clear();
    SA_s.shrink_to_fit();

    pi_.clear();
    pi_.shrink_to_fit();
}

template <typename uint_t>
void move_r<uint_t>::construction::build_rsl_() {
    if (log) {
        time = now();
        std::cout << "building RS_L'" << std::flush;
    }
    
    idx._RS_L_ = std::move(string_rank_select_support<uint_t>([this](uint_t i){return idx.L_(i);},0,r_-1,p));

    if (log) {
        if (mf_idx != NULL) *mf_idx << " time_build_rsl_=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}