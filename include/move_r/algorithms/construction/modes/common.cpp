#include <ips4o.hpp>

template <typename uint_t>
void move_r<uint_t>::construction::preprocess_t(bool in_memory, bool map_t, std::ifstream* t_file) {
    if (log) std::cout << "preprocessing T" << std::flush;

    // contains_uchar_thr[i_p][c] = true <=> thread i_p found the character c in its section of T[0..n-2].
    std::vector<std::vector<uint8_t>> contains_uchar_thr(p,std::vector<uint8_t>(256,0));

    if (in_memory) {
        // Iterate over T[0..n-2] and report the occurence of each found character in T[0..n-2] in contains_uchar_thr.
        #pragma omp parallel for num_threads(p)
        for (uint64_t i=0; i<n-1; i++) {
            contains_uchar_thr[omp_get_thread_num()][char_to_uchar(T[i])] = 1;
        }
    } else {
        uint_t t_buffer_size = std::max((uint_t)1,n/500);
        std::string T_buffer;
        no_init_resize(T_buffer,t_buffer_size);
        uint_t current_buffer_size;
        uint_t n_ = n-1;

        while (n_ > 0) {
            current_buffer_size = std::min(n_,t_buffer_size);
            read_from_file(*t_file,T_buffer.c_str(),current_buffer_size);
            
            // Iterate over T[0..n-2] and report the occurence of each found character in T[0..n-2] in contains_uchar_thr.
            #pragma omp parallel for num_threads(p)
            for (uint64_t i=0; i<current_buffer_size; i++) {
                contains_uchar_thr[omp_get_thread_num()][char_to_uchar(T_buffer[i])] = 1;
            }

            n_ -= current_buffer_size;
        }
    }

    // contains_uchar[c] = 1 <=> c in T[0..n-2].
    contains_uchar.resize(256,0);

    /* Combine the results of each thread's sections in contains_uchar_thr[0..i_p-1][0..255]
    into contains_uchar[0..255]. */
    for (uint16_t i=0; i<256; i++) {
        for (uint16_t i_p=0; i_p<p; i_p++) {
            if (contains_uchar_thr[i_p][i] == 1) {
                contains_uchar[i] = 1;
                break;
            }
        }
    }

    contains_uchar_thr.clear();
    contains_uchar_thr.shrink_to_fit();

    // The number of distinct characters in T[0..n-2].
    sigma = 1;

    // Count the number of ones in contains_uchar[0..255] in sigma.
    for (uint16_t i=0; i<256; i++) {
        if (contains_uchar[i] == 1) {
            sigma++;
        }
    }
    
    idx.sigma = sigma;
    bool contains_invalid_char = false;

    for (uint8_t i=0; i<min_valid_char; i++) {
        if (contains_uchar[i] == 1) {
            contains_invalid_char = true;
            break;
        }
    }

    // If an invalid character, we have to remap the characters in T[0..n-2].
    if (contains_invalid_char) {
        if (sigma > 253) {
            /* If T[0..n-2] contains more than 253 distinct characters, we cannot remap them into the 
            range [0..255] without using 0 or 1, hence we cannot build an index for T. */
            std::cout << "Error: the input contains more than 253 distinct characters" << std::endl;
        }

        idx.chars_remapped = true;
        
        // build the mapping function map_char, that remaps the characters of T, so it does
        // not contain 0 or 1; also build its inverse function unmap_char

        idx.map_char.resize(256,0);
        idx.unmap_char.resize(256,0);

        /* To preserve the order among characters in T[0..n-2], we start by mapping smallest
           character in T[0..n-2] to 2, the second smallest to 3, ... . */

        // The character, to map the currently next largest character in T[0..n-2] to.
        uint16_t j = min_valid_char;

        for (uint16_t i=0; i<256; i++) {
            if (contains_uchar[i] == 1) {
                idx.map_char[i] = j;
                idx.unmap_char[j] = i;
                j++;
            }
        }

        // unmark the invalid chars
        for (uint8_t i=0; i<min_valid_char; i++) {
            contains_uchar[i] = 0;
        }

        // Since the largest character has been mapped to j-1, T contains c > min_valid_char exactly if c < j
        for (uint16_t i=min_valid_char; i<256; i++) {
            contains_uchar[i] = i < j;
        }

        // Apply map_char to T.
        if (map_t) {
            #pragma omp parallel for num_threads(p)
            for (uint64_t i=0; i<n-1; i++) {
                T[i] = idx.map_to_internal(T[i]);
            }
        }
    }

    contains_uchar[terminator] = 1;
    
    if (!build_count_support) {
        contains_uchar.clear();
        contains_uchar.shrink_to_fit();
    }

    if (log) {
        if (measurement_file_index != NULL) *measurement_file_index << " time_preprocess_t=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::build_rsl_() {
    if (log) {
        time = now();
        std::cout << "building RS_L'" << std::flush;
    }

    
    std::vector<char> chars;

    if (!build_from_sa_and_l) {
        for (uint16_t i=0; i<256; i++) {
            if (contains_uchar[i] == 1) {
                chars.emplace_back(uchar_to_char((uint8_t)i));
            }
        }
    }

    contains_uchar.clear();
    contains_uchar.shrink_to_fit();
    
    idx.RS_L_ = std::move(string_rank_select_support<uint_t>([this](uint_t i){return idx.L_<char>(i);},0,r_-1,chars,p));

    chars.clear();
    chars.shrink_to_fit();

    if (log) {
        if (measurement_file_index != NULL) *measurement_file_index << " time_build_rsl_=" << time_diff_ns(time,now());
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::process_c_array() {
    C[p].resize(256,0);

    /* Now, C[i_p][c] is the number of occurences of c in L[b..e], where [b..e] is the range of the
    thread i_p in [0..p-1]. Also, we have C[p][0..255] = 0. */

    /* We want to have C[i_p][c] = rank(L,c,b-1), where b is the iteration range start position of
    thread i_p in [0..p-1]. Also, we want C[p][0..255] to be the C-array, that is C[p][c] stores
    the number of occurences of all smaller characters c' < c in L[0..n-1], for c in [0..255]. */

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
void move_r<uint_t>::construction::build_mlf() {
    if (log) {
        if (measurement_file_move_data_structures != NULL) {
            *measurement_file_move_data_structures << "RESULT"
                << " type=build_mlf"
                << " text=" << name_textfile
                << " p=" << p
                << " a=" << a;
        }
        time = now();
        std::cout << std::endl << "building M_LF" << std::flush;
    }

    idx.M_LF = std::move(move_data_structure_lf<uint_t>(I_LF,n,p,a,log,measurement_file_move_data_structures));
    r_ = idx.M_LF.num_intervals();
    idx.r_ = r_;

    if (log) {
        if (measurement_file_move_data_structures != NULL) *measurement_file_move_data_structures << std::endl;
        if (measurement_file_index != NULL) {
            *measurement_file_index << " time_build_mlf=" << time_diff_ns(time,now())
                << " r_=" << r_;
        }
        std::cout << std::endl;
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
        if (measurement_file_move_data_structures != NULL) {
            *measurement_file_move_data_structures << "RESULT"
                    << " type=build_mphi"
                    << " text=" << name_textfile
                    << " p=" << p
                    << " a=" << a;
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

    idx.M_Phi = std::move(move_data_structure_phi<uint_t>(I_Phi,n,p,a,log,measurement_file_move_data_structures));
    r__ = idx.M_Phi.num_intervals();
    idx.r__ = r__;

    if (log) {
        if (measurement_file_move_data_structures != NULL) *measurement_file_move_data_structures << std::endl;
        if (measurement_file_index != NULL) {
            *measurement_file_index << " time_build_mphi=" << time_diff_ns(time,now());
            *measurement_file_index << " r__=" << r__;
        }
        std::cout << std::endl;
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::build_saidxoffs(uint_t r_) {
    time = now();
    if (log) std::cout << "building SA_offs and SA_idx" << std::flush;

    (*reinterpret_cast<std::vector<no_init<uint_t>>*>(&pi_)).resize(this->r_);

    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<this->r_; i++) {
        pi_[i] = i;
    }

    auto comp_pi_ = [this](uint_t i, uint_t j){return SA_s[i] < SA_s[j];};
    if (p > 1) {
        ips4o::parallel::sort(pi_.begin(),pi_.end(),comp_pi_);
    } else {
        ips4o::sort(pi_.begin(),pi_.end(),comp_pi_);
    }

    omega_idx = std::max((uint8_t)8,(uint8_t)(std::ceil(std::log2(r__)/(double)8)*8));
    omega_offs = idx.M_Phi.width_offs();
    
    idx.omega_idx = omega_idx;
    idx.omega_offs = omega_offs;

    idx.SA_idxoffs = std::move(interleaved_vectors<uint_t>({(uint8_t)(omega_idx/8),(uint8_t)(omega_offs/8)},this->r_,false));

    /* Now we will divide the range [0..n-1] up into p non-overlapping sub-ranges [s[i_p]..s[i_p+1]-1],
    for each i_p in [0..p-1], with 0 = s[0] < s[1] < ... < s[p] = n, where
    s[i_p] = min {s' in [0,n-1], s.t. x[i_p] + u[i_p] - 2 >= i_p * lfloor (r'+r'')/p rfloor, where 
                    x[i_p] = min {x' in [0,r''-1], s.t. M_Phi.p(x') >= s'} and
                    u[i_p] = min {u' in [0,r'-1], s.t. SA_s[u'] >= s'}
    }.
    By doing so, we ensure, that the number of the input intervals of M_Phi starting in the range
    [s[i_p]..s[i_p+1]-1] plus the number of suffix array samples in SA_s lying in the range
    [s[i_p]..s[i_p+1]-1] is lfloor (r'+r'')/p rfloor +- 1. This property is useful, because it
    ensures that if with each thread i_p, we simultaneously iterate over those, then each thread
    iterates over almost exactly the same number lfloor (r'+r'')/p rfloor +- 1 of entries in M_Phi
    and SA_s combined. This way, we can acheive good load-balancing. Because we do not have to access
    s[0..p] later, we will not store those values in an array. */

    /* [0..p], x[i_p] = min {x' in [0,r''-1], s.t. M_Phi.p(x') >= s'} stores the number of input
    intervals in M_Phi starting before s[i_p]. */
    std::vector<uint_t> x(p+1);
    x[0] = 0;
    x[p] = r__;

    /* [0..p], u[i_p] = min {u' in [0,r'-1], s.t. SA_s[u'] >= s'} stores the number of suffix array
    samples in SA_s that are smaller than s[i_p]. */
    std::vector<uint_t> u(p+1);
    u[0] = 0;
    u[p] = r_;

    // Compute s[1..p-1], x[1..p-1] and u[1..p-1].
    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // The optimal value i_p * lfloor (r'+r'')/p rfloor for s[i_p].
        uint_t o = i_p*((r_+r__)/p);

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

            // Find the minimum x' in [0,r''-1], s.t. M_Phi.p(x') >= m_s.
            l_x = 0;
            r_x = r__-1;
            while (l_x != r_x) {
                m_x = l_x+(r_x-l_x)/2;
                if (idx.M_Phi.p(m_x) < m_s) {
                    l_x = m_x+1;
                } else {
                    r_x = m_x;
                }
            }

            // Find the minimum u' in [0,r'-1], s.t. SA_s[pi'[u']] >= m_s.
            l_u = 0;
            r_u = r_-1;
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

        /* Check if the range [u[i_p]..u[i_p+1]-1], over which the thread i_p
        has to iterate in SA_s, is empty. */
        if (u[i_p+1] > u[i_p]) {
            // Iteration range start position in M_Phi.
            uint_t i = x[i_p];
            // Iteration range start position in SA_s.
            uint_t j = u[i_p];

            // Iteration range end position in M_Phi.
            uint_t i_ = x[i_p+1]-1;
            // Iteration range end position in SA_s.
            uint_t j_ = u[i_p+1]-1;

            /* Check if the first suffix array sample lies before x[i_p]-th 
            input interval of M_Phi. */
            if (SA_s[pi_[j]] < idx.M_Phi.p(i)) {
                i--;
            }

            /* Iterate until one of the iteration end positions i_ and j_ has
            been reached. */
            while (i <= i_ && j <= j_) {

                /* Iterate over the suffix array samples that lie in the current
                i-th input interval of M_Phi. */
                while (j <= j_ && SA_s[pi_[j]] < idx.M_Phi.p(i+1)) {
                    
                    /* Because each of those j-th largest suffix array samples lie in the i-th
                    input interval of M_Phi, we can set SA_idx[pi'[j]] = i for each of them. */
                    idx.set_SA_idx(pi_[j],i);
                    idx.set_SA_offs(pi_[j],SA_s[pi_[j]]-idx.M_Phi.p(i));

                    j++;
                }

                i++;
            };
        }
    }

    x.clear();
    x.shrink_to_fit();

    u.clear();
    u.shrink_to_fit();

    SA_s.clear();
    SA_s.shrink_to_fit();

    if (log) {
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_r<uint_t>::construction::build_de(uint_t r_) {
    idx.D_e.resize(p_r-1);
    
    #pragma omp parallel for num_threads(p)
    for (uint16_t i=0; i<p_r-1; i++) {
        uint_t opt = (i+1)*((n-1)/p_r);

        // Left interval limit of the binary search.
        uint_t b = 0;
        // Left interval limit of the binary search.
        uint_t e = r_-1;
        // Candidate position for the binary search.
        uint_t c;

        while (b != e) {
            c = b+(e-b)/2;
            if (opt <= (((int64_t)idx.SA_s(pi_[c]))-1)%n) {
                e = c;
            } else {
                b = c+1;
            }
        }

        idx.D_e[i] = std::make_pair(pi_[b],(uint_t)((((int64_t)idx.SA_s(pi_[b]))-1)%n));
    }

    pi_.clear();
    pi_.shrink_to_fit();
}