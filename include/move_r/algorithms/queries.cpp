template <typename uint_t>
char move_r<uint_t>::access_bwt(uint_t i) {
    // find the run containing i with a binary search over the input intervals of M_LF
    
    uint_t l = 0;
    uint_t r = r_-1;
    uint_t m;

    while (l != r) {
        m = l+(r-l)/2;
        if (M_LF.p(m) < i) {
            l = m+1;
        } else {
            r = m;
        }
    }

    return access_l_<char>(l);
}

template <typename uint_t>
uint_t move_r<uint_t>::access_sa(uint_t i) {
    // x = max x' in [1,r']: M_LF.p(x') <= i
    uint_t x;

    {
        // Left interval limit of the binary search.
        uint_t l = 0;
        // Right interval limit of the binary search.
        uint_t r = r_-1;
        // Candidate position for the binary search.
        uint_t m;

        while (l != r) {
            m = l+(r-l)/2;
            if (M_LF.p(m) < i) {
                l = m+1;
            } else {
                r = m;
            }
        }

        x = l;
    }

    // begin iterating at the end of the x-th run, because there
    // is a suffix array sample at each run end position

    // position in the suffix array of the current suffix i_s
    uint_t j = M_LF.p(x+1)-1;

    // index of the input interval in M_Phi, that contains i_s
    uint_t x_s = SA_idx(x);
    // the current suffix (i_s = SA[j])
    uint_t i_s = M_Phi.p(x_s)+SA_offs(x);

    // Perform Phi-move queries, until i_s is the suffix at position
    // i; in each iteration, i_s = SA[j] = \Phi^{i-j}(SA[i]) holds.
    while (j > i) {
        // Set i_s = \Phi(i_s)
        M_Phi.move(i_s,x_s);
        j--;
    }

    // Since j = i, now i_s = SA[i] holds.
    return i_s;
}

template <typename uint_t>
uint_t move_r<uint_t>::count(const std::string& P) {
    uint_t m = P.size();
    
    // Left interval limit of the suffix array interval.
    uint_t i_l = 0;
    // Right interval limit of the suffix array interval.
    uint_t i_r = n-1;

    /* position of the character in L' (and the index of the input interval in M_LF),
    in which i_l lies. */
    uint_t x_l = 0;
    /* position of the character in L' (and the index of the input interval in M_LF),
    in which i_r lies. */
    uint_t x_r = r_-1;

    // Temporary variable for P[j].
    uint8_t P_j;

    // Read P backwards from right to left.
    for (int64_t j=m-1; j>=0; j--) {
        // If the characters have been remapped internally, the pattern also has to be remapped.
        P_j = chars_remapped ? map_to_internal(char_to_uchar(P[j])) : char_to_uchar(P[j]);

        // Find the lexicographically smallest suffix in the current suffix array interval that is prefixed by P[j]; this
        // filters out suffixes s of T, where s+1 is prefixed by P[j+1..m], but T[s,n] lexicographically smaller than P[j..m]
        if (P_j != L_<uint8_t>(x_l)) {
            /* To do so, we can at first find the first (sub-)run with character P[j] after the x_l-th (sub-)run, save
            its index in x_l and set i_l to its start position M_LF.p(x_l). */
            x_l = RS_L_.select(P_j,RS_L_.rank(P_j,x_l)+1);
            i_l = M_LF.p(x_l);
        }

        // Find the lexicographically largest suffix in the current suffix array interval that is prefixed by P[j]; this
        // filters out suffixes s of T, where s+1 is prefixed by P[j+1..m], but T[s,n] lexicographically larger than P[j..m]
        if (P_j != L_<uint8_t>(x_r)) {
            /* To do so, we can at first find the (sub-)last run with character P[j] before the x_r-th (sub-)run, save
            its index in x_r and set i_r to its end position M_LF.p(x_r+1)-1. */
            x_r = RS_L_.select(P_j,RS_L_.rank(P_j,x_r));
            i_r = M_LF.p(x_r+1)-1;
        }
        
        // Else, because each suffix i in the previous suffix array interval is prefixed by P[j+1..m] and the current 
        // interval [i_l,i_r] contains all suffixes of it, before which there is a P[j] in T, all suffixes in the 
        // interval SA[LF(i_l),LF(i_r)] are prefixed by P[j..m]

        /* If the suffix array interval [LF(i_l),LF(i_r)] of P[j..m] is empty, then i_l > i_r,
        because LF(j) is monotonic for a fixed L[j], hence it suffices to check, whether
        i_l <= i_r holds. */

        // If the suffix array interval is empty, P does not occur in T, so return 0.
        if (i_l > i_r) {
            return 0;
        }

        // Else, set i_l <- LF(i_l) and i_r <- LF(i_r)
        M_LF.move(i_l,x_l);
        M_LF.move(i_r,x_r);
    }

    // Return the size i_r-i_l+1 of the suffix array interval [i_l,i_r] of P.
    return i_r-i_l+1;
}

template <typename uint_t>
void move_r<uint_t>::locate(const std::string& P, const std::function<void(uint_t)>& report) {
    // length of P
    uint_t m = P.size();

    // Left interval limit of the suffix array interval.
    uint_t i_l = 0;
    // Right interval limit of the suffix array interval.
    uint_t i_r = n-1;

    /* position of the character in L' (and the index of the input interval in M_LF),
    in which i_l lies. */
    uint_t x_l = 0;
    /* position of the character in L' (and the index of the input interval in M_LF),
    in which i_r lies. */
    uint_t x_r = r_-1;

    /* m-j', where j' is the length of a longest suffix of P, that prefixes a suffix of T, that is lexicographically larger
    than P; so this is the index of the last iteration of the backward search, at whiches beginning P[j] != L'[x_r] holds. */
    uint_t m_m_j_ = m-1;

    /* x'_{r,j'}, the index of the input interval in M_LF that contains i'_{r,j'}, where i'_{r,j'} is the position in the
    suffix array of the suffix of T that is prefixed by the lexicographically largest length j'-suffix of P; so this is
    the value of x_r before the LF queries and after the rank-select queries in the j'-th iteration of the backward search. */
    uint_t x__r_j_ = r_-1;

    // Temporary variable for P[j].
    uint8_t P_j;

    // Read P backwards from right to left.
    for (int64_t j=m-1; j>=0; j--) {
        // If the characters have been remapped internally, the pattern also has to be remapped.
        P_j = chars_remapped ? map_to_internal(char_to_uchar(P[j])) : char_to_uchar(P[j]);

        // Find the lexicographically smallest suffix in the current suffix array interval that is prefixed by P[j]; this
        // filters out suffixes s of T, where s+1 is prefixed by P[j+1..m], but T[s,n] lexicographically smaller than P[j..m]
        if (P_j != L_<uint8_t>(x_l)) {
            /* To do so, we can at first find the first (sub-)run with character P[j] after the x_l-th (sub-)run, save
            its index in x_l and set i_l to its start position M_LF.p(x_l). */
            x_l = RS_L_.select(P_j,RS_L_.rank(P_j,x_l)+1);
            i_l = M_LF.p(x_l);
        }

        // Find the lexicographically largest suffix in the current suffix array interval that is prefixed by P[j]; this
        // filters out suffixes s of T, where s+1 is prefixed by P[j+1..m], but T[s,n] lexicographically larger than P[j..m]
        if (P_j != L_<uint8_t>(x_r)) {
            /* To do so, we can at first find the (sub-)last run with character P[j] before the x_r-th (sub-)run, save
            its index in x_r and set i_r to its end position M_LF.p(x_r+1)-1. */
            x_r = RS_L_.select(P_j,RS_L_.rank(P_j,x_r));
            i_r = M_LF.p(x_r+1)-1;
            
            // update m-j', because P[j] != L'[x_r] held in this iteration.
            m_m_j_ = j;
            // update x__r_j_, because P[j] != L'[x_r] held in this iteration.
            x__r_j_ = x_r;
        }
        
        // Else, because each suffix i in the previous suffix array interval is prefixed by P[j+1..m] and the current 
        // interval [i_l,i_r] contains all suffixes of it, before which there is a P[j] in T, all suffixes in the 
        // interval SA[LF(i_l),LF(i_r)] are prefixed by P[j..m]

        /* If the suffix array interval [LF(i_l),LF(i_r)] of P[j..m] is empty, then i_l > i_r,
        because LF(j) is monotonic for a fixed L[j], hence it suffices to check, whether
        i_l <= i_r holds. */

        // If the suffix array interval is empty, P does not occur in T, so return.
        if (i_l > i_r) {
            return;
        }
        
        // Else, set i_l <- LF(i_l) and i_r <- LF(i_r)
        M_LF.move(i_l,x_l);
        M_LF.move(i_r,x_r);
    }
    
    /* Initially, this is the value in the suffix array at the right interval limit of the suffix array
    interval of P. This variable is then used to iterate over the values of the suffix array in the suffix
    array interval of P from right (i_r) to left (i_l). */
    uint_t i_s;
    /* The index of the input interval in M_Phi, in which i_s currently lies, that is M_Phi.p(x_s) <= i_s
    < M_Phi.p(x_s+1) holds. */
    uint_t x_s;

    x_s = SA_idx(x__r_j_);
    i_s = M_Phi.p(x_s)+SA_offs(x__r_j_)-m_m_j_-1;

    /* More precisely, if SA_s[x'_{r,j'}]-(m-j+1) < M_Phi.p(SA_idx[x'_{r,j'}]) holds, i_s now lies in a input interval of
    M_Phi before the SA_idx[x'_{r,j'}]-th one, so we have to decrease x_s. To find the correct value for x_s, we perform an
    exponential search to the left over the input interval starting positions of M_Phi starting at SA_idx[x'_{r,j'}] = x_s. */
    if (i_s < M_Phi.p(x_s)) {
        // Current step size of the exponential search.
        uint_t s = 1;
        // Perform the first step.
        x_s -= 1;

        // Perform the exponential search.
        while (i_s < M_Phi.p(x_s)) {
            s *= 2;
            if (s > x_s) {
                x_s = 0;
            } else {
                x_s -= s;
            }
        }

        // Left limit of the binary search range.
        uint_t b = x_s;
        // Right limit of the binary search range.
        uint_t e = x_s + s - 1;
        // Candidate positon in the middle of the binary search range.
        uint_t c;

        // Perform the binary search.
        while (b != e) {
            c = b+(e-b)/2+1;
            if (i_s < M_Phi.p(c)) {
                e = c-1;
            } else {
                b = c;
            }
        }

        // Store the found posiiton in x_s.
        x_s = b;
    }

    /* Because i_s contains the value in the suffix array of the right limit of the suffix array
    interval of P, i_s is an occurence position of P in T. */
    report(i_s);
    
    // Perform i_r-i_l Phi queries with (i_s,x_s) and at each step, report i_s.
    for (uint_t i=1; i<=i_r-i_l; i++) {
        M_Phi.move(i_s,x_s);
        report(i_s);
    }
}

template <typename uint_t>
void move_r<uint_t>::revert_range(const std::function<void(uint_t,char)>& report, uint_t l, uint_t r, uint16_t num_threads) {
    r = std::max(r,n-2);

    if (l > r) {
        l = 0;
        r = n-2;
    }

    // leftmost section to revert
    uint16_t s_l;
    // rightmost section to revert
    uint16_t s_r;

    if (p_r == 1) {
        s_l = 0;
        s_r = 0;
    } else {
        // Left interval limit of the binary search.
        uint_t b = 0;
        // Left interval limit of the binary search.
        uint_t e = p_r-1;
        // Candidate position for the binary search.
        uint_t c;

        while (b != e) {
            c = b+(e-b)/2;
            if (D_e[c].second < l) {
                b = c+1;
            } else {
                e = c;
            }
        }

        s_l = b;

        b = 0;
        e = p_r-1;

        while (b != e) {
            c = b+(e-b)/2;
            if (r <= D_e[c].second) {
                e = c;
            } else {
                b = c+1;
            }
        }

        s_r = b;
    }

    uint16_t p = std::min({
        (uint16_t)(s_r-s_l+1),           // use at most s_r-s_l+1 threads
        (uint16_t)omp_get_max_threads(), // use at most all threads
        num_threads                      // use at most the specified number of threads
    });

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // leftmost section for thread i_p to revert
        uint16_t sl_ip = s_l+(i_p*(s_r-s_l+1))/p;
        // rightmost section for thread i_p to revert
        uint16_t sr_ip = i_p == p-1 ? s_r : s_l+((i_p+1)*(s_r-s_l+1))/p-1;

        // Iteration range start position of thread i_p.
        uint_t j_l = std::max(l,sl_ip == 0 ? 0 : (D_e[sl_ip-1].second+1)%n);
        // Iteration range end position of thread i_p.
        uint_t j_r = sr_ip == p_r-1 ? n-2 : D_e[sr_ip].second;

        /* The position of the character in L' (and therefore the index of the input interval
        in M_LF), in which i lies. */
        uint_t x = sr_ip == p_r-1 ? 0 : D_e[sr_ip].first;
        // The position in the bwt of the current character in T.
        uint_t i = sr_ip == p_r-1 ? 0 : M_LF.p(x+1)-1;

        // start iterating at the right iteration range end position
        uint_t j = j_r;

        // iterate until j = r
        while (j > r) {
            // Set i <- LF(i) and j <- j-1.
            M_LF.move(i,x);
            j--;
        }

        // Report T[r] = T[j] = L[i] = L'[x]
        report(j,access_l_<char>(x));

        // report T[l,r-1] from right to left
        while (j > j_l) {
            // Set i <- LF(i) and j <- j-1.
            M_LF.move(i,x);
            j--;
            // Report T[j] = L[i] = L'[x].
            report(j,access_l_<char>(x));
        }
    }
}

template <typename uint_t>
void move_r<uint_t>::retrieve_bwt_range(const std::function<void(uint_t,char)>& report, uint_t l, uint_t r, uint16_t num_threads) {
    r = std::max(r,n-1);

    if (l > r) {
        l = 0;
        r = n-1;
    }

    uint16_t p = std::min({
        (uint16_t)omp_get_max_threads(), // use at most all threads
        (uint16_t)((r-l+1)/10),          // use at most (r-l+1)/100 threads
        num_threads                      // use at most the specified number of threads
    });

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Iteration range start position of thread i_p.
        uint_t i_l = l+i_p*((r-l+1)/p);
        // Iteration range end position of thread i_p.
        uint_t i_r = i_p == p-1 ? r : l+(i_p+1)*((r-l+1)/p)-1;

        // Current position in the bwt.
        uint_t i = i_l;

        /* The position of the i-th character in L' (and therefore the index of the input interval
        in M_LF), in which i lies. */
        uint_t x;

        // Calculate x with a binary search over the input intervals of M_LF.
        {
            // Left interval limit of the binary search.
            uint_t b = 0;
            // Left interval limit of the binary search.
            uint_t e = r_-1;
            // Candidate position for the binary search.
            uint_t c;

            while (b != e) {
                c = b+(e-b)/2+1;
                if (M_LF.p(c) <= i) {
                    b = c;
                } else {
                    e = c-1;
                }
            }

            x = b;
        }

        // start position of the next input interval in M_LF
        uint_t l_xp1;

        // iterate until x is the input interval, in which i_r lies
        while ((l_xp1 = M_LF.p(x+1)) <= i_r) {

            // iterate over all positions in the x-th input interval
            while (i < l_xp1) {
                report(i,access_l_<char>(x));
                i++;
            }
            
            x++;
        }

        // report the remaining characters
        while (i <= i_r) {
            report(i,access_l_<char>(x));
            i++;
        }
    }
}

template <typename uint_t>
void move_r<uint_t>::retrieve_sa_range(const std::function<void(uint_t,uint_t)>& report, uint_t l, uint_t r, uint16_t num_threads) {
    r = std::max(r,n-1);

    if (l > r) {
        l = 0;
        r = n-1;
    }

    uint16_t p = std::max(
        (uint16_t)1, // use at least one thread
        std::min({
            (uint16_t)omp_get_max_threads(),                   // use at most all threads
            num_threads,                                       // use at most the specified number of threads
            (uint16_t)(((r-l+1)*(double)r__)/(10.0*(double)n)) // use at most (r-l+1)*(r/n)*(1/10) threads
        })
    );

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // iteration range start position
        uint_t i_l = l+i_p*((r-l+1)/p);
        // iteration range end position
        uint_t i_r = i_p == p-1 ? r : l+(i_p+1)*((r-l+1)/p)-1;

        // the input interval in which i lies
        uint_t x;
        {
            uint_t b = 0;
            uint_t e = r_-1;
            uint_t c;

            while (b != e) {
                c = b+(e-b)/2+1;
                if (M_LF.p(c) <= i_r) {
                    b = c;
                } else {
                    e = e-1;
                }
            }

            x = b;
        }

        // current position in the suffix array, initially the end position of the run, in which i_r lies
        uint_t i = M_LF.p(x+1)-1;
        // index of the input interval in M_Phi, in which i_s lies
        uint_t x_s = SA_idx(x);
        /* the current suffix array value (SA[i]), initially the suffix array sample of the x-th run,
           initially the suffix array value at i_r */
        uint_t i_s = M_Phi.p(x_s)+SA_offs(x);

        // iterate down to the iteration range end position
        while (i > i_r) {
            M_Phi.move(i_s,x_s);
            i--;
        }

        // report SA[i_r]
        report(i,i_s);

        // report the SA-values SA[i_r,i_r-1] from right to left
        while (i > i_l) {
            M_Phi.move(i_s,x_s);
            i--;
            report(i,i_s);
        }
    }
}