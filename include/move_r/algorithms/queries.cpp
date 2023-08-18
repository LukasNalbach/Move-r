template <typename uint_t>
void move_r<uint_t>::setup_phi_move_pair(uint_t x, uint_t& s, uint_t& s_) {
    // the index of the pair in M_Phi, that creates the output interval with start position s = SA[l'_{x+1}-1]
    uint_t x_s_ = SA_idx_phi[x];

    // set s_ to the index of the input interval in M_Phi, that contains s
    s_ = M_Phi.idx(x_s_);
    
    // compute s
    s = M_Phi.p(s_)+M_Phi.offs(x_s_);
}

template <typename uint_t>
char move_r<uint_t>::access_bwt(uint_t i) {
    // find the run containing i with a binary search over the input intervals of M_LF
    return access_l_(bin_search_max_leq<uint_t>(i,0,r_-1,[this](uint_t p){return M_LF.p(p);}));
}

template <typename uint_t>
uint_t move_r<uint_t>::access_sa(uint_t i) {
    // x = max x' in [1,r']: M_LF.p(x') <= i
    uint_t x = bin_search_max_leq<uint_t>(i,0,r_-1,[this](uint_t x_){return M_LF.p(x_);});

    if (i == M_LF.p(x)) return M_Phi.p(SA_idx_phi[x]);

    // begin iterating at the end of the x-th run, because there
    // is a suffix array sample at each run end position

    // position in the suffix array of the current suffix s
    uint_t j = M_LF.p(x+1)-1;

    // index of the input interval in M_Phi, that contains s
    uint_t s_;
    // the current suffix (s = SA[j])
    uint_t s;

    setup_phi_move_pair(x,s,s_);

    // Perform Phi-move queries, until s is the suffix at position
    // i; in each iteration, s = SA[j] = \Phi^{i-j}(SA[i]) holds.
    while (j > i) {
        // Set s = \Phi(s)
        M_Phi.move(s,s_);
        j--;
    }

    // Since j = i, now s = SA[i] holds.
    return s;
}

template <typename uint_t>
uint_t move_r<uint_t>::count(const std::string& P) {
    uint_t m = P.size();
    
    // Left interval limit of the suffix array interval.
    uint_t b = 0;
    // Right interval limit of the suffix array interval.
    uint_t e = n-1;

    // index of the input interval in M_LF containing b */
    uint_t b_ = 0;
    // index of the input interval in M_LF containing e */
    uint_t e_ = r_-1;

    // Temporary variable for P[j].
    char P_j;

    // Read P backwards from right to left.
    for (int64_t j=m-1; j>=0; j--) {
        // If the characters have been remapped internally, the pattern also has to be remapped.
        P_j = chars_remapped ? map_to_internal(P[j]) : P[j];

        if (!RS_L_.contains_character(P_j)) return 0;

        // Find the lexicographically smallest suffix in the current suffix array interval that is prefixed by P[j]
        if (P_j != L_(b_)) {
            /* To do so, we can at first find the first (sub-)run with character P[j] after the b_-th (sub-)run, save
            its index in b_ and set b to its start position M_LF.p(b_). */
            b_ = RS_L_.rank(P_j,b_);
            if (b_ == RS_L_.num_occurrences(P_j)) return 0;
            b_ = RS_L_.select(P_j,b_+1);
            b = M_LF.p(b_);
        }

        // Find the lexicographically largest suffix in the current suffix array interval that is prefixed by P[j]
        if (P_j != L_(e_)) {
            /* To do so, we can at first find the (sub-)last run with character P[j] before the e_-th (sub-)run, save
            its index in e_ and set e to its end position M_LF.p(e_+1)-1. */
            e_ = RS_L_.rank(P_j,e_);
            if (e_ == 0) return 0;
            e_ = RS_L_.select(P_j,e_);
            e = M_LF.p(e_+1)-1;
        }
        
        // Else, because each suffix i in the previous suffix array interval starts with P[j+1..m] and the current 
        // interval [b,e] contains all suffixes of it, before which there is a P[j] in T, all suffixes in the 
        // interval SA[LF(b),LF(e)] start with P[j..m]

        /* If the suffix array interval [LF(b),LF(e)] of P[j..m] is empty, then b > e,
        because LF(j) is monotonic for a fixed L[j], hence it suffices to check, whether
        b <= e holds. */

        // If the suffix array interval is empty, P does not occur in T, so return 0.
        if (b > e) return 0;
        
        // Else, set b <- LF(b) and e <- LF(e)
        M_LF.move(b,b_);
        M_LF.move(e,e_);
    }

    // Return the size e-b+1 of the suffix array interval [b,e] of P.
    return e-b+1;
}

template <typename uint_t>
void move_r<uint_t>::locate(const std::string& P, const std::function<void(uint_t)>& report) {
    // length of P
    uint_t m = P.size();

    // Left interval limit of the suffix array interval.
    uint_t b = 0;
    // Right interval limit of the suffix array interval.
    uint_t e = n-1;

    /* index of the input interval in M_LF containing b. */
    uint_t b_ = 0;
    /* index of the input interval in M_LF containing e. */
    uint_t e_ = r_-1;

    uint_t y; // y
    uint_t hat_e_ap_y; // \hat{e}'_y

    // Temporary variable for P[j].
    char P_j;

    // Read P backwards from right to left.
    for (int64_t j=m-1; j>=0; j--) {
        // If the characters have been remapped internally, the pattern also has to be remapped.
        P_j = chars_remapped ? map_to_internal(P[j]) : P[j];

        if (!RS_L_.contains_character(P_j)) return;

        // Find the lexicographically smallest suffix in the current suffix array interval that is prefixed by P[j]
        if (P_j != L_(b_)) {
            /* To do so, we can at first find the first (sub-)run with character P[j] after the b_-th (sub-)run, save
            its index in b_ and set b to its start position M_LF.p(b_). */
            b_ = RS_L_.rank(P_j,b_);
            if (b_ == RS_L_.num_occurrences(P_j)) return;
            b_ = RS_L_.select(P_j,b_+1);
            b = M_LF.p(b_);
        }

        // Find the lexicographically largest suffix in the current suffix array interval that is prefixed by P[j]
        if (P_j != L_(e_)) {
            /* To do so, we can at first find the (sub-)last run with character P[j] before the e_-th (sub-)run, save
            its index in e_ and set e to its end position M_LF.p(e_+1)-1. */
            e_ = RS_L_.rank(P_j,e_);
            if (e_ == 0) return;
            e_ = RS_L_.select(P_j,e_);
            e = M_LF.p(e_+1)-1;
            
            // update y.
            y = j;
            // update \hat{e}'_y.
            hat_e_ap_y = e_;
        }
        
        // Else, because each suffix i in the previous suffix array interval starts with P[j+1..m] and the current 
        // interval [b,e] contains all suffixes of it, before which there is a P[j] in T, all suffixes in the 
        // interval SA[LF(b),LF(e)] start with P[j..m]

        /* If the suffix array interval [LF(b),LF(e)] of P[j..m] is empty, then b > e,
        because LF(j) is monotonic for a fixed L[j], hence it suffices to check, whether
        b <= e holds. */

        // If the suffix array interval is empty, P does not occur in T, so return.
        if (b > e) return;
        
        // Else, set b <- LF(b) and e <- LF(e)
        M_LF.move(b,b_);
        M_LF.move(e,e_);
    }
    
    /* Initially, this is the value in the suffix array at the right interval limit of the suffix array
    interval of P. This variable is then used to iterate over the values of the suffix array in the suffix
    array interval of P from right (e) to left (b). */
    uint_t s;
    // The index of the input interval in M_Phi containing s.
    uint_t s_;

    setup_phi_move_pair(hat_e_ap_y,s,s_);
    s -= y+1;

    /* Because s contains the value in the suffix array of the right limit of the suffix array
    interval of P, s is an occurrence position of P in T. */
    report(s);

    if (b < e) {
        /* If s_1 < M_Phi.p[\hat{s}'_y] now an input interval of M_Phi before the \hat{s}'_y-th one contains s,
        so we have to decrease s_. To find the correct value for s_, we perform an exponential search to the
        left over the input interval starting positions of M_Phi starting at \hat{s}'_y = s_. */
        if (s < M_Phi.p(s_)) {
            s_ = exp_search_max_leq<uint_t,LEFT>(s,s_-y-1,s_,[this](uint_t x){return M_Phi.p(x);});
        }
        
        // Perform e-b Phi queries with (s,s_) and at each step, report s.
        for (uint_t i=b; i<e; i++) {
            M_Phi.move(s,s_);
            report(s);
        }
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
        s_l = bin_search_min_gt<uint_t>(l,0,p_r-1,[this](uint_t x){return D_e[x].second;});
        s_r = bin_search_min_geq<uint_t>(r,0,p_r-1,[this](uint_t x){return D_e[x].second;});
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

        /* index of the input interval in M_LF containing i. */
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
        report(j,access_l_(x));

        // report T[l,r-1] from right to left
        while (j > j_l) {
            // Set i <- LF(i) and j <- j-1.
            M_LF.move(i,x);
            j--;
            // Report T[j] = L[i] = L'[x].
            report(j,access_l_(x));
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
        uint_t b = l+i_p*((r-l+1)/p);
        // Iteration range end position of thread i_p.
        uint_t e = i_p == p-1 ? r : l+(i_p+1)*((r-l+1)/p)-1;

        // Current position in the bwt.
        uint_t i = b;

        /* index of the input interval in M_LF containing i. */
        uint_t x = bin_search_max_leq<uint_t>(i,0,r_-1,[this](uint_t x_){return M_LF.p(x_);});

        // start position of the next input interval in M_LF
        uint_t l_xp1;

        // iterate until x is the input interval containing e
        while ((l_xp1 = M_LF.p(x+1)) <= e) {

            // iterate over all positions in the x-th input interval
            while (i < l_xp1) {
                report(i,access_l_(x));
                i++;
            }
            
            x++;
        }

        // report the remaining characters
        while (i <= e) {
            report(i,access_l_(x));
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
        uint_t b = l+i_p*((r-l+1)/p);
        // iteration range end position
        uint_t e = i_p == p-1 ? r : l+(i_p+1)*((r-l+1)/p)-1;

        // the input interval of M_LF containing i
        uint_t x = bin_search_max_leq<uint_t>(e,0,r_-1,[this](uint_t x_){return M_LF.p(x_);});
        // current position in the suffix array, initially the end position of the input interval of M_LF containing e
        uint_t i = M_LF.p(x+1)-1;
        // index of the input interval in M_Phi containing s
        uint_t s_;
        /* the current suffix array value (SA[i]), initially the suffix array sample of the x-th run,
           initially the suffix array value at e */
        uint_t s;

        setup_phi_move_pair(x,s,s_);

        // iterate down to the iteration range end position
        while (i > e) {
            M_Phi.move(s,s_);
            i--;
        }

        // report SA[e]
        report(i,s);

        // report the SA-values SA[e,e-1] from right to left
        while (i > b) {
            M_Phi.move(s,s_);
            i--;
            report(i,s);
        }
    }
}