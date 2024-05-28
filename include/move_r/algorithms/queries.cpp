template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
void move_r<locate_support,sym_t,pos_t>::setup_phi_move_pair(pos_t& x, pos_t& s, pos_t& s_) const {
    // the index of the pair in M_Phi creating the output interval with start position s = SA[M_LF.p[x+1]-1]
    pos_t x_s_ = SA_Phi(x);

    // set s_ to the index of the input interval in M_Phi containing s
    s_ = M_Phi().idx(x_s_);
    
    // compute s
    s = M_Phi().p(s_)+M_Phi().offs(x_s_);
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
sym_t move_r<locate_support,sym_t,pos_t>::BWT(pos_t i) const {
    // find the index of the input interval in M_LF containing i with a binary search.
    return unmap_symbol(L_(bin_search_max_leq<pos_t>(i,0,r_-1,[this](pos_t x){return M_LF().p(x);})));
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
pos_t move_r<locate_support,sym_t,pos_t>::SA(pos_t i) const {
    if constexpr (locate_support == _rlzdsa) {
        pos_t x_p,x_lp,x_cp,x_r,s_np;

        // variable to store SA[i] in
        pos_t s;

        // initialize the rlzdsa context to position i
        init_rlzdsa(i,s,x_p,x_lp,x_cp,x_r,s_np);
        
        // compute SA[i]
        next_rlzdsa(i,s,x_p,x_lp,x_cp,x_r,s_np);

        return s;
    } else {
        // index of the input interval in M_LF containing i.
        pos_t x = bin_search_max_leq<pos_t>(i,0,r_-1,[this](pos_t x_){return M_LF().p(x_);});

        /* if i is a run start position (i = M_LF.p(x)), then
            SA[i] = Phi^{-1}(SA_s'[(x-1) mod r'])
                  = Phi^{-1}(M_Phi.q(SA_Phi[(x-1) mod r']))
                  =          M_Phi.p(SA_Phi[(x-1) mod r'])
        */
        if (i == M_LF().p(x)) {
            if (x == 0) {
                return M_Phi().p(SA_Phi(r_-1));
            } else if (SA_Phi(x-1) != r__) {
                return M_Phi().p(SA_Phi(x-1));
            }
        }

        // increment x until the end position of the x-th input interval of M_LF is an end position of a bwt run
        while (SA_Phi(x) == r__) {
            x++;
        }

        // begin iterating at the end of the x-th run, because there is a
        // suffix array sample at the end position of the x-th input interval

        // position in the suffix array of the current suffix s
        pos_t j = M_LF().p(x+1)-1;

        // index of the input interval in M_Phi containing s
        pos_t s_;
        // the current suffix (s = SA[j])
        pos_t s;

        setup_phi_move_pair(x,s,s_);

        // Perform Phi-move queries until s is the suffix at position
        // i; in each iteration, s = SA[j] = \Phi^{i-j}(SA[i]) holds.
        while (j > i) {
            // Set s = \Phi(s)
            M_Phi().move(s,s_);
            j--;
        }

        // Since j = i, now s = SA[i] holds.
        return s;
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
pos_t move_r<locate_support,sym_t,pos_t>::query_context::next_occ() {
    if constexpr (locate_support == _rlzdsa) {
        if (i == b) {
            // compute the suffix array value at b
            s = idx->SA_s(hat_b_ap_z)-(z+1);
            i++;
            
            // check if there is more than one occurrence
            if (b < e) {
                idx->init_rlzdsa(i,x_p,x_lp,x_cp,x_r,s_np);
            }

            return s;
        } else {
            idx->next_rlzdsa(i,s,x_p,x_lp,x_cp,x_r,s_np);
            return s;
        }
    } else {
        if (i == b) {
            // compute the suffix array value at b
            idx->init_phi(b,e,s,s_,hat_e_ap_y,y);
            i++;
            return s;
        } else {
            idx->M_Phi().move(s,s_);
            i++;
            return s;
        }
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
void move_r<locate_support,sym_t,pos_t>::query_context::locate(std::vector<pos_t>& Occ) {
    if constexpr (locate_support == _rlzdsa) {
        if (i == b) {
            // compute the suffix array value at b
            s = idx->SA_s(hat_b_ap_z)-(z+1);
            Occ.emplace_back(s);
            i++;
            
            // check if there is more than one occurrence
            if (b < e) {
                idx->init_rlzdsa(i,x_p,x_lp,x_cp,x_r,s_np);
            }
        }

        // compute the remaining occurrences SA(b,e]
        if (i <= e) {
            locate_rlzdsa(i,e,s,x_p,x_lp,x_cp,x_r,s_np,Occ);
        }
    } else {
        // compute the suffix array value at b
        if (i == b) {
            idx->init_phi(i,b,e,s,s_,hat_e_ap_y,y);
            Occ.emplace_back(s);
            i++;
        }
        
        // compute the remaining occurrences SA(b,e]
        while (i <= e) {
            M_Phi().move(s,s_);
            Occ.emplace_back(s);
            i++;
        }
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
void move_r<locate_support,sym_t,pos_t>::backward_search_step(
    sym_t sym,
    pos_t& b, pos_t& e,
    pos_t& b_, pos_t& e_,
    int64_t& y, pos_t& hat_e_ap_y,
    int64_t& z, pos_t& hat_b_ap_z
) const {
    // If the characters have been remapped internally, the pattern also has to be remapped.
    sym = map_symbol(sym);

    // If sym does not occur in L', then P[i..m] does not occur in T
    if (sym == 0 || !RS_L_().contains(sym)) {e = 0; b = 1; return;}

    // Find the lexicographically smallest suffix in the current suffix array interval that is prefixed by P[i]
    if (sym != L_(b_)) {
        /* To do so, we can at first find the first (sub-)run with character P[i] after the b_-th (sub-)run, save
        its index in b_ and set b to its start position M_LF.p(b_). */
        b_ = RS_L_().rank(sym,b_);
        if (b_ == RS_L_().frequency(sym)) {e = 0; b = 1; return;}
        b_ = RS_L_().select(sym,b_+1);
        b = M_LF().p(b_);
        
        // update z (Case 1).
        z = 0;
        // update \hat{b}'_z.
        hat_b_ap_z = b_;
    } else {
        z++;
    }

    // Find the lexicographically largest suffix in the current suffix array interval that is prefixed by P[i]
    if (sym != L_(e_)) {
        /* To do so, we can at first find the (sub-)last run with character P[i] before the e_-th (sub-)run, save
        its index in e_ and set e to its end position M_LF.p(e_+1)-1. */
        e_ = RS_L_().rank(sym,e_);
        if (e_ == 0) {e = 0; b = 1; return;}
        e_ = RS_L_().select(sym,e_);
        e = M_LF().p(e_+1)-1;
        
        // update y (Case 1).
        y = 0;
        // update \hat{e}'_y.
        hat_e_ap_y = e_;
    } else {
        y++;
    }
    /* Else Case 2 holds, hence y(i) = y(i+1) and therefore, \hat{e}'_{y(i)} = \hat{e}'_{y(i+1)}, so
        don't change y and \hat{e}'_y */
    
    // Else, because each suffix i in the previous suffix array interval starts with P[i+1..m] and the current 
    // interval [b,e] contains all suffixes of it, before which there is a P[i] in T, all suffixes in the 
    // interval SA[LF(b),LF(e)] start with P[i..m]

    /* If the suffix array interval [LF(b),LF(e)] of P[i..m] is empty, then b > e,
    because LF(i) is monotonic for a fixed L[i], hence it suffices to check, whether
    b <= e holds. */

    // If the suffix array interval is empty, P does not occur in T, so return.
    if (b > e) {
        return;
    }
    
    /* Else, set b <- LF(b) and e <- LF(e). The following two optimizations increase query throughput slightly
        if there are only few occurrences */
    if (b_ == e_) {
        if (b == e) {
            /* If \hat{b'}_i == \hat{e'}_i and b'_i = e'_i, then computing
            (e_i,\hat{e}_i) <- M_LF.move(e'_i,\hat{e'}_i) is redundant */
            M_LF().move(b,b_);
            e = b;
            e_ = b_;
        } else {
            /* If \hat{b'}_i == \hat{e'}_i, but b'_i != e'_i, then e_i = b_i + e'_i - b'_i and therefore
            \hat{b'}_i < \hat{e'}_i, hence we can compute \hat{e'}_i by setting e_ <- \hat{b'}_i = b_ and
            incrementing e_ until e < M_LF.p[e_+1] holds; This takes O(a) time because of the a-balancedness property */
            pos_t diff_eb = e - b;
            M_LF().move(b,b_);
            e = b + diff_eb;
            e_ = b_;
            
            while (e >= M_LF().p(e_+1)) {
                e_++;
            }
        }
    } else {
        M_LF().move(b,b_);
        M_LF().move(e,e_);
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
void move_r<locate_support,sym_t,pos_t>::init_phi(
    pos_t& b, pos_t& e,
    pos_t& s, pos_t& s_,
    pos_t& hat_e_ap_y, int64_t& y
) const {
    setup_phi_move_pair(hat_e_ap_y,s,s_);
    s -= y+1;

    // If there is more than one occurrence and s < M_Phi.p[s_], now an input interval of M_Phi before 
    // the s_-th one contains s, so we have to decrease s_. To find the correct value for s_, we perform
    // an exponential search to the left over the input interval starting positions of M_Phi starting at s_.
    if (b < e && s < M_Phi().p(s_)) {
        s_ = exp_search_max_leq<pos_t,LEFT>(s,0,s_,[this](pos_t x){return M_Phi().p(x);});
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
void move_r<locate_support,sym_t,pos_t>::init_rlzdsa(
    pos_t& i,
    pos_t& x_p, pos_t& x_lp, pos_t& x_cp, pos_t& x_r, pos_t& s_np
) const {
    // copy-phrase index of the last copy-phrase
    pos_t x_cp_lcp = SCP().rank_1(i+1);
    // phrase-index of the last copy-phrase (possibly -1, if there is none before or at i)
    int64_t x_p_lcp = x_cp_lcp == 0 ? int64_t{-1} : PT().select_0(x_cp_lcp);
    // number of literal phrases between the current and the next copy-phrase
    pos_t n_lp = PT().select_0(x_cp_lcp+1)-x_p_lcp-1;
    // starting position of the next copy-phrase
    pos_t s_ncp = SCP().select_1(x_cp_lcp+1);

    if (i >= s_ncp-n_lp) {
        // there is a literal phrase at position i
        s_np = i+1;
        x_cp = x_cp_lcp;
        x_r = SR(x_cp);
        x_p = x_p_lcp+1+(i-(s_ncp-n_lp));
    } else {
        // i lies within a copy-phrase
        x_cp = x_cp_lcp-1;
        x_p = x_p_lcp;
        x_r = SR(x_cp)+(i-SCP().select_1(x_cp_lcp));
        s_np = s_ncp-n_lp;
    }

    x_lp = x_p-x_cp;
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
void move_r<locate_support,sym_t,pos_t>::init_rlzdsa(
    pos_t& i, pos_t& s,
    pos_t& x_p, pos_t& x_lp, pos_t& x_cp, pos_t& x_r, pos_t& s_np
) const {
    // copy-phrase index of the last copy-phrase
    pos_t x_cp_lcp = SCP().rank_1(i+1);
    // phrase-index of the last copy-phrase (possibly -1, if there is none before or at i)
    int64_t x_p_lcp = x_cp_lcp == 0 ? int64_t{-1} : PT().select_0(x_cp_lcp);
    // number of literal phrases between the current and the next copy-phrase
    pos_t n_lp = PT().select_0(x_cp_lcp+1)-x_p_lcp-1;
    // starting position of the next copy-phrase
    pos_t s_ncp = SCP().select_1(x_cp_lcp+1);

    if (i >= s_ncp-n_lp) {
        // there is a literal phrase at i
        s_np = i+1;
        x_cp = x_cp_lcp;
        x_r = SR(x_cp);
        x_p = x_p_lcp+1+(i-(s_ncp-n_lp));
        x_lp = x_p-x_cp;
    } else {
        // i lies within a copy-phrase, so find the last literal phrase before the copy-phrase contianing i

        // literal phrase-index of the last copy-phrase before i
        pos_t x_lp_llp = PT().rank_1(x_p_lcp)-1;

        // phrase-index of the last literal phrase before i
        pos_t x_p_llp = PT().select_1(x_lp_llp+1);

        // copy-phrase-index of the first copy-phrase after the last literal phrase before i
        pos_t x_cp_fcpallp = PT().rank_0(x_p_llp);

        // start iterating at the last literal phrase before i
        pos_t j = SCP().select_1(x_cp_fcpallp+1)-1;
        x_cp = x_cp_fcpallp;
        x_r = SR(x_cp);
        x_p = x_p_llp;
        s_np = j+1;
        x_lp = x_p-x_cp;

        // get the suffix array value at j
        s = LP(x_lp);

        // prepare to decode the first copy-phrase after the last literal phrase before i
        j++;
        x_p++;
        x_lp++;
        s_np = SCP().select_1(x_cp+2)-(PT().select_0(x_cp+2)-x_p)+1;

        // decode all copy-phrases after the last literal phrase before i and up to the copy-phrase containing i
        while (j < i) {
            s += R(x_r);
            s -= n;
            j++;
            x_r++;

            // there is a new copy-phrase starting at j
            if (s_np <= j) {
                x_p++;
                x_cp++;
                x_r = SR(x_cp);
                s_np = SCP().select_1(x_cp+2)-(PT().select_0(x_cp+2)-x_p)+1;
            }
        }
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
void move_r<locate_support,sym_t,pos_t>::next_rlzdsa(
    pos_t& i, pos_t& s,
    pos_t& x_p, pos_t& x_lp, pos_t& x_cp, pos_t& x_r, pos_t& s_np
) const {
    if (PT(x_p) == 1) {
        // literal phrase
        s = LP(x_lp);
        i++;
        x_p++;
        x_lp++;

        if (i < n) {
            if (PT(x_p) == 1) {
                // the next phrase is a literal phrase
                s_np++;
            } else {
                // the next phrase is a copy-phrase
                s_np = SCP().select_1(x_cp+2)-(PT().select_0(x_cp+2)-x_p)+1;
            }
        }
    } else {
        // copy-prhase
        s += R(x_r);
        s -= n;
        i++;
        x_r++;

        // there is a new phrase starting at i
        if (s_np <= i && i < n) {
            x_p++;
            x_cp++;
            x_r = SR(x_cp);

            if (PT(x_p) == 1) {
                // the next phrase is a literal phrase
                s_np++;
            } else {
                // the next phrase is a copy-phrase
                s_np = SCP().select_1(x_cp+2)-(PT().select_0(x_cp+2)-x_p)+1;
            }
        }
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
void move_r<locate_support,sym_t,pos_t>::locate_rlzdsa(
    pos_t& i, pos_t& e, pos_t& s,
    pos_t& x_p, pos_t& x_lp, pos_t& x_cp, pos_t& x_r, pos_t& s_np,
    std::vector<pos_t>& Occ
) const {
    while (i <= e) {
        // decode all copy phrases before the next literal phrase
        while (!PT(x_p) && i <= e) {
            // decode the x_cp-th copy phrase
            while (i < s_np && i <= e) {
                s += R(x_r);
                s -= n;
                Occ.emplace_back(s);
                i++;
                x_r++;
            }

            x_p++;
            x_cp++;
            x_r = SR(x_cp);
            s_np = SCP().select_1(x_cp+2)-(PT().select_0(x_cp+2)-x_p)+1;
        }

        // decode all literal phrases before the next copy phrase
        while (PT(x_p) && i <= e) {
            // decode the x_lp-th literal phrase
            s = LP(x_lp);
            Occ.emplace_back(s);
            i++;
            x_p++;
            x_lp++;
        }

        // set s_np to the starting position of the next (the x_lp-th)
        // literal phrase after the current (the x_cp-th) copy phrase
        s_np = SCP().select_1(x_cp+2)-(PT().select_0(x_cp+2)-x_p)+1;
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
pos_t move_r<locate_support,sym_t,pos_t>::count(const input_t& P) const {
    pos_t l,b,e,b_,e_,hat_e_ap_y,hat_b_ap_z;
    int64_t y,z;

    init_backward_search(b,e,b_,e_,y,hat_e_ap_y,z,hat_b_ap_z);

    for (int64_t i=P.size()-1; i>=0; i--) {
        backward_search_step(P[i],b,e,b_,e_,y,hat_e_ap_y,z,hat_b_ap_z);

        if (b > e) {
            return 0;
        }
    }

    return e-b+1;
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
void move_r<locate_support,sym_t,pos_t>::locate(const input_t& P, std::vector<pos_t>& Occ) const {
    pos_t l,b,e,b_,e_,hat_e_ap_y,hat_b_ap_z;
    int64_t y,z;

    init_backward_search(b,e,b_,e_,y,hat_e_ap_y,z,hat_b_ap_z);

    for (int64_t i=P.size()-1; i>=0; i--) {
        backward_search_step(P[i],b,e,b_,e_,y,hat_e_ap_y,z,hat_b_ap_z);

        if (b > e) {
            return;
        }
    }

    Occ.reserve(e-b+1);
    
    if constexpr (locate_support == _rlzdsa) {
        pos_t s = SA_s(hat_b_ap_z)-(z+1);
        Occ.emplace_back(s);

        if (b < e) {
            pos_t i = b+1;
            pos_t x_p,x_lp,x_cp,x_r,s_np;

            init_rlzdsa(i,x_p,x_lp,x_cp,x_r,s_np);
            locate_rlzdsa(i,e,s,x_p,x_lp,x_cp,x_r,s_np,Occ);
        }
    } else {
        pos_t s,s_;
        init_phi(b,e,s,s_,hat_e_ap_y,y);
        Occ.emplace_back(s);

        if (b < e) {
            pos_t i = b+1;
            
            while (i <= e) {
                M_Phi().move(s,s_);
                Occ.emplace_back(s);
                i++;
            }
        }
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
void move_r<locate_support,sym_t,pos_t>::revert(const std::function<void(pos_t,sym_t)>& report, retrieve_params params) const {
    adjust_retrieve_params(params,n-2);

    pos_t l = params.l;
    pos_t r = params.r;

    // leftmost section to revert
    uint16_t s_l;
    // rightmost section to revert
    uint16_t s_r;

    if (p_r == 1) {
        s_l = 0;
        s_r = 0;
    } else {
        s_l = bin_search_min_gt<pos_t>(l,0,p_r-1,[this](pos_t x){return _D_e[x].second;});
        s_r = bin_search_min_geq<pos_t>(r,0,p_r-1,[this](pos_t x){return _D_e[x].second;});
    }

    uint16_t p = std::max(
        (uint16_t)1,                         // use at least one thread
        std::min({
            (uint16_t)(s_r-s_l+1),           // use at most s_r-s_l+1 threads
            (uint16_t)omp_get_max_threads(), // use at most all threads
            params.num_threads               // use at most the specified number of threads
        })
    );

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // leftmost section for thread i_p to revert
        uint16_t sl_ip = s_l+(i_p*(s_r-s_l+1))/p;
        // rightmost section for thread i_p to revert
        uint16_t sr_ip = i_p == p-1 ? s_r : s_l+((i_p+1)*(s_r-s_l+1))/p-1;

        // Iteration range start position of thread i_p.
        pos_t j_l = std::max(l,sl_ip == 0 ? 0 : (_D_e[sl_ip-1].second+1)%n);
        // Iteration range end position of thread i_p.
        pos_t j_r = sr_ip == p_r-1 ? n-2 : _D_e[sr_ip].second;

        // index of the input interval in M_LF containing i.
        pos_t x = sr_ip == p_r-1 ? 0 : _D_e[sr_ip].first;
        // The position in the bwt of the current character in T.
        pos_t i = sr_ip == p_r-1 ? 0 : M_LF().p(x+1)-1;

        // start iterating at the right iteration range end position
        pos_t j = j_r;

        // iterate until j = r
        while (j > r) {
            // Set i <- LF(i) and j <- j-1.
            M_LF().move(i,x);
            j--;
        }

        // Report T[r] = T[j] = L[i] = L'[x]
        report(j,unmap_symbol(L_(x)));

        // report T[l,r-1] from right to left
        while (j > j_l) {
            // Set i <- LF(i) and j <- j-1.
            M_LF().move(i,x);
            j--;
            // Report T[j] = L[i] = L'[x].
            report(j,unmap_symbol(L_(x)));
        }
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
void move_r<locate_support,sym_t,pos_t>::BWT(const std::function<void(pos_t,sym_t)>& report, retrieve_params params) const {
    adjust_retrieve_params(params,n-1);

    pos_t l = params.l;
    pos_t r = params.r;

    uint16_t p = std::max(
        (uint16_t)1,                         // use at least one thread
        std::min({
            (uint16_t)omp_get_max_threads(), // use at most all threads
            (uint16_t)((r-l+1)/10),          // use at most (r-l+1)/100 threads
            params.num_threads               // use at most the specified number of threads
        })
    );

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Iteration range start position of thread i_p.
        pos_t b = l+i_p*((r-l+1)/p);
        // Iteration range end position of thread i_p.
        pos_t e = i_p == p-1 ? r : l+(i_p+1)*((r-l+1)/p)-1;

        // Current position in the bwt.
        pos_t i = b;

        // index of the input interval in M_LF containing i.
        pos_t x = bin_search_max_leq<pos_t>(i,0,r_-1,[this](pos_t x_){return M_LF().p(x_);});

        // start position of the next input interval in M_LF
        pos_t l_xp1;

        // iterate until x is the input interval containing e
        while ((l_xp1 = M_LF().p(x+1)) <= e) {

            // iterate over all positions in the x-th input interval
            while (i < l_xp1) {
                report(i,unmap_symbol(L_(x)));
                i++;
            }
            
            x++;
        }

        // report the remaining characters
        while (i <= e) {
            report(i,unmap_symbol(L_(x)));
            i++;
        }
    }
}

template <move_r_locate_supp locate_support, typename sym_t, typename pos_t>
void move_r<locate_support,sym_t,pos_t>::SA(const std::function<void(pos_t,pos_t)>& report, retrieve_params params) const {
    adjust_retrieve_params(params,n-1);

    pos_t l = params.l;
    pos_t r = params.r;

    uint16_t p = std::max(
        (uint16_t)1,                                           // use at least one thread
        std::min({
            (uint16_t)omp_get_max_threads(),                   // use at most all threads
            params.num_threads,                                // use at most the specified number of threads
            (uint16_t)(((r-l+1)*(double)r__)/(10.0*(double)n)) // use at most (r-l+1)*(r/n)*(1/10) threads
        })
    );

    if constexpr (locate_support == _rlzdsa) {
        #pragma omp parallel num_threads(p)
        {
            // Index in [0..p-1] of the current thread.
            uint16_t i_p = omp_get_thread_num();

            // iteration range start position
            pos_t b = l+i_p*((r-l+1)/p);
            // iteration range end position
            pos_t e = i_p == p-1 ? r : l+(i_p+1)*((r-l+1)/p)-1;

            pos_t x_p,x_lp,x_cp,x_r,s_np;

            // current position in the suffix array
            pos_t i = b;
            // current suffix array value
            pos_t s;

            // initialize the rlzdsa context to position i = b
            init_rlzdsa(i,s,x_p,x_lp,x_cp,x_r,s_np);

            // decode and report SA[b..e]
            while (i <= e) {
                next_rlzdsa(i,s,x_p,x_lp,x_cp,x_r,s_np);
                report(i-1,s);
            }
        }
    } else {
        #pragma omp parallel num_threads(p)
        {
            // Index in [0..p-1] of the current thread.
            uint16_t i_p = omp_get_thread_num();

            // iteration range start position
            pos_t b = l+i_p*((r-l+1)/p);
            // iteration range end position
            pos_t e = i_p == p-1 ? r : l+(i_p+1)*((r-l+1)/p)-1;

            // the input interval of M_LF containing i
            pos_t x = bin_search_max_leq<pos_t>(e,0,r_-1,[this](pos_t x_){return M_LF().p(x_);});

            // increment x until the end position of the x-th input interval of M_LF is an end position of a bwt run
            while (SA_Phi(x) == r__) {
                x++;
            }

            // current position in the suffix array, initially the end position of the x-th interval of M_LF
            pos_t i = M_LF().p(x+1)-1;

            // index of the input interval in M_Phi containing s
            pos_t s_;
            /* the current suffix array value (SA[i]), initially the suffix array sample of the x-th run,
            initially the suffix array value at e */
            pos_t s;

            setup_phi_move_pair(x,s,s_);

            // iterate down to the iteration range end position
            while (i > e) {
                M_Phi().move(s,s_);
                i--;
            }

            // report SA[e]
            report(i,s);

            // report the SA-values SA[b,e-1] from right to left
            while (i > b) {
                M_Phi().move(s,s_);
                i--;
                report(i,s);
            }
        }
    }
}