#pragma once

#include <move_r/move_r.hpp>

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::setup_phi_m1_move_pair(pos_t& x, pos_t& s, pos_t& s_) const
    requires(support == _locate_move)
{
    // the index of the pair in M_Phi^{-1} creating the output interval with starting position s = SA[M_LF.p[x]]
    pos_t x_s_ = SA_Phi_m1(x);

    // set s_ to the index of the input interval in M_Phi^{-1} containing s
    s_ = M_Phi_m1().idx(x_s_);

    // compute s
    s = M_Phi_m1().p(s_) + M_Phi_m1().offs(x_s_);
}

template <move_r_support support, typename sym_t, typename pos_t>
sym_t move_r<support, sym_t, pos_t>::BWT(pos_t i) const
{
    // find the index of the input interval in M_LF containing i with a binary search.
    return unmap_symbol(L_(bin_search_max_leq<pos_t>(i, 0, r_ - 1, [this](pos_t x) { return M_LF().p(x); })));
}

template <move_r_support support, typename sym_t, typename pos_t>
pos_t move_r<support, sym_t, pos_t>::SA(pos_t i) const
    requires(supports_multiple_locate)
{
    if constexpr (support == _locate_rlzdsa) {
        pos_t x_p, x_lp, x_cp, x_r, s_np;

        // variable to store SA[i] in
        pos_t s;

        // initialize the rlzdsa context to position i
        init_rlzdsa(i, s, x_p, x_lp, x_cp, x_r, s_np);

        // compute SA[i]
        next_rlzdsa(i, s, x_p, x_lp, x_cp, x_r, s_np);

        return s;
    } else {
        // index of the input interval in M_LF containing i.
        pos_t x = bin_search_max_leq<pos_t>(i, 0, r_ - 1, [this](pos_t x_) { return M_LF().p(x_); });

        /* if i is a bwt run end position (i = M_LF.p(x+1)-1) and SA_Phi^{-1}[x+1] != r'', then
            SA[i] = Phi(SA_s[(x+1) mod r'])
                  = Phi(M_Phi^{-1}.q(SA_Phi^{-1}[(x+1) mod r']))
                  =     M_Phi^{-1}.p(SA_Phi^{-1}[(x+1) mod r'])
        */
        if (i == M_LF().p(x + 1) - 1) [[unlikely]] {
            pos_t xp1 = (x + 1) == r_ ? 0 : (x + 1);

            if (SA_Phi_m1(xp1) != r__) [[unlikely]] {
                return M_Phi_m1().p(SA_Phi_m1(xp1));
            }
        }

        // decrement x until the starting position of the x-th input interval of M_LF is a starting position of a bwt run
        while (SA_Phi_m1(x) == r__) {
            x--;
        }

        // begin iterating at the start of the x-th run, because there is a
        // suffix array sample at the end position of the x-th input interval

        // position in the suffix array of the current suffix s
        pos_t j = M_LF().p(x);

        // index of the input interval in M_Phi^{-1} containing s
        pos_t s_;
        // the current suffix (s = SA[j])
        pos_t s;

        setup_phi_m1_move_pair(x, s, s_);

        // Perform Phi-move queries until s is the suffix at position
        // i; in each iteration, s = SA[j] = \Phi^{i-j}(SA[i]) holds.
        while (j < i) {
            // Set s = \Phi(s)
            M_Phi_m1().move(s, s_);
            j++;
        }

        // Since j = i, now s = SA[i] holds.
        return s;
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
bool move_r<support, sym_t, pos_t>::query_context::prepend(sym_t sym)
{
    pos_t b_tmp = b;
    pos_t e_tmp = e;
    pos_t b__tmp = b_;
    pos_t e__tmp = e_;
    pos_t hat_b_ap_y_tmp = hat_b_ap_y;
    int64_t y_tmp = y;
    pos_t hat_e_ap_z_tmp = hat_e_ap_z;
    int64_t z_tmp = z;

    if (idx->backward_search_step(sym, b, e, b_, e_, hat_b_ap_y, y, hat_e_ap_z, z)) {
        l++;
        i = b;

        return true;
    } else {
        b = b_tmp;
        e = e_tmp;
        b_ = b__tmp;
        e_ = e__tmp;
        hat_b_ap_y = hat_b_ap_y_tmp;
        y = y_tmp;
        hat_e_ap_z = hat_e_ap_z_tmp;
        z = z_tmp;

        return false;
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
pos_t move_r<support, sym_t, pos_t>::query_context::next_occ()
    requires(supports_multiple_locate)
{
    if constexpr (support == _locate_rlzdsa) {
        if (i == b) [[unlikely]] {
            // compute the suffix array value at b
            s = idx->SA_s(hat_b_ap_y) - (y + 1);
            i++;

            // check if there is more than one occurrence
            if (b < e) {
                idx->init_rlzdsa(i, x_p, x_lp, x_cp, x_r, s_np);
            }

            return s;
        } else {
            idx->next_rlzdsa(i, s, x_p, x_lp, x_cp, x_r, s_np);
            return s;
        }
    } else {
        if (i == b) [[unlikely]] {
            // compute the suffix array value at b
            idx->init_phi_m1(b, e, s, s_, hat_b_ap_y, y);
            i++;
            return s;
        } else {
            idx->M_Phi_m1().move(s, s_);
            i++;
            return s;
        }
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
pos_t move_r<support, sym_t, pos_t>::query_context::one_occ() const
    requires(supports_locate)
{
    return idx->SA_s(hat_b_ap_y) - (y + 1);
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::query_context::locate(std::vector<pos_t>& Occ)
    requires(supports_multiple_locate)
{
    Occ.reserve(Occ.size() + num_occ_rem());

    if constexpr (support == _locate_rlzdsa) {
        if (i == b) [[unlikely]] {
            // compute the suffix array value at b
            s = idx->SA_s(hat_b_ap_y) - (y + 1);
            Occ.emplace_back(s);
            i++;

            // check if there is more than one occurrence
            if (b < e) [[likely]] {
                idx->init_rlzdsa(i, x_p, x_lp, x_cp, x_r, s_np);
            }
        }

        // compute the remaining occurrences SA(b,e]
        if (i <= e) [[likely]] {
            pos_t o = Occ.size();
            no_init_resize(Occ, o + num_occ_rem());
            idx->write_rlzdsa_right(i, e, s, x_p, x_lp, x_cp, x_r, s_np, Occ, o);
        }
    } else {
        // compute the suffix array value at b
        if (i == b) [[unlikely]] {
            idx->init_phi_m1(b, e, s, s_, hat_b_ap_y, y);
            Occ.emplace_back(s);
            i++;
        }

        // compute the remaining occurrences SA(b,e]
        while (i <= e) {
            idx->M_Phi_m1().move(s, s_);
            Occ.emplace_back(s);
            i++;
        }
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
bool move_r<support, sym_t, pos_t>::backward_search_step(
    sym_t sym,
    pos_t& b, pos_t& e,
    pos_t& b_, pos_t& e_,
    pos_t& hat_b_ap_y, int64_t& y,
    pos_t& hat_e_ap_z, int64_t& z) const
{
    // If the characters have been remapped internally, the pattern also has to be remapped.
    i_sym_t i_sym = map_symbol(sym);

    // If sym does not occur in L', then P[i..m] does not occur in T
    if constexpr (byte_alphabet) {
        if (!RS_L_().contains(i_sym)) [[unlikely]] {
            return false;
        }
    } else {
        if (i_sym == 0) [[unlikely]] {
            return false;
        }
    }

    // Find the lexicographically smallest suffix in the current suffix array interval that is prefixed by P[i]
    if (i_sym != L_(b_)) {
        /* To do so, we can at first find the first (sub-)run with character P[i] after the b_-th (sub-)run, save
        its index in b_ and set b to its start position M_LF.p(b_). */

        if constexpr (byte_alphabet) {
            pos_t max_b_ = std::min<pos_t>(b_ + max_scan_l_, e_);

            while (b_ < max_b_ && i_sym != L_(b_)) {
                b_++;
            }
        }

        if (int_alphabet || (i_sym != L_(b_) && b_ < e_)) {
            b_ = RS_L_().rank(i_sym, b_);

            if (b_ == RS_L_().frequency(i_sym)) [[unlikely]] {
                return false;
            }

            b_ = RS_L_().select(i_sym, b_ + 1);

            if (b_ > e_) [[unlikely]] {
                return false;
            }
        }

        b = M_LF().p(b_);

        // update y (Case 1).
        y = 0;
        // update \hat{b}'_y.
        hat_b_ap_y = b_;
    } else {
        y++;
    }

    // Find the lexicographically largest suffix in the current suffix array interval that is prefixed by P[i]
    if (i_sym != L_(e_)) {
        /* To do so, we can at first find the (sub-)last run with character P[i] before the e_-th (sub-)run, save
        its index in e_ and set e to its end position M_LF.p(e_+1)-1. */

        if constexpr (byte_alphabet) {
            pos_t min_e_ = std::max<pos_t>(e_ <= max_scan_l_ ? 0 : e_ - max_scan_l_, b_);

            while (e_ > min_e_ && i_sym != L_(e_)) {
                e_--;
            }
        }

        if (int_alphabet || (i_sym != L_(e_) && e_ > b_)) {
            e_ = RS_L_().select(i_sym, RS_L_().rank(i_sym, e_));
        }

        e = M_LF().p(e_ + 1) - 1;

        // update z (Case 1).
        z = 0;
        // update \hat{e}'_z.
        hat_e_ap_z = e_;
    } else {
        z++;
    }

    // Else, because each suffix i in the previous suffix array interval starts with P[i+1..m] and the current
    // interval [b,e] contains all suffixes of it, before which there is a P[i] in T, all suffixes in the
    // interval SA[LF(b),LF(e)] start with P[i..m]

    /* If the suffix array interval [LF(b),LF(e)] of P[i..m] is empty, then b > e,
    because LF(i) is monotonic for a fixed L[i], hence it suffices to check, whether
    b <= e holds. */

    // If the suffix array interval is empty, P does not occur in T, so return false.
    if (b > e) [[unlikely]] {
        return false;
    }

    /* Else, set b <- LF(b) and e <- LF(e). The following two optimizations increase query throughput slightly
        if there are only few occurrences */
    if (b_ == e_) {
        if (b == e) {
            /* If \hat{b'}_i == \hat{e'}_i and b'_i = e'_i, then computing
            (e_i,\hat{e}_i) <- M_LF.move(e'_i,\hat{e'}_i) is redundant */
            M_LF().move(b, b_);
            e = b;
            e_ = b_;
        } else {
            /* If \hat{b'}_i == \hat{e'}_i, but b'_i != e'_i, then e_i = b_i + e'_i - b'_i and therefore
            \hat{b'}_i < \hat{e'}_i, hence we can compute \hat{e'}_i by setting e_ <- \hat{b'}_i = b_ and
            incrementing e_ until e < M_LF.p[e_+1] holds; This takes O(a) time because of the a-balancedness property */
            pos_t diff_eb = e - b;
            M_LF().move(b, b_);
            e = b + diff_eb;
            e_ = b_;

            while (e >= M_LF().p(e_ + 1)) {
                e_++;
            }
        }
    } else {
        M_LF().move(b, b_);
        M_LF().move(e, e_);
    }

    return true;
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::init_phi_m1(
    pos_t& b, pos_t& e,
    pos_t& s, pos_t& s_,
    pos_t& hat_b_ap_y, int64_t& y) const
    requires(support == _locate_move)
{
    setup_phi_m1_move_pair(hat_b_ap_y, s, s_);
    s -= y + 1;

    // If there is more than one occurrence and s < M_Phi^{-1}.p[s_], now an input interval of M_Phi^{-1} before
    // the s_-th one contains s, so we have to decrease s_. To find the correct value for s_, we perform
    // an exponential search to the left over the input interval starting positions of M_Phi^{-1} starting at s_.
    if (b < e && s < M_Phi_m1().p(s_)) {
        s_ = exp_search_max_leq<pos_t, LEFT>(s, 0, s_, [this](pos_t x) { return M_Phi_m1().p(x); });
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::init_rlzdsa(
    pos_t& i,
    pos_t& x_p, pos_t& x_lp, pos_t& x_cp, pos_t& x_r, pos_t& s_np) const
    requires(support == _locate_rlzdsa)
{
    // index in SCP_S of the last sampled copy phrase starting before or at i
    pos_t x_scps = SCP_S().rank_1(i + 1);

    if (x_scps == 0) [[unlikely]] {
        // i lies before the first copy phrase
        s_np = i + 1;
        x_p = i;
        x_lp = i;
        x_cp = 0;
        x_r = SR(0);
    } else {
        // copy-phrase index of the last copy-phrase starting before or at i
        pos_t x_cp_lcp = x_scps == 1 ? 0 : (x_scps - 1) * sr_scp;
        // starting position of the last copy-phrase starting before or at i
        pos_t s_lcp = SCP_S().select_1(x_scps);
        // phrase index of the last copy-phrase starting before or at i
        pos_t x_p_lcp = PT().select_0(x_cp_lcp + 1);
        // phrase index of the x_cp_lcp+1-th copy phrase
        pos_t x_p_ncp = PT().select_0(x_cp_lcp + 2);
        // number of literal phrases between the current and the next copy-phrase
        pos_t n_lp = x_p_ncp - x_p_lcp - 1;
        // starting position of the next copy-phrase
        pos_t s_ncp = s_lcp + CPL(x_cp_lcp) + n_lp;

        // find the last copy-phrase starting before or at i
        while (s_ncp <= i) {
            x_cp_lcp++;
            s_lcp = s_ncp;
            x_p_lcp = x_p_ncp;
            x_p_ncp = PT().select_0(x_cp_lcp + 2);
            n_lp = x_p_ncp - x_p_lcp - 1;
            s_ncp += CPL(x_cp_lcp) + n_lp;
        }

        if (i >= s_ncp - n_lp) {
            // there is a literal phrase at position i
            s_np = i + 1;
            x_cp = x_cp_lcp + 1;
            x_r = SR(x_cp);
            x_p = x_p_lcp + 1 + (i - (s_ncp - n_lp));
        } else {
            // i lies within a copy-phrase
            x_cp = x_cp_lcp;
            x_p = x_p_lcp;
            x_r = SR(x_cp) + (i - s_lcp);
            s_np = s_ncp - n_lp;
        }

        x_lp = x_p - x_cp;
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::init_rlzdsa(
    pos_t& i, pos_t& s,
    pos_t& x_p, pos_t& x_lp, pos_t& x_cp, pos_t& x_r, pos_t& s_np) const
    requires(support == _locate_rlzdsa)
{
    // index of the input interval in M_LF containing i.
    pos_t x = bin_search_max_leq<pos_t>(i, 0, r_ - 1, [this](pos_t x_) { return M_LF().p(x_); });

    while (SA_s(x) == n) {
        x--;
    }

    s = SA_s(x);
    pos_t j = M_LF().p(x);

    if (j == i) [[unlikely]] {
        init_rlzdsa(i, x_p, x_lp, x_cp, x_r, s_np);

        if (!PT(x_p)) {
            pos_t sad_i = R(x_r);

            // set s <- SA[i-1] = SA[i]-SA^d[i]
            if (sad_i < n) {
                s += n - sad_i;
            } else {
                s -= sad_i - n;
            }
        }
    } else {
        j++;
        init_rlzdsa(j, x_p, x_lp, x_cp, x_r, s_np);
        skip_rlzdsa_right(j, i, s, x_p, x_lp, x_cp, x_r, s_np);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::skip_rlzdsa_right(
    pos_t& i, pos_t& e, pos_t& s,
    pos_t& x_p, pos_t& x_lp, pos_t& x_cp, pos_t& x_r, pos_t& s_np) const
    requires(support == _locate_rlzdsa)
{
    while (i < e) {
        // decode all copy-phrases before the next literal phrase
        while (i < e && !PT(x_p)) {
            // decode the x_cp-th copy-phrase
            while (i < s_np && i < e) {
                s += R(x_r);
                s -= n;
                i++;
                x_r++;
            }

            if (i < e || i == s_np) [[likely]] {
                x_p++;
                x_cp++;
                x_r = SR(x_cp);
                s_np += PT(x_p) ? 1 : CPL(x_cp);
            }
        }

        // decode all literal phrases before the next copy-phrase
        while (i < e && PT(x_p)) {
            // decode the x_lp-th literal phrase
            s = LP(x_lp);
            i++;
            x_p++;
            x_lp++;
            s_np++;
        }

        if (i < e || PT(x_p)) [[likely]] {
            // set s_np to the starting position of the next (the x_lp-th)
            // literal phrase after the current (the x_cp-th) copy-phrase
            s_np += CPL(x_cp) - 1;
        }
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::next_rlzdsa(
    pos_t& i, pos_t& s,
    pos_t& x_p, pos_t& x_lp, pos_t& x_cp, pos_t& x_r, pos_t& s_np) const
    requires(support == _locate_rlzdsa)
{
    if (PT(x_p)) {
        // literal phrase
        s = LP(x_lp);
        i++;
        x_p++;
        x_lp++;

        if (i < n) [[likely]] {
            if (PT(x_p)) {
                // the next phrase is a literal phrase
                s_np++;
            } else {
                // the next phrase is a copy-phrase
                s_np += CPL(x_cp);
            }
        }
    } else {
        // copy-prhase
        s += R(x_r);
        s -= n;
        i++;
        x_r++;

        // there is a new phrase starting at i
        if (s_np <= i && i < n) [[unlikely]] {
            x_p++;
            x_cp++;
            x_r = SR(x_cp);

            if (PT(x_p)) {
                // the next phrase is a literal phrase
                s_np++;
            } else {
                // the next phrase is a copy-phrase
                s_np += CPL(x_cp);
            }
        }
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::write_rlzdsa_right(
    pos_t& i, pos_t& e, pos_t& s,
    pos_t& x_p, pos_t& x_lp, pos_t& x_cp, pos_t& x_r, pos_t& s_np,
    std::vector<pos_t>& vec, pos_t o) const
    requires(support == _locate_rlzdsa)
{
    while (true) {
        // decode all copy-phrases before the next literal phrase
        while (!PT(x_p)) {
            // decode the x_cp-th copy-phrase
            while (i < s_np) {
                s += R(x_r);
                s -= n;
                vec[o] = s;

                if (i == e) [[unlikely]] {
                    return;
                }
                
                o++;
                i++;
                x_r++;
            }

            x_p++;
            x_cp++;
            x_r = SR(x_cp);
            s_np += PT(x_p) ? 1 : CPL(x_cp);
        }

        // decode all literal phrases before the next copy-phrase
        while (PT(x_p)) {
            // decode the x_lp-th literal phrase
            s = LP(x_lp);
            vec[o] = s;

            if (i == e) [[unlikely]] {
                return;
            }

            o++;
            i++;
            x_p++;
            x_lp++;
            s_np++;
        }

        // set s_np to the starting position of the next (the x_lp-th)
        // literal phrase after the current (the x_cp-th) copy-phrase
        s_np += CPL(x_cp) - 1;
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::report_rlzdsa_right(
    pos_t& i, pos_t& e, pos_t& s,
    pos_t& x_p, pos_t& x_lp, pos_t& x_cp, pos_t& x_r, pos_t& s_np,
    const std::function<void(pos_t, pos_t)>& report) const
    requires(support == _locate_rlzdsa)
{
    while (true) {
        // decode all copy-phrases before the next literal phrase
        while (!PT(x_p)) {
            // decode the x_cp-th copy-phrase
            while (i < s_np) {
                s += R(x_r);
                s -= n;
                report(i, s);

                if (i == e) [[unlikely]] {
                    return;
                }

                i++;
                x_r++;
            }

            x_p++;
            x_cp++;
            x_r = SR(x_cp);
            s_np += PT(x_p) ? 1 : CPL(x_cp);
        }

        // decode all literal phrases before the next copy-phrase
        while (PT(x_p)) {
            // decode the x_lp-th literal phrase
            s = LP(x_lp);
            report(i, s);

            if (i == e) [[unlikely]] {
                return;
            }

            i++;
            x_p++;
            x_lp++;
            s_np++;
        }

        // set s_np to the starting position of the next (the x_lp-th)
        // literal phrase after the current (the x_cp-th) copy-phrase
        s_np += CPL(x_cp) - 1;
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
pos_t move_r<support, sym_t, pos_t>::count(const inp_t& P) const
{
    pos_t b, e, b_, e_, hat_b_ap_y, hat_e_ap_z;
    int64_t y, z;

    init_backward_search(b, e, b_, e_, hat_b_ap_y, y, hat_e_ap_z, z);

    for (int64_t i = P.size() - 1; i >= 0; i--) {
        if (!backward_search_step(P[i], b, e, b_, e_, hat_b_ap_y, y, hat_e_ap_z, z)) {
            return 0;
        }
    }

    return e - b + 1;
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::locate(const inp_t& P, std::vector<pos_t>& Occ) const
    requires(supports_multiple_locate)
{
    pos_t b, e, b_, e_, hat_b_ap_y, hat_e_ap_z;
    int64_t y, z;

    init_backward_search(b, e, b_, e_, hat_b_ap_y, y, hat_e_ap_z, z);

    for (int64_t i = P.size() - 1; i >= 0; i--) {
        if (!backward_search_step(P[i], b, e, b_, e_, hat_b_ap_y, y, hat_e_ap_z, z)) {
            return;
        }
    }

    if constexpr (support == _locate_rlzdsa) {
        pos_t o = Occ.size();
        no_init_resize(Occ, o + e - b + 1);
        pos_t s = SA_s(hat_b_ap_y) - (y + 1);
        Occ[o] = s;

        if (b < e) {
            pos_t i = b + 1;
            pos_t x_p, x_lp, x_cp, x_r, s_np;

            init_rlzdsa(i, x_p, x_lp, x_cp, x_r, s_np);
            write_rlzdsa_right(i, e, s, x_p, x_lp, x_cp, x_r, s_np, Occ, o + 1);
        }
    } else {
        pos_t s, s_;
        init_phi_m1(b, e, s, s_, hat_b_ap_y, y);
        Occ.emplace_back(s);

        if (b < e) {
            pos_t i = b + 1;

            while (i <= e) {
                M_Phi_m1().move(s, s_);
                Occ.emplace_back(s);
                i++;
            }
        }
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::revert(const std::function<void(pos_t, sym_t)>& report, retrieve_params params) const
{
    adjust_retrieve_params(params, n - 2);

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
        s_l = bin_search_min_gt<pos_t>(l, 0, p_r - 1, [this](pos_t x) { return _D_e[x].second; });
        s_r = bin_search_min_geq<pos_t>(r, 0, p_r - 1, [this](pos_t x) { return _D_e[x].second; });
    }

    uint16_t p = std::max(
        (uint16_t)1, // use at least one thread
        std::min({
            (uint16_t)(s_r - s_l + 1), // use at most s_r-s_l+1 threads
            (uint16_t)omp_get_max_threads(), // use at most all threads
            params.num_threads // use at most the specified number of threads
        }));

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // leftmost section for thread i_p to revert
        uint16_t sl_ip = s_l + (i_p * (s_r - s_l + 1)) / p;
        // rightmost section for thread i_p to revert
        uint16_t sr_ip = i_p == p - 1 ? s_r : s_l + ((i_p + 1) * (s_r - s_l + 1)) / p - 1;

        // Iteration range start position of thread i_p.
        pos_t j_l = std::max(l, sl_ip == 0 ? 0 : (_D_e[sl_ip - 1].second + 1) % n);
        // Iteration range end position of thread i_p.
        pos_t j_r = sr_ip == p_r - 1 ? n - 2 : _D_e[sr_ip].second;

        // index of the input interval in M_LF containing i.
        pos_t x = sr_ip == p_r - 1 ? 0 : _D_e[sr_ip].first;
        // The position in the bwt of the current character in T.
        pos_t i = sr_ip == p_r - 1 ? 0 : M_LF().p(x);

        // start iterating at the right iteration range end position
        pos_t j = j_r;

        // iterate until j = r
        while (j > r) {
            // Set i <- LF(i) and j <- j-1.
            M_LF().move(i, x);
            j--;
        }

        // Report T[r] = T[j] = L[i] = L'[x]
        report(j, unmap_symbol(L_(x)));

        // report T[l,r-1] from right to left
        while (j > j_l) {
            // Set i <- LF(i) and j <- j-1.
            M_LF().move(i, x);
            j--;
            // Report T[j] = L[i] = L'[x].
            report(j, unmap_symbol(L_(x)));
        }
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::BWT(const std::function<void(pos_t, sym_t)>& report, retrieve_params params) const
{
    adjust_retrieve_params(params, n - 1);

    pos_t l = params.l;
    pos_t r = params.r;

    uint16_t p = std::max(
        (uint16_t)1, // use at least one thread
        std::min({
            (uint16_t)omp_get_max_threads(), // use at most all threads
            (uint16_t)((r - l + 1) / 10), // use at most (r-l+1)/100 threads
            params.num_threads // use at most the specified number of threads
        }));

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Iteration range start position of thread i_p.
        pos_t b = l + i_p * ((r - l + 1) / p);
        // Iteration range end position of thread i_p.
        pos_t e = i_p == p - 1 ? r : l + (i_p + 1) * ((r - l + 1) / p) - 1;

        // Current position in the bwt.
        pos_t i = b;

        // index of the input interval in M_LF containing i.
        pos_t x = bin_search_max_leq<pos_t>(i, 0, r_ - 1, [this](pos_t x_) { return M_LF().p(x_); });

        // start position of the next input interval in M_LF
        pos_t l_xp1;

        // iterate until x is the input interval containing e
        while ((l_xp1 = M_LF().p(x + 1)) <= e) {

            // iterate over all positions in the x-th input interval
            while (i < l_xp1) {
                report(i, unmap_symbol(L_(x)));
                i++;
            }

            x++;
        }

        // report the remaining characters
        while (i <= e) {
            report(i, unmap_symbol(L_(x)));
            i++;
        }
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::SA(const std::function<void(pos_t, pos_t)>& report, retrieve_params params) const
    requires(support == _locate_move)
{
    adjust_retrieve_params(params, n - 1);

    pos_t l = params.l;
    pos_t r = params.r;

    uint16_t p = std::max(
        (uint16_t)1, // use at least one thread
        std::min({
            (uint16_t)omp_get_max_threads(), // use at most all threads
            params.num_threads, // use at most the specified number of threads
            (uint16_t)(((r - l + 1) * (double)r__) / (10.0 * (double)n)) // use at most (r-l+1)*(r/n)*(1/10) threads
        }));

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // iteration range start position
        pos_t b = l + i_p * ((r - l + 1) / p);
        // iteration range end position
        pos_t e = i_p == p - 1 ? r : l + (i_p + 1) * ((r - l + 1) / p) - 1;

        // the input interval of M_LF containing i
        pos_t x = bin_search_max_leq<pos_t>(b, 0, r_ - 1, [this](pos_t x_) { return M_LF().p(x_); });

        // decrement x until the starting position of the x-th input interval of M_LF is a starting position of a bwt run
        while (SA_Phi_m1(x) == r__) {
            x--;
        }

        // current position in the suffix array, initially the starting position of the x-th interval of M_LF
        pos_t i = M_LF().p(x);

        // index of the input interval in M_Phi^{-1} containing s
        pos_t s_;
        /* the current suffix array value (SA[i]), initially the suffix array sample of the x-th run,
        initially the suffix array value at b */
        pos_t s;

        setup_phi_m1_move_pair(x, s, s_);

        // iterate up to the iteration range starting position
        while (i < b) {
            M_Phi_m1().move(s, s_);
            i++;
        }

        // report SA[b]
        report(i, s);

        // report the SA-values SA[b+1,e] from left to right
        while (i < e) {
            M_Phi_m1().move(s, s_);
            i++;
            report(i, s);
        }
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::SA(std::vector<pos_t>& SA_range, retrieve_params params) const
    requires(support == _locate_move)
{
    SA_range.reserve(SA_range.size() + (params.r - params.l) + 1);
    SA([&](pos_t i, pos_t s) { SA_range[i] = s; }, params);
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::SA(const std::function<void(pos_t, pos_t)>& report, retrieve_params params) const
    requires(support == _locate_rlzdsa)
{
    adjust_retrieve_params(params, n - 1);

    pos_t l = params.l;
    pos_t r = params.r;

    uint16_t p = std::max(
        (uint16_t)1, // use at least one thread
        std::min({
            (uint16_t)omp_get_max_threads(), // use at most all threads
            params.num_threads, // use at most the specified number of threads
            (uint16_t)(((r - l + 1) * (double)r__) / (10.0 * (double)n)) // use at most (r-l+1)*(r/n)*(1/10) threads
        }));

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // iteration range start position
        pos_t b = l + i_p * ((r - l + 1) / p);
        // iteration range end position
        pos_t e = i_p == p - 1 ? r : l + (i_p + 1) * ((r - l + 1) / p) - 1;

        pos_t x_p, x_lp, x_cp, x_r, s_np;

        // current suffix array value
        pos_t s;

        // initialize the rlzdsa context to position b
        init_rlzdsa(b, s, x_p, x_lp, x_cp, x_r, s_np);

        // decode and report SA(b..e]
        report_rlzdsa_right(b, e, s, x_p, x_lp, x_cp, x_r, s_np, report);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
void move_r<support, sym_t, pos_t>::SA(std::vector<pos_t>& SA_range, retrieve_params params) const
    requires(support == _locate_rlzdsa)
{
    adjust_retrieve_params(params, n - 1);

    pos_t l = params.l;
    pos_t r = params.r;

    uint16_t p = std::max(
        (uint16_t)1, // use at least one thread
        std::min({
            (uint16_t)omp_get_max_threads(), // use at most all threads
            params.num_threads, // use at most the specified number of threads
            (uint16_t)(((r - l + 1) * (double)r__) / (10.0 * (double)n)) // use at most (r-l+1)*(r/n)*(1/10) threads
        }));

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // iteration range start position
        pos_t b = l + i_p * ((r - l + 1) / p);
        // iteration range end position
        pos_t e = i_p == p - 1 ? r : l + (i_p + 1) * ((r - l + 1) / p) - 1;

        pos_t x_p, x_lp, x_cp, x_r, s_np;

        // current suffix array value
        pos_t s;

        // initialize the rlzdsa context to position b
        init_rlzdsa(b, s, x_p, x_lp, x_cp, x_r, s_np);

        // decode and append SA[b..e] to SA_range
        write_rlzdsa_right(b, e, s, x_p, x_lp, x_cp, x_r, s_np, SA_range, b - l);
    }
}

template <move_r_support support, typename sym_t, typename pos_t>
template <typename output_t, bool output_reversed>
void move_r<support, sym_t, pos_t>::retrieve_range(
    void (move_r<support, sym_t, pos_t>::*retrieve_method)(
        const std::function<void(pos_t, output_t)>&, move_r<support, sym_t, pos_t>::retrieve_params) const,
    std::string file_name, move_r<support, sym_t, pos_t>::retrieve_params params) const
{
    pos_t l = params.l;
    pos_t r = params.r;
    uint16_t num_threads = params.num_threads;
    uint64_t buffer_size_per_thread = std::max<uint64_t>(1024,
        params.max_bytes_alloc != -1 ? params.max_bytes_alloc / num_threads : ((r - l + 1) * sizeof(output_t)) / (num_threads * 500));

    std::filesystem::resize_file(std::filesystem::current_path() / file_name, (r - l + 1) * sizeof(output_t));
    std::vector<sdsl::int_vector_buffer<sizeof(output_t) * 8>> file_bufs;

    for (uint16_t i = 0; i < num_threads; i++) {
        file_bufs.emplace_back(sdsl::int_vector_buffer<sizeof(output_t) * 8>(file_name, std::ios::in, buffer_size_per_thread, sizeof(output_t) * 8, true));
    }

    (this->*retrieve_method)([&](pos_t pos, output_t val) {
        file_bufs[omp_get_thread_num()][pos] = *reinterpret_cast<uint64_t*>(&val);
    }, params);
}