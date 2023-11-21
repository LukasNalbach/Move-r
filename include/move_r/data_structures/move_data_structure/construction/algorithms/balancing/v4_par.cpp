template <typename uint_t>
inline typename move_data_structure<uint_t>::construction::tout_node_t_v2v3v4* move_data_structure<uint_t>::construction::balance_upto_v4_par(
    lin_node_t_v2v3v4 *ln_IpA,
    tout_node_t_v2v3v4 *tn_J,
    tout_node_t_v2v3v4* tn_J_,
    uint_t q_u,
    uint_t p_cur,
    uint_t *i_
) {
    // Index in [0..p-1] of the current thread.
    uint16_t i_p = omp_get_thread_num();

    uint_t p_j = tn_J->v.v.first;
    uint_t q_j = tn_J->v.v.second;
    uint_t d_j = tn_J_->v.v.second - q_j;

    // d = p_{i+2a} - q_j is the maximum integer, so that [q_j, q_j + d) has a incoming edges in the permutation graph.
    uint_t d = ln_IpA->v.first - q_j;

    // Create the pair (p_j + d, q_j + d), which creates two new input intervals [p_j, p_j + d) and [p_j + d, p_j + d_j).
    tout_node_t_v2v3v4 *tn_NEW = new_nodes_2v3v4[i_p].emplace_back(tout_node_t_v2v3v4(lin_node_t_v2v3v4(pair_t{p_j + d, q_j + d})));
    T_out_v2v3v4[i_p].insert_hint(tn_NEW,tn_J);

    if (!(s[i_p] <= p_j + d && p_j + d < s[i_p+1])) {
        // If the new pair must be inserted in L_in_v2v3v4[i_p_] of another thread i_p_ != i_p, find i_p_ with a binary search.
        uint16_t i_p_ = bin_search_max_leq<uint_t>(p_j+d,0,p-1,[this](uint_t x){return s[x];});
        Q_v4[i_p_][i_p].enqueue(q_node_t_v34{&tn_NEW->v,&tn_J->v});
    } else {
        // Else insert it in L_in_v2v3v4[i_p].
        L_in_v2v3v4[i_p].insert_after_node(&tn_NEW->v,&tn_J->v);

        if (p_j + d < q_u) {
            if (p_j + d < q_j || q_j + d_j <= p_j + d) {
                tout_node_t_v2v3v4 *tn_Y = T_out_v2v3v4[i_p].max_leq(lin_node_t_v2v3v4(pair_t{0,p_j + d}));
                uint_t q_y = tn_Y->v.v.second;

                // find the output interval starting after [q_y, q_y + d_y)
                tout_node_t_v2v3v4 *tn_Yp1 = tn_Y->nxt();

                // find the first input interval [p_z, p_z + d_z) connected to [q_y, q_y + d_y) in the permutation graph
                lin_node_t_v2v3v4 *ln_Z = &tn_NEW->v;
                uint_t i__ = 1;
                while (ln_Z->pr != NULL && ln_Z->pr->v.first >= q_y) {
                    ln_Z = ln_Z->pr;
                    i__++;
                }
                ln_Z = &tn_NEW->v;

                lin_node_t_v2v3v4 *ln_ZpA = is_a_heavy_v2v3v4(&ln_Z,&i__,tn_Y,tn_Yp1);
                if (ln_ZpA != NULL) {
                    balance_upto_v4_par(ln_ZpA,tn_Y,tn_Yp1,q_u,p_cur,i_);
                }
            }
        } else if (p_j + d < p_cur) {
            (*i_)++;
        }
    }
    return tn_NEW;
}

/**
 * @brief balances the disjoint interval sequence in L_in_v2v3v4[0..p-1] and T_out_v2v3v4[0..p-1] in parallel
 */
template <typename uint_t>
void move_data_structure<uint_t>::construction::balance_v4_par() {
    if (log) log_message("balancing (phase 1)");

    Q_v4.resize(p);

    /** @brief [0..p-1] D_done[i] stores whether thread i is done */
    std::vector<bool> D_done(p,false);

    for (uint16_t i=0; i<p; i++) {
        Q_v4[i] = std::vector<moodycamel::ConcurrentQueue<q_node_t_v34>>(p);
        D_done[i] = false;
    }

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // points to to the pair (p_i,q_i).
        lin_node_t_v2v3v4 *ln_I = L_in_v2v3v4[i_p].head();
        // points to the pair (p_j,q_j).
        typename tout_t_v2v3v4::avl_it tn_J = T_out_v2v3v4[i_p].iterator();
        // points to the pair (p_{j'},q_{j'}), where q_j + d_j = q_{j'}.
        typename tout_t_v2v3v4::avl_it tn_Jp1 = T_out_v2v3v4[i_p].iterator(T_out_v2v3v4[i_p].second_smallest());

        // temporary variables
        lin_node_t_v2v3v4 *ln_IpA;
        uint_t i_ = 1;

        // At the start of each iteration, [p_i, p_i + d_i) is the first input interval connected to [q_j, q_j + d_j) in the permutation graph
        bool stop = false;
        while (!stop) {
            ln_IpA = is_a_heavy_v2v3v4(&ln_I,&i_,tn_J.current(),tn_Jp1.current());

            // If [q_j, q_j + d_j) is a-heavy, balance it and all output intervals that start before it and become a-heavy in the process.
            if (ln_IpA != NULL) {
                tn_J.set(balance_upto_v4_par(ln_IpA,tn_J.current(),tn_Jp1.current(),tn_J.current()->v.v.second,ln_I->v.first,&i_));
                continue;
            }

            // Find the next output interval with an incoming edge in the permutation graph and the first input interval connected to it.
            do {
                if (!tn_Jp1.has_next()) {stop = true; break;}
                tn_J.set(tn_Jp1.current());
                tn_Jp1.next();
                while (ln_I->v.first < tn_J.current()->v.v.second) {
                    if (ln_I->sc == NULL) {stop = true; break;}
                    ln_I = ln_I->sc;
                }
            } while (!stop && ln_I->v.first >= tn_Jp1.current()->v.v.second);
            i_ = 1;
        }
    }

    if (log) {
        if (mf != NULL) *mf << " time_balance_phase_1=" << time_diff_ns(time);
        time = log_runtime(time);
        log_message("balancing (phase 2)");
    }

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // temporary variables
        q_node_t_v34 pair_ins;
        lin_node_t_v2v3v4 *ln_I,*ln_Im1,*ln_Z,*ln_ZpA;
        tout_node_t_v2v3v4 *tn_Y,*tn_Yp1;
        bool done_this,done_other;
        uint_t i_ = 1;
        uint_t q_y;

        done_this = false;
        done_other = false;

        while (!done_other || !done_this) {
            done_this = true;

            // Consider all ((p_i,q_i),(p_{i-1},q_{i-1})), where (p_i,q_i) must be inserted into L_in[i_p] after (p_{i-1},q_{i-1})
            for (uint16_t i=0; i<p; i++) {
                while (Q_v4[i_p][i].try_dequeue(pair_ins)) {
                    if (done_this) {
                        done_this = false;
                        if (D_done[i_p]) {
                            D_done[i_p] = false;
                        }
                    }

                    ln_I = pair_ins.first;
                    ln_Im1 = pair_ins.second;

                    // Insert (p_i,q_i) into L_in[i_p] after (p_{i-1},q_{i-1})
                    L_in_v2v3v4[i_p].insert_after_node(ln_I,ln_Im1);

                    // check if an output interval could have become a-heavy by inserting (p_i,q_i)
                    tn_Y = T_out_v2v3v4[i_p].max_leq(lin_node_t_v2v3v4(pair_t{0,ln_I->v.first}));
                    q_y = tn_Y->v.v.second;

                    // find the output interval starting after [q_y, q_y + d_y)
                    tn_Yp1 = tn_Y->nxt();

                    // find the first input interval [p_z, p_z + d_z) connected to [q_y, q_y + d_y) in the permutation graph
                    ln_Z = ln_I;
                    i_ = 1;
                    while (ln_Z->pr != NULL && ln_Z->pr->v.first >= q_y) {
                        ln_Z = ln_Z->pr;
                        i_++;
                    }
                    ln_Z = ln_I;

                    /* check if [q_y, q_y + d_y) is a-heavy and if yes, balance it and all output intervals starting before it,
                       that might get a-heavy in the process. */
                    ln_ZpA = is_a_heavy_v2v3v4(&ln_Z,&i_,tn_Y,tn_Yp1);
                    if (ln_ZpA != NULL) {
                        balance_upto_v4_par(ln_ZpA,tn_Y,tn_Yp1,s[i_p+1],s[i_p+1],&i_);
                    }
                }
            }

            if (done_this) {
                done_other = true;
                for (uint16_t i=0; i<p; i++) {
                    if (!D_done[i] || Q_v4[i][i_p].size_approx() != 0) {
                        done_other = false;
                    }
                }

                if (!D_done[i_p]) {
                    D_done[i_p] = true;
                }
            }
        }
    }
    
    Q_v4.clear();
    Q_v4.shrink_to_fit();

    if (log) {
        if (mf != NULL) *mf << " time_balance_phase_2=" << time_diff_ns(time);
        time = log_runtime(time);
    }
}