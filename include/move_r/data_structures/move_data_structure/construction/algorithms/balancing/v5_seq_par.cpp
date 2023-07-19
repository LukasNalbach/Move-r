template <typename uint_t>
inline uint_t move_data_structure_phi<uint_t>::construction::is_a_heavy_v5_seq_par(tin_it_t_v5& tn_I, uint_t q_J_) {
    // current number of input interval starting in [q_j, q_j + d_j)
    uint16_t e = 1;

    // count the number of input interval starting in [q_j, q_j + d_j) up to a
    while (e <= a && (*tn_I).first < q_J_) {
        tn_I++;
        e++;
    }

    // if there are less than a input interval starting in [q_j, q_j + d_j), then it is a-balanced, so return 0
    if ((e <= a) || q_J_ <= (*tn_I).first) {
        return 0;
    } else {
        // else tn_I points to the a+1-th input interval starting in [q_j, q_j + d_j)
        // and since q_j + d = p_{i+a} we can now set qj_pd
        uint_t qj_pd = (*tn_I).first;

        // count the number of input interval starting in [q_j, q_j + d_j) up to 2a
        while (e < two_a && (*tn_I).first < q_J_) {
            tn_I++;
            e++;
        }

        // if there are less than 2a input interval starting in [q_j, q_j + d_j), then it is a-balanced, so return 0
        if (e < two_a || q_J_ <= (*tn_I).first) {
            return 0;
        } else {
            // else return q_j + d
            return qj_pd;
        }
    }
}

template <typename uint_t>
inline move_data_structure_phi<uint_t>::construction::tout_it_t_v5 move_data_structure_phi<uint_t>::construction::balance_upto_v5_seq_par(tout_it_t_v5& tn_J_, uint_t qj_pd, uint_t q_u) {
    // Index in [0..p-1] of the current thread.
    uint16_t i_p = omp_get_thread_num();

    uint_t q_J_ = (*tn_J_).second; // q_j'
    tout_it_t_v5 tn_J = tn_J_; // iterator pointing to the pair (p_j',q_j') in T_out_v5[i_p]
    tn_J--;
    uint_t q_j = (*tn_J).second; // iterator pointing to the pair (p_j,q_j) in T_out_v5[i_p]
    uint_t pj_pd = (*tn_J).first + (qj_pd - q_j); // p_j + d
    pair_t pr_new{pj_pd,qj_pd}; // the newly created pair

    // insert the newly created pair into the current thread's tree in T_out_v5
    tout_it_t_v5 tout_n_new = T_out_v5[i_p].emplace_hint(tn_J_,pr_new);

    // check, if the newly created has to be inserted into a tree in T_in_v5 of another thread
    if (p != 1 && !(s[i_p] <= pj_pd && pj_pd < s[i_p+1])) {
        uint16_t l = 0;
        uint16_t r = p-1;
        uint16_t m;

        // calculate the index i_p' of the thread, into whiches tree in T_in_v5, the newly created pair has to be inserted
        while (l != r) {
            m = l+(r-l)/2+1;
            if (pj_pd >= s[m]) {
                l = m;
            } else {
                r = m-1;
            }
        }

        // store the newly created pair in Q_v5[i_p'][i_p]
        Q_v5[l][i_p].emplace_back(pr_new);
    } else {
        // else the newly created pair can be inserted into the current thread's tree in T_in_v5
        tin_it_t_v5 tin_n_new = node(T_in_v5[i_p].emplace(pr_new));

        // check whether the output interval containing p_j + d can have become unbalanced by splitting [p_j, p_j + d_j)
        if (pj_pd < q_u && (pj_pd < q_j || q_J_ <= pj_pd)) {
            // if yes, find the node in T_out_v5[i_p], that creates the pair (p_y,q_y), where p_j + d in [q_j, q_y + d_y)
            tout_it_t_v5 tn_Y = T_out_v5[i_p].lower_bound(pair_t{0,pj_pd});

            if (pj_pd < (*tn_Y).second) {
                tn_Y--;
            }

            // find the node in T_in_v5[i_p], that creates the first input interval starting in [q_j, q_y + d_y)
            while ((*tn_Y).second < (*tin_n_new).first) {
                tin_n_new--;
            }

            if ((*tin_n_new).first < (*tn_Y).second) {
                tin_n_new++;
            }

            // iterate one step with tn_Y, s.t. it now points to the output interval starting direclty after [q_y, q_y + d_y)
            tn_Y++;
            uint_t qy_pd_; // q_y + d', where d' = p_{z+a} - q_y

            // check if [q_y, q_y + d_y) is a-heavy and balanced
            if ((qy_pd_ = is_a_heavy_v5_seq_par(tin_n_new,(*tn_Y).second))) {
                // if yes, balance it and all a-heavy output intervals in [s[i_p],s[i_p+1]) starting before q_u that
                // become a-heavy in the process
                balance_upto_v5_seq_par(tn_Y,qy_pd_,q_u);
                
                // because we inserted another pair into T_out_v5[i_p] in the recursive call of balance_upto_v5_seq_par, tout_n_new
                // may not point to (p_j + d, q_j + d) anymore, so set it to T_out_v5[i_p].end() (which is constant)
                tout_n_new = T_out_v5[i_p].end();
            }
        }
    }

    // return an iterator pointing to the newly created pair (p_j + d, q_j + d) in T_out_v5, if no recursive call has been
    // made; else an iterator pointing to T_out_v5[i_p].end()
    return tout_n_new;
}

template <typename uint_t>
void move_data_structure_phi<uint_t>::construction::balance_v5_seq_par() {
    if (log) {
        std::string msg = "balancing";
        if (p > 1) msg.append(" (phase 1)");
        log_message(msg);
    }

    if (p > 1) {
        Q_v5.resize(p,std::vector<pair_arr_t>(p));
        Q_v5_.resize(p,std::vector<pair_arr_t>(p));

        for (uint16_t i=0; i<p; i++) {
            for (uint16_t j=0; j<p; j++) {
                Q_v5[i][j].reserve(k/(32*a*p*p));
                Q_v5_[i][j].reserve(k/(32*a*p*p));
            }
        }
    }

    // first phase of the balancing algorithm
    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        tin_it_t_v5 tn_I = T_in_v5[i_p].begin(); // iterator pointing to the pair in T_in_v5[i_p], that creates (p_i,q_i)
        tout_it_t_v5 tn_J = T_out_v5[i_p].begin(); // iterator pointing to the pair in T_out_v5[i_p], that creates (p_j,q_j)
        /* iterator pointing to the pair in T_out_v5[i_p], that creates (p_j',q_j'), where [q_j', q_j' + d_j') is the output interval,
           that starts direclty after [q_j, q_j + d_j) */
        tout_it_t_v5 tn_J_ = tn_J; 
        tn_J_++;
        uint_t qj_pd; // q_j + d (temporary variable)
        bool stop = false;

        while (!stop) {
            // check if the current output interval [q_j, q_j + d_j) is a-heavy
            if ((qj_pd = is_a_heavy_v5_seq_par(tn_I,(*tn_J_).second))) {
                // if yes, balance it and all a-heavy output intervals in [s[i_p],s[i_p+1]) starting before q_j + d that
                // become a-heavy in the process
                tn_J = balance_upto_v5_seq_par(tn_J_,qj_pd,qj_pd);
                
                // because we inserted a pair into T_in_v5[i_p], the iterator tn_I may now be invalid, so reset it
                tn_I = T_in_v5[i_p].find(pair_t{qj_pd,0});

                // if tn_J points to T_out_v5[i_p].end(), balance_upto_v5_seq_par made a recursive call, so reset tn_J
                if (tn_J == T_out_v5[i_p].end()) {
                    tn_J = T_out_v5[i_p].find(pair_t{0,qj_pd});
                }

                // iterate to the next output interval
                tn_J_ = tn_J;
                tn_J_++;

                // now tn_I points to (p_{i+a}, q_{i+a}) and tn_J points to (p_j + d, q_j + d), so we can directly start
                // the next iteration
                continue;
            }

            // else, find the next output interval that contains an input interval
            do {
                tn_J = tn_J_;
                tn_J_++;
                if (tn_J_ == T_out_v5[i_p].end()) {stop = true; break;}
                while (!stop && (*tn_I).first < (*tn_J).second) {
                    tn_I++;
                    if (tn_I == T_in_v5[i_p].end()) {stop = true; break;}
                }
            } while (!stop && (*tn_I).first >= (*tn_J_).second);
        }
    }

    if (log) {
        if (mf != NULL) {
            *mf << " time_balance_phase_1=" << time_diff_ns(time);
            if (p == 1) *mf << " time_balance_phase_2=0" << time_diff_ns(time);
        }
        time = log_runtime(time);
    }

    // second phase of the balancing algorithm
    if (p > 1) {
        if (log) log_message("balancing (phase 2)");

        bool done;

        #pragma omp parallel num_threads(p)
        {
            // Index in [0..p-1] of the current thread.
            uint16_t i_p = omp_get_thread_num();

            uint_t qy_pd; // q_y + d
            tin_it_t_v5 tn_new = T_in_v5[i_p].end(); // iterator pointing to the newly created pair in T_out_v5[i_p]
            tout_it_t_v5 tn_Y = T_out_v5[i_p].end(); // iterator pointing to the pair (p_y, q_y) in T_out_v5[i_p]

            while (true) {
                #pragma omp barrier

                #pragma omp single
                {
                    // the pairs that have been inserted into Q_v5 in the first phase (last iteration of the second phase)
                    // now have to be inserted into T_in_v5[0..p-1], so swap Q_v5 with Q_v5_
                    std::swap(Q_v5,Q_v5_);
                    done = true;
                }

                #pragma omp barrier

                // check whether Q_v5_ is empty
                for (uint16_t i_p_=0; i_p_<p; i_p_++) {
                    if (!Q_v5_[i_p][i_p_].empty()) {
                        done = false;
                        break;
                    }
                }

                #pragma omp barrier
                
                // if Q_v5_ is empty, there are no a-heavy output intervals, so break
                if (done) {
                    break;
                }

                // iterate over all pairs to insert into T_in_v5[i_p]
                for (pair_arr_t& vec : Q_v5_[i_p]) {
                    for (pair_t& pr_new : vec) { // a newly created pair (p_j + d, q_j + d)
                        // insert (p_j + d, q_j + d) into T_in_v5[i_p]
                        tn_new = node(T_in_v5[i_p].emplace(pr_new));

                        // find the pair in T_out_v5[i_p] containing (p_j + d, q_j + d)
                        tn_Y = T_out_v5[i_p].lower_bound(pair_t{0,pr_new.first});

                        if ((*tn_new).first < (*tn_Y).second) {
                            tn_Y--;
                        }

                        // find the node in T_in_v5[i_p], that creates the first input interval starting in [q_j, q_y + d_y)
                        while ((*tn_Y).second < (*tn_new).first) {
                            tn_new--;
                        }

                        if ((*tn_new).first < (*tn_Y).second) {
                            tn_new++;
                        }

                        // iterate one step with tn_Y, s.t. it now points to the output interval starting direclty after [q_y, q_y + d_y)
                        tn_Y++;

                        // check if [q_y, q_y + d_y) is a-heavy and balanced
                        if ((qy_pd = is_a_heavy_v5_seq_par(tn_new,(*tn_Y).second))) {
                            // if yes, balance it and all a-heavy output intervals in [s[i_p],s[i_p+1]), that become a-heavy in the process
                            balance_upto_v5_seq_par(tn_Y,qy_pd,n);
                        }
                    }

                    vec.clear();
                }
            }
        }
        
        Q_v5.clear();
        Q_v5.shrink_to_fit();

        Q_v5_.clear();
        Q_v5_.shrink_to_fit();

        if (log) {
            if (mf != NULL) *mf << " time_balance_phase_2=" << time_diff_ns(time);
            time = log_runtime(time);
        }
    }

    T_out_v5.clear();
    T_out_v5.shrink_to_fit();

    s.clear();
    s.shrink_to_fit();
}