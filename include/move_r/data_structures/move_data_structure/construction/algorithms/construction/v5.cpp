template <typename uint_t>
void move_data_structure_phi<uint_t>::construction::build_tin_tout_v5() {
    if (log) log_message("building pi");
        
    // build pi to construct T_out faster
    build_pi_for_I();
    calculate_seperation_positions_for_I();

    if (log) {
        if (mf != NULL) *mf << " time_build_pi=" << time_diff_ns(time);
        time = log_runtime(time);
        log_message("building T_out");
    }

    T_out_v5.resize(p);

    // build T_out_v5[0..p-1]
    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        uint_t b = u[i_p];
        uint_t e = u[i_p+1];

        for (uint_t i=b; i<e; i++) {
            T_out_v5[i_p].emplace_hint(T_out_v5[i_p].end(),I[pi[i]]);
        }
    }

    u.clear();
    u.shrink_to_fit();

    pi.clear();
    pi.shrink_to_fit();

    if (log) {
        if (mf != NULL) *mf << " time_build_tout=" << time_diff_ns(time);
        time = log_runtime(time);
        log_message("building T_in");
    }

    T_in_v5.resize(p);

    // build T_in_v5[0..p-1]
    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        if (x[i_p] < x[i_p+1]) {
            T_in_v5[i_p].insert(&I[x[i_p]],&I[x[i_p+1]]);
        }
    }

    x.clear();
    x.shrink_to_fit();

    if (delete_i) {
        // Now, we do not need I anymore.
        I.clear();
        I.shrink_to_fit();
    }

    // make sure there is an input interval starting at s[0], s[1], ... s[p-1]
    for (uint16_t i=1; i<p; i++) {
        if (T_in_v5[i].empty() || (*T_in_v5[i].begin()).first != s[i]) {
            pair_t pr_split = *T_in_v5[i-1].rbegin();
            pair_t pr_new{s[i], pr_split.second + (s[i] - pr_split.first)};

            uint16_t l = 0;
            uint16_t r = p-1;
            uint16_t m;

            while (l != r) {
                m = l+(r-l)/2+1;
                if (s[m] > pr_new.second) {
                    r = m-1;
                } else {
                    l = m;
                }
            }

            T_in_v5[i].emplace(pr_new);
            T_out_v5[l].emplace(pr_new);
        }
    }

    // make sure there is an output interval starting at s[0], s[1], ... s[p-1]
    for (uint16_t i=1; i<p; i++) {
        if (T_out_v5[i].empty() || (*T_out_v5[i].begin()).second != s[i]) {
            pair_t pr_split = *T_out_v5[i-1].rbegin();
            pair_t pr_new{pr_split.first + (s[i] - pr_split.second), s[i]};

            uint16_t l = 0;
            uint16_t r = p-1;
            uint16_t m;

            while (l != r) {
                m = l+(r-l)/2+1;
                if (s[m] > pr_new.first) {
                    r = m-1;
                } else {
                    l = m;
                }
            }

            T_in_v5[l].emplace(pr_new);
            T_out_v5[i].emplace(pr_new);
        }
    }

    // add dummy pairs (s[i+1],s[i+1]) to the i-th trees, for i in [0..p-1]
    // this enables us to calculate the length of the last input- and output intervals in each section
    for (uint16_t i=0; i<p; i++) {
        T_in_v5[i].emplace(pair_t{s[i+1],s[i+1]});
        T_out_v5[i].emplace(pair_t{s[i+1],s[i+1]});
    }

    if (log) {
        if (mf != NULL) *mf << " time_build_tin=" << time_diff_ns(time);
        time = log_runtime(time);
        log_message("splitting too long intervals");
    }

    T_out_temp_v5.resize(p,std::vector<tout_t_v5>(p));

    // iterate over all trees in T_in_v5 and split every input interval that is longer than l_max
    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        tin_it_t_v5 tn_I = T_in_v5[i_p].begin();
        pair_t pr_last = *tn_I;
        tn_I++;
        uint_t l,m,r;
        
        while (tn_I != T_in_v5[i_p].end()) {
            while (((*tn_I).first - pr_last.first) > l_max) {
                tn_I = T_in_v5[i_p].emplace_hint(tn_I,pair_t{pr_last.first + l_max, pr_last.second + l_max});
                pr_last = *tn_I;
                tn_I++;

                l = 0;
                r = p-1;

                while (l != r) {
                    m = (l+r)/2+1;
                    if (s[m] <= pr_last.second) {
                        l = m;
                    } else {
                        r = m-1;
                    }
                }

                // store each new pair in its corresponding tree in T_out_temp_v5
                T_out_temp_v5[l][i_p].emplace(pr_last);
            }

            pr_last = *tn_I;
            tn_I++;
        }
    }

    // merge T_out_v5 with T_out_temp_v5
    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        for (uint16_t i=0; i<p; i++) {
            T_out_v5[i_p].merge(T_out_temp_v5[i_p][i]);
        }
    }

    T_out_temp_v5.clear();
    T_out_temp_v5.shrink_to_fit();

    if (log) {
        if (mf != NULL) *mf << " time_split_too_long_input_intervals=" << time_diff_ns(time);
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_data_structure_phi<uint_t>::construction::build_dp_dq_v5() {
    // recalculate the indices x of the seperation positions in the interval sequence
    x.resize(p+1);
    x[0] = 0;

    for (uint16_t i=0; i<p; i++) {
        T_in_v5[i].erase(*T_in_v5[i].rbegin());
        x[i+1] = x[i] + T_in_v5[i].size();
    }

    k_ = x[p];

    // resize the interleaved vectors in the move data structure
    mds.resize(n,k_);
    
    if (log) {
        float k__k = std::round(100.0*k_/k)/100.0;
        if (mf != NULL) {
            *mf << " k=" << k;
            *mf << " k_=" << k_;
        }
        std::cout << "k' = " << k_ << ", k'/k = " << k__k << std::endl;
        log_message("building D_p and D_q");
    }

    no_init_resize(D_q,k_+1);
    D_q[k_] = n;

    // write the input interval starting positions to D_p (in the move data structure) and
    // write the output interval starting positions to D_q
    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        uint_t b = x[i_p];
        uint_t e = x[i_p+1];

        tin_it_t_v5 tin_it = T_in_v5[i_p].begin();

        for (uint_t i=b; i<e; i++) {
            mds.set_p(i,(*tin_it).first);
            D_q[i] = (*tin_it).second;
            tin_it++;
        }
    }

    x.clear();
    x.shrink_to_fit();

    T_in_v5.clear();
    T_in_v5.shrink_to_fit();

    if (log) {
        if (mf != NULL) *mf << " time_build_dp_dq=" << time_diff_ns(time);
        time = log_runtime(time);
    }
}