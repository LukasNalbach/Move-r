#include <ips4o.hpp>

template <typename uint_t>
inline uint_t move_data_structure_phi<uint_t>::construction::interval_length_v2v3_seq(lin_node_t_v2v3v4 *ln) {
    return (ln->sc != NULL ? ln->sc->v.first : n) - ln->v.first;
}

template <typename uint_t>
inline typename move_data_structure_phi<uint_t>::construction::lin_node_t_v2v3v4* move_data_structure_phi<uint_t>::construction::is_a_heavy_v2v3v4(
    lin_node_t_v2v3v4 **ln_IpI_,
    uint_t* i_,
    tout_node_t_v2v3v4 *tn_J,
    tout_node_t_v2v3v4 *tn_J_
) {
    // [l,r] = [q_j, q_j + d_j)

    // r + 1
    uint_t rp1 = tn_J_ == NULL ?
        tn_J->v.v.second + interval_length_v2v3_seq(&tn_J->v)
        : tn_J_->v.v.second;

    /* If |[l,r]| < 2a, there cannot be at least 2a input intervals
    connected to [l,r] in the permutation graph. */
    if (rp1-tn_J->v.v.second < two_a) return NULL;

    uint_t i_start = *i_;

    /* Count the number i of input intervals connected to [l,r]
    in the permutation graph and stop as soon as i > a. */
    while (true) {
        if (*i_ > a) break;
        if ((*ln_IpI_)->sc == NULL) return NULL;
        *ln_IpI_ = (*ln_IpI_)->sc;
        (*i_)++;
        if ((*ln_IpI_)->v.first >= rp1) return NULL;
    }

    lin_node_t_v2v3v4 *ln_IpA = *ln_IpI_;

    // Count further and stop as soon as i >= 2a.
    while (true) {
        if (*i_ >= two_a) {
            // if i_ > a+1 held in the beginning, correct ln_IpA
            while (i_start > a+1) {
                ln_IpA = ln_IpA->pr;
                i_start--;
            }
            *i_ = a;
            return ln_IpA;
        }
        if ((*ln_IpI_)->sc == NULL) return NULL;
        *ln_IpI_ = (*ln_IpI_)->sc;
        (*i_)++;
        if ((*ln_IpI_)->v.first >= rp1) return NULL;
    }
}

template <typename uint_t>
void move_data_structure_phi<uint_t>::construction::build_lin_tout_v2v3v4() {
    L_in_v2v3v4.resize(p);
    T_out_v2v3v4.resize(p);

    if (log) log_message("building pi");

    build_pi_for_I();

    if (log) {
        if (mf != NULL) *mf << " time_build_pi=" << time_diff_ns(time);
        time = log_runtime(time);
        log_message("building L_in");
    }

    calculate_seperation_positions_for_I();

    // write I into nodes_v2v3v4 and build L_in_v2v3v4[0..p-1]
    (*reinterpret_cast<std::vector<typename avl_tree<typename doubly_linked_list<std::pair<no_init<uint_t>,no_init<uint_t>>>::doubly_linked_list_node>::avl_node>*>(&nodes_v2v3v4)).resize(k);

    std::vector<std::queue<lin_node_t_v2v3v4*>> Q_i(2*p);

    #pragma omp parallel num_threads(p)
    {
        #pragma omp single
        {
            for (uint16_t i_p=0; i_p<2*p; i_p++) {
                uint_t l2 = x[i_p/2];
                uint_t r2 = x[i_p/2+1];

                if (r2 > l2) {
                    uint_t m2 = l2+(r2-l2)/2;

                    uint_t l = i_p%2 == 0 ? l2 : (m2+1);
                    uint_t r = i_p%2 == 0 ? m2 : r2-1;

                    #pragma omp task
                    {
                        for (uint_t i=l; i<=r; i++) {
                            nodes_v2v3v4[i].v.v = I[i];
                            nodes_v2v3v4[i].v.pr = &nodes_v2v3v4[i-1].v;
                            nodes_v2v3v4[i].v.sc = &nodes_v2v3v4[i+1].v;

                            if (I[i+1].first - I[i].first > l_max) {
                                Q_i[i_p].emplace(&nodes_v2v3v4[i].v);
                            }
                        }
                    }
                }
            }

            #pragma omp taskwait
        }
    }
    
    // Adjust the sizes of the lists in L_in_v2v3v4, and set their head and tail nodes.
    for (uint16_t i=0; i<p; i++) {
        L_in_v2v3v4[i].set_size(x[i+1]-x[i]);

        if (!L_in_v2v3v4[i].empty()) {
            L_in_v2v3v4[i].set_head(&nodes_v2v3v4[x[i]].v);
            L_in_v2v3v4[i].set_tail(&nodes_v2v3v4[x[i+1]-1].v);
            L_in_v2v3v4[i].head()->pr = NULL;
            L_in_v2v3v4[i].tail()->sc = NULL;
        }
    }

    x.clear();
    x.shrink_to_fit();

    // Now, we do not need I anymore.
    I.clear();
    I.shrink_to_fit();

    if (log) {
        if (mf != NULL) *mf << " time_build_lin=" << time_diff_ns(time);
        time = log_runtime(time);
        log_message("building T_out");
    }

    /* This function returns for the value i in [0,k-1] the node in nodes_v2v3v4[0..k-1],
    that stores the pair, which creates the i-th output interval. */
    std::function<tout_node_t_v2v3v4*(uint_t)> at = [this](uint_t i){
        return &nodes_v2v3v4[pi[i]];
    };

    // build T_out_v2v3v4[0..p-1] from nodes_v2v3v4
    #pragma omp parallel num_threads(p)
    {
        #pragma omp single
        {
            for (uint16_t i=0; i<p; i++) {
                
                if (u[i+1] > u[i]) {
                    #pragma omp task
                    {
                        /* Build T_out_v2v3v4[i] out of the pairs, that create the output starting
                        in the range[s[i]..s[i+1]-1], which are located at the positions
                        at[u[i]], at[u[i+1]], ..., at[u[i+1]-1]. */
                        T_out_v2v3v4[i].insert_array(u[i],u[i+1]-1,at,2);
                    }
                }
            }

            #pragma omp taskwait
        }
    }

    u.clear();
    u.shrink_to_fit();

    if (log) {
        if (mf != NULL) *mf << " time_build_tout=" << time_diff_ns(time);
        time = log_runtime(time);
        log_message("splitting too long intervals");
    }

    pi.clear();
    pi.shrink_to_fit();

    // build new_nodes_2v3v4[0..p-1]
    new_nodes_2v3v4.reserve(p);

    for (uint16_t i=0; i<p; i++) {
        new_nodes_2v3v4.emplace_back(dynamic_insert_only_no_copy<tout_node_t_v2v3v4>(k/(double)(16*p*(a-1))));
    }

    /* make sure each avl tree T_out_v2v3v4[i], with i in [0..p-1], contains a pair creating
    an output interval starting at s[i] */
    for (uint16_t i=1; i<p; i++) {
        if (T_out_v2v3v4[i].empty() || T_out_v2v3v4[i].min()->v.v.second != s[i]) {
            lin_node_t_v2v3v4 *ln = &T_out_v2v3v4[i-1].max()->v;

            tout_node_t_v2v3v4 *tn = new_nodes_2v3v4[i].emplace_back(
                tout_node_t_v2v3v4(lin_node_t_v2v3v4(pair_t{
                    ln->v.first+s[i]-ln->v.second,s[i]
                }))
            );

            // find i_ in [0,p-1], so that s[i_] <= tn->v.first < s[i_+1]
            uint16_t i_ = bin_search_max_leq<uint_t>(ln->v.first,0,p-1,[this](uint_t x){return s[x];});

            L_in_v2v3v4[i_].insert_after_node(&tn->v,ln);
            T_out_v2v3v4[i].insert_node(tn);

            if (tn->v.sc == NULL) {
                if (i_ < p-1) {
                    uint16_t i__ = i_+1;

                    while (i__ <= p-1 && L_in_v2v3v4[i__].empty()) {
                        i__++;
                    }

                    if (i__ <= p-1 && !L_in_v2v3v4[i__].empty() && L_in_v2v3v4[i__].head()->v.first - tn->v.v.first > l_max) {
                        Q_i[2*i_].emplace(&tn->v);
                    }
                }
            } else if (tn->v.sc->v.first - tn->v.v.first > l_max) {
                Q_i[2*i_].emplace(&tn->v);
            }
        }
    }

    /* make sure each list L_in_v2v3v4[i], with i in [0..p-1], contains a pair creating an
    input interval starting at s[i] */
    for (uint16_t i=1; i<p; i++) {
        if (L_in_v2v3v4[i].empty() || L_in_v2v3v4[i].head()->v.first != s[i]) {
            lin_node_t_v2v3v4 *ln = L_in_v2v3v4[i-1].tail();

            tout_node_t_v2v3v4 *tn = new_nodes_2v3v4[i].emplace_back(
                tout_node_t_v2v3v4(lin_node_t_v2v3v4(pair_t{
                    s[i],ln->v.second+s[i]-ln->v.first
                }))
            );

            // find i_ in [0,p-1], so that s[i_] <= tn->v.v.second < s[i_+1]
            uint16_t i_ = bin_search_max_leq<uint_t>(tn->v.v.second,0,p-1,[this](uint_t x){return s[x];});

            T_out_v2v3v4[i_].insert_node(tn);
            L_in_v2v3v4[i].push_front_node(&tn->v);

            if (tn->v.sc == NULL) {
                if (i < p-1) {
                    uint16_t i__ = i+1;

                    while (i__ <= p-1 && L_in_v2v3v4[i__].empty()) {
                        i__++;
                    }

                    if (i__ <= p-1 && !L_in_v2v3v4[i__].empty() && L_in_v2v3v4[i__].head()->v.first - tn->v.v.first > l_max) {
                        Q_i[2*i].emplace(&tn->v);
                    }
                }
            } else if (tn->v.sc->v.first - tn->v.v.first > l_max) {
                Q_i[2*i].emplace(&tn->v);
            }
        }
    }

    // insert the pair (s[i+1],s[i+1]) into each avl tree T_out_v2v3v4[i], with i in [0..p-1]
    for (uint16_t i=0; i<p; i++) {
        tout_node_t_v2v3v4* tn = new_nodes_2v3v4[i].emplace_back(
            tout_node_t_v2v3v4(lin_node_t_v2v3v4(pair_t{
                s[i+1],s[i+1]
            }))
        );
        
        L_in_v2v3v4[i].push_back_node(&tn->v);
        T_out_v2v3v4[i].insert_node(tn);
    }

    std::vector<std::vector<std::queue<tout_node_t_v2v3v4*>>> Q_o(p,std::vector<std::queue<tout_node_t_v2v3v4*>>(p));

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        lin_node_t_v2v3v4 *ln_cur,*ln_last;
        tout_node_t_v2v3v4 *tn_tmp;

        for (uint16_t i=2*i_p; i<2*(i_p+1); i++) {
            while (!Q_i[i].empty()) {
                ln_cur = Q_i[i].front();
                Q_i[i].pop();
                
                if (ln_cur->sc != NULL && ln_cur->sc->v.first - ln_cur->v.first > l_max) {
                    // find i_p' in [0,p-1], so that s[i_p'] <= tn->v.v.second < s[i_p'+1]
                    uint16_t i_p_ = bin_search_max_leq<uint_t>(ln_cur->v.second,0,p-1,[this](uint_t x){return s[x];});

                    do {
                        ln_last = ln_cur;

                        tn_tmp = new_nodes_2v3v4[i_p].emplace_back(
                            tout_node_t_v2v3v4(lin_node_t_v2v3v4(pair_t{
                                ln_last->v.first+l_max,
                                ln_last->v.second+l_max
                            }))
                        );

                        ln_cur = &tn_tmp->v;
                        L_in_v2v3v4[i_p].insert_after_node(ln_cur,ln_last);

                        while (ln_cur->v.second >= s[i_p_+1]) {
                            i_p_++;
                        }
                        
                        Q_o[i_p_][i_p].emplace(tn_tmp);
                    } while (ln_cur->sc != NULL && ln_cur->sc->v.first - ln_cur->v.first > l_max);
                }
            }
        }
    }

    for (uint16_t i=0; i<p; i++) {
        L_in_v2v3v4[i].remove_node(L_in_v2v3v4[i].tail());
    }

    Q_i.clear();

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        for (uint16_t i=0; i<p; i++) {
            while (!Q_o[i_p][i].empty()) {
                T_out_v2v3v4[i_p].insert_node(Q_o[i_p][i].front());
                Q_o[i_p][i].pop();
            }
        }
    }

    Q_o.clear();

    if (log) {
        if (mf != NULL) *mf << " time_split_too_long_input_intervals=" << time_diff_ns(time);
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_data_structure_phi<uint_t>::construction::build_dp_dq_v2v3v4() {
    s.clear();
    s.shrink_to_fit();
    
    // [0..p-1], x[i] stores the number of input intervals in I starting before s[i]
    x.resize(p+1);
    x[0] = 0;

    for (uint16_t i=0; i<p; i++) {
        // remove the pairs (s[i+1],s[i+1]) from T_out_v2v3v4[i], for each i in [0..p-1]
        T_out_v2v3v4[i].remove_node(T_out_v2v3v4[i].max());

        // calculate x
        x[i+1] = x[i] + L_in_v2v3v4[i].size();
    }

    k_ = x[p];

    if (log) {
        float k__k = std::round(100.0*k_/k)/100.0;
        if (mf != NULL) {
            *mf << " k=" << k;
            *mf << " k_=" << k_;
        }
        std::cout << "k' = " << k_ << ", k'/k = " << k__k << std::endl;
        log_message("building D_p and D_q");
    }

    mds.resize(n,k_);

    no_init_resize(D_q,k_+1);
    D_q[k_] = n;

    #pragma omp parallel num_threads(p)
    {
        #pragma omp single
        {
            for (uint16_t i_p=0; i_p<p; i_p++) {
                uint_t l = x[i_p];
                uint_t r = x[i_p+1];

                auto ln_cur = L_in_v2v3v4[i_p].head();

                for (uint_t i=l; i<r; i++) {
                    mds.set_p(i,ln_cur->v.first);
                    D_q[i] = ln_cur->v.second;
                    ln_cur = ln_cur->sc;
                }
            }

            #pragma omp taskwait
        }
    }

    // Deconstruct the additional data structures.
    x.clear();
    x.shrink_to_fit();

    for (uint16_t i=0; i<p; i++) {
        L_in_v2v3v4[i].disconnect_nodes();
        T_out_v2v3v4[i].disconnect_nodes();
        new_nodes_2v3v4[i].clear();
    }
    
    L_in_v2v3v4.clear();
    L_in_v2v3v4.shrink_to_fit();

    T_out_v2v3v4.clear();
    T_out_v2v3v4.shrink_to_fit();

    new_nodes_2v3v4.clear();
    new_nodes_2v3v4.shrink_to_fit();

    nodes_v2v3v4.clear();
    nodes_v2v3v4.shrink_to_fit();

    if (log) {
        if (mf != NULL) *mf << " time_build_dp_dq=" << time_diff_ns(time);
        time = log_runtime(time);
    }
}