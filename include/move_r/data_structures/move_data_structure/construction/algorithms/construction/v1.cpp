template <typename uint_t>
void move_data_structure_phi<uint_t>::construction::build_tin_tout_v1() {
    if (log) log_message("building T_in");

    T_in_v1 = avl_tree<std::pair<uint_t,uint_t>>(
        [](auto n1, auto n2){return n1.first < n2.first;},
        [](auto n1, auto n2){return n1.first > n2.first;},
        [](auto n1, auto n2){return n1.first == n2.first;}
    );

    T_out_v1 = avl_tree<std::pair<uint_t,uint_t>>(
        [](auto n1, auto n2){return n1.second < n2.second;},
        [](auto n1, auto n2){return n1.second > n2.second;},
        [](auto n1, auto n2){return n1.second == n2.second;}
    );

    T_e_v1 = avl_tree<std::pair<uint_t,uint_t>>(
        [](auto n1, auto n2){return n1.first < n2.first;},
        [](auto n1, auto n2){return n1.first > n2.first;},
        [](auto n1, auto n2){return n1.first == n2.first;}
    );

    // build T_in_v1 and T_out_v1
    for (uint_t i=0; i<k; i++) {
        T_in_v1.insert_or_update(I[i]);
    }

    T_in_v1.insert_or_update(std::make_pair(n,n));

    if (log) {
        if (mf != NULL) *mf << " time_build_tin=" << time_diff_ns(time);
        time = log_runtime(time);
        log_message("building T_out");
    }

    for (uint_t i=0; i<k; i++) {
        T_out_v1.insert_or_update(I[i]);
    }

    T_out_v1.insert_or_update(std::make_pair(n,n));

    avl_node<std::pair<uint_t,uint_t>> *node_cur = T_in_v1.minimum();

    while (node_cur != T_in_v1.maximum()) {
        while (node_cur->nxt()->v.first - node_cur->v.first > l_max) {
            T_out_v1.insert_or_update(std::make_pair(node_cur->v.first+l_max,node_cur->v.second+l_max));
            node_cur = T_in_v1.insert_or_update(std::make_pair(node_cur->v.first+l_max,node_cur->v.second+l_max));
        }
        node_cur = node_cur->nxt();
    }

    // Now, we do not need I anymore.
    I.clear();
    I.shrink_to_fit();

    if (log) {
        if (mf != NULL) *mf << " time_build_tout=" << time_diff_ns(time);
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_data_structure_phi<uint_t>::construction::build_dp_dq_v1() {
    if (log) log_message("building D_p and D_q");

    mds.resize(n,k_);
    
    (*reinterpret_cast<std::vector<no_init<uint_t>>*>(&D_q)).resize(k_+1);
    D_q[k_] = n;

    auto it = T_in_v1.iterator();

    for (uint_t i=0; i<=k_; i++) {
        mds.set_p(i,it.current()->v.first);
        D_q[i] = it.current()->v.second;
        it.next();
    }

    T_in_v1.delete_nodes();

    if (log) {
        if (mf != NULL) *mf << " time_write_dp_dq=" << time_diff_ns(time);
        time = log_runtime(time);
    }
}

template <typename uint_t>
void move_data_structure_phi<uint_t>::construction::build_didx_doffs_v1() {
    if (log) log_message("building D_idx and D_offs");

    uint_t l;
    uint_t r;
    uint_t m;
    
    for (uint_t j=0; j<k_; j++) {
        /* For each output interval [q_j, q_j + d_j), find the input interval [p_i, p_i + d_i),
        q_j is in and set D_idx[j] = i. Find the maximum integer i in [0,k'-1], so that p_i <= q_j
        with a binary search over D_pair. */

        l = 0;
        r = k_-1;

        while (l != r) {
            m = (l+r)/2+1;
            
            if (mds.p(m) > D_q[j]) {
                r = m-1;
            } else {
                l = m;
            }
        }

        mds.set_idx(j,l);
        mds.set_offs(j,D_q[j]-mds.p(l));
    }

    if (log) {
        if (mf != NULL) *mf << " time_build_didx_doffs=" << time_diff_ns(time);
        time = log_runtime(time);
    }
}