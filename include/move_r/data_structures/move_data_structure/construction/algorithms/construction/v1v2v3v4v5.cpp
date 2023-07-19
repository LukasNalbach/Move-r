#include <ips4o.hpp>

template <typename uint_t>
void move_data_structure_phi<uint_t>::construction::build_pi_for_I() {
    (*reinterpret_cast<std::vector<no_init<uint_t>>*>(&pi)).resize(k);
    
    // write the identity permutation of [0..k-1] to pi
    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<k; i++) {
        pi[i] = i;
    }

    // sort pi by the output interval starting positions in I
    auto comp_pi = [this](uint_t i, uint_t j){return I[i].second < I[j].second;};
    if (p > 1) {
        ips4o::parallel::sort(&pi[0],&pi[k],comp_pi);
    } else {
        ips4o::sort(&pi[0],&pi[k],comp_pi);
    }
}

template <typename uint_t>
void move_data_structure_phi<uint_t>::construction::build_pi_for_dq() {
    (*reinterpret_cast<std::vector<no_init<uint_t>>*>(&pi)).resize(k_+1);

    // write the identity permutation of [0..k'] to pi
    #pragma omp parallel for num_threads(p)
    for (uint64_t i=0; i<=k_; i++) {
        pi[i] = i;
    }

    // sort pi by D_q
    auto comp_pi = [this](uint_t i, uint_t j){return D_q[i] < D_q[j];};
    if (p > 1) {
        ips4o::parallel::sort(&pi[0],&pi[k_],comp_pi);
    } else {
        ips4o::sort(&pi[0],&pi[k_],comp_pi);
    }
}

template <typename uint_t>
void move_data_structure_phi<uint_t>::construction::calculate_seperation_positions_for_I() {
    s.resize(p+1);
    s[0] = 0;
    s[p] = n;

    x.resize(p+1);
    x[0] = 0;
    x[p] = k;

    u.resize(p+1);
    u[0] = 0;
    u[p] = k;

    // calculate seperation positions
    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // The optimal value i_p * lfloor 2k/p rfloor for s[i_p].
        uint_t o = i_p*((2*k)/p);

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
            /* Set the Candidate position for s[i_p] to the position
            in the middle between l_s and r_s. */
            m_s = l_s+(r_s-l_s)/2;

            // Find the minimum x' in [0,k-1], s.t. p_{x'} >= m_s.
            l_x = 0;
            r_x = k-1;
            while (l_x != r_x) {
                m_x = l_x+(r_x-l_x)/2;
                if (I[m_x].first < m_s) {
                    l_x = m_x+1;
                } else {
                    r_x = m_x;
                }
            }

            // Find the minimum u' in [0,k-1], s.t. q_{pi[u']} >= m_s.
            l_u = 0;
            r_u = k-1;
            while (l_u != r_u) {
                m_u = l_u+(r_u-l_u)/2;
                if (I[pi[m_u]].second < m_s) {
                    l_u = m_u+1;
                } else {
                    r_u = m_u;
                }
            }

            /* If l_s = r_s, then l_s is an optimal value for s[i_p] and l_x
            and l_u are valid values for x[i_p] and u[i_p], respectively. */
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

        // Store l_s, l_x and l_u in s[i_p], x[i_p] and u[i_p], respectively.
        s[i_p] = l_s;
        x[i_p] = l_x;
        u[i_p] = l_u;
    }
}

template <typename uint_t>
void move_data_structure_phi<uint_t>::construction::calculate_seperation_positions_for_dq_and_mds() {
    x.resize(p+1);
    x[0] = 0;
    x[p] = k_;
    
    u.resize(p+1);
    u[0] = 0;
    u[p] = k_;

    // Compute s[1..p-1], x[1..p-1] and u[1..p-1].
    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // The optimal value i_p * lfloor (r'+r'')/p rfloor for s[i_p].
        uint_t o = i_p*((2*k_)/p);

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
            r_x = k_-1;
            while (l_x != r_x) {
                m_x = l_x+(r_x-l_x)/2;
                if (mds.p(m_x) < m_s) {
                    l_x = m_x+1;
                } else {
                    r_x = m_x;
                }
            }

            // Find the minimum u' in [0,r'-1], s.t. SA_s[pi'[u']] >= m_s.
            l_u = 0;
            r_u = k_-1;
            while (l_u != r_u) {
                m_u = l_u+(r_u-l_u)/2;
                if (D_q[pi[m_u]] < m_s) {
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
    }
}

/**
 * @brief builds D_offs and D_idx
 */
template <typename uint_t>
void move_data_structure_phi<uint_t>::construction::build_didx_doffs_v2v3v4v5() {
    if (log) log_message("building D_offs and D_idx");

    build_pi_for_dq();

    calculate_seperation_positions_for_dq_and_mds();
        
    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        /* Check if the range [u[i_p]..u[i_p+1]-1], over which the thread i_p
        has to iterate in D_q, is empty. */
        if (u[i_p+1] > u[i_p]) {
            // Iteration range start position in D_p.
            uint_t i = x[i_p];
            // Iteration range start position in D_q.
            uint_t j = u[i_p];

            // Iteration range end position in D_p.
            uint_t i_ = x[i_p+1]-1;
            // Iteration range end position in D_q.
            uint_t j_ = u[i_p+1]-1;

            /* Check if the first value in D_q lies before x[i_p]-th 
            input interval. */
            if (D_q[pi[j]] < mds.p(i)) {
                i--;
            }

            /* Iterate until one of the iteration end positions i_ and j_ has
            been reached. */
            while (i <= i_ && j <= j_) {

                /* Iterate over the values in D_q that lie in the current
                i-th input interval. */
                while (j <= j_ && D_q[pi[j]] < mds.p(i+1)) {
                    
                    /* Because each of those j-th largest value in D_q lie in the i-th
                    input interval, we can set D_idx[pi[j]] = i for each of them. */
                    mds.set_idx(pi[j],i);
                    mds.set_offs(pi[j],D_q[pi[j]]-mds.p(i));

                    j++;
                }

                i++;
            }
        }
    }

    x.clear();
    x.shrink_to_fit();

    u.clear();
    u.shrink_to_fit();

    if (log) {
        if (mf != NULL) *mf << " time_build_didx_doffs=" << time_diff_ns(time);
        time = log_runtime(time);
    }
}