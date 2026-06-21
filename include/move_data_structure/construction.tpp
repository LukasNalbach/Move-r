/**
 * part of LukasNalbach/Move-r
 *
 * MIT License
 *
 * Copyright (c) Lukas Nalbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#pragma once

#include <ips4o.hpp>
#include <move_data_structure/move_data_structure.hpp>

// ############################# COMMON METHODS #############################

template <typename pos_t>
void move_data_structure<pos_t>::construction::build_pi_for_I()
{
    no_init_resize(pi, k);

    // write the identity permutation of [0..k-1] to pi
    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < k; i++) {
        pi[i] = i;
    }

    // sort pi by the output interval starting positions in I
    auto comp_pi = [&](pos_t i, pos_t j) { return I[i].second < I[j].second; };
    if (p > 1) {
        ips4o::parallel::sort(pi.begin(), pi.end(), comp_pi);
    } else {
        ips4o::sort(pi.begin(), pi.end(), comp_pi);
    }
}

template <typename pos_t>
void move_data_structure<pos_t>::construction::build_pi_for_dq()
{
    no_init_resize(pi, k_ + 1);

    // write the identity permutation of [0..k'] to pi
    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i <= k_; i++) {
        pi[i] = i;
    }

    // sort pi by D_q
    auto comp_pi = [&](pos_t i, pos_t j) { return D_q[i] < D_q[j]; };
    if (p > 1) {
        ips4o::parallel::sort(pi.begin(), pi.end(), comp_pi);
    } else {
        ips4o::sort(pi.begin(), pi.end(), comp_pi);
    }
}

template <typename pos_t>
void move_data_structure<pos_t>::construction::calculate_seperation_positions_for_I()
{
    s.resize(p + 1);
    s[0] = 0;
    s[p] = n;

    x.resize(p + 1);
    x[0] = 0;
    x[p] = k;

    u.resize(p + 1);
    u[0] = 0;
    u[p] = k;

    // calculate seperation positions
    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // The optimal value i_p * lfloor 2k/p rfloor for s[i_p].
        pos_t o = i_p * ((2 * k) / p);

        // Search range of the binary search for s[i_p], initialized to [0..n-1].
        pos_t l_s = 0;
        pos_t r_s = n - 1;
        // The values x[i_p] and u[i_p] corresponding to the current s[i_p].
        pos_t l_x;
        pos_t l_u;

        // Perform a binary search over [0..n-1], to find s[i_p].
        while (true) {
            /* Set the candidate position for s[i_p] to the position
            in the middle between l_s and r_s. */
            pos_t m_s = l_s + (r_s - l_s) / 2;

            // Find the minimum x' in [0,k-1], s.t. p_{x'} >= m_s.
            l_x = bin_search_min_geq<pos_t>(m_s, 0, k - 1, [&](pos_t y) { return I[y].first; });

            // Find the minimum u' in [0,k-1], s.t. q_{pi[u']} >= m_s.
            l_u = bin_search_min_geq<pos_t>(m_s, 0, k - 1, [&](pos_t y) { return I[pi[y]].second; });

            /* If l_s = r_s, then l_s is an optimal value for s[i_p] and l_x
            and l_u are valid values for x[i_p] and u[i_p], respectively. */
            if (l_s == r_s) {
                break;
            }

            // Else, adjust the range for the binary search over [0..n-1].
            if (l_x + l_u < o) {
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

template <typename pos_t>
void move_data_structure<pos_t>::construction::calculate_seperation_positions_for_dq_and_mds()
{
    x.resize(p + 1);
    x[0] = 0;
    x[p] = k_;

    u.resize(p + 1);
    u[0] = 0;
    u[p] = k_;

    // Compute s[1..p-1], x[1..p-1] and u[1..p-1].
    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // The optimal value i_p * lfloor (r'+r'')/p rfloor for s[i_p].
        pos_t o = i_p * ((2 * k_) / p);

        // Search range of the binary search for s[i_p], initialized to [0..n-1].
        pos_t l_s = 0;
        pos_t r_s = n - 1;
        // The values x[i_p] and u[i_p] corresponding to the current s[i_p].
        pos_t l_x;
        pos_t l_u;

        // Perform a binary search over [0..n-1], to find s[i_p].
        while (true) {
            /* Set the candidate position for s[i_p] to the position in the middle
            between l_s and r_s. */
            pos_t m_s = l_s + (r_s - l_s) / 2;

            // Find the minimum x' in [0,r''-1], s.t. D_p[x'] >= m_s.
            l_x = bin_search_min_geq<pos_t>(m_s, 0, k_ - 1, [&](pos_t y) { return mds.p(y); });

            // Find the minimum u' in [0,r'-1], s.t. D_q[pi'[u']] >= m_s.
            l_u = bin_search_min_geq<pos_t>(m_s, 0, k_ - 1, [&](pos_t y) { return D_q[pi[y]]; });

            /* If l_s = r_s, then l_s is an optimal value for s[i_p] and l_x and
            l_u are valid values for x[i_p] and u[i_p], respectively. */
            if (l_s == r_s) {
                break;
            }

            // Else, adjust the range for the binary search over [0..n-1].
            if (l_x + l_u < o) {
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
template <typename pos_t>
void move_data_structure<pos_t>::construction::build_didx_doffs()
{
    log_message(log, "building D_offs and D_idx");

    build_pi_for_dq();

    calculate_seperation_positions_for_dq_and_mds();

    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        // Check if thread i_p's section D_q[u[i_p]..u[i_p+1]-1] is empty.
        if (u[i_p] < u[i_p + 1]) {
            // Iteration range start position in D_p.
            pos_t i = x[i_p];
            // Iteration range start position in D_q.
            pos_t j = u[i_p];
            // Iteration range end position in D_q + 1.
            pos_t j_ = u[i_p + 1];

            // Check if the first value D_q[pi[j]] lies before the x[i_p]-th input interval.
            while (D_q[pi[j]] < mds.p(i)) {
                i--;
            }

            // Check if the first value D_q[pi[j]] lies after the x[i_p]-th input interval.
            while (mds.p(i + 1) <= D_q[pi[j]]) {
                i++;
            }

            // Iterate until one of the iteration end positions i_ and j_ has been reached.
            while (j < j_) {

                // Iterate over the values in D_q that lie in the current i-th input interval.
                while (j < j_ && D_q[pi[j]] < mds.p(i + 1)) {

                    /* Because each of those j-th largest values in D_q lie in the i-th
                    input interval, we can set D_idx[pi[j]] = i for each of them. */
                    mds.set_idx(pi[j], i);
                    mds.set_offs(pi[j], D_q[pi[j]] - mds.p(i));

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

    log_phase_end(log, time, mf, "time_build_didx_doffs");
}
// ############################# INTERVAL SEQUENCE CONSTRUCTION #############################

template <typename pos_t>
void move_data_structure<pos_t>::construction::build_tin_tout()
{
    log_message(log, "building pi");

    // build pi to construct T_out faster
    build_pi_for_I();
    calculate_seperation_positions_for_I();

    log_phase_end(log, time, mf, "time_build_pi");
    log_message(log, "building T_out");

    T_out.resize(p);

    // build T_out[0..p-1]
    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        pos_t b = u[i_p];
        pos_t e = u[i_p + 1];

        for (pos_t i = b; i < e; i++) {
            T_out[i_p].emplace_hint(T_out[i_p].end(), I[pi[i]]);
        }
    }

    u.clear();
    u.shrink_to_fit();

    pi.clear();
    pi.shrink_to_fit();

    log_phase_end(log, time, mf, "time_build_tout");
    log_message(log, "building T_in");

    T_in.resize(p);

    // build T_in[0..p-1]
    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        if (x[i_p] < x[i_p + 1]) {
            T_in[i_p].insert(I.begin() + x[i_p], I.begin() + x[i_p + 1]);
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
    for (uint16_t i = 1; i < p; i++) {
        if (T_in[i].empty() || (*T_in[i].begin()).first != s[i]) {
            pair_t pr_split = *T_in[i - 1].rbegin();
            pair_t pr_new { s[i], pr_split.second + (s[i] - pr_split.first) };
            uint16_t i_p_ = bin_search_max_leq<pos_t>(pr_new.second, 0, p - 1, [&](pos_t x) { return s[x]; });
            T_in[i].emplace(pr_new);
            T_out[i_p_].emplace(pr_new);
        }
    }

    // make sure there is an output interval starting at s[0], s[1], ... s[p-1]
    for (uint16_t i = 1; i < p; i++) {
        if (T_out[i].empty() || (*T_out[i].begin()).second != s[i]) {
            pair_t pr_split = *T_out[i - 1].rbegin();
            pair_t pr_new { pr_split.first + (s[i] - pr_split.second), s[i] };
            uint16_t i_p_ = bin_search_max_leq<pos_t>(pr_new.first, 0, p - 1, [&](pos_t x) { return s[x]; });
            T_in[i_p_].emplace(pr_new);
            T_out[i].emplace(pr_new);
        }
    }

    // add dummy pairs (s[i+1],s[i+1]) to the i-th trees, for i in [0..p-1]
    // this enables us to calculate the length of the last input- and output intervals in each section
    for (uint16_t i = 0; i < p; i++) {
        T_in[i].emplace(pair_t { s[i + 1], s[i + 1] });
        T_out[i].emplace(pair_t { s[i + 1], s[i + 1] });
    }

    log_phase_end(log, time, mf, "time_build_tin");
    log_message(log, "splitting too long intervals");

    T_out_temp.resize(p, std::vector<tout_t>(p));

    // iterate over all trees in T_in and split every input interval that is longer than l_max
    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        tin_it_t tn_I = T_in[i_p].begin();
        pair_t pr_Im1 = *tn_I;
        tn_I++;

        while (tn_I != T_in[i_p].end()) {
            while (((*tn_I).first - pr_Im1.first) > l_max) {
                tn_I = T_in[i_p].emplace_hint(tn_I, pair_t { pr_Im1.first + l_max, pr_Im1.second + l_max });
                pr_Im1 = *tn_I;
                tn_I++;

                uint16_t i_p_ = bin_search_max_leq<pos_t>(pr_Im1.second, 0, p - 1, [&](pos_t x) { return s[x]; });

                // store each new pair in its corresponding tree in T_out_temp
                T_out_temp[i_p_][i_p].emplace(pr_Im1);
            }

            pr_Im1 = *tn_I;
            tn_I++;
        }
    }

    // merge T_out with T_out_temp
    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        for (uint16_t i = 0; i < p; i++) {
            T_out[i_p].merge(T_out_temp[i_p][i]);
        }
    }

    T_out_temp.clear();
    T_out_temp.shrink_to_fit();

    log_phase_end(log, time, mf, "time_split_too_long_input_intervals");
}

template <typename pos_t>
void move_data_structure<pos_t>::construction::build_dp_dq()
{
    // recalculate the indices x of the seperation positions in the interval sequence
    x.resize(p + 1);
    x[0] = 0;

    for (uint16_t i = 0; i < p; i++) {
        T_in[i].erase(*T_in[i].rbegin());
        x[i + 1] = x[i] + T_in[i].size();
    }

    k_ = x[p];

    // resize the interleaved vectors in the move data structure
    mds.resize(n, k_, width_l_);

    if (log) {
        float k__k = std::round(100.0 * k_ / k) / 100.0;
        if (mf != nullptr) {
            *mf << " k=" << k;
            *mf << " k_=" << k_;
        }
        std::cout << "k' = " << k_ << ", k'/k = " << k__k << std::endl;
        log_message("building D_p and D_q");
    }

    // D_q stays byte-aligned (transient build buffer); its values are output positions in [0,n]
    D_q = interleaved_byte_aligned_vectors<pos_t, pos_t>({ byte_width(n) });
    D_q.resize_no_init(k_ + 1);
    D_q.template set_parallel<0, pos_t>(k_, n);

    // write the input interval starting positions to D_p (in the move data structure) and
    // write the output interval starting positions to D_q
    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        pos_t b = x[i_p];
        pos_t e = x[i_p + 1];

        tin_it_t tin_it = T_in[i_p].begin();

        for (pos_t i = b; i < e; i++) {
            mds.set_p(i, (*tin_it).first);
            D_q.template set_parallel<0, pos_t>(i, (*tin_it).second);
            tin_it++;
        }
    }

    x.clear();
    x.shrink_to_fit();

    T_in.clear();
    T_in.shrink_to_fit();

    log_phase_end(log, time, mf, "time_build_dp_dq");
}
// ############################# BALANCING #############################

template <typename pos_t>
inline pos_t move_data_structure<pos_t>::construction::is_a_heavy(tin_it_t& tn_I, pos_t q_J_)
{
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
        pos_t qj_pd = (*tn_I).first;

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

template <typename pos_t>
inline typename move_data_structure<pos_t>::construction::tout_it_t move_data_structure<pos_t>::construction::balance_upto(tout_it_t& tn_J_, pos_t qj_pd, pos_t q_u)
{
    // Index in [0..p-1] of the current thread.
    uint16_t i_p = omp_get_thread_num();

    pos_t q_J_ = (*tn_J_).second; // q_j'
    tout_it_t tn_J = tn_J_; // iterator pointing to the pair (p_j',q_j') in T_out[i_p]
    tn_J--;
    pos_t q_j = (*tn_J).second; // iterator pointing to the pair (p_j,q_j) in T_out[i_p]
    pos_t pj_pd = (*tn_J).first + (qj_pd - q_j); // p_j + d
    pair_t pr_new { pj_pd, qj_pd }; // the newly created pair

    // insert the newly created pair into the current thread's tree in T_out
    tout_it_t tout_n_new = T_out[i_p].emplace_hint(tn_J_, pr_new);

    // check, if the newly created has to be inserted into a tree in T_in of another thread
    if (p != 1 && !(s[i_p] <= pj_pd && pj_pd < s[i_p + 1])) {
        // calculate the index i_p' of the thread, into whiches tree in T_in, the newly created pair has to be inserted
        uint16_t i_p_ = bin_search_max_leq<pos_t>(pj_pd, 0, p - 1, [&](pos_t x) { return s[x]; });

        // store the newly created pair in Q[i_p'][i_p]
        Q[i_p_][i_p].emplace_back(pr_new);
    } else {
        // else the newly created pair can be inserted into the current thread's tree in T_in
        tin_it_t tin_n_new = node(T_in[i_p].emplace(pr_new));

        // check whether the output interval containing p_j + d can have become unbalanced by splitting [p_j, p_j + d_j)
        if (pj_pd < q_u && (pj_pd < q_j || q_J_ <= pj_pd)) {
            // if yes, find the node in T_out[i_p] creating the pair (p_y,q_y), where p_j + d in [q_j, q_y + d_y)
            tout_it_t tn_Y = T_out[i_p].lower_bound(pair_t { 0, pj_pd });

            if (pj_pd < (*tn_Y).second) {
                tn_Y--;
            }

            // find the node in T_in[i_p] creating the first input interval starting in [q_j, q_y + d_y)
            while ((*tn_Y).second < (*tin_n_new).first) {
                tin_n_new--;
            }

            if ((*tin_n_new).first < (*tn_Y).second) {
                tin_n_new++;
            }

            // iterate one step with tn_Y, s.t. it now points to the output interval starting direclty after [q_y, q_y + d_y)
            tn_Y++;
            pos_t qy_pd_; // q_y + d', where d' = p_{z+a} - q_y

            // check if [q_y, q_y + d_y) is a-heavy and balanced
            if ((qy_pd_ = is_a_heavy(tin_n_new, (*tn_Y).second))) {
                // if yes, balance it and all a-heavy output intervals in [s[i_p],s[i_p+1]) starting before q_u that
                // become a-heavy in the process
                balance_upto(tn_Y, qy_pd_, q_u);

                // because we inserted another pair into T_out[i_p] in the recursive call of balance_upto, tout_n_new
                // may not point to (p_j + d, q_j + d) anymore, so return T_out[i_p].end() (which is constant)
                return T_out[i_p].end();
            }
        }
    }

    // return an iterator pointing to the newly created pair (p_j + d, q_j + d) in T_out, since no recursive call has been made
    return tout_n_new;
}

template <typename pos_t>
void move_data_structure<pos_t>::construction::balance()
{
    log_message(log, p > 1 ? "balancing (phase 1)" : "balancing");

    if (p > 1) {
        Q.resize(p, std::vector<pair_arr_t>(p));
        Q_.resize(p, std::vector<pair_arr_t>(p));

        for (uint16_t i = 0; i < p; i++) {
            for (uint16_t j = 0; j < p; j++) {
                Q[i][j].reserve(k / (32 * a * p * p));
                Q_[i][j].reserve(k / (32 * a * p * p));
            }
        }
    }

    // first phase of the balancing algorithm
    #pragma omp parallel num_threads(p)
    {
        // Index in [0..p-1] of the current thread.
        uint16_t i_p = omp_get_thread_num();

        tin_it_t tn_I = T_in[i_p].begin(); // iterator pointing to the pair in T_in[i_p] creating (p_i,q_i)
        tout_it_t tn_J = T_out[i_p].begin(); // iterator pointing to the pair in T_out[i_p] creating (p_j,q_j)
        /* iterator pointing to the pair in T_out[i_p] creating (p_j',q_j'), where [q_j', q_j' + d_j') is the output interval,
           that starts direclty after [q_j, q_j + d_j) */
        tout_it_t tn_J_ = tn_J;
        tn_J_++;
        pos_t qj_pd; // q_j + d (temporary variable)
        bool stop = false;

        while (!stop) {
            // check if the current output interval [q_j, q_j + d_j) is a-heavy
            if ((qj_pd = is_a_heavy(tn_I, (*tn_J_).second))) {
                // if yes, balance it and all a-heavy output intervals in [s[i_p],s[i_p+1]) starting before q_j + d that
                // become a-heavy in the process
                tn_J = balance_upto(tn_J_, qj_pd, qj_pd);

                // because we inserted a pair into T_in[i_p], the iterator tn_I may now be invalid, so reset it
                tn_I = T_in[i_p].find(pair_t { qj_pd, 0 });

                // if tn_J points to T_out[i_p].end(), balance_upto made a recursive call, so reset tn_J
                if (tn_J == T_out[i_p].end()) {
                    tn_J = T_out[i_p].find(pair_t { 0, qj_pd });
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

                if (tn_J_ == T_out[i_p].end()) {
                    stop = true;
                    break;
                }

                while (!stop && (*tn_I).first < (*tn_J).second) {
                    tn_I++;

                    if (tn_I == T_in[i_p].end()) {
                        stop = true;
                        break;
                    }
                }
            } while (!stop && (*tn_I).first >= (*tn_J_).second);
        }
    }

    if (log) {
        if (mf != nullptr) {
            *mf << " time_balance_phase_1=" << time_diff_ns(time);
            if (p == 1)
                *mf << " time_balance_phase_2=0" << time_diff_ns(time);
        }
        time = log_runtime(time);
    }

    // second phase of the balancing algorithm
    if (p > 1) {
        log_message(log, "balancing (phase 2)");

        bool done;

        #pragma omp parallel num_threads(p)
        {
            // Index in [0..p-1] of the current thread.
            uint16_t i_p = omp_get_thread_num();

            pos_t qy_pd; // q_y + d
            tin_it_t tn_new = T_in[i_p].end(); // iterator pointing to the newly created pair in T_out[i_p]
            tout_it_t tn_Y = T_out[i_p].end(); // iterator pointing to the pair (p_y, q_y) in T_out[i_p]

            while (true) {
                #pragma omp barrier

                #pragma omp single
                {
                    // the pairs that have been inserted into Q in the first phase (last iteration of the second phase)
                    // now have to be inserted into T_in[0..p-1], so swap Q with Q_
                    std::swap(Q, Q_);
                    done = true;
                }

                #pragma omp barrier

                // check whether Q_ is empty
                for (uint16_t i_p_ = 0; i_p_ < p; i_p_++) {
                    if (!Q_[i_p][i_p_].empty()) {
                        done = false;
                        break;
                    }
                }

                #pragma omp barrier

                // if Q_ is empty, there are no a-heavy output intervals, so break
                if (done) {
                    break;
                }

                // iterate over all pairs to insert into T_in[i_p]
                for (pair_arr_t& vec : Q_[i_p]) {
                    for (pair_t& pr_new : vec) { // a newly created pair (p_j + d, q_j + d)
                        // insert (p_j + d, q_j + d) into T_in[i_p]
                        tn_new = node(T_in[i_p].emplace(pr_new));

                        // find the pair in T_out[i_p] containing (p_j + d, q_j + d)
                        tn_Y = T_out[i_p].lower_bound(pair_t { 0, pr_new.first });

                        if ((*tn_new).first < (*tn_Y).second) {
                            tn_Y--;
                        }

                        // find the node in T_in[i_p] creating the first input interval starting in [q_j, q_y + d_y)
                        while ((*tn_Y).second < (*tn_new).first) {
                            tn_new--;
                        }

                        if ((*tn_new).first < (*tn_Y).second) {
                            tn_new++;
                        }

                        // iterate one step with tn_Y, s.t. it now points to the output interval starting direclty after [q_y, q_y + d_y)
                        tn_Y++;

                        // check if [q_y, q_y + d_y) is a-heavy
                        if ((qy_pd = is_a_heavy(tn_new, (*tn_Y).second))) {
                            // if yes, balance it and all a-heavy output intervals in [s[i_p],s[i_p+1]) becoming a-heavy in the process
                            balance_upto(tn_Y, qy_pd, n);
                        }
                    }

                    vec.clear();
                }
            }
        }

        Q.clear();
        Q.shrink_to_fit();

        Q_.clear();
        Q_.shrink_to_fit();

        log_phase_end(log, time, mf, "time_balance_phase_2");
    }

    T_out.clear();
    T_out.shrink_to_fit();

    s.clear();
    s.shrink_to_fit();
}
// ############################# CORRECTNESS VERIFICATION #############################

template <typename pos_t>
void move_data_structure<pos_t>::construction::verify_correctness()
{
    std::cout << "verifying correctness of the interval sequence:" << std::endl;
    bool correct = true;

    // check if the input interval starting positions ascend
    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < k_; i++) {
        if (!(mds.p(i) < mds.p(i + 1))) {
            #pragma omp critical
            {
                std::cout << "input interval starting positions do not ascend:" << std::endl;
                std::cout << "i = " << i << std::endl;
                std::cout << "p_i = " << mds.p(i) << std::endl;
                std::cout << "p_{i+1} = " << mds.p(i + 1) << std::endl
                          << std::endl;
                correct = false;
            }
        }
    }

    // check if an input interval is too long
    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < k_; i++) {
        if (mds.p(i + 1) - mds.p(i) > l_max) {
            #pragma omp critical
            {
                std::cout << "too long input interval (> l_max = " << l_max << "):" << std::endl;
                std::cout << "i = " << i << std::endl;
                std::cout << "p_i = " << mds.p(i) << std::endl;
                std::cout << "p_{i+1} = " << mds.p(i + 1) << std::endl
                          << std::endl;
                correct = false;
            }
        }
    }

    // check if the output interval lengths do not match the input interval lengths
    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < k_; i++) {
        if (D_q[pi[i + 1]] - D_q[pi[i]] != mds.p(pi[i] + 1) - mds.p(pi[i])) {
            #pragma omp critical
            {
                std::cout << "input- and output interval lengths do not match:" << std::endl;
                std::cout << "i = " << i << std::endl;
                std::cout << "p_{pi[i]} = " << mds.p(pi[i]) << std::endl;
                std::cout << "p_{pi[i]+1} = " << mds.p(pi[i] + 1) << std::endl;
                std::cout << "q_{pi[i]} = " << D_q[pi[i]] << std::endl;
                std::cout << "q_{pi[i+1]} = " << D_q[pi[i + 1]] << std::endl
                          << std::endl;
                correct = false;
            }
        }
    }

    pos_t i = 0;
    pos_t j = 0;
    pos_t e;

    // check if there is an a-heavy output interval by iterating over the input- and output intervals with pi
    while (i < k && j < k) {
        while (i < k && mds.p(i) < D_q[pi[j]]) {
            i++;
        }

        if (mds.p(i) < D_q[pi[j + 1]]) {
            e = 1;

            while (i + 1 < k && mds.p(i + 1) < D_q[pi[j + 1]]) {
                i++;
                e++;

                if (e >= two_a) {
                    std::cout << "a-heavy output interval:" << std::endl;
                    std::cout << "i = " << i - two_a + 1 << std::endl;
                    std::cout << "j = " << j << std::endl;
                    std::cout << "q_{pi[j]} = " << D_q[pi[j]] << std::endl;
                    std::cout << "p_{i} = " << mds.p(i - two_a + 1) << std::endl;
                    std::cout << "p_{i+2a} = " << mds.p(i) << std::endl;
                    std::cout << "q_{pi[j+1]} = " << D_q[pi[j + 1]] << std::endl
                              << std::endl;
                    correct = false;
                    j++;
                    break;
                }
            }

            i++;
        } else {
            j++;
        }
    }

    // check whether D_idx has been calculated correctly
    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < k_; i++) {
        if (!(mds.p(mds.idx(i)) <= D_q[i] && D_q[i] < mds.p(mds.idx(i) + 1))) {
            #pragma omp critical
            {
                std::cout << "wrong value:" << std::endl;
                std::cout << "i = " << i << std::endl;
                std::cout << "D_p[D_idx[i]] = " << mds.p(mds.idx(i)) << std::endl;
                std::cout << "D_q(i) = " << D_q[i] << std::endl;
                std::cout << "D_p[D_idx[i]+1] = " << mds.p(mds.idx(i) + 1) << std::endl
                          << std::endl;
                correct = false;
            }
        }
    }

    // check if D_offs has been calculated correctly (by checking if we can recalculate D_q together with D_idx and D_p)
    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < k_; i++) {
        if (mds.q(i) != D_q[i]) {
            #pragma omp critical
            {
                std::cout << "wrong value:" << std::endl;
                std::cout << "i = " << i << std::endl;
                std::cout << "D_idx[i] = " << mds.idx(i) << std::endl;
                std::cout << "D_p[D_idx[i]] = " << mds.p(mds.idx(i)) << std::endl;
                std::cout << "D_offs[i] = " << mds.offs(i) << std::endl;
                std::cout << "mds.q(i) = D_p[D_idx[i]] + D_offs[i] = " << mds.q(i);
                std::cout << " != D_q[i] = " << D_q[i] << std::endl
                          << std::endl;
                correct = false;
            }
        }
    }

    std::cout << "the interval sequence has " << (correct ? "" : "not ") << "been built correctly" << std::endl;
}