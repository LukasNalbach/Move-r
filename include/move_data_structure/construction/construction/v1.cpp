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

#include <move_data_structure/move_data_structure.hpp>

template <typename pos_t>
void move_data_structure<pos_t>::construction::build_tin_tout_v1()
{
    if (log) log_message("building T_in");

    // build T_in_v1 and T_out_v1
    for (pos_t i = 0; i < k; i++) {
        T_in_v1.insert(I[i]);
    }

    T_in_v1.insert(pair_t { n, n });

    if (log) {
        if (mf != NULL)
            *mf << " time_build_tin=" << time_diff_ns(time);
        time = log_runtime(time);
        log_message("building T_out");
    }

    for (pos_t i = 0; i < k; i++) {
        T_out_v1.insert(I[i]);
    }

    T_out_v1.insert(pair_t { n, n });
    tin_node_t_v1* node_cur = T_in_v1.min();

    while (node_cur != T_in_v1.max()) {
        while (node_cur->nxt()->v.first - node_cur->v.first > l_max) {
            T_out_v1.insert(pair_t { node_cur->v.first + l_max, node_cur->v.second + l_max });
            node_cur = T_in_v1.insert(pair_t { node_cur->v.first + l_max, node_cur->v.second + l_max });
        }
        
        node_cur = node_cur->nxt();
    }

    if (delete_i) {
        // Now, we do not need I anymore.
        I.clear();
        I.shrink_to_fit();
    }

    if (log) {
        if (mf != NULL)
            *mf << " time_build_tout=" << time_diff_ns(time);
        time = log_runtime(time);
    }
}

template <typename pos_t>
void move_data_structure<pos_t>::construction::build_dp_dq_v1()
{
    if (log)
        log_message("building D_p and D_q");

    mds.resize(n, k_, width_l_);
    D_q = interleaved_byte_aligned_vectors<pos_t, pos_t>({ (uint8_t)(mds.omega_p / 8) });
    D_q.resize_no_init(k_);
    D_q.template set_parallel<0, pos_t>(k_, n);

    auto it = T_in_v1.iterator();

    for (pos_t i = 0; i <= k_; i++) {
        mds.set_p(i, it.current()->v.first);
        D_q.template set_parallel<0, pos_t>(i, it.current()->v.second);
        it.next();
    }

    T_in_v1.delete_nodes();

    if (log) {
        if (mf != NULL)
            *mf << " time_write_dp_dq=" << time_diff_ns(time);
        time = log_runtime(time);
    }
}

template <typename pos_t>
void move_data_structure<pos_t>::construction::build_didx_doffs_v1()
{
    if (log)
        log_message("building D_idx and D_offs");

    for (pos_t j = 0; j < k_; j++) {
        /* For each output interval [q_j, q_j + d_j), find the input interval [p_i, p_i + d_i)
        containing q_j and set D_idx[j] = i. Find the maximum integer i in [0,k'-1], s.t.
        p_i <= q_j with a binary search over D_pair. */
        pos_t i = bin_search_max_leq<pos_t>(D_q[j], 0, k_ - 1, [&](pos_t x) { return mds.p(x); });

        mds.set_idx(j, i);
        mds.set_offs(j, D_q[j] - mds.p(i));
    }

    if (log) {
        if (mf != NULL)
            *mf << " time_build_didx_doffs=" << time_diff_ns(time);
        time = log_runtime(time);
    }
}