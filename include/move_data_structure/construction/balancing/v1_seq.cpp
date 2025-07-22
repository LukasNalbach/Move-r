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
void move_data_structure<pos_t>::construction::balance_v1_seq()
{
    if (log) log_message("building T_e");

    // build T_e_v1
    pos_t q_i, q_next;
    tin_node_t_v1* node_cur = T_in_v1.min();
    tin_node_t_v1* node_cur_2;
    uint32_t e;
    
    while (node_cur != T_in_v1.max()) {
        /* For each output interval [q_i, q_i + d_i), find the first input interval connected
        to it in the permutation graph. */
        q_i = node_cur->v.second;
        q_next = q_i + node_cur->nxt()->v.first - node_cur->v.first;
        node_cur_2 = T_in_v1.min_geq(pair_t { q_i, 0 });

        // Count the number of input intervals connected to it in the permutation graph.
        e = 0;

        while (node_cur_2 != nullptr) {
            if (node_cur_2->v.first < q_next) {
                e++;

                if (e == two_a) {
                    // If there are at least 2a, insert it's corresponding pair into T_e_v1.
                    T_e_v1.insert(node_cur->v);
                    break;
                }
            } else {
                break;
            }

            node_cur_2 = node_cur_2->nxt();
        }

        node_cur = node_cur->nxt();
    }

    if (log) {
        if (mf != nullptr)
            *mf << " time_build_te=" << time_diff_ns(time);
        time = log_runtime(time);
        log_message("balancing");
    }

    // balance the disjoint interval sequence
    pos_t d, q_j, p_j, q_y, d_j, d_y;
    pair_t pair_NEW, pair_Y;
    tin_node_t_v1* node_Ipa;
    tout_te_node_t_v1 *node_NEW, *node_Y, *min;
    std::vector<std::tuple<pos_t, pos_t, pair_t>> intervals_to_check;
    while (!T_e_v1.empty()) {
        /* Find the pair creating the first a-heavy output interval [q_j, q_j + d_j)
        and remove it from T_e_v1. */
        min = T_e_v1.min();
        p_j = min->v.first;
        q_j = min->v.second;
        delete T_e_v1.remove(min->v);

        // Find the a+1-st input interval in [q_j, q_j + d_j) and set d = p_{i+a} - q_j.
        /* d is the smallest integer, so that [q_j, q_j + d) has a incoming edges in the
        permutation graph. */
        node_Ipa = T_in_v1.min_geq(pair_t { q_j, 0 });

        for (uint16_t i = 0; i < a; i++) {
            node_Ipa = node_Ipa->nxt();
        }

        d = node_Ipa->v.first - q_j;

        // Create the new pair (p_j + d, q_j + d) and insert it into T_in_v1 and T_out_v1.
        pair_NEW = pair_t { p_j + d, q_j + d };
        T_in_v1.insert(pair_NEW);
        node_NEW = T_out_v1.insert(pair_NEW);

        /* Find the output interval [q_y, q_y + d_y), [p_j + d, p_j + d_j) is connected
        to in the permutation graph. */
        node_Y = T_out_v1.max_leq(pair_t { 0, p_j + d });
        pair_Y = node_Y->v;
        q_y = pair_Y.second;

        /* The number of input intervals connected to [q_j + d, q_j + d_j) and [q_y, q_y + d_y)
        may have changed. For each, check if it has at least 2a incoming edges in the permutation
        graph and insert it into T_e_v1, if it has. */
        d_j = node_NEW->nxt()->v.second - q_j;
        d_y = node_Y->nxt()->v.second - q_y;

        intervals_to_check = {
            std::make_tuple(q_j + d, q_j + d_j - 1, pair_NEW),
            std::make_tuple(q_y, q_y + d_y - 1, pair_Y)
        };

        for (auto tup : intervals_to_check) {
            /* Find the first input interval connected to the output interval in the
            permutation graph. */
            node_cur = T_in_v1.min_geq(pair_t { std::get<0>(tup), 0 });

            // Count the number of input intervals connected to it in the permutation graph.
            e = 0;
            while (node_cur != nullptr) {
                if (node_cur->v.first <= std::get<1>(tup)) {
                    e++;

                    // If there are at least 2a, insert it's corresponding pair into T_e_v1.
                    if (e == two_a) {
                        T_e_v1.insert(std::get<2>(tup));
                        break;
                    }
                } else {
                    break;
                }

                node_cur = node_cur->nxt();
            }
        }
    }
    // Because T_e_v1 is empty, there are no a-heavy output intervals.

    k_ = T_in_v1.size() - 1;
    T_out_v1.delete_nodes();

    if (log) {
        time = log_runtime(time);
        float k__k = std::round(100.0 * k_ / k) / 100.0;

        if (mf != nullptr) {
            *mf << " k=" << k;
            *mf << " k_=" << k_;
            *mf << " time_balance=" << time_diff_ns(time);
        }

        std::cout << "k' = " << k_ << ", k'/k = " << k__k << std::endl;
    }
}