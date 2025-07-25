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
void move_data_structure<pos_t>::construction::balance_v2_seq()
{
    if (log) log_message("building T_e");

    std::vector<te_node_t_v2> nodes_te = std::vector<te_node_t_v2>();
    nodes_te.reserve(k / (two_a));

    // points to to the pair (p_i,q_i).
    lin_node_t_v2v3v4* ln_I = L_in_v2v3v4[0].head();
    // points to the pair (p_j,q_j).
    typename tout_t_v2v3v4::avl_it it_J = T_out_v2v3v4[0].iterator();
    // points to the pair (p_{j'},q_{j'}), where q_j + d_j = q_{j'}.
    typename tout_t_v2v3v4::avl_it it_Jp1 = T_out_v2v3v4[0].iterator(T_out_v2v3v4[0].second_smallest());

    // temporary variables
    lin_node_t_v2v3v4* ln_IpA;
    pos_t i_ = 1;

    /* At the start of each iteration, [p_i, p_i + d_i) is the first input interval connected
    to [q_j, q_j + d_j) in the permutation graph */
    bool stop = false;
    while (!stop) {
        ln_IpA = is_a_heavy_v2v3v4(&ln_I, &i_, it_J.current(), it_Jp1.current());
        i_ = 1;

        /* If [q_j, q_j + d_j) is a-heavy, balance it and all output intervals starting before
        it that might get a-heavy in the process. */
        if (ln_IpA != nullptr) {
            nodes_te.push_back(te_node_t_v2(te_pair_t_v2 { ln_IpA, it_J.current() }));
        }

        /* Find the next output interval with an incoming edge in the permutation graph and the
        first input interval connected to it. */
        do {
            if (!it_Jp1.has_next()) {
                stop = true;
                break;
            }
            
            it_J.set(it_Jp1.current());
            it_Jp1.next();

            while (ln_I->v.first < it_J.current()->v.v.second) {
                if (ln_I->sc == nullptr) {
                    stop = true;
                    break;
                }

                ln_I = ln_I->sc;
            }
        } while (!stop && ln_I->v.first >= it_Jp1.current()->v.v.second);
    }

    // Build T_e_v2 from nodes_te.
    if (!nodes_te.empty()) {
        T_e_v2.insert_array((pos_t)0, (pos_t)nodes_te.size() - 1, [&](pos_t i) { return &nodes_te[i]; });
    }

    if (log) {
        if (mf != nullptr)
            *mf << " time_build_te=" << time_diff_ns(time);
        time = log_runtime(time);
        log_message("balancing");
    }

    // temporary variables
    pos_t p_j, q_j, d_j, d, q_y;
    tout_node_t_v2v3v4 *tn_J, *tn_NEW, *tn_Y;
    lin_node_t_v2v3v4 *ln_Ip2A, *ln_Z, *ln_ZpA;

    while (!T_e_v2.empty()) {
        /* Find the first a-heavy output interval and the a+1-st input interval connected to
        it in the permutation graph. */
        te_node_t_v2* min = T_e_v2.min();
        ln_IpA = min->v.first;
        tn_J = min->v.second;

        p_j = tn_J->v.v.first;
        q_j = tn_J->v.v.second;
        d_j = interval_length_v2v3_seq(&tn_J->v);

        // d is the smalles integer, so that [q_j, q_j + d) has 2 incoming edges in the permutation graph.
        d = ln_IpA->v.first - q_j;

        /* Create the pair (p_j + d, q_j + d), which creates two new input intervals [p_j, p_j + d) and
        [p_j + d, p_j + d_j). */
        tn_NEW = new_nodes_2v3v4[0].emplace_back(tout_node_t_v2v3v4(lin_node_t_v2v3v4(pair_t { p_j + d, q_j + d })));
        L_in_v2v3v4[0].insert_after_node(&tn_NEW->v, &tn_J->v);
        T_out_v2v3v4[0].insert_hint(tn_NEW, tn_J);

        // Check if [q_j + d, q_j + d_j) is a-heavy.
        i_ = 1;
        ln_Ip2A = is_a_heavy_v2v3v4(&ln_IpA, &i_, tn_NEW);

        ln_ZpA = nullptr;
        /* If p_j + d \in [q_j, q_j + d) or [q_j + d, q_j + d_j), [q_j + d, q_j + d_j) is the only
        possibly new a-heavy output interval. */
        if (p_j + d < q_j || q_j + d_j <= p_j + d) {
            /* Else find the output interval [q_y, q_y + d_y), to which [p_j + d, p_j + d_j) is connected
            in the permutation graph. */
            tn_Y = T_out_v2v3v4[0].max_leq(lin_node_t_v2v3v4(pair_t { 0, p_j + d }));
            q_y = tn_Y->v.v.second;

            /*
                Find the first input interval [p_z, p_z + d_z) connected to [q_y, q_y + d_y) in the permutation graph.
                Case 1: [p_z, p_z + d_z) ends before q_j.
                    Then [q_y, q_y + d_y) does as well and because [q_y, q_y + d_y) was not a-heavy before [q_j + d, q_j + d_j) has been
                    created, [p_z, p_z + d_z) can be reached by iterating backwards in L_in_v2v3v4 at most 2a-1 steps starting from (p_j + d, q_j + d).
                Case 2: [p_z, p_z + d_z) starts after q_j + d_j.
                    Case 2.1: [p_z, p_z + d_z) is found by iterating backwards at most 2a-1 steps in L_in_v2v3v4 starting from (p_j + d, q_j + d).
                        Then check if [q_y, q_y + d_y) is a-heavy and insert ((p_{z+a},q_{z+a}),(p_y,q_y)) into T_e if it is.
                    Case 2.2: [p_z, p_z + d_z) is not found by iterating backwards at most 2a-1 steps in L_in_v2v3v4 starting from (p_j + d, q_j + d).
                        Then there are at least 2a input intervals connected to [q_y, q_y + d_y) in the permutation graph, hence there exists a pair
                        ((p_{x+a},q_{x+a}),(p_y,q_y)) in T_e, where [p_{x+a}, p_{x+a} + d_{x+a}) was the first input interval connected to
                        [q_y, q_y + d_y) in the permutation graph before [q_j + d, q_j + d_j) has been created. It still is, because if it
                        was not, [p_z, p_z + d_z) would have been found in a < 2a-1 steps, hence ((p_{x+a},q_{x+a}),(p_y,q_y)) still is a valid pair in T_e.
            */

            // find (p_z,q_z)
            ln_Z = &tn_NEW->v;
            i_ = 1;
            while (i_ < two_a && ln_Z->pr != nullptr && ln_Z->pr->v.first >= q_y) {
                ln_Z = ln_Z->pr;
                i_++;
            }

            // check if case 2.1 holds
            if (p_j + d < q_j || ln_Z->pr == nullptr || ln_Z->pr->v.first < q_y) {
                ln_Z = &tn_NEW->v;
                pos_t i__ = i_;

                // check if [q_y, q_y + d_y) is a-heavy
                ln_ZpA = is_a_heavy_v2v3v4(&ln_Z, &i_, tn_Y);

                if (ln_ZpA != nullptr && i__ > a + 1 && ln_Z->sc != nullptr &&
                    ln_Z->sc->v.first < q_y + interval_length_v2v3_seq(&tn_Y->v)
                ) {
                    ln_ZpA = nullptr;
                }
            }
        }

        if (ln_ZpA != nullptr) {
            if (ln_Ip2A != nullptr) {
                // [q_j + d, q_j + d_j) and [q_y, q_y + d_y) are both new a-heavy output intervals
                if (ln_ZpA->v.first < ln_Ip2A->v.first) {
                    // and [q_y, q_y + d_y) ends before [q_j + d, q_j + d_j)
                    min->v = te_pair_t_v2 { ln_Ip2A, tn_NEW };
                    T_e_v2.emplace_hint(te_pair_t_v2 { ln_ZpA, tn_Y }, min);
                } else {
                    // and [q_y, q_y + d_y) starts after [q_j + d, q_j + d_j)
                    if (T_e_v2.size() == 1 || tn_Y->v.v.second < T_e_v2.second_smallest()->v.second->v.v.second) {
                        // and is the second a-heavy output interval
                        min->v = te_pair_t_v2 { ln_ZpA, tn_Y };
                        T_e_v2.emplace_hint(te_pair_t_v2 { ln_Ip2A, tn_NEW }, min);
                    } else {
                        // and is not the second a-heavy output interval
                        min->v = te_pair_t_v2 { ln_Ip2A, tn_NEW };
                        T_e_v2.insert(te_pair_t_v2 { ln_ZpA, tn_Y });
                    }
                }
            } else {
                // [q_y, q_y + d_y) is the only new a-heavy output interval
                if (T_e_v2.size() == 1 || tn_Y->v.v.second < T_e_v2.second_smallest()->v.second->v.v.second) {
                    // and is the first a-heavy output interval
                    min->v = te_pair_t_v2 { ln_ZpA, tn_Y };
                } else {
                    // and is not the first a-heavy output interval
                    T_e_v2.remove_node(min);

                    if (min < &nodes_te.front() || &nodes_te.back() < min) {
                        delete min;
                    }

                    T_e_v2.insert(te_pair_t_v2 { ln_ZpA, tn_Y });
                }
            }
        } else {
            if (ln_Ip2A != nullptr) {
                // [q_j + d, q_j + d_j) is the only new a-heavy output interval
                min->v = te_pair_t_v2 { ln_Ip2A, tn_NEW };
            } else {
                // there is no new a-heavy output interval
                T_e_v2.remove_node(min);

                if (min < &nodes_te.front() || &nodes_te.back() < min) {
                    delete min;
                }
            }
        }
    }

    if (log) {
        if (mf != nullptr) {
            *mf << " time_balance_phase_1=" << time_diff_ns(time)
                << " time_balance_phase_2=" << 0;
        }
        time = log_runtime(time);
    }
}