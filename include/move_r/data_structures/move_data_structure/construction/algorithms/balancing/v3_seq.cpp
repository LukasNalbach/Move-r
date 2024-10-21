#pragma once

#include <move_r/data_structures/move_data_structure/move_data_structure.hpp>

template <typename pos_t>
inline typename move_data_structure<pos_t>::construction::tout_node_t_v2v3v4* move_data_structure<pos_t>::construction::balance_upto_v3_seq(
    lin_node_t_v2v3v4* ln_IpA,
    tout_node_t_v2v3v4* tn_J,
    pos_t q_u,
    pos_t p_cur,
    pos_t* i_)
{
    pos_t p_j = tn_J->v.v.first;
    pos_t q_j = tn_J->v.v.second;
    pos_t d_j = interval_length_v2v3_seq(&tn_J->v);

    // d = p_{i+2a} - q_j is the maximum integer, so that [q_j, q_j + d) has a incoming edges in the permutation graph.
    pos_t d = ln_IpA->v.first - q_j;

    // Create the pair (p_j + d, q_j + d), which creates two new input intervals [p_j, p_j + d) and [p_j + d, p_j + d_j).
    tout_node_t_v2v3v4* tn_NEW = new_nodes_2v3v4[0].emplace_back(tout_node_t_v2v3v4(lin_node_t_v2v3v4(pair_t { p_j + d, q_j + d })));
    T_out_v2v3v4[0].insert_hint(tn_NEW, tn_J);
    L_in_v2v3v4[0].insert_after_node(&tn_NEW->v, &tn_J->v);

    /*
        Case 1: pj + d ≥ qu + du
            It is irrelevant, if [qy, qy + qy) is a-heavy, since it starts after qu.
        Case 2: pj + d ∈ [qu, qu + du)
            Then [qy, qy +dy) = [qu, qu+du), so there are a+1 input intervals starting in [qu, qu+du),
            hence, it and all output intervals starting before qu are balanced.
        Case 3: pj + d < qu
            Case 3.1: pj + d ∈ [qj , qj + dj)
                Then [qy, qy + dy) equals either [qj , qj + d) or [qj + d, qj + dj). Since before inserting
                (pj + d, qj + d), there were < 2a input intervals starting in [qj , qj + dj), there are now
                ≤ a + 1 input intervals starting in [qj , qj + d) and [qj + d, qj + dj) each, hence, they
                are both balanced.
            Case 3.2: pj + d /∈ [qj , qj + dj)
                [qy, qy + dy) is possibly a-heavy.
    */

    // check if case 3.2 holds
    if (p_j + d < q_u) {
        if (p_j + d < q_j || q_j + d_j <= p_j + d) {
            // find [q_y, q_y + d_y)
            tout_node_t_v2v3v4* tn_Y = T_out_v2v3v4[0].max_leq(lin_node_t_v2v3v4(pair_t { 0, p_j + d }));

            // find the first input interval [p_z, p_z + d_z) connected to [q_y, q_y + d_y) in the permutation graph
            lin_node_t_v2v3v4* ln_Z = &tn_NEW->v;
            pos_t i__ = 1;

            while (ln_Z->pr != NULL && ln_Z->pr->v.first >= tn_Y->v.v.second) {
                ln_Z = ln_Z->pr;
                i__++;
            }

            ln_Z = &tn_NEW->v;

            /* check if [q_y, q_y + d_y) is a-heavy and if yes, balance it and all output intervals starting before it,
               that might get a-heavy in the process. */
            lin_node_t_v2v3v4* ln_ZpA = is_a_heavy_v2v3v4(&ln_Z, &i__, tn_Y);

            if (ln_ZpA != NULL) {
                balance_upto_v3_seq(ln_ZpA, tn_Y, q_u, p_cur, i_);
            }
        }
    } else if (p_j + d < p_cur) {
        (*i_)++;
    }

    return tn_NEW;
}

/**
 * @brief balances the disjoint interval sequence in L_in_v2v3v4[0] and T_out_v2v3v4[0] sequentially
 */
template <typename pos_t>
void move_data_structure<pos_t>::construction::balance_v3_seq()
{
    if (log)
        log_message("balancing");

    // points to to the pair (p_i,q_i).
    lin_node_t_v2v3v4* ln_I = L_in_v2v3v4[0].head();
    // points to the pair (p_j,q_j).
    typename tout_t_v2v3v4::avl_it tn_J = T_out_v2v3v4[0].iterator();
    // points to the pair (p_{j'},q_{j'}), where q_j + d_j = q_{j'}.
    typename tout_t_v2v3v4::avl_it tn_Jp1 = T_out_v2v3v4[0].iterator(T_out_v2v3v4[0].second_smallest());

    // temporary variables
    lin_node_t_v2v3v4* ln_IpA;
    pos_t i_ = 1;

    // At the start of each iteration, [p_i, p_i + d_i) is the first input interval connected to [q_j, q_j + d_j) in the permutation graph
    bool stop = false;
    
    while (!stop) {
        ln_IpA = is_a_heavy_v2v3v4(&ln_I, &i_, tn_J.current(), tn_Jp1.current());

        // If [q_j, q_j + d_j) is a-heavy, balance it and all output intervals starting before it becoming might get a-heavy in the process.
        if (ln_IpA != NULL) {
            tn_J.set(balance_upto_v3_seq(ln_IpA, tn_J.current(), tn_J.current()->v.v.second, ln_I->v.first, &i_));
            continue;
        }

        // Find the next output interval with an incoming edge in the permutation graph and the first input interval connected to it.
        do {
            if (!tn_Jp1.has_next()) {
                stop = true;
                break;
            }
            tn_J.set(tn_Jp1.current());
            tn_Jp1.next();

            while (ln_I->v.first < tn_J.current()->v.v.second) {
                if (ln_I->sc == NULL) {
                    stop = true;
                    break;
                }

                ln_I = ln_I->sc;
            }
        } while (!stop && ln_I->v.first >= tn_Jp1.current()->v.v.second);

        i_ = 1;
    }

    if (log) {
        if (mf != NULL) {
            *mf << " time_balance_phase_1=" << time_diff_ns(time)
                << " time_balance_phase_2=" << 0;
        }
        time = log_runtime(time);
    }
}