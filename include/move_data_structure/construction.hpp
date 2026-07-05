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

#include <misc/files.hpp>
#include <misc/search.hpp>
#include <misc/log.hpp>

#include <gtl/btree.hpp>
#include <omp.h>

/**
 * @brief constructs a move data structure out of a disjoint interval sequence
 * @tparam pos_t unsigned integer type of the interval starting positions
 */
template <typename pos_t>
class move_data_structure<pos_t>::construction {
    static_assert(std::is_same_v<pos_t, uint32_t> || std::is_same_v<pos_t, uint64_t>);

public:
    construction() = delete;
    construction(construction&&) = delete;
    construction(const construction&) = delete;
    construction& operator=(construction&&) = delete;
    construction& operator=(const construction&) = delete;
    ~construction() { }

    // ############################# TYPES #############################

    /**
     * @brief comparator that orders the pairs (p_i,q_i) by their input interval starting position p_i
     */
    struct in_cmp {
        bool operator()(const pair_t& p1, const pair_t& p2) const
        {
            return p1.first < p2.first;
        }
    };

    /**
     * @brief comparator that orders the pairs (p_i,q_i) by their output interval starting position q_i
     */
    struct out_cmp {
        bool operator()(const pair_t& p1, const pair_t& p2) const
        {
            return p1.second < p2.second;
        }
    };

    using tin_t = gtl::btree_set<pair_t, in_cmp>;
    using tout_t = gtl::btree_set<pair_t, out_cmp>;

    using tin_it_t = typename tin_t::iterator;
    using tout_it_t = typename tout_t::iterator;

    using q_t = std::vector<std::vector<pair_arr_t>>;

    // ############################# VARIABLES #############################

    /* 1 + epsilon is the maximum factor, by which the number of intervals can increase in the
     * process of splitting too long intervals*/
    static constexpr double epsilon = 0.125;
    move_data_structure<pos_t>& mds; // the move data structure to construct
    pair_arr_t& I; // the disjoint interval sequence to construct the move data structure out of
    pos_t n; // maximum value, n = p_{k-1} + d_{k-1}, k <= n
    pos_t k; // number of intervals in the (possibly a-heavy) inteval sequence I, 0 < k
    pos_t k_; // number of intervals in the a-balanced inteval sequence B_a(I), 0 < k <= k'
    uint16_t a; // balancing parameter, restricts size increase to the factor (a/(a-1)), 2 <= a
    uint16_t two_a; // 2*a
    uint16_t p; // number of threads to use
    bool log; // toggles log messages
    bool delete_i; // controls whether I should be deleted when not needed anymore
    std::ostream* mf; // measurement file
    pos_t l_max; // maximum interval length
    uint64_t baseline_mem_usage; // baseline memory allocation in bytes
    std::chrono::steady_clock::time_point time; // time point of the start of the last construction phase
    std::chrono::steady_clock::time_point time_start; // time point of the start of the entire construction
    interleaved_byte_aligned_vectors<pos_t, pos_t> D_q; // [0..k'-1] output interval starting positions (ordered by the input interval starting positions)
    std::vector<pos_t> pi; // [0..k'-1] permutation storing the order of the output interval starting postions
    uint8_t width_l_; // width of L_

    /**
     * @brief [0..p-1] section start positions in the range [0..n], 0 = s[0] < s[1] < ... < s[p-1] = n.
     *        Before building T_out, s is chosen so that |T_in[0]| + |T_out[0]|
     *         ~ |T_in[1]| + |T_out[1]| ~ ... ~ |T_in[p-1]| + |T_out[p-1]|, that
     *        is s.t. s[i_p] = min {s' in [0,n-1], s.t. x[i_p] + u[i_p] - 2 >= i_p * lfloor 2k/p rfloor, where
                                x[i_p] = min {x' in [0,k-1], , s.t. p_{x'} >= s'} and
                                u[i_p] = min {u' in [0,k-1], s.t. q_{pi[u']} >= s'}
              } holds.
     */
    std::vector<pos_t> s;

    // [0..p], x[i] stores the number of input intervals in I starting before s[i]
    std::vector<pos_t> x;

    // [0..p], u[i] stores the number of output intervals in I starting before s[i]
    std::vector<pos_t> u;

    /**
     * @brief [0..p-1] b-trees; T_in[i_p] stores the pairs (p_i,q_i) in ascending order of p_i,
     *        where s[i_p] <= p_i < s[i_p+1] and i_p in [0..p-1]. T_in[0]T_in[1]...T_in[p-1] = I.
     */
    std::vector<tin_t> T_in;

    /**
     * @brief [0..p-1] b-trees; T_out[i_p] stores the pairs (p_i,q_i) in ascending order of q_i,
     *        where s[i_p] <= q_i < s[i_p+1] and i_p in [0..p-1].
     */
    std::vector<tout_t> T_out;

    /**
     * @brief [0..p-1][0..p-1] b-trees temporarily storing the pairs that have been inserted into
     *        T_in[0..p-1] in order to split the too long intervals; T_out_temp[i][j] stores the
     *        pairs that have already been inserted into T_in[j] and have to be inserted into T_out[i]
     */
    std::vector<std::vector<tout_t>> T_out_temp;

    /** @brief [0..p-1][0..p-1] stores vectors of pairs;
     *         Q[i] stores the pairs to insert into thread i's section [s[i]..s[i+1]) */
    q_t Q;

    /** @brief swap variable for Q */
    q_t Q_;

    // ############################# METHODS #############################

    /**
     * @brief builds the move data structure mds
     * @param mds the move data structure to build
     * @param I disjoint interval sequence
     * @param n n = p_{k-1} + d_{k-1}, k <= n
     * @param delete_i controls whether I should be deleted when not needed anymore
     * @param width_l_ bit width of the L_ values stored in the move data structure
     * @param params construction parameters
     * @param pi_mphi vector to move pi into after the construction
     */
    construction(
        move_data_structure<pos_t>& mds,
        pair_arr_t& I,
        pos_t n,
        bool delete_i,
        uint8_t width_l_,
        mds_params params,
        std::vector<pos_t>* pi_mphi = nullptr)
        : mds(mds)
        , I(I)
    {
        this->n = n;
        this->k = I.size();
        this->a = params.a;
        this->p = params.num_threads;
        this->delete_i = delete_i;
        this->log = params.log;
        this->mf = params.mf;
        this->width_l_ = width_l_;

        if (log) {
            time = now();
            time_start = time;
            std::cout << std::endl;
        }

        mds.a = a;
        mds.k = k;
        two_a = 2 * a;

        if (p > 1 && 1000 * p > k) {
            p = std::max<pos_t>(1, k / 1000);
            if (log)
                std::cout << "warning: p > k/1000, setting p to k/1000 ~ " << p << std::endl;
        }

        omp_set_num_threads(p);
        // tight bit-width for D_offs (its values are limited to l_max = 2^omega_offs - 1 below)
        mds.omega_offs = std::bit_width((uint64_t)((8.0 * n) / k));

        /* in order to store D_offs with at most omega_offs bits, we have to limit the interval length
         * to l_max = 2^omega_offs */
        l_max = (pos_t{1} << mds.omega_offs) - 1;

        // build the move data structure
        build();

        // verify the correctness of the construction, if in debug mode and printing log messages
        #ifndef NDEBUG
        if (log) {
            verify_correctness();
        }
        #endif

        if (pi_mphi == nullptr) {
            pi.clear();
            pi.shrink_to_fit();
        } else {
            *pi_mphi = std::move(pi);
        }

        D_q.clear();
        D_q.shrink_to_fit();

        if (log) {
            log_message("move data structure built");
            log_runtime(time_start);
            if (mf != nullptr)
                *mf << " time_total=" << time_diff_ns(time_start, time);
        }

        this->mf = nullptr;
    }

    /**
     * @brief builds the move data structure mds from the disjoint interval sequence I
     */
    void build()
    {
        // build T_in and T_out
        build_tin_tout();

        // balance the disjoint interval sequence stored in T_in and T_out
        balance();

        // build D_p and D_q
        build_dp_dq();

        // build D_offs and D_idx
        build_didx_doffs();
    }

    // ############################# COMMON METHODS #############################

    /**
     * @brief builds the permutation pi for the output interval starting positions stored in I
     */
    void build_pi_for_I();

    /**
     * @brief builds the permutation pi for the output interval starting positions stored in D_q
     */
    void build_pi_for_dq();

    /**
     * @brief calculates s, x and u for I
     */
    void calculate_seperation_positions_for_I();

    /**
     * @brief calculates s, x and u for D_p in mds and D_q
     */
    void calculate_seperation_positions_for_dq_and_mds();

    /**
     * @brief builds D_offs and D_idx
     */
    void build_didx_doffs();

    /**
     * @brief verifies the correctness of the resulting interval sequence
     */
    void verify_correctness();

    // ############################# INTERVAL SEQUENCE CONSTRUCTION #############################

    /**
     * @brief returns the node of an insert result into a b-tree in T_in
     * @param insert_result the result of an emplace/insert into a b-tree in T_in
     * @return an iterator pointing to the inserted (or already present) pair
     */
    tin_it_t node(std::pair<tin_it_t, bool> insert_result)
    {
        return insert_result.first;
    }

    /**
     * @brief builds T_in[0..p-1] and T_out[0..p-1] out of the disjoint interval sequence I
     */
    void build_tin_tout();

    /**
     * @brief builds D_p in mds and D_q out of the disjoint interval sequence stored in T_in and T_out
     */
    void build_dp_dq();

    // ############################# BALANCING #############################

    /**
     * @brief checks wheter the output interval [q_j, q_j + d_j) that starts directly before q_j'
     *        is a-heavy; if [q_j, q_j + d_j) is a-heavy, it returns q_j + d (where d = p_{i+a})
     *        and iterates tn_I further to at most the pair that creates the 2a-th input interval
     *        starting in [q_j, q_j+ d_j); else it iterates tn_I further to the first input interval
     *        starting after q_j' and returns 0
     * @param tn_I an iterator pointing to the pair (p_i,q_i)
     * @param q_J_ q_j' the starting position of the output interval that starts directly after [q_j, q_j+ d_j)
     * @return q_j + d (where d = p_{i+a}); the position at which [q_j, q_j + d_j) has to be split, if
     *         [q_j, q_j + d_j) is a-heavy, or 0 if it is a-balanced
     */
    inline pos_t is_a_heavy(tin_it_t& tn_I, pos_t q_J_);

    /**
     * @brief balances the output interval [q_j, q_j + d_j) and all a-heavy output intervals in [s[i_p],s[i_p+1])
     *        starting before q_u that have become a-heavy in the process by inserting the newly created pair into
     *        T_out[i_p] and Q[0..p-1][i_p]
     * @param tn_J_ iterator to the pair (p_j',q_j') in T_out, [q_j', q_j' + d_j') must be the first
     *              output interval starting after [q_j, q_j + d_j)
     * @param qj_pd the position to split [q_j, q_j + d_j) at; [q_j, q_j + d_j) must be the first a-heavy output
     *              interval starting in [s[i_p],s[i_p+1]) and before q_u
     * @param q_u starting position of an output interval starting in [s[i_p]..s[i_p+1]] (q_j + d <= q_u)
     * @return an iterator pointing to the newly created pair (p_j + d, q_j + d) in T_out, if no recursive call
     *         has been made in this call of balance_upto, else returns an iterator pointing to T_out[i_p].end()
     */
    inline tout_it_t balance_upto(tout_it_t& tn_J_, pos_t qj_pd, pos_t q_u);

    /**
     * @brief balances the disjoint interval sequence in T_in[0..p-1] and T_out[0..p-1] sequentially or in parallel
     */
    void balance();
};

#include "construction.tpp"
