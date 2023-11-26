#pragma once

#include <list>
#include <set>
#include <queue>
#include <omp.h>
#include <concurrentqueue.h>
#include <move_r/data_structures/avl_tree.hpp>
#include <move_r/data_structures/doubly_linked_list.hpp>
#include <move_r/data_structures/dynamic_insert_only_no_copy.hpp>
#include <absl/container/btree_set.h>

/**
 * @brief constructs a move data structure out of a disjoint interval sequence
 * @tparam uint_t unsigned integer type of the interval starting positions
 */
template <typename uint_t>
class move_data_structure<uint_t>::construction {
    static_assert(std::is_same<uint_t,uint32_t>::value || std::is_same<uint_t,uint64_t>::value);

    public:
    construction() = delete;
    construction(construction&&) = delete;
    construction(const construction&) = delete;
    construction& operator=(construction&&) = delete;
    construction& operator=(const construction&) = delete;
    ~construction() {}

    // ############################# COMMON TYPES #############################

    using pair_t = std::pair<uint_t,uint_t>;
    using pair_arr_t = std::vector<pair_t>;

    /**
     * @brief comparator for the pairs in T_in_v1 and T_in_v5
     */
    struct in_cmp_v1v5 {
        bool operator()(const pair_t &p1, const pair_t &p2) const {
            return p1.first < p2.first;
        }
    };

    /**
     * @brief comparator for the pairs in T_out_v1 and T_out_v5
     */
    struct out_cmp_v1v5 {
        bool operator()(const pair_t &p1, const pair_t &p2) const {
            return p1.second < p2.second;
        }
    };

    // ############################# COMMON VARIABLES #############################

    static constexpr uint8_t v = 5; // construction method
    /* 1 + epsilon is the maximum factor, by which the number of intervals can increase in the 
     * process of splitting too long intervals*/
    static constexpr double epsilon = 0.125;
    move_data_structure<uint_t>& mds; // the move data structure to construct
    pair_arr_t& I; // the disjoint interval sequence to construct the move data structure out of
    uint_t n; // maximum value, n = p_{k-1} + d_{k-1}, k <= n
    uint_t k; // number of intervals in the (possibly a-heavy) inteval sequence I, 0 < k
    uint_t k_; // number of intervals in the a-balanced inteval sequence B_a(I), 0 < k <= k'
    uint16_t a; // balancing parameter, restricts size increase to the factor (a/(a-1)), 2 <= a
    uint16_t two_a; // 2*a
    uint16_t p; // number of threads to use
    bool log; // toggles log messages
    bool delete_i;
    std::ostream* mf; // measurement file
    uint_t l_max; // maximum interval length
    uint64_t baseline_memory_allocation; // baseline memory allocation in bytes
    std::chrono::steady_clock::time_point time; // time point of the start of the last construction phase
    std::chrono::steady_clock::time_point time_start; // time point of the start of the entire construction
    interleaved_vectors<uint_t> D_q; // [0..k'-1] output interval starting positions (ordered by the input interval starting positions)
    std::vector<uint_t> pi; // [0..k'-1] permutation storing the order of the output interval starting postions

    /**
     * @brief [0..p-1] section start positions in the range [0..n], 0 = s[0] < s[1] < ... < s[p-1] = n.
     *        Before building T_out, s is chosen so that |L_in[0]| + |T_out[0]|
     *         ~ |L_in[1]| + |T_out[1]| ~ ... ~ |L_in[p-1]| + |T_out[p-1]|, that
     *        is s.t. s[i_p] = min {s' in [0,n-1], s.t. x[i_p] + u[i_p] - 2 >= i_p * lfloor 2k/p rfloor, where
                                x[i_p] = min {x' in [0,k-1], , s.t. p_{x'} >= s'} and
                                u[i_p] = min {u' in [0,k-1], s.t. q_{pi[u']} >= s'}
              } holds.
     */
    std::vector<uint_t> s;

    // [0..p], x[i] stores the number of input intervals in I starting before s[i]
    std::vector<uint_t> x;

    // [0..p], u[i] stores the number of output intervals in I starting before s[i]
    std::vector<uint_t> u;

    // ############################# COMMON METHODS #############################

    /**
     * @brief builds the move data structure mds
     * @param mds the move data structure to build
     * @param I disjoint interval sequence
     * @param n n = p_{k-1} + d_{k-1}, k <= n
     * @param k k = |I|
     * @param p number of threads to use
     * @param a balancing parameter, restricts size increase to the factor
     *          (1+1/(a-1)) and restricts move query runtime to 2a, 2 <= a
     * @param delete_i controls whether I should be deleted when not needed anymore
     * @param pi_mphi vector to move pi_mphi into after the construction
     * @param log enables log messages during build process
     * @param mf output stream to write runtime and space usage to if log is enabled
     */
    construction(
        move_data_structure<uint_t>& mds,
        std::vector<std::pair<uint_t,uint_t>>& I,
        uint_t n,
        uint16_t p = omp_get_max_threads(),
        uint16_t a = 8,
        bool delete_i = false,
        std::vector<uint_t>* pi_mphi = NULL,
        bool log = false,
        std::ostream* mf = NULL
    ) : mds(mds), I(I) {
        if (log) {
            time = now();
            time_start = time;
            std::cout << std::endl;
        }

        this->n = n;
        this->k = I.size();
        this->a = a;
        this->delete_i = delete_i;
        this->log = log;
        this->mf = mf;

        mds.a = a;
        mds.k = k;
        two_a = 2*a;

        if (p > 1 && 1000*p > k) {
            p = std::max((uint_t)1,k/1000);
            if (log) std::cout << "warning: p > k/1000, setting p to k/1000 ~ " << p << std::endl;
        }

        if (v < 3 && p > 1) {
            p = 1;
        }

        this->p = p;
        omp_set_num_threads(p);

        /* set omega_offs <- min {omega in {8,16,24,32,40} | n/(k2^omega) <= epsilon}, which ensures
         * k' <= k*(1+epsilon)*a/(a-1) */
        for (uint8_t omega=8; omega<=40; omega+=8) {
            if (n/((pow(2,omega)-1)*(double)k) <= epsilon) {
                mds.omega_offs = omega;
                break;
            }
        }

        /* in order to store D_offs with at most omega_offs bits, we have to limit the interval length
         * to l_max = 2^omega_offs */
        l_max = pow(2,mds.omega_offs)-1;

        // choose the correct construction method
        if constexpr (v == 1) {
            v1();
        } else if constexpr (2 <= v && v <= 4) {
            v2v3v4();
        } else {
            v5();
        }

        // verify the correctness of the construction, if in debug mode and printing log messages
        #ifndef NDEBUG
        if (log) {
            if constexpr (v == 1) build_pi_for_dq();
            verify_correctness();
        }
        #endif

        if (pi_mphi == NULL) {
            pi.clear();
            pi.shrink_to_fit();
        } else {
            *pi_mphi = std::move(pi);
            pi_mphi = NULL;
        }

        D_q.clear();
        D_q.shrink_to_fit();

        if (log) {
            log_message("move data structure built");
            log_runtime(time_start);
            if (mf != NULL) *mf << " time_total=" << time_diff_ns(time_start,time);
        }

        this->mf = NULL;
    }

    /**
     * @brief builds the permutation pi for the output interval starting positions stored in I
     */
    void build_pi_for_I();

    /**
     * @brief builds the permutation pi for the output interval starting positions stored in D_q
     */
    void build_pi_for_dq();

    /**
     * @brief calcualtes s, x and u for I
     */
    void calculate_seperation_positions_for_I();

    /**
     * @brief calculates s, x and u for D_p in mds and D_q
     */
    void calculate_seperation_positions_for_dq_and_mds();

    /**
     * @brief verifies the correctness of the resulting interval sequence
     */
    void verify_correctness();

    // ############################# V1 #############################
    
    using tin_t_v1 = avl_tree<pair_t,in_cmp_v1v5>;
    using tout_te_t_v1 = avl_tree<pair_t,out_cmp_v1v5>;
    using tin_node_t_v1 = typename tin_t_v1::avl_node;
    using tout_te_node_t_v1 = typename tout_te_t_v1::avl_node;
    
    // stores the pairs in I sorted by p_i
    tin_t_v1 T_in_v1;

    // stores the pairs in I sorted by q_i
    tout_te_t_v1 T_out_v1;

    /* stores the pairs in I sorted by p_i creating output intervals with at least 4 incoming edges in
    the permutation graph */
    tout_te_t_v1 T_e_v1;

    /**
     * @brief builds the move data structure mds using the construction method v1
     */
    void v1() {
        // Build T_in_v1 and T_out_v1
        build_tin_tout_v1();

        balance_v1_seq();

        // Build D_p and D_q
        build_dp_dq_v1();

        // Build D_offs and D_idx
        build_didx_doffs_v1();
    };

    /**
     * @brief builds T_in_v1 and T_out_v1
     */
    void build_tin_tout_v1();

    /**
     * @brief balances the disjoint interval sequence stored in T_in_v1 and T_out_v1
     */
    void balance_v1_seq();

    /**
     * @brief builds D_p in mds and D_q out of the disjoint interval sequence stored in T_in_v1 and T_out_v1 
     */
    void build_dp_dq_v1();

    /**
     * @brief builds D_idx and D_offs in mds using binary searches over D_p in mds
     */
    void build_didx_doffs_v1();

    // ############################# V2/V3/V4 SEQUENTIAL/PARALLEL #############################

    using lin_t_v2v3v4 = doubly_linked_list<pair_t>;
    using lin_node_t_v2v3v4 = typename lin_t_v2v3v4::doubly_linked_list_node;

    /**
     * @brief comparator for the list nodes in T_out_v2, T_e_v2, T_out_v3 and T_out_v5
     */
    struct tout_cmp_v2v3v4 {
        bool operator()(const lin_node_t_v2v3v4 &n1, const lin_node_t_v2v3v4 &n2) const {
            return n1.v.second < n2.v.second;
        }
    };
    
    using tout_t_v2v3v4 = avl_tree<lin_node_t_v2v3v4,tout_cmp_v2v3v4>;
    using tout_node_t_v2v3v4 = typename tout_t_v2v3v4::avl_node;

    /** 
     * @brief [0..p-1] doubly linked lists; L_in_v2v3v4[i_p] stores the pairs (p_i,q_i) in ascending order of p_i,
     *        where s[i_p] <= p_i < s[i_p+1] and i_p in [0..p-1]. L_in_v2v3v4[0]L_in_v2v3v4[1]...L_in_v2v3v4[p-1] = I.
     */
    std::vector<lin_t_v2v3v4> L_in_v2v3v4;
    /**
     * @brief [0..p-1] avl trees; T_out_v2v3v4[i_p] stores nodes of lists in L_in_v2v3v4 in ascending order of q_i,
     *        for each pair (p_i,q_i), s[i_p] <= q_i < s[i_p+1] holds, with i_p in [0..p-1].
     */
    std::vector<tout_t_v2v3v4> T_out_v2v3v4;
    /**
     * @brief stores the nodes in L_in_v2v3v4[0..p-1] and T_out_v2v3v4[0..p-1] which they were initially created with
     */
    std::vector<tout_node_t_v2v3v4> nodes_v2v3v4;
    /**
     * @brief [0..p-1] new_nodes_2v3v4[i_p] stores the newly created nodes in L_in_v2v3v4[0..p-1] and T_out_v2v3v4[0..p-1],
     *        which were created by thread i_p.
     */
    std::vector<dynamic_insert_only_no_copy<tout_node_t_v2v3v4>> new_nodes_2v3v4;

    /**
     * @brief builds the move data structure mds using the construction method v2, v3 or v4
     */
    void v2v3v4() {
        // Build L_in_v2v3v4 and T_out_v2v3v4
        build_lin_tout_v2v3v4();

        // Choose the correct balancing algorithm
        if (p == 1) {
            if constexpr (v == 2) {
                balance_v2_seq();
            } else {
                balance_v3_seq();
            }
        } else {
            if constexpr (v == 3) {
                balance_v3_par();
            } else {
                balance_v4_par();
            }
        }

        // Build D_p and D_q
        build_dp_dq_v2v3v4();

        // Build D_offs and D_idx
        build_didx_doffs_v2v3v4v5();
    }

    /**
     * @brief builds L_in_v2v3v4[0..p-1] and T_out_v2v3v4[0..p-1] out of the disjoint interval sequence I
     */
    void build_lin_tout_v2v3v4();

    /**
     * @brief inserts the pairs in L_in_v2v3v4[0..p-1] into D_pair and writes D_idx
     */
    void build_dp_dq_v2v3v4();

    /**
     * @brief builds D_offs and D_idx
     */
    void build_didx_doffs_v2v3v4v5();

    /**
     * @brief returns the length of the input/output interval starting at p_i/q_i
     * @param ln (p_i,q_i)
     * @return length of the input/output interval starting at p_i/q_i
     */
    inline uint_t interval_length_v2v3_seq(lin_node_t_v2v3v4 *ln);

    /**
     * @brief checks if the output interval [q_j, q_j + d_j) is a-heavy and iterates
     *        ln_IpI_ to the last output interval connected to it in the permutation graph
     * @param ln_IpI_ (p_{i+i_},q_{i+i_}), [p_i, p_i + d_i) must be the first input
     *                 interval connected to [q_j, q_j + d_j) in the permutation graph
     * @param i_ 1 <= i_ <= 2a
     * @param tn_J (p_j,q_j)
     * @param tn_J_ (p_{j'},q_{j'}), with q_j + d_j = q_{j'}
     * @return (p_{i+a},q_{i+a}) if [q_j, q_j + d_j) is a-heavy, else NULL
     */
    inline lin_node_t_v2v3v4* is_a_heavy_v2v3v4(
        lin_node_t_v2v3v4 **ln_IpI_,
        uint_t* i_,
        tout_node_t_v2v3v4 *tn_J,
        tout_node_t_v2v3v4 *tn_J_ = NULL
    );

    // ############################# V2 #############################
    
    using te_pair_t_v2 = std::pair<lin_node_t_v2v3v4*,tout_node_t_v2v3v4*>;

    /**
     * @brief comparator for the list nodes in T_out_v2, T_e_v2, T_out_v3 and T_out_v5
     */
    struct te_cmp_v2 {
        bool operator()(const te_pair_t_v2 &p1, const te_pair_t_v2 &p2) const {
            return p1.second->v.v.second < p2.second->v.v.second;
        }
    };

    using te_t_v2 = avl_tree<te_pair_t_v2,te_cmp_v2>;
    using te_node_t_v2 = typename te_t_v2::avl_node;

    /* 
        contains pairs (lin_node_t_v2v3v4 *p1, lin_node_t_v2v3v4 *p2),
        where p2 is associated with an output interval with at least 2a incoming edges in the permutation graph
        and p1 is the pair associated with the a+1-st input interval in the output interval associated with p2.
        The pairs are ordered by the starting position of the a-heavy output intervals associated with p2.
    */
    te_t_v2 T_e_v2;

    /**
     * @brief balances the disjoint interval sequence in L_in_v2v3v4[0] and T_out_v2v3v4[0] sequentially
     */
    void balance_v2_seq();

    // ############################# V3 SEQUENTIAL #############################

    /**
     * @brief balances the output interval [q_j, q_j + d_j) and all a-heavy output intervals
     *        starting before q_u that have become a-heavy in the process
     * @param ln_IpA (p_{i+a},q_{i+a}), [p_i, p_i + d_i) must be the first input interval connected
     *                to [q_j, q_j + d_j) in the permutation graph
     * @param tn_J (p_j,q_j), [q_j, q_j + d_j) must be the first a-heavy output interval
     * @param q_u starting position of an output interval (q_j + d <= q_u)
     * @param p_cur starting position of the current input interval
     * @param i_ the current input interval is the i_-th input interval in [q_j, q_j + d_j)
     * @return the newly created pair (p_j+d,q_j+d)
     */
    inline tout_node_t_v2v3v4* balance_upto_v3_seq(
        lin_node_t_v2v3v4 *ln_IpA,
        tout_node_t_v2v3v4 *tn_J,
        uint_t q_u,
        uint_t p_cur,
        uint_t* i_
    );

    /**
     * @brief balances the disjoint interval sequence in L_in_v2v3v4[0] and T_out_v2v3v4[0] sequentially
     */
    void balance_v3_seq();

    // ############################# V3/V4 PARALLEL #############################

    using q_node_t_v34 = std::pair<lin_node_t_v2v3v4*,lin_node_t_v2v3v4*>;

    // ############################# V3 PARALLEL #############################

    using q_t_v3 = std::vector<std::vector<std::queue<q_node_t_v34>>>;

    /** @brief [0..p-1][0..p-1] stores queues with tuples (*p1,*p2);
     *         Q_v3[i] stores the tuples to insert into thread i's section [s[i]..s[i+1]) */
    q_t_v3 Q_v3;
    
    /** @brief swap variable for Q_v3 */
    q_t_v3 Q_v3_;

    /**
     * @brief balances the output interval [q_j, q_j + d_j) and all a-heavy output intervals starting
     *        in [s[i_p],s[i_p+1]) and before q_u that have become a-heavy in the process by inserting
     *        the newly created pair into T_out_v2v3v4[i_p] and either Q_v3[0..p-1][i_p] or L_in_v2v3v4[i_p]
     * @param ln_IpA (p_{i+a},q_{i+a}), [p_i, p_i + d_i) must be the first input interval connected
     *               to [q_j, q_j + d_j) in the permutation graph
     * @param tn_J (p_j,q_j), [q_j, q_j + d_j) must be the first a-heavy output interval starting 
     *             in [s[i_p]..q_u]
     * @param tn_J_ (p_j',q_j'), [q_j', q_j' + d_j') must be the first output interval starting
     *              after [q_j, q_j + d_j)
     * @param q_u starting position of an output interval starting in [s[i_p]..s[i_p+1]] (q_j + d <= q_u)
     * @param p_cur starting position of the current input interval
     * @param i_ the current input interval is the i_-th input interval in [q_j, q_j + d_j)
     * @return the newly created pair (p_j+d,q_j+d)
     */
    inline tout_node_t_v2v3v4* balance_upto_v3_par(
        lin_node_t_v2v3v4* ln_IpA,
        tout_node_t_v2v3v4* tn_J,
        tout_node_t_v2v3v4* tn_J_,
        uint_t q_u,
        uint_t p_cur,
        uint_t* i_
    );

    /**
     * @brief balances the disjoint interval sequence in L_in_v2v3v4[0..p-1] and T_out_v2v3v4[0..p-1] in parallel
     */
    void balance_v3_par();

    // ############################# V4 PARALLEL #############################

    using q_t_v4 = std::vector<std::vector<moodycamel::ConcurrentQueue<q_node_t_v34>>>;
    
    /** @brief [0..p-1][0..p-1] stores queues with tuples (*p1,*p2);
     *         Q_v4[i] stores the tuples to insert into thread i's section [s[i]..s[i+1]) */
    q_t_v4 Q_v4;

    /**
     * @brief balances the output interval [q_j, q_j + d_j) and all a-heavy output intervals starting
     *        in [s[i_p],s[i_p+1]) and before q_u that have become a-heavy in the process by inserting
     *        the newly created pair into T_out_v2v3v4[i_p] and either Q_v4[0..p-1] or L_in_v2v3v4[i_p]
     * @param ln_IpA (p_{i+a},q_{i+a}), [p_i, p_i + d_i) must be the first input interval connected
     *               to [q_j, q_j + d_j) in the permutation graph
     * @param tn_J (p_j,q_j), [q_j, q_j + d_j) must be the first a-heavy output interval starting
     *             in [s[i_p]..q_u]
     * @param tn_J_ (p_j',q_j'), [q_j', q_j' + d_j') must be the first output interval starting
     *              after [q_j, q_j + d_j)
     * @param q_u starting position of an output interval starting in [s[i_p]..s[i_p+1]] (q_j + d <= q_u)
     * @param p_cur starting position of the current input interval
     * @param i_ the current input interval is the i_-th input interval in [q_j, q_j + d_j)
     * @return the newly created pair (p_j+d,q_j+d)
     */
    inline tout_node_t_v2v3v4* balance_upto_v4_par(
        lin_node_t_v2v3v4 *ln_IpA,
        tout_node_t_v2v3v4 *tn_J,
        tout_node_t_v2v3v4* tn_J_,
        uint_t q_u,
        uint_t p_cur,
        uint_t* i_
    );

    /**
     * @brief balances the disjoint interval sequence in L_in_v2v3v4[0..p-1] and T_out_v2v3v4[0..p-1] in parallel
     */
    void balance_v4_par();

    // ############################# V5 SEQUENTIAL/PARALLEL #############################

    using tin_t_v5 = absl::btree_set<pair_t,in_cmp_v1v5>;
    using tout_t_v5 = absl::btree_set<pair_t,out_cmp_v1v5>;

    using tin_it_t_v5 = typename tin_t_v5::iterator;
    using tout_it_t_v5 = typename tout_t_v5::iterator;

    /** 
     * @brief [0..p-1] b-trees; T_in_v5[i_p] stores the pairs (p_i,q_i) in ascending order of p_i,
     *        where s[i_p] <= p_i < s[i_p+1] and i_p in [0..p-1]. T_in_v5[0]T_in_v5[1]...T_in_v5[p-1] = I.
     */
    std::vector<tin_t_v5> T_in_v5;

    /**
     * @brief [0..p-1] b-trees; T_out_v5[i_p] stores the pairs (p_i,q_i) in ascending order of q_i,
     *        where s[i_p] <= q_i < s[i_p+1] and i_p in [0..p-1].
     */
    std::vector<tout_t_v5> T_out_v5;

    /**
     * @brief [0..p-1][0..p-1] b-trees temporarily storing the pairs that have been inserted into
     *        T_in_v5[0..p-1] in order to split the too long intervals; T_out_temp_v5[i][j] stores the
     *        pairs that have already been inserted into T_in_v5[j] and have to be inserted into T_out_v5[i]
     */
    std::vector<std::vector<tout_t_v5>> T_out_temp_v5;

    using q_t_v5 = std::vector<std::vector<pair_arr_t>>;

    /**
     * @brief [0..p-1][0..p-1] stores vectors of pairs;
     *        Q_v5[i] stores the pairs to insert into thread i's section [s[i]..s[i+1]) */
    q_t_v5 Q_v5;

    /** @brief swap variable for Q_v5 */
    q_t_v5 Q_v5_;

    // returns the node of an insert result into a b-tree in T_in_v5
    tin_it_t_v5 node(std::pair<tin_it_t_v5,bool> insert_result) {
        return insert_result.first;
    }

    /**
     * @brief builds the move data structure mds using the construction method v5
     */
    void v5() {
        // build T_in_v5 and T_out_v5
        build_tin_tout_v5();
        
        // balance the disjoint interval sequence stored in T_in_v5 and T_out_v5
        balance_v5_seq_par();

        // Build D_p and D_q
        build_dp_dq_v5();

        // Build D_offs and D_idx
        build_didx_doffs_v2v3v4v5();
    }

    /**
     * @brief builds T_in_v5 and T_out_v5
     */
    void build_tin_tout_v5();

    /**
     * @brief builds D_p and D_q
     */
    void build_dp_dq_v5();

    /**
     * @brief checks wheter the output interval [q_j, q_j + d_j) that starts directly before q_j'
     *        is a-heavy; if [q_j, q_j + d_j) is a-heavy, it returns q_j + d (where d = p_{i+a})
     *        and iterates tn_I further to at most the pair that creates the 2a-th input interval
     *        starting in [q_j, q_j+ d_j); else it iterates tn_I further to the first input interval
     *        starting after q_j' and returns 0
     * @param tn_I an iterator pointing to the pair (p_i,q_i)
     * @param q_J_ q_j' the starting position of the output interval that starts directly after [q_j, q_j+ d_j)
     * @return q_j + d (where d = p_{i+a}); the position at which [q_j, q_j + d_j) has to be split
     */
    inline uint_t is_a_heavy_v5_seq_par(tin_it_t_v5& tn_I, uint_t q_J_);

    /**
     * @brief balances the output interval [q_j, q_j + d_j) and all a-heavy output intervals in [s[i_p],s[i_p+1])
     *        starting before q_u that have become a-heavy in the process by inserting the newly created pair into
     *        T_out_v5[i_p] and Q_v5[0..p-1][i_p]
     * @param tn_J_ iterator to the pair (p_j',q_j') in T_out_v5, [q_j', q_j' + d_j') must be the first
     *              output interval starting after [q_j, q_j + d_j)
     * @param qj_pd the position to split [q_j, q_j + d_j) at; [q_j, q_j + d_j) must be the first a-heavy output
     *              interval starting in [s[i_p],s[i_p+1]) and before q_u
     * @param q_u starting position of an output interval starting in [s[i_p]..s[i_p+1]] (q_j + d <= q_u)
     * @return an iterator pointing to the newly created pair (p_j + d, q_j + d) in T_out_v5, if no recursive call
     *         has been made in this call of balance_upto_v5_seq_par, else returns an iterator pointing to T_out_v5[i_p].end()
     */
    inline tout_it_t_v5 balance_upto_v5_seq_par(tout_it_t_v5& tn_J_, uint_t qj_pd, uint_t q_u);

    /**
     * @brief balances the disjoint interval sequence in L_in_v5[0..p-1] and T_out_v5[0..p-1] sequentially or in parallel
     */
    void balance_v5_seq_par();
};

#include "algorithms/construction/v1.cpp"
#include "algorithms/construction/v1v2v3v4v5.cpp"
#include "algorithms/construction/v2v3v4.cpp"
#include "algorithms/construction/v5.cpp"

#include "algorithms/balancing/v1_seq.cpp"
#include "algorithms/balancing/v2_seq.cpp"
#include "algorithms/balancing/v3_seq.cpp"
#include "algorithms/balancing/v3_par.cpp"
#include "algorithms/balancing/v4_par.cpp"
#include "algorithms/balancing/v5_seq_par.cpp"

#include "algorithms/misc/verify_correctness.cpp"