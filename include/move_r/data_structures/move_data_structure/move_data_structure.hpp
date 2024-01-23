#pragma once

#include <cmath>
#include <omp.h>
#include <move_r/misc/utils.hpp>
#include <move_r/data_structures/interleaved_vectors.hpp>

struct mds_params {
    uint16_t num_threads = omp_get_max_threads(); // the number of threads to use during the construction
    uint16_t a = 8; // balancing parameter, restricts the number of intervals in the resulting move data structure to k*(a/(a-1))
    bool log = false; // controls whether to print log messages during the construction
    std::ostream* mf = NULL; // measurement file to write runtime data to
};

/**
 * @brief move data structure
 * @tparam uint_t unsigned integer type of the interval starting positions
 */
template <typename uint_t = uint32_t>
class move_data_structure {
    static_assert(std::is_same<uint_t,uint32_t>::value || std::is_same<uint_t,uint64_t>::value);

    protected:
    class construction;

    using pair_t = std::pair<uint_t,uint_t>;
    using pair_arr_t = std::vector<pair_t>;

    bool is_move_data_structure_str = false; // true <=> this move data structure is of type move_data_structure_str
    uint_t n = 0; // n = p_{k_'-1} + d_{k_'-1}, k_' <= n
    uint_t k = 0; // k, number of intervals in the original disjoint inteval sequence I
    uint_t k_ = 0; // k', number of intervals in the balanced disjoint inteval sequence B_a(I), k <= k_'
    uint16_t a = 0; // balancing parameter, restricts the number of intervals in the resulting move data structure to k*(a/(a-1))
    uint8_t omega_p = 0; // word width of D_p
    uint8_t omega_idx = 0; // word width of D_idx
    uint8_t omega_offs = 0; // word width of D_offs
    interleaved_vectors<uint_t> data; // interleaved vectors storing D_p, D_idx and D_offs (and L', for M_LF)
    
    public:
    move_data_structure() = default;
    move_data_structure(move_data_structure&& other) = default;
    move_data_structure(const move_data_structure& other) = default;
    move_data_structure& operator=(move_data_structure&& other) = default;
    move_data_structure& operator=(const move_data_structure& other) = default;

    /**
     * @brief Constructs a new move data structure from a disjoint interval sequence
     * @param I a disjoint interval sequence
     * @param n n = p_k + d_j
     * @param params construction parameters
     * @param pi_mphi vector to move pi into after the construction
     */
    move_data_structure(pair_arr_t&& I, uint_t n, mds_params params = {}, std::vector<uint_t>* pi_mphi = NULL) {
        construction(*this,I,n,true,params,pi_mphi);
    }

    /**
     * @brief Constructs a new move data structure from a disjoint interval sequence
     * @param I a disjoint interval sequence
     * @param n n = p_k + d_j
     * @param params construction parameters
     * @param pi_mphi vector to move pi into after the construction
     */
    move_data_structure(pair_arr_t& I, uint_t n, mds_params params = {}, std::vector<uint_t>* pi_mphi = NULL) {
        construction(*this,I,n,false,params,pi_mphi);
    }

    /**
     * @brief returns the size of the data structure in bytes
     * @return size of the data structure in bytes
     */
    uint64_t size_in_bytes() {
        return
            1+2*sizeof(uint_t)+3+ // variables
            data.size_in_bytes(); // data
    }

    /**
     * @brief returns the maximum value n = p_k + d_k of the stored disjoint interval sequence
     * @return n = p_k + d_k
     */
    inline uint_t max_value() {
        return n;
    }

    /**
     * @brief returns the number k' of intervals in the move data structure
     * @return k'
     */
    inline uint_t num_intervals() {
        return k_;
    }

    /**
     * @brief returns a
     * @return a 
     */
    inline uint16_t balancing_parameter() {
        return a;
    }

    /**
     * @brief returns the number omega_p of bits used by one entry in D_p (word width of D_p)
     * @return omega_p 
     */
    inline uint8_t width_p() {
        return omega_p;
    }

    /**
     * @brief returns the number omega_p of bits used by one entry in D_idx (word width of D_idx)
     * @return omega_idx
     */
    inline uint8_t width_idx() {
        return omega_idx;
    }

    /**
     * @brief returns the number omega_p of bits used by one entry in D_offs (word width of D_offs)
     * @return omega_offs
     */
    inline uint8_t width_offs() {
        return omega_offs;
    }

    protected:
    /**
     * @brief resizes the move data structure to size k_
     * @param n maximum value
     * @param k_ size
     */
    void resize(uint_t n, uint_t k_) {
        this->n = n;
        this->k_ = k_;

        omega_p = std::max((uint8_t)8,(uint8_t)(std::ceil(std::log2(n+1)/(double)8)*8));
        omega_idx = std::max((uint8_t)8,(uint8_t)(std::ceil(std::log2(k_+1)/(double)8)*8));

        if (is_move_data_structure_str) {
            data = std::move(interleaved_vectors<uint_t>({
                (uint8_t)(omega_p/8),
                (uint8_t)(omega_idx/8),
                (uint8_t)(omega_offs/8),
                (uint8_t)1
            }));
        } else {
            data = std::move(interleaved_vectors<uint_t>({
                (uint8_t)(omega_p/8),
                (uint8_t)(omega_idx/8),
                (uint8_t)(omega_offs/8)
            }));
        }

        data.resize_no_init(k_+1);
        set_p(k_,n);
        set_idx(k_,k_);
        set_offs(k_,0);
    }

    /**
     * @brief sets D_p[x] to p
     * @param x [0..k'-1] interval index
     * @param p [0..n-1] value to set D_p[x] to
     */
    inline void set_p(uint_t x, uint_t p) {
        data.template set<0>(x,p);
    }

    /**
     * @brief sets D_idx[x] to p
     * @param x [0..k'-1] interval index
     * @param idx [0..k'-1] value to set D_idx[x] to
     */
    inline void set_idx(uint_t x, uint_t idx) {
        data.template set<1>(x,idx);
    }

    /**
     * @brief sets D_offs[x] to offs
     * @param x [0..k'-1] interval index
     * @param offs [0..2^omega_offs-1] value to set D_offs[x] to
     */
    inline void set_offs(uint_t x, uint_t offs) {
        data.template set<2>(x,offs);
    }

    public:
    /**
     * @brief returns D_p[x]
     * @param x [0..k_']
     * @return D_p[x]
     */
    inline uint_t p(uint_t x) {
        return data.template get<0>(x);
    }

    /**
     * @brief returns q_x
     * @param x [0..k_'-1]
     * @return q_x
     */
    inline uint_t q(uint_t x) {
        return p(idx(x))+offs(x);
    }

    /**
     * @brief returns D_idx[x]
     * @param x [0..k_'-1]
     * @return D_idx[x]
     */
    inline uint_t idx(uint_t x) {
        return data.template get<1>(x);
    }

    /**
     * @brief returns D_offs[x]
     * @param x [0..k_'-1]
     * @return D_offs[x]
     */
    inline uint_t offs(uint_t x) {
        return data.template get<2>(x);
    }

    /**
     * @brief calculates the move query Move(I,i,x) = (i',x') by changing (i,x)
     *        to (i',x'), with i' = f_I(i) and i' in [p_x', p_x' + d_x')
     * @param i [0..n-1]
     * @param x [0..k_'-1], where i in [p_x, p_x + d_x)
     */
    inline void move(uint_t& i, uint_t& x) {
        i = q(x)+(i-p(x));
        x = idx(x);
        while (i >= p(x+1)) {
            x++;
        }
    }

    /**
     * @brief performs the move query Move(I,i,x) = (i',x') by changing ix = (i,x) to
     *        (i',x'), with i' = f_I(i) and i' in [p_x', p_x' + d_x') and returns (i',x')
     * @param i [0..n-1]
     * @param x [0..k_'-1], where i in [p_x, p_x + d_x)
     * @returns
     */
    inline pair_t move(pair_t ix) {
        move(ix.first,ix.second);
        return ix;
    }

    /**
     * @brief serializes the move data structure to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) {
        out.write((char*)&n,sizeof(uint_t));
        out.write((char*)&k,sizeof(uint_t));
        out.write((char*)&k_,sizeof(uint_t));
        out.write((char*)&a,sizeof(uint16_t));
        out.write((char*)&omega_p,1);
        out.write((char*)&omega_idx,1);
        out.write((char*)&omega_offs,1);
        data.serialize(out);
    }

    /**
     * @brief loads the move data structure from an input stream
     * @param in input stream
     */
    void load(std::istream& in) {
        in.read((char*)&n,sizeof(uint_t));
        in.read((char*)&k,sizeof(uint_t));
        in.read((char*)&k_,sizeof(uint_t));
        in.read((char*)&a,sizeof(uint16_t));
        in.read((char*)&omega_p,1);
        in.read((char*)&omega_idx,1);
        in.read((char*)&omega_offs,1);
        data.load(in);
    }
};

#include "construction/construction.hpp"