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
 * @tparam pos_t unsigned integer type of the interval starting positions
 */
template <typename pos_t = uint32_t>
class move_data_structure {
    static_assert(std::is_same<pos_t,uint32_t>::value || std::is_same<pos_t,uint64_t>::value);

    protected:
    class construction;

    using pair_t = std::pair<pos_t,pos_t>;
    using pair_arr_t = std::vector<pair_t>;

    pos_t n = 0; // n = p_{k_'-1} + d_{k_'-1}, k_' <= n
    pos_t k = 0; // k, number of intervals in the original disjoint inteval sequence I
    pos_t k_ = 0; // k', number of intervals in the balanced disjoint inteval sequence B_a(I), k <= k_'
    uint16_t a = 0; // balancing parameter, restricts the number of intervals in the resulting move data structure to k*(a/(a-1))
    uint8_t omega_p = 0; // word width of D_p
    uint8_t omega_idx = 0; // word width of D_idx
    uint8_t omega_offs = 0; // word width of D_offs
    uint8_t omega_l_ = 0; // word width of L_
    interleaved_vectors<pos_t,pos_t> data; // interleaved vectors storing D_p, D_idx and D_offs (and L', for M_LF)
    
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
    move_data_structure(pair_arr_t&& I, pos_t n, mds_params params = {}, std::vector<pos_t>* pi_mphi = NULL) {
        construction(*this,I,n,true,0,params,pi_mphi);
    }

    /**
     * @brief Constructs a new move data structure from a disjoint interval sequence
     * @param I a disjoint interval sequence
     * @param n n = p_k + d_j
     * @param params construction parameters
     * @param pi_mphi vector to move pi into after the construction
     */
    move_data_structure(pair_arr_t& I, pos_t n, mds_params params = {}, std::vector<pos_t>* pi_mphi = NULL) {
        construction(*this,I,n,false,0,params,pi_mphi);
    }

    /**
     * @brief returns the size of the data structure in bytes
     * @return size of the data structure in bytes
     */
    uint64_t size_in_bytes() const {
        return
            1+2*sizeof(pos_t)+3+ // variables
            data.size_in_bytes(); // data
    }

    /**
     * @brief returns the maximum value n = p_k + d_k of the stored disjoint interval sequence
     * @return n = p_k + d_k
     */
    inline pos_t max_value() const {
        return n;
    }

    /**
     * @brief returns the number k' of intervals in the move data structure
     * @return k'
     */
    inline pos_t num_intervals() const {
        return k_;
    }

    /**
     * @brief returns a
     * @return a 
     */
    inline uint16_t balancing_parameter() const {
        return a;
    }

    /**
     * @brief returns the number omega_p of bits used by one entry in D_p (word width of D_p)
     * @return omega_p 
     */
    inline uint8_t width_p() const {
        return omega_p;
    }

    /**
     * @brief returns the number omega_p of bits used by one entry in D_idx (word width of D_idx)
     * @return omega_idx
     */
    inline uint8_t width_idx() const {
        return omega_idx;
    }

    /**
     * @brief returns the number omega_p of bits used by one entry in D_offs (word width of D_offs)
     * @return omega_offs
     */
    inline uint8_t width_offs() const {
        return omega_offs;
    }

    protected:
    /**
     * @brief resizes the move data structure to size k_
     * @param n maximum value
     * @param k_ size
     */
    void resize(pos_t n, pos_t k_, uint8_t omega_l_) {
        this->n = n;
        this->k_ = k_;

        omega_p = std::max((uint8_t)8,(uint8_t)(std::ceil(std::log2(n+1)/(double)8)*8));
        omega_idx = std::max((uint8_t)8,(uint8_t)(std::ceil(std::log2(k_+1)/(double)8)*8));
        this->omega_l_ = omega_l_;

        if (omega_l_ > 0) {
            data = std::move(interleaved_vectors<pos_t,pos_t>({
                (uint8_t)(omega_p/8),
                (uint8_t)(omega_idx/8),
                (uint8_t)(omega_offs/8),
                (uint8_t)(omega_l_/8)
            }));
        } else {
            data = std::move(interleaved_vectors<pos_t,pos_t>({
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
    inline void set_p(pos_t x, pos_t p) {
        data.template set<0,pos_t>(x,p);
    }

    /**
     * @brief sets D_idx[x] to p
     * @param x [0..k'-1] interval index
     * @param idx [0..k'-1] value to set D_idx[x] to
     */
    inline void set_idx(pos_t x, pos_t idx) {
        data.template set<1,pos_t>(x,idx);
    }

    /**
     * @brief sets D_offs[x] to offs
     * @param x [0..k'-1] interval index
     * @param offs [0..2^omega_offs-1] value to set D_offs[x] to
     */
    inline void set_offs(pos_t x, pos_t offs) {
        data.template set<2,pos_t>(x,offs);
    }

    public:
    /**
     * @brief returns D_p[x]
     * @param x [0..k_']
     * @return D_p[x]
     */
    inline pos_t p(pos_t x) const {
        return data.template get<0,pos_t>(x);
    }

    /**
     * @brief returns q_x
     * @param x [0..k_'-1]
     * @return q_x
     */
    inline pos_t q(pos_t x) const {
        return p(idx(x))+offs(x);
    }

    /**
     * @brief returns D_idx[x]
     * @param x [0..k_'-1]
     * @return D_idx[x]
     */
    inline pos_t idx(pos_t x) const {
        return data.template get<1,pos_t>(x);
    }

    /**
     * @brief returns D_offs[x]
     * @param x [0..k_'-1]
     * @return D_offs[x]
     */
    inline pos_t offs(pos_t x) const {
        return data.template get<2,pos_t>(x);
    }

    /**
     * @brief calculates the move query Move(I,i,x) = (i',x') by changing (i,x)
     *        to (i',x'), with i' = f_I(i) and i' in [p_x', p_x' + d_x')
     * @param i [0..n-1]
     * @param x [0..k_'-1], where i in [p_x, p_x + d_x)
     */
    inline void move(pos_t& i, pos_t& x) const {
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
    inline pair_t move(pair_t ix) const {
        move(ix.first,ix.second);
        return ix;
    }

    /**
     * @brief serializes the move data structure to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) const {
        out.write((char*)&n,sizeof(pos_t));
        out.write((char*)&k,sizeof(pos_t));
        out.write((char*)&k_,sizeof(pos_t));
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
        in.read((char*)&n,sizeof(pos_t));
        in.read((char*)&k,sizeof(pos_t));
        in.read((char*)&k_,sizeof(pos_t));
        in.read((char*)&a,sizeof(uint16_t));
        in.read((char*)&omega_p,1);
        in.read((char*)&omega_idx,1);
        in.read((char*)&omega_offs,1);
        data.load(in);
    }

    std::ostream& operator>>(std::ostream& os) const {
        serialize(os);
        return os;
    }

    std::istream& operator<<(std::istream& is) {
        load(is);
        return is;
    }
};

#include "construction/construction.hpp"