#pragma once

#include <cmath>
#include <omp.h>
#include <move_r/data_structures/interleaved_vectors.hpp>

/**
 * @brief move data structure
 * @tparam uint_t unsigned integer type of the interval starting positions
 */
template <typename uint_t = uint32_t>
class move_data_structure_phi {
	static_assert(std::is_same<uint_t,uint32_t>::value || std::is_same<uint_t,uint64_t>::value);

    protected:
    class construction;

    bool is_move_data_structure_lf = false; // true <=> this move data structure is of type move_data_structure_lf
    uint_t n = 0; // n = p_{k_'-1} + d_{k_'-1}, k_' <= n
    uint_t k_ = 0; // k', number of intervals in the balanced disjoint inteval sequence B_a(I), k <= k_'
    uint8_t omega_p = 0; // word width of D_p
    uint8_t omega_idx = 0; // word width of D_idx
    uint8_t omega_offs = 0; // word width of D_offs
    interleaved_vectors<uint_t> data; // interleaved vectors storing D_p, D_idx and D_offs (and L', for M_LF)
    
    public:
    move_data_structure_phi() = default;
    move_data_structure_phi(move_data_structure_phi&& other) = default;
    move_data_structure_phi(const move_data_structure_phi& other) = default;
    move_data_structure_phi& operator=(move_data_structure_phi&& other) = default;
    move_data_structure_phi& operator=(const move_data_structure_phi& other) = default;

    /**
     * @brief Constructs a new move data structure from a disjoint interval sequence
     * @param I a disjoint interval sequence
     * @param n n = p_k + d_j
     * @param p the number of threads to use during the construction
     * @param a balancing parameter, restricts the number of intervals in the resulting move data structure to k*(a/(a-1))
     * @param log controls whether to print log messages during the construction
     * @param mf measurement file to write runtime data to
     */
    move_data_structure_phi(
        std::vector<std::pair<uint_t,uint_t>>& I,
        uint_t n,
        uint16_t p = omp_get_max_threads(),
        uint16_t a = 8,
        bool log = false,
        std::ostream* mf = NULL
    ) : move_data_structure_phi(std::move(I),n,p,a,log,mf) {}

    /**
     * @brief Constructs a new move data structure from a disjoint interval sequence
     * @param I a disjoint interval sequence
     * @param n n = p_k + d_j
     * @param p the number of threads to use during the construction
     * @param a balancing parameter, restricts the number of intervals in the resulting move data structure to k*(a/(a-1))
     * @param log controls whether to print log messages during the construction
     * @param mf measurement file to write runtime data to
     */
    move_data_structure_phi(
        std::vector<std::pair<uint_t,uint_t>>&& I,
        uint_t n,
        uint16_t p = omp_get_max_threads(),
        uint16_t a = 8,
        bool log = false,
        std::ostream* mf = NULL
    ) {
        construction mdsc(*this,std::move(I),n,p,a,log,mf);
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
    inline uint_t sequence_width() {
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
     * @brief returns the number omega_p of bytes used by one entry in D_p (word width of D_p)
     * @return omega_p 
     */
    inline uint8_t width_p() {
        return omega_p;
    }

    /**
     * @brief returns the number omega_p of bytes used by one entry in D_idx (word width of D_idx)
     * @return omega_idx
     */
    inline uint8_t width_idx() {
        return omega_idx;
    }

    /**
     * @brief returns the number omega_p of bytes used by one entry in D_offs (word width of D_offs)
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

        omega_p = std::max((uint8_t)1,(uint8_t)(ceil(log2(n)/(double)8)*8));
        omega_idx = std::max((uint8_t)1,(uint8_t)(ceil(log2(k_)/(double)8)*8));

        if (is_move_data_structure_lf) {
            data = std::move(interleaved_vectors<uint_t>({
                (uint8_t)(omega_p/8),
                (uint8_t)(omega_idx/8),
                (uint8_t)(omega_offs/8),
                (uint8_t)1
            },k_+1,false));
        } else {
            data = std::move(interleaved_vectors<uint_t>({
                (uint8_t)(omega_p/8),
                (uint8_t)(omega_idx/8),
                (uint8_t)(omega_offs/8)
            },k_+1,false));
        }

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
     * @brief returns D_idx[x]
     * @param x [0..k_'-1]
     * @return D_idx[x]
     */
    inline uint_t offs(uint_t x) {
        return data.template get<2>(x);
    }

    /**
     * @brief performs the move query Move(I,i,x) = (i',x') by changing ix = (i,x) to
     *        (i',x'), with i' = f_I(i) and i' in [p_x', p_x' + d_x') and returns (i',x')
     * @param i [0..n-1]
     * @param x [0..k_'-1], where i in [p_x, p_x + d_x)
     * @returns
     */
    inline std::pair<uint_t,uint_t> move(std::pair<uint_t,uint_t> ix) {
        ix.first = q(ix.second)+(ix.first-p(ix.second));
        ix.second = idx(ix.second);
        while (ix.first >= p(ix.second+1)) {
            ix.second++;
        }
        return ix;
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
     * @brief serializes the move data structure to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) {
        out.write((char*)&n,sizeof(uint_t));
        out.write((char*)&k_,sizeof(uint_t));
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
        in.read((char*)&k_,sizeof(uint_t));
        in.read((char*)&omega_p,1);
        in.read((char*)&omega_idx,1);
        in.read((char*)&omega_offs,1);
        data.load(in);
    }
};

#include "construction/construction.hpp"