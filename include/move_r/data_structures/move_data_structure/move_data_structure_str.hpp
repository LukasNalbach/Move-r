#pragma once

#include "move_data_structure.hpp"

/**
 * @brief move data structure, additionally stores a string interleaved with the arrays needed for move queries
 * @tparam uint_t unsigned integer type of the interval starting positions
 */
template <typename uint_t = uint32_t>
class move_data_structure_str : public move_data_structure<uint_t> {
    static_assert(std::is_same<uint_t,uint32_t>::value || std::is_same<uint_t,uint64_t>::value);

    using pair_t = typename move_data_structure<uint_t>::pair_t;
    using pair_arr_t = typename move_data_structure<uint_t>::pair_arr_t;

    /**
     * @brief Constructs a new move data structure from a disjoint interval sequence
     * @param I a disjoint interval sequence
     * @param n n = p_k + d_j
     * @param delete_i whether I can be deleted during the construction
     * @param params construction parameters
     * @param pi_mphi vector to move pi into after the construction
     */
    void build(pair_arr_t& I, uint_t n, bool delete_i, mds_params params, std::vector<uint_t>* pi_mphi) {
        move_data_structure<uint_t>::is_move_data_structure_str = true;
        typename move_data_structure<uint_t>::construction(
            *reinterpret_cast<move_data_structure<uint_t>*>(this),I,n,delete_i,params,pi_mphi
        );
        set_character(move_data_structure<uint_t>::k_,0);
    }

    public:
    move_data_structure_str() = default;

    /**
     * @brief Constructs a new move data structure from a disjoint interval sequence
     * @param I a disjoint interval sequence
     * @param n n = p_k + d_j
     * @param params construction parameters
     * @param pi_mphi vector to move pi into after the construction
     */
    move_data_structure_str(pair_arr_t&& I, uint_t n, mds_params params = {}, std::vector<uint_t>* pi_mphi = NULL) {
        build(I,n,true,params,pi_mphi);
    }

    /**
     * @brief Constructs a new move data structure from a disjoint interval sequence
     * @param I a disjoint interval sequence
     * @param n n = p_k + d_j
     * @param params construction parameters
     * @param pi_mphi vector to move pi into after the construction
     */
    move_data_structure_str(pair_arr_t& I, uint_t n, mds_params params = {}, std::vector<uint_t>* pi_mphi = NULL) {
        build(I,n,false,params,pi_mphi);
    }

    /**
     * @brief returns the character at position x
     * @param x index in [0..k_'-1]
     * @return the character at position x
     */
    inline char character(uint_t x) {
        return move_data_structure<uint_t>::data.template get_unsafe<3,char>(x);
    }

    /**
     * @brief sets the character at position x to c
     * @param x index in [0..k_'-1]
     * @param c a character
     */
    inline void set_character(uint_t x, char c) {
        move_data_structure<uint_t>::data.template set_unsafe<3,char>(x,c);
    }
};