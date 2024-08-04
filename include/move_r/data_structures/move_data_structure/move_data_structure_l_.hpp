#pragma once

#include "move_data_structure.hpp"

/**
 * @brief move data structure, additionally stores a string interleaved with the arrays needed for move queries
 * @tparam pos_t unsigned integer type of the interval starting positions
 */
template <typename pos_t = uint32_t, typename l_t = char>
class move_data_structure_l_ : public move_data_structure<pos_t> {
    static_assert(std::is_same_v<pos_t,uint32_t> || std::is_same_v<pos_t,uint64_t>);

    using pair_t = typename move_data_structure<pos_t>::pair_t; // pair type
    using pair_arr_t = typename move_data_structure<pos_t>::pair_arr_t; // pair array type

    /**
     * @brief Constructs a new move data structure from a disjoint interval sequence
     * @param I a disjoint interval sequence
     * @param n n = p_k + d_j
     * @param delete_i whether I can be deleted during the construction
     * @param params construction parameters
     */
    void build(pair_arr_t& I, pos_t n, bool delete_i, mds_params params, uint8_t omega_l_) {
        typename move_data_structure<pos_t>::construction(
            *reinterpret_cast<move_data_structure<pos_t>*>(this),I,n,delete_i,omega_l_,params,NULL
        );
        set_L_(move_data_structure<pos_t>::k_,0);
    }

    public:
    move_data_structure_l_() = default;

    /**
     * @brief Constructs a new move data structure from a disjoint interval sequence
     * @param I a disjoint interval sequence
     * @param n n = p_k + d_j
     * @param params construction parameters
     */
    move_data_structure_l_(pair_arr_t&& I, pos_t n, mds_params params = {}, uint8_t omega_l_ = sizeof(l_t)) {
        build(I,n,true,params,omega_l_);
    }

    /**
     * @brief Constructs a new move data structure from a disjoint interval sequence
     * @param I a disjoint interval sequence
     * @param n n = p_k + d_j
     * @param params construction parameters
     */
    move_data_structure_l_(pair_arr_t& I, pos_t n, mds_params params = {}, uint8_t omega_l_ = sizeof(l_t)) {
        build(I,n,false,params,omega_l_);
    }

    /**
     * @brief returns the number omega_l_ of bits used by one entry in L' (word width of L')
     * @return omega_l_ 
     */
    inline uint8_t width_l_() const {
        return move_data_structure<pos_t>::omega_l_;
    }

    /**
     * @brief returns the value in L_ at position x
     * @param x index in [0..k_'-1]
     * @return the value in L_ at position x
     */
    inline l_t L_(pos_t x) const {
        if constexpr (sizeof(l_t) == 1) {
            return move_data_structure<pos_t>::data.template get_unsafe<3,l_t>(x);
        } else {
            return move_data_structure<pos_t>::data.template get<3,l_t>(x);
        }
    }

    /**
     * @brief sets the value in L_ at position x to c
     * @param x index in [0..k_'-1]
     * @param v a value
     */
    inline void set_L_(pos_t x, l_t v) {
        if constexpr (sizeof(l_t) == 1) {
            move_data_structure<pos_t>::data.template set_unsafe<3,l_t>(x,v);
        } else {
            move_data_structure<pos_t>::data.template set<3,l_t>(x,v);
        }
    }
};