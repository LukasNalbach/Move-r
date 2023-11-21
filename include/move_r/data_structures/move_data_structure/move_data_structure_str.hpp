#pragma once

#include "move_data_structure.hpp"

/**
 * @brief move data structure, additionally stores a string interleaved with the arrays needed for move queries
 * @tparam uint_t unsigned integer type of the interval starting positions
 */
template <typename uint_t = uint32_t>
class move_data_structure_str : public move_data_structure<uint_t> {
    static_assert(std::is_same<uint_t,uint32_t>::value || std::is_same<uint_t,uint64_t>::value);

    /**
     * @brief builds the move_data_structure_str
     * @param I a disjoint interval sequence
     * @param n n = p_k + d_j
     * @param p the number of threads to use during the construction
     * @param a balancing parameter, restricts the number of intervals in the resulting move data structure to k*(a/(a-1))
     * @param log controls whether to print log messages during the construction
     * @param mf measurement file to write runtime data to
     * @param delete_i whether I can be deleted during the construction
     */
    void build(
        std::vector<std::pair<uint_t,uint_t>>& I,
        uint_t n,
        uint16_t p,
        uint16_t a,
        bool log,
        std::ostream* mf,
        bool delete_i
    ) {
        move_data_structure<uint_t>::is_move_data_structure_str = true;
        typename move_data_structure<uint_t>::construction mdsc(*reinterpret_cast<move_data_structure<uint_t>*>(this),I,n,p,a,delete_i,NULL,log,mf);
        set_character(move_data_structure<uint_t>::k_,0);
    }

    public:
    move_data_structure_str() = default;

    /**
     * @brief Constructs a new move data structure from a disjoint interval sequence
     * @param I a disjoint interval sequence
     * @param n n = p_k + d_j
     * @param p the number of threads to use during the construction
     * @param a balancing parameter, restricts the number of intervals in the resulting move data structure to k*(a/(a-1))
     * @param log controls whether to print log messages during the construction
     * @param mf measurement file to write runtime data to
     */
    move_data_structure_str(
        std::vector<std::pair<uint_t,uint_t>>& I,
        uint_t n,
        uint16_t p = omp_get_max_threads(),
        uint16_t a = 8,
        bool log = false,
        std::ostream* mf = NULL
    ) {
        build(I,n,p,a,log,mf,false);
    }

    /**
     * @brief Constructs a new move data structure from a disjoint interval sequence
     * @param I a disjoint interval sequence
     * @param n n = p_k + d_j
     * @param p the number of threads to use during the construction
     * @param a balancing parameter, restricts the number of intervals in the resulting move data structure to k*(a/(a-1))
     * @param log controls whether to print log messages during the construction
     * @param mf measurement file to write runtime data to
     */
    move_data_structure_str(
        std::vector<std::pair<uint_t,uint_t>>&& I,
        uint_t n,
        uint16_t p = omp_get_max_threads(),
        uint16_t a = 8,
        bool log = false,
        std::ostream* mf = NULL
    ) {
        build(I,n,p,a,log,mf,true);
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