#pragma once

#include <vector>
#include <iostream>
#include <sdsl/bit_vectors.hpp>

/**
 * @brief wrapper class for the sd_vector from sdsl
 * @tparam pos_t unsigned integer type
 */
template <typename pos_t = uint32_t>
class sd_array {
    static_assert(std::is_same_v<pos_t,uint32_t> || std::is_same_v<pos_t,uint64_t>);

    protected:
    sdsl::sd_vector<> sd_vector; // the sd_vector
    sdsl::sd_vector<>::rank_0_type rank_0_support; // rank_0 support for sd_vector
    sdsl::sd_vector<>::rank_1_type rank_1_support; // rank_1 support for sd_vector
    sdsl::sd_vector<>::select_0_type select_0_support; // select_1 support for sd_vector
    sdsl::sd_vector<>::select_1_type select_1_support; // select_1 support for sd_vector

    pos_t zeros = 0;
    pos_t ones = 0;

    /**
     * @brief copies another sd_array object into this object
     * @param other another sd_array object
     */
    void copy_from_other(const sd_array& other) {
        sd_vector = other.sd_vector;
        setup();
    }

    /**
     * @brief moves another sd_array object into this object
     * @param other another sd_array object
     */
    void move_from_other(sd_array&& other) {
        other.reset();
        sd_vector = std::move(other.sd_vector);
        setup();
    }

    /**
     * @brief sets rank_1-, select_0- and select_1-support to sd_vector
     */
    void setup() {
        rank_0_support.set_vector(&sd_vector);
        rank_1_support.set_vector(&sd_vector);
        select_0_support.set_vector(&sd_vector);
        select_1_support.set_vector(&sd_vector);

        if (size() > 0) {
            ones = rank_1(size());
            zeros = size()-ones;
        }
    }

    /**
     * @brief resets rank_1-, select_0- and select_1-support
     */
    void reset() {
        rank_0_support.set_vector(NULL);
        rank_1_support.set_vector(NULL);
        select_0_support.set_vector(NULL);
        select_1_support.set_vector(NULL);

        ones = 0;
        zeros = 0;
    }

    public:
    sd_array() {setup();}
    sd_array(sd_array&& other) {move_from_other(std::move(other));}
    sd_array(const sd_array& other) {copy_from_other(other);}
    sd_array& operator=(sd_array&& other) {move_from_other(std::move(other));return *this;}
    sd_array& operator=(const sd_array& other) {copy_from_other(other);return *this;}

    ~sd_array() {
        reset();
    }

    /**
     * @brief constructs a new sd_array from a bit vector
     * @param bit_vector a bit vector
     */
    sd_array(const sdsl::bit_vector& bit_vector) {
        sd_vector = std::move(sdsl::sd_vector<>(bit_vector));
        setup();
    }

    /**
     * @brief constructs a new sd_array from an sd_vector
     * @param bit_vector an sd_vector
     */
    sd_array(const sdsl::sd_vector<>& sd_vector) {
        this->sd_vector = sd_vector;
        setup();
    }

    /**
     * @brief constructs a new sd_array for an sd_vector and
     *        moves the sd_vector into the sd_array
     * @param bit_vector an sd_vector
     */
    sd_array(sdsl::sd_vector<>&& sd_vector) {
        this->sd_vector = std::move(sd_vector);
        setup();
    }

    /**
     * @brief returns the size of the bit vector
     * @return the size of the bit vector 
     */
    inline pos_t size() const {
        return sd_vector.size();
    }

    /**
     * @brief returns whether the input bit vector is empty
     * @return whether the bit vector is empty
     */
    inline bool empty() const {
        return size() == 0;
    }

    /**
     * @brief returns the number of ones in the bit vector
     * @return the number of ones in the bit vector 
     */
    inline pos_t num_ones() const {
        return ones;
    }

    /**
     * @brief returns the number of zeros in the bit vector
     * @return the number of ones in the bit vector 
     */
    inline pos_t num_zeros() const {
        return zeros;
    }

    /**
     * @brief returns the number of ones before index i
     * @param i [0..size]
     * @return the number of ones before index i 
     */
    inline pos_t rank_1(pos_t i) const {
        return rank_1_support.rank(i);
    }

    /**
     * @brief returns the index of the i-th one
     * @param i [1..number of ones]
     * @return the index of the i-th one 
     */
    inline pos_t select_1(pos_t i) const {
        return select_1_support.select(i);
    }

    /**
     * @brief returns the number of zeros before index i
     * @param i [0..size]
     * @return the number of zeros before index i 
     */
    inline pos_t rank_0(pos_t i) const {
        return rank_0_support.rank(i);
    }

    /**
     * @brief returns the index of the i-th zero
     * @param i [1..number of zeros]
     * @return the index of the i-th zero 
     */
    inline pos_t select_0(pos_t i) const {
        return select_0_support.select(i);
    }

    /**
     * @brief returns the index of the next one after index i
     * @param i [1..size-1]
     * @return the index of the next one after index i
     */
    inline pos_t next_1(pos_t i) const {
        return select_1(rank_1(i+1)+1);
    }

    /**
     * @brief returns the index of the previous one before index i
     * @param i [1..size-1]
     * @return the index of the previous one before index i
     */
    inline pos_t previous_1(pos_t i) const {
        return select_1(rank_1(i));
    }

    /**
     * @brief returns the index of the next zero after index i
     * @param i [1..size-1]
     * @return the index of the next zero after index i
     */
    inline pos_t next_0(pos_t i) const {
        return select_0(rank_0(i+1)+1);
    }

    /**
     * @brief returns the index of the previous zero before index i
     * @param i [1..size-1]
     * @return the index of the previous zero before index i
     */
    inline pos_t previous_0(pos_t i) const {
        return select_0(rank_0(i));
    }

    /**
     * @brief returns whether there is a one at index i
     * @param i [0..size-1]
     * @return whether there is a one at index i
     */
    inline bool operator[](pos_t i) const {
        return sd_vector[i];
    }

    /**
     * @brief returns the size of the data structure in bytes
     * @return size of the data structure in bytes
     */
    uint64_t size_in_bytes() const {
        return sdsl::size_in_bytes(sd_vector);
    }

    /**
     * @brief serializes the sd_array to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) const {
        sd_vector.serialize(out);
    }

    /**
     * @brief loads the sd_array from an input stream
     * @param in input stream
     */
    void load(std::istream& in) {
        sd_vector.load(in);
        setup();
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