#pragma once

#include <vector>
#include <iostream>
#include <sdsl/bit_vectors.hpp>

/**
 * @brief wrapper class for the bit_vector from sdsl
 * @tparam uint_t unsigned integer type
 */
template <typename uint_t = uint32_t, bool build_rank_support = false, bool build_select_0_support = false, bool build_select_1_support = false>
class plain_bit_vector {
    static_assert(std::is_same<uint_t,uint32_t>::value || std::is_same<uint_t,uint64_t>::value);

    protected:
    sdsl::bit_vector vec; // the bit_vector
    sdsl::bit_vector::rank_1_type rank_1_support; // rank_1 support for vec
    sdsl::bit_vector::select_1_type select_1_support; // select_1 support for vec
    sdsl::bit_vector::select_0_type select_0_support; // select_1 support for vec

    uint_t zeros = 0;
    uint_t ones = 0;

    /**
     * @brief copies another plain_bit_vector object into this object
     * @param other another plain_bit_vector object
     */
    void copy_from_other(const plain_bit_vector& other) {
        vec = other.vec;
        rank_1_support = other.rank_1_support;
        select_0_support = other.select_0_support;
        select_1_support = other.select_1_support;

        rank_1_support.set_vector(&vec);
        select_0_support.set_vector(&vec);
        select_1_support.set_vector(&vec);
    }

    /**
     * @brief moves another plain_bit_vector object into this object
     * @param other another plain_bit_vector object
     */
    void move_from_other(plain_bit_vector&& other) {
        vec = std::move(other.vec);
        rank_1_support = std::move(other.rank_1_support);
        select_0_support = std::move(other.select_0_support);
        select_1_support = std::move(other.select_1_support);

        rank_1_support.set_vector(&vec);
        select_0_support.set_vector(&vec);
        select_1_support.set_vector(&vec);

        ones = other.ones;
        zeros = other.zeros;

        other.reset();
    }

    /**
     * @brief builds rank_1-, select_0- and select_1-support for vec
     */
    void setup() {
        if constexpr (build_rank_support) rank_1_support = sdsl::bit_vector::rank_1_type(&vec);
        if constexpr (build_select_0_support) select_0_support = sdsl::bit_vector::select_0_type(&vec);
        if constexpr (build_select_1_support) select_1_support = sdsl::bit_vector::select_1_type(&vec);

        if (size() > 0) {
            if constexpr (build_rank_support) {
                ones = rank_1(size());
            } else {
                ones = sdsl::bit_vector::rank_1_type(&vec).rank(size());
            }
            
            zeros = size()-ones;
        }
    }

    /**
     * @brief resets rank_1-, select_0- and select_1-support
     */
    void reset() {
        rank_1_support.set_vector(NULL);
        select_0_support.set_vector(NULL);
        select_1_support.set_vector(NULL);

        ones = 0;
        zeros = 0;
    }

    public:
    plain_bit_vector() {setup();}
    plain_bit_vector(plain_bit_vector&& other) {move_from_other(std::move(other));}
    plain_bit_vector(const plain_bit_vector& other) {copy_from_other(other);}
    plain_bit_vector& operator=(plain_bit_vector&& other) {move_from_other(std::move(other));return *this;}
    plain_bit_vector& operator=(const plain_bit_vector& other) {copy_from_other(other);return *this;}

    ~plain_bit_vector() {
        reset();
    }

    /**
     * @brief constructs a new plain_bit_vector from a bit vector
     * @param bit_vector a bit vector
     */
    plain_bit_vector(const sdsl::bit_vector& vec) {
        this->vec = vec;
        setup();
    }

    /**
     * @brief constructs a new plain_bit_vector from a bit vector
     * @param bit_vector a bit vector
     */
    plain_bit_vector(sdsl::bit_vector&& vec) {
        this->vec = std::move(vec);
        setup();
    }

    /**
     * @brief returns the size of the bit vector
     * @return the size of the bit vector 
     */
    uint_t size() {
        return vec.size();
    }

    /**
     * @brief returns the number of ones in the bit vector
     * @return the number of ones in the bit vector 
     */
    uint_t num_ones() {
        return ones;
    }

    /**
     * @brief returns the number of zeros in the bit vector
     * @return the number of ones in the bit vector 
     */
    uint_t num_zeros() {
        return zeros;
    }

    /**
     * @brief returns the number of ones before index i
     * @param i [0..size]
     * @return the number of ones before index i 
     */
    inline uint_t rank_1(uint_t i) {
        return rank_1_support.rank(i);
    }

    /**
     * @brief returns the index of the i-th one
     * @param i [1..number of ones]
     * @return the index of the i-th one 
     */
    inline uint_t select_1(uint_t i) {
        return select_1_support.select(i);
    }

    /**
     * @brief returns the number of zeros before index i
     * @param i [0..size]
     * @return the number of zeros before index i 
     */
    inline uint_t rank_0(uint_t i) {
        return i-(rank_1(i));
    }

    /**
     * @brief returns the index of the i-th zero
     * @param i [1..number of zeros]
     * @return the index of the i-th zero 
     */
    inline uint_t select_0(uint_t i) {
        return select_0_support.select(i);
    }

    /**
     * @brief returns the index of the next one after index i
     * @param i [1..size-1]
     * @return the index of the next one after index i
     */
    inline uint_t next_1(uint_t i) {
        return select_1(rank_1(i+1)+1);
    }

    /**
     * @brief returns the index of the previous one before index i
     * @param i [1..size-1]
     * @return the index of the previous one before index i
     */
    inline uint_t previous_1(uint_t i) {
        return select_1(rank_1(i));
    }

    /**
     * @brief returns the index of the next zero after index i
     * @param i [1..size-1]
     * @return the index of the next zero after index i
     */
    inline uint_t next_0(uint_t i) {
        return select_0(rank_0(i+1)+1);
    }

    /**
     * @brief returns the index of the previous zero before index i
     * @param i [1..size-1]
     * @return the index of the previous zero before index i
     */
    inline uint_t previous_0(uint_t i) {
        return select_0(rank_0(i));
    }

    /**
     * @brief returns whether there is a one at index i
     * @param i [0..size-1]
     * @return whether there is a one at index i
     */
    inline bool operator[](uint_t i) {
        return vec[i];
    }

    /**
     * @brief returns the size of the data structure in bytes
     * @return size of the data structure in bytes
     */
    uint64_t size_in_bytes() {
        return sdsl::size_in_bytes(vec)+
               sdsl::size_in_bytes(rank_1_support)+
               sdsl::size_in_bytes(select_0_support)+
               sdsl::size_in_bytes(select_1_support);
    }

    /**
     * @brief serializes the plain_bit_vector to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) {
        out.write((char*)&ones,sizeof(uint_t));
        out.write((char*)&zeros,sizeof(uint_t));

        vec.serialize(out);

        rank_1_support.serialize(out);
        select_0_support.serialize(out);
        select_1_support.serialize(out);
    }

    /**
     * @brief loads the plain_bit_vector from an input stream
     * @param in input stream
     */
    void load(std::istream& in) {
        in.read((char*)&ones,sizeof(uint_t));
        in.read((char*)&zeros,sizeof(uint_t));
        
        vec.load(in);

        rank_1_support.load(in);
        select_0_support.load(in);
        select_1_support.load(in);
        
        rank_1_support.set_vector(&vec);
        select_0_support.set_vector(&vec);
        select_1_support.set_vector(&vec);
    }
};