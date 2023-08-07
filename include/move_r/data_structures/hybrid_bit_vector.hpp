#pragma once

#include <vector>
#include <iostream>
#include <optional>
#include "sd_array.hpp"
#include "plain_bit_vector.hpp"

/**
 * @brief hybrid bit vector (either an sd_array or a plain bit vector)
 * @tparam uint_t unsigned integer type
 */
template <typename uint_t = uint32_t, bool build_rank_support = false, bool build_select_0_support = false, bool build_select_1_support = false>
class hybrid_bit_vector {
    static_assert(std::is_same<uint_t,uint32_t>::value || std::is_same<uint_t,uint64_t>::value);

    using plain_bv_t = plain_bit_vector<uint_t,build_rank_support,build_select_0_support,build_select_1_support>;

    /* if there are less than (compression_threshold * size of the bitvector)
     * ones in the bit vector, then it should be compressed */
    static constexpr double compression_threshold = 0.1;

    protected:
    std::optional<sd_array<uint_t>> sd_arr; // the sd_array
    std::optional<plain_bv_t> plain_bit_vec; // the plain bit vector

    /**
     * @brief returns, whether plain_bit_vec should be compressed
     * @param plain_bit_vec a bit vector 
     * @return whether plain_bit_vec should be compressed
     */
    bool should_be_compressed(const sdsl::bit_vector& plain_bit_vec) {
        return sdsl::bit_vector::rank_1_type(&plain_bit_vec).rank(plain_bit_vec.size()) < plain_bit_vec.size()*compression_threshold;
    }

    public:
    hybrid_bit_vector() = default;

    /**
     * @brief constructs a new hybrid_bit_vector from a bit vector
     * @param bit_vector a bit vector
     */
    hybrid_bit_vector(const sdsl::bit_vector& plain_bit_vec) {
        if (should_be_compressed(plain_bit_vec)) {
            sd_arr = std::move(sd_array<uint_t>(plain_bit_vec));
        } else {
            this->plain_bit_vec = std::move(plain_bv_t(plain_bit_vec));
        }
    }

    /**
     * @brief constructs a new hybrid_bit_vector from a bit vector
     * @param bit_vector a bit vector
     */
    hybrid_bit_vector(sdsl::bit_vector&& plain_bit_vec) {
        if (should_be_compressed(plain_bit_vec)) {
            sd_arr = std::move(sd_array<uint_t>(plain_bit_vec));
        } else {
            this->plain_bit_vec = std::move(plain_bv_t(std::move(plain_bit_vec)));
        }
    }

    inline bool is_initialized() {
        return sd_arr.has_value() || plain_bit_vec.has_value();
    }

    /**
     * @brief returns, whether the bit vector is compressed
     * @return whether the bit vector is compressed
     */
    inline bool is_compressed() {
        return sd_arr.has_value();
    }

    /**
     * @brief returns the size of the bit vector
     * @return the size of the bit vector 
     */
    inline uint_t size() {
        if (is_compressed()) {
            return sd_arr.value().size();
        } else {
            return plain_bit_vec.value().size();
        }
    }

    /**
     * @brief returns the number of ones in the bit vector
     * @return the number of ones in the bit vector 
     */
    inline uint_t num_ones() {
        if (is_compressed()) {
            return sd_arr.value().num_ones();
        } else {
            return plain_bit_vec.value().num_ones();
        }
    }

    /**
     * @brief returns the number of zeros in the bit vector
     * @return the number of ones in the bit vector 
     */
    inline uint_t num_zeros() {
        if (is_compressed()) {
            return sd_arr.value().num_zeros();
        } else {
            return plain_bit_vec.value().num_zeros();
        }
    }

    /**
     * @brief returns the number of ones before index i
     * @param i [0..size]
     * @return the number of ones before index i 
     */
    inline uint_t rank_1(uint_t i) {
        static_assert(build_rank_support);
        if (is_compressed()) {
            return sd_arr.value().rank_1(i);
        } else {
            return plain_bit_vec.value().rank_1(i);
        }
    }

    /**
     * @brief returns the index of the i-th one
     * @param i [1..number of ones]
     * @return the index of the i-th one 
     */
    inline uint_t select_1(uint_t i) {
        static_assert(build_select_1_support);
        if (is_compressed()) {
            return sd_arr.value().select_1(i);
        } else {
            return plain_bit_vec.value().select_1(i);
        }
    }

    /**
     * @brief returns the number of zeros before index i
     * @param i [0..size]
     * @return the number of zeros before index i 
     */
    inline uint_t rank_0(uint_t i) {
        static_assert(build_rank_support);
        if (is_compressed()) {
            return sd_arr.value().rank_0(i);
        } else {
            return plain_bit_vec.value().rank_0(i);
        }
    }

    /**
     * @brief returns the index of the i-th zero
     * @param i [1..number of zeros]
     * @return the index of the i-th zero 
     */
    inline uint_t select_0(uint_t i) {
        static_assert(build_select_0_support);
        if (is_compressed()) {
            return sd_arr.value().select_0(i);
        } else {
            return plain_bit_vec.value().select_0(i);
        }
    }

    /**
     * @brief returns the index of the next one after index i
     * @param i [1..size-1]
     * @return the index of the next one after index i
     */
    inline uint_t next_1(uint_t i) {
        if (is_compressed()) {
            return sd_arr.value().next_1(i);
        } else {
            return plain_bit_vec.value().next_1(i);
        }
    }

    /**
     * @brief returns the index of the previous one before index i
     * @param i [1..size-1]
     * @return the index of the previous one before index i
     */
    inline uint_t previous_1(uint_t i) {
        if (is_compressed()) {
            return sd_arr.value().previous_1(i);
        } else {
            return plain_bit_vec.value().previous_1(i);
        }
    }

    /**
     * @brief returns the index of the next zero after index i
     * @param i [1..size-1]
     * @return the index of the next zero after index i
     */
    inline uint_t next_0(uint_t i) {
        if (is_compressed()) {
            return sd_arr.value().next_0(i);
        } else {
            return plain_bit_vec.value().next_0(i);
        }
    }

    /**
     * @brief returns the index of the previous zero before index i
     * @param i [1..size-1]
     * @return the index of the previous zero before index i
     */
    inline uint_t previous_0(uint_t i) {
        if (is_compressed()) {
            return sd_arr.value().previous_0(i);
        } else {
            return plain_bit_vec.value().previous_0(i);
        }
    }

    /**
     * @brief returns whether there is a one at index i
     * @param i [0..size-1]
     * @return whether there is a one at index i
     */
    inline bool operator[](uint_t i) {
        if (is_compressed()) {
            return sd_arr.value()[i];
        } else {
            return plain_bit_vec.value()[i];
        }
    }

    /**
     * @brief returns the size of the data structure in bytes
     * @return size of the data structure in bytes
     */
    inline uint64_t size_in_bytes() {
        if (!is_initialized()) return 0;

        if (is_compressed()) {
            return sd_arr.value().size_in_bytes();
        } else {
            return plain_bit_vec.value().size_in_bytes();
        }
    }

    /**
     * @brief serializes the hybrid_bit_vector to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) {
        bool is_init = is_initialized();
        out.write((char*)&is_init,1);
        if (!is_init) return;
        bool compressed = is_compressed();
        out.write((char*)&compressed,1);

        if (compressed ) {
            sd_arr.value().serialize(out);
        } else {
            plain_bit_vec.value().serialize(out);
        }
    }

    /**
     * @brief loads the hybrid_bit_vector from an input stream
     * @param in input stream
     */
    void load(std::istream& in) {
        bool is_init;
        in.read((char*)&is_init,1);
        if (!is_init) return;
        bool compressed;
        in.read((char*)&compressed,1);

        if (compressed) {
            sd_arr = std::move(sd_array<uint_t>());
            sd_arr.value().load(in);
        } else {
            plain_bit_vec = std::move(plain_bv_t());
            plain_bit_vec.value().load(in);
        }
    }
};