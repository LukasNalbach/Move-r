#pragma once

#include <vector>
#include <iostream>
#include <optional>
#include "sd_array.hpp"
#include "plain_bit_vector.hpp"

/**
 * @brief hybrid bit vector (either an sd_array or a plain bit vector)
 * @tparam pos_t unsigned integer type
 */
template <typename pos_t = uint32_t, bool build_rank_support = false, bool build_select_0_support = false, bool build_select_1_support = false>
class hybrid_bit_vector {
    static_assert(std::is_same_v<pos_t,uint32_t> || std::is_same_v<pos_t,uint64_t>);

    public:
    using plain_bv_t = plain_bit_vector<pos_t,build_rank_support,build_select_0_support,build_select_1_support>;

    /* if there are less than (compression_threshold * size of the bitvector)
     * ones in the bit vector, then it should be compressed */
    static constexpr double compression_threshold = 0.1;

    protected:
    std::optional<sd_array<pos_t>> sd_arr; // the sd_array
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
            sd_arr = std::move(sd_array<pos_t>(plain_bit_vec));
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
            sd_arr = std::move(sd_array<pos_t>(plain_bit_vec));
        } else {
            this->plain_bit_vec = std::move(plain_bv_t(std::move(plain_bit_vec)));
        }
    }

    /**
     * @brief constructs a new hybrid_bit_vector from an sd_vector
     * @param sd_vec an sd_vector
     */
    hybrid_bit_vector(const sdsl::sd_vector<>& sd_vec) {
        sd_arr = std::move(sd_array<pos_t>(sd_vec));
    }

    /**
     * @brief constructs a new hybrid_bit_vector from an sd_vector
     * @param sd_vec an sd_vector
     */
    hybrid_bit_vector(sdsl::sd_vector<>&& sd_vec) {
        sd_arr = std::move(sd_array<pos_t>(std::move(sd_vec)));
    }

    inline bool is_initialized() const {
        return sd_arr.has_value() || plain_bit_vec.has_value();
    }

    /**
     * @brief returns, whether the bit vector is compressed
     * @return whether the bit vector is compressed
     */
    inline bool is_compressed() const {
        return sd_arr.has_value();
    }

    /**
     * @brief returns the size of the bit vector
     * @return the size of the bit vector 
     */
    inline pos_t size() const {
        if (is_compressed()) {
            return sd_arr.value().size();
        } else {
            return plain_bit_vec.value().size();
        }
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
    inline pos_t num_zeros() const {
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
    inline pos_t rank_1(pos_t i) const {
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
    inline pos_t select_1(pos_t i) const {
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
    inline pos_t rank_0(pos_t i) const {
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
    inline pos_t select_0(pos_t i) const {
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
    inline pos_t next_1(pos_t i) const {
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
    inline pos_t previous_1(pos_t i) const {
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
    inline pos_t next_0(pos_t i) const {
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
    inline pos_t previous_0(pos_t i) const {
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
    inline bool operator[](pos_t i) const {
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
    inline uint64_t size_in_bytes() const {
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
    void serialize(std::ostream& out) const {
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
            sd_arr = std::move(sd_array<pos_t>());
            sd_arr.value().load(in);
        } else {
            plain_bit_vec = std::move(plain_bv_t());
            plain_bit_vec.value().load(in);
        }
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