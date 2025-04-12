// Copyright (c) 2017, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 * sparse_sd_vector: a wrapper on sd_vector of the sdsl library, with support for rank/select1
 */

//============================================================================

#pragma once

#include <vector>
#include <sdsl/bit_vectors.hpp>

namespace custom_r_index {

class sparse_sd_vector {

protected:
    // bitvector length
    uint64_t u = 0;

    sdsl::sd_vector<> sdv;
    sdsl::sd_vector<>::rank_1_type rank1;
    sdsl::sd_vector<>::select_1_type select1;

public:
    /*
     * empty constructor. Initialize bitvector with length 0.
     */
    sparse_sd_vector() { }

    /*
     * constructor. build bitvector given a vector of bools
     */
    sparse_sd_vector(std::vector<bool>& b)
    {
        if (b.size() == 0)
            return;

        u = b.size();
        sdsl::bit_vector bv(b.size());

        for (uint64_t i = 0; i < b.size(); ++i)
            bv[i] = b[i];

        sdv = sdsl::sd_vector<>(bv);
        rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        select1 = sdsl::sd_vector<>::select_1_type(&sdv);
    }

    /*
     * constructor. build bitvector given a sdsl::bit_vector
     */
    sparse_sd_vector(sdsl::bit_vector& bv)
    {
        sdv = sdsl::sd_vector<>(bv);
        rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        select1 = sdsl::sd_vector<>::select_1_type(&sdv);
    }

    sparse_sd_vector& operator=(const sparse_sd_vector& other)
    {

        u = other.sdv.size();
        sdv = sdsl::sd_vector<>(other.sdv);
        rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        select1 = sdsl::sd_vector<>::select_1_type(&sdv);

        return *this;
    }

    /*
     * argument: position i in the bitvector
     * returns: bit in position i
     * only access! the bitvector is static.
     */
    bool operator[](uint64_t i) const
    {
        assert(i < size());
        return sdv[i];
    }

    bool at(uint64_t i) const
    {
        return operator[](i);
    }

    /*
     * argument: position i in the bitvector
     * returns: number of bits equal to 1 before position i excluded
     */
    uint64_t rank(uint64_t i) const
    {
        assert(i <= size());
        return rank1(i);
    }

    /*
     * input: position 0<=i<=n
     * output: predecessor of i (position i excluded) in
     * rank space (apply select to get bitvector space)
     */
    uint64_t predecessor_rank(uint64_t i) const
    {
        /*
         * i must have a predecessor
         */
        assert(rank(i) > 0);

        return rank(i) - 1;
    }

    /*
     * input: position 0<=i<=n
     * output: predecessor of i (i excluded) in
     * bitvector space
     */
    uint64_t predecessor(uint64_t i) const
    {
        /*
         * i must have a predecessor
         */
        assert(rank(i) > 0);

        return select(rank(i) - 1);
    }

    /*
     * input: position 0<=i<=n
     * output: rank of predecessor of i (i excluded) in
     * bitvector space. If i does not have a predecessor,
     * return rank of the last bit set in the bitvector
     */
    uint64_t predecessor_rank_circular(uint64_t i) const
    {
        return rank(i) == 0 ? number_of_1() - 1 : rank(i) - 1;
    }

    /*
     * retrieve length of the i-th gap (i>=0). gap length includes the leading 1
     * \param i<number_of_1()
     *
     */
    uint64_t gapAt(uint64_t i) const
    {
        assert(i < number_of_1());

        if (i == 0)
            return select(0) + 1;

        return select(i) - select(i - 1);
    }

    /*
     * argument: integer i>0
     * returns: position of the i-th one in the bitvector. i starts from 0!
     */
    uint64_t select(uint64_t i) const
    {
        assert(i < number_of_1());
        return select1(i + 1); // in sd_vector, i starts from 1
    }

    /*
     * returns: size of the bitvector
     */
    uint64_t size() const { return u; }

    /*
     * returns: number of 1s in the bitvector
     */
    uint64_t number_of_1() const { return rank1(size()); }

    uint64_t size_in_bytes() const
    {
        return sdsl::size_in_bytes(sdv);
    }

    /* serialize the structure to the ostream
     * \param out	 the ostream
     */
    void serialize(std::ostream& out) const
    {
        out.write((char*) &u, sizeof(u));

        if (u == 0)
            return;

        sdv.serialize(out);
    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream& in)
    {
        in.read((char*) &u, sizeof(u));

        if (u == 0)
            return;

        sdv.load(in);
        rank1 = sdsl::sd_vector<>::rank_1_type(&sdv);
        select1 = sdsl::sd_vector<>::select_1_type(&sdv);
    }
};

}