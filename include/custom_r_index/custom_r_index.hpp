// Copyright (c) 2017, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 * r_index.hpp
 *
 *  Created on: Apr 13, 2017
 *      Author: nico
 *
 * Small version of the r-index: O(r) words of space, O(log(n/r)) locate time per occurrence
 *
 */

#pragma once

#include "rle_string.hpp"
#include "sparse_sd_vector.hpp"
#include <misc/utils.hpp>

namespace custom_r_index {

/*
 * sparse RLBWT: r (log sigma + (1+epsilon) * log (n/r)) (1+o(1)) bits
 */
class index {

protected:
    // F column of the BWT (vector of 256 elements)
    std::vector<uint64_t> F;

    // L column of the BWT, run-length compressed
    rle_string bwt;

    // vector containing the sufffix array samples at the run start/end positions
    sdsl::int_vector<> samples;

    uint64_t r = 0; // number of BWT runs

public:
    index() = default;

    uint64_t size_in_bytes() const
    {
        return sizeof(this) +
            sdsl::size_in_bytes(samples) +
            bwt.size_in_bytes() +
            sizeof(uint64_t) * F.size(); // F
    }

    /*
     * Build index
     */
    template <typename int_t>
    index(const std::string& BWT, const std::vector<int_t>& SA, bool sa_samples = false)
    {
        bwt = rle_string(BWT);
        uint64_t n = BWT.size();

        // build F column
        F = std::vector<uint64_t>(256, 0);
        for (unsigned char c : BWT)
            F[c]++;

        for (uint64_t i = 255; i > 0; --i)
            F[i] = F[i - 1];

        F[0] = 0;

        for (uint64_t i = 1; i < 256; ++i)
            F[i] += F[i - 1];

        r = bwt.number_of_runs();

        if (!sa_samples) return;

        samples.width(std::bit_width(n));
        samples.resize(r);
        uint64_t run = 0;

        for (uint64_t i = 1; i < n; i++) {
            if (BWT[i] != BWT[i - 1]) {
                samples[run] = SA[i - 1];
                run++;
            }
        }

        samples[r - 1] = SA[n - 1];
    }

    uint64_t num_bwt_runs() const
    {
        return bwt.number_of_runs();
    }

    /*
     * get full BWT range
     */
    std::pair<uint64_t, uint64_t> full_range() const
    {
        // inclusive range
        return { 0, bwt_size() - 1 };
    }

    unsigned char operator[](uint64_t i) const
    {
        return bwt[i];
    }

    uint64_t sample(uint64_t i) const
    {
        return samples[i];
    }

    uint64_t sample_pos(uint64_t i) const
    {
        return bwt.run_range(i).second;
    }

    /*
     * \param r inclusive range of a string w
     * \param c character
     * \return inclusive range of cw
     */
    std::pair<uint64_t, uint64_t> LF(std::pair<uint64_t, uint64_t> rn, unsigned char c) const
    {
        // if character does not appear in the text, return empty pair
        if ((c == 255 and F[c] == bwt_size()) || F[c] >= F[c + 1])
            return { 1, 0 };

        // number of c before the interval
        uint64_t c_before = bwt.rank(rn.first, c);

        // number of c inside the interval rn
        uint64_t c_inside = bwt.rank(rn.second + 1, c) - c_before;

        // if there are no c in the interval, return empty range
        if (c_inside == 0)
            return { 1, 0 };

        uint64_t l = F[c] + c_before;

        return { l, l + c_inside - 1 };
    }

    // backward navigation of the BWT
    uint64_t LF(uint64_t i) const
    {
        auto c = bwt[i];
        return F[c] + bwt.rank(i, c);
    }

    // forward navigation of the BWT
    uint64_t FL(uint64_t i) const
    {
        // i-th character in first BWT column
        auto c = F_at(i);

        // this c is the j-th (counting from 0)
        uint64_t j = i - F[c];

        return bwt.select(j, (unsigned char) c);
    }

    // forward navigation of the BWT, where for efficiency we give c=F[i] as input
    uint64_t FL(uint64_t i, unsigned char c) const
    {
        // i-th character in first BWT column
        assert(c == F_at(i));

        // this c is the j-th (counting from 0)
        uint64_t j = i - F[c];

        return bwt.select(j, (unsigned char) c);
    }

    /*
     * access column F at position i
     */
    unsigned char F_at(uint64_t i) const
    {
        uint64_t c = (upper_bound(F.begin(), F.end(), i) - F.begin()) - 1;
        assert(c < 256);
        assert(i >= F[c]);

        return (unsigned char) c;
    }

    /*
     * Return BWT range of character c
     */
    std::pair<uint64_t, uint64_t> get_char_range(unsigned char c) const
    {
        // if character does not appear in the text, return empty pair
        if ((c == 255 and F[c] == bwt_size()) || F[c] >= F[c + 1])
            return { 1, 0 };

        uint64_t l = F[c];
        uint64_t r = bwt_size() - 1;

        if (c < 255)
			r = F[c + 1] - 1;

        return { l, r };
    }

    /*
     * Return BWT range of pattern P
     */
    std::pair<uint64_t, uint64_t> count(const std::string& P) const
    {
        auto range = full_range();
        uint64_t m = P.size();

        for (uint64_t i = 0; i < m and range.second >= range.first; ++i)
            range = LF(range, P[m - i - 1]);

        return range;
    }

    /*
     * Return number of occurrences of P in the text
     */
    uint64_t occ(std::string& P) const
    {
        auto rn = count(P);
        return rn.second >= rn.first ? (rn.second - rn.first) + 1 : 0;
    }

    /*
     * get number of runs in the BWT (terminator character included)
     */
    uint64_t number_of_runs() const
    {
        return bwt.number_of_runs();
    }

    /*
     * get BWT in string format
     */
    const std::string get_bwt() const
    {
        return bwt.toString();
    }

    /* serialize the structure to the ostream
     * \param out	 the ostream
     */
    void serialize(std::ostream& out) const
    {
        assert(F.size() > 0);
        assert(bwt.size() > 0);

        out.write((char*)F.data(), 256 * sizeof(uint64_t));
        bwt.serialize(out);

        samples.serialize(out);
    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream& in)
    {
        F = std::vector<uint64_t>(256);
        in.read((char*)F.data(), 256 * sizeof(uint64_t));

        bwt.load(in);
        r = bwt.number_of_runs();

        samples.load(in);
    }

    /*
     * save the structure to the path specified.
     * \param path_prefix prefix of the index files. suffix ".ri" will be automatically added
     */
    void save_to_file(std::string path_prefix) const
    {
        std::string path = std::string(path_prefix).append(".ri");

        std::ofstream out(path);
        serialize(out);
        out.close();
    }

    /*
     * load the structure from the path specified.
     * \param path: full file name
     */
    void load_from_file(std::string path)
    {
        std::ifstream in(path);
        load(in);
        in.close();
    }

    uint64_t text_size() const
    {
        return bwt.size() - 1;
    }

    uint64_t bwt_size() const
    {
        return bwt.size();
    }

    /*
     * returns <<l,r>, SA[l/r] >, where l,r are the inclusive ranges of the pattern P. If P does not occur, then l>r
     *
     * returns <range, j,k>
     *
     */    
    std::tuple<uint64_t, uint64_t, uint64_t> count_and_get_occ(const std::string& P) const
    {
        assert(samples.size() != 0);
        std::pair<uint64_t, uint64_t> range = full_range();

        // k = SA[r]
        uint64_t k = samples[r - 1];

        std::pair<uint64_t, uint64_t> range1;
        uint64_t m = P.size();

        for (uint64_t i = 0; i < m and range.second >= range.first; ++i) {
            unsigned char c = P[m - i - 1];
            range1 = LF(range, c);

            // if suffix can be left-extended with char
            if (range1.first <= range1.second) {
                if (bwt[range.second] == c) {
                    // last c is at the start&end of range. Then, we have this sample by induction!
                    assert(k > 0);
                    k--;
                } else {
                    // find last c in range (there must be one because range1 is not empty)
                    // and get its sample (must be sampled because it is at the start/end of a run)
                    uint64_t rnk = bwt.rank(range.second, c);

                    // this is the rank of the first/last c
                    rnk--;

                    // jump to the corresponding BWT position
                    uint64_t j = bwt.select(rnk, c);

                    // the c must be in the range
                    assert(j >= range.first and j <= range.second);

                    // run of position j
                    uint64_t run_of_j = bwt.run_of_position(j);

                    k = samples[run_of_j];

                    if (k == 0) {
                        k = bwt.size() - 1;
                    } else {
                        k--;
                    }
                }
            }

            range = range1;
        }

        return { range.first, range.second, k };
    }
};

}