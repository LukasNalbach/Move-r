// Copyright (c) 2017, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 * huff_string.hpp
 *
 *  Created on: May 18, 2015
 *      Author: nicola
 *
 *  Huffman-compressed string with access/rank/select. The class is a wrapper on sdsl::wt_huff, with a simpler constructor
 */

#pragma once

#include <misc/utils.hpp>
#include <sdsl/wavelet_trees.hpp>

namespace custom_r_index {

class huff_string {

protected:
    sdsl::wt_huff<> wt;

public:
    huff_string() { }

    huff_string(std::string& s)
    {
        s.push_back(0);
        construct_im(wt, s.c_str(), 1);

        assert(wt.size() == s.size() - 1);
    }

    unsigned char operator[](uint64_t i) const
    {
        assert(i < wt.size());
        return wt[i];
    }

    uint64_t size() const
    {
        return wt.size();
    }

    uint64_t rank(uint64_t i, unsigned char c) const
    {
        assert(i <= wt.size());
        return wt.rank(i, c);
    }

    /*
     * position of i-th character c. i starts from 0!
     */
    uint64_t select(uint64_t i, unsigned char c) const
    {
        return wt.select(i + 1, c);
    }

    uint64_t size_in_bytes() const
    {
        return sizeof(this) + sdsl::size_in_bytes(wt);
    }

    /* serialize the structure to the ostream
     * \param out	 the ostream
     */
    uint64_t serialize(std::ostream& out) const
    {
        return wt.serialize(out);
    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream& in)
    {
        wt.load(in);
    }
};

}