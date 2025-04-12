// Copyright (c) 2017, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 * rle_string.hpp
 *
 *  Created on: May 18, 2015
 *      Author: nicola
 *
 *  A run-length encoded string with rank/access functionalities.
 *
 *
 *  space of the structure: R * (H0 + log(n/R) + log(n/R)/B ) (1+o(1)) bits, n being text length,
 *  R number of runs, B block length, and H0 zero-order entropy of the run heads.
 *
 *  Time for all operations: O( B*(log(n/R)+H0) )
 *
 *  From the paper
 *
 *  Djamal Belazzougui, Fabio Cunial, Travis Gagie, Nicola Prezza and Mathieu Raffinot.
 *  Flexible Indexing of Repetitive Collections. Computability in Europe (CiE) 2017)
 *
 */

#pragma once

#include "huff_string.hpp"
#include "sparse_sd_vector.hpp"

namespace custom_r_index {

class rle_string {

protected:
    // block size: bitvector 'runs' has R/B bits set (R being number of runs)
    uint64_t B = 0;

    sparse_sd_vector runs;

    // for each letter, its runs stored contiguously
    std::vector<sparse_sd_vector> runs_per_letter;

    // store run heads in a compressed string supporting access/rank
    huff_string run_heads;

    // text length and number of runs
    uint64_t n = 0;
    uint64_t R = 0;

public:
    rle_string() { }

    /*
     * constructor: build structure on the input string
     * \param input the input string without 0x0 bytes in it.
     * \param B block size. The main sparse bitvector has R/B bits set (R being number of runs)
     *
     */
    rle_string(const std::string& input, uint64_t B = 2)
    {
        assert(not contains0(input));

        this->B = B;
        n = input.size();
        R = 0;
        auto runs_per_letter_bv = std::vector<std::vector<bool>>(256);

        // runs in main bitvector
        std::vector<bool> runs_bv;
        std::string run_heads_s;
        unsigned char last_c = input[0];

        for (uint64_t i = 1; i < input.size(); ++i) {
            if ((unsigned char) input[i] != last_c) {
                run_heads_s.push_back(last_c);
                runs_per_letter_bv[last_c].push_back(true);
                last_c = input[i];

                // push back a bit set only at the end of a block
                runs_bv.push_back(R % B == B - 1);
                R++;
            } else {
                runs_bv.push_back(false);
                runs_per_letter_bv[last_c].push_back(false);
            }
        }

        run_heads_s.push_back(last_c);
        runs_per_letter_bv[last_c].push_back(true);
        runs_bv.push_back(false);
        R++;

        assert(run_heads_s.size() == R);
        assert(R == count_runs(input));

        // std::cout << "runs in BWT(input) = " << count_runs(input) << std::endl;
        // std::cout << "runs in rle bwt = " << R << std::endl << std::endl;
        // now compact structures

        assert(runs_bv.size() == input.size());
        uint64_t t = 0;
        for (uint64_t i = 0; i < 256; ++i)
            t += runs_per_letter_bv[i].size();
        assert(t == input.size());
        runs = sparse_sd_vector(runs_bv);

        // a fast direct array: char -> bitvector.
        runs_per_letter = std::vector<sparse_sd_vector>(256);
		
        for (uint64_t i = 0; i < 256; ++i)
            runs_per_letter[i] = sparse_sd_vector(runs_per_letter_bv[i]);

        run_heads = huff_string(run_heads_s);
        assert(run_heads.size() == R);
    }

    unsigned char operator[](uint64_t i) const
    {
        assert(i < n);
        return run_heads[run_of(i).first];
    }

    /*
     * position of i-th character c. i starts from 0!
     */
    uint64_t select(uint64_t i, unsigned char c) const
    {
        assert(i < runs_per_letter[c].size());

        // i-th c is inside j-th c-run (j starts from 0)
        assert(i < runs_per_letter[c].size());
        uint64_t j = runs_per_letter[c].rank(i);

        // starting position of i-th c inside its run
        assert(j == 0 || i >= runs_per_letter[c].select(j - 1) + 1);
        uint64_t before = (j == 0 ? i : i - (runs_per_letter[c].select(j - 1) + 1));

        // position in run_heads
        uint64_t r = run_heads.select(j, c);

        // k = number of bits before position of interest in the main string
        // here, k is initialized looking at the sampled runs
        assert(r / B == 0 || r / B - 1 < runs.number_of_1());
        uint64_t k = (r / B == 0 ? 0 : runs.select(r / B - 1) + 1);

        // now add remaining run lengths to k
        for (uint64_t t = (r / B) * B; t < r; ++t) {
            k += run_at(t);
        }

        return k + before;
    }

    /*
     * number of c before position i
     */
    uint64_t rank(uint64_t i, unsigned char c) const
    {
        assert(i <= n);

        // letter does not exist in the text
        if (runs_per_letter[c].size() == 0)
            return 0;

        if (i == n)
            return runs_per_letter[c].size();

        uint64_t last_block = runs.rank(i);
        uint64_t current_run = last_block * B;

        // current position in the string: the first of a block
        uint64_t pos = 0;
        if (last_block > 0)
            pos = runs.select(last_block - 1) + 1;

        assert(pos <= i);

        uint64_t dist = i - pos;

        // otherwise, scan at most B runs
        while (pos < i) {
            pos += run_at(current_run);
            current_run++;

            if (pos <= i)
                dist = i - pos;
        }

        if (pos > i)
            current_run--;

        // position i is inside run current_run
        assert(current_run < R);

        // number of c runs before the current run
        uint64_t rk = run_heads.rank(current_run, c);

        // number of c before i in the current run
        uint64_t tail = (run_heads[current_run] == c) * dist;

        // in this case, either there are no c before position i
        // or the current run is the first of the kind ccc...cc
        if (rk == 0)
            return tail;

        return runs_per_letter[c].select(rk - 1) + 1 + tail;
    }

    /*
     * text position i is inside this run
     */
    uint64_t run_of_position(uint64_t i) const
    {
        assert(i < n);

        uint64_t last_block = runs.rank(i);
        uint64_t current_run = last_block * B;

        // current position in the string: the first of a block
        uint64_t pos = 0;
        if (last_block > 0)
            pos = runs.select(last_block - 1) + 1;

        assert(pos <= i);

        // otherwise, scan at most B runs
        while (pos < i) {
            pos += run_at(current_run);
            current_run++;
        }

        if (pos > i)
            current_run--;

        // position i is inside run current_run
        assert(current_run < R);

        return current_run;
    }

    // break range: given a range <l',r'> on the string and a character c, this function
    // breaks <l',r'> in maximal sub-ranges containing character c.
    // for simplicity and efficiency, we assume that characters at range extremities are both 'c'
    // thanks to the encoding (run-length), this function is quite efficient: O(|result|) ranks and selects
    std::vector<std::pair<uint64_t, uint64_t>> break_range(std::pair<uint64_t, uint64_t> rn, unsigned char c) const
    {
        auto l = rn.first;
        auto r = rn.second;

        assert(l <= r);
        assert(r < size());

        assert(operator[](l) == c);
        assert(operator[](r) == c);

        // retrieve runs that contain positions l and r
        auto run_l = run_of(l);
        auto run_r = run_of(r);

        // in this case rn contains only character c: do not break
        if (run_l.first == run_r.first)
            return { rn };

        std::vector<std::pair<uint64_t, uint64_t>> result;

        // first range: from l to the end of the run containing position l
        result.push_back({ l, run_l.second });

        // rank of c's of interest in run_heads
        uint64_t rank_l = run_heads.rank(run_l.first, c);
        uint64_t rank_r = run_heads.rank(run_r.first, c);

        // now retrieve run bounds of all c-runs of interest
        for (uint64_t j = rank_l + 1; j < rank_r; ++j) {
            result.push_back(run_range(run_heads.select(j, c)));
        }

        // now last (possibly incomplete) run
        auto range = run_range(run_heads.select(rank_r, c));
        result.push_back({ range.first, r });

        return result;
    }

    uint64_t size() const { return n; }

    /*
     * return inclusive range of j-th run in the string
     */
    std::pair<uint64_t, uint64_t> run_range(uint64_t j) const
    {
        assert(j < run_heads.size());

        uint64_t this_block = j / B;
        uint64_t current_run = this_block * B;
        uint64_t pos = (this_block == 0 ? 0 : runs.select(this_block - 1) + 1);

        while (current_run < j) {
            pos += run_at(current_run);
            current_run++;
        }

        assert(current_run == j);

        return { pos, pos + run_at(j) - 1 };
    }

    // length of i-th run
    uint64_t run_at(uint64_t i) const
    {
        assert(i < R);
        unsigned char c = run_heads[i];

        return runs_per_letter[c].gapAt(run_heads.rank(i, c));
    }

    uint64_t number_of_runs() const { return R; }

    uint64_t size_in_bytes() const
    {
        uint64_t bytes = sizeof(this) + runs.size_in_bytes() + run_heads.size_in_bytes();

        for (auto& vec : runs_per_letter) {
            bytes += vec.size_in_bytes();
        }

        return bytes;
    }

    /* serialize the structure to the ostream
     * \param out	 the ostream
     */
    void serialize(std::ostream& out) const
    {
        out.write((char*) &n, sizeof(n));
        out.write((char*) &R, sizeof(R));
        out.write((char*) &B, sizeof(B));

        if (n == 0)
            return;

        runs.serialize(out);

        for (uint64_t i = 0; i < 256; ++i)
            runs_per_letter[i].serialize(out);

        run_heads.serialize(out);
    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream& in)
    {
        in.read((char*) &n, sizeof(n));
        in.read((char*) &R, sizeof(R));
        in.read((char*) &B, sizeof(B));

        if (n == 0)
            return;

        runs.load(in);
        runs_per_letter = std::vector<sparse_sd_vector>(256);

        for (uint64_t i = 0; i < 256; ++i)
            runs_per_letter[i].load(in);

        run_heads.load(in);
    }

    std::string toString() const
    {
        std::string s;

        for (uint64_t i = 0; i < size(); ++i)
            s.push_back(operator[](i));

        return s;
    }

    /*
     * input: inclusive range rn, character c
     *
     * return the position j that is closest to rn.first,
     * such that character in position j is c and that is
     * adjacent to a position j' inside rn that contains a
     * character != c
     *
     * rn must contain c and at least another character d!=c
     *
     */
    uint64_t closest_run_break(std::pair<uint64_t, uint64_t> rn, unsigned char c) const
    {
        /*
         * case 1: range begins with a c-run: return last position of the run
         */
        if (operator[](rn.first) == c) {
            uint64_t i = run_of_position(rn.first);
            uint64_t j = run_range(i).second;

            // j must be inside rn, i.e. rn must not contain only c
            // j must not be last position of rn: this would imply
            // that rn contain only c
            assert(j < rn.second);

            return j;
        } else {
            // case 2: first c-run starts in the middle of the range

            // rank i of first c in the range
            uint64_t i = rank(rn.first, c);
            assert(i < rank(size(), c));

            // map from rank space to string position:
            // i now is the first position inside the range that contains c
            i = select(i, c);
            assert(operator[](i) == c);
            assert(i <= rn.second);

            return i;
        }
    }

private:
    static uint64_t count_runs(const std::string& s)
    {
        uint64_t runs = 1;

        for (uint64_t i = 1; i < s.size(); ++i) {
            if (s[i] != s[i - 1])
                runs++;
        }

        return runs;
    }

    //<j=run of position i, last position of j-th run>
    std::pair<uint64_t, uint64_t> run_of(uint64_t i) const
    {
        uint64_t last_block = runs.rank(i);
        uint64_t current_run = last_block * B;

        // current position in the string: the first of a block
        uint64_t pos = 0;
        if (last_block > 0)
            pos = runs.select(last_block - 1) + 1;

        assert(pos <= i);

        while (pos < i) {
            pos += run_at(current_run);
            current_run++;
        }

        assert(pos >= i);

        if (pos > i) {
            current_run--;
        } else { // pos==i
            pos += run_at(current_run);
        }

        assert(pos > 0);
        assert(current_run < R);

        return { current_run, pos - 1 };
    }

    static bool contains0(const std::string& s)
    {
        for (auto c : s)
            if (c == 0)
                return true;

        return false;
    }
};

}