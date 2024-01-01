#ifndef BWT_SCAN_H
#define BWT_SCAN_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include "sdsl/wavelet_trees.hpp"

using std::cerr;
using std::cout;
using std::endl;

/* note: pfbwt uses 0x0 for end of string, 0x1 for end of word */
/* r-index uses 1 for end of string (we can change this, nbd?) */

/* scans bwt and returns suffix array samples (run-starts & run-ends), F array, and terminator position
* NOTE: any char c s.t. c < terminator is treated as terminator (ie. if c==0 and terminator == 1, then c is treated as 1) 
*
* TODO: make this less ugly.
*/
std::tuple<uint64_t, uint64_t> bwt_scan_ssa(std::string bwt_fname, std::vector<std::pair<uint64_t, uint64_t>>& samples_first, std::vector<uint64_t>& samples_last, std::vector<uint64_t>& F, uint64_t* term_pos=nullptr) {
    static const uint8_t TERMINATOR = 1; // assuming this is some very low number
    
    sdsl::wt_huff<> bwt;
    sdsl::construct(bwt, bwt_fname, 1);
    uint64_t n = bwt.size();
    uint64_t r = 0;

    sdsl::bit_vector run_starts(n,0);
    sdsl::bit_vector run_ends(n, 0);
    F = std::vector<ulint>(256,0);
    // populate F array
    uint8_t lastchar=bwt[0], nextchar=bwt[0], nextnextchar=bwt[1];
    uint64_t i, j, k;
    for (i = 0; i < bwt.size()-1; ++i) {
        nextnextchar = bwt[i+1];
        if (nextchar <= TERMINATOR) {
            F[TERMINATOR] += 1;
            if (term_pos) *term_pos = i;
        } else 
            F[nextchar] += 1;
        if (nextchar != lastchar || i == 0) { // i is runstart or first char of bwt
            run_starts[i] = 1;
            r += 1;
        }
        if (nextchar != nextnextchar) 
            run_ends[i] = 1;
        lastchar = nextchar;
        nextchar = nextnextchar;
    }
    // deal with i==n-1 case (nextchar == bwt[n-2+1] == bwt[n-1], lastchar == bwt[n-2])
    if (nextchar <= TERMINATOR) 
        F[TERMINATOR] += 1; 
    else 
        F[nextchar] += 1; 
    if (nextchar != lastchar || i == 0)  {
        r += 1;
        run_starts[i] = 1;
    }
    run_ends[i] = 1; // i==n-1 aka last character in bwt
    // make sure only one occurence of terminator exists
    if (F[TERMINATOR] > 1) {
        cerr << "error: bwt can only have one occurence of TERMINATOR! (found " << F[TERMINATOR] << " occurences" << endl;
        exit(1);
    }
    // properly finish populating F
    uint64_t running_sum = 1;
    F[TERMINATOR] = 0;
    for (i = TERMINATOR+1; i < F.size(); i++) {
        j = F[i];
        F[i] = running_sum;
        running_sum += j;
    }

    // rank-ify run_starts and run_ends
    sdsl::bit_vector::rank_1_type rs_rank(&run_starts);
    sdsl::bit_vector::rank_1_type re_rank(&run_ends);
    assert(run_starts.size() == run_ends.size());
    assert(sdsl::util::cnt_one_bits(run_starts) == sdsl::util::cnt_one_bits(run_ends));
    uint64_t nruns = sdsl::util::cnt_one_bits(run_starts);
    samples_first.resize(nruns);
    samples_last.assign(nruns, 0);
    k = 0; // k is bwt position
    j = n-1; // j is sa sample
    uint8_t c;
    // find suffix array samples using LF mapping and induction
    // S[LF(k)] = S[k]-1
    for (i = 0; i < n; ++i) {
        if (run_starts[k]) {
            auto idx = rs_rank(k);
            samples_first[idx] = std::pair<uint64_t, uint64_t>(j ? j-1 : n-1, idx); 
        } 
        if (run_ends[k]) {
            auto idx = re_rank(k);
            samples_last[idx] = j ? j-1 : n-1;
        }
        // prepare for next iteration
        j--;
        k = F[bwt[k]] + bwt.rank(k, bwt[k]); // LF mapping, ie. next bwt position
    }

    /* DEBUG:
    std::ofstream saf("samples_first_bwtscan"); std::ofstream sal("samples_last_bwtscan");
    for (auto i: samples_first) saf << std::get<0>(i) << " " << std::get<1>(i) << endl;
    for (auto i : samples_last) sal << i << endl;
    saf.close(); sal.close();
    std::ofstream ffs("F");
    for (auto i = 0; i < F.size(); ++i) {
        ffs << i << " " << F[i] << endl;
    }
    */
    return std::tuple<uint64_t,uint64_t>(n,r);
}

#endif
