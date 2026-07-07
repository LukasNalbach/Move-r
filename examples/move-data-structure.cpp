/**
 * part of LukasNalbach/Move-r
 *
 * MIT License
 *
 * Copyright (c) Lukas Nalbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <move_data_structure/move_data_structure.hpp>
#include <misc/utils.hpp>
#include <misc/log.hpp>

int main()
{
    // Build a move data structure from the disjoint interval
    // sequence I = (0,1),(1,0) with n = 2
    move_data_structure<uint32_t, POS> mds({ { 0, 1 }, { 1, 0 } }, 2);

    // create a pair to perform move queries with
    std::pair<uint32_t, uint32_t> ix { 0, 0 };

    // perform some move queries
    std::cout << to_string<>(ix = mds.move(ix)) << std::endl;
    std::cout << to_string<>(ix = mds.move(ix)) << std::endl;

    // build a move data structure that additionally interleaves an extra row (here a char column, intended
    // e.g. for the characters of the bwt (sub-)runs) via the row_ts... template parameter; the extra-row
    // widths (here { 8 } bits) are passed after the construction parameters

    // use at most 4 threads and set a = 2
    move_data_structure<uint32_t, POS, char> mds_str({
        { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 }, { 4, 0 } }, 8,
        { .num_threads = 4, .a = 2 }, { 8 });

    // this disjoint interval sequence is not 2-balanced, because the output
    // interval [0,3] contains 4 >= 2a = 4 input intervals

    // print the pairs of the resulting disjoint interval sequence
    for (uint32_t i = 0; i < mds_str.num_intervals(); i++) {
        std::cout << to_string<>({ mds_str.p(i), mds_str.q(i) });
    }

    // the balancing algorithm has added the pair (6, 2) to balance the sequence
}