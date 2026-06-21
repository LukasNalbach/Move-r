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

#include <move_rb/move_rb.hpp>

int main()
{
    // build a bi-directional move-r index (move_rb); in addition to count and
    // locate, it supports extending the currently matched pattern to both the
    // left and the right, as well as approximate pattern matching
    move_rb<> index("This is a test string");

    // build a bi-directional index that uses a relative lempel-ziv encoded
    // differential suffix array (rlzsa) for the locate support
    move_rb<_locate_rlzsa> index_rlzsa("This is a test string");

    // print the index size
    std::cout << format_size(index.size_in_bytes()) << std::endl;

    // ############################# BI-DIRECTIONAL SEARCH #############################

    // get an empty (locate-) search context for the index
    auto ctx = index.empty_context<LOCATE>();

    // search the pattern "test" by extending it in both directions: start with
    // "e", extend to the left to "te", then to the right to "tes" and "test";
    // extend(...) returns the extended context and whether the extended pattern
    // occurs in the input
    auto [c1, ok1] = ctx.extend('e', LEFT);  // "e"
    auto [c2, ok2] = c1.extend('t', LEFT);   // "te"
    auto [c3, ok3] = c2.extend('s', RIGHT);  // "tes"
    auto [c4, ok4] = c3.extend('t', RIGHT);  // "test"
    ctx = c4;

    // print the number of occurrences of "test"
    std::cout << ctx.num_occ() << std::endl;

    // print the suffix-array interval [b,e] of "test" in the forward text
    auto [b, e] = ctx.forward_sa_interval();
    std::cout << "b = " << b << ", e = " << e << std::endl;

    // enter the locate phase for the matched pattern and store all of its
    // occurrences in a vector
    auto loc_ctx = ctx.locate_phase();
    auto Occ = loc_ctx.locate();
    for (auto o : Occ) std::cout << o << ", ";
    std::cout << std::endl;

    // if only the number of occurrences is needed, a count-context avoids
    // maintaining the information required for locating
    auto count_ctx = index.empty_context<COUNT>();
    auto [d1, found] = count_ctx.extend('i', RIGHT); // "i"
    auto [d2, _]     = d1.extend('s', RIGHT);        // "is"
    std::cout << d2.num_occ() << std::endl;

    // ############################# APPROXIMATE PATTERN MATCHING #############################

    // count the occurrences of "test" with at most 1 mismatch (w.r.t. the
    // hamming distance), using a pigeon-hole search scheme
    search_scheme_t scheme = pigeon_hole_scheme(1);
    std::cout << index.count_hamming_dist("test", scheme) << std::endl;

    // locate the occurrences of "best" with at most 1 error (w.r.t. the edit
    // distance); each occurrence stores its position, length and error count
    std::vector<aprx_occ_t<uint32_t>> aprx_Occ =
        index.locate<EDIT_DISTANCE>("best", min_u_scheme(1));
    for (auto occ : aprx_Occ) {
        std::cout << "(pos=" << occ.pos << " len=" << occ.len
                  << " err=" << occ.err << ") ";
    }
    std::cout << std::endl;
}
