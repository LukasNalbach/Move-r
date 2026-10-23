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


// naive references of the index tests: suffix array intervals

#ifndef MOVE_R_TESTS_TEST_INDEX_REFERENCES_HPP
#define MOVE_R_TESTS_TEST_INDEX_REFERENCES_HPP

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>


/**
 * @brief returns the suffix array interval [b, e] of the pattern P in text; b > e <=> P does not occur in text
 * @tparam sym_t symbol type of the text and of the pattern
 * @param SA [0..n] the suffix array of text with the sentinel (n = |text|)
 * @param text the text, in the alphabet the suffix array is sorted by
 * @param P the pattern
 */
template <typename sym_t>
inline std::pair<int64_t, int64_t> reference_sa_interval(
    const std::vector<int64_t>& SA, const std::vector<sym_t>& text, const std::vector<sym_t>& P)
{
    const int64_t n = text.size();

    // the suffix of text starting at sa (terminated by the sentinel) is lexicographically smaller than P
    auto suffix_less = [&](int64_t sa, const std::vector<sym_t>& P) {
        for (uint64_t i = 0; i < P.size(); i++) {
            if (sa + (int64_t) i >= n) return true; // the sentinel is smaller than every symbol of P
            if (text[sa + i] != P[i]) return text[sa + i] < P[i];
        }

        return false;
    };

    // P is lexicographically smaller than the suffix of text starting at sa, and not a prefix of it
    auto pattern_less = [&](const std::vector<sym_t>& P, int64_t sa) {
        for (uint64_t i = 0; i < P.size(); i++) {
            if (sa + (int64_t) i >= n) return false;
            if (text[sa + i] != P[i]) return P[i] < text[sa + i];
        }

        return false;
    };

    return {
        std::lower_bound(SA.begin(), SA.end(), P, suffix_less) - SA.begin(),
        std::upper_bound(SA.begin(), SA.end(), P, pattern_less) - SA.begin() - 1
    };
}

#endif // MOVE_R_TESTS_TEST_INDEX_REFERENCES_HPP
