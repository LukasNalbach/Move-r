/**
 * part of LukasNalbach/Move-r
 *
 * MIT License
 *
 * Copyright (c) Patrick Dinklage, Jan Zumbrink, Lukas Nalbach
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

#pragma once

#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <memory>
#include <vector>

#include <libsais.h>
#include <libsais64.h>
#include <ips4o.hpp>

#include <misc/utils.hpp>
#include <ordered/btree/map.hpp>
#include <rmq/rmq.hpp>

template <typename int_t = int32_t>
struct lzend_phr_t {
    int_t lnk; // phrase id, where the source ends (a new phrase can extend multiple phrases)
    int_t len; // len of phrase (including the extension)
    int_t ext; // extension
};

template <typename int_t = int32_t>
struct input_val_idx_pair_t { // used for constructing the transformed reverse delta suffix array
    int_t value;
    int_t index;

    bool operator<(const input_val_idx_pair_t& other) const
    {
        return value < other.value;
    }
};

template <typename int_t = int32_t>
std::vector<lzend_phr_t<int_t>> construct_lzend_of_reverse(std::vector<int_t>& dsa, int_t h = -1, bool log = false)
{
    auto time = now();
    int_t n = dsa.size();

    // Now reverse the given DSA and transform every value into a non-negative value, such that
    // the order of the values is preserved, but the distances between the values are reduced to one.
    // This transformation makes the construction of the suffix array of the reversed DSA faster,
    // because the suffix array construction time rises linearly with the size of the given alphabet.
    if (log) std::cout << "reversing and compactifying alphabet of DSA" << std::flush;

    std::vector<input_val_idx_pair_t<int_t>> value_index_pairs;
    value_index_pairs.reserve(n);

    for (int_t i = 0; i < n; i++) {
        value_index_pairs.emplace_back(input_val_idx_pair_t<int_t> { dsa[i], i });
    }

    // sort all value index pairs by the original DSA value
    ips4o::parallel::sort(value_index_pairs.begin(), value_index_pairs.end());
    
    // use this sorted vector of value index pairs to construct the transformed reverse DSA
    std::vector<int_t> rev_dsa;
    no_init_resize(rev_dsa, n);

    int_t last_value = value_index_pairs[0].value;
    int_t new_value = 0;
    rev_dsa[n - value_index_pairs[0].index - 1] = new_value;

    for (int_t i = 1; i < n; i++) {
        if (value_index_pairs[i].value != last_value) {
            // if the current value is unequal to the last value, the reversed DSA should also have an unequal value
            // if the current value is equal to the last value, the rev_dsa should have the same value for both indices
            last_value = value_index_pairs[i].value;
            new_value++;
        }

        rev_dsa[n - value_index_pairs[i].index - 1] = new_value;
    }

    // discard value_index_pairs
    value_index_pairs.clear();
    value_index_pairs.shrink_to_fit();

    // return the new alphabet size, because it is important to know for the sa construction
    int_t alphabet_size = new_value + 1;

    if (log) time = log_runtime(time);

    // construct suffix array of the reversed DSA
    if (log) std::cout << "building SA of rev(DSA)" << std::flush;
    std::vector<int_t> sa;
    no_init_resize(sa, n);
    
    if constexpr (std::is_same_v<int_t, int32_t>) {
        libsais_int(rev_dsa.data(), sa.data(), n, alphabet_size, 0);
    } else {
        libsais64_long(rev_dsa.data(), sa.data(), n, alphabet_size, 0);
    }

    if (log) time = log_runtime(time);

    // construct PLCP of rev(DSA)
    if (log) std::cout << "building PLCP of rev(DSA)" << std::flush;
    std::vector<int_t> isa;
    no_init_resize(isa, n);
    auto& plcp = isa;
    
    if constexpr (std::is_same_v<int_t, int32_t>) {
        libsais_plcp_int(rev_dsa.data(), sa.data(), plcp.data(), n);
    } else {
        libsais64_plcp_int(rev_dsa.data(), sa.data(), plcp.data(), n);
    }

    if (log) time = log_runtime(time);

    // discard rev(DSA)
    rev_dsa.clear();
    rev_dsa.shrink_to_fit();

    // construct LCP of rev(DSA)
    if (log) std::cout << "building LCP of rev(DSA)" << std::flush;
    std::vector<int_t> lcp;
    no_init_resize(lcp, n);

    if constexpr (std::is_same_v<int_t, int32_t>) {
        libsais_lcp(plcp.data(), sa.data(), lcp.data(), n);
    } else {
        libsais64_lcp(plcp.data(), sa.data(), lcp.data(), n);
    }

    if (log) time = log_runtime(time);

    // construct RMQ data structure
    if (log) std::cout << "building RMQ of LCP" << std::flush;
    rmq::RMQ<int_t> rmq(lcp.data(), n);
    if (log) time = log_runtime(time);

    // construct permuted ISA (PISA) of rev(DSA)
    if (log) std::cout << "building permuted ISA (PISA) of rev(DSA)" << std::flush;
    auto& pisa = isa;

    for (int_t i = 0; i < n; i++) {
        pisa[n - sa[i] - 1] = i;
    }

    if (log) time = log_runtime(time);

    // discard SA of rev(DSA)
    sa.clear();
    sa.shrink_to_fit();

    // construct_lzend_of_reverse
    if (log) std::cout << "building LZ-End parsing of rev(DSA)" << std::flush;

    // initialize predecessor/successor
    ordered::btree::Map<int_t, int_t> marked;

    // helpers
    struct Candidate {
        int_t lex_pos;
        int_t lnk;
        int_t len;
    };

    auto lex_smaller_phrase = [&](int_t const x) {
        auto const r = marked.predecessor(x - 1);
        return r.exists
            ? Candidate { r.key, r.value, lcp[rmq(r.key + 1, x)] }
            : Candidate { 0, 0, 0 };
    };

    auto lex_greater_phrase = [&](int_t const x) {
        auto const r = marked.successor(x + 1);
        return r.exists
            ? Candidate { r.key, r.value, lcp[rmq(x + 1, r.key)] }
            : Candidate { 0, 0, 0 };
    };

    std::vector<lzend_phr_t<int_t>> parsing;
    parsing.push_back({ 0, 1, dsa[0] }); // initial empty phrase
    int_t z = 0; // index of latest phrase (number of phrases would be z+1)

    for (int_t i = 1; i < n; i++) {
        int_t const len1 = parsing[z].len;
        int_t const len2 = len1 + (z > 0 ? parsing[z - 1].len : 0);
        int_t const isa_last = isa[i - 1];

        // find source phrase candidates
        int_t p1 = -1, p2 = -1;

        auto find_copy_source = [&](std::function<Candidate(int_t)> f) {
            auto c = f(isa_last);

            if (c.len >= len1) {
                p1 = c.lnk; // only mark the small enough phrases omits the phrase length check here
                
                // only check for merge candidates if the new phrase would be small enough (or the max phrase size is unbounded)
                if (i > len1 && (h == -1 || len2 < h)) { 
                    if (c.lnk == z - 1) c = f(c.lex_pos);
                    if (c.len >= len2) p2 = c.lnk;
                }
            }
        };

        if (h == -1 || len1 < h) {
            find_copy_source(lex_smaller_phrase);
            
            if (p1 == -1 || p2 == -1) {
                find_copy_source(lex_greater_phrase);
            }
        }

        // case distinction according to Lemma 1
        if (p2 != -1) {
            // merge last two phrases
            marked.erase(isa[i - 1 - len1]);
            parsing.pop_back();
            --z;
            parsing.back() = lzend_phr_t<int_t> { p2, len2 + 1, dsa[i] };
        } else if (p1 != -1 && (h == -1 || len1 < h)) {
            // extend last phrase
            parsing.back() = lzend_phr_t<int_t> { p1, len1 + 1, dsa[i] };
        } else {
            // lazily mark previous phrase
            marked.insert(isa_last, z);

            // begin new phrase
            parsing.push_back(lzend_phr_t<int_t> { 0, 1, dsa[i] });
            ++z;
        }
    }

    if (log) time = log_runtime(time);

    return parsing;
}