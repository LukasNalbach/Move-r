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

#ifndef MOVE_R_DATA_STRUCTURES_SYMBOL_GROUPS_HPP
#define MOVE_R_DATA_STRUCTURES_SYMBOL_GROUPS_HPP

#include <algorithm>
#include <bit>
#include <cstdint>
#include <istream>
#include <limits>
#include <omp.h>
#include <ostream>
#include <vector>

#include <data_structures/rank_select_set.hpp>
#include <misc/utils.hpp>

/**
 * @brief for some symbols of a sequence, the groups of 2^k consecutive positions that contain them: for each of these
 *        symbols a rank_select_set of its groups, which gives the next and the previous group of a position that
 *        contains the symbol
 */
class symbol_groups {
public:
    // ############################# TYPES AND PARAMETERS #############################

    // returned by next and previous if there is no such group
    static constexpr uint64_t none = std::numeric_limits<uint64_t>::max();

protected:
    // the entry of set_of_symbol for a symbol whose groups are not stored
    static constexpr uint32_t no_set = std::numeric_limits<uint32_t>::max();

    // ############################# VARIABLES #############################

    uint8_t log_group_size = 0; // log2 of the number of positions of a group

    // the set of each symbol: [0..sigma-1] set_of_symbol[c] is the index in sets of the groups of c, or no_set
    no_init_vector<uint32_t> set_of_symbol;

    /* the groups of the symbols whose groups are stored: [0..sets.size()-1] sets[set_of_symbol[c]] contains the groups
       h (in [0..ceil(n / 2^log_group_size) - 1]) whose positions h * 2^log_group_size..(h + 1) * 2^log_group_size - 1
       contain c */
    std::vector<rank_select_set> sets;

public:
    // ############################# CONSTRUCTORS #############################

    symbol_groups() = default;

    /**
     * @brief builds the groups of the symbols c with stored[c] that occur in input[0..n-1]
     * @tparam read_t type of read
     * @param read read(i) returns the symbol at position i
     * @param n the length of the input
     * @param group_size the number of positions of a group (a power of two)
     * @param stored [0..sigma-1] whether the groups of each symbol are stored
     * @param num_threads the number of threads
     */
    template <typename read_t>
    symbol_groups(const read_t& read, uint64_t n, uint64_t group_size, const std::vector<bool>& stored,
        uint64_t num_threads)
    {
        log_group_size = uint8_t(std::bit_width(group_size) - 1);
        const uint64_t sigma = stored.size();
        const uint64_t num_groups = div_ceil<uint64_t>(n, group_size);
        // the symbols whose groups are stored, numbered in increasing order: stored_rank[c] (no_set if !stored[c])
        std::vector<uint32_t> stored_rank(sigma, no_set);
        uint32_t num_stored = 0;
        for (uint64_t c = 0; c < sigma; c++) if (stored[c]) stored_rank[c] = num_stored++;
        const uint64_t threads = std::clamp<uint64_t>(num_threads, 1, std::max<uint64_t>(1, num_groups));
        // groups[t][s] = the groups of the s-th symbol in the groups of thread t, in increasing order
        std::vector<std::vector<std::vector<uint64_t>>> groups(threads, std::vector<std::vector<uint64_t>>(num_stored));

        // each thread collects the groups of a range of groups
        #pragma omp parallel for num_threads(threads) schedule(static, 1)
        for (uint64_t t = 0; t < threads; t++) {
            for (uint64_t h = t * num_groups / threads; h < (t + 1) * num_groups / threads; h++) {
                const uint64_t end = std::min(n, (h + 1) << log_group_size);

                for (uint64_t i = h << log_group_size; i < end; i++) {
                    const uint64_t c = uint64_t(read(i));
                    if (c >= sigma || stored_rank[c] == no_set) continue;
                    std::vector<uint64_t>& of_c = groups[t][stored_rank[c]];
                    if (of_c.empty() || of_c.back() != h) of_c.push_back(h);
                }
            }
        }

        // concatenate the groups of each symbol over the threads and build its set
        set_of_symbol.assign(sigma, no_set);

        for (uint64_t c = 0; c < sigma; c++) {
            if (stored_rank[c] == no_set) continue;
            std::vector<uint64_t> values;

            for (uint64_t t = 0; t < threads; t++) {
                std::vector<uint64_t>& of_c = groups[t][stored_rank[c]];
                values.insert(values.end(), of_c.begin(), of_c.end());
                std::vector<uint64_t>().swap(of_c);
            }

            if (values.empty()) continue;
            set_of_symbol[c] = uint32_t(sets.size());
            sets.emplace_back(values, num_groups, num_threads);
        }
    }

    // ############################# MISC PUBLIC METHODS #############################

    // returns whether no groups are stored
    inline bool empty() const { return sets.empty(); }

    // returns whether the groups of the symbol c are stored
    inline bool stores(uint64_t c) const { return c < set_of_symbol.size() && set_of_symbol[c] != no_set; }

    // returns the number of positions of a group
    inline uint64_t group_size() const { return uint64_t { 1 } << log_group_size; }

    // returns the size in bytes
    uint64_t size_in_bytes() const
    {
        uint64_t bytes = set_of_symbol.size() * sizeof(uint32_t);
        for (const rank_select_set& set : sets) bytes += set.size_in_bytes();
        return bytes;
    }

    // ############################# QUERY METHODS #############################

    // returns the group of position i
    inline uint64_t group_of(uint64_t i) const { return i >> log_group_size; }

    /**
     * @brief returns the first group from h on that contains the symbol c, or none
     * @param c a symbol whose groups are stored
     * @param h a group
     */
    inline uint64_t next(uint64_t c, uint64_t h) const
    {
        const rank_select_set& set = sets[set_of_symbol[c]];
        const uint64_t rank = set.rank_1(std::min(h, set.size()));
        return rank == set.num_ones() ? none : set.select_1(rank + 1);
    }

    /**
     * @brief returns the last group up to h that contains the symbol c, or none
     * @param c a symbol whose groups are stored
     * @param h a group
     */
    inline uint64_t previous(uint64_t c, uint64_t h) const
    {
        const rank_select_set& set = sets[set_of_symbol[c]];
        const uint64_t rank = set.rank_1(std::min(h + 1, set.size()));
        return rank == 0 ? none : set.select_1(rank);
    }

    // ############################# SERIALIZATION METHODS #############################

    /**
     * @brief writes the groups to an output stream
     * @param out the output stream
     */
    void serialize(std::ostream& out) const
    {
        const uint64_t sigma = set_of_symbol.size();
        const uint64_t num_sets = sets.size();
        out.write((char*) &log_group_size, sizeof(uint8_t));
        out.write((char*) &sigma, sizeof(uint64_t));
        out.write((char*) set_of_symbol.data(), sigma * sizeof(uint32_t));
        out.write((char*) &num_sets, sizeof(uint64_t));
        for (const rank_select_set& set : sets) set.serialize(out);
    }

    /**
     * @brief reads groups written by serialize
     * @param in the input stream
     */
    void load(std::istream& in)
    {
        uint64_t sigma;
        uint64_t num_sets;
        in.read((char*) &log_group_size, sizeof(uint8_t));
        in.read((char*) &sigma, sizeof(uint64_t));
        resize_no_init(set_of_symbol, sigma);
        in.read((char*) set_of_symbol.data(), sigma * sizeof(uint32_t));
        in.read((char*) &num_sets, sizeof(uint64_t));
        sets.resize(num_sets);
        for (rank_select_set& set : sets) set.load(in);
    }
};

#endif // MOVE_R_DATA_STRUCTURES_SYMBOL_GROUPS_HPP
