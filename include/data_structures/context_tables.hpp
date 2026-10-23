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

#ifndef MOVE_R_DATA_STRUCTURES_CONTEXT_TABLES_HPP
#define MOVE_R_DATA_STRUCTURES_CONTEXT_TABLES_HPP

#include <algorithm>
#include <cstdint>
#include <functional>
#include <istream>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

#include <data_structures/blocked_positions.hpp>
#include <data_structures/elias_fano_bit_vector.hpp>
#include <data_structures/interleaved_bit_aligned_vectors.hpp>
#include <data_structures/rank_select_set.hpp>
#include <stdexcept>
#include <misc/files.hpp>
#include <misc/log_structures.hpp>
#include <misc/utils.hpp>

// order of the entries of a level: by the intervals in T (forward) or by the intervals of the reversed patterns in T^R (reverse)
enum context_order_t : uint8_t { _forward_order = 0, _reverse_order = 1 };

// direction of an extension: primary = the extension direction every index supports, secondary = the other one
enum context_extension_t : uint8_t { _primary_extension = 0, _secondary_extension = 1 };

/**
 * @brief returns whether the keys of an extension direction in an order are computed as
 *        symbol * (number of parents) + parent (true) or as parent * sigma + symbol (false)
 * @param order order of the entries
 * @param key extension direction
 */
inline constexpr bool symbol_major_keys(context_order_t order, context_extension_t key)
{ return (order == _forward_order) == (key == _primary_extension); }

/**
 * @brief stores, for every pattern of length <= max_length() that occurs in T (or for every pattern of the one
 *        length single_length()), its suffix array interval in T and the interval of the reversed pattern in T^R,
 *        and the keys to go from a pattern to its extensions by one symbol
 *
 * A level stores the patterns of one length in two orders: sorted by the starting position of their suffix
 * array interval in T (forward order) and sorted by the starting position of the interval of their reverse
 * in T^R (reverse order). The entry of a pattern is its rank in an order. For each order, a level stores
 * 1. the starting positions of the intervals; the intervals of the patterns of one length are adjacent, so
 *    an interval ends at the position before the starting position of the next interval, and
 * 2. for each extension direction the keys of the entries; the key of an entry is computed from the added
 *    symbol and the entry of its parent (the pattern without the added symbol). The keys are increasing,
 *    so the entry of a pattern is the rank of its key.
 * The fields of the samples of a level are defined by the caller (see move_rb::build_context_tables).
 */
class context_tables {
public:
    // ############################# TYPES AND PARAMETERS #############################

    // invalid entry (returned if a pattern does not occur)
    static constexpr uint32_t no_entry = std::numeric_limits<uint32_t>::max();

    // identifies the format of serialized context tables
    static constexpr uint64_t format = 0x637478746162'0002;

    // the patterns of one length in one order
    struct order_t {
        /* the starting positions of the suffix array intervals of the entries: [0..num_entries] beg[e] is the
           starting position of the interval of entry e; beg[num_entries] = n; increasing */
        blocked_positions beg;

        /* (entry, end position of its interval) for each entry whose interval does not end at the position before
           the next interval (the suffixes in between are shorter than the pattern length); sorted by entry, at
           most (pattern length - 1) pairs */
        std::vector<std::pair<uint32_t, uint64_t>> end_exceptions;

        /* the keys of the entries for an extension in primary direction: [0..num_entries-1] the keys in increasing
           order, the key of entry e has rank e (see symbol_major_keys); in a single-length level, the
           patterns packed by packed_key */
        rank_select_set primary_keys;

        // the keys of the entries for an extension in secondary direction (as primary_keys)
        rank_select_set secondary_keys;

        // returns the size of the order in bytes
        uint64_t size_in_bytes() const
        {
            return beg.size_in_bytes() + primary_keys.size_in_bytes() + secondary_keys.size_in_bytes()
                + end_exceptions.size() * sizeof(std::pair<uint32_t, uint64_t>);
        }

        // writes the order to out
        void serialize(std::ostream& out) const
        {
            beg.serialize(out);
            primary_keys.serialize(out);
            secondary_keys.serialize(out);
            uint64_t size = end_exceptions.size();
            out.write((char*) &size, sizeof(uint64_t));

            // write each pair as two 64-bit words (the padding bytes of a pair are not initialized)
            static_assert(sizeof(std::pair<uint32_t, uint64_t>) == 2 * sizeof(uint64_t));
            std::vector<uint64_t> words;
            words.reserve(2 * size);

            for (const auto& [exc_entry, exc_end] : end_exceptions) {
                words.push_back(exc_entry);
                words.push_back(exc_end);
            }

            if (size != 0) out.write((const char*) words.data(), size * sizeof(std::pair<uint32_t, uint64_t>));
        }

        // reads an order written by serialize
        void load(std::istream& in)
        {
            beg.load(in);
            primary_keys.load(in);
            secondary_keys.load(in);
            uint64_t size;
            in.read((char*) &size, sizeof(uint64_t));
            end_exceptions.resize(size);
            if (size != 0) in.read((char*) end_exceptions.data(),
                size * sizeof(std::pair<uint32_t, uint64_t>));
        }
    };

    // the patterns of one length in both orders
    struct level_t {
        uint64_t num_entries = 0; // number of patterns (entries) of the level

        // orders[_forward_order] and orders[_reverse_order]; the latter is empty if the index has no reverse index
        order_t orders[2];

        // marks the entries (in forward order) with a stored suffix array sample, i.e., those whose interval
        // contains no run boundary of M_LF (see move_rb::sample_at_run_boundary)
        elias_fano_bit_vector sample_marks;

        // the stored samples: [0..number of marked entries - 1] the sample of the marked entry with rank i in sample_marks
        interleaved_bit_aligned_vectors samples;

        // returns the size of the level in bytes
        uint64_t size_in_bytes() const
        {
            return orders[0].size_in_bytes() + orders[1].size_in_bytes()
                + sample_marks.size_in_bytes() + samples.size_in_bytes() - sizeof(samples);
        }

        // writes the level to out
        void serialize(std::ostream& out) const
        {
            out.write((char*) &num_entries, sizeof(uint64_t));
            samples.serialize(out);
            sample_marks.serialize(out);
            for (const order_t& o : orders) o.serialize(out);
        }

        // reads a level written by serialize
        void load(std::istream& in)
        {
            in.read((char*) &num_entries, sizeof(uint64_t));
            samples.load(in);
            sample_marks.load(in);
            for (order_t& o : orders) o.load(in);
        }
    };

protected:
    // ############################# VARIABLES #############################

    uint64_t sigma = 0; // size of the alphabet the keys are computed with

    // the levels: [0..levels.size()-1] levels[l - 1] stores the patterns of length l (see single_len)
    std::vector<level_t> levels;

    // 0 <=> the levels of all lengths in [1..levels.size()] are stored; else levels.back() stores the patterns of
    // length single_len with their packed symbols as primary keys (see packed_key), and the levels before
    // it store the patterns of the lengths in [1..levels.size() - 1]
    uint64_t single_len = 0;

    // the entries of the patterns of length 1: [0..sigma-1] the entry of the pattern c (the same in both orders), or no_entry
    std::vector<uint32_t> symbol_entries;

    // ############################# INTERNAL METHODS #############################

    // returns the level storing the patterns of length len (see single_len)
    inline const level_t& level_of(uint64_t len) const
    { return single_len != 0 && len == single_len ? levels.back() : levels[len - 1]; }

    // returns the key of the extension of the pattern of length len with entry entry by sym (see symbol_major_keys)
    template <context_extension_t key, context_order_t order>
    inline uint64_t extension_key(uint64_t len, uint32_t entry, uint64_t sym) const
    {
        return symbol_major_keys(order, key)
            ? sym * num_entries(len) + entry
            : uint64_t(entry) * sigma + sym;
    }

public:
    // ############################# CONSTRUCTORS #############################

    context_tables() = default;

    // ############################# MISC PUBLIC METHODS #############################

    // returns whether no level is stored
    inline bool empty() const { return levels.empty(); }

    // returns the maximum length l, s.t. the patterns of all lengths in [1..l] are stored (0 if empty)
    inline uint64_t max_length() const { return single_len == 0 ? levels.size() : levels.size() - 1; }

    // returns single_len
    inline uint64_t single_length() const { return single_len; }

    // returns the alphabet size (the largest symbol + 1)
    inline uint64_t alphabet_size() const { return sigma; }

    // returns the number of patterns of length len (1 for the empty pattern)
    inline uint64_t num_entries(uint64_t len) const { return len == 0 ? 1 : levels[len - 1].num_entries; }

    // returns the number of stored patterns (of all levels)
    inline uint64_t num_patterns() const
    {
        uint64_t sum = 0;
        for (const level_t& l : levels) sum += l.num_entries;
        return sum;
    }

    // returns the average size of a stored pattern in bits
    inline double bits_per_pattern() const
    {
        uint64_t patterns = num_patterns();
        return patterns == 0 ? 0.0 : (8.0 * double(size_in_bytes())) / double(patterns);
    }

    // returns the size of the data structure in bytes
    uint64_t size_in_bytes() const
    {
        uint64_t size = symbol_entries.size() * sizeof(uint32_t);
        for (const level_t& l : levels) size += l.size_in_bytes();
        return size;
    }

    /**
     * @brief logs the alphabet size, the lengths and the numbers of patterns, one per line, the patterns of the
     *        lengths 1 to max_length() as a tree (each pattern below the pattern its last primary extension
     *        extended), and the patterns of the single-length level as a table; a pattern is logged with its suffix
     *        array interval in T, the interval of its reverse in T^R (if the reverse order is stored) and whether a
     *        sample is stored for it
     * @param symbol returns a symbol (in the internal effective alphabet) as it is logged
     * @param primary_prepends true <=> an extension in primary direction adds a symbol on the left
     * @param max_depth the patterns up to this length are logged in the tree
     * @param max_entries at most this many patterns are logged in the tree and in the table
     */
    void log_data_structures(const std::function<std::string(uint64_t)>& symbol, bool primary_prepends,
        uint64_t max_depth = default_logged_tree_depth, uint64_t max_entries = default_logged_entries) const;

    // ############################# ACCESS METHODS #############################

    // returns the stored samples of the patterns of length len, indexed by sample_rank()
    inline const interleaved_bit_aligned_vectors& samples(uint64_t len) const { return level_of(len).samples; }

    /* returns the keys of the extensions of the patterns of length len, i.e., the keys of the patterns of length
       len + 1, which are stored in levels[len] */
    template <context_extension_t key, context_order_t order>
    inline const rank_select_set& keys(uint64_t len) const
    {
        return key == _primary_extension ? levels[len].orders[order].primary_keys : levels[len].orders[order].secondary_keys;
    }

    /**
     * @brief returns the suffix array interval of a stored pattern
     * @tparam order _forward_order: the interval in T, _reverse_order: the interval of the reverse pattern in T^R
     * @param len length of the pattern
     * @param entry entry of the pattern in the given order
     * @return [b, e], the interval
     */
    template <context_order_t order>
    inline std::pair<uint64_t, uint64_t> interval(uint64_t len, uint32_t entry) const
    {
        const order_t& o = level_of(len).orders[order];
        // starting positions of the interval and of the next interval
        auto [b, next] = o.beg.pair_at(entry);

        for (const auto& [exc_entry, exc_end] : o.end_exceptions) {
            if (exc_entry == entry) return { b, exc_end };
        }

        return { b, next - 1 };
    }

    // returns the starting position of the suffix array interval of a stored pattern (see interval)
    template <context_order_t order>
    inline uint64_t interval_beg(uint64_t len, uint32_t entry) const
    { return level_of(len).orders[order].beg[entry]; }

    // ############################# QUERY METHODS #############################

    /**
     * @brief returns the entry of the pattern of length len whose suffix array interval starts at b
     * @tparam order order of the interval and of the returned entry
     * @param len length of the pattern (a level of this length must be stored)
     * @param b starting position of the interval
     * @return the entry, or no_entry if no interval of a pattern of length len starts at b
     */
    template <context_order_t order>
    inline uint32_t entry_of(uint64_t len, uint64_t b) const
    {
        const level_t& l = level_of(len);
        const uint64_t i = l.orders[order].beg.find(b);
        if (i >= l.num_entries) [[unlikely]] return no_entry;
        return uint32_t(i);
    }

    // returns whether the keys of the extensions of the patterns of length len are stored
    template <context_extension_t key, context_order_t order>
    inline bool has_keys(uint64_t len) const
    { return len < levels.size() && !keys<key, order>(len).empty(); }

    // returns whether the extensions of each pattern of length len are consecutive entries in the given order (see extension_range)
    template <context_extension_t key, context_order_t order>
    inline bool has_extension_range(uint64_t len) const
    { return !symbol_major_keys(order, key) && has_keys<key, order>(len); }

    /**
     * @brief computes the range of the entries of the extensions of a pattern (sorted by the added symbol);
     *        has_extension_range<key, order>(len) must hold
     * @param len length of the pattern to extend
     * @param entry entry of the pattern to extend (0 if len = 0)
     * @param start output: first entry of the range
     * @param end output: last entry of the range + 1
     */
    template <context_extension_t key, context_order_t order>
    inline void extension_range(uint64_t len, uint32_t entry, uint32_t& start, uint32_t& end) const
    {
        static_assert(!symbol_major_keys(order, key));
        // the keys of the extensions of entry lie in [entry * sigma, (entry + 1) * sigma)
        const rank_select_set& k = keys<key, order>(len);
        start = uint32_t(k.rank_1(uint64_t(entry) * sigma));
        end = uint32_t(k.rank_1((uint64_t(entry) + 1) * sigma));
    }

    /**
     * @brief returns the entry of a pattern extended by a symbol
     * @tparam key extension direction
     * @tparam order order of the given and the returned entry
     * @param len length of the pattern to extend
     * @param entry entry of the pattern to extend (0 if len = 0)
     * @param sym symbol to extend with (in the internal effective alphabet)
     * @return entry of the extended pattern, or no_entry if it does not occur
     */
    template <context_extension_t key, context_order_t order>
    inline uint32_t extend(uint64_t len, uint32_t entry, uint64_t sym) const
    {
        // the extensions of the empty pattern are stored in symbol_entries
        if (len == 0) return symbol_entries[sym];

        const rank_select_set& k = keys<key, order>(len);
        const uint64_t rank = k.index_of(extension_key<key, order>(len, entry, sym));
        if (rank == k.num_ones()) [[unlikely]] return no_entry;
        return uint32_t(rank);
    }

    // returns the entry of a pattern extended by a symbol (see extend); the extended pattern must occur
    template <context_extension_t key, context_order_t order>
    inline uint32_t extend_occurring(uint64_t len, uint32_t entry, uint64_t sym) const
    {
        if (len == 0) return symbol_entries[sym];
        return uint32_t(keys<key, order>(len).rank_1(extension_key<key, order>(len, entry, sym)));
    }

    /**
     * @brief returns the index of the stored sample of a pattern (see level_t::samples)
     * @param len length of the pattern
     * @param entry entry of the pattern in forward order
     * @return index in samples(len), or no_entry if no sample is stored for the pattern
     */
    inline uint32_t sample_rank(uint64_t len, uint32_t entry) const
    {
        const elias_fano_bit_vector& s = level_of(len).sample_marks;
        uint64_t rank = s.rank_1(entry);

        // check whether entry is marked
        if (rank >= s.num_ones() || s.select_1(rank + 1) != entry) return no_entry;
        return uint32_t(rank);
    }

    // returns whether the keys of the patterns of length len are stored (see parent)
    template <context_extension_t key, context_order_t order>
    inline bool has_parent(uint64_t len) const
    { return len >= 1 && len <= max_length() && !keys<key, order>(len - 1).empty(); }

    /**
     * @brief returns the entry of a pattern without the symbol added by an extension;
     *        has_parent<key, order>(len) must hold
     * @tparam key direction of the extension that added the symbol to remove
     * @tparam order order of the given and the returned entry
     * @param len length of the pattern
     * @param entry entry of the pattern
     * @return entry of the pattern of length len - 1 (0 if len = 1)
     */
    template <context_extension_t key, context_order_t order>
    inline uint32_t parent(uint64_t len, uint32_t entry) const
    {
        const rank_select_set& k = keys<key, order>(len - 1);

        if constexpr (symbol_major_keys(order, key)) {
            return uint32_t(k.select_1(entry + 1) % num_entries(len - 1));
        } else {
            return uint32_t(k.select_1(entry + 1) / sigma);
        }
    }

    /**
     * @brief returns the key of a pattern in a single-length table, i.e., its symbols packed to the base
     *        sigma with sym(0) as the most significant symbol
     * @tparam sym_fnc_t type of sym
     * @param sym function returning the symbols (in the internal effective alphabet); for the forward order,
     *        sym(i) is the i-th symbol of the pattern, for the reverse order the i-th symbol of the reverse pattern
     * @param len length of the pattern (must be single_length())
     */
    template <typename sym_fnc_t>
    inline uint64_t packed_key(sym_fnc_t sym, uint64_t len) const
    {
        uint64_t key = 0;
        for (uint64_t i = 0; i < len; i++) key = key * sigma + sym(i);
        return key;
    }

    /**
     * @brief returns the entry of a pattern in a single-length table
     * @tparam order order of the returned entry (the key must be packed for this order, see packed_key)
     * @param key packed key of the pattern
     * @return the entry, or no_entry if the pattern does not occur
     */
    template <context_order_t order>
    inline uint32_t lookup(uint64_t key) const
    {
        const rank_select_set& k = levels.back().orders[order].primary_keys;
        const uint64_t rank = k.index_of(key);
        if (rank == k.num_ones()) [[unlikely]] return no_entry;
        return uint32_t(rank);
    }

    // ############################# CONSTRUCTION METHODS #############################

    // sets the alphabet size (must be set before the first level is appended)
    void set_alphabet_size(uint64_t size) { sigma = size; }

    // sets single_len
    void set_single_length(uint64_t len) { single_len = len; }

    // appends a level (the levels must be appended in increasing order of the pattern length)
    void push_level(level_t&& level)
    {
        levels.emplace_back(std::move(level));
    }

    /**
     * @brief sets the entries of the patterns of length 1 (see symbol_entries)
     * @param entries [0..sigma-1] entry of each symbol (no_entry if it does not occur)
     */
    void set_symbol_entries(std::vector<uint32_t>&& entries)
    {
        symbol_entries = std::move(entries);
        advise_huge_pages(symbol_entries.data(), symbol_entries.size() * sizeof(uint32_t));
    }

    // removes all levels
    void clear()
    {
        sigma = 0;
        single_len = 0;
        levels.clear();
        levels.shrink_to_fit();
        symbol_entries.clear();
        symbol_entries.shrink_to_fit();
    }

    // ############################# SERIALIZATION METHODS #############################

    /**
     * @brief writes the context tables to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) const
    {
        out.write((const char*) &format, sizeof(uint64_t));
        uint64_t num_levels = levels.size();
        out.write((char*) &num_levels, sizeof(uint64_t));
        out.write((char*) &sigma, sizeof(uint64_t));
        out.write((char*) &single_len, sizeof(uint64_t));
        write_vector(out, symbol_entries);

        for (const level_t& l : levels) l.serialize(out);
    }

    /**
     * @brief reads context tables written by serialize
     * @param in input stream
     */
    void load(std::istream& in)
    {
        uint64_t stored_format;
        in.read((char*) &stored_format, sizeof(uint64_t));

        if (stored_format != format)
            throw std::runtime_error("the index stores context tables of an older format; build it again");

        uint64_t num_levels;
        in.read((char*) &num_levels, sizeof(uint64_t));
        in.read((char*) &sigma, sizeof(uint64_t));
        in.read((char*) &single_len, sizeof(uint64_t));
        read_vector(in, symbol_entries);
        levels.resize(num_levels);

        for (level_t& l : levels) l.load(in);
    }
};

#include <data_structures/context_tables.tpp>

#endif // MOVE_R_DATA_STRUCTURES_CONTEXT_TABLES_HPP
