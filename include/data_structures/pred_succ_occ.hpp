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

#ifndef MOVE_R_DATA_STRUCTURES_PRED_SUCC_OCC_HPP
#define MOVE_R_DATA_STRUCTURES_PRED_SUCC_OCC_HPP

#include <bit>
#include <cstring>
#include <functional>
#include <limits>
#include <optional>
#include <vector>

#include <data_structures/interleaved_bit_aligned_vectors.hpp>
#include <data_structures/symbol_groups.hpp>
#include <misc/utils.hpp>

/**
 * @brief character predecessor/successor queries over an input read through a stored function read(i),
 *        using the stored occurrences of each symbol
 *
 * For each symbol, the input positions are divided into blocks whose size is a power of two that is chosen
 * s.t. a block contains about occs_per_block occurrences of the symbol. The occurrences of a symbol are
 * stored in increasing order as offsets within their blocks, and for each block the index of the first
 * occurrence in or after it is stored. A query reads the entry of a block and scans the occurrences of the
 * symbol in that block. The occurrences of a symbol with an average gap <= max_gap_scanned are not stored; instead,
 * the groups of max_gap_scanned positions that contain it are stored (see symbol_groups). A query for such a symbol
 * scans the input up to the end of the group after the one of the query position, and then the next group that
 * contains the symbol.
 *
 * @tparam read_fnc_t type of the read function (read(i) returns the symbol at index i)
 */
template <typename read_fnc_t>
class pred_succ_occ {
protected:
    // ############################# TYPES AND PARAMETERS #############################

    // default upper bound of succ (no bound)
    static constexpr uint64_t unbounded = std::numeric_limits<uint64_t>::max();

    // maximum average gap between the occurrences of a symbol whose occurrences are not stored
    static constexpr uint64_t max_gap_scanned = 16;
    static_assert(std::has_single_bit(max_gap_scanned));

    // the fields stored for each symbol in symbol_fields
    enum symbol_field_t : uint8_t {
        _num_occs = 0, // number of occurrences (0 <=> the occurrences of the symbol are not stored)
        _block_offset = 1, // bit position in first_occ_of_block of the entry of the first block of the symbol
        _occ_offset = 2, // bit position in occ_offsets of the first occurrence of the symbol
        _log_block_size = 3 // log2 of the block size of the symbol
    };

    // the stored fields of a symbol and the values computed from them
    struct symbol_t {
        uint64_t block_offset;
        uint64_t occ_offset;
        uint64_t num_occs;
        uint64_t num_blocks; // div_ceil(input_size, block size) + 1
        uint8_t width_block; // bits per block entry (std::bit_width(num_occs))
        uint8_t width_occ; // bits per occurrence (log_block_size)
        uint8_t log_block_size;
    };

    // ############################# VARIABLES #############################

    mutable std::optional<read_fnc_t> read; // the stored read function (rebindable with set_read)

    // average number of occurrences of a symbol per block
    uint64_t occs_per_block = 1;

    uint64_t input_size = 0; // number of input positions
    uint64_t sigma = 0; // alphabet size

    // the fields of the symbols: [0..sigma-1] the fields of symbol c (see symbol_field_t)
    interleaved_bit_aligned_vectors symbol_fields;

    /* for each symbol and block, the index of the first occurrence of the symbol in or after the block, bit-packed:
       the entries of symbol c start at bit _block_offset(c); entry b in [0..num_blocks(c)-1] (bit_width(num_occs(c))
       bits) is the index (among the occurrences of c) of the first occurrence of c at a position >= b * 2^log_block_size(c),
       or num_occs(c) if there is none */
    no_init_vector<uint64_t> first_occ_of_block;

    /* the occurrences of each symbol as offsets within their blocks, bit-packed: the occurrences of symbol c start at
       bit _occ_offset(c); entry i in [0..num_occs(c)-1] (log_block_size(c) bits) is (position of the i-th occurrence of
       c) mod 2^log_block_size(c) */
    no_init_vector<uint64_t> occ_offsets;

    // the groups of max_gap_scanned positions of the symbols whose occurrences are not stored
    symbol_groups frequent_symbol_groups;

    // ############################# INTERNAL METHODS #############################

    // returns the fields of a symbol
    inline symbol_t fields_of(uint64_t sym) const
    {
        symbol_t s;
        s.num_occs = symbol_fields.template get<_num_occs, uint64_t>(sym);
        s.block_offset = symbol_fields.template get<_block_offset, uint64_t>(sym);
        s.occ_offset = symbol_fields.template get<_occ_offset, uint64_t>(sym);
        s.log_block_size = uint8_t(symbol_fields.template get<_log_block_size, uint64_t>(sym));
        const uint64_t block_size = uint64_t { 1 } << s.log_block_size;
        s.num_blocks = ((input_size + block_size - 1) >> s.log_block_size) + 1;
        s.width_block = std::bit_width(s.num_occs);
        s.width_occ = s.log_block_size;
        return s;
    }

    /**
     * @brief returns the index (among the occurrences of a symbol) of its first occurrence in or after a block
     * @param s fields of the symbol
     * @param block [0..s.num_blocks-1] block index
     */
    inline uint64_t block_first_occ(const symbol_t& s, uint64_t block) const
    {
        return read_bits(first_occ_of_block, s.block_offset + block * s.width_block, s.width_block);
    }

    /**
     * @brief returns the input position of an occurrence of a symbol
     * @param s fields of the symbol
     * @param block index of the block containing the occurrence
     * @param i [0..s.num_occs-1] index of the occurrence among the occurrences of the symbol
     */
    inline uint64_t occ_at(const symbol_t& s, uint64_t block, uint64_t i) const
    {
        return (block << s.log_block_size)
            | read_bits(occ_offsets, s.occ_offset + i * s.width_occ, s.width_occ);
    }

    /**
     * @brief returns the smallest index in [x, max] whose symbol is sym (a value > max if none), with the groups of sym
     * @param sym a symbol whose occurrences are not stored
     * @param x lower bound (the first position of a group)
     * @param max upper bound
     */
    [[gnu::noinline]] uint64_t succ_in_groups(uint64_t sym, uint64_t x, uint64_t max) const
    {
        const read_fnc_t& rd = *read;
        const uint64_t h = frequent_symbol_groups.stores(sym)
            ? frequent_symbol_groups.next(sym, frequent_symbol_groups.group_of(x)) : symbol_groups::none;
        if (h == symbol_groups::none || h * frequent_symbol_groups.group_size() > max) return max + 1;
        for (x = h * frequent_symbol_groups.group_size(); x <= max && sym != symbol_bit_pattern(rd(x)); x++) {}
        return x;
    }

    /**
     * @brief returns the largest index in [min, x) whose symbol is sym, with the groups of sym
     * @param sym a symbol whose occurrences are not stored and that occurs in [min, x)
     * @param x upper bound (the first position of a group h > 0)
     * @param min lower bound
     */
    [[gnu::noinline]] uint64_t pred_in_groups(uint64_t sym, uint64_t x, uint64_t min) const
    {
        const read_fnc_t& rd = *read;
        const uint64_t h = frequent_symbol_groups.stores(sym)
            ? frequent_symbol_groups.previous(sym, frequent_symbol_groups.group_of(x) - 1) : symbol_groups::none;
        if (h == symbol_groups::none) return min;
        for (x = (h + 1) * frequent_symbol_groups.group_size() - 1; x > min && sym != symbol_bit_pattern(rd(x)); x--) {}
        return x;
    }

public:
    // ############################# CONSTRUCTORS #############################

    pred_succ_occ() = default;

    /**
     * @brief builds the data structure over input[0..input_size-1] read through the function read
     * @param read function returning the symbol at index i
     * @param input_size number of input positions
     * @param alphabet_size alphabet size
     * @param occs_per_block average number of occurrences of a symbol per block
     * @param num_threads number of threads to use for counting the occurrences
     * @param scan_frequent true <=> the occurrences of a symbol with an average gap <= max_gap_scanned are not stored
     *        (see pred_succ_occ); else a query takes O(occs_per_block) time for every symbol
     */
    pred_succ_occ(read_fnc_t read, uint64_t input_size, uint64_t alphabet_size,
        uint64_t occs_per_block = 1, uint64_t num_threads = 1, bool scan_frequent = true)
    {
        this->input_size = input_size;
        this->sigma = alphabet_size;
        this->occs_per_block = std::max<uint64_t>(1, occs_per_block);
        this->read.emplace(read);
        const read_fnc_t& rd = *this->read;

        // count[c] = number of occurrences of c
        interleaved_bit_aligned_vectors count({ std::max<uint8_t>(1, std::bit_width(input_size)) });
        count.resize(sigma);

        // count the occurrences (in parallel: with one counter array per thread)
        {
            const uint64_t chunks = num_threads <= 1 ? 1 : num_threads;
            std::vector<std::vector<uint64_t>> count_chunk(
                chunks == 1 ? 0 : chunks, std::vector<uint64_t>(sigma, 0));

            if (chunks == 1) {
                for (uint64_t i = 0; i < input_size; i++) {
                    uint64_t c = symbol_bit_pattern(rd(i));
                    count.template set<0, uint64_t>(c, count.template get<0, uint64_t>(c) + 1);
                }
            } else {
                #pragma omp parallel for num_threads(num_threads) schedule(static)
                for (uint64_t t = 0; t < chunks; t++) {
                    std::vector<uint64_t>& chunk_count = count_chunk[t];

                    for (uint64_t i = (t * input_size) / chunks;
                         i < ((t + 1) * input_size) / chunks; i++) {
                        chunk_count[symbol_bit_pattern(rd(i))]++;
                    }
                }

                #pragma omp parallel for num_threads(num_threads)
                for (uint64_t c = 0; c < sigma; c++) {
                    uint64_t total = 0;
                    for (uint64_t t = 0; t < chunks; t++) total += count_chunk[t][c];
                    count.template set_parallel<0, uint64_t>(c, total);
                }
            }
        }

        // returns log2 of the block size of symbol c, or 0 if its occurrences are not stored
        auto log_block_size_of = [&](uint64_t c) -> uint64_t {
            uint64_t occurrences = count.template get<0, uint64_t>(c);
            if (occurrences == 0) return 0;

            uint64_t gap = input_size / occurrences;
            if (scan_frequent && gap <= max_gap_scanned) return 0;

            // about occs_per_block occurrences per block, rounded down to a power of two
            return std::max<uint64_t>(1, std::bit_width(gap * this->occs_per_block) - 1);
        };

        auto is_stored = [&](uint64_t c) { return log_block_size_of(c) != 0; };

        // compute the sizes of first_occ_of_block and occ_offsets in bits and the maximum block size
        uint64_t block_bits = 0;
        uint64_t occ_bits = 0;
        uint64_t max_log_block_size = 0;

        for (uint64_t c = 0; c < sigma; c++) {
            if (!is_stored(c)) continue;
            const uint64_t log_block_size = log_block_size_of(c);
            const uint64_t num_occs = count.template get<0, uint64_t>(c);
            block_bits += (div_ceil<uint64_t>(input_size, uint64_t { 1 } << log_block_size) + 1)
                * std::bit_width(num_occs);
            occ_bits += num_occs * log_block_size;
            max_log_block_size = std::max(max_log_block_size, log_block_size);
        }

        // the fields of the symbols, in the order of symbol_field_t
        symbol_fields = interleaved_bit_aligned_vectors({
            std::max<uint8_t>(1, std::bit_width(input_size)),        // _num_occs
            std::max<uint8_t>(1, std::bit_width(block_bits)),        // _block_offset
            std::max<uint8_t>(1, std::bit_width(occ_bits)),          // _occ_offset
            std::max<uint8_t>(1, std::bit_width(max_log_block_size)) }); // _log_block_size
        symbol_fields.resize(sigma);

        {
            uint64_t block_offset = 0;
            uint64_t occ_offset = 0;

            for (uint64_t c = 0; c < sigma; c++) {
                const uint64_t log_block_size = log_block_size_of(c);
                const uint64_t num_occs = log_block_size == 0
                    ? 0 : count.template get<0, uint64_t>(c);

                symbol_fields.template set<_block_offset, uint64_t>(c, block_offset);
                symbol_fields.template set<_occ_offset, uint64_t>(c, occ_offset);
                symbol_fields.template set<_num_occs, uint64_t>(c, num_occs);
                symbol_fields.template set<_log_block_size, uint64_t>(c, log_block_size);
                if (num_occs == 0) continue;

                block_offset += (div_ceil<uint64_t>(input_size, uint64_t { 1 } << log_block_size) + 1)
                    * std::bit_width(num_occs);
                occ_offset += num_occs * log_block_size;
            }
        }

        first_occ_of_block.assign(block_bits / 64 + 2, 0);
        occ_offsets.assign(occ_bits / 64 + 2, 0);
        advise_huge_pages(first_occ_of_block.data(), first_occ_of_block.size() * sizeof(uint64_t));
        advise_huge_pages(occ_offsets.data(), occ_offsets.size() * sizeof(uint64_t));

        /* fill first_occ_of_block and occ_offsets in one pass over the input; cur_occ[c] = number of occurrences of c
           before the current position, next_block[c] = index of the next block of c whose entry is written */
        uint8_t counter_width = std::max<uint8_t>(1, std::bit_width(input_size));
        interleaved_bit_aligned_vectors cur_occ({ counter_width });
        interleaved_bit_aligned_vectors next_block({ counter_width });
        cur_occ.resize(sigma);
        next_block.resize(sigma);

        for (uint64_t i = 0; i < input_size; i++) {
            uint64_t c = symbol_bit_pattern(rd(i));
            const symbol_t s = fields_of(c);
            if (s.num_occs == 0) continue;

            while (next_block.template get<0, uint64_t>(c) < s.num_blocks
                && (next_block.template get<0, uint64_t>(c) << s.log_block_size) <= i) {
                write_bits(first_occ_of_block, s.block_offset + next_block.template get<0, uint64_t>(c) * s.width_block,
                    s.width_block, cur_occ.template get<0, uint64_t>(c));
                next_block.template set<0, uint64_t>(c, next_block.template get<0, uint64_t>(c) + 1);
            }

            uint64_t mask = (uint64_t { 1 } << s.log_block_size) - 1;
            write_bits(occ_offsets, s.occ_offset + cur_occ.template get<0, uint64_t>(c) * s.width_occ, s.width_occ, i & mask);
            cur_occ.template set<0, uint64_t>(c, cur_occ.template get<0, uint64_t>(c) + 1);
        }

        // the entries of the blocks after the last occurrence of a symbol are its number of occurrences
        for (uint64_t c = 0; c < sigma; c++) {
            const symbol_t s = fields_of(c);
            if (s.num_occs == 0) continue;

            while (next_block.template get<0, uint64_t>(c) < s.num_blocks) {
                write_bits(first_occ_of_block, s.block_offset + next_block.template get<0, uint64_t>(c) * s.width_block,
                    s.width_block, cur_occ.template get<0, uint64_t>(c));
                next_block.template set<0, uint64_t>(c, next_block.template get<0, uint64_t>(c) + 1);
            }
        }

        // the groups of the symbols that occur and whose occurrences are not stored
        std::vector<bool> scanned(sigma);
        for (uint64_t c = 0; c < sigma; c++) scanned[c] = count.template get<0, uint64_t>(c) != 0 && !is_stored(c);
        frequent_symbol_groups = symbol_groups([&rd](uint64_t i) { return symbol_bit_pattern(rd(i)); }, input_size,
            max_gap_scanned, scanned, num_threads);
    }

    // ############################# MISC PUBLIC METHODS #############################

    // (re-)binds the read function; call after loading, or after the object read references has moved
    inline void set_read(read_fnc_t read) const { this->read.emplace(read); }

    // returns whether the data structure is empty (i.e., has not been built)
    inline bool empty() const { return symbol_fields.empty(); }

    // returns whether the occurrences of a symbol are stored (else its groups are, see pred_succ_occ)
    inline bool stores_occurrences(uint64_t sym) const { return symbol_fields.template get<_num_occs, uint64_t>(sym) != 0; }

    // returns the groups of the symbols whose occurrences are not stored
    inline const symbol_groups& frequent_groups() const { return frequent_symbol_groups; }

    // returns the size of the data structure in bytes
    uint64_t size_in_bytes() const
    {
        return symbol_fields.size_in_bytes() + first_occ_of_block.size() * sizeof(uint64_t)
            + occ_offsets.size() * sizeof(uint64_t) + frequent_symbol_groups.size_in_bytes();
    }

    // ############################# QUERY METHODS #############################

    /**
     * @brief returns the smallest index in [x, max] whose symbol is sym (a value > max if none)
     * @param sym a symbol
     * @param x lower bound
     * @param max upper bound
     */
    inline uint64_t succ(uint64_t sym, uint64_t x, uint64_t max = unbounded) const
    {
        const read_fnc_t& rd = *read;
        const symbol_t s = fields_of(sym);
        // max + 1 is returned if there is no such index, so max is limited to input_size - 1 (unbounded + 1 = 0)
        max = std::min<uint64_t>(max, input_size - 1);

        if (s.num_occs == 0) [[unlikely]] {
            // the occurrences of sym are not stored: scan the input up to the end of the group after the one of x
            const uint64_t scan_end = std::min(max, (frequent_symbol_groups.group_of(x) + 2) * frequent_symbol_groups.group_size() - 1);
            while (x <= scan_end && sym != symbol_bit_pattern(rd(x))) x++;
            return x <= scan_end || scan_end == max ? x : succ_in_groups(sym, x, max);
        }

        uint64_t block = x >> s.log_block_size;
        uint64_t offset = x & ((uint64_t { 1 } << s.log_block_size) - 1);
        uint64_t i = block_first_occ(s, block);
        uint64_t end = block + 1 < s.num_blocks ? block_first_occ(s, block + 1) : uint64_t(s.num_occs);

        // search the occurrences in the block of x
        for (; i < end; i++) {
            uint64_t occ = read_bits(occ_offsets, s.occ_offset + i * s.width_occ, s.width_occ);
            if (occ >= offset) {
                uint64_t pos = (block << s.log_block_size) | occ;
                return pos <= max ? pos : max + 1;
            }
        }

        if (i == s.num_occs) [[unlikely]] return max + 1;

        // the result is the first occurrence after the block of x: find the block containing occurrence i
        while (++block < s.num_blocks && block_first_occ(s, block) == i) {}
        uint64_t pos = occ_at(s, block - 1, i);
        return pos <= max ? pos : max + 1;
    }

    /**
     * @brief returns the largest index in [min, x] whose symbol is sym (sym must occur in [min, x])
     * @param sym a symbol
     * @param x upper bound
     * @param min lower bound
     */
    inline uint64_t pred(uint64_t sym, uint64_t x, uint64_t min = 0) const
    {
        const read_fnc_t& rd = *read;
        const symbol_t s = fields_of(sym);

        if (s.num_occs == 0) [[unlikely]] {
            // the occurrences of sym are not stored: scan the input down to the start of the group before the one of x
            const uint64_t h_x = frequent_symbol_groups.group_of(x);
            const uint64_t scan_begin = std::max(min, h_x == 0 ? 0 : (h_x - 1) * frequent_symbol_groups.group_size());
            while (x > scan_begin && sym != symbol_bit_pattern(rd(x))) x--;
            return x > scan_begin || x == min || sym == symbol_bit_pattern(rd(x)) ? x : pred_in_groups(sym, x, min);
        }

        uint64_t block = x >> s.log_block_size;
        uint64_t offset = x & ((uint64_t { 1 } << s.log_block_size) - 1);
        uint64_t begin = block_first_occ(s, block);
        uint64_t end = block + 1 < s.num_blocks ? block_first_occ(s, block + 1) : uint64_t(s.num_occs);

        // search the occurrences in the block of x from right to left
        for (uint64_t i = end; i > begin; i--) {
            uint64_t occ = read_bits(occ_offsets, s.occ_offset + (i - 1) * s.width_occ, s.width_occ);
            if (occ <= offset) return (block << s.log_block_size) | occ;
        }

        // the result is the last occurrence before the block of x: find the block containing occurrence begin - 1
        while (block > 0 && block_first_occ(s, block) == begin) block--;
        return occ_at(s, block, begin - 1);
    }

    // ############################# SERIALIZATION METHODS #############################

    /**
     * @brief writes the data structure to an output stream (the read function is not serialized)
     * @param out output stream
     */
    void serialize(std::ostream& out) const
    {
        out.write((char*) &input_size, sizeof(uint64_t));
        out.write((char*) &sigma, sizeof(uint64_t));
        out.write((char*) &occs_per_block, sizeof(uint64_t));
        symbol_fields.serialize(out);
        uint64_t num_blocks = first_occ_of_block.size();
        uint64_t num_occs = occ_offsets.size();
        out.write((char*) &num_blocks, sizeof(uint64_t));
        out.write((char*) &num_occs, sizeof(uint64_t));
        out.write((char*) first_occ_of_block.data(), num_blocks * sizeof(uint64_t));
        out.write((char*) occ_offsets.data(), num_occs * sizeof(uint64_t));
        frequent_symbol_groups.serialize(out);
    }

    // reads the data structure from an input stream; set_read has to be called afterwards
    void load(std::istream& in)
    {
        in.read((char*) &input_size, sizeof(uint64_t));
        in.read((char*) &sigma, sizeof(uint64_t));
        in.read((char*) &occs_per_block, sizeof(uint64_t));
        symbol_fields.load(in);
        uint64_t num_blocks;
        uint64_t num_occs;
        in.read((char*) &num_blocks, sizeof(uint64_t));
        in.read((char*) &num_occs, sizeof(uint64_t));
        resize_no_init(first_occ_of_block, num_blocks);
        resize_no_init(occ_offsets, num_occs);
        advise_huge_pages(first_occ_of_block.data(), first_occ_of_block.size() * sizeof(uint64_t));
        advise_huge_pages(occ_offsets.data(), occ_offsets.size() * sizeof(uint64_t));
        in.read((char*) first_occ_of_block.data(), num_blocks * sizeof(uint64_t));
        in.read((char*) occ_offsets.data(), num_occs * sizeof(uint64_t));
        frequent_symbol_groups.load(in);
    }
};

#endif // MOVE_R_DATA_STRUCTURES_PRED_SUCC_OCC_HPP
