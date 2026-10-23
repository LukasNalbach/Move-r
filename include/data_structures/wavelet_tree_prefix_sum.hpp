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

#ifndef MOVE_R_DATA_STRUCTURES_WAVELET_TREE_PREFIX_SUM_HPP
#define MOVE_R_DATA_STRUCTURES_WAVELET_TREE_PREFIX_SUM_HPP

#include <bit>
#include <cassert>
#include <cstdint>
#include <istream>
#include <ostream>
#include <utility>
#include <vector>

#include <misc/utils.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v.hpp>
#include <sdsl/util.hpp>

#include <data_structures/elias_fano_bit_vector.hpp>
#include <data_structures/interleaved_bit_aligned_vectors.hpp>

/**
 * @brief a wavelet tree over a sequence of symbols with weights, supporting the query smaller(i, j, c): the
 *        sum of the weights of the entries in [i, j) whose symbol is smaller than c
 *
 * A node represents a range of the alphabet and contains the entries whose symbol lies in it. It splits its
 * range at the symbol that divides its entries most evenly; an internal node is identified by this symbol,
 * and the leaves are in symbol order. A level stores the blocks (the entries) of its internal nodes
 * concatenated in symbol order: a bit vector marking the entries that belong to a right child, and the prefix
 * sums of the weights of the entries after each block has been partitioned stably into the entries of the
 * left and of the right child. The entries of a leaf are not stored. A query descends from the root to a
 * leaf and adds, wherever it descends to a right child, the weights of the entries of the left child.
 */
class wavelet_tree_prefix_sum {
public:
    // ############################# TYPES AND PARAMETERS #############################

    // a child that is a leaf (no alphabet range is split at symbol 0)
    static constexpr uint32_t leaf = 0;

protected:
    // the fields of _nodes
    static constexpr uint8_t field_left = 0;
    static constexpr uint8_t field_right = 1;
    static constexpr uint8_t field_start = 2;
    static constexpr uint8_t field_ones = 3;

    // ############################# VARIABLES #############################

    uint64_t num_entries = 0; // length of the sequence
    uint64_t total = 0; // sum of all weights
    uint64_t sigma = 0; // the largest symbol + 1
    uint32_t root = leaf; // the symbol at which the root splits its alphabet range
    uint32_t height = 0; // number of levels

    /* the bit vectors of the levels: [0..height-1] _bits[l][i] = 1 <=> the entry at position i of level l (the blocks of
       the internal nodes of level l concatenated in symbol order) belongs to the right child of its node */
    std::vector<sdsl::bit_vector> _bits;

    std::vector<sdsl::rank_support_v<>> _rank; // [0..height-1] _rank[l] is the rank support of _bits[l]

    /* the prefix sums of the weights of the levels: [0..height-1] _weight_prefix_sums[l] stores PS_l[0..m_l] as an
       Elias-Fano bit vector (select_1(k + 1) = PS_l[k]), where m_l is the number of entries of level l and PS_l[k] is
       the sum of the weights of the first k entries of level l after each block has been partitioned */
    std::vector<elias_fano_bit_vector> _weight_prefix_sums;

    // [1..sigma-1] the internal nodes, indexed by the symbol at which they split their range; fields: the
    // left and the right child (leaf <=> the child is a leaf), the starting position of the block of the
    // node in its level and the number of ones in its level before the block
    interleaved_bit_aligned_vectors _nodes;

    // ############################# INTERNAL METHODS #############################

    // sets the bit vector of the rank support of each level (required after copying or moving)
    void rebind()
    {
        for (uint32_t l = 0; l < height; l++) _rank[l].set_vector(&_bits[l]);
        advise_levels();
    }

    // advises huge pages for the bit vectors of all levels
    void advise_levels()
    {
        for (uint32_t l = 0; l < height; l++)
            advise_huge_pages((void*) _bits[l].data(), (_bits[l].size() + 7) / 8);
    }

    // returns the sum of the weights of the first k entries of a level (after its blocks have been partitioned)
    inline uint64_t weight_prefix(uint32_t level, uint64_t k) const { return _weight_prefix_sums[level].select_1(k + 1); }

    /**
     * @brief builds the subtree of the alphabet range [lo, hi)
     * @param cum [0..sigma] cum[c] = number of entries with a symbol < c
     * @param lo first symbol of the range
     * @param hi last symbol of the range + 1
     * @return the symbol at which the range is split, or leaf if it contains only one symbol
     */
    uint32_t build_subtree(const std::vector<uint64_t>& cum, uint64_t lo, uint64_t hi)
    {
        if (hi - lo <= 1) return leaf;

        uint64_t split;

        if (cum[hi] == cum[lo]) {
            // no entry has a symbol in the range: split it in the middle
            split = lo + (hi - lo) / 2;
        } else {
            // binary search for the first symbol with at least half of the entries of the range before it
            uint64_t half = cum[lo] + (cum[hi] - cum[lo]) / 2;
            uint64_t search_lo = lo + 1;
            uint64_t search_hi = hi - 1;

            while (search_lo < search_hi) {
                uint64_t mid = search_lo + (search_hi - search_lo) / 2;
                if (cum[mid] < half) search_lo = mid + 1; else search_hi = mid;
            }

            split = search_lo;
            // choose the previous symbol, if it divides the entries more evenly
            if (split > lo + 1 && half - cum[split - 1] < cum[split] - half) split--;

            // one side would contain no entry: split the range in the middle
            if (cum[split] == cum[lo] || cum[split] == cum[hi]) split = lo + (hi - lo) / 2;
        }

        _nodes.set<field_left>(split, build_subtree(cum, lo, split));
        _nodes.set<field_right>(split, build_subtree(cum, split, hi));
        return uint32_t(split);
    }

public:
    // ############################# CONSTRUCTORS #############################

    wavelet_tree_prefix_sum() = default;

    wavelet_tree_prefix_sum(const wavelet_tree_prefix_sum& other) { *this = other; }
    wavelet_tree_prefix_sum(wavelet_tree_prefix_sum&& other) noexcept { *this = std::move(other); }

    wavelet_tree_prefix_sum& operator=(const wavelet_tree_prefix_sum& other)
    {
        if (this == &other) return *this;
        num_entries = other.num_entries; total = other.total; sigma = other.sigma;
        root = other.root; height = other.height;
        _bits = other._bits; _rank = other._rank; _weight_prefix_sums = other._weight_prefix_sums; _nodes = other._nodes;
        rebind(); // the copied rank supports point to the bit vectors of other
        return *this;
    }

    wavelet_tree_prefix_sum& operator=(wavelet_tree_prefix_sum&& other) noexcept
    {
        if (this == &other) return *this;
        num_entries = other.num_entries; total = other.total; sigma = other.sigma;
        root = other.root; height = other.height;
        _bits = std::move(other._bits); _rank = std::move(other._rank);
        _weight_prefix_sums = std::move(other._weight_prefix_sums); _nodes = std::move(other._nodes);
        rebind();
        return *this;
    }

    /**
     * @brief builds the wavelet tree over a sequence of symbols with weights
     * @param symbols symbols of the entries (overwritten during the construction)
     * @param weights weights of the entries (overwritten during the construction)
     * @param sigma the largest symbol + 1
     * @param num_threads number of threads to use for building the prefix sums
     */
    wavelet_tree_prefix_sum(no_init_vector<uint32_t>&& symbols, no_init_vector<uint64_t>&& weights,
        uint64_t sigma, uint64_t num_threads = 1)
        : num_entries(symbols.size()), sigma(sigma)
    {
        if (num_entries == 0 || sigma <= 1) return;

        // cum[c] = number of entries with a symbol < c
        std::vector<uint64_t> cum(sigma + 1, 0);

        for (uint64_t i = 0; i < num_entries; i++) {
            total += weights[i];
            cum[symbols[i] + 1]++;
        }

        for (uint64_t c = 0; c < sigma; c++) cum[c + 1] += cum[c];
        uint8_t width_sym = uint8_t(std::bit_width(sigma));
        uint8_t width_pos = uint8_t(std::bit_width(num_entries));
        _nodes = interleaved_bit_aligned_vectors(
            { width_sym, width_sym, width_pos, width_pos, 0, 0, 0, 0 });
        _nodes.resize(sigma);
        root = build_subtree(cum, 0, sigma);
        if (root == leaf) return; // the alphabet range is not split

        // symbols, weights and prefix sums of the current level after partitioning (written by each level
        // before they are read)
        no_init_vector<uint32_t> partitioned_symbols;
        no_init_vector<uint64_t> partitioned_weights;
        no_init_vector<uint64_t> prefix_sums;
        resize_no_init(partitioned_symbols, num_entries);
        resize_no_init(partitioned_weights, num_entries);
        resize_no_init(prefix_sums, num_entries + 1);
        // (split symbol, starting position, size) of the blocks of the internal nodes of the current
        // level, in symbol order
        std::vector<std::tuple<uint32_t, uint64_t, uint64_t>> blocks { { root, 0, num_entries } };
        std::vector<std::tuple<uint32_t, uint64_t, uint64_t>> next_blocks;
        // block_zeros[s] = number of entries in the block of node s that belong to its left child
        std::vector<uint64_t> block_zeros(sigma, 0);
        uint64_t level_size = num_entries;

        // build the levels from top to bottom
        while (!blocks.empty()) {
            sdsl::bit_vector bits(level_size, 0);
            uint64_t num_ones = 0;

            for (const auto& [split, start, len] : blocks) {
                _nodes.set<field_start>(split, start);
                _nodes.set<field_ones>(split, num_ones);
                uint64_t zeros = 0;

                for (uint64_t i = start; i < start + len; i++) {
                    if (symbols[i] < split) zeros++;
                }

                block_zeros[split] = zeros;
                uint64_t z = start;
                uint64_t o = start + zeros;

                // stably partition the entries of the block into those of the left and of the right child
                for (uint64_t i = start; i < start + len; i++) {
                    uint64_t pos;

                    if (symbols[i] < split) {
                        pos = z++;
                    } else {
                        pos = o++;
                        bits[i] = 1;
                        num_ones++;
                    }

                    partitioned_symbols[pos] = symbols[i];
                    partitioned_weights[pos] = weights[i];
                }
            }

            _bits.emplace_back(std::move(bits));
            _rank.emplace_back(sdsl::rank_support_v<>(&_bits.back()));

            // prefix sums of the weights of the level (the universe is the sum of the weights of the level + 1)
            prefix_sums[0] = 0;
            for (uint64_t k = 0; k < level_size; k++) prefix_sums[k + 1] = prefix_sums[k] + partitioned_weights[k];
            _weight_prefix_sums.emplace_back(elias_fano_bit_vector(
                [&](uint64_t k) { return prefix_sums[k]; }, level_size + 1, prefix_sums[level_size] + 1,
                num_threads));

            // the next level consists of the blocks of the children that are not leaves, in the same order
            next_blocks.clear();
            uint64_t write_pos = 0;

            for (const auto& [split, start, len] : blocks) {
                uint32_t left = uint32_t(_nodes.get<field_left>(split));
                uint32_t right = uint32_t(_nodes.get<field_right>(split));
                uint64_t zeros = block_zeros[split];

                if (left != leaf) {
                    next_blocks.emplace_back(left, write_pos, zeros);
                    for (uint64_t i = start; i < start + zeros; i++) {
                        symbols[write_pos] = partitioned_symbols[i];
                        weights[write_pos] = partitioned_weights[i];
                        write_pos++;
                    }
                }

                if (right != leaf) {
                    next_blocks.emplace_back(right, write_pos, len - zeros);
                    for (uint64_t i = start + zeros; i < start + len; i++) {
                        symbols[write_pos] = partitioned_symbols[i];
                        weights[write_pos] = partitioned_weights[i];
                        write_pos++;
                    }
                }
            }

            level_size = write_pos;
            blocks.swap(next_blocks);
        }

        height = uint32_t(_bits.size());
        rebind();
    }

    // ############################# MISC PUBLIC METHODS #############################

    // returns whether no level is stored (smaller() returns 0 then)
    inline bool empty() const { return height == 0; }

    // returns the number of levels
    inline uint32_t num_levels() const { return height; }

    // returns the size of the bit vectors and their rank supports in bytes
    uint64_t size_bits_in_bytes() const
    {
        uint64_t size = 0;
        for (uint32_t l = 0; l < height; l++)
            size += sdsl::size_in_bytes(_bits[l]) + sdsl::size_in_bytes(_rank[l]);
        return size;
    }

    // returns the size of the prefix sums in bytes
    uint64_t size_prefix_sums_in_bytes() const
    {
        uint64_t size = 0;
        for (uint32_t l = 0; l < height; l++) size += _weight_prefix_sums[l].size_in_bytes();
        return size;
    }

    // returns the size of the nodes in bytes
    uint64_t size_tree_in_bytes() const { return _nodes.size_in_bytes(); }

    // returns the size of the data structure in bytes
    uint64_t size_in_bytes() const
    {
        return sizeof(*this) + size_bits_in_bytes() + size_prefix_sums_in_bytes() + size_tree_in_bytes();
    }

    // ############################# QUERY METHODS #############################

    /**
     * @brief returns the sum of the weights of the entries in [i, j) whose symbol is smaller than c
     * @param i first entry of the range
     * @param j last entry of the range + 1
     * @param c [0..sigma-1] a symbol
     */
    uint64_t smaller(uint64_t i, uint64_t j, uint64_t c) const
    {
        if (i >= j || c == 0 || height == 0) [[unlikely]] return 0;
        assert(c < sigma);

        uint64_t sum = 0;
        uint32_t split = root;

        for (uint32_t level = 0; level < height; level++) {
            uint64_t start = _nodes.get<field_start>(split);
            uint64_t ones_before = _nodes.get<field_ones>(split);
            // number of zeros in the block before i resp. before j
            uint64_t z_i = (i - start) - (_rank[level](i) - ones_before);
            uint64_t z_j = (j - start) - (_rank[level](j) - ones_before);

            if (c >= split) {
                // the entries of the left child have a symbol < split <= c: add their weights
                sum += weight_prefix(level, start + z_j) - weight_prefix(level, start + z_i);
                uint32_t right = uint32_t(_nodes.get<field_right>(split));
                if (right == leaf) break;

                // map i and j to the block of the right child
                uint64_t child_start = _nodes.get<field_start>(right);
                i = child_start + ((i - start) - z_i);
                j = child_start + ((j - start) - z_j);
                split = right;
            } else {
                // map i and j to the block of the left child
                uint32_t left = uint32_t(_nodes.get<field_left>(split));
                if (left == leaf) break;
                uint64_t child_start = _nodes.get<field_start>(left);
                i = child_start + z_i;
                j = child_start + z_j;
                split = left;
            }
        }

        return sum;
    }

    // ############################# SERIALIZATION METHODS #############################

    /**
     * @brief writes the data structure to an output stream
     * @param out output stream
     */
    void serialize(std::ostream& out) const
    {
        out.write((char*) &num_entries, sizeof(uint64_t));
        out.write((char*) &total, sizeof(uint64_t));
        out.write((char*) &sigma, sizeof(uint64_t));
        out.write((char*) &root, sizeof(uint32_t));
        out.write((char*) &height, sizeof(uint32_t));

        _nodes.serialize(out);

        for (uint32_t l = 0; l < height; l++) {
            _bits[l].serialize(out);
            _rank[l].serialize(out);
            _weight_prefix_sums[l].serialize(out);
        }
    }

    /**
     * @brief reads the data structure from an input stream
     * @param in input stream
     */
    void load(std::istream& in)
    {
        in.read((char*) &num_entries, sizeof(uint64_t));
        in.read((char*) &total, sizeof(uint64_t));
        in.read((char*) &sigma, sizeof(uint64_t));
        in.read((char*) &root, sizeof(uint32_t));
        in.read((char*) &height, sizeof(uint32_t));

        _nodes.load(in);

        _bits.resize(height);
        _rank.resize(height);
        _weight_prefix_sums.resize(height);

        for (uint32_t l = 0; l < height; l++) {
            _bits[l].load(in);
            _rank[l].load(in, &_bits[l]);
            _weight_prefix_sums[l].load(in);
        }

        advise_levels();
    }
};

#endif // MOVE_R_DATA_STRUCTURES_WAVELET_TREE_PREFIX_SUM_HPP
