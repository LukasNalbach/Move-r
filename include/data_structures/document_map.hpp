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

/**
 * @file document_map.hpp
 * @brief the documents of a text as consecutive ranges of positions
 */

#ifndef MOVE_R_DATA_STRUCTURES_DOCUMENT_MAP_HPP
#define MOVE_R_DATA_STRUCTURES_DOCUMENT_MAP_HPP

#include <algorithm>
#include <bit>
#include <cstdint>
#include <istream>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <misc/files.hpp>

/**
 * @brief the documents of a text: the text is split into consecutive ranges of positions, and each range
 *        belongs to a document (a document can consist of several ranges); the range of a position is found by a
 *        binary search over the ranges that start in its bucket of 2^_bucket_shift positions
 */
class document_map {
public:
    // ############################# CONSTANTS #############################

    // the largest number of documents and of ranges: their ids are stored as uint32_t
    static constexpr uint64_t max_num_ids = std::numeric_limits<uint32_t>::max();

    /**
     * @brief throws std::length_error if more than max_num_ids documents or ranges would be created
     * @param num the number of documents or ranges
     * @param what what is counted, e.g. "documents"
     */
    static void check_num_ids(uint64_t num, const std::string& what)
    {
        if (num > max_num_ids) {
            throw std::length_error("too many " + what + ": " + std::to_string(num) + ", at most 2^"
                + std::to_string(8 * sizeof(uint32_t)) + " - 1 are supported");
        }
    }

protected:
    // ############################# VARIABLES #############################

    // the starting positions of the ranges: [0..num_ranges()-1] _starts[k] is the first position of range k; _starts[0] = 0, strictly increasing
    std::vector<uint64_t> _starts;

    // the documents of the ranges: [0..num_ranges()-1] _documents[k] is the document of range k
    std::vector<uint32_t> _documents;

    // the names of the documents: [0..num_documents()-1]
    std::vector<std::string> _names;

    // the text positions are split into buckets of 2^_bucket_shift positions, about as many as there are ranges
    uint8_t _bucket_shift = 0;

    /* the first range of each bucket: [0..number of buckets] _bucket_first[x] is the index of the range containing
       the first position of bucket x; the last entry is the index of the last range */
    std::vector<uint32_t> _bucket_first;

    // ############################# INTERNAL METHODS #############################

    // builds _bucket_shift and _bucket_first
    void build_buckets()
    {
        const uint64_t universe = _starts.back() + 1;
        _bucket_shift = std::bit_width(std::max<uint64_t>(1, universe / _starts.size())) - 1;
        const uint64_t num_buckets = ((universe - 1) >> _bucket_shift) + 1;
        _bucket_first.resize(num_buckets + 1);
        uint64_t k = 0;

        for (uint64_t x = 0; x < num_buckets; x++) {
            while (k + 1 < _starts.size() && _starts[k + 1] <= (x << _bucket_shift)) k++;
            _bucket_first[x] = uint32_t(k);
        }

        _bucket_first[num_buckets] = uint32_t(_starts.size() - 1);
    }

public:
    // ############################# CONSTRUCTORS #############################

    document_map() = default;

    /**
     * @brief builds the map
     * @param starts the starting positions of the ranges in the text (strictly increasing, starts[0] = 0)
     * @param documents documents[k] is the document of the k-th range
     * @param names the names of the documents
     */
    document_map(std::vector<uint64_t> starts, std::vector<uint32_t> documents, std::vector<std::string> names)
        : _starts(std::move(starts)), _documents(std::move(documents)), _names(std::move(names))
    {
        if (_starts.empty() || _starts[0] != 0 || _starts.size() != _documents.size())
            throw std::invalid_argument("document_map: the ranges must start at position 0, one document per range");
        check_num_ids(_names.size(), "documents");
        check_num_ids(_starts.size(), "document ranges");

        for (uint64_t k = 0; k < _starts.size(); k++) {
            if (k != 0 && _starts[k] <= _starts[k - 1])
                throw std::invalid_argument("document_map: the starting positions must be strictly increasing");
            if (_documents[k] >= _names.size()) throw std::invalid_argument("document_map: a range has no document");
        }

        build_buckets();
    }

    // ############################# MISC PUBLIC METHODS #############################

    // returns the number of documents
    inline uint64_t num_documents() const { return _names.size(); }

    // returns whether the map is empty
    inline bool empty() const { return _names.empty(); }

    // returns the number of ranges
    inline uint64_t num_ranges() const { return _starts.size(); }

    // returns the size in bytes
    uint64_t size_in_bytes() const
    {
        uint64_t size = _starts.size() * sizeof(uint64_t) + _documents.size() * sizeof(uint32_t)
            + _bucket_first.size() * sizeof(uint32_t);
        for (const std::string& name : _names) size += name.size();
        return size;
    }

    // ############################# ACCESS METHODS #############################

    // returns the name of a document
    inline const std::string& name(uint64_t document) const { return _names[document]; }

    // returns the first position of a range
    inline uint64_t range_start(uint64_t range) const { return _starts[range]; }

    // returns the document of a range
    inline uint64_t range_document(uint64_t range) const { return _documents[range]; }

    // ############################# QUERY METHODS #############################

    /**
     * @brief returns the document of a text position
     * @param position a position in the text
     */
    inline uint64_t document_of(uint64_t position) const { return _documents[range_of(position)]; }

    /**
     * @brief returns the index of the range containing a text position: the range begins in the bucket of the
     *        position or contains the first position of the bucket (see _bucket_first)
     * @param position a position in the text
     */
    inline uint64_t range_of(uint64_t position) const
    {
        const uint64_t bucket = std::min<uint64_t>(position >> _bucket_shift, _bucket_first.size() - 2);
        const auto first = _starts.begin() + _bucket_first[bucket];
        const auto last = _starts.begin() + _bucket_first[bucket + 1] + 1;
        return std::upper_bound(first, last, position) - _starts.begin() - 1;
    }

    // ############################# SERIALIZATION METHODS #############################

    // writes the map to an output stream
    void serialize(std::ostream& out) const
    {
        const uint64_t num_ranges = _starts.size();
        const uint64_t num_names = _names.size();
        out.write((const char*) &num_ranges, sizeof(uint64_t));
        out.write((const char*) _starts.data(), num_ranges * sizeof(uint64_t));
        out.write((const char*) _documents.data(), num_ranges * sizeof(uint32_t));
        out.write((const char*) &num_names, sizeof(uint64_t));

        for (const std::string& name : _names) write_string(out, name);
    }

    // reads a map written by serialize from an input stream
    void load(std::istream& in)
    {
        uint64_t num_ranges = 0;
        uint64_t num_names = 0;
        in.read((char*) &num_ranges, sizeof(uint64_t));
        _starts.resize(num_ranges);
        _documents.resize(num_ranges);
        in.read((char*) _starts.data(), num_ranges * sizeof(uint64_t));
        in.read((char*) _documents.data(), num_ranges * sizeof(uint32_t));
        in.read((char*) &num_names, sizeof(uint64_t));
        _names.resize(num_names);

        for (std::string& name : _names) read_string(in, name);

        if (!_starts.empty()) build_buckets();
    }
};

#endif // MOVE_R_DATA_STRUCTURES_DOCUMENT_MAP_HPP
