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
 * @file document_array.hpp
 * @brief the documents of the text of an index
 */

#ifndef MOVE_R_DATA_STRUCTURES_DOCUMENT_ARRAY_HPP
#define MOVE_R_DATA_STRUCTURES_DOCUMENT_ARRAY_HPP

#include <cstdint>
#include <data_structures/document_map.hpp>
#include <istream>
#include <ostream>
#include <utility>

/**
 * @brief the documents of the text of an index that locates all occurrences (see document_map); the documents of a
 *        pattern are computed from its located occurrences (see document_query)
 */
class document_array {
protected:
    document_map _documents; // the documents of the text

public:
    document_array() = default;

    /**
     * @brief stores the documents of a text
     * @param documents the documents of the text
     */
    explicit document_array(document_map documents) : _documents(std::move(documents)) {}

    // returns whether no documents are stored
    inline bool empty() const { return _documents.empty(); }

    // returns the documents of the text
    inline const document_map& documents() const { return _documents; }

    // returns the number of documents
    inline uint64_t num_documents() const { return _documents.num_documents(); }

    // returns the size in bytes
    inline uint64_t size_in_bytes() const { return _documents.size_in_bytes(); }

    // writes the documents to an output stream
    void serialize(std::ostream& out) const { _documents.serialize(out); }

    // reads documents written by serialize from an input stream
    void load(std::istream& in) { _documents.load(in); }
};

#endif // MOVE_R_DATA_STRUCTURES_DOCUMENT_ARRAY_HPP
