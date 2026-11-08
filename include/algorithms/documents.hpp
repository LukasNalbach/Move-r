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
 * @file documents.hpp
 * @brief document listing: the documents containing a pattern, from its located occurrences
 */

#ifndef MOVE_R_ALGORITHMS_DOCUMENTS_HPP
#define MOVE_R_ALGORITHMS_DOCUMENTS_HPP

#include <algorithm>
#include <cstdint>
#include <data_structures/document_array.hpp>
#include <utility>
#include <vector>

/**
 * @brief finds the documents containing a pattern from its occurrences: from a source of the occurrences, a function
 *        called as source(report), which calls report(occ) for the text position occ of every occurrence (e.g., with
 *        the occurrences located by a query context of the index)
 */
class document_query {
protected:
    const document_array& docs; // the documents of the text
    // the documents of the occurrences of the last query: [0..number of occurrences - 1] in increasing order, without
    // duplicates after decode_unique
    std::vector<uint64_t> decoded_docs;
    // the documents of top_k: [0..number of documents - 1] (frequency, document); the first k are sorted by decreasing
    // frequency (ties by increasing document)
    std::vector<std::pair<uint64_t, uint64_t>> frequencies;

    // writes the documents of the occurrences of a source into decoded_docs and sorts them
    template <typename source_t>
    void decode_sorted(source_t& source)
    {
        decoded_docs.clear();
        source([&](uint64_t occ) { decoded_docs.push_back(docs.documents().document_of(occ)); });
        std::sort(decoded_docs.begin(), decoded_docs.end());
    }

    // writes the documents of the occurrences of a source into decoded_docs and removes the duplicates
    template <typename source_t>
    void decode_unique(source_t& source)
    {
        decode_sorted(source);
        decoded_docs.erase(std::unique(decoded_docs.begin(), decoded_docs.end()), decoded_docs.end());
    }

public:
    /**
     * @brief prepares the queries
     * @param docs the documents of the text
     */
    explicit document_query(const document_array& docs) : docs(docs) {}

    /**
     * @brief reports the documents containing a pattern, each once, in increasing order
     * @tparam source_t type of source
     * @tparam fnc_t type of report
     * @param source the source of the occurrences of the pattern (see document_query)
     * @param report function called as report(document)
     */
    template <typename source_t, typename fnc_t>
    void list(source_t source, fnc_t report)
    {
        decode_unique(source);
        for (uint64_t document : decoded_docs) report(document);
    }

    /**
     * @brief returns the number of documents containing a pattern
     * @tparam source_t type of source
     * @param source the source of the occurrences of the pattern (see document_query)
     */
    template <typename source_t>
    uint64_t count(source_t source)
    {
        decode_unique(source);
        return decoded_docs.size();
    }

    /**
     * @brief reports the documents containing a pattern with the numbers of occurrences of the pattern in them, in
     *        increasing order of the documents
     * @tparam source_t type of source
     * @tparam fnc_t type of report
     * @param source the source of the occurrences of the pattern (see document_query)
     * @param report function called as report(document, frequency)
     */
    template <typename source_t, typename fnc_t>
    void list_with_frequencies(source_t source, fnc_t report)
    {
        decode_sorted(source);

        for (uint64_t i = 0; i < decoded_docs.size();) {
            uint64_t j = i + 1;
            while (j < decoded_docs.size() && decoded_docs[j] == decoded_docs[i]) j++;
            report(decoded_docs[i], j - i);
            i = j;
        }
    }

    /**
     * @brief reports the k documents with the most occurrences of a pattern (ties are broken by the smaller
     *        document), in decreasing order of their frequencies
     * @tparam source_t type of source
     * @tparam fnc_t type of report
     * @param source the source of the occurrences of the pattern (see document_query)
     * @param k the number of documents
     * @param report function called as report(document, frequency)
     */
    template <typename source_t, typename fnc_t>
    void top_k(source_t source, uint64_t k, fnc_t report)
    {
        frequencies.clear();
        list_with_frequencies(source, [&](uint64_t document, uint64_t frequency) {
            frequencies.emplace_back(frequency, document); });
        k = std::min<uint64_t>(k, frequencies.size());

        std::partial_sort(frequencies.begin(), frequencies.begin() + k, frequencies.end(),
            [](const auto& x, const auto& y) { return x.first != y.first ? x.first > y.first : x.second < y.second; });

        for (uint64_t i = 0; i < k; i++) report(frequencies[i].second, frequencies[i].first);
    }
};

#endif // MOVE_R_ALGORITHMS_DOCUMENTS_HPP
