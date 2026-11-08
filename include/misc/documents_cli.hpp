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
 * @file documents_cli.hpp
 * @brief the option -docs of move-rb build: the documents of a FASTA input (see document_map)
 */

#ifndef MOVE_R_MISC_DOCUMENTS_CLI_HPP
#define MOVE_R_MISC_DOCUMENTS_CLI_HPP

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <data_structures/document_map.hpp>
#include <data_structures/fasta_sequence_data.hpp>
#include <misc/log.hpp>

/**
 * @brief returns the documents of a FASTA input whose records are assigned to documents; record k covers the
 *        text positions from the separator before it to the position before the separator of record k + 1
 * @param data the sequences (records) of the input
 * @param record_document the document of each record: [0..data.num_sequences()-1] record_document[k] is the document
 *        of record k
 * @param names the names of the documents
 */
inline document_map document_map_of_records(const fasta_sequence_data& data,
    const std::vector<uint32_t>& record_document, std::vector<std::string> names)
{
    // the ranges of consecutive records of the same document: starts[x] is the text position at which range x
    // starts (starts[0] = 0, increasing), documents[x] is its document
    std::vector<uint64_t> starts;
    std::vector<uint32_t> documents;

    for (uint64_t k = 0; k < data.num_sequences(); k++) {
        if (!documents.empty() && documents.back() == record_document[k]) continue;
        starts.push_back(k == 0 ? 0 : data.sequence_start(k) - 1);
        documents.push_back(record_document[k]);
    }

    return document_map(std::move(starts), std::move(documents), std::move(names));
}

/**
 * @brief returns the documents of a FASTA input given with the option -docs
 * @param source "records" (one document per record), "files" (one document per input file) or the path of a file
 *        with one line "<record name> <document name>" per record (the documents are numbered in the order of their
 *        first lines)
 * @param data the sequences (records) of the input
 * @param input_files the input files, in the order in which they are concatenated
 */
inline document_map document_map_of_cli(const std::string& source, const fasta_sequence_data& data,
    const std::vector<std::string>& input_files)
{
    const uint64_t num_records = data.num_sequences();
    if (num_records == 0) throw std::invalid_argument("-docs needs a FASTA input with at least one record");
    // the document of each record: [0..num_records-1] record_document[k] is the document of record k
    std::vector<uint32_t> record_document(num_records);
    // [0..d-1] the name of document x
    std::vector<std::string> names;

    if (source == "records") {
        document_map::check_num_ids(num_records, "documents");

        for (uint64_t k = 0; k < num_records; k++) {
            record_document[k] = uint32_t(k);
            names.push_back(data.sequence_name(k));
        }
    } else if (source == "files") {
        document_map::check_num_ids(input_files.size(), "documents");
        uint64_t k = 0;

        for (const std::string& file : input_files) {
            std::ifstream in(file);
            std::string line;
            const uint32_t document = uint32_t(names.size());
            names.push_back(std::filesystem::path(file).filename().string());

            while (std::getline(in, line)) {
                if (line.empty() || line[0] != '>') continue;
                if (k == num_records) throw std::invalid_argument("-docs files: the files have more records than the text");
                record_document[k++] = document;
            }
        }

        if (k != num_records) throw std::invalid_argument("-docs files: the files have fewer records than the text");
    } else {
        std::ifstream in(source);
        if (!in.good()) throw std::runtime_error("cannot read " + source);
        std::unordered_map<std::string, uint32_t> document_of_record;
        std::unordered_map<std::string, uint32_t> document_of_name;
        std::string line;

        while (std::getline(in, line)) {
            std::istringstream fields(line);
            std::string record;
            std::string document;
            if (!(fields >> record)) continue;
            if (!(fields >> document)) throw std::invalid_argument("-docs: the line '" + line + "' has no document");
            auto [it, added] = document_of_name.try_emplace(document, uint32_t(names.size()));

            if (added) {
                document_map::check_num_ids(names.size() + 1, "documents");
                names.push_back(document);
            }

            document_of_record[record] = it->second;
        }

        for (uint64_t k = 0; k < num_records; k++) {
            auto it = document_of_record.find(data.sequence_name(k));
            if (it == document_of_record.end())
                throw std::invalid_argument("-docs: the record " + data.sequence_name(k) + " has no document");
            record_document[k] = it->second;
        }
    }

    return document_map_of_records(data, record_document, std::move(names));
}

/**
 * @brief returns the documents of a FASTA input given with the option -docs (see document_map_of_cli) and logs their
 *        number
 * @param source the value of -docs
 * @param data the sequences (records) of the input
 * @param input_files the input files, in the order in which they are concatenated
 */
inline document_map read_documents_option(const std::string& source, const fasta_sequence_data& data,
    const std::vector<std::string>& input_files)
{
    document_map documents = document_map_of_cli(source, data, input_files);
    log_out() << format_count(documents.num_documents(), "document") << " in "
              << format_count(documents.num_ranges(), "range") << std::endl;
    return documents;
}

#endif // MOVE_R_MISC_DOCUMENTS_CLI_HPP
