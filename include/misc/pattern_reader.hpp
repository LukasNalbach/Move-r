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
 * @file pattern_reader.hpp
 * @brief reads the patterns of a command-line tool from a file or from stdin, and parses the header of a pizza&chili
 *        patterns file
 */

#ifndef MOVE_R_MISC_PATTERN_READER_HPP
#define MOVE_R_MISC_PATTERN_READER_HPP

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <misc/files.hpp>

// throws the error about a malformed header of a pizza&chili patterns file
[[noreturn]] inline void throw_header_error()
{
    throw std::runtime_error("malformed header in the patterns file (for the format, see "
        "http://pizzachili.dcc.uchile.cl/experiments.html)");
}

/**
 * @brief parses the integer value following a key in a header line (e.g. "number=" in "number=100 length=10");
 *        the value consists of the decimal digits directly after the key
 * @param str the header line
 * @param key the key to look for
 * @return the value following the key, or -1 if the key is not found, has no value or the value is larger
 *         than 2^63 - 1
 */
inline int64_t header_value(const std::string& str, const std::string& key)
{
    uint64_t start_pos = str.find(key);

    if (start_pos == std::string::npos || start_pos + key.size() >= str.size()) {
        return -1;
    }

    start_pos += key.size();
    uint64_t value = 0;
    const std::from_chars_result result = std::from_chars(str.data() + start_pos, str.data() + str.size(), value);
    if (result.ec != std::errc() || value > uint64_t(std::numeric_limits<int64_t>::max())) return -1;
    return int64_t(value);
}

/**
 * @brief reads the number of patterns from a patterns-file header
 * @param header the header line
 * @return the number of patterns (throws std::runtime_error if the header is malformed)
 */
inline uint64_t num_patterns_from_header(std::string header)
{
    int64_t num = header_value(header, "number=");

    if (num < 1) {
        throw_header_error();
    }

    return num;
}

/**
 * @brief reads the pattern length from a patterns-file header
 * @param header the header line
 * @return the pattern length (throws std::runtime_error if the header is malformed)
 */
inline uint64_t pattern_length_from_header(std::string header)
{
    int64_t len = header_value(header, "length=");

    if (len < 1) {
        throw_header_error();
    }

    return len;
}

/**
 * @brief reads and parses the header line of a pizza&chili patterns file
 * @param patterns_file a patterns file positioned at its header
 * @return the number of patterns and their length
 */
inline std::pair<uint64_t, uint64_t> read_patterns_header(std::istream& patterns_file)
{
    std::string header;
    std::getline(patterns_file, header);
    return { num_patterns_from_header(header), pattern_length_from_header(header) };
}

// the description of a patterns file in the help texts of the command-line tools (see pattern_reader)
inline const std::string patterns_file_help = "the patterns, in pizza&chili, FASTA or FASTQ format or one per line "
    "(detected from the input; a file of one pattern per line whose first pattern starts with '>' is read as FASTA)";

/**
 * @brief reads patterns in one of four formats, which is detected from the input: pizza&chili (a first line
 *        "# number=<k> length=<m> ..." followed by k patterns of length m, which a file has to hold), FASTQ (a first
 *        record of four lines: '@' and the name, the symbols, '+', as many qualities as symbols), FASTA (the input
 *        starts with '>') or one pattern per line (otherwise; empty lines are skipped). The lines read to detect the
 *        format are buffered, so the input may be stdin. A pattern without a name in its format is named pat<i>,
 *        where i is its index.
 */
class pattern_reader {
public:
    // format of the input
    enum format_t : uint8_t {
        _pizza_chili,
        _fasta,
        _fastq,
        _lines
    };

    // a pattern
    struct pattern_t {
        std::string name; // the name (up to the first whitespace of the header in FASTA and FASTQ)
        std::string seq; // the symbols
        std::string quality; // the qualities (FASTQ only)
    };

    /**
     * @brief opens the input
     * @param path the file to read from; "" or "-" reads from stdin
     */
    explicit pattern_reader(const std::string& path)
    {
        if (path.empty() || path == "-") {
            in = &std::cin;
        } else {
            file = std::make_unique<std::ifstream>(path, std::ios::binary);
            failed = !file->good();
            in = file.get();
            std::error_code ec;
            if (!failed) size = std::filesystem::file_size(path, ec);
            if (ec) size = 0;
        }

        if (failed) return;
        const int first = in->peek();

        if (first == '#') {
            std::string header;
            std::getline(*in, header);
            const int64_t num = header_value(header, "number=");
            const int64_t len = header_value(header, "length=");
            // the bytes the header line and the patterns take in a file
            const uint64_t pattern_bytes = header.size() + 1
                + uint64_t(std::max<int64_t>(num, 0)) * uint64_t(std::max<int64_t>(len, 0));

            if (num >= 1 && len >= 1 && (file == nullptr || size >= pattern_bytes)) {
                fmt = _pizza_chili;
                count = num;
                length = len;
            } else {
                buffered.push_back(header);
            }
        } else if (first == '>') {
            fmt = _fasta;
        } else if (first == '@') {
            // the four lines of a FASTQ record: header, symbols, separator, qualities
            constexpr uint64_t record_lines = 4;
            std::string l;
            while (buffered.size() < record_lines && std::getline(*in, l)) buffered.push_back(l);
            for (std::string& b : buffered) chomp(b);
            const bool fastq = buffered.size() == record_lines && !buffered[2].empty() && buffered[2][0] == '+'
                && buffered[1].size() == buffered[3].size();
            if (fastq) fmt = _fastq;
        }
    }

    // returns true <=> the input could be opened
    inline bool good() const { return in != nullptr && !failed; }

    // returns the format of the input
    inline format_t format() const { return fmt; }

    // returns the number of patterns, if the format states it (pizza&chili), else 0
    inline uint64_t num_patterns() const { return count; }

    // returns the length of the patterns, if the format states it (pizza&chili), else 0
    inline uint64_t pattern_length() const { return length; }

    // returns the size of the input in bytes, if it is a file, else 0
    inline uint64_t file_bytes() const { return size; }

    // returns the number of bytes of the input read so far, if it is a file, else 0
    inline uint64_t bytes_read()
    {
        if (file == nullptr) return 0;
        const std::streamoff pos = file->tellg();
        return pos < 0 ? size : uint64_t(pos);
    }

    /**
     * @brief reads the next pattern
     * @param p is set to the pattern
     * @return false <=> there is no pattern left
     */
    bool next(pattern_t& p)
    {
        p.quality.clear();

        // returns the name in a FASTA or FASTQ header (without its first symbol, up to the first whitespace)
        auto name_of = [](const std::string& header) {
            const uint64_t end = header.find_first_of(" \t\r", 1);
            return header.substr(1, end == std::string::npos ? std::string::npos : end - 1);
        };

        if (fmt == _pizza_chili) {
            if (index == count) return false;
            p.seq.resize(length);
            in->read(p.seq.data(), length);
            if (uint64_t(in->gcount()) != length) return false;
            p.name = "pat" + std::to_string(index++);
            return true;
        }

        if (fmt == _fasta) {
            if (line.empty() && !read_line(line)) return false;
            chomp(line);
            p.name = name_of(line);
            p.seq.clear();
            line.clear();

            while (read_line(line)) {
                chomp(line);
                if (!line.empty() && line[0] == '>') break;
                p.seq += line;
                line.clear();
            }

            index++;
            return true;
        }

        if (fmt == _fastq) {
            std::string header, plus;
            if (!read_line(header)) return false;
            chomp(header);
            if (header.empty()) return false;
            p.name = name_of(header);
            read_line(p.seq);
            read_line(plus);
            read_line(p.quality);
            chomp(p.seq);
            chomp(p.quality);
            index++;
            return true;
        }

        while (read_line(p.seq)) {
            chomp(p.seq);
            if (p.seq.empty()) continue;
            p.name = "pat" + std::to_string(index++);
            return true;
        }

        return false;
    }

protected:
    // removes a trailing carriage return
    static void chomp(std::string& s)
    {
        if (!s.empty() && s.back() == '\r') s.pop_back();
    }

    // reads the next line: first the buffered lines, then from the input; returns false at the end of the input
    bool read_line(std::string& l)
    {
        if (next_buffered < buffered.size()) {
            l = std::move(buffered[next_buffered++]);
            return true;
        }

        return bool(std::getline(*in, l));
    }

    std::unique_ptr<std::ifstream> file; // the file, unless the input is stdin
    std::istream* in = nullptr; // the input
    format_t fmt = _lines;
    bool failed = false; // true <=> the input could not be opened
    uint64_t count = 0; // the number of patterns (pizza&chili)
    uint64_t length = 0; // the length of the patterns (pizza&chili)
    uint64_t size = 0; // the size of the input in bytes, if it is a file, else 0
    uint64_t index = 0; // the index of the next pattern
    std::string line; // the next line of a FASTA file (its header), if it has been read already
    std::vector<std::string> buffered; // the lines read to detect the format, which are read again by next()
    uint64_t next_buffered = 0; // the index in buffered of the next line read by next()
};

/**
 * @brief reads all patterns of a file (see pattern_reader)
 * @param path the file to read from; "" or "-" reads from stdin
 * @param patterns is set to the patterns
 * @return false <=> the file could not be read
 */
inline bool read_all_patterns(const std::string& path, std::vector<pattern_reader::pattern_t>& patterns)
{
    pattern_reader reader(path);
    if (!reader.good()) return false;
    patterns.clear();
    pattern_reader::pattern_t pattern;
    while (reader.next(pattern)) patterns.push_back(pattern);
    return true;
}

#endif // MOVE_R_MISC_PATTERN_READER_HPP
