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
 * @file index_cli.hpp
 * @brief functions shared by the command-line tools and benchmarks of the indexes: loading an index with a progress
 *        display, printing the information about an index, closing output files and formatting query throughputs;
 *        includes the measurement file and RESULT lines (misc/measurement.hpp)
 */

#ifndef MOVE_R_MISC_INDEX_CLI_HPP
#define MOVE_R_MISC_INDEX_CLI_HPP

#include <cstdint>
#include <iostream>
#include <streambuf>
#include <string>

#include <misc/files.hpp>
#include <misc/index_facts.hpp>
#include <misc/log.hpp>
#include <misc/measurement.hpp>

/**
 * @brief read-only stream buffer that reads from another stream buffer and counts the bytes read as the
 *        progress of the running phase (see phase_progress)
 */
class index_load_progress_streambuf : public std::streambuf {
    std::streambuf* source; // the stream buffer to read from
    phase_progress progress; // number of bytes read

protected:
    // reads count bytes from the underlying buffer and advances the progress by the number of bytes read
    std::streamsize xsgetn(char* data, std::streamsize count) override
    {
        std::streamsize read = source->sgetn(data, count);
        if (read > 0) progress.advance(static_cast<uint64_t>(read));
        return read;
    }

    int_type underflow() override { return source->sgetc(); }

    // reads one character and advances the progress by 1
    int_type uflow() override
    {
        int_type value = source->sbumpc();
        if (!traits_type::eq_int_type(value, traits_type::eof())) progress.advance();
        return value;
    }

    std::streamsize showmanyc() override { return source->in_avail(); }

    // seeks in the underlying buffer; the progress is not changed
    pos_type seekoff(
        off_type offset, std::ios_base::seekdir direction,
        std::ios_base::openmode mode = std::ios_base::in | std::ios_base::out) override
    {
        return source->pubseekoff(offset, direction, mode);
    }

    // seeks in the underlying buffer; the progress is not changed
    pos_type seekpos(
        pos_type position,
        std::ios_base::openmode mode = std::ios_base::in | std::ios_base::out) override
    {
        return source->pubseekpos(position, mode);
    }

public:
    // total: total number of bytes that are read from source
    index_load_progress_streambuf(std::streambuf* source, uint64_t total)
        : source(source), progress(true, total)
    {
    }

    // ends the progress output (see phase_progress::finish)
    void finish() { progress.finish(); }
};

/**
 * @brief loads an index as a logged phase: prints "loading the <name>", the progress of reading its bytes with the
 *        elapsed and the remaining time, and the time it took
 * @tparam index_t type of the index
 * @param index the index to load
 * @param in input stream, positioned at the beginning of the serialized index
 * @param name name of the index, e.g. "index"
 * @return the size of the serialized index in bytes (from the position of in to the end of the stream)
 */
template <typename index_t>
uint64_t load_index_with_progress(index_t& index, std::istream& in, const std::string& name = "index")
{
    std::streampos begin = in.tellg();
    in.seekg(0, std::ios::end);
    std::streampos end = in.tellg();
    in.seekg(begin);

    // number of bytes of the serialized index
    uint64_t total = begin >= 0 && end >= begin
        ? static_cast<uint64_t>(end - begin) : 0;

    auto time = now();
    log_phase_start(true, time, "loading the " + name);
    index_load_progress_streambuf buffer(in.rdbuf(), total);
    std::istream monitored(&buffer);

    try {
        index.load(monitored);
    } catch (...) {
        buffer.finish();
        end_phase_line();
        throw;
    }

    buffer.finish();
    log_phase_end(true, time);
    return total;
}

/**
 * @brief returns the description of the alphabet of an index for the facts of the info modes
 * @param byte_alphabet whether the symbols of the text are bytes
 * @param input_order whether the symbols are ordered as in the input (else by decreasing frequency)
 */
inline std::string alphabet_description(bool byte_alphabet, bool input_order)
{
    return std::string("alphabet: ") + (byte_alphabet ? "bytes" : "integers")
        + (input_order ? ", in the order of the input" : ", by decreasing frequency");
}

/**
 * @brief prints the information about an index (the output of the info modes) after its loading line: the facts,
 *        one per line, the size of the index ("<name> size: ..."), the sizes of its data structures and the modes that
 *        support it
 * @tparam log_sizes_fnc_t type of log_sizes
 * @param facts the facts
 * @param name the name of the index in the size line, e.g. "index"
 * @param index_bytes the size of the serialized index in bytes
 * @param input_bytes the size of the input in bytes (with the terminator)
 * @param log_sizes log_sizes() logs the sizes of the data structures
 * @param modes the modes, separated by commas
 */
template <typename log_sizes_fnc_t>
void print_index_info(const index_facts& facts, const std::string& name, uint64_t index_bytes, uint64_t input_bytes,
    log_sizes_fnc_t log_sizes, const std::string& modes)
{
    facts.print(std::cout);
    log_output_size(name, index_bytes, input_bytes);
    log_sizes();
    std::cout << "modes: " << modes << std::endl;
}

/**
 * @brief closes an output file and drops it from the page cache
 * @param file the stream to close (nothing happens if it is not open)
 * @param path the path the stream was opened with
 */
inline void close_output_file(direct_ofstream& file, const std::string& path)
{
    if (!file.is_open()) return;
    file.close();
    drop_from_page_cache(path, true);
}

/**
 * @brief formats a query throughput (queries over a duration) with an appropriate unit
 * @param num_queries the number of queries
 * @param ns the duration in nanoseconds
 * @return the formatted query throughput
 */
inline std::string format_query_throughput(uint64_t num_queries, uint64_t ns)
{
    std::string str;
    double queries_per_ns = num_queries / (double)ns;

    if (queries_per_ns < 0.000001) {
        str = std::to_string(queries_per_ns * 1000000000) + " queries/s";
    } else if (queries_per_ns < 0.001) {
        str = std::to_string(queries_per_ns * 1000000) + " queries/ms";
    } else if (queries_per_ns < 1) {
        str = std::to_string(queries_per_ns * 1000) + " queries/us";
    } else {
        str = std::to_string(queries_per_ns) + " queries/ns";
    }

    return str;
}

#endif // MOVE_R_MISC_INDEX_CLI_HPP
