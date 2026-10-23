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
 * @file parallel_patterns.hpp
 * @brief processes the patterns of a file with several threads: the patterns are read in blocks, the threads take
 *        chunks of consecutive patterns of a block, and the output of the patterns is written in their order; while
 *        the threads process a block, one of them first writes the output of the previous block and reads the next
 *        one; prints the summary of the queries of all query modes
 */

#ifndef MOVE_R_MISC_PARALLEL_PATTERNS_HPP
#define MOVE_R_MISC_PARALLEL_PATTERNS_HPP

#include <algorithm>
#include <array>
#include <atomic>
#include <charconv>
#include <concepts>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <omp.h>
#include <ostream>
#include <string>
#include <string_view>
#include <vector>

#include <misc/fasta.hpp>
#include <misc/log.hpp>
#include <misc/pattern_reader.hpp>
#include <misc/threads.hpp>
#include <misc/utils.hpp>

/**
 * @brief returns the default number of threads of the query modes: the number of CPUs the threads of an OpenMP region
 *        may run on (see taskset, OMP_PLACES and OMP_PROC_BIND), at most omp_get_max_threads(); it runs a parallel
 *        region, so it is called after the command line is parsed
 */
inline uint64_t default_query_threads()
{
    const uint64_t max_threads = std::max<int>(1, omp_get_max_threads());
    platform::cpu_mask all;
    for (const platform::cpu_mask& cpus : omp_thread_cpus(max_threads)) all |= cpus;
    return std::clamp<uint64_t>(all.count(), 1, max_threads);
}

/** @brief the options of process_patterns */
struct pattern_processing_t {
    uint64_t threads = 1; // the number of threads
    bool reverse_complements = false; // true <=> each pattern is followed by its reverse complement (named <name>_rc)
    std::vector<std::ostream*> outputs; // the output streams (nullptr: not written); each pattern has one buffer per stream
    bool log = true; // true <=> the search is logged as a phase with its progress (see process_patterns)
    uint64_t chunk_symbols = uint64_t{1} << 20; // the number of symbols from which on a chunk is cut
    // the number of symbols from which on a block after the first is cut (with few threads, its chunks are longer than
    // chunk_symbols)
    uint64_t block_symbols = uint64_t{1} << 24;
    uint64_t chunks_per_thread = 2; // the number of chunks per thread a block is cut into
};

/** @brief the totals of process_patterns */
struct pattern_totals_t {
    uint64_t patterns = 0; // the number of patterns (with their reverse complements)
    uint64_t symbols = 0; // the number of symbols of the patterns
    uint64_t time_wall = 0; // the time (in ns) from reading the first pattern to writing the last output
};

/** @brief consecutive patterns of a block, processed by one thread */
class pattern_chunk {
public:
    pattern_chunk(const pattern_reader::pattern_t* patterns, std::string* buffers, uint64_t size, uint64_t num_outputs)
        : patterns(patterns), buffers(buffers), num(size), num_outputs(num_outputs) {}

    // returns the number of patterns of the chunk
    uint64_t size() const { return num; }

    // returns pattern j of the chunk
    const pattern_reader::pattern_t& pattern(uint64_t j) const { return patterns[j]; }

    // returns the output of pattern j of the chunk to the output stream s
    std::string& out(uint64_t j, uint64_t s = 0) const { return buffers[j * num_outputs + s]; }

private:
    const pattern_reader::pattern_t* patterns; // the patterns of the chunk
    std::string* buffers; // the outputs of the patterns of the chunk, num_outputs per pattern
    uint64_t num; // the number of patterns
    uint64_t num_outputs; // the number of output streams
};

/** @brief appends text and integers to a string, in the way an output stream writes them (symbols as characters) */
class string_appender {
public:
    explicit string_appender(std::string& s) : s(&s) {}

    string_appender& operator<<(char c) { s->push_back(c); return *this; }
    string_appender& operator<<(signed char c) { s->push_back(char(c)); return *this; }
    string_appender& operator<<(unsigned char c) { s->push_back(char(c)); return *this; }
    string_appender& operator<<(std::string_view v) { s->append(v); return *this; }
    string_appender& operator<<(const char* v) { s->append(v); return *this; }
    string_appender& operator<<(const std::string& v) { s->append(v); return *this; }

    template <std::integral int_t>
        requires(!std::same_as<int_t, bool> && !std::same_as<int_t, char> && !std::same_as<int_t, signed char> &&
            !std::same_as<int_t, unsigned char>)
    string_appender& operator<<(int_t v)
    {
        char buf[std::numeric_limits<int_t>::digits10 + 2]; // the digits and the sign
        s->append(buf, std::to_chars(buf, buf + sizeof(buf), v).ptr);
        return *this;
    }

private:
    std::string* s; // the string
};

/** @brief a block of patterns of process_patterns and the output of its patterns */
struct pattern_block_t {
    std::vector<pattern_reader::pattern_t> patterns; // the patterns (the first size ones belong to the block)
    std::vector<std::string> buffers; // the outputs of the patterns, one per pattern and output stream
    std::vector<uint64_t> chunk_starts; // the first pattern of each chunk, followed by size
    uint64_t size = 0; // the number of patterns of the block
    uint64_t symbols = 0; // the number of symbols of the patterns of the block
};

/**
 * @brief reads the next block of patterns and cuts it into chunks
 * @param reader the reader of the patterns
 * @param options the options of process_patterns
 * @param block the block
 * @param first true <=> the block is the first one, which is read before the threads start (it has at least
 *        options.chunk_symbols symbols per chunk, the others at least options.block_symbols in total)
 */
inline void read_pattern_block(pattern_reader& reader, const pattern_processing_t& options, pattern_block_t& block,
    bool first)
{
    const uint64_t chunks = options.threads <= 1 ? 1 : options.chunks_per_thread * options.threads;
    const uint64_t min_symbols = first ? chunks * options.chunk_symbols
        : std::max(options.block_symbols, chunks * options.chunk_symbols);
    block.size = 0;
    block.symbols = 0;

    while (block.size < chunks || block.symbols < min_symbols) {
        if (block.patterns.size() < block.size + 2) block.patterns.resize(block.size + 2);
        pattern_reader::pattern_t& p = block.patterns[block.size];
        if (!reader.next(p)) break;
        block.size++;
        block.symbols += p.seq.size();

        if (options.reverse_complements) {
            pattern_reader::pattern_t& rc = block.patterns[block.size++];
            rc.name = p.name + "_rc";
            rc.seq = reverse_complement(p.seq);
            rc.quality.assign(p.quality.rbegin(), p.quality.rend());
            block.symbols += rc.seq.size();
        }
    }

    block.buffers.resize(block.size * options.outputs.size());
    block.chunk_starts.clear();
    block.chunk_starts.push_back(0);
    const uint64_t symbols_per_chunk = std::max<uint64_t>(1, block.symbols / chunks);
    uint64_t chunk_syms = 0;

    for (uint64_t j = 0; j < block.size; j++) {
        chunk_syms += block.patterns[j].seq.size();

        if (chunk_syms >= symbols_per_chunk) {
            block.chunk_starts.push_back(j + 1);
            chunk_syms = 0;
        }
    }

    if (block.chunk_starts.back() != block.size) block.chunk_starts.push_back(block.size);
}

/**
 * @brief writes the outputs of the patterns of a block in their order and frees them
 * @param options the options of process_patterns
 * @param block the block
 * @param progress the progress of the search phase, which is advanced by the number of patterns of the block
 */
inline void write_pattern_block(const pattern_processing_t& options, pattern_block_t& block, phase_progress& progress)
{
    const uint64_t num_outputs = options.outputs.size();

    for (uint64_t j = 0; j < block.size; j++) {
        for (uint64_t s = 0; s < num_outputs; s++) {
            std::string& buffer = block.buffers[j * num_outputs + s];
            if (options.outputs[s] != nullptr && !buffer.empty()) options.outputs[s]->write(buffer.data(), buffer.size());
            std::string().swap(buffer);
        }
    }

    progress.advance(block.size);
}

/**
 * @brief processes the patterns of a reader with several threads (see the file comment); the outputs are written in
 *        the order of the patterns; with options.log, an empty line and the phase "searching <number> patterns using
 *        <p> threads" are printed, with the progress if the reader knows the number of patterns
 * @tparam process_fnc_t type of process
 * @param reader the reader of the patterns
 * @param options the options
 * @param process process(chunk, thread) processes the patterns of a chunk on thread thread (in [0, options.threads))
 *        and appends their output to chunk.out(j, s)
 * @return the totals
 */
template <typename process_fnc_t>
pattern_totals_t process_patterns(pattern_reader& reader, const pattern_processing_t& options, process_fnc_t process)
{
    pattern_totals_t totals;
    const uint64_t threads = std::max<uint64_t>(1, options.threads);
    const uint64_t num_outputs = options.outputs.size();
    // the number of patterns, with their reverse complements (0 if the reader does not know it)
    const uint64_t num_patterns = reader.num_patterns() * (options.reverse_complements ? 2 : 1);
    auto time = now();
    if (options.log) log_out() << std::endl;
    log_phase_start(options.log, time, "searching " + (num_patterns == 0 ? std::string("the") : format_number(num_patterns))
        + " patterns using " + format_threads(threads));
    phase_progress progress(options.log, num_patterns);
    // block b is held in blocks[b % 3]: while block b is processed, block b - 1 is written and block b + 1 is read
    std::array<pattern_block_t, 3> blocks;
    auto start = now();

    auto read = [&](pattern_block_t& block, bool first) {
        read_pattern_block(reader, options, block, first);
        totals.patterns += block.size;
        totals.symbols += block.symbols;
    };

    read(blocks[0], true);
    uint64_t b = 0;

    for (; blocks[b % 3].size != 0; b++) {
        pattern_block_t& block = blocks[b % 3];
        const uint64_t num_chunks = block.chunk_starts.size() - 1;
        std::atomic<uint64_t> next_chunk = 0;

        run_on_threads(threads, [&](uint64_t thread, uint64_t) {
            if (thread == 0) {
                write_pattern_block(options, blocks[(b + 2) % 3], progress);
                blocks[(b + 2) % 3].size = 0;
                read(blocks[(b + 1) % 3], false);
            }

            for (uint64_t c = next_chunk++; c < num_chunks; c = next_chunk++) {
                const uint64_t first = block.chunk_starts[c];
                process(pattern_chunk(block.patterns.data() + first, block.buffers.data() + first * num_outputs,
                    block.chunk_starts[c + 1] - first, num_outputs), thread);
            }
        });
    }

    write_pattern_block(options, blocks[(b + 2) % 3], progress);
    progress.finish();
    totals.time_wall = time_diff_ns(start, now());
    log_phase_end(options.log, time);
    return totals;
}

// ############################# QUERY SUMMARY #############################

/** @brief the summary of the queries of a patterns file (see print_query_summary) */
struct query_summary_t {
    std::string pattern_details; // appended to the line "patterns: ..." after a comma (empty: nothing)
    std::vector<std::string> results; // the lines after the line "patterns: ..." (see format_count_line)
    uint64_t time_query = 0; // the query time in nanoseconds, summed over the threads
    std::string result_name; // the name of a result whose time is printed, e.g. "occurrence" (empty: none)
    uint64_t num_results = 0; // the number of results
    uint64_t additional_memory = 0; // the memory allocated during the search in bytes (0: not printed)
};

/**
 * @brief formats a line of the results of the queries: "<name>: <count> (<count / patterns> per pattern)"
 * @param name the name of the result, e.g. "occurrences"
 * @param count the number of results
 * @param num_patterns the number of patterns
 * @param details appended to the parentheses after a comma (empty: nothing)
 */
inline std::string format_count_line(const std::string& name, uint64_t count, uint64_t num_patterns,
    const std::string& details = "")
{
    return name + ": " + format_number(count) + " (" + format_ratio(num_patterns == 0 ? 0.0 :
        double(count) / double(num_patterns)) + " per pattern" + (details.empty() ? "" : ", " + details) + ")";
}

/**
 * @brief prints the summary of the queries of a patterns file:
 *        patterns: 200 (average length 20 symbols), <pattern_details>
 *        <results>
 *        query time (summed over 4 threads): 4883 us = 24 us per pattern, 1220 ns per symbol, 379 ns per occurrence
 *        wall time: 2627 us
 *        additional memory: 348 KB
 * @param totals the totals of process_patterns
 * @param threads the number of threads
 * @param summary the results and times of the queries
 */
inline void print_query_summary(const pattern_totals_t& totals, uint64_t threads, const query_summary_t& summary)
{
    const uint64_t num_patterns = totals.patterns;
    std::cout << "patterns: " << format_number(num_patterns);
    if (num_patterns != 0) std::cout << " (average length " << format_ratio(double(totals.symbols) / num_patterns)
        << " symbols)";
    if (!summary.pattern_details.empty()) std::cout << ", " << summary.pattern_details;
    std::cout << '\n';
    for (const std::string& line : summary.results) std::cout << line << '\n';
    std::cout << "query time" << (threads > 1 ? " (summed over " + format_threads(threads) + ")" : std::string())
              << ": " << format_time(summary.time_query);

    if (num_patterns != 0) {
        std::cout << " = " << format_time(summary.time_query / num_patterns) << " per pattern";
        if (totals.symbols != 0) std::cout << ", " << format_time(summary.time_query / totals.symbols) << " per symbol";
        if (!summary.result_name.empty() && summary.num_results != 0)
            std::cout << ", " << format_time(summary.time_query / summary.num_results) << " per " << summary.result_name;
    }

    std::cout << "\nwall time: " << format_time(totals.time_wall) << '\n';
    if (summary.additional_memory != 0) std::cout << "additional memory: " << format_size(summary.additional_memory) << '\n';
    std::cout << std::flush;
}

#endif // MOVE_R_MISC_PARALLEL_PATTERNS_HPP
