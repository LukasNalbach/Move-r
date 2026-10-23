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

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <misc/index_cli.hpp>
#include <move_r/move_r.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "move-r-bench-ops.tpp"

std::string index_path;
direct_ofstream measurement_file;
move_r_bench_ops::config_t config;
// default of max_pattern_length
constexpr uint64_t default_max_pattern_length = uint64_t { 1 } << 22;
// the largest pattern length measured (-M)
uint64_t max_pattern_length = default_max_pattern_length;
// SA(i) and BWT(i) are measured at max(min_access_queries, access_queries_per_sample * -q) positions (measurement
// setup)
constexpr uint64_t min_access_queries = 1 << 16;
constexpr uint64_t access_queries_per_sample = 128;
// time in ns after which the SA(i) queries are stopped (measurement setup)
constexpr uint64_t sa_time_budget_ns = 2'000'000'000;

// the widths of the columns of the printed tables: pattern length, operation, locate; operation, queries, ns/query
constexpr int length_column_width = 10;
constexpr int operation_column_width = 11;
constexpr int locate_column_width = 12;
constexpr int access_name_column_width = 12;
constexpr int access_queries_column_width = 12;
constexpr int access_time_column_width = 14;

[[noreturn]] void help(const std::string& msg = "")
{
    if (!msg.empty()) std::cout << "error: " << msg << "\n\n";
    std::cout
        << "move-r-bench-ops: measures the operations of a move-r index (prepend, one-occ, next-occ,\n"
        << "                  locate, SA(i) and BWT(i)).\n\n"
        << "usage: move-r-bench-ops [options] <index_file>\n"
        << "   -q <integer>   number of patterns per pattern length m (default: "
        << bench_ops::default_num_patterns << ")\n"
        << "   -r <integer>   seed for sampling the patterns\n"
        << "   -M <integer>   largest pattern length m (default: " << default_max_pattern_length << ")\n"
        << "   -m <file>      append machine-readable RESULT lines to <file>\n"
        << "   <index_file>   index built by move-r build\n\n"
        << "For m = 1, 2, 4, ... up to -M, the patterns are sampled from the index (no input file is\n"
        << "needed): a pattern is extended by prepending L[i] for a uniformly random position i in\n"
        << "its suffix array interval.\n\n"
        << "prepend, one-occ and next-occ are run once per pattern, in random order over all\n"
        << "operations, pattern lengths and patterns; each query is timed separately and the time\n"
        << "of the two clock reads is subtracted. They are reported in ns/query per m.\n\n"
        << "locate is reported in ns/occurrence per m (measured on the patterns with at most "
        << bench_ops::max_occurrences_per_call << "\n"
        << "occurrences). SA(i) (only with random access to the suffix array) and BWT(i) are reported in\n"
        << "ns/query for uniformly random positions i.\n";
    std::exit(msg.empty() ? 0 : 1);
}

// prints one row per pattern length m and one column per operation
void print_grid(std::string_view profile, const move_r_bench_ops::grid_t& grid)
{
    std::cout << '\n' << profile
              << " | ns per query, one query per pattern, up to " << grid.samples
              << " patterns per m | clock overhead of " << std::fixed << std::setprecision(1)
              << grid.clock_overhead_ns << " ns per query subtracted\n";

    std::cout << std::setw(length_column_width) << "m";
    for (auto kind : grid.kinds)
        std::cout << std::setw(operation_column_width) << move_r_bench_ops::query_kind_name(kind);
    std::cout << std::setw(locate_column_width) << "locate*" << '\n';

    for (uint64_t length = 0; length < grid.lengths.size(); length++) {
        std::cout << std::setw(length_column_width) << grid.lengths[length];

        for (uint64_t kind = 0; kind < grid.kinds.size(); kind++) {
            // print "-" if no query of this operation was measured for this m
            const auto& cell = grid.at(kind, length);
            if (cell.queries == 0) { std::cout << std::setw(operation_column_width) << "-"; continue; }

            std::cout << std::setw(operation_column_width)
                      << std::setprecision(grid.ns_per_query(kind, length) < bench_ops::one_decimal_below_ns ? 1 : 0)
                      << grid.ns_per_query(kind, length);
        }

        const auto& locate = grid.locate[length];
        if (locate.occurrences == 0) std::cout << std::setw(locate_column_width) << "-";
        else std::cout << std::setw(locate_column_width) << std::setprecision(1)
                       << double(locate.time_ns) / double(locate.occurrences);

        std::cout << '\n';
    }

    std::cout << "  prepend: one backward-search step (a count query for a pattern of length m performs m steps)\n";

    auto has = [&](move_r_bench_ops::query_kind_t kind) {
        return std::find(grid.kinds.begin(), grid.kinds.end(), kind) != grid.kinds.end();
    };

    if (has(move_r_bench_ops::query_kind_t::next_occ)) {
        std::cout << "  one-occ: computes one occurrence of the matched pattern; next-occ: computes the next\n"
                  << "  occurrence (measured for the second occurrence of the pattern)\n";
    }

    std::cout << "  * locate in ns per occurrence, the other columns in ns per query\n";
    std::cout << std::defaultfloat;
}

// prints the results of SA(i) and BWT(i)
void print_access(std::string_view profile, const move_r_bench_ops::access_grid_t& access)
{
    std::cout << "\n" << profile << " | SA(i) and BWT(i) at uniformly random positions i, in ns per query\n"
              << std::setw(access_name_column_width) << "operation" << std::setw(access_queries_column_width)
              << "queries" << std::setw(access_time_column_width) << "ns/query" << '\n';

    auto row = [](const char* name, const move_r_bench_ops::access_cell_t& cell) {
        std::cout << std::setw(access_name_column_width) << name << std::setw(access_queries_column_width)
                  << cell.queries << std::fixed << std::setprecision(1) << std::setw(access_time_column_width)
                  << cell.ns_per_query() << std::defaultfloat << '\n';
    };

    if (access.has_sa) row("SA", access.sa);
    row("BWT", access.bwt);

    if (access.has_sa && access.sa.queries < access.bwt.queries)
        std::cout << "  SA was stopped after " << access.sa.queries << " of the "
                  << access.bwt.queries << " positions (time budget reached)\n";
}

// writes one RESULT line per operation and pattern length to the measurement file
void emit_grid(std::string_view profile, const move_r_bench_ops::grid_t& grid)
{
    if (!measurement_file.is_open()) return;

    for (uint64_t length = 0; length < grid.lengths.size(); length++) {
        for (uint64_t kind = 0; kind < grid.kinds.size(); kind++) {
            const auto& cell = grid.at(kind, length);

            result_line(measurement_file)
                ("algo", "move_r_bench_ops")
                ("profile", profile)
                ("query", move_r_bench_ops::query_kind_name(grid.kinds[kind]))
                ("m", grid.lengths[length])
                ("samples", grid.samples)
                ("queries", cell.queries)
                ("time_ns", cell.time_ns)
                ("time_per_query_ns", grid.ns_per_query(kind, length))
                ("clock_overhead_ns", grid.clock_overhead_ns)
                ("checksum", grid.checksum);
        }

        const auto& locate = grid.locate[length];
        if (locate.calls == 0) continue;

        result_line(measurement_file)
            ("algo", "move_r_bench_ops")
            ("profile", profile)
            ("query", "locate")
            ("m", grid.lengths[length])
            ("samples", grid.samples)
            ("queries", locate.calls)
            ("occurrences", locate.occurrences)
            ("time_ns", locate.time_ns)
            ("time_per_occ_ns", locate.occurrences == 0 ? 0.0
                                                        : locate.time_ns / double(locate.occurrences))
            ("checksum", locate.checksum);
    }
}

// writes one RESULT line for SA(i) and one for BWT(i) to the measurement file
void emit_access(std::string_view profile, const move_r_bench_ops::access_grid_t& access)
{
    if (!measurement_file.is_open()) return;

    auto emit = [&](const char* query, const move_r_bench_ops::access_cell_t& cell) {
        if (cell.queries == 0) return;

        result_line(measurement_file)
            ("algo", "move_r_bench_ops")
            ("profile", profile)
            ("query", query)
            ("queries", cell.queries)
            ("time_ns", cell.time_ns)
            ("time_per_query_ns", cell.ns_per_query())
            ("clock_overhead_ns", 0)
            ("checksum", cell.checksum);
    };

    emit("SA", access.sa);
    emit("BWT", access.bwt);
}

// loads the index and measures its operations
template <move_r_support support>
void benchmark_index(direct_ifstream& file, std::string_view profile)
{
    move_r<support> index;
    load_index_with_progress(index, file, std::string(profile) + " index");
    file.close();
    drop_from_page_cache(index_path, true);

    // length of the input (without the terminator symbol)
    uint64_t indexed_text_length = index.input_size();
    if (indexed_text_length < bench_ops::min_indexed_text_length) help("the indexed text is too short");

    uint64_t max_length =
        std::min<uint64_t>(max_pattern_length, indexed_text_length - bench_ops::pattern_length_margin);
    std::vector<uint64_t> lengths;
    for (uint64_t length = 1; length <= max_length; length *= 2) lengths.push_back(length);

    auto report_progress = [](uint64_t length) {
        std::cout << "\rsampling patterns of length " << length << " ... " << std::flush;
    };

    try {
        std::cout << "sampling patterns from the index ...\n";
        auto grid = move_r_bench_ops::benchmark(index, lengths, config, report_progress);
        std::cout << "\r" << std::string(bench_ops::progress_line_width, ' ') << "\r";
        print_grid(profile, grid);
        emit_grid(profile, grid);

        // SA(i) and BWT(i) are measured after the other operations
        const uint64_t access_queries = std::max<uint64_t>(
            min_access_queries, config.num_samples * access_queries_per_sample);
        auto access = move_r_bench_ops::benchmark_access(index, config, access_queries, sa_time_budget_ns);
        print_access(profile, access);
        emit_access(profile, access);
    } catch (const std::exception& error) {
        help(std::string("could not run the benchmark: ") + error.what());
    }
}

int main(int argc, char** argv)
{
    index_path = bench_ops::parse_options(argc, argv, config, max_pattern_length, measurement_file, help);
    direct_ifstream file(index_path);
    if (!file.good()) help("could not open the index file");

    // the locate support the index was built with
    move_r_support support = move_r<>::peek_support(file);
    file.seekg(0);

    switch (support) {
        case _count:         benchmark_index<_count>(file, "count"); break;
        case _locate_one:    benchmark_index<_locate_one>(file, "locate-one"); break;
        case _locate_move:   benchmark_index<_locate_move>(file, "locate-move"); break;
        case _locate_rlzsa:  benchmark_index<_locate_rlzsa>(file, "locate-rlzsa"); break;
        default: help("<index_file> is not an index built by move-r build");
    }

    return 0;
}
