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
#include <fstream>
#include <iomanip>
#include <iostream>
#include <misc/index_cli.hpp>
#include <misc/move_rb_cli.hpp>
#include <move_rb/move_rb.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "move-rb-bench-ops.tpp"

std::string index_path;
direct_ofstream measurement_file;
move_rb_bench_ops::config_t config;
// default of max_pattern_length
constexpr uint64_t default_max_pattern_length = uint64_t { 1 } << 18;
// the largest pattern length measured (-M)
uint64_t max_pattern_length = default_max_pattern_length;
// the widths of the columns of the printed tables: pattern length, operation, locate
constexpr int length_column_width = 8;
constexpr int operation_column_width = 11;
constexpr int locate_column_width = 12;

[[noreturn]] void help(const std::string& msg = "")
{
    if (!msg.empty()) std::cout << "error: " << msg << "\n\n";
    std::cout
        << "move-rb-bench-ops: measures the operations of a move-rb index (extension and locate).\n\n"
        << "usage: move-rb-bench-ops [options] <index_file>\n"
        << "   -q <integer>   number of patterns sampled per pattern length m; each operation is\n"
        << "                  measured once on each pattern (default: " << bench_ops::default_num_patterns << ")\n"
        << "   -r <integer>   seed for sampling the patterns\n"
        << "   -M <integer>   largest pattern length m (default: " << default_max_pattern_length << ")\n"
        << "   -m <file>      appends RESULT lines to the file\n"
        << "   <index_file>   index built by move-rb build\n\n"
        << "The operations are measured for the pattern lengths m = 1, 2, 3, 4, 6, 8, ..., 64, 96,\n"
        << "128, 256, ... up to -M. The patterns are sampled from the index, so no input file is needed:\n"
        << "a pattern is extended symbol by symbol, where each symbol is chosen with a probability\n"
        << "proportional to the number of occurrences of the extended pattern.\n\n"
        << "All queries are executed in random order and measured individually. The median time in\n"
        << "ns per query is printed for each operation and m, once with COUNT contexts and once with\n"
        << "LOCATE contexts; for locate, the time in ns per occurrence is printed.\n";
    std::exit(msg.empty() ? 0 : 1);
}

// returns the column header of an operation
std::string_view query_kind_column(move_rb_bench_ops::query_kind_t kind)
{
    using move_rb_bench_ops::query_kind_t;

    switch (kind) {
        case query_kind_t::extend_left:           return "ext-L";
        case query_kind_t::extend_left_switch:    return "ext-L/sw";
        case query_kind_t::extend_right:          return "ext-R";
        case query_kind_t::extend_right_switch:   return "ext-R/sw";
        case query_kind_t::extend_left_all:       return "ext-L-all";
        case query_kind_t::extend_right_all:      return "ext-R-all";
        case query_kind_t::extend_left_chain:     return "ext-L-ch";
        case query_kind_t::extend_right_chain:    return "ext-R-ch";
    }

    return "?";
}

// prints the measurements as a table with one row per pattern length and one column per operation
void print_grid(const move_rb_bench_ops::grid_t& grid)
{
    std::cout << "\n" << (grid.locate_contexts ? "LOCATE contexts" : "COUNT contexts")
              << " | median time in ns per query, " << grid.samples << " patterns per m"
              << " (mean, 90th percentile and maximum: see the -m file) | clock overhead of "
              << std::fixed << std::setprecision(1)
              << grid.clock_overhead_ns << " ns per query subtracted\n";

    std::cout << std::setw(length_column_width) << "m";
    for (auto kind : grid.kinds) std::cout << std::setw(operation_column_width) << query_kind_column(kind);
    if (grid.locate_measured) std::cout << std::setw(locate_column_width) << "locate*";
    std::cout << '\n';

    for (uint64_t length = 0; length < grid.lengths.size(); length++) {
        std::cout << std::setw(length_column_width) << grid.lengths[length];

        for (uint64_t kind = 0; kind < grid.kinds.size(); kind++) {
            std::cout << std::setw(operation_column_width)
                      << std::setprecision(grid.median_ns(kind, length) < bench_ops::one_decimal_below_ns ? 1 : 0)
                      << grid.median_ns(kind, length);
        }

        if (grid.locate_measured) {
            const auto& locate = grid.locate[length];
            std::cout << std::setw(locate_column_width) << std::setprecision(1)
                      << (locate.occurrences == 0 ? 0.0 : double(locate.time_ns) / double(locate.occurrences));
        }

        std::cout << '\n';
    }

    // legend
    auto has = [&](move_rb_bench_ops::query_kind_t kind) {
        return std::find(grid.kinds.begin(), grid.kinds.end(), kind) != grid.kinds.end();
    };

    if (has(move_rb_bench_ops::query_kind_t::extend_left_all)) {
        std::cout << "  ext-L: extension of the pattern by one symbol to the left; ext-L-all: enumeration of all\n"
                  << "  extensions to the left, in ns per extension.\n";
    }

    if (has(move_rb_bench_ops::query_kind_t::extend_right_all)) {
        std::cout << "  ext-R: extension of the pattern by one symbol to the right; ext-R-all: enumeration of all\n"
                  << "  extensions to the right, in ns per extension.\n";
    }

    if (has(move_rb_bench_ops::query_kind_t::extend_left_switch)) {
        std::cout << "  /sw: the last symbol of the pattern has been added in the other direction, i.e., the\n"
                  << "  operation first updates the input intervals of its direction.\n";
    }

    std::cout << "  -ch: " << move_rb_bench_ops::chain_length << " consecutive extensions in the same direction,"
              << " in ns per extension.\n";

    if (grid.locate_measured) {
        std::cout << "  * locate in ns per occurrence; all other columns in ns per query.\n";
    }

    std::cout << std::defaultfloat;
}

// writes the measurements to the measurement file (one RESULT line per operation and pattern length)
void emit_grid(std::string_view support, const move_rb_bench_ops::grid_t& grid)
{
    if (!measurement_file.is_open()) return;
    // the file name of the index (the part of its path after the last '/' or '\\')
    const std::string index = index_path.substr(index_path.find_last_of("/\\") + 1);

    for (uint64_t length = 0; length < grid.lengths.size(); length++) {
        for (uint64_t kind = 0; kind < grid.kinds.size(); kind++) {
            const auto& cell = grid.at(kind, length);

            result_line(measurement_file)
                ("algo", "move_rb_bench_ops")
                ("index", index)
                ("support", support)
                ("query", move_rb_bench_ops::query_kind_name(grid.kinds[kind]))
                ("context", grid.locate_contexts ? "locate" : "count")
                ("m", grid.lengths[length])
                ("samples", grid.samples)
                ("queries", cell.queries)
                ("extensions", cell.extensions)
                ("time_ns", cell.time_ns)
                ("time_per_query_ns", grid.ns_per_query(kind, length))
                ("median_ns", grid.median_ns(kind, length))
                ("p90_ns", grid.p90_ns(kind, length))
                ("max_ns", grid.max_ns(kind, length))
                ("clock_overhead_ns", grid.clock_overhead_ns)
                ("checksum", grid.checksum);
        }

        if (!grid.locate_measured) continue;
        const auto& locate = grid.locate[length];

        result_line(measurement_file)
            ("algo", "move_rb_bench_ops")
            ("index", index)
            ("support", support)
            ("query", "locate")
            ("context", "locate")
            ("m", grid.lengths[length])
            ("samples", grid.samples)
            ("queries", locate.calls)
            ("occurrences", locate.occurrences)
            ("time_ns", locate.time_ns)
            ("time_per_occ_ns", locate.occurrences == 0 ? 0.0 :
                locate.time_ns / double(locate.occurrences))
            ("checksum", locate.checksum);
    }
}

// loads the index and measures its operations on COUNT contexts and (if it supports locate) on LOCATE contexts
template <typename index_t>
void benchmark_index(direct_ifstream& file, std::string_view support)
{
    index_t index;
    load_index_with_progress(index, file);
    file.close();
    drop_from_page_cache(index_path, true);
    std::cout << "stored contexts up to length " << index.context_tables_max_length() << '\n';

    // length of the input without the terminator
    uint64_t indexed_text_length = index.forward_index().input_size();
    if (indexed_text_length < bench_ops::min_indexed_text_length)
        help("the indexed text is too short");

    uint64_t max_length =
        std::min<uint64_t>(max_pattern_length, indexed_text_length - bench_ops::pattern_length_margin);
    std::vector<uint64_t> lengths = bench_ops::pattern_lengths(max_length);

    auto report_progress = [](uint64_t length) {
        std::cout << "\rsampling contexts of length " << length << " ... " << std::flush;
    };

    try {
        std::cout << "sampling COUNT contexts from the index ...\n";
        auto count_grid = move_rb_bench_ops::benchmark<COUNT>(
            index, lengths, config, report_progress);
        print_grid(count_grid);
        emit_grid(support, count_grid);

        if constexpr (index_t::supports_locate) {
            std::cout << "\nsampling LOCATE contexts from the index ...\n";
            auto locate_grid = move_rb_bench_ops::benchmark<LOCATE>(
                index, lengths, config, report_progress);
            print_grid(locate_grid);
            emit_grid(support, locate_grid);
        } else {
            std::cout << "count-only index: skipping LOCATE-context and locate workloads\n";
        }
    } catch (const std::exception& error) {
        help(std::string("could not prepare or run operation workload: ") +
             error.what());
    }
}

int main(int argc, char** argv)
{
    index_path = bench_ops::parse_options(argc, argv, config, max_pattern_length, measurement_file, help);
    direct_ifstream file(index_path);
    if (!file.good()) help("could not open the index file");

    if (!dispatch_move_rb_index(read_move_rb_support(file), [&]<move_r_support support>() {
        benchmark_index<move_rb<support>>(file, move_r_support_name(support));
    })) help("the file is not an index built with move-rb build");
}
