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
 * @file bench_ops.hpp
 * @brief parts shared by move-r-bench-ops and move-rb-bench-ops
 */

#ifndef MOVE_R_BENCH_BENCH_OPS_HPP
#define MOVE_R_BENCH_BENCH_OPS_HPP

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include <misc/log.hpp>
#include <misc/measurement.hpp>

namespace bench_ops {

// the default of -seed: the seed for sampling the patterns and for the order of the queries (part of the
// measurement setup: changing it changes the measured queries)
inline constexpr uint64_t default_seed = 0x6a09e667f3bcc909ULL;
// the default of -q of move-r-bench-ops and move-rb-bench-ops: the number of patterns sampled per pattern length
// (measurement setup)
inline constexpr uint64_t default_num_patterns = 512;
// xored into the seed for the order of the queries (move-rb-bench-ops with locate)
inline constexpr uint64_t schedule_seed_salt = 0x3c6ef372fe94f82bULL;
// xored into the seed for the positions of the access queries (move-r-bench-ops)
inline constexpr uint64_t access_seed_salt = 0x5851f42d4c957f2dULL;

// locate is measured only on patterns with at most this many occurrences (measurement setup)
inline constexpr uint64_t max_occurrences_per_call = 1000;

// the printed tables show a time in ns with one decimal below this value, else without decimals
inline constexpr double one_decimal_below_ns = 100;
// the number of characters of the progress line that is overwritten after the sampling
inline constexpr uint64_t progress_line_width = 48;
// the indexed text must have at least min_indexed_text_length symbols; the longest measured pattern
// is pattern_length_margin symbols shorter than the text (measurement setup)
inline constexpr uint64_t min_indexed_text_length = 4;
inline constexpr uint64_t pattern_length_margin = 2;

// the number of pairs of clock reads measure_clock_overhead averages
inline constexpr uint64_t clock_overhead_reps = uint64_t { 1 } << 18;

// the measured pattern lengths up to 96, then the powers of two from first_doubled_pattern_length
// (measurement setup, see pattern_lengths)
inline constexpr uint64_t small_pattern_lengths[] = { 1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64,
    96 };
inline constexpr uint64_t first_doubled_pattern_length = 128;

// the multiplier of the xorshift64* generator (Vigna, "An experimental exploration of Marsaglia's xorshift
// generators, scrambled", the variant with the shifts 12, 25, 27)
inline constexpr uint64_t xorshift64_star_multiplier = 0x2545f4914f6cdd1dULL;

// measurements of locating all occurrences of the patterns of one length
struct locate_cell_t {
    uint64_t calls = 0; // number of measured locate queries
    uint64_t occurrences = 0; // number of reported occurrences
    uint64_t time_ns = 0; // total time in ns
    uint64_t checksum = 0; // hash of the reported occurrences
};

// a measured query: the index of the operation, of the pattern length and of the sampled pattern
struct schedule_entry_t {
    uint32_t kind;
    uint32_t length;
    uint32_t sample;
};

/**
 * @brief returns the next value of a xorshift* random number generator
 * @param state state of the generator (not 0); it is advanced
 */
inline uint64_t next_random(uint64_t& state)
{
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;
    return state * xorshift64_star_multiplier;
}

// returns the average time in ns of the two clock reads around a measured query
inline double measure_clock_overhead()
{
    constexpr uint64_t reps = clock_overhead_reps;
    volatile uint64_t sink = 0;
    uint64_t total = 0;

    for (uint64_t i = 0; i < reps; i++) {
        auto begin = now();
        sink = sink + 1;
        total += time_diff_ns(begin, now());
    }

    return double(total) / double(reps);
}

/**
 * @brief returns the measured pattern lengths 1, 2, 3, 4, 6, 8, ..., 64, 96, 128, 256, 512, ...
 * @param max_length the largest length
 */
inline std::vector<uint64_t> pattern_lengths(uint64_t max_length)
{
    std::vector<uint64_t> lengths;

    for (uint64_t m : small_pattern_lengths)
        if (m <= max_length) lengths.push_back(m);

    for (uint64_t m = first_doubled_pattern_length; m <= max_length; m *= 2)
        lengths.push_back(m);

    return lengths;
}

/**
 * @brief returns the average over the operations of the total time measured at one pattern length; this is the
 *        time budget for measuring locate at that length
 * @tparam grid_t results of a benchmark, with kinds (the operations) and at(kind, length) (cell with time_ns)
 * @param grid the results
 * @param length index of the pattern length
 */
template <typename grid_t>
uint64_t mean_operation_time(const grid_t& grid, uint64_t length)
{
    uint64_t total_time = 0;

    for (uint64_t kind = 0; kind < grid.kinds.size(); kind++)
        total_time += grid.at(kind, length).time_ns;

    return total_time / std::max<uint64_t>(1, grid.kinds.size());
}

/**
 * @brief parses an unsigned integer option value
 * @param value the value
 * @param option name of the option (for the error message)
 * @param help prints an error message and exits
 */
template <typename help_fnc_t>
uint64_t parse_u64(const char* value, const char* option, help_fnc_t help)
{
    try {
        size_t parsed = 0;
        uint64_t result = std::stoull(value, &parsed);
        if (value[parsed] == '\0') return result;
    } catch (...) { }

    help(std::string("invalid integer for ") + option);
    return 0;
}

/**
 * @brief parses the command line of move-r-bench-ops and move-rb-bench-ops: -h, -q <integer>, -r <integer>,
 *        -M <integer>, -m <file> and <index_file>
 * @param argc number of arguments
 * @param argv the arguments
 * @param config num_samples (-q) and seed (-r) are set in it
 * @param max_pattern_length is set to the value of -M
 * @param measurement_file is opened with the file of -m
 * @param help prints an error message (if not empty) and the usage information and exits
 * @return the path of the index file
 */
template <typename config_t, typename help_fnc_t>
std::string parse_options(int argc, char** argv, config_t& config, uint64_t& max_pattern_length,
    direct_ofstream& measurement_file, help_fnc_t help)
{
    std::vector<std::string> positional;

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") help("");

        if (arg == "-q") {
            if (++i >= argc) help("missing value after -q");
            config.num_samples = parse_u64(argv[i], "-q", help);
            if (config.num_samples == 0) help("-q must be positive");
        } else if (arg == "-r") {
            if (++i >= argc) help("missing value after -r");
            config.seed = parse_u64(argv[i], "-r", help);
        } else if (arg == "-M") {
            if (++i >= argc) help("missing value after -M");
            max_pattern_length = parse_u64(argv[i], "-M", help);
            if (max_pattern_length == 0) help("-M must be positive");
        } else if (arg == "-m") {
            if (++i >= argc) help("missing value after -m");
            if (!open_measurement_file(measurement_file, argv[i])) help("could not open the measurement file");
        } else if (!arg.empty() && arg[0] == '-') {
            help("unknown option '" + arg + "'");
        } else {
            positional.emplace_back(std::move(arg));
        }
    }

    if (positional.size() != 1) help("expected <index_file>");
    if (!std::filesystem::is_regular_file(positional[0])) help("the index file does not exist");
    return positional[0];
}

} // namespace bench_ops

#endif // MOVE_R_BENCH_BENCH_OPS_HPP
