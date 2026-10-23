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
 * @file extract_modes_cli.hpp
 * @brief the extract modes of move-r and move-rb (extract-sa, extract-bwt, extract-rlbwt): options,
 *        help, loading the index and writing the values to the output file
 */

#ifndef MOVE_R_MISC_EXTRACT_MODES_CLI_HPP
#define MOVE_R_MISC_EXTRACT_MODES_CLI_HPP

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <omp.h>

#include <misc/cli_help.hpp>
#include <misc/files.hpp>
#include <misc/index_cli.hpp>
#include <misc/log.hpp>
#include <misc/value_writer.hpp>
#include <move_r/move_r.hpp>

// the format of a file of n values in [0..n-1] (see write_positions)
inline constexpr std::string_view positions_format
    = "input size + 1 values of 40 bits each (48 bits each for an input of 2^40 or more symbols), little endian";

// the format of a file with the BWT (see write_bwt)
inline constexpr std::string_view bwt_format = "input size + 1 symbols, one byte each for a byte alphabet, else the "
    "fewest bytes each that hold the largest symbol, little endian; the terminator is written as 0";

/** @brief the options of an extract mode */
struct extract_cli_t {
    std::string path_index_file; // <index_file>
    std::string path_output_file; // <output_file>
    bool reversed = false; // -rev: the output is written for the reversed text
    uint64_t num_threads = omp_get_max_threads(); // -p
};

/** @brief what an extract mode writes and what it reads (for its help) */
struct extract_mode_t {
    std::string tool; // move-r or move-rb
    std::string name; // the name of the mode, e.g. extract-sa
    std::string description; // what the mode writes (also the line of the mode in the list of modes)
    std::string output; // the format of the output file and the indexes that provide it
    std::string index_explanation; // the explanation of <index_file>
    std::string reversed_explanation; // the explanation of -rev; empty if the mode has no option -rev
};

/**
 * @brief prints the help of an extract mode
 * @param out the stream to print to
 * @param mode the mode
 */
inline void print_extract_help(std::ostream& out, const extract_mode_t& mode)
{
    help_printer h(out);
    h.intro(mode.tool + " " + mode.name + ": " + mode.description + ".",
        mode.tool + " " + mode.name + " [options] <index_file> <output_file>");
    h.group("arguments");
    h.option("<index_file>", "", mode.index_explanation);
    h.option("<output_file>", "", "the file the values are written to");

    if (!mode.reversed_explanation.empty()) {
        h.group("what is written");
        h.option("-rev", "", mode.reversed_explanation);
    }

    h.group("execution");
    h.option_threads();
    h.text(mode.output);
    h.text("The suffixes are sorted by the internal order of the symbols, which is by decreasing frequency unless the "
        "index was built with -keep-order (see " + mode.tool + " build and " + mode.tool + " info).");
}

/**
 * @brief runs an extract mode: reads the options, opens the index file and calls dispatch(index_file, cli), which reads
 *        the type of the index from the file and calls extract_from_index with it
 * @tparam dispatch_fnc_t type of dispatch
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments; argv[0] is the name of the mode
 * @param mode the mode
 * @param dispatch dispatch(index_file, cli) loads the index and writes the output
 * @return the exit code
 */
template <typename dispatch_fnc_t>
int run_extract_mode(int argc, char** argv, const extract_mode_t& mode, dispatch_fnc_t dispatch)
{
    return run_cli(argc, argv, [&](std::ostream& out) { print_extract_help(out, mode); }, [&] {
        extract_cli_t cli;
        cli_arguments args(argc, argv);

        while (args.at_option()) {
            const std::string s = args.next();

            if (s == "-p") {
                cli.num_threads = args.integer("-p", 1);
            } else if (s == "-rev" && !mode.reversed_explanation.empty()) {
                cli.reversed = true;
            } else {
                throw unknown_option(s);
            }
        }

        if (args.remaining() != 2) throw usage_error("expected <index_file> and <output_file>");
        cli.path_index_file = args.next();
        cli.path_output_file = args.next();
        direct_ifstream index_file(cli.path_index_file);
        require_readable(index_file, cli.path_index_file);
        dispatch(index_file, cli);
    });
}

/**
 * @brief checks whether an index type provides the output of an extract mode, loads the index and writes the output in
 *        the logged phase "writing <output_file>"
 * @tparam index_t type of the index stored in the index file
 * @tparam missing_fnc_t type of missing
 * @tparam write_fnc_t type of write
 * @param index_file the index file, positioned at the beginning of the index
 * @param cli the options
 * @param index_name the name of the index in the loading line (see load_index_with_progress)
 * @param missing missing.template operator()<index_t>(cli) returns the error message if index_t does not provide the
 *        output, else an empty string
 * @param write write(index, cli) writes the output
 * @throws std::runtime_error the message of missing
 */
template <typename index_t, typename missing_fnc_t, typename write_fnc_t>
void extract_from_index(direct_ifstream& index_file, const extract_cli_t& cli, const std::string& index_name,
    missing_fnc_t missing, write_fnc_t write)
{
    const std::string error = missing.template operator()<index_t>(cli);
    if (!error.empty()) throw std::runtime_error(error);
    index_t index;
    load_index_with_progress(index, index_file, index_name);
    index_file.close();
    drop_from_page_cache(cli.path_index_file, true);
    auto time = now();
    log_phase_start(true, time, "writing " + cli.path_output_file);
    write(index, cli);
    log_phase_end(true, time);
}

// the check of extract_from_index for the modes whose output every index provides
inline constexpr auto every_index_provides = []<typename index_t>(const extract_cli_t&) { return std::string(); };

/**
 * @brief calls fnc with the move_r index of the text or, with -rev, with the one of the reversed text
 * @tparam index_t type of the move_rb index
 * @tparam fnc_t type of fnc
 * @param index the move_rb index
 * @param cli the options
 * @param fnc a generic callable taking the move_r index
 */
template <typename index_t, typename fnc_t>
void with_selected_index(const index_t& index, const extract_cli_t& cli, fnc_t fnc)
{
    if (!cli.reversed) fnc(index.forward_index()); else fnc(index.backward_index());
}

/**
 * @brief writes n values to the output file; they are reported by retrieve(report) as report(w, i, value) for every i
 *        in [0..n-1], by at most cli.num_threads threads, with the slot w in [0..cli.num_threads * slots_per_thread -
 *        1]; each slot writes through its own buffer, and the reports of a slot that follow each other are buffered if
 *        their positions increase by one; the buffers of a thread take about one huge page; the
 *        progress over the values is shown in the line of the running phase
 * @tparam bytes_per_value the number of bytes per value
 * @tparam retrieve_fnc_t type of retrieve
 * @param cli the options
 * @param n the number of values
 * @param slots_per_thread the number of slots per thread
 * @param retrieve reports the values
 */
template <uint8_t bytes_per_value, typename retrieve_fnc_t>
void write_range(const extract_cli_t& cli, uint64_t n, uint64_t slots_per_thread, retrieve_fnc_t retrieve)
{
    direct_positioned_file file(cli.path_output_file, n * bytes_per_value);
    const uint64_t num_slots = cli.num_threads * slots_per_thread;
    const uint64_t buffer_values = div_ceil<uint64_t>(huge_page_bytes / bytes_per_value, slots_per_thread);
    // [0..num_slots-1] the writer of each slot
    std::vector<value_writer<bytes_per_value>> writers;
    writers.reserve(num_slots);
    for (uint64_t w = 0; w < num_slots; w++) writers.emplace_back(file, 0, 0, buffer_values);

    phase_progress progress(true, n);

    retrieve([&](uint64_t w, uint64_t i, uint64_t value) {
        writers[w].put(i, value);
        progress.at(i);
    });

    writers.clear();
    if (!file.close_file()) throw std::runtime_error("cannot write " + cli.path_output_file);
}

/**
 * @brief writes n values in [0..n-1] (e.g. suffix array values) with write_range: 40 bits each, or 48 bits each if they
 *        do not fit 40 bits (see wide_values)
 * @tparam retrieve_fnc_t type of retrieve
 * @param cli the options
 * @param n the number of values
 * @param slots_per_thread the number of slots per thread (see write_range)
 * @param retrieve reports the values (see write_range)
 */
template <typename retrieve_fnc_t>
void write_positions(const extract_cli_t& cli, uint64_t n, uint64_t slots_per_thread, retrieve_fnc_t retrieve)
{
    if (wide_values(n - 1)) write_range<uint48_t::width>(cli, n, slots_per_thread, retrieve);
    else write_range<uint40_t::width>(cli, n, slots_per_thread, retrieve);
}

/**
 * @brief writes the suffix array of a move_r index (see positions_format), decoded in pieces (see
 *        move_r::SA_range_in_pieces)
 * @tparam move_r_t type of the move_r index (with move_r_t::supports_multiple_locate)
 * @param index the move_r index
 * @param cli the options
 */
template <typename move_r_t>
void write_suffix_array(const move_r_t& index, const extract_cli_t& cli)
{
    const uint64_t n = index.input_size() + 1;
    const typename move_r_t::retrieve_params params { .l = 0, .r = n - 1, .num_threads = cli.num_threads };
    write_positions(cli, n, params.num_strands, [&](auto report) { index.SA_range_in_pieces(report, params); });
}

/**
 * @brief writes the BWT of a move_r index (see bwt_format)
 * @tparam move_r_t type of the move_r index
 * @param index the move_r index
 * @param cli the options
 */
template <typename move_r_t>
void write_bwt(const move_r_t& index, const extract_cli_t& cli)
{
    const uint64_t n = index.input_size() + 1;

    auto retrieve = [&](auto report) {
        index.BWT_range([&](uint64_t i, uint64_t c) { report(omp_get_thread_num(), i, terminator_as_zero(c)); },
            typename move_r_t::retrieve_params { .l = 0, .r = n - 1, .num_threads = cli.num_threads });
    };

    // the fewest bytes that hold the largest symbol
    const uint8_t symbol_bytes = index.byte_alphabet ? sizeof(uint8_t) : byte_width(index.max_input_symbol());

    move_r_with_symbol_type(symbol_bytes, [&](auto type) {
        write_range<sizeof(typename decltype(type)::type)>(cli, n, 1, retrieve);
    });
}

#endif // MOVE_R_MISC_EXTRACT_MODES_CLI_HPP
