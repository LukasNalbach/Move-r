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

/**
 * @file revert_cli.hpp
 * @brief the mode revert of move-r and move-rb: options, help, loading the index and writing
 *        the text
 */

#ifndef MOVE_R_MISC_REVERT_CLI_HPP
#define MOVE_R_MISC_REVERT_CLI_HPP

#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>

#include <misc/cli_help.hpp>
#include <misc/files.hpp>
#include <misc/index_cli.hpp>
#include <misc/log.hpp>
#include <misc/memory.hpp>

/** @brief the options of a revert mode */
struct revert_cli_t {
    bool in_memory = false; // -im: the whole text is kept in memory while reverting, then written to the output file
    bool reversed = false; // -rev: the reversed text is written
    uint64_t num_threads = omp_get_max_threads(); // -p
    std::string path_index_file; // <index_file>
    std::string path_output_file; // <output_file>
    std::string text_name; // -m: <name>, the name of the text in the measurement file
    direct_ofstream mf; // -m: the measurement file
};

/** @brief what a revert mode reverts (for its help) */
struct revert_mode_t {
    std::string tool; // move-r or move-rb
    std::string description; // what the mode does (also the line of the mode in the list of modes)
    std::string index_explanation; // the explanation of <index_file>
    std::string reversed_explanation; // the explanation of -rev; empty if the mode has no option -rev
};

/**
 * @brief prints the help of a revert mode
 * @param out the stream to print to
 * @param mode the mode
 */
inline void print_revert_help(std::ostream& out, const revert_mode_t& mode)
{
    help_printer h(out);
    h.intro(mode.tool + " revert: " + mode.description + ".",
        mode.tool + " revert [options] <index_file> <output_file>");
    h.group("arguments");
    h.option("<index_file>", "", mode.index_explanation);
    h.option("<output_file>", "", "the file the text is written to");

    if (!mode.reversed_explanation.empty()) {
        h.group("what is written");
        h.option("-rev", "", mode.reversed_explanation);
    }

    h.group("execution");
    h.option_threads();
    h.option("-im", "", "keep the whole text in memory while reverting (n bytes)");
    h.group("measurement");
    h.option_measurement();
}

/**
 * @brief calls fnc(index) with the move_r index that reverts the text: the index itself, or the index of the text
 *        (the reversed text, with -rev) of a move_rb
 * @tparam index_t type of the index
 * @param index the index
 * @param reversed true <=> the reversed text is reverted
 * @param fnc the function
 */
template <typename index_t, typename fnc_t>
inline void with_revert_index(const index_t& index, bool reversed, fnc_t fnc)
{
    if constexpr (requires { index.forward_index(); }) {
        if (!reversed) fnc(index.forward_index()); else fnc(index.backward_index());
    } else {
        fnc(index);
    }
}

/**
 * @brief writes the text (the reversed text, with -rev) of an index to a file or returns it
 * @tparam index_t type of the index
 * @param index the index
 * @param cli the options
 * @param text if not nullptr, the text is stored in it instead of writing it to cli.path_output_file
 */
template <typename index_t>
void revert_text(const index_t& index, const revert_cli_t& cli, no_init_vector<uint8_t>* text)
{
    with_revert_index(index, cli.reversed, [&](const auto& idx) {
        if (text != nullptr) *text = idx.revert(cli.num_threads);
        else idx.revert(cli.path_output_file, cli.num_threads);
    });
}

/**
 * @brief loads an index, reverts it as a logged phase and writes the measurement data
 * @tparam index_t type of the index stored in the index file
 * @param index_file the index file, positioned at the beginning of the index
 * @param output_file the output file (opened)
 * @param cli the options
 * @param index_name the name of the index in the loading line (see load_index_with_progress)
 * @param algo the value of the key algo in the measurement file
 */
template <typename index_t>
void revert_index(direct_ifstream& index_file, direct_ofstream& output_file, revert_cli_t& cli,
    const std::string& index_name, const std::string& algo)
{
    index_t index;
    load_index_with_progress(index, index_file, index_name);
    index_file.close();
    drop_from_page_cache(cli.path_index_file, true);
    const uint64_t length = index.input_size();
    auto time = now();
    const auto time_start = time;
    const std::string with_threads = " with " + format_threads(cli.num_threads);

    if (cli.in_memory) {
        log_phase_start(true, time, "reverting the index in memory" + with_threads);
        no_init_vector<uint8_t> text; // the text: [0..length-1] text[i] = T[i]
        revert_text(index, cli, &text);
        log_phase_end(true, time);
        log_phase_start(true, time, "writing " + cli.path_output_file);

        {
            phase_progress progress(true, length);
            write_in_chunks(output_file, (const char*) text.data(), length, &progress);
        }

        log_phase_end(true, time);
    } else {
        output_file.close();
        log_phase_start(true, time, "reverting the index to " + cli.path_output_file + with_threads);
        revert_text(index, cli, nullptr);
        log_phase_end(true, time);
    }

    const uint64_t time_revert = time_diff_ns(time_start, now());
    if (output_file.is_open()) output_file.close();
    drop_from_page_cache(cli.path_output_file, true);

    if (cli.mf.is_open()) {
        {
            result_line line(cli.mf);
            line("algo", algo)
                ("text", cli.text_name)
                ("n", length)
                ("num_threads", cli.num_threads)
                ("reversed", int(cli.reversed))
                ("time_revert", time_revert);
            index.log_data_structure_sizes(cli.mf);
        }

        cli.mf.close();
    }
}

/**
 * @brief runs a revert mode: reads the options, opens the index file and the output file, and calls
 *        dispatch(index_file, output_file, cli), which reads the type of the index from the file and calls
 *        revert_index with it
 * @tparam dispatch_fnc_t type of dispatch
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments; argv[0] is the name of the mode
 * @param mode the mode
 * @param dispatch dispatch(index_file, output_file, cli) loads and reverts the index
 * @return the exit code
 */
template <typename dispatch_fnc_t>
int run_revert_mode(int argc, char** argv, const revert_mode_t& mode, dispatch_fnc_t dispatch)
{
    return run_cli(argc, argv, [&](std::ostream& out) { print_revert_help(out, mode); }, [&] {
        revert_cli_t cli;
        cli_arguments args(argc, argv);

        while (args.at_option()) {
            const std::string s = args.next();

            if (s == "-m") {
                read_measurement_option(args, cli.mf, cli.text_name);
            } else if (s == "-im") {
                cli.in_memory = true;
            } else if (s == "-p") {
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
        direct_ofstream output_file(cli.path_output_file);
        require_created(output_file, cli.path_output_file);
        dispatch(index_file, output_file, cli);
    });
}

#endif // MOVE_R_MISC_REVERT_CLI_HPP
