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
 * @file query_cli.hpp
 * @brief the options, the SAM header and the interactive mode shared by the query modes of move-r and move-rb
 */

#ifndef MOVE_R_MISC_QUERY_CLI_HPP
#define MOVE_R_MISC_QUERY_CLI_HPP

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include <data_structures/fasta_sequence_data.hpp>
#include <misc/cli_help.hpp>
#include <misc/files.hpp>
#include <misc/log.hpp>
#include <misc/measurement.hpp>
#include <misc/parallel_patterns.hpp>

// ############################# OPTIONS #############################

/** @brief the options shared by the query modes (each mode accepts a subset, see query_option_t) */
struct query_options_t {
    uint64_t threads = 0; // -p (0: default_query_threads())
    std::string path_output_file; // -o
    std::string text_name; // the <name> of -m
    direct_ofstream measurement_file; // the <file> of -m
    bool reverse_complements = false; // the reverse complement of each pattern is searched, too
};

/** @brief the options of query_options_t, as bits of a mask of the options a mode accepts */
enum query_option_t : uint8_t {
    _option_threads = 1, // -p <threads>
    _option_output = 2, // -o <file>
    _option_measurement = 4 // -m <file> <name>
};

/**
 * @brief reads the option, if it is one of the options of query_options_t that the mode accepts
 * @param args the arguments, positioned after the option
 * @param option the option
 * @param options the options, the value is stored in
 * @param accepted the options the mode accepts (a mask of query_option_t)
 * @return whether option is one of the accepted options
 */
inline bool read_query_option(cli_arguments& args, const std::string& option, query_options_t& options,
    uint8_t accepted)
{
    if (option == "-p" && (accepted & _option_threads)) {
        options.threads = args.integer(option, 1);
    } else if (option == "-o" && (accepted & _option_output)) {
        options.path_output_file = args.value(option);
    } else if (option == "-m" && (accepted & _option_measurement)) {
        read_measurement_option(args, options.measurement_file, options.text_name);
    } else {
        return false;
    }

    return true;
}

// sets the number of threads to default_query_threads() if -p was not given (call after the options are read)
inline void set_default_query_threads(query_options_t& options)
{
    if (options.threads == 0) options.threads = default_query_threads();
}

// ############################# SAM #############################

/**
 * @brief writes the header of a SAM file: the version, one line @SQ per sequence and the program line
 * @param out the SAM file
 * @param sequences the sequences of the text
 * @param tool the name of the tool, e.g. "move-r"
 * @param command_line the command line
 */
inline void write_sam_header(std::ostream& out, const fasta_sequence_data& sequences, const std::string& tool,
    const std::string& command_line)
{
    out << "@HD\tVN:1.6\tSO:unsorted\n";
    for (uint64_t i = 0; i < sequences.num_sequences(); i++)
        out << "@SQ\tSN:" << sequences.sequence_name(i) << "\tLN:" << sequences.sequence_length(i) << '\n';
    out << "@PG\tID:" << tool << "\tPN:" << tool << "\tCL:" << command_line << '\n';
}

/**
 * @brief returns the command line as one string (for the program line of a SAM file)
 * @param tool the name of the tool
 * @param argc the number of command-line arguments of the mode
 * @param argv the command-line arguments of the mode; argv[0] is the name of the mode
 */
inline std::string command_line_of(const std::string& tool, int argc, char** argv)
{
    std::string command_line = tool;
    for (int i = 0; i < argc; i++) command_line += ' ' + std::string(argv[i]);
    return command_line;
}

// ############################# INTERACTIVE MODE #############################

// the number of occurrences printed per pattern in the interactive mode, unless it is changed with "limit <n>"
inline constexpr uint64_t default_max_printed_occurrences = 100;

/**
 * @brief prints at most max_printed occurrences, each after a blank, and " ... (<x> more)" if there are more
 * @tparam occurrences_t type of a random-access range of occurrences
 * @param out the stream to print to
 * @param occurrences the occurrences
 * @param max_printed the maximum number of occurrences printed (0: all)
 */
template <typename occurrences_t>
void print_occurrences(std::ostream& out, const occurrences_t& occurrences, uint64_t max_printed)
{
    const uint64_t num = occurrences.size();
    const uint64_t printed = max_printed == 0 ? num : std::min<uint64_t>(num, max_printed);
    for (uint64_t i = 0; i < printed; i++) out << ' ' << uint64_t(occurrences[i]);
    if (printed < num) out << " ... (" << format_number(num - printed) << " more)";
}

/**
 * @brief reads lines from stdin after the prompt "> " and passes each to answer, until a line is "quit" or "exit" or
 *        the input ends; the line "limit <n>" sets the number of occurrences answer prints per pattern (0: all)
 * @tparam answer_fnc_t type of answer
 * @param intro printed before the first prompt, e.g. "enter one pattern per line ('quit' exits)"
 * @param answer answer(line, max_printed) answers a nonempty line
 */
template <typename answer_fnc_t>
void run_interactive(const std::string& intro, answer_fnc_t answer)
{
    std::cout << intro << std::endl;
    uint64_t max_printed = default_max_printed_occurrences;
    std::string line;

    while (true) {
        std::cout << "> " << std::flush;

        if (!std::getline(std::cin, line)) {
            std::cout << std::endl;
            break;
        }

        if (line == "quit" || line == "exit") break;
        if (line.empty()) continue;

        if (line.rfind("limit ", 0) == 0) {
            try {
                max_printed = parse_integer_option("limit", line.substr(6));
                std::cout << (max_printed == 0 ? std::string("all occurrences are printed")
                    : "at most " + format_number(max_printed) + " occurrences are printed per pattern") << std::endl;
            } catch (const usage_error& error) {
                std::cout << "error: " << error.what() << " (0: all)" << std::endl;
            }

            continue;
        }

        answer(line, max_printed);
    }
}

#endif // MOVE_R_MISC_QUERY_CLI_HPP
