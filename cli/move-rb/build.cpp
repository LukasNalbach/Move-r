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
 * @file build.cpp
 * @brief the mode build of move-rb
 */

#include <iostream>
#include <optional>
#include <string>

#include <misc/build_cli.hpp>
#include <misc/cli_help.hpp>
#include <misc/move_r_cli.hpp>
#include <misc/move_rb_cli.hpp>
#include <move_rb/move_rb.hpp>

#include "modes.hpp"

namespace {

// the options of move-rb build
struct move_rb_build_cli_t : build_cli_t {
    move_r_support support = _locate_move; // -s
    // -context-table: the maximum size of the context tables as a share of the size of the index file
    double context_table_share = move_r_params().context_table_share;
    // -context-table-one-length (see move_r_params::context_tables_single_length)
    bool context_tables_single_length = false;
} cli;

/**
 * @brief prints the help of move-rb build
 * @param out the stream to print to
 */
void print_help(std::ostream& out)
{
    help_printer h(out);
    h.intro("move-rb build: builds a bidirectional index of a text, which extends a pattern in both directions.",
        "move-rb build [options] <input_file>...");
    h.group("arguments");
    h.option("<input_file>...", "", "the text; in the FASTA mode one or more FASTA files, which are concatenated");
    h.group("input");
    print_fasta_options(h, "FASTA mode: strip the header lines, replace each run of symbols outside the alphabet of "
        "-A with one separator symbol, which cannot be matched, and store the names and boundaries of the sequences "
        "(for SAM and PAF output)");
    h.option("-keep-order", "", "map the symbols in their order in the input (default: by decreasing frequency); "
        "then the extract modes write the suffix arrays and BWTs in the lexicographic order of the input");
    h.group("what is stored");
    h.option("-s <support>", "count | locate-one | locate-move | locate-rlzsa (default: locate-move)",
        "count: counting\n"
        "locate-one: one occurrence of a pattern\n"
        "locate-move: all occurrences, with the suffix array in move data structures\n"
        "locate-rlzsa: all occurrences, with the suffix array compressed with relative Lempel-Ziv");
    h.option("-context-table <share>", "number >= 0 (default: " + help_value(move_r_params().context_table_share) + ")",
        "the maximum size of the table of the suffix-array intervals of short patterns, as a share of the size of "
        "the index file; 0: no table");
    h.option("-context-table-one-length", "", "store the table for one pattern length only");
    h.group("output");
    h.option("-o <base_name>", "a path (default: <input_file>)", "the index is written to "
        "<base_name>.move-rb-<support>, with <support> = count, locate-one, locate (for -s locate-move) or rlzsa");
    h.group("execution");
    print_execution_options(h, "the BWTs are built with the prefix-free parsing (pfp-rlbwt) and stored in "
        "<input_file>.bwt-runs and <input_file>.bwt-runs-rev, which later builds reuse without checking them "
        "against the input");
    h.group("measurement");
    print_build_measurement_options(h, true);
}

/**
 * @brief reads an option of move-rb build
 * @param args the arguments, positioned after the option
 * @param option the option
 * @throws usage_error if the option is unknown or its value is missing or invalid
 */
void read_option(cli_arguments& args, const std::string& option)
{
    if (read_build_option(args, option, cli, true)) return;

    if (option == "-s") {
        if (!parse_move_r_support(args.value(option), cli.support))
            throw usage_error("-s expects count, locate-one, locate-move or locate-rlzsa");
    } else if (option == "-context-table") {
        cli.context_table_share = read_nonnegative_number(args, option);
    } else if (option == "-context-table-one-length") {
        cli.context_tables_single_length = true;
    } else {
        throw unknown_option(option);
    }
}

/**
 * @brief builds the index and writes it to disk
 * @tparam support the support type of the index
 */
template <move_r_support support>
void build()
{
    fasta_cli_t fasta(cli);
    fasta.scan(true);

    move_rb<support> index = build_with_fallback(cli.mode, cli.construction_mode_given, fasta.input_length(),
        [&](move_r_construction_mode construction_mode) {
        // true <=> the construction reads the FASTA file itself (only with a single file)
        const bool fasta_direct = cli.fasta_mode && cli.path_input_files.size() == 1;
        std::string build_file = fasta.construction_input_file(fasta_direct);

        return move_rb<support>(build_file, {
            .input_file = build_file,
            .mode = construction_mode,
            .num_threads = cli.num_threads,
            .a = cli.balancing_parameter,
            .reorder_alphabet = cli.reorder_alphabet,
            .log = true,
            .measurement_file = cli.measurement(),
            .mf_mds = cli.measurement_mds(),
            .text_name = cli.text_name,
            .fasta = fasta.fasta_input(fasta_direct),
            .context_table_share = cli.context_table_share,
            .context_tables_single_length = cli.context_tables_single_length
        });
    });

    fasta.finish_scan();

    if (cli.fasta_mode) {
        // stores the names and the boundaries of the FASTA sequences (for SAM, BED and PAF output)
        index.set_fasta_sequence_data(std::move(fasta.data));
        fasta.remove_preprocessed_copy();
    }

    write_index(index, cli);
}

} // namespace

int move_rb_build(int argc, char** argv)
{
    return run_cli(argc, argv, print_help, [&] {
        cli_arguments args(argc, argv);
        while (args.at_option()) read_option(args, args.next());
        read_input_files(args, cli);

        cli.path_index_file = move_rb_index_file_name(cli.index_base_name, cli.support);
        begin_build(cli, "the move-rb index", "s = " + move_r_support_option_name(cli.support) + ", balance = "
            + std::to_string(cli.balancing_parameter));

        // the RESULT line of the measurement file; the construction appends its measurements to it
        std::optional<result_line> measurement;

        if (cli.measurement_file.is_open()) {
            measurement.emplace(cli.measurement_file);
            (*measurement)
                ("algo", std::string("build_move_rb_") + std::string(move_r_support_name(cli.support)))
                ("text", cli.text_name)
                ("num_threads", cli.num_threads)
                ("a", cli.balancing_parameter);
        }

        dispatch_move_rb_index(cli.support, []<move_r_support support>() { build<support>(); });
        // ends the RESULT line
        measurement.reset();
    });
}
