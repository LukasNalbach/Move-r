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
 * @file fasta-to-text.cpp
 * @brief the tool fasta-to-text: converts a FASTA file to the text that move-r build indexes in the FASTA mode
 */

#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <system_error>

#include <misc/cli_help.hpp>
#include <misc/fasta.hpp>
#include <misc/log.hpp>

namespace {

/** @brief the options of fasta-to-text */
struct fasta_to_text_cli_t {
    std::string alphabet = "ACGT"; // -A: the symbols that are kept as they are
    std::string path_fasta_file; // <fasta_file>
    std::string path_output_file; // -o (default: <fasta_file>.txt)
};

/**
 * @brief prints the help of fasta-to-text
 * @param out the stream to print to
 */
void print_help(std::ostream& out)
{
    const fasta_to_text_cli_t defaults;
    help_printer h(out);
    h.intro("fasta-to-text: converts a FASTA file to the text that move-r build indexes in the FASTA mode: the line "
        "breaks are removed, the letters are upper-cased, and each header line and each run of symbols outside the "
        "alphabet is replaced with one separator symbol.", "fasta-to-text [options] <fasta_file>");
    h.group("arguments");
    h.option("<fasta_file>", "", "the FASTA file");
    h.group("input");
    h.option("-A <alphabet>", "a string (default: " + defaults.alphabet + ")", "the symbols that are kept as they are");
    h.group("output");
    h.option("-o <file>", "a path (default: <fasta_file>.txt)", "the file the text is written to");
}

} // namespace

int main(int argc, char** argv)
{
    return run_cli(argc, argv, print_help, [&] {
        fasta_to_text_cli_t cli;
        cli_arguments args(argc, argv);

        while (args.at_option()) {
            const std::string s = args.next();

            if (s == "-A") {
                cli.alphabet = args.value("-A");
            } else if (s == "-o") {
                cli.path_output_file = args.value("-o");
            } else {
                throw unknown_option(s);
            }
        }

        if (args.empty()) throw usage_error("expected <fasta_file>");
        cli.path_fasta_file = args.next();
        if (!args.empty()) throw usage_error("expected exactly one <fasta_file>");
        if (cli.path_output_file.empty()) cli.path_output_file = cli.path_fasta_file + ".txt";

        std::error_code error;
        const uint64_t fasta_size = std::filesystem::file_size(cli.path_fasta_file, error);
        if (error) throw std::runtime_error("cannot read " + cli.path_fasta_file);

        auto time = now();
        log_phase_start(true, time, "converting " + cli.path_fasta_file + " (" + format_size(fasta_size) + ")");
        const fasta_sequence_data data = process_fasta(cli.path_fasta_file, cli.path_output_file, cli.alphabet);
        log_phase_end(true, time);

        log_out() << "wrote " << cli.path_output_file << " ("
                  << format_size(std::filesystem::file_size(cli.path_output_file)) << ", "
                  << format_number(data.num_sequences()) << " sequences)" << std::endl;
    });
}
