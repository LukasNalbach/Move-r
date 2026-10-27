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
 * @file gen-patterns.cpp
 * @brief the tool gen-patterns: writes random substrings of a file as patterns in the Pizza&Chili format
 */

#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

#include <misc/cli_help.hpp>
#include <misc/fasta.hpp>
#include <misc/files.hpp>
#include <misc/log.hpp>
#include <misc/memory.hpp>
#include <misc/temp_files.hpp>

namespace {

/** @brief the options of gen-patterns */
struct gen_patterns_cli_t {
    std::string path_input_file; // <file>
    uint64_t pattern_length = 0; // <length>
    uint64_t num_patterns = 0; // <number>
    std::string path_patterns_file; // <patterns_file>
    std::string forbidden; // <forbidden>: the bytes no pattern contains
    uint64_t seed = 0; // -seed: the seed of the random number generator (0: a random seed)
    bool fasta = false; // -fasta: whether the patterns are drawn from the text of a FASTA file
    std::string fasta_alphabet = std::string(fasta_default_alphabet); // -A
};

/**
 * @brief prints the help of gen-patterns
 * @param out the stream to print to
 */
void print_help(std::ostream& out)
{
    const gen_patterns_cli_t defaults;
    help_printer h(out);
    h.intro("gen-patterns: generates patterns in the Pizza&Chili format: random substrings of a file that contain no "
        "byte of <forbidden>.", "gen-patterns [options] <file> <length> <number> <patterns_file> [<forbidden>]");
    h.group("arguments");
    h.option("<file>", "", "the file the patterns are taken from");
    h.option("<length>", "integer in [1, size of <file>]", "the length of each pattern in bytes");
    h.option("<number>", "integer >= 1", "the number of patterns");
    h.option("<patterns_file>", "", "the file the patterns are written to");
    h.option("<forbidden>", "a string (default: empty)", "the bytes no pattern contains");
    h.group("input");
    h.option("-fasta <mode>", "on | off (default: off)", "on: the patterns are drawn from the text that the build "
        "modes index in their FASTA mode (the header lines removed, each run of symbols outside the alphabet of -A "
        "replaced by one separator symbol) and contain no separator symbol");
    h.option_fasta_alphabet();
    h.group("what is computed");
    h.option("-seed <seed>", "integer in [0, " + std::to_string(std::numeric_limits<uint32_t>::max()) + "] (default: "
        + std::to_string(defaults.seed) + ")", "the seed of the random number generator (0: a random seed)");
    h.text("The patterns file starts with the line \"# number=<number> length=<length> file=<file> "
        "forbidden=<forbidden>\" (<file> without its directory), followed by the patterns without separators.");
}

} // namespace

int main(int argc, char** argv)
{
    return run_cli(argc, argv, print_help, [&] {
        gen_patterns_cli_t cli;
        cli_arguments args(argc, argv);

        while (args.at_option()) {
            const std::string s = args.next();
            if (s == "-seed") {
                cli.seed = args.integer("-seed", 0, std::numeric_limits<uint32_t>::max());
            } else if (s == "-fasta") {
                const std::string mode = args.value(s);
                if (mode != "on" && mode != "off") throw usage_error("-fasta expects on or off");
                cli.fasta = mode == "on";
            } else if (s == "-A") {
                cli.fasta_alphabet = args.value(s);
                if (cli.fasta_alphabet.find(fasta_separator_symbol) != std::string::npos)
                    throw usage_error("-A must not contain the separator symbol");
            } else {
                throw unknown_option(s);
            }
        }

        if (args.remaining() != 4 && args.remaining() != 5)
            throw usage_error("expected <file>, <length>, <number>, <patterns_file> and optionally <forbidden>");

        cli.path_input_file = args.next();
        cli.pattern_length = parse_integer_option("<length>", args.next(), 1);
        cli.num_patterns = parse_integer_option("<number>", args.next(), 1);
        cli.path_patterns_file = args.next();
        if (!args.empty()) cli.forbidden = args.next();

        // in the FASTA mode, the patterns are drawn from the text of the FASTA file, written to a temporary file
        temporary_files temporary;
        std::string path_text = cli.path_input_file;

        if (cli.fasta) {
            path_text = temporary_file_path("gen-patterns-fasta-");
            temporary.paths.push_back(path_text);
            auto time = now();
            log_phase_start(true, time, "storing the FASTA text");
            process_fasta({ cli.path_input_file }, path_text, cli.fasta_alphabet);
            log_phase_end(true, time);
        }

        direct_ifstream input_file(path_text);
        require_readable(input_file, path_text);
        input_file.seekg(0, std::ios::end);
        const uint64_t input_size = input_file.tellg();

        if (cli.pattern_length > input_size) {
            throw std::runtime_error("<length> must be at most the size of " + std::string(cli.fasta ? "the text of " : "")
                + cli.path_input_file + " (" + format_number(input_size) + " bytes)");
        }

        direct_ofstream patterns_file(cli.path_patterns_file);
        require_created(patterns_file, cli.path_patterns_file);

        // is_forbidden[c] = true <=> the byte c occurs in <forbidden> or is the separator symbol of the FASTA mode
        std::array<bool, std::numeric_limits<uint8_t>::max() + 1> is_forbidden {};
        for (char c : cli.forbidden) is_forbidden[uint8_t(c)] = true;
        if (cli.fasta) is_forbidden[uint8_t(fasta_separator_symbol)] = true;

        std::string file_name = cli.path_input_file;
        file_name = file_name.substr(file_name.find_last_of("/\\") + 1);
        patterns_file << "# number=" << cli.num_patterns << " length=" << cli.pattern_length << " file=" << file_name
                      << " forbidden=" << cli.forbidden << "\n";
        input_file.seekg(0, std::ios::beg);

        const uint64_t seed = cli.seed != 0 ? cli.seed : std::random_device()();
        auto time = now();
        log_phase_start(true, time, "generating " + format_number(cli.num_patterns) + " patterns of length "
            + format_number(cli.pattern_length) + " (seed " + std::to_string(seed) + ")");
        std::mt19937 random(seed);
        // the starting positions of the substrings of length <length>
        std::uniform_int_distribution<uint64_t> pos_distribution(0, input_size - cli.pattern_length);
        std::string pattern;
        resize_no_init(pattern, cli.pattern_length);
        phase_progress progress(true, cli.num_patterns);

        for (uint64_t i = 0; i < cli.num_patterns; i++) {
            bool has_forbidden;

            // draw substrings until one contains no forbidden byte
            do {
                input_file.seekg(pos_distribution(random), std::ios::beg);
                read_in_chunks(input_file, pattern.data(), cli.pattern_length);
                has_forbidden = false;

                if (!cli.forbidden.empty() || cli.fasta) {
                    for (char c : pattern) {
                        if (is_forbidden[uint8_t(c)]) {
                            has_forbidden = true;
                            break;
                        }
                    }
                }
            } while (has_forbidden);

            patterns_file.write(pattern.data(), cli.pattern_length);
            progress.at(i);
        }

        progress.finish();
        log_phase_end(true, time);
        input_file.close();
        patterns_file.close();
        if (!cli.fasta) drop_from_page_cache(cli.path_input_file, true);
        drop_from_page_cache(cli.path_patterns_file, true);
    });
}
