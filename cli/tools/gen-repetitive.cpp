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
 * @file gen-repetitive.cpp
 * @brief the tool gen-repetitive: generates a repetitive text (see generate_repetitive_input)
 */

#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <misc/cli_help.hpp>
#include <misc/files.hpp>
#include <misc/log.hpp>
#include <misc/repetitive_input.hpp>

namespace {

// the smallest symbol of the generated texts; the alphabet is [first_symbol..first_symbol+sigma-1]
constexpr uint8_t first_symbol = 'a';

// the smallest alphabet size: a mutation replaces a symbol by another one
constexpr uint64_t min_sigma = 2;

// the largest alphabet size: all symbols fit in a byte
constexpr uint64_t max_sigma = uint64_t { std::numeric_limits<uint8_t>::max() } + 1 - first_symbol;

/** @brief the options of gen-repetitive */
struct gen_repetitive_cli_t {
    // the parameters of the text (-kind, -n, -l, -b, -r, -seed)
    repetitive_params_t params { .kind = _repetitive_versioned, .size = uint64_t { 8 } << 20, .base_length = 65536,
        .block_length = 1024, .mutation_rate = repetitive_params_t().mutation_rate, .seed = 1 };
    uint64_t sigma = 4; // -sigma: the alphabet size
    std::string path_output_file; // <output_file>
};

// returns the names of the kinds of repetitive texts, separated by " | "
std::string kind_names()
{
    std::string names;

    for (uint8_t i = 0; i < _repetitive_num_kinds; i++) {
        if (i != 0) names += " | ";
        names += repetitiveness_kind_name(repetitive_input_kind_t(i));
    }

    return names;
}

/**
 * @brief prints the help of gen-repetitive
 * @param out the stream to print to
 */
void print_help(std::ostream& out)
{
    const gen_repetitive_cli_t defaults;
    const repetitive_params_t& params = defaults.params;
    help_printer h(out);
    h.intro("gen-repetitive: generates a repetitive text.", "gen-repetitive [options] <output_file>");
    h.group("arguments");
    h.option("<output_file>", "", "the file the text is written to");
    h.group("text");
    h.option("-kind <kind>", kind_names() + " (default: " + std::string(repetitiveness_kind_name(params.kind)) + ")",
        "the kind of repetitiveness");
    h.option("-n <length>", "integer >= 1 (default: " + std::to_string(params.size) + ")",
        "the length of the text in bytes");
    h.option("-sigma <sigma>", "integer in [" + std::to_string(min_sigma) + ", " + std::to_string(max_sigma)
        + "] (default: " + std::to_string(defaults.sigma) + ")", "the alphabet size; the symbols are the bytes a, b, c, ...");
    h.option("-l <length>", "integer >= 1 (default: " + std::to_string(params.base_length) + ")",
        "the length of the repeated base string in bytes; n/r is roughly (n/l)/(1-1/sigma)");
    h.option("-b <length>", "integer >= 1 (default: " + std::to_string(params.block_length) + ")",
        "the length of the blocks permuted by blockmove, in bytes");
    h.option("-r <rate>", "number in [0, 1] (default: " + help_value(params.mutation_rate) + ")",
        "the mutation rate per position");
    h.option("-seed <seed>", "integer >= 0 (default: " + std::to_string(params.seed) + ")",
        "the seed of the random number generator (0: a random seed)");
}

/**
 * @brief returns the kind of repetitive text with a name
 * @param name the name (see repetitiveness_kind_name)
 * @throws usage_error if there is no kind with the name
 */
repetitive_input_kind_t parse_kind(const std::string& name)
{
    for (uint8_t i = 0; i < _repetitive_num_kinds; i++)
        if (name == repetitiveness_kind_name(repetitive_input_kind_t(i))) return repetitive_input_kind_t(i);

    throw usage_error("-kind expects " + kind_names());
}

} // namespace

int main(int argc, char** argv)
{
    return run_cli(argc, argv, print_help, [&] {
        gen_repetitive_cli_t cli;
        repetitive_params_t& params = cli.params;
        cli_arguments args(argc, argv);

        while (args.at_option()) {
            const std::string s = args.next();

            if (s == "-kind") {
                params.kind = parse_kind(args.value("-kind"));
            } else if (s == "-n") {
                params.size = args.integer("-n", 1);
            } else if (s == "-sigma") {
                cli.sigma = args.integer("-sigma", min_sigma, max_sigma);
            } else if (s == "-l") {
                params.base_length = args.integer("-l", 1);
            } else if (s == "-b") {
                params.block_length = args.integer("-b", 1);
            } else if (s == "-r") {
                params.mutation_rate = args.number("-r", 0, 1);
            } else if (s == "-seed") {
                params.seed = args.integer("-seed");
            } else {
                throw unknown_option(s);
            }
        }

        if (args.remaining() != 1) throw usage_error("expected <output_file>");
        cli.path_output_file = args.next();
        direct_ofstream out(cli.path_output_file, std::ios::binary);
        require_created(out, cli.path_output_file);

        const std::vector<uint8_t> text = generate_repetitive_input<std::vector<uint8_t>>(
            params, first_symbol, uint8_t(first_symbol + cli.sigma - 1));
        out.write((const char*) text.data(), text.size());
        out.close();

        log_out() << "wrote " << cli.path_output_file << " (" << format_number(text.size()) << " bytes): kind "
                  << repetitiveness_kind_name(params.kind) << ", sigma = " << cli.sigma << ", base length "
                  << format_number(params.base_length) << ", mutation rate " << help_value(params.mutation_rate)
                  << ", seed " << params.seed << std::endl;
    });
}
