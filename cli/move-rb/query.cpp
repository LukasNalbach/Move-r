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
 * @file query.cpp
 * @brief the modes count and locate of move-rb
 */

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <ips2ra.hpp>

#include <misc/apm_query_cli.hpp>
#include <misc/cli_help.hpp>
#include <misc/fasta.hpp>
#include <misc/files.hpp>
#include <misc/index_cli.hpp>
#include <misc/move_rb_cli.hpp>
#include <misc/query_cli.hpp>
#include <move_rb/move_rb.hpp>

#include "modes.hpp"

namespace {

// ############################# OPTIONS #############################

/** @brief the options of move-rb count and move-rb locate */
struct query_cli_t {
    bool locate = true; // move-rb locate (else move-rb count)
    query_options_t options; // -p, -o, -m
    apm_query_cli_t apm; // the options of the approximate queries (see read_apm_query_option)
    std::vector<std::string> options_given; // the options given, in their order (for the warning without patterns)
    std::string path_index_file;
    std::string path_patterns_file; // empty: the queries are read from stdin
    std::ofstream output_file; // -o
} cli;

/**
 * @brief prints the help of move-rb count or move-rb locate (see query_cli_t::locate)
 * @param out the stream to print to
 */
void print_help(std::ostream& out)
{
    help_printer h(out);
    const bool locate = cli.locate;

    if (locate) {
        h.intro("move-rb locate: locates the occurrences of patterns with at most k errors (Hamming or edit distance).",
            "move-rb locate [options] <index_file> [<patterns_file>]");
    } else {
        h.intro("move-rb count: counts the occurrences of patterns with at most k errors (Hamming distance).",
            "move-rb count [options] <index_file> [<patterns_file>]");
    }

    h.group("arguments");
    h.option("<index_file>", "", locate ? "an index built with move-rb build (not with -s count or locate-one)"
        : "an index built with move-rb build");
    h.option("<patterns_file>", "", patterns_file_help
        + (locate ? " (with -sam, the name of a pattern becomes QNAME and its FASTQ qualities become QUAL)" : "")
        + "; without it, the queries are read from stdin (see below)");

    if (locate) {
        h.group("input");
        print_apm_input_options(h);
    }

    h.group("what is computed");
    print_apm_search_options(h, locate);
    print_apm_model_options(h, locate);
    h.group("output");
    h.option("-o <file>", "a path", std::string("write one line \"<name>\\t<count>\" per pattern") +
        (locate ? ", followed by \"\\t<position>:<length>:<errors>\" per occurrence" : ""));
    print_apm_output_options(h, locate, "move-rb build in the FASTA mode");
    h.group("execution");
    h.option_threads();
    print_apm_execution_options(h);
    h.group("measurement");
    h.option_measurement();
    h.group("queries from stdin (one per line)");
    h.option("count <pattern>", "", "the number of exact occurrences");
    h.option("locate <pattern>", "", "the exact occurrences (sorted)");
    h.option("hamming-count <k> <pattern>", "", "the number of occurrences with at most k errors (Hamming distance)");
    h.option("hamming-locate <k> <pattern>", "", "the occurrences with at most k errors (Hamming distance)");
    h.option("edit-locate <k> <pattern>", "", "the occurrences with at most k errors (edit distance, k <= "
        + std::to_string(max_edit_errors) + ")");
    h.option("limit <n>", "", "print at most n occurrences per query (0: all; default: "
        + std::to_string(default_max_printed_occurrences) + ")");
    h.option("help", "", "print the commands");
    h.option("quit, exit", "", "leave the program");
    h.text("Each approximate occurrence is printed with its CIGAR and its alignment (the pattern above the matched "
        "text). For an index built in the FASTA mode, a position is printed as <sequence name>:<1-based "
        "position in the sequence>. The approximate queries are answered with the search scheme "
        + std::string(default_scheme_name) + "; <pattern> is the rest of the line and may contain spaces.");
}

/**
 * @brief reads an option of move-rb count or move-rb locate
 * @param args the arguments, positioned after the option
 * @param option the option
 * @throws usage_error if the option is unknown or its value is invalid
 */
void read_option(cli_arguments& args, const std::string& option)
{
    cli.options_given.push_back(option);
    if (read_query_option(args, option, cli.options, _option_threads | _option_output | _option_measurement)) return;
    if (!read_apm_query_option(args, option, cli.apm)) throw unknown_option(option);
}

// ############################# INTERACTIVE MODE #############################

// the commands of the interactive mode (besides 'limit <n>', 'quit' and 'exit', see run_interactive)
constexpr std::string_view interactive_commands = "count <pattern>, locate <pattern>, hamming-count <k> <pattern>, "
    "hamming-locate <k> <pattern>, edit-locate <k> <pattern>, help";

// returns the description of the commands of the interactive mode
std::string commands_description()
{
    return "commands: " + std::string(interactive_commands) + " ('quit' exits, 'limit <n>' prints at most n "
        "occurrences per pattern); the approximate queries are answered with the search scheme "
        + std::string(default_scheme_name);
}

/**
 * @brief writes the position of an occurrence: <sequence name>:<1-based position in the sequence> for an index built
 *        from a FASTA input, else the position in the text
 * @tparam index_t the index type
 * @param out the stream to write to
 * @param index the index
 * @param pos the position of the occurrence in the text
 * @param last_sequence index of the sequence of the previous position (sequence_index searches from it; the positions
 *        must be passed in increasing order)
 */
template <typename index_t>
void print_position(std::ostream& out, const index_t& index, uint64_t pos, uint64_t& last_sequence)
{
    const fasta_sequence_data& sequences = index.seq_data();

    if (!sequences.has_sequences()) {
        out << pos;
        return;
    }

    const uint64_t seq_idx = sequences.sequence_index(pos, last_sequence);
    out << sequences.sequence_name(seq_idx) << ':' << (pos - sequences.sequence_start(seq_idx) + 1);
}

/**
 * @brief writes the alignment of an occurrence: its CIGAR and the two rows of the alignment (the pattern above the
 *        matched text, with '-' for a gap); the matched text is computed from the pattern and the CIGAR (see
 *        apply_cigar)
 * @param out the stream to write to
 * @param P the pattern
 * @param cigar the alignment of P against the matched text
 */
void print_alignment(std::ostream& out, const std::string& P, const cigar_t<uint8_t>& cigar)
{
    out << "      cigar ";
    append_cigar(out, cigar);
    out << '\n';
    const std::string text = apply_cigar(P, cigar);
    std::string row_p; // the upper row: the pattern
    std::string row_t; // the lower row: the matched text
    uint64_t i = 0; // position in the pattern
    uint64_t j = 0; // position in the matched text

    for (const cigar_run_t<uint8_t>& run : cigar) {
        switch (run.op) {
            case cigar_op_t::MATCH:
            case cigar_op_t::MISMATCH:
                for (uint32_t r = 0; r < run.len; r++) {
                    row_p += i < P.size() ? P[i++] : '?';
                    row_t += j < text.size() ? text[j++] : '?';
                }

                break;

            // a symbol of the pattern that is missing in the text
            case cigar_op_t::INS:
                row_p += i < P.size() ? P[i++] : '?';
                row_t += '-';
                break;

            // a symbol of the text that is missing in the pattern
            case cigar_op_t::DEL:
                row_p += '-';
                row_t += j < text.size() ? text[j++] : '?';
                break;
        }
    }

    out << "      P     " << row_p << '\n';
    out << "      T     " << row_t << std::endl;
}

/**
 * @brief splits the rest of a query line after its command into a number of leading words and the pattern (the rest
 *        of the line, which may contain spaces)
 * @param line the rest of the query line
 * @param num_words the number of leading words
 * @param words is set to the leading words
 * @param pattern is set to the pattern
 * @return whether num_words words and a nonempty pattern were found
 */
bool split_query(const std::string& line, uint64_t num_words, std::vector<std::string>& words, std::string& pattern)
{
    uint64_t i = 0;
    words.clear();

    for (uint64_t w = 0; w < num_words; w++) {
        while (i < line.size() && line[i] == ' ') i++;
        const uint64_t beg = i;
        while (i < line.size() && line[i] != ' ') i++;
        if (beg == i) return false;
        words.push_back(line.substr(beg, i - beg));
    }

    while (i < line.size() && line[i] == ' ') i++;
    pattern = line.substr(i);
    return !pattern.empty();
}

/**
 * @brief returns the number of exact occurrences of P (extending P to the right symbol by symbol)
 * @tparam support the locate support of the index
 * @param index the index
 * @param P the pattern
 */
template <move_r_support support>
uint64_t exact_count(const move_rb<support>& index, const std::string& P)
{
    auto ctx = index.template empty_context<COUNT>();

    for (char c : P) {
        auto [next, found] = ctx.extend(c, RIGHT);
        if (!found) return 0;
        ctx = next;
    }

    return ctx.num_occ();
}

/**
 * @brief prints the number of occurrences and the query time of an answer of the interactive mode, e.g. "3
 *        occurrences (12 us)"
 * @param count the number of occurrences
 * @param start the time the query started
 */
void print_answer_count(uint64_t count, std::chrono::steady_clock::time_point start)
{
    const uint64_t time = time_diff_ns(start, now());
    std::cout << format_number(count) << (count == 1 ? " occurrence" : " occurrences") << " (" << format_time(time)
              << ")";
}

/**
 * @brief answers a query of the interactive mode (see interactive_commands)
 * @tparam support the locate support of the index
 * @param index the index
 * @param line the query line
 * @param max_printed the maximum number of occurrences printed (0: all)
 */
template <move_r_support support>
void answer_query(const move_rb<support>& index, const std::string& line, uint64_t max_printed)
{
    // true <=> the index supports locating all occurrences of a pattern
    constexpr bool can_locate = move_rb<support>::supports_multiple_locate;
    const uint64_t space = line.find(' ');
    const std::string command = line.substr(0, space);
    const std::string rest = space == std::string::npos ? "" : line.substr(space + 1);

    if (command == "help") {
        std::cout << commands_description() << std::endl;
        return;
    }

    // exact queries: <command> <pattern>
    if (command == "count" || command == "locate") {
        if (rest.empty()) {
            std::cout << "error: expected " << command << " <pattern>" << std::endl;
            return;
        }

        auto start = now();

        if (command == "count") {
            print_answer_count(exact_count(index, rest), start);
            std::cout << std::endl;
        } else if constexpr (can_locate) {
            auto ctx = index.template empty_context<LOCATE>();
            bool found = true;

            for (char c : rest) {
                auto [next, extended] = ctx.extend(c, RIGHT);

                if (!extended) {
                    found = false;
                    break;
                }

                ctx = next;
            }

            std::vector<uint48_t> occ;
            if (found) occ = ctx.locate_phase().locate();
            ips2ra::sort(occ.begin(), occ.end(), key_of_value {});
            const uint64_t num_printed = max_printed == 0 ? occ.size() : std::min<uint64_t>(occ.size(), max_printed);
            print_answer_count(occ.size(), start);
            if (!occ.empty()) std::cout << ':';
            uint64_t last_sequence = 0;

            for (uint64_t i = 0; i < num_printed; i++) {
                std::cout << ' ';
                print_position(std::cout, index, occ[i], last_sequence);
            }

            if (num_printed < occ.size()) std::cout << " ... (" << format_number(occ.size() - num_printed) << " more)";
            std::cout << std::endl;
        } else {
            std::cout << "error: this index cannot locate (it was built with -s count or locate-one)" << std::endl;
        }

        return;
    }

    // approximate queries: <command> <k> <pattern>
    if (command == "hamming-count" || command == "hamming-locate" || command == "edit-locate") {
        std::vector<std::string> words;
        std::string pattern;

        if (!split_query(rest, 1, words, pattern)) {
            std::cout << "error: expected " << command << " <k> <pattern>" << std::endl;
            return;
        }

        uint64_t max_errors = 0;

        try {
            max_errors = parse_integer_option("<k>", words[0], 0,
                command == "edit-locate" ? max_edit_errors : max_scheme_errors);
        } catch (const usage_error& error) {
            std::cout << "error: " << error.what() << std::endl;
            return;
        }

        if (command == "edit-locate" && !index.forward_index().byte_alphabet) {
            std::cout << "error: the index " << edit_distance_needs_byte_alphabet << std::endl;
            return;
        }

        const search_scheme_t scheme = min_u_scheme(max_errors);

        if (pattern.size() < scheme.p) {
            std::cout << "error: the pattern has fewer symbols than the search scheme has parts (" << scheme.p << ")"
                      << std::endl;
            return;
        }

        auto start = now();

        if (command == "hamming-count") {
            print_answer_count(index.count_hamming_dist(pattern, scheme), start);
            std::cout << std::endl;
        } else if constexpr (can_locate) {
            // the occurrences with their CIGARs (see print_alignment)
            std::vector<aprx_occ_t<CIGAR>> occ;
            const distance_metric_t metric = command == "hamming-locate" ? HAMMING_DISTANCE : EDIT_DISTANCE;
            const std::vector<cigar_t<uint8_t>> cigars = index.template locate<CIGAR>(pattern, scheme, metric, occ);
            if (metric == EDIT_DISTANCE) filter_edit_distance_occurrences(occ, max_errors);
            ips2ra::sort(occ.begin(), occ.end(), key_of_pos {});
            const uint64_t num_printed = max_printed == 0 ? occ.size() : std::min<uint64_t>(occ.size(), max_printed);
            print_answer_count(occ.size(), start);
            std::cout << (occ.empty() ? "" : ":") << std::endl;
            uint64_t last_sequence = 0;

            for (uint64_t i = 0; i < num_printed; i++) {
                const aprx_occ_t<CIGAR>& o = occ[i];
                std::cout << "   ";
                print_position(std::cout, index, o.pos(), last_sequence);
                std::cout << ": length " << o.len() << ", " << o.err() << (o.err() == 1 ? " error" : " errors") << std::endl;
                print_alignment(std::cout, pattern, cigars[o.cig_idx]);
            }

            if (num_printed < occ.size()) std::cout << "   ... (" << format_number(occ.size() - num_printed) << " more)"
                << std::endl;
        } else {
            std::cout << "error: this index cannot locate (it was built with -s count or locate-one)" << std::endl;
        }

        return;
    }

    std::cout << "error: unknown command " << command << '\n' << commands_description() << std::endl;
}

// ############################# RUNNING THE QUERIES #############################

/**
 * @brief loads the index and runs the queries of the patterns file or of stdin
 * @tparam support the locate support of the index
 * @param index_file the index file, positioned at the beginning of the index
 */
template <move_r_support support>
void run_queries(direct_ifstream& index_file)
{
    using index_t = move_rb<support>;
    index_t index;
    load_index_with_progress(index, index_file, "move-rb");
    index_file.close();
    drop_from_page_cache(cli.path_index_file, true);

    if (!cli.path_patterns_file.empty()) require_apm_alphabet(cli.apm, index, cli.path_index_file);

    if (cli.path_patterns_file.empty()) {
        const auto& fwd = index.forward_index();
        run_interactive("n = " + format_number(fwd.input_size()) + ", r = " + format_number(fwd.num_bwt_runs())
            + ", index size " + format_size(index.size_in_bytes())
            + (index_t::supports_multiple_locate ? "" : " (it cannot locate)") + "; " + commands_description(),
            [&](const std::string& line, uint64_t max_printed) { answer_query(index, line, max_printed); });
        return;
    }

    const apm_query_run_t run { .options = cli.options, .output_file = cli.output_file,
        .path_index_file = cli.path_index_file, .path_patterns_file = cli.path_patterns_file,
        .built_with = "move-rb build in the FASTA mode",
        .measurement_name = "move_rb_" + std::string(move_r_support_name(support)) };
    auto log_sizes = [&](std::ostream& out) { index.log_data_structure_sizes(out); };

    if (!cli.locate) {
        run_apm_count_queries(index, cli.apm, run, log_sizes);
    } else if constexpr (index_t::supports_multiple_locate) {
        if (apm_locate_with_cigar(cli.apm)) run_apm_locate_queries<CIGAR>(index, cli.apm, run, log_sizes);
        else run_apm_locate_queries<NO_CIGAR>(index, cli.apm, run, log_sizes);
    }
}

/**
 * @brief runs move-rb count or move-rb locate
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments; argv[0] is the name of the mode
 * @param locate true <=> move-rb locate
 * @return the exit code
 */
int run_query(int argc, char** argv, bool locate)
{
    cli.locate = locate;
    cli.apm.locate = locate;
    cli.apm.tool = "move-rb";

    return run_cli(argc, argv, print_help, [&] {
        cli.apm.command_line = command_line_of("move-rb", argc, argv);
        cli_arguments args(argc, argv);
        while (args.at_option()) read_option(args, args.next());
        if (args.remaining() < 1 || args.remaining() > 2)
            throw usage_error("expected <index_file> and at most one <patterns_file>");
        cli.path_index_file = args.next();
        if (!args.empty()) cli.path_patterns_file = args.next();

        if (!cli.path_patterns_file.empty()) {
            check_apm_query_options(cli.apm);
        } else if (!cli.options_given.empty()) {
            std::string options;
            for (const std::string& option : cli.options_given) options += (options.empty() ? "" : ", ") + option;
            std::cerr << "warning: without <patterns_file>, " << (cli.options_given.size() == 1 ? "the option "
                + options + " is ignored" : "the options " + options + " are ignored") << std::endl;
        }

        set_default_query_threads(cli.options);
        direct_ifstream index_file(cli.path_index_file);
        require_readable(index_file, cli.path_index_file);
        const move_r_support support = read_move_rb_support(index_file);

        if ((support == _count || support == _locate_one) && locate && !cli.path_patterns_file.empty())
            throw std::runtime_error(cli.path_index_file + " cannot locate (it was built with -s count or locate-one)");

        if (!cli.path_patterns_file.empty()) {
            open_output_file(cli.output_file, cli.options.path_output_file);
            open_apm_query_files(cli.apm);
        }

        const bool known = dispatch_move_rb_index(support,
            [&]<move_r_support index_support>() { run_queries<index_support>(index_file); });

        if (!known) throw std::runtime_error(cli.path_index_file + " is not an index built with move-rb build");
        close_output_file(cli.output_file, cli.options.path_output_file);
        close_apm_query_files(cli.apm);
    });
}

} // namespace

int move_rb_count(int argc, char** argv) { return run_query(argc, argv, false); }

int move_rb_locate(int argc, char** argv) { return run_query(argc, argv, true); }
