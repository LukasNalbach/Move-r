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
 * @brief the modes count and locate of move-r
 */

#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <ips2ra.hpp>
#include <misc/alignment_output_cli.hpp>
#include <misc/cli_help.hpp>
#include <misc/fasta.hpp>
#include <misc/files.hpp>
#include <misc/index_cli.hpp>
#include <misc/memory.hpp>
#include <misc/move_r_cli.hpp>
#include <misc/parallel_patterns.hpp>
#include <misc/pattern_reader.hpp>
#include <misc/query_cli.hpp>
#include <move_r/move_r.hpp>

#include "modes.hpp"

namespace {

/** @brief the options and files of move-r count and move-r locate */
struct query_cli_t {
    bool locate = false; // true <=> move-r locate
    std::string query_name; // count or locate
    query_options_t options; // -p, -o, -m
    exact_alignment_cli_t exact; // -sam, -bed, -rc
    std::string path_index_file;
    std::string path_patterns_file; // empty <=> the patterns are read from stdin interactively
    std::string command_line; // the command line, for the program line of the SAM header
    direct_ifstream index_file;
    direct_ofstream output_file;
} cli;

/**
 * @brief prints the help of move-r count or move-r locate (cli.locate)
 * @param out the stream to print to
 */
void print_help(std::ostream& out)
{
    help_printer h(out);

    if (cli.locate) {
        h.intro("move-r locate: locates the occurrences of patterns (with an index built with -s locate-one: one "
            "occurrence of each pattern).", "move-r locate [options] <index_file> [<patterns_file>]");
    } else {
        h.intro("move-r count: counts the occurrences of patterns.",
            "move-r count [options] <index_file> [<patterns_file>]");
    }

    h.group("arguments");
    h.option("<index_file>", "", "an index built with move-r build");
    h.option("<patterns_file>", "", patterns_file_help + "; without it, "
        "the patterns are read from stdin (one per line) and the results are printed");

    if (cli.locate) {
        h.group("input");
        print_exact_alignment_input_option(h);
    }

    h.group("output");
    h.option("-o <file>", "a path", std::string("write one line \"<name>\\t<count>\" per pattern") +
        (cli.locate ? ", followed by the occurrences in ascending order" : ""));

    if (cli.locate) print_exact_alignment_output_options(h, "move-r build in the FASTA mode");

    h.group("execution");
    h.option_threads();
    h.group("measurement");
    h.option_measurement();
}

/**
 * @brief returns the occurrences of a pattern: all of them, or one if the index supports locating only one
 * @tparam support the support type of the index
 * @param index the index
 * @param pattern the pattern
 * @param occurrences is set to the occurrences
 * @return the number of occurrences of the pattern
 */
template <move_r_support support>
uint64_t locate_pattern(const move_r<support>& index, const std::string& pattern, std::vector<uint48_t>& occurrences)
{
    occurrences.clear();

    if constexpr (move_r<support>::supports_multiple_locate) {
        occurrences = index.locate(pattern);
        return occurrences.size();
    } else {
        auto ctx = index.query();
        for (uint64_t j = pattern.size(); j > 0; j--) if (!ctx.prepend(pattern[j - 1])) return 0;
        if (pattern.empty()) return index.input_size();
        occurrences.push_back(ctx.one_occ());
        return ctx.num_occ();
    }
}

/**
 * @brief counts or locates the patterns of the patterns file, writes the results and prints the summary
 * @tparam support the support type of the index
 * @param index the index
 */
template <move_r_support support>
void run_queries(const move_r<support>& index)
{
    const bool sam_output = cli.exact.sam_file.is_open();
    const bool bed_output = cli.exact.bed_file.is_open();
    const fasta_sequence_data& sequences = index.seq_data();

    if ((sam_output || bed_output) && !sequences.has_sequences())
        throw std::runtime_error("-sam and -bed need an index built with move-r build in the FASTA mode");

    pattern_reader patterns(cli.path_patterns_file);
    require_readable(patterns, cli.path_patterns_file);
    if (sam_output) write_sam_header(cli.exact.sam_file, sequences, "move-r", cli.command_line);
    const bool search_reverse_complements = (sam_output || bed_output) && cli.exact.reverse_complements;

    // the counters and buffers of a thread
    struct alignas(cache_line_bytes) thread_state_t {
        uint64_t last_sequence = 0; // index of the sequence of the previous occurrence (sequence_index searches from it)
        uint64_t num_occurrences = 0;
        uint64_t time_query = 0;
        // text positions of the occurrences of the current pattern, sorted before they are written
        std::vector<uint48_t> occurrences;
        // text positions of the occurrences of the reverse complement of the current pattern, sorted before they are
        // written
        std::vector<uint48_t> rc_occurrences;
    };

    std::vector<thread_state_t> states(cli.options.threads);
    enum output_t : uint64_t { occurrences_out, sam_out, bed_out }; // the order of processing.outputs
    pattern_processing_t processing { .threads = cli.options.threads, .outputs = {
        cli.output_file.is_open() ? &cli.output_file : nullptr,
        sam_output ? &cli.exact.sam_file : nullptr,
        bed_output ? &cli.exact.bed_file : nullptr } };

    const pattern_totals_t totals = process_patterns(patterns, processing, [&](const pattern_chunk& chunk, uint64_t thread) {
        thread_state_t& st = states[thread];

        for (uint64_t j = 0; j < chunk.size(); j++) {
            const pattern_reader::pattern_t& pattern = chunk.pattern(j);
            const uint64_t m = pattern.seq.size();
            uint64_t count = 0;
            auto start = now();

            if (cli.locate) {
                if constexpr (move_r<support>::supports_locate) {
                    count = locate_pattern(index, pattern.seq, st.occurrences);

                    if (search_reverse_complements)
                        count += locate_pattern(index, reverse_complement(pattern.seq), st.rc_occurrences);
                }
            } else {
                count = index.count(pattern.seq);
            }

            st.time_query += time_diff_ns(start, now());
            st.num_occurrences += count;

            if (!st.occurrences.empty()) ips2ra::sort(st.occurrences.begin(), st.occurrences.end(), key_of_value {});
            if (!st.rc_occurrences.empty()) ips2ra::sort(st.rc_occurrences.begin(), st.rc_occurrences.end(), key_of_value {});

            if (cli.output_file.is_open()) {
                string_appender out(chunk.out(j, occurrences_out));
                out << pattern.name << '\t' << count;
                for (uint64_t occ : st.occurrences) out << '\t' << occ;
                out << '\n';
            }

            if (sam_output || bed_output) {
                string_appender sam(chunk.out(j, sam_out));
                string_appender bed(chunk.out(j, bed_out));
                write_exact_records(sam_output ? &sam : nullptr, bed_output ? &bed : nullptr, sequences,
                    st.last_sequence, pattern.name, m, count, st.occurrences, st.rc_occurrences);
            }

            st.occurrences.clear();
            st.rc_occurrences.clear();
        }
    });

    const uint64_t num_patterns = totals.patterns;
    const uint64_t num_symbols = totals.symbols;
    uint64_t num_occurrences = 0;
    uint64_t time_query = 0;

    for (const thread_state_t& st : states) {
        num_occurrences += st.num_occurrences;
        time_query += st.time_query;
    }

    if (num_patterns == 0) throw std::runtime_error(cli.path_patterns_file + " contains no patterns");

    print_query_summary(totals, cli.options.threads, {
        .results = { format_count_line("occurrences", num_occurrences, num_patterns) },
        .time_query = time_query,
        .result_name = cli.locate ? "occurrence" : "",
        .num_results = num_occurrences });

    direct_ofstream& measurement_file = cli.options.measurement_file;

    if (measurement_file.is_open()) {
        {
            result_line line(measurement_file);
            line("algo", cli.query_name + "_move_r_" + std::string(move_r_support_name(support)))
                ("text", cli.options.text_name)
                ("n", index.input_size())
                ("pattern_length", num_symbols / num_patterns)
                ("num_patterns", num_patterns)
                ("num_occurrences", num_occurrences)
                ("time_query", time_query)
                ("num_threads", cli.options.threads)
                ("time_wall", totals.time_wall);
            index.log_data_structure_sizes(measurement_file);
        }

        measurement_file.close();
    }
}

/**
 * @brief counts or locates the patterns read from stdin (one per line) and prints the results
 * @tparam support the support type of the index
 * @param index the index
 */
template <move_r_support support>
void run_interactive_queries(const move_r<support>& index)
{
    const std::string intro = "n = " + format_number(index.input_size()) + ", r = "
        + format_number(index.num_bwt_runs()) + ", index size " + format_size(index.size_in_bytes())
        + "; enter one pattern per line ('quit' exits"
        + (cli.locate ? ", 'limit <n>' prints at most n occurrences per pattern)" : ")");
    std::vector<uint48_t> occurrences;

    // returns "<count> occurrence(s)"
    auto occurrence_count = [](uint64_t count) {
        return format_number(count) + (count == 1 ? " occurrence" : " occurrences");
    };

    run_interactive(intro, [&](const std::string& pattern, uint64_t max_printed) {
        auto start = now();

        if (cli.locate) {
            if constexpr (move_r<support>::supports_locate) {
                const uint64_t count = locate_pattern(index, pattern, occurrences);
                ips2ra::sort(occurrences.begin(), occurrences.end(), key_of_value {});
                std::cout << occurrence_count(count) << " (" << format_time(time_diff_ns(start, now())) << ")";
                if (!occurrences.empty()) std::cout << (move_r<support>::supports_multiple_locate ? ":" : ", e.g.:");
                print_occurrences(std::cout, occurrences, max_printed);
                std::cout << std::endl;
            }
        } else {
            const uint64_t count = index.count(pattern);
            std::cout << occurrence_count(count) << " (" << format_time(time_diff_ns(start, now())) << ")" << std::endl;
        }
    });
}

/**
 * @brief runs move-r count or move-r locate
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments
 * @param locate true <=> move-r locate
 * @return the exit code
 */
int run_query(int argc, char** argv, bool locate)
{
    cli.locate = locate;
    cli.query_name = locate ? "locate" : "count";

    return run_cli(argc, argv, print_help, [&] {
        cli_arguments args(argc, argv);

        while (args.at_option()) {
            const std::string s = args.next();
            if (read_query_option(args, s, cli.options, _option_threads | _option_output | _option_measurement)) continue;
            if (!locate || !read_exact_alignment_option(args, s, cli.exact)) throw unknown_option(s);
        }

        if (args.empty() || args.remaining() > 2) throw usage_error("expected <index_file> and at most one <patterns_file>");
        set_default_query_threads(cli.options);
        cli.command_line = command_line_of("move-r", argc, argv);
        cli.path_index_file = args.next();
        if (!args.empty()) cli.path_patterns_file = args.next();
        cli.index_file.open(cli.path_index_file);
        require_readable(cli.index_file, cli.path_index_file);

        if (!cli.path_patterns_file.empty() && cli.path_patterns_file != "-")
            require_readable(std::ifstream(cli.path_patterns_file), cli.path_patterns_file);

        open_output_file(cli.output_file, cli.options.path_output_file);
        open_exact_alignment_files(cli.exact);
        const move_r_support support = move_r<>::peek_support(cli.index_file);

        const bool known = dispatch_move_r_index(support, [&]<move_r_support supp>() {
            if constexpr (!move_r<supp>::supports_locate) {
                if (locate) throw std::runtime_error(cli.path_index_file + " does not support locate (see move-r build -s)");
            }

            move_r<supp> index;
            load_index_with_progress(index, cli.index_file);
            cli.index_file.close();
            drop_from_page_cache(cli.path_index_file, true);
            if (cli.path_patterns_file.empty()) run_interactive_queries(index);
            else run_queries(index);
        });

        if (!known) throw std::runtime_error(cli.path_index_file + " is not an index built with move-r build");
        close_output_file(cli.output_file, cli.options.path_output_file);
        close_exact_alignment_files(cli.exact);
    });
}

} // namespace

int move_r_count(int argc, char** argv) { return run_query(argc, argv, false); }

int move_r_locate(int argc, char** argv) { return run_query(argc, argv, true); }
