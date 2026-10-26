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
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <ips2ra.hpp>

#include <algorithms/move_rb/apm/model.hpp>
#include <misc/alignment_output_cli.hpp>
#include <misc/cli_help.hpp>
#include <misc/fasta.hpp>
#include <misc/files.hpp>
#include <misc/index_cli.hpp>
#include <misc/memory.hpp>
#include <misc/move_rb_cli.hpp>
#include <misc/parallel_patterns.hpp>
#include <misc/pattern_reader.hpp>
#include <misc/query_cli.hpp>
#include <move_rb/move_rb.hpp>

#include "modes.hpp"

namespace {

// ############################# OPTIONS #############################

// the largest number of errors of a search scheme (-k)
constexpr uint64_t max_scheme_errors = std::numeric_limits<decltype(search_scheme_t::k)>::max();

// the largest number of errors of the edit distance (-k with -d edit)
constexpr uint64_t max_edit_errors = edit_distance_max_errors_128;

// the search scheme without -scheme and of the approximate queries of the interactive mode
constexpr std::string_view default_scheme_name = "min-u";

/** @brief the options of move-rb count and move-rb locate */
struct query_cli_t {
    bool locate = true; // move-rb locate (else move-rb count)
    query_options_t options; // -p, -o, -m
    uint64_t max_errors = 0; // -k
    distance_metric_t metric = HAMMING_DISTANCE; // -d
    std::string scheme_name = std::string(default_scheme_name); // -scheme: the name of a built-in scheme or a file
    bool bestmap = false; // -bestmap: the number of errors is increased from 0 up to -k until an occurrence is found
    uint64_t bestmap_strata = 0; // -x: with -bestmap, the number of further error counts that are reported
    apm_params_t apm_params; // -partition, -no-seed
    apm_model_t apm_model; // -cost, -max, -exact, -class, -iupac, -qual, -tau, -cost-dp
    std::string path_hist_file; // -hist
    std::string path_sam_file; // -sam
    std::string path_paf_file; // -paf
    bool write_seq = false; // -seq: SEQ and QUAL are written in the primary SAM record of a read
    bool reverse_complements = true; // -rc: with -sam or -paf, the reverse complement of each pattern is searched too
    bool best_only = false; // -best: only the primary SAM record of a read is written
    bool no_cigar = false; // -no-cigar: no CIGAR is computed or written
    std::vector<std::string> options_given; // the options given, in their order (for the warning without patterns)
    std::string path_index_file;
    std::string path_patterns_file; // empty: the queries are read from stdin
    std::string command_line; // the command line (for the program line of a SAM file)
    // the search schemes of the queries of the patterns file: [0..e] schemes_by_errors[e'] = the scheme for e'
    // errors, with e = max_errors, if -bestmap or an error model is given; else only the scheme for max_errors errors
    std::vector<search_scheme_t> schemes_by_errors;
    direct_ofstream output_file; // -o
    direct_ofstream hist_file; // -hist
    direct_ofstream sam_file; // -sam
    direct_ofstream paf_file; // -paf
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
    h.option("<patterns_file>", "", std::string("the patterns, in pizza&chili, FASTA or FASTQ format or one per line")
        + (locate ? " (with -sam, the name of a pattern becomes QNAME and its FASTQ qualities become QUAL)" : "")
        + "; without it, the queries are read from stdin (see below)");

    if (locate) {
        h.group("input");
        print_reverse_complements_option(h, true, "with -sam or -paf: also search the reverse complement of each "
            "pattern; its occurrences are written with the SAM flag 0x10 and the PAF strand -");
    }

    h.group("what is computed");
    h.option("-k <errors>", "integer in [0, " + std::to_string(max_scheme_errors) + "] (default: "
        + std::to_string(query_cli_t().max_errors) + ")", "the maximum number of errors");
    if (locate) h.option("-d <distance>", "hamming | edit (default: hamming)", "edit: k <= "
        + std::to_string(max_edit_errors));
    h.option("-scheme <scheme>", "min-u | pigeon-hole | suffix-filter | 01 | <file> (default: "
        + std::string(default_scheme_name) + ")", "the search scheme; a file also sets k");
    h.option("-bestmap", "", "increase the number of errors from 0 up to k until an occurrence is found and report "
        "only the occurrences with the fewest errors");
    h.option("-x <strata>", "integer >= 0 (default: " + std::to_string(query_cli_t().bestmap_strata) + ")",
        "with -bestmap: also report the occurrences with up to this many more errors");

    h.group("error model (by default, every error costs 1; -k bounds the total cost of an alignment)");
    h.option("-cost <rule>", "<ops>=<cost>[@<range>][/<from>,<to>], repeatable", "sets the cost of the operations "
        "<ops> in <range> (for substitutions of a pattern symbol in <from> by a text symbol in <to>; empty: all); "
        "a later rule overrides an earlier one");
    h.option("-max <counter>", "<ops>=<n>[@<range>], n <= " + std::to_string(apm_max_counter_bound) + ", at most "
        + std::to_string(apm_max_counters), "at most n of the operations <ops> in <range>");
    h.option("-exact <range>", "", "no errors in <range> (as -cost sid=inf@<range>)");
    h.option("-class <class>", "<from>,<to>[@<range>], repeatable", "a pattern symbol in <from> also matches the text "
        "symbols in <to>");
    h.option("-iupac", "", "an IUPAC code in a pattern matches each of its bases");
    h.option("-qual", "", "a substitution costs the Phred quality of the pattern symbol (FASTQ patterns), before -cost");
    h.option("-tau <cost>", "auto | integer >= 1 (default: auto)", "a substitution cheaper than this is no error for "
        "the search scheme (its cost still counts); auto: chosen per pattern, such that the part of a search without "
        "errors has the most symbols without free alternatives");

    if (locate) {
        h.option("-cost-dp <mode>", "auto | on | prune | off (default: auto)", "with -d edit: whether the search "
            "bounds the cost by a dynamic program next to its matrix of errors (on), and also prunes out its contexts "
            "by their costs if not every error costs the same and gap openings and counters do not both apply "
            "(prune); auto: prune if it can prune and neither gap openings nor counters apply, else on if the costs "
            "allow at most one insertion or deletion but the search more errors, else off");
    }

    h.text(std::string("<ops>: letters s (substitution), i (insertion: a pattern symbol absent from the text), d "
        "(deletion: a text symbol absent from the pattern), o (opening a gap of insertions or deletions); <cost>: an "
        "integer, or inf (the operation is not allowed). <range>: <a>:<b>, <a>:, :<b> or <i>, positions of the "
        "pattern") + (locate ? " (of the read with -sam or -paf)" : "") + "; a negative value counts from the end "
        "(-5: the last 5 symbols). With an error model, the error count of an occurrence is its cost"
        + (locate ? "; with -d edit, every valid alignment has a reported occurrence at most 2e + 1 positions away "
        "from its start whose cost is not higher, where e is the largest number of errors of a valid alignment." : "."));

    h.group("output");
    h.option("-o <file>", "a path", std::string("write one line \"<name>\\t<count>\" per pattern") +
        (locate ? ", followed by \"\\t<position>:<length>:<errors>\" per occurrence" : ""));

    if (locate) {
        h.option("-sam <file>", "a path", "write the occurrences in SAM format, with CIGAR and the tags NM and MD "
            "(index built with move-rb build in the FASTA mode)");
        h.option("-paf <file>", "a path", "write the occurrences in PAF format, one line per occurrence (index built "
            "with move-rb build in the FASTA mode); can be combined with -sam");
        h.option("-seq", "", "with -sam: write the read (SEQ) and its qualities (QUAL) in the primary record");
        h.option("-best", "", "with -sam: write only the primary record of each read (an occurrence with the fewest "
            "errors)");
        h.option("-no-cigar", "", "with -sam or -paf: no CIGAR is computed or written");
    } else {
        h.option("-hist <file>", "a path", "write one line \"<name>\\t<c_0>\\t<c_1>...\\t<c_k>\" per pattern, c_e = "
            "the number of occurrences with e errors");
    }

    h.group("execution");
    h.option_threads();
    h.option("-partition <partitioning>", "dynamic | uniform (default: dynamic)",
        "how a pattern is partitioned into the parts of the search scheme");
    h.option("-no-seed", "", "match the first part of each search symbol by symbol instead of starting at its search "
        "context computed by the partitioning");
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
        "position in the sequence>. The approximate queries use the search scheme " + std::string(default_scheme_name)
        + "; <pattern> is the rest of the line and may contain spaces.");
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

    // reads the value of an option of the error model with parse
    auto read_model_option = [&](auto parse) {
        const std::string value = args.value(option);

        try {
            parse(value);
        } catch (const std::exception& error) {
            throw usage_error(option + ": " + error.what());
        }
    };

    if (option == "-d" && cli.locate) {
        const std::string metric = args.value(option);
        if (metric == "hamming") cli.metric = HAMMING_DISTANCE;
        else if (metric == "edit") cli.metric = EDIT_DISTANCE;
        else throw usage_error("-d expects hamming or edit");
    } else if (option == "-scheme") {
        cli.scheme_name = args.value(option);
    } else if (option == "-k") {
        cli.max_errors = args.integer(option, 0, max_scheme_errors);
    } else if (option == "-partition") {
        const std::string partitioning = args.value(option);
        if (partitioning == "uniform") cli.apm_params.partitioning = _partition_uniform;
        else if (partitioning == "dynamic") cli.apm_params.partitioning = _partition_dynamic;
        else throw usage_error("-partition expects dynamic or uniform");
    } else if (option == "-no-seed") {
        cli.apm_params.start_at_part_contexts = false;
    } else if (option == "-bestmap") {
        cli.bestmap = true;
    } else if (option == "-x") {
        cli.bestmap_strata = args.integer(option);
    } else if (option == "-hist" && !cli.locate) {
        cli.path_hist_file = args.value(option);
    } else if (option == "-sam" && cli.locate) {
        cli.path_sam_file = args.value(option);
    } else if (option == "-paf" && cli.locate) {
        cli.path_paf_file = args.value(option);
    } else if (option == "-seq" && cli.locate) {
        cli.write_seq = true;
    } else if (option == "-rc" && cli.locate) {
        cli.reverse_complements = read_on_off(args, option);
    } else if (option == "-best" && cli.locate) {
        cli.best_only = true;
    } else if (option == "-no-cigar" && cli.locate) {
        cli.no_cigar = true;
    } else if (option == "-cost") {
        read_model_option([](const std::string& value) { cli.apm_model.costs.emplace_back(parse_apm_cost_rule(value)); });
    } else if (option == "-max") {
        read_model_option([](const std::string& value) { cli.apm_model.counters.emplace_back(parse_apm_counter(value)); });
    } else if (option == "-class") {
        read_model_option([](const std::string& value) { cli.apm_model.classes.emplace_back(parse_apm_class(value)); });
    } else if (option == "-tau") {
        read_model_option([](const std::string& value) {
            cli.apm_model.tau = value == "auto" ? 0 : uint16_t(parse_apm_value(value, "tau")); });
    } else if (option == "-exact") {
        read_model_option([](const std::string& value) {
            cli.apm_model.costs.emplace_back(apm_cost_rule_t{.ops = APM_SUB | APM_INS | APM_DEL,
                .cost = apm_forbidden, .range = parse_apm_range(value)}); });
    } else if (option == "-cost-dp" && cli.locate) {
        const std::string mode = args.value(option);
        if (mode == "auto") cli.apm_model.cost_dp = APM_COST_DP_AUTO;
        else if (mode == "on") cli.apm_model.cost_dp = APM_COST_DP_ON;
        else if (mode == "prune") cli.apm_model.cost_dp = APM_COST_DP_PRUNE;
        else if (mode == "off") cli.apm_model.cost_dp = APM_COST_DP_OFF;
        else throw usage_error("-cost-dp expects auto, on, prune or off");
    } else if (option == "-iupac") {
        cli.apm_model.iupac = true;
    } else if (option == "-qual") {
        cli.apm_model.quality = true;
    } else {
        throw unknown_option(option);
    }
}

// ############################# SEARCH SCHEMES #############################

/**
 * @brief builds the search schemes of the queries of the patterns file (see query_cli_t::schemes_by_errors) from
 *        -scheme and -k; a scheme file sets max_errors
 * @throws usage_error if -scheme names neither a built-in scheme nor a file, or with -bestmap and a scheme file
 */
void build_search_schemes()
{
    const std::string& name = cli.scheme_name;

    if (name != "min-u" && name != "pigeon-hole" && name != "suffix-filter" && name != "01") {
        if (!std::filesystem::exists(name)) throw usage_error("-scheme expects min-u, pigeon-hole, suffix-filter, 01 "
            "or a file");
        std::string file_content;
        const uint64_t file_size = std::filesystem::file_size(name);
        resize_no_init(file_content, file_size);
        direct_ifstream scheme_file(name);
        require_readable(scheme_file, name);
        scheme_file.read(file_content.data(), file_size);
        const search_scheme_t scheme = parse_search_scheme(file_content);
        cli.max_errors = scheme.k;
        if (cli.bestmap) throw usage_error("-bestmap cannot be combined with a search scheme file");
        cli.schemes_by_errors.push_back(scheme);
        return;
    }

    // returns the built-in search scheme for the given number of errors
    auto scheme_of = [&](uint64_t errors) {
        if (name == "pigeon-hole") return pigeon_hole_scheme(errors);
        else if (name == "suffix-filter") return suffix_filter_scheme(errors);
        else if (name == "min-u") return min_u_scheme(errors);
        else return zero_one_scheme(errors);
    };

    if (cli.bestmap || !cli.apm_model.is_default()) {
        for (uint64_t e = 0; e < cli.max_errors; e++) cli.schemes_by_errors.push_back(scheme_of(e));
    }

    cli.schemes_by_errors.push_back(scheme_of(cli.max_errors));
}

// returns the search scheme for at most the given number of errors (see query_cli_t::schemes_by_errors)
const search_scheme_t& scheme_for_errors(uint64_t errors)
{
    return cli.schemes_by_errors[std::min<uint64_t>(errors, cli.schemes_by_errors.size() - 1)];
}

/**
 * @brief returns the search scheme of a query with at most the given number of errors, or of the error model compiled
 *        for its pattern
 * @param errors the maximum number of errors
 * @param costs the error model compiled for the pattern (nullptr: none)
 * @throws std::runtime_error if the scheme allows fewer errors than the error model
 */
const search_scheme_t& scheme_for_query(uint64_t errors, const apm_costs_t* costs)
{
    if (!costs) return scheme_for_errors(errors);
    const search_scheme_t& scheme = scheme_for_errors(costs->max_errors());
    if (scheme.k < costs->max_errors()) throw std::runtime_error("the search scheme allows fewer errors than the "
        "error model");
    return scheme;
}

/**
 * @brief returns the query parameters for a pattern: with an error model, the model compiled for the pattern (it is
 *        valid until the next call on the same thread)
 * @tparam index_t type of the index
 * @param index the index
 * @param pattern the pattern as it is searched
 * @param qual the qualities of the read the pattern was taken from
 * @param reverse true <=> the pattern is the reverse complement of the read
 * @param bound the maximum cost of an alignment
 * @param metric the distance metric
 */
template <typename index_t>
apm_params_t query_params(const index_t& index, const std::string& pattern, const std::string& qual, bool reverse,
    uint64_t bound, distance_metric_t metric)
{
    apm_params_t params = cli.apm_params;
    if (cli.apm_model.is_default()) return params;
    static thread_local apm_costs_t costs;
    costs.compile(cli.apm_model, pattern, qual, reverse, uint16_t(bound), metric, index.forward_index().map_int());
    params.costs = &costs;
    return params;
}

// returns the details of the patterns line of the query summary, e.g. "search scheme min-u, at most 1 error (edit
// distance)"
std::string search_description()
{
    return "search scheme " + cli.scheme_name + ", at most " + format_number(cli.max_errors)
        + (cli.max_errors == 1 ? " error" : " errors")
        + (cli.metric == HAMMING_DISTANCE ? " (Hamming distance)" : " (edit distance)");
}

// returns the memory allocated during the queries in bytes (0 if malloc_count is not linked)
uint64_t additional_memory(uint64_t baseline)
{
    return is_malloc_count_linked ? malloc_count_peak() - baseline : 0;
}

// ############################# COUNT #############################

/**
 * @brief counts the approximate occurrences of the patterns of the patterns file and prints the summary
 * @tparam support the locate support of the index
 * @param index the index
 */
template <move_r_support support>
void run_count_queries(const move_rb<support>& index)
{
    pattern_reader patterns(cli.path_patterns_file);
    require_readable(patterns, cli.path_patterns_file);
    const bool write_hist = cli.hist_file.is_open();
    const uint64_t threads = cli.options.threads;
    const uint64_t baseline_alloc = malloc_count_current();
    malloc_count_reset_peak();

    // the counters and buffers of a thread
    struct alignas(cache_line_bytes) thread_state_t {
        uint64_t num_occurrences = 0;
        uint64_t time_count = 0;
        // [0..k] hist[e] = the number of occurrences of the current pattern with e errors (-hist)
        std::vector<uint64_t> hist;
    };

    std::vector<thread_state_t> states(threads);
    enum output_t : uint64_t { occurrences_out, hist_out }; // the order of processing.outputs
    pattern_processing_t processing { .threads = threads, .outputs = {
        cli.output_file.is_open() ? &cli.output_file : nullptr,
        write_hist ? &cli.hist_file : nullptr } };

    const pattern_totals_t totals = process_patterns(patterns, processing,
        [&](const pattern_chunk& chunk, uint64_t thread) {
            thread_state_t& st = states[thread];

            for (uint64_t j = 0; j < chunk.size(); j++) {
                const pattern_reader::pattern_t& p = chunk.pattern(j);
                const std::string& pattern = p.seq;
                uint64_t count = 0;
                auto start = now();
                apm_params_t params = query_params(index, pattern, p.quality, false, cli.max_errors,
                    HAMMING_DISTANCE);
                const apm_costs_t* costs = params.costs;
                const search_scheme_t& scheme = scheme_for_query(cli.max_errors, costs);

                if (pattern.size() < scheme.p || (costs && pattern.size() == scheme.p)) {
                    throw std::runtime_error("the pattern " + p.name + " has fewer symbols than the search scheme "
                        "has parts");
                }

                // with -hist, the count is the sum of the histogram
                if (write_hist) {
                    st.hist = index.count_hamming_dist_histogram(pattern, scheme, params);
                    for (uint64_t c : st.hist) count += c;
                } else {
                    count = index.count_hamming_dist(pattern, scheme, params);
                }

                st.num_occurrences += count;
                st.time_count += time_diff_ns(start, now());

                if (cli.output_file.is_open()) {
                    string_appender out(chunk.out(j, occurrences_out));
                    out << p.name << '\t' << count << '\n';
                }

                if (write_hist) {
                    string_appender out(chunk.out(j, hist_out));
                    out << p.name;
                    for (uint64_t c : st.hist) out << '\t' << c;
                    out << '\n';
                }
            }
        });

    const uint64_t num_patterns = totals.patterns;
    uint64_t num_occurrences = 0;
    uint64_t time_count = 0;

    for (const thread_state_t& st : states) {
        num_occurrences += st.num_occurrences;
        time_count += st.time_count;
    }

    if (num_patterns == 0) throw std::runtime_error(cli.path_patterns_file + " contains no patterns");
    const uint64_t pattern_length = totals.symbols / num_patterns;
    print_query_summary(totals, threads, query_summary_t {
        .pattern_details = search_description(),
        .results = { format_count_line("occurrences", num_occurrences, num_patterns) },
        .time_query = time_count,
        .additional_memory = additional_memory(baseline_alloc) });

    if (cli.options.measurement_file.is_open()) {
        {
            result_line line(cli.options.measurement_file);
            line("algo", std::string("count_move_rb_") + std::string(move_r_support_name(support)))
                ("text", cli.options.text_name)
                ("n", index.forward_index().input_size())
                ("pattern_length", pattern_length)
                ("num_patterns", num_patterns)
                ("max_mismatches", cli.max_errors)
                ("num_occurrences", num_occurrences)
                ("time_query", time_count)
                ("num_threads", threads)
                ("time_wall", totals.time_wall);
            index.log_data_structure_sizes(cli.options.measurement_file);
        }

        cli.options.measurement_file.close();
    }
}

// ############################# LOCATE #############################

/**
 * @brief locates the occurrences of a read on one strand (with -sam or -paf), sorts them by position and removes the
 *        edit-distance occurrences that overlap a better one (see filter_edit_distance_occurrences); a read with at
 *        most p symbols (p: the number of parts of the search scheme) has no occurrences
 * @tparam cm whether the occurrences carry a CIGAR
 * @tparam index_t type of the index
 * @param index the index
 * @param read the read as it is searched
 * @param qual the qualities of the read
 * @param reverse whether read is the reverse complement of the read
 * @param strand is set to the occurrences and their CIGARs
 */
template <cigar_mode_t cm, typename index_t>
void locate_strand(const index_t& index, const std::string& read, const std::string& qual, bool reverse,
    strand_occurrences_t<cm>& strand)
{
    strand.clear();

    // locates the occurrences with at most the given number of errors (or the cost bound of the error model) into occ
    auto locate_with_bound = [&](uint64_t errors, std::vector<aprx_occ_t<cm>>& occ) {
        apm_params_t params = query_params(index, read, qual, reverse, errors, cli.metric);
        const apm_costs_t* costs = params.costs;
        const search_scheme_t& scheme = scheme_for_query(errors, costs);
        if ((errors > 0 || costs) && read.size() <= scheme.p) return;
        strand.cigars = index.template locate<cm>(read, scheme, cli.metric, occ, params);
        // the edit-distance occurrences of an error model are already sorted and selected (see apm_params_t)
        if (costs != nullptr && cli.metric == EDIT_DISTANCE) return;
        ips2ra::sort(occ.begin(), occ.end(), key_of_pos {});
        filter_edit_distance_occurrences<cm>(occ, costs ? costs->max_errors() : errors);
    };

    if (cli.bestmap) locate_best<cm>(cli.max_errors, cli.bestmap_strata, strand.occurrences, locate_with_bound);
    else locate_with_bound(cli.max_errors, strand.occurrences);
}

/**
 * @brief locates the approximate occurrences of the patterns of the patterns file and prints the summary
 * @tparam support the locate support of the index
 * @tparam cm whether the occurrences carry a CIGAR (with -sam or -paf without -no-cigar)
 * @param index the index
 */
template <move_r_support support, cigar_mode_t cm>
void run_locate_queries(const move_rb<support>& index)
{
    pattern_reader patterns(cli.path_patterns_file);
    require_readable(patterns, cli.path_patterns_file);
    const bool write_sam = cli.sam_file.is_open();
    const bool write_paf = cli.paf_file.is_open();
    const uint64_t threads = cli.options.threads;
    const fasta_sequence_data& sequences = index.seq_data();

    if ((write_sam || write_paf) && !sequences.has_sequences()) {
        throw std::runtime_error(std::string(write_sam ? "-sam" : "-paf") + " needs an index built with move-rb build "
            "in the FASTA mode (" + cli.path_index_file + " stores no sequences)");
    }

    if (write_sam) write_sam_header(cli.sam_file, sequences, "move-rb", cli.command_line);
    const uint64_t baseline_alloc = malloc_count_current();
    malloc_count_reset_peak();

    // the counters and buffers of a thread
    struct alignas(cache_line_bytes) thread_state_t {
        std::vector<aprx_occ_t<>> occurrences; // the occurrences of the current pattern (without -sam and -paf)
        strand_occurrences_t<cm> fwd; // with -sam or -paf: the occurrences of the current read
        strand_occurrences_t<cm> rc; // with -sam or -paf: the occurrences of the reverse complement of the read
        uint64_t last_sequence_sam = 0; // index of the sequence of the previous SAM record (see sequence_index)
        uint64_t last_sequence_paf = 0; // index of the sequence of the previous PAF record (see sequence_index)
        uint64_t num_occurrences = 0;
        uint64_t time_locate = 0;
        uint64_t checksum = 0;
    };

    std::vector<thread_state_t> states(threads);
    enum output_t : uint64_t { occurrences_out, sam_out, paf_out }; // the order of processing.outputs
    pattern_processing_t processing { .threads = threads, .outputs = {
        cli.output_file.is_open() ? &cli.output_file : nullptr,
        write_sam ? &cli.sam_file : nullptr,
        write_paf ? &cli.paf_file : nullptr } };

    const pattern_totals_t totals = process_patterns(patterns, processing,
        [&](const pattern_chunk& chunk, uint64_t thread) {
            thread_state_t& st = states[thread];
            std::vector<aprx_occ_t<>>& occurrences = st.occurrences;

            for (uint64_t j = 0; j < chunk.size(); j++) {
                const pattern_reader::pattern_t& p = chunk.pattern(j);
                const std::string& qname = p.name;
                const std::string& qual = p.quality;
                const std::string& pattern = p.seq;

                // -sam or -paf: the read is aligned on both strands
                if (write_sam || write_paf) {
                    auto start = now();
                    locate_strand<cm>(index, pattern, qual, false, st.fwd);
                    std::string rc_seq;

                    if (cli.reverse_complements) {
                        rc_seq = reverse_complement(pattern);
                        locate_strand<cm>(index, rc_seq, qual, true, st.rc);
                    } else {
                        st.rc.clear();
                    }

                    st.time_locate += time_diff_ns(start, now());
                    const uint64_t num_read_occurrences = st.fwd.occurrences.size() + st.rc.occurrences.size();
                    st.num_occurrences += num_read_occurrences;
                    st.checksum += num_read_occurrences;

                    if (write_sam) {
                        string_appender sam(chunk.out(j, sam_out));
                        write_sam_records<cm>(sam, sequences, st.last_sequence_sam, qname, pattern, rc_seq, qual,
                            st.fwd, st.rc, cli.write_seq, cli.best_only);
                    }

                    if (write_paf) {
                        string_appender paf(chunk.out(j, paf_out));
                        write_paf_records<cm>(paf, sequences, st.last_sequence_paf, qname, pattern.size(), st.fwd,
                            false);
                        write_paf_records<cm>(paf, sequences, st.last_sequence_paf, qname, pattern.size(), st.rc,
                            true);
                    }

                    continue;
                }

                auto start = now();

                // locates the occurrences with at most the given number of errors (or the cost bound of the error
                // model) into occ
                auto locate_with_bound = [&](uint64_t errors, std::vector<aprx_occ_t<>>& occ) {
                    apm_params_t params = query_params(index, pattern, qual, false, errors, cli.metric);
                    const apm_costs_t* costs = params.costs;
                    const search_scheme_t& scheme = scheme_for_query(errors, costs);
                    if (costs && pattern.size() <= scheme.p) return;
                    index.locate(pattern, scheme, cli.metric, occ, params);

                    // with an error model, the edit-distance occurrences are already sorted and selected (see
                    // apm_params_t)
                    if (cli.metric == EDIT_DISTANCE && costs == nullptr) {
                        ips2ra::sort(occ.begin(), occ.end(), key_of_pos {});
                        filter_edit_distance_occurrences<>(occ, errors);
                    }
                };

                if (cli.bestmap) locate_best<>(cli.max_errors, cli.bestmap_strata, occurrences, locate_with_bound);
                else locate_with_bound(cli.max_errors, occurrences);
                st.num_occurrences += occurrences.size();
                st.time_locate += time_diff_ns(start, now());

                if (cli.output_file.is_open()) {
                    if (cli.metric == HAMMING_DISTANCE) ips2ra::sort(occurrences.begin(), occurrences.end(), key_of_pos {});
                    string_appender out(chunk.out(j, occurrences_out));
                    out << qname << '\t' << occurrences.size();
                    for (const auto& occ : occurrences) out << '\t' << occ.pos() << ':' << occ.len() << ':' << occ.err();
                    out << '\n';
                }

                st.checksum += occurrences.size();
                for (const auto& occ : occurrences) st.checksum += occ.pos() + occ.err();
                occurrences.clear();
            }
        });

    const uint64_t num_patterns = totals.patterns;
    uint64_t num_occurrences = 0;
    uint64_t time_locate = 0;
    uint64_t checksum = 0;

    for (const thread_state_t& st : states) {
        num_occurrences += st.num_occurrences;
        time_locate += st.time_locate;
        checksum += st.checksum;
    }

    const uint64_t pattern_length = num_patterns == 0 ? 0 : totals.symbols / num_patterns;
    print_query_summary(totals, threads, query_summary_t {
        .pattern_details = search_description(),
        .results = { format_count_line("occurrences", num_occurrences, num_patterns) },
        .time_query = time_locate,
        .result_name = "occurrence",
        .num_results = num_occurrences,
        .additional_memory = additional_memory(baseline_alloc) });

    if (cli.options.measurement_file.is_open()) {
        {
            result_line line(cli.options.measurement_file);
            line("algo", std::string("locate_move_rb_") + std::string(move_r_support_name(support)))
                ("dist_metric", (cli.metric == HAMMING_DISTANCE ? "hamming" : "edit"))
                ("text", cli.options.text_name)
                ("n", index.forward_index().input_size())
                ("pattern_length", pattern_length)
                ("num_patterns", num_patterns)
                ("max_mismatches", cli.max_errors)
                ("num_occurrences", num_occurrences)
                ("time_query", time_locate)
                ("num_threads", threads)
                ("time_wall", totals.time_wall)
                ("checksum", checksum);
            index.log_data_structure_sizes(cli.options.measurement_file);
        }

        cli.options.measurement_file.close();
    }

    if (write_sam) std::cout << "wrote SAM output to " << cli.path_sam_file << std::endl;
    if (write_paf) std::cout << "wrote PAF output to " << cli.path_paf_file << std::endl;
}

// ############################# INTERACTIVE MODE #############################

// the commands of the interactive mode (besides 'limit <n>', 'quit' and 'exit', see run_interactive)
constexpr std::string_view interactive_commands = "count <pattern>, locate <pattern>, hamming-count <k> <pattern>, "
    "hamming-locate <k> <pattern>, edit-locate <k> <pattern>, help";

// returns the description of the commands of the interactive mode
std::string commands_description()
{
    return "commands: " + std::string(interactive_commands) + " ('quit' exits, 'limit <n>' prints at most n "
        "occurrences per pattern); the approximate queries use the search scheme " + std::string(default_scheme_name);
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
    // true <=> the index can report all occurrences of a pattern
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
    load_index_with_progress(index, index_file, "bidirectional index");
    index_file.close();
    drop_from_page_cache(cli.path_index_file, true);

    if (cli.path_patterns_file.empty()) {
        const auto& fwd = index.forward_index();
        run_interactive("n = " + format_number(fwd.input_size()) + ", r = " + format_number(fwd.num_bwt_runs())
            + ", index size " + format_size(index.size_in_bytes())
            + (index_t::supports_multiple_locate ? "" : " (it cannot locate)") + "; " + commands_description(),
            [&](const std::string& line, uint64_t max_printed) { answer_query(index, line, max_printed); });
        return;
    }

    if (!cli.locate) {
        run_count_queries(index);
    } else if constexpr (index_t::supports_multiple_locate) {
        if ((cli.sam_file.is_open() || cli.paf_file.is_open()) && cli.no_cigar) {
            run_locate_queries<support, NO_CIGAR>(index);
        } else {
            run_locate_queries<support, CIGAR>(index);
        }
    }
}

/**
 * @brief checks the options that only apply to the queries of a patterns file and builds the search schemes
 * @throws usage_error if the options do not fit together
 */
void check_pattern_file_options()
{
    build_search_schemes();

    if (cli.metric == EDIT_DISTANCE && cli.max_errors > max_edit_errors)
        throw usage_error("the edit distance allows at most " + std::to_string(max_edit_errors) + " errors");
    if (cli.best_only && cli.path_sam_file.empty()) throw usage_error("-best needs -sam");
}

/**
 * @brief opens an output file, if its path is not empty
 * @param file the stream to open
 * @param path the path of the file
 * @throws std::runtime_error if the file cannot be created
 */
void open_output_file(direct_ofstream& file, const std::string& path)
{
    if (path.empty()) return;
    file.open(path);
    require_created(file, path);
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

    return run_cli(argc, argv, print_help, [&] {
        cli.command_line = command_line_of("move-rb", argc, argv);
        cli_arguments args(argc, argv);
        while (args.at_option()) read_option(args, args.next());
        if (args.remaining() < 1 || args.remaining() > 2)
            throw usage_error("expected <index_file> and at most one <patterns_file>");
        cli.path_index_file = args.next();
        if (!args.empty()) cli.path_patterns_file = args.next();

        if (!cli.path_patterns_file.empty()) {
            check_pattern_file_options();
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
            open_output_file(cli.hist_file, cli.path_hist_file);
            open_output_file(cli.sam_file, cli.path_sam_file);
            open_output_file(cli.paf_file, cli.path_paf_file);
        }

        const bool known = dispatch_move_rb_index(support,
            [&]<move_r_support index_support>() { run_queries<index_support>(index_file); });

        if (!known) throw std::runtime_error(cli.path_index_file + " is not an index built with move-rb build");
        close_output_file(cli.output_file, cli.options.path_output_file);
        close_output_file(cli.hist_file, cli.path_hist_file);
        close_output_file(cli.sam_file, cli.path_sam_file);
        close_output_file(cli.paf_file, cli.path_paf_file);
    });
}

} // namespace

int move_rb_count(int argc, char** argv) { return run_query(argc, argv, false); }

int move_rb_locate(int argc, char** argv) { return run_query(argc, argv, true); }
