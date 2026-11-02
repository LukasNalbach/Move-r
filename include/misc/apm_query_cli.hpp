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
 * @file apm_query_cli.hpp
 * @brief the approximate queries of the modes count and locate of move-rb: their options, their help
 *        and the queries of a patterns file
 */

#ifndef MOVE_R_MISC_APM_QUERY_CLI_HPP
#define MOVE_R_MISC_APM_QUERY_CLI_HPP

#include <algorithm>
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
#include <misc/parallel_patterns.hpp>
#include <misc/pattern_reader.hpp>
#include <misc/query_cli.hpp>

// ############################# OPTIONS #############################

// the largest number of errors of a search scheme (-k)
inline constexpr uint64_t max_scheme_errors = std::numeric_limits<decltype(search_scheme_t::k)>::max();

// the largest number of errors of the edit distance (-k with -d edit)
inline constexpr uint64_t max_edit_errors = edit_distance_max_errors_128;

// the search scheme without -scheme and of the approximate queries of the interactive mode of move-rb
inline constexpr std::string_view default_scheme_name = "min-u";

// the error of the edit distance and of the error models for an index without a byte alphabet (more distinct bytes
// than internal symbols of one byte, see move_r::byte_alphabet)
inline const std::string edit_distance_needs_byte_alphabet = "has more distinct bytes than internal symbols of one "
    "byte; only the Hamming distance without an error model is supported";

/** @brief the options of the approximate queries of count and locate */
struct apm_query_cli_t {
    bool locate = true; // the mode locate (else count)
    bool given = false; // true <=> one of the options was given
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
    std::string path_bed_file; // -bed
    bool write_seq = false; // -seq: SEQ and QUAL are written in the primary SAM record of a read
    bool reverse_complements = true; // -rc: with -sam, -paf or -bed, the reverse complement of each pattern is searched too
    bool best_only = false; // -best: only the primary SAM record of a read is written
    bool no_cigar = false; // -no-cigar: no CIGAR is computed or written
    std::string tool; // the name of the tool, e.g. move-rb (for the program line of a SAM file)
    std::string command_line; // the command line (for the program line of a SAM file)
    // the search schemes of the queries of the patterns file: [0..e] schemes_by_errors[e'] = the scheme for e'
    // errors, with e = max_errors, if -bestmap or an error model is given; else only the scheme for max_errors errors
    std::vector<search_scheme_t> schemes_by_errors;
    direct_ofstream hist_file; // -hist
    direct_ofstream sam_file; // -sam
    direct_ofstream paf_file; // -paf
    direct_ofstream bed_file; // -bed
};

// the files and the names of the approximate queries of a patterns file
struct apm_query_run_t {
    query_options_t& options; // -p, -o, -m
    direct_ofstream& output_file; // -o
    const std::string& path_index_file;
    const std::string& path_patterns_file;
    std::string built_with; // how an index with sequences is built, e.g. "move-rb build in the FASTA mode"
    std::string measurement_name; // the name of the index in the measurement file (algo = <mode>_<measurement_name>)
};

// ############################# HELP #############################

/**
 * @brief prints the option -rc of locate
 * @param h the help printer
 */
inline void print_apm_input_options(help_printer& h)
{
    print_reverse_complements_option(h, true, "with -sam, -paf or -bed: also search the reverse complement of each "
        "pattern; its occurrences are written with the SAM flag 0x10 and the PAF and BED strand -");
}

/**
 * @brief prints the options -k, -d (locate), -scheme, -bestmap and -x
 * @param h the help printer
 * @param locate true <=> the help of locate
 */
inline void print_apm_search_options(help_printer& h, bool locate)
{
    h.option("-k <errors>", "integer in [0, " + std::to_string(max_scheme_errors) + "] (default: "
        + std::to_string(apm_query_cli_t().max_errors) + ")", "the maximum number of errors");
    if (locate) h.option("-d <distance>", "hamming | edit (default: hamming)", "edit: k <= "
        + std::to_string(max_edit_errors));
    h.option("-scheme <scheme>", "min-u | pigeon-hole | suffix-filter | 01 | <file> (default: "
        + std::string(default_scheme_name) + ")", "the search scheme; a file also sets k");
    h.option("-bestmap", "", "increase the number of errors from 0 up to k until an occurrence is found and report "
        "only the occurrences with the fewest errors");
    h.option("-x <strata>", "integer >= 0 (default: " + std::to_string(apm_query_cli_t().bestmap_strata) + ")",
        "with -bestmap: also report the occurrences with up to this many more errors");
}

/**
 * @brief prints the group of the options of the error model
 * @param h the help printer
 * @param locate true <=> the help of locate
 */
inline void print_apm_model_options(help_printer& h, bool locate)
{
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
}

/**
 * @brief prints the options -sam, -paf, -bed, -seq, -best and -no-cigar (locate) or -hist (count)
 * @param h the help printer
 * @param locate true <=> the help of locate
 * @param built_with how an index with sequences is built, e.g. "move-rb build in the FASTA mode"
 */
inline void print_apm_output_options(help_printer& h, bool locate, const std::string& built_with)
{
    if (locate) {
        h.option("-sam <file>", "a path", "write the occurrences in SAM format, with CIGAR and the tags NM and MD (index "
            "built with " + built_with + ")");
        h.option("-paf <file>", "a path", "write the occurrences in PAF format, one line per occurrence (index built "
            "with " + built_with + "); can be combined with -sam");
        h.option("-bed <file>", "a path", "write the occurrences in BED format (0-based half-open intervals; index "
            "built with " + built_with + "); can be combined with -sam and -paf");
        h.option("-seq", "", "with -sam: write the read (SEQ) and its qualities (QUAL) in the primary record");
        h.option("-best", "", "with -sam: write only the primary record of each read (an occurrence with the fewest "
            "errors)");
        h.option("-no-cigar", "", "with -sam or -paf: no CIGAR is computed or written");
    } else {
        h.option("-hist <file>", "a path", "write one line \"<name>\\t<c_0>\\t<c_1>...\\t<c_k>\" per pattern, c_e = "
            "the number of occurrences with e errors");
    }
}

/**
 * @brief prints the options -partition and -no-seed
 * @param h the help printer
 */
inline void print_apm_execution_options(help_printer& h)
{
    h.option("-partition <partitioning>", "dynamic | uniform (default: dynamic)",
        "how a pattern is partitioned into the parts of the search scheme");
    h.option("-no-seed", "", "match the first part of each search symbol by symbol instead of starting at its search "
        "context computed by the partitioning");
}

/**
 * @brief reads an option of the approximate queries
 * @param args the arguments, positioned after the option
 * @param option the option
 * @param apm the options of the approximate queries
 * @return whether option is one of them
 * @throws usage_error if the value of the option is invalid
 */
inline bool read_apm_query_option(cli_arguments& args, const std::string& option, apm_query_cli_t& apm)
{
    // reads the value of an option of the error model with parse
    auto read_model_option = [&](auto parse) {
        const std::string value = args.value(option);

        try {
            parse(value);
        } catch (const std::exception& error) {
            throw usage_error(option + ": " + error.what());
        }
    };

    apm_model_t& model = apm.apm_model;

    if (option == "-d" && apm.locate) {
        const std::string metric = args.value(option);
        if (metric == "hamming") apm.metric = HAMMING_DISTANCE;
        else if (metric == "edit") apm.metric = EDIT_DISTANCE;
        else throw usage_error("-d expects hamming or edit");
    } else if (option == "-scheme") {
        apm.scheme_name = args.value(option);
    } else if (option == "-k") {
        apm.max_errors = args.integer(option, 0, max_scheme_errors);
    } else if (option == "-partition") {
        const std::string partitioning = args.value(option);
        if (partitioning == "uniform") apm.apm_params.partitioning = _partition_uniform;
        else if (partitioning == "dynamic") apm.apm_params.partitioning = _partition_dynamic;
        else throw usage_error("-partition expects dynamic or uniform");
    } else if (option == "-no-seed") {
        apm.apm_params.start_at_part_contexts = false;
    } else if (option == "-bestmap") {
        apm.bestmap = true;
    } else if (option == "-x") {
        apm.bestmap_strata = args.integer(option);
    } else if (option == "-hist" && !apm.locate) {
        apm.path_hist_file = args.value(option);
    } else if (option == "-sam" && apm.locate) {
        apm.path_sam_file = args.value(option);
    } else if (option == "-paf" && apm.locate) {
        apm.path_paf_file = args.value(option);
    } else if (option == "-bed" && apm.locate) {
        apm.path_bed_file = args.value(option);
    } else if (option == "-seq" && apm.locate) {
        apm.write_seq = true;
    } else if (option == "-rc" && apm.locate) {
        apm.reverse_complements = read_on_off(args, option);
    } else if (option == "-best" && apm.locate) {
        apm.best_only = true;
    } else if (option == "-no-cigar" && apm.locate) {
        apm.no_cigar = true;
    } else if (option == "-cost") {
        read_model_option([&](const std::string& value) { model.costs.emplace_back(parse_apm_cost_rule(value)); });
    } else if (option == "-max") {
        read_model_option([&](const std::string& value) { model.counters.emplace_back(parse_apm_counter(value)); });
    } else if (option == "-class") {
        read_model_option([&](const std::string& value) { model.classes.emplace_back(parse_apm_class(value)); });
    } else if (option == "-tau") {
        read_model_option([&](const std::string& value) {
            model.tau = value == "auto" ? 0 : uint16_t(parse_apm_value(value, "tau")); });
    } else if (option == "-exact") {
        read_model_option([&](const std::string& value) {
            model.costs.emplace_back(apm_cost_rule_t{.ops = APM_SUB | APM_INS | APM_DEL,
                .cost = apm_forbidden, .range = parse_apm_range(value)}); });
    } else if (option == "-cost-dp" && apm.locate) {
        const std::string mode = args.value(option);
        if (mode == "auto") model.cost_dp = APM_COST_DP_AUTO;
        else if (mode == "on") model.cost_dp = APM_COST_DP_ON;
        else if (mode == "prune") model.cost_dp = APM_COST_DP_PRUNE;
        else if (mode == "off") model.cost_dp = APM_COST_DP_OFF;
        else throw usage_error("-cost-dp expects auto, on, prune or off");
    } else if (option == "-iupac") {
        model.iupac = true;
    } else if (option == "-qual") {
        model.quality = true;
    } else {
        return false;
    }

    apm.given = true;
    return true;
}

// ############################# SEARCH SCHEMES #############################

/**
 * @brief builds the search schemes of the queries of the patterns file (see apm_query_cli_t::schemes_by_errors) from
 *        -scheme and -k; a scheme file sets max_errors
 * @param apm the options of the approximate queries
 * @throws usage_error if -scheme names neither a built-in scheme nor a file, or with -bestmap and a scheme file
 */
inline void build_search_schemes(apm_query_cli_t& apm)
{
    const std::string& name = apm.scheme_name;

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
        apm.max_errors = scheme.k;
        if (apm.bestmap) throw usage_error("-bestmap cannot be combined with a search scheme file");
        apm.schemes_by_errors.push_back(scheme);
        return;
    }

    // returns the built-in search scheme for the given number of errors
    auto scheme_of = [&](uint64_t errors) {
        if (name == "pigeon-hole") return pigeon_hole_scheme(errors);
        else if (name == "suffix-filter") return suffix_filter_scheme(errors);
        else if (name == "min-u") return min_u_scheme(errors);
        else return zero_one_scheme(errors);
    };

    if (apm.bestmap || !apm.apm_model.is_default()) {
        for (uint64_t e = 0; e < apm.max_errors; e++) apm.schemes_by_errors.push_back(scheme_of(e));
    }

    apm.schemes_by_errors.push_back(scheme_of(apm.max_errors));
}

// returns the search scheme for at most the given number of errors (see apm_query_cli_t::schemes_by_errors)
inline const search_scheme_t& scheme_for_errors(const apm_query_cli_t& apm, uint64_t errors)
{
    return apm.schemes_by_errors[std::min<uint64_t>(errors, apm.schemes_by_errors.size() - 1)];
}

/**
 * @brief returns the search scheme of a query with at most the given number of errors, or of the error model compiled
 *        for its pattern
 * @param apm the options of the approximate queries
 * @param errors the maximum number of errors
 * @param costs the error model compiled for the pattern (nullptr: none)
 * @throws std::runtime_error if the scheme allows fewer errors than the error model
 */
inline const search_scheme_t& scheme_for_query(const apm_query_cli_t& apm, uint64_t errors, const apm_costs_t* costs)
{
    if (!costs) return scheme_for_errors(apm, errors);
    const search_scheme_t& scheme = scheme_for_errors(apm, costs->max_errors());
    if (scheme.k < costs->max_errors()) throw std::runtime_error("the search scheme allows fewer errors than the "
        "error model");
    return scheme;
}

/**
 * @brief returns the query parameters for a pattern: with an error model, the model compiled for the pattern (it is
 *        valid until the next call on the same thread)
 * @tparam index_t type of the index
 * @param apm the options of the approximate queries
 * @param index the index
 * @param pattern the pattern as it is searched
 * @param qual the qualities of the read the pattern was taken from
 * @param reverse true <=> the pattern is the reverse complement of the read
 * @param bound the maximum cost of an alignment
 * @param metric the distance metric
 */
template <typename index_t>
apm_params_t query_params(const apm_query_cli_t& apm, const index_t& index, const std::string& pattern,
    const std::string& qual, bool reverse, uint64_t bound, distance_metric_t metric)
{
    apm_params_t params = apm.apm_params;
    if (apm.apm_model.is_default()) return params;
    static thread_local apm_costs_t costs;
    costs.compile(apm.apm_model, pattern, qual, reverse, uint16_t(bound), metric, index.forward_index().map_int());
    params.costs = &costs;
    return params;
}

// returns the details of the patterns line of the query summary, e.g. "search scheme min-u, at most 1 error (edit
// distance)"
inline std::string search_description(const apm_query_cli_t& apm)
{
    return "search scheme " + apm.scheme_name + ", at most " + format_number(apm.max_errors)
        + (apm.max_errors == 1 ? " error" : " errors")
        + (apm.metric == HAMMING_DISTANCE ? " (Hamming distance)" : " (edit distance)");
}

// returns the memory allocated during the queries in bytes (0 if malloc_count is not linked)
inline uint64_t additional_memory(uint64_t baseline)
{
    return is_malloc_count_linked ? malloc_count_peak() - baseline : 0;
}

/**
 * @brief checks the options that only apply to the queries of a patterns file and builds the search schemes
 * @param apm the options of the approximate queries
 * @throws usage_error if the options do not fit together
 */
inline void check_apm_query_options(apm_query_cli_t& apm)
{
    build_search_schemes(apm);

    if (apm.metric == EDIT_DISTANCE && apm.max_errors > max_edit_errors)
        throw usage_error("the edit distance allows at most " + std::to_string(max_edit_errors) + " errors");
    if (apm.best_only && apm.path_sam_file.empty()) throw usage_error("-best needs -sam");
}

// opens the files of -hist, -sam, -paf and -bed (see open_output_file)
inline void open_apm_query_files(apm_query_cli_t& apm)
{
    open_output_file(apm.hist_file, apm.path_hist_file);
    open_output_file(apm.sam_file, apm.path_sam_file);
    open_output_file(apm.paf_file, apm.path_paf_file);
    open_output_file(apm.bed_file, apm.path_bed_file);
}

// closes the files of -hist, -sam, -paf and -bed (see close_output_file)
inline void close_apm_query_files(apm_query_cli_t& apm)
{
    close_output_file(apm.hist_file, apm.path_hist_file);
    close_output_file(apm.sam_file, apm.path_sam_file);
    close_output_file(apm.paf_file, apm.path_paf_file);
    close_output_file(apm.bed_file, apm.path_bed_file);
}

// returns whether locate computes the CIGARs of the occurrences: with -sam or -paf, without -no-cigar
inline bool apm_locate_with_cigar(const apm_query_cli_t& apm)
{
    return (apm.sam_file.is_open() || apm.paf_file.is_open()) && !apm.no_cigar;
}

/**
 * @brief throws if an index cannot answer the approximate queries of a patterns file: the edit distance and the error
 *        models rank the symbols with the map of the bytes (see move_r::map_int)
 * @tparam index_t type of the index
 * @param apm the options of the approximate queries
 * @param index the index
 * @param path_index_file the index file
 * @throws std::runtime_error if the index has no byte alphabet but the queries need one
 */
template <typename index_t>
void require_apm_alphabet(const apm_query_cli_t& apm, const index_t& index, const std::string& path_index_file)
{
    if (!index.forward_index().byte_alphabet && (apm.metric == EDIT_DISTANCE || !apm.apm_model.is_default()))
        throw std::runtime_error(path_index_file + " " + edit_distance_needs_byte_alphabet);
}

// ############################# COUNT #############################

/**
 * @brief counts the approximate occurrences of the patterns of the patterns file and prints the summary
 * @tparam index_t type of the index (a move_rb that supports approximate pattern matching)
 * @tparam log_sizes_fnc_t type of log_sizes
 * @param index the index
 * @param apm the options of the approximate queries
 * @param run the files and the names of the queries
 * @param log_sizes called as log_sizes(measurement_file) to append the sizes of the data structures to the measurement
 *        line
 */
template <typename index_t, typename log_sizes_fnc_t>
void run_apm_count_queries(const index_t& index, apm_query_cli_t& apm, const apm_query_run_t& run,
    log_sizes_fnc_t log_sizes)
{
    pattern_reader patterns(run.path_patterns_file);
    require_readable(patterns, run.path_patterns_file);
    const bool write_hist = apm.hist_file.is_open();
    const uint64_t threads = run.options.threads;
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
        run.output_file.is_open() ? &run.output_file : nullptr,
        write_hist ? &apm.hist_file : nullptr } };

    const pattern_totals_t totals = process_patterns(patterns, processing,
        [&](const pattern_chunk& chunk, uint64_t thread) {
            thread_state_t& st = states[thread];

            for (uint64_t j = 0; j < chunk.size(); j++) {
                const pattern_reader::pattern_t& p = chunk.pattern(j);
                const std::string& pattern = p.seq;
                uint64_t count = 0;
                auto start = now();
                apm_params_t params = query_params(apm, index, pattern, p.quality, false, apm.max_errors,
                    HAMMING_DISTANCE);
                const apm_costs_t* costs = params.costs;
                const search_scheme_t& scheme = scheme_for_query(apm, apm.max_errors, costs);

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

                if (run.output_file.is_open()) {
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

    if (num_patterns == 0) throw std::runtime_error(run.path_patterns_file + " contains no patterns");
    const uint64_t pattern_length = totals.symbols / num_patterns;
    print_query_summary(totals, threads, query_summary_t {
        .pattern_details = search_description(apm),
        .results = { format_count_line("occurrences", num_occurrences, num_patterns) },
        .time_query = time_count,
        .additional_memory = additional_memory(baseline_alloc) });

    if (run.options.measurement_file.is_open()) {
        {
            result_line line(run.options.measurement_file);
            line("algo", "count_" + run.measurement_name)
                ("text", run.options.text_name)
                ("n", index.forward_index().input_size())
                ("pattern_length", pattern_length)
                ("num_patterns", num_patterns)
                ("max_mismatches", apm.max_errors)
                ("num_occurrences", num_occurrences)
                ("time_query", time_count)
                ("num_threads", threads)
                ("time_wall", totals.time_wall);
            log_sizes(run.options.measurement_file);
        }

        run.options.measurement_file.close();
    }
}

// ############################# LOCATE #############################

/**
 * @brief locates the occurrences of a read on one strand (with -sam, -paf or -bed), sorts them by position and removes the
 *        edit-distance occurrences that overlap a better one (see filter_edit_distance_occurrences); a read with at
 *        most p symbols (p: the number of parts of the search scheme) has no occurrences
 * @tparam cm whether the occurrences carry a CIGAR
 * @tparam index_t type of the index
 * @param apm the options of the approximate queries
 * @param index the index
 * @param read the read as it is searched
 * @param qual the qualities of the read
 * @param reverse whether read is the reverse complement of the read
 * @param strand is set to the occurrences and their CIGARs
 */
template <cigar_mode_t cm, typename index_t>
void locate_strand(const apm_query_cli_t& apm, const index_t& index, const std::string& read, const std::string& qual,
    bool reverse, strand_occurrences_t<cm>& strand)
{
    strand.clear();

    // locates the occurrences with at most the given number of errors (or the cost bound of the error model) into occ
    auto locate_with_bound = [&](uint64_t errors, std::vector<aprx_occ_t<cm>>& occ) {
        apm_params_t params = query_params(apm, index, read, qual, reverse, errors, apm.metric);
        const apm_costs_t* costs = params.costs;
        const search_scheme_t& scheme = scheme_for_query(apm, errors, costs);
        if ((errors > 0 || costs) && read.size() <= scheme.p) return;
        strand.cigars = index.template locate<cm>(read, scheme, apm.metric, occ, params);
        // the edit-distance occurrences of an error model are already sorted and selected (see apm_params_t)
        if (costs != nullptr && apm.metric == EDIT_DISTANCE) return;
        ips2ra::sort(occ.begin(), occ.end(), key_of_pos {});
        filter_edit_distance_occurrences<cm>(occ, costs ? costs->max_errors() : errors);
    };

    if (apm.bestmap) locate_best<cm>(apm.max_errors, apm.bestmap_strata, strand.occurrences, locate_with_bound);
    else locate_with_bound(apm.max_errors, strand.occurrences);
}

/**
 * @brief locates the approximate occurrences of the patterns of the patterns file and prints the summary
 * @tparam cm whether the occurrences carry a CIGAR (with -sam or -paf without -no-cigar)
 * @tparam index_t type of the index (a move_rb that supports approximate pattern matching and locates all occurrences)
 * @tparam log_sizes_fnc_t type of log_sizes
 * @param index the index
 * @param apm the options of the approximate queries
 * @param run the files and the names of the queries
 * @param log_sizes called as log_sizes(measurement_file) to append the sizes of the data structures to the measurement
 *        line
 */
template <cigar_mode_t cm, typename index_t, typename log_sizes_fnc_t>
void run_apm_locate_queries(const index_t& index, apm_query_cli_t& apm, const apm_query_run_t& run,
    log_sizes_fnc_t log_sizes)
{
    pattern_reader patterns(run.path_patterns_file);
    require_readable(patterns, run.path_patterns_file);
    const bool write_sam = apm.sam_file.is_open();
    const bool write_paf = apm.paf_file.is_open();
    const bool write_bed = apm.bed_file.is_open();
    // true <=> the patterns are reads aligned on both strands
    const bool write_alignments = write_sam || write_paf || write_bed;
    const uint64_t threads = run.options.threads;
    const fasta_sequence_data& sequences = index.seq_data();

    if (write_alignments && !sequences.has_sequences()) {
        throw std::runtime_error(std::string(write_sam ? "-sam" : write_paf ? "-paf" : "-bed") + " needs an index "
            "built with " + run.built_with + " (" + run.path_index_file + " stores no sequences)");
    }

    if (write_sam) write_sam_header(apm.sam_file, sequences, apm.tool, apm.command_line);
    const uint64_t baseline_alloc = malloc_count_current();
    malloc_count_reset_peak();

    // the counters and buffers of a thread
    struct alignas(cache_line_bytes) thread_state_t {
        std::vector<aprx_occ_t<>> occurrences; // the occurrences of the current pattern (without -sam, -paf and -bed)
        strand_occurrences_t<cm> fwd; // with -sam, -paf or -bed: the occurrences of the current read
        strand_occurrences_t<cm> rc; // with -sam, -paf or -bed: the occurrences of the reverse complement of the read
        uint64_t last_sequence_sam = 0; // index of the sequence of the previous SAM record (see sequence_index)
        uint64_t last_sequence_paf = 0; // index of the sequence of the previous PAF record (see sequence_index)
        uint64_t last_sequence_bed = 0; // index of the sequence of the previous BED record (see sequence_index)
        uint64_t num_occurrences = 0;
        uint64_t time_locate = 0;
        uint64_t checksum = 0;
    };

    std::vector<thread_state_t> states(threads);
    enum output_t : uint64_t { occurrences_out, sam_out, paf_out, bed_out }; // the order of processing.outputs
    pattern_processing_t processing { .threads = threads, .outputs = {
        run.output_file.is_open() ? &run.output_file : nullptr,
        write_sam ? &apm.sam_file : nullptr,
        write_paf ? &apm.paf_file : nullptr,
        write_bed ? &apm.bed_file : nullptr } };

    const pattern_totals_t totals = process_patterns(patterns, processing,
        [&](const pattern_chunk& chunk, uint64_t thread) {
            thread_state_t& st = states[thread];
            std::vector<aprx_occ_t<>>& occurrences = st.occurrences;

            for (uint64_t j = 0; j < chunk.size(); j++) {
                const pattern_reader::pattern_t& p = chunk.pattern(j);
                const std::string& qname = p.name;
                const std::string& qual = p.quality;
                const std::string& pattern = p.seq;

                // -sam, -paf or -bed: the read is aligned on both strands
                if (write_alignments) {
                    auto start = now();
                    locate_strand<cm>(apm, index, pattern, qual, false, st.fwd);
                    std::string rc_seq;

                    if (apm.reverse_complements) {
                        rc_seq = reverse_complement(pattern);
                        locate_strand<cm>(apm, index, rc_seq, qual, true, st.rc);
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
                            st.fwd, st.rc, apm.write_seq, apm.best_only);
                    }

                    if (write_paf) {
                        string_appender paf(chunk.out(j, paf_out));
                        write_paf_records<cm>(paf, sequences, st.last_sequence_paf, qname, pattern.size(), st.fwd,
                            false);
                        write_paf_records<cm>(paf, sequences, st.last_sequence_paf, qname, pattern.size(), st.rc,
                            true);
                    }

                    if (write_bed) {
                        string_appender bed(chunk.out(j, bed_out));
                        write_bed_records<cm>(bed, sequences, st.last_sequence_bed, qname, st.fwd, false);
                        write_bed_records<cm>(bed, sequences, st.last_sequence_bed, qname, st.rc, true);
                    }

                    continue;
                }

                auto start = now();

                // locates the occurrences with at most the given number of errors (or the cost bound of the error
                // model) into occ
                auto locate_with_bound = [&](uint64_t errors, std::vector<aprx_occ_t<>>& occ) {
                    apm_params_t params = query_params(apm, index, pattern, qual, false, errors, apm.metric);
                    const apm_costs_t* costs = params.costs;
                    const search_scheme_t& scheme = scheme_for_query(apm, errors, costs);
                    if (costs && pattern.size() <= scheme.p) return;
                    index.locate(pattern, scheme, apm.metric, occ, params);

                    // with an error model, the edit-distance occurrences are already sorted and selected (see
                    // apm_params_t)
                    if (apm.metric == EDIT_DISTANCE && costs == nullptr) {
                        ips2ra::sort(occ.begin(), occ.end(), key_of_pos {});
                        filter_edit_distance_occurrences<>(occ, errors);
                    }
                };

                if (apm.bestmap) locate_best<>(apm.max_errors, apm.bestmap_strata, occurrences, locate_with_bound);
                else locate_with_bound(apm.max_errors, occurrences);
                st.num_occurrences += occurrences.size();
                st.time_locate += time_diff_ns(start, now());

                if (run.output_file.is_open()) {
                    if (apm.metric == HAMMING_DISTANCE) ips2ra::sort(occurrences.begin(), occurrences.end(), key_of_pos {});
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
        .pattern_details = search_description(apm),
        .results = { format_count_line("occurrences", num_occurrences, num_patterns) },
        .time_query = time_locate,
        .result_name = "occurrence",
        .num_results = num_occurrences,
        .additional_memory = additional_memory(baseline_alloc) });

    if (run.options.measurement_file.is_open()) {
        {
            result_line line(run.options.measurement_file);
            line("algo", "locate_" + run.measurement_name)
                ("dist_metric", (apm.metric == HAMMING_DISTANCE ? "hamming" : "edit"))
                ("text", run.options.text_name)
                ("n", index.forward_index().input_size())
                ("pattern_length", pattern_length)
                ("num_patterns", num_patterns)
                ("max_mismatches", apm.max_errors)
                ("num_occurrences", num_occurrences)
                ("time_query", time_locate)
                ("num_threads", threads)
                ("time_wall", totals.time_wall)
                ("checksum", checksum);
            log_sizes(run.options.measurement_file);
        }

        run.options.measurement_file.close();
    }

    if (write_sam) std::cout << "wrote SAM output to " << apm.path_sam_file << std::endl;
    if (write_paf) std::cout << "wrote PAF output to " << apm.path_paf_file << std::endl;
    if (write_bed) std::cout << "wrote BED output to " << apm.path_bed_file << std::endl;
}

#endif // MOVE_R_MISC_APM_QUERY_CLI_HPP
