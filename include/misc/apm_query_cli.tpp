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
 * @file apm_query_cli.tpp
 * @brief the approximate queries of move-rb: printing and reading their options, building the search schemes and
 *        running the queries of a patterns file
 */

#ifndef MOVE_R_MISC_APM_QUERY_CLI_TPP
#define MOVE_R_MISC_APM_QUERY_CLI_TPP

#include <misc/apm_query_cli.hpp>

inline void print_apm_search_options(help_printer& h, bool locate)
{
    h.option("-k <errors>", "integer in [0, " + std::to_string(max_scheme_errors) + "] (default: "
        + std::to_string(apm_query_cli_t().max_errors) + ")", "the maximum number of errors");
    if (locate) h.option("-d <distance>", "hamming | edit (default: hamming)", "edit: k <= "
        + std::to_string(max_edit_errors) + "; an alignment does not start or end with a deletion");
    h.option("-scheme <scheme>", "min-u | pigeon-hole | suffix-filter | 01 | <file> (default: "
        + std::string(default_scheme_name) + ")", "the search scheme (with a file, k is read from the file, and the "
        "scheme must be complete as defined in docs/move-rb.md, section Approximate pattern matching)");
    h.option("-bestmap", "", "increase the number of errors from 0 up to k until an occurrence is found and report "
        "only the occurrences with the fewest errors");
    h.option("-x <strata>", "integer >= 0 (default: " + std::to_string(apm_query_cli_t().bestmap_strata) + ")",
        "with -bestmap: also report the occurrences with up to this many more errors");
}

inline void print_apm_model_options(help_printer& h, bool locate)
{
    h.group("error model (by default, each error has cost 1; an alignment is valid if its total cost is at most k)");
    h.option("-cost <rule>", "<ops>=<cost>[@<range>][/<from>,<to>], repeatable", "set the cost of the operations "
        "<ops> in <range> (for substitutions of a pattern symbol in <from> by a text symbol in <to>; empty: all); "
        "where rules overlap, the later one applies");
    h.option("-max <counter>", "<ops>=<n>[@<range>], n <= " + std::to_string(apm_max_counter_bound) + ", at most "
        + std::to_string(apm_max_counters), "at most n of the operations <ops> in <range>");
    h.option("-exact <range>", "", "no errors in <range> (as -cost sid=inf@<range>)");
    h.option("-class <class>", "<from>,<to>[@<range>], repeatable", "a pattern symbol in <from> also matches the text "
        "symbols in <to>");
    h.option("-models <file>", "a path", "a file with one line \"<name> <rules>\" per pattern: the rules (-cost, -max, "
        "-exact and -class with their values, separated by white space) apply to the pattern <name> after those of the "
        "command line");
    h.option("-iupac", "", "an IUPAC code in a pattern matches each of its bases");
    h.option("-qual", "", "FASTQ patterns: the cost of substituting a pattern symbol is its Phred quality Q from the "
        "quality line of the FASTQ record (the character Q + 33), unless a rule of -cost applies to it");
    h.option("-tau <cost>", "auto | integer >= 1 (default: auto)", "a substitution with a lower cost is not counted "
        "as an error in the search scheme, but its cost is added to the cost of the alignment\n"
        "auto: chosen per pattern such that the error-free part of a search has, on average, as many positions as "
        "possible at which every other text symbol counts as an error");

    h.text(std::string("<ops>: letters s (substitution), i (insertion: a pattern symbol absent from the text), d "
        "(deletion: a text symbol absent from the pattern), o (opening a gap of insertions or deletions); <cost>: an "
        "integer, or inf (the operation is not allowed). <range>: <a>:<b>, <a>:, :<b> or <i>, positions of the "
        "pattern") + (locate ? " (of the read with -sam or -paf)" : "") + "; negative values are counted from the end "
        "(-5: the last 5 symbols). With an error model, the cost of an occurrence is reported as its number of errors"
        + (locate ? ". With -d edit, for every valid alignment, an occurrence with at most the same cost is reported "
        "whose start is at most 2e + 1 positions away from the start of the alignment, where e is the largest number "
        "of errors in a valid alignment." : "."));
}

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
        h.option("-no-cigar", "", "with -sam or -paf: do not compute or write CIGARs");
        h.option("-strings <file>", "a path", "write one line per alignment of a read on a strand (index built with "
            + built_with + "): the name of the read, the strand, the text string, the number of errors, the CIGAR, "
            "the number of occurrences and the number of sequences that contain one; occurrences with the same CIGAR "
            "share a line, and a text symbol that matches through a class (-class, -iupac) is written as the read "
            "symbol");
        h.option("-strings-seqs", "", "with -strings: also write the names of these sequences, separated by commas");
    } else {
        h.option("-hist <file>", "a path", "write one line \"<name>\\t<c_0>\\t<c_1>...\\t<c_k>\" per pattern, c_e = "
            "the number of occurrences with e errors");
    }
}

// returns the cost rule of -exact <range>
inline apm_cost_rule_t exact_rule(const std::string& range)
{ return apm_cost_rule_t{.ops = APM_SUB | APM_INS | APM_DEL, .cost = apm_forbidden, .range = parse_apm_range(range)}; }

// reads the file of -models (see apm_query_cli_t::pattern_models)
inline void read_pattern_models(const std::string& path, apm_query_cli_t& apm)
{
    std::ifstream in(path);
    if (!in) throw usage_error("-models: cannot read " + path);
    std::string line;

    for (uint64_t line_number = 1; std::getline(in, line); line_number++) {
        std::istringstream words(line);
        std::string name, option, value;
        if (!(words >> name)) continue;
        apm_model_t& rules = apm.pattern_models[name];

        try {
            while (words >> option) {
                if (!(words >> value)) throw std::runtime_error(option + " has no value");
                if (option == "-cost") rules.costs.emplace_back(parse_apm_cost_rule(value));
                else if (option == "-max") rules.counters.emplace_back(parse_apm_counter(value));
                else if (option == "-exact") rules.costs.emplace_back(exact_rule(value));
                else if (option == "-class") rules.classes.emplace_back(parse_apm_class(value));
                else throw std::runtime_error("unknown rule " + option);
            }
        } catch (const std::exception& error) {
            throw usage_error("-models: " + path + ", line " + std::to_string(line_number) + ": " + error.what());
        }
    }
}

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
    } else if (option == "-strings" && apm.locate) {
        apm.path_strings_file = args.value(option);
    } else if (option == "-strings-seqs" && apm.locate) {
        apm.strings_sequences = true;
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
        read_model_option([&](const std::string& value) { model.costs.emplace_back(exact_rule(value)); });
    } else if (option == "-models") {
        read_pattern_models(args.value(option), apm);
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

    if (apm.bestmap || has_error_model(apm)) {
        for (uint64_t e = 0; e < apm.max_errors; e++) apm.schemes_by_errors.push_back(scheme_of(e));
    }

    apm.schemes_by_errors.push_back(scheme_of(apm.max_errors));
}

inline const search_scheme_t& scheme_for_query(const apm_query_cli_t& apm, uint64_t errors, const apm_costs_t* costs)
{
    if (!costs) return scheme_for_errors(apm, errors);
    const search_scheme_t& scheme = scheme_for_errors(apm, costs->max_errors());
    if (scheme.k < costs->max_errors()) throw std::runtime_error("the search scheme has fewer errors than the "
        "largest number of errors of the error model");
    return scheme;
}

template <typename index_t>
inline apm_params_t query_params(const apm_query_cli_t& apm, const index_t& index, const std::string& name,
    const std::string& pattern, const std::string& qual, bool reverse, uint64_t bound, distance_metric_t metric)
{
    apm_params_t params;
    const auto rules = apm.pattern_models.find(name);
    if (rules == apm.pattern_models.end() && apm.apm_model.is_default()) return params;
    static thread_local apm_costs_t costs;
    static thread_local apm_model_t model;
    model = apm.apm_model;

    // the rules of the pattern in -models are applied after those of the command line
    if (rules != apm.pattern_models.end()) {
        const apm_model_t& pattern_rules = rules->second;
        model.costs.insert(model.costs.end(), pattern_rules.costs.begin(), pattern_rules.costs.end());
        model.counters.insert(model.counters.end(), pattern_rules.counters.begin(), pattern_rules.counters.end());
        model.classes.insert(model.classes.end(), pattern_rules.classes.begin(), pattern_rules.classes.end());
        if (model.counters.size() > apm_max_counters) throw std::runtime_error("the pattern " + name + " has more "
            "than " + std::to_string(apm_max_counters) + " counters (-max)");
    }

    costs.compile(model, pattern, qual, reverse, uint16_t(bound), metric, index.forward_index().map_int());
    params.costs = &costs;
    return params;
}

inline void check_apm_query_options(apm_query_cli_t& apm)
{
    build_search_schemes(apm);

    if (apm.metric == EDIT_DISTANCE && apm.max_errors > max_edit_errors)
        throw usage_error("with -d edit, -k is at most " + std::to_string(max_edit_errors));
    if (apm.best_only && apm.path_sam_file.empty()) throw usage_error("-best needs -sam");
    if (apm.strings_sequences && apm.path_strings_file.empty()) throw usage_error("-strings-seqs needs -strings");
    if (apm.no_cigar && !apm.path_strings_file.empty()) throw usage_error("-strings needs the CIGARs (no -no-cigar)");
}

template <typename index_t, typename log_sizes_fnc_t>
inline void run_apm_count_queries(const index_t& index, apm_query_cli_t& apm, const apm_query_run_t& run, log_sizes_fnc_t log_sizes)
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
                const pattern_reader::pattern_t& record = chunk.pattern(j);
                const std::string& pattern = record.seq;
                uint64_t count = 0;
                auto start = now();
                apm_params_t params = query_params(apm, index, record.name, pattern, record.quality, false,
                    apm.max_errors, HAMMING_DISTANCE);
                const apm_costs_t* costs = params.costs;
                const search_scheme_t& scheme = scheme_for_query(apm, apm.max_errors, costs);

                if ((apm.max_errors > 0 || costs) && !fits_search_scheme(pattern.size(), scheme)) {
                    throw std::runtime_error(too_short_for_scheme("the pattern " + record.name, pattern.size(), scheme));
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
                    out << record.name << '\t' << count << '\n';
                }

                if (write_hist) {
                    string_appender out(chunk.out(j, hist_out));
                    out << record.name;
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

    const uint64_t pattern_length = num_patterns == 0 ? 0 : totals.symbols / num_patterns;
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

template <cigar_mode_t cm, typename index_t>
inline void locate_strand(const apm_query_cli_t& apm, const index_t& index, const std::string& name,
    const std::string& read, const std::string& qual, bool reverse, strand_occurrences_t<cm>& strand)
{
    strand.clear();

    // locates the occurrences with at most the given number of errors (or the cost bound of the error model) into occ
    auto locate_with_bound = [&](uint64_t errors, std::vector<aprx_occ_t<cm>>& occ) {
        apm_params_t params = query_params(apm, index, name, read, qual, reverse, errors, apm.metric);
        const apm_costs_t* costs = params.costs;
        const search_scheme_t& scheme = scheme_for_query(apm, errors, costs);
        if ((errors > 0 || costs) && !fits_search_scheme(read.size(), scheme)) {
            throw std::runtime_error(too_short_for_scheme("the read " + name, read.size(), scheme));
        }

        strand.cigars = index.template locate<cm>(read, scheme, apm.metric, occ, params);
        // the edit-distance occurrences of an error model are already sorted and selected (see apm_params_t)
        if (costs != nullptr && apm.metric == EDIT_DISTANCE) return;
        ips2ra::sort(occ.begin(), occ.end(), key_of_pos {});
        filter_edit_distance_occurrences<cm>(occ, costs ? costs->max_errors() : errors);
    };

    if (apm.bestmap) locate_best<cm>(apm.max_errors, apm.bestmap_strata, strand.occurrences, locate_with_bound);
    else locate_with_bound(apm.max_errors, strand.occurrences);
}

template <cigar_mode_t cm, typename index_t, typename log_sizes_fnc_t>
inline void run_apm_locate_queries(const index_t& index, apm_query_cli_t& apm, const apm_query_run_t& run, log_sizes_fnc_t log_sizes)
{
    pattern_reader patterns(run.path_patterns_file);
    require_readable(patterns, run.path_patterns_file);
    const bool write_sam = apm.sam_file.is_open();
    const bool write_paf = apm.paf_file.is_open();
    const bool write_bed = apm.bed_file.is_open();
    const bool write_strings = apm.strings_file.is_open();
    // true <=> the patterns are reads whose occurrences are written as alignments
    const bool write_alignments = write_sam || write_paf || write_bed || write_strings;
    apm.apm_model.writes_occurrences = write_alignments || run.output_file.is_open();
    const uint64_t threads = run.options.threads;
    const fasta_sequence_data& sequences = index.seq_data();

    if (write_alignments && !sequences.has_sequences()) {
        throw std::runtime_error(std::string(write_sam ? "-sam" : write_paf ? "-paf" : write_bed ? "-bed" : "-strings")
            + " needs an index built with " + run.built_with + " (" + run.path_index_file + " stores no sequences)");
    }

    if (write_sam) write_sam_header(apm.sam_file, sequences, apm.tool, apm.command_line);
    const uint64_t baseline_alloc = malloc_count_current();
    malloc_count_reset_peak();

    // the counters and buffers of a thread
    struct alignas(cache_line_bytes) thread_state_t {
        std::vector<aprx_occ_t<>> occurrences; // the occurrences of the current pattern (without -sam, -paf, -bed and -strings)
        strand_occurrences_t<cm> fwd; // with -sam, -paf, -bed or -strings: the occurrences of the current read
        strand_occurrences_t<cm> rc; // with -sam, -paf, -bed or -strings: the occurrences of the reverse complement of the read
        uint64_t last_sequence_sam = 0; // index of the sequence of the previous SAM record (see sequence_index)
        uint64_t last_sequence_paf = 0; // index of the sequence of the previous PAF record (see sequence_index)
        uint64_t last_sequence_bed = 0; // index of the sequence of the previous BED record (see sequence_index)
        uint64_t last_sequence_strings = 0; // index of the sequence of the previous -strings occurrence
        string_groups_t string_groups; // the buffers of -strings
        uint64_t num_occurrences = 0;
        uint64_t time_locate = 0;
        uint64_t checksum = 0;
    };

    std::vector<thread_state_t> states(threads);
    // the order of processing.outputs
    enum output_t : uint64_t { occurrences_out, sam_out, paf_out, bed_out, strings_out };
    pattern_processing_t processing { .threads = threads, .outputs = {
        run.output_file.is_open() ? &run.output_file : nullptr,
        write_sam ? &apm.sam_file : nullptr,
        write_paf ? &apm.paf_file : nullptr,
        write_bed ? &apm.bed_file : nullptr,
        write_strings ? &apm.strings_file : nullptr } };

    const pattern_totals_t totals = process_patterns(patterns, processing,
        [&](const pattern_chunk& chunk, uint64_t thread) {
            thread_state_t& st = states[thread];
            std::vector<aprx_occ_t<>>& occurrences = st.occurrences;

            for (uint64_t j = 0; j < chunk.size(); j++) {
                const pattern_reader::pattern_t& record = chunk.pattern(j);
                const std::string& qname = record.name;
                const std::string& qual = record.quality;
                const std::string& pattern = record.seq;

                // -sam, -paf, -bed or -strings: the read is aligned on both strands (with -rc)
                if (write_alignments) {
                    auto start = now();
                    locate_strand<cm>(apm, index, qname, pattern, qual, false, st.fwd);
                    std::string rc_seq;

                    if (apm.reverse_complements) {
                        rc_seq = reverse_complement(pattern);
                        locate_strand<cm>(apm, index, qname, rc_seq, qual, true, st.rc);
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

                    if constexpr (cm == CIGAR) if (write_strings) {
                        string_appender strings(chunk.out(j, strings_out));
                        write_string_records(strings, sequences, st.last_sequence_strings, qname, pattern, st.fwd,
                            false, apm.strings_sequences, st.string_groups);
                        write_string_records(strings, sequences, st.last_sequence_strings, qname, rc_seq, st.rc, true,
                            apm.strings_sequences, st.string_groups);
                    }

                    chunk.write(j);
                    continue;
                }

                auto start = now();

                /* locates the occurrences with at most the given number of errors (or the cost bound of the error
                   model) into occ */
                auto locate_with_bound = [&](uint64_t errors, std::vector<aprx_occ_t<>>& occ) {
                    apm_params_t params = query_params(apm, index, qname, pattern, qual, false, errors, apm.metric);
                    const apm_costs_t* costs = params.costs;
                    const search_scheme_t& scheme = scheme_for_query(apm, errors, costs);
                    if ((errors > 0 || costs) && !fits_search_scheme(pattern.size(), scheme)) {
                        throw std::runtime_error(too_short_for_scheme("the pattern " + qname, pattern.size(), scheme));
                    }

                    index.locate(pattern, scheme, apm.metric, occ, params);

                    // with an error model, the edit-distance occurrences are already sorted and selected (see apm_params_t)
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
                    chunk.write(j);
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
    if (write_strings) std::cout << "wrote the strings to " << apm.path_strings_file << std::endl;
}

#endif // MOVE_R_MISC_APM_QUERY_CLI_TPP
