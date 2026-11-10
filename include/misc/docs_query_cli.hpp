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
 * @file docs_query_cli.hpp
 * @brief the mode docs of move-rb: its options, their help and the document queries of a patterns file
 */

#ifndef MOVE_R_MISC_DOCS_QUERY_CLI_HPP
#define MOVE_R_MISC_DOCS_QUERY_CLI_HPP

#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <algorithms/documents.hpp>
#include <misc/cli_help.hpp>
#include <misc/index_cli.hpp>
#include <misc/memory.hpp>
#include <misc/parallel_patterns.hpp>
#include <misc/pattern_reader.hpp>
#include <misc/query_cli.hpp>

// ############################# OPTIONS #############################

// what is reported for each pattern (-report)
enum class docs_report_t { list, count, frequencies, top_k };

// the names of docs_report_t in the RESULT line: [0..3] the name of docs_report_t(x)
inline constexpr const char* docs_report_names[] = { "list", "count", "frequencies", "top_k" };

/** @brief the options of the mode docs */
struct docs_query_cli_t {
    docs_report_t report = docs_report_t::list; // -report
    uint64_t top_k = 0; // -top (0: not given)
};

// the files and the names of the document queries of a patterns file
struct docs_query_run_t {
    query_options_t& options; // -p, -o, -m
    std::ofstream& output_file; // -o
    const std::string& path_index_file;
    const std::string& path_patterns_file;
    std::string measurement_name; // the name of the index in the measurement file (algo = docs_<report>_<name>)
};

// ############################# HELP #############################

/**
 * @brief prints the options -report and -top
 * @param h the help printer
 */
inline void print_docs_query_options(help_printer& h)
{
    h.option("-report <what>", "list | count | frequencies | top (default: list)",
        "list: the documents containing the pattern\n"
        "count: their number\n"
        "frequencies: the documents with the number of occurrences of the pattern in each\n"
        "top: the -top documents with the most occurrences");
    h.option("-top <k>", "integer >= 1", "with -report top (and only with it): the number of documents");
}

/**
 * @brief prints the option -o of the mode docs
 * @param h the help printer
 */
inline void print_docs_output_option(help_printer& h)
{
    h.option("-o <file>", "a path", "write one line \"<name>\\t<occurrences>\\t<documents>\\t<names>\" per pattern, "
        "where <names> are the names of the documents separated by commas (with -report frequencies and top as "
        "<name>:<occurrences>)");
}

// ############################# READING THE OPTIONS #############################

/**
 * @brief reads an option of the mode docs
 * @param args the arguments, positioned after the option
 * @param option the option
 * @param docs the options of the mode docs
 * @return whether option is one of them
 * @throws usage_error if the value of the option is invalid
 */
inline bool read_docs_query_option(cli_arguments& args, const std::string& option, docs_query_cli_t& docs)
{
    if (option == "-report") {
        const std::string report = args.value(option);
        if (report == "list") docs.report = docs_report_t::list;
        else if (report == "count") docs.report = docs_report_t::count;
        else if (report == "frequencies") docs.report = docs_report_t::frequencies;
        else if (report == "top") docs.report = docs_report_t::top_k;
        else throw usage_error("-report expects list, count, frequencies or top");
    } else if (option == "-top") {
        docs.top_k = args.integer(option, 1);
    } else {
        return false;
    }

    return true;
}

/**
 * @brief checks that the options of the mode docs fit together
 * @param docs the options of the mode docs
 * @throws usage_error if they do not
 */
inline void check_docs_query_options(const docs_query_cli_t& docs)
{
    if (docs.report == docs_report_t::top_k && docs.top_k == 0) throw usage_error("-report top needs -top <k>");
    if (docs.report != docs_report_t::top_k && docs.top_k != 0) throw usage_error("-top needs -report top");
}

// ############################# QUERIES #############################

/**
 * @brief finds the documents that contain the patterns of the patterns file, writes them and prints the summary
 * @tparam search_fnc_t type of search
 * @tparam log_sizes_fnc_t type of log_sizes
 * @param docs the documents of the index
 * @param n the length of the text
 * @param cli the options of the mode docs
 * @param run the files and the names of the queries
 * @param search called as search(P, found) for a nonempty pattern P; if P occurs, it calls found(b, e, locate) with
 *        the suffix array interval [b, e] of P and a function locate, called as locate(report), that calls
 *        report(occ) for the text position occ of every occurrence of P
 * @param log_sizes called as log_sizes(measurement_file) to append the sizes of the data structures to the measurement
 *        line
 */
template <typename search_fnc_t, typename log_sizes_fnc_t>
void run_docs_queries(const document_array& docs, uint64_t n, const docs_query_cli_t& cli, const docs_query_run_t& run,
    search_fnc_t search, log_sizes_fnc_t log_sizes)
{
    pattern_reader reader(run.path_patterns_file);
    require_readable(reader, run.path_patterns_file);
    const uint64_t threads = run.options.threads;
    const docs_report_t report = cli.report;
    const uint64_t baseline_alloc = malloc_count_current();
    malloc_count_reset_peak();

    // the counters and buffers of a thread
    struct alignas(cache_line_bytes) thread_state_t {
        document_query query;
        uint64_t num_occurrences = 0;
        uint64_t num_documents = 0;
        uint64_t checksum = 0;
        uint64_t time_search = 0;
        uint64_t time_documents = 0;
        // (document, frequency) of the documents that contain the current pattern, in the order document_query
        // reports them (-report top: by decreasing frequency); frequency = 0 with -report list
        std::vector<std::pair<uint64_t, uint64_t>> found;
    };

    std::vector<thread_state_t> states;
    states.reserve(threads);

    for (uint64_t t = 0; t < threads; t++)
        states.push_back({ .query = document_query(docs) });

    pattern_processing_t processing { .threads = threads,
        .outputs = { run.output_file.is_open() ? &run.output_file : nullptr } };

    const pattern_totals_t totals = process_patterns(reader, processing,
        [&](const pattern_chunk& chunk, uint64_t thread) {
            thread_state_t& st = states[thread];
            document_query& query = st.query;
            std::vector<std::pair<uint64_t, uint64_t>>& found = st.found;

            for (uint64_t j = 0; j < chunk.size(); j++) {
                const pattern_reader::pattern_t& pattern = chunk.pattern(j);
                const std::string& P = pattern.seq;
                found.clear();
                uint64_t occ = 0; // the number of occurrences of the pattern
                uint64_t num = 0; // the number of documents containing the pattern
                auto start = now();
                auto searched = start;
                bool occurs = false;

                // computes the documents of the pattern with the suffix array interval [b, e]
                auto documents_of = [&](uint64_t b, uint64_t e, auto&& locate) {
                    searched = now();
                    occurs = true;
                    occ = e - b + 1;

                    if (report == docs_report_t::count) {
                        num = query.count(locate);
                    } else if (report == docs_report_t::list) {
                        query.list(locate, [&](uint64_t document) { found.emplace_back(document, 0); });
                    } else if (report == docs_report_t::frequencies) {
                        query.list_with_frequencies(locate, [&](uint64_t document, uint64_t frequency) {
                            found.emplace_back(document, frequency); });
                    } else {
                        query.top_k(locate, cli.top_k, [&](uint64_t document, uint64_t frequency) {
                            found.emplace_back(document, frequency); });
                    }

                    if (report != docs_report_t::count) num = found.size();
                };

                if (!P.empty()) search(P, documents_of);
                const auto end = now();
                if (!occurs) searched = end;
                st.time_search += time_diff_ns(start, searched);
                st.time_documents += time_diff_ns(searched, end);
                st.num_occurrences += occ;
                st.num_documents += num;
                st.checksum += num;

                if (report != docs_report_t::count) {
                    for (const auto& [document, frequency] : found) st.checksum += document + frequency;
                }

                if (run.output_file.is_open()) {
                    string_appender out(chunk.out(j));
                    out << pattern.name << '\t' << occ << '\t' << num << '\t';

                    if (report != docs_report_t::count) {
                        for (uint64_t i = 0; i < found.size(); i++) {
                            if (i != 0) out << ',';
                            out << docs.documents().name(found[i].first);
                            if (report == docs_report_t::frequencies || report == docs_report_t::top_k) out << ':' << found[i].second;
                        }
                    }

                    out << '\n';
                    chunk.write(j);
                }
            }
        });

    if (totals.patterns == 0) throw std::runtime_error(run.path_patterns_file + " contains no patterns");
    const uint64_t num_patterns = totals.patterns;
    uint64_t num_occurrences = 0;
    uint64_t num_documents = 0;
    uint64_t checksum = 0;
    uint64_t time_search = 0;
    uint64_t time_documents = 0;

    for (const thread_state_t& st : states) {
        num_occurrences += st.num_occurrences;
        num_documents += st.num_documents;
        checksum += st.checksum;
        time_search += st.time_search;
        time_documents += st.time_documents;
    }

    const std::string summed = threads > 1 ? " (both summed over " + format_threads(threads) + ")" : "";
    print_query_summary(totals, threads, query_summary_t {
        .pattern_details = "report " + (report == docs_report_t::top_k ? "top " + format_number(cli.top_k)
            : std::string(docs_report_names[uint64_t(report)])),
        .results = {
            format_count_line("occurrences", num_occurrences, num_patterns),
            format_count_line("documents", num_documents, num_patterns),
            "search time: " + format_time(time_search) + ", document time: " + format_time(time_documents)
                + (num_documents == 0 ? "" : " = " + format_time(time_documents / num_documents) + " per document")
                + summed },
        .time_query = time_search + time_documents,
        .additional_memory = is_malloc_count_linked ? malloc_count_peak() - baseline_alloc : 0 });

    if (run.options.measurement_file.is_open()) {
        {
            result_line line(run.options.measurement_file);
            line("algo", std::string("docs_") + docs_report_names[uint64_t(report)] + "_" + run.measurement_name)
                ("text", run.options.text_name)
                ("n", n)
                ("num_documents_text", docs.num_documents())
                ("num_patterns", num_patterns)
                ("pattern_length", totals.symbols / num_patterns)
                ("num_occurrences", num_occurrences)
                ("num_documents", num_documents)
                ("time_search", time_search)
                ("time_documents", time_documents)
                ("num_threads", threads)
                ("time_wall", totals.time_wall)
                ("checksum", checksum);
            log_sizes(run.options.measurement_file);
        }

        run.options.measurement_file.close();
    }
}

#endif // MOVE_R_MISC_DOCS_QUERY_CLI_HPP
