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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY type, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

// move-rb-gen-apm-queries: generates the pattern files of move-rb-bench-apm. For each k and m it writes the files
//   <text_name>.patterns-count-hamming-k<k>-m<m>
//   <text_name>.patterns-locate-hamming-k<k>-m<m>
//   <text_name>.patterns-locate-edit-k<k>-m<m>
// of substrings of the text at random positions. Patterns are added to a file until their queries on a move_rb index
// take the target time (--time) in total.

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <clocale>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <poll.h>
#include <sys/wait.h>
#include <unistd.h>

#include <misc/index_cli.hpp>
#include <misc/move_rb_cli.hpp>
#include <move_rb/move_rb.hpp>

#include "bench_common.hpp"
#include "timeout_runner.hpp"

#include <ips4o.hpp>
#include <ips2ra.hpp>

// ############################# CONFIGURATION #############################

std::string path_input_file;         // the text the index has been built for
std::string path_index_file;         // a move_rb index of the text
double      target_seconds = default_target_seconds; // --time: total time in seconds of the queries of a pattern file

std::string out_dir = ".";           // -o: directory of the pattern files
std::string text_name;               // -n: name of the text in the file names (default: file name of the text)
std::string scheme_str = "min_u";    // -s: search scheme (see search_scheme_of)
std::string forbidden;               // -x: symbols a pattern must not contain (e.g. a FASTA separator)
bool gen_count = true;               // --only: generate the count pattern files
bool gen_locate = true;              // --only: generate the locate pattern files
uint64_t min_patterns = default_min_patterns; // --min: minimum number of patterns per file
uint64_t seed = 0;                   // --seed: seed for the positions of the patterns (0: random)
double   timeout_seconds = default_timeout_seconds; // --timeout: timeout of a query in seconds (0: none)
// --max-skips: a file is given up after this many timeouts in a row
uint64_t max_consecutive_skips = default_max_consecutive_skips;

// the default of -k (measurement setup)
const std::vector<int64_t> default_error_counts = { 4, 7, 10, 13 };
std::vector<int64_t>  ks = default_error_counts; // -k: the numbers of errors
std::vector<uint64_t> ms = default_pattern_lengths; // -m: the pattern lengths

direct_ifstream index_file;

// an operation of a pattern file
struct query_type_t {
    const char* op; // count or locate
    const char* metric; // hamming or edit
    bool is_count; // true <=> op is count
    bool is_edit; // true <=> metric is edit
};

// ############################# HELPERS #############################

[[noreturn]] void help(const std::string& msg = "")
{
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "move-rb-gen-apm-queries: generate the query (pattern) sets for move-rb-bench-apm, auto-calibrating" << std::endl;
    std::cout << "                 the number of patterns per file (against a move_rb index) so that" << std::endl;
    std::cout << "                 running that file's operation takes ~<time> seconds (see --time)." << std::endl << std::endl;
    std::cout << "                 For every (k, m) it writes three files, one per benchmark operation:" << std::endl;
    std::cout << "                 <text_name>.patterns-count-hamming, -locate-hamming and -locate-edit (-k<k>-m<m>)" << std::endl << std::endl;
    std::cout << "usage: move-rb-gen-apm-queries [options] <input_file> <move_rb_index>" << std::endl;
    std::cout << "   -o <dir>          output directory for the pattern files (default: current directory)" << std::endl;
    std::cout << "   -n <text_name>    text name used in the file names (default: basename of <input_file>)" << std::endl;
    std::cout << "   -x <chars>        characters a generated pattern must not contain (e.g. a FASTA separator or" << std::endl;
    std::cout << "                     the terminator); such patterns are skipped (default: none)" << std::endl;
    std::cout << "   -k <list>         comma-separated list of error counts k (default: " << join_csv(default_error_counts)
              << ")" << std::endl;
    std::cout << "   -m <list>         comma-separated list of pattern lengths (default: "
              << join_csv(default_pattern_lengths) << ")" << std::endl;
    std::cout << "   -s <scheme>       search scheme: pigeon_hole, suffix_filter, min_u or 01 (default: min_u)" << std::endl;
    std::cout << "   --time <s>        target count/locate duration per pattern file, in seconds (default: "
              << default_target_seconds << "s)" << std::endl;
    std::cout << "   --only <op>       restrict to count or locate files: count, locate or both (default: both)" << std::endl;
    std::cout << "   --min <N>         minimum number of patterns per file (default: " << default_min_patterns << ")"
              << std::endl;
    std::cout << "   --seed <s>        seed for the random pattern positions (default: random)" << std::endl;
    std::cout << "   --timeout <s>     per-pattern query timeout in seconds; a query runs in a forked worker that is" << std::endl;
    std::cout << "                     killed if it exceeds the timeout, then the pattern is discarded (0 = off, default: "
              << default_timeout_seconds << ")" << std::endl;
    std::cout << "   --max-skips <N>   give up on a (k,m) after this many consecutive timeouts (default: "
              << default_max_consecutive_skips << ")" << std::endl;
    std::cout << "   <input_file>      the text the index was built for" << std::endl;
    std::cout << "   <move_rb_index>   a move_rb index of <input_file> that supports locate" << std::endl;
    std::exit(0);
}

// ############################# GENERATION #############################

/**
 * @brief generates and writes all pattern files; patterns are added to a file until their queries on the index take
 *        target_seconds in total
 * @tparam support the locate support of the index
 */
template <move_r_support support>
void generate()
{
    using idx_t = move_rb<support>;
    idx_t loaded;
    idx_t* index = &loaded;
    load_index_with_progress(*index, index_file, "move_rb index");
    index_file.close();
    drop_from_page_cache(path_index_file, true);

    direct_ifstream tf(path_input_file, std::ios::binary);
    tf.seekg(0, std::ios::end);
    uint64_t text_size = tf.tellg();

    const uint64_t target_ns = (uint64_t) (target_seconds * 1e9);
    const uint64_t timeout_ns = (uint64_t) (timeout_seconds * 1e9);
    std::mt19937_64 rng(seed != 0 ? seed : std::random_device{}());

    std::vector<query_type_t> types;
    if (gen_count)  { types.push_back({"count",  "hamming", true,  false}); }
    if (gen_locate) { types.push_back({"locate", "hamming", false, false});
                      types.push_back({"locate", "edit",    false, true }); }

    std::cout << "calibrating to ~" << target_seconds << "s per file\n" << std::endl;

    for (const query_type_t& type : types) {
        for (int64_t k : ks) {
            search_scheme_t scheme = search_scheme_of(scheme_str, k);

            for (uint64_t m : ms) {
                if (m < scheme.p) {
                    std::cerr << "skipping " << type.op << "-" << type.metric << " k" << k << " m" << m
                              << ": pattern length < " << (int32_t) scheme.p << " parts in the search scheme" << std::endl;
                    continue;
                }

                std::string pattern;
                resize_no_init(pattern, m);
                std::string buffer; // the concatenated patterns to write
                std::uniform_int_distribution<uint64_t> pos_distrib(0, text_size - m - 1);
                uint64_t n = 0, total_ns = 0, total_occ = 0, skipped = 0, consecutive_skips = 0;

                auto run_query = [index, scheme, type](const std::string& pat) -> uint64_t {
                    if (type.is_count) {
                        return index->count_hamming_dist(pat, scheme);
                    }
                    std::vector<aprx_occ_t<>> occ;
                    if (type.is_edit) {
                        index->locate(pat, scheme, EDIT_DISTANCE, occ);
                        ips2ra::sort(occ.begin(), occ.end(), key_of_pos {});
                        filter_edit_distance_occurrences<>(occ, (uint64_t) scheme.k);
                    } else {
                        index->locate(pat, scheme, HAMMING_DISTANCE, occ);
                    }
                    return occ.size();
                };

                // runs the query of a pattern and measures its time
                auto timed_query = [&run_query](const std::string& pat) {
                    auto t1 = now();
                    uint64_t c = run_query(pat);
                    return timed_result_t { c, time_diff_ns(t1) };
                };

                std::optional<timeout_runner<decltype(timed_query)>> runner;
                if (timeout_ns != 0) runner.emplace(timed_query, m);

                // add patterns until their queries take target_ns in total (but at least min_patterns patterns)
                while (n < min_patterns || total_ns < target_ns) {
                    tf.seekg(pos_distrib(rng), std::ios::beg);
                    tf.read(pattern.data(), m);

                    // skip a pattern that contains a forbidden symbol
                    if (!forbidden.empty() && pattern.find_first_of(forbidden) != std::string::npos) continue;

                    std::optional<timed_result_t> r = runner ? runner->run(pattern, timeout_ns) : timed_query(pattern);

                    // the query timed out: discard the pattern; give up after max_consecutive_skips timeouts in a row
                    if (!r.has_value()) {
                        skipped++;
                        if (++consecutive_skips > max_consecutive_skips) {
                            std::cerr << "giving up on " << type.op << "-" << type.metric
                                      << " k" << k << " m" << m << " after " << consecutive_skips
                                      << " consecutive timeouts" << std::endl;
                            break;
                        }
                        continue;
                    }
                    consecutive_skips = 0;

                    total_occ += r->count;
                    total_ns += r->elapsed_ns;

                    buffer.append(pattern);
                    n++;
                }

                if (n < min_patterns) {
                    std::cerr << "warning: not enough patterns for " << type.op << "-" << type.metric
                              << " k" << k << " m" << m << " completed within the timeout; "
                              << "skipping" << std::endl;
                    continue;
                }

                std::string fname = (std::filesystem::path(out_dir) / (text_name + ".patterns-" +
                    type.op + "-" + type.metric + "-k" + std::to_string(k) + "-m" + std::to_string(m))).string();
                write_pattern_file(fname, buffer, n, m, text_name, forbidden);

                std::cout << type.op << "-" << type.metric << " k" << k << " m" << m
                          << ": " << n << " patterns";
                if (skipped > 0) std::cout << " (" << skipped << " skipped)";
                std::cout << ", " << (total_ns / 1e9) << "s";
                if (n > 0) std::cout << " (" << (total_ns / n) / 1e3 << " us/pattern, "
                                     << (total_occ / n) << " occ/pattern)";
                std::cout << std::endl;
            }
        }
    }

    tf.close();
    drop_from_page_cache(path_input_file, true);

    std::cout << "\ndone; pattern files written to " << out_dir << std::endl;
}

// ############################# DRIVER #############################

int main(int argc, char** argv)
{
    std::setlocale(LC_ALL, "C");
    std::signal(SIGPIPE, SIG_IGN);

    int32_t i = 1;
    while (i < argc && argv[i][0] == '-') {
        std::string o = argv[i++];
        auto need = [&](const char* opt) { if (i >= argc) help(std::string("error: missing argument after ") + opt); return std::string(argv[i++]); };
        if      (o == "-o") out_dir = need("-o");
        else if (o == "-n") text_name = need("-n");
        else if (o == "-s") scheme_str = need("-s");
        else if (o == "-x") forbidden = need("-x");
        else if (o == "--min")  min_patterns = std::stoull(need("--min"));
        else if (o == "--seed") seed = std::stoull(need("--seed"));
        else if (o == "--time")      { if (!parse_double(need("--time"), target_seconds)) help("error: --time expects a number"); }
        else if (o == "--timeout")   { if (!parse_double(need("--timeout"), timeout_seconds)) help("error: --timeout expects a number"); }
        else if (o == "--max-skips") max_consecutive_skips = std::stoull(need("--max-skips"));
        else if (o == "-k") { ks.clear(); for (auto& v : split_csv(need("-k"))) ks.push_back(std::stoll(v)); }
        else if (o == "-m") { ms.clear(); for (auto& v : split_csv(need("-m"))) ms.push_back(std::stoull(v)); }
        else if (o == "--only") {
            std::string ops = need("--only");
            if      (ops == "count")  { gen_count = true;  gen_locate = false; }
            else if (ops == "locate") { gen_count = false; gen_locate = true; }
            else if (ops == "both")   { gen_count = true;  gen_locate = true; }
            else help("error: --only must be count, locate or both");
        } else {
            help("error: unrecognized option '" + o + "'");
        }
    }

    if (argc - i != 2) help();
    path_input_file = argv[i++];
    path_index_file = argv[i++];

    if (target_seconds <= 0) help("error: --time must be > 0");
    if (timeout_seconds < 0) help("error: --timeout must be >= 0");
    if (min_patterns < 1) min_patterns = 1;
    if (!is_search_scheme_name(scheme_str)) help("error: -s must be pigeon_hole, suffix_filter, min_u or 01");
    if (!std::filesystem::exists(path_input_file)) help("error: <input_file> does not exist");
    if (text_name.empty()) text_name = std::filesystem::path(path_input_file).filename().string();
    if (!std::filesystem::is_directory(out_dir)) help("error: output directory '" + out_dir + "' does not exist");

    index_file.open(path_index_file, std::ios::binary);
    if (!index_file.good()) help("error: cannot open <move_rb_index>");

    const move_r_support support = read_move_rb_support(index_file);
    if (support == _count)
        help("error: this index only supports count, not locate; rebuild with locate support");

    if (!dispatch_move_rb_index<_require_all_occs>(support,
        [&]<move_r_support index_support>() { generate<index_support>(); }))
        help("error: unsupported support type in <move_rb_index>");

    return 0;
}
