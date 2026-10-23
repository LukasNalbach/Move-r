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

// move-rb-gen-ext-queries: generates the pattern files of move-rb-bench-ext. For each m it writes the files
// <text_name>.patterns-count-m<m> and <text_name>.patterns-locate-m<m> of substrings of the text at random positions.
// Patterns are added to a file until their queries on a move_rb index take the target time (--time) in total. The
// extension plans are not stored: move-rb-bench-ext draws the same plans (see ext_query.hpp).

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <clocale>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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

#include "bench_common.hpp"
#include "ext_query.hpp"
#include "timeout_runner.hpp"
#include <misc/index_cli.hpp>
#include <misc/move_rb_cli.hpp>
#include <move_rb/move_rb.hpp>

// ############################# CONFIGURATION #############################

std::string path_input_file;         // the text the index has been built for
std::string path_index_file;         // a move_rb index of the text
double      target_seconds = default_target_seconds; // --time: total time in seconds of the queries of a pattern file

std::string out_dir = ".";           // -o: directory of the pattern files
std::string text_name;               // -n: name of the text in the file names (default: file name of the text)
bool gen_count = true;               // --only: generate the count pattern files
bool gen_locate = true;              // --only: generate the locate pattern files
std::string forbidden;               // -x: symbols a pattern must not contain (e.g. a FASTA separator)
uint64_t min_patterns = default_min_patterns; // --min: minimum number of patterns per file
double   timeout_seconds = default_timeout_seconds; // --timeout: timeout of a query in seconds (0: none)
// --max-skips: a file is given up after this many timeouts in a row
uint64_t max_consecutive_skips = default_max_consecutive_skips;

std::vector<uint64_t> ms = default_pattern_lengths; // -m: the pattern lengths

// an operation of a pattern file
struct ext_op_t {
    const char* op; // count or locate
    bool is_locate; // true <=> op is locate
};

/* a request to the child of timeout_runner: [0..m-1] the pattern, [m..m+7] plan.start, [m+8..2m+6] the directions of
   plan.order (0 = LEFT, 1 = RIGHT) */

// returns the number of bytes of a request for a pattern of length m
inline uint64_t request_size(uint64_t m) { return m + sizeof(uint64_t) + (m == 0 ? 0 : m - 1); }

// encodes a pattern and its extension plan as a request
inline void encode_request(const std::string& pat, const ext_plan_t& plan, std::string& request)
{
    const uint64_t m = pat.size();
    request.resize(request_size(m));
    std::memcpy(request.data(), pat.data(), m);
    std::memcpy(request.data() + m, &plan.start, sizeof(uint64_t));
    for (uint64_t i = 0; i + 1 < m; i++) request[m + sizeof(uint64_t) + i] = plan.order[i] == LEFT ? 0 : 1;
}

// decodes a request for a pattern of length m into the pattern and its extension plan
inline void decode_request(const std::string& request, uint64_t m, std::string& pat, ext_plan_t& plan)
{
    pat.assign(request.data(), m);
    std::memcpy(&plan.start, request.data() + m, sizeof(uint64_t));
    plan.order.resize(m == 0 ? 0 : m - 1);
    for (uint64_t i = 0; i + 1 < m; i++) plan.order[i] = request[m + sizeof(uint64_t) + i] == 0 ? LEFT : RIGHT;
}

direct_ifstream index_file;

// ############################# HELPERS #############################

[[noreturn]] void help(const std::string& msg = "")
{
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "move-rb-gen-ext-queries: generate the query (pattern) sets for move-rb-bench-ext, auto-calibrating" << std::endl;
    std::cout << "                 the number of patterns per file (against a move_rb index) so that running that" << std::endl;
    std::cout << "                 file's locate operation takes ~<time> seconds (see --time)." << std::endl << std::endl;
    std::cout << "                 For every pattern length m it writes two files, one per operation:" << std::endl;
    std::cout << "                 <text_name>.patterns-count-m<m> and <text_name>.patterns-locate-m<m>" << std::endl << std::endl;
    std::cout << "usage: move-rb-gen-ext-queries [options] <input_file> <move_rb_index>" << std::endl;
    std::cout << "   -o <dir>          output directory for the pattern files (default: current directory)" << std::endl;
    std::cout << "   -n <text_name>    text name used in the file names (default: basename of <input_file>)" << std::endl;
    std::cout << "   -x <chars>        characters a generated pattern must not contain (e.g. a FASTA separator or" << std::endl;
    std::cout << "                     the terminator); such patterns are skipped (default: none)" << std::endl;
    std::cout << "   -m <list>         comma-separated list of pattern lengths (default: "
              << join_csv(default_pattern_lengths) << ")" << std::endl;
    std::cout << "   --time <s>        target query duration per pattern file, in seconds (default: "
              << default_target_seconds << "s)" << std::endl;
    std::cout << "   --only <op>       restrict to count or locate files: count, locate or both (default: both)" << std::endl;
    std::cout << "   --min <N>         minimum number of patterns per file (default: " << default_min_patterns << ")"
              << std::endl;
    std::cout << "   --timeout <s>     per-pattern query timeout in seconds; a query runs in a forked worker that is" << std::endl;
    std::cout << "                     killed if it exceeds the timeout, then the pattern is discarded (0 = off, default: "
              << default_timeout_seconds << ")" << std::endl;
    std::cout << "   --max-skips <N>   give up on a set after this many consecutive timeouts (default: "
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

    direct_ifstream tf(path_input_file, std::ios::binary);
    tf.seekg(0, std::ios::end);
    uint64_t text_size = tf.tellg();

    const uint64_t target_ns = (uint64_t) (target_seconds * 1e9);
    const uint64_t timeout_ns = (uint64_t) (timeout_seconds * 1e9);

    std::vector<ext_op_t> ops;
    if (gen_count)  ops.push_back({"count",  false});
    if (gen_locate) ops.push_back({"locate", true });

    std::cout << "calibrating to ~" << target_seconds << "s per file\n" << std::endl;

    std::mt19937_64 pos_rng(std::random_device{}());

    for (const ext_op_t& op : ops) {
        for (uint64_t m : ms) {
            if (m < 1) continue;
            if (m + 1 >= text_size) {
                std::cerr << "skipping " << op.op << " m" << m << ": pattern length >= text size" << std::endl;
                continue;
            }

            std::string pattern;
            resize_no_init(pattern, m);
            std::string buffer; // the concatenated patterns to write
            std::uniform_int_distribution<uint64_t> pos_distrib(0, text_size - m - 1);

            // the generator of the extension plans (move-rb-bench-ext draws the same plans)
            std::mt19937_64 plan_rng(ext_plan_seed);

            uint64_t n = 0, total_ns = 0, total_occ = 0, skipped = 0, consecutive_skips = 0;

            // the query (count or locate) of a pattern with its extension plan
            auto run_query = [index, is_locate = op.is_locate](const std::string& pat, const ext_plan_t& plan) -> uint64_t {
                if (is_locate) {
                    uint64_t occ = 0;
                    ext_locate<idx_t>(*index, pat, plan, [&](uint64_t){ occ++; });
                    return occ;
                }
                return ext_count<idx_t>(*index, pat, plan);
            };

            // runs the query of a pattern and measures its time
            auto timed_query = [&run_query](const std::string& pat, const ext_plan_t& plan) {
                auto t1 = now();
                uint64_t c = run_query(pat, plan);
                return timed_result_t { c, time_diff_ns(t1) };
            };

            // the query of a request in the child of timeout_runner (the request is decoded before the measurement)
            auto timed_request = [&timed_query, m, pat = std::string(), plan = ext_plan_t()](
                const std::string& request) mutable {
                decode_request(request, m, pat, plan);
                return timed_query(pat, plan);
            };

            std::optional<timeout_runner<decltype(timed_request)>> runner;
            if (timeout_ns != 0) runner.emplace(timed_request, request_size(m));
            std::string request;

            // add patterns until their queries take target_ns in total (but at least min_patterns patterns)
            while (n < min_patterns || total_ns < target_ns) {
                tf.seekg(pos_distrib(pos_rng), std::ios::beg);
                tf.read(pattern.data(), m);

                // skip a pattern that contains a forbidden symbol (before its plan is drawn)
                if (!forbidden.empty() && pattern.find_first_of(forbidden) != std::string::npos) continue;

                // draw the plan (the generator is restored if the pattern is discarded)
                const std::mt19937_64 plan_rng_saved = plan_rng;
                const ext_plan_t plan = make_ext_plan(m, plan_rng);

                std::optional<timed_result_t> r;

                if (runner) {
                    encode_request(pattern, plan, request);
                    r = runner->run(request, timeout_ns);
                } else {
                    r = timed_query(pattern, plan);
                }

                // the query timed out: discard the pattern; give up after max_consecutive_skips timeouts in a row
                if (!r.has_value()) {
                    plan_rng = plan_rng_saved;
                    skipped++;
                    if (++consecutive_skips > max_consecutive_skips) {
                        std::cerr << "giving up on " << op.op << " m" << m << " after " << consecutive_skips
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
                std::cerr << "warning: not enough patterns for " << op.op << " m" << m
                          << " completed within the timeout; skipping" << std::endl;
                continue;
            }

            std::string fname = (std::filesystem::path(out_dir) /
                (text_name + ".patterns-" + op.op + "-m" + std::to_string(m))).string();
            write_pattern_file(fname, buffer, n, m, text_name, forbidden);

            std::cout << op.op << " m" << m << ": " << n << " patterns";
            if (skipped > 0) std::cout << " (" << skipped << " skipped)";
            std::cout << ", " << (total_ns / 1e9) << "s";
            if (n > 0) std::cout << " (" << (total_ns / n) / 1e3 << " us/pattern, "
                                 << (total_occ / n) << " occ/pattern)";
            std::cout << std::endl;
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
        else if (o == "-x") forbidden = need("-x");
        else if (o == "--min")  min_patterns = std::stoull(need("--min"));
        else if (o == "--time")      { if (!parse_double(need("--time"), target_seconds)) help("error: --time expects a number"); }
        else if (o == "--timeout")   { if (!parse_double(need("--timeout"), timeout_seconds)) help("error: --timeout expects a number"); }
        else if (o == "--max-skips") max_consecutive_skips = std::stoull(need("--max-skips"));
        else if (o == "-m") { ms.clear(); for (auto& v : split_csv(need("-m"))) ms.push_back(std::stoull(v)); }
        else if (o == "--only") {
            std::string ops = need("--only");
            if      (ops == "count")  { gen_count = true;  gen_locate = false; }
            else if (ops == "locate") { gen_count = false; gen_locate = true; }
            else if (ops == "both")   { gen_count = true;  gen_locate = true; }
            else help("error: --only must be count, locate or both");
        }
        else {
            help("error: unrecognized option '" + o + "'");
        }
    }

    if (argc - i != 2) help();
    path_input_file = argv[i++];
    path_index_file = argv[i++];

    if (target_seconds <= 0) help("error: --time must be > 0");
    if (timeout_seconds < 0) help("error: --timeout must be >= 0");
    if (min_patterns < 1) min_patterns = 1;
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
