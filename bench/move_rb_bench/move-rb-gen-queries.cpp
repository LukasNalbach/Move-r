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

// move-rb-gen-queries: generates the query (pattern) sets for move-rb-bench. For every
// (k, m) combination it writes three pattern files - one per benchmark operation:
//   <text_name>.patterns-count-hamming-k<k>-m<m>
//   <text_name>.patterns-locate-hamming-k<k>-m<m>
//   <text_name>.patterns-locate-edit-k<k>-m<m>
// The number of patterns in each file is auto-calibrated (using a move_rb index)
// so that running that file's operation takes approximately a target duration.

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <clocale>
#include <csignal>
#include <cstdio>
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

#include <move_rb/move_rb.hpp>

#include <ips4o.hpp>

// ############################# CONFIGURATION #############################

std::string path_input_file;         // the text the index was built for
std::string path_index_file;         // a move_rb index for that text
double      target_seconds = 1.0;    // target count/locate duration per pattern file (set via --time)

std::string out_dir = ".";           // where the pattern files are written
std::string text_name;               // name used in the file names (default: basename of input)
std::string scheme_str = "min_u";    // search scheme (k taken from the grid)
bool gen_count = true;               // generate the count-pattern files
bool gen_locate = true;              // generate the locate-pattern files
uint64_t min_patterns = 1;           // at least this many patterns per file
uint64_t seed = 0;                   // RNG seed (0 => random)
double   timeout_seconds = 10.0;     // per-pattern query timeout (0 => disabled)
uint64_t max_consecutive_skips = 10; // give up on a (k, m) after this many timeouts in a row

std::vector<int64_t>  ks = {4, 7, 10, 13};
std::vector<uint64_t> ms = {10, 20, 40, 80, 160, 320, 640, 1280};

std::ifstream index_file;

// the four benchmark operations a pattern file can target
struct query_type_t { const char* op; const char* metric; bool is_count; bool is_edit; };

// ############################# HELPERS #############################

void help(std::string msg)
{
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "move-rb-gen-queries: generate the query (pattern) sets for move-rb-bench, auto-calibrating" << std::endl;
    std::cout << "                 the number of patterns per file (against a move_rb index) so that" << std::endl;
    std::cout << "                 running that file's operation takes ~<time> seconds (see --time)." << std::endl << std::endl;
    std::cout << "                 For every (k, m) it writes three files, one per benchmark operation:" << std::endl;
    std::cout << "                 <text_name>.patterns-count-hamming, -locate-hamming and -locate-edit (-k<k>-m<m>)" << std::endl << std::endl;
    std::cout << "usage: move-rb-gen-queries [options] <input_file> <move_rb_index>" << std::endl;
    std::cout << "   -o <dir>          output directory for the pattern files (default: current directory)" << std::endl;
    std::cout << "   -n <text_name>    text name used in the file names (default: basename of <input_file>)" << std::endl;
    std::cout << "   -k <list>         comma-separated list of error counts k (default: 4,7,10,13)" << std::endl;
    std::cout << "   -m <list>         comma-separated list of pattern lengths (default: 10,20,40,80,160,320,640,1280)" << std::endl;
    std::cout << "   -s <scheme>       search scheme: pigeon_hole, suffix_filter, min_u or 01 (default: min_u)" << std::endl;
    std::cout << "   --time <s>        target count/locate duration per pattern file, in seconds (default: 1s)" << std::endl;
    std::cout << "   --only <ops>      restrict to count or locate files: count, locate or count,locate (default: both)" << std::endl;
    std::cout << "   --min <N>         minimum number of patterns per file (default: 1)" << std::endl;
    std::cout << "   --seed <s>        seed for the random pattern positions (default: random)" << std::endl;
    std::cout << "   --timeout <s>     per-pattern query timeout in seconds; a query runs in a forked worker that is" << std::endl;
    std::cout << "                     killed if it exceeds the timeout, then the pattern is discarded (0 = off, default: 10)" << std::endl;
    std::cout << "   --max-skips <N>   give up on a (k,m) after this many consecutive timeouts (default: 10)" << std::endl;
    std::cout << "   <input_file>      the text the index was built for" << std::endl;
    std::cout << "   <move_rb_index>   a move_rb index (.move-rb or .move-rb-rlzsa) of <input_file>" << std::endl;
    exit(0);
}

std::vector<std::string> split_csv(const std::string& s)
{
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) if (!item.empty()) out.push_back(item);
    return out;
}

search_scheme_t make_scheme(int64_t k)
{
    if (scheme_str == "suffix_filter") return suffix_filter_scheme(k);
    if (scheme_str == "min_u")         return min_u_scheme(k);
    if (scheme_str == "01")            return zero_one_scheme(k);
    return min_u_scheme(k); // default: min_u
}

// the outcome of a timed query: its result count and the in-process time the query itself took
struct timed_result_t { uint64_t count; uint64_t elapsed_ns; };

// reads exactly n bytes from fd into buf; false on EOF or error (e.g. the peer closed its pipe end)
static bool read_full(int fd, void* buf, size_t n)
{
    char* p = static_cast<char*>(buf);
    while (n > 0) {
        ssize_t r = read(fd, p, n);
        if (r > 0) { p += r; n -= (size_t) r; }
        else if (r == 0) return false;     // EOF: the peer is gone
        else if (errno == EINTR) continue; // interrupted: retry
        else return false;                 // error
    }
    return true;
}

// writes exactly n bytes from buf to fd; false if the peer is gone (EPIPE) or on error
static bool write_full(int fd, const void* buf, size_t n)
{
    const char* p = static_cast<const char*>(buf);
    while (n > 0) {
        ssize_t w = write(fd, p, n);
        if (w > 0) { p += w; n -= (size_t) w; }
        else if (w < 0 && errno == EINTR) continue;
        else return false;
    }
    return true;
}

// waits up to timeout_ns for fd to become readable; false on timeout or error
static bool wait_readable(int fd, uint64_t timeout_ns)
{
    struct pollfd pfd { fd, POLLIN, 0 };
    struct timespec ts { (time_t) (timeout_ns / 1000000000ull), (long) (timeout_ns % 1000000000ull) };
    int r;
    do { r = ppoll(&pfd, 1, &ts, nullptr); } while (r < 0 && errno == EINTR);
    return r > 0 && (pfd.revents & POLLIN);
}

/**
 * @brief runs queries on a forked worker PROCESS so each one can be hard-timed-out.
 * @tparam query_fnc_t query function: (const std::string& pattern) -> uint64_t occurrence count
 */
template <typename query_fnc_t>
class timeout_runner {
    query_fnc_t query;
    uint64_t m;          // bytes per pattern (the child reads exactly this many per query)
    pid_t child = -1;
    int to_child = -1;   // parent -> child: the next pattern
    int from_child = -1; // child -> parent: the timed_result_t

    void spawn()
    {
        int down[2], up[2]; // down: parent->child, up: child->parent
        if (pipe(down) != 0 || pipe(up) != 0) { perror("pipe"); std::exit(1); }
        pid_t pid = fork();
        if (pid < 0) { perror("fork"); std::exit(1); }

        if (pid == 0) {
            // child: inherits the loaded index (COW). Runs queries until the parent closes the pipe (or
            // it is killed on timeout), timing each one itself so the measurement excludes the IPC.
            close(down[1]); close(up[0]);
            std::string pat;
            no_init_resize(pat, m);
            while (read_full(down[0], pat.data(), m)) {
                auto t1 = now();
                uint64_t c = query(pat);
                timed_result_t r { c, time_diff_ns(t1) };
                if (!write_full(up[1], &r, sizeof(r))) break;
            }
            _exit(0);
        }

        close(down[0]); close(up[1]);
        child = pid; to_child = down[1]; from_child = up[0];
    }

    // SIGKILLs the worker (freeing all its memory at once) and reaps it
    void reap()
    {
        if (child <= 0) return;
        kill(child, SIGKILL);
        int st;
        while (waitpid(child, &st, 0) < 0 && errno == EINTR) {}
        close(to_child); close(from_child);
        child = -1; to_child = from_child = -1;
    }

public:
    timeout_runner(query_fnc_t q, uint64_t pattern_len) : query(std::move(q)), m(pattern_len) { spawn(); }
    ~timeout_runner() { reap(); }

    // runs one query on the worker, waiting up to timeout_ns; nullopt on timeout (worker killed + replaced)
    std::optional<timed_result_t> run(const std::string& pattern, uint64_t timeout_ns)
    {
        if (write_full(to_child, pattern.data(), m) && wait_readable(from_child, timeout_ns)) {
            timed_result_t r;
            if (read_full(from_child, &r, sizeof(r))) return r;
        }

        reap();
        spawn();
        return std::nullopt;
    }
};

// ############################# GENERATION #############################

/**
 * @brief generates and writes all pattern files for the loaded index, calibrating the number
 *        of patterns in each so that its operation takes ~target_seconds
 * @tparam pos_t the index integer type
 * @tparam support the move_rb support (_locate_move or _locate_rlzsa)
 */
template <typename pos_t, move_r_support support>
void generate()
{
    using idx_t = move_rb<support, char, pos_t>;
    std::cout << "loading the move_rb index ..." << std::endl;
    auto index = std::make_shared<idx_t>();
    index->load(index_file);
    index_file.close();

    std::ifstream tf(path_input_file, std::ios::binary);
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
            search_scheme_t scheme = make_scheme(k);

            for (uint64_t m : ms) {
                if (m < scheme.p) {
                    std::cerr << "skipping " << type.op << "-" << type.metric << " k" << k << " m" << m
                              << ": pattern length < " << (int) scheme.p << " parts in the search scheme" << std::endl;
                    continue;
                }

                std::string pattern;
                no_init_resize(pattern, m);
                std::string buffer; // the concatenated patterns to write
                std::uniform_int_distribution<uint64_t> pos_distrib(0, text_size - m - 1);
                uint64_t n = 0, total_ns = 0, total_occ = 0, skipped = 0, consecutive_skips = 0;

                auto run_query = [index, scheme, type](const std::string& pat) -> uint64_t {
                    if (type.is_count) {
                        return index->count_hamming_dist(pat, scheme);
                    }
                    std::vector<aprx_occ_t<pos_t>> occ;
                    auto report = [&](aprx_occ_t<pos_t> o) { occ.emplace_back(o); };
                    if (type.is_edit) {
                        index->template locate<EDIT_DISTANCE>(pat, scheme, report);
                        ips4o::sort(occ.begin(), occ.end());
                        filter_aprx_occurrences<pos_t>(occ, (pos_t) scheme.k);
                    } else {
                        index->template locate<HAMMING_DISTANCE>(pat, scheme, report);
                    }
                    return occ.size();
                };

                std::optional<timeout_runner<decltype(run_query)>> runner;
                if (timeout_ns != 0) runner.emplace(run_query, m);

                // add patterns until the accumulated query time reaches the target
                // (a single very expensive pattern still yields one pattern)
                while (n < min_patterns || total_ns < target_ns) {
                    tf.seekg(pos_distrib(rng), std::ios::beg);
                    tf.read(pattern.data(), m);

                    std::optional<timed_result_t> r;
                    if (runner) {
                        r = runner->run(pattern, timeout_ns);
                    } else {
                        auto t1 = now();
                        uint64_t c = run_query(pattern);
                        r = timed_result_t { c, time_diff_ns(t1) };
                    }

                    // the query exceeded the timeout: discard this pattern and try another, unless
                    // too many in a row time out (then this (k, m) is just too slow to benchmark)
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
                std::ofstream of(fname, std::ios::binary);
                of << "# number=" << n << " length=" << m << " file=" << text_name << " forbidden=\n";
                of.write(buffer.data(), buffer.size());
                of.close();

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

    std::cout << "\ndone; pattern files written to " << out_dir << std::endl;
}

// ############################# DRIVER #############################

int main(int argc, char** argv)
{
    std::setlocale(LC_ALL, "C");
    std::signal(SIGPIPE, SIG_IGN);

    int i = 1;
    while (i < argc && argv[i][0] == '-') {
        std::string o = argv[i++];
        auto need = [&](const char* opt) { if (i >= argc) help(std::string("error: missing argument after ") + opt); return std::string(argv[i++]); };
        if      (o == "-o") out_dir = need("-o");
        else if (o == "-n") text_name = need("-n");
        else if (o == "-s") scheme_str = need("-s");
        else if (o == "--min")  min_patterns = std::stoull(need("--min"));
        else if (o == "--seed") seed = std::stoull(need("--seed"));
        else if (o == "--time")      target_seconds = std::stod(need("--time"));
        else if (o == "--timeout")   timeout_seconds = std::stod(need("--timeout"));
        else if (o == "--max-skips") max_consecutive_skips = std::stoull(need("--max-skips"));
        else if (o == "-k") { ks.clear(); for (auto& v : split_csv(need("-k"))) ks.push_back(std::stoll(v)); }
        else if (o == "-m") { ms.clear(); for (auto& v : split_csv(need("-m"))) ms.push_back(std::stoull(v)); }
        else if (o == "--only") {
            std::string ops = need("--only");
            gen_count = ops.find("count") != std::string::npos;
            gen_locate = ops.find("locate") != std::string::npos;
            if (!gen_count && !gen_locate) help("error: --only must contain count and/or locate");
        } else {
            help("error: unrecognized option '" + o + "'");
        }
    }

    if (argc - i != 2) help("");
    path_input_file = argv[i++];
    path_index_file = argv[i++];

    if (target_seconds <= 0) help("error: --time must be > 0");
    if (timeout_seconds < 0) help("error: --timeout must be >= 0");
    if (min_patterns < 1) min_patterns = 1;
    if (scheme_str != "pigeon_hole" && scheme_str != "suffix_filter" &&
        scheme_str != "min_u" && scheme_str != "01")
        help("error: -s must be pigeon_hole, suffix_filter, min_u or 01");
    if (!std::filesystem::exists(path_input_file)) help("error: <input_file> does not exist");
    if (text_name.empty()) text_name = std::filesystem::path(path_input_file).filename().string();
    if (!std::filesystem::is_directory(out_dir)) help("error: output directory '" + out_dir + "' does not exist");

    index_file.open(path_index_file, std::ios::binary);
    if (!index_file.good()) help("error: cannot open <move_rb_index>");

    bool is_64_bit;
    move_r_support support;
    index_file.read((char*) &is_64_bit, 1);
    index_file.read((char*) &support, sizeof(move_r_support));
    index_file.seekg(0, std::ios::beg);

    if (support == _count) {
        help("error: this index only supports count, not locate; rebuild with locate support");
    } else if (support == _locate_move) {
        if (is_64_bit) generate<uint64_t, _locate_move>();
        else           generate<uint32_t, _locate_move>();
    } else if (support == _locate_rlzsa) {
        if (is_64_bit) generate<uint64_t, _locate_rlzsa>();
        else           generate<uint32_t, _locate_rlzsa>();
    } else {
        help("error: unknown move_rb support in the index file");
    }

    return 0;
}
