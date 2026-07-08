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

// move-rb-gen-ext-queries: generates the pattern sets for move-rb-bench-ext. For every length m it writes a count
// and a locate file (<text_name>.patterns-{count,locate}-m<m>) of sampled text substrings, each with its pattern
// count calibrated against a move_rb index to run for ~<time> seconds. The extension plans are not stored: both
// the position and plan RNGs are reset to a fixed seed per pattern set (see ext_query.hpp).

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <clocale>
#include <csignal>
#include <cstdio>
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

#include "ext_query.hpp" // ext_plan_t / make_ext_plan / ext_locate (the shared query model)
#include <move_rb/move_rb.hpp>

// ############################# CONFIGURATION #############################

std::string path_input_file;         // the text the index was built for
std::string path_index_file;         // a move_rb index for that text
double      target_seconds = 1.0;    // target locate duration per pattern file (set via --time)

std::string out_dir = ".";           // where the pattern files are written
std::string text_name;               // name used in the file names (default: basename of input)
bool gen_count = true;               // generate the count-pattern files
bool gen_locate = true;              // generate the locate-pattern files
std::string forbidden;               // characters a generated pattern must not contain (e.g. a FASTA separator)
uint64_t min_patterns = 1;           // at least this many patterns per file
uint64_t seed = 0;                   // position-sampling RNG seed (0 => random); reset per pattern set, not recorded
double   timeout_seconds = 10.0;     // per-pattern query timeout (0 => disabled)
uint64_t max_consecutive_skips = 10; // give up on a set after this many timeouts in a row

std::vector<uint64_t> ms = {10, 20, 40, 80, 160, 320, 640, 1280};

// the two benchmark operations a pattern file targets; count and locate get separate files and calibrations
struct ext_op_t { const char* op; bool is_locate; };

std::ifstream index_file;

// ############################# HELPERS #############################

void help(std::string msg)
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
    std::cout << "   -m <list>         comma-separated list of pattern lengths (default: 10,20,40,80,160,320,640,1280)" << std::endl;
    std::cout << "   --time <s>        target query duration per pattern file, in seconds (default: 1s)" << std::endl;
    std::cout << "   --only <op>       restrict to count or locate files: count, locate or both (default: both)" << std::endl;
    std::cout << "   --min <N>         minimum number of patterns per file (default: 1)" << std::endl;
    std::cout << "   --seed <s>        seed for the random pattern positions (default: random); reset per pattern set so" << std::endl;
    std::cout << "                     generation is reproducible. The extension plans use a fixed seed and are not recorded." << std::endl;
    std::cout << "   --timeout <s>     per-pattern query timeout in seconds; a query runs in a forked worker that is" << std::endl;
    std::cout << "                     killed if it exceeds the timeout, then the pattern is discarded (0 = off, default: 10)" << std::endl;
    std::cout << "   --max-skips <N>   give up on a set after this many consecutive timeouts (default: 10)" << std::endl;
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
 * @brief runs extension-locate queries on a forked worker PROCESS so each one can be hard-timed-out. Each query
 *        request is a fixed-size byte blob (per pattern of length m): the m pattern bytes, then the 8-byte start
 *        position, then the m-1 order bytes (0 = LEFT, 1 = RIGHT).
 * @tparam query_fnc_t query function: (const std::string& pattern, const ext_plan_t& plan) -> uint64_t occ count
 */
template <typename query_fnc_t>
class timeout_runner {
    query_fnc_t query;
    uint64_t m;                 // pattern length
    uint64_t req_size;          // bytes per request (m + 8 + (m-1))
    pid_t child = -1;
    int to_child = -1;   // parent -> child: the next request
    int from_child = -1; // child -> parent: the timed_result_t

    // decodes a request blob into (pattern, plan)
    static void decode(const std::vector<char>& req, uint64_t m, std::string& pat, ext_plan_t& plan)
    {
        pat.assign(req.data(), m);
        std::memcpy(&plan.start, req.data() + m, sizeof(uint64_t));
        plan.order.resize(m == 0 ? 0 : m - 1);
        for (uint64_t i = 0; i + 1 < m; i++) plan.order[i] = req[m + 8 + i] == 0 ? LEFT : RIGHT;
    }

    void spawn()
    {
        int down[2], up[2]; // down: parent->child, up: child->parent
        if (pipe(down) != 0 || pipe(up) != 0) { perror("pipe"); std::exit(1); }
        pid_t pid = fork();
        if (pid < 0) { perror("fork"); std::exit(1); }

        if (pid == 0) {
            // child: inherits the loaded index (COW). Runs queries until the parent closes the pipe (or it is
            // killed on timeout), timing each one itself so the measurement excludes the IPC.
            close(down[1]); close(up[0]);
            std::vector<char> req(req_size);
            std::string pat;
            ext_plan_t plan;
            while (read_full(down[0], req.data(), req_size)) {
                decode(req, m, pat, plan);
                auto t1 = now();
                uint64_t c = query(pat, plan);
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

    // encodes (pattern, plan) into a request blob
    void encode(const std::string& pat, const ext_plan_t& plan, std::vector<char>& req) const
    {
        req.resize(req_size);
        std::memcpy(req.data(), pat.data(), m);
        std::memcpy(req.data() + m, &plan.start, sizeof(uint64_t));
        for (uint64_t i = 0; i + 1 < m; i++) req[m + 8 + i] = plan.order[i] == LEFT ? 0 : 1;
    }

public:
    timeout_runner(query_fnc_t q, uint64_t pattern_len)
        : query(std::move(q)), m(pattern_len), req_size(pattern_len + 8 + (pattern_len == 0 ? 0 : pattern_len - 1))
    { spawn(); }
    ~timeout_runner() { reap(); }

    // runs one query on the worker, waiting up to timeout_ns; nullopt on timeout (worker killed + replaced)
    std::optional<timed_result_t> run(const std::string& pattern, const ext_plan_t& plan, uint64_t timeout_ns)
    {
        std::vector<char> req;
        encode(pattern, plan, req);
        if (write_full(to_child, req.data(), req_size) && wait_readable(from_child, timeout_ns)) {
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
 * @brief generates and writes all pattern files for the loaded index, calibrating the number of patterns in each
 *        so that its locate operation takes ~target_seconds
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

    std::vector<ext_op_t> ops;
    if (gen_count)  ops.push_back({"count",  false});
    if (gen_locate) ops.push_back({"locate", true });

    std::cout << "calibrating to ~" << target_seconds << "s per file (position seed " << seed << ")\n" << std::endl;

    for (const ext_op_t& op : ops) {
        for (uint64_t m : ms) {
            if (m < 1) continue;
            if (m + 1 >= text_size) {
                std::cerr << "skipping " << op.op << " m" << m << ": pattern length >= text size" << std::endl;
                continue;
            }

            std::string pattern;
            no_init_resize(pattern, m);
            std::string buffer; // the concatenated patterns to write
            std::uniform_int_distribution<uint64_t> pos_distrib(0, text_size - m - 1);

            // both RNGs are reset per pattern set: the position RNG (--seed) so generation is reproducible, the
            // plan RNG (ext_plan_seed) so move-rb-bench-ext replays the same plans with nothing stored in the file
            std::mt19937_64 pos_rng(seed);
            std::mt19937_64 plan_rng(ext_plan_seed);

            uint64_t n = 0, total_ns = 0, total_occ = 0, skipped = 0, consecutive_skips = 0;

            // the ext query (count or locate) for one pattern and its extension plan
            auto run_query = [index, is_locate = op.is_locate](const std::string& pat, const ext_plan_t& plan) -> uint64_t {
                if (is_locate) {
                    uint64_t occ = 0;
                    ext_locate<idx_t>(*index, pat, plan, [&](pos_t){ occ++; });
                    return occ;
                }
                return ext_count<idx_t>(*index, pat, plan);
            };

            std::optional<timeout_runner<decltype(run_query)>> runner;
            if (timeout_ns != 0) runner.emplace(run_query, m);

            // add patterns until the accumulated query time reaches the target (a single very expensive pattern
            // still yields one pattern)
            while (n < min_patterns || total_ns < target_ns) {
                tf.seekg(pos_distrib(pos_rng), std::ios::beg);
                tf.read(pattern.data(), m);

                // skip a forbidden pattern before drawing a plan, so the plan RNG stays aligned with accepted ones
                if (!forbidden.empty() && pattern.find_first_of(forbidden) != std::string::npos) continue;

                // draw the plan, saving the RNG so it can be restored if the pattern is discarded
                const std::mt19937_64 plan_rng_saved = plan_rng;
                const ext_plan_t plan = make_ext_plan(m, plan_rng);

                std::optional<timed_result_t> r;
                if (runner) {
                    r = runner->run(pattern, plan, timeout_ns);
                } else {
                    auto t1 = now();
                    uint64_t c = run_query(pattern, plan);
                    r = timed_result_t { c, time_diff_ns(t1) };
                }

                // the query exceeded the timeout: discard this pattern and try another, unless too many in a row
                // time out (then this set is just too slow to benchmark)
                if (!r.has_value()) {
                    plan_rng = plan_rng_saved; // restore: a discarded pattern must not consume a plan
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
            std::ofstream of(fname, std::ios::binary);
            of << "# number=" << n << " length=" << m << " file=" << text_name << " forbidden=" << forbidden << "\n";
            of.write(buffer.data(), buffer.size());
            of.close();

            std::cout << op.op << " m" << m << ": " << n << " patterns";
            if (skipped > 0) std::cout << " (" << skipped << " skipped)";
            std::cout << ", " << (total_ns / 1e9) << "s";
            if (n > 0) std::cout << " (" << (total_ns / n) / 1e3 << " us/pattern, "
                                 << (total_occ / n) << " occ/pattern)";
            std::cout << std::endl;
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
        else if (o == "-x") forbidden = need("-x");
        else if (o == "--min")  min_patterns = std::stoull(need("--min"));
        else if (o == "--seed") seed = std::stoull(need("--seed"));
        else if (o == "--time")      target_seconds = std::stod(need("--time"));
        else if (o == "--timeout")   timeout_seconds = std::stod(need("--timeout"));
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

    if (argc - i != 2) help("");
    path_input_file = argv[i++];
    path_index_file = argv[i++];

    if (target_seconds <= 0) help("error: --time must be > 0");
    if (timeout_seconds < 0) help("error: --timeout must be >= 0");
    if (min_patterns < 1) min_patterns = 1;
    if (!std::filesystem::exists(path_input_file)) help("error: <input_file> does not exist");
    if (text_name.empty()) text_name = std::filesystem::path(path_input_file).filename().string();
    if (!std::filesystem::is_directory(out_dir)) help("error: output directory '" + out_dir + "' does not exist");

    // resolve a concrete position seed (reset per pattern set); the plans use the fixed ext_plan_seed instead
    if (seed == 0) seed = std::random_device{}();

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
