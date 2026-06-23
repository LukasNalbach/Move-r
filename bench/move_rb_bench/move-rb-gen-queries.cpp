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
#include <chrono>
#include <clocale>
#include <condition_variable>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <move_rb/move_rb.hpp>

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
    std::cout << "   --timeout <s>     per-pattern query timeout in seconds; a query is run in a thread and, if it" << std::endl;
    std::cout << "                     exceeds the timeout, the pattern is discarded (0 = off, default: 10)" << std::endl;
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

/**
 * @brief runs queries on a worker thread
 * @tparam query_fnc_t query function
 */
template <typename query_fnc_t>
class timeout_runner {
    // shared between the main thread and the (possibly abandoned) worker; held via shared_ptr so it
    // outlives an abandoned worker
    struct slot_t {
        std::mutex m;
        std::condition_variable cv;
        std::string pattern;
        bool go = false;   // a request is pending
        bool done = false; // the result is ready
        bool quit = false; // the worker should exit (shutdown or abandonment)
        timed_result_t result {};
    };

    query_fnc_t query;
    std::shared_ptr<slot_t> s;

    void spawn()
    {
        s = std::make_shared<slot_t>();
        std::thread([sl = s, q = query]() {
            std::unique_lock<std::mutex> lk(sl->m);
            while (true) {
                sl->cv.wait(lk, [&] { return sl->go || sl->quit; });
                if (sl->quit) return;
                std::string p = std::move(sl->pattern);
                sl->go = false;
                lk.unlock();

                auto t1 = now();
                uint64_t c = q(p);
                uint64_t elapsed = time_diff_ns(t1);

                lk.lock();
                sl->result = { c, elapsed };
                sl->done = true;
                bool abandoned = sl->quit; // abandoned while computing -> exit after signaling
                sl->cv.notify_one();
                if (abandoned) return;
            }
        }).detach();
    }

public:
    explicit timeout_runner(query_fnc_t q) : query(std::move(q)) { spawn(); }

    // tells the current worker to exit once finished (the detached thread keeps its slot alive)
    ~timeout_runner() { std::lock_guard<std::mutex> lk(s->m); s->quit = true; s->cv.notify_one(); }

    // runs the query on the worker, waiting up to timeout_ns; nullopt on timeout (worker abandoned)
    std::optional<timed_result_t> run(const std::string& pattern, uint64_t timeout_ns)
    {
        {
            std::lock_guard<std::mutex> lk(s->m);
            s->pattern = pattern;
            s->go = true;
            s->done = false;
        }
        s->cv.notify_one();

        std::unique_lock<std::mutex> lk(s->m);
        if (s->cv.wait_for(lk, std::chrono::nanoseconds(timeout_ns), [&] { return s->done; }))
            return s->result;

        // timeout: abandon this worker (it exits after finishing its slow query) and start a fresh one
        s->quit = true;
        lk.unlock();
        s->cv.notify_one();
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
    // shared so that a query thread abandoned on timeout can keep reading the index after generate()
    // returns: each query closure captures this shared_ptr, so the index stays alive until the last
    // thread (the main thread or an abandoned worker) releases it, then it is freed - no leak
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

                // runs one query and returns its result count. Captures the index (shared_ptr) and
                // scheme/type/k by value, and uses its OWN occ vector, so it carries no shared mutable
                // state and keeps the index alive - that makes it safe to run on a detached timeout
                // thread (and concurrently with an abandoned one). Timing is done by the caller.
                auto run_query = [index, scheme, type, k](const std::string& pat) -> uint64_t {
                    if (type.is_count) {
                        return index->count_hamming_dist(pat, scheme);
                    }
                    std::vector<aprx_occ_t<pos_t>> occ;
                    if (type.is_edit) {
                        index->template locate<EDIT_DISTANCE>(pat, scheme,
                            [&](aprx_occ_t<pos_t> o) { occ.emplace_back(o); });
                        ips4o::sort(occ.begin(), occ.end());
                        filter_aprx_occurrences<pos_t>(occ, (pos_t) k);
                    } else {
                        index->template locate<HAMMING_DISTANCE>(pat, scheme,
                            [&](aprx_occ_t<pos_t> o) { occ.emplace_back(o); });
                    }
                    return occ.size();
                };

                // with a timeout, queries run on a persistent worker thread that enforces it;
                // without one they run in-process. The runner is created here so it exists only while
                // this (operation, k, m) is being calibrated.
                std::optional<timeout_runner<decltype(run_query)>> runner;
                if (timeout_ns != 0) runner.emplace(run_query);

                // add patterns until the accumulated query time reaches the target
                // (a single very expensive pattern still yields one pattern)
                while (n < min_patterns || total_ns < target_ns) {
                    tf.seekg(pos_distrib(rng), std::ios::beg);
                    tf.read(pattern.data(), m);

                    std::optional<timed_result_t> r;
                    if (runner) {
                        r = runner->run(pattern, timeout_ns); // hard timeout on a worker thread
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
