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

#include "indexes/b_move_adapter.hpp"
#include "indexes/br_index_adapter.hpp"
#include <move_rb/move_rb.hpp>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <ips4o.hpp>

// ############################# CONFIGURATION #############################

// the maximum number of errors the edit-distance algorithm supports
static constexpr int EDIT_K_LIMIT = 20;

// the operation a pattern file is to be measured with
enum bench_op_t : uint8_t {
    OP_COUNT = 0, // approximate count
    OP_LOCATE = 1 // approximate locate
};

// a single pattern file to be measured; the (op, metric) pair is one of the four
// benchmark operations (count/locate x hamming/edit) encoded in the file name
struct bench_job_t {
    bench_op_t op; // count or locate
    distance_metric_t metric; // HAMMING_DISTANCE or EDIT_DISTANCE
    int64_t k; // maximum number of errors (from the file name)
    uint64_t m; // pattern length (from the file name)
    std::string path; // path to the pattern file
};

// the configuration of one benchmark run
struct bench_config_t {
    std::string bmove_dir; // directory containing the b-move index files (empty => b-move not measured)
    std::string br_index_path; // path to the br-index (.bri) file (empty => br-index not measured)
    std::string move_rb_move_path; // path to the move_rb index using move data structures
    std::string move_rb_rlzsa_path; // path to the move_rb index using an rlzsa
    std::string text_name; // name of the original text (used in the output and pattern file names)
    std::string patterns_path; // directory containing the pattern files
    std::vector<bench_job_t> jobs; // the pattern files to measure

    // search scheme to use (like move-rb-locate's -s): the name of a built-in scheme
    // (pigeon_hole, suffix_filter, min_u, 01) or the path to a scheme file
    std::string scheme = "min_u";
    bool scheme_from_file = false; // true <=> scheme is a file; then file_scheme holds it (k is fixed by the file)
    search_scheme_t file_scheme; // the scheme parsed from the file (only used if scheme_from_file)

    // optional filters (used to drive fine-grained, externally isolated runs):
    std::string only;          // if non-empty, measure only this index (bmove|br_index|move_rb_move|move_rb_rlzsa)
    std::string metric = "all"; // which distance metric to run: all|hamming|edit
};

/**
 * @brief builds the search scheme to use for a pattern file with k errors: a built-in
 *        scheme is instantiated with k, a file scheme is returned as-is (its k is fixed)
 */
inline search_scheme_t make_scheme(const bench_config_t& cfg, int64_t k)
{
    if (cfg.scheme_from_file)           return cfg.file_scheme;
    if (cfg.scheme == "suffix_filter")  return suffix_filter_scheme(k);
    if (cfg.scheme == "min_u")          return min_u_scheme(k);
    if (cfg.scheme == "01")             return zero_one_scheme(k);
    return min_u_scheme(k); // default
}

// ############################# RESULT OUTPUT #############################

/**
 * @brief writes a single RESULT line to std::cout
 */
inline void print_result(
    const std::string& algo, const std::string& text, uint64_t n,
    const std::string& dist_metr, const std::string& search_scheme,
    uint64_t pattern_length, uint64_t num_patterns,
    int64_t max_mismatches, uint64_t num_occurrences,
    const std::string& time_key, uint64_t time_value)
{
    std::cout << "RESULT"
              << " algo=" << algo
              << " text=" << text
              << " n=" << n
              << " dist_metr=" << dist_metr
              << " search_scheme=" << search_scheme
              << " pattern_length=" << pattern_length
              << " num_patterns=" << num_patterns
              << " max_mismatches=" << max_mismatches
              << " num_occurrences=" << num_occurrences
              << " " << time_key << "=" << time_value
              << std::endl;
}

// ############################# GENERIC MEASUREMENT #############################

/**
 * @brief opens a pattern file and reads its (pizza&chili-style) header
 * @return whether the file could be opened (and the header parsed)
 */
inline bool open_patterns(const bench_job_t& job, std::ifstream& pf, uint64_t& num_patterns, uint64_t& m)
{
    pf.open(job.path, std::ios::binary);

    if (!pf.good()) {
        std::cerr << "warning: cannot open pattern file " << job.path << std::endl;
        return false;
    }

    std::string header;
    std::getline(pf, header);
    num_patterns = number_of_patterns(header);
    m = patterns_length(header);
    return true;
}

/**
 * @brief measures approximate count (w.r.t. hamming distance) of all patterns in job.path
 *        for the given index and writes a RESULT line
 * @tparam idx_t the index type (move_rb, b_move_adapter or br_index_adapter)
 */
template <typename idx_t>
void run_count(idx_t& index, uint64_t n, const bench_config_t& cfg, const bench_job_t& job,
               const search_scheme_t& scheme, const std::string& index_name)
{
    using pos_t = typename idx_t::pos_type;

    std::ifstream pf;
    uint64_t num_patterns, m;
    if (!open_patterns(job, pf, num_patterns, m)) return;

    if (m < scheme.p) {
        std::cerr << "warning: pattern length " << m << " < number of parts in the search scheme; "
                  << "skipping " << job.path << std::endl;
        return;
    }

    std::string pattern;
    no_init_resize(pattern, m);
    uint64_t num_occurrences = 0;
    uint64_t time_count = 0;

    for (uint64_t i = 0; i < num_patterns; i++) {
        pf.read(pattern.data(), m);
        auto t1 = now();
        pos_t c = index.count_hamming_dist(pattern, scheme);
        time_count += time_diff_ns(t1);
        num_occurrences += c;
    }

    print_result("count_" + index_name, cfg.text_name, n, "hamming", cfg.scheme,
        m, num_patterns, scheme.k, num_occurrences, "time_count", time_count);
}

/**
 * @brief measures approximate locate (w.r.t. dist_metr) of all patterns in job.path
 *        for the given index and writes a RESULT line
 * @tparam idx_t the index type (move_rb, b_move_adapter or br_index_adapter)
 * @tparam dist_metr the distance metric (HAMMING_DISTANCE or EDIT_DISTANCE)
 */
template <typename idx_t, distance_metric_t dist_metr>
void run_locate(idx_t& index, uint64_t n, const bench_config_t& cfg, const bench_job_t& job,
                const search_scheme_t& scheme, const std::string& index_name)
{
    using pos_t = typename idx_t::pos_type;

    std::ifstream pf;
    uint64_t num_patterns, m;
    if (!open_patterns(job, pf, num_patterns, m)) return;

    if (m < scheme.p) {
        std::cerr << "warning: pattern length " << m << " < number of parts in the search scheme; "
                  << "skipping " << job.path << std::endl;
        return;
    }

    std::string pattern;
    no_init_resize(pattern, m);
    std::vector<aprx_occ_t<pos_t>> occurrences;
    uint64_t num_occurrences = 0;
    uint64_t time_locate = 0;

    for (uint64_t i = 0; i < num_patterns; i++) {
        pf.read(pattern.data(), m);
        occurrences.clear();
        auto t1 = now();
        index.template locate<dist_metr>(pattern, scheme, [&](aprx_occ_t<pos_t> occ) { occurrences.emplace_back(occ); });

        // for edit distance the same occurrence may be reported by several search contexts;
        // deduplicate (as the move-rb-locate tool does) before counting
        if constexpr (dist_metr == EDIT_DISTANCE) {
            ips4o::sort(occurrences.begin(), occurrences.end());
            filter_aprx_occurrences<pos_t>(occurrences, (pos_t) scheme.k);
        }

        time_locate += time_diff_ns(t1);
        num_occurrences += occurrences.size();
    }

    const std::string dist_str = dist_metr == HAMMING_DISTANCE ? "hamming" : "edit";
    print_result("locate_" + index_name, cfg.text_name, n, dist_str, cfg.scheme,
        m, num_patterns, scheme.k, num_occurrences, "time_locate", time_locate);
}

/**
 * @brief measures every job for the given (already loaded) index. Each pattern file encodes
 *        exactly one of the four benchmark operations (count/locate x hamming/edit) and is
 *        measured with that operation
 * @tparam idx_t the index type (move_rb, b_move_adapter or br_index_adapter)
 */
template <typename idx_t>
void measure_all_jobs(idx_t& index, uint64_t n, const bench_config_t& cfg, const std::string& index_name)
{
    const bool do_ham  = cfg.metric == "all" || cfg.metric == "hamming";
    const bool do_edit = cfg.metric == "all" || cfg.metric == "edit";

    for (const bench_job_t& job : cfg.jobs) {
        if (job.metric == HAMMING_DISTANCE && !do_ham) continue;
        if (job.metric == EDIT_DISTANCE && !do_edit) continue;

        if (job.metric == EDIT_DISTANCE && job.k > EDIT_K_LIMIT) {
            std::cerr << "warning: k = " << job.k << " > " << EDIT_K_LIMIT
                      << " is not supported for edit distance; skipping " << job.path << std::endl;
            continue;
        }

        if (job.op == OP_COUNT && job.metric == EDIT_DISTANCE) {
            std::cerr << "warning: edit-distance count is not supported; skipping " << job.path << std::endl;
            continue;
        }

        search_scheme_t scheme = make_scheme(cfg, job.k);
        std::cerr << "measuring " << index_name << " on " << job.path << " ..." << std::endl;

        if (job.op == OP_COUNT) {
            run_count<idx_t>(index, n, cfg, job, scheme, index_name);
        } else {
            if (job.metric == HAMMING_DISTANCE) run_locate<idx_t, HAMMING_DISTANCE>(index, n, cfg, job, scheme, index_name);
            else                                run_locate<idx_t, EDIT_DISTANCE>(index, n, cfg, job, scheme, index_name);
        }
    }
}

// ############################# PER-INDEX MEASUREMENT #############################

void measure_bmove(const bench_config_t& cfg)
{
    // b-move stores its index in several files sharing the base name <bmove_dir>/<text_name>
    std::string base = (std::filesystem::path(cfg.bmove_dir) / cfg.text_name).string();

    // BMove's loader does not fail gracefully on missing files (it aborts), so check
    // for the core index files up front and skip the b-move measurement if they are absent
    for (const char* ext : {".cct", ".plcp", ".smpf", ".smpl"}) {
        if (!std::filesystem::exists(base + ext)) {
            throw std::runtime_error("b-move index file " + base + ext + " not found");
        }
    }

    std::cerr << "loading the b-move index from " << base << " ..." << std::endl;
    BMove bmove_index(base, false);
    b_move_adapter index(bmove_index);

    measure_all_jobs<b_move_adapter>(index, index.input_size(), cfg, "bmove");
}

void measure_br_index(const bench_config_t& cfg)
{
    std::cerr << "loading the br-index from " << cfg.br_index_path << " ..." << std::endl;
    br_index_adapter::br_index_t br_index;
    std::ifstream in(cfg.br_index_path);

    if (!in.good()) {
        std::cerr << "error: cannot open br-index file " << cfg.br_index_path << std::endl;
        return;
    }

    br_index.load(in);
    in.close();

    br_index_adapter index(br_index);
    measure_all_jobs<br_index_adapter>(index, index.input_size(), cfg, "br_index");
}

namespace {

/**
 * @brief loads a move_rb index from path and measures all jobs
 * @tparam support _locate_move or _locate_rlzsa
 * @tparam pos_t index integer type (uint32_t / uint64_t)
 */
template <move_r_support support, typename pos_t>
void measure_move_rb_impl(const bench_config_t& cfg, const std::string& path, const std::string& index_name)
{
    using idx_t = move_rb<support, char, pos_t>;
    idx_t index;

    std::ifstream in(path);
    if (!in.good()) {
        std::cerr << "error: cannot open move_rb index file " << path << std::endl;
        return;
    }
    index.load(in);
    in.close();

    measure_all_jobs<idx_t>(index, index.forward_index().input_size(), cfg, index_name);
}

/**
 * @brief reads the bit-width flag of a serialized move_rb index
 * @return true <=> the index was built with a 64-bit position type
 */
bool index_is_64_bit(const std::string& path)
{
    std::ifstream in(path);
    bool is_64_bit = false;
    in.read((char*) &is_64_bit, 1);
    in.close();
    return is_64_bit;
}

} // namespace

void measure_move_rb_move(const bench_config_t& cfg)
{
    std::cerr << "loading the move_rb (move) index from " << cfg.move_rb_move_path << " ..." << std::endl;

    if (index_is_64_bit(cfg.move_rb_move_path)) {
        measure_move_rb_impl<_locate_move, uint64_t>(cfg, cfg.move_rb_move_path, "move_rb_move");
    } else {
        measure_move_rb_impl<_locate_move, uint32_t>(cfg, cfg.move_rb_move_path, "move_rb_move");
    }
}

void measure_move_rb_rlzsa(const bench_config_t& cfg)
{
    std::cerr << "loading the move_rb (rlzsa) index from " << cfg.move_rb_rlzsa_path << " ..." << std::endl;

    if (index_is_64_bit(cfg.move_rb_rlzsa_path)) {
        measure_move_rb_impl<_locate_rlzsa, uint64_t>(cfg, cfg.move_rb_rlzsa_path, "move_rb_rlzsa");
    } else {
        measure_move_rb_impl<_locate_rlzsa, uint32_t>(cfg, cfg.move_rb_rlzsa_path, "move_rb_rlzsa");
    }
}

// ############################# DRIVER #############################

void help(const std::string& msg)
{
    if (!msg.empty()) std::cout << msg << std::endl;
    std::cout << "move-rb-bench: benchmarks approximate count- and locate-performance of" << std::endl;
    std::cout << "               b-move, br-index and move_rb (move & rlzsa). The indexes must be" << std::endl;
    std::cout << "               built beforehand (each with its own build tool); the pattern sets" << std::endl;
    std::cout << "               are generated by move-rb-queries." << std::endl << std::endl;
    std::cout << "usage: move-rb-bench [-s <scheme>] [--bmove <dir>] [--br-index <file>] <move_rb_move_file> <move_rb_rlzsa_file> <text_name> <patterns_path>" << std::endl;
    std::cout << "   -s <scheme>          search scheme to use (default: min_u): one of" << std::endl;
    std::cout << "                        pigeon_hole, suffix_filter, min_u, 01, or a path to a search-scheme file." << std::endl;
    std::cout << "                        For the built-in schemes the number of errors k is taken from each pattern" << std::endl;
    std::cout << "                        file name; a scheme file fixes k itself." << std::endl;
    std::cout << "   --bmove <dir>        (optional) directory containing the b-move index files (base name <dir>/<text_name>);" << std::endl;
    std::cout << "                        if omitted, b-move is not measured (b-move only supports the ACGT alphabet)" << std::endl;
    std::cout << "   --br-index <file>    (optional) path to the br-index (.bri) file; if omitted, br-index is not measured" << std::endl;
    std::cout << "   --only <index>       (optional) measure only one index: bmove|br_index|move_rb_move|move_rb_rlzsa" << std::endl;
    std::cout << "   --metric <metric>    (optional) run only one distance metric: all|hamming|edit (default: all)" << std::endl;
    std::cout << "   <move_rb_move_file>  path to the move_rb index built with move data structures" << std::endl;
    std::cout << "   <move_rb_rlzsa_file> path to the move_rb index built with an rlzsa" << std::endl;
    std::cout << "   <text_name>          name of the original text (used in the output and pattern file names)" << std::endl;
    std::cout << "   <patterns_path>      directory containing the pattern files (from move-rb-queries) of the form" << std::endl;
    std::cout << "                        <text_name>.patterns-{count,locate}-{hamming,edit}-k<value>-m<value>" << std::endl;
    exit(0);
}

/**
 * @brief resolves cfg.scheme: a built-in scheme name is kept as-is; anything else is
 *        treated as a path to a search-scheme file, which is parsed into cfg.file_scheme
 */
void resolve_scheme(bench_config_t& cfg)
{
    if (cfg.scheme == "pigeon_hole" ||
        cfg.scheme == "suffix_filter" ||
        cfg.scheme == "min_u" ||
        cfg.scheme == "01"
    ) return;

    if (!std::filesystem::exists(cfg.scheme)) {
        help("error: -s must be pigeon_hole, suffix_filter, min_u, 01 or a path to a search-scheme file");
    }

    uint64_t size = std::filesystem::file_size(cfg.scheme);
    std::string content;
    no_init_resize(content, size);
    std::ifstream f(cfg.scheme);
    f.read(content.data(), size);
    cfg.file_scheme = parse_search_scheme(content);
    cfg.scheme_from_file = true;
}

/**
 * @brief tries to parse a pattern file name of the form
 *        <text_name>.patterns-{count,locate}-{hamming,edit}-k<value>-m<value>
 * @return whether name matched (and job was filled)
 */
bool parse_pattern_file_name(const std::string& name, const std::string& text_name, const std::string& path, bench_job_t& job)
{
    const std::string prefix = text_name + ".patterns-";
    if (name.compare(0, prefix.size(), prefix) != 0) return false;
    std::string rest = name.substr(prefix.size());

    if (rest.compare(0, 6, "count-") == 0) {
        job.op = OP_COUNT;
        rest = rest.substr(6);
    } else if (rest.compare(0, 7, "locate-") == 0) {
        job.op = OP_LOCATE;
        rest = rest.substr(7);
    } else {
        return false;
    }

    if (rest.compare(0, 8, "hamming-") == 0) {
        job.metric = HAMMING_DISTANCE;
        rest = rest.substr(8);
    } else if (rest.compare(0, 5, "edit-") == 0) {
        job.metric = EDIT_DISTANCE;
        rest = rest.substr(5);
    } else {
        return false;
    }

    // rest must now be "k<value>-m<value>"
    if (rest.empty() || rest[0] != 'k') return false;
    size_t dash = rest.find("-m");
    if (dash == std::string::npos) return false;

    std::string k_str = rest.substr(1, dash - 1);
    std::string m_str = rest.substr(dash + 2);

    if (k_str.empty() || m_str.empty()) return false;
    for (char c : k_str) if (!std::isdigit((unsigned char) c)) return false;
    for (char c : m_str) if (!std::isdigit((unsigned char) c)) return false;

    job.k = std::stoll(k_str);
    job.m = std::stoull(m_str);
    job.path = path;
    return true;
}

int main(int argc, char** argv)
{
    bench_config_t cfg;
    int arg_idx = 1;

    // optional flags
    while (arg_idx < argc && argv[arg_idx][0] == '-') {
        std::string opt = argv[arg_idx++];
        if (opt == "-s") {
            if (arg_idx >= argc) help("error: missing scheme after -s");
            cfg.scheme = argv[arg_idx++];
        } else if (opt == "--bmove") {
            if (arg_idx >= argc) help("error: missing directory after --bmove");
            cfg.bmove_dir = argv[arg_idx++];
        } else if (opt == "--br-index") {
            if (arg_idx >= argc) help("error: missing file after --br-index");
            cfg.br_index_path = argv[arg_idx++];
        } else if (opt == "--only") {
            if (arg_idx >= argc) help("error: missing index name after --only");
            cfg.only = argv[arg_idx++];
        } else if (opt == "--metric") {
            if (arg_idx >= argc) help("error: missing metric after --metric");
            cfg.metric = argv[arg_idx++];
        } else {
            help("error: unrecognized option '" + opt + "'");
        }
    }

    // the four positional arguments
    if (argc - arg_idx != 4) help("");
    cfg.move_rb_move_path = argv[arg_idx++];
    cfg.move_rb_rlzsa_path = argv[arg_idx++];
    cfg.text_name = argv[arg_idx++];
    cfg.patterns_path = argv[arg_idx++];

    resolve_scheme(cfg);

    if (!std::filesystem::is_directory(cfg.patterns_path)) {
        help("error: <patterns_path> is not a directory");
    }

    // collect all matching pattern files
    for (const auto& entry : std::filesystem::directory_iterator(cfg.patterns_path)) {
        if (!entry.is_regular_file()) continue;
        bench_job_t job;
        if (parse_pattern_file_name(entry.path().filename().string(), cfg.text_name, entry.path().string(), job)) {
            cfg.jobs.emplace_back(job);
        }
    }

    if (cfg.jobs.empty()) {
        help("error: no pattern files of the form <text_name>.patterns-{count,locate}-{hamming,edit}-k<value>-m<value> found in <patterns_path>");
    }

    // deterministic order: count before locate, then by k, then by m
    std::sort(cfg.jobs.begin(), cfg.jobs.end(), [](const bench_job_t& a, const bench_job_t& b) {
        if (a.op != b.op) return a.op < b.op;
        if (a.k != b.k) return a.k < b.k;
        return a.m < b.m;
    });

    std::cerr << "found " << cfg.jobs.size() << " pattern file(s) for text '" << cfg.text_name << "'" << std::endl;

    // measure each index independently; a failure to load one index (e.g. a missing
    // file) only skips that index instead of aborting the whole benchmark
    auto guarded = [](const char* name, void (*fnc)(const bench_config_t&), const bench_config_t& cfg) {
        try {
            fnc(cfg);
        } catch (const std::exception& e) {
            std::cerr << "warning: skipping " << name << " (" << e.what() << ")" << std::endl;
        }
    };

    // optionally restrict to a single index (--only)
    auto want = [&](const char* name) { return cfg.only.empty() || cfg.only == name; };

    // b-move and br-index are optional; only measure them if their paths were given
    if (!cfg.bmove_dir.empty()     && want("bmove"))     guarded("b-move", measure_bmove, cfg);
    if (!cfg.br_index_path.empty() && want("br_index"))  guarded("br-index", measure_br_index, cfg);
    if (want("move_rb_move"))  guarded("move_rb (move)", measure_move_rb_move, cfg);
    if (want("move_rb_rlzsa")) guarded("move_rb (rlzsa)", measure_move_rb_rlzsa, cfg);

    return 0;
}
