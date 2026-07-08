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

// move-rb-bench-ext: benchmarks the bidirectional-extension primitive (see ext_query.hpp) of move_rb (move &
// rlzsa), both columba flavors and br-index. Count and locate use separate pattern files (from
// move-rb-gen-ext-queries): a count file measures the search alone, a locate file the search plus enumerating all
// occurrences. Each index is passed by its own flag and measured only if given.

#include "ext_query.hpp" // ext_plan_t / make_ext_plan / ext_count / ext_locate (the shared query model)
#include "indexes/columba_capi.hpp" // columba native index (both flavors) behind a C ABI, via shared plugins
#include "indexes/br_index_adapter.hpp" // br-index: move-r-apm adapter (exposes the shared bidirectional interface)
#include <move_rb/move_rb.hpp>

#include <algorithm>
#include <clocale>
#include <cstdint>
#include <cstdlib>
#include <dlfcn.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <unistd.h>
#include <vector>

#include <malloc_count.h>

// ############################# CONFIGURATION #############################

// the operation a pattern file is measured with (count and locate get separate files)
enum ext_op_t : uint8_t {
    OP_COUNT = 0,  // the bidirectional search alone (its result is the occurrence count)
    OP_LOCATE = 1  // the bidirectional search followed by enumerating all occurrences
};

// a single pattern file to be measured (one operation, one pattern length m)
struct bench_job_t {
    ext_op_t op; // count or locate (from the file name)
    uint64_t m; // pattern length (from the file name / header)
    std::string path; // path to the pattern file
};

// the configuration of one benchmark run
struct bench_config_t {
    std::string text_name; // name of the original text (used in the output and pattern file names)
    std::string patterns_path; // directory containing the pattern files
    // per-index locations, each set from its own optional flag; an index is measured iff its path is given
    std::string move_rb_move_path;  // --move-rb: path to the move_rb index using move data structures
    std::string move_rb_rlzsa_path; // --move-rb-rlzsa: path to the move_rb index using an rlzsa
    std::string br_index_path;      // --bri: path to the br-index (.bri) file
    std::string columba_path;       // --columba: base name of the columba FM-index
    std::string columba_rlc_path;   // --columba-rlc: base name of the columba-rlc (b-move) index
    std::vector<bench_job_t> jobs;  // the pattern files to measure

    // each pattern file's query set is repeated until at least this many seconds have been measured for it
    // (>= one full pass), so short-running sets still accumulate a stable timing
    double min_time_seconds = 10.0;

    // which operations to measure (from --only): count files, locate files or both
    bool do_count = true;
    bool do_locate = true;

    // if non-empty, only pattern files with these lengths m (--m) are measured
    std::vector<uint64_t> m_filter;
};

// ############################# RESULT OUTPUT #############################

// per-pattern-set memory measurement (bytes); peak_mem is the peak additional heap allocated while measuring the
// set, occ_mem the peak footprint of a single pattern's occurrence vector. Both 0 for count / when malloc_count
// is disabled.
struct mem_stats_t {
    uint64_t peak_mem = 0;
    uint64_t occ_mem = 0;
};

// in-memory footprint (bytes) of the index currently being measured, captured as the malloc_count heap growth
// while the index is loaded. Set once per index (loads are sequential) and emitted on every RESULT line of that
// index. 0 if malloc_count is disabled.
static uint64_t g_index_mem = 0;

// measures the heap growth caused by running load_fn() and stores it in g_index_mem
template <typename load_fnc_t>
inline void load_measured(load_fnc_t load_fn)
{
    const uint64_t before = malloc_count_current();
    load_fn();
    const uint64_t after = malloc_count_current();
    g_index_mem = after > before ? after - before : 0;
}

/**
 * @brief writes a single RESULT line to std::cout. The schema matches move-rb-bench-apm (so the same parsers work);
 *        the fields that do not apply to exact extension are fixed: dist_metr=exact, search_scheme=none,
 *        max_mismatches=0, cigar=0
 */
inline void print_result(
    const std::string& algo, const std::string& text, uint64_t n,
    uint64_t pattern_length, uint64_t num_patterns, uint64_t num_occurrences,
    const std::string& time_key, uint64_t time_value, const mem_stats_t& mem = {})
{
    std::cout << "RESULT"
              << " algo=" << algo
              << " text=" << text
              << " n=" << n
              << " dist_metr=exact"
              << " cigar=0"
              << " search_scheme=none"
              << " pattern_length=" << pattern_length
              << " num_patterns=" << num_patterns
              << " max_mismatches=0"
              << " num_occurrences=" << num_occurrences
              << " index_mem=" << g_index_mem
              << " peak_mem=" << mem.peak_mem
              << " occ_mem=" << mem.occ_mem
              << " cigar_mem=0"
              << " " << time_key << "=" << time_value
              << std::endl;
}

// ############################# PATTERN FILE READING #############################

// reads a whole pattern file (plain pizza&chili header + the fixed-length patterns); returns false (with a
// warning) if the file cannot be opened, is empty or is truncated
inline bool read_patterns(const bench_job_t& job, std::vector<std::string>& patterns, uint64_t& m)
{
    std::ifstream pf(job.path, std::ios::binary);
    if (!pf.good()) { std::cerr << "warning: cannot open pattern file " << job.path << std::endl; return false; }

    std::string header;
    std::getline(pf, header);
    const uint64_t num_patterns = number_of_patterns(header);
    m = patterns_length(header);

    if (num_patterns == 0) { std::cerr << "warning: no patterns in " << job.path << "; skipping" << std::endl; return false; }
    patterns.resize(num_patterns);
    for (uint64_t i = 0; i < num_patterns; i++) { no_init_resize(patterns[i], m); pf.read(patterns[i].data(), m); }
    if (!pf) { std::cerr << "warning: " << job.path << " is truncated; skipping" << std::endl; return false; }
    return true;
}

// builds the extension plan of each pattern in a set (from the plan RNG reset to ext_plan_seed)
inline std::vector<ext_plan_t> build_plans(uint64_t m, uint64_t num_patterns)
{
    std::vector<ext_plan_t> plans(num_patterns);
    std::mt19937_64 plan_rng(ext_plan_seed);
    for (uint64_t i = 0; i < num_patterns; i++) plans[i] = make_ext_plan(m, plan_rng);
    return plans;
}

// ############################# GENERIC MEASUREMENT #############################

/**
 * @brief measures one pattern file (one operation) for a directly-linked index (move_rb or the br-index adapter):
 *        a count file with the bidirectional search alone, a locate file with the search plus enumerating all
 *        occurrences, replayed until at least cfg.min_time_seconds have been measured, and writes the RESULT line
 * @tparam idx_t the index type (move_rb or br_index_adapter)
 */
template <typename idx_t>
void run_ext(idx_t& index, uint64_t n, const bench_config_t& cfg, const bench_job_t& job, const std::string& index_name)
{
    using pos_t = typename idx_t::pos_type;

    std::vector<std::string> patterns;
    uint64_t m;
    if (!read_patterns(job, patterns, m)) return;
    const uint64_t num_patterns = patterns.size();
    const std::vector<ext_plan_t> plans = build_plans(m, num_patterns);
    const uint64_t min_time_ns = (uint64_t) (cfg.min_time_seconds * 1e9);

    if (job.op == OP_COUNT) {
        // counting: the bidirectional search alone (the result is the occurrence count)
        uint64_t num_occurrences = 0, time_count = 0, passes = 0;
        do {
            const uint64_t before = time_count;
            for (uint64_t i = 0; i < num_patterns; i++) {
                auto t1 = now();
                pos_t c = ext_count<idx_t>(index, patterns[i], plans[i]);
                time_count += time_diff_ns(t1);
                num_occurrences += c;
            }
            passes++;
            if (time_count == before) break;
        } while (time_count < min_time_ns);

        print_result("count_" + index_name, cfg.text_name, n, m, num_patterns,
            num_occurrences / passes, "time_count", time_count / passes);
    } else {
        // locating: the search followed by enumerating all occurrences
        std::vector<pos_t> occ;
        uint64_t num_occurrences = 0, time_locate = 0, passes = 0, peak_occ = 0;
        const uint64_t mem_baseline = malloc_count_current();
        malloc_count_reset_peak();
        do {
            const uint64_t before = time_locate;
            for (uint64_t i = 0; i < num_patterns; i++) {
                occ.clear();
                auto t1 = now();
                ext_locate<idx_t>(index, patterns[i], plans[i], [&](pos_t o) { occ.emplace_back(o); });
                time_locate += time_diff_ns(t1);
                num_occurrences += occ.size();
                peak_occ = std::max<uint64_t>(peak_occ, occ.size());
            }
            passes++;
            if (time_locate == before) break;
        } while (time_locate < min_time_ns);

        mem_stats_t mem;
        mem.peak_mem = malloc_count_peak() - mem_baseline;
        mem.occ_mem = peak_occ * sizeof(pos_t);
        print_result("locate_" + index_name, cfg.text_name, n, m, num_patterns,
            num_occurrences / passes, "time_locate", time_locate / passes, mem);
    }
}

/**
 * @brief measures every job for the given (already loaded) directly-linked index
 * @tparam idx_t the index type (move_rb or br_index_adapter)
 */
template <typename idx_t>
void measure_all_jobs(idx_t& index, uint64_t n, const bench_config_t& cfg, const std::string& index_name)
{
    for (const bench_job_t& job : cfg.jobs) {
        std::cerr << "measuring " << index_name << " on " << job.path << " ..." << std::endl;
        run_ext<idx_t>(index, n, cfg, job, index_name);
    }
}

// ############################# PER-INDEX MEASUREMENT: move_rb & br-index #############################

void measure_br_index(const bench_config_t& cfg)
{
    std::cerr << "loading the br-index from " << cfg.br_index_path << " ..." << std::endl;
    br_index_adapter::br_index_t br_index;
    std::ifstream in(cfg.br_index_path);
    if (!in.good()) throw std::runtime_error("br-index file " + cfg.br_index_path + " not found");
    load_measured([&]{ br_index.load(in); });
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
    if (!in.good()) throw std::runtime_error("move_rb index file " + path + " not found");
    load_measured([&]{ index.load(in); });
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
    if (index_is_64_bit(cfg.move_rb_move_path))
        measure_move_rb_impl<_locate_move, uint64_t>(cfg, cfg.move_rb_move_path, "move_rb_move");
    else
        measure_move_rb_impl<_locate_move, uint32_t>(cfg, cfg.move_rb_move_path, "move_rb_move");
}

void measure_move_rb_rlzsa(const bench_config_t& cfg)
{
    std::cerr << "loading the move_rb (rlzsa) index from " << cfg.move_rb_rlzsa_path << " ..." << std::endl;
    if (index_is_64_bit(cfg.move_rb_rlzsa_path))
        measure_move_rb_impl<_locate_rlzsa, uint64_t>(cfg, cfg.move_rb_rlzsa_path, "move_rb_rlzsa");
    else
        measure_move_rb_impl<_locate_rlzsa, uint32_t>(cfg, cfg.move_rb_rlzsa_path, "move_rb_rlzsa");
}

// ############################# PER-INDEX MEASUREMENT: columba (via plugin) #############################

// directory of the running executable (the columba plugin .so files sit next to it in build/bench/)
inline std::string executable_dir()
{
    char buf[4096];
    ssize_t len = ::readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    if (len <= 0) return ".";
    buf[len] = '\0';
    std::string path(buf);
    size_t slash = path.find_last_of('/');
    return slash == std::string::npos ? std::string(".") : path.substr(0, slash);
}

// measures one columba flavor's bidirectional-extension performance through its plugin's C ABI. Each pattern file
// is one operation (count or locate); the extension plans are rebuilt from the fixed ext_plan_seed exactly as for
// the directly-linked indexes and passed to the plugin in the C-ABI's 0=LEFT/1=RIGHT encoding
void measure_columba_ext_api(const columba_api_t& api, const bench_config_t& cfg, const std::string& base)
{
    const std::string flavor = api.flavor;
    std::cerr << "loading the " << flavor << " (native) index from " << base << " ..." << std::endl;

    void* handle = nullptr;
    // the search scheme is irrelevant to exact extension; load with a valid built-in name
    load_measured([&]{ handle = api.load(base.c_str(), "columba"); });
    if (handle == nullptr) {
        std::cerr << "warning: skipping " << flavor << " (index files not found at " << base << ")" << std::endl;
        return;
    }
    const uint64_t n = api.input_size(handle);
    void* ext = api.ext_make(handle);

    for (const bench_job_t& job : cfg.jobs) {
        std::vector<std::string> patterns;
        uint64_t m;
        std::cerr << "measuring " << flavor << " (ext) on " << job.path << " ..." << std::endl;
        if (!read_patterns(job, patterns, m)) continue;
        const uint64_t num_patterns = patterns.size();
        const std::vector<ext_plan_t> plans = build_plans(m, num_patterns);

        // the plugin's C ABI takes the order as 0=LEFT/1=RIGHT bytes; precompute them (out of the timed region)
        std::vector<std::vector<uint8_t>> orders(num_patterns);
        for (uint64_t i = 0; i < num_patterns; i++) {
            orders[i].resize(plans[i].order.size());
            for (size_t j = 0; j < plans[i].order.size(); j++) orders[i][j] = plans[i].order[j] == LEFT ? 0 : 1;
        }

        const uint64_t min_time_ns = (uint64_t) (cfg.min_time_seconds * 1e9);
        const int want_locate = job.op == OP_LOCATE ? 1 : 0;

        uint64_t num_occurrences = 0, time = 0, passes = 0, occ_mem = 0;
        const uint64_t mem_baseline = malloc_count_current();
        malloc_count_reset_peak();
        do {
            const uint64_t before = time;
            for (uint64_t i = 0; i < num_patterns; i++) {
                uint64_t ob = 0;
                auto t1 = now();
                uint64_t c = api.ext_run(ext, patterns[i].data(), m, plans[i].start, orders[i].data(), want_locate, &ob);
                time += time_diff_ns(t1);
                num_occurrences += c;
                occ_mem = std::max(occ_mem, ob);
            }
            passes++;
            if (time == before) break;
        } while (time < min_time_ns);

        if (job.op == OP_COUNT) {
            print_result("count_" + flavor, cfg.text_name, n, m, num_patterns,
                num_occurrences / passes, "time_count", time / passes);
        } else {
            mem_stats_t mem;
            mem.peak_mem = malloc_count_peak() - mem_baseline;
            mem.occ_mem = occ_mem;
            print_result("locate_" + flavor, cfg.text_name, n, m, num_patterns,
                num_occurrences / passes, "time_locate", time / passes, mem);
        }
    }

    api.ext_free(ext);
    api.destroy(handle);
}

// dlopen()s one columba flavor plugin (RTLD_LOCAL keeps its columba symbols private, so both flavors coexist in
// one process) and measures its bidirectional extension through the api getter @p api_symbol
void measure_columba_plugin(const std::string& so_name, const std::string& api_symbol,
                            const bench_config_t& cfg, const std::string& base)
{
    const std::string so_path = executable_dir() + "/" + so_name;
    void* lib = dlopen(so_path.c_str(), RTLD_NOW | RTLD_LOCAL);
    if (lib == nullptr) {
        std::cerr << "warning: cannot load columba plugin " << so_path << " (" << dlerror() << ")" << std::endl;
        return;
    }
    using api_getter_t = const columba_api_t* (*)();
    api_getter_t getter = reinterpret_cast<api_getter_t>(dlsym(lib, api_symbol.c_str()));
    if (getter == nullptr) {
        std::cerr << "warning: " << api_symbol << " not found in " << so_path << std::endl;
        dlclose(lib);
        return;
    }
    measure_columba_ext_api(*getter(), cfg, base);
    // the plugin is intentionally left mapped (not dlclose'd) until process exit
}

// reads the index-position width (sizeof(length_t) in bytes: 4 for 32-bit, 8 for 64-bit) that the columba index
// at base name @p base was built with, from the second line of its <base>.meta. Returns 0 if the meta is missing.
int columba_meta_width_bytes(const std::string& base)
{
    std::ifstream ifs(base + ".meta");
    long tag = 0, width_bytes = 0;
    if (ifs && (ifs >> tag >> width_bytes)) return (int) width_bytes;
    return 0;
}

// measures one columba flavor: picks the 32- or 64-bit plugin from the width recorded in the index's meta file
void measure_columba_flavor(const bench_config_t& cfg, const std::string& base,
                            const std::string& plugin_base, const char* api_symbol)
{
    const int width = columba_meta_width_bytes(base);
    if (width == 0)
        std::cerr << "warning: cannot read " << base << ".meta; defaulting to the 64-bit plugin" << std::endl;
    const std::string so = "lib" + plugin_base + (width == 4 ? "_32" : "_64") + ".so";
    std::cerr << "columba index at " << base << " is " << (width ? width * 8 : 64) << "-bit; loading " << so << std::endl;
    measure_columba_plugin(so, api_symbol, cfg, base);
}

// measures the given columba flavors (run-length-compressed = b-move, and the FM-index), each via its plugin
void measure_columba_native(const bench_config_t& cfg)
{
    if (!cfg.columba_rlc_path.empty())
        measure_columba_flavor(cfg, cfg.columba_rlc_path, "columba_rlc_plugin", "columba_rlc_api");
    if (!cfg.columba_path.empty())
        measure_columba_flavor(cfg, cfg.columba_path, "columba_fm_plugin", "columba_fm_api");
}

// ############################# DRIVER #############################

// splits a comma-separated string into its non-empty parts
inline std::vector<std::string> split_csv(const std::string& list)
{
    std::vector<std::string> parts;
    for (size_t start = 0; start <= list.size();) {
        size_t comma = list.find(',', start);
        std::string part = list.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
        if (!part.empty()) parts.push_back(part);
        if (comma == std::string::npos) break;
        start = comma + 1;
    }
    return parts;
}

void help(const std::string& msg)
{
    if (!msg.empty()) std::cout << msg << std::endl;
    std::cout << "move-rb-bench-ext: benchmarks the bidirectional-extension primitive (exact count and subsequent" << std::endl;
    std::cout << "               locate) shared by move_rb (move & rlzsa), both columba flavors and br-index. Each" << std::endl;
    std::cout << "               query starts at one fixed random position inside the pattern and extends it left" << std::endl;
    std::cout << "               and right in a random order until the whole pattern is matched. The indexes must be" << std::endl;
    std::cout << "               built beforehand; each is measured only if its path is given below. The pattern sets" << std::endl;
    std::cout << "               are generated by move-rb-gen-ext-queries." << std::endl << std::endl;
    std::cout << "usage: move-rb-bench-ext [...] <index flags...> <text_name> <patterns_path>" << std::endl;
    std::cout << "   --only <op>          (optional) measure only count or locate files: count, locate or both (default: both)" << std::endl;
    std::cout << "   --m <list>           (optional) measure only pattern files with these pattern lengths m (comma-separated)" << std::endl;
    std::cout << "   --time <s>           (optional) replay each pattern file's query set until at least <s> seconds" << std::endl;
    std::cout << "                        have been measured for it (>= one pass); the reported time_* is then the" << std::endl;
    std::cout << "                        average time for one pass over the set. 0 = a single pass; default: 10" << std::endl;
    std::cout << "   index flags          (all optional, all independent) the location of each index to measure;" << std::endl;
    std::cout << "                        an index is measured only if its flag is given. At least one is required:" << std::endl;
    std::cout << "     --move-rb <file>         move_rb index using move data structures" << std::endl;
    std::cout << "     --move-rb-rlzsa <file>   move_rb index using an rlzsa" << std::endl;
    std::cout << "     --bri <file>             br-index (.bri)" << std::endl;
    std::cout << "     --columba <base>         columba FM-index base name (columba is ACGT-only)" << std::endl;
    std::cout << "     --columba-rlc <base>     columba-rlc (b-move) index base name (ACGT-only)" << std::endl;
    std::cout << "   <text_name>          name of the original text (used in the output and pattern file names)" << std::endl;
    std::cout << "   <patterns_path>      directory containing the pattern files (from move-rb-gen-ext-queries) of the" << std::endl;
    std::cout << "                        form <text_name>.patterns-{count,locate}-m<value>" << std::endl;
    exit(0);
}

/**
 * @brief tries to parse a pattern file name of the form <text_name>.patterns-{count,locate}-m<value>
 * @return whether name matched (and job was filled)
 */
bool parse_pattern_file_name(const std::string& name, const std::string& text_name, const std::string& path, bench_job_t& job)
{
    const std::string prefix = text_name + ".patterns-";
    if (name.compare(0, prefix.size(), prefix) != 0) return false;
    std::string rest = name.substr(prefix.size());

    if (rest.compare(0, 7, "count-m") == 0) {
        job.op = OP_COUNT;
        rest = rest.substr(7);
    } else if (rest.compare(0, 8, "locate-m") == 0) {
        job.op = OP_LOCATE;
        rest = rest.substr(8);
    } else {
        return false;
    }

    if (rest.empty()) return false;
    for (char c : rest) if (!std::isdigit((unsigned char) c)) return false;
    job.m = std::stoull(rest);
    job.path = path;
    return true;
}

int main(int argc, char** argv)
{
    std::setlocale(LC_ALL, "C");
    bench_config_t cfg;
    int arg_idx = 1;

    // optional flags
    while (arg_idx < argc && argv[arg_idx][0] == '-') {
        std::string opt = argv[arg_idx++];
        if (opt == "--time") {
            if (arg_idx >= argc) help("error: missing seconds after --time");
            cfg.min_time_seconds = std::stod(argv[arg_idx++]);
            if (cfg.min_time_seconds < 0) help("error: --time must be >= 0");
        } else if (opt == "--only") {
            if (arg_idx >= argc) help("error: missing operation after --only");
            std::string op = argv[arg_idx++];
            if      (op == "count")  { cfg.do_count = true;  cfg.do_locate = false; }
            else if (op == "locate") { cfg.do_count = false; cfg.do_locate = true; }
            else if (op == "both")   { cfg.do_count = true;  cfg.do_locate = true; }
            else help("error: --only must be count, locate or both");
        } else if (opt == "--m") {
            if (arg_idx >= argc) help("error: missing m list after --m");
            for (const std::string& s : split_csv(argv[arg_idx++])) cfg.m_filter.push_back(std::stoull(s));
        } else if (opt == "--move-rb") {
            if (arg_idx >= argc) help("error: missing file after --move-rb");
            cfg.move_rb_move_path = argv[arg_idx++];
        } else if (opt == "--move-rb-rlzsa") {
            if (arg_idx >= argc) help("error: missing file after --move-rb-rlzsa");
            cfg.move_rb_rlzsa_path = argv[arg_idx++];
        } else if (opt == "--bri") {
            if (arg_idx >= argc) help("error: missing file after --bri");
            cfg.br_index_path = argv[arg_idx++];
        } else if (opt == "--columba") {
            if (arg_idx >= argc) help("error: missing base name after --columba");
            cfg.columba_path = argv[arg_idx++];
        } else if (opt == "--columba-rlc") {
            if (arg_idx >= argc) help("error: missing base name after --columba-rlc");
            cfg.columba_rlc_path = argv[arg_idx++];
        } else {
            help("error: unrecognized option '" + opt + "'");
        }
    }

    // the two positional arguments
    if (argc - arg_idx != 2) help("");
    cfg.text_name = argv[arg_idx++];
    cfg.patterns_path = argv[arg_idx++];

    // at least one index must be given (each is measured only if its path was passed)
    if (cfg.move_rb_move_path.empty() && cfg.move_rb_rlzsa_path.empty() && cfg.br_index_path.empty() &&
        cfg.columba_path.empty() && cfg.columba_rlc_path.empty()) {
        help("error: no index given; pass at least one of --move-rb, --move-rb-rlzsa, --bri, --columba, --columba-rlc");
    }

    if (!std::filesystem::is_directory(cfg.patterns_path)) {
        help("error: <patterns_path> is not a directory");
    }

    // collect all matching pattern files (restricted to the --m values if given)
    auto in_filter = [](const auto& filter, auto value) {
        return filter.empty() || std::find(filter.begin(), filter.end(), value) != filter.end();
    };
    for (const auto& entry : std::filesystem::directory_iterator(cfg.patterns_path)) {
        if (!entry.is_regular_file()) continue;
        bench_job_t job;
        if (parse_pattern_file_name(entry.path().filename().string(), cfg.text_name, entry.path().string(), job)) {
            const bool op_ok = job.op == OP_COUNT ? cfg.do_count : cfg.do_locate;
            if (op_ok && in_filter(cfg.m_filter, job.m)) cfg.jobs.emplace_back(job);
        }
    }

    if (cfg.jobs.empty()) {
        help("error: no pattern files of the form <text_name>.patterns-{count,locate}-m<value> found in <patterns_path>");
    }

    // deterministic order: count before locate, then by pattern length
    std::sort(cfg.jobs.begin(), cfg.jobs.end(), [](const bench_job_t& a, const bench_job_t& b) {
        if (a.op != b.op) return a.op < b.op;
        return a.m < b.m;
    });

    std::cerr << "found " << cfg.jobs.size() << " pattern file(s) for text '" << cfg.text_name << "'" << std::endl;

    // measure each index independently; a failure to load one index only reports that index as skipped
    auto guarded = [](const char* name, void (*fnc)(const bench_config_t&), const bench_config_t& cfg) {
        try {
            fnc(cfg);
        } catch (const std::exception& e) {
            std::cerr << "warning: skipping " << name << " (" << e.what() << ")" << std::endl;
        }
    };

    if (!cfg.br_index_path.empty())     guarded("br-index", measure_br_index, cfg);
    if (!cfg.columba_rlc_path.empty() ||
        !cfg.columba_path.empty())      guarded("columba (native)", measure_columba_native, cfg);
    if (!cfg.move_rb_move_path.empty()) guarded("move_rb (move)", measure_move_rb_move, cfg);
    if (!cfg.move_rb_rlzsa_path.empty()) guarded("move_rb (rlzsa)", measure_move_rb_rlzsa, cfg);

    return 0;
}
