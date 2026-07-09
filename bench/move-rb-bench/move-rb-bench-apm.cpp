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

#include "indexes/columba_capi.hpp" // columba native algorithm (both flavors) behind a C ABI, via shared plugins
#include "indexes/br_index_adapter.hpp" // br-index: move-r-apm adapter + native algorithm (br_index_native)
#include <move_rb/move_rb.hpp>

#include <algorithm>
#include <clocale>
#include <cstdint>
#include <dlfcn.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <unistd.h>
#include <vector>

#include <ips4o.hpp>
#include <ips2ra.hpp>
#include <malloc_count.h>

// ############################# CONFIGURATION #############################

// the maximum number of errors the edit-distance algorithm supports
static constexpr int EDIT_K_LIMIT = 20;

// the operation a pattern file is to be measured with
enum bench_op_t : uint8_t {
    OP_COUNT = 0, // approximate count
    OP_LOCATE = 1 // approximate locate
};

// which algorithm(s) to run for each index
enum algo_sel_t : uint8_t { ALGO_NATIVE = 0, ALGO_NON_NATIVE = 1, ALGO_BOTH = 2 };

// which CIGAR variant of locate to measure
enum cigar_sel_t : uint8_t { CIGAR_OFF = 0, CIGAR_BOTH = 1, CIGAR_ONLY = 2 };

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
    std::string text_name; // name of the original text (used in the output and pattern file names)
    std::string patterns_path; // directory containing the pattern files
    // per-index locations, each set from its own optional flag; an index is measured iff its path is given. The
    // indexes are fully independent (they may live in different directories / use different base names), so the
    // two columba flavors no longer have to share a base name (which would collide on .cct/.pos/.sna/.fsid/.meta)
    std::string move_rb_move_path;  // --move-rb: path to the move_rb index using move data structures
    std::string move_rb_rlzsa_path; // --move-rb-rlzsa: path to the move_rb index using an rlzsa
    std::string br_index_path;      // --bri: path to the br-index (.bri) file
    std::string columba_path;       // --columba: base name of the columba FM-index
    std::string columba_rlc_path;   // --columba-rlc: base name of the columba-rlc (b-move) index
    std::vector<bench_job_t> jobs; // the pattern files to measure

    // search scheme to use (like move-rb-locate's -s): the name of a built-in scheme
    // (pigeon_hole, suffix_filter, min_u, 01) or the path to a scheme file
    std::string scheme = "min_u";
    bool scheme_from_file = false; // true <=> scheme is a file; then file_scheme holds it (k is fixed by the file)
    search_scheme_t file_scheme; // the scheme parsed from the file (only used if scheme_from_file)

    std::string metric = "all"; // which distance metric to run: all|hamming|edit

    // which operations to measure (from --only): count files, locate files or both
    bool do_count = true;
    bool do_locate = true;

    // each pattern file's query set is repeated until at least this many seconds have been
    // measured for it (>= one full pass), so short-running sets still accumulate a stable timing
    double min_time_seconds = 10.0;

    algo_sel_t algo = ALGO_BOTH; // --algo: which algorithm(s) to run (native / non-native / both)
    cigar_sel_t cigar_mode = CIGAR_OFF; // --cigar: which CIGAR variant of locate to measure (off / both / only)

    // whether the native / non-native algorithms are to be run (from --algo)
    bool run_native() const { return algo == ALGO_NATIVE || algo == ALGO_BOTH; }
    bool run_non_native() const { return algo == ALGO_NON_NATIVE || algo == ALGO_BOTH; }
    // whether the plain / CIGAR-producing locate is to be measured (from --cigar)
    bool locate_plain() const { return cigar_mode == CIGAR_OFF || cigar_mode == CIGAR_BOTH; }
    bool locate_cigar() const { return cigar_mode == CIGAR_BOTH || cigar_mode == CIGAR_ONLY; }

    // if non-empty, only pattern files with these error counts k (--k) resp. lengths m (--m) are measured
    std::vector<int64_t> k_filter;
    std::vector<uint64_t> m_filter;
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

// per-pattern-set memory measurement (bytes). peak_mem is the peak additional heap allocated while measuring the
// set (malloc_count peak minus the baseline before it); occ_mem is the peak footprint of the occurrence vector
// (peak occurrence count x element size); cigar_mem is the peak footprint of the occurrences' CIGAR data. All 0
// when not measured (e.g. count, or malloc_count disabled).
struct mem_stats_t {
    uint64_t peak_mem = 0;
    uint64_t occ_mem = 0;
    uint64_t cigar_mem = 0;
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
 * @brief writes a single RESULT line to std::cout
 */
inline void print_result(
    const std::string& algo, const std::string& text, uint64_t n,
    const std::string& dist_metr, const std::string& search_scheme,
    uint64_t pattern_length, uint64_t num_patterns,
    int64_t max_mismatches, uint64_t num_occurrences,
    const std::string& time_key, uint64_t time_value, bool cigar = false,
    const mem_stats_t& mem = {})
{
    std::cout << "RESULT"
              << " algo=" << algo
              << " text=" << text
              << " n=" << n
              << " dist_metr=" << dist_metr
              << " cigar=" << (cigar ? 1 : 0)
              << " search_scheme=" << search_scheme
              << " pattern_length=" << pattern_length
              << " num_patterns=" << num_patterns
              << " max_mismatches=" << max_mismatches
              << " num_occurrences=" << num_occurrences
              << " index_mem=" << g_index_mem
              << " peak_mem=" << mem.peak_mem
              << " occ_mem=" << mem.occ_mem
              << " cigar_mem=" << mem.cigar_mem
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

    if (num_patterns == 0) {
        std::cerr << "warning: no patterns in " << job.path << "; skipping" << std::endl;
        return;
    }

    std::vector<std::string> patterns(num_patterns);
    for (uint64_t i = 0; i < num_patterns; i++) {
        no_init_resize(patterns[i], m);
        pf.read(patterns[i].data(), m);
    }
    if (!pf) {
        std::cerr << "warning: " << job.path << " is truncated (expected " << num_patterns
                  << " x " << m << " bytes); skipping" << std::endl;
        return;
    }

    const uint64_t min_time_ns = (uint64_t) (cfg.min_time_seconds * 1e9);
    uint64_t num_occurrences = 0;
    uint64_t time_count = 0;
    uint64_t passes = 0;

    do {
        const uint64_t time_before = time_count;
        for (uint64_t i = 0; i < num_patterns; i++) {
            auto t1 = now();
            pos_t c = index.count_hamming_dist(patterns[i], scheme);
            time_count += time_diff_ns(t1);
            num_occurrences += c;
        }
        passes++;
        if (time_count == time_before) break;
    } while (time_count < min_time_ns);

    print_result("count_" + index_name, cfg.text_name, n, "hamming", cfg.scheme,
        m, num_patterns, scheme.k, num_occurrences / passes, "time_count", time_count / passes);
}

/**
 * @brief measures approximate locate (w.r.t. dist_metr) of all patterns in job.path
 *        for the given index and writes a RESULT line
 * @tparam idx_t the index type (move_rb, b_move_adapter or br_index_adapter)
 * @tparam dist_metr the distance metric (HAMMING_DISTANCE or EDIT_DISTANCE)
 * @tparam mode whether to additionally produce a per-occurrence CIGAR alignment (measures the alignment cost)
 */
template <typename idx_t, distance_metric_t dist_metr, cigar_mode_t mode = NO_CIGAR>
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

    if (num_patterns == 0) {
        std::cerr << "warning: no patterns in " << job.path << "; skipping" << std::endl;
        return;
    }

    std::vector<std::string> patterns(num_patterns);
    for (uint64_t i = 0; i < num_patterns; i++) {
        no_init_resize(patterns[i], m);
        pf.read(patterns[i].data(), m);
    }
    if (!pf) {
        std::cerr << "warning: " << job.path << " is truncated (expected " << num_patterns
                  << " x " << m << " bytes); skipping" << std::endl;
        return;
    }

    std::vector<aprx_occ_t<pos_t, mode>> occurrences;
    const uint64_t min_time_ns = (uint64_t) (cfg.min_time_seconds * 1e9);
    uint64_t num_occurrences = 0;
    uint64_t time_locate = 0;
    uint64_t passes = 0;

    // peak occurrence count and peak CIGAR-data bytes over the (post-filter) result of any single pattern
    uint64_t peak_occ = 0, peak_cigar = 0;
    const uint64_t mem_baseline = malloc_count_current(); // heap already used before this set (mostly the index)
    malloc_count_reset_peak();

    do {
        const uint64_t time_before = time_locate;
        for (uint64_t i = 0; i < num_patterns; i++) {
            occurrences.clear();
            auto t1 = now();
            auto collect = [&](aprx_occ_t<pos_t, mode> occ) { occurrences.emplace_back(std::move(occ)); };
            // the CIGAR overload takes an explicit mode template argument; the plain path keeps the single-argument
            // form so the adapters (whose locate has no mode parameter) still compile
            if constexpr (mode == CIGAR) index.template locate<dist_metr, CIGAR>(patterns[i], scheme, collect);
            else                         index.template locate<dist_metr>(patterns[i], scheme, collect);

            // for edit distance the same occurrence may be reported by several search contexts;
            // deduplicate (as the move-rb-locate tool does) before counting
            if constexpr (dist_metr == EDIT_DISTANCE) {
                ips2ra::sort(occurrences.begin(), occurrences.end(), [](const auto& o){ return o.pos; });
                filter_edit_distance_occurrences<pos_t, mode>(occurrences, (pos_t) scheme.k);
            }

            time_locate += time_diff_ns(t1);
            num_occurrences += occurrences.size();

            // memory bookkeeping (outside the timed region): track the largest single-pattern result
            peak_occ = std::max<uint64_t>(peak_occ, occurrences.size());
            if constexpr (mode == CIGAR) {
                uint64_t cig = 0;
                for (const auto& o : occurrences) cig += (uint64_t) o.cigar.size() * sizeof(cigar_run_t);
                peak_cigar = std::max(peak_cigar, cig);
            }
        }
        passes++;
        if (time_locate == time_before) break;
    } while (time_locate < min_time_ns);

    mem_stats_t mem;
    mem.peak_mem = malloc_count_peak() - mem_baseline;
    mem.occ_mem = peak_occ * sizeof(aprx_occ_t<pos_t, mode>);
    mem.cigar_mem = peak_cigar;

    const std::string dist_str = dist_metr == HAMMING_DISTANCE ? "hamming" : "edit";
    print_result("locate_" + index_name, cfg.text_name, n, dist_str, cfg.scheme,
        m, num_patterns, scheme.k, num_occurrences / passes, "time_locate", time_locate / passes, mode == CIGAR, mem);
}

/**
 * @brief runs run_locate in CIGAR mode, but only for indexes whose locate can actually emit a per-occurrence
 *        CIGAR alignment (i.e. move_rb). For indexes that cannot (the b-move / br-index adapters, whose locate has
 *        no CIGAR overload) this is a no-op with a note, so --cigar never fails a run
 * @tparam idx_t the index type
 * @tparam dist_metr the distance metric (HAMMING_DISTANCE or EDIT_DISTANCE)
 */
template <typename idx_t, distance_metric_t dist_metr>
void run_locate_cigar(idx_t& index, uint64_t n, const bench_config_t& cfg, const bench_job_t& job,
                      const search_scheme_t& scheme, const std::string& index_name)
{
    using pos_t = typename idx_t::pos_type;

    if constexpr (requires(idx_t& idx, const std::string& P, const search_scheme_t& s) {
        idx.template locate<dist_metr, CIGAR>(P, s, [](aprx_occ_t<pos_t, CIGAR>){});
    }) {
        run_locate<idx_t, dist_metr, CIGAR>(index, n, cfg, job, scheme, index_name);
    } else {
        std::cerr << "note: " << index_name << " cannot produce CIGARs; skipping its CIGAR locate for "
                  << job.path << std::endl;
    }
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
        } else if (job.metric == HAMMING_DISTANCE) {
            if (cfg.locate_plain()) run_locate<idx_t, HAMMING_DISTANCE>(index, n, cfg, job, scheme, index_name);
            if (cfg.locate_cigar()) run_locate_cigar<idx_t, HAMMING_DISTANCE>(index, n, cfg, job, scheme, index_name);
        } else {
            if (cfg.locate_plain()) run_locate<idx_t, EDIT_DISTANCE>(index, n, cfg, job, scheme, index_name);
            if (cfg.locate_cigar()) run_locate_cigar<idx_t, EDIT_DISTANCE>(index, n, cfg, job, scheme, index_name);
        }
    }
}

// ############################# PER-INDEX MEASUREMENT #############################

void measure_br_index(const bench_config_t& cfg)
{
    std::cerr << "loading the br-index from " << cfg.br_index_path << " ..." << std::endl;
    br_index_adapter::br_index_t br_index;
    std::ifstream in(cfg.br_index_path);

    if (!in.good()) {
        throw std::runtime_error("br-index file " + cfg.br_index_path + " not found");
    }

    load_measured([&]{ br_index.load(in); });
    in.close();

    br_index_adapter index(br_index);
    measure_all_jobs<br_index_adapter>(index, index.input_size(), cfg, "br_index");
}

// ############################# NATIVE COMPETITOR ALGORITHMS #############################
// In addition to the adapter measurements above (move-r's apm algorithm driven on each structure), b-move and
// br-index are measured running their OWN approximate-matching algorithms, so the benchmark covers every algorithm.

namespace {

// reads a whole pattern file (pizza&chili header + the fixed-length patterns)
inline bool read_patterns(const bench_job_t& job, std::vector<std::string>& patterns, uint64_t& m)
{
    std::ifstream pf;
    uint64_t num_patterns;
    if (!open_patterns(job, pf, num_patterns, m)) return false;
    if (num_patterns == 0) { std::cerr << "warning: no patterns in " << job.path << "; skipping" << std::endl; return false; }
    patterns.resize(num_patterns);
    for (uint64_t i = 0; i < num_patterns; i++) { no_init_resize(patterns[i], m); pf.read(patterns[i].data(), m); }
    if (!pf) { std::cerr << "warning: " << job.path << " is truncated; skipping" << std::endl; return false; }
    return true;
}

// replays a pattern set until at least cfg.min_time_seconds have been measured (>= one pass) and writes the
// averaged RESULT line; measure_one(i, time) measures pattern i, adds its time and returns its occurrence count.
// peak_mem is measured over the whole replay; occ_mem / cigar_mem (if provided) are read after the loop -- the
// caller's measure_one populates them as the peak single-pattern occurrence / CIGAR footprint.
template <typename measure_t>
void replay(const bench_config_t& cfg, const bench_job_t& job, uint64_t num_patterns, uint64_t m,
            uint64_t n, const std::string& algo, const std::string& dist_metr,
            const std::string& scheme_str, const std::string& time_key, measure_t measure_one,
            bool cigar = false, const uint64_t* occ_mem = nullptr, const uint64_t* cigar_mem = nullptr)
{
    const uint64_t min_time_ns = (uint64_t) (cfg.min_time_seconds * 1e9);
    uint64_t num_occurrences = 0, time = 0, passes = 0;
    const uint64_t mem_baseline = malloc_count_current();
    malloc_count_reset_peak();
    do {
        const uint64_t time_before = time;
        for (uint64_t i = 0; i < num_patterns; i++) num_occurrences += measure_one(i, time);
        passes++;
        if (time == time_before) break;
    } while (time < min_time_ns);
    mem_stats_t mem;
    mem.peak_mem = malloc_count_peak() - mem_baseline;
    mem.occ_mem = occ_mem ? *occ_mem : 0;
    mem.cigar_mem = cigar_mem ? *cigar_mem : 0;
    print_result(algo, cfg.text_name, n, dist_metr, scheme_str, m, num_patterns,
        job.k, num_occurrences / passes, time_key, time / passes, cigar, mem);
}

// maps a move-r search-scheme name to the closest columba built-in scheme name
inline std::string columba_scheme_name(const std::string& mr_scheme)
{
    if (mr_scheme == "min_u")       return "minU";
    if (mr_scheme == "pigeon_hole") return "pigeon";
    if (mr_scheme == "01")          return "01*0";
    return "minU"; // suffix_filter / file schemes have no columba equivalent -> closest comparable
}

} // namespace

// measures one columba flavor's NATIVE approximate-locate algorithm (locate only; hamming & edit) through its
// plugin's C ABI (@p api). algo is locate_columba_rlc (run-length-compressed = b-move) or locate_columba
// (FM-index). Under --cigar every locate job is additionally measured with CIGAR generation (columba computes
// real CIGARs -- from the matched string in RLC, from the located text in FM).
void measure_columba_api(const columba_api_t& api, const bench_config_t& cfg, const std::string& base)
{
    const std::string scheme_str = columba_scheme_name(cfg.scheme);
    const std::string flavor = api.flavor;
    std::cerr << "loading the " << flavor << " (native) index from " << base << " ..." << std::endl;

    // load inside load_measured so the heap growth of the load is captured as index_mem
    void* handle = nullptr;
    load_measured([&]{ handle = api.load(base.c_str(), scheme_str.c_str()); });
    if (handle == nullptr) {
        std::cerr << "warning: skipping " << flavor << " (index files not found at " << base << ")" << std::endl;
        return;
    }
    const uint64_t n = api.input_size(handle);

    const bool do_ham  = cfg.metric == "all" || cfg.metric == "hamming";
    const bool do_edit = cfg.metric == "all" || cfg.metric == "edit";

    for (const bench_job_t& job : cfg.jobs) {
        if (job.metric == HAMMING_DISTANCE && !do_ham) continue;
        if (job.metric == EDIT_DISTANCE && !do_edit) continue;

        if (job.op == OP_COUNT) {
            if (!cfg.run_non_native()) continue; // columba's only count here is move-rb's (non-native) apm count
            if (job.metric != HAMMING_DISTANCE) continue; // edit-distance count is unsupported
            if (job.k > api.max_k_hamming) {
                std::cerr << "warning: k=" << job.k << " > " << api.max_k_hamming << " unsupported by columba; skipping "
                          << job.path << std::endl;
                continue;
            }
            std::vector<std::string> patterns; uint64_t m;
            std::cerr << "measuring " << flavor << " (apm count) on " << job.path << " ..." << std::endl;
            if (!read_patterns(job, patterns, m)) continue;
            const uint64_t num_patterns = patterns.size();
            replay(cfg, job, num_patterns, m, n, "count_" + flavor + "_apm", "hamming",
                cfg.scheme, "time_count", [&](uint64_t i, uint64_t& time) {
                    auto t1 = now();
                    uint64_t c = api.count_apm(handle, patterns[i].data(), m, (int) job.k);
                    time += time_diff_ns(t1);
                    return c;
                });
            continue;
        }

        if (job.op != OP_LOCATE) continue; // (any other op is not measured for columba)
        const bool is_edit = job.metric == EDIT_DISTANCE;
        const int k_limit = is_edit ? api.max_k_edit : api.max_k_hamming;
        if (job.k > k_limit) {
            std::cerr << "warning: k=" << job.k << " > " << k_limit << " unsupported by columba; skipping "
                      << job.path << std::endl;
            continue;
        }

        const bool want_native = cfg.run_native() && (cfg.locate_plain() || cfg.locate_cigar());
        const bool want_apm    = cfg.run_non_native() && cfg.locate_plain();
        if (!want_native && !want_apm) continue;

        std::vector<std::string> patterns; uint64_t m;
        if (!read_patterns(job, patterns, m)) continue;
        const uint64_t num_patterns = patterns.size();
        const int metric_edit = is_edit ? 1 : 0;
        const std::string dist_str = is_edit ? "edit" : "hamming";

        const std::string job_scheme = job.k <= 7 ? scheme_str : "columba";
        void* bundles = nullptr;
        if (want_native) {
            api.set_scheme(handle, job_scheme.c_str());
            std::cerr << "measuring " << flavor << " (native, " << job_scheme << ") on " << job.path << " ..." << std::endl;

            // precompute the read bundles inside the plugin (kept out of the per-pattern locate time)
            std::vector<const char*> ptrs(num_patterns);
            std::vector<uint64_t> lens(num_patterns);
            for (uint64_t i = 0; i < num_patterns; i++) { ptrs[i] = patterns[i].data(); lens[i] = patterns[i].size(); }
            bundles = api.make_bundles(handle, ptrs.data(), lens.data(), num_patterns);

            // plain locate (no CIGAR); occ_mem = peak single-pattern occurrence-vector footprint (reported by the plugin)
            if (cfg.locate_plain()) {
                uint64_t occ_mem = 0;
                replay(cfg, job, num_patterns, m, n, "locate_" + flavor, dist_str,
                    job_scheme, "time_locate", [&](uint64_t i, uint64_t& time) {
                        uint64_t ob = 0, cb = 0;
                        auto t1 = now();
                        uint64_t c = api.locate(handle, bundles, i, metric_edit, (int) job.k, /* cigar */ 0, &ob, &cb);
                        time += time_diff_ns(t1);
                        occ_mem = std::max(occ_mem, ob);
                        return c;
                    }, /* cigar */ false, &occ_mem);
            }

            // additionally, CIGAR-producing locate (cigar=1) -- columba computes a real alignment per occurrence
            if (cfg.locate_cigar()) {
                uint64_t occ_mem_c = 0, cigar_mem = 0;
                replay(cfg, job, num_patterns, m, n, "locate_" + flavor, dist_str,
                    job_scheme, "time_locate", [&](uint64_t i, uint64_t& time) {
                        uint64_t ob = 0, cb = 0;
                        auto t1 = now();
                        uint64_t c = api.locate(handle, bundles, i, metric_edit, (int) job.k, /* cigar */ 1, &ob, &cb);
                        time += time_diff_ns(t1);
                        occ_mem_c = std::max(occ_mem_c, ob);
                        cigar_mem = std::max(cigar_mem, cb);
                        return c;
                    }, /* cigar */ true, &occ_mem_c, &cigar_mem);
            }
        }

        // move-rb's OWN apm algorithm run on the columba index (same algorithm, different index), reported as
        // locate_<flavor>_apm and deduplicated like move-rb. It has no CIGAR variant (hence skipped under --cigar only)
        if (want_apm) {
            std::cerr << "measuring " << flavor << " (apm) on " << job.path << " ..." << std::endl;
            uint64_t occ_mem_apm = 0;
            replay(cfg, job, num_patterns, m, n, "locate_" + flavor + "_apm", dist_str,
                cfg.scheme, "time_locate", [&](uint64_t i, uint64_t& time) {
                    uint64_t ob = 0;
                    auto t1 = now();
                    uint64_t c = api.locate_apm(handle, patterns[i].data(), m, metric_edit, (int) job.k, &ob);
                    time += time_diff_ns(t1);
                    occ_mem_apm = std::max(occ_mem_apm, ob);
                    return c;
                }, false, &occ_mem_apm);
        }

        if (bundles) api.free_bundles(bundles);
    }
    api.destroy(handle);
}

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

// dlopen()s one columba flavor plugin (RTLD_LOCAL keeps its columba symbols private, so both flavors coexist in
// one process) and measures it through the api getter @p api_symbol
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
    measure_columba_api(*getter(), cfg, base);
    // the plugin is intentionally left mapped (not dlclose'd) until process exit
}

// reads the index-position width (sizeof(length_t), in bytes: 4 for 32-bit, 8 for 64-bit) that the columba index at
// base name @p base was built with, from the second line of its <base>.meta (line 1 = build tag, line 3 = flavor).
// Returns 0 if the meta file is missing or unreadable.
int columba_meta_width_bytes(const std::string& base)
{
    std::ifstream ifs(base + ".meta");
    long tag = 0, width_bytes = 0;
    if (ifs && (ifs >> tag >> width_bytes)) return (int) width_bytes;
    return 0;
}

// measures one columba flavor. Its length_t width is fixed at compile time and must match the provided index (columba's
// loader throws otherwise), so the 32- or 64-bit plugin (lib<plugin_base>_{32,64}.so) is picked from the width recorded
// in the index's meta file. @p base is the flavor's index base name (given explicitly on the command line via
// --columba / --columba-rlc). A 32-bit index caps inputs at ~4 GB but is smaller/faster.
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

// measures the given columba flavors (run-length-compressed = b-move, and the FM-index), each via its plugin,
// from its own explicitly given base name
void measure_columba_native(const bench_config_t& cfg)
{
    if (!cfg.columba_rlc_path.empty())
        measure_columba_flavor(cfg, cfg.columba_rlc_path, "columba_rlc_plugin", "columba_rlc_api");
    if (!cfg.columba_path.empty())
        measure_columba_flavor(cfg, cfg.columba_path, "columba_fm_plugin", "columba_fm_api");
}

// measures br-index's NATIVE approximate count & locate algorithm (hamming distance only)
void measure_br_index_native(const bench_config_t& cfg)
{
    std::cerr << "loading the br-index (native) from " << cfg.br_index_path << " ..." << std::endl;
    std::optional<br_index_native> index_opt;
    load_measured([&]{ index_opt.emplace(cfg.br_index_path); });
    br_index_native& index = *index_opt;
    const uint64_t n = index.bwt_size();

    for (const bench_job_t& job : cfg.jobs) {
        if (job.metric != HAMMING_DISTANCE) continue; // br-index: hamming distance only
        std::vector<std::string> patterns; uint64_t m;
        std::cerr << "measuring br-index (native) on " << job.path << " ..." << std::endl;
        if (!read_patterns(job, patterns, m)) continue;
        const uint64_t num_patterns = patterns.size();
        if (job.op == OP_COUNT)
            replay(cfg, job, num_patterns, m, n, "count_br_index_native", "hamming", "none", "time_count",
                [&](uint64_t i, uint64_t& time) {
                    auto t1 = now(); uint64_t c = index.count_hamming_dist(patterns[i], job.k); time += time_diff_ns(t1); return c;
                });
        else if (cfg.locate_plain())
            replay(cfg, job, num_patterns, m, n, "locate_br_index_native", "hamming", "none", "time_locate",
                [&](uint64_t i, uint64_t& time) {
                    auto t1 = now(); uint64_t c = index.locate_hamming_dist(patterns[i], job.k); time += time_diff_ns(t1); return c;
                });
    }
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
        throw std::runtime_error("move_rb index file " + path + " not found");
    }
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
    std::cout << "move-rb-bench-apm: benchmarks approximate count- and locate-performance of" << std::endl;
    std::cout << "               b-move, br-index and move_rb (move & rlzsa). The indexes must be" << std::endl;
    std::cout << "               built beforehand (each with its own build tool); each is measured only if" << std::endl;
    std::cout << "               its path is given below. The pattern sets are generated by move-rb-gen-apm-queries." << std::endl << std::endl;
    std::cout << "usage: move-rb-bench-apm [...] <index flags...> <text_name> <patterns_path>" << std::endl;
    std::cout << "   -s <scheme>          search scheme to use (default: min_u): one of" << std::endl;
    std::cout << "                        pigeon_hole, suffix_filter, min_u, 01, or a path to a search-scheme file." << std::endl;
    std::cout << "                        For the built-in schemes the number of errors k is taken from each pattern" << std::endl;
    std::cout << "                        file name; a scheme file fixes k itself." << std::endl;
    std::cout << "   --metric <metric>    (optional) run only one distance metric: all|hamming|edit (default: all)" << std::endl;
    std::cout << "   --algo <selection>   (optional) which algorithm(s) to run: native|non-native|both (default: both)." << std::endl;
    std::cout << "                        native     = each index's own algorithm: move_rb (move/rlzsa) and the" << std::endl;
    std::cout << "                                     columba / b-move / br-index native locate." << std::endl;
    std::cout << "                        non-native = move-rb's apm algorithm on a foreign structure: the columba" << std::endl;
    std::cout << "                                     *_apm and br-index-adapter measurements." << std::endl;
    std::cout << "   --only <op>          (optional) measure only count or locate files: count, locate or both (default: both)" << std::endl;
    std::cout << "   --k <list>           (optional) measure only pattern files with these error counts k (comma-separated)" << std::endl;
    std::cout << "   --m <list>           (optional) measure only pattern files with these pattern lengths m (comma-separated)" << std::endl;
    std::cout << "   --cigar [mode]       (optional) which CIGAR variant of locate to measure: off|both|only" << std::endl;
    std::cout << "                        (bare --cigar = both; default: off). off = only the plain locate; both = the" << std::endl;
    std::cout << "                        plain locate and additionally a CIGAR-producing locate (RESULT field cigar=1);" << std::endl;
    std::cout << "                        only = just the CIGAR locate. Applies to locate only (count has no CIGAR) and" << std::endl;
    std::cout << "                        to the indexes that emit a CIGAR (move_rb and columba-native); the others are" << std::endl;
    std::cout << "                        skipped for locate under 'only'." << std::endl;
    std::cout << "   --time <s>           (optional) replay each pattern file's query set until at least <s> seconds" << std::endl;
    std::cout << "                        have been measured for it (>= one pass); the reported time_* is then the" << std::endl;
    std::cout << "                        average time for one pass over the set (num_patterns and num_occurrences stay" << std::endl;
    std::cout << "                        per-pass). 0 = a single pass; default: 10" << std::endl;
    std::cout << "   index flags          (all optional, all independent) the location of each index to measure;" << std::endl;
    std::cout << "                        an index is measured only if its flag is given. At least one is required:" << std::endl;
    std::cout << "     --move-rb <file>         move_rb index using move data structures" << std::endl;
    std::cout << "     --move-rb-rlzsa <file>   move_rb index using an rlzsa" << std::endl;
    std::cout << "     --bri <file>             br-index (.bri); drives both br_index and br_index_native" << std::endl;
    std::cout << "     --columba <base>         columba FM-index base name (columba is ACGT-only)" << std::endl;
    std::cout << "     --columba-rlc <base>     columba-rlc (b-move) index base name (ACGT-only)" << std::endl;
    std::cout << "   <text_name>          name of the original text (used in the output and pattern file names)" << std::endl;
    std::cout << "   <patterns_path>      directory containing the pattern files (from move-rb-gen-apm-queries) of the form" << std::endl;
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
    std::setlocale(LC_ALL, "C");
    bench_config_t cfg;
    int arg_idx = 1;

    // optional flags
    while (arg_idx < argc && argv[arg_idx][0] == '-') {
        std::string opt = argv[arg_idx++];
        if (opt == "-s") {
            if (arg_idx >= argc) help("error: missing scheme after -s");
            cfg.scheme = argv[arg_idx++];
        } else if (opt == "--metric") {
            if (arg_idx >= argc) help("error: missing metric after --metric");
            cfg.metric = argv[arg_idx++];
        } else if (opt == "--only") {
            if (arg_idx >= argc) help("error: missing operation after --only");
            std::string op = argv[arg_idx++];
            if      (op == "count")  { cfg.do_count = true;  cfg.do_locate = false; }
            else if (op == "locate") { cfg.do_count = false; cfg.do_locate = true; }
            else if (op == "both")   { cfg.do_count = true;  cfg.do_locate = true; }
            else help("error: --only must be count, locate or both");
        } else if (opt == "--time") {
            if (arg_idx >= argc) help("error: missing seconds after --time");
            cfg.min_time_seconds = std::stod(argv[arg_idx++]);
            if (cfg.min_time_seconds < 0) help("error: --time must be >= 0");
        } else if (opt == "--algo") {
            if (arg_idx >= argc) help("error: missing selection after --algo");
            std::string a = argv[arg_idx++];
            if      (a == "native")     cfg.algo = ALGO_NATIVE;
            else if (a == "non-native") cfg.algo = ALGO_NON_NATIVE;
            else if (a == "both")       cfg.algo = ALGO_BOTH;
            else help("error: --algo must be native, non-native or both");
        } else if (opt == "--cigar") {
            if (arg_idx < argc) {
                std::string mode = argv[arg_idx];
                if      (mode == "off")  { cfg.cigar_mode = CIGAR_OFF;  arg_idx++; }
                else if (mode == "both") { cfg.cigar_mode = CIGAR_BOTH; arg_idx++; }
                else if (mode == "only") { cfg.cigar_mode = CIGAR_ONLY; arg_idx++; }
                else cfg.cigar_mode = CIGAR_BOTH;
            } else {
                cfg.cigar_mode = CIGAR_BOTH;
            }
        } else if (opt == "--k") {
            if (arg_idx >= argc) help("error: missing k list after --k");
            for (const std::string& s : split_csv(argv[arg_idx++])) cfg.k_filter.push_back(std::stoll(s));
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

    resolve_scheme(cfg);

    if (!std::filesystem::is_directory(cfg.patterns_path)) {
        help("error: <patterns_path> is not a directory");
    }

    // collect all matching pattern files (restricted to the --k / --m values if given)
    auto in_filter = [](const auto& filter, auto value) {
        return filter.empty() || std::find(filter.begin(), filter.end(), value) != filter.end();
    };
    for (const auto& entry : std::filesystem::directory_iterator(cfg.patterns_path)) {
        if (!entry.is_regular_file()) continue;
        bench_job_t job;
        if (parse_pattern_file_name(entry.path().filename().string(), cfg.text_name, entry.path().string(), job)) {
            const bool op_ok = job.op == OP_COUNT ? cfg.do_count : cfg.do_locate;
            if (op_ok && in_filter(cfg.k_filter, job.k) && in_filter(cfg.m_filter, job.m)) cfg.jobs.emplace_back(job);
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
    // file) only reports that index as skipped instead of aborting the whole benchmark
    auto guarded = [](const char* name, void (*fnc)(const bench_config_t&), const bench_config_t& cfg) {
        try {
            fnc(cfg);
        } catch (const std::exception& e) {
            std::cerr << "warning: skipping " << name << " (" << e.what() << ")" << std::endl;
        }
    };

    // measure every index whose path was given; a failure to load one is reported and skipped (guarded). columba
    // (both flavors) is handled per-flavor inside measure_columba_native. --bri drives both the br-index adapter
    // and its native algorithm; columba is gated on either flavor's path
    if (!cfg.br_index_path.empty() && cfg.run_non_native())
        guarded("br-index (apm)", measure_br_index, cfg);
    if (!cfg.columba_rlc_path.empty() || !cfg.columba_path.empty())
        guarded("columba", measure_columba_native, cfg);
    if (!cfg.br_index_path.empty() && cfg.run_native())
        guarded("br-index (native)", measure_br_index_native, cfg);
    if (!cfg.move_rb_move_path.empty() && cfg.run_native())
        guarded("move_rb (move)", measure_move_rb_move, cfg);
    if (!cfg.move_rb_rlzsa_path.empty() && cfg.run_native())
        guarded("move_rb (rlzsa)", measure_move_rb_rlzsa, cfg);

    return 0;
}
