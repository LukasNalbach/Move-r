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

// move-rb-bench-ext: measures the exact bidirectional search of ext_query.hpp (count, and count followed by locate)
// on move_rb (move & rlzsa), both columba variants and the br-index. The count and locate pattern files are generated
// by move-rb-gen-ext-queries.

#include "bench_common.hpp"
#include "ext_query.hpp"
#include "indexes/br_index_adapter.hpp"
#include <misc/index_cli.hpp>
#include <misc/pattern_reader.hpp>
#include <move_rb/move_rb.hpp>

#include <algorithm>
#include <clocale>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <random>
#include <string>
#include <string_view>
#include <vector>

#include <malloc_count.h>

// ############################# CONFIGURATION #############################

// the operation of a pattern file
enum ext_op_t : uint8_t {
    OP_COUNT = 0, // the bidirectional search (its result is the number of occurrences)
    OP_LOCATE = 1 // the bidirectional search followed by reporting all occurrences
};

// a pattern file
struct bench_job_t {
    ext_op_t op; // the operation (from the file name)
    uint64_t m; // the pattern length (from the file name)
    std::string path; // path of the pattern file
};

// the configuration of the benchmark
struct bench_config_t : bench_options_t {
    std::vector<bench_job_t> jobs; // the pattern files to measure
};

// ############################# RESULT OUTPUT #############################

/**
 * @brief prints a result and writes a RESULT line (see write_result_line; dist_metr = exact, k = 0)
 * @param algo name of the measured algorithm and index
 * @param text name of the text
 * @param n length of the text
 * @param pattern_length the pattern length
 * @param num_patterns number of patterns
 * @param num_occurrences number of occurrences of all patterns
 * @param time_key time_count or time_locate
 * @param time_value time of one pass over the patterns in ns
 * @param mem the memory measurements
 */
inline void print_result(
    const std::string& algo, const std::string& text, uint64_t n,
    uint64_t pattern_length, uint64_t num_patterns, uint64_t num_occurrences,
    const std::string& time_key, uint64_t time_value, const mem_stats_t& mem = {})
{
    std::cout << algo << " | exact "
              << (time_key == "time_count" ? "count" : "locate")
              << " | m=" << pattern_length
              << ": " << format_time(
                    num_patterns == 0 ? 0 : time_value / num_patterns)
              << "/query"
              << " (" << num_occurrences << " occurrences)" << std::endl;

    write_result_line(algo, text, n, "exact", false, "none", pattern_length, num_patterns, 0, num_occurrences,
        mem, nullptr, 0, time_key, time_value);
}

// returns the extension plans of a pattern set (drawn from a generator with the seed ext_plan_seed)
inline std::vector<ext_plan_t> build_plans(uint64_t m, uint64_t num_patterns)
{
    std::vector<ext_plan_t> plans(num_patterns);
    std::mt19937_64 plan_rng(ext_plan_seed);
    for (uint64_t i = 0; i < num_patterns; i++) plans[i] = make_ext_plan(m, plan_rng);
    return plans;
}

// ############################# MEASUREMENT #############################

/**
 * @brief measures a pattern file on move_rb or the br-index adapter and writes the RESULT line
 * @tparam idx_t move_rb or br_index_adapter
 * @param index the index
 * @param n length of the text
 * @param cfg the configuration
 * @param job the pattern file
 * @param index_name name of the index in the RESULT line
 */
template <typename idx_t>
void run_ext(idx_t& index, uint64_t n, const bench_config_t& cfg, const bench_job_t& job, const std::string& index_name)
{
    std::vector<std::string> patterns;
    uint64_t m;
    if (!read_pattern_file(job.path, 0, patterns, m)) return;
    const uint64_t num_patterns = patterns.size();
    const std::vector<ext_plan_t> plans = build_plans(m, num_patterns);

    if (job.op == OP_COUNT) {
        replay_result_t result = replay_pattern_set(cfg.min_time_seconds, num_patterns, [&](uint64_t i, uint64_t& time) {
            auto t1 = now();
            uint64_t c = ext_count<idx_t>(index, patterns[i], plans[i]);
            time += time_diff_ns(t1);
            return c;
        });

        print_result("count_" + index_name, cfg.text_name, n, m, num_patterns,
            result.num_occurrences, "time_count", result.time_ns);
    } else {
        std::vector<uint64_t> occ;
        uint64_t peak_occ = 0; // largest number of occurrences of a pattern
        const uint64_t mem_baseline = malloc_count_current();
        malloc_count_reset_peak();

        replay_result_t result = replay_pattern_set(cfg.min_time_seconds, num_patterns, [&](uint64_t i, uint64_t& time) {
            occ.clear();
            auto t1 = now();
            ext_locate<idx_t>(index, patterns[i], plans[i], [&](uint64_t o) { occ.emplace_back(o); });
            time += time_diff_ns(t1);
            peak_occ = std::max<uint64_t>(peak_occ, occ.size());
            return occ.size();
        });

        mem_stats_t mem;
        mem.peak_mem = malloc_count_peak() - mem_baseline;
        mem.occ_mem = peak_occ * sizeof(uint64_t);
        print_result("locate_" + index_name, cfg.text_name, n, m, num_patterns,
            result.num_occurrences, "time_locate", result.time_ns, mem);
    }
}

// measures all pattern files on move_rb or the br-index adapter
template <typename idx_t>
void measure_all_jobs(idx_t& index, uint64_t n, const bench_config_t& cfg, const std::string& index_name)
{
    for (const bench_job_t& job : cfg.jobs) {
        std::cerr << "measuring " << index_name << " on " << job.path << " ..." << std::endl;
        run_ext<idx_t>(index, n, cfg, job, index_name);
    }
}

// loads the br-index and measures all pattern files
void measure_br_index(const bench_config_t& cfg)
{
    std::cerr << "loading the br-index from " << cfg.br_index_path << " ..." << std::endl;
    br_index_adapter::br_index_t br_index;
    load_index_measured(br_index, cfg.br_index_path, "br-index");
    br_index_adapter index(br_index);
    measure_all_jobs<br_index_adapter>(index, index.input_size(), cfg, "br_index");
}

/**
 * @brief loads a move_rb index and measures all pattern files
 * @tparam support _locate_move or _locate_rlzsa
 */
template <move_r_support support>
void measure_move_rb(const bench_config_t& cfg, const std::string& path, const std::string& index_name)
{
    move_rb<support> index;
    load_index_measured(index, path, "move_rb index");
    measure_all_jobs(index, index.forward_index().input_size(), cfg, index_name);
}

void measure_move_rb_move(const bench_config_t& cfg)
{
    std::cerr << "loading the move_rb (move) index from " << cfg.move_rb_move_path << " ..." << std::endl;
    measure_move_rb<_locate_move>(cfg, cfg.move_rb_move_path, "move_rb_move");
}

void measure_move_rb_rlzsa(const bench_config_t& cfg)
{
    std::cerr << "loading the move_rb (rlzsa) index from " << cfg.move_rb_rlzsa_path << " ..." << std::endl;
    measure_move_rb<_locate_rlzsa>(cfg, cfg.move_rb_rlzsa_path, "move_rb_rlzsa");
}

/**
 * @brief measures all pattern files on a columba index through the C interface of its plugin; the extension plans
 *        are passed as bytes (0 = LEFT, 1 = RIGHT)
 * @param api the C interface of the plugin
 * @param cfg the configuration
 * @param base base name of the columba index
 */
void measure_columba_ext_api(const columba_api_t& api, const bench_config_t& cfg, const std::string& base)
{
    const std::string flavor = api.flavor;
    std::cerr << "loading the " << flavor << " (native) index from " << base << " ..." << std::endl;

    // the search scheme does not matter for exact matching
    void* handle = nullptr;
    auto load_time = now();
    log_phase_start(true, load_time, "loading the " + flavor + " index");
    load_measured([&] { handle = api.load(base.c_str(), "columba"); });
    log_phase_end(true, load_time);

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
        if (!read_pattern_file(job.path, 0, patterns, m)) continue;
        const uint64_t num_patterns = patterns.size();
        const std::vector<ext_plan_t> plans = build_plans(m, num_patterns);

        // the extension plans as bytes: [0..num_patterns-1] orders[i][j] = 0 (LEFT) or 1 (RIGHT) for plans[i].order[j]
        std::vector<std::vector<uint8_t>> orders(num_patterns);

        for (uint64_t i = 0; i < num_patterns; i++) {
            orders[i].resize(plans[i].order.size());
            for (uint64_t j = 0; j < plans[i].order.size(); j++) orders[i][j] = plans[i].order[j] == LEFT ? 0 : 1;
        }

        const int32_t with_locate = job.op == OP_LOCATE ? 1 : 0;
        uint64_t occ_mem = 0;
        const uint64_t mem_baseline = malloc_count_current();
        malloc_count_reset_peak();

        replay_result_t result = replay_pattern_set(cfg.min_time_seconds, num_patterns, [&](uint64_t i, uint64_t& time) {
            uint64_t occ_bytes = 0;
            auto t1 = now();
            uint64_t c = api.ext_run(ext, patterns[i].data(), m, plans[i].start, orders[i].data(), with_locate, &occ_bytes);
            time += time_diff_ns(t1);
            occ_mem = std::max(occ_mem, occ_bytes);
            return c;
        });

        if (job.op == OP_COUNT) {
            print_result("count_" + flavor, cfg.text_name, n, m, num_patterns,
                result.num_occurrences, "time_count", result.time_ns);
        } else {
            mem_stats_t mem;
            mem.peak_mem = malloc_count_peak() - mem_baseline;
            mem.occ_mem = occ_mem;
            print_result("locate_" + flavor, cfg.text_name, n, m, num_patterns,
                result.num_occurrences, "time_locate", result.time_ns, mem);
        }
    }

    api.ext_free(ext);
    api.destroy(handle);
}

// measures the given columba indexes (b-move and FM-index), each through its plugin
void measure_columba_native(const bench_config_t& cfg)
{
    if (!cfg.columba_rlc_path.empty()) {
        with_columba_plugin(cfg.columba_rlc_path, "columba-rlc-plugin", "columba_rlc_api",
            [&](const columba_api_t& api) { measure_columba_ext_api(api, cfg, cfg.columba_rlc_path); });
    }

    if (!cfg.columba_path.empty()) {
        with_columba_plugin(cfg.columba_path, "columba-fm-plugin", "columba_fm_api",
            [&](const columba_api_t& api) { measure_columba_ext_api(api, cfg, cfg.columba_path); });
    }
}

// ############################# DRIVER #############################

[[noreturn]] void help(const std::string& msg = "")
{
    if (!msg.empty()) std::cout << msg << std::endl;
    std::cout << "move-rb-bench-ext: measures the exact bidirectional search (count, and count followed by locate) of" << std::endl;
    std::cout << "               move_rb (move & rlzsa), both columba variants and the br-index. Each query starts" << std::endl;
    std::cout << "               at a random position of the pattern and adds the other symbols to the left and to" << std::endl;
    std::cout << "               the right in a random order. Each index is measured if its path is given. The" << std::endl;
    std::cout << "               pattern files are generated by move-rb-gen-ext-queries." << std::endl << std::endl;
    std::cout << "usage: move-rb-bench-ext [...] <index flags...> <text_name> <patterns_path>" << std::endl;
    std::cout << "   --only <op>          (optional) measure only count or locate files: count, locate or both (default: both)" << std::endl;
    std::cout << "   --m <list>           (optional) measure only pattern files with these pattern lengths m (comma-separated)" << std::endl;
    std::cout << "   -m <file>            append machine-readable RESULT lines to <file>" << std::endl;
    std::cout << "   --time <s>           (optional) repeat each pattern set until at least <s> seconds have been measured" << std::endl;
    std::cout << "                        (at least once); the times in the RESULT lines are for one pass over the set." << std::endl;
    std::cout << "                        0 = one pass; default: " << default_min_time_seconds << std::endl;
    std::cout << "   index flags          (optional) the indexes to measure; at least one is required:" << std::endl;
    std::cout << "     --move-rb <file>         move_rb index using move data structures" << std::endl;
    std::cout << "     --move-rb-rlzsa <file>   move_rb index using an rlzsa" << std::endl;
    std::cout << "     --bri <file>             br-index (.bri)" << std::endl;
    std::cout << "     --columba <base>         base name of a columba FM-index (columba supports only ACGT)" << std::endl;
    std::cout << "     --columba-rlc <base>     base name of a columba b-move index (only ACGT)" << std::endl;
    std::cout << "   <text_name>          name of the text (used in the RESULT lines and in the names of the pattern files)" << std::endl;
    std::cout << "   <patterns_path>      directory of the pattern files (from move-rb-gen-ext-queries), named" << std::endl;
    std::cout << "                        <text_name>.patterns-{count,locate}-m<value>" << std::endl;
    std::exit(0);
}

/**
 * @brief parses a pattern file name of the form <text_name>.patterns-{count,locate}-m<value>
 * @return whether name has this form (then job is set)
 */
bool parse_pattern_file_name(const std::string& name, const std::string& text_name, const std::string& path, bench_job_t& job)
{
    const std::string prefix = text_name + ".patterns-";
    if (name.compare(0, prefix.size(), prefix) != 0) return false;
    std::string rest = name.substr(prefix.size());

    // removes prefix from the front of rest, if rest starts with it
    auto consume = [&rest](std::string_view prefix) {
        if (!rest.starts_with(prefix)) return false;
        rest = rest.substr(prefix.size());
        return true;
    };

    if (consume("count-m")) job.op = OP_COUNT;
    else if (consume("locate-m")) job.op = OP_LOCATE;
    else return false;

    if (rest.empty()) return false;
    for (char c : rest) if (!std::isdigit((uint8_t) c)) return false;
    job.m = std::stoull(rest);
    job.path = path;
    return true;
}

int main(int argc, char** argv)
{
    std::setlocale(LC_ALL, "C");
    bench_config_t cfg;
    int32_t arg_idx = 1;

    while (arg_idx < argc && argv[arg_idx][0] == '-') {
        std::string opt = argv[arg_idx++];
        if (!cfg.parse_option(opt, argc, argv, arg_idx, help)) help("error: unrecognized option '" + opt + "'");
    }

    if (argc - arg_idx != 2) help();
    cfg.text_name = argv[arg_idx++];
    cfg.patterns_path = argv[arg_idx++];

    if (!cfg.has_index())
        help("error: no index given; pass at least one of --move-rb, --move-rb-rlzsa, --bri, --columba, --columba-rlc");

    if (!std::filesystem::is_directory(cfg.patterns_path)) help("error: <patterns_path> is not a directory");

    // collect the pattern files (only those with the --m values, if given)
    for (const auto& entry : std::filesystem::directory_iterator(cfg.patterns_path)) {
        if (!entry.is_regular_file()) continue;
        bench_job_t job;

        if (parse_pattern_file_name(entry.path().filename().string(), cfg.text_name, entry.path().string(), job)) {
            const bool op_ok = job.op == OP_COUNT ? cfg.do_count : cfg.do_locate;
            const bool m_ok = cfg.m_filter.empty() ||
                std::find(cfg.m_filter.begin(), cfg.m_filter.end(), job.m) != cfg.m_filter.end();
            if (op_ok && m_ok) cfg.jobs.emplace_back(job);
        }
    }

    if (cfg.jobs.empty())
        help("error: no pattern files of the form <text_name>.patterns-{count,locate}-m<value> found in <patterns_path>");

    // count before locate, then by pattern length
    std::sort(cfg.jobs.begin(), cfg.jobs.end(), [](const bench_job_t& a, const bench_job_t& b) {
        if (a.op != b.op) return a.op < b.op;
        return a.m < b.m;
    });

    std::cerr << "found " << cfg.jobs.size() << " pattern file(s) for text '" << cfg.text_name << "'" << std::endl;

    if (!cfg.br_index_path.empty()) run_guarded("br-index", measure_br_index, cfg);
    if (!cfg.columba_rlc_path.empty() || !cfg.columba_path.empty()) run_guarded("columba (native)", measure_columba_native, cfg);
    if (!cfg.move_rb_move_path.empty()) run_guarded("move_rb (move)", measure_move_rb_move, cfg);
    if (!cfg.move_rb_rlzsa_path.empty()) run_guarded("move_rb (rlzsa)", measure_move_rb_rlzsa, cfg);

    return 0;
}
