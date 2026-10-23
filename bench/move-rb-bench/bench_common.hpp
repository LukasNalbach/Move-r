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
 * @file bench_common.hpp
 * @brief parts shared by move-rb-bench-apm, move-rb-bench-ext, move-rb-gen-apm-queries and move-rb-gen-ext-queries:
 *        reading pattern files, repeating a pattern set until a minimum time is measured, RESULT lines, loading the
 *        indexes and the columba plugins, and the shared options
 */

#ifndef MOVE_R_BENCH_MOVE_RB_BENCH_BENCH_COMMON_HPP
#define MOVE_R_BENCH_MOVE_RB_BENCH_BENCH_COMMON_HPP

#include <climits>
#include <cstdint>
#include <dlfcn.h>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unistd.h>
#include <vector>

#include <malloc_count.h>
#include <misc/index_cli.hpp>
#include <misc/measurement.hpp>
#include <misc/parse.hpp>
#include <misc/pattern_reader.hpp>
#include <misc/search_schemes.hpp>

#include "indexes/columba_capi.hpp"

// ############################# OPTIONS #############################

// the default of --time of move-rb-bench-apm and move-rb-bench-ext: each pattern set is repeated until at least this
// many seconds have been measured (measurement setup)
inline constexpr double default_min_time_seconds = 10.0;

// the defaults of move-rb-gen-apm-queries and move-rb-gen-ext-queries (measurement setup): --time, --min, --timeout,
// --max-skips and -m
inline constexpr double default_target_seconds = 1.0;
inline constexpr uint64_t default_min_patterns = 1;
inline constexpr double default_timeout_seconds = 10.0;
inline constexpr uint64_t default_max_consecutive_skips = 10;
inline const std::vector<uint64_t> default_pattern_lengths = { 10, 20, 40, 80, 160, 320, 640, 1280 };

// returns the values of a list separated by commas
template <typename value_t>
inline std::string join_csv(const std::vector<value_t>& values)
{
    std::string out;
    for (const value_t& value : values) out += (out.empty() ? "" : ",") + std::to_string(value);
    return out;
}

/**
 * @brief splits a comma-separated list into its non-empty parts
 * @param list the list
 */
inline std::vector<std::string> split_csv(const std::string& list)
{
    std::vector<std::string> parts;

    for (uint64_t start = 0; start <= list.size();) {
        uint64_t comma = list.find(',', start);
        std::string part = list.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
        if (!part.empty()) parts.push_back(part);
        if (comma == std::string::npos) break;
        start = comma + 1;
    }

    return parts;
}

// returns whether name is the name of a built-in search scheme (see search_scheme_of)
inline bool is_search_scheme_name(const std::string& name)
{
    return name == "pigeon_hole" || name == "suffix_filter" || name == "min_u" || name == "01";
}

/**
 * @brief returns a built-in search scheme for k errors
 * @param name pigeon_hole, suffix_filter, 01 or min_u (any other name: min_u)
 * @param k the number of errors
 */
inline search_scheme_t search_scheme_of(const std::string& name, int64_t k)
{
    if (name == "pigeon_hole")   return pigeon_hole_scheme(k);
    if (name == "suffix_filter") return suffix_filter_scheme(k);
    if (name == "01")            return zero_one_scheme(k);
    return min_u_scheme(k);
}

// the measurement file (-m)
inline direct_ofstream g_measurement_file;

// the options shared by move-rb-bench-apm and move-rb-bench-ext
struct bench_options_t {
    std::string text_name; // name of the text (used in the RESULT lines and in the names of the pattern files)
    std::string patterns_path; // directory of the pattern files
    // the indexes; an index is measured if its path is given
    std::string move_rb_move_path; // --move-rb: move_rb index with move data structures (locate-move)
    std::string move_rb_rlzsa_path; // --move-rb-rlzsa: move_rb index with an rlzsa (locate-rlzsa)
    std::string br_index_path; // --bri: br-index (.bri)
    std::string columba_path; // --columba: base name of the columba FM-index
    std::string columba_rlc_path; // --columba-rlc: base name of the columba b-move index
    // --time: each pattern set is repeated until at least this many seconds have been measured (at least once)
    double min_time_seconds = default_min_time_seconds;
    bool do_count = true; // --only: measure the count pattern files
    bool do_locate = true; // --only: measure the locate pattern files
    std::vector<uint64_t> m_filter; // --m: if not empty, only the pattern files with these pattern lengths are measured

    // returns whether at least one index is given
    bool has_index() const
    {
        return !move_rb_move_path.empty() || !move_rb_rlzsa_path.empty() || !br_index_path.empty() ||
            !columba_path.empty() || !columba_rlc_path.empty();
    }

    /**
     * @brief parses an option if it is one of the shared options
     * @param opt the option
     * @param argc number of arguments
     * @param argv the arguments
     * @param arg_idx index of the argument after opt; it is advanced past the value of the option
     * @param help prints an error message and the usage information and exits
     * @return whether opt is a shared option
     */
    template <typename help_fnc_t>
    bool parse_option(const std::string& opt, int argc, char** argv, int32_t& arg_idx, help_fnc_t help)
    {
        // returns the value of the option
        auto value = [&](const char* what) -> std::string {
            if (arg_idx >= argc) help("error: missing " + std::string(what) + " after " + opt);
            return argv[arg_idx++];
        };

        if (opt == "--time") {
            if (!parse_double(value("seconds"), min_time_seconds) || min_time_seconds < 0)
                help("error: --time must be a number >= 0");
        } else if (opt == "--only") {
            std::string op = value("operation");
            if      (op == "count")  { do_count = true;  do_locate = false; }
            else if (op == "locate") { do_count = false; do_locate = true; }
            else if (op == "both")   { do_count = true;  do_locate = true; }
            else help("error: --only must be count, locate or both");
        } else if (opt == "--m") {
            for (const std::string& m : split_csv(value("m list"))) m_filter.push_back(std::stoull(m));
        } else if (opt == "-m") {
            if (!open_measurement_file(g_measurement_file, value("measurement file")))
                help("error: cannot open or create the measurement file");
        } else if (opt == "--move-rb") {
            move_rb_move_path = value("file");
        } else if (opt == "--move-rb-rlzsa") {
            move_rb_rlzsa_path = value("file");
        } else if (opt == "--bri") {
            br_index_path = value("file");
        } else if (opt == "--columba") {
            columba_path = value("base name");
        } else if (opt == "--columba-rlc") {
            columba_rlc_path = value("base name");
        } else {
            return false;
        }

        return true;
    }
};

// ############################# RESULT LINES #############################

// memory measured for one pattern set, in bytes (0 if not measured)
struct mem_stats_t {
    uint64_t peak_mem = 0; // peak of the additional heap memory while the pattern set is measured
    uint64_t occ_mem = 0; // the largest size of the occurrences of one pattern
    uint64_t cigar_mem = 0; // the largest size of the CIGAR strings of the occurrences of one pattern
};

// heap memory in bytes allocated while the index currently measured has been loaded (see load_measured)
inline uint64_t g_index_mem = 0;

// stores the heap memory allocated while load_fnc() runs in g_index_mem
template <typename load_fnc_t>
inline void load_measured(load_fnc_t load_fnc)
{
    const uint64_t before = malloc_count_current();
    load_fnc();
    const uint64_t after = malloc_count_current();
    g_index_mem = after > before ? after - before : 0;
}

/**
 * @brief writes a RESULT line to the measurement file (if -m is given)
 * @param algo name of the measured algorithm and index
 * @param text name of the text
 * @param n length of the text
 * @param dist_metr hamming, edit or exact
 * @param cigar whether CIGAR strings have been computed
 * @param search_scheme name of the search scheme
 * @param pattern_length the pattern length
 * @param num_patterns number of patterns
 * @param max_mismatches k, the maximum number of errors
 * @param num_occurrences number of occurrences (of all patterns, in one pass)
 * @param mem the memory measurements
 * @param extra_key key of an additional value (written before the time if not nullptr)
 * @param extra_value the additional value
 * @param time_key time_count or time_locate
 * @param time_value time of one pass over the patterns in ns
 */
inline void write_result_line(const std::string& algo, const std::string& text, uint64_t n,
    const std::string& dist_metr, bool cigar, const std::string& search_scheme, uint64_t pattern_length,
    uint64_t num_patterns, int64_t max_mismatches, uint64_t num_occurrences, const mem_stats_t& mem,
    const char* extra_key, uint64_t extra_value, const std::string& time_key, uint64_t time_value)
{
    if (!g_measurement_file.is_open()) return;

    result_line line(g_measurement_file);
    line("algo", algo)
        ("text", text)
        ("n", n)
        ("dist_metr", dist_metr)
        ("cigar", cigar ? 1 : 0)
        ("search_scheme", search_scheme)
        ("pattern_length", pattern_length)
        ("num_patterns", num_patterns)
        ("max_mismatches", max_mismatches)
        ("num_occurrences", num_occurrences)
        ("index_mem", g_index_mem)
        ("peak_mem", mem.peak_mem)
        ("occ_mem", mem.occ_mem)
        ("cigar_mem", mem.cigar_mem);
    if (extra_key != nullptr) line(extra_key, extra_value);
    line(time_key.c_str(), time_value);
}

// ############################# PATTERN SETS #############################

/**
 * @brief reads a pattern file (a pizza&chili header and the patterns, all of the same length) and removes it from
 *        the page cache; prints a warning if it cannot be read, is empty or is truncated
 * @param path path of the pattern file
 * @param max_patterns if not 0, only the first max_patterns patterns are returned
 * @param patterns is set to the patterns
 * @param m is set to the pattern length
 * @return whether the patterns have been read
 */
inline bool read_pattern_file(const std::string& path, uint64_t max_patterns, std::vector<std::string>& patterns,
    uint64_t& m)
{
    direct_ifstream file(path, std::ios::binary);

    if (!file.good()) {
        std::cerr << "warning: cannot open pattern file " << path << std::endl;
        return false;
    }

    uint64_t num_patterns;
    std::tie(num_patterns, m) = read_patterns_header(file);

    if (num_patterns == 0) {
        std::cerr << "warning: no patterns in " << path << "; skipping" << std::endl;
        return false;
    }

    patterns.resize(num_patterns);

    for (uint64_t i = 0; i < num_patterns; i++) {
        resize_no_init(patterns[i], m);
        file.read(patterns[i].data(), m);
    }

    const bool complete = bool(file);
    file.close();
    drop_from_page_cache(path, true);

    if (!complete) {
        std::cerr << "warning: " << path << " is truncated (expected " << num_patterns << " x " << m
                  << " bytes); skipping" << std::endl;
        return false;
    }

    if (max_patterns != 0 && patterns.size() > max_patterns) patterns.resize(max_patterns);
    return true;
}

// the result of repeating a pattern set (see replay_pattern_set)
struct replay_result_t {
    uint64_t num_occurrences = 0; // number of occurrences of all patterns in one pass
    uint64_t time_ns = 0; // average time of one pass in ns
};

/**
 * @brief measures the patterns of a set one after the other and repeats the set until at least min_time_seconds have
 *        been measured (at least once, and not again if a pass measures no time)
 * @param min_time_seconds the minimum time in seconds
 * @param num_patterns number of patterns
 * @param measure_one called as measure_one(i, time_ns) for the i-th pattern; it adds the measured time to time_ns and
 *        returns the number of occurrences
 */
template <typename measure_fnc_t>
replay_result_t replay_pattern_set(double min_time_seconds, uint64_t num_patterns, measure_fnc_t measure_one)
{
    const uint64_t min_time_ns = (uint64_t) (min_time_seconds * 1e9); // seconds to ns
    uint64_t num_occurrences = 0, time_ns = 0, passes = 0;

    do {
        const uint64_t time_before = time_ns;
        for (uint64_t i = 0; i < num_patterns; i++) num_occurrences += measure_one(i, time_ns);
        passes++;
        if (time_ns == time_before) break;
    } while (time_ns < min_time_ns);

    return { num_occurrences / passes, time_ns / passes };
}

// ############################# INDEXES #############################

/**
 * @brief loads an index (and stores the heap memory it takes in g_index_mem) and removes its file from the page cache
 * @param index the index
 * @param path path of the index file
 * @param name name of the index in the messages
 */
template <typename index_t>
void load_index_measured(index_t& index, const std::string& path, const std::string& name)
{
    direct_ifstream in(path);
    if (!in.good()) throw std::runtime_error(name + " file " + path + " not found");
    load_measured([&] { load_index_with_progress(index, in, name); });
    in.close();
    drop_from_page_cache(path, true);
}

// returns the directory of the running executable (the columba plugins are next to it)
inline std::string executable_dir()
{
    char buf[PATH_MAX];
    int64_t len = ::readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    if (len <= 0) return ".";
    buf[len] = '\0';
    std::string path(buf);
    uint64_t slash = path.find_last_of('/');
    return slash == std::string::npos ? std::string(".") : path.substr(0, slash);
}

/**
 * @brief returns the width in bytes (4 or 8) of the positions of a columba index (sizeof(length_t)), read from the
 *        second line of <base>.meta, or 0 if it cannot be read
 * @param base base name of the columba index
 */
inline int32_t columba_meta_width_bytes(const std::string& base)
{
    direct_ifstream ifs(base + ".meta");
    int64_t tag = 0, width_bytes = 0;
    if (ifs && (ifs >> tag >> width_bytes)) return (int32_t) width_bytes;
    return 0;
}

/**
 * @brief loads the columba plugin lib<plugin_base>-{32,64}.so that matches the width of the positions of a columba
 *        index (see columba_meta_width_bytes) and calls measure(api) with its C interface; the plugin is not unloaded
 * @param base base name of the columba index
 * @param plugin_base columba-rlc-plugin or columba-fm-plugin
 * @param api_symbol columba_rlc_api or columba_fm_api
 * @param measure called with the columba_api_t of the plugin
 */
template <typename measure_fnc_t>
void with_columba_plugin(const std::string& base, const std::string& plugin_base, const char* api_symbol,
    measure_fnc_t measure)
{
    const int32_t width = columba_meta_width_bytes(base);
    if (width == 0) std::cerr << "warning: cannot read " << base << ".meta; defaulting to the 64-bit plugin" << std::endl;
    const bool narrow = width == int32_t(sizeof(uint32_t));
    const std::string so_name = "lib" + plugin_base + (narrow ? "-32" : "-64") + ".so";
    std::cerr << "columba index at " << base << " is "
              << (width ? width * std::numeric_limits<uint8_t>::digits : std::numeric_limits<uint64_t>::digits)
              << "-bit; loading " << so_name << std::endl;

    // RTLD_LOCAL keeps the symbols of the plugin private, so that the plugins of both columba variants can be loaded
    const std::string so_path = executable_dir() + "/" + so_name;
    void* lib = dlopen(so_path.c_str(), RTLD_NOW | RTLD_LOCAL);

    if (lib == nullptr) {
        std::cerr << "warning: cannot load columba plugin " << so_path << " (" << dlerror() << ")" << std::endl;
        return;
    }

    using api_getter_t = const columba_api_t* (*)();
    api_getter_t getter = reinterpret_cast<api_getter_t>(dlsym(lib, api_symbol));

    if (getter == nullptr) {
        std::cerr << "warning: " << api_symbol << " not found in " << so_path << std::endl;
        dlclose(lib);
        return;
    }

    measure(*getter());
}

/**
 * @brief calls measure(cfg) and prints a warning instead of exiting if it throws (e.g., if an index file is missing)
 * @param name name of the measured index in the warning
 */
template <typename config_t>
void run_guarded(const char* name, void (*measure)(const config_t&), const config_t& cfg)
{
    try {
        measure(cfg);
    } catch (const std::exception& e) {
        std::cerr << "warning: skipping " << name << " (" << e.what() << ")" << std::endl;
    }
}

#endif // MOVE_R_BENCH_MOVE_RB_BENCH_BENCH_COMMON_HPP
