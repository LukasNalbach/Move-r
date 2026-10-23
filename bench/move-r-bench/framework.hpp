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

#ifndef MOVE_R_BENCH_MOVE_R_BENCH_FRAMEWORK_HPP
#define MOVE_R_BENCH_MOVE_R_BENCH_FRAMEWORK_HPP

#define MOVE_R_BENCH 1

#include <bit>
#include <limits>
#include <string_view>
#include <misc/files.hpp>
#include <misc/index_cli.hpp>
#include <misc/parse.hpp>
#include <misc/pattern_reader.hpp>
#include <move_r/move_r.hpp>
#include <omp.h>
#include <string>
#include <thread>
#include <vector>

// global state of the benchmark (used by framework.hpp, indexes.hpp and move-r-bench.cpp)
inline bool chars_remapped = false; // true <=> the symbols of the input have been remapped (see preprocess_input)
// the remapping of the symbols: [0..255] map_char[c] is the symbol c is replaced with (only if chars_remapped)
inline std::vector<uint8_t> map_char;
inline uint64_t max_num_threads; // number of threads for preprocessing the input
inline std::string input; // the input text T (without a terminator)
inline uint64_t input_size; // n, the length of T plus one (for the terminator)
inline std::string path_mf; // path of the measurement file (-m)
inline bool check_correctness = false; // -c: check the reported occurrences
inline uint64_t external_peak_memory_usage; // peak memory in bytes of the external tools (see update_peak_memory_usage)
inline direct_ofstream mf; // the measurement file (-m)
inline direct_ifstream input_file;
inline direct_ifstream patterns_file_1;
inline direct_ifstream patterns_file_2;
inline std::string path_input_file;
inline std::string path_patterns_file_1;
inline std::string path_patterns_file_2;
inline std::string text_name; // file name of the input
inline int32_t arg_idx = 1; // index of the next command-line argument

// the number of byte values
inline constexpr uint64_t num_byte_values = uint64_t { std::numeric_limits<uint8_t>::max() } + 1;
// the symbols below this value (and '\n') must not occur in the input of the compared indexes (see preprocess_input)
inline constexpr uint8_t min_valid_uchar = 3;
// the queries are run this many times without time measurement before the measured pass (measurement setup)
inline constexpr uint64_t num_warmup_passes = 2;

/**
 * @brief returns the largest peak memory usage in bytes in a log file of malloc_count (the values in KiB after
 *        ", peak")
 * @param log_file input stream of the log file
 */
inline uint64_t malloc_count_peak_memory_usage(std::istream& log_file)
{
    std::string log_file_content;
    log_file.seekg(0, std::ios::end);
    resize_no_init(log_file_content, log_file.tellg());
    log_file.seekg(0, std::ios::beg);
    log_file.read((char*) &log_file_content[0], log_file_content.size());
    int64_t size = (int64_t) log_file_content.size();
    int64_t pos = 0;
    double cur_peak_kib = 0.0;

    auto is_digit = [&](int64_t i) { return '0' <= log_file_content[i] && log_file_content[i] <= '9'; };

    constexpr std::string_view peak_marker = ", peak";

    while ((pos = log_file_content.find(peak_marker, pos)) != (int64_t) std::string::npos) {
        pos += peak_marker.size();
        while (pos < size && !is_digit(pos)) pos++;
        std::string str_cur_peak;
        while (pos < size && (is_digit(pos) || log_file_content[pos] == '.'))
            str_cur_peak.push_back(log_file_content[pos++]);
        double value;
        if (parse_double(str_cur_peak, value)) cur_peak_kib = std::max(cur_peak_kib, value);
    }

    return (uint64_t) (cur_peak_kib * 1024.0);
}

/**
 * @brief replaces the symbols of the input by min_valid_uchar, min_valid_uchar + 1, ... (without '\n') in their
 *        order, if the input contains a symbol below min_valid_uchar or '\n'
 */
void preprocess_input()
{
    std::vector<std::vector<uint8_t>> contains_uchar_thr(max_num_threads, std::vector<uint8_t>(num_byte_values, 0));

    #pragma omp parallel for num_threads(max_num_threads)
    for (uint64_t pos = 0; pos < input_size - 1; pos++) {
        contains_uchar_thr[omp_get_thread_num()][std::bit_cast<uint8_t>(input[pos])] = 1;
    }

    std::vector<uint8_t> contains_uchar(num_byte_values, 0);

    for (uint64_t cur_char = 0; cur_char < num_byte_values; cur_char++) {
        for (uint64_t cur_thread = 0; cur_thread < max_num_threads; cur_thread++) {
            if (contains_uchar_thr[cur_thread][cur_char] == 1) {
                contains_uchar[cur_char] = 1;
                break;
            }
        }
    }

    contains_uchar_thr.clear();
    contains_uchar_thr.shrink_to_fit();
    uint64_t alphabet_size = 1;

    for (uint64_t cur_char = 0; cur_char < num_byte_values; cur_char++) {
        if (contains_uchar[cur_char] == 1) {
            alphabet_size++;
        }
    }

    bool contains_invalid_char = false;

    for (uint8_t cur_char = 0; cur_char < min_valid_uchar; cur_char++) {
        if (contains_uchar[cur_char] == 1) {
            contains_invalid_char = true;
            break;
        }
    }

    // make sure the input does not contain '\n'
    if (contains_uchar[uint8_t('\n')] == 1) {
        contains_invalid_char = true;
    }

    if (contains_invalid_char) {
        // the valid symbols: from min_valid_uchar on, without '\n'
        constexpr uint64_t max_distinct_chars = num_byte_values - min_valid_uchar - 1;

        if (alphabet_size - 1 > max_distinct_chars) {
            std::cout << "Error: the input contains more than " << max_distinct_chars << " distinct characters"
                      << std::endl;
        }

        chars_remapped = true;

        map_char.resize(num_byte_values, 0);
        uint64_t next_valid_char = min_valid_uchar;

        for (uint64_t cur_char = 0; cur_char < num_byte_values; cur_char++) {
            if (contains_uchar[cur_char] == 1) {
                map_char[cur_char] = next_valid_char;
                next_valid_char++;

                // make sure the input does not contain '\n'
                if (next_valid_char == uint8_t('\n')) {
                    next_valid_char++;
                }
            }
        }

        #pragma omp parallel for num_threads(max_num_threads)
        for (uint64_t pos = 0; pos < input_size - 1; pos++) {
            input[pos] = std::bit_cast<char>(map_char[std::bit_cast<uint8_t>(input[pos])]);
        }
    }

    contains_uchar.clear();
    contains_uchar.shrink_to_fit();
}

// sets external_peak_memory_usage to the maximum of its value and the peak memory usage in a malloc_count log file
void update_peak_memory_usage(std::istream& log_file)
{
    external_peak_memory_usage = std::max(external_peak_memory_usage, malloc_count_peak_memory_usage(log_file));
}

// builds an index of the input (specialized for each index in indexes.hpp)
template <typename idx_t>
void build_index(idx_t& index, uint64_t num_threads);

// returns the number of occurrences of a pattern (specialized for each index in indexes.hpp)
template <typename idx_t>
uint64_t count_pattern(idx_t& index, std::string& pattern);

// appends the occurrences of a pattern to occurrences (specialized for each index in indexes.hpp)
template <typename idx_t>
void locate_pattern(idx_t& index, std::string& pattern, std::vector<uint64_t>& occurrences);

// measurements of the construction of an index
struct build_result {
    uint64_t num_threads;
    uint64_t time_build; // in ns
    uint64_t peak_memory_usage; // in bytes
    uint64_t size_index; // in bytes
};

// measurements of a set of queries
struct query_result {
    uint64_t num_queries;
    uint64_t pattern_length;
    uint64_t num_occurrences; // of all patterns
    uint64_t time_query; // total time in ns
};

// replaces the symbols of a string as preprocess_input replaced the symbols of the input
void map_string(std::string& str)
{
    for (uint64_t pos = 0; pos < str.size(); pos++) {
        str[pos] = std::bit_cast<char>(map_char[std::bit_cast<uint8_t>(str[pos])]);
    }
}

// counts the patterns of the first pattern file: num_warmup_passes passes without and one pass with time measurement
template <typename idx_t>
query_result count_patterns(idx_t& index)
{
    patterns_file_1.seekg(0);
    auto [num_queries, pattern_length] = read_patterns_header(patterns_file_1);
    uint64_t num_occurrences = 0;
    uint64_t time_query = 0;
    std::string pattern;
    resize_no_init(pattern, pattern_length);
    std::chrono::steady_clock::time_point t1, t2;

    for (uint64_t pass = 0; pass < num_warmup_passes; pass++) {
        for (uint64_t cur_query = 0; cur_query < num_queries; cur_query++) {
            patterns_file_1.read((char*) &pattern[0], pattern_length);
            if (chars_remapped) map_string(pattern);
            count_pattern<idx_t>(index, pattern);
        }

        patterns_file_1.seekg(0);
        read_patterns_header(patterns_file_1); // skip the header
    }

    for (uint64_t cur_query = 0; cur_query < num_queries; cur_query++) {
        patterns_file_1.read((char*) &pattern[0], pattern_length);
        if (chars_remapped) map_string(pattern);
        t1 = now();
        num_occurrences += count_pattern<idx_t>(index, pattern);
        t2 = now();
        time_query += time_diff_ns(t1, t2);
    }

    return query_result { num_queries, pattern_length, num_occurrences, time_query };
}

// locates the patterns of a pattern file: num_warmup_passes passes without and one pass with time measurement (with
// -c, the occurrences are checked against the input and their number against count)
template <typename idx_t>
query_result locate_patterns(idx_t& index, std::istream& patterns_file)
{
    patterns_file.seekg(0);
    auto [num_queries, pattern_length] = read_patterns_header(patterns_file);
    uint64_t num_occurrences = 0;
    std::vector<uint64_t> occurrences;
    uint64_t time_query = 0;
    std::string pattern;
    resize_no_init(pattern, pattern_length);
    std::chrono::steady_clock::time_point t1, t2;
    bool correct = true;

    for (uint64_t pass = 0; pass < num_warmup_passes; pass++) {
        for (uint64_t cur_query = 0; cur_query < num_queries; cur_query++) {
            patterns_file.read((char*) &pattern[0], pattern_length);
            if (chars_remapped) map_string(pattern);
            locate_pattern<idx_t>(index, pattern, occurrences);
            occurrences.clear();
        }

        patterns_file.seekg(0);
        read_patterns_header(patterns_file); // skip the header
    }

    for (uint64_t cur_query = 0; cur_query < num_queries; cur_query++) {
        patterns_file.read((char*) &pattern[0], pattern_length);
        if (chars_remapped) map_string(pattern);
        t1 = now();
        locate_pattern<idx_t>(index, pattern, occurrences);
        t2 = now();
        time_query += time_diff_ns(t1, t2);
        num_occurrences += occurrences.size();

        if (check_correctness) {
            for (uint64_t occurrence : occurrences) {
                for (uint64_t pos = 0; pos < pattern_length; pos++) {
                    if (input[occurrence + pos] != pattern[pos]) {
                        correct = false;
                        break;
                    }
                }

                if (!correct) break;
            }

            if (!correct || occurrences.size() != count_pattern<idx_t>(index, pattern)) {
                correct = false;
                break;
            }
        }

        occurrences.clear();
    }

    if (check_correctness) {
        if (correct) std::cout << " (no wrong occurrences)";
        else std::cout << " (wrong occurrences)";
    }

    return query_result { num_queries, pattern_length, num_occurrences, time_query };
}

// writes the RESULT lines of an index (construction, count, locate on both pattern files)
void write_measurement_data(
    bool bench_count,
    bool bench_locate,
    std::string index_log_name,
    uint64_t size_index,
    build_result result_build = build_result {},
    query_result result_count = query_result {},
    query_result result_locate_1 = query_result {},
    query_result result_locate_2 = query_result {})
{
    {
        result_line line(mf);
        line("type", "comparison_build")
            ("implementation", index_log_name)
            ("text", text_name)
            ("text_length", input_size)
            ("num_threads", result_build.num_threads)
            ("time_construction", result_build.time_build)
            ("peak_memory_usage", result_build.peak_memory_usage)
            ("size_index", size_index);
    }

    if (bench_count) {
        result_line line(mf);
        line("type", "comparison_count")
            ("implementation", index_log_name)
            ("text", text_name)
            ("text_length", input_size)
            ("num_queries", result_count.num_queries)
            ("pattern_length", result_count.pattern_length)
            ("num_occurrences", result_count.num_occurrences)
            ("time_query", result_count.time_query)
            ("size_index", size_index);
    }

    if (bench_locate) {
        std::vector<query_result> locate_results = { result_locate_1, result_locate_2 };

        for (query_result& res : locate_results) {
            result_line line(mf);
            line("type", "comparison_locate")
                ("implementation", index_log_name)
                ("text", text_name)
                ("text_length", input_size)
                ("num_queries", res.num_queries)
                ("pattern_length", res.pattern_length)
                ("num_occurrences", res.num_occurrences)
                ("time_query", res.time_query)
                ("size_index", result_build.size_index);
        }
    }
}

/**
 * @brief builds an index of the input (one thread), measures count and locate and writes the RESULT lines
 * @tparam bench_count whether count is measured
 * @tparam bench_locate whether locate is measured
 * @param index_name name of the index in the output
 * @param index_log_name name of the index in the RESULT lines
 */
template <typename idx_t, bool bench_count, bool bench_locate>
void bench_index(std::string index_name, std::string index_log_name)
{
    idx_t index;
    std::chrono::steady_clock::time_point t1, t2;
    uint64_t m1, m2;
    build_result result_build;
    query_result result_count, result_locate_1, result_locate_2;
    external_peak_memory_usage = 0;

    std::cout << "############## benchmarking " << index_name << " ##############" << std::endl << std::endl;
    std::this_thread::sleep_for(std::chrono::seconds(1));
    std::cout << "building " << index_name << std::flush;

    malloc_count_reset_peak();
    m1 = malloc_count_current();
    t1 = now();
    build_index<idx_t>(index, 1);
    t2 = now();
    m2 = malloc_count_current();

    result_build = build_result { 1, time_diff_ns(t1, t2), std::max({ malloc_count_peak() - m1, external_peak_memory_usage, m2 - m1 }), m2 - m1 };

    std::cout << std::endl;
    std::cout << "build time: " << format_time(result_build.time_build) << std::endl;
    std::cout << "peak memory usage: " << format_size(result_build.peak_memory_usage) << std::endl;
    std::cout << "index size: " << format_size(result_build.size_index) << std::endl << std::endl;

    std::this_thread::sleep_for(std::chrono::seconds(1));

    if constexpr (bench_count) {
        std::cout << "counting the first set of patterns" << std::flush;
        result_count = count_patterns<idx_t>(index);
        if (!check_correctness) std::cout << ": " << format_query_throughput(result_count.num_queries, result_count.time_query);
        std::cout << std::endl << "total number of occurrences: " << result_count.num_occurrences << std::endl;
        std::this_thread::sleep_for(std::chrono::seconds(1));
    }

    if constexpr (bench_locate) {
        std::cout << "locating the first set of patterns" << std::flush;
        result_locate_1 = locate_patterns<idx_t>(index, patterns_file_1);
        if (!check_correctness) std::cout << ": " << format_query_throughput(result_locate_1.num_queries, result_locate_1.time_query);
        std::cout << std::endl << "total number of occurrences: " << result_locate_1.num_occurrences << std::endl;

        std::this_thread::sleep_for(std::chrono::seconds(1));

        std::cout << "locating the second set of patterns" << std::flush;
        result_locate_2 = locate_patterns<idx_t>(index, patterns_file_2);
        if (!check_correctness) std::cout << ": " << format_query_throughput(result_locate_2.num_queries, result_locate_2.time_query);
        std::cout << std::endl << "total number of occurrences: " << result_locate_2.num_occurrences << std::endl;
    }

    if constexpr (bench_count || bench_locate) std::cout << std::endl;

    if (mf.is_open()) {
        write_measurement_data(
            bench_count,
            bench_locate,
            index_log_name,
            result_build.size_index,
            result_build,
            result_count,
            result_locate_1,
            result_locate_2);
    }
}

#endif // MOVE_R_BENCH_MOVE_R_BENCH_FRAMEWORK_HPP
