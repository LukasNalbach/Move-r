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

#include <filesystem>
#include <iostream>
#include <move_rb/move_rb.hpp>
#include <misc/progress.hpp>

static constexpr int min_args = 4;
int arg_idx = 1;
int64_t k = -1;
std::string scheme_str;
search_scheme_t search_scheme;
std::ofstream mf;
std::string path_index_file;
std::string path_patterns_file;
std::ifstream index_file;
std::ifstream patterns_file;
std::string name_text_file;
std::string path_output_file;
std::ofstream output_file;
bool output_occurrences;

/**
 * @brief prints the usage information and exits
 * @param msg an optional error message printed before the usage information
 */
void help(std::string msg)
{
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "move-rb-count: count all approximate occurrences of the input patterns w.r.t. hamming distance." << std::endl << std::endl;
    std::cout << "usage: move-rb-count [...] -s <scheme> [-k <mismatches>] <index_file> <patterns_file>" << std::endl;
    std::cout << "   -m <m_file> <text_name>    m_file is the file to write measurement data to," << std::endl;
    std::cout << "                              text_name should be the name of the original file" << std::endl;
    std::cout << "   <scheme>                   search scheme to use (pigeon_hole, suffix_filter, min_u, 01 or path to a file)" << std::endl;
    std::cout << "   <mismatches>               maximum number of allowed mismatches; applies only to" << std::endl;
    std::cout << "                              pigeon_hole, suffix_filter, min_u and 01 search schemes" << std::endl;
    std::cout << "   -o <output_file>           write pattern counts to this file (in ASCII format; one line per pattern)" << std::endl;
    std::cout << "   <index_file>               index file (with extension .move-r)" << std::endl;
    std::cout << "   <patterns_file>            file in move-rb-patterns format containing the patterns." << std::endl;
    exit(0);
}

/**
 * @brief parses the next command-line argument(s)
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments
 * @return whether there are further options to parse
 */
bool parse_args(char** argv, int argc)
{
    if (arg_idx >= argc - min_args) return false;
    std::string s = argv[arg_idx];
    if (s == "-s") return false;
    arg_idx++;

    if (s == "-m") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -o option.");
        std::string path_m_file = argv[arg_idx++];
        mf.open(path_m_file, std::filesystem::exists(path_m_file) ? std::ios::app : std::ios::out);
        if (!mf.good()) help("error: cannot open nor create <m_file>");
        name_text_file = argv[arg_idx++];
    } else if (s == "-o") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -o option.");
        output_occurrences = true;
        path_output_file = argv[arg_idx++];
    } else {
        help("error: unrecognized '" + s + "' option");
    }

    return true;
}

/**
 * @brief returns the next command-line argument, exiting with an error if there is none
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments
 * @param what a description of the expected argument (used in the error message)
 * @return the next argument
 */
const char* next_arg(int argc, char** argv, const std::string& what)
{
    if (arg_idx >= argc) help("error: missing " + what);
    return argv[arg_idx++];
}

/**
 * @brief loads the index and benchmarks counting the input patterns
 * @tparam pos_t index integer type
 * @tparam support the move-r locate-support type
 */
template <typename pos_t, move_r_support support>
void measure_count()
{
    std::cout << std::setprecision(4);
    auto time = now();
    log_phase_start(true, time, "loading the index");
    using idx_t = move_rb<support, char, pos_t>;
    idx_t index;
    index.load(index_file);
    index_file.close();
    log_phase_end(true, time);
    index.log_data_structure_sizes();
    
    std::cout << std::endl << "searching patterns ... " << std::endl;
    std::string header;
    std::getline(patterns_file, header);
    uint64_t num_patterns = number_of_patterns(header);
    uint64_t pattern_length = patterns_length(header);

    if (pattern_length < search_scheme.p)
        help("error: pattern length < number of parts defined in the search scheme");

    uint64_t num_occurrences = 0;
    uint64_t time_count = 0;
    uint64_t count = 0;
    std::chrono::steady_clock::time_point t2, t3;
    std::string pattern;
    no_init_resize(pattern, pattern_length);
    uint64_t baseline_alloc = malloc_count_current();
    malloc_count_reset_peak();
    progress_meter meter(num_patterns);

    for (uint64_t i = 0; i < num_patterns; i++) {
        patterns_file.read(pattern.data(), pattern_length);
        t2 = now();

        count = index.count_hamming_dist(pattern, search_scheme);
        num_occurrences += count;

        t3 = now();
        time_count += time_diff_ns(t2, t3);

        if (output_occurrences) {
            output_file << count << std::endl;
        }

        meter.step();
    }

    meter.finish();
    patterns_file.close();
    std::cout << "average occurrences per pattern: " << (num_occurrences / num_patterns) << std::endl;
    std::cout << "additional memory consumption during the search phase: " << format_size(malloc_count_peak() - baseline_alloc) << std::endl;
    std::cout << "number of patterns: " << num_patterns << std::endl;
    std::cout << "pattern length: " << pattern_length << std::endl;
    std::cout << "maximum number of mismatches: " << k << std::endl;
    std::cout << "search scheme: " << scheme_str << std::endl;
    std::cout << "total number of occurrences: " << num_occurrences << std::endl;
    std::cout << "count time: " << format_time(time_count) << std::endl;
    std::cout << "            " << format_time(time_count / num_patterns) << "/pattern" << std::endl;
    std::cout << "            " << format_time(time_count / (num_patterns * pattern_length)) << "/character" << std::endl;

    if (mf.is_open()) {
        mf << "RESULT";
        mf << " algo=count_move_rb_" << move_r_support_suffix(support);
        mf << " text=" << name_text_file;
        mf << " n=" << index.forward_index().input_size();
        mf << " pattern_length=" << pattern_length;
        mf << " num_patterns=" << num_patterns;
        mf << " max_mismatches=" << k;
        mf << " num_occurrences=" << num_occurrences;
        mf << " time_count=" << time_count;
        index.log_data_structure_sizes(mf);
        mf << std::endl;
        mf.close();
    }
}

/**
 * @brief program entry point
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments
 * @return the exit code
 */
int main(int argc, char** argv)
{
    if (argc - 1 < min_args) help("");
    while (parse_args(argv, argc));

    std::string arg = next_arg(argc, argv, "'-s <scheme>'");
    if (arg != "-s") help("error: expected '-s <scheme>'");
    scheme_str = next_arg(argc, argv, "<scheme>");
    bool is_default_scheme =
        scheme_str == "pigeon_hole" ||
        scheme_str == "suffix_filter" ||
        scheme_str == "min_u" ||
        scheme_str == "01";
    if (!is_default_scheme) {
        if (std::filesystem::exists(scheme_str)) {
            std::string file_content;
            uint64_t file_size = std::filesystem::file_size(scheme_str);
            no_init_resize(file_content, file_size);
            std::ifstream ifile(scheme_str);
            ifile.read(file_content.data(), file_size);
            search_scheme = parse_search_scheme(file_content);
        } else help("error: invalid option after -s");
    }
    
    if (is_default_scheme) {
        arg = next_arg(argc, argv, "'-k <mismatches>'");
        if (arg != "-k") help("error: expected '-k <mismatches>'");
        k = atoi(next_arg(argc, argv, "<mismatches>"));
        if (k < 0) help("error: invalid k value");
        if      (scheme_str == "pigeon_hole")   search_scheme = pigeon_hole_scheme(k);
        else if (scheme_str == "suffix_filter") search_scheme = suffix_filter_scheme(k);
        else if (scheme_str == "min_u")         search_scheme = min_u_scheme(k);
        else if (scheme_str == "01")            search_scheme = zero_one_scheme(k);
    } else {
        k = search_scheme.k;
    }

    if (arg_idx + 2 > argc) help("error: missing <index_file> and/or <patterns_file>");
    path_index_file = argv[arg_idx];
    path_patterns_file = argv[arg_idx + 1];

    index_file.open(path_index_file);
    patterns_file.open(path_patterns_file);

    if (!index_file.good()) help("error: could not read <index_file>");
    if (!patterns_file.good()) help("error: could not read <patterns_file>");

    if (output_occurrences) {
        output_file.open(path_output_file);
        if (!output_file.good()) help("error: could not create <output_file>");
    }

    bool is_64_bit;
    index_file.read((char*) &is_64_bit, 1);
    move_r_support _support;
    index_file.read((char*) &_support, sizeof(move_r_support));
    index_file.seekg(0, std::ios::beg);

    if (_support == _count) {
        if (is_64_bit) measure_count<uint64_t, _count>();
        else           measure_count<uint32_t, _count>();
    } else if (_support == _locate_move) {
        if (is_64_bit) measure_count<uint64_t, _locate_move>();
        else           measure_count<uint32_t, _locate_move>();
    } else if (_support == _locate_rlzsa) {
        if (is_64_bit) measure_count<uint64_t, _locate_rlzsa>();
        else           measure_count<uint32_t, _locate_rlzsa>();
    }
    
    return 0;
}