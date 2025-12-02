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
#include <move_r/move_rb.hpp>

static constexpr int min_args = 8;
int arg_idx = 1;
int64_t k = -1;
distance_metric_t dist_metr = NO_METRIC;
search_scheme_t search_scheme;
std::ofstream mf;
std::string path_index_file;
std::string path_patterns_file;
std::ifstream index_file;
std::ifstream patterns_file;
std::string name_text_file;

void help(std::string msg)
{
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "move-rb-count: count all approximate occurrences of the input patterns." << std::endl << std::endl;
    std::cout << "usage: move-rb-count [options] -k <mismatches> -d <metric> -s <scheme> <index_file> <patterns_file>" << std::endl;
    std::cout << "   -m <m_file> <text_name>    m_file is the file to write measurement data to," << std::endl;
    std::cout << "                              text_name should be the name of the original file" << std::endl;
    std::cout << "   <mismatches>               maximum number of allowed mismatches" << std::endl;
    std::cout << "   <metric>                   distance metric to use (hamming or edit)" << std::endl;
    std::cout << "   <scheme>                   search scheme to use (pigeon_hole, suffix_filter, 01 or path to a file)" << std::endl;
    std::cout << "   <index_file>               index file (with extension .move-r)" << std::endl;
    std::cout << "   <patterns_file>            file in move-rb-patterns format containing the patterns." << std::endl;
    exit(0);
}

void parse_args(char** argv, int argc)
{
    std::string s = argv[arg_idx];
    arg_idx++;

    if (s == "-m") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -o option.");
        std::string path_m_file = argv[arg_idx++];
        mf.open(path_m_file, std::filesystem::exists(path_m_file) ? std::ios::app : std::ios::out);
        if (!mf.good()) help("error: cannot open nor create <m_file>");
        name_text_file = argv[arg_idx++];
    } else {
        help("error: unrecognized '" + s + "' option");
    }
}

template <typename pos_t, move_r_support support>
void measure_count()
{
    std::cout << std::setprecision(4);
    std::cout << "loading the index" << std::flush;
    auto time = now();
    using idx_t = move_rb<support, char, pos_t>;
    idx_t index;
    index.load(index_file);
    index_file.close();
    time = log_runtime(time);
    std::cout << std::endl;
    index.log_data_structure_sizes();
    
    std::cout << std::endl << "searching patterns ... " << std::endl;
    std::string header;
    std::getline(patterns_file, header);
    uint64_t num_patterns = number_of_patterns(header);
    uint64_t pattern_length = patterns_length(header);
    uint64_t perc;
    uint64_t last_perc = 0;
    uint64_t num_occurrences = 0;
    uint64_t time_count = 0;
    std::chrono::steady_clock::time_point t2, t3;
    std::string pattern;
    no_init_resize(pattern, pattern_length);
    uint64_t baseline_alloc = malloc_count_current();
    malloc_count_reset_peak();

    for (uint64_t i = 0; i < num_patterns; i++) {
        perc = (100 * i) / num_patterns;

        if (perc > last_perc) {
            std::cout << perc << "% done .." << std::endl;
            last_perc = perc;
        }

        patterns_file.read(pattern.data(), pattern_length);
        t2 = now();
        num_occurrences += index.count_with_mismatches(pattern, search_scheme);
        t3 = now();
        time_count += time_diff_ns(t2, t3);
    }

    patterns_file.close();
    std::cout << "average occurrences per pattern: " << (num_occurrences / num_patterns) << std::endl;
    std::cout << "additional memory consumption during the search phase: " << format_size(malloc_count_peak() - baseline_alloc) << std::endl;
    std::cout << "number of patterns: " << num_patterns << std::endl;
    std::cout << "pattern length: " << pattern_length << std::endl;
    std::cout << "total number of occurrences: " << num_occurrences << std::endl;
    std::cout << "count time: " << format_time(time_count) << std::endl;
    std::cout << "            " << format_time(time_count / num_patterns) << "/pattern" << std::endl;
    std::cout << "            " << format_time(time_count / (num_patterns * pattern_length)) << "/character" << std::endl;

    if (mf.is_open()) {
        mf << "RESULT";
        mf << " algo=count_move_rb_" << move_r_support_suffix(support);
        mf << " text=" << name_text_file;
        mf << " a=" << index.forward_index().balancing_parameter();
        mf << " n=" << index.forward_index().input_size();
        mf << " sigma=" << std::to_string(index.forward_index().alphabet_size());
        mf << " r=" << index.forward_index().num_bwt_runs();
        mf << " r_rev=" << index.backward_index().num_bwt_runs();
        mf << " r_=" << index.forward_index().M_LF().num_intervals();
        mf << " r_rev_=" << index.backward_index().M_LF().num_intervals();

        if constexpr (idx_t::supports_multiple_locate) {
            if constexpr (idx_t::has_locate_move) {
                mf << " r__=" << index.forward_index().M_Phi_m1().num_intervals();
                mf << " r___=" << index.forward_index().M_Phi().num_intervals();
            } else if constexpr (idx_t::has_rlzsa) {
                mf << " z=" << index.forward_index().num_phrases_rlzsa();
                mf << " z_l=" << index.forward_index().num_literal_phrases_rlzsa();
                mf << " z_c=" << index.forward_index().num_copy_phrases_rlzsa();
            }
        }

        mf << " pattern_length=" << pattern_length;
        index.log_data_structure_sizes(mf);
        mf << " num_patterns=" << num_patterns;
        mf << " num_switches=" << num_occurrences;
        mf << " num_occurrences=" << num_occurrences;
        mf << " time_count=" << time_count;
        mf << std::endl;
        mf.close();
    }
}

int main(int argc, char** argv)
{
    if (argc - 1 < min_args) help("");
    while (arg_idx < argc - min_args) parse_args(argv, argc);

    std::string arg = argv[arg_idx++];
    if (arg != "-k") help("");
    k = atoi(argv[arg_idx++]);
    if (k < 0) help("error: invalid k value");

    arg = argv[arg_idx++];
    if (arg != "-d") help("");
    std::string str = argv[arg_idx++];
    if      (str == "hamming") dist_metr = HAMMING_DISTANCE;
    else if (str == "edit")    dist_metr = EDIT_DISTANCE;
    else help("error: invalid option after -d");

    arg = argv[arg_idx++];
    if (arg != "-s") help("");
    str = argv[arg_idx++];
    if      (str == "pigeon_hole")   search_scheme = pigeon_hole_scheme(k, dist_metr);
    else if (str == "suffix_filter") search_scheme = suffix_filter_scheme(k, dist_metr);
    else if (str == "01")            search_scheme = zero_one_scheme(k, dist_metr);
    else if (std::filesystem::exists(str)) {
        std::string file_content;
        uint64_t file_size = std::filesystem::file_size(str);
        no_init_resize(file_content, file_size);
        std::ifstream ifile(str);
        ifile.read(file_content.data(), file_size);
        search_scheme = parse_search_scheme(file_content, dist_metr);
    } else help("error: invalid option after -s");

    if (k != search_scheme.k_max) help("error: provided search scheme and k value are not compatible");

    path_index_file = argv[arg_idx];
    path_patterns_file = argv[arg_idx + 1];

    index_file.open(path_index_file);
    patterns_file.open(path_patterns_file);

    if (!index_file.good()) help("error: could not read <index_file>");
    if (!patterns_file.good()) help("error: could not read <patterns_file>");

    bool is_64_bit;
    index_file.read((char*) &is_64_bit, 1);
    move_r_support _support;
    index_file.read((char*) &_support, sizeof(move_r_support));
    index_file.seekg(0, std::ios::beg);

    if (_support == _count) {
        if (is_64_bit) {
            measure_count<uint64_t, _count>();
        } else {
            measure_count<uint32_t, _count>();
        }
    } else if (_support == _locate_move) {
        if (is_64_bit) {
            measure_count<uint64_t, _locate_move>();
        } else {
            measure_count<uint32_t, _locate_move>();
        }
    } else if (_support == _locate_rlzsa) {
        if (is_64_bit) {
            measure_count<uint64_t, _locate_rlzsa>();
        } else {
            measure_count<uint32_t, _locate_rlzsa>();
        }
    }
}