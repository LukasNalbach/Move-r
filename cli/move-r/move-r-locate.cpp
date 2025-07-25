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
#include <move_r/move_r.hpp>

int ptr = 1;
bool output_occurrences = false;
bool check_correctness = false;
std::string input;
std::ofstream mf;
std::string path_index_file;
std::string path_patterns_file;
std::string path_output_file;
std::ifstream index_file;
std::ifstream patterns_file;
std::ifstream input_file;
std::ofstream output_file;
std::string name_text_file;
std::string path_input_file;

void help(std::string msg)
{
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "move-r-locate: locate all occurrences of the input patterns." << std::endl << std::endl;
    std::cout << "usage: move-r-locate [options] <index_file> <patterns>" << std::endl;
    std::cout << "   -m <m_file> <text_name>    m_file is the file to write measurement data to," << std::endl;
    std::cout << "                              text_name should be the name of the original file" << std::endl;
    std::cout << "   -i <input_file>            input_file must be the file the index was built for" << std::endl;
    std::cout << "                              (required for locate_rlzsa_bin_search and the -c option)" << std::endl;
    std::cout << "   -c                         checks correctness of each pattern occurrence on <input_file>" << std::endl;
    std::cout << "   -o <output_file>           write pattern occurrences to this file (ASCII)" << std::endl;
    std::cout << "   <index_file>               index file (with extension .move-r)" << std::endl;
    std::cout << "   <patterns_file>            file in pizza&chili format containing the patterns" << std::endl;
    exit(0);
}

void parse_args(char** argv, int argc, int& ptr)
{
    std::string s = argv[ptr];
    ptr++;

    if (s == "-c") {
        check_correctness = true;
    } else if (s == "-m") {
        if (ptr >= argc - 1) help("error: missing parameter after -o option.");
        std::string path_m_file = argv[ptr++];
        mf.open(path_m_file, std::filesystem::exists(path_m_file) ? std::ios::app : std::ios::out);
        if (!mf.good()) help("error: cannot open nor create <m_file>");
        name_text_file = argv[ptr++];
    } else if (s == "-i") {
        if (ptr >= argc - 1) help("error: missing parameter after -i option.");
        path_input_file = argv[ptr++];
        input_file.open(path_input_file);
        if (!input_file.good()) help("error: cannot open <input_file>");
    } else if (s == "-o") {
        if (ptr >= argc - 1) help("error: missing parameter after -o option.");
        output_occurrences = true;
        path_output_file = argv[ptr++];
    } else {
        help("error: unrecognized '" + s + "' option");
    }
}

template <typename pos_t, move_r_support support>
void measure_locate()
{
    std::cout << std::setprecision(4);
    std::cout << "loading the index" << std::flush;
    auto time = now();
    using idx_t = move_r<support, char, pos_t>;
    idx_t index;
    index.load(index_file);
    time = log_runtime(time);
    index_file.close();
    std::cout << std::endl;
    index.log_data_structure_sizes();

    if (support == _locate_rlzsa_bin_search || check_correctness) {
        if (path_input_file == "") help("error: <input_file> not provided");
        std::cout << std::endl << "loading input file" << std::flush;
        no_init_resize(input, index.input_size());
        read_from_file(input_file, input.data(), input.size());
        if constexpr (support == _locate_rlzsa_bin_search) index.set_input(input);
        input_file.close();
        time = log_runtime(time);
    }

    std::cout << std::endl << "searching patterns ... " << std::endl;
    std::string header;
    std::getline(patterns_file, header);
    uint64_t num_patterns = number_of_patterns(header);
    uint64_t pattern_length = patterns_length(header);
    uint64_t perc;
    uint64_t last_perc = 0;
    uint64_t num_occurrences = 0;
    uint64_t time_locate = 0;
    std::chrono::steady_clock::time_point t2, t3;
    std::string pattern;
    no_init_resize(pattern, pattern_length);
    std::vector<pos_t> occurrences;
    bool is_sorted, equal;
    pos_t count;

    for (uint64_t i = 0; i < num_patterns; i++) {
        perc = (100 * i) / num_patterns;
        is_sorted = false;

        if (perc > last_perc) {
            std::cout << perc << "% done .." << std::endl;
            last_perc = perc;
        }

        patterns_file.read(pattern.data(), pattern_length);
        t2 = now();
        occurrences = index.locate(pattern);
        t3 = now();
        time_locate += time_diff_ns(t2, t3);
        num_occurrences += occurrences.size();

        if (check_correctness) {
            for (pos_t occ : occurrences) {
                if (input.substr(occ, pattern_length) != pattern) {
                    std::cout << "error: wrong occurrence: " << occ << " of pattern '" << pattern << "'" << std::endl;
                    exit(-1);
                }
            }
        }

        if (output_occurrences) {
            ips4o::sort(occurrences.begin(), occurrences.end());
            for (pos_t occ : occurrences) output_file << occ << " ";
            output_file << std::endl;
        }

        occurrences.clear();
    }

    std::cout << "average occurrences per pattern: " << (num_occurrences / num_patterns) << std::endl;
    std::cout << "number of patterns: " << num_patterns << std::endl;
    std::cout << "pattern length: " << pattern_length << std::endl;
    std::cout << "total number of occurrences: " << num_occurrences << std::endl;
    std::cout << "locate time: " << format_time(time_locate) << std::endl;
    std::cout << "             " << format_time(time_locate / num_patterns) << "/pattern" << std::endl;
    std::cout << "             " << format_time(time_locate / num_occurrences) << "/occurrence" << std::endl;

    if (mf.is_open()) {
        mf << "RESULT";
        mf << " algo=locate_move_r_" << move_r_support_suffix(support);
        mf << " text=" << name_text_file;
        mf << " a=" << index.balancing_parameter();
        mf << " n=" << index.input_size();
        mf << " sigma=" << std::to_string(index.alphabet_size());
        mf << " r=" << index.num_bwt_runs();
        mf << " r_=" << index.M_LF().num_intervals();

        if constexpr (idx_t::supports_multiple_locate) {
            if constexpr (idx_t::has_locate_move) {
                mf << " r__=" << index.M_Phi_m1().num_intervals();
            } else if constexpr (idx_t::has_rlzsa) {
                mf << " z=" << index.num_phrases_rlzsa();
                mf << " z_l=" << index.num_literal_phrases_rlzsa();
                mf << " z_c=" << index.num_copy_phrases_rlzsa();
            } if constexpr (idx_t::has_lzendsa) {
                mf << " z=" << index.num_phrases_lzendsa();
            }
        }

        mf << " pattern_length=" << pattern_length;
        index.log_data_structure_sizes(mf);
        mf << " num_patterns=" << num_patterns;
        mf << " num_occurrences=" << num_occurrences;
        mf << " time_locate=" << time_locate;
        mf << std::endl;
        mf.close();
    }
}

int main(int argc, char** argv)
{
    if (argc < 3) help("");
    while (ptr < argc - 2) parse_args(argv, argc, ptr);

    path_index_file = argv[ptr];
    path_patterns_file = argv[ptr + 1];

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

    if (_support == _count || _support == _locate_one) {
        std::cout << "error: this index does not support locate" << std::endl;
        exit(0);
    } else if (_support == _locate_move) {
        if (is_64_bit) {
            measure_locate<uint64_t, _locate_move>();
        } else {
            measure_locate<uint32_t, _locate_move>();
        }
    } else if (_support == _locate_rlzsa) {
        if (is_64_bit) {
            measure_locate<uint64_t, _locate_rlzsa>();
        } else {
            measure_locate<uint32_t, _locate_rlzsa>();
        }
    } else if (_support == _locate_rlzsa_bin_search) {
        if (is_64_bit) {
            measure_locate<uint64_t, _locate_rlzsa_bin_search>();
        } else {
            measure_locate<uint32_t, _locate_rlzsa_bin_search>();
        }
    } else if (_support == _locate_lzendsa) {
        if (is_64_bit) {
            measure_locate<uint64_t, _locate_lzendsa>();
        } else {
            measure_locate<uint32_t, _locate_lzendsa>();
        }
    }
}