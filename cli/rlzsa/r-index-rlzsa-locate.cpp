/**
 * part of LukasNalbach/Move-r
 *
 * MIT License
 *
 * Copyright (c) Jan Zumbrink, Lukas Nalbach
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

#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <string>

#include "malloc_count.h"
#include <misc/cli.hpp>
#include <rlzsa/r_index_rlzsa.hpp>

void help()
{
    std::cout << "r-index-rlzsa-locate: locates all occurences of the provided patterns." << std::endl << std::endl;
    std::cout << "Usage: r-index-rlzsa-locate [options] <index file> <pattern file>" << std::endl;
    std::cout << "\t<index file>    path to the computed index (file with extension .r-index-rlzsa)" << std::endl;
    std::cout << "\t<pattern>       path to the pattern file in the pizza&chili format" << std::endl;
    std::cout << "\t-c <input file> checks correctness of the results on <input file>" << std::endl;
}

template <typename int_t>
void locate(std::ifstream& index_file, std::ifstream& pattern_in, std::string file_name, std::string input_file)
{
    size_t pre_load_memory = malloc_count_current();
    std::cout << "Loading r-index-rlzsa index" << std::flush;
    r_index_rlzsa<int_t> index;
    index.load(index_file);
    std::cout << " done" << std::endl;
    std::string input;

    if (!input_file.empty()) {
        std::cout << "Loading <input file>" << std::flush;
        std::ifstream input_ifile(input_file);
        uint64_t input_size = std::filesystem::file_size(input_file);
        no_init_resize(input, input_size);
        read_from_file(input_ifile, input.data(), input_size);
        std::cout << " done" << std::endl;
    }

    std::cout << "Loading patterns" << std::flush;
    std::string pattern_header;
    std::getline(pattern_in, pattern_header);
    uint64_t pattern_length = get_pattern_length(pattern_header);
    uint64_t pattern_count = get_pattern_count(pattern_header);
    std::vector<std::string> patterns = load_patterns(pattern_in, pattern_length, pattern_count);
    std::cout << " found " << pattern_count << " patterns of length " << pattern_length << "." << std::endl;

    std::cout << "Locate: " << std::flush;
    int64_t occ_total = 0;
    auto t1 = now();

    for (std::string& pattern : patterns) {
        std::vector<int_t> occurrences = index.template locate<int_t>(pattern);
        occ_total += occurrences.size();
        
        if (!input_file.empty()) {
            for (int_t occ : occurrences) {
                if (input.substr(occ, pattern.size()) != pattern) {
                    std::cout << "error: wrong occurrence: " << occ << " of pattern '" << pattern << "'" << std::endl;
                    exit(-1);
                }
            }
        }
    }

    auto t2 = now();
    uint64_t time_ns = time_diff_ns(t1, t2);
    std::cout << "Located " << patterns.size() << " patterns (with " << occ_total << " occurences) in " << format_time(time_ns) << std::endl;
    uint64_t size_index = index.size_in_bytes();

    std::cout << "RESULT"
        << " algo=r_index_rlzsa_locate"
        << " time_locate=" << time_ns
        << " size_index=" << size_index
        << " num_occurrences=" << occ_total
        << " text=" << file_name
        << " pattern_length=" << pattern_length
        << " num_patterns=" << patterns.size()
        << " z=" << index.sa_encoding().num_phrases()
        << std::endl;
}

int main(int argc, char** argv)
{
    std::set<std::string> allowed_value_options;
    std::set<std::string> allowed_literal_options;
    allowed_value_options.insert("-c");

    CommandLineArguments a = parse_args(argc, argv, allowed_value_options, allowed_literal_options, 2);

    if (!a.success) {
        help();
        return -1;
    }

    std::string file_name = a.last_parameter.at(0);
    file_name = file_name.substr(file_name.find_last_of("/\\") + 1);
    std::string input_file;

    for (Option value_option : a.value_options) {
        if (value_option.name == "-c") {
            input_file = value_option.value;
        }
    }

    std::ifstream index_file(a.last_parameter.at(0));
    std::ifstream patterns_file(a.last_parameter.at(1));
    uint8_t long_integer_flag;
    index_file.read((char*) &long_integer_flag, sizeof(uint8_t));

    if (long_integer_flag == 0) {
        locate<int32_t>(index_file, patterns_file, file_name, input_file);
    } else {
        locate<int64_t>(index_file, patterns_file, file_name, input_file);
    }
}