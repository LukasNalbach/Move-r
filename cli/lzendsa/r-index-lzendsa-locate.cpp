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
#include <lzendsa/r_index_lzendsa.hpp>

void help()
{
    std::cout << "r-index-lzendsa-locate: locates all occurences of the provided patterns." << std::endl << std::endl;
    std::cout << "Usage: r-index-lzendsa-locate [options] <index file> <pattern file>" << std::endl;
    std::cout << "\t<index file>    path to the computed index (file with extension .r-index-lzendsa)" << std::endl;
    std::cout << "\t<pattern>       path to the pattern file in the pizza&chili format" << std::endl;
    std::cout << "\t-c <input file> checks correctness of the results on <input file>" << std::endl;
}

template <typename int_t>
void locate(std::ifstream& index_file, std::ifstream& patterns_file, std::string file_name, std::string input_file)
{
    size_t pre_load_memory = malloc_count_current();
    std::cout << "Loading r-index-lzendsa index" << std::flush;
    r_index_lzendsa<int_t> index;
    index.load(index_file);
    std::cout << " done" << std::endl;
    std::string input;

    if (!input_file.empty()) {
        std::cout << "Loading <input file>" << std::flush;
        uint64_t input_size = std::filesystem::file_size(input_file);
        std::ifstream input_ifile(input_file);
        no_init_resize(input, input_size);
        read_from_file(input_ifile, input.data(), input_size);
        std::cout << " done" << std::endl;
    }

    std::string pattern_header;
    std::getline(patterns_file, pattern_header);
    uint64_t pattern_length = get_pattern_length(pattern_header);
    uint64_t pattern_count = get_pattern_count(pattern_header);
    std::cout << "Found " << pattern_count << " patterns of length " << pattern_length << "." << std::endl;
    std::cout << "Locate: " << std::flush;
    int64_t occ_total = 0;
    uint64_t time_ns = 0;
    std::string pattern;
    no_init_resize(pattern, pattern_length);

    for (uint64_t i = 0; i < pattern_count; i++) {
        patterns_file.read(pattern.data(), pattern_length);

        auto t1 = now();
        std::vector<int_t> occurrences = index.locate(pattern);
        auto t2 = now();

        time_ns += time_diff_ns(t1, t2);
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

    std::cout << "Located " << pattern_count << " patterns (with " << occ_total << " occurences) in " << format_time(time_ns) << std::endl;
    uint64_t size_index = index.size_in_bytes();

    std::cout << "RESULT"
        << " algo=r_index_lzendsa_locate"
        << " time_locate=" << time_ns
        << " size_index=" << size_index
        << " num_occurrences=" << occ_total
        << " text=" << file_name
        << " pattern_length=" << pattern_length
        << " num_patterns=" << pattern_count
        << " z=" << index.sa_encoding().num_phrases()
        << " h=" << index.sa_encoding().maximum_phrase_length()
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