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
    std::cout << "\t(-h             sets h only for the RESULT line)" << std::endl;
    std::cout << "\t-filename       sets the filename only for the RESULT line" << std::endl;
}

template <typename int_t>
void locate(std::ifstream& index_file, std::ifstream& pattern_in, std::string& filename, int32_t h)
{
    size_t pre_load_memory = malloc_count_current();
    std::cout << "Loading r-index-lzendsa index" << std::flush;
    r_index_lzendsa<int_t> index;
    index.load(index_file);
    std::cout << " done." << std::endl;

    std::cout << "Loading patterns" << std::flush;
    std::string pattern_header;
    std::getline(pattern_in, pattern_header);
    uint64_t pattern_length = get_pattern_length(pattern_header);
    uint64_t pattern_count = get_pattern_count(pattern_header);
    std::vector<std::string> patterns = load_patterns(pattern_in, pattern_length, pattern_count);
    std::cout << " found " << pattern_count << " patterns of length " << pattern_length << "." << std::endl;

    std::cout << "Locate" << std::flush;
    int64_t occ_total = 0;
    auto t1 = now();

    for (std::string& pattern : patterns) {
        std::vector<int_t> x = index.locate(pattern);
        occ_total += x.size();
    }

    auto t2 = now();
    uint64_t time_ns = time_diff_ns(t1, t2);
    std::cout << "Located " << patterns.size() << " patterns (with " << occ_total << " occurences) in " << format_time(time_ns) << std::endl;
    uint64_t index_size = index.size_in_bytes();

    std::cout << "RESULT"
        << " algo=r_index_lzendsa_locate"
        << " time_ns=" << time_ns
        << " index_size=" << index_size
        << " occ_total=" << occ_total
        << " file=" << filename
        << " m=" << pattern_length
        << " z=" << index.encoding().num_phrases()
        << " h=" << h
        << std::endl;
}

int main(int argc, char** argv)
{
    std::set<std::string> allowed_value_options;
    std::set<std::string> allowed_literal_options;

    allowed_value_options.insert("-h");
    allowed_value_options.insert("-filename");

    CommandLineArguments a = parse_args(argc, argv, allowed_value_options, allowed_literal_options, 2);

    if (!a.success) {
        help();
        return -1;
    }

    int32_t h = 0;
    std::string filename = a.last_parameter.at(0);
    filename = filename.substr(filename.find_last_of("/\\") + 1);

    for (Option value_option : a.value_options) {
        if (value_option.name == "-h") {
            h = std::stoi(value_option.value);
        }

        if (value_option.name == "-filename") {
            filename = value_option.value;
        }
    }

    std::ifstream index_file(a.last_parameter.at(0));
    std::ifstream patterns_file(a.last_parameter.at(1));
    uint8_t long_integer_flag;
    index_file.read((char*) &long_integer_flag, sizeof(uint8_t));

    if (long_integer_flag == 0) {
        locate<int32_t>(index_file, patterns_file, filename, h);
    } else {
        locate<int64_t>(index_file, patterns_file, filename, h);
    }
}