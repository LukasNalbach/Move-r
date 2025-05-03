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
#include <sstream>
#include <string>
#include <vector>

#include <misc/cli.hpp>
#include <lzendsa/lzendsa.hpp>

void help()
{
    std::cout << "lzendsa-count: count all occurences of the input patterns." << std::endl << std::endl;
    std::cout << "Usage: lzendsa-count [options] <lzendsa file> <text file> <pattern file>" << std::endl;
    std::cout << "\t<lzendsa file>  path to lzendsa file (should the binary representation of the lzendsa construction)" << std::endl;
    std::cout << "\t<text file>     path to text file (should contain text)" << std::endl;
    std::cout << "\t<pattern file>  path to file containing the pattern in pizza&chili format." << std::endl;
}

template <typename int_t>
void count(std::string& input, std::ifstream& index_file, std::ifstream& patterns_file, std::string output_filename, std::string file_name)
{
    std::cout << "Loading lzendsa" << std::flush;
    lzendsa<int_t> index;
    index.load(index_file);
    index.set_input(input);
    std::cout << " done" << std::endl;

    std::string pattern_header;
    std::getline(patterns_file, pattern_header);
    uint64_t pattern_length = get_pattern_length(pattern_header);
    uint64_t pattern_count = get_pattern_count(pattern_header);
    std::cout << "Found " << pattern_count << " patterns of length " << pattern_length << "." << std::endl;
    std::cout << "Count: " << std::flush;
    int64_t occ_total = 0;
    uint64_t time_ns = 0;
    std::string pattern;
    no_init_resize(pattern, pattern_length);

    for (uint64_t i = 0; i < pattern_count; i++) {
        patterns_file.read(pattern.data(), pattern_length);

        auto t1 = now();
        auto [beg, end] = index.template count(pattern);
        auto t2 = now();

        time_ns += time_diff_ns(t1, t2);
        occ_total += end >= beg ? end - beg + 1 : 0;
    }

    std::cout << "Counted " << pattern_count << " patterns (with " << occ_total
        << " occurences) in " << format_time(time_ns) << std::endl;

    std::cout << "RESULT"
        << " algo=lzendsa_count"
        << " time_count=" << time_ns
        << " size_index=" << index.size_in_bytes()
        << " num_occurrences=" << occ_total
        << " text=" << file_name
        << " pattern_length=" << pattern_length
        << " n=" << input.length()
        << " d=" << index.delta()
        << " h=" << index.sa_encoding().maximum_phrase_length()
        << " size_index=" << index.size_in_bytes()
        << std::endl;
}

int main(int argc, char** argv)
{
    std::set<std::string> allowed_value_options;
    std::set<std::string> allowed_literal_options;
    CommandLineArguments a = parse_args(argc, argv, allowed_value_options, allowed_literal_options, 3);

    if (!a.success) {
        help();
        return -1;
    }

    std::string pattern_out_file = "";
    std::string file_name = a.last_parameter.at(0);
    file_name = file_name.substr(file_name.find_last_of("/\\") + 1);

    std::string lzendsa_file = a.last_parameter.at(0);
    std::string text_file = a.last_parameter.at(1);
    std::string pattern_file = a.last_parameter.at(2);

    std::ifstream in(lzendsa_file);
    std::ifstream text_in(text_file);
    std::ifstream patterns_file(pattern_file);

    std::cout << "Loading text file" << std::flush;
    std::string input;
    uint64_t input_size = std::filesystem::file_size(text_file);
    no_init_resize(input, input_size);
    read_from_file(text_in, input.data(), input_size);
    std::cout << " (" << format_size(input_size) << ")" << std::endl;

    uint8_t long_integer_flag;
    in.read((char*) &long_integer_flag, sizeof(long_integer_flag));

    if (long_integer_flag == 0) {
        count<int32_t>(input, in, patterns_file, pattern_out_file, file_name);
    } else {
        count<int64_t>(input, in, patterns_file, pattern_out_file, file_name);
    }
}