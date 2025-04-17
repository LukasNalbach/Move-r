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
#include <vector>

#include <misc/cli.hpp>
#include <rlzsa/rlzsa.hpp>

void help()
{
    std::cout << "rlzsa-locate: locates all occurences in the suffix array intervals." << std::endl << std::endl;
    std::cout << "Usage: rlzsa-locate <rlzsa file> <text file> <pattern file>" << std::endl;
    std::cout << "\t<rlzsa file>     path to rlzsa file (should the binary representation of the rlzsa construction)" << std::endl;
    std::cout << "\t<text file>        path to text file (should contain text)" << std::endl;
    std::cout << "\t<pattern file>     path to file containing the pattern in pizza&chili format." << std::endl;
    std::cout << "\t--c                checks correctness of the results on <text file>" << std::endl;
    std::cout << "\t(-filename         sets the filename only for the RESULT line)" << std::endl;
    std::cout << "\t(-m                sets m only for the RESULT line)" << std::endl;
}

template <typename int_t>
void locate(std::string& input, std::ifstream& index_file, std::ifstream& patterns_file, std::string filename, bool check_correctness, int64_t m)
{
    std::cout << "Loading rlzsa" << std::flush;
    rlzsa<int_t> index;
    index.load(index_file);
    index.set_input(input);
    std::cout << " done." << std::endl;

    std::cout << "Loading patterns" << std::flush;
    std::string pattern_header;
    std::getline(patterns_file, pattern_header);
    uint64_t pattern_length = get_pattern_length(pattern_header);
    uint64_t pattern_count = get_pattern_count(pattern_header);
    std::vector<std::string> patterns = load_patterns(patterns_file, pattern_length, pattern_count);
    std::cout << " found " << pattern_count << " patterns of length " << pattern_length << "." << std::endl;
    std::cout << "Locate" << std::flush;
    int64_t occ_total = 0;
    auto t1 = now();

    for (std::string& pattern : patterns) {
        std::vector<int_t> occurrences = index.template locate<int_t>(pattern);
        occ_total += occurrences.size();
        
        if (check_correctness) {
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
    int_t n = index.input_size();
    std::cout << "Located " << pattern_count << " patterns (with " << occ_total
        << " occurences) in " << format_time(time_ns) << std::endl;
    
    std::cout << "RESULT"
        << " algo=rlzsa_locate"
        << " time_ns=" << time_ns
        << " occ_total=" << occ_total
        << " file=" << filename
        << " m=" << m
        << " n=" << n
        << " d=" << index.delta()
        << " index_size=" << index.size_in_bytes()
        << std::endl;
}

int main(int argc, char** argv)
{
    std::set<std::string> allowed_value_options;
    std::set<std::string> allowed_literal_options;

    allowed_value_options.insert("-filename");
    allowed_value_options.insert("-m");
    allowed_value_options.insert("-t");
    allowed_literal_options.insert("--c");

    CommandLineArguments a = parse_args(argc, argv, allowed_value_options, allowed_literal_options, 3);

    if (!a.success) {
        help();
        return -1;
    }

    std::string filename = a.last_parameter.at(0);
    filename = filename.substr(filename.find_last_of("/\\") + 1);
    bool check_correctness = false;
    int64_t m = 0;

    for (Option value_option : a.value_options) {
        if (value_option.name == "-filename") {
            filename = value_option.value;
        }

        if (value_option.name == "-m") {
            m = std::stol(value_option.value);
        }
    }

    if (a.literal_options.contains("--c")) {
        check_correctness = true;
    }

    std::string rlzsa_file = a.last_parameter.at(0);
    std::string text_file = a.last_parameter.at(1);
    std::string pattern_file = a.last_parameter.at(2);

    std::ifstream in(rlzsa_file);
    std::ifstream text_in(text_file);
    std::ifstream patterns_file(pattern_file);

    std::cout << "Loading text file" << std::flush;
    std::string input = std::string(std::istreambuf_iterator<char>(text_in), {});
    std::cout << " done (n = " << input.length() << ")." << std::endl;

    uint8_t long_integer_flag;
    in.read((char*) &long_integer_flag, sizeof(long_integer_flag));

    if (long_integer_flag == 0) {
        locate<int32_t>(input, in, patterns_file, filename, check_correctness, m);
    } else {
        locate<int64_t>(input, in, patterns_file, filename, check_correctness, m);
    }
}