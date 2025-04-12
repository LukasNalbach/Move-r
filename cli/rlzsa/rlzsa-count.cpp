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
#include <rlzsa/rlzsa.hpp>

void help()
{
    std::cout << "rlzsa-count: count all occurences of the input patterns." << std::endl << std::endl;
    std::cout << "Usage: rlzsa-count [options] <rlzsa file> <text file> <pattern file>" << std::endl;
    std::cout << "\t<rlzsa file>     path to rlzsa file (should the binary representation of the rlzsa construction)" << std::endl;
    std::cout << "\t<text file>        path to text file (should contain text)" << std::endl;
    std::cout << "\t<pattern file>     path to file containing the pattern in pizza&chili format." << std::endl;
    std::cout << "\t-o                 file where the intervals will be written to (if no value is provided, then the results are not written onto the disk)." << std::endl;
    std::cout << "\t(-filename         sets the filename only for the RESULT line)" << std::endl;
    std::cout << "\t(-h                sets h only for the RESULT line)" << std::endl;
}

template <typename int_t>
void count(std::string& input, std::ifstream& index_file, std::ifstream& patterns_file, std::string output_filename, std::string filename, int32_t h)
{
    std::cout << "Loading rlzsa" << std::flush;
    rlzsa<int_t> index;
    index.load(index_file);
    std::cout << " done." << std::endl;

    std::cout << "Loading patterns" << std::flush;
    std::string pattern_header;
    std::getline(patterns_file, pattern_header);
    uint64_t pattern_length = get_pattern_length(pattern_header);
    uint64_t pattern_count = get_pattern_count(pattern_header);
    std::vector<std::string> patterns_str = load_patterns(patterns_file, pattern_length, pattern_count);
    std::vector<std::vector<uint8_t>> patterns = strToUint8Vec(patterns_str);
    std::cout << " found " << pattern_count << " patterns of length " << pattern_length << "." << std::endl;
    std::cout << "Count" << std::flush;

    struct Interval {
        int64_t beg;
        int64_t end;
    };

    auto t1 = now();
    std::vector<Interval> intervals;
    intervals.reserve(pattern_count);

    for (std::vector<uint8_t> pattern : patterns) {
        Interval interval;// = index.count(input, pattern);
        intervals.push_back(interval);
    }

    auto t2 = now();
    int64_t occ_total = 0;

    for (Interval interval : intervals) {
        occ_total += interval.end - interval.beg + 1;
    }
    
    uint64_t time_ns = time_diff_ns(t1, t2);
    std::cout << std::endl << "Found all patterns in " << format_time(time_ns)
        << " ms (occ_total=" << occ_total << ")." << std::endl;

    if (output_filename.length() > 0) {
        std::ofstream out(output_filename);

        for (Interval interval : intervals) {
            out << interval.beg << " " << interval.end << std::endl;
        }

        out.close();
    }

    std::cout << "RESULT"
        << " algo=rlzsa_count"
        << " time_ns=" << time_ns
        << " index_size=" << index.size_in_bytes()
        << " occ_total=" << occ_total
        << " file=" << filename
        << " m=" << pattern_length
        << " n=" << input.length()
        << " d=" << index.delta()
        << " h=" << h
        << " index_size=" << index.size_in_bytes()
        << std::endl;
}

int main(int argc, char** argv)
{
    std::set<std::string> allowed_value_options;
    std::set<std::string> allowed_literal_options;

    allowed_value_options.insert("-o");
    allowed_value_options.insert("-h");
    allowed_value_options.insert("-filename");

    CommandLineArguments a = parse_args(argc, argv, allowed_value_options, allowed_literal_options, 3);

    if (!a.success) {
        help();
        return -1;
    }

    std::string pattern_out_file = "";
    std::string filename = a.last_parameter.at(0);
    filename = filename.substr(filename.find_last_of("/\\") + 1);
    int32_t h = -1;

    for (Option value_option : a.value_options) {
        if (value_option.name == "-o") {
            pattern_out_file = value_option.value.append(".intervals");
        }

        if (value_option.name == "-filename") {
            filename = value_option.value;
        }

        if (value_option.name == "-h") {
            h = std::stoi(value_option.value);
        }
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
        count<int32_t>(input, in, patterns_file, pattern_out_file, filename, h);
    } else {
        count<int64_t>(input, in, patterns_file, pattern_out_file, filename, h);
    }
}