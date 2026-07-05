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

/**
 * @brief prints the usage information and exits
 */
void help()
{
    std::cout << "rlzsa-count: count all occurences of the input patterns." << std::endl << std::endl;
    std::cout << "Usage: rlzsa-count [...] <rlzsa file> <text file> <pattern file>" << std::endl;
    std::cout << "   <rlzsa file>     path to rlzsa file (should the binary representation of the rlzsa construction)" << std::endl;
    std::cout << "   <text file>      path to text file (should contain text)" << std::endl;
    std::cout << "   <pattern file>   path to file containing the pattern in pizza&chili format." << std::endl;
    std::cout << "   -o               file where the intervals will be written to (if no value is provided, then the results are not written onto the disk)." << std::endl;
}

template <typename int_t>
void count(std::string& input, std::ifstream& index_file, std::ifstream& patterns_file, std::string /*output_filename*/, std::string file_name)
{
    std::cout << "Loading rlzsa" << std::flush;
    rlzsa<int_t> index;
    index.load(index_file);
    index.set_input(input);
    std::cout << " done" << std::endl;

    query_stats stats = benchmark_count(patterns_file,
        [&](std::string& p) { return index.count(p); });

    std::cout << "RESULT"
        << " algo=rlzsa_count"
        << " time_count=" << stats.time_ns
        << " size_index=" << index.size_in_bytes()
        << " num_occurrences=" << stats.occ_total
        << " text=" << file_name
        << " pattern_length=" << stats.pattern_length
        << " n=" << input.length()
        << " d=" << index.delta()
        << " size_index=" << index.size_in_bytes()
        << std::endl;
}

/**
 * @brief program entry point
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments
 * @return the exit code
 */
int main(int argc, char** argv)
{
    std::set<std::string> allowed_value_options;
    std::set<std::string> allowed_literal_options;
    allowed_value_options.insert("-o");
    CommandLineArguments a = parse_args(argc, argv,
        allowed_value_options, allowed_literal_options, 3);

    if (!a.success) {
        help();
        return -1;
    }

    std::string pattern_out_file = "";
    std::string file_name = a.last_param.at(0);
    file_name = file_name.substr(file_name.find_last_of("/\\") + 1);

    for (Option value_option : a.value_options) {
        if (value_option.name == "-o") {
            pattern_out_file = value_option.value.append(".intervals");
        }
    }

    std::string rlzsa_file = a.last_param.at(0);
    std::string text_file = a.last_param.at(1);
    std::string pattern_file = a.last_param.at(2);

    require_file(rlzsa_file);
    require_file(text_file);
    require_file(pattern_file);

    std::ifstream in(rlzsa_file);
    std::ifstream text_in(text_file);
    std::ifstream patterns_file(pattern_file);

    std::cout << "Loading text file" << std::flush;
    std::string input;
    uint64_t input_size = std::filesystem::file_size(text_file);
    no_init_resize(input, input_size);
    read_from_file(text_in, input.data(), input_size);
    std::cout << " (" << format_size(input_size) << ")" << std::endl;

    bool is_64_bit;
    in.read((char*) &is_64_bit, sizeof(is_64_bit));

    if (!is_64_bit) count<int32_t>(input, in, patterns_file, pattern_out_file, file_name);
    else            count<int64_t>(input, in, patterns_file, pattern_out_file, file_name);
    
    return 0;
}