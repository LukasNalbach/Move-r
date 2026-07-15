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

#include <malloc_count.h>
#include <misc/cli.hpp>
#include <rlzsa/r_index_rlzsa.hpp>

/**
 * @brief prints the usage information and exits
 */
void help()
{
    std::cout << "r-index-rlzsa-locate: locates all occurences of the provided patterns." << std::endl << std::endl;
    std::cout << "Usage: r-index-rlzsa-locate [...] <index file> <pattern file>" << std::endl;
    std::cout << "   <index file>    path to the computed index (file with extension .r-index-rlzsa)" << std::endl;
    std::cout << "   <pattern>       path to the pattern file in the pizza&chili format" << std::endl;
    std::cout << "   -c <input_file> checks correctness of the results on <input_file>" << std::endl;
}

template <typename int_t>
void locate(std::ifstream& index_file, std::ifstream& patterns_file, std::string file_name, std::string input_file)
{
    std::cerr << "Loading r-index-rlzsa index" << std::flush;
    r_index_rlzsa<int_t> index;
    index.load(index_file);
    index_file.close();
    std::cerr << " done" << std::endl;
    std::string input;

    if (!input_file.empty()) {
        std::cerr << "Loading <input_file>" << std::flush;
        std::ifstream input_ifile(input_file);
        uint64_t input_size = std::filesystem::file_size(input_file);
        no_init_resize(input, input_size);
        read_from_file(input_ifile, input.data(), input_size);
        std::cerr << " done" << std::endl;
    }

    query_stats stats = benchmark_locate(patterns_file,
        [&](std::string& p) { return index.template locate<int_t>(p); },
        input_file.empty() ? nullptr : &input);

    std::cout << "RESULT"
        << " algo=r_index_rlzsa_locate"
        << " time_locate=" << stats.time_ns
        << " size_index=" << index.size_in_bytes()
        << " num_occurrences=" << stats.occ_total
        << " text=" << file_name
        << " pattern_length=" << stats.pattern_length
        << " num_patterns=" << stats.pattern_count
        << " z=" << index.sa_encoding().num_phrases()
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
    allowed_value_options.insert("-c");

    CommandLineArguments a = parse_args(argc, argv,
        allowed_value_options, allowed_literal_options, 2);

    if (!a.success) {
        help();
        return -1;
    }

    std::string file_name = a.last_param.at(0);
    file_name = file_name.substr(file_name.find_last_of("/\\") + 1);
    std::string input_file;

    for (Option value_option : a.value_options) {
        if (value_option.name == "-c") {
            input_file = value_option.value;
        }
    }

    require_file(a.last_param.at(0));
    require_file(a.last_param.at(1));
    if (!input_file.empty()) require_file(input_file);
    std::ifstream index_file(a.last_param.at(0));
    std::ifstream patterns_file(a.last_param.at(1));
    bool is_64_bit;
    index_file.read((char*) &is_64_bit, sizeof(uint8_t));

    if (!is_64_bit) locate<int32_t>(index_file, patterns_file, file_name, input_file);
    else            locate<int64_t>(index_file, patterns_file, file_name, input_file);
    
    return 0;
}