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

#include <misc/cli.hpp>
#include <rlzsa/rlzsa.hpp>

void help()
{
    std::cout << "rlzsa-random-access: measures the random access time of the rlzsa construction." << std::endl << std::endl;
    std::cout << "Usage: rlzsa-random-access [options] <rlzsa file> <integer file>" << std::endl;
    std::cout << "\t<rlzsa file>     path to rlzsa file (should the binary representation of the rlzsa construction)." << std::endl;
    std::cout << "\t<integer file>     path to file which contains an integer value at each line." << std::endl;
    std::cout << "\t-l                 interval length" << std::endl;
    std::cout << "\t(-filename         sets the filename, only for the RESULT line)" << std::endl;
}

template <typename int_t>
void random_access(std::ifstream& index_file, std::vector<int_t>& indices, std::string filename, int64_t interval_length)
{
    std::cout << "Loading rlzsa" << std::flush;
    rlzsa<int_t> index;
    index.load(index_file);
    std::cout << " done." << std::endl;

    std::cout << "Start random access" << std::flush;
    auto t1 = now();

    for (int_t i = 0; i < indices.size(); i++) {
        std::vector<int_t> _x = index.sa_values(indices[i], interval_length);
    }

    auto t2 = now();
    uint64_t time_ns = time_diff_ns(t1, t2);
    std::cout << " done." << std::endl;
    std::cout << "Measured time random access: " << format_time(time_ns) << " (with " << indices.size() << " Iterations)" << std::endl;

    std::cout << "RESULT"
        << " algo=rlzsa_random_access"
        << " time_ns=" << time_ns
        << " iterations=" << indices.size()
        << " file=" << filename
        << " interval_length=" << interval_length
        << " d=" << index.delta()
        << " n=" << index.input_size()
        << " index_size=" << index.size_in_bytes()
        << std::endl;
}

int main(int argc, char** argv)
{
    std::set<std::string> allowed_value_options;
    std::set<std::string> allowed_literal_options;

    allowed_value_options.insert("-l");
    allowed_value_options.insert("-filename");
    CommandLineArguments a = parse_args(argc, argv, allowed_value_options, allowed_literal_options, 2);

    if (!a.success) {
        help();
        return -1;
    }

    int64_t l = 1;
    std::string filename = a.last_parameter.at(0);
    filename = filename.substr(filename.find_last_of("/\\") + 1);

    for (Option value_option : a.value_options) {
        if (value_option.name == "-l") {
            l = std::stol(value_option.value);
        }

        if (value_option.name == "-filename") {
            filename = value_option.value;
        }
    }

    std::string rlzsa_file = a.last_parameter.at(0);
    std::string arbitrary_indices_file = a.last_parameter.at(1);

    std::ifstream index_file(rlzsa_file);
    std::ifstream arbitrary_indices_in(arbitrary_indices_file);

    uint8_t long_integer_flag;
    index_file.read((char*) &long_integer_flag, sizeof(long_integer_flag));

    std::cout << "Loading Indices" << std::flush;
    std::string last_index;

    if (long_integer_flag == 0) {
        std::vector<int32_t> indices;

        while (arbitrary_indices_in >> last_index) {
            indices.push_back(std::stoi(last_index));
        }

        std::cout << " done." << std::endl << std::flush;
        random_access<int32_t>(index_file, indices, filename, l);
    } else {
        std::vector<int64_t> indices;

        while (arbitrary_indices_in >> last_index) {
            indices.push_back(std::stol(last_index));
        }

        std::cout << " done." << std::endl << std::flush;
        random_access<int64_t>(index_file, indices, filename, l);
    }
}