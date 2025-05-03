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
#include "malloc_count.h"

#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <string>

#include <misc/cli.hpp>
#include <rlzsa/rlzsa.hpp>

void help()
{
    std::cout << "rlzsa-build: builds the rlzsa-index from the input file." << std::endl << std::endl;
    std::cout << "Usage: rlzsa-build [options] <text file>" << std::endl;
    std::cout << "\t<text file>     path to the input file (should contain text)" << std::endl;
    std::cout << "\t-o              path to the desired output file (the extension .rlzsa will be added automatically)" << std::endl;
    std::cout << "\t-d              delta, if not provided the sample will be about 10\% of the index size" << std::endl;
    std::cout << "\t-ref-size       reference size (in relation to the input size, default = 0.1 (10\%))" << std::endl;
    std::cout << "\t--bigbwt        use Big-BWT instead of libsais" << std::endl;
    std::cout << "\t--f64           explicitly use 64-bit-integers regardless of the file size" << std::endl;
}

int main(int argc, char** argv)
{
    std::set<std::string> allowed_value_options;
    std::set<std::string> allowed_literal_options;

    allowed_value_options.insert("-o");
    allowed_value_options.insert("-d");
    allowed_value_options.insert("-ref-size");
    allowed_literal_options.insert("--f64");
    allowed_literal_options.insert("--bigbwt");

    CommandLineArguments parsed_args = parse_args(argc, argv, allowed_value_options, allowed_literal_options, 1);

    if (!parsed_args.success) {
        help();
        return -1;
    }

    std::string o = parsed_args.last_parameter.at(0);
    o.append(".rlzsa");
    std::string filepath = parsed_args.last_parameter.at(0);
    std::string file_name = filepath.substr(filepath.find_last_of("/\\") + 1);
    bool use_bigbwt = false;
    bool use64 = false;
    int64_t d = -1;
    double relative_reference_size = 0.1;

    if (parsed_args.literal_options.contains("--f64")) {
        use64 = true;
    }

    if (parsed_args.literal_options.contains("--bigbwt")) {
        use_bigbwt = true;
    }

    for (Option value_option : parsed_args.value_options) {
        if (value_option.name == "-o") {
            o = value_option.value.append(".rlzsa");
        }

        if (value_option.name == "-d") {
            d = std::stol(value_option.value);
        }

        if (value_option.name == "-ref-size") {
            relative_reference_size = std::stod(value_option.value);
        }
    }

    std::string input;

    if (use_bigbwt) {
        input = filepath;
    } else {
        std::ifstream ifs(parsed_args.last_parameter.at(0));
        input = std::string(std::istreambuf_iterator<char>(ifs), {});
    }

    int64_t n = use_bigbwt ? std::filesystem::file_size(input) : input.length();

    auto t1 = now();
    rlzsa<int32_t> rlzsa_32;
    rlzsa<int64_t> rlzsa_64;

    if (n <= INT32_MAX && !use64) {
        std::cout << "Constructing using 32-bit integers" << std::endl;
        rlzsa_32 = rlzsa<int32_t>(input, relative_reference_size, d, use_bigbwt, true);
    } else {
        std::cout << "Constructing using 64-bit integers" << std::endl;
        rlzsa_64 = rlzsa<int64_t>(input, relative_reference_size, d, use_bigbwt, true);
    }

    auto t2 = now();
    uint64_t z;

    if (n <= INT32_MAX && !use64) {
        z = rlzsa_32.num_phrases();
        rlzsa_32.size_in_bytes();
    } else {
        z = rlzsa_64.num_phrases();
        rlzsa_64.size_in_bytes();
    }

    std::cout << "rlzsa index successfully constructed (z=" << z << ")" << std::endl;
    uint64_t size_index = sizeof(uint8_t);
    std::ofstream out(o);

    if (n <= INT32_MAX && !use64) {
        uint8_t long_integer_flag = 0;
        out.write((char*) &long_integer_flag, sizeof(uint8_t));
        rlzsa_32.serialize(out);
        size_index += rlzsa_32.size_in_bytes();
    } else {
        uint8_t long_integer_flag = 1;
        out.write((char*) &long_integer_flag, sizeof(uint8_t));
        rlzsa_64.serialize(out);
        size_index += rlzsa_64.size_in_bytes();
    }

    out.close();
    std::cout << "Wrote " << format_size(size_index) << " bytes to disk." << std::endl;
    uint64_t memory_peak = malloc_count_peak();
    uint64_t time_ns = time_diff_ns(t1, t2);

    std::cout << "RESULT"
        << " algo=rlzsa_build"
        << " time_ms=" << time_ns
        << " peak_mem_usage=" << memory_peak
        << " text=" << file_name
        << " size_index=" << size_index
        << " n=" << n
        << " z=" << z
        << " d=" << d
        << std::endl;
}