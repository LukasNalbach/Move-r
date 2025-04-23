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
#include <lzendsa/r_index_lzendsa.hpp>

void help()
{
    std::cout << "r-index-lzendsa-build: builds the r-index-lzendsa-Index from the input file." << std::endl << std::endl;
    std::cout << "Usage: r-index-lzendsa-build [options] <text file>" << std::endl;
    std::cout << "\t<text file>     path to the input file (should contain text)" << std::endl;
    std::cout << "\t-o              path to the desired output file (the extension .r-index-lzendsa will be added automatically)" << std::endl;
    std::cout << "\t-h              longest phrase length, leave blank or put -1 for unbounded phrase length (default: 8192)" << std::endl;
    std::cout << "\t--bigbwt        use Big-BWT instead of libsais" << std::endl;
    std::cout << "\t--f64           explicitly use 64-bit-integers regardless of the file size" << std::endl;
    std::cout << "\t--lzend-samples use SA-samples at end positions of LZ-End phrases instead of the SA-samples in the r-index" << std::endl;
    std::cout << "\t-filename       sets the filename only for the RESULT line" << std::endl;
}

int main(int argc, char** argv)
{
    std::set<std::string> allowed_value_options;
    std::set<std::string> allowed_literal_options;

    allowed_value_options.insert("-o");
    allowed_value_options.insert("-h");
    allowed_value_options.insert("-filename");
    allowed_literal_options.insert("--bigbwt");
    allowed_literal_options.insert("--f64");
    allowed_literal_options.insert("--lzend-samples");

    CommandLineArguments parsed_args = parse_args(argc, argv, allowed_value_options, allowed_literal_options, 1);

    if (!parsed_args.success) {
        help();
        return -1;
    }

    std::string o = parsed_args.last_parameter.at(0);
    o.append(".r-index-lzendsa");
    std::string filepath = parsed_args.last_parameter.at(0);
    std::string filename = filepath.substr(filepath.find_last_of("/\\") + 1);
    bool use_bigbwt = false;
    bool use64 = false;
    bool use_rindex_samples = true;
    int64_t h = 8192;

    if (parsed_args.literal_options.contains("--f64")) {
        use64 = true;
    }

    if (parsed_args.literal_options.contains("--bigbwt")) {
        use_bigbwt = true;
    }

    if (parsed_args.literal_options.contains("--lzend-samples")) {
        use_rindex_samples = false;
    }

    for (Option value_option : parsed_args.value_options) {
        if (value_option.name == "-o") {
            o = value_option.value.append(".r-index-lzendsa");
        }

        if (value_option.name == "-filename") {
            filename = value_option.value;
        }

        if (value_option.name == "-h") {
            h = std::stol(value_option.value);
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
    std::cout << "File " << parsed_args.last_parameter.at(0) << " successfully loaded (" << format_size(n) << ")" << std::endl;

    auto t1 = now();
    r_index_lzendsa<int32_t> index_32;
    r_index_lzendsa<int64_t> index_64;

    if (n <= INT32_MAX && !use64) {
        std::cout << "Using 32-bit-integers" << std::endl;
        index_32 = r_index_lzendsa<int32_t>(input, h, use_bigbwt, use_rindex_samples, true);
    } else {
        std::cout << "Using 64-bit-integers" << std::endl;
        index_64 = r_index_lzendsa<int64_t>(input, h, use_bigbwt, use_rindex_samples, true);
    }

    auto t2 = now();
    uint64_t size_index = sizeof(uint8_t);
    std::ofstream out(o);
    uint64_t z;

    if (n <= INT32_MAX && !use64) {
        z = index_32.sa_encoding().num_phrases();
        uint8_t long_integer_flag = 0;
        out.write((char*) &long_integer_flag, sizeof(uint8_t));
        index_32.serialize(out);
        size_index += index_32.size_in_bytes();
    } else {
        z = index_64.sa_encoding().num_phrases();
        uint8_t long_integer_flag = 1;
        out.write((char*) &long_integer_flag, sizeof(uint8_t));
        index_64.serialize(out);
        size_index += index_64.size_in_bytes();
    }

    out.close();
    std::cout << "Wrote " << format_size(size_index) << " bytes to disk." << std::endl;
    uint64_t memory_peak = malloc_count_peak();
    uint64_t time_ns = time_diff_ns(t1, t2);

    std::cout << "RESULT"
        << " algo=r_index_lzendsa_build"
        << " time_construction=" << time_ns
        << " peak_mem_usage=" << memory_peak
        << " text=" << filename
        << " size_index=" << size_index
        << " n=" << n
        << " z=" << z
        << " h=" << h
        << std::endl;
}