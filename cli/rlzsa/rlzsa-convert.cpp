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
#include <omp.h>

#include <misc/cli.hpp>
#include <rlzsa/r_index_rlzsa.hpp>
#include <rlzsa/rlzsa.hpp>

std::ifstream input_file;
std::ofstream output_file;
int64_t d = -1;

void help()
{
    std::cout << "rlzsa-convert: converts an r-index-rlzsa index to an rlzsa index." << std::endl << std::endl;
    std::cout << "Usage: rlzsa-convert [options] <input index> <output index>" << std::endl;
    std::cout << "\t<input index>    path to the input index (file with extension .r-index-rlzsa)" << std::endl;
    std::cout << "\t<output index>   path to the output index (file with extension .rlzsa)" << std::endl;
    std::cout << "\t-d               delta, if not provided the sample will be about 10\% of the index size" << std::endl;
}

template <typename int_t>
void convert()
{
    std::cout << "Loading r-index-rlzsa index" << std::flush;
    r_index_rlzsa<int_t> input_index;
    input_index.load(input_file);
    input_file.close();
    std::cout << " done" << std::endl;
    
    std::cout << "extracting SA" << std::flush;
    uint64_t n = input_index.input_size();
    uint16_t p = omp_get_max_threads();
    std::vector<int_t> sa;
    no_init_resize(sa, n);

    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();

        uint64_t b = i_p * (n / p);
        uint64_t e = i_p == p - 1 ? n - 1 : ((i_p + 1) * (n / p) - 1);

        input_index.sa_encoding().template extract<int_t>(b, e,
            [&](uint64_t pos, uint64_t val){sa[pos] = val;});
    }
    
    std::cout << " done" << std::endl;
    rlzsa<int_t> output_index(input_index.sa_encoding(), sa, d, true);
    output_index.serialize(output_file);
}

int main(int argc, char** argv)
{
    std::set<std::string> allowed_value_options;
    std::set<std::string> allowed_literal_options;
    allowed_value_options.insert("-d");
    CommandLineArguments parsed_args = parse_args(argc, argv, allowed_value_options, allowed_literal_options, 2);

    if (!parsed_args.success) {
        help();
        return -1;
    }

    for (Option value_option : parsed_args.value_options) {
        if (value_option.name == "-d") {
            d = std::stol(value_option.value);
        }
    }

    input_file.open(parsed_args.last_parameter.at(0));
    output_file.open(parsed_args.last_parameter.at(1));
    uint8_t long_integer_flag;
    input_file.read((char*) &long_integer_flag, sizeof(uint8_t));
    output_file.write((char*) &long_integer_flag, sizeof(uint8_t));

    if (long_integer_flag == 0) {
        convert<int32_t>();
    } else {
        convert<int64_t>();
    }
}