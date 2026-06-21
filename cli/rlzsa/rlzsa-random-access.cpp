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
#include <random>

#include <misc/cli.hpp>
#include <rlzsa/rlzsa.hpp>

/**
 * @brief prints the usage information and exits
 */
void help()
{
    std::cout << "rlzsa-random-access: measures the random access time of the rlzsa construction." << std::endl << std::endl;
    std::cout << "Usage: rlzsa-random-access [...] <rlzsa file>" << std::endl;
    std::cout << "   <rlzsa file>     path to rlzsa file (should the binary representation of the rlzsa construction)." << std::endl;
    std::cout << "   -s               seed used to calculate the interval starting positions (default: 42)" << std::endl;
    std::cout << "   -q               number of queries (default: 300000)" << std::endl;
    std::cout << "   -l               interval length (default: 7)" << std::endl;
}

/**
 * @brief benchmarks random-access queries on the index
 * @tparam int_t template parameter
 * @param index_file the input stream of the index file
 * @param file_name the name of the original text file
 * @param len the length of each random-access query range
 * @param num_queries the number of queries to perform
 * @param seed the seed for the random number generator
 */
template <typename int_t>
void random_access(std::ifstream& index_file, std::string file_name, uint64_t len, uint64_t num_queries, uint64_t seed)
{
    std::cout << "Loading rlzsa" << std::flush;
    rlzsa<int_t> index;
    index.load(index_file);
    std::cout << ", done" << std::endl;

    std::cout << "Start random access" << std::flush;
    std::mt19937 gen(seed);
    std::uniform_int_distribution<uint64_t> index_distrib(0, index.input_size() - len);
    uint64_t time_ns = 0;
    uint64_t checksum = 0;

    for (uint64_t i = 0; i < num_queries; i++) {
        uint64_t beg = index_distrib(gen);
        uint64_t end = beg + len - 1;
        
        auto t1 = now();
        std::vector<int_t> res = index.template sa_values<int_t>(beg, end);
        auto t2 = now();

        for (int_t val : res) checksum += val;
        time_ns += time_diff_ns(t1, t2);
    }

    std::cout << ", done (checksum: " << checksum << ")" << std::endl;
    std::cout << "Measured time random access: " << format_time(time_ns) <<
        " (with " << num_queries << " queries)" << std::endl;

    std::cout << "RESULT"
        << " algo=rlzsa_random_access"
        << " seed=" << seed
        << " time_access=" << time_ns
        << " num_queries=" << num_queries
        << " text=" << file_name
        << " len=" << len
        << " d=" << index.delta()
        << " n=" << index.input_size()
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

    allowed_value_options.insert("-l");
    allowed_value_options.insert("-q");
    allowed_value_options.insert("-s");

    CommandLineArguments a = parse_args(argc, argv,
        allowed_value_options, allowed_literal_options, 1);

    if (!a.success) {
        help();
        return -1;
    }

    uint64_t len = 7;
    uint64_t seed = 42;
    uint64_t num_queries = 300000;
    std::string file_name = a.last_param.at(0);
    file_name = file_name.substr(file_name.find_last_of("/\\") + 1);

    for (Option value_option : a.value_options) {
        if (value_option.name == "-l") {
            len = parse_int_arg(value_option.value, "-l");
        }

        if (value_option.name == "-q") {
            num_queries = parse_int_arg(value_option.value, "-q");
        }

        if (value_option.name == "-s") {
            seed = parse_int_arg(value_option.value, "-s");
        }
    }

    std::string rlzsa_file = a.last_param.at(0);
    require_file(rlzsa_file);
    std::ifstream index_file(rlzsa_file);

    bool is_64_bit;
    index_file.read((char*) &is_64_bit, sizeof(is_64_bit));

    if (!is_64_bit) random_access<int32_t>(index_file, file_name, len, num_queries, seed);
    else            random_access<int64_t>(index_file, file_name, len, num_queries, seed);
    
    return 0;
}