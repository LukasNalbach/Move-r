/**
 * part of LukasNalbach/Move-r
 *
 * MIT License
 *
 * Copyright (c) Lukas Nalbach
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

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string_view>
#include <vector>

#include <misc/files.hpp>

#include "framework.hpp"
#include "indexes.hpp"

[[noreturn]] void help(const std::string& msg = "")
{
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "move-r-bench: benchmarks construction- and query performance of move-r, block-rlbwt-2, block-rlbwt-v," << std::endl;
    std::cout << "              block-rlbwt-r, r-index, r-index-f, rcomp-glfig and online-rlbwt;" << std::endl;
    std::cout << "              has to be executed from the base folder." << std::endl << std::endl;
    std::cout << "usage: move-r-bench [...] <input_file> <patterns_file_1> <patterns_file_2>" << std::endl;
    std::cout << "   -c                 check for correctness if possible; disables the -m option; will not print" << std::endl;
    std::cout << "                      runtime data if the runtime could be affected by checking for correctness" << std::endl;
    std::cout << "   -m <m_file>        writes measurement data to m_file" << std::endl;
    std::cout << "   <input_file>       input file" << std::endl;
    std::cout << "   <patterns_file_1>  file containing patterns (pattern length ~ number of occurrences) from <input_file>" << std::endl;
    std::cout << "                      to count and locate" << std::endl;
    std::cout << "   <patterns_file_2>  file containing patterns (pattern length << number of occurrences) from <input_file>" << std::endl;
    std::cout << "                      to locate" << std::endl;
    std::exit(0);
}

// measures all indexes
void bench_indexes()
{
    bench_index<move_r<_locate_move>, true, true>("move-r", "move_r");
    bench_index<move_r<_locate_rlzsa>, true, true>("move-rlzsa", "move_rlzsa");

    block_rlbwt_data bd = prepare_blockrlbwt();
    measure_blockrlbwt("block_rlbwt_2", bd);
    measure_blockrlbwt("block_rlbwt_v", bd);
    measure_blockrlbwt("block_rlbwt_r", bd);
    cleanup_grlbwt(bd);

    bench_index<r_index, false, false>("r-index", "r_index");

    bench_index<r_index_f<>, true, false>("r-index-f", "r_index_f");
    bench_index<rcomp_glfig, true, true>("rcomp-glfig", "rcomp_glfig");
    bench_index<online_rlbwt, true, true>("online-rlbwt", "online_rlbwt");
}

// parses the option at argv[arg_idx]
void parse_args(char** argv, int argc)
{
    std::string s = argv[arg_idx];
    arg_idx++;

    if (s == "-m") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -m option");
        path_mf = argv[arg_idx++];
        mf.open(path_mf, std::filesystem::exists(path_mf) ? std::ios::app : std::ios::out);
        if (!mf.good()) help("error: cannot open nor create measurement file");
    } else if (s == "-c") {
        check_correctness = true;
    } else {
        help("error: unrecognized '" + s + "' option");
    }
}

// the positional arguments: <input_file>, <patterns_file_1> and <patterns_file_2>
constexpr int num_positional_args = 3;
// the number of significant digits of the printed numbers
constexpr int printed_significant_digits = 4;

int main_bench_indexes(int argc, char** argv)
{
    if (argc < 1 + num_positional_args) help();
    while (arg_idx < argc - num_positional_args) parse_args(argv, argc);

    path_input_file = argv[arg_idx];
    path_patterns_file_1 = argv[arg_idx + 1];
    path_patterns_file_2 = argv[arg_idx + num_positional_args - 1];
    max_num_threads = 1;

    input_file.open(path_input_file);
    patterns_file_1.open(path_patterns_file_1);
    patterns_file_2.open(path_patterns_file_2);

    if (!input_file.good()) help("error: invalid input, could not read <input_file>");
    if (!patterns_file_1.good()) help("error: invalid input, could not read <patterns_file_1>");
    if (!patterns_file_2.good()) help("error: invalid input, could not read <patterns_file_2>");
    if (mf.is_open() && check_correctness) help("error: cannot output measurement data when checking for correctness");

    system("chmod +x build/external/pfp-thresholds/*");
    system("chmod +x build/external/block_RLBWT/*");
    system("chmod +x build/external/grlBWT/*");
    std::cout << std::setprecision(printed_significant_digits);
    text_name = path_input_file.substr(path_input_file.find_last_of("/\\") + 1);
    std::cout << "benchmarking " << path_input_file << std::flush;

    input_file.seekg(0, std::ios::end);
    input_size = (uint64_t) input_file.tellg() + 1;
    input_file.seekg(0, std::ios::beg);
    input.reserve(input_size);
    resize_no_init(input, input_size - 1);
    read_in_chunks(input_file, input.data(), input_size - 1);
    input_file.close();
    drop_from_page_cache(path_input_file, true);
    preprocess_input();

    std::cout << " (" << format_size(input_size - 1) << ")" << std::endl;
    if (check_correctness) std::cout << "correctness will be checked if possible" << std::endl;
    std::cout << std::endl;

    bench_indexes();

    patterns_file_1.close();
    patterns_file_2.close();
    drop_from_page_cache(path_patterns_file_1, true);
    drop_from_page_cache(path_patterns_file_2, true);

    return 0;
}

int main(int argc, char** argv)
{
    if (argc == 1) help();
    return main_bench_indexes(argc, argv);
}