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

#include <climits>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <misc/utils.hpp>
#include <vector>

void help(std::string msg)
{
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "patterns-to-reads: converts patterns in the Pizza&Chili" << std::endl;
    std::cout << "format to reads in the fasta format." << std::endl << std::endl;
    std::cout << "usage: patterns-to-reads <patterns> <reads>" << std::endl;
    std::cout << "       <patterns> is the file containing the input patterns" << std::endl;
    std::cout << "       <reads> is the output file to create" << std::endl;
    exit(0);
}

int main(int argc, char* argv[])
{
    if (argc != 3) help("invalid input: wrong number of arguments");
    auto time = now();

    std::ifstream patterns_file(argv[1]);
    if (!patterns_file.good()) help("invalid input: could not read <patterns>");

    std::ofstream reads_file(argv[2]);
    if (!reads_file.is_open()) help("invalid input: could not create <reads>");
    
    std::string header;
    std::getline(patterns_file, header);
    uint64_t num_patterns = number_of_patterns(header);
    uint64_t pattern_length = patterns_length(header);
    std::cout << "converting " << num_patterns << " patterns of length " << pattern_length << std::flush;

    std::string pattern;
    no_init_resize(pattern, pattern_length);
    
    for (uint64_t i = 0; i < num_patterns; i++) {
        patterns_file.read(pattern.data(), pattern_length);
        reads_file << ">" << i << std::endl;
        reads_file.write(pattern.data(), pattern_length);
        reads_file << std::endl;
    }

    patterns_file.close();
    reads_file.close();
    log_runtime(time);

    return 0;
}