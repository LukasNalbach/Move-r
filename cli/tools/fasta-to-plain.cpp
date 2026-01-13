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

#include <filesystem>
#include <climits>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <misc/utils.hpp>
#include <misc/log.hpp>
#include <misc/files.hpp>

void help(std::string msg)
{
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "fasta-to-plain: converts a fasta file to a plain text file" << std::endl;
    std::cout << "usage: fasta-to-plain <input_file> <output_file>" << std::endl;
    exit(0);
}

int main(int argc, char* argv[])
{
    if (argc != 3) help("error: wrong number of arguments");
    auto time = now();

    std::string path_ifile = argv[1];
    std::string path_ofile = argv[2];

    std::ifstream ifile(path_ifile);
    if (!ifile.good()) help("error: could not read <input_file>");

    std::ofstream ofile(path_ofile);
    if (!ofile.good()) help("error: could not create <output_file>");

    uint64_t input_size = std::filesystem::file_size(path_ifile);
    std::cout << "converting " << path_ifile << " (" <<
        format_size(input_size) << ") to a plain file" << std::flush;

    time = now();
    std::string line;

    while (!std::getline(ifile, line).eof()) {
        if (line[0] != '>') [[likely]] {
            ofile.write(line.data(), line.size());
        }
    }

    ifile.close();
    ofile.close();

    uint64_t output_size = std::filesystem::file_size(path_ofile);
    std::cout << " (" << format_size(output_size) << ")";
    log_runtime(time);

    return 0;
}