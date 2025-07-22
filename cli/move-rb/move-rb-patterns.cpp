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

#include <ips4o.hpp>

void help(std::string msg)
{
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "move-rb-patterns: generate patterns from a file." << std::endl << std::endl;
    std::cout << "usage: move-rb-patterns <file> <length> <number> <switches> <patterns file> <forbidden>" << std::endl;
    std::cout << "       randomly extracts <number> substrings of length <length> from <file>," << std::endl;
    std::cout << "       avoiding substrings containing characters in <forbidden>, where" << std::endl;
    std::cout << "       the extension direction changes <switches> times." << std::endl;
    std::cout << "       The output file <patterns file> has a first line of the form:" << std::endl;
    std::cout << "       # number=<number> length=<length> file=<file> forbidden=<forbidden>" << std::endl;
    std::cout << "       and then <number> tuples of the form <pattern, directions> come" << std::endl;
    std::cout << "       successively without any separator, where directions are <length>" << std::endl;
    std::cout << "       consecutive uint8_t values, where 1 means left-extension and 2 means right-extension" << std::endl;
    exit(0);
}

int main(int argc, char* argv[])
{
    if (argc < 6 || 7 < argc)
        help("invalid input: wrong number of arguments");

    std::ifstream input_file(argv[1]);
    if (!input_file.good())
        help("invalid input: could not read <file>");

    std::cout << std::setprecision(4);
    input_file.seekg(0, std::ios::end);
    int64_t input_size = input_file.tellg();
    int64_t num_patterns = atoi(argv[3]);
    int64_t num_switches = atoi(argv[4]);
    int64_t pattern_length = atoi(argv[2]);

    if (pattern_length < 0 || pattern_length >= input_size)
        help("Error: length must be >= 1 and <= file length");
    if (num_patterns < 0) help("Error: number of patterns must be >= 1");
    if (num_switches < 0 || num_switches >= pattern_length)
        help("Error: switches must be >= 0 and < length");

    std::ofstream output_file(argv[5]);
    if (!output_file.is_open()) help("invalid input: could not create <patterns file>");

    std::string forbidden = "";
    if (argc == 7) forbidden = argv[6];
    std::vector<uint8_t> is_forbidden;

    if (!forbidden.empty()) {
        is_forbidden.resize(256, 0);
        for (uint64_t i = 0; i < forbidden.size(); i++) is_forbidden[forbidden[i]] = 1;
    }

    std::string basename = argv[1];
    basename = basename.substr(basename.find_last_of("/\\") + 1);
    output_file << "# number=" << num_patterns << " length=" << pattern_length <<
        " switches=" << num_switches << " file=" << basename << " forbidden=\n";
    input_file.seekg(0, std::ios::beg);

    std::cout << "generating " << num_patterns << " pattern queries of length " << pattern_length
        << " with " << num_switches << " random direction switches, each" << std::flush;
    std::uniform_int_distribution<uint64_t> pattern_pos_distrib(0, input_size - pattern_length - 1);
    std::uniform_int_distribution<uint32_t> start_pos_distrib(0, pattern_length - 1);
    std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);
    std::random_device rd;
    std::mt19937 gen(rd());
    uint64_t pattern_pos;
    std::string pattern;
    no_init_resize(pattern, pattern_length);
    std::vector<uint64_t> pattern_positions;
    no_init_resize(pattern_positions, pattern_length - 1);
    std::iota(pattern_positions.begin(), pattern_positions.end(), 1);
    std::vector<uint64_t> switch_positions;
    no_init_resize(switch_positions, num_switches + 1);
    switch_positions[num_switches] = pattern_length;
    std::vector<direction> dirs;
    dirs.reserve(pattern_length);
    uint64_t start_pos;
    bool found_forbidden;
    auto time = now();

    for (int64_t i = 0; i < num_patterns; i++) {
        do {
            pattern_pos = pattern_pos_distrib(gen);
            input_file.seekg(pattern_pos, std::ios::beg);
            read_from_file(input_file, pattern.c_str(), pattern_length);
            found_forbidden = false;

            if (!forbidden.empty()) {
                for (int64_t i = 0; i < pattern_length; i++) {
                    if (is_forbidden[pattern[i]] == 1) {
                        found_forbidden = true;
                        break;
                    }
                }
            }
        } while (found_forbidden);

        std::shuffle(pattern_positions.begin(), pattern_positions.end(), gen);
        std::copy(pattern_positions.begin(), pattern_positions.begin() + num_switches, switch_positions.begin());
        ips4o::sort(switch_positions.begin(), switch_positions.end());
        direction dir = prob_distrib(gen) < 0.5 ? LEFT : RIGHT;
        if (num_switches >= 1) for (uint32_t i = 0; i < switch_positions.front(); i++) dirs.emplace_back(dir);
        else for (uint32_t i = 0; i < pattern_length; i++) dirs.emplace_back(dir);
        
        for (uint32_t i = 1; i <= num_switches; i++) {
            dir = flip(dir);
            for (uint32_t j = switch_positions[i - 1]; j < switch_positions[i]; j++) dirs.emplace_back(dir);
        }
        
        start_pos = std::count(dirs.begin(), dirs.end(), LEFT) - 1;

        output_file.write(pattern.c_str(), pattern_length);
        output_file.write((char*) &start_pos, sizeof(uint64_t));
        output_file.write((char*) &dirs[0], pattern_length * sizeof(direction));

        dirs.clear();
    }

    input_file.close();
    output_file.close();
    time = log_runtime(time);
}