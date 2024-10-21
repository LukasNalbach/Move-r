#include <climits>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <move_r/misc/utils.hpp>
#include <vector>

void help(std::string msg)
{
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "move-r-patterns: generate patterns from a file." << std::endl << std::endl;
    std::cout << "usage: move-r-patterns <file> <length> <number> <patterns file> <forbidden>" << std::endl;
    std::cout << "       randomly extracts <number> substrings of length <length> from <file>," << std::endl;
    std::cout << "       avoiding substrings containing characters in <forbidden>." << std::endl;
    std::cout << "       The output file, <patterns file> has a first line of the form:" << std::endl;
    std::cout << "       # number=<number> length=<length> file=<file> forbidden=<forbidden>" << std::endl;
    std::cout << "       and then the <number> patterns come successively without any separator" << std::endl;
    exit(0);
}

int main(int argc, char* argv[])
{
    if (argc < 5 || 6 < argc)
        help("invalid input: wrong number of arguments");

    std::ifstream input_file(argv[1]);
    if (!input_file.good())
        help("invalid input: could not read <file>");

    std::cout << std::setprecision(4);
    input_file.seekg(0, std::ios::end);
    int64_t input_size = input_file.tellg();
    int64_t num_patterns = atoi(argv[3]);
    int64_t pattern_length = atoi(argv[2]);

    if (pattern_length < 0 || pattern_length >= input_size)
        help("Error: length must be >= 1 and <= file length");
    if (num_patterns < 0)
        help("Error: number of patterns must be >= 1");

    std::ofstream output_file(argv[4]);
    if (!output_file.is_open())
        help("invalid input: could not create <patterns file>");

    std::string forbidden = "";
    if (argc == 6) forbidden = argv[5];
    std::vector<uint8_t> is_forbidden;

    if (!forbidden.empty()) {
        is_forbidden.resize(256, 0);
        for (uint64_t i = 0; i < forbidden.size(); i++)
            is_forbidden[forbidden[i]] = 1;
    }

    std::string basename = argv[1];
    basename = basename.substr(basename.find_last_of("/\\") + 1);
    output_file << "# number=" << num_patterns << " length=" << pattern_length << " file=" << basename << " forbidden=\n";
    input_file.seekg(0, std::ios::beg);

    std::cout << "generating " << num_patterns << " petterns of length " << pattern_length << std::flush;
    uint64_t pos_random;
    std::string pattern;
    no_init_resize(pattern, pattern_length);
    bool found_forbidden = false;
    auto time = now();

    for (int64_t i = 0; i < num_patterns; i++) {
        do {
            pos_random = std::rand() % (input_size - pattern_length);
            input_file.seekg(pos_random, std::ios::beg);
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

        output_file.write(pattern.c_str(), pattern_length);
    }

    input_file.close();
    output_file.close();
    time = log_runtime(time);
}