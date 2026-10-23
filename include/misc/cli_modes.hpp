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

/**
 * @file cli_modes.hpp
 * @brief the modes of a command-line tool (e.g. move-r build, move-r count), each implemented in its own translation
 *        unit
 */

#ifndef MOVE_R_MISC_CLI_MODES_HPP
#define MOVE_R_MISC_CLI_MODES_HPP

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

/** @brief a mode of a command-line tool */
struct cli_mode_t {
    const char* name; // name of the mode, the first argument of the tool
    const char* description; // one line describing the mode
    int (*run)(int argc, char** argv); // runs the mode; argv[0] is the name of the mode
};

/**
 * @brief prints the modes of a tool
 * @param out the stream to print to
 * @param tool name of the tool
 * @param description one line describing the tool
 * @param modes the modes
 */
inline void print_cli_modes(std::ostream& out, const char* tool, const char* description,
    const std::vector<cli_mode_t>& modes)
{
    out << tool << ": " << description << "\n\n";
    out << "usage: " << tool << " <mode> [...]\n";

    // the column in which the descriptions of the modes start
    uint64_t width = 0;
    for (const cli_mode_t& mode : modes) width = std::max<uint64_t>(width, std::strlen(mode.name));

    for (const cli_mode_t& mode : modes)
        out << "   " << std::left << std::setw(width + 2) << mode.name << mode.description << '\n';

    out << "\nrun '" << tool << " <mode>' without further arguments for the options of a mode" << std::endl;
}

/**
 * @brief runs the mode named by the first argument; without arguments or with -h or --help, the modes are printed
 *        (exit code 0), for an unknown mode an error and the modes are printed to stderr (exit code 1)
 * @param tool name of the tool
 * @param description one line describing the tool
 * @param modes the modes
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments
 * @return the exit code of the mode
 */
inline int run_cli_mode(const char* tool, const char* description, const std::vector<cli_mode_t>& modes, int argc,
    char** argv)
{
    if (argc < 2 || std::strcmp(argv[1], "-h") == 0 || std::strcmp(argv[1], "--help") == 0) {
        print_cli_modes(std::cout, tool, description, modes);
        return 0;
    }

    for (const cli_mode_t& mode : modes)
        if (std::strcmp(argv[1], mode.name) == 0) return mode.run(argc - 1, argv + 1);

    std::cerr << "error: unknown mode " << argv[1] << "\n\n";
    print_cli_modes(std::cerr, tool, description, modes);
    return 1;
}

#endif // MOVE_R_MISC_CLI_MODES_HPP
