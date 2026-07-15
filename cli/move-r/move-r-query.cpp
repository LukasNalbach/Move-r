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
#include <iostream>
#include <string>
#include <ips2ra.hpp>
#include <move_r/move_r.hpp>

int32_t arg_idx = 1;
std::string path_index_file;
std::ifstream index_file;

/**
 * @brief prints the usage information and exits
 * @param msg an optional error message printed before the usage information
 */
void help(std::string msg)
{
    if (msg != "") std::cout << msg << std::endl << std::endl;
    std::cout << "move-r-query: loads a move-r index and answers count/locate queries typed interactively." << std::endl << std::endl;
    std::cout << "usage: move-r-query [options] <index_file>" << std::endl;
    std::cout << "   <index_file>      index file (with extension .move-r)" << std::endl << std::endl;
    std::cout << "interactive commands (read from stdin, one per line):" << std::endl;
    std::cout << "   count <pattern>   prints the number of occurrences of <pattern>" << std::endl;
    std::cout << "   locate <pattern>  prints the (sorted) occurrences of <pattern>" << std::endl;
    std::cout << "   help              prints the list of commands" << std::endl;
    std::cout << "   quit / exit       leaves the program" << std::endl;
    std::cout << "(the pattern is the remainder of the line and may contain spaces)" << std::endl;
    exit(msg == "" ? 0 : 1);
}

/**
 * @brief prints the list of interactive commands
 */
void print_commands() { std::cout << "commands: count <pattern> | locate <pattern> | help | quit" << std::endl; }

/**
 * @brief loads the index and runs the interactive query loop
 * @tparam pos_t index integer type
 * @tparam support the move-r locate-support type
 */
template <typename pos_t, move_r_support support>
void run()
{
    constexpr bool can_locate = support != _count && support != _locate_one;

    using idx_t = move_r<support, char, pos_t>;
    idx_t index;
    index.load(index_file);
    index_file.close();

    std::cout << "loaded index for an input of " << index.input_size() << " symbols ("
        << index.num_bwt_runs() << " BWT runs, " << format_size(index.size_in_bytes()) << ")." << std::endl;
    if constexpr (!can_locate)
        std::cout << "this index only supports count." << std::endl;
    print_commands();

    std::string line;
    while (true) {
        std::cout << "> " << std::flush;
        if (!std::getline(std::cin, line)) break;
        if (line.empty()) continue;

        uint64_t sep = line.find(' ');
        std::string command = line.substr(0, sep);
        std::string pattern = sep == std::string::npos ? "" : line.substr(sep + 1);

        if (command == "quit" || command == "exit") {
            break;
        } else if (command == "help") {
            print_commands();
        } else if (command == "count") {
            if (pattern.empty()) { std::cout << "error: missing pattern" << std::endl; continue; }
            auto t = now();
            pos_t c = index.count(pattern);
            std::cout << c << " occurrence(s) (" << format_time(time_diff_ns(t, now())) << ")" << std::endl;
        } else if (command == "locate") {
            if constexpr (can_locate) {
                if (pattern.empty()) { std::cout << "error: missing pattern" << std::endl; continue; }
                auto t = now();
                std::vector<pos_t> occ = index.locate(pattern);
                ips2ra::sort(occ.begin(), occ.end());
                std::cout << occ.size() << " occurrence(s) (" << format_time(time_diff_ns(t, now())) << ")";
                if (!occ.empty()) std::cout << ":";
                for (pos_t o : occ) std::cout << " " << o;
                std::cout << std::endl;
            } else {
                std::cout << "error: this index does not support locate" << std::endl;
            }
        } else {
            std::cout << "error: unknown command '" << command << "'" << std::endl;
            print_commands();
        }
    }
}

/**
 * @brief parses the next command-line argument(s)
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments
 */
void parse_args(char** argv, int32_t /*argc*/)
{
    std::string s = argv[arg_idx++];
    help("error: unrecognized '" + s + "' option");
}

/**
 * @brief program entry point
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments
 * @return the exit code
 */
int main(int argc, char** argv)
{
    if (argc < 2) help("");
    while (arg_idx < argc - 1) parse_args(argv, argc);

    path_index_file = argv[arg_idx];
    index_file.open(path_index_file);
    if (!index_file.good()) help("error: could not read <index_file>");

    bool is_64_bit;
    index_file.read((char*) &is_64_bit, 1);
    move_r_support support;
    index_file.read((char*) &support, sizeof(move_r_support));
    index_file.seekg(0, std::ios::beg);

    if (support == _count) {
        if (is_64_bit) run<uint64_t, _count>();
        else           run<uint32_t, _count>();
    } else if (support == _locate_one) {
        if (is_64_bit) run<uint64_t, _locate_one>();
        else           run<uint32_t, _locate_one>();
    } else if (support == _locate_move) {
        if (is_64_bit) run<uint64_t, _locate_move>();
        else           run<uint32_t, _locate_move>();
    } else if (support == _locate_rlzsa) {
        if (is_64_bit) run<uint64_t, _locate_rlzsa>();
        else           run<uint32_t, _locate_rlzsa>();
    } else {
        help("error: unknown support type in the index file");
    }

    return 0;
}
