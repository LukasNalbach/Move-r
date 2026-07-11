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
#include <move_rb/move_rb.hpp>

int arg_idx = 1;
std::string path_index_file;
std::ifstream index_file;

// maximum number of occurrences printed per query (0 = unlimited); changeable with the "limit <n>" command
uint64_t print_limit = 100;

/**
 * @brief prints the usage information and exits
 * @param msg an optional error message printed before the usage information
 */
void help(std::string msg)
{
    if (msg != "") std::cout << msg << std::endl << std::endl;
    std::cout << "move-rb-query: loads a move-rb index and answers exact and approximate pattern-matching" << std::endl;
    std::cout << "               queries typed interactively." << std::endl << std::endl;
    std::cout << "usage: move-rb-query <index_file>" << std::endl;
    std::cout << "   <index_file>   index file (with extension .move-rb)" << std::endl << std::endl;
    std::cout << "interactive commands (read from stdin, one per line):" << std::endl;
    std::cout << "   count <pattern>            exact number of occurrences of <pattern>" << std::endl;
    std::cout << "   locate <pattern>           exact (sorted) occurrences of <pattern>" << std::endl;
    std::cout << "   hamming-count <k> <pat>    number of occurrences with <= k mismatches (hamming)" << std::endl;
    std::cout << "   hamming-locate <k> <pat>   occurrences with <= k mismatches (hamming)" << std::endl;
    std::cout << "   edit-locate <k> <pat>      occurrences with <= k errors (edit distance, k <= 20)" << std::endl;
    std::cout << "   limit <n>                  print at most <n> occurrences per query (0 = unlimited)" << std::endl;
    std::cout << "   help                       prints the list of commands" << std::endl;
    std::cout << "   quit / exit                leaves the program" << std::endl << std::endl;
    std::cout << "   approximate queries use the min_u search scheme." << std::endl;
    std::cout << "   <pattern> / <pat> is the remainder of the line and may contain spaces" << std::endl;
    exit(msg == "" ? 0 : 1);
}

/**
 * @brief prints the list of interactive commands
 */
void print_commands()
{
    std::cout << "commands: count <pat> | locate <pat> | hamming-count <k> <pat> |" << std::endl;
    std::cout << "          hamming-locate <k> <pat> | edit-locate <k> <pat> | limit <n> | help | quit" << std::endl;
    std::cout << "approximate queries use the min_u search scheme." << std::endl;
}

/**
 * @brief splits a query line into a fixed number of leading tokens followed by the pattern (the
 *        remainder of the line, which may contain spaces)
 * @param line the query line (without the leading command word)
 * @param num_tokens the number of leading tokens to extract
 * @param tokens output vector receiving the extracted tokens
 * @param pattern output string receiving the remainder of the line
 * @return whether num_tokens tokens and a non-empty pattern could be extracted
 */
bool split_query(const std::string& line, int num_tokens, std::vector<std::string>& tokens, std::string& pattern)
{
    uint64_t i = 0;
    tokens.clear();

    for (int t = 0; t < num_tokens; t++) {
        while (i < line.size() && line[i] == ' ') i++;
        uint64_t start = i;
        while (i < line.size() && line[i] != ' ') i++;
        if (start == i) return false;
        tokens.push_back(line.substr(start, i - start));
    }

    while (i < line.size() && line[i] == ' ') i++;
    pattern = line.substr(i);
    return !pattern.empty();
}

/**
 * @brief computes the exact number of occurrences of P (extending P to the right symbol by symbol)
 * @tparam pos_t index integer type
 * @tparam support the move-rb locate-support type
 * @param index the index to query
 * @param P the pattern to count
 * @return the number of occurrences of P in the input
 */
template <typename pos_t, move_r_support support>
pos_t exact_count(const move_rb<support, char, pos_t>& index, const std::string& P)
{
    auto ctx = index.template empty_context<COUNT>();
    for (char c : P) {
        auto [next, found] = ctx.extend(c, RIGHT);
        if (!found) return 0;
        ctx = next;
    }
    return ctx.num_occ();
}

/**
 * @brief loads the index and runs the interactive query loop
 * @tparam pos_t index integer type
 * @tparam support the move-rb locate-support type
 */
template <typename pos_t, move_r_support support>
void run()
{
    constexpr bool can_locate = support != _count;

    using idx_t = move_rb<support, char, pos_t>;
    idx_t index;
    index.load(index_file);
    index_file.close();

    std::cout << "loaded index for an input of " << index.forward_index().input_size() << " symbols ("
        << format_size(index.size_in_bytes()) << ")." << std::endl;
    if constexpr (!can_locate)
        std::cout << "this index only supports count (exact count and hamming-count)." << std::endl;
    print_commands();

    std::string line, pattern;
    std::vector<std::string> tokens;

    while (true) {
        std::cout << "> " << std::flush;
        if (!std::getline(std::cin, line)) break;
        if (line.empty()) continue;

        uint64_t sep = line.find(' ');
        std::string command = line.substr(0, sep);
        std::string rest = sep == std::string::npos ? "" : line.substr(sep + 1);

        if (command == "quit" || command == "exit") {
            break;
        } else if (command == "help") {
            print_commands();
            continue;
        } else if (command == "limit") {
            try {
                print_limit = std::stoull(rest);
            } catch (...) {
                std::cout << "error: usage: limit <n>  (0 = unlimited)" << std::endl;
                continue;
            }
            if (print_limit == 0) std::cout << "occurrence output is now unlimited" << std::endl;
            else std::cout << "occurrence output limited to " << print_limit << " per query" << std::endl;
            continue;
        }

        // exact queries: command <pattern>
        if (command == "count" || command == "locate") {
            pattern = rest;
            if (pattern.empty()) { std::cout << "error: missing pattern" << std::endl; continue; }
            auto t = now();

            if (command == "count") {
                pos_t c = exact_count<pos_t, support>(index, pattern);
                std::cout << c << " occurrence(s) (" << format_time(time_diff_ns(t, now())) << ")" << std::endl;
            } else {
                if constexpr (can_locate) {
                    auto ctx = index.template empty_context<LOCATE>();
                    bool found = true;
                    for (char c : pattern) {
                        auto [next, ok] = ctx.extend(c, RIGHT);
                        if (!ok) { found = false; break; }
                        ctx = next;
                    }
                    std::vector<pos_t> occ;
                    if (found) occ = ctx.locate_phase().locate();
                    ips2ra::sort(occ.begin(), occ.end());
                    uint64_t to_print = print_limit == 0 ? occ.size() : std::min<uint64_t>(occ.size(), print_limit);
                    std::cout << occ.size() << " occurrence(s) (" << format_time(time_diff_ns(t, now())) << ")";
                    if (!occ.empty()) std::cout << ":";
                    for (uint64_t i = 0; i < to_print; i++) std::cout << " " << occ[i];
                    if (to_print < occ.size()) std::cout << " ... (" << (occ.size() - to_print) << " more)";
                    std::cout << std::endl;
                } else {
                    std::cout << "error: this index does not support locate" << std::endl;
                }
            }
            continue;
        }

        // approximate queries: command <k> <scheme> <pattern>
        if (command == "hamming-count" || command == "hamming-locate" || command == "edit-locate") {
            if (!split_query(rest, 1, tokens, pattern)) {
                std::cout << "error: usage: " << command << " <k> <pattern>" << std::endl;
                continue;
            }

            int64_t k;
            try {
                k = std::stoll(tokens[0]);
            } catch (...) {
                std::cout << "error: <k> must be an integer" << std::endl;
                continue;
            }
            if (k < 0) { std::cout << "error: <k> must be >= 0" << std::endl; continue; }
            if (command == "edit-locate" && k > 20) { std::cout << "error: k > 20 is not supported for edit distance" << std::endl; continue; }
            if (k >= 256) { std::cout << "error: k must be < 256" << std::endl; continue; }

            search_scheme_t scheme = min_u_scheme(k);
            if (pattern.size() < scheme.p) {
                std::cout << "error: pattern length (" << pattern.size() << ") < number of parts in the search scheme ("
                    << (int) scheme.p << ")" << std::endl;
                continue;
            }

            auto t = now();

            if (command == "hamming-count") {
                pos_t c = index.count_hamming_dist(pattern, scheme);
                std::cout << c << " occurrence(s) (" << format_time(time_diff_ns(t, now())) << ")" << std::endl;
            } else if constexpr (can_locate) {
                std::vector<aprx_occ_t<pos_t>> occ;
                if (command == "hamming-locate") {
                    index.template locate<HAMMING_DISTANCE>(pattern, scheme, [&](aprx_occ_t<pos_t> o){ occ.emplace_back(o); });
                } else {
                    index.template locate<EDIT_DISTANCE>(pattern, scheme, [&](aprx_occ_t<pos_t> o){ occ.emplace_back(o); });
                    filter_edit_distance_occurrences<pos_t>(occ, (pos_t) k);
                }
                ips2ra::sort(occ.begin(), occ.end(), [](const auto& o){ return o.pos; });
                uint64_t to_print = print_limit == 0 ? occ.size() : std::min<uint64_t>(occ.size(), print_limit);
                std::cout << occ.size() << " occurrence(s) (" << format_time(time_diff_ns(t, now())) << ")";
                if (!occ.empty()) std::cout << ":";
                for (uint64_t i = 0; i < to_print; i++) {
                    const aprx_occ_t<pos_t>& o = occ[i];
                    std::cout << " (pos=" << o.pos << " len=" << o.len << " err=" << o.err << ")";
                }
                if (to_print < occ.size()) std::cout << " ... (" << (occ.size() - to_print) << " more)";
                std::cout << std::endl;
            } else {
                std::cout << "error: this index does not support locate" << std::endl;
            }
            continue;
        }

        std::cout << "error: unknown command '" << command << "'" << std::endl;
        print_commands();
    }
}

/**
 * @brief program entry point
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments
 * @return the exit code
 */
int main(int argc, char** argv)
{
    if (argc != 2) help("");
    std::string first = argv[1];
    if (first == "-h" || first == "--help") help("");

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
