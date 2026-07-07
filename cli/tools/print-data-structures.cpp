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

#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <lzendsa/lzendsa.hpp>
#include <move_data_structure/move_data_structure.hpp>
#include <move_r/move_r.hpp>
#include <move_rb/move_rb.hpp>
#include <rlzsa/rlzsa.hpp>

std::string structure;                  // the data structure to build
std::string input;                      // the input (text, or interval sequence for mds)
bool input_provided = false;            // whether an input was given via -t, -f or stdin
int64_t n_mds = -1;                     // the universe size n (only for mds)
int64_t d = -1;                         // sampling parameter delta (rlzsa / lzendsa; -1 => auto)
int64_t h = 8192;                       // maximum phrase length (lzendsa)
uint16_t a = 8;                         // balancing parameter (move_r / move_rb / mds)
uint16_t p = omp_get_max_threads();     // number of threads (move_r / move_rb)
bool verbose = false;                   // whether to print construction log messages
move_r_support support = _locate_move;  // the support type (move_r / move_rb)

/**
 * @brief prints the usage information and exits
 * @param msg an optional error message printed before the usage information
 */
void help(std::string msg)
{
    if (msg != "") std::cout << msg << std::endl << std::endl;
    std::cout << "print-data-structures: builds a data structure from a (small) input and prints its contents." << std::endl << std::endl;
    std::cout << "usage: print-data-structures <structure> [options]" << std::endl << std::endl;
    std::cout << "   <structure>      the data structure to build; one of:" << std::endl;
    std::cout << "                      rlzsa      the old rlzsa" << std::endl;
    std::cout << "                      rlzsa_opt  the new rlzsa (optimized)" << std::endl;
    std::cout << "                      lzendsa    the lzendsa" << std::endl;
    std::cout << "                      mds        a move data structure" << std::endl;
    std::cout << "                      move_r     a move-r index" << std::endl;
    std::cout << "                      move_rb    a bidirectional move-r index" << std::endl << std::endl;
    std::cout << "   input (if none of -t/-f is given, the input is read from stdin):" << std::endl;
    std::cout << "   -t <input>       use <input> as the input; for mds this is a disjoint interval" << std::endl;
    std::cout << "                    sequence such as \"(0,1),(1,0)\" (any non-digit separators are allowed)" << std::endl;
    std::cout << "   -f <file>        read the input from <file>" << std::endl << std::endl;
    std::cout << "   options:" << std::endl;
    std::cout << "   -s <support>     support type (move_r / move_rb); default: locate_move" << std::endl;
    std::cout << "                      move_r:  count, locate_one, locate_move, locate_rlzsa" << std::endl;
    std::cout << "                      move_rb: count, locate_move, locate_rlzsa" << std::endl;
    std::cout << "   -n <integer>     universe size n (required for mds)" << std::endl;
    std::cout << "   -a <integer>     balancing parameter, a >= 2 (move_r / move_rb / mds); default: 8" << std::endl;
    std::cout << "   -d <integer>     sampling parameter delta (rlzsa / lzendsa); default: auto (~10%)" << std::endl;
    std::cout << "   -h <integer>     maximum phrase length (lzendsa); default: 8192" << std::endl;
    std::cout << "   -p <integer>     number of threads (move_r / move_rb); default: all" << std::endl;
    std::cout << "   -v               print construction log messages" << std::endl;
    exit(msg == "" ? 0 : 1);
}

/**
 * @brief parses the support type given with the -s option
 * @param s the support string
 */
void parse_support(const std::string& s)
{
    if      (s == "count")                   support = _count;
    else if (s == "locate_one")              support = _locate_one;
    else if (s == "locate_move")             support = _locate_move;
    else if (s == "locate_rlzsa")            support = _locate_rlzsa;
    else help("error: unknown support '" + s + "' provided with -s");
}

/**
 * @brief parses a disjoint interval sequence from str by extracting every non-negative integer and
 *        grouping the integers into consecutive pairs (p_0, q_0), (p_1, q_1), ...
 * @tparam pos_t the index integer type
 * @param str the string to parse
 * @return the disjoint interval sequence
 */
template <typename pos_t>
std::vector<std::pair<pos_t, pos_t>> parse_intervals(const std::string& str)
{
    std::vector<pos_t> values;
    uint64_t i = 0;

    while (i < str.size()) {
        if (std::isdigit((unsigned char) str[i])) {
            uint64_t value = 0;
            while (i < str.size() && std::isdigit((unsigned char) str[i])) {
                value = value * 10 + (str[i] - '0');
                i++;
            }
            values.push_back((pos_t) value);
        } else {
            i++;
        }
    }

    if (values.size() % 2 != 0)
        help("error: the interval sequence must contain an even number of integers (p, q pairs)");

    std::vector<std::pair<pos_t, pos_t>> intervals;
    intervals.reserve(values.size() / 2);
    for (uint64_t j = 0; j + 1 < values.size(); j += 2)
        intervals.push_back({ values[j], values[j + 1] });

    return intervals;
}

/**
 * @brief builds the old rlzsa from the input text and prints its contents (the decoded suffix array)
 * @tparam int_t the rlzsa integer type
 */
template <typename int_t>
void build_rlzsa()
{
    rlzsa<int_t> index(input, d, false, verbose);
    uint64_t n = index.input_size();

    std::cout << std::endl << "n = " << n << ", z = " << index.num_phrases()
        << ", size = " << format_size(index.size_in_bytes()) << std::endl;
    std::cout << "(the old rlzsa has no contents-logging facility; printing the decoded suffix array)" << std::endl;

    log_indexed("SA:", index.template sa_values<uint64_t>(0, n - 1));
}

/**
 * @brief builds the new rlzsa_opt (via a move_r<_locate_rlzsa> index) from the input text and prints
 *        its contents and data-structure sizes
 * @tparam pos_t the index integer type
 */
template <typename pos_t>
void build_rlzsa_opt()
{
    move_r<_locate_rlzsa, char, pos_t> index(input, { .num_threads = p, .a = a, .log = verbose });

    std::cout << "n = " << index.input_size()
        << ", z = " << index.num_phrases_rlzsa() << std::endl << std::endl;
    index.rlzsa().log_data_structures();
}

/**
 * @brief builds the lzendsa from the input text and prints its contents and data-structure sizes
 * @tparam int_t the lzendsa integer type
 */
template <typename int_t>
void build_lzendsa()
{
    lzendsa<int_t> index(input, d, h, false, verbose);

    std::cout << "n = " << index.input_size() << ", z = " << index.num_phrases()
        << ", h = " << h << ", size = " << format_size(index.size_in_bytes()) << std::endl << std::endl;
    index.sa_encoding().log_data_structures();
}

/**
 * @brief builds a move data structure from the parsed interval sequence and prints its contents
 * @tparam pos_t the index integer type
 */
template <typename pos_t>
void build_mds()
{
    std::vector<std::pair<pos_t, pos_t>> intervals = parse_intervals<pos_t>(input);

    if (intervals.empty()) help("error: the interval sequence is empty");
    if (n_mds < 0) help("error: <n> (the universe size) is required for mds; provide it with -n");

    move_data_structure<pos_t, POS> mds(std::move(intervals), (pos_t) n_mds,
        { .num_threads = p, .a = a });

    std::cout << "n = " << n_mds << ", k' = " << mds.num_intervals()
        << ", a = " << mds.balancing_parameter() << std::endl << std::endl;
    mds.log_data_structures();
}

/**
 * @brief builds a move_r index from the input text and prints its contents and data-structure sizes
 * @tparam pos_t the index integer type
 * @tparam supp the move-r support type
 */
template <typename pos_t, move_r_support supp>
void build_move_r()
{
    move_r<supp, char, pos_t> index(input, { .num_threads = p, .a = a, .log = verbose });
    index.log_data_structures();
}

/**
 * @brief builds a move_rb index from the input text and prints its contents and data-structure sizes
 * @tparam pos_t the index integer type
 * @tparam supp the move-r support type
 */
template <typename pos_t, move_r_support supp>
void build_move_rb()
{
    move_rb<supp, char, pos_t> index(input, { .num_threads = p, .a = a, .log = verbose });
    index.log_data_structures();
}

/**
 * @brief dispatches a move_r build to the right (pos_t, support) instantiation
 * @param use_64_bit whether to use 64-bit positions
 */
void dispatch_move_r(bool use_64_bit)
{
    auto run = [&]<move_r_support supp>() {
        if (use_64_bit) build_move_r<uint64_t, supp>();
        else            build_move_r<uint32_t, supp>();
    };

    switch (support) {
        case _count:                   run.template operator()<_count>(); break;
        case _locate_one:              run.template operator()<_locate_one>(); break;
        case _locate_move:             run.template operator()<_locate_move>(); break;
        case _locate_rlzsa:            run.template operator()<_locate_rlzsa>(); break;
        default:                       help("error: unsupported support type for move_r");
    }
}

/**
 * @brief dispatches a move_rb build to the right (pos_t, support) instantiation; move_rb only
 *        supports count, locate_move and locate_rlzsa
 * @param use_64_bit whether to use 64-bit positions
 */
void dispatch_move_rb(bool use_64_bit)
{
    auto run = [&]<move_r_support supp>() {
        if (use_64_bit) build_move_rb<uint64_t, supp>();
        else            build_move_rb<uint32_t, supp>();
    };

    switch (support) {
        case _count:        run.template operator()<_count>(); break;
        case _locate_move:  run.template operator()<_locate_move>(); break;
        case _locate_rlzsa: run.template operator()<_locate_rlzsa>(); break;
        default:            help("error: move_rb only supports count, locate_move and locate_rlzsa");
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
    if (argc < 2) help("");
    structure = argv[1];
    if (structure == "-h" || structure == "--help") help("");

    int i = 2;
    while (i < argc) {
        std::string s = argv[i++];
        auto value = [&](const std::string& opt) -> std::string {
            if (i >= argc) help("error: missing parameter after " + opt);
            return argv[i++];
        };

        if (s == "-t") {
            input = value("-t");
            input_provided = true;
        } else if (s == "-f") {
            std::string path = value("-f");
            std::ifstream ifs(path);
            if (!ifs.good()) help("error: cannot open input file '" + path + "'");
            input = std::string(std::istreambuf_iterator<char>(ifs), {});
            input_provided = true;
        } else if (s == "-s") {
            parse_support(value("-s"));
        } else if (s == "-n") {
            n_mds = std::stoll(value("-n"));
        } else if (s == "-a") {
            a = std::stoi(value("-a"));
            if (a < 2) help("error: a < 2");
        } else if (s == "-d") {
            d = std::stoll(value("-d"));
        } else if (s == "-h") {
            h = std::stoll(value("-h"));
        } else if (s == "-p") {
            p = std::stoi(value("-p"));
            if (p < 1) help("error: p < 1");
            if (p > omp_get_max_threads()) help("error: p > maximum number of threads");
        } else if (s == "-v") {
            verbose = true;
        } else {
            help("error: unrecognized option '" + s + "'");
        }
    }

    // if no input was provided via -t or -f, read it from stdin
    if (!input_provided) {
        input = std::string(std::istreambuf_iterator<char>(std::cin), {});
        // strip a single trailing newline that interactive input usually adds
        if (!input.empty() && input.back() == '\n') input.pop_back();
    }

    std::cout << std::setprecision(4);

    if (structure == "mds") {
        bool use_64_bit = n_mds >= UINT_MAX;
        if (use_64_bit) build_mds<uint64_t>();
        else            build_mds<uint32_t>();
        return 0;
    }

    if (input.empty()) help("error: the input is empty");

    // for small inputs the automatic SA-sampling rate (rlzsa / lzendsa) can round down to 0; for an
    // inspection tool on small inputs, sampling every SA value (d = 1) is the sensible default
    if (d < 0 && input.size() < 10000) d = 1;

    if (structure == "rlzsa") {
        bool use_64_bit = input.size() > INT_MAX;
        if (use_64_bit) build_rlzsa<int64_t>();
        else            build_rlzsa<int32_t>();
    } else if (structure == "lzendsa") {
        bool use_64_bit = input.size() > INT_MAX;
        if (use_64_bit) build_lzendsa<int64_t>();
        else            build_lzendsa<int32_t>();
    } else if (structure == "rlzsa_opt") {
        bool use_64_bit = (input.size() + 1) >= UINT_MAX;
        if (use_64_bit) build_rlzsa_opt<uint64_t>();
        else            build_rlzsa_opt<uint32_t>();
    } else if (structure == "move_r") {
        dispatch_move_r((input.size() + 1) >= UINT_MAX);
    } else if (structure == "move_rb") {
        dispatch_move_rb((input.size() + 1) >= UINT_MAX);
    } else {
        help("error: unknown structure '" + structure + "'");
    }

    return 0;
}
