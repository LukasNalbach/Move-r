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
 * @file print-data-structures.cpp
 * @brief the tool print-data-structures: builds a data structure from a small input and prints its contents
 */

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <omp.h>

#include <misc/cli_help.hpp>
#include <misc/log_structures.hpp>
#include <misc/move_r_cli.hpp>
#include <misc/parse.hpp>
#include <move_data_structure/move_data_structure.hpp>
#include <move_r/move_r.hpp>
#include <move_rb/move_rb.hpp>

namespace {

// the number of significant digits of the printed numbers
constexpr int printed_significant_digits = 4;

// the data structures of <data_structure>, in the order of the help
const std::vector<std::string> data_structure_names = { "mds", "sa_rlz", "move-r", "move-rb", "context_tables" };

/** @brief the options of print-data-structures */
struct print_cli_t {
    std::string data_structure; // <data_structure> (see data_structure_names)
    std::string input; // the input: a text, or for mds a disjoint interval sequence (-t, -f or stdin)
    bool has_input = false; // whether the input was given with -t or -f
    uint64_t mds_n = 0; // -n: the universe size n of the move data structure (0: not given)
    uint64_t a = min_balancing_parameter; // -balance: the balancing parameter
    uint64_t num_threads = omp_get_max_threads(); // -p
    uint64_t rlz_max_reference_len = move_r_params {}.rlz_max_reference_len; // -rlz-ref
    move_r_support support = _locate_move; // -s: the support type (move-r and move-rb)
    uint64_t max_depth = default_logged_tree_depth; // -depth: the depth up to which trees are printed
    uint64_t max_entries = default_logged_entries; // -entries: the maximum number of printed entries per table or tree
    // -context-table: the maximum size of the context tables of move-rb and context_tables as a share of the index
    double context_table_share = 1;
    bool context_tables_single_length = false; // -context-table-one-length
} cli;

/**
 * @brief prints the help of print-data-structures
 * @param out the stream to print to
 */
void print_help(std::ostream& out)
{
    const print_cli_t defaults;
    help_printer h(out);
    h.intro("print-data-structures: builds a data structure from a small input and prints its contents.",
        "print-data-structures <data_structure> [options]");
    h.group("arguments");
    h.option("<data_structure>", "mds | sa_rlz | move-r | move-rb | context_tables",
        "mds: a move data structure of a disjoint interval sequence\n"
        "sa_rlz: the RLZ encoding of the differential suffix array SA^d (of move-r -s locate-rlzsa)\n"
        "move-r, move-rb: a move-r index or a bidirectional move-r index (move-rb: with its context tables)\n"
        "context_tables: the context tables of a move-rb index");
    h.group("input");
    h.option("-t <s>", "a string", "the input; for mds the pairs (p_0, q_0), (p_1, q_1), ... of a disjoint interval "
        "sequence, e.g. \"(0,1),(1,0)\" (all symbols other than digits separate the integers)");
    h.option("-f <file>", "a path", "read the input from this file");
    h.group("data structure");
    h.option("-s <support>", "count | locate-one | locate-move | locate-rlzsa (default: "
        + move_r_support_option_name(defaults.support) + "); move-rb: count, locate-move or locate-rlzsa",
        "the support type (move-r and move-rb)");
    h.option("-n <n>", "integer >= 1 (required for mds)", "the universe size of the move data structure: its "
        "intervals partition [0..n-1]");
    h.option("-balance <a>", "integer >= " + std::to_string(min_balancing_parameter) + " (default: "
        + std::to_string(defaults.a) + ")", "the balancing parameter of the move data structures");
    h.option("-rlz-ref <len>", "integer >= 0 (default: " + std::to_string(defaults.rlz_max_reference_len)
        + ", chosen by the construction)", "sa_rlz: the number of values of SA^d in the reference");
    h.option("-context-table <share>", "real >= 0 (default: " + help_value(defaults.context_table_share) + ")",
        "move-rb and context_tables: the maximum size of the context tables as a share of the size of the index");
    h.option("-context-table-one-length", "", "move-rb and context_tables: store the context tables for one pattern "
        "length only");
    h.group("output");
    h.option("-depth <d>", "integer >= 0 (default: " + std::to_string(defaults.max_depth) + ")",
        "the depth up to which trees are printed (the roots have depth 0)");
    h.option("-entries <k>", "integer >= 1 (default: " + std::to_string(defaults.max_entries) + ")",
        "the maximum number of printed values or nodes per table or tree");
    h.group("execution");
    h.option_threads();
    h.text("Without -t and -f, the input is read from stdin (without a final line break).");
}

/**
 * @brief returns the non-negative integers of a string in their order; all symbols other than digits separate them
 * @param str the string
 */
std::vector<uint64_t> parse_integers(const std::string& str)
{
    std::vector<uint64_t> values;
    constexpr uint64_t decimal_base = 10;
    uint64_t i = 0;

    while (i < str.size()) {
        if (std::isdigit((uint8_t) str[i])) {
            uint64_t value = 0;

            while (i < str.size() && std::isdigit((uint8_t) str[i])) {
                value = value * decimal_base + (str[i] - '0');
                i++;
            }

            values.push_back(value);
        } else {
            i++;
        }
    }

    return values;
}

/**
 * @brief parses a disjoint interval sequence: groups the integers of str (see parse_integers) into consecutive pairs
 *        (p_0, q_0), (p_1, q_1), ...
 * @param str the string to parse
 * @return the pairs (p_i, q_i)
 * @throws std::runtime_error if the number of integers is odd
 */
std::vector<std::pair<uint64_t, uint64_t>> parse_intervals(const std::string& str)
{
    // the integers of str in their order (p_0, q_0, p_1, q_1, ...)
    const std::vector<uint64_t> values = parse_integers(str);

    if (values.size() % 2 != 0)
        throw std::runtime_error("the interval sequence must contain an even number of integers (pairs (p_i, q_i))");

    std::vector<std::pair<uint64_t, uint64_t>> intervals;
    intervals.reserve(values.size() / 2);
    for (uint64_t j = 0; j + 1 < values.size(); j += 2) intervals.emplace_back(values[j], values[j + 1]);
    return intervals;
}

/**
 * @brief reads the value of an option as a decimal number >= 0 (see parse_double)
 * @param args the arguments, positioned after the option
 * @param option the option (for the error message)
 */
double read_share(cli_arguments& args, const std::string& option)
{
    double value = 0;
    if (!parse_double(args.value(option), value) || !(value >= 0)) throw usage_error(option + " expects a number >= 0");
    return value;
}

// the construction parameters of the indexes
move_r_params index_params()
{
    return { .num_threads = cli.num_threads, .a = cli.a, .rlz_max_reference_len = cli.rlz_max_reference_len };
}

// prints a line "<name> = <value>"
void print_fact(const std::string& name, uint64_t value) { std::cout << name << " = " << value << std::endl; }

// builds a move_r<_locate_rlzsa> of the input and prints the RLZ encoding of SA^d
void print_sa_rlz()
{
    move_r<_locate_rlzsa> index(cli.input, index_params());
    print_fact("n", index.input_size() + 1);
    index.SA_rlz().log_data_structures();
}

// builds a move data structure of the interval sequence of the input and prints it
void print_mds()
{
    const std::vector<std::pair<uint64_t, uint64_t>> intervals = parse_intervals(cli.input);
    if (intervals.empty()) throw std::runtime_error("the interval sequence is empty");
    if (cli.mds_n == 0) throw usage_error("mds requires -n <n>");
    uint64_t max_value = cli.mds_n;
    for (const auto& [p_i, q_i] : intervals) max_value = std::max({ max_value, p_i, q_i });
    mds_pairs pairs(intervals.size(), max_value);

    for (uint64_t i = 0; i < intervals.size(); i++) {
        pairs.set<0>(i, intervals[i].first);
        pairs.set<1>(i, intervals[i].second);
    }

    move_data_structure<POS> mds(std::move(pairs), cli.mds_n, { .num_threads = cli.num_threads, .a = cli.a });
    print_fact("n", cli.mds_n);
    print_fact("k'", mds.num_intervals());
    print_fact("a", mds.balancing_parameter());
    std::cout << std::endl;
    mds.log_data_structures();
}

/**
 * @brief builds a move_r of the input and prints its contents
 * @tparam support the support type
 */
template <move_r_support support>
void print_move_r()
{
    move_r<support> index(cli.input, index_params());
    index.log_data_structures();
}

/**
 * @brief builds a move_rb of the input and prints its contents
 * @tparam support the support type
 */
template <move_r_support support>
void print_move_rb()
{
    move_r_params params = index_params();
    params.context_table_share = cli.context_table_share;
    params.context_tables_single_length = cli.context_tables_single_length;
    move_rb<support> index(cli.input, params);

    if (cli.data_structure == "context_tables") {
        index.log_context_tables(cli.max_depth, cli.max_entries);
        return;
    }

    index.log_data_structures();
    std::cout << std::endl << "context tables:" << std::endl;
    index.log_context_tables(cli.max_depth, cli.max_entries);
}

// builds a move_r with the support type of -s and prints it
void dispatch_move_r()
{
    switch (cli.support) {
        case _count:        print_move_r<_count>(); break;
        case _locate_one:   print_move_r<_locate_one>(); break;
        case _locate_move:  print_move_r<_locate_move>(); break;
        case _locate_rlzsa: print_move_r<_locate_rlzsa>(); break;
        default: throw usage_error("move-r supports -s count, locate-one, locate-move and locate-rlzsa");
    }
}

// builds a move_rb with the support type of -s and prints it
void dispatch_move_rb()
{
    switch (cli.support) {
        case _count:        print_move_rb<_count>(); break;
        case _locate_move:  print_move_rb<_locate_move>(); break;
        case _locate_rlzsa: print_move_rb<_locate_rlzsa>(); break;
        default: throw usage_error("move-rb supports -s count, locate-move and locate-rlzsa");
    }
}

/**
 * @brief reads the arguments into cli
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments
 */
void read_arguments(int argc, char** argv)
{
    cli_arguments args(argc, argv);
    if (args.at_option()) throw usage_error("expected <data_structure> before the options");
    cli.data_structure = args.next();

    while (!args.empty()) {
        const std::string s = args.next();

        if (s == "-t") {
            cli.input = args.value("-t");
            cli.has_input = true;
        } else if (s == "-f") {
            const std::string path = args.value("-f");
            direct_ifstream file(path);
            require_readable(file, path);
            cli.input = std::string(std::istreambuf_iterator<char>(file), {});
            file.close();
            drop_from_page_cache(path, true);
            cli.has_input = true;
        } else if (s == "-s") {
            const std::string name = args.value("-s");
            if (!parse_move_r_support(name, cli.support)) throw usage_error("unknown support type " + name);
        } else if (s == "-n") {
            cli.mds_n = args.integer("-n", 1);
        } else if (s == "-balance") {
            cli.a = args.integer("-balance", min_balancing_parameter);
        } else if (s == "-p") {
            cli.num_threads = args.integer("-p", 1, omp_get_max_threads());
        } else if (s == "-rlz-ref") {
            cli.rlz_max_reference_len = args.integer("-rlz-ref");
        } else if (s == "-depth") {
            cli.max_depth = args.integer("-depth");
        } else if (s == "-entries") {
            cli.max_entries = args.integer("-entries", 1);
        } else if (s == "-context-table") {
            cli.context_table_share = read_share(args, "-context-table");
        } else if (s == "-context-table-one-length") {
            cli.context_tables_single_length = true;
        } else {
            throw unknown_option(s);
        }
    }

    const std::string& name = cli.data_structure;
    if (std::find(data_structure_names.begin(), data_structure_names.end(), name) == data_structure_names.end())
        throw usage_error("unknown data structure " + name);
}

} // namespace

int main(int argc, char** argv)
{
    return run_cli(argc, argv, print_help, [&] {
        read_arguments(argc, argv);

        if (!cli.has_input) {
            cli.input = std::string(std::istreambuf_iterator<char>(std::cin), {});
            if (!cli.input.empty() && cli.input.back() == '\n') cli.input.pop_back();
        }

        std::cout << std::setprecision(printed_significant_digits);

        if (cli.data_structure == "mds") {
            print_mds();
            return;
        }

        if (cli.input.empty()) throw std::runtime_error("the input is empty");
        const std::string& name = cli.data_structure;
        if (name == "sa_rlz") print_sa_rlz();
        else if (name == "move-r") dispatch_move_r();
        else dispatch_move_rb();
    });
}
