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
 * @file docs.cpp
 * @brief the mode docs of move-rb
 */

#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

#include <misc/cli_help.hpp>
#include <misc/docs_query_cli.hpp>
#include <misc/index_cli.hpp>
#include <misc/move_rb_cli.hpp>
#include <misc/query_cli.hpp>
#include <move_rb/move_rb.hpp>

#include "modes.hpp"

namespace {

/** @brief the options of move-rb docs */
struct docs_cli_t {
    query_options_t options; // -p, -o, -m
    docs_query_cli_t docs; // -report, -top
    std::string path_index_file;
    std::string path_patterns_file;
    direct_ofstream output_file; // -o
} cli;

/**
 * @brief prints the help of move-rb docs
 * @param out the stream to print to
 */
void print_help(std::ostream& out)
{
    help_printer h(out);
    h.intro("move-rb docs: finds the documents that contain each pattern.",
        "move-rb docs [options] <index_file> <patterns_file>");
    h.group("arguments");
    h.option("<index_file>", "", "an index built with move-rb build -docs (with -s locate-move or locate-rlzsa)");
    h.option("<patterns_file>", "", patterns_file_help);
    h.group("what is computed");
    print_docs_query_options(h);
    h.group("output");
    print_docs_output_option(h);
    h.group("execution");
    h.option_threads();
    h.group("measurement");
    h.option_measurement();
}

/**
 * @brief reads an option of move-rb docs
 * @param args the arguments, positioned after the option
 * @param option the option
 * @throws usage_error if the option is unknown or its value is invalid
 */
void read_option(cli_arguments& args, const std::string& option)
{
    if (read_query_option(args, option, cli.options, _option_threads | _option_output | _option_measurement)) return;
    if (!read_docs_query_option(args, option, cli.docs)) throw unknown_option(option);
}

/**
 * @brief loads the index, finds the documents of the patterns of the patterns file and prints the summary
 * @tparam support the locate support of the index
 * @param index_file the index file, positioned at the beginning of the index
 */
template <move_r_support support>
void run_queries(direct_ifstream& index_file)
{
    move_rb<support> index;
    load_index_with_progress(index, index_file, "move-rb");
    index_file.close();
    drop_from_page_cache(cli.path_index_file, true);

    if (!index.has_documents())
        throw std::runtime_error(cli.path_index_file + " has no documents (build it with move-rb build -docs)");

    const auto& fwd = index.forward_index();
    const docs_query_run_t run { .options = cli.options, .output_file = cli.output_file,
        .path_index_file = cli.path_index_file, .path_patterns_file = cli.path_patterns_file,
        .measurement_name = "move_rb_" + std::string(move_r_support_name(support)) };

    run_docs_queries(index.documents(), fwd.input_size(), cli.docs, run,
        [&](const std::string& P, auto&& found) {
            auto ctx = fwd.query();
            bool occurs = true;
            for (uint64_t i = P.size(); occurs && i-- > 0;) occurs = ctx.prepend(P[i]);
            if (!occurs) return;
            const auto [b, e] = ctx.sa_interval();
            found(b, e, [&](auto&& report) { ctx.locate([&](uint64_t, uint64_t occ) { report(occ); }); });
        },
        [&](std::ostream& out) { index.log_data_structure_sizes(out); });
}

} // namespace

int move_rb_docs(int argc, char** argv)
{
    return run_cli(argc, argv, print_help, [&] {
        cli_arguments args(argc, argv);
        while (args.at_option()) read_option(args, args.next());
        if (args.remaining() != 2) throw usage_error("expected <index_file> and <patterns_file>");
        check_docs_query_options(cli.docs);
        cli.path_index_file = args.next();
        cli.path_patterns_file = args.next();
        set_default_query_threads(cli.options);
        direct_ifstream index_file(cli.path_index_file);
        require_readable(index_file, cli.path_index_file);
        const move_r_support support = read_move_rb_support(index_file);

        if (support == _count || support == _locate_one)
            throw std::runtime_error(cli.path_index_file + " has no documents (it was built with -s count or locate-one)");

        open_output_file(cli.output_file, cli.options.path_output_file);
        const bool known = dispatch_move_rb_index<_require_all_occs>(support,
            [&]<move_r_support index_support>() { run_queries<index_support>(index_file); });

        if (!known) throw std::runtime_error(cli.path_index_file + " is not an index built with move-rb build");
        close_output_file(cli.output_file, cli.options.path_output_file);
    });
}
