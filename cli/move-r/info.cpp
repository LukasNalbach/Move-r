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
 * @file info.cpp
 * @brief the mode info of move-r
 */

#include <iostream>
#include <string>

#include <misc/cli_help.hpp>
#include <misc/files.hpp>
#include <misc/index_cli.hpp>
#include <misc/move_r_cli.hpp>
#include <move_r/move_r.hpp>

#include "modes.hpp"

namespace {

/**
 * @brief prints the help of move-r info
 * @param out the stream to print to
 */
void print_help(std::ostream& out)
{
    help_printer h(out);
    h.intro("move-r info: prints the support type, the parameters, the sizes of the data structures and the modes of "
        "an index.", "move-r info <index_file>");
    h.group("arguments");
    h.option("<index_file>", "", "an index built with move-r build");
}

/**
 * @brief prints the information about an index
 * @tparam support the support type of the index
 * @param index the index
 * @param index_bytes the size of the serialized index in bytes
 */
template <move_r_support support>
void print_info(const move_r<support>& index, uint64_t index_bytes)
{
    using index_t = move_r<support>;
    std::string modes = "count";
    if constexpr (index_t::supports_locate) modes += index_t::supports_multiple_locate ? ", locate" : ", locate (one occurrence)";
    modes += ", revert";
    if constexpr (index_t::supports_multiple_locate) modes += ", extract-sa";
    modes += ", extract-bwt, extract-rlbwt";

    index_facts facts;
    facts.add("support: " + move_r_support_option_name(support));
    index.add_parameters(facts);
    facts.add(alphabet_description(index.byte_alphabet, index.alphabet_in_input_order()));
    print_index_info(facts, "index", index_bytes, index.input_size(), [&] { index.log_data_structure_sizes(); },
        modes);
}

} // namespace

int move_r_info(int argc, char** argv)
{
    return run_cli(argc, argv, print_help, [&] {
        if (argc != 2) throw usage_error("expected exactly one <index_file>");
        const std::string path_index_file = argv[1];
        direct_ifstream index_file(path_index_file);
        require_readable(index_file, path_index_file);
        const move_r_support support = move_r<>::peek_support(index_file);

        const bool known = dispatch_move_r_index(support, [&]<move_r_support supp>() {
            move_r<supp> index;
            print_info(index, load_index_with_progress(index, index_file));
        });

        if (!known) throw std::runtime_error(path_index_file + " is not an index built with move-r build");
    });
}
