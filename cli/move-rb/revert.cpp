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
 * @file revert.cpp
 * @brief the mode revert of move-rb
 */

#include <stdexcept>
#include <string>

#include <misc/move_rb_cli.hpp>
#include <misc/revert_cli.hpp>
#include <move_rb/move_rb.hpp>

#include "modes.hpp"

int move_rb_revert(int argc, char** argv)
{
    const revert_mode_t mode { .tool = "move-rb", .description = move_rb_revert_description,
        .index_explanation = "an index built with move-rb build",
        .reversed_explanation = "write the reversed text" };

    return run_revert_mode(argc, argv, mode,
        [](direct_ifstream& index_file, direct_ofstream& output_file, revert_cli_t& cli) {
            const bool known = dispatch_move_rb_index(read_move_rb_support(index_file), [&]<move_r_support support>() {
                revert_index<move_rb<support>>(index_file, output_file, cli, "bidirectional index",
                    "revert_move_rb_" + std::string(move_r_support_name(support)));
            });

            if (!known) throw std::runtime_error(cli.path_index_file + " is not an index built with move-rb build");
        });
}
