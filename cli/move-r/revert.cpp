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
 * @brief the mode revert of move-r
 */

#include <stdexcept>
#include <string>

#include <misc/move_r_cli.hpp>
#include <misc/revert_cli.hpp>
#include <move_r/move_r.hpp>

#include "modes.hpp"

int move_r_revert(int argc, char** argv)
{
    const revert_mode_t mode { .tool = "move-r", .description = move_r_revert_description,
        .index_explanation = "an index built with move-r build" };

    return run_revert_mode(argc, argv, mode,
        [](direct_ifstream& index_file, direct_ofstream& output_file, revert_cli_t& cli) {
            const move_r_support support = move_r<>::peek_support(index_file);

            const bool known = dispatch_move_r_index(support, [&]<move_r_support supp>() {
                revert_index<move_r<supp>>(index_file, output_file, cli, "index",
                    "revert_move_r_" + std::string(move_r_support_name(supp)));
            });

            if (!known) throw std::runtime_error(cli.path_index_file + " is not an index built with move-r build");
        });
}
