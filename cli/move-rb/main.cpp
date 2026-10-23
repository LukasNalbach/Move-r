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
 * @file main.cpp
 * @brief the command-line tool move-rb: runs the mode named by the first argument
 */

#include <misc/cli_modes.hpp>

#include "modes.hpp"

int main(int argc, char** argv)
{
    return run_cli_mode("move-rb", move_rb_description, {
        { "build", move_rb_build_description, move_rb_build },
        { "count", move_rb_count_description, move_rb_count },
        { "locate", move_rb_locate_description, move_rb_locate },
        { "revert", move_rb_revert_description, move_rb_revert },
        { "extract-sa", move_rb_extract_sa_description, move_rb_extract_sa },
        { "extract-bwt", move_rb_extract_bwt_description, move_rb_extract_bwt },
        { "extract-rlbwt", move_rb_extract_rlbwt_description, move_rb_extract_rlbwt },
        { "info", move_rb_info_description, move_rb_info }
    }, argc, argv);
}
