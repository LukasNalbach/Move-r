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
 * @file temporary_input.hpp
 * @brief removes an input file written by the construction and its caches
 */

#ifndef MOVE_R_MOVE_RB_CONSTRUCTION_TEMPORARY_INPUT_HPP
#define MOVE_R_MOVE_RB_CONSTRUCTION_TEMPORARY_INPUT_HPP

#include <filesystem>
#include <string>
#include <system_error>

#include <pfp_rlbwt/pfp_rlbwt.hpp>

/**
 * @brief removes an input file written by the construction, together with the caches of the RLBWT and its samples
 *        and the files of the prefix-free parsing that are named after it
 * @param input_file path of the input file
 */
inline void remove_temporary_input(const std::string& input_file)
{
    std::error_code ec;
    if (input_file.empty()) return;

    // the input file and its caches (see move_r::pfp_cache_path)
    const std::string cache = input_file + pfp_rlbwt::ext_runs;
    for (const std::string& file : { input_file, cache, cache + "-rev" }) std::filesystem::remove(file, ec);
}

#endif // MOVE_R_MOVE_RB_CONSTRUCTION_TEMPORARY_INPUT_HPP
