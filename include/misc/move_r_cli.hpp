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
 * @file move_r_cli.hpp
 * @brief the support types of a move-r for the command-line tool move-r
 */

#ifndef MOVE_R_MISC_MOVE_R_CLI_HPP
#define MOVE_R_MISC_MOVE_R_CLI_HPP

#include <algorithm>
#include <move_r/move_r.hpp>
#include <string>

/**
 * @brief calls a callback with the support type of a serialized move-r as template parameter
 * @tparam callback_t type of callback
 * @param support the support type read from the serialized index (see move_r::peek_support)
 * @param callback a generic callable with the template parameter <support>
 * @return false <=> the support type is unknown
 */
template <typename callback_t>
bool dispatch_move_r_index(move_r_support support, callback_t&& callback)
{
    switch (support) {
        case _count:          callback.template operator()<_count>(); return true;
        case _locate_one:     callback.template operator()<_locate_one>(); return true;
        case _locate_move:    callback.template operator()<_locate_move>(); return true;
        case _locate_rlzsa:   callback.template operator()<_locate_rlzsa>(); return true;
        default:              return false;
    }
}

/**
 * @brief returns the name of a support type as a value of the option -s of move-r build: the name of the support type
 *        with '-' instead of '_' (see move_r_support_name), e.g. locate-move
 * @param support the support type
 */
inline std::string move_r_support_option_name(move_r_support support)
{
    std::string name(move_r_support_name(support));
    std::replace(name.begin(), name.end(), '_', '-');
    return name;
}

/**
 * @brief returns the support type with a name of the option -s of move-r build
 * @param name the name (count, locate-one, locate-move or locate-rlzsa)
 * @param support is set to the support type
 * @return false <=> there is no support type with the name
 */
inline bool parse_move_r_support(const std::string& name, move_r_support& support)
{
    for (move_r_support candidate : { _count, _locate_one, _locate_move, _locate_rlzsa }) {
        if (name == move_r_support_option_name(candidate)) {
            support = candidate;
            return true;
        }
    }

    return false;
}

#endif // MOVE_R_MISC_MOVE_R_CLI_HPP
