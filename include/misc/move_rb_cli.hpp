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
 * @file move_rb_cli.hpp
 * @brief dispatches a serialized move-rb to the template parameters of its index type (for the modes of move-rb)
 */

#ifndef MOVE_R_MISC_MOVE_RB_CLI_HPP
#define MOVE_R_MISC_MOVE_RB_CLI_HPP

#include <cstdint>
#include <move_rb/move_rb.hpp>

/** @brief the locate support a command-line tool requires of an index (see move_r_support) */
enum move_rb_locate_requirement_t : uint8_t {
    _require_count = 0, // none (the tool only counts)
    _require_one_occ = 1, // one occurrence of a pattern (see locate_context_t::sampled_occ)
    _require_all_occs = 2 // all occurrences of a pattern (see move_rb::supports_multiple_locate)
};

/**
 * @brief calls a callback with the type of locate support of a serialized move-rb as template parameter, if the
 *        index supports the required locate queries
 * @tparam required_locate the required locate support
 * @tparam callback_t type of callback
 * @param support the type of locate support read from the serialized index (see read_move_rb_support)
 * @param callback a generic callable with the template parameter <support>
 * @return whether support is the locate support of a move-rb with the required locate queries
 */
template <move_rb_locate_requirement_t required_locate = _require_count, typename callback_t>
bool dispatch_move_rb_index(move_r_support support, callback_t&& callback)
{
    if (support == _count) {
        if constexpr (required_locate == _require_count) {
            callback.template operator()<_count>();
            return true;
        }
    } else if (support == _locate_one) {
        if constexpr (required_locate != _require_all_occs) {
            callback.template operator()<_locate_one>();
            return true;
        }
    } else if (support == _locate_move) {
        callback.template operator()<_locate_move>();
        return true;
    } else if (support == _locate_rlzsa) {
        callback.template operator()<_locate_rlzsa>();
        return true;
    }

    return false;
}

#endif // MOVE_R_MISC_MOVE_RB_CLI_HPP
