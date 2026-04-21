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

#pragma once

#include <cstdint>

#include "utils.hpp"

template <typename pos_t, typename fnc_t>
inline static pos_t bin_search_max_leq(pos_t value, pos_t left, pos_t right, fnc_t value_at)
{
    pos_t middle;

    while (left != right) {
        middle = left + (right - left) / 2 + 1;

        if (value_at(middle) <= value) {
            left = middle;
        } else {
            right = middle - 1;
        }
    }

    return left;
}

template <typename pos_t, typename fnc_t>
inline static pos_t bin_search_min_geq(pos_t value, pos_t left, pos_t right, fnc_t value_at)
{
    pos_t middle;

    while (left != right) {
        middle = left + (right - left) / 2;

        if (value <= value_at(middle)) {
            right = middle;
        } else {
            left = middle + 1;
        }
    }

    return left;
}

template <typename pos_t, typename fnc_t>
inline static pos_t bin_search_max_lt(pos_t value, pos_t left, pos_t right, fnc_t value_at)
{
    pos_t middle;

    while (left != right) {
        middle = left + (right - left) / 2 + 1;

        if (value_at(middle) < value) {
            left = middle;
        } else {
            right = middle - 1;
        }
    }

    return left;
}

template <typename pos_t, typename fnc_t>
inline static pos_t bin_search_min_gt(pos_t value, pos_t left, pos_t right, fnc_t value_at)
{
    pos_t middle;

    while (left != right) {
        middle = left + (right - left) / 2;

        if (value < value_at(middle)) {
            right = middle;
        } else {
            left = middle + 1;
        }
    }

    return left;
}

template <typename pos_t, direction_t search_dir, typename fnc_t>
inline static pos_t exp_search_max_leq(pos_t value, pos_t left, pos_t right, fnc_t value_at)
{
    if (left == right) [[unlikely]] return left;
    pos_t cur_step_size = 1;

    if constexpr (search_dir == LEFT) {
        right -= cur_step_size;

        while (value < value_at(right)) {
            cur_step_size *= 2;

            if (right < left + cur_step_size) {
                cur_step_size = right;
                right = 0;
                break;
            }

            right -= cur_step_size;
        }

        return bin_search_max_leq<pos_t>(value, right, right + cur_step_size - 1, value_at);
    } else {
        left += cur_step_size;

        while (value_at(left) < value) {
            cur_step_size *= 2;

            if (right < cur_step_size || right - cur_step_size < left) {
                cur_step_size = right - left;
                left = right;
                break;
            }

            left += cur_step_size;
        }

        return bin_search_max_leq<pos_t>(value, left - cur_step_size, left, value_at);
    }
}