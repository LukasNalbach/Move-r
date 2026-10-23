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

#ifndef MOVE_R_MISC_TEXT_LCE_HPP
#define MOVE_R_MISC_TEXT_LCE_HPP

#include <bit>
#include <cstdint>
#include <cstring>

#include <misc/utils.hpp>

/**
 * @brief returns the length of the longest common extension of two positions of a text in memory, i.e., the
 *        maximum l <= max_len with T[i..i+l-1] = T[j..j+l-1] (dir = RIGHT) or T[i-l+1..i] = T[j-l+1..j]
 *        (dir = LEFT); the text is compared 8 bytes at a time
 * @tparam dir direction of the comparison (LEFT: towards the start of the text, RIGHT: towards its end)
 * @param text the text
 * @param sym_width number of bytes per symbol
 * @param i position in the text (in symbols)
 * @param j position in the text (in symbols)
 * @param max_len upper bound for the result
 * @return the lce (in symbols)
 */
template <direction_t dir>
inline uint64_t text_lce(const char* text, uint8_t sym_width, uint64_t i, uint64_t j, uint64_t max_len)
{
    constexpr bool left = dir == LEFT;
    // positions of the first (RIGHT) or the last (LEFT) byte of T[i] and T[j]
    const uint64_t a = (i + left) * sym_width - left;
    const uint64_t b = (j + left) * sym_width - left;
    const uint64_t max_bytes = max_len * sym_width;
    uint64_t num_bytes = 0; // number of equal bytes

    while (num_bytes + 8 <= max_bytes) {
        // the 8 bytes starting (RIGHT) or ending (LEFT) at the next bytes to compare
        uint64_t word_a, word_b;
        std::memcpy(&word_a, text + (left ? a - num_bytes - 7 : a + num_bytes), 8);
        std::memcpy(&word_b, text + (left ? b - num_bytes - 7 : b + num_bytes), 8);
        const uint64_t diff = word_a ^ word_b;
        if (diff != 0) return (num_bytes + (left ? std::countl_zero(diff) : std::countr_zero(diff)) / 8) / sym_width;
        num_bytes += 8;
    }

    while (num_bytes < max_bytes &&
           text[left ? a - num_bytes : a + num_bytes] == text[left ? b - num_bytes : b + num_bytes]) num_bytes++;

    return num_bytes / sym_width;
}

#endif // MOVE_R_MISC_TEXT_LCE_HPP
