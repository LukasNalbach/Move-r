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
 * @file directional_substring.hpp
 * @brief a view of a substring that is read left-to-right or right-to-left
 */

#ifndef MOVE_R_MISC_DIRECTIONAL_SUBSTRING_HPP
#define MOVE_R_MISC_DIRECTIONAL_SUBSTRING_HPP

#include <cstdint>
#include <string>

#include <misc/utils.hpp>

/**
 * @brief a view of the substring string[beg..end] that can be read either left-to-right or right-to-left; the
 *        substring is empty if end < beg (beg and end are signed: end can be -1)
 * @tparam inp_t the underlying sequence type (e.g. std::string or std::vector<uint8_t>)
 */
template <typename inp_t = std::string>
class directional_substring {
  public:
    using sym_t = typename inp_t::value_type;

  protected:
    const inp_t& string; // the underlying sequence
    int64_t beg = 1; // position of the first symbol of the substring in string
    int64_t end = 0; // position of the last symbol of the substring in string
    direction_t dir = NO_DIR; // direction in which the substring is read

  public:
    directional_substring() = default;

    /**
     * @brief constructs a view of string[beg..end] that is read in direction dir
     * @param string the underlying sequence
     * @param beg position of the first symbol of the substring
     * @param end position of the last symbol of the substring
     * @param dir direction in which the substring is read (default: RIGHT)
     */
    directional_substring(const inp_t& string, uint64_t beg, uint64_t end, direction_t dir = RIGHT)
        : string(string), beg(beg), end(end), dir(dir) {}

    /**
     * @brief returns the i-th symbol of the substring in the reading direction: counted from beg if it is read
     *        left-to-right (RIGHT), from end if it is read right-to-left (LEFT)
     * @param i [0..size()-1] index in the substring
     * @return the i-th symbol of the substring
     */
    inline const sym_t& operator[](uint64_t i) const
    {
        return dir == LEFT ? string[end - i] : string[beg + i];
    }

    /**
     * @brief returns the position of the i-th symbol (in reading direction) in the underlying sequence
     * @param i [0..size()-1] index in the substring
     * @return the position of the i-th symbol in the underlying sequence
     */
    inline uint64_t pos(uint64_t i) const { return dir == LEFT ? end - i : beg + i; }

    // returns the underlying sequence
    inline const inp_t& underlying() const { return string; }

    // returns the length of the substring
    inline uint64_t size() const
    {
        if (empty()) return 0;
        return end - beg + 1;
    }

    // returns whether the substring is empty
    inline bool empty() const { return end < beg; }

    // sets the direction in which the substring is read
    void set_direction(direction_t dir) { this->dir = dir; }

    // returns the direction in which the substring is read
    inline direction_t direction() const { return dir; }
};

#endif // MOVE_R_MISC_DIRECTIONAL_SUBSTRING_HPP
