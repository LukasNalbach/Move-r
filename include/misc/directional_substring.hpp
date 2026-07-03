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
#include <string>

#include <misc/utils.hpp>

/**
 * @brief a view of the substring string[start..end], that can be read either left-to-right or right-to-left
 * @tparam pos_t unsigned integer position type
 * @tparam inp_t the underlying sequence type (e.g. std::string or std::vector<uint8_t>)
 */
template <typename pos_t, typename inp_t = std::string>
class directional_substring {
  public:
    using sym_t = typename inp_t::value_type;

  protected:
    const inp_t& string;
    pos_t start = 1;
    pos_t end = 0;
    direction_t dir = NO_DIR;

  public:
    directional_substring() = default;

    /**
     * @brief constructs a view of string[start..end] that is read in direction dir
     * @param string the underlying sequence
     * @param start start index of the substring (inclusive)
     * @param end end index of the substring (inclusive)
     * @param dir direction in which the substring is read (default: RIGHT)
     */
    directional_substring(const inp_t& string, pos_t start, pos_t end, direction_t dir = RIGHT)
        : string(string), start(start), end(end), dir(dir) {}

    /**
     * @brief returns the i-th symbol of the substring in the current reading direction: counted from the start
     *        when read left-to-right (RIGHT), from the end when read right-to-left (LEFT)
     * @param i index into the substring
     * @return the i-th symbol of the substring
     */
    inline const sym_t& operator[](pos_t i) const
    {
        return dir == LEFT ? string[end - i] : string[start + i];
    }

    /**
     * @brief returns the index of the i-th symbol (in reading direction) within the underlying sequence
     * @param i index into the substring
     * @return the position of the i-th symbol in the underlying sequence
     */
    inline pos_t pos(pos_t i) const
    {
        return dir == LEFT ? end - i : start + i;
    }

    /** @brief the underlying sequence this is a substring of */
    inline const inp_t& underlying() const
    {
        return string;
    }

    /**
     * @brief returns the length of the substring
     * @return the length of the substring
     */
    inline pos_t size() const
    {
        if (empty()) return 0;
        return end - start + 1;
    }

    /**
     * @brief returns whether the substring is empty
     * @return whether the substring is empty
     */
    inline bool empty() const
    {
        return end < start;
    }

    /**
     * @brief sets the direction in which the substring is read
     * @param dir the new reading direction
     */
    void set_direction(direction_t dir)
    {
        this->dir = dir;
    }

    /**
     * @brief returns the direction in which the substring is read
     * @return the current reading direction
     */
    inline direction_t direction() const
    {
        return dir;
    }
};