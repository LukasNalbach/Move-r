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

#include "utils.hpp"

template <typename pos_t>
class sub_string {
  private:
    const std::string& string;
    pos_t start = 1;
    pos_t end = 0;
    direction_t dir = NO_DIR;

    using access_fnc_t = const char&(sub_string::*)(pos_t) const;
    access_fnc_t access;

    inline const char& right_access(pos_t i) const
    {
        return string[start + i];
    }

    inline const char& left_access(pos_t i) const
    {
        return string[end - i];
    }

  public:
    sub_string() = default;
    
    sub_string(const std::string& string, pos_t start, pos_t end, direction_t dir = RIGHT)
        : string(string), start(start), end(end), dir(dir)
    {
        set_direction(dir);
    }

    inline const char& operator[](pos_t i) const
    {
        return (this->*access)(i);
    }

    inline pos_t size() const
    {
        if (empty()) return 0;
        return end - start + 1;
    }

    inline bool empty() const
    {
        return end < start;
    }

    void set_direction(direction_t dir)
    {
        this->dir = dir;
        access = dir == LEFT ? &sub_string<pos_t>::left_access
                             : &sub_string<pos_t>::right_access;
    }

    inline direction_t direction() const
    {
        return dir;
    }
};