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

#ifndef MOVE_R_MISC_SERIALIZED_SIZE_HPP
#define MOVE_R_MISC_SERIALIZED_SIZE_HPP

#include <cstdint>
#include <ostream>
#include <streambuf>

// a stream buffer that counts the bytes written to it and discards them
class byte_counting_streambuf : public std::streambuf {
public:
    uint64_t bytes = 0; // number of bytes written so far

protected:
    std::streamsize xsputn(const char*, std::streamsize count) override
    {
        bytes += uint64_t(count);
        return count;
    }

    int_type overflow(int_type c) override
    {
        if (!traits_type::eq_int_type(c, traits_type::eof())) bytes++;
        return traits_type::not_eof(c);
    }
};

/**
 * @brief returns the number of bytes x.serialize writes; the bytes are not stored
 * @param x an object with a method serialize(std::ostream&) const
 */
template <typename T>
uint64_t serialized_size_in_bytes(const T& x)
{
    byte_counting_streambuf buffer;
    std::ostream out(&buffer);
    x.serialize(out);
    return buffer.bytes;
}

#endif // MOVE_R_MISC_SERIALIZED_SIZE_HPP
