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

#include <algorithm>
#include <cstdint>
#include <string>
#include <iostream>
#include <fstream>
#include <climits>
#include <vector>
#include <omp.h>
#include <sdsl/int_vector_buffer.hpp>
#include <misc/utils.hpp>

static constexpr uint64_t byte_file_block_size = 128 * 1024;

/**
 * @brief applies f to every byte T[0..len-1] of the plain byte file `file`, in place, using up to p threads
 * @tparam fnc_t type of the byte-mapping function (uint8_t -> uint8_t)
 * @param file path of the plain byte file to transform in place
 * @param len number of bytes of `file` to transform (T[0..len-1])
 * @param p number of threads to use
 * @param f byte-mapping function; each byte c is replaced with f(c)
 */
template <typename fnc_t>
inline void transform_byte_file_parallel(const std::string& file, uint64_t len, uint16_t p, fnc_t f)
{
    static constexpr uint64_t block = byte_file_block_size;
    std::vector<sdsl::int_vector_buffer<8>> buf;
    buf.reserve(p);

    for (uint16_t i = 0; i < p; i++) {
        buf.emplace_back(sdsl::int_vector_buffer<8>(file, std::ios::in, block, 8, true));
    }

    uint64_t chunk = div_ceil<uint64_t>(div_ceil<uint64_t>(len, uint64_t(p)), block) * block;

    #pragma omp parallel num_threads(p)
    {
        uint16_t i_p = omp_get_thread_num();
        uint64_t b = std::min<uint64_t>(len, uint64_t(i_p) * chunk);
        uint64_t e = std::min<uint64_t>(len, b + chunk);

        for (uint64_t i = b; i < e; i++) {
            buf[i_p][i] = f(buf[i_p][i]);
        }
    }
}

/**
 * @brief prints an error message about a malformed patterns-file header and exits
 */
static void print_header_error()
{
    std::cout << "Error: malformed header in patterns file" << std::endl;
    std::cout << "Take a look here for more info on the file format: http://pizzachili.dcc.uchile.cl/experiments.html" << std::endl;
    exit(0);
}

/**
 * @brief parses the integer value following a key in a header line (e.g. "number=" in "number=100 length=10")
 * @param str the header line
 * @param key the key to look for
 * @return the value following the key, or -1 if the key is not found or has no value
 */
static int64_t value_from_key(std::string str, std::string key)
{
    uint64_t start_pos = str.find(key);

    if (start_pos == std::string::npos || start_pos + key.size() >= str.size()) {
        return -1;
    }

    start_pos += key.size();
    uint64_t end_pos = std::min(str.size(), str.substr(start_pos).find(" "));

    if (end_pos == std::string::npos) {
        return -1;
    }

    return std::atoi(str.substr(start_pos).substr(0, end_pos).c_str());
}

/**
 * @brief reads the number of patterns from a patterns-file header
 * @param header the header line
 * @return the number of patterns (exits with an error if the header is malformed)
 */
inline uint64_t number_of_patterns(std::string header)
{
    int64_t num = value_from_key(header, "number=");

    if (num < 1) {
        print_header_error();
    }

    return num;
}

/**
 * @brief reads the pattern length from a patterns-file header
 * @param header the header line
 * @return the pattern length (exits with an error if the header is malformed)
 */
inline uint64_t patterns_length(std::string header)
{
    int64_t len = value_from_key(header, "length=");

    if (len < 1) {
        print_header_error();
    }

    return len;
}

/**
 * @brief reads size bytes from the input stream into data, in chunks that fit into a 32-bit signed integer
 * @param in the input stream
 * @param data the destination buffer (must hold at least size bytes)
 * @param size the number of bytes to read
 */
inline void read_from_file(std::istream& in, const char* data, uint64_t size)
{
    uint64_t size_left = size;
    uint64_t bytes_to_read;

    while (size_left > 0) {
        bytes_to_read = std::min<uint64_t>(size_left, INT_MAX);
        in.read((char*) &data[size - size_left], bytes_to_read);
        size_left -= bytes_to_read;
    }
}

/**
 * @brief writes size bytes from data to the output stream, in chunks that fit into a 32-bit signed integer
 * @param out the output stream
 * @param data the source buffer (must hold at least size bytes)
 * @param size the number of bytes to write
 */
inline void write_to_file(std::ostream& out, const char* data, uint64_t size)
{
    uint64_t size_left = size;
    uint64_t bytes_to_write;

    while (size_left > 0) {
        bytes_to_write = std::min<uint64_t>(size_left, INT_MAX);
        out.write((char*) &data[size - size_left], bytes_to_write);
        size_left -= bytes_to_write;
    }
}