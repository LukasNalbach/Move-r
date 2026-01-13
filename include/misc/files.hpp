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
#include <iostream>
#include <fstream>
#include <limits>

static void print_header_error()
{
    std::cout << "Error: malformed header in patterns file" << std::endl;
    std::cout << "Take a look here for more info on the file format: http://pizzachili.dcc.uchile.cl/experiments.html" << std::endl;
    exit(0);
}

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

static uint64_t number_of_patterns(std::string header)
{
    int64_t num = value_from_key(header, "number=");

    if (num == -1) {
        print_header_error();
    }

    return num;
}

static uint64_t patterns_length(std::string header)
{
    int64_t len = value_from_key(header, "length=");

    if (len == -1) {
        print_header_error();
    }

    return len;
}

static void read_from_file(std::istream& in, const char* data, uint64_t size)
{
    uint64_t size_left = size;
    uint64_t bytes_to_read;

    while (size_left > 0) {
        bytes_to_read = std::min<uint64_t>(size_left, std::numeric_limits<int32_t>::max());
        in.read((char*) &data[size - size_left], bytes_to_read);
        size_left -= bytes_to_read;
    }
}

static void write_to_file(std::ostream& out, const char* data, uint64_t size)
{
    uint64_t size_left = size;
    uint64_t bytes_to_write;

    while (size_left > 0) {
        bytes_to_write = std::min<uint64_t>(size_left, std::numeric_limits<int32_t>::max());
        out.write((char*) &data[size - size_left], bytes_to_write);
        size_left -= bytes_to_write;
    }
}