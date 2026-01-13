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
#include <chrono>

template <typename pos_t = uint32_t>
inline static std::string to_string(std::pair<pos_t, pos_t> pair)
{
    std::string str = "(";
    str.append(std::to_string(pair.first));
    str.push_back(',');
    str.append(std::to_string(pair.second));
    str.push_back(')');
    return str;
}

static std::chrono::steady_clock::time_point now()
{
    return std::chrono::steady_clock::now();
}

static std::string format_time(uint64_t ns)
{
    std::string time_str;

    if (ns > 10000000000) {
        time_str = std::to_string(ns / 1000000000) + " s";
    } else if (ns > 10000000) {
        time_str = std::to_string(ns / 1000000) + " ms";
    } else if (ns > 10000) {
        time_str = std::to_string(ns / 1000) + " us";
    } else {
        time_str = std::to_string(ns) + " ns";
    }

    return time_str;
}

static std::string format_query_throughput(uint64_t num_queries, uint64_t ns)
{
    std::string str;
    double queries_per_ns = num_queries / (double)ns;

    if (queries_per_ns < 0.000001) {
        str = std::to_string(queries_per_ns * 1000000000) + " queries/s";
    } else if (queries_per_ns < 0.001) {
        str = std::to_string(queries_per_ns * 1000000) + " queries/ms";
    } else if (queries_per_ns < 1) {
        str = std::to_string(queries_per_ns * 1000) + " queries/us";
    } else {
        str = std::to_string(queries_per_ns) + " queries/ns";
    }

    return str;
}

static std::string format_size(uint64_t B)
{
    std::string size_str;

    if (B > 10000000000) {
        size_str = std::to_string(B / 1000000000) + " GB";
    } else if (B > 10000000) {
        size_str = std::to_string(B / 1000000) + " MB";
    } else if (B > 10000) {
        size_str = std::to_string(B / 1000) + " KB";
    } else {
        size_str = std::to_string(B) + " B";
    }

    return size_str;
}

static std::string format_threads(uint16_t p)
{
    if (p == 1) {
        return "1 thread";
    } else {
        return std::to_string(p) + " threads";
    }
}

static uint64_t time_diff_min(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2)
{
    return std::chrono::duration_cast<std::chrono::minutes>(t2 - t1).count();
}

static uint64_t time_diff_ns(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2)
{
    return std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
}

static uint64_t time_diff_min(std::chrono::steady_clock::time_point t)
{
    return time_diff_min(t, std::chrono::steady_clock::now());
}

static uint64_t time_diff_ns(std::chrono::steady_clock::time_point t)
{
    return time_diff_ns(t, std::chrono::steady_clock::now());
}

static std::chrono::steady_clock::time_point log_runtime(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2)
{
    std::cout << ", in ~ " << format_time(time_diff_ns(t1, t2)) << std::endl;
    return std::chrono::steady_clock::now();
}

static std::chrono::steady_clock::time_point log_runtime(std::chrono::steady_clock::time_point t)
{
    return log_runtime(t, std::chrono::steady_clock::now());
}

static void log_message(std::string message)
{
    std::cout << message << std::flush;
}