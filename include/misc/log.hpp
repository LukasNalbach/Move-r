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
#include <iostream>
#include <iomanip>
#include <sstream>

/**
 * @brief formats a pair as the string "(first,second)"
 * @tparam pos_t type of the pair elements
 * @param pair a pair
 * @return the string "(first,second)"
 */
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

/**
 * @brief returns the current time point of the steady clock
 * @return the current time point
 */
static std::chrono::steady_clock::time_point now()
{
    return std::chrono::steady_clock::now();
}

/**
 * @brief formats a duration in nanoseconds with an appropriate unit (ns, us, ms or s)
 * @param ns a duration in nanoseconds
 * @return the formatted duration
 */
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

/**
 * @brief formats a query throughput (queries over a duration) with an appropriate unit
 * @param num_queries the number of queries
 * @param ns the duration in nanoseconds
 * @return the formatted query throughput
 */
inline std::string format_query_throughput(uint64_t num_queries, uint64_t ns)
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

/**
 * @brief computes an index-construction throughput in MB/s (megabytes = 10^6 bytes, matching format_size)
 * @param n_bytes the number of input bytes processed (the text length)
 * @param time_ns the construction duration in nanoseconds
 * @return the throughput in MB/s (0 if time_ns is 0)
 */
inline double construction_throughput_mb_per_s(uint64_t n_bytes, uint64_t time_ns)
{
    return time_ns == 0 ? 0.0 : (n_bytes * 1000.0) / time_ns; // (n_bytes / 1e6) / (time_ns / 1e9)
}

/**
 * @brief formats an index-construction throughput in MB/s (megabytes = 10^6 bytes, matching format_size)
 * @param n_bytes the number of input bytes processed (the text length)
 * @param time_ns the construction duration in nanoseconds
 * @return the formatted throughput, e.g. "123.45 MB/s"
 */
inline std::string format_construction_throughput(uint64_t n_bytes, uint64_t time_ns)
{
    std::ostringstream str;
    str << std::fixed << std::setprecision(2) << construction_throughput_mb_per_s(n_bytes, time_ns) << " MB/s";
    return str.str();
}

/**
 * @brief formats a size in bytes with an appropriate unit (B, KB, MB or GB)
 * @param B a size in bytes
 * @return the formatted size
 */
inline std::string format_size(uint64_t B)
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

/**
 * @brief formats a thread count as "1 thread" or "p threads"
 * @param p a number of threads
 * @return the formatted thread count
 */
inline std::string format_threads(uint16_t p)
{
    if (p == 1) {
        return "1 thread";
    } else {
        return std::to_string(p) + " threads";
    }
}

/**
 * @brief returns the time difference (in minutes) between two time points
 * @param t1 the earlier time point
 * @param t2 the later time point
 * @return the time difference t2 - t1 in minutes
 */
static uint64_t time_diff_min(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2)
{
    return std::chrono::duration_cast<std::chrono::minutes>(t2 - t1).count();
}

/**
 * @brief returns the time difference (in nanoseconds) between two time points
 * @param t1 the earlier time point
 * @param t2 the later time point
 * @return the time difference t2 - t1 in nanoseconds
 */
static uint64_t time_diff_ns(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2)
{
    return std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
}

/**
 * @brief returns the time difference (in minutes) between a time point and now
 * @param t the earlier time point
 * @return the time difference now - t in minutes
 */
inline uint64_t time_diff_min(std::chrono::steady_clock::time_point t)
{
    return time_diff_min(t, std::chrono::steady_clock::now());
}

/**
 * @brief returns the time difference (in nanoseconds) between a time point and now
 * @param t the earlier time point
 * @return the time difference now - t in nanoseconds
 */
inline uint64_t time_diff_ns(std::chrono::steady_clock::time_point t)
{
    return time_diff_ns(t, std::chrono::steady_clock::now());
}

/**
 * @brief logs the runtime between two time points and returns the current time point
 * @param t1 the earlier time point
 * @param t2 the later time point
 * @return the current time point
 */
static std::chrono::steady_clock::time_point log_runtime(std::chrono::steady_clock::time_point t1, std::chrono::steady_clock::time_point t2)
{
    std::cout << ", in ~ " << format_time(time_diff_ns(t1, t2)) << std::endl;
    return std::chrono::steady_clock::now();
}

/**
 * @brief logs the runtime between a time point and now and returns the current time point
 * @param t the earlier time point
 * @return the current time point
 */
static std::chrono::steady_clock::time_point log_runtime(std::chrono::steady_clock::time_point t)
{
    return log_runtime(t, std::chrono::steady_clock::now());
}

/**
 * @brief starts a logged phase: prints msg and (re)sets the phase timer to now()
 * @param log whether logging is enabled
 * @param time phase timer; set to now() (only when logging)
 * @param msg message describing the phase
 */
inline void log_phase_start(bool log, std::chrono::steady_clock::time_point& time, const std::string& msg)
{
    if (log) {
        time = now();
        std::cout << msg << std::flush;
    }
}

/**
 * @brief ends a logged phase: optionally writes the elapsed time to the measurement
 *        stream mf under the key mf_key, then logs the runtime and resets the timer
 * @param log whether logging is enabled
 * @param time phase timer; read for the elapsed time, then reset by log_runtime
 * @param mf measurement-file stream (nullptr => nothing written)
 * @param mf_key measurement-file key for the elapsed time (empty => nothing written)
 */
inline void log_phase_end(bool log, std::chrono::steady_clock::time_point& time, std::ostream* mf = nullptr, const std::string& mf_key = "")
{
    if (log) {
        if (!mf_key.empty() && mf != nullptr)
            *mf << " " << mf_key << "=" << time_diff_ns(time, now());
        time = log_runtime(time);
    }
}

/**
 * @brief prints a message to std::cout
 * @param message the message to print
 */
inline void log_message(std::string message)
{
    std::cout << message << std::flush;
}

/**
 * @brief prints message (without resetting any phase timer), if logging is enabled;
 *        use for phase-start messages that must not reset the timer and for section headers
 * @param log whether logging is enabled
 * @param message message to print (include a trailing "\n" for a header line)
 */
inline void log_message(bool log, const std::string& message)
{
    if (log) std::cout << message << std::flush;
}