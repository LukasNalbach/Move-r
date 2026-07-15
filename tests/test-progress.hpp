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

#include <atomic>
#include <chrono>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <string>
#include <vector>
#include <unistd.h>

#include <omp.h>

/**
 * @brief the number of random inputs to validate in a test, overridable via the MOVE_R_TEST_ITERATIONS
 *        environment variable (lower it for a quick smoke test, raise it for a soak run)
 * @param default_iterations the default iteration count if the environment variable is unset/invalid
 */
inline uint64_t fuzz_iterations(uint64_t default_iterations = 1000)
{
    if (const char* env = std::getenv("MOVE_R_TEST_ITERATIONS")) {
        uint64_t value = std::strtoull(env, nullptr, 10);
        if (value > 0) return value;
    }

    return default_iterations;
}

/**
 * @brief one named functionality (e.g. count, locate, random-access, serialize) to fuzz; run(iteration) is
 *        expected to build a fresh random structure and verify exactly this one functionality
 *
 * parallel controls how the threads are used: structures whose construction and queries already parallelise
 * internally (e.g. the indexes) leave it false, so their iterations run one at a time and each iteration uses
 * all threads; structures that are built and queried sequentially per instance (e.g. the bit vectors) set it
 * true, so run_fuzz spreads the iterations themselves across the threads instead. A parallel run's callback
 * must be thread-safe (in particular it must use a thread-local random number generator).
 */
struct fuzz_functionality {
    std::string name;
    std::function<void(uint64_t)> run;
    bool parallel = false;
};

/**
 * @brief fuzzes a structure/index by running each of its functionalities for a fixed number of iterations
 *        (each iteration builds a random data structure and verifies exactly that one functionality); a
 *        googletest-style (green-bracketed) progress bar labelled "structure/functionality" tracks progress
 *        and reports, per functionality, how many inputs were validated and how long it took
 * @param structure the name of the structure/index being fuzzed (shown in the progress bar)
 * @param functionalities the functionalities to fuzz (the total iteration count is split evenly over them)
 * @param iterations the total number of iterations to run for the whole test (defaults to fuzz_iterations())
 */
inline void run_fuzz(
    const std::string& structure, const std::vector<fuzz_functionality>& functionalities,
    uint64_t iterations = fuzz_iterations())
{
    using clock = std::chrono::steady_clock;
    const bool tty = isatty(fileno(stderr)) != 0;
    const char* green = tty ? "\033[0;32m" : "";
    const char* reset = tty ? "\033[m" : "";

    const uint64_t total_iterations = std::max<uint64_t>(1, iterations);
    const uint64_t num_functionalities = std::max<uint64_t>(1, functionalities.size());

    auto seconds_since = [](clock::time_point t) {
        return std::chrono::duration<double>(clock::now() - t).count();
    };

    int32_t last_pct = -1;
    auto draw = [&](const std::string& label, uint64_t done, double elapsed) {
        constexpr int32_t width = 10; // 10 columns inside the brackets, like googletest
        const double fraction = (double) done / total_iterations;
        int32_t filled = (int32_t) (fraction * width + 0.5);
        if (filled > width) filled = width;
        std::string bar(filled, '=');
        bar.resize(width, ' ');
        // show the total iterations completed so far out of all iterations (rather than a percentage)
        // on a terminal redraw the same line (\r); otherwise emit a new line per milestone
        std::fprintf(stderr, "%s%s[%s]%s %-28s %7" PRIu64 "/%-7" PRIu64 " %5.1fs   ",
            tty ? "\r" : "\n", green, bar.c_str(), reset, (structure + "/" + label).c_str(),
            (uint64_t) done, (uint64_t) total_iterations, elapsed);
        std::fflush(stderr);
    };

    // split the total iteration count evenly over the functionalities (the first few absorb the remainder)
    uint64_t completed = 0; // iterations of already-finished functionalities
    for (uint64_t f = 0; f < functionalities.size(); f++) {
        const fuzz_functionality& functionality = functionalities[f];
        const uint64_t slice_iterations = total_iterations / num_functionalities
            + (f < total_iterations % num_functionalities ? 1 : 0);
        const auto slice_start = clock::now();

        if (functionality.parallel) {
            // the per-instance work is sequential, so spread the iterations over the threads; the callback
            // must be thread-safe (thread-local RNG). progress is tracked via an atomic and drawn by thread 0
            std::atomic<uint64_t> finished { 0 };

            #pragma omp parallel for schedule(dynamic)
            for (int64_t iteration = 0; iteration < (int64_t) slice_iterations; iteration++) {
                functionality.run((uint64_t) iteration);
                const uint64_t done = completed + finished.fetch_add(1) + 1;
                if (omp_get_thread_num() == 0) {
                    const int32_t pct = (int32_t) (100.0 * done / total_iterations + 0.5);
                    if (pct != last_pct && (tty || pct % 10 == 0)) {
                        last_pct = pct;
                        draw(functionality.name, done, seconds_since(slice_start));
                    }
                }
            }
        } else {
            for (uint64_t iteration = 0; iteration < slice_iterations; iteration++) {
                functionality.run(iteration);
                const uint64_t done = completed + iteration + 1;
                const int32_t pct = (int32_t) (100.0 * done / total_iterations + 0.5);
                if (pct != last_pct && (tty || pct % 10 == 0)) {
                    last_pct = pct;
                    draw(functionality.name, done, seconds_since(slice_start));
                }
            }
        }

        // leave the completed functionality (with the running iteration total and its wall time) on its own line
        completed += slice_iterations;
        draw(functionality.name, completed, seconds_since(slice_start));
        std::fputc('\n', stderr);
        last_pct = -1;
    }
}
