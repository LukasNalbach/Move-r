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
#include <iostream>
#include <string>
#include <utility>

/**
 * @brief a percentage indicator that redraws itself in place (via a carriage
 *        return) instead of printing a new line for every update, so the
 *        terminal shows a single, live progress line while a loop runs
 */
class progress_meter {
    std::string prefix;
    uint64_t total;
    uint64_t done = 0;
    int64_t last_perc = -1;

    /**
     * @brief redraws the progress line for the given percentage
     * @param perc the percentage to display
     */
    void draw(int64_t perc)
    {
        std::cout << '\r' << prefix << perc << "% done" << std::flush;
    }

public:
    /**
     * @brief constructs the meter and draws the initial 0% state
     * @param total total number of steps (a value of 0 disables all output)
     * @param prefix text printed in front of the percentage
     */
    progress_meter(uint64_t total, std::string prefix = "")
        : prefix(std::move(prefix)), total(total)
    {
        if (total != 0) draw(0);
    }

    /**
     * @brief advances the meter by one step, redrawing only when the integer
     *        percentage actually changes
     */
    void step()
    {
        if (total == 0) return;
        done++;
        int64_t perc = (100 * done) / total;

        if (perc != last_perc) {
            last_perc = perc;
            draw(perc);
        }
    }

    /**
     * @brief completes the meter by drawing 100% and ending the line
     */
    void finish()
    {
        if (total == 0) return;

        if (last_perc != 100) {
            last_perc = 100;
            std::cout << '\r' << prefix << "100% done";
        }

        std::cout << std::endl;
    }
};
