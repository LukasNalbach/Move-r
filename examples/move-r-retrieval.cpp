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

#include <move_r/move_r.hpp>

int main()
{
    // build an index
    move_r<> index("This is a test string");

    // retrieve the range [8,17] of the original text and store
    // it in a string using at most 2 threads
    std::string reverted_range = index.revert_range(
        { .l = 8, .r = 17, .num_threads = 2 });
    for (auto c : reverted_range) std::cout << c;
    std::cout << std::endl;

    // print the original text from right to left without storing it
    // using 1 thread
    index.revert_range([](auto, auto c) { std::cout << c; }, { .num_threads = 1 });
    std::cout << std::endl;

    // retrieve the suffix array values in the range [2,6] using at
    // most 4 threads and store them in a vector
    std::vector<uint32_t> SA_range = index.SA_range(
        { .l = 2, .r = 6, .num_threads = 4 });
    for (auto s : SA_range) std::cout << s << ", ";
    std::cout << std::endl;

    // print SA[1]
    std::cout << index.SA(1) << std::endl;

    // retrieve the BWT in the range [7,14] from left to right
    // using 1 thread
    index.BWT_range([](auto, auto s) { std::cout << s << ", "; },
        { .l = 7, .r = 14, .num_threads = 1 });
    std::cout << std::endl;

    // print BWT[16]
    std::cout << index.BWT(16) << std::endl;
}