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

/**
 * @file strands.hpp
 * @brief processes independent pieces of work with several threads, each thread advancing several pieces alternately,
 *        one step per piece and turn, so that the memory accesses of the pieces of a thread overlap
 */

#ifndef MOVE_R_MISC_STRANDS_HPP
#define MOVE_R_MISC_STRANDS_HPP

#include <atomic>
#include <cstdint>
#include <utility>
#include <vector>

#include <omp.h>

#include <misc/memory.hpp>

/**
 * @brief processes num_pieces pieces with up to num_threads threads; each thread advances num_strands pieces (its
 *        strands) alternately, one step per piece and turn; a finished piece is replaced by the next piece not taken
 *        by any thread; with one strand, each thread processes its pieces one after another
 * @tparam strand_t type of a strand; it has a member w, the strand slot in [0..num_threads * num_strands - 1], which
 *         is unique among all threads and kept when the strand takes the next piece
 * @tparam start_fnc_t type of start
 * @tparam step_fnc_t type of step
 * @param num_pieces the number of pieces
 * @param num_threads the number of threads
 * @param num_strands the number of pieces each thread advances alternately
 * @param start start(k, s) sets the strand s to the beginning of the k-th piece, k in [0..num_pieces-1]
 * @param step step(s) processes the next step of the piece of s; returns false if the piece is finished
 */
template <typename strand_t, typename start_fnc_t, typename step_fnc_t>
void run_strands(uint64_t num_pieces, uint64_t num_threads, uint64_t num_strands, start_fnc_t start,
    step_fnc_t step)
{
    if (num_strands == 1) {
        #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 1)
        for (uint64_t k = 0; k < num_pieces; k++) {
            strand_t s;
            start(k, s);
            s.w = omp_get_thread_num();
            while (step(s)) { }
        }

        return;
    }

    std::atomic<uint64_t> next_piece = 0; // the next piece not taken by any thread

    #pragma omp parallel num_threads(num_threads)
    {
        // [0..num_strands-1] the strands of this thread; [0..active-1] are not finished
        std::vector<strand_t, cache_aligned_no_init_allocator<strand_t>> strands(num_strands);
        uint64_t active = 0;

        // sets s to the next piece not taken by any thread; returns false if there is none
        auto take_piece = [&](strand_t& s) {
            const uint64_t k = next_piece.fetch_add(1, std::memory_order_relaxed);
            if (k >= num_pieces) return false;
            start(k, s);
            return true;
        };

        while (active < num_strands && take_piece(strands[active])) {
            strands[active].w = omp_get_thread_num() * num_strands + active;
            active++;
        }

        while (active > 0) {
            for (uint64_t t = 0; t < active;) {
                strand_t& s = strands[t];

                if (!step(s) && !take_piece(s)) {
                    std::swap(s, strands[--active]);
                    continue;
                }

                t++;
            }
        }
    }
}

#endif // MOVE_R_MISC_STRANDS_HPP
