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

// Shared query model for move-rb-bench-ext and move-rb-gen-ext-queries: a bidirectional exact search that starts
// at a fixed position inside the pattern and extends it left/right in a given order. As every pattern is a
// substring of the text, every extension succeeds. Both tools draw the extension plan from a plan RNG seeded with
// ext_plan_seed and reset per pattern set, so every index replays the same searches without a seed in the file.

#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include <misc/apm.hpp> // direction_t (LEFT/RIGHT), query_support_t (COUNT/LOCATE)

// the seed the plan RNG is reset to at the start of each pattern set
inline constexpr uint64_t ext_plan_seed = 0x9E3779B97F4A7C15ull;

// the order in which a pattern is assembled during a bidirectional-extension query
struct ext_plan_t {
    uint64_t start = 0; // position in [0, m) the search starts at
    // the m-1 extension steps after the first character (LEFT = prepend, RIGHT = append)
    std::vector<direction_t> order;
};

/**
 * @brief draws the next extension plan for a pattern of length @p m from @p rng, advancing it
 * @param m the pattern length
 * @param rng the plan RNG
 * @return the start position and the L/R extension order
 */
inline ext_plan_t make_ext_plan(uint64_t m, std::mt19937_64& rng)
{
    ext_plan_t plan;
    plan.start = m == 0 ? 0 : std::uniform_int_distribution<uint64_t>(0, m - 1)(rng);

    uint64_t left_rem = plan.start;          // characters still to be prepended (positions start-1 .. 0)
    uint64_t right_rem = m - 1 - plan.start; // characters still to be appended (positions start+1 .. m-1)
    plan.order.reserve(m == 0 ? 0 : m - 1);

    while (left_rem > 0 || right_rem > 0) {
        bool go_left;
        if (left_rem == 0) go_left = false;
        else if (right_rem == 0) go_left = true;
        // pick a side weighted by how many characters remain on it, so the interleaving is unbiased
        else go_left = std::uniform_int_distribution<uint64_t>(1, left_rem + right_rem)(rng) <= left_rem;

        if (go_left) { plan.order.push_back(LEFT); left_rem--; }
        else         { plan.order.push_back(RIGHT); right_rem--; }
    }

    return plan;
}

/**
 * @brief performs the bidirectional exact search of @p P for the given plan and returns the occurrence count
 * @tparam idx_t an index exposing the shared bidirectional query interface (move_rb or an index adapter)
 * @param index the index to query
 * @param P the pattern (a substring of the indexed text)
 * @param plan the start position and L/R extension order
 * @return the number of occurrences of P (the width of the final SA-interval)
 */
template <typename idx_t>
inline typename idx_t::pos_type ext_count(const idx_t& index, const std::string& P, const ext_plan_t& plan)
{
    using pos_t = typename idx_t::pos_type;
    const uint64_t m = P.size();
    if (m == 0) return 0;

    auto ctx = index.template empty_context<COUNT>();
    { auto [nxt, ok] = ctx.extend(P[plan.start], RIGHT); if (!ok) return 0; ctx = nxt; }

    int64_t lp = (int64_t) plan.start - 1; // next character to prepend
    int64_t rp = (int64_t) plan.start + 1; // next character to append
    for (direction_t dir : plan.order) {
        if (dir == LEFT) { auto [nxt, ok] = ctx.extend(P[lp--], LEFT);  if (!ok) return 0; ctx = nxt; }
        else             { auto [nxt, ok] = ctx.extend(P[rp++], RIGHT); if (!ok) return 0; ctx = nxt; }
    }

    return (pos_t) ctx.num_occ();
}

/**
 * @brief performs the bidirectional exact search of @p P for the given plan and then locates all occurrences,
 *        calling @p report for each one
 * @tparam idx_t an index exposing the shared bidirectional query interface (move_rb or an index adapter)
 * @tparam report_fnc_t a callable invoked with each occurrence position
 * @param index the index to query
 * @param P the pattern (a substring of the indexed text)
 * @param plan the start position and L/R extension order
 * @param report called with every occurrence position of P
 */
template <typename idx_t, typename report_fnc_t>
inline void ext_locate(const idx_t& index, const std::string& P, const ext_plan_t& plan, report_fnc_t report)
{
    const uint64_t m = P.size();
    if (m == 0) return;

    auto ctx = index.template empty_context<LOCATE>();
    { auto [nxt, ok] = ctx.extend(P[plan.start], RIGHT); if (!ok) return; ctx = nxt; }

    int64_t lp = (int64_t) plan.start - 1;
    int64_t rp = (int64_t) plan.start + 1;
    for (direction_t dir : plan.order) {
        if (dir == LEFT) { auto [nxt, ok] = ctx.extend(P[lp--], LEFT);  if (!ok) return; ctx = nxt; }
        else             { auto [nxt, ok] = ctx.extend(P[rp++], RIGHT); if (!ok) return; ctx = nxt; }
    }

    ctx.locate_phase().locate(report);
}
