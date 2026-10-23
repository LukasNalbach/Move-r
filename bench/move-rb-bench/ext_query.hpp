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
 * @file ext_query.hpp
 * @brief the queries of move-rb-bench-ext and move-rb-gen-ext-queries: an exact bidirectional search that starts at
 *        a random position of the pattern and adds the other symbols to the left and to the right in a random order
 *        (the extension plan); the plans of a pattern set are drawn from a generator with the seed ext_plan_seed,
 *        so both tools compute the same plans
 */

#ifndef MOVE_R_BENCH_MOVE_RB_BENCH_EXT_QUERY_HPP
#define MOVE_R_BENCH_MOVE_RB_BENCH_EXT_QUERY_HPP

#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include <misc/apm.hpp>
#include <misc/hashing.hpp>

// the seed of the generator the extension plans of a pattern set are drawn from (measurement setup: the value of
// golden_ratio_64)
inline constexpr uint64_t ext_plan_seed = golden_ratio_64;

// the order in which the symbols of a pattern P of length m are added
struct ext_plan_t {
    uint64_t start = 0; // [0..m-1] position of the first symbol
    // the directions of the other symbols: [0..m-2] direction of the (i+2)-th added symbol (LEFT: the next symbol
    // before the matched part of P is prepended, RIGHT: the next symbol after it is appended)
    std::vector<direction_t> order;
};

/**
 * @brief draws the extension plan of a pattern of length m: a uniformly random start position, then each of the
 *        remaining directions with a probability proportional to the number of symbols left on that side
 * @param m the pattern length
 * @param rng the generator
 */
inline ext_plan_t make_ext_plan(uint64_t m, std::mt19937_64& rng)
{
    ext_plan_t plan;
    plan.start = m == 0 ? 0 : std::uniform_int_distribution<uint64_t>(0, m - 1)(rng);

    uint64_t left_rem = plan.start; // number of symbols left to prepend
    uint64_t right_rem = m - 1 - plan.start; // number of symbols left to append
    plan.order.reserve(m == 0 ? 0 : m - 1);

    while (left_rem > 0 || right_rem > 0) {
        bool go_left;
        if (left_rem == 0) go_left = false;
        else if (right_rem == 0) go_left = true;
        else go_left = std::uniform_int_distribution<uint64_t>(1, left_rem + right_rem)(rng) <= left_rem;

        if (go_left) { plan.order.push_back(LEFT); left_rem--; }
        else         { plan.order.push_back(RIGHT); right_rem--; }
    }

    return plan;
}

/**
 * @brief matches P in the order of an extension plan
 * @param P the pattern (not empty)
 * @param plan the extension plan
 * @param ctx an empty search context; it is set to the search context of P
 * @return whether P occurs in the text
 */
template <typename ctx_t>
inline bool ext_search(const std::string& P, const ext_plan_t& plan, ctx_t& ctx)
{
    { auto [nxt, ok] = ctx.extend(P[plan.start], RIGHT); if (!ok) return false; ctx = nxt; }

    int64_t lp = (int64_t) plan.start - 1; // position of the next symbol to prepend
    int64_t rp = (int64_t) plan.start + 1; // position of the next symbol to append

    for (direction_t dir : plan.order) {
        if (dir == LEFT) { auto [nxt, ok] = ctx.extend(P[lp--], LEFT);  if (!ok) return false; ctx = nxt; }
        else             { auto [nxt, ok] = ctx.extend(P[rp++], RIGHT); if (!ok) return false; ctx = nxt; }
    }

    return true;
}

/**
 * @brief returns the number of occurrences of P, matched in the order of an extension plan
 * @param index move_rb or an index adapter (see indexes/index_adapter.hpp)
 * @param P the pattern
 * @param plan the extension plan
 */
template <typename idx_t>
inline uint64_t ext_count(const idx_t& index, const std::string& P, const ext_plan_t& plan)
{
    if (P.empty()) return 0;
    auto ctx = index.template empty_context<COUNT>();
    return ext_search(P, plan, ctx) ? (uint64_t) ctx.num_occ() : 0;
}

/**
 * @brief matches P in the order of an extension plan and reports its occurrences
 * @param index move_rb or an index adapter (see indexes/index_adapter.hpp)
 * @param P the pattern
 * @param plan the extension plan
 * @param report called with each occurrence
 */
template <typename idx_t, typename report_fnc_t>
inline void ext_locate(const idx_t& index, const std::string& P, const ext_plan_t& plan, report_fnc_t report)
{
    if (P.empty()) return;
    auto ctx = index.template empty_context<LOCATE>();
    if (ext_search(P, plan, ctx)) ctx.locate_phase().locate(report);
}

#endif // MOVE_R_BENCH_MOVE_RB_BENCH_EXT_QUERY_HPP
