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

#include <random>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <move_r/move_r.hpp>
#include <rlzsa_opt/rlzsa_opt.hpp>
#include <misc/strings.hpp>
#include <misc/utils.hpp>
#include "test-index-common.hpp"
#include "test-progress.hpp"

std::random_device rd;
std::mt19937 gen(rd());
uint16_t max_num_threads = omp_get_max_threads();

/**
 * @brief builds a move_r index (with an rlzsa) over a random repetitive input and verifies every query of its
 * rlzsa decode context (rlzsa_opt::decode_context_t) in parallel against the reference suffix array: decode(),
 * init_left()/init_right(), turn_left()/turn_right(), next()/prev(), skip_left()/skip_right(),
 * report_left()/report_right() and the pos()/value() + set_pos()/set_value() accessors. The decode context is
 * generic over the index integer type pos_t only (it decodes suffix-array positions, independent of the input
 * symbol type), so the input value type is fixed and only pos_t is varied.
 * @tparam pos_t index integer type (uint32_t or uint64_t)
 */
enum decode_op { decode_report, decode_step, decode_skip, decode_turn };

// each of the following verifies one group of rlzsa decode-context operations for one query (1 <= l < c <= r
// <= n-1), against the reference suffix array sa; verify_decode_context below just dispatches to the selected one

// range reporting: report_right decodes SA[l..r] left-to-right, report_left decodes SA[c-1..l] right-to-left
template <typename pos_t>
static void verify_decode_report(const rlzsa_opt<pos_t>& rlz, const std::vector<int64_t>& sa, pos_t l, pos_t c, pos_t r)
{
    {
        auto dec = rlz.decode();
        dec.init_right(l);
        dec.set_value(sa[l - 1]);
        EXPECT_EQ(dec.pos(), l);

        pos_t cnt = 0;
        dec.report_right(r, [&](pos_t p, pos_t v) {
            EXPECT_EQ(v, sa[p]);
            cnt++;
        });
        EXPECT_EQ(cnt, r - l + 1);
    }

    {
        auto dec = rlz.decode();
        dec.init_left(c);
        dec.set_value(sa[c]);
        EXPECT_EQ(dec.pos(), c);
        EXPECT_EQ(dec.value(), sa[c]);

        std::vector<pos_t> res;
        dec.report_left(l, [&](pos_t v) { res.push_back(v); });
        EXPECT_EQ(res.size(), c - l);
        for (pos_t k = 0; k < res.size() && k < c - l; k++)
            EXPECT_EQ(res[k], sa[c - 1 - k]);
    }
}

// single-stepping: next() decodes SA[l..r] one value forward, prev() decodes SA[c-1..l] one value backward,
// plus the set_pos()/set_value() accessor round-trip
template <typename pos_t>
static void verify_decode_step(const rlzsa_opt<pos_t>& rlz, const std::vector<int64_t>& sa, pos_t l, pos_t c, pos_t r)
{
    {
        auto dec = rlz.decode();
        dec.init_right(l);
        dec.set_value(sa[l - 1]);
        for (pos_t p = l; p <= r; p++)
            EXPECT_EQ(dec.next(), sa[p]);
    }

    {
        auto dec = rlz.decode();
        dec.init_left(c);
        dec.set_value(sa[c]);
        for (pos_t p = c; p > l; p--)
            EXPECT_EQ(dec.prev(), sa[p - 1]);
    }

    {
        auto dec = rlz.decode();
        dec.set_pos(c);
        EXPECT_EQ(dec.pos(), c);
        dec.set_value(sa[c]);
        EXPECT_EQ(dec.value(), sa[c]);
    }
}

// skipping: skip_right jumps l -> c then decodes SA[c..r]; skip_left jumps r -> c then decodes SA[c-1..l]
template <typename pos_t>
static void verify_decode_skip(const rlzsa_opt<pos_t>& rlz, const std::vector<int64_t>& sa, pos_t l, pos_t c, pos_t r)
{
    {
        auto dec = rlz.decode();
        dec.init_right(l);
        dec.set_value(sa[l - 1]);
        dec.skip_right(c);
        EXPECT_EQ(dec.pos(), c);
        EXPECT_EQ(dec.value(), sa[c - 1]);

        pos_t cnt = 0;
        dec.report_right(r, [&](pos_t v) {
            EXPECT_EQ(v, sa[c + cnt]);
            cnt++;
        });
        EXPECT_EQ(cnt, r - c + 1);
    }

    {
        auto dec = rlz.decode();
        dec.init_left(r);
        dec.set_value(sa[r]);
        dec.skip_left(c);
        EXPECT_EQ(dec.pos(), c);
        EXPECT_EQ(dec.value(), sa[c]);

        std::vector<pos_t> res;
        dec.report_left(l, [&](pos_t v) { res.push_back(v); });
        EXPECT_EQ(res.size(), c - l);
        for (pos_t k = 0; k < res.size() && k < c - l; k++)
            EXPECT_EQ(res[k], sa[c - 1 - k]);
    }
}

// direction turning: turn_left makes a right-ready context left-ready (then prev()), turn_right the converse
// (then next()), and turn_left+turn_right is the identity
template <typename pos_t>
static void verify_decode_turn(const rlzsa_opt<pos_t>& rlz, const std::vector<int64_t>& sa, pos_t l, pos_t c, pos_t r)
{
    {
        auto dec = rlz.decode();
        dec.init_right(c);
        dec.turn_left();
        EXPECT_EQ(dec.pos(), c);
        dec.set_value(sa[c]);
        for (pos_t p = c; p > l; p--)
            EXPECT_EQ(dec.prev(), sa[p - 1]);
    }

    {
        auto dec = rlz.decode();
        dec.init_left(c);
        dec.turn_right();
        EXPECT_EQ(dec.pos(), c);
        dec.set_value(sa[c - 1]);
        for (pos_t p = c; p <= r; p++)
            EXPECT_EQ(dec.next(), sa[p]);
    }

    {
        auto dec = rlz.decode();
        dec.init_right(c);
        dec.turn_left();
        dec.turn_right();
        dec.set_value(sa[c - 1]);
        for (pos_t p = c; p <= r; p++)
            EXPECT_EQ(dec.next(), sa[p]);
    }
}

// builds a move_r index (with an rlzsa) over a random repetitive input and verifies one group of decode-context
// operations in parallel against the reference suffix array
template <typename pos_t>
static void verify_decode_context(decode_op op)
{
    char min_sym = std::numeric_limits<char>::min();
    char max_sym = std::numeric_limits<char>::max() - 2;

    // a repetitive input over a small, random alphabet, so the rlzsa contains both literal and copy phrases
    std::string input = random_repetitive_input<std::string>(500, 165000, min_sym, max_sym);
    uint64_t input_size = input.size();

    // build from a copy (the constructor appends a sentinel) and keep input for the reference suffix array
    move_r<_locate_rlzsa, char, pos_t> index(input, { .mode = _suffix_array });

    auto sa = compute_reference_sa<char, std::string>(input, input_size, max_num_threads);
    const rlzsa_opt<pos_t>& rlz = index.rlzsa();
    ASSERT_FALSE(rlz.empty());

    pos_t n = rlz.input_size();
    ASSERT_EQ(n, sa.size());

    uint64_t total_queries = std::min<uint64_t>(2000, std::max<uint64_t>(100, n / 20));
    uint64_t queries_per_thread = std::max<uint64_t>(1, total_queries / max_num_threads);

    #pragma omp parallel num_threads(max_num_threads)
    {
        std::random_device rd_thr;
        std::mt19937 gen_thr(rd_thr());

        for (uint64_t q = 0; q < queries_per_thread; q++) {
            pos_t l = std::uniform_int_distribution<pos_t>(1, n - 2)(gen_thr);
            pos_t r = std::uniform_int_distribution<pos_t>(l + 1, n - 1)(gen_thr);
            pos_t c = std::uniform_int_distribution<pos_t>(l + 1, r)(gen_thr);

            switch (op) {
                case decode_report: verify_decode_report<pos_t>(rlz, sa, l, c, r); break;
                case decode_step:   verify_decode_step<pos_t>(rlz, sa, l, c, r);   break;
                case decode_skip:   verify_decode_skip<pos_t>(rlz, sa, l, c, r);   break;
                case decode_turn:   verify_decode_turn<pos_t>(rlz, sa, l, c, r);   break;
            }
        }
    }
}

static void rlzsa_functionality(uint64_t i, decode_op op)
{
    if (i % 2 == 0) verify_decode_context<uint32_t>(op);
    else            verify_decode_context<uint64_t>(op);
}

TEST(test_rlzsa_opt, report) { run_fuzz("rlzsa-opt", { { "report", [](uint64_t i) { rlzsa_functionality(i, decode_report); } } }); }
TEST(test_rlzsa_opt, step)   { run_fuzz("rlzsa-opt", { { "step",   [](uint64_t i) { rlzsa_functionality(i, decode_step); } } }); }
TEST(test_rlzsa_opt, skip)   { run_fuzz("rlzsa-opt", { { "skip",   [](uint64_t i) { rlzsa_functionality(i, decode_skip); } } }); }
TEST(test_rlzsa_opt, turn)   { run_fuzz("rlzsa-opt", { { "turn",   [](uint64_t i) { rlzsa_functionality(i, decode_turn); } } }); }
