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

#include <gtest/gtest.h>
#include <gtl/phmap.hpp>
#include <move_rb/move_rb.hpp>
#include <misc/apm.hpp>
#include <ips2ra.hpp>
#include "test-index-common.hpp"
#include "test-progress.hpp"

std::random_device rd;
std::mt19937 gen(rd());
uint16_t max_num_threads = omp_get_max_threads();
uint64_t min_input_size = 1;
uint64_t max_input_size = 3000;
uint64_t min_bigbwt_input_size = 60;
uint64_t max_pattern_length = 100;
uint64_t mismatches_limit = 10;

std::uniform_real_distribution<double> prob_distrib(0.0, 1.0);
std::lognormal_distribution<double> a_distrib(2.0, 3.0);

// the individual functionalities a move_rb index is fuzzed for
enum mrb_func { mrb_extract, mrb_count_hamming, mrb_locate_hamming, mrb_locate_edit, mrb_locate_edit_cigar };

// extraction: forward and backward text
template <typename sym_t, typename inp_t, typename pos_t, move_r_support support>
static void verify_move_rb_extract(const move_rb<support, sym_t, pos_t>& index, const inp_t& input)
{
    pos_t n = input.size();
    inp_t input_reversed(input.rbegin(), input.rend());

    auto verify_revert = [&](const auto& idx, const inp_t& expected) {
        inp_t reverted = idx.revert_range({ .num_threads = max_num_threads });
        ASSERT_EQ((pos_t) reverted.size(), n);
        #pragma omp parallel for num_threads(max_num_threads)
        for (pos_t i = 0; i < n; i++) EXPECT_EQ(reverted[i], expected[i]);
    };

    verify_revert(index.forward_index(), input);
    verify_revert(index.backward_index(), input_reversed);
}

// hamming-distance count: the reported occurrence count must equal the naive count
template <typename sym_t, typename inp_t, typename pos_t, move_r_support support>
static void verify_count_hamming(const move_rb<support, sym_t, pos_t>& index, const inp_t& pattern,
    std::span<const sym_t> pattern_span, const search_scheme_t& scheme, pos_t k, std::span<const sym_t> input_span)
{
    auto correct = locate<pos_t, HAMMING_DISTANCE>(input_span, pattern_span, k);
    EXPECT_EQ(index.count_hamming_dist(pattern, scheme), correct.size());
}

// hamming-distance locate: the reported occurrences must equal the naive occurrences exactly
template <typename sym_t, typename inp_t, typename pos_t, move_r_support support>
static void verify_locate_hamming(const move_rb<support, sym_t, pos_t>& index, const inp_t& pattern,
    std::span<const sym_t> pattern_span, const search_scheme_t& scheme, pos_t k, std::span<const sym_t> input_span)
{
    auto correct = locate<pos_t, HAMMING_DISTANCE>(input_span, pattern_span, k);
    auto occ = index.template locate<HAMMING_DISTANCE>(pattern, scheme);
    ips2ra::sort(occ.begin(), occ.end(), [](const auto& o){ return o.pos; });
    EXPECT_EQ(occ.size(), correct.size());
    EXPECT_EQ(occ, correct);
}

// edit-distance locate: every reported occurrence has the claimed edit distance (<= k),
// and the reported set covers every naive occurrence both before and after filtering
template <typename sym_t, typename inp_t, typename pos_t, move_r_support support>
static void verify_locate_edit(const move_rb<support, sym_t, pos_t>& index, const inp_t& pattern,
    std::span<const sym_t> pattern_span, const search_scheme_t& scheme, pos_t k, std::span<const sym_t> input_span)
{
    pos_t input_size = input_span.size();
    auto correct = locate<pos_t, EDIT_DISTANCE>(input_span, pattern_span, k);
    auto occ = index.template locate<EDIT_DISTANCE>(pattern, scheme);
    ips2ra::sort(occ.begin(), occ.end(), [](const auto& o){ return o.pos; });

    gtl::flat_hash_map<std::span<const sym_t>, pos_t, span_hash, span_eq> dist_cache;
    for (const auto& o : occ) {
        if (o.pos + o.len > input_size) { ADD_FAILURE(); continue; }
        std::span<const sym_t> match = input_span.subspan(o.pos, o.len);
        auto [it, inserted] = dist_cache.try_emplace(match, 0);
        if (inserted) it->second = edit_dist_bounded<pos_t>(pattern_span, match, k);
        EXPECT_TRUE(it->second == o.err && o.err <= k);
    }

    EXPECT_TRUE(verify_edit_distance_coverage(occ, correct, k));
    filter_edit_distance_occurrences(occ, k);
    EXPECT_TRUE(verify_edit_distance_coverage(occ, correct, k));
}

// edit-distance locate with CIGARs: additionally every occurrence's CIGAR must validly align
// the pattern against the matched text with exactly the claimed number of edits
template <typename sym_t, typename inp_t, typename pos_t, move_r_support support>
static void verify_locate_edit_cigar(const move_rb<support, sym_t, pos_t>& index, const inp_t& pattern,
    std::span<const sym_t> pattern_span, const search_scheme_t& scheme, pos_t k, std::span<const sym_t> input_span)
{
    pos_t input_size = input_span.size();
    auto correct = locate<pos_t, EDIT_DISTANCE>(input_span, pattern_span, k);

    std::vector<aprx_occ_t<pos_t, CIGAR>> cigar_results;
    // in CIGAR mode locate returns the CIGARs (one per SA-interval); each occurrence's o.cig_idx indexes into them
    std::vector<cigar_t<sym_t>> cigars = index.template locate<EDIT_DISTANCE, CIGAR>(pattern, scheme, [&](aprx_occ_t<pos_t, CIGAR> o){
        cigar_results.emplace_back(std::move(o));
    });

    std::vector<aprx_occ_t<pos_t>> cigar_occurrences;
    gtl::flat_hash_map<std::span<const sym_t>, pos_t, span_hash, span_eq> dist_cache;
    for (const auto& o : cigar_results) {
        cigar_occurrences.push_back({o.pos, o.len, o.err});
        if (o.pos + o.len > input_size) { ADD_FAILURE(); continue; }
        std::span<const sym_t> match = input_span.subspan(o.pos, o.len);
        EXPECT_EQ(num_cigar_edits(pattern_span.data(), pattern_span.size(), match.data(), match.size(), cigars[o.cig_idx]), o.err);
        auto [it, inserted] = dist_cache.try_emplace(match, 0);
        if (inserted) it->second = edit_dist_bounded<pos_t>(pattern_span, match, k);
        EXPECT_TRUE(it->second == o.err && o.err <= k);
    }

    ips2ra::sort(cigar_occurrences.begin(), cigar_occurrences.end(), [](const auto& o){ return o.pos; });
    EXPECT_TRUE(verify_edit_distance_coverage(cigar_occurrences, correct, k));
    filter_edit_distance_occurrences(cigar_occurrences, k);
    EXPECT_TRUE(verify_edit_distance_coverage(cigar_occurrences, correct, k));
}

// builds one random move_rb index and verifies exactly one functionality of it: either its extraction queries,
// or (in a parallel query loop over random patterns) one approximate-match query type
template <typename sym_t, typename inp_t, typename pos_t, move_r_support support>
void test_move_rb(mrb_func func)
{
    // reserve the three largest symbol values (needed by Big-BWT and move_r construction); a byte alphabet then
    // has at most 256 - 3 == 253 distinct symbols, matching move_r's distinct-symbol limit (256 - min_valid_char)
    sym_t min_sym = std::numeric_limits<sym_t>::min();
    sym_t max_sym = std::numeric_limits<sym_t>::max() - 3;

    inp_t input = random_repetitive_input<inp_t>(min_input_size, max_input_size, min_sym, max_sym);
    pos_t input_size = input.size();

    // build move-r with a random number of threads and balancing parameter, using Big-BWT half the time for
    // large enough string inputs (falling back to the in-memory suffix-array construction if it is unavailable)
    auto make = [&](move_r_construction_mode mode) {
        return move_rb<support, sym_t, pos_t>(input, {
            .mode = mode,
            .num_threads = max_num_threads,
            .a = std::min<uint16_t>(2 + a_distrib(gen), 32767)
        });
    };

    auto build_index = [&]() {
        if constexpr (std::is_same_v<inp_t, std::string>) {
            std::uniform_int_distribution<int> use_bigbwt(0, 15);
            if (input.size() >= min_bigbwt_input_size && use_bigbwt(gen) == 0) {
                inp_t input_backup = input;
                try {
                    auto index = make(_bigbwt);
                    input = input_backup;
                    return index;
                } catch (const std::exception&) {
                    input = input_backup;
                }
            }
        }
        return make(_suffix_array);
    };
    move_rb<support, sym_t, pos_t> index = build_index();

    if (func == mrb_extract) {
        verify_move_rb_extract(index, input);
        return;
    }

    // a read-only view for the naive reference algorithms; created AFTER the build, which reallocates input's
    // buffer (its contents are preserved) -- a view captured earlier would dangle
    std::span<const sym_t> input_span(input);

    std::uniform_int_distribution<pos_t> pattern_pos_distrib(0, input_size - 1);
    std::uniform_int_distribution<pos_t> pattern_length_distrib(1, std::min<pos_t>(max_pattern_length, input_size));
    // enough queries (spread over the threads) that the approximate-match queries -- not the construction --
    // dominate the per-iteration time
    pos_t num_queries = std::min<pos_t>(1'500, std::max<pos_t>(300, input_size / 4));

    #pragma omp parallel for num_threads(max_num_threads), schedule(dynamic,1)
    for (pos_t cur_query = 0; cur_query < num_queries; cur_query++) {
        std::random_device rd_thr;
        std::mt19937 gen_thr(rd_thr());
        pos_t pattern_pos = pattern_pos_distrib(gen_thr);
        pos_t pattern_length = std::min<pos_t>(input_size - pattern_pos, pattern_length_distrib(gen_thr));
        inp_t pattern;
        no_init_resize(pattern, pattern_length);
        for (pos_t i = 0; i < pattern_length; i++)
            pattern[i] = input[pattern_pos + i];

        mutate_pattern<pos_t>(pattern, input, gen_thr);
        std::span<const sym_t> pattern_span(pattern);
        // bound the error budget by the pattern length: edit distance is bounded by pattern_length / 3 (the
        // edit search cost grows steeply with k), hamming distance by the laxer pattern_length / 2
        const bool is_edit = (func == mrb_locate_edit || func == mrb_locate_edit_cigar);
        const int32_t max_k = std::min<int32_t>(
            is_edit ? int32_t(pattern.size()) / 3 : int32_t(pattern.size()) / 2 - 1, mismatches_limit);
        pos_t k = std::uniform_int_distribution<pos_t>(0, std::max<int32_t>(max_k, 0))(gen_thr);
        const search_scheme_t scheme = min_u_scheme(k);

        switch (func) {
            case mrb_count_hamming:     verify_count_hamming    (index, pattern, pattern_span, scheme, k, input_span); break;
            case mrb_locate_hamming:    verify_locate_hamming   (index, pattern, pattern_span, scheme, k, input_span); break;
            case mrb_locate_edit:       verify_locate_edit      (index, pattern, pattern_span, scheme, k, input_span); break;
            case mrb_locate_edit_cigar: verify_locate_edit_cigar(index, pattern, pattern_span, scheme, k, input_span); break;
            case mrb_extract: break; // handled above
        }
    }
}

template <typename sym_t, typename inp_t>
static void move_rb_combination(uint64_t i, mrb_func func)
{
    auto run = [&]<move_r_support support>() {
        if (i % 2 == 0) test_move_rb<sym_t, inp_t, uint32_t, support>(func);
        else            test_move_rb<sym_t, inp_t, uint64_t, support>(func);
    };

    if ((i / 2) % 2 == 0) run.template operator()<_locate_move>();
    else                  run.template operator()<_locate_rlzsa>();
}

static void move_rb_functionality(uint64_t i, mrb_func func)
{
    switch ((i / 4) % 3) {
        case 0: move_rb_combination<char,    std::string>(i, func);          break;
        case 1: move_rb_combination<uint8_t, std::vector<uint8_t>>(i, func); break;
        case 2: move_rb_combination<int8_t,  std::vector<int8_t>>(i, func);  break;
    }
}

TEST(test_move_rb, extract)           { run_fuzz("move-rb", { { "extract",           [](uint64_t i) { move_rb_functionality(i, mrb_extract); } } }); }
TEST(test_move_rb, count_hamming)     { run_fuzz("move-rb", { { "count-hamming",     [](uint64_t i) { move_rb_functionality(i, mrb_count_hamming); } } }); }
TEST(test_move_rb, locate_hamming)    { run_fuzz("move-rb", { { "locate-hamming",    [](uint64_t i) { move_rb_functionality(i, mrb_locate_hamming); } } }); }
TEST(test_move_rb, locate_edit)       { run_fuzz("move-rb", { { "locate-edit",       [](uint64_t i) { move_rb_functionality(i, mrb_locate_edit); } } }); }
TEST(test_move_rb, locate_edit_cigar) { run_fuzz("move-rb", { { "locate-edit-cigar", [](uint64_t i) { move_rb_functionality(i, mrb_locate_edit_cigar); } } }); }
