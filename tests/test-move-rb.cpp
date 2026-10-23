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

#include <algorithm>
#include <numeric>

#include <data_structures/fasta_sequence_data.hpp>

#include "test-move-rb-common.hpp"

// the largest input of the tests with error models (their reference aligns the pattern at every starting position)
inline uint64_t max_model_input_size = 400;
// the smallest value run_fuzz reduces max_model_input_size to
inline constexpr uint64_t min_model_fuzz_input_size = 64;
// the largest length of the repeated block of the inputs of apm_model_repetitive
inline uint64_t max_model_block_length = 24;
// the number of patterns apm_model_repetitive queries per index
inline uint64_t model_queries_per_index = 20;
// apm_model_repetitive replaces each symbol by a random one with probability 1 / mutation_inverse_probability
inline constexpr uint64_t mutation_inverse_probability = 10;
// apm_separator_not_largest: the largest number of errors k of a query, and the smallest number of records (Z occurs
// in every second one)
inline constexpr uint64_t separator_test_max_k = 2;
inline constexpr uint64_t separator_test_min_records = 2;
// apm_text_ends: the input length is in [min_text_end_input_size, min_text_end_input_size + max_text_end_input_size)
// (run_fuzz reduces max_text_end_input_size down to min_text_end_input_size), the number of patterns per index, and a
// pattern is cut from the start, from the end or from a random position of the input with probability
// 1 / text_end_pattern_inverse_probability each
inline uint64_t max_text_end_input_size = 200;
inline constexpr uint64_t min_text_end_input_size = 32;
inline constexpr uint64_t text_end_queries_per_index = 20;
inline constexpr uint64_t text_end_pattern_inverse_probability = 3;
// apm_text_ends: the largest number of errors k of a query
inline constexpr uint64_t text_end_max_k = 3;
// large_alphabet_extensions: the internal symbols of the marker and of the separator, and the smallest continuation
// symbol (above both)
inline constexpr uint32_t large_alphabet_marker = 3;
inline constexpr uint32_t large_alphabet_separator = 4;
inline constexpr uint32_t continuation_base = 10000;

// the index operations, for every support (see cycle_supports)
TEST(test_move_rb, index)
{
    run_fuzz("move-rb", { { "index",
        [](uint64_t i) { cycle_supports(i, _test_index_operations); } } },
        fuzz_iterations(), &max_input_size, min_fuzz_input_size);
}

// approximate matching
TEST(test_move_rb, hamming)
{
    run_fuzz("move-rb", {
        { "count-hamming",  [](uint64_t i) { cycle_multiple_locate_supports(i, _test_count_hamming); } },
        { "locate-hamming", [](uint64_t i) { cycle_multiple_locate_supports(i, _test_locate_hamming); } },
    }, fuzz_iterations(), &max_input_size, min_fuzz_input_size);
}

TEST(test_move_rb, edit)
{
    run_fuzz("move-rb", {
        { "locate-edit",       [](uint64_t i) { cycle_multiple_locate_supports(i, _test_locate_edit); } },
        { "locate-edit-cigar", [](uint64_t i) { cycle_multiple_locate_supports(i, _test_locate_edit_cigar); } },
    }, fuzz_iterations(), &max_input_size, min_fuzz_input_size);
}

TEST(test_move_rb, apm_model)
{
    auto run = [](move_rb_test_t test) {
        return [test](uint64_t i) { cycle_multiple_locate_supports(i, test, min_input_size, max_model_input_size); };
    };

    run_fuzz("move-rb", {
        { "count-hamming-model",  run(_test_count_hamming_model) },
        { "locate-hamming-model", run(_test_locate_hamming_model) },
        { "locate-edit-model",    run(_test_locate_edit_model) },
    }, fuzz_iterations(), &max_model_input_size, min_model_fuzz_input_size);
}

// the approximate-matching algorithms with random error models on short repetitive inputs over 2 to 4 symbols: copies
// of a random block, each symbol replaced by a random one with probability 1 / mutation_inverse_probability
TEST(test_move_rb, apm_model_repetitive)
{
    using idx_t = move_rb<_locate_move>;
    const std::string bases = "ACGT";

    run_fuzz("move-rb", { { "apm-model-repetitive", [&](uint64_t) {
        std::mt19937_64 gen_test(gen());
        const std::string syms = bases.substr(0,
            min_random_alphabet_size + gen_test() % (bases.size() - min_random_alphabet_size + 1));
        const uint64_t n = max_model_block_length + gen_test() % max_model_input_size;
        std::string block(1 + gen_test() % max_model_block_length, 'A');
        for (char& c : block) c = syms[gen_test() % syms.size()];
        std::string input(n, 'A');

        for (uint64_t i = 0; i < n; i++) {
            input[i] = gen_test() % mutation_inverse_probability == 0
                ? syms[gen_test() % syms.size()] : block[i % block.size()];
        }

        move_r_params params;
        params.mode = _suffix_array;
        params.num_threads = 1;
        idx_t index(input, params);
        std::span<const uint8_t> input_bytes = byte_span(input);

        for (uint64_t q = 0; q < model_queries_per_index; q++) {
            const uint64_t m = std::min<uint64_t>(n, 1 + gen_test() % max_model_block_length);
            std::string pattern = input.substr(gen_test() % (n - m + 1), m);
            apply_random_edits(pattern, input, gen_test);
            std::span<const uint8_t> pattern_bytes = byte_span(pattern);
            const bool edit = gen_test() % 2 == 0;
            std::string quality;
            auto costs = random_query_model(index, input_bytes, pattern_bytes, edit ? EDIT_DISTANCE : HAMMING_DISTANCE,
                max_errors_for(pattern.size(), edit), gen_test, quality);
            if (!costs) continue;

            const apm_params_t apm_params = random_apm_params(gen_test);

            if (edit) {
                verify_edit_model(index, pattern_bytes, input_bytes, *costs, apm_params, gen_test);
            } else {
                const bool locate = gen_test() % 2 == 0;
                verify_hamming_model(index, pattern_bytes, input_bytes, *costs, apm_params, locate, gen_test);
            }
        }
    } } }, fuzz_iterations());
}

/* locate with the Hamming and the edit distance (with and without CIGARs) with each search scheme for patterns at the
   text ends, which the edit-distance search and its text-end pass both find: each position is reported once, and the
   occurrences are those of the naive reference (see verify_locate_hamming and verify_locate_edit) */
TEST(test_move_rb, apm_text_ends)
{
    using idx_t = move_rb<_locate_move>;
    const std::string bases = "ACGT";

    run_fuzz("move-rb", { { "apm-text-ends", [&](uint64_t) {
        std::mt19937_64 gen_test(gen());
        const std::string syms = bases.substr(0,
            min_random_alphabet_size + gen_test() % (bases.size() - min_random_alphabet_size + 1));
        const uint64_t n = min_text_end_input_size + gen_test() % max_text_end_input_size;
        std::string input(n, 'A');
        for (char& c : input) c = syms[gen_test() % syms.size()];

        move_r_params params;
        params.mode = _suffix_array;
        params.num_threads = 1;
        idx_t index(input, params);
        std::span<const char> input_span(input);

        for (uint64_t q = 0; q < text_end_queries_per_index; q++) {
            const uint64_t m = std::min<uint64_t>(n, 1 + gen_test() % max_pattern_length);
            const uint64_t where = gen_test() % text_end_pattern_inverse_probability;
            const uint64_t pos = where == 0 ? 0 : (where == 1 ? n - m : gen_test() % (n - m + 1));
            std::string pattern = input.substr(pos, m);
            apply_random_edits(pattern, input, gen_test);
            std::span<const char> pattern_span(pattern);
            const apm_params_t apm_params = random_apm_params(gen_test);

            for (uint64_t k = 1; k <= text_end_max_k; k++) {
                for (const search_scheme_t& scheme : { min_u_scheme(k), pigeon_hole_scheme(k), suffix_filter_scheme(k),
                    zero_one_scheme(k) }) {
                    if (scheme.p >= pattern.size()) continue;
                    if (k <= max_errors_for(pattern.size(), false)) {
                        verify_locate_hamming(index, pattern, pattern_span, scheme, k, input_span, apm_params);
                    }
                    if (k > max_errors_for(pattern.size(), true)) continue;
                    verify_locate_edit<false>(index, pattern, pattern_span, scheme, k, input_span, apm_params);
                    verify_locate_edit<true>(index, pattern, pattern_span, scheme, k, input_span, apm_params);
                }
            }
        }
    } } }, fuzz_iterations(), &max_text_end_input_size, min_text_end_input_size);
}

// the parsing of the rules of an error model and their compilation for a pattern
TEST(test_move_rb, apm_model_compile)
{
    std::vector<uint8_t> alph_map(256, 0);
    alph_map['A'] = 1;
    alph_map['C'] = 2;
    alph_map['G'] = 3;
    alph_map['T'] = 4;
    auto rank = [&](char c) { return uint64_t(alph_map[uint8_t(c)] - 1); };

    EXPECT_EQ(parse_apm_range("-5:").absolute(10), std::make_pair(uint64_t(5), uint64_t(10)));
    EXPECT_EQ(parse_apm_range(":3").absolute(10), std::make_pair(uint64_t(0), uint64_t(3)));
    EXPECT_EQ(parse_apm_range("2").absolute(10), std::make_pair(uint64_t(2), uint64_t(3)));
    EXPECT_EQ(parse_apm_range("-1").absolute(10), std::make_pair(uint64_t(9), uint64_t(10)));
    EXPECT_EQ(parse_apm_range("3:-3").absolute(10), std::make_pair(uint64_t(3), uint64_t(7)));

    apm_cost_rule_t rule = parse_apm_cost_rule("si=2@-5:/C,T");
    EXPECT_EQ(rule.ops, APM_SUB | APM_INS);
    EXPECT_EQ(rule.cost, 2);
    EXPECT_EQ(rule.from, "C");
    EXPECT_EQ(rule.to, "T");
    EXPECT_EQ(parse_apm_cost_rule("d=inf").cost, apm_forbidden);
    EXPECT_EQ(parse_apm_counter("io=1@2:4").ops, APM_INS | APM_OPEN);
    EXPECT_THROW(parse_apm_cost_rule("x=1"), std::invalid_argument);
    EXPECT_THROW(parse_apm_counter("s=200"), std::invalid_argument);

    // P = ACGTRCGTAC: positions 2..4 are exact, T also matches C, R is an IUPAC code
    const std::string P = "ACGTRCGTAC";
    apm_model_t model;
    model.costs.emplace_back(parse_apm_cost_rule("d=2"));
    model.costs.emplace_back(apm_cost_rule_t{.ops = APM_SUB | APM_INS | APM_DEL, .cost = apm_forbidden,
        .range = parse_apm_range("2:5")});
    model.costs.emplace_back(parse_apm_cost_rule("o=3@-2:"));
    model.classes.emplace_back(parse_apm_class("T,C"));
    model.counters.emplace_back(parse_apm_counter("s=0@0:1"));
    model.iupac = true;
    apm_costs_t costs(model, P, "", false, 6, EDIT_DISTANCE, alph_map);

    EXPECT_EQ(costs.sub(2, rank('A')), apm_forbidden);
    EXPECT_EQ(costs.ins(4), apm_forbidden);
    EXPECT_EQ(costs.del(3), apm_forbidden);
    EXPECT_EQ(costs.del(4), apm_forbidden);
    EXPECT_EQ(costs.del(2), 2);
    EXPECT_EQ(costs.del(5), 2);
    EXPECT_EQ(costs.del(0), apm_forbidden);
    EXPECT_EQ(costs.del(10), apm_forbidden);
    EXPECT_EQ(costs.sub(0, rank('C')), apm_forbidden);
    EXPECT_TRUE(costs.matches(7, rank('C')));
    EXPECT_FALSE(costs.matches(5, rank('T')));
    EXPECT_TRUE(costs.matches(4, rank('A')) && costs.matches(4, rank('G')) && !costs.matches(4, rank('C')));
    EXPECT_EQ(costs.ins_open(8), 3);
    EXPECT_EQ(costs.ins_open(7), 0);
    EXPECT_EQ(costs.del_open(9), 3);
    EXPECT_TRUE(costs.is_affine());
    EXPECT_TRUE(costs.is_literal(2) && !costs.is_literal(3) && !costs.is_literal(4));
    EXPECT_EQ(costs.max_errors(2, 4), 0u);
    EXPECT_EQ(costs.max_errors(), 6u);

    // on the reverse strand, the range [0, 1) of the read is the last position of the pattern and T,C is A,G
    apm_costs_t reverse(model, std::string("GTACGYACGT"), "", true, 6, EDIT_DISTANCE, alph_map);
    EXPECT_EQ(reverse.sub(9, rank('A')), apm_forbidden);
    EXPECT_EQ(reverse.sub(0, rank('A')), 1);
    EXPECT_TRUE(reverse.matches(2, rank('G')));

    // a substitution costs the quality of its position
    apm_model_t quality_model;
    quality_model.quality = true;
    apm_costs_t quality(quality_model, std::string("ACG"), "!#%", false, 4, HAMMING_DISTANCE, alph_map);
    EXPECT_EQ(quality.sub(0, rank('C')), 0);
    EXPECT_EQ(quality.sub(1, rank('A')), 2);
    EXPECT_EQ(quality.sub(2, rank('A')), 4);
    EXPECT_EQ(quality.ins(0), apm_forbidden);
    EXPECT_EQ(quality.tau(), 1);

    // two of ten positions have quality 2: with tau = 3 they are free and a search allows 1 error instead of 2
    apm_costs_t low(quality_model, std::string("ACGTACGTAC"), "II#IIII#II", false, 40, HAMMING_DISTANCE, alph_map);
    EXPECT_EQ(low.tau(), 3);
    EXPECT_EQ(low.max_errors(), 1u);
    EXPECT_TRUE(low.has_free(2) && !low.has_free(3));

    // the cheapest alignment of P against a prefix of a window
    apm_model_t unit;
    apm_costs_t unit_costs(unit, std::string("ACGT"), "", false, 2, EDIT_DISTANCE, alph_map);
    cigar_t<uint8_t> cigar;
    EXPECT_EQ(unit_costs.align_prefix(std::string("ACGTT"), &cigar), std::make_pair(uint64_t(0), uint64_t(4)));
    EXPECT_EQ(unit_costs.align_prefix(std::string("AGTA"), &cigar), std::make_pair(uint64_t(1), uint64_t(3)));
    EXPECT_EQ(unit_costs.align_prefix(std::string("TTTT")).first, 3u);
}

/* A FASTA text whose separator is not the largest internal symbol: Z occurs once in every second record, so it is
   rarer than the separator and is mapped after it. A Hamming search must not extend with the separator, but must
   still extend with the symbols after it in the internal order. */
TEST(test_move_rb, apm_separator_not_largest)
{
    constexpr uint64_t max_k = separator_test_max_k;

    using idx_t = move_rb<_locate_move>;
    std::mt19937_64 gen_test(gen());
    const uint8_t separator = fasta_separator_symbol;
    const uint64_t num_records = separator_test_min_records + gen_test() % (max_input_size / max_pattern_length);
    std::string input;
    std::vector<uint64_t> starts;
    std::vector<std::string> names;

    for (uint64_t r = 0; r < num_records; r++) {
        input += char(separator);
        starts.push_back(input.size());
        names.push_back("r" + std::to_string(r));
        const uint64_t length = 1 + gen_test() % max_pattern_length;
        const uint64_t z = r % 2 == 1 ? gen_test() % length : length;
        for (uint64_t i = 0; i < length; i++) input += i == z ? 'Z' : (gen_test() % 2 == 0 ? 'A' : 'C');
    }

    move_r_params params;
    params.mode = _suffix_array;
    params.num_threads = 1;
    params.reorder_alphabet = true;
    idx_t index(input, params);
    fasta_sequence_data sequences;
    sequences.set_sequences(input.size() + 1, starts, names, separator);
    index.set_fasta_sequence_data(std::move(sequences));
    ASSERT_LT(index.forward_index().map_symbol(char(separator)), index.forward_index().map_symbol('Z'))
        << "the separator is not mapped before Z";

    for (uint64_t q = 0; q < max_pattern_length; q++) {
        const uint64_t r = gen_test() % num_records;
        const uint64_t record_end = r + 1 < num_records ? starts[r + 1] - 1 : input.size();
        const uint64_t length = record_end - starts[r];
        const uint64_t k = 1 + gen_test() % max_k;
        const search_scheme_t scheme = min_u_scheme(k);
        // a pattern has more symbols than the search scheme has parts
        const uint64_t min_m = scheme.p + 1;
        if (length < min_m) continue;
        const uint64_t m = min_m + gen_test() % (length - min_m + 1);
        const uint64_t pos = starts[r] + gen_test() % (length - m + 1);
        std::string pattern = input.substr(pos, m);
        const std::string replacements = "ACZ";
        pattern[gen_test() % m] = replacements[gen_test() % replacements.size()];

        // the windows without a separator with at most k mismatches
        std::vector<std::pair<uint64_t, uint64_t>> correct;

        for (uint64_t p = 0; p + m <= input.size(); p++) {
            uint64_t mismatches = 0;
            bool crosses = false;

            for (uint64_t i = 0; i < m; i++) {
                crosses |= uint8_t(input[p + i]) == separator;
                mismatches += input[p + i] != pattern[i];
            }

            if (!crosses && mismatches <= k) correct.emplace_back(p, mismatches);
        }

        std::vector<std::pair<uint64_t, uint64_t>> reported;
        for (const auto& occ : index.locate(pattern, scheme, HAMMING_DISTANCE)) reported.emplace_back(occ.pos(), occ.err());
        std::sort(reported.begin(), reported.end());
        EXPECT_EQ(reported, correct) << "pattern " << pattern << ", k = " << k;
    }
}

/* The input has more than move_r_params::min_sigma_L_prev_min distinct symbols, so L_prev_min is built. The suffixes
   starting with marker are ordered by their distinct continuation and their preceding symbols alternate. The left
   extensions of marker span thousands of runs of L' and report two symbols, the right extensions report every
   continuation. */
TEST(test_move_rb, large_alphabet_extensions)
{
    using idx_t = move_rb<_locate_move>;
    constexpr uint32_t marker = large_alphabet_marker;
    constexpr uint32_t separator = large_alphabet_separator;
    const uint32_t num_continuations = uint32_t(move_r_params {}.min_sigma_L_prev_min);
    // the symbols of a group: a preceding symbol, the marker, a continuation and the separator
    constexpr uint64_t symbols_per_group = 4;

    std::vector<uint32_t> input;
    input.reserve(symbols_per_group * num_continuations);

    for (uint32_t i = 0; i < num_continuations; i++) {
        input.emplace_back(1 + (i & 1));
        input.emplace_back(marker);
        input.emplace_back(continuation_base + i);
        input.emplace_back(separator);
    }

    move_r_params params;
    params.mode = _suffix_array;
    params.num_threads = 1;
    params.a = min_balancing_parameter;
    params.reorder_alphabet = false;
    params.context_table_share = 0.0;
    idx_t index(input, params);

    auto verify = [&](const idx_t& queried, const std::string& label) {
        EXPECT_FALSE(queried.forward_index().L_prev_min().empty()) << label;
        EXPECT_FALSE(queried.backward_index().L_prev_min().empty()) << label;

        auto marker_context = [&]() {
            auto root = queried.template empty_context<COUNT>();
            auto [ctx, occurs] = root.template extend<LEFT>(marker);
            EXPECT_TRUE(occurs) << label;
            return ctx;
        };

        std::vector<uint64_t> left;
        auto left_ctx = marker_context();
        left_ctx.template for_each_extension<LEFT, NO_TERM>(
            [&](uint64_t sym, const auto&) { left.emplace_back(sym); return true; });
        EXPECT_EQ(left, (std::vector<uint64_t> { 1, 2 })) << label;

        std::vector<uint64_t> right;
        auto right_ctx = marker_context();
        right_ctx.template for_each_extension<RIGHT, NO_TERM>(
            [&](uint64_t sym, const auto&) { right.emplace_back(sym); return true; });

        std::vector<uint64_t> expected_right(num_continuations);
        std::iota(expected_right.begin(), expected_right.end(), continuation_base);
        EXPECT_EQ(right, expected_right) << label;
    };

    verify(index, "built");

    verify(serialize_and_load(index), "loaded");

    params.mode = _suffix_array_space;
    idx_t space_index(input, params);
    verify(space_index, "suffix-array-space");
}

// gives access to the SA-sample of a search context (see move_rb::sample_info_pack_t)
struct sample_info_access : move_rb<> {
    using pack_t = sample_info_pack_t;
    static constexpr sample_t sample_at_run_end = END;
};

// the fields of an SA-sample hold their largest values: o, i and x positions up to 2^48 - 1
TEST(test_move_rb, sample_field_widths)
{
    sample_info_access::pack_t sample {};
    sample.set_d(BOTH_DIRS);
    sample.set_t(sample_info_access::sample_at_run_end);
    sample.set_o(uint48_t::max);
    sample.set_i(uint48_t::max);
    sample.set_x(uint48_t::max);
    EXPECT_EQ(sample.d(), BOTH_DIRS);
    EXPECT_EQ(sample.t(), sample_info_access::sample_at_run_end);
    EXPECT_EQ(uint64_t(sample.o()), uint48_t::max);
    EXPECT_EQ(uint64_t(sample.i()), uint48_t::max);
    EXPECT_EQ(uint64_t(sample.x()), uint48_t::max);
}
