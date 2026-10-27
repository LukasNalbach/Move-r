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

// fuzz tests of move_rb (see test-move-rb.cpp)

#ifndef MOVE_R_TESTS_TEST_MOVE_RB_COMMON_HPP
#define MOVE_R_TESTS_TEST_MOVE_RB_COMMON_HPP

#include <filesystem>
#include <fstream>
#include <map>
#include <numeric>
#include <set>

#include <gtest/gtest.h>
#include <gtl/phmap.hpp>
#include <move_rb/move_rb.hpp>
#include <misc/apm.hpp>
#include <algorithms/move_rb/apm/model.hpp>
#include <ips2ra.hpp>
#include "test-references.hpp"
#include "test-patterns.hpp"
#include "test-fuzz.hpp"
#include "test-inputs.hpp"
#include "test-index-references.hpp"
#include "test-parameters.hpp"
#include "test-revert.hpp"

inline std::mt19937_64 gen(test_seed());
inline uint64_t max_num_threads = omp_get_max_threads();

// the smallest and the largest length of a random input (run_fuzz reduces max_input_size)
inline uint64_t min_input_size = 1;
inline uint64_t max_input_size = 3000;
// the smallest value run_fuzz reduces max_input_size to
inline constexpr uint64_t min_fuzz_input_size = 256;
// the largest number of distinct symbols of a random input of 32-bit symbols
inline constexpr uint64_t max_integer_alphabet_size = uint64_t(std::numeric_limits<int32_t>::max()) + 1;
// the balancing parameter a of a random index is min_balancing_parameter + a value drawn from a_distrib, at most
// max_balancing_parameter
inline std::lognormal_distribution<double> a_distrib(balancing_parameter_log_mean, balancing_parameter_log_deviation);

// the largest length of a random pattern
inline uint64_t max_pattern_length = 100;
// the inverse of the probability that a random pattern contains a symbol that does not occur in the input
inline constexpr uint64_t absent_symbol_inverse_probability = 8;
// the largest number of errors k of an approximate-matching query
inline uint64_t max_query_errors = 10;

// the number of consecutive fuzz iterations with the same support, and with the same input type (see cycle_supports
// and cycle_input_types)
inline constexpr uint64_t iterations_per_support = 2;
inline constexpr uint64_t iterations_per_input_type = 4;

// the largest cost of a cost rule of a random error model (other than apm_forbidden)
inline uint64_t max_model_cost = 3;
// the largest bound of a counter of a random error model
inline uint64_t max_model_counter = 3;
// the largest number of cost rules, counters and classes of a random error model
inline uint64_t max_model_rules = 3;
// the largest Phred quality of a random pattern quality
inline uint64_t max_model_quality = 5;
// a random cost rule of an error model forbids its operations with probability 1 / forbidden_rule_inverse_probability
inline constexpr uint64_t forbidden_rule_inverse_probability = 6;
// the costs of the substitutions of a random error model are the qualities of the positions with probability
// 1 / quality_model_inverse_probability
inline constexpr uint64_t quality_model_inverse_probability = 4;
// the largest tau of a random error model (0: chosen per pattern)
inline constexpr uint64_t max_model_tau = 2;
// the range of a random cost rule is the whole pattern with probability 1 / whole_range_inverse_probability
inline constexpr uint64_t whole_range_inverse_probability = 3;

// the largest number of errors of a query with a pattern of length m is m / min_symbols_per_edit_error for the edit
// distance and m / min_symbols_per_hamming_part - 1 for the Hamming distance (see max_errors_for)
inline constexpr uint64_t min_symbols_per_edit_error = 3;
inline constexpr uint64_t min_symbols_per_hamming_part = 2;

// verify_context_tables checks context_check_num_patterns patterns whose lengths are in [len, max(len,
// min(context_check_max_length_factor * len, context_check_max_length))] for the longest stored context length len
inline constexpr uint64_t context_check_num_patterns = 64;
inline constexpr uint64_t context_check_max_length_factor = 3;
inline constexpr uint64_t context_check_max_length = 40;

// verify_extend runs n / input_symbols_per_sequence sequences, in [min_num_sequences, max_num_sequences] (not
// thorough: num_sequences_not_thorough), of num_ops_per_sequence extensions each
inline constexpr uint64_t input_symbols_per_sequence = 8;
inline constexpr uint64_t min_num_sequences = 4;
inline constexpr uint64_t max_num_sequences = 200;
inline constexpr uint64_t num_sequences_not_thorough = 8;
inline constexpr uint64_t num_ops_per_sequence = 32;

// verify_algorithm runs n / input_symbols_per_query queries, in [min_num_queries, max_num_queries] (not thorough:
// num_queries_not_thorough)
inline constexpr uint64_t input_symbols_per_query = 4;
inline constexpr uint64_t min_num_queries = 16;
inline constexpr uint64_t max_num_queries = 1'500;
inline constexpr uint64_t num_queries_not_thorough = 8;

// symbol type of an input container
template <typename inp_t>
using sym_of = typename inp_t::value_type;

// the part of move_rb one fuzz iteration checks (see build_and_verify_move_rb)
enum move_rb_test_t {
    // the construction, extend, for_each_extension, prepare_extend_all and the context tables
    _test_index_operations,
    // count_hamming_dist
    _test_count_hamming,
    // locate with the Hamming distance
    _test_locate_hamming,
    // locate with the edit distance
    _test_locate_edit,
    // locate<CIGAR> with the edit distance
    _test_locate_edit_cigar,
    // count_hamming_dist and count_hamming_dist_histogram with a random error model (see random_apm_model)
    _test_count_hamming_model,
    // locate with the Hamming distance and a random error model
    _test_locate_hamming_model,
    // locate with the edit distance and a random error model
    _test_locate_edit_model
};

// true <=> test is an approximate-matching test (these are run with byte alphabets only)
inline constexpr bool is_apm_test(move_rb_test_t test)
{
    return test == _test_count_hamming || test == _test_locate_hamming
        || test == _test_locate_edit || test == _test_locate_edit_cigar
        || test == _test_count_hamming_model || test == _test_locate_hamming_model || test == _test_locate_edit_model;
}

// true <=> test is an approximate-matching test with the edit distance
inline constexpr bool is_edit_distance_test(move_rb_test_t test)
{
    return test == _test_locate_edit || test == _test_locate_edit_cigar || test == _test_locate_edit_model;
}

// true <=> test is an approximate-matching test with a random error model
inline constexpr bool is_model_test(move_rb_test_t test)
{
    return test == _test_count_hamming_model || test == _test_locate_hamming_model || test == _test_locate_edit_model;
}

/**
 * @brief returns the largest number of errors k of an approximate-matching query with a pattern of length m:
 *        m / min_symbols_per_edit_error for the edit distance, m / min_symbols_per_hamming_part - 1 for the Hamming
 *        distance (then the k + 1 parts of the search scheme min_u_scheme(k) have at least min_symbols_per_hamming_part
 *        symbols each), at most max_query_errors
 */
inline uint64_t max_errors_for(uint64_t m, bool edit_distance)
{
    const uint64_t max_k = edit_distance ? m / min_symbols_per_edit_error
        : (m < min_symbols_per_hamming_part ? 0 : m / min_symbols_per_hamming_part - 1);
    return std::min(max_k, max_query_errors);
}

// returns the symbols of a container of 1-byte symbols as bytes
template <typename inp_t>
inline std::span<const uint8_t> byte_span(const inp_t& container)
{
    static_assert(sizeof(sym_of<inp_t>) == 1);
    return { reinterpret_cast<const uint8_t*>(container.data()), container.size() };
}

/**
 * @brief returns a random query pattern: a substring of the input of length at most max_pattern_length, with random
 *        edits (see apply_random_edits); with probability 1 / absent_symbol_inverse_probability, one of its symbols is
 *        replaced by absent
 * @param absent a symbol that does not occur in the input
 */
template <typename inp_t>
inline inp_t random_query_pattern(const inp_t& input, sym_of<inp_t> absent, std::mt19937_64& gen_query)
{
    const uint64_t n = input.size();
    const uint64_t pattern_pos = std::uniform_int_distribution<uint64_t>(0, n - 1)(gen_query);
    const uint64_t max_length = std::min<uint64_t>(max_pattern_length, n);
    const uint64_t pattern_length =
        std::min<uint64_t>(n - pattern_pos, std::uniform_int_distribution<uint64_t>(1, max_length)(gen_query));
    inp_t pattern(input.begin() + pattern_pos, input.begin() + pattern_pos + pattern_length);
    apply_random_edits(pattern, input, gen_query);
    if (gen_query() % absent_symbol_inverse_probability == 0) pattern[gen_query() % pattern.size()] = absent;
    return pattern;
}

// returns random apm_params_t (the partitioning and the seeding must not change the result)
inline apm_params_t random_apm_params(std::mt19937_64& gen_query)
{
    return {
        .partitioning = partitioning_t(gen_query() % 2),
        .start_at_part_contexts = gen_query() % 2 == 0
    };
}

// returns a copy of index made by serializing and loading it
template <typename idx_t>
inline idx_t serialize_and_load(const idx_t& index)
{
    std::stringstream stream;
    index.serialize(stream);
    idx_t loaded;
    loaded.load(stream);
    return loaded;
}

/* the naive references of one index; the text is stored in the internal alphabet of the index, such that its
   suffixes are ordered as in the index */
struct naive_reference {
    std::vector<int32_t> text; // [0..n-1] the input in the internal alphabet of the index
    std::vector<int32_t> text_rev; // [0..n-1] text_rev[i] = text[n-1-i]
    std::vector<int64_t> SA; // [0..n] the suffix array of text (including the terminator)
    std::vector<int64_t> SA_rev; // [0..n] the suffix array of text_rev
    // [0..alphabet_size()+1] input_symbol[s] = the input symbol of the internal symbol s (see move_r::unmap_symbol)
    std::vector<uint64_t> input_symbol;
    // the internal symbol of each input symbol that occurs in the input
    std::map<uint64_t, uint64_t> internal_symbol;
};

// returns the naive references for an index built from input
template <typename idx_t, typename inp_t>
naive_reference build_naive_reference(const idx_t& index, const inp_t& input)
{
    uint64_t n = input.size();
    naive_reference ref;
    resize_for_overwrite(ref.text, n);

    #pragma omp parallel for num_threads(max_num_threads)
    for (uint64_t i = 0; i < n; i++) ref.text[i] = int32_t(index.map_symbol(input[i]));

    ref.SA = compute_reference_sa<int32_t, std::vector<int32_t>>(ref.text, n, max_num_threads);

    ref.text_rev.assign(ref.text.rbegin(), ref.text.rend());
    ref.SA_rev = compute_reference_sa<int32_t, std::vector<int32_t>>(ref.text_rev, n, max_num_threads);

    // the mappings between the internal symbols and the input symbols
    ref.input_symbol.resize(index.alphabet_size() + 2);
    for (uint64_t s = 0; s < ref.input_symbol.size(); s++) ref.input_symbol[s] = index.unmap_symbol(s);
    for (uint64_t i = 0; i < n; i++) ref.internal_symbol[ref.input_symbol[uint64_t(ref.text[i])]] = uint64_t(ref.text[i]);

    return ref;
}

// ############################# CONSTRUCTION #############################

// verifies the data structures of an index, the reverted text and the suffix array against the naive references
template <typename inp_t, move_r_support support>
void verify_construction(
    const move_rb<support>& index, const inp_t& input, const naive_reference& ref, const std::string& label)
{
    using idx_t = move_rb<support>;
    uint64_t n = input.size();

    // the index of T^R stores no M_Phi, no M_Phi^{-1} and no rlz encoding of SA^d
    using bwd_t = std::remove_cvref_t<decltype(index.backward_index())>;
    static_assert(!bwd_t::has_mphi);
    static_assert(!bwd_t::has_locate_move);
    static_assert(!bwd_t::has_SA_rlz);

    // a locating index keeps the SA samples of the index of T^R
    if constexpr (idx_t::supports_locate) {
        const auto& bwd = index.backward_index();
        uint64_t num_intervals = bwd.M_LF().num_intervals();

        // the SA samples are stored at the run boundaries only (see move_r::run_starts_at)
        for (uint64_t x = 0; x < num_intervals; x++) {
            if (bwd.run_starts_at(x)) EXPECT_LE(bwd.SA_run_start(x), n);
            if (bwd.run_ends_at(x)) EXPECT_LE(bwd.SA_run_end(x), n);
        }
    }

    // the index stores PSB_L' for every alphabet size
    EXPECT_TRUE(index.has_PSB_L_()) << "sigma = " << index.alphabet_size() + 1;

    // PSO_L' is always built
    EXPECT_FALSE(index.forward_index().PSO_L_().empty());

    // container type the text is reverted to; the width of its symbols must be the width of the input symbols
    // (see move_r::symbol_idx)
    using out_t = std::vector<std::make_unsigned_t<sym_of<inp_t>>>;

    auto verify_revert = [&](const auto& idx, bool reversed) {
        using usym_t = std::make_unsigned_t<sym_of<inp_t>>;

        // the index reverts to unsigned symbols, so a signed input symbol is reinterpreted as unsigned
        auto expected = [&](uint64_t i) {
            sym_of<inp_t> symbol = reversed ? input[n - 1 - i] : input[i];
            usym_t as_unsigned;
            std::memcpy(&as_unsigned, &symbol, sizeof(symbol));
            return uint64_t(as_unsigned);
        };

        verify_parallel_revert<usym_t>(n, expected,
            [&](uint64_t p, uint64_t k) { return idx.template revert<out_t>(p, k); },
            [&](const std::string& file, uint64_t p, uint64_t k) { idx.template revert<usym_t>(file, p, k); },
            std::string(reversed ? "T^R, " : "T, ") + label);
    };

    verify_revert(index.forward_index(), false);
    verify_revert(index.backward_index(), true);

    // the suffix array
    if constexpr (idx_t::supports_multiple_locate) {
        auto sa = index.forward_index().SA_range({ .num_threads = max_num_threads });
        ASSERT_EQ((uint64_t) sa.size(), n + 1);

        for (uint64_t i = 0; i <= n; i++)
            EXPECT_EQ(uint64_t(sa[i]), uint64_t(ref.SA[i])) << "SA[" << i << "], " << label;
    }
}

// ############################# INDEX OPERATIONS #############################

/* returns the input symbols the pattern with the suffix array interval [interval.first, interval.second] and the length
   len can be extended with in the direction dir, each with its number of occurrences (the terminator is not counted) */
inline std::map<uint64_t, uint64_t> reference_extensions(
    const naive_reference& ref, direction_t dir, std::pair<int64_t, int64_t> interval, uint64_t len)
{
    uint64_t n = ref.text.size();
    std::map<uint64_t, uint64_t> extensions;

    for (int64_t i = interval.first; i <= interval.second; i++) {
        uint64_t sa = uint64_t(ref.SA[i]);

        if (dir == LEFT) {
            if (sa == 0) continue; // the terminator is not reported
            extensions[ref.input_symbol[uint64_t(ref.text[sa - 1])]]++;
        } else {
            if (sa + len >= n) continue;
            extensions[ref.input_symbol[uint64_t(ref.text[sa + len])]]++;
        }
    }

    return extensions;
}

// verifies the enumeration of the symbols a pattern can be extended with in the direction dir, with
// for_each_extension() and with prepare_extend_all()/extend_next(); returns false <=> a check failed
template <direction_t dir, typename ctx_t>
bool verify_extensions(
    const naive_reference& ref, ctx_t& ctx, std::pair<int64_t, int64_t> interval, uint64_t len,
    const std::string& op, const std::string& label)
{
    std::map<uint64_t, uint64_t> correct = reference_extensions(ref, dir, interval, len);
    // the symbols reported by for_each_extension with their numbers of occurrences
    std::map<uint64_t, uint64_t> reported;
    // the symbols in the order for_each_extension reports them
    std::vector<uint64_t> order;

    ctx.template for_each_extension<dir, NO_TERM>([&](uint64_t sym, const auto& extended) {
        reported[sym] = extended.num_occ();
        order.emplace_back(sym);
        return true;
    });

    if (reported != correct) {
        ADD_FAILURE() << op << ": for_each_extension<" << (dir == LEFT ? "LEFT" : "RIGHT")
            << "> reported " << reported.size() << " symbols, " << correct.size() << " occur"
            << " (m = " << len << ", " << label << ")";
        return false;
    }

    // the same enumeration with prepare_extend_all() and extend_next()
    std::map<uint64_t, uint64_t> prepared;
    auto extending = ctx.extend_phase();
    extending.template prepare_extend_all<NO_TERM>(dir);

    while (extending.can_extend()) {
        uint64_t sym = extending.next_symbol();
        prepared[sym] = extending.extend_next().num_occ();
    }

    if (prepared != correct) {
        ADD_FAILURE() << op << ": the prepared enumeration in direction "
            << (dir == LEFT ? "LEFT" : "RIGHT") << " differs from for_each_extension"
            << " (m = " << len << ", " << label << ")";
        return false;
    }

    // the symbols must be reported in ascending order of the internal alphabet
    for (uint64_t i = 1; i < order.size(); i++) {
        auto rank_of = [&](uint64_t sym) {
            auto it = ref.internal_symbol.find(sym);
            return it == ref.internal_symbol.end() ? 0 : it->second;
        };

        if (rank_of(order[i - 1]) >= rank_of(order[i])) {
            ADD_FAILURE() << op << ": the extensions are not in the order of the internal alphabet"
                << " (m = " << len << ", " << label << ")";
            return false;
        }
    }

    return true;
}

// verifies the context tables: search() and context_of() must return the contexts that are reached by
// extending symbol by symbol
template <typename inp_t, move_r_support support>
void verify_context_tables(
    const move_rb<support>& index, const inp_t& input, const std::string& label)
{
    uint64_t n = input.size();
    // the length of the longest stored context
    uint64_t len = index.context_tables_max_length();
    // no context tables are stored, or the input is too short
    if (!index.has_context_tables() || len == 0 || n < len + 2) return;

    std::uniform_int_distribution<uint64_t> length_distrib(len,
        std::max<uint64_t>(len, std::min<uint64_t>(context_check_max_length_factor * len, context_check_max_length)));

    for (uint64_t q = 0; q < context_check_num_patterns; q++) {
        uint64_t m = std::min<uint64_t>(length_distrib(gen), n);
        std::uniform_int_distribution<uint64_t> pos_distrib(0, n - m);
        uint64_t pos = pos_distrib(gen);
        inp_t P(input.begin() + pos, input.begin() + pos + m);

        // the context of P, computed by extending symbol by symbol
        auto slow = index.template empty_context<COUNT>();
        bool occurs = true;

        for (uint64_t i = P.size(); i > 0 && occurs; i--) {
            auto [next, ok] = slow.template extend<LEFT>(P[i - 1]);
            occurs = ok;
            if (ok) slow = std::move(next);
        }

        ASSERT_TRUE(occurs) << "a substring of the input does not occur in it, " << label;

        // the context of P, computed with search() (uses the context tables, see move_rb::search)
        bool found = false;
        auto fast = index.template search<COUNT>(P, found);
        ASSERT_TRUE(found) << "search() missed a substring of the input, " << label;
        ASSERT_EQ(fast.num_occ(), slow.num_occ()) << "search() returned the context of another pattern";
        ASSERT_EQ(fast.pattern_length(), P.size());

        // context_of() returns the context of the longest stored suffix of P
        uint64_t matched = 0;
        auto entry = index.template context_of<COUNT>(P, matched);
        ASSERT_LE(matched, P.size()) << "the table matched more symbols than the pattern has";

        if (matched != 0) {
            auto extended = index.template empty_context<COUNT>();

            for (uint64_t i = P.size(); i > P.size() - matched; i--) {
                auto [next, ok] = extended.template extend<LEFT>(P[i - 1]);
                ASSERT_TRUE(ok);
                extended = std::move(next);
            }

            ASSERT_EQ(entry.num_occ(), extended.num_occ()) << "the context of the table belongs to another pattern";
            ASSERT_EQ(entry.pattern_length(), matched);
        }
    }
}

/**
 * @brief performs random sequences of extensions to the left and to the right over the substrings of the input; after
 *        every extension, the suffix array intervals of the matched pattern in T and T^R, its occurrences and the
 *        symbols it can be extended with are verified against the naive references
 * @param absent a symbol that does not occur in the input
 * @param thorough false <=> the index was loaded from a stream or built from a cache (fewer sequences are performed)
 */
template <typename inp_t, move_r_support support>
void verify_extend(const move_rb<support>& index, const inp_t& input,
    const naive_reference& ref, sym_of<inp_t> absent, const std::string& label, bool thorough)
{
    // the operations of a sequence
    enum operation_t : uint64_t { _extend_left, _extend_right };

    using idx_t = move_rb<support>;
    static constexpr query_support_t query_support = idx_t::supports_locate ? LOCATE : COUNT;
    uint64_t n = input.size();

    uint64_t num_sequences = thorough
        ? std::clamp<uint64_t>(n / input_symbols_per_sequence, min_num_sequences, max_num_sequences)
        : num_sequences_not_thorough;

    // [0..num_sequences-1] the seed of the random number generator of each sequence
    std::vector<uint64_t> seeds(num_sequences);
    for (uint64_t& seed : seeds) seed = gen();

    #pragma omp parallel for num_threads(max_num_threads), schedule(dynamic,1)
    for (uint64_t cur_sequence = 0; cur_sequence < num_sequences; cur_sequence++) {
        std::mt19937_64 gen_sequence(seeds[cur_sequence]);
        std::uniform_int_distribution<uint64_t> pos_distrib(0, n - 1);
        std::uniform_int_distribution<uint64_t> op_distrib(_extend_left, _extend_right);

        // the currently matched pattern is always the substring input[pos, pos + len), so every extension
        // succeeds
        uint64_t pos = pos_distrib(gen_sequence);
        uint64_t len = 0;
        auto ctx = index.template empty_context<query_support>();

        // verifies the context of the matched pattern; returns false <=> a check failed
        auto verify = [&](const std::string& op) {
            std::vector<int32_t> P(ref.text.begin() + pos, ref.text.begin() + pos + len);
            std::vector<int32_t> P_rev(P.rbegin(), P.rend());
            auto [b, e] = ctx.forward_sa_interval();
            auto [b_R, e_R] = ctx.backward_sa_interval();
            auto [b_ref, e_ref] = reference_sa_interval(ref.SA, ref.text, P);
            auto [b_R_ref, e_R_ref] = reference_sa_interval(ref.SA_rev, ref.text_rev, P_rev);

            bool correct_intervals =
                (int64_t) ctx.pattern_length() == (int64_t) len &&
                (int64_t) b == b_ref && (int64_t) e == e_ref &&
                (int64_t) b_R == b_R_ref && (int64_t) e_R == e_R_ref;
            EXPECT_TRUE(correct_intervals) << op << ": [" << b << "," << e << "] vs ["
                << b_ref << "," << e_ref << "], reverse [" << b_R << "," << e_R << "] vs ["
                << b_R_ref << "," << e_R_ref << "]"
                << " (" << label << ")";
            if (!correct_intervals) return false;

            // the symbols the pattern can be extended with
            if (!verify_extensions<LEFT>(ref, ctx, { b_ref, e_ref }, len, op, label)) return false;
            if (!verify_extensions<RIGHT>(ref, ctx, { b_ref, e_ref }, len, op, label)) return false;

            // the pattern cannot be extended with a symbol that does not occur in the input
            if (!ref.internal_symbol.contains(symbol_bit_pattern(absent))) {
                auto [ctx_absent, occurs_absent] = ctx.template extend<LEFT>(absent);
                EXPECT_FALSE(occurs_absent) << op << ": extend<LEFT> by a symbol that does not occur, "
                    << label;
                if (occurs_absent) return false;

                auto [ctx_r, occurs_r] = ctx.template extend<RIGHT>(absent);
                EXPECT_FALSE(occurs_r) << op << ": extend<RIGHT> by a symbol that does not occur, "
                    << label;
                if (occurs_r) return false;
            }

            // the occurrences (the empty pattern occurs at every position and is skipped)
            if constexpr (query_support == LOCATE) {
                if (len != 0) {
                    if constexpr (idx_t::supports_multiple_locate) {
                        auto reported = ctx.locate_phase().locate();
                        std::vector<uint64_t> occ(reported.begin(), reported.end());
                        std::vector<uint64_t> correct;

                        for (int64_t i = b_ref; i <= e_ref; i++) correct.emplace_back(ref.SA[i]);
                        std::sort(occ.begin(), occ.end());
                        std::sort(correct.begin(), correct.end());
                        EXPECT_EQ(occ, correct) << op;
                        if (occ != correct) return false;
                    } else {
                        // sampled_occ() must report one of the occurrences
                        uint64_t occ = ctx.locate_phase().sampled_occ();
                        bool is_occurrence = false;

                        for (int64_t i = b_ref; i <= e_ref; i++)
                            if (uint64_t(ref.SA[i]) == occ) is_occurrence = true;

                        EXPECT_TRUE(is_occurrence) << op << ": " << occ
                            << " is not an occurrence of the matched pattern, " << label;
                        if (!is_occurrence) return false;
                    }
                }
            }

            return true;
        };

        for (uint64_t cur_op = 0; cur_op < num_ops_per_sequence; cur_op++) {
            bool valid = true;

            switch (operation_t(op_distrib(gen_sequence))) {
                case _extend_left: {
                    if (pos == 0) continue;
                    pos--; len++;
                    auto [ctx_ext, occurs] = ctx.template extend<LEFT>(input[pos]);
                    EXPECT_TRUE(occurs) << "extend<LEFT> at " << pos << ", " << label;
                    if (!occurs) { valid = false; break; }
                    ctx = ctx_ext;
                    valid = verify("extend<LEFT>");
                    break;
                }
                case _extend_right: {
                    if (pos + len == n) continue;
                    len++;
                    auto [ctx_ext, occurs] = ctx.template extend<RIGHT>(input[pos + len - 1]);
                    EXPECT_TRUE(occurs) << "extend<RIGHT> at " << (pos + len - 1) << ", " << label;
                    if (!occurs) { valid = false; break; }
                    ctx = ctx_ext;
                    valid = verify("extend<RIGHT>");
                    break;
                }
            }

            if (!valid) break;
        }
    }

    if (thorough) verify_context_tables(index, input, label);
}

// ############################# ALGORITHMS #############################

// count_hamming_dist must return the number of naive occurrences with Hamming distance at most k
template <typename inp_t, move_r_support support>
void verify_count_hamming(const move_rb<support>& index, const inp_t& pattern,
    std::span<const sym_of<inp_t>> pattern_span, const search_scheme_t& scheme, uint64_t k,
    std::span<const sym_of<inp_t>> input_span, apm_params_t params)
{
    auto correct = locate_naive(input_span, pattern_span, k, HAMMING_DISTANCE);
    EXPECT_EQ(index.count_hamming_dist(pattern, scheme, params), correct.size());
}

// locate with the Hamming distance must report the naive occurrences with Hamming distance at most k
template <typename inp_t, move_r_support support>
void verify_locate_hamming(const move_rb<support>& index, const inp_t& pattern,
    std::span<const sym_of<inp_t>> pattern_span, const search_scheme_t& scheme, uint64_t k,
    std::span<const sym_of<inp_t>> input_span, apm_params_t params)
{
    auto correct = locate_naive(input_span, pattern_span, k, HAMMING_DISTANCE);
    auto occ = index.locate(pattern, scheme, HAMMING_DISTANCE, params);
    ips2ra::sort(occ.begin(), occ.end(), key_of_pos {});
    EXPECT_EQ(occ.size(), correct.size());
    EXPECT_EQ(occ, correct);
}

/**
 * @brief checks locate with the edit distance: no position is reported twice, every reported occurrence has the
 *        reported edit distance (at most k), and the occurrences cover every naive occurrence before and after
 *        filtering them (see filter_edit_distance_occurrences)
 * @tparam with_cigar true <=> locate<CIGAR> is called, and the CIGAR of each occurrence must be an alignment with the
 *         reported number of edits
 */
template <bool with_cigar, typename inp_t, move_r_support support>
void verify_locate_edit(const move_rb<support>& index, const inp_t& pattern,
    std::span<const sym_of<inp_t>> pattern_span, const search_scheme_t& scheme, uint64_t k,
    std::span<const sym_of<inp_t>> input_span, apm_params_t params)
{
    using sym_t = sym_of<inp_t>;
    uint64_t input_size = input_span.size();
    auto correct = locate_naive(input_span, pattern_span, k, EDIT_DISTANCE);
    std::vector<aprx_occ_t<>> occ;

    if constexpr (with_cigar) {
        std::vector<aprx_occ_t<CIGAR>> occ_cigar;
        // the CIGARs (one per suffix array interval); o.cig_idx is the index of the CIGAR of the occurrence o
        auto cigars = index.template locate<CIGAR>(pattern, scheme, EDIT_DISTANCE, occ_cigar, params);

        for (const auto& o : occ_cigar) {
            occ.push_back({o.pos(), o.len(), o.err()});
            if (o.pos() + o.len() > input_size) continue; // reported below
            EXPECT_EQ(num_cigar_edits(pattern_span.data(), pattern_span.size(), input_span.data() + o.pos(), o.len(),
                cigars[o.cig_idx]), o.err());
        }
    } else {
        occ = index.locate(pattern, scheme, EDIT_DISTANCE, params);
    }

    ips2ra::sort(occ.begin(), occ.end(), key_of_pos {});

    for (uint64_t x = 1; x < occ.size(); x++) {
        EXPECT_NE(occ[x - 1].pos(), occ[x].pos()) << "position " << occ[x].pos() << " is reported twice";
    }

    // the edit distance of the pattern to each distinct matched substring
    gtl::flat_hash_map<std::span<const sym_t>, uint64_t, span_hash, span_eq> dist_cache;

    for (const auto& o : occ) {
        if (o.pos() + o.len() > input_size) { ADD_FAILURE(); continue; }
        std::span<const sym_t> match = input_span.subspan(o.pos(), o.len());
        auto [it, inserted] = dist_cache.try_emplace(match, 0);
        if (inserted) it->second = edit_dist_bounded(pattern_span, match, k);
        EXPECT_TRUE(it->second == o.err() && o.err() <= k);
    }

    EXPECT_TRUE(verify_edit_distance_coverage(occ, correct, k));
    filter_edit_distance_occurrences(occ, k);
    EXPECT_TRUE(verify_edit_distance_coverage(occ, correct, k));
}

// returns a random range of positions of a pattern of length m (see apm_range_t)
inline apm_range_t random_apm_range(uint64_t m, std::mt19937_64& gen_query)
{
    if (gen_query() % whole_range_inverse_probability == 0) return {};
    int64_t beg = int64_t(gen_query() % (m + 1));
    int64_t end = beg + int64_t(gen_query() % (m + 1 - beg));
    if (gen_query() % 2 == 0) beg -= int64_t(m);
    if (gen_query() % 2 == 0) end -= int64_t(m);
    return {beg, end};
}

// returns a random non-empty subset of the symbols syms
inline std::string random_symbols(const std::string& syms, std::mt19937_64& gen_query)
{
    std::string out;
    for (char c : syms) if (gen_query() % 2 == 0) out += c;
    if (out.empty()) out += syms[gen_query() % syms.size()];
    return out;
}

// returns a random error model over the symbols syms of the text for a pattern of length m
inline apm_model_t random_apm_model(const std::string& syms, uint64_t m, std::mt19937_64& gen_query)
{
    // the largest cost_dp of a random model
    constexpr uint64_t max_model_cost_dp = APM_COST_DP_OFF;

    apm_model_t model;
    const uint8_t ops_all = APM_SUB | APM_INS | APM_DEL | APM_OPEN;

    for (uint64_t r = gen_query() % (max_model_rules + 1); r > 0; r--) {
        apm_cost_rule_t rule;
        rule.ops = uint8_t(1 + gen_query() % ops_all);
        rule.cost = gen_query() % forbidden_rule_inverse_probability == 0
            ? apm_forbidden : uint16_t(1 + gen_query() % max_model_cost);
        rule.range = random_apm_range(m, gen_query);

        if (rule.ops == APM_SUB && gen_query() % 2 == 0) {
            rule.from = random_symbols(syms, gen_query);
            rule.to = random_symbols(syms, gen_query);
        }

        model.costs.emplace_back(rule);
    }

    for (uint64_t r = gen_query() % (max_model_rules + 1); r > 0; r--) {
        model.counters.emplace_back(apm_counter_t{.ops = uint8_t(1 + gen_query() % ops_all),
            .max = uint8_t(gen_query() % (max_model_counter + 1)), .range = random_apm_range(m, gen_query)});
    }

    for (uint64_t r = gen_query() % (max_model_rules + 1); r > 0; r--) {
        model.classes.emplace_back(apm_class_t{.from = random_symbols(syms, gen_query),
            .to = random_symbols(syms, gen_query), .range = random_apm_range(m, gen_query)});
    }

    model.quality = gen_query() % quality_model_inverse_probability == 0;
    model.tau = uint16_t(gen_query() % (max_model_tau + 1));
    model.cost_dp = apm_cost_dp_t(gen_query() % (max_model_cost_dp + 1));
    return model;
}

// returns a random Phred+33 quality string of length m
inline std::string random_quality(uint64_t m, std::mt19937_64& gen_query)
{
    std::string qual(m, '!');
    for (char& c : qual) c = char('!' + gen_query() % (max_model_quality + 1));
    return qual;
}

// true <=> cnt[j] <= the bound of counter j of the model, for every counter j
inline bool within_counters(const apm_costs_t& costs, const std::vector<uint64_t>& cnt)
{
    for (uint64_t j = 0; j < costs.num_counters(); j++) if (cnt[j] > costs.counter_bound(j)) return false;
    return true;
}

// adds the counter increments packed into inc (the increment of counter j in byte j) to the counters cnt
inline void add_counts(std::vector<uint64_t>& cnt, uint64_t inc)
{
    const uint64_t lane_bits = std::numeric_limits<uint8_t>::digits;
    for (uint64_t j = 0; j < cnt.size(); j++) cnt[j] += (inc >> (lane_bits * j)) & std::numeric_limits<uint8_t>::max();
}

// the last operation of an alignment (see reference_model_costs)
enum last_operation_t : uint64_t { _last_match, _last_insertion, _last_deletion };

/**
 * @brief returns the cheapest costs of aligning P to the substrings of the text under an error model:
 *        [0..n-1][0..len_max] result[s][len] = the cost of the cheapest valid alignment of P to text[s..s+len-1]
 *        (len >= 1), or bound + 1 if there is none or len = 0 (len_max = min(n - s, m + max_errors())); DP over (row,
 *        column, last operation, counter values)
 */
inline std::vector<std::vector<uint64_t>> reference_model_costs(
    std::span<const uint8_t> text, std::span<const uint8_t> P, const apm_costs_t& costs)
{
    const uint64_t n = text.size(), m = P.size();
    // the largest difference of the lengths of P and the aligned substring
    const uint64_t kb = costs.metric() == HAMMING_DISTANCE ? 0 : costs.max_errors();
    const uint64_t none = uint64_t(costs.bound()) + 1;
    std::vector<std::vector<uint64_t>> result(n);

    // an alignment state: the last operation (see last_operation_t) and the values of the counters
    using key_t = std::pair<uint64_t, std::vector<uint64_t>>;
    // the cheapest cost of each alignment state of a DP cell
    using cell_t = std::vector<std::pair<key_t, uint64_t>>;

    for (uint64_t s = 0; s < n; s++) {
        const uint64_t len_max = std::min(n - s, m + kb);
        /* the DP matrix: [0..m][0..len_max] cell(i, j) = the cheapest costs of the alignment states of P[0..i-1]
           aligned to text[s..s+j-1] */
        std::vector<cell_t> cells((m + 1) * (len_max + 1));
        auto cell = [&](uint64_t i, uint64_t j) -> cell_t& { return cells[i * (len_max + 1) + j]; };

        // stores the cost of the state (state, cnt) in the cell c, if it is cheaper than the stored one and valid
        auto relax = [&](cell_t& c, uint64_t state, std::vector<uint64_t> cnt, uint64_t cost) {
            if (cost > costs.bound() || !within_counters(costs, cnt)) return;
            for (auto& [key, val] : c) if (key.first == state && key.second == cnt) { val = std::min(val, cost); return; }
            c.emplace_back(key_t{state, std::move(cnt)}, cost);
        };

        relax(cell(0, 0), _last_match, std::vector<uint64_t>(costs.num_counters(), 0), 0);

        for (uint64_t i = 0; i <= m; i++) {
            for (uint64_t j = 0; j <= len_max; j++) {
                if (j + kb < i || j > i + kb) continue;

                if (i > 0 && j > 0) {
                    uint64_t r = costs.rank(text[s + j - 1]);
                    bool match = costs.matches(i - 1, r);
                    uint64_t c = match ? 0 : costs.sub(i - 1, r);

                    if (c != apm_forbidden) {
                        for (const auto& [key, val] : cell(i - 1, j - 1)) {
                            std::vector<uint64_t> cnt = key.second;
                            if (!match) add_counts(cnt, costs.sub_count(i - 1));
                            relax(cell(i, j), _last_match, cnt, val + c);
                        }
                    }
                }

                if (i > 0 && costs.ins(i - 1) != apm_forbidden) {
                    for (const auto& [key, val] : cell(i - 1, j)) {
                        std::vector<uint64_t> cnt = key.second;
                        bool opens = key.first != _last_insertion;
                        add_counts(cnt, costs.ins_count(i - 1));
                        if (opens) add_counts(cnt, costs.ins_open_count(i - 1));
                        relax(cell(i, j), _last_insertion, cnt, val + costs.ins(i - 1) + (opens ? costs.ins_open(i - 1) : 0));
                    }
                }

                if (j > 0 && costs.del(i) != apm_forbidden) {
                    for (const auto& [key, val] : cell(i, j - 1)) {
                        std::vector<uint64_t> cnt = key.second;
                        bool opens = key.first != _last_deletion;
                        add_counts(cnt, costs.del_count(i));
                        if (opens) add_counts(cnt, costs.del_open_count(i));
                        relax(cell(i, j), _last_deletion, cnt, val + costs.del(i) + (opens ? costs.del_open(i) : 0));
                    }
                }
            }
        }

        result[s].assign(len_max + 1, none);

        for (uint64_t j = 1; j <= len_max; j++) {
            for (const auto& [key, val] : cell(m, j)) result[s][j] = std::min(result[s][j], val);
        }
    }

    return result;
}

// returns the cost of the alignment of P to text[pos..pos+len-1] given by a CIGAR under a model, or bound + 1 if the
// CIGAR is no valid alignment of them
inline uint64_t model_cigar_cost(std::span<const uint8_t> text, std::span<const uint8_t> P, uint64_t pos, uint64_t len,
    const cigar_t<uint8_t>& cigar, const apm_costs_t& costs)
{
    const uint64_t invalid = uint64_t(costs.bound()) + 1;
    uint64_t i = 0, j = 0, cost = 0;
    last_operation_t state = _last_match;
    std::vector<uint64_t> cnt(costs.num_counters(), 0);

    for (const cigar_run_t<uint8_t>& run : cigar) {
        for (uint64_t x = 0; x < run.len; x++) {
            if (run.op == cigar_op_t::MATCH || run.op == cigar_op_t::MISMATCH) {
                if (i >= P.size() || j >= len) return invalid;
                uint8_t t = text[pos + j];
                if (run.op == cigar_op_t::MATCH && t != P[i]) return invalid;
                if (run.op == cigar_op_t::MISMATCH && (t != run.sym || t == P[i])) return invalid;
                uint64_t r = costs.rank(t);

                if (!costs.matches(i, r)) {
                    if (costs.sub(i, r) == apm_forbidden) return invalid;
                    cost += costs.sub(i, r);
                    add_counts(cnt, costs.sub_count(i));
                }

                i++; j++; state = _last_match;
            } else if (run.op == cigar_op_t::INS) {
                if (i >= P.size() || P[i] != run.sym || costs.ins(i) == apm_forbidden) return invalid;
                cost += costs.ins(i) + (state != _last_insertion ? costs.ins_open(i) : 0);
                add_counts(cnt, costs.ins_count(i));
                if (state != _last_insertion) add_counts(cnt, costs.ins_open_count(i));
                i++; state = _last_insertion;
            } else {
                if (j >= len || text[pos + j] != run.sym || costs.del(i) == apm_forbidden) return invalid;
                cost += costs.del(i) + (state != _last_deletion ? costs.del_open(i) : 0);
                add_counts(cnt, costs.del_count(i));
                if (state != _last_deletion) add_counts(cnt, costs.del_open_count(i));
                j++; state = _last_deletion;
            }
        }
    }

    if (i != P.size() || j != len || cost > costs.bound() || !within_counters(costs, cnt)) return invalid;
    return cost;
}

/**
 * @brief returns a random error model compiled for the pattern (see apm_costs_t), or none if its search would allow
 *        more than max_errors errors or have at least as many parts as the pattern has symbols
 * @param quality output: the random quality string of the pattern
 */
template <typename idx_t>
std::optional<apm_costs_t> random_query_model(const idx_t& index, std::span<const uint8_t> input,
    std::span<const uint8_t> pattern, distance_metric_t metric, uint64_t max_errors, std::mt19937_64& gen_query,
    std::string& quality)
{
    std::string syms;
    for (uint8_t c : input) if (syms.find(char(c)) == std::string::npos) syms += char(c);
    apm_model_t model = random_apm_model(syms, pattern.size(), gen_query);
    quality = random_quality(pattern.size(), gen_query);
    uint16_t bound = uint16_t(gen_query() % (max_errors * max_model_cost + 1));
    std::optional<apm_costs_t> costs;
    costs.emplace(model, pattern, quality, gen_query() % 2 == 0, bound, metric, index.forward_index().map_int());
    if (costs->max_errors() > max_errors || min_u_scheme(costs->max_errors()).p >= pattern.size()) return std::nullopt;
    return costs;
}

/* count_hamming_dist, count_hamming_dist_histogram (locate = false) or locate (locate = true) with the Hamming distance
   and an error model: the count, the histogram and the occurrences must equal the reference, and every CIGAR must be an
   alignment with the reported cost */
template <typename idx_t>
void verify_hamming_model(const idx_t& index, std::span<const uint8_t> pattern, std::span<const uint8_t> input,
    const apm_costs_t& costs, apm_params_t params, bool locate, std::mt19937_64& gen_query)
{
    params.costs = &costs;
    const search_scheme_t scheme = min_u_scheme(costs.max_errors());
    const uint64_t m = pattern.size();
    std::vector<uint8_t> P(pattern.begin(), pattern.end());
    // the occurrences of the reference, in increasing order of their positions
    std::vector<aprx_occ_t<>> correct;
    // [0..bound] histogram[c] = the number of occurrences with cost c
    std::vector<uint64_t> histogram(costs.bound() + 1, 0);
    auto ref = reference_model_costs(input, pattern, costs);

    for (uint64_t s = 0; s + m <= input.size(); s++) {
        if (ref[s][m] > costs.bound()) continue;
        correct.emplace_back(aprx_occ_t<>(s, m, ref[s][m]));
        histogram[ref[s][m]]++;
    }

    if (!locate) {
        EXPECT_EQ(index.count_hamming_dist(P, scheme, params), correct.size());
        EXPECT_EQ(index.count_hamming_dist_histogram(P, scheme, params), histogram);
        return;
    }

    if constexpr (idx_t::supports_multiple_locate) {
        if (gen_query() % 2 == 0) {
            auto occ = index.locate(P, scheme, HAMMING_DISTANCE, params);
            ips2ra::sort(occ.begin(), occ.end(), key_of_pos {});
            EXPECT_EQ(occ, correct);
        } else {
            std::vector<aprx_occ_t<CIGAR>> occ_cigar;
            auto cigars = index.template locate<CIGAR>(P, scheme, HAMMING_DISTANCE, occ_cigar, params);
            std::vector<aprx_occ_t<>> occ;

            for (const auto& o : occ_cigar) {
                occ.push_back({o.pos(), o.len(), o.err()});
                EXPECT_EQ(model_cigar_cost(input, pattern, o.pos(), o.len(), cigars[o.cig_idx], costs), o.err());
            }

            ips2ra::sort(occ.begin(), occ.end(), key_of_pos {});
            EXPECT_EQ(occ, correct);
        }
    }
}

// edit-distance locate with an error model: every reported occurrence is a valid alignment with the reported cost (and
// CIGAR), each position is reported at most once, and every valid alignment starting at p has a reported occurrence in
// [p - 2 * max_errors() - 1, p + 2 * max_errors() + 1] whose cost is at most its own
template <typename idx_t>
void verify_edit_model(const idx_t& index, std::span<const uint8_t> pattern, std::span<const uint8_t> input,
    const apm_costs_t& costs, apm_params_t params, std::mt19937_64& gen_query)
{
    params.costs = &costs;
    const search_scheme_t scheme = min_u_scheme(costs.max_errors());
    const uint64_t max_shift = 2 * costs.max_errors() + 1;
    std::vector<uint8_t> P(pattern.begin(), pattern.end());
    auto ref = reference_model_costs(input, pattern, costs);
    std::vector<aprx_occ_t<>> occ;
    std::vector<aprx_occ_t<CIGAR>> occ_cigar;
    std::vector<cigar_t<uint8_t>> cigars;

    if (gen_query() % 2 == 0) {
        cigars = index.template locate<CIGAR>(P, scheme, EDIT_DISTANCE, occ_cigar, params);
        for (const auto& o : occ_cigar) occ.push_back({o.pos(), o.len(), o.err()});
    } else {
        occ = index.locate(P, scheme, EDIT_DISTANCE, params);
    }

    for (const auto& o : occ_cigar) {
        if (o.pos() + o.len() > input.size()) { ADD_FAILURE(); continue; }
        EXPECT_EQ(model_cigar_cost(input, pattern, o.pos(), o.len(), cigars[o.cig_idx], costs), o.err())
            << "pos " << o.pos() << ", len " << o.len();
    }

    ips2ra::sort(occ.begin(), occ.end(), key_of_pos {});

    for (uint64_t x = 0; x < occ.size(); x++) {
        const aprx_occ_t<>& o = occ[x];
        EXPECT_FALSE(x > 0 && occ[x - 1].pos() == o.pos()) << "pos " << o.pos() << " reported twice";
        EXPECT_TRUE(o.pos() < input.size() && o.len() < ref[o.pos()].size() && ref[o.pos()][o.len()] == o.err())
            << "pos " << o.pos() << ", len " << o.len() << ", cost " << o.err() << " is no valid alignment";
    }

    // the reported occurrences in [s - max_shift, s + max_shift] as a window over occ
    uint64_t next = 0;

    for (uint64_t s = 0; s < input.size(); s++) {
        uint64_t best = *std::min_element(ref[s].begin(), ref[s].end());
        if (best > costs.bound()) continue;
        while (next < occ.size() && occ[next].pos() + max_shift < s) next++;
        bool covered = false;

        for (uint64_t x = next; x < occ.size() && occ[x].pos() <= s + max_shift && !covered; x++) {
            covered = occ[x].err() <= best;
        }

        EXPECT_TRUE(covered) << "pos " << s << " (cost " << best << ") has no reported occurrence at most " << max_shift
            << " positions away with at most its cost";
    }
}

/**
 * @brief verifies one approximate-matching algorithm with random patterns drawn from the input
 * @param test the algorithm to verify (not _test_index_operations)
 * @param absent a symbol that does not occur in the input
 * @param thorough false <=> the index was loaded from a stream or built from a cache (fewer queries are performed)
 */
template <typename inp_t, move_r_support support>
void verify_algorithm(const move_rb<support>& index, const inp_t& input,
    move_rb_test_t test, sym_of<inp_t> absent, const std::string& label, bool thorough)
{
    using idx_t = move_rb<support>;
    // true <=> the index and the input type support the approximate-matching algorithms (byte alphabets only)
    static constexpr bool supports_apm = sizeof(sym_of<inp_t>) == 1;
    // the edit distance and the error models need an index with a byte alphabet (see move_r::byte_alphabet)
    if ((is_edit_distance_test(test) || is_model_test(test)) && !index.forward_index().byte_alphabet) return;
    uint64_t n = input.size();
    std::span<const sym_of<inp_t>> input_span(input);

    uint64_t num_queries = thorough
        ? std::clamp<uint64_t>(n / input_symbols_per_query, min_num_queries, max_num_queries)
        : num_queries_not_thorough;

    // [0..num_queries-1] the seed of the random number generator of each query
    std::vector<uint64_t> seeds(num_queries);
    for (uint64_t& seed : seeds) seed = gen();

    #pragma omp parallel for num_threads(max_num_threads), schedule(dynamic,1)
    for (uint64_t cur_query = 0; cur_query < num_queries; cur_query++) {
        std::mt19937_64 gen_query(seeds[cur_query]);
        inp_t pattern = random_query_pattern(input, absent, gen_query);
        std::span<const sym_of<inp_t>> pattern_span(pattern);
        const bool is_edit = is_edit_distance_test(test);
        const uint64_t max_k = max_errors_for(pattern.size(), is_edit);
        uint64_t k = std::uniform_int_distribution<uint64_t>(0, max_k)(gen_query);
        const search_scheme_t scheme = min_u_scheme(k);
        const apm_params_t params = random_apm_params(gen_query);

        if (is_model_test(test)) {
            if constexpr (supports_apm) {
                std::span<const uint8_t> pattern_bytes = byte_span(pattern);
                std::span<const uint8_t> input_bytes = byte_span(input);
                std::string quality;
                auto costs = random_query_model(index, input_bytes, pattern_bytes,
                    is_edit ? EDIT_DISTANCE : HAMMING_DISTANCE, max_k, gen_query, quality);
                if (!costs) continue;

                if (!is_edit) {
                    verify_hamming_model(index, pattern_bytes, input_bytes, *costs, params,
                        test == _test_locate_hamming_model, gen_query);
                } else if constexpr (idx_t::supports_multiple_locate) {
                    verify_edit_model(index, pattern_bytes, input_bytes, *costs, params, gen_query);
                }
            }

            continue;
        }

        switch (test) {
            case _test_count_hamming:
                if constexpr (supports_apm)
                    verify_count_hamming(index, pattern, pattern_span, scheme, k, input_span, params);
                break;
            case _test_locate_hamming:
                if constexpr (supports_apm && idx_t::supports_multiple_locate)
                    verify_locate_hamming(index, pattern, pattern_span, scheme, k, input_span, params);
                break;
            case _test_locate_edit:
                if constexpr (supports_apm && idx_t::supports_multiple_locate)
                    verify_locate_edit<false>(index, pattern, pattern_span, scheme, k, input_span, params);
                break;
            case _test_locate_edit_cigar:
                if constexpr (supports_apm && idx_t::supports_multiple_locate)
                    verify_locate_edit<true>(index, pattern, pattern_span, scheme, k, input_span, params);
                break;
            default: break;
        }
    }
}

// ############################# TEST DRIVERS #############################

/**
 * @brief builds one random move_rb index and verifies one part of it (and of the indexes built from the pfp cache and
 *        loaded from a stream); the input, the alphabet size, the construction mode, the input source (file or
 *        memory), the balancing parameter, the number of threads and the context tables are chosen randomly
 * @tparam inp_t the container type of the input
 * @param test the part to verify
 * @param min_n the smallest length of the input
 * @param max_n the largest length of the input
 */
template <typename inp_t, move_r_support support>
void build_and_verify_move_rb(move_rb_test_t test, uint64_t min_n, uint64_t max_n)
{
    using sym_t = sym_of<inp_t>;
    using idx_t = move_rb<support>;
    static constexpr bool byte_input = sizeof(sym_t) == 1;

    uint64_t input_seed = std::uniform_int_distribution<uint64_t>(1, ~uint64_t(0))(gen);
    inp_t input;
    uint64_t span; // the number of distinct symbols the input may contain
    // a symbol that does not occur in the input (one that occurs, if a byte input may contain all 256 byte values)
    sym_t absent;

    if constexpr (byte_input) {
        const random_input_alphabet alphabet(gen);
        span = alphabet.bytes.size();
        input = alphabet.template input<inp_t>(min_n, max_n, input_seed);
        absent = sym_t(alphabet.absent());
    } else {
        // the number of distinct symbols (log-uniform); the symbol range starts at 1 or at the smallest value of sym_t
        // (negative for a signed type)
        span = random_log_uniform_size(min_random_alphabet_size, max_integer_alphabet_size, gen);
        const sym_t min_sym = gen() % 2 == 0 ? sym_t(1) : std::numeric_limits<sym_t>::min();
        const sym_t max_sym = sym_t(int64_t(min_sym) + int64_t(span) - 1);
        input = random_repetitive_input<inp_t>(min_n, max_n, min_sym, max_sym, input_seed);
        absent = sym_t(int64_t(max_sym) + 1);
    }

    uint64_t n = input.size();

    double context_table_share = context_table_shares[
        std::uniform_int_distribution<uint64_t>(0, context_table_shares.size() - 1)(gen)];

    move_r_params params {
        .mode = _suffix_array,
        .num_threads = std::uniform_int_distribution<uint64_t>(1, max_num_threads)(gen),
        .a = min_balancing_parameter
            + uint64_t(std::min(a_distrib(gen), double(max_balancing_parameter - min_balancing_parameter))),
        .reorder_alphabet = gen() % 2 == 0,
        .revert_start_samples = revert_start_sample_counts[gen() % revert_start_sample_counts.size()],
        .context_table_share = context_table_share,
        .context_tables_single_length = gen() % 2 == 0
    };

    // use the prefix-free parsing for large enough inputs (see pfp_inverse_probability), else one of the two in-memory
    // modes
    const bool use_pfp = n >= min_pfp_input_size && gen() % pfp_inverse_probability == 0;
    if (!use_pfp) params.mode = gen() % 2 == 0 ? _suffix_array : _suffix_array_space;
    else params.mode = _pfp;

    // half of the byte inputs are written to a file and the index is built from the file (the pfp cache
    // is written only for an input file, see move_r::pfp_cache_path)
    std::string input_file;
    if (byte_input && gen() % 2 == 0) {
        input_file = temporary_file_path("move-rb-test-");
        direct_ofstream out(input_file);
        ASSERT_TRUE(out.good()) << "could not write the input to " << input_file;
        out.write(reinterpret_cast<const char*>(input.data()), n * sizeof(sym_t));
        out.close();
        params.input_file = input_file;
    }

    // empty container that passes the symbol type to the constructor if the input is read from a file
    inp_t symbol_type;

    auto make = [&]() {
        return idx_t(input_file.empty() ? input : symbol_type, params, false);
    };

    idx_t index = make();
    naive_reference ref = build_naive_reference<idx_t>(index, input);

    // parameters of the index, appended to the failure messages
    std::string label = std::string("support ") + std::to_string(int(support))
        + ", " + std::to_string(8 * sizeof(sym_t)) + "-bit symbols"
        + ", n = " + std::to_string(n)
        + ", sigma = " + std::to_string(index.alphabet_size() + 1)
        + ", a = " + std::to_string(params.a)
        + ", mode " + (params.mode == _pfp ? "pfp"
            : (params.mode == _suffix_array ? "suffix-array" : "suffix-array-space"))
        + (input_file.empty() ? ", from memory" : ", from a file")
        + ", " + std::to_string(params.num_threads) + " threads"
        + (params.reorder_alphabet ? ", reordered" : ", not reordered")
        + ", revert starts " + std::to_string(params.revert_start_samples)
        + ", context share " + std::to_string(params.context_table_share)
        + (params.context_tables_single_length ? ", single length" : ", chained");


    // verifies the part test of an index built from the input (thorough: see verify_extend)
    auto verify_functionality = [&](const idx_t& idx, bool thorough) {
        if (test == _test_index_operations) verify_extend(idx, input, ref, absent, label, thorough);
        else verify_algorithm(idx, input, test, absent, label, thorough);
    };

    verify_construction(index, input, ref, label);
    verify_functionality(index, true);

    // build a second index from the pfp cache (the RLBWT and the SA samples) written by the first
    // construction and verify it against the references
    if (use_pfp && !input_file.empty()) {
        EXPECT_TRUE(std::filesystem::exists(input_file + ".bwt-runs")) << "the cache was not written";
        EXPECT_TRUE(std::filesystem::exists(input_file + ".bwt-runs-rev"))
            << "the reverse cache was not written";

        idx_t from_cache = make();
        verify_construction(from_cache, input, ref, label + ", read back from the cache");
        verify_functionality(from_cache, false);
    }

    // verify the index after a serialize/load round-trip
    idx_t loaded = serialize_and_load(index);
    verify_construction(loaded, input, ref, label + ", loaded from a stream");
    verify_functionality(loaded, false);

    // remove the input file, its caches and the temporary files of the prefix-free parsing
    if (!input_file.empty()) remove_temporary_input(input_file);
}

/* runs build_and_verify_move_rb in the fuzz iteration i with the container type number
   (i / iterations_per_input_type) mod (number of types) of the input: std::string, std::vector<uint8_t>,
   std::vector<int8_t> (byte alphabets), std::vector<uint32_t>, std::vector<int32_t> (integer alphabets); the
   approximate-matching tests use the three byte types only */
template <move_r_support supp>
void cycle_input_types(uint64_t i, move_rb_test_t test,
    uint64_t min_n = min_input_size, uint64_t max_n = max_input_size)
{
    // the number of input types, and of input types with a byte alphabet
    constexpr uint64_t num_input_types = 5;
    constexpr uint64_t num_byte_input_types = 3;

    switch ((i / iterations_per_input_type) % (is_apm_test(test) ? num_byte_input_types : num_input_types)) {
        case 0: build_and_verify_move_rb<std::string, supp>(test, min_n, max_n); break;
        case 1: build_and_verify_move_rb<std::vector<uint8_t>, supp>(test, min_n, max_n); break;
        case 2: build_and_verify_move_rb<std::vector<int8_t>, supp>(test, min_n, max_n); break;
        case 3: build_and_verify_move_rb<std::vector<uint32_t>, supp>(test, min_n, max_n); break;
        case 4: build_and_verify_move_rb<std::vector<int32_t>, supp>(test, min_n, max_n); break;
    }
}

/* the instantiations of cycle_input_types, declared (prefix extern) or defined (empty prefix);
   tests/test-move-rb-instantiations.cpp is compiled once per support */
#define TEST_MOVE_RB_CYCLE_INPUT_TYPES(prefix, supp) \
    prefix template void cycle_input_types<supp>(uint64_t, move_rb_test_t, uint64_t, uint64_t);
TEST_MOVE_RB_CYCLE_INPUT_TYPES(extern, _locate_move)
TEST_MOVE_RB_CYCLE_INPUT_TYPES(extern, _locate_rlzsa)
TEST_MOVE_RB_CYCLE_INPUT_TYPES(extern, _count)
TEST_MOVE_RB_CYCLE_INPUT_TYPES(extern, _locate_one)

// runs cycle_input_types in the fuzz iteration i with one of the supports that locate all occurrences, which changes
// every iterations_per_support iterations
inline void cycle_multiple_locate_supports(uint64_t i, move_rb_test_t test,
    uint64_t min_n = min_input_size, uint64_t max_n = max_input_size)
{
    if ((i / iterations_per_support) % 2 == 0) cycle_input_types<_locate_move>(i, test, min_n, max_n);
    else cycle_input_types<_locate_rlzsa>(i, test, min_n, max_n);
}

// runs cycle_input_types in the fuzz iteration i with one of the four supports, which changes every
// iterations_per_support iterations
inline void cycle_supports(uint64_t i, move_rb_test_t test,
    uint64_t min_n = min_input_size, uint64_t max_n = max_input_size)
{
    constexpr uint64_t num_supports = 4;

    switch ((i / iterations_per_support) % num_supports) {
        case 0: cycle_input_types<_locate_move>(i, test, min_n, max_n); break;
        case 1: cycle_input_types<_locate_rlzsa>(i, test, min_n, max_n); break;
        case 2: cycle_input_types<_count>(i, test, min_n, max_n); break;
        case 3: cycle_input_types<_locate_one>(i, test, min_n, max_n); break;
    }
}

#endif // MOVE_R_TESTS_TEST_MOVE_RB_COMMON_HPP
