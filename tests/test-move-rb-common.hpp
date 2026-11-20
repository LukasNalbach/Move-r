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
 * @file test-move-rb-common.hpp
 * @brief fuzz tests of move_rb (see test-move-rb.cpp)
 */

#ifndef MOVE_R_TESTS_TEST_MOVE_RB_COMMON_HPP
#define MOVE_R_TESTS_TEST_MOVE_RB_COMMON_HPP

#include <filesystem>
#include <fstream>
#include <map>
#include <numeric>
#include <optional>
#include <set>
#include <span>

#include <gtest/gtest.h>
#include <gtl/phmap.hpp>
#include <move_rb/move_rb.hpp>
#include <algorithms/documents.hpp>
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

// one generator per thread (run_fuzz runs parallel functionalities on all threads)
inline thread_local std::mt19937_64 gen(test_seed() + omp_get_thread_num());
inline uint64_t max_num_threads = omp_get_max_threads();

// returns the number of threads of a construction: 1 in a parallel functionality (run_fuzz runs it on all threads),
// else a random number in [1..max_num_threads]
inline uint64_t random_num_threads()
{
    return omp_in_parallel() ? 1 : std::uniform_int_distribution<uint64_t>(1, max_num_threads)(gen);
}

// the smallest and the largest length of a random input (run_fuzz reduces max_input_size)
inline uint64_t min_input_size = 1;
inline uint64_t max_input_size = 3000;
// the smallest value run_fuzz reduces max_input_size to
inline constexpr uint64_t min_fuzz_input_size = 256;
// the balancing parameter a of a random index is min_balancing_parameter + a value drawn from a_distrib, at most
// max_balancing_parameter
inline thread_local std::lognormal_distribution<double> a_distrib(balancing_parameter_log_mean,
    balancing_parameter_log_deviation);

// a byte input is split into sequences (see split_into_sequences) with probability 1 / sequences_inverse_probability
inline constexpr uint64_t sequences_inverse_probability = 4;
/* the input of an index-operations check has an integer alphabet of at least move_r_params::min_sigma_L_prev_min
   symbols (L_prev_min is built) and large_alphabet_symbol_occurrences times as many symbols, with probability 1 /
   large_alphabet_inverse_probability */
inline constexpr uint64_t large_alphabet_inverse_probability = 16;
inline constexpr uint64_t large_alphabet_symbol_occurrences = 4;
/* the input of an approximate-matching check without an error model has between long_text_min_length and
   long_text_max_length symbols (so that a position in a pattern can exceed 16 bits) with probability 1 /
   long_text_inverse_probability */
inline constexpr uint64_t long_text_inverse_probability = 32;
inline constexpr uint64_t long_text_min_length = uint64_t(1) << std::numeric_limits<uint16_t>::digits;
inline constexpr uint64_t long_text_max_length = 2 * long_text_min_length;
// TEST(test_move_rb, parallel_construction): the inputs have input_symbols_per_thread symbols per thread of the
// construction, which then uses all of its threads (see move_r::construction::limit_num_threads)
inline constexpr uint64_t input_symbols_per_thread = 2 * move_r_min_symbols_per_thread;

// the largest length of a random pattern
inline constexpr uint64_t max_pattern_length = 100;
/* a pattern of an approximate-matching check without an error model is nearly the whole text with probability 1 /
   whole_text_inverse_probability (see random_pattern) if n * min(n, n / r) <= max_whole_text_pattern_steps for a text
   with n symbols and r BWT runs: the search follows each path with an error along a repeat of the text, which has about
   n / r symbols, so it visits about m * min(m, n / r) contexts for m = n */
inline constexpr uint64_t whole_text_inverse_probability = 8;
inline constexpr uint64_t max_whole_text_pattern_steps = uint64_t(1) << 18;
// the largest number of errors k of an approximate-matching query
inline constexpr uint64_t max_query_errors = 10;
/* on a text with more than max_exhaustive_query_input_size symbols, an approximate-matching query has at most one
   error per pattern_symbols_per_query_error symbols of its pattern (with an error model, a cost bound of at most one
   per pattern_symbols_per_query_error symbols); with more errors, nearly every substring of the text is an occurrence,
   and the search visits all of them */
inline constexpr uint64_t pattern_symbols_per_query_error = 8;
inline constexpr uint64_t max_exhaustive_query_input_size = 128;
// a random option of an error model or a search scheme is changed (see randomly_changed_value) with probability
// 1 / changed_option_inverse_probability
inline constexpr uint64_t changed_option_inverse_probability = 8;
// the k of the search scheme of a query with an error model is smaller than the one the model needs (see
// apm_costs_t::max_errors) with probability 1 / too_few_errors_inverse_probability
inline constexpr uint64_t too_few_errors_inverse_probability = 8;

// the number of consecutive fuzz iterations with the same support, and with the same input type (see cycle_profiles
// and cycle_input_types)
inline constexpr uint64_t iterations_per_support = 2;
inline constexpr uint64_t iterations_per_input_type = 4;

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

/* verify_algorithm runs n / input_symbols_per_query queries, at most max_reference_symbols / n (the naive references
   read about n symbols per query), in [min_num_queries, max_num_queries] (not thorough: num_queries_not_thorough) */
inline constexpr uint64_t input_symbols_per_query = 4;
inline constexpr uint64_t max_reference_symbols = uint64_t(1) << 22;
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
    // random documents of the input and the document queries (see document_query)
    _test_documents,
    // count_hamming_dist, count_hamming_dist_histogram, locate and locate<CIGAR> with the Hamming distance
    _test_apm_hamming,
    // locate and locate<CIGAR> with the edit distance
    _test_apm_edit,
    // _test_apm_hamming with a random error model (see random_error_model)
    _test_apm_hamming_model,
    // _test_apm_edit with a random error model
    _test_apm_edit_model
};

// true <=> test is an approximate-matching test (these are run with byte alphabets only)
inline constexpr bool is_apm_test(move_rb_test_t test)
{
    return test == _test_apm_hamming || test == _test_apm_edit || test == _test_apm_hamming_model
        || test == _test_apm_edit_model;
}

// returns the symbols of a container of 1-byte symbols as bytes
template <typename inp_t>
inline std::span<const uint8_t> byte_span(const inp_t& container)
{
    static_assert(sizeof(sym_of<inp_t>) == 1);
    return { reinterpret_cast<const uint8_t*>(container.data()), container.size() };
}

/* returns random apm_params_t (the partitioning, the seeding, the lowering of the bounds and the anchored search must
   not change the result, with the edit distance its coverage of the naive occurrences, see verify_edit_occurrences) */
inline apm_params_t random_apm_params(std::mt19937_64& gen_query)
{
    return {
        .partitioning = partitioning_t(gen_query() % (_partition_optimized + 1)),
        .start_at_part_contexts = gen_query() % 2 == 0,
        .lower_upper_bounds = gen_query() % 2 == 0,
        .anchored_search = gen_query() % 2 == 0
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
    // [0..alphabet_size()] input_symbol[s] = the input symbol of the internal symbol s (see move_r::unmap_symbol)
    std::vector<uint64_t> input_symbol;
    // the internal symbol of each input symbol that occurs in the input
    std::map<uint64_t, uint64_t> internal_symbol;
    // the positions of the separators of the sequences in text, in increasing order (empty: one sequence)
    std::vector<uint64_t> separators;
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
    ref.input_symbol.resize(index.alphabet_size() + 1);
    for (uint64_t s = 0; s < ref.input_symbol.size(); s++) ref.input_symbol[s] = index.unmap_symbol(s);
    for (uint64_t i = 0; i < n; i++) ref.internal_symbol[ref.input_symbol[uint64_t(ref.text[i])]] = uint64_t(ref.text[i]);

    return ref;
}

// ############################# CONSTRUCTION #############################

/* returns the numbers of the distinct substrings without a separator of a text with n symbols for each length L in
   [0..n] (1 for the empty string), from its suffix array SA ([0..n], with the terminator) and the positions of its
   separators: the suffixes whose first L symbols contain no separator and that share fewer than L symbols with the
   suffix before them in SA */
inline std::vector<uint64_t> naive_distinct_substrings(const std::vector<int32_t>& text, const std::vector<int64_t>& SA,
    const std::vector<uint64_t>& separators)
{
    const uint64_t n = text.size();
    // [0..n] the number of positions from j to the next separator or to the end of the text
    std::vector<uint64_t> to_separator(n + 1, 0);
    std::vector<uint8_t> is_separator(n, 0);
    for (uint64_t s : separators) is_separator[s] = 1;
    for (uint64_t j = n; j-- > 0;) to_separator[j] = is_separator[j] ? 0 : to_separator[j + 1] + 1;

    std::vector<uint64_t> rank(n + 1), lcp(n + 1, 0);
    for (uint64_t i = 0; i <= n; i++) rank[SA[i]] = i;

    // the longest common prefix of each suffix with the one before it in SA (Kasai et al.)
    for (uint64_t j = 0, h = 0; j < n; j++) {
        if (rank[j] == 0) {
            h = 0;
            continue;
        }

        const uint64_t k = SA[rank[j] - 1];
        while (j + h < n && k + h < n && text[j + h] == text[k + h]) h++;
        lcp[rank[j]] = h;
        if (h > 0) h--;
    }

    // the suffix SA[i] starts a distinct substring of each length in [lcp[i] + 1, to_separator[SA[i]]]
    std::vector<int64_t> diff(n + 2, 0);

    for (uint64_t i = 0; i <= n; i++) {
        const uint64_t len = to_separator[SA[i]];

        if (lcp[i] < len) {
            diff[lcp[i] + 1]++;
            diff[len + 1]--;
        }
    }

    std::vector<uint64_t> distinct(n + 1);
    distinct[0] = 1;
    int64_t sum = 0;

    for (uint64_t L = 1; L <= n; L++) {
        sum += diff[L];
        distinct[L] = uint64_t(sum);
    }

    return distinct;
}

// verifies the data structures of an index, the reverted text and the suffix array against the naive references
template <typename inp_t, move_r_support support>
void verify_construction(
    const move_rb<support>& index, const inp_t& input, const naive_reference& ref, const std::string& label)
{
    using idx_t = move_rb<support>;
    uint64_t n = input.size();

    // the index of T^R stores no M_Phi, no M_Phi_inv and no rlz encoding of SA^d
    using bwd_t = std::remove_cvref_t<decltype(index.backward_index())>;
    static_assert(!bwd_t::has_mphi);
    static_assert(!bwd_t::has_locate_move);
    static_assert(!bwd_t::has_SA_rlz);

    // a locating index stores the SA samples of the index of T^R
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

    // the numbers of the distinct substrings of each counted length (see move_rb::distinct_substrings)
    EXPECT_TRUE(index.has_distinct_substrings()) << "the numbers of the distinct substrings are missing, " << label;

    if (index.has_distinct_substrings()) {
        const std::vector<uint64_t> distinct = naive_distinct_substrings(ref.text, ref.SA, ref.separators);
        const uint64_t max_len = std::min<uint64_t>(n, std::remove_cvref_t<decltype(index.forward_index())>::max_counted_lcp);
        uint64_t L = 0;
        while (L <= max_len && index.distinct_substrings(L) == distinct[L]) L++;
        EXPECT_GT(L, max_len) << "distinct substrings of length " << L << ": " << index.distinct_substrings(L)
            << " instead of " << distinct[L] << ", " << label;
    }

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

// verifies the context tables: context_of() and table_context_of() must return the contexts that are reached by
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

        // the context of P, computed with context_of() (uses the context tables, see move_rb::context_of)
        bool found = false;
        auto fast = index.template context_of<COUNT>(P, found);
        ASSERT_TRUE(found) << "context_of() missed a substring of the input, " << label;
        ASSERT_EQ(fast.num_occ(), slow.num_occ()) << "context_of() returned the context of another pattern";
        ASSERT_EQ(fast.pattern_length(), P.size());

        // table_context_of() returns the context of the longest stored suffix of P
        uint64_t matched = 0;
        auto entry = index.template table_context_of<COUNT>(P, matched);
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

// ############################# APPROXIMATE PATTERN MATCHING #############################

/**
 * @brief checks the text form of a search scheme: parse_search_scheme returns the scheme; with probability 1 /
 *        changed_option_inverse_probability, also a randomly changed scheme (see randomly_changed_scheme), which is
 *        parsed iff it satisfies the rules of parse_search_scheme, else std::runtime_error is thrown
 * @param scheme the search scheme
 * @param gen_query random number generator
 */
inline void verify_search_scheme_text(const search_scheme_t& scheme, std::mt19937_64& gen_query)
{
    const std::string text = search_scheme_text(scheme);
    EXPECT_TRUE(equal_search_schemes(parse_search_scheme(text), scheme)) << "parse_search_scheme of\n" << text;
    if (gen_query() % changed_option_inverse_probability != 0) return;
    const search_scheme_t changed = randomly_changed_scheme(scheme, gen_query);
    const std::string changed_text = search_scheme_text(changed);

    if (satisfies_scheme_rules(changed)) {
        EXPECT_TRUE(equal_search_schemes(parse_search_scheme(changed_text), changed))
            << "parse_search_scheme of\n" << changed_text;
    } else {
        EXPECT_THROW(parse_search_scheme(changed_text), std::runtime_error) << "parse_search_scheme of\n" << changed_text;
    }
}

/**
 * @brief checks that the rules of an error model are parsed into the model they were written from (see
 *        random_model_t::options and parse_model_options)
 * @param parsed the parsed model
 * @param model the model
 * @param m the pattern length
 */
inline void verify_parsed_model(const apm_model_t& parsed, const random_model_t& model, uint64_t m)
{
    // the number of rules of each kind parsed so far
    uint64_t num_costs = 0, num_counters = 0, num_classes = 0;
    EXPECT_EQ(parsed.iupac, model.iupac);
    EXPECT_EQ(parsed.quality, model.quality);
    EXPECT_EQ(parsed.tau, model.tau);

    for (const model_rule_t& rule : model.rules) {
        const std::string what = rule.option() + " " + rule.text();
        apm_range_t range;

        if (rule.kind == model_rule_t::_cost) {
            ASSERT_LT(num_costs, parsed.costs.size()) << what;
            const apm_cost_rule_t& cost = parsed.costs[num_costs++];
            EXPECT_TRUE(cost.ops == rule.ops && cost.cost == rule.value && cost.from == rule.from && cost.to == rule.to)
                << what;
            range = cost.range;
        } else if (rule.kind == model_rule_t::_counter) {
            ASSERT_LT(num_counters, parsed.counters.size()) << what;
            const apm_counter_t& counter = parsed.counters[num_counters++];
            EXPECT_TRUE(counter.ops == rule.ops && counter.max == rule.value) << what;
            range = counter.range;
        } else {
            ASSERT_LT(num_classes, parsed.classes.size()) << what;
            const apm_class_t& symbol_class = parsed.classes[num_classes++];
            EXPECT_TRUE(symbol_class.from == rule.from && symbol_class.to == rule.to) << what;
            range = symbol_class.range;
        }

        EXPECT_EQ(range.absolute(m), rule.range.positions(m)) << what << ", m = " << m;
    }

    EXPECT_EQ(num_costs, parsed.costs.size());
    EXPECT_EQ(num_counters, parsed.counters.size());
    EXPECT_EQ(num_classes, parsed.classes.size());
}

/**
 * @brief checks the approximate-matching queries of one pattern against the naive references: with the Hamming
 *        distance count_hamming_dist, count_hamming_dist_histogram, locate and locate<CIGAR> (the occurrences and every
 *        CIGAR), with the edit distance locate or locate<CIGAR> (see verify_edit_occurrences, or with an error model:
 *        every occurrence is a valid alignment with the reported cost and lies in one sequence, no position is
 *        reported twice, and every valid alignment that starts at p has an occurrence with at most its cost that starts
 *        at most 2 * max_errors() + 1 positions away); the queries must throw std::invalid_argument if the pattern is
 *        too short for the search scheme (see fits_search_scheme; without an error model only for k > 0, the exact
 *        search does not use the scheme), if the k of the scheme is smaller than the one the error model needs, or if
 *        the index has no byte alphabet and the query uses the edit distance or a model
 *
 * Let e = max_query_errors, or min(max_query_errors, m / pattern_symbols_per_query_error) on a text with more than
 * max_exhaustive_query_input_size symbols. Without an error model, k is drawn from [0, min(m, e)]. With an error model
 * (see random_error_model), its options are parsed (see verify_parsed_model; with probability 1 /
 * changed_option_inverse_probability, one value is changed instead, and parsing must throw std::invalid_argument iff the
 * value is invalid, see is_valid_model_value), compiled for the pattern on a random strand with a random quality and a
 * random cost bound in [0, e * max_model_cost], or in [0, e] on a text with more than max_exhaustive_query_input_size
 * symbols (compiling must throw std::invalid_argument iff a rule has symbols for operations that do not allow them), and
 * k is the largest number of errors of the model (with probability 1 / too_few_errors_inverse_probability less). The
 * search scheme is drawn for k (see random_search_scheme, verify_search_scheme_text).
 *
 * @param index the index
 * @param text the text of the index
 * @param P the pattern
 * @param metric the distance metric
 * @param with_model whether the query has a random error model
 * @param gen_query random number generator
 * @param label the parameters of the index, for the failure messages
 */
template <typename idx_t>
void verify_apm_query(const idx_t& index, const apm_text_t& text, std::span<const uint8_t> P, distance_metric_t metric,
    bool with_model, std::mt19937_64& gen_query, const std::string& label)
{
    auto uniform = [&](uint64_t min, uint64_t max) { return std::uniform_int_distribution<uint64_t>(min, max)(gen_query); };
    const uint64_t m = P.size();
    const std::vector<uint8_t> pattern(P.begin(), P.end());
    apm_params_t params = random_apm_params(gen_query);
    const bool byte_alphabet = index.forward_index().byte_alphabet;
    std::optional<apm_costs_t> costs;
    naive_model_costs model_costs;
    // the largest number of errors of the query (see pattern_symbols_per_query_error)
    const bool exhaustive_allowed = text.T.size() <= max_exhaustive_query_input_size;
    const uint64_t max_k = exhaustive_allowed ? max_query_errors : std::min(max_query_errors, m / pattern_symbols_per_query_error);
    // the k of the search scheme, and the smallest one the query accepts
    uint64_t k = uniform(0, std::min(m, max_k));
    uint64_t required_k = 0;

    if (with_model) {
        // the symbols of the symbol sets of the rules: the symbols of the text and of P, without the syntax of the rules
        std::string symbols;
        for (uint8_t c : distinct_symbols(text.T)) symbols += char(c);
        for (uint8_t c : distinct_symbols(P)) if (symbols.find(char(c)) == std::string::npos) symbols += char(c);
        std::erase_if(symbols, [](char c) { return model_syntax_symbols.find(c) != std::string_view::npos; });
        const random_model_t model = random_error_model(symbols, m, gen_query);
        std::vector<std::pair<std::string, std::string>> options = model.options();

        if (gen_query() % changed_option_inverse_probability == 0) {
            auto& [option, value] = options[uniform(0, options.size() - 1)];
            value = randomly_changed_value(value, gen_query);

            if (is_valid_model_value(option, value)) {
                EXPECT_NO_THROW(parse_model_options(options)) << option << " " << value;
            } else {
                EXPECT_THROW(parse_model_options(options), std::invalid_argument) << option << " " << value;
            }

            return;
        }

        apm_model_t parsed = parse_model_options(options);
        ASSERT_NO_FATAL_FAILURE(verify_parsed_model(parsed, model, m));
        parsed.cost_dp = apm_cost_dp_t(uniform(0, APM_COST_DP_PRUNE));
        const std::string quality = random_quality(m, gen_query);
        const bool reverse_strand = gen_query() % 2 == 0;
        const uint64_t bound = uniform(0, exhaustive_allowed ? max_k * max_model_cost : max_k);

        if (!byte_alphabet) {
            EXPECT_THROW(apm_costs_t(parsed, pattern, quality, reverse_strand, bound, metric,
                index.forward_index().map_int()), std::invalid_argument) << label;
            return;
        }

        if (model.has_disallowed_symbols()) {
            EXPECT_THROW(apm_costs_t(parsed, pattern, quality, reverse_strand, bound, metric,
                index.forward_index().map_int()), std::invalid_argument) << label;
            return;
        }

        costs.emplace(parsed, pattern, quality, reverse_strand, bound, metric, index.forward_index().map_int());
        // the searches of a model with more errors take too long
        if (costs->max_errors() > max_query_errors) return;
        params.costs = &*costs;
        model_costs = naive_model(model, P, quality, reverse_strand, bound, metric);
        required_k = costs->max_errors();
        k = required_k > 0 && gen_query() % too_few_errors_inverse_probability == 0 ? uniform(0, required_k - 1)
            : required_k;
    }

    const search_scheme_t scheme = random_search_scheme(k, gen_query);
    ASSERT_NO_FATAL_FAILURE(verify_search_scheme_text(scheme, gen_query));
    const bool throws = (!fits_search_scheme(m, scheme) && (k > 0 || with_model)) || k < required_k
        || (metric == EDIT_DISTANCE && !byte_alphabet);
    const std::string what = std::string(metric == HAMMING_DISTANCE ? "Hamming" : "edit") + " distance, m = "
        + std::to_string(m) + ", k = " + std::to_string(k) + ", p = " + std::to_string(scheme.p)
        + (with_model ? ", bound " + std::to_string(model_costs.bound) : "") + ", " + label;

    // checks that the CIGAR of each occurrence in one sequence is an alignment with its cost (reports the first that is not)
    auto verify_cigars = [&](const std::vector<aprx_occ_t<CIGAR>>& occ, const std::vector<cigar_t<uint8_t>>& cigars) {
        for (const auto& o : occ) {
            if (!text.in_one_sequence(o.pos(), o.len())) continue; // reported with the occurrences
            const uint64_t cost = with_model
                ? model_cigar_cost(text, P, o.pos(), o.len(), cigars[o.cig_idx], model_costs)
                : num_cigar_edits(P.data(), m, text.T.data() + o.pos(), o.len(), cigars[o.cig_idx]);

            if (cost != o.err()) {
                ADD_FAILURE() << "the CIGAR of the occurrence at " << o.pos() << " of length " << o.len() << " has cost "
                    << cost << ", not " << o.err() << ", " << what;
                return;
            }
        }
    };

    // checks that query() throws std::invalid_argument iff throws, and returns whether it does not throw
    auto runs = [&](auto query) {
        if (!throws) return true;
        EXPECT_THROW(query(), std::invalid_argument) << what;
        return false;
    };

    if (metric == HAMMING_DISTANCE) {
        // the occurrences of the reference, in increasing order of their positions
        std::vector<aprx_occ_t<>> correct;

        if (!with_model) {
            correct = reference_hamming_occurrences(text, P, k);
        } else if (!throws) {
            const auto reference = reference_model_costs(text, P, model_costs, metric);

            for (uint64_t s = 0; s + m <= text.T.size(); s++) {
                if (m < reference[s].size() && reference[s][m] <= model_costs.bound) correct.emplace_back(s, m, reference[s][m]);
            }
        }

        // [0..k] (with a model [0..bound]) histogram[e] = the number of occurrences with e errors
        std::vector<uint64_t> histogram((with_model ? model_costs.bound : k) + 1, 0);
        for (const aprx_occ_t<>& o : correct) histogram[o.err()]++;

        if (runs([&] { index.count_hamming_dist(pattern, scheme, params); }))
            EXPECT_EQ(index.count_hamming_dist(pattern, scheme, params), correct.size()) << what;
        if (runs([&] { index.count_hamming_dist_histogram(pattern, scheme, params); }))
            EXPECT_EQ(index.count_hamming_dist_histogram(pattern, scheme, params), histogram) << what;

        if constexpr (idx_t::supports_multiple_locate) {
            if (runs([&] { index.locate(pattern, scheme, HAMMING_DISTANCE, params); })) {
                auto occ = index.locate(pattern, scheme, HAMMING_DISTANCE, params);
                std::sort(occ.begin(), occ.end());
                EXPECT_EQ(occ, correct) << what;
            }

            std::vector<aprx_occ_t<CIGAR>> occ_cigar;

            if (runs([&] { index.template locate<CIGAR>(pattern, scheme, HAMMING_DISTANCE, occ_cigar, params); })) {
                const auto cigars = index.template locate<CIGAR>(pattern, scheme, HAMMING_DISTANCE, occ_cigar, params);
                verify_cigars(occ_cigar, cigars);
                std::vector<aprx_occ_t<>> occ;
                for (const auto& o : occ_cigar) occ.emplace_back(o.pos(), o.len(), o.err());
                std::sort(occ.begin(), occ.end());
                EXPECT_EQ(occ, correct) << "with CIGARs, " << what;
            }
        }

        return;
    }

    if constexpr (idx_t::supports_multiple_locate) {
        std::vector<aprx_occ_t<>> occ;
        std::vector<aprx_occ_t<CIGAR>> occ_cigar;
        std::vector<cigar_t<uint8_t>> cigars;
        const bool with_cigar = gen_query() % 2 == 0;

        if (with_cigar) {
            if (!runs([&] { index.template locate<CIGAR>(pattern, scheme, EDIT_DISTANCE, occ_cigar, params); })) return;
            cigars = index.template locate<CIGAR>(pattern, scheme, EDIT_DISTANCE, occ_cigar, params);
            for (const auto& o : occ_cigar) occ.emplace_back(o.pos(), o.len(), o.err());
        } else {
            if (!runs([&] { index.locate(pattern, scheme, EDIT_DISTANCE, params); })) return;
            occ = index.locate(pattern, scheme, EDIT_DISTANCE, params);
        }

        verify_cigars(occ_cigar, cigars);
        std::sort(occ.begin(), occ.end());

        if (!with_model) {
            verify_edit_occurrences(occ, reference_edit_occurrences(text, P, k), P, text, k);
            return;
        }

        const auto reference = reference_model_costs(text, P, model_costs, metric);
        const uint64_t max_shift = 2 * costs->max_errors() + 1;

        for (uint64_t x = 0; x < occ.size(); x++) {
            const aprx_occ_t<>& o = occ[x];
            EXPECT_FALSE(x > 0 && occ[x - 1].pos() == o.pos()) << "position " << o.pos() << " is reported twice, " << what;
            EXPECT_TRUE(o.pos() < text.T.size() && o.len() < reference[o.pos()].size()
                && reference[o.pos()][o.len()] == o.err())
                << "the occurrence at " << o.pos() << " of length " << o.len() << " with cost " << o.err()
                << " is no valid alignment, " << what;
        }

        // the reported occurrences in [s - max_shift, s + max_shift] as a window over occ
        uint64_t next = 0;

        for (uint64_t s = 0; s < text.T.size(); s++) {
            const uint64_t best = *std::min_element(reference[s].begin(), reference[s].end());
            if (best > model_costs.bound) continue;
            while (next < occ.size() && occ[next].pos() + max_shift < s) next++;
            bool covered = false;
            for (uint64_t x = next; x < occ.size() && occ[x].pos() <= s + max_shift && !covered; x++) covered = occ[x].err() <= best;
            EXPECT_TRUE(covered) << "the alignment at " << s << " with cost " << best << " has no occurrence at most "
                << max_shift << " positions away with at most its cost, " << what;
        }
    }
}

/**
 * @brief verifies one approximate-matching algorithm with random patterns drawn from the input (see random_pattern)
 * @param test the algorithm to verify (an approximate-matching test)
 * @param symbols the symbols the edits of the patterns draw from (see edit_symbols)
 * @param sequence_starts the starting positions of the sequences of the input (empty: the input is one sequence)
 * @param thorough false <=> the index was loaded from a stream or built from a cache (fewer queries are performed)
 */
template <typename inp_t, move_r_support support>
void verify_algorithm(const move_rb<support>& index, const inp_t& input, move_rb_test_t test,
    const std::vector<sym_of<inp_t>>& symbols, const std::vector<uint64_t>& sequence_starts, const std::string& label,
    bool thorough)
{
    uint64_t n = input.size();

    uint64_t num_queries = thorough
        ? std::clamp<uint64_t>(std::min(n / input_symbols_per_query, max_reference_symbols / n), min_num_queries,
            max_num_queries)
        : num_queries_not_thorough;

    // [0..num_queries-1] the seed of the random number generator of each query
    std::vector<uint64_t> seeds(num_queries);
    for (uint64_t& seed : seeds) seed = gen();

    #pragma omp parallel for num_threads(max_num_threads), schedule(dynamic,1)
    for (uint64_t cur_query = 0; cur_query < num_queries; cur_query++) {
        std::mt19937_64 gen_query(seeds[cur_query]);

        if constexpr (sizeof(sym_of<inp_t>) == 1) {
            if (is_apm_test(test)) {
                const apm_text_t text(byte_span(input), !sequence_starts.empty(), uint8_t(fasta_separator_symbol));
                const bool with_model = test == _test_apm_hamming_model || test == _test_apm_edit_model;
                const bool whole_text_allowed = !with_model
                    && n * std::min(n, n / index.num_bwt_runs()) <= max_whole_text_pattern_steps;
                inp_t pattern = random_pattern(input, max_pattern_length, whole_text_allowed ? whole_text_inverse_probability
                    : 0, symbols, sequence_starts, gen_query);
                // a pattern of a text with sequences does not contain their separator (like a read of a FASTA file)
                if (!sequence_starts.empty()) std::erase(pattern, sym_of<inp_t>(fasta_separator_symbol));
                if (pattern.empty()) continue;
                verify_apm_query(index, text, byte_span(pattern), test == _test_apm_hamming || test == _test_apm_hamming_model
                    ? HAMMING_DISTANCE : EDIT_DISTANCE, with_model, gen_query, label);
            }
        }
    }
}

/**
 * @brief verifies the documents of an index and the document queries against a naive search
 * @param map the documents the index was built with
 * @param symbols the symbols the edits of the patterns draw from (see edit_symbols)
 */
template <typename inp_t, move_r_support support>
void verify_documents(const move_rb<support>& index, const inp_t& input, const document_map& map,
    const std::vector<sym_of<inp_t>>& symbols, const std::string& label) requires(move_rb<support>::supports_multiple_locate)
{
    const uint64_t n = input.size();
    ASSERT_TRUE(index.has_documents()) << "no documents, " << label;
    const document_array& docs = index.documents();
    const document_map& stored = docs.documents();
    ASSERT_EQ(docs.num_documents(), map.num_documents()) << label;
    ASSERT_EQ(stored.num_ranges(), map.num_ranges()) << label;

    for (uint64_t k = 0; k < map.num_ranges(); k++) {
        ASSERT_EQ(stored.range_start(k), map.range_start(k)) << "start of range " << k << ", " << label;
        ASSERT_EQ(stored.range_document(k), map.range_document(k)) << "document of range " << k << ", " << label;
    }

    for (uint64_t d = 0; d < map.num_documents(); d++)
        ASSERT_EQ(stored.name(d), map.name(d)) << "name of document " << d << ", " << label;

    std::mt19937_64 gen_query(gen());
    document_query query(docs);

    for (uint64_t q = 0; q < num_document_queries; q++) {
        const inp_t pattern = random_pattern(input, max_pattern_length, 0, symbols, {}, gen_query);
        const uint64_t m = pattern.size();

        // the number of occurrences of the pattern in each document that contains it
        std::map<uint64_t, uint64_t> frequency;

        for (uint64_t p = 0; p + m <= n; p++) {
            if (std::equal(pattern.begin(), pattern.end(), input.begin() + p)) frequency[map.document_of(p)]++;
        }

        auto ctx = index.forward_index().query();
        bool occurs = true;
        for (uint64_t i = m; occurs && i-- > 0;) occurs = ctx.prepend(pattern[i]);
        ASSERT_EQ(occurs, !frequency.empty()) << "pattern occurrence, " << label;
        if (!occurs) continue;

        // the occurrences of the pattern; locate reports the occurrences not reported yet, so each pass locates with a
        // copy of the context
        auto source = [&](auto&& report) {
            auto pass = ctx;
            pass.locate([&](uint64_t, uint64_t occ) { report(occ); });
        };

        // the documents that contain the pattern, in increasing order
        std::vector<uint64_t> expected;
        for (const auto& [document, freq] : frequency) expected.push_back(document);
        std::vector<uint64_t> listed;
        query.list(source, [&](uint64_t document) { listed.push_back(document); });
        EXPECT_EQ(listed, expected) << "list, " << label;
        EXPECT_EQ(query.count(source), expected.size()) << "count, " << label;

        std::map<uint64_t, uint64_t> listed_frequency;
        query.list_with_frequencies(source, [&](uint64_t document, uint64_t freq) {
            EXPECT_TRUE(listed_frequency.empty() || listed_frequency.rbegin()->first < document)
                << "list_with_frequencies: not in increasing order, " << label;
            listed_frequency[document] = freq;
        });
        EXPECT_EQ(listed_frequency, frequency) << "list_with_frequencies, " << label;

        // the top k by decreasing frequency, ties by the smaller document
        std::vector<std::pair<uint64_t, uint64_t>> ranked;
        for (const auto& [document, freq] : frequency) ranked.emplace_back(freq, document);
        std::sort(ranked.begin(), ranked.end(), [](const auto& x, const auto& y) {
            return x.first != y.first ? x.first > y.first : x.second < y.second; });
        const uint64_t k = 1 + gen_query() % (frequency.size() + 1);
        std::vector<std::pair<uint64_t, uint64_t>> top;
        query.top_k(source, k, [&](uint64_t document, uint64_t freq) { top.emplace_back(freq, document); });
        ranked.resize(std::min<uint64_t>(k, ranked.size()));
        EXPECT_EQ(top, ranked) << "top_k, " << label;
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

    inp_t input;
    // the starting positions of the sequences of a byte input (empty: the input is one sequence)
    std::vector<uint64_t> sequence_starts;
    const bool with_sequences = byte_input && gen() % sequences_inverse_probability == 0;

    if constexpr (byte_input) {
        random_input_alphabet alphabet(gen);
        if (with_sequences) alphabet.remove(uint8_t(fasta_separator_symbol));

        if ((test == _test_apm_hamming || test == _test_apm_edit) && gen() % long_text_inverse_probability == 0) {
            min_n = long_text_min_length;
            max_n = long_text_max_length;
        }

        input = alphabet.template input<inp_t>(min_n, max_n, gen);
        if (with_sequences) sequence_starts = split_into_sequences(input, gen);
    } else {
        uint64_t min_sigma = 1;
        uint64_t max_sigma = max_n;

        if (test == _test_index_operations && gen() % large_alphabet_inverse_probability == 0) {
            min_sigma = max_sigma = move_r_params {}.min_sigma_L_prev_min;
            min_n = max_n = large_alphabet_symbol_occurrences * min_sigma;
        }

        const integer_alphabet<sym_t> alphabet = random_integer_alphabet<sym_t>(min_sigma, max_sigma, gen);
        input = random_text<inp_t>(min_n, max_n, alphabet.sigma, alphabet, gen);
    }

    uint64_t n = input.size();
    // the symbols the edits of the patterns draw from; the last one does not occur in the input (see edit_symbols)
    const std::vector<sym_t> symbols = edit_symbols(input);
    // a symbol that does not occur in the input (one that occurs, if a byte input contains all 256 byte values)
    const sym_t absent = symbols.back();
    // the sequences of the index
    fasta_sequence_data sequence_data;
    std::vector<std::string> sequence_names(sequence_starts.size());
    if (with_sequences) sequence_data.set_sequences(n + 1, sequence_starts, sequence_names, uint8_t(fasta_separator_symbol));

    double context_table_share = context_table_shares[
        std::uniform_int_distribution<uint64_t>(0, context_table_shares.size() - 1)(gen)];

    move_r_params params {
        .mode = _suffix_array,
        .num_threads = random_num_threads(),
        .a = min_balancing_parameter
            + uint64_t(std::min(a_distrib(gen), double(max_balancing_parameter - min_balancing_parameter))),
        .reorder_alphabet = gen() % 2 == 0,
        .revert_start_samples = revert_start_sample_counts[gen() % revert_start_sample_counts.size()],
        .context_table_share = context_table_share,
        .context_tables_single_length = gen() % 2 == 0
    };

    // the documents of the input (_test_documents only)
    document_map documents;

    if (test == _test_documents) {
        documents = random_documents(n + 1, max_document_ranges, gen);
        params.documents = &documents;
    }

    /* the sequences passed to the construction, which counts the distinct substrings without a separator (see
       move_r_params::counts_distinct_substrings): a random part of the separators between the sequences is passed as
       the separators of masked runs */
    fasta_sequence_data construction_sequences;
    no_init_vector<uint64_t> masked_separators;

    if (with_sequences) {
        std::vector<uint64_t> starts { sequence_starts[0] };

        for (uint64_t j = 1; j < sequence_starts.size(); j++) {
            if (gen() % 2 == 0) masked_separators.emplace_back(sequence_starts[j] - 1);
            else starts.emplace_back(sequence_starts[j]);
        }

        construction_sequences.set_sequences(n + 1, starts, std::vector<std::string>(starts.size()),
            uint8_t(fasta_separator_symbol));
        params.fasta.sequences = &construction_sequences;
        params.fasta.masked_separators = &masked_separators;
    }

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

    // empty container with which the symbol type is passed to the constructor if the input is read from a file
    inp_t symbol_type;

    auto make = [&]() {
        idx_t built(input_file.empty() ? input : symbol_type, params, false);
        if (with_sequences) built.set_fasta_sequence_data(sequence_data);
        return built;
    };

    idx_t index = make();
    naive_reference ref = build_naive_reference<idx_t>(index, input);

    if (with_sequences) {
        for (uint64_t s : sequence_starts) if (s > 0) ref.separators.emplace_back(s - 1);
    }

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
        + (params.context_tables_single_length ? ", single length" : ", chained")
        + (with_sequences ? ", " + std::to_string(sequence_starts.size()) + " sequences" : "")
        + (test == _test_documents ? ", " + std::to_string(documents.num_documents()) + " documents in "
            + std::to_string(documents.num_ranges()) + " ranges" : "");

    // verifies the part test of an index built from the input (thorough: see verify_extend)
    auto verify_functionality = [&](const idx_t& idx, bool thorough) {
        if (test == _test_index_operations) {
            verify_extend(idx, input, ref, absent, label, thorough);
        } else if (test == _test_documents) {
            if constexpr (idx_t::supports_multiple_locate) verify_documents(idx, input, documents, symbols, label);
        } else {
            verify_algorithm(idx, input, test, symbols, sequence_starts, label, thorough);
        }
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

// the supports of the indexes of a test (see cycle_profiles): the four supports, or the two that locate all occurrences
enum support_set_t { _all_supports, _multiple_locate_supports };

/**
 * @brief runs cycle_input_types in the fuzz iteration i with one of the supports of a set, which changes every
 *        iterations_per_support iterations
 * @param i the iteration
 * @param test the part to verify
 * @param supports the supports
 * @param min_n the smallest length of the input
 * @param max_n the largest length of the input
 */
inline void cycle_profiles(uint64_t i, move_rb_test_t test, support_set_t supports, uint64_t min_n = min_input_size,
    uint64_t max_n = max_input_size)
{
    // the number of supports of the switch below, of which the first num_multiple_locate_supports locate all occurrences
    constexpr uint64_t num_supports = 4;
    constexpr uint64_t num_multiple_locate_supports = 2;
    const uint64_t num = supports == _all_supports ? num_supports : num_multiple_locate_supports;

    switch ((i / iterations_per_support) % num) {
        case 0: cycle_input_types<_locate_move>(i, test, min_n, max_n); break;
        case 1: cycle_input_types<_locate_rlzsa>(i, test, min_n, max_n); break;
        case 2: cycle_input_types<_count>(i, test, min_n, max_n); break;
        case 3: cycle_input_types<_locate_one>(i, test, min_n, max_n); break;
    }
}

#endif // MOVE_R_TESTS_TEST_MOVE_RB_COMMON_HPP
