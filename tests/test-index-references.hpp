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
 * @file test-index-references.hpp
 * @brief naive references of the index tests: suffix array intervals, random documents and the approximate occurrences
 *        of a pattern without and with an error model
 */

#ifndef MOVE_R_TESTS_TEST_INDEX_REFERENCES_HPP
#define MOVE_R_TESTS_TEST_INDEX_REFERENCES_HPP

#include <algorithm>
#include <array>
#include <cstdint>
#include <random>
#include <set>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include <gtl/phmap.hpp>

#include <algorithms/move_rb/apm/model.hpp>
#include <data_structures/document_map.hpp>
#include <misc/apm.hpp>
#include <misc/fasta.hpp>
#include <misc/hash.hpp>
#include "test-patterns.hpp"

/**
 * @brief returns the suffix array interval [b, e] of the pattern P in text; b > e <=> P does not occur in text
 * @tparam sym_t symbol type of the text and of the pattern
 * @param SA [0..n] the suffix array of text with the sentinel (n = |text|)
 * @param text the text, in the alphabet the suffix array is sorted by
 * @param P the pattern
 */
template <typename sym_t>
inline std::pair<int64_t, int64_t> reference_sa_interval(
    const std::vector<int64_t>& SA, const std::vector<sym_t>& text, const std::vector<sym_t>& P)
{
    const int64_t n = text.size();

    // the suffix of text starting at sa (terminated by the sentinel) is lexicographically smaller than P
    auto suffix_less = [&](int64_t sa, const std::vector<sym_t>& P) {
        for (uint64_t i = 0; i < P.size(); i++) {
            if (sa + (int64_t) i >= n) return true; // the sentinel is smaller than every symbol of P
            if (text[sa + i] != P[i]) return text[sa + i] < P[i];
        }

        return false;
    };

    // P is lexicographically smaller than the suffix of text starting at sa, and not a prefix of it
    auto pattern_less = [&](const std::vector<sym_t>& P, int64_t sa) {
        for (uint64_t i = 0; i < P.size(); i++) {
            if (sa + (int64_t) i >= n) return false;
            if (text[sa + i] != P[i]) return P[i] < text[sa + i];
        }

        return false;
    };

    return {
        std::lower_bound(SA.begin(), SA.end(), P, suffix_less) - SA.begin(),
        std::upper_bound(SA.begin(), SA.end(), P, pattern_less) - SA.begin() - 1
    };
}

/**
 * @brief returns random documents of a text of length n (including the sentinel): the positions are split into
 *        ranges of random lengths, and each range is assigned to one of at most as many documents as there are ranges;
 *        document d is named "d<d>"
 * @param n the number of positions (>= 1)
 * @param max_num_ranges the maximum number of ranges (>= 1)
 * @param gen random number generator
 */
template <typename gen_t>
inline document_map random_documents(uint64_t n, uint64_t max_num_ranges, gen_t& gen)
{
    const uint64_t num_ranges = std::uniform_int_distribution<uint64_t>(1, std::min(n, max_num_ranges))(gen);

    // the starting positions of the ranges, starting with 0
    std::set<uint64_t> start_set { 0 };
    while (start_set.size() < num_ranges) start_set.insert(std::uniform_int_distribution<uint64_t>(1, n - 1)(gen));
    std::vector<uint64_t> starts(start_set.begin(), start_set.end());

    // [0..num_ranges-1] documents[x] = document of the x-th range
    const uint64_t num_documents = std::uniform_int_distribution<uint64_t>(1, num_ranges)(gen);
    std::vector<uint32_t> documents(num_ranges);
    for (uint32_t& document : documents) document = uint32_t(gen() % num_documents);

    std::vector<std::string> names;
    for (uint64_t d = 0; d < num_documents; d++) names.push_back("d" + std::to_string(d));
    return document_map(std::move(starts), std::move(documents), std::move(names));
}

// ############################# APPROXIMATE PATTERN MATCHING #############################

/**
 * @brief the text of an approximate-matching query: an occurrence lies in one sequence of the text, i.e., it contains
 *        no separator of the sequences (if the text has sequences, see fasta_sequence_data)
 */
struct apm_text_t {
    std::span<const uint8_t> T; // the text
    // [0..n] separators_before[i] = the number of separators in T[0..i-1] (all 0 if the text has no sequences)
    std::vector<uint64_t> separators_before;

    /**
     * @param T the text
     * @param has_sequences whether the text has sequences
     * @param separator the separator of the sequences
     */
    apm_text_t(std::span<const uint8_t> T, bool has_sequences, uint8_t separator) : T(T), separators_before(T.size() + 1, 0)
    {
        for (uint64_t i = 0; i < T.size(); i++)
            separators_before[i + 1] = separators_before[i] + (has_sequences && T[i] == separator ? 1 : 0);
    }

    // returns whether T[pos..pos+len-1] lies in one sequence
    bool in_one_sequence(uint64_t pos, uint64_t len) const
    {
        return pos + len <= T.size() && separators_before[pos + len] == separators_before[pos];
    }

    // returns the largest len such that T[pos..pos+len-1] lies in one sequence
    uint64_t sequence_rest(uint64_t pos) const
    {
        uint64_t end = pos;
        while (end < T.size() && separators_before[end + 1] == separators_before[pos]) end++;
        return end - pos;
    }
};

/**
 * @brief returns the occurrences of P with at most k mismatches (Hamming distance) that lie in one sequence of the
 *        text, in increasing order of their positions
 */
inline std::vector<aprx_occ_t<>> reference_hamming_occurrences(const apm_text_t& text, std::span<const uint8_t> P,
    uint64_t k)
{
    std::vector<aprx_occ_t<>> occurrences;
    const uint64_t m = P.size();

    for (uint64_t pos = 0; pos + m <= text.T.size(); pos++) {
        if (!text.in_one_sequence(pos, m)) continue;
        uint64_t mismatches = 0;
        for (uint64_t i = 0; i < m && mismatches <= k; i++) mismatches += text.T[pos + i] != P[i];
        if (mismatches <= k) occurrences.emplace_back(pos, m, mismatches);
    }

    return occurrences;
}

/**
 * @brief returns the naive occurrences of P with at most k errors (edit distance) in each sequence of the text (see
 *        locate_naive), in increasing order of their positions
 */
inline std::vector<aprx_occ_t<>> reference_edit_occurrences(const apm_text_t& text, std::span<const uint8_t> P,
    uint64_t k)
{
    std::vector<aprx_occ_t<>> occurrences;

    for (uint64_t beg = 0; beg <= text.T.size();) {
        const uint64_t len = text.sequence_rest(beg);

        for (const aprx_occ_t<>& o : locate_naive(text.T.subspan(beg, len), P, k, EDIT_DISTANCE))
            occurrences.emplace_back(beg + o.pos(), o.len(), o.err());

        beg += len + 1;
    }

    return occurrences;
}

/**
 * @brief checks the occurrences of a pattern located with the edit distance: no position is reported twice, every
 *        occurrence lies in one sequence of the text and has the reported edit distance (at most k), and the
 *        occurrences cover the naive occurrences before and after filtering them (see verify_edit_distance_coverage
 *        and filter_edit_distance_occurrences)
 * @param occ the located occurrences, sorted by their positions
 * @param correct the naive occurrences (see reference_edit_occurrences)
 * @param pattern the pattern
 * @param text the text
 * @param k the maximum number of errors
 */
inline void verify_edit_occurrences(std::vector<aprx_occ_t<>> occ, const std::vector<aprx_occ_t<>>& correct,
    std::span<const uint8_t> pattern, const apm_text_t& text, uint64_t k)
{
    for (uint64_t x = 1; x < occ.size(); x++) {
        EXPECT_NE(occ[x - 1].pos(), occ[x].pos()) << "position " << occ[x].pos() << " is reported twice";
    }

    // the edit distance of the pattern to each distinct matched substring
    gtl::flat_hash_map<std::span<const uint8_t>, uint64_t, span_hash, span_eq> dist_cache;

    for (const auto& o : occ) {
        if (!text.in_one_sequence(o.pos(), o.len())) {
            ADD_FAILURE() << "the occurrence at " << o.pos() << " of length " << o.len() << " is not in one sequence";
            continue;
        }

        std::span<const uint8_t> match = text.T.subspan(o.pos(), o.len());
        auto [it, inserted] = dist_cache.try_emplace(match, 0);
        if (inserted) it->second = edit_dist_bounded(pattern, match, k);
        EXPECT_TRUE(it->second == o.err() && o.err() <= k) << "the occurrence at " << o.pos() << " of length "
            << o.len() << " has " << it->second << " errors, not " << o.err();
    }

    EXPECT_TRUE(verify_edit_distance_coverage(occ, correct, k));
    filter_edit_distance_occurrences(occ, k);
    EXPECT_TRUE(verify_edit_distance_coverage(occ, correct, k));
}

/**
 * @brief an error model compiled naively for a pattern P from its rules (see random_model_t and apm_model_t): the
 *        reference of apm_costs_t; a text symbol is a byte
 */
struct naive_model_costs {
    // the operations an unbounded counter counts: [0..m-1] at the positions of P, [0..m] at the boundaries
    struct counter_t {
        uint64_t bound = 0;
        std::vector<bool> sub, ins, ins_open; // [0..m-1]
        std::vector<bool> del, del_open; // [0..m]
    };

    uint64_t m = 0; // m = |P|
    uint64_t bound = 0; // the largest cost of a valid alignment
    // [0..m-1][0..255] matches[i][c] <=> P[i] matches the text symbol c
    std::vector<std::array<bool, num_byte_values>> matches;
    // [0..m-1][0..255] sub[i][c] = the cost of aligning P[i] to the text symbol c (0 if it matches)
    std::vector<std::array<uint16_t, num_byte_values>> sub;
    std::vector<uint16_t> ins, ins_open; // [0..m-1] the costs of inserting P[i] and of opening a gap with it
    // [0..m] the costs of a deletion at the boundary b between P[b - 1] and P[b] and of opening a gap with it
    std::vector<uint16_t> del, del_open;
    std::vector<counter_t> counters; // the counters with a bound of at least 1
};

/**
 * @brief compiles an error model naively for a pattern: by default, a symbol matches itself (and its bases with
 *        -iupac), a substitution costs 1 (with -qual the Phred quality of its position), an insertion and a deletion
 *        cost 1 and opening a gap costs 0; then the classes add matches, the cost rules set costs in their order, a
 *        counter with bound 0 forbids its operations, a match costs 0, there is no deletion at the ends of P and, with
 *        the Hamming distance, no insertion and no deletion; on the reverse strand, the ranges are mirrored, the symbol
 *        sets are complemented (see reverse_complement) and the qualities are reversed
 * @param model the model
 * @param P the pattern as it is searched
 * @param quality the qualities of the read (one per symbol)
 * @param reverse_strand whether P is the reverse complement of the read
 * @param bound the largest cost of a valid alignment
 * @param metric the distance metric
 */
inline naive_model_costs naive_model(const random_model_t& model, std::span<const uint8_t> P,
    const std::string& quality, bool reverse_strand, uint64_t bound, distance_metric_t metric)
{
    const uint64_t m = P.size();
    naive_model_costs costs;
    costs.m = m;
    costs.bound = bound;
    costs.matches.assign(m, {});
    costs.sub.assign(m, {});
    costs.ins.assign(m, 1);
    costs.ins_open.assign(m, 0);
    costs.del.assign(m + 1, 1);
    costs.del_open.assign(m + 1, 0);

    // returns the positions [b, e) of P in the range of a rule
    auto positions_of = [&](const model_rule_t& rule) -> std::pair<uint64_t, uint64_t> {
        auto [b, e] = rule.range.positions(m);
        return reverse_strand ? std::make_pair(m - e, m - b) : std::make_pair(b, e);
    };

    // returns whether a symbol set of a rule contains the symbol c (an empty set contains every symbol)
    auto contains = [&](const std::string& set, uint8_t c) {
        const std::string symbols = reverse_strand ? reverse_complement(set) : set;
        return set.empty() || symbols.find(char(c)) != std::string::npos;
    };

    for (uint64_t i = 0; i < m; i++) {
        costs.matches[i][P[i]] = true;
        if (model.iupac) for (char base : iupac_bases(P[i])) costs.matches[i][uint8_t(base)] = true;
        const uint8_t q = model.quality ? uint8_t(quality[reverse_strand ? m - 1 - i : i]) : 0;
        costs.sub[i].fill(!model.quality ? 1 : (q > '!' ? q - '!' : 0));
    }

    for (const model_rule_t& rule : model.rules) {
        if (rule.kind != model_rule_t::_class) continue;
        auto [b, e] = positions_of(rule);

        for (uint64_t i = b; i < e; i++) {
            if (!contains(rule.from, P[i])) continue;
            for (uint64_t c = 0; c < num_byte_values; c++) if (!rule.to.empty() && contains(rule.to, c)) costs.matches[i][c] = true;
        }
    }

    for (const model_rule_t& rule : model.rules) {
        if (rule.kind != model_rule_t::_cost) continue;
        auto [b, e] = positions_of(rule);
        const uint16_t cost = uint16_t(rule.value);

        for (uint64_t i = b; i < e; i++) {
            if (!contains(rule.from, P[i])) continue;
            if (rule.ops & APM_SUB) for (uint64_t c = 0; c < num_byte_values; c++) if (contains(rule.to, c)) costs.sub[i][c] = cost;
            if (rule.ops & APM_INS) costs.ins[i] = cost;
            if (rule.ops & APM_OPEN) costs.ins_open[i] = cost;
        }

        for (uint64_t x = b + 1; x < e; x++) {
            if (rule.ops & APM_DEL) costs.del[x] = cost;
            if (rule.ops & APM_OPEN) costs.del_open[x] = cost;
        }
    }

    for (const model_rule_t& rule : model.rules) {
        if (rule.kind != model_rule_t::_counter) continue;
        auto [b, e] = positions_of(rule);
        naive_model_costs::counter_t counter;
        counter.bound = rule.value;
        counter.sub.assign(m, false);
        counter.ins.assign(m, false);
        counter.ins_open.assign(m, false);
        counter.del.assign(m + 1, false);
        counter.del_open.assign(m + 1, false);

        for (uint64_t i = b; i < e; i++) {
            counter.sub[i] = rule.ops & APM_SUB;
            counter.ins[i] = rule.ops & APM_INS;
            counter.ins_open[i] = rule.ops & APM_OPEN;
        }

        for (uint64_t x = b + 1; x < e; x++) {
            counter.del[x] = rule.ops & APM_DEL;
            counter.del_open[x] = rule.ops & APM_OPEN;
        }

        if (counter.bound > 0) {
            costs.counters.push_back(std::move(counter));
            continue;
        }

        // a counter with bound 0 forbids its operations
        for (uint64_t i = 0; i < m; i++) {
            if (counter.sub[i]) costs.sub[i].fill(apm_forbidden);
            if (counter.ins[i]) costs.ins[i] = apm_forbidden;
            if (counter.ins_open[i]) costs.ins_open[i] = apm_forbidden;
        }

        for (uint64_t x = 0; x <= m; x++) {
            if (counter.del[x]) costs.del[x] = apm_forbidden;
            if (counter.del_open[x]) costs.del_open[x] = apm_forbidden;
        }
    }

    for (uint64_t i = 0; i < m; i++) for (uint64_t c = 0; c < num_byte_values; c++) if (costs.matches[i][c]) costs.sub[i][c] = 0;
    costs.del[0] = costs.del[m] = apm_forbidden;

    if (metric == HAMMING_DISTANCE) {
        std::fill(costs.ins.begin(), costs.ins.end(), apm_forbidden);
        std::fill(costs.del.begin(), costs.del.end(), apm_forbidden);
    }

    return costs;
}

// the last operation of an alignment (see reference_model_costs)
enum last_operation_t : uint64_t { _last_match, _last_insertion, _last_deletion };

/**
 * @brief an alignment state of the reference of an error model: the cost of an alignment and the values of the
 *        counters, after an operation
 */
struct model_alignment_t {
    uint64_t cost = 0;
    std::vector<uint64_t> counts; // [0..|counters|-1] the value of each counter

    /**
     * @brief adds an operation with a cost to the alignment; returns false if the operation is forbidden, the cost
     *        exceeds the bound or a counter exceeds its bound
     * @param costs the model
     * @param cost the cost of the operation
     * @param counted counted(counter) returns whether the counter counts the operation
     */
    template <typename counted_fnc_t>
    bool add(const naive_model_costs& costs, uint64_t cost, counted_fnc_t counted)
    {
        if (cost == apm_forbidden) return false;
        this->cost += cost;

        for (uint64_t j = 0; j < costs.counters.size(); j++) {
            if (counted(costs.counters[j])) counts[j]++;
            if (counts[j] > costs.counters[j].bound) return false;
        }

        return this->cost <= costs.bound;
    }
};

/**
 * @brief returns the cheapest costs of aligning P to the substrings of the text that lie in one sequence, under an
 *        error model: [0..n-1][0..len_max] result[s][len] = the cost of the cheapest valid alignment of P to
 *        text[s..s+len-1] (len >= 1), or bound + 1 if there is none or len = 0 (len_max = min(the rest of the sequence
 *        from s, m + bound)); DP over (row, column, last operation, counter values)
 * @param text the text
 * @param P the pattern
 * @param costs the model compiled for P (see naive_model)
 * @param metric the distance metric
 */
inline std::vector<std::vector<uint64_t>> reference_model_costs(const apm_text_t& text, std::span<const uint8_t> P,
    const naive_model_costs& costs, distance_metric_t metric)
{
    const uint64_t n = text.T.size(), m = P.size();
    // the largest difference of the lengths of P and the aligned substring (an insertion and a deletion cost at least 1)
    const uint64_t kb = metric == HAMMING_DISTANCE ? 0 : costs.bound;
    const uint64_t none = costs.bound + 1;
    std::vector<std::vector<uint64_t>> result(n);

    // an alignment state: the last operation (see last_operation_t) and the values of the counters
    using key_t = std::pair<uint64_t, std::vector<uint64_t>>;
    // the cheapest cost of each alignment state of a DP cell
    using cell_t = std::vector<std::pair<key_t, uint64_t>>;

    for (uint64_t s = 0; s < n; s++) {
        const uint64_t len_max = std::min(text.sequence_rest(s), m + kb);
        /* the DP matrix: [0..m][0..len_max] cell(i, j) = the cheapest costs of the alignment states of P[0..i-1]
           aligned to text[s..s+j-1] */
        std::vector<cell_t> cells((m + 1) * (len_max + 1));
        auto cell = [&](uint64_t i, uint64_t j) -> cell_t& { return cells[i * (len_max + 1) + j]; };

        // stores the cost of an alignment state in the cell c, if it is cheaper than the stored one
        auto relax = [&](cell_t& c, uint64_t state, const model_alignment_t& alignment) {
            for (auto& [key, val] : c) {
                if (key.first == state && key.second == alignment.counts) {
                    val = std::min(val, alignment.cost);
                    return;
                }
            }

            c.emplace_back(key_t { state, alignment.counts }, alignment.cost);
        };

        relax(cell(0, 0), _last_match, model_alignment_t { 0, std::vector<uint64_t>(costs.counters.size(), 0) });

        for (uint64_t i = 0; i <= m; i++) {
            for (uint64_t j = 0; j <= len_max; j++) {
                if (j + kb < i || j > i + kb) continue;

                if (i > 0 && j > 0) {
                    const uint8_t t = text.T[s + j - 1];
                    const bool match = costs.matches[i - 1][t];

                    for (const auto& [key, val] : cell(i - 1, j - 1)) {
                        model_alignment_t alignment { val, key.second };
                        if (alignment.add(costs, costs.sub[i - 1][t], [&](const auto& c) { return !match && c.sub[i - 1]; }))
                            relax(cell(i, j), _last_match, alignment);
                    }
                }

                if (i > 0) {
                    for (const auto& [key, val] : cell(i - 1, j)) {
                        const bool opens = key.first != _last_insertion;
                        model_alignment_t alignment { val, key.second };
                        const bool valid = alignment.add(costs, costs.ins[i - 1], [&](const auto& c) { return c.ins[i - 1]; })
                            && (!opens || alignment.add(costs, costs.ins_open[i - 1], [&](const auto& c) { return c.ins_open[i - 1]; }));
                        if (valid) relax(cell(i, j), _last_insertion, alignment);
                    }
                }

                if (j > 0) {
                    for (const auto& [key, val] : cell(i, j - 1)) {
                        const bool opens = key.first != _last_deletion;
                        model_alignment_t alignment { val, key.second };
                        const bool valid = alignment.add(costs, costs.del[i], [&](const auto& c) { return c.del[i]; })
                            && (!opens || alignment.add(costs, costs.del_open[i], [&](const auto& c) { return c.del_open[i]; }));
                        if (valid) relax(cell(i, j), _last_deletion, alignment);
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

/**
 * @brief returns the cost of the alignment of P to text[pos..pos+len-1] given by a CIGAR under an error model, or
 *        bound + 1 if the CIGAR is no valid alignment of them
 */
inline uint64_t model_cigar_cost(const apm_text_t& text, std::span<const uint8_t> P, uint64_t pos, uint64_t len,
    const cigar_t<uint8_t>& cigar, const naive_model_costs& costs)
{
    const uint64_t invalid = costs.bound + 1;
    if (!text.in_one_sequence(pos, len)) return invalid;
    uint64_t i = 0, j = 0;
    last_operation_t state = _last_match;
    model_alignment_t alignment { 0, std::vector<uint64_t>(costs.counters.size(), 0) };

    for (const cigar_run_t<uint8_t>& run : cigar) {
        for (uint64_t x = 0; x < run.len; x++) {
            bool valid = true;

            if (run.op == cigar_op_t::MATCH || run.op == cigar_op_t::MISMATCH) {
                if (i >= P.size() || j >= len) return invalid;
                const uint8_t t = text.T[pos + j];
                const bool match = costs.matches[i][t];
                // a match is a pair of symbols the model matches (e.g., with symbol classes), not only equal symbols
                if (run.op == cigar_op_t::MATCH && !match) return invalid;
                if (run.op == cigar_op_t::MISMATCH && (t != run.sym || match)) return invalid;
                valid = alignment.add(costs, costs.sub[i][t], [&](const auto& c) { return !match && c.sub[i]; });
                i++; j++; state = _last_match;
            } else if (run.op == cigar_op_t::INS) {
                if (i >= P.size() || P[i] != run.sym) return invalid;
                valid = alignment.add(costs, costs.ins[i], [&](const auto& c) { return c.ins[i]; }) && (state == _last_insertion
                    || alignment.add(costs, costs.ins_open[i], [&](const auto& c) { return c.ins_open[i]; }));
                i++; state = _last_insertion;
            } else {
                if (j >= len || text.T[pos + j] != run.sym) return invalid;
                valid = alignment.add(costs, costs.del[i], [&](const auto& c) { return c.del[i]; }) && (state == _last_deletion
                    || alignment.add(costs, costs.del_open[i], [&](const auto& c) { return c.del_open[i]; }));
                j++; state = _last_deletion;
            }

            if (!valid) return invalid;
        }
    }

    return i == P.size() && j == len ? alignment.cost : invalid;
}

#endif // MOVE_R_TESTS_TEST_INDEX_REFERENCES_HPP
