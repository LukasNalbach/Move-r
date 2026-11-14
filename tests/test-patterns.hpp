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
 * @file test-patterns.hpp
 * @brief random patterns, search schemes and error models for the tests
 */

#ifndef MOVE_R_TESTS_TEST_PATTERNS_HPP
#define MOVE_R_TESTS_TEST_PATTERNS_HPP

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <limits>
#include <random>
#include <regex>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <algorithms/move_rb/apm/model.hpp>
#include <misc/search_schemes.hpp>
#include "test-inputs.hpp"

// ############################# PATTERNS #############################

// apply_random_edits draws max_edits from [0, |pattern| / pattern_symbols_per_max_edit]
inline constexpr uint64_t pattern_symbols_per_max_edit = 2;

/**
 * @brief applies k random edits (substitution, insertion or deletion of one symbol) to a pattern, with k drawn from
 *        [0, max_edits] and max_edits drawn from [0, |pattern| / pattern_symbols_per_max_edit]; the inserted and
 *        substituted symbols are drawn from the positions of input or, if symbols is not empty, with probability 1/2
 *        uniformly from symbols
 * @tparam inp_t container of the pattern and of the input (e.g. std::string or std::vector<int32_t>)
 * @tparam gen_t random number generator type
 * @param pattern the pattern, edited in place (a non-empty pattern remains non-empty)
 * @param input the input the symbols are drawn from
 * @param gen random number generator
 * @param symbols symbols an edit can insert or substitute (e.g. the distinct symbols of input and symbols that do not
 *        occur in it)
 */
template <typename inp_t, typename gen_t>
inline void apply_random_edits(inp_t& pattern, const inp_t& input, gen_t& gen,
    const std::vector<typename inp_t::value_type>& symbols = {})
{
    uint64_t input_size = input.size();
    uint64_t max_edits = std::uniform_int_distribution<uint64_t>(0, pattern.size() / pattern_symbols_per_max_edit)(gen);
    uint64_t num_edits = std::uniform_int_distribution<uint64_t>(0, max_edits)(gen);

    for (uint64_t k = 0; k < num_edits; k++) {
        uint64_t pos = std::uniform_int_distribution<uint64_t>(0, pattern.size() - 1)(gen);
        typename inp_t::value_type sym = symbols.empty() || std::bernoulli_distribution()(gen)
            ? input[std::uniform_int_distribution<uint64_t>(0, input_size - 1)(gen)]
            : symbols[std::uniform_int_distribution<uint64_t>(0, symbols.size() - 1)(gen)];
        // 0: deletion, 1: insertion, 2: substitution (also if the deletion or insertion is not possible)
        uint64_t edit_type = gen() % 3;

        if (edit_type == 0 && pattern.size() > 1)               pattern.erase(pattern.begin() + pos);
        else if (edit_type == 1 && pattern.size() < input_size) pattern.insert(pattern.begin() + pos, sym);
        else                                                    pattern[pos] = sym;
    }
}

/**
 * @brief returns the starting positions of the occurrences of a pattern in an input, by comparing it at every position
 * @tparam inp_t container of the input and of the pattern (e.g. std::string or std::vector<int32_t>)
 * @param input the input
 * @param pattern the pattern
 * @return the starting positions of the occurrences of pattern in input, in increasing order
 */
template <typename inp_t>
inline std::vector<uint64_t> locate_naive(const inp_t& input, const inp_t& pattern)
{
    std::vector<uint64_t> occurrences;
    uint64_t input_size = input.size();
    uint64_t pattern_length = pattern.size();

    for (uint64_t i = 0; i + pattern_length <= input_size; i++) {
        bool match = true;

        for (uint64_t j = 0; j < pattern_length; j++) {
            if (input[i + j] != pattern[j]) {
                match = false;
                break;
            }
        }

        if (match) occurrences.emplace_back(i);
    }

    return occurrences;
}

// a random pattern gets a symbol that does not occur in the text with probability 1 / absent_symbol_inverse_probability
// (see random_pattern)
inline constexpr uint64_t absent_symbol_inverse_probability = 8;

// the places of a random pattern (see random_pattern)
enum pattern_place_t : uint64_t {
    _place_random,         // a random starting position
    _place_text_start,     // the start of the text
    _place_text_end,       // the end of the text
    _place_sequence_start, // the start of a random sequence
    _place_sequence_end,   // the end of a random sequence
    _num_places
};

/**
 * @brief returns the distinct symbols of a text in increasing order
 * @tparam inp_t container of the text
 */
template <typename inp_t>
inline std::vector<typename inp_t::value_type> distinct_symbols(const inp_t& text)
{
    std::vector<typename inp_t::value_type> symbols(text.begin(), text.end());
    std::sort(symbols.begin(), symbols.end());
    symbols.erase(std::unique(symbols.begin(), symbols.end()), symbols.end());
    return symbols;
}

/**
 * @brief returns the symbols the edits of the random patterns of a text draw from (see apply_random_edits): the distinct
 *        symbols of the text in increasing order and the smallest value of the symbol type that does not occur in the
 *        text (if there is one)
 * @tparam inp_t container of the text
 */
template <typename inp_t>
inline std::vector<typename inp_t::value_type> edit_symbols(const inp_t& text)
{
    using sym_t = typename inp_t::value_type;
    std::vector<sym_t> symbols = distinct_symbols(text);
    sym_t absent = std::numeric_limits<sym_t>::min();

    for (sym_t c : symbols) {
        if (c != absent) break;
        if (absent == std::numeric_limits<sym_t>::max()) return symbols;
        absent++;
    }

    symbols.push_back(absent);
    return symbols;
}

/**
 * @brief returns a random pattern drawn from a text: a substring of a log-uniform length in [1, max_length] (with
 *        probability 1 / whole_text_inverse_probability, the text length minus a log-uniform number in [0, max_length -
 *        1]) at a random place (see pattern_place_t, each with the same probability; a substring at the start or at
 *        the end of a non-empty sequence lies in that sequence), with random edits (see apply_random_edits) with
 *        probability 1/2; with probability 1 / absent_symbol_inverse_probability, a random symbol of the pattern is
 *        replaced by the last of symbols
 * @tparam inp_t container of the text and of the pattern
 * @tparam gen_t random number generator type
 * @param text the text (|text| >= 1)
 * @param max_length the largest length of a substring, and the largest difference of the lengths of the text and of a
 *        substring that is drawn as nearly the whole text
 * @param whole_text_inverse_probability see above (0: no substring is drawn as nearly the whole text)
 * @param symbols the symbols an edit draws from with probability 1/2 (see apply_random_edits); the last one does not
 *        occur in the text (see edit_symbols)
 * @param sequence_starts the starting positions of the sequences of the text in increasing order (empty: the text is
 *        one sequence)
 * @param gen random number generator
 */
template <typename inp_t, typename gen_t>
inline inp_t random_pattern(const inp_t& text, uint64_t max_length, uint64_t whole_text_inverse_probability,
    const std::vector<typename inp_t::value_type>& symbols, const std::vector<uint64_t>& sequence_starts, gen_t& gen)
{
    auto uniform = [&](uint64_t min, uint64_t max) { return std::uniform_int_distribution<uint64_t>(min, max)(gen); };
    const uint64_t n = text.size();
    const uint64_t length_limit = std::max<uint64_t>(1, std::min(max_length, n));
    const bool whole_text = whole_text_inverse_probability != 0 && gen() % whole_text_inverse_probability == 0;
    uint64_t m = whole_text ? n + 1 - random_log_uniform_size(1, length_limit, gen)
        : random_log_uniform_size(1, length_limit, gen);
    const pattern_place_t place = pattern_place_t(uniform(0, _num_places - 1));
    // the sequence [beg, end) of a pattern at the start or at the end of a sequence
    uint64_t beg = 0;
    uint64_t end = n;

    if ((place == _place_sequence_start || place == _place_sequence_end) && !sequence_starts.empty()) {
        const uint64_t j = uniform(0, sequence_starts.size() - 1);
        beg = sequence_starts[j];
        end = j + 1 < sequence_starts.size() ? sequence_starts[j + 1] - 1 : n;
        if (end > beg) m = std::min(m, end - beg);
    }

    uint64_t pos = 0;

    switch (place) {
        case _place_random: pos = uniform(0, n - m); break;
        case _place_text_start: pos = 0; break;
        case _place_text_end: pos = n - m; break;
        case _place_sequence_start: pos = std::min(beg, n - m); break;
        default: pos = std::min(std::max(end, m) - m, n - m); break;
    }

    inp_t pattern(text.begin() + pos, text.begin() + pos + m);
    if (std::bernoulli_distribution()(gen)) apply_random_edits(pattern, text, gen, symbols);
    if (!symbols.empty() && gen() % absent_symbol_inverse_probability == 0) pattern[uniform(0, pattern.size() - 1)] = symbols.back();
    return pattern;
}

/**
 * @brief draws a random pattern from an input (see random_pattern, without sequences and without patterns that are
 *        nearly the whole input) and returns it together with the starting positions of its occurrences in the input
 * @tparam inp_t container of the input and of the pattern (e.g. std::string or std::vector<int32_t>)
 * @tparam gen_t random number generator type
 * @param input the input (|input| >= 1)
 * @param max_pattern_length the maximum length of the substring
 * @param symbols the symbols an edit draws from with probability 1/2 (see random_pattern)
 * @param gen random number generator
 * @return the pattern and the starting positions of its occurrences in input, in increasing order
 */
template <typename inp_t, typename gen_t>
inline std::pair<inp_t, std::vector<uint64_t>> random_pattern_and_occurrences(
    const inp_t& input, uint64_t max_pattern_length, const std::vector<typename inp_t::value_type>& symbols, gen_t& gen)
{
    inp_t pattern = random_pattern(input, max_pattern_length, 0, symbols, {}, gen);
    std::vector<uint64_t> occurrences = locate_naive(input, pattern);
    return { std::move(pattern), std::move(occurrences) };
}

// ############################# SEARCH SCHEMES #############################

// the largest number of errors of a random complete search scheme (see random_complete_scheme)
inline constexpr uint64_t max_random_scheme_errors = 4;

/**
 * @brief returns the text form of a search scheme (see parse_search_scheme): the header "p=<p> k=<k>" and one line
 *        "{pi} {L} {U}" per search
 * @param scheme the search scheme
 */
inline std::string search_scheme_text(const search_scheme_t& scheme)
{
    std::string text = "p=" + std::to_string(scheme.p) + " k=" + std::to_string(scheme.k) + "\n";

    for (const search_t& S : scheme.S) {
        for (uint8_t search_step_t::* field : { &search_step_t::part, &search_step_t::k_min, &search_step_t::k_max }) {
            text += "{";
            for (uint64_t j = 0; j < S.size(); j++) text += (j == 0 ? "" : ",") + std::to_string(S[j].*field);
            text += "} ";
        }

        text += "\n";
    }

    return text;
}

// returns whether two search schemes have the same k, p and searches
inline bool equal_search_schemes(const search_scheme_t& a, const search_scheme_t& b)
{
    if (a.k != b.k || a.p != b.p || a.S.size() != b.S.size()) return false;

    for (uint64_t s = 0; s < a.S.size(); s++) {
        if (a.S[s].size() != b.S[s].size()) return false;

        for (uint64_t j = 0; j < a.S[s].size(); j++) {
            if (a.S[s][j].part != b.S[s][j].part || a.S[s][j].k_min != b.S[s][j].k_min
                || a.S[s][j].k_max != b.S[s][j].k_max) return false;
        }
    }

    return true;
}

/**
 * @brief returns whether a search scheme satisfies the rules of parse_search_scheme: 1 <= p <= 256, k <= 255, every
 *        search has p steps with parts in [0, p - 1] and bounds k_min <= k_max <= k, each part after the first of a
 *        search is adjacent to the parts before it, and the bounds of a step are at least those of the previous step
 * @param scheme the search scheme
 */
inline bool satisfies_scheme_rules(const search_scheme_t& scheme)
{
    if (scheme.p < 1 || scheme.p > num_byte_values) return false;

    for (const search_t& S : scheme.S) {
        if (S.size() != scheme.p) return false;
        uint64_t lo = S[0].part;
        uint64_t hi = S[0].part;

        for (uint64_t j = 0; j < S.size(); j++) {
            if (S[j].part >= scheme.p || S[j].k_min > S[j].k_max || S[j].k_max > scheme.k) return false;
            if (j == 0) continue;
            if (S[j].k_min < S[j - 1].k_min || S[j].k_max < S[j - 1].k_max) return false;

            if (S[j].part + 1 == lo) lo--;
            else if (S[j].part == hi + 1) hi++;
            else return false;
        }
    }

    return true;
}

/**
 * @brief returns whether a search accepts a distribution of errors over the parts: after each step, the number of
 *        errors of the parts matched so far lies within the bounds of the step
 * @param S the search
 * @param errors [0..p-1] errors[i] is the number of errors in part i
 */
inline bool search_accepts(const search_t& S, const std::vector<uint64_t>& errors)
{
    uint64_t sum = 0;

    for (const search_step_t& step : S) {
        sum += errors[step.part];
        if (sum < step.k_min || sum > step.k_max) return false;
    }

    return true;
}

/**
 * @brief calls report(errors) for every distribution errors[0..n-1] of at most k errors over n places
 * @param n the number of places
 * @param k the largest number of errors
 * @param report called as report(errors)
 */
template <typename fnc_t>
inline void for_each_error_distribution(uint64_t n, uint64_t k, fnc_t report)
{
    std::vector<uint64_t> errors(n, 0);

    // assigns the errors of the places from place i on, with at most remaining errors left
    std::function<void(uint64_t, uint64_t)> assign = [&](uint64_t i, uint64_t remaining) {
        if (i == n) {
            report(errors);
            return;
        }

        for (errors[i] = 0; errors[i] <= remaining; errors[i]++) assign(i + 1, remaining - errors[i]);
        errors[i] = 0;
    };

    assign(0, k);
}

/**
 * @brief returns the numbers of errors per part that a search counts for an alignment whose errors are distributed over
 *        the p parts and the p - 1 boundaries of adjacent parts (see parse_search_scheme): an error inside a part is
 *        counted in this part, a deletion at the boundary of two parts in the part that the search matches later
 * @param S the search
 * @param errors [0..2p-2] errors[i] is the number of errors inside part i for i < p, errors[p + j] the number of
 *        deletions at the boundary of the parts j and j + 1 for j < p - 1
 * @return [0..p-1] the number of errors counted in part i
 */
inline std::vector<uint64_t> counted_errors(const search_t& S, const std::vector<uint64_t>& errors)
{
    const uint64_t p = S.size();
    // [0..p-1] step_of[i] is the step of S that matches part i
    std::vector<uint64_t> step_of(p);
    for (uint64_t x = 0; x < p; x++) step_of[S[x].part] = x;
    std::vector<uint64_t> counted(errors.begin(), errors.begin() + p);
    for (uint64_t j = 0; j + 1 < p; j++) counted[step_of[j] > step_of[j + 1] ? j : j + 1] += errors[p + j];
    return counted;
}

/**
 * @brief returns a random search scheme for k errors that is complete for the edit distance (see parse_search_scheme):
 *        p is drawn from [1, zero_one_scheme(k).p], then between 1 and p random searches (each part after the first
 *        adjacent to the parts before it, non-decreasing bounds) are drawn, and for every distribution of at most k
 *        errors over the parts and their boundaries that no search accepts (see counted_errors), a search with a random
 *        order of the parts is added whose upper bounds are the numbers of errors it counts for that distribution
 * @param k the number of errors
 * @param gen random number generator
 */
template <typename gen_t>
inline search_scheme_t random_complete_scheme(uint64_t k, gen_t& gen)
{
    auto uniform = [&](uint64_t min, uint64_t max) { return std::uniform_int_distribution<uint64_t>(min, max)(gen); };
    search_scheme_t scheme { .k = uint8_t(k), .p = uniform(1, zero_one_scheme(k).p), .S = {} };
    const uint64_t p = scheme.p;

    // returns a random order of the parts, in which each part after the first is adjacent to the parts before it
    auto random_order = [&] {
        std::vector<uint64_t> order { uniform(0, p - 1) };
        uint64_t lo = order[0];
        uint64_t hi = order[0];

        while (order.size() < p) {
            const bool left = lo > 0 && (hi + 1 == p || std::bernoulli_distribution()(gen));
            order.push_back(left ? --lo : ++hi);
        }

        return order;
    };

    for (uint64_t s = uniform(1, p); s > 0; s--) {
        search_t S;
        uint64_t k_min = 0;
        uint64_t k_max = 0;

        for (uint64_t part : random_order()) {
            k_max = uniform(k_max, k);
            if (std::bernoulli_distribution()(gen)) k_min = uniform(k_min, k_max);
            S.push_back({ uint8_t(part), uint8_t(k_min), uint8_t(k_max) });
        }

        scheme.S.push_back(std::move(S));
    }

    for_each_error_distribution(2 * p - 1, k, [&](const std::vector<uint64_t>& errors) {
        for (const search_t& S : scheme.S) if (search_accepts(S, counted_errors(S, errors))) return;
        search_t S;
        for (uint64_t part : random_order()) S.push_back({ uint8_t(part), 0, 0 });
        const std::vector<uint64_t> counted = counted_errors(S, errors);
        uint64_t sum = 0;

        for (search_step_t& step : S) {
            sum += counted[step.part];
            step.k_max = uint8_t(sum);
        }

        scheme.S.push_back(std::move(S));
    });

    return scheme;
}

/**
 * @brief draws a search scheme for k errors: min-u, pigeonhole, suffix filter, 01 or (for k <= max_random_scheme_errors)
 *        a random complete scheme (see random_complete_scheme), each with the same probability, and its mirror (see
 *        mirrored_scheme) with probability 1/2
 * @param k the number of errors
 * @param gen random number generator
 */
template <typename gen_t>
inline search_scheme_t random_search_scheme(uint64_t k, gen_t& gen)
{
    enum scheme_kind_t : uint64_t { min_u, pigeon_hole, suffix_filter, zero_one, random_complete, num_kinds };
    const uint64_t max_kind = k <= max_random_scheme_errors ? num_kinds - 1 : random_complete - 1;
    search_scheme_t scheme;

    switch (scheme_kind_t(std::uniform_int_distribution<uint64_t>(0, max_kind)(gen))) {
        case min_u: scheme = min_u_scheme(k); break;
        case pigeon_hole: scheme = pigeon_hole_scheme(k); break;
        case suffix_filter: scheme = suffix_filter_scheme(k); break;
        case zero_one: scheme = zero_one_scheme(k); break;
        default: scheme = random_complete_scheme(k, gen); break;
    }

    return std::bernoulli_distribution()(gen) ? mirrored_scheme(scheme) : scheme;
}

/**
 * @brief returns a search scheme with one search changed in one of three ways, each with probability 1/3: the order of
 *        its parts is replaced by a random permutation, the bounds of one step are swapped, or the lower or the upper
 *        bound of one step is set to a random value in [0, k]; it may violate the rules of parse_search_scheme (see
 *        satisfies_scheme_rules)
 * @param scheme the search scheme (with at least one search)
 * @param gen random number generator
 */
template <typename gen_t>
inline search_scheme_t randomly_changed_scheme(search_scheme_t scheme, gen_t& gen)
{
    enum change_t : uint64_t { permuted_parts, swapped_bounds, random_bound, num_changes };
    auto uniform = [&](uint64_t min, uint64_t max) { return std::uniform_int_distribution<uint64_t>(min, max)(gen); };
    search_t& S = scheme.S[uniform(0, scheme.S.size() - 1)];
    search_step_t& step = S[uniform(0, S.size() - 1)];

    switch (change_t(uniform(0, num_changes - 1))) {
        case permuted_parts: {
            std::vector<uint8_t> parts;
            for (const search_step_t& s : S) parts.push_back(s.part);
            std::shuffle(parts.begin(), parts.end(), gen);
            for (uint64_t j = 0; j < S.size(); j++) S[j].part = parts[j];
            break;
        }
        case swapped_bounds: std::swap(step.k_min, step.k_max); break;
        default: (std::bernoulli_distribution()(gen) ? step.k_min : step.k_max) = uint8_t(uniform(0, scheme.k)); break;
    }

    return scheme;
}

// ############################# ERROR MODELS #############################

// the largest cost of a cost rule of a random error model (other than inf)
inline constexpr uint64_t max_model_cost = 3;
// the largest bound of a counter of a random error model
inline constexpr uint64_t max_model_counter = 3;
// the largest number of cost rules, of counters and of classes of a random error model
inline constexpr uint64_t max_model_rules = 3;
// the largest Phred quality of a random pattern quality
inline constexpr uint64_t max_model_quality = 5;
// the largest tau of a random error model (0: chosen per pattern)
inline constexpr uint64_t max_model_tau = 2;
// the operations of a random cost rule are forbidden with probability 1 / forbidden_rule_inverse_probability
inline constexpr uint64_t forbidden_rule_inverse_probability = 6;
// a random error model has the options -iupac and -qual with probability 1 / model_flag_inverse_probability each
inline constexpr uint64_t model_flag_inverse_probability = 4;
// a cost rule has symbols for operations that do not allow them with probability 1 / invalid_symbols_inverse_probability
inline constexpr uint64_t invalid_symbols_inverse_probability = 16;

// the characters with a meaning in the rules of an error model (see parse_apm_cost_rule): they are not drawn as symbols
inline constexpr std::string_view model_syntax_symbols = ",@/";
// the characters of the rules of an error model, from which changed rules draw half of their characters
inline constexpr std::string_view model_rule_characters = "sido=@/,:-0123456789inf";

// the letters of the operations in the rules of an error model, with their operations (see parse_apm_cost_rule)
inline constexpr std::array<std::pair<char, apm_op_t>, 4> model_op_letters {
    std::pair { 's', APM_SUB }, std::pair { 'i', APM_INS }, std::pair { 'd', APM_DEL }, std::pair { 'o', APM_OPEN } };

/**
 * @brief a range of positions of P in the form of a rule of an error model (see parse_apm_range): the whole pattern
 *        (no range), a position, or an interval with an optional start and an optional end
 */
struct model_range_t {
    enum form_t : uint64_t { _whole, _position, _interval, _num_forms };

    form_t form = _whole;
    int64_t beg = 0; // the position (_position) or the start (_interval)
    int64_t end = 0; // the end (_interval)
    bool has_beg = false; // whether an interval has a start
    bool has_end = false; // whether an interval has an end

    // returns the range as it follows '@' in a rule (empty for _whole)
    std::string text() const
    {
        if (form == _whole) return "";
        if (form == _position) return std::to_string(beg);
        return (has_beg ? std::to_string(beg) : "") + ":" + (has_end ? std::to_string(end) : "");
    }

    /**
     * @brief returns the positions [b, e) of a pattern of length m that lie in the range: a negative value x is the
     *        position m + x, a position i is the interval [i, i + 1), and b and e are clamped to [0, m] with b <= e
     * @param m the pattern length
     */
    std::pair<uint64_t, uint64_t> positions(uint64_t m) const
    {
        auto clamped = [&](int64_t x) { return uint64_t(std::clamp<int64_t>(x, 0, int64_t(m))); };
        auto absolute = [&](int64_t x) { return clamped(x < 0 ? int64_t(m) + x : x); };
        if (form == _whole) return { 0, m };
        const uint64_t b = form == _position || has_beg ? absolute(beg) : 0;
        uint64_t e = m;
        if (form == _position) e = beg < 0 ? clamped(int64_t(m) + beg + 1) : clamped(beg + 1);
        else if (has_end) e = absolute(end);
        return { b, std::max(b, e) };
    }
};

/**
 * @brief a rule of an error model as the value of an option of the command line (see apm_query_cli_t::apm_model): a
 *        cost rule (-cost), a counter (-max) or a symbol class (-class)
 */
struct model_rule_t {
    enum kind_t : uint64_t { _cost, _counter, _class, _num_kinds };

    kind_t kind = _cost;
    uint8_t ops = 0; // the operations of a cost rule or a counter (a mask of apm_op_t)
    uint64_t value = 0; // the cost of a cost rule (apm_forbidden: inf) or the bound of a counter
    model_range_t range;
    std::string from; // the symbols of P of a cost rule or a class
    std::string to; // the text symbols of a cost rule or a class

    // returns the option of the rule
    std::string option() const { return kind == _cost ? "-cost" : (kind == _counter ? "-max" : "-class"); }

    // returns the value of the option (see parse_apm_cost_rule, parse_apm_counter and parse_apm_class)
    std::string text() const
    {
        const std::string at = range.form == model_range_t::_whole ? "" : "@" + range.text();
        if (kind == _class) return from + "," + to + at;
        std::string letters;
        for (auto [letter, op] : model_op_letters) if (ops & op) letters += letter;
        const std::string symbols = kind == _cost && (!from.empty() || !to.empty()) ? "/" + from + "," + to : "";
        return letters + "=" + (value == apm_forbidden ? std::string("inf") : std::to_string(value)) + at + symbols;
    }

    // returns whether the rule has symbols for operations that do not allow them (see apm_costs_t::compile)
    bool has_disallowed_symbols() const
    {
        return kind == _cost && ((!to.empty() && (ops & ~APM_SUB)) || (!from.empty() && (ops & (APM_DEL | APM_OPEN))));
    }
};

/**
 * @brief a random error model: its rules in the order of the command line, and the options -iupac, -qual and -tau
 *        (see apm_model_t)
 */
struct random_model_t {
    std::vector<model_rule_t> rules;
    bool iupac = false;
    bool quality = false;
    uint64_t tau = 0; // 0: auto

    // returns the options of the command line as pairs of an option and its value (empty for -iupac and -qual)
    std::vector<std::pair<std::string, std::string>> options() const
    {
        std::vector<std::pair<std::string, std::string>> result;
        for (const model_rule_t& rule : rules) result.emplace_back(rule.option(), rule.text());
        if (iupac) result.emplace_back("-iupac", "");
        if (quality) result.emplace_back("-qual", "");
        result.emplace_back("-tau", tau == 0 ? std::string("auto") : std::to_string(tau));
        return result;
    }

    // returns whether a rule has symbols for operations that do not allow them
    bool has_disallowed_symbols() const
    {
        return std::any_of(rules.begin(), rules.end(), [](const model_rule_t& rule) { return rule.has_disallowed_symbols(); });
    }
};

/**
 * @brief returns a random range of positions of a pattern of length m (see model_range_t): the whole pattern with
 *        probability 1/3, else a position or an interval with values in [-(m + 1), m + 1]
 */
template <typename gen_t>
inline model_range_t random_model_range(uint64_t m, gen_t& gen)
{
    auto value = [&] {
        const int64_t x = int64_t(std::uniform_int_distribution<uint64_t>(0, m + 1)(gen));
        return std::bernoulli_distribution()(gen) ? -x : x;
    };

    model_range_t range;
    range.form = model_range_t::form_t(std::uniform_int_distribution<uint64_t>(0, model_range_t::_num_forms - 1)(gen));
    range.beg = value();
    range.end = value();
    range.has_beg = std::bernoulli_distribution()(gen);
    range.has_end = std::bernoulli_distribution()(gen);
    return range;
}

/**
 * @brief returns a random subset of symbols (empty with probability 1/2 if allow_empty)
 */
template <typename gen_t>
inline std::string random_symbol_set(const std::string& symbols, bool allow_empty, gen_t& gen)
{
    std::string set;
    if (symbols.empty() || (allow_empty && std::bernoulli_distribution()(gen))) return set;
    for (char c : symbols) if (std::bernoulli_distribution()(gen)) set += c;
    if (set.empty()) set += symbols[std::uniform_int_distribution<uint64_t>(0, symbols.size() - 1)(gen)];
    return set;
}

/**
 * @brief returns a random error model for a pattern of length m: up to max_model_rules cost rules, counters and classes
 *        in random order with random operations, costs, bounds and ranges (see random_model_range) and symbol sets
 *        drawn from symbols; the symbols of a cost rule are drawn only for the operations that allow them, except
 *        with probability 1 / invalid_symbols_inverse_probability; -iupac and -qual with probability 1 /
 *        model_flag_inverse_probability each and a tau in [0, max_model_tau]
 * @param symbols the symbols the symbol sets are drawn from (without the characters of model_syntax_symbols)
 * @param m the pattern length
 * @param gen random number generator
 */
template <typename gen_t>
inline random_model_t random_error_model(const std::string& symbols, uint64_t m, gen_t& gen)
{
    auto uniform = [&](uint64_t min, uint64_t max) { return std::uniform_int_distribution<uint64_t>(min, max)(gen); };
    const uint64_t all_ops = APM_SUB | APM_INS | APM_DEL | APM_OPEN;
    random_model_t model;

    for (uint64_t kind = 0; kind < model_rule_t::_num_kinds; kind++) {
        for (uint64_t r = uniform(0, max_model_rules); r > 0; r--) {
            model_rule_t rule;
            rule.kind = model_rule_t::kind_t(kind);
            rule.ops = uint8_t(uniform(1, all_ops));
            rule.range = random_model_range(m, gen);

            if (rule.kind == model_rule_t::_cost) {
                // an insertion and a deletion cost at least 1
                const uint64_t min_cost = rule.ops & (APM_INS | APM_DEL) ? 1 : 0;
                rule.value = gen() % forbidden_rule_inverse_probability == 0 ? apm_forbidden : uniform(min_cost, max_model_cost);
                const bool any_symbols = gen() % invalid_symbols_inverse_probability == 0;
                if (any_symbols || !(rule.ops & (APM_DEL | APM_OPEN))) rule.from = random_symbol_set(symbols, true, gen);
                if (any_symbols || rule.ops == APM_SUB) rule.to = random_symbol_set(symbols, true, gen);
            } else if (rule.kind == model_rule_t::_counter) {
                rule.value = uniform(0, max_model_counter);
            } else {
                rule.from = random_symbol_set(symbols, true, gen);
                rule.to = random_symbol_set(symbols, false, gen);
            }

            model.rules.push_back(std::move(rule));
        }
    }

    std::shuffle(model.rules.begin(), model.rules.end(), gen);
    model.iupac = gen() % model_flag_inverse_probability == 0;
    model.quality = gen() % model_flag_inverse_probability == 0;
    model.tau = uniform(0, max_model_tau);
    return model;
}

/**
 * @brief returns a value of an option with one random change: a character replaced, removed or inserted (each with
 *        probability 1/3); a new character is one of model_rule_characters or (with probability 1/2) a random byte
 */
template <typename gen_t>
inline std::string randomly_changed_value(std::string value, gen_t& gen)
{
    auto uniform = [&](uint64_t min, uint64_t max) { return std::uniform_int_distribution<uint64_t>(min, max)(gen); };
    const char c = std::bernoulli_distribution()(gen) ? char(uniform(0, num_byte_values - 1))
        : model_rule_characters[uniform(0, model_rule_characters.size() - 1)];
    const uint64_t change = value.empty() ? 2 : uniform(0, 2); // 0: replace, 1: remove, 2: insert

    if (change == 0) value[uniform(0, value.size() - 1)] = c;
    else if (change == 1) value.erase(uniform(0, value.size() - 1), 1);
    else value.insert(value.begin() + uniform(0, value.size()), c);
    return value;
}

/**
 * @brief returns whether the value of an option of an error model is valid by the grammar of the rules (see
 *        parse_apm_cost_rule, parse_apm_counter, parse_apm_class and parse_apm_value): -cost "<ops>=<cost>[@<range>]
 *        [/<from>,<to>]" with a cost of at most apm_forbidden, -max "<ops>=<bound>[@<range>]" with a bound of at most
 *        apm_max_counter_bound, -class "<from>,<to>[@<range>]" and -tau "auto", "inf" or a number
 * @param option the option
 * @param value its value
 */
inline bool is_valid_model_value(const std::string& option, const std::string& value)
{
    static const std::string range = "(-?[0-9]+|(-?[0-9]+)?:(-?[0-9]+)?)";
    static const std::regex cost_rule("([sido]+)=(inf|[0-9]+)(@" + range + ")?(/[^,]*,[\\s\\S]*)?");
    static const std::regex counter("([sido]+)=(inf|[0-9]+)(@" + range + ")?");
    static const std::regex symbol_class("[^,@]*,[^@]*(@" + range + ")?");
    static const std::regex tau("auto|inf|[0-9]+");

    // returns the value of the second group of a match (inf: apm_forbidden)
    auto number = [](const std::smatch& match) {
        return match[2] == "inf" ? uint64_t(apm_forbidden) : std::stoull(match[2].str());
    };

    std::smatch match;
    if (option == "-cost") return std::regex_match(value, match, cost_rule) && number(match) <= apm_forbidden;
    if (option == "-max") return std::regex_match(value, match, counter) && number(match) <= apm_max_counter_bound;
    if (option == "-class") return std::regex_match(value, symbol_class);
    if (option == "-tau") return std::regex_match(value, tau);
    return true;
}

/**
 * @brief parses the options of an error model as move-rb reads them (see read_apm_query_option)
 * @param options pairs of an option and its value
 * @throws std::invalid_argument if a value is invalid
 */
inline apm_model_t parse_model_options(const std::vector<std::pair<std::string, std::string>>& options)
{
    apm_model_t model;

    for (const auto& [option, value] : options) {
        if (option == "-cost") model.costs.push_back(parse_apm_cost_rule(value));
        else if (option == "-max") model.counters.push_back(parse_apm_counter(value));
        else if (option == "-class") model.classes.push_back(parse_apm_class(value));
        else if (option == "-tau") model.tau = value == "auto" ? 0 : uint16_t(parse_apm_value(value, "tau"));
        else if (option == "-iupac") model.iupac = true;
        else if (option == "-qual") model.quality = true;
    }

    return model;
}

/**
 * @brief returns a random Phred+33 quality string of length m with qualities in [0, max_model_quality]
 */
template <typename gen_t>
inline std::string random_quality(uint64_t m, gen_t& gen)
{
    std::string quality(m, '!');
    for (char& c : quality) c = char('!' + std::uniform_int_distribution<uint64_t>(0, max_model_quality)(gen));
    return quality;
}

#endif // MOVE_R_TESTS_TEST_PATTERNS_HPP
