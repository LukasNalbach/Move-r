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

#ifndef MOVE_R_DATA_STRUCTURES_CONTEXT_TABLES_TPP
#define MOVE_R_DATA_STRUCTURES_CONTEXT_TABLES_TPP

#include <map>

#include <data_structures/context_tables.hpp>

inline void context_tables::log_data_structures(const std::function<std::string(uint64_t)>& symbol,
    bool primary_prepends, uint64_t max_depth, uint64_t max_entries) const
{
    using pattern_t = std::vector<uint64_t>; // a pattern in the internal effective alphabet
    const uint64_t num_levels = max_length();
    const bool has_reverse = !levels.empty() && !levels[0].orders[_reverse_order].beg.empty();
    std::cout << "sigma = " << sigma << std::endl;
    std::cout << "maximum length = " << num_levels << std::endl;
    for (uint64_t len = 1; len <= num_levels; len++)
        std::cout << "patterns of length " << len << " = " << format_number(num_entries(len)) << std::endl;
    if (single_len != 0) {
        std::cout << "patterns of length " << single_len << " = " << format_number(levels.back().num_entries)
                  << std::endl;
    }

    std::cout << "reverse order = " << (has_reverse ? "yes" : "no") << std::endl;
    std::cout << "size = " << format_size(size_in_bytes()) << std::endl;
    if (levels.empty()) return;

    // returns a pattern as it is logged: its symbols in double quotes, if each is one character, else separated by
    // spaces in parentheses
    auto pattern_string = [&](const pattern_t& pattern) {
        std::vector<std::string> symbols;
        for (uint64_t sym : pattern) symbols.push_back(symbol(sym));
        const bool characters = std::all_of(symbols.begin(), symbols.end(),
            [](const std::string& s) { return s.size() == 1; });
        std::string text;
        for (const std::string& s : symbols) text += (characters || text.empty() ? "" : " ") + s;
        return characters ? "\"" + text + "\"" : "(" + text + ")";
    };

    // returns the pattern of a parent extended by sym in primary direction
    auto extended = [&](const pattern_t& parent, uint64_t sym) {
        pattern_t pattern = parent;
        if (primary_prepends) pattern.insert(pattern.begin(), sym);
        else pattern.push_back(sym);
        return pattern;
    };

    /* the patterns of the entries of the lengths 1 to num_levels in an order, computed from the keys of their
       parents: [0..num_levels] patterns[len][e] is the pattern of entry e of length len; the levels whose keys are
       not stored are left empty */
    auto patterns_in = [&]<context_order_t order>() {
        std::vector<std::vector<pattern_t>> patterns(num_levels + 1);
        patterns[0].emplace_back();
        patterns[1].resize(num_entries(1));

        for (uint64_t sym = 0; sym < symbol_entries.size(); sym++)
            if (symbol_entries[sym] != no_entry) patterns[1][symbol_entries[sym]] = { sym };

        for (uint64_t len = 2; len <= num_levels && has_parent<_primary_extension, order>(len); len++) {
            const rank_select_set& k = keys<_primary_extension, order>(len - 1);
            const uint64_t num_parents = num_entries(len - 1);
            patterns[len].resize(num_entries(len));

            for (uint64_t e = 0; e < num_entries(len); e++) {
                const uint64_t key = k.select_1(e + 1);
                const bool symbol_major = symbol_major_keys(order, _primary_extension);
                const uint64_t sym = symbol_major ? key / num_parents : key % sigma;
                const uint64_t parent_entry = symbol_major ? key % num_parents : key / sigma;
                patterns[len][e] = extended(patterns[len - 1][parent_entry], sym);
            }
        }

        return patterns;
    };

    const std::vector<std::vector<pattern_t>> patterns = patterns_in.template operator()<_forward_order>();
    // the entry in reverse order of each pattern
    std::map<pattern_t, uint32_t> reverse_entry;

    if (has_reverse) {
        const std::vector<std::vector<pattern_t>> patterns_R = patterns_in.template operator()<_reverse_order>();

        for (uint64_t len = 1; len <= num_levels; len++)
            for (uint64_t e = 0; e < patterns_R[len].size(); e++) reverse_entry[patterns_R[len][e]] = uint32_t(e);
    }

    // returns the logged intervals of the entry of a pattern in forward order (see log_data_structures)
    auto intervals_string = [&](uint64_t len, uint32_t entry, uint32_t entry_R) {
        const auto [b, e] = interval<_forward_order>(len, entry);
        std::string text = "[" + std::to_string(b) + ", " + std::to_string(e) + "]";

        if (entry_R != no_entry) {
            const auto [b_R, e_R] = interval<_reverse_order>(len, entry_R);
            text += ", reversed: [" + std::to_string(b_R) + ", " + std::to_string(e_R) + "]";
        }

        return text + (sample_rank(len, entry) != no_entry ? ", sample" : "");
    };

    // [0..num_levels-1] children[len][parent] are the entries of length len + 1 whose parent is the entry parent
    std::vector<std::vector<std::vector<uint32_t>>> children(num_levels);

    for (uint64_t len = 1; len <= num_levels && !patterns[len].empty(); len++) {
        children[len - 1].resize(num_entries(len - 1));

        for (uint64_t e = 0; e < num_entries(len); e++) {
            const uint32_t parent_entry = len == 1 ? 0 : parent<_primary_extension, _forward_order>(len, uint32_t(e));
            children[len - 1][parent_entry].push_back(uint32_t(e));
        }
    }

    std::cout << std::endl << "patterns of the lengths 1 to " << num_levels << ":" << std::endl;
    tree_log tree(max_depth, max_entries);

    // logs the pattern of an entry of length len and the patterns below it
    auto log_pattern = [&](auto& self, uint64_t len, uint32_t entry) -> void {
        const auto it = reverse_entry.find(patterns[len][entry]);
        const uint32_t entry_R = it == reverse_entry.end() ? no_entry : it->second;
        const std::string text = pattern_string(patterns[len][entry]) + ": " + intervals_string(len, entry, entry_R);
        if (!tree.node(len - 1, text)) return;
        if (len == num_levels || children[len].empty()) return;
        for (uint32_t child : children[len][entry]) self(self, len + 1, child);
    };

    if (!children[0].empty()) for (uint32_t entry : children[0][0]) log_pattern(log_pattern, 1, entry);
    tree.limits();
    if (single_len == 0) return;

    // the patterns of the single-length level, decoded from their packed keys (see packed_key)
    const level_t& single = levels.back();
    std::vector<pattern_t> single_patterns(std::min<uint64_t>(single.num_entries, max_entries));

    for (uint64_t e = 0; e < single_patterns.size(); e++) {
        uint64_t key = single.orders[_forward_order].primary_keys.select_1(e + 1);
        single_patterns[e].resize(single_len);

        for (uint64_t i = single_len; i-- > 0;) {
            single_patterns[e][i] = key % sigma;
            key /= sigma;
        }
    }

    // returns the entry in reverse order of a pattern of the single-length level
    auto single_entry_R = [&](const pattern_t& pattern) {
        if (!has_reverse) return no_entry;
        const uint64_t key_R = packed_key([&](uint64_t i) { return pattern[single_len - 1 - i]; }, single_len);
        return lookup<_reverse_order>(key_R);
    };

    std::cout << std::endl << "patterns of length " << single_len << ":" << std::endl;

    log_limited_table(single.num_entries, max_entries, [&](aligned_table& table, uint64_t num) {
        table.add_row("pattern:", num, [&](uint64_t e) { return pattern_string(single_patterns[e]); });
        table.add_row("b:", num, [&](uint64_t e) { return interval<_forward_order>(single_len, uint32_t(e)).first; });
        table.add_row("e:", num, [&](uint64_t e) { return interval<_forward_order>(single_len, uint32_t(e)).second; });

        if (has_reverse) {
            table.add_row("b_R:", num, [&](uint64_t e) {
                return interval<_reverse_order>(single_len, single_entry_R(single_patterns[e])).first; });
            table.add_row("e_R:", num, [&](uint64_t e) {
                return interval<_reverse_order>(single_len, single_entry_R(single_patterns[e])).second; });
        }

        table.add_row("sample:", num, [&](uint64_t e) {
            return sample_rank(single_len, uint32_t(e)) != no_entry ? "yes" : "no"; });
    });
}

#endif // MOVE_R_DATA_STRUCTURES_CONTEXT_TABLES_TPP
