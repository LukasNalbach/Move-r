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

#ifndef MOVE_R_MOVE_RB_CONSTRUCTION_CONTEXT_TABLES_TPP
#define MOVE_R_MOVE_RB_CONSTRUCTION_CONTEXT_TABLES_TPP

#include <misc/sort.hpp>

#include <data_structures/blocked_positions.hpp>
#include <data_structures/elias_fano_bit_vector.hpp>
#include <data_structures/rank_select_set.hpp>

#include <move_rb/move_rb.hpp>

template <move_r_support support>
template <direction_t dir, typename start_t, typename report_t, typename end_t>
void move_rb<support>::extend_level_parents(uint64_t len, uint64_t num_parents,
    const no_init_vector<uint32_t>& parent_perm, uint64_t round_parents, uint64_t p, start_t start_chunk,
    report_t report, end_t end_round) const
{
    constexpr query_support_t qs = supports_locate ? LOCATE : COUNT;
    bool next_round = true;

    for (uint64_t first = 0; first < num_parents && next_round; first += round_parents) {
        const uint64_t last = std::min<uint64_t>(num_parents, first + round_parents);
        const uint64_t chunks = std::min<uint64_t>(p, last - first);

        #pragma omp parallel for num_threads(p) schedule(dynamic, 1)
        for (uint64_t c = 0; c < chunks; c++) {
            start_chunk(c);
            const uint64_t beg = first + (c * (last - first)) / chunks;
            const uint64_t end = first + ((c + 1) * (last - first)) / chunks;

            for (uint64_t parent = beg; parent < end; parent++) {
                search_context_t<qs> parent_ctx;

                if (len == 1) parent_ctx = empty_context<qs>();
                else load_context<qs>(parent_ctx, len - 1, uint32_t(parent), parent_perm[parent]);

                auto ec = parent_ctx.extend_phase();
                ec.template prepare_extend_all<dir, NO_TERM>();

                while (ec.can_extend()) {
                    uint64_t sym = ec.next_internal_symbol();
                    report(c, parent, sym, ec.template extend_next<dir>());
                }
            }
        }

        next_round = end_round(chunks);
    }
}

template <move_r_support support>
void move_rb<support>::build_level_samples(::context_tables::level_t& lvl,
    const std::vector<sample_info_pack_t>& samples)
{
    uint64_t max_o = 0;
    uint64_t max_i = 0;
    uint64_t max_x = 0;

    for (const sample_info_pack_t& smpl : samples) {
        max_o = std::max<uint64_t>(max_o, smpl.o());
        max_i = std::max<uint64_t>(max_i, smpl.i());
        max_x = std::max<uint64_t>(max_x, smpl.x());
    }

    lvl.samples = interleaved_bit_aligned_vectors({
        2, 1,
        std::max<uint8_t>(1, std::bit_width(max_o)),
        std::max<uint8_t>(1, std::bit_width(max_i)),
        std::max<uint8_t>(1, std::bit_width(max_x)),
        0, 0, 0 });
    lvl.samples.resize_no_init(samples.size());

    for (uint64_t i = 0; i < samples.size(); i++) {
        lvl.samples.template set<_smpl_d>(i, uint64_t(samples[i].d()));
        lvl.samples.template set<_smpl_t>(i, uint64_t(samples[i].t()));
        lvl.samples.template set<_smpl_o>(i, samples[i].o());
        lvl.samples.template set<_smpl_i>(i, samples[i].i());
        lvl.samples.template set<_smpl_x>(i, samples[i].x());
    }
}

template <move_r_support support>
void move_rb<support>::build_reverse_order(::context_tables::level_t& lvl,
    std::vector<uint64_t>& beg_R, std::vector<uint64_t>& width, const std::vector<uint64_t>& primary_keys,
    const std::vector<uint64_t>& secondary_keys, const no_init_vector<uint32_t>& parent_perm, uint64_t num_parents,
    const std::vector<uint64_t>& packed_rev, no_init_vector<uint32_t>& perm, std::vector<uint64_t>& packed_rev_sorted,
    uint64_t p) const
{
    const uint64_t num_entries = lvl.num_entries;
    auto& rev = lvl.orders[_reverse_order];

    // [0..num_entries-1] (beg_R[i], i) sorted by beg_R[i]
    std::vector<std::pair<uint64_t, uint32_t>> by_beg_R;
    by_beg_R.reserve(num_entries);
    for (uint64_t i = 0; i < num_entries; i++)
        by_beg_R.emplace_back(beg_R[i], uint32_t(i));

    radix_sort(by_beg_R.begin(), by_beg_R.end(), key_of_first {}, p);
    beg_R = std::vector<uint64_t>();

    std::vector<uint64_t> beg_R_sorted; // [0..num_entries] the left limits in reverse order, then n
    std::vector<uint64_t> keys_rev; // [0..num_entries-1] the keys of the entries in reverse order
    beg_R_sorted.reserve(num_entries + 1);
    keys_rev.reserve(num_entries);
    if (!packed_rev.empty()) packed_rev_sorted.reserve(num_entries);
    resize_no_init(perm, num_entries);

    // the right limit of an interval is stored only if the next interval does not start at the next position
    uint64_t prev_e_R = 0;
    bool has_prev_R = false;

    for (uint64_t j = 0; j < num_entries; j++) {
        const uint64_t i = by_beg_R[j].second;
        const uint64_t b_R = by_beg_R[j].first;
        perm[i] = uint32_t(j);

        if (has_prev_R && prev_e_R + 1 != b_R)
            rev.end_exceptions.emplace_back(uint32_t(j - 1), prev_e_R);

        prev_e_R = b_R + width[i];
        has_prev_R = true;
        beg_R_sorted.emplace_back(b_R);

        // in reverse order, the key of the primary direction is parent_perm[parent] * sigma + sym
        const uint64_t sym = primary_keys[i] / num_parents;
        const uint64_t parent = primary_keys[i] - sym * num_parents;
        keys_rev.emplace_back(uint64_t(parent_perm[parent]) * sigma + sym);
        if (!packed_rev.empty()) packed_rev_sorted.emplace_back(packed_rev[i]);
    }

    if (has_prev_R && prev_e_R + 1 != n) rev.end_exceptions.emplace_back(uint32_t(num_entries - 1), prev_e_R);
    beg_R_sorted.emplace_back(n);
    width = std::vector<uint64_t>();

    rev.beg = blocked_positions(beg_R_sorted);
    beg_R_sorted = std::vector<uint64_t>();
    rev.primary_keys = rank_select_set(keys_rev, sigma * num_parents + 1, p);

    // in reverse order, the key of the secondary direction is sym * num_parents + parent_perm[parent]
    if (!secondary_keys.empty()) {
        for (uint64_t j = 0; j < num_entries; j++) {
            const uint64_t i = by_beg_R[j].second;
            const uint64_t parent = secondary_keys[i] / sigma;
            const uint64_t sym = secondary_keys[i] - parent * sigma;
            keys_rev[j] = sym * num_parents + parent_perm[parent];
        }

        rev.secondary_keys = rank_select_set(keys_rev, sigma * num_parents + 1, p);
    }
}

template <move_r_support support>
void move_rb<support>::build_context_tables(double share, uint64_t base,
    bool single_length, uint64_t p)
{
    _context_tables.clear();
    if (share <= 0) return;

    constexpr query_support_t qs = supports_locate ? LOCATE : COUNT;
    constexpr direction_t primary_dir = ctx_primary_dir;
    constexpr direction_t secondary_dir = ctx_secondary_dir;
    using ctx_t = search_context_t<qs>;

    // upper bound for the size of an entry of a built level in bytes
    const uint64_t bytes_per_entry = 2 * (1 + 1 + 1) + 1;
    // memory used per entry during the construction of a level in bytes
    constexpr uint64_t bytes_per_row = sizeof(context_table_row_t) + 9 * sizeof(uint64_t);
    // number of parents per thread, that are extended in one round (see extend_level_parents)
    constexpr uint64_t round_parents_per_thread = 4096;

    // maximum size of the tables in bytes
    uint64_t budget = uint64_t(share * double(base == 0 ? size_in_bytes() : base));
    uint64_t used = 0; // size of the built levels in bytes
    double last_bytes_per_entry = 0; // size per entry of the last built level in bytes
    std::vector<uint32_t> entry_of_symbol; // [0..sigma-1] the entry of the pattern c in level 1
    _context_tables.set_alphabet_size(sigma);

    /* single-length mode: only the levels 1..num_short_levels and the last level that fits into the budget are
       stored; the last level is keyed by the packed patterns (see context_tables::single_len), so sigma^len < 2^64
       has to hold. A pattern P[0..len-1] is packed as sum_i P[i] * sigma^(len-1-i), where P[0] is the symbol of the
       last primary extension. */
    typename ::context_tables::level_t single_length_level; // the last level that fits into the budget
    std::vector<typename ::context_tables::level_t> short_levels; // the levels 1..num_short_levels
    static constexpr uint64_t num_short_levels = 2;
    std::vector<uint64_t> single_packed; // the packed patterns of single_length_level, in forward order
    std::vector<uint64_t> single_packed_rev; // the packed reversed patterns of single_length_level, in reverse order
    uint64_t single_len = 0; // length of the patterns of single_length_level
    std::vector<uint64_t> parent_packed(1, 0); // the packed patterns of the previous level, in forward order
    std::vector<uint64_t> parent_packed_rev(1, 0); // the packed reversed patterns of the previous level, in forward order
    // [0..num_parents-1] parent_perm[i] = entry in reverse order of the i-th entry in forward order of the previous level
    no_init_vector<uint32_t> parent_perm(1, 0);
    constexpr uint64_t word_bits = 8 * sizeof(uint64_t);
    const uint64_t max_single_len = sigma < 2 ? word_bits
        : uint64_t(double(word_bits - 1) / std::log2(double(sigma)));

    // sigma^(len-1), the factor of the symbol of a primary extension in a packed pattern
    uint64_t sigma_pow = 1;

    for (uint64_t len = 1;; len++) {
        uint64_t num_parents = _context_tables.num_entries(len - 1);
        if (used >= budget) break;

        // maximum number of entries of the level, limited by the budget and by the memory used to build it
        uint64_t max_entries = std::min<uint64_t>({
            (single_length ? budget : budget - used) / bytes_per_entry,
            (16 * budget) / bytes_per_row,
            uint64_t(no_ctx_entry) - 1 });

        // the level has at least num_parents - 1 entries; stop if they are estimated to exceed the budget
        if (last_bytes_per_entry > 0 && 0.95 * double(num_parents - 1) * last_bytes_per_entry >
                double(single_length ? budget : budget - used)) break;

        // 1. extend the parents in the primary direction

        /* bucket[c] stores the patterns obtained by a primary extension with c, sorted by parent; the
           concatenation of the buckets is the level in forward order */
        std::vector<std::vector<context_table_row_t>> bucket(sigma);
        uint64_t num_entries = 0;
        bool fits = true;

        // row_chunk[c][sym] stores the patterns of the c-th chunk of a round obtained by an extension with sym
        const uint64_t round_parents = std::max<uint64_t>(1, p * round_parents_per_thread);
        std::vector<std::vector<std::vector<context_table_row_t>>> row_chunk(p);
        for (auto& per_thread : row_chunk) per_thread.resize(sigma);

        extend_level_parents<primary_dir>(len, num_parents, parent_perm, round_parents, p,
            [&](uint64_t c) { for (auto& per_sym : row_chunk[c]) per_sym.clear(); },
            [&](uint64_t c, uint64_t parent, uint64_t sym, ctx_t child) {
                context_table_row_t& row = row_chunk[c][sym].emplace_back();
                row.b = child.b;
                row.e = child.e;
                row.parent = uint32_t(parent);
                row.b_R = child.b_R;

                // the SA-sample is stored only if the interval contains no SA-sample of the index
                if constexpr (supports_locate && qs == LOCATE) {
                    child.template update_input_intervals<LEFT>();
                    const auto& [cb, ce, cb_R, ce_R, cb_, ce_, cb_R_, ce_R_, cs] =
                        child.template const_vars<LEFT>();

                    row.stores_sample = !ctx_t::sample_inside(*this, cb, ce, cb_);
                    if (row.stores_sample) [[unlikely]] row.s = cs;
                } else {
                    row.stores_sample = false;
                }
            },
            [&](uint64_t chunks) {
                // append the chunks to the buckets in chunk order (the buckets remain sorted by parent)
                for (uint64_t c = 0; c < chunks; c++) {
                    for (uint64_t sym = 0; sym < sigma; sym++) {
                        std::vector<context_table_row_t>& rows = row_chunk[c][sym];
                        num_entries += rows.size();

                        if (num_entries > max_entries) [[unlikely]] { fits = false; break; }
                        bucket[sym].insert(bucket[sym].end(), rows.begin(), rows.end());
                    }

                    if (!fits) break;
                }

                return fits;
            });

        row_chunk = std::vector<std::vector<std::vector<context_table_row_t>>>();
        if (!fits || num_entries == 0) break;

        // 2. the data of the entries in forward order

        typename ::context_tables::level_t lvl;
        lvl.num_entries = num_entries;
        auto& fwd = lvl.orders[_forward_order];

        std::vector<uint64_t> beg_vals; // [0..num_entries] the left limits of the intervals in forward order, then n
        // [0..num_entries-1] the key sym * num_parents + parent of the primary direction of each entry in forward order
        std::vector<uint64_t> primary_keys;
        // [0..num_entries-1] the key parent * sigma + sym of the secondary direction of each entry in forward order
        std::vector<uint64_t> secondary_keys;
        std::vector<uint64_t> packed; // [0..num_entries-1] the packed patterns in forward order (single-length mode)
        std::vector<uint64_t> packed_rev; // [0..num_entries-1] the packed reversed patterns in forward order (single-length mode)
        std::vector<uint64_t> sample_marks; // the entries (in forward order) whose interval contains no SA-sample
        std::vector<uint64_t> beg_R_vals; // [0..num_entries-1] the left limits of the intervals in the suffix array of T^R
        std::vector<uint64_t> width; // [0..num_entries-1] e - b of each interval
        std::vector<sample_info_pack_t> kept_smpl; // the SA-samples of the entries in sample_marks
        beg_vals.reserve(num_entries + 1);
        primary_keys.reserve(num_entries);
        if (single_length) packed.reserve(num_entries);

        beg_R_vals.reserve(num_entries);
        width.reserve(num_entries);
        if (single_length) packed_rev.reserve(num_entries);

        /* the right limit of an interval is stored only if the next interval does not start at the next
           position (see order_t::end_exceptions) */
        uint64_t prev_e = 0; // right limit of the previous interval
        bool has_prev = false;

        for (uint64_t sym = 0; sym < sigma; sym++) {
            for (const context_table_row_t& row : bucket[sym]) {
                if (has_prev && prev_e + 1 != row.b)
                    fwd.end_exceptions.emplace_back(uint32_t(beg_vals.size() - 1), prev_e);

                prev_e = row.e;
                has_prev = true;

                beg_R_vals.emplace_back(row.b_R);
                width.emplace_back(row.e - row.b);
                if (single_length) packed_rev.emplace_back(parent_packed_rev[row.parent] * sigma + sym);

                if constexpr (supports_locate && qs == LOCATE) {
                    if (row.stores_sample) [[unlikely]] {
                        sample_marks.emplace_back(beg_vals.size());
                        kept_smpl.emplace_back(row.s);
                    }
                }

                beg_vals.emplace_back(row.b);
                primary_keys.emplace_back(sym * num_parents + row.parent);
                if (single_length) packed.emplace_back(sym * sigma_pow + parent_packed[row.parent]);
            }

            bucket[sym] = std::vector<context_table_row_t>();
        }

        if (has_prev && prev_e + 1 != n) fwd.end_exceptions.emplace_back(uint32_t(num_entries - 1), prev_e);
        beg_vals.emplace_back(n);
        if (!sample_marks.empty()) lvl.sample_marks = elias_fano_bit_vector(sample_marks, num_entries, p);
        sample_marks = std::vector<uint64_t>();

        // store the SA-samples with the bit widths of their largest values
        if constexpr (supports_locate && qs == LOCATE) {
            build_level_samples(lvl, kept_smpl);
            kept_smpl = std::vector<sample_info_pack_t>();
        }

        fwd.primary_keys = rank_select_set(primary_keys, sigma * num_parents + 1, p);

        // 3. the keys of the secondary direction

        /* extending the parents in forward order with each symbol in the secondary direction yields the entries of
           the level in forward order; if the interval of an extension differs from the interval of the entry, no
           keys of the secondary direction are stored */
        if (!single_length) {
            secondary_keys.reserve(num_entries);
            bool ordered = true;

            // the keys and the left limits of the extensions of chunk c of a round
            std::vector<std::vector<uint64_t>> sec_chunk_keys(p);
            std::vector<std::vector<uint64_t>> sec_chunk_begs(p);

            extend_level_parents<secondary_dir>(len, num_parents, parent_perm, round_parents, p,
                [&](uint64_t c) { sec_chunk_keys[c].clear(); sec_chunk_begs[c].clear(); },
                [&](uint64_t c, uint64_t parent, uint64_t sym, ctx_t child) {
                    sec_chunk_keys[c].emplace_back(parent * sigma + sym);
                    sec_chunk_begs[c].emplace_back(child.b);
                },
                [&](uint64_t chunks) {
                    // compare the extensions with the entries in forward order
                    for (uint64_t c = 0; c < chunks && ordered; c++) {
                        for (uint64_t k = 0; k < sec_chunk_keys[c].size(); k++) {
                            if (secondary_keys.size() >= num_entries || beg_vals[secondary_keys.size()] != sec_chunk_begs[c][k]) [[unlikely]] {
                                ordered = false;
                                break;
                            }

                            secondary_keys.emplace_back(sec_chunk_keys[c][k]);
                        }
                    }

                    return ordered;
                });

            sec_chunk_keys = std::vector<std::vector<uint64_t>>();
            sec_chunk_begs = std::vector<std::vector<uint64_t>>();
            if (!ordered || secondary_keys.size() != num_entries) secondary_keys = std::vector<uint64_t>();
            else fwd.secondary_keys = rank_select_set(secondary_keys, sigma * num_parents + 1, p);
        }

        fwd.beg = blocked_positions(beg_vals);
        beg_vals = std::vector<uint64_t>();

        // 4. the level in reverse order, i.e., sorted by the left limits of the intervals in the suffix array of T^R
        no_init_vector<uint32_t> perm; // [0..num_entries-1] perm[i] = entry in reverse order of the i-th entry in forward order
        std::vector<uint64_t> packed_rev_sorted; // packed_rev in reverse order

        build_reverse_order(lvl, beg_R_vals, width, primary_keys, secondary_keys, parent_perm, num_parents,
            packed_rev, perm, packed_rev_sorted, p);

        secondary_keys = std::vector<uint64_t>();
        primary_keys = std::vector<uint64_t>();

        // 5. add the level
        if (len == 1) {
            entry_of_symbol.resize(sigma, no_ctx_entry);

            // in level 1, the key of the primary direction is the symbol
            for (uint32_t i = 0; i < num_entries; i++)
                entry_of_symbol[fwd.primary_keys.select_1(i + 1)] = i;

            used += entry_of_symbol.size() * sizeof(uint32_t);
        }

        const uint64_t lvl_bytes = lvl.size_in_bytes();
        if (!single_length && used + lvl_bytes > budget) break;
        used += lvl_bytes;
        last_bytes_per_entry = double(lvl_bytes) / double(num_entries);
        sigma_pow *= sigma;

        if (single_length) {
            // each level with len > num_short_levels replaces the previous single_length_level
            if (len > num_short_levels) {
                single_length_level = lvl;
                single_packed = packed;
                single_packed_rev = packed_rev_sorted;
                single_len = len;
            }

            parent_packed = std::move(packed);
            parent_packed_rev = std::move(packed_rev);
            parent_perm = std::move(perm);
            if (len <= num_short_levels) short_levels.emplace_back(lvl);
            _context_tables.push_level(std::move(lvl));
            if (len >= max_single_len) break;
            continue;
        }

        parent_perm = std::move(perm);
        _context_tables.push_level(std::move(lvl));
    }

    if (single_length) {
        // the tables of the single-length mode: the levels 1..num_short_levels and single_length_level
        _context_tables.clear();
        _context_tables.set_alphabet_size(sigma);

        if (single_len == 0 || short_levels.empty()) return;

        for (auto& l : short_levels) _context_tables.push_level(std::move(l));
        _context_tables.set_symbol_entries(std::move(entry_of_symbol));

        // replace the keys of single_length_level by the packed patterns (see context_tables::single_len)
        single_length_level.orders[_forward_order].secondary_keys = rank_select_set();
        single_length_level.orders[_forward_order].primary_keys = rank_select_set(single_packed, single_packed.back() + 1, p);
        single_length_level.orders[_reverse_order].secondary_keys = rank_select_set();

        single_length_level.orders[_reverse_order].primary_keys =
            rank_select_set(single_packed_rev, single_packed_rev.back() + 1, p);

        _context_tables.set_single_length(single_len);
        _context_tables.push_level(std::move(single_length_level));
        return;
    }

    if (_context_tables.empty()) return;
    _context_tables.set_symbol_entries(std::move(entry_of_symbol));
}

#endif // MOVE_R_MOVE_RB_CONSTRUCTION_CONTEXT_TABLES_TPP
