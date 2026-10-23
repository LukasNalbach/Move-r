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

#ifndef MOVE_R_MOVE_RB_QUERIES_TPP
#define MOVE_R_MOVE_RB_QUERIES_TPP

#include <move_rb/move_rb.hpp>

template <move_r_support support>
template <typename fnc_t>
inline uint64_t move_rb<support>::input_interval(uint64_t i, const elias_fano_bit_vector& sd_arr, fnc_t interval_start)
{
    uint64_t x = (std::max<uint64_t>(1, sd_arr.rank_1(i)) - 1) * input_interval_sample_rate;

    while (i >= interval_start(x + 1)) {
        x++;
    }

    return x;
}

template <move_r_support support>
template <direction_t dir>
inline uint64_t move_rb<support>::input_interval_near(uint64_t i, uint64_t x) const
{
    const auto& M_LF = index<dir>().M_LF();
    auto p = [&](uint64_t y) { return M_LF.p(y); };

    // a bidirectional M_LF searches near x itself (see move_data_structure::interval_near)
    if constexpr (std::remove_cvref_t<decltype(M_LF)>::has_inverse) {
        uint64_t near = M_LF.interval_near(i, x);
        if (near != M_LF.num_intervals()) [[likely]] return near;
        return mlf_interval<dir>(i);
    }

    // exponential search over at most max_scanned_intervals input intervals to the right or to the left of x
    if (p(x) <= i) {
        if (i < p(x + 1)) [[likely]] return x;

        uint64_t hi = std::min<uint64_t>(x + max_scanned_intervals, M_LF.num_intervals() - 1);
        if (i < p(hi + 1)) return exp_search_max_leq<RIGHT>(i, x, hi, p);
    } else {
        uint64_t lo = x > max_scanned_intervals ? x - max_scanned_intervals : 0;
        if (p(lo) <= i) return exp_search_max_leq<LEFT>(i, lo, x, p);
    }

    return mlf_interval<dir>(i);
}

template <move_r_support support>
template <query_support_t query_support, typename sym_t>
auto move_rb<support>::search(move_r_pattern_view<sym_t> P, bool& occurs) const
    -> search_context_t<query_support>
{
    uint64_t matched = 0;
    auto ctx = context_of<query_support>(P, matched);
    occurs = true;

    for (uint64_t i = P.size() - matched; i > 0; i--) {
        auto [next, ok] = ctx.template extend<ctx_primary_dir>(P[i - 1]);
        if (!ok) { occurs = false; return ctx; }
        ctx = std::move(next);
    }

    return ctx;
}

template <move_r_support support>
template <query_support_t query_support, typename sym_t>
auto move_rb<support>::context_of(move_r_pattern_view<sym_t> P, uint64_t& matched) const
    -> search_context_t<query_support>
{
    matched = 0;
    auto ctx = empty_context<query_support>();
    const uint64_t len = _context_tables.single_length();
    if (len == 0 || P.size() < len) return ctx;

    // look up the keys of P[first..first+len-1] in the forward and in the reverse order; the first resp.
    // the last symbol is the most significant symbol of the key (see context_tables::packed_key)
    const uint64_t first = P.size() - len;
    uint64_t key = _context_tables.packed_key(
        [&](uint64_t i) { return map_symbol(P[first + i]); }, len);
    uint32_t entry = _context_tables.template lookup<_forward_order>(key);
    if (entry == ::context_tables::no_entry) return ctx;
    uint32_t entry_R = ::context_tables::no_entry;

    uint64_t key_R = _context_tables.packed_key(
        [&](uint64_t i) { return map_symbol(P[first + len - 1 - i]); }, len);
    entry_R = _context_tables.template lookup<_reverse_order>(key_R);

    load_context<query_support>(ctx, len, entry, entry_R);
    matched = len;
    return ctx;
}

template <move_r_support support>
template <typename sym_t>
uint64_t move_rb<support>::count_hamming_dist(
    move_r_pattern_view<sym_t> P, const search_scheme_t& scheme, apm_params_t params) const
    {
        std::vector<uint8_t> P_(P.begin(), P.end());

        return apm_hamming{*this, params}.count(P_, scheme);
    }

template <move_r_support support>
template <typename sym_t>
std::vector<uint64_t> move_rb<support>::count_hamming_dist_histogram(
    move_r_pattern_view<sym_t> P, const search_scheme_t& scheme, apm_params_t params) const
    {
        std::vector<uint8_t> P_(P.begin(), P.end());

        return apm_hamming{*this, params}.count_histogram(P_, scheme);
    }

template <move_r_support support>
template <cigar_mode_t mode>
std::vector<cigar_t<uint8_t>> move_rb<support>::locate_apm(
    const std::vector<uint8_t>& P, const search_scheme_t& scheme, distance_metric_t dist_metr,
    std::vector<aprx_occ_t<mode>>& occ, apm_params_t params) const
    requires(supports_multiple_locate)
    {
        assert(dist_metr == HAMMING_DISTANCE || dist_metr == EDIT_DISTANCE);

        auto collect = [&occ](aprx_occ_t<mode> o) { occ.emplace_back(o); };

        if (dist_metr == HAMMING_DISTANCE)
            return apm_hamming<move_rb>{*this, params}.template locate<mode>(P, scheme, collect);

        return apm_edit<move_rb>{*this, params}.template locate<mode>(P, scheme, collect);
    }

#endif // MOVE_R_MOVE_RB_QUERIES_TPP
