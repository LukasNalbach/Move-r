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

#ifndef MOVE_R_MISC_APM_HPP
#define MOVE_R_MISC_APM_HPP

#include <array>
#include <cstdint>
#include <cstring>
#include <ranges>
#include <span>
#include <string>
#include <type_traits>
#include <variant>
#include <vector>

#include <tsl/sparse_set.h>

#include <misc/utils.hpp>
#include <misc/portable_layout.hpp>
#include <misc/hash.hpp>
#include <misc/log.hpp>
#include <misc/files.hpp>
#include <misc/search.hpp>
#include <misc/search_schemes.hpp>
#include <data_structures/edit_distance_matrix.hpp>
#include <misc/cigar.hpp>

/** @brief whether an approximate-pattern-matching query only counts occurrences or also locates them */
enum query_support_t : uint8_t {
    COUNT = 0,
    LOCATE = 1
};

/**
 * @brief whether an enumeration of the right extensions of a pattern reports the terminator symbol; a left
 *        extension with the terminator is never reported
 */
enum extension_scope_t : uint8_t {
    NO_TERM = 0, // all symbols of the text except the terminator
    TERM = 1     // all symbols of the text and the terminator (as the smallest symbol, reported as TERMINATOR_SYMBOL)
};

// bit-field widths of a packed aprx_occ_t: position, length and error count (10 bytes together)
static constexpr uint32_t aprx_occ_pos_bits = 48; // positions up to 2^48 - 1 (see uint48_t)
static constexpr uint32_t aprx_occ_len_bits = 24; // occurrence lengths up to 2^24 - 1 (i.e. m + k)
static constexpr uint32_t aprx_occ_err_bits = 8;  // error counts up to 2^8 - 1 (i.e. k)

/**
 * @brief an approximate occurrence: its starting position, length and error count in bit-fields; with a CIGAR, also
 *        the index cig_idx of its CIGAR in the vector of CIGARs returned by locate (one CIGAR per SA interval, which
 *        all occurrences of the SA interval share); 10 bytes without and 14 bytes with a CIGAR index. The bit-fields
 *        have the layout of a packed struct of GCC (see packed_bit_fields).
 * @tparam mode whether the occurrence also carries a CIGAR alignment
 */
template <cigar_mode_t mode = NO_CIGAR>
struct __attribute__((packed)) aprx_occ_t {
    using pos_field = next_bit_field_t<void, aprx_occ_pos_bits>;
    using len_field = next_bit_field_t<pos_field, aprx_occ_len_bits>;
    using err_field = next_bit_field_t<len_field, aprx_occ_err_bits>;

    packed_bit_fields<bytes_of_bit_fields<err_field>> bits; // pos, len and err
    STRING_UTILS_NO_UNIQUE_ADDRESS std::conditional_t<mode == CIGAR, uint32_t, std::monostate> cig_idx = {};

    aprx_occ_t() = default;

    aprx_occ_t(uint64_t pos, uint64_t len, uint64_t err)
    { bits.template assign<pos_field, len_field, err_field>(pos, len, err); }

    aprx_occ_t(uint64_t pos, uint64_t len, uint64_t err, uint32_t cigar_idx) requires(mode == CIGAR)
        : aprx_occ_t(pos, len, err)
    { cig_idx = cigar_idx; }

    uint64_t pos() const { return bits.template get<pos_field>(); } // the starting position
    uint64_t len() const { return bits.template get<len_field>(); } // the length
    uint64_t err() const { return bits.template get<err_field>(); } // the error count

    void set_pos(uint64_t value) { bits.template set<pos_field>(value); }
    void set_len(uint64_t value) { bits.template set<len_field>(value); }
    void set_err(uint64_t value) { bits.template set<err_field>(value); }

    /** @brief orders occurrences by position, then by error count, then by length (the CIGAR is ignored) */
    bool operator<(const aprx_occ_t& other) const
    {
        return pos() != other.pos() ? pos() < other.pos() :
              (err() != other.err() ? err() < other.err() :
                                      len() < other.len());
    }

    bool operator==(const aprx_occ_t& other) const
    {
        return pos() == other.pos() &&
               len() == other.len() &&
               err() == other.err();
    }

    /** @brief writes an occurrence to an output stream as "(pos, len, err)" */
    friend std::ostream& operator<<(std::ostream& os, const aprx_occ_t& occ) {
        return os << "(" << occ.pos() << ", " << occ.len() << ", " << occ.err() << ")";
    }
};

// IntelliSense (EDG) lays out empty [[no_unique_address]] members of packed structs differently than GCC and clang, so
// the size check is only done by the compiler
#ifndef __INTELLISENSE__
static_assert(sizeof(aprx_occ_t<NO_CIGAR>) == (aprx_occ_pos_bits + aprx_occ_len_bits + aprx_occ_err_bits) / 8);
#endif
static_assert(sizeof(aprx_occ_t<CIGAR>) == sizeof(aprx_occ_t<NO_CIGAR>) + sizeof(uint32_t));

// key function for radix_sort: returns the starting position of an occurrence (see aprx_occ_t)
struct key_of_pos {
    template <typename occ_t>
    inline auto operator()(const occ_t& occurrence) const { return occurrence.pos(); }
};

/**
 * @brief computes the hamming distance of two equally long sequences, stopping early once it exceeds k
 * @tparam inp_t contiguous range type of the sequences (e.g. std::string, std::vector, std::span)
 * @param str_1 the first sequence
 * @param str_2 the second sequence (must have the same length as str_1)
 * @param k the distance threshold
 * @return the hamming distance of str_1 and str_2 (or a value > k if it exceeds k)
 */
template <std::ranges::contiguous_range inp_t>
static uint64_t hamming_dist_bounded(const inp_t& str_1, const inp_t& str_2, uint64_t k)
{
    using sym_t = std::ranges::range_value_t<inp_t>;
    assert(str_1.size() == str_2.size());
    const uint64_t n = str_1.size();
    uint64_t dist = 0;

    /* for a byte alphabet, compare eight symbols per step: XOR two 8-byte words, then count the mismatching (nonzero)
       bytes by OR-reducing the bits of each byte into its low bit and popcounting */
    if constexpr (sizeof(sym_t) == 1) {
        const auto* b1 = reinterpret_cast<const uint8_t*>(std::ranges::data(str_1));
        const auto* b2 = reinterpret_cast<const uint8_t*>(std::ranges::data(str_2));
        uint64_t i = 0;

        for (; i + 8 <= n && dist <= k; i += 8) {
            uint64_t w1, w2;
            std::memcpy(&w1, b1 + i, 8);
            std::memcpy(&w2, b2 + i, 8);
            uint64_t x = w1 ^ w2;
            x |= x >> 1;
            x |= x >> 2;
            x |= x >> 4;
            x &= 0x0101010101010101ULL; // low bit of each byte set iff that byte differed
            dist += ::popcount(x);
        }

        for (; i < n && dist <= k; i++) dist += b1[i] != b2[i];
    } else {
        const auto* p1 = std::ranges::data(str_1);
        const auto* p2 = std::ranges::data(str_2);

        for (uint64_t i = 0; i < n && dist <= k; i++) dist += p1[i] != p2[i];
    }

    return dist;
}

/**
 * @brief naively locates all occurrences of P in T with at most k mismatches (w.r.t. hamming distance)
 * @tparam inp_t contiguous range type of the text and pattern
 * @param T the text
 * @param P the pattern
 * @param k the maximum number of mismatches
 * @return the occurrences of P in T with at most k mismatches
 */
template <std::ranges::contiguous_range inp_t>
static std::vector<aprx_occ_t<>> locate_hamming_dist(const inp_t& T, const inp_t& P, uint64_t k)
{
    using sym_t = std::ranges::range_value_t<inp_t>;
    std::span<const sym_t> Ts(T), Ps(P); // read-only views for the windowed distance calls
    uint64_t m = Ps.size();
    uint64_t n = Ts.size();

    std::vector<aprx_occ_t<>> Occ;

    for (uint64_t i = 0; i + m <= n; i++) {
        uint64_t dist = hamming_dist_bounded(Ts.subspan(i, m), Ps, k);
        if (dist <= k) Occ.emplace_back(aprx_occ_t<>(i, m, dist));
    }

    return Occ;
}

/**
 * @brief computes the edit distance of two sequences with a banded bit-parallel dynamic program (Myers), capped at
 *        k + 1. Requires k <= edit_distance_max_errors_128 and the longer sequence to be at most edit_distance_max_length long.
 * @tparam inp_t contiguous range type of the sequences (a byte alphabet, since the matrix indexes symbols by byte)
 * @param T the first sequence
 * @param P the second sequence
 * @param k the distance threshold
 * @return the edit distance of T and P, or k + 1 if it exceeds k
 */
template <std::ranges::contiguous_range inp_t>
static int64_t edit_dist_bounded(const inp_t& T, const inp_t& P, int64_t k)
{
    const int64_t nT = T.size(), nP = P.size();
    if (std::abs(nT - nP) > k) return k + 1;
    const bool p_longer = nP > nT;
    const int64_t nR = p_longer ? nP : nT, nQ = p_longer ? nT : nP;
    if (nQ == 0) return nR;
    assert(k <= edit_distance_max_errors_128);
    assert(nR <= edit_distance_max_length);

    auto run = [&]<typename word_t>() -> int64_t {
        edit_distance_matrix<word_t> mat;
        if (p_longer) setup_edit_matrix<word_t>(mat, P, nP, T, nT, k);
        else          setup_edit_matrix<word_t>(mat, T, nT, P, nP, k);
        return std::min<int64_t>(mat(nQ, nR), k + 1);
    };

    return k <= edit_distance_max_errors_64 ?
        run.template operator()<uint64_t>() :
        run.template operator()<uint128_t>();
}

/**
 * @brief over all prefixes M[0..len) of M (len in the band [m-k, m+k]), finds the one minimizing ed(P, M[0..len))
 * @tparam word_t word type of the edit-distance matrix (chosen by the caller to fit k)
 * @tparam inp_t contiguous range type of the matched string (a byte alphabet, since the matrix indexes symbols by byte)
 * @param mat the banded edit-distance matrix, with the pattern P already set as its input (columns)
 * @param M the matched string whose best pattern-length prefix is sought
 * @param m the pattern length (= mat's number of columns minus one)
 * @param k the error threshold
 * @return (minimum error, shortest length achieving it), or (k + 1, 0) if no prefix stays within k errors
 */
template <typename word_t, std::ranges::contiguous_range inp_t>
std::pair<uint64_t, uint64_t> edit_dist_prefix(
    edit_distance_matrix<word_t>& mat, const inp_t& M, uint64_t m, uint64_t k)
{
    mat.init(k);
    const uint64_t n2 = M.size();
    const uint64_t l_min = m > k ? m - k : 0;
    const uint64_t i_max = std::min<uint64_t>(n2, m + k);
    const uint16_t final_col = mat.num_cols() - 1;
    uint64_t best_err = k + 1, best_len = 0;

    for (uint64_t i = 1; i <= i_max; i++) {
        bool alive = mat.compute_row(i, M[i - 1]);
        if (i >= l_min) {
            uint64_t err = mat(i, final_col);
            if (err < best_err) { best_err = err; best_len = i; }
        }
        if (!alive) break;
    }

    return {best_err, best_len};
}

/**
 * @brief locates all occurrences of P in T with at most k errors (w.r.t. edit distance) by running the banded
 *        bit-parallel matrix (Myers; columns = P, set up once) over every starting position of T, keeping for
 *        each starting position the shortest minimum-error occurrence
 * @tparam inp_t contiguous range type of the text and pattern (a byte alphabet, since the matrix indexes symbols by byte)
 * @param T the text
 * @param P the pattern (1 <= |P| <= edit_distance_max_length)
 * @param k the maximum number of errors (k <= edit_distance_max_errors_128)
 * @return the occurrences of P in T with at most k errors
 */
template <std::ranges::contiguous_range inp_t>
static std::vector<aprx_occ_t<>> locate_edit_dist(const inp_t& T, const inp_t& P, uint64_t k)
{
    static_assert(sizeof(std::ranges::range_value_t<inp_t>) == 1);
    const int64_t n = T.size(), m = P.size();
    assert(m > 0 && m <= edit_distance_max_length && int64_t(k) <= edit_distance_max_errors_128);
    const int64_t l_min = std::max<int64_t>(1, m - int64_t(k)); // occurrences have positive length (no empty window)
    const int64_t l_max = m + k;

    auto run = [&]<typename word_t>() {
        edit_distance_matrix<word_t> mat;
        std::vector<uint8_t> alph_map;
        uint16_t sigma = build_alphabet_map(P, m, alph_map);
        directional_substring<inp_t> P_col(P, 0, m - 1, RIGHT);
        mat.set_input(P_col, sigma, alph_map);

        std::vector<aprx_occ_t<>> Occ;

        for (int64_t i = 0; i < n; i++) {
            mat.init(k);
            aprx_occ_t<> occ(uint64_t(i), 0, uint64_t(k + 1));
            const int64_t j_max = std::min(l_max, n - i);

            for (int64_t j = 1; j <= j_max; j++) {
                // once the whole band exceeds k, no longer window can come back under the budget
                if (!mat.compute_row(j, T[i + j - 1])) [[unlikely]] break;

                // a final-column cell holds ed(T[i..i+j), P) and is exact whenever it is within the budget
                if (j >= l_min && mat.is_in_final_column(j)) {
                    int64_t err = mat(j, m);

                    if (err <= int64_t(k) && (err < int64_t(occ.err()) ||
                       (err == int64_t(occ.err()) && j < int64_t(occ.len())))
                    ) {
                        occ.set_err(err);
                        occ.set_len(j);
                    }
                }
            }

            if (occ.err() <= k) Occ.emplace_back(occ);
        }

        return Occ;
    };

    return int64_t(k) <= edit_distance_max_errors_64 ?
        run.template operator()<uint64_t>() :
        run.template operator()<uint128_t>();
}

/**
 * @brief removes redundant approximate occurrences in place: an occurrence is dropped if the previous kept occurrence
 *        and the next occurrence both have an error count at most its own and their positions differ by at most
 *        2k + 3
 * @tparam mode whether the occurrences carry a CIGAR alignment (moved along with the kept occurrences)
 * @param Occ the occurrences (sorted by position); filtered in place
 * @param k the maximum number of errors
 */
template <cigar_mode_t mode = NO_CIGAR>
static void filter_edit_distance_occurrences(std::vector<aprx_occ_t<mode>>& Occ, uint64_t k)
{
    if (Occ.size() <= 2) return;
    uint64_t write_idx = 0;

    for (uint64_t read_idx = 0; read_idx < Occ.size(); read_idx++) {
        auto occ = std::move(Occ[read_idx]);

        // drop the last kept occurrence (mid) while it lies between two occurrences that cover it
        while (write_idx >= 2) {
            const auto& mid = Occ[write_idx - 1];
            const auto& prev = Occ[write_idx - 2];

            if (occ.pos() - prev.pos() <= 2 * k + 3 &&
                prev.err() <= mid.err() &&
                occ.err() <= mid.err()
            ) {
                write_idx--;
            } else {
                break;
            }
        }

        Occ[write_idx++] = std::move(occ);
    }

    Occ.resize(write_idx);
}

/**
 * @brief verifies that every occurrence in occ_all is covered by an occurrence in occ_filtered with an equal or
 *        smaller error count within 2k+1 positions (i.e. filtering did not drop any required occurrence)
 * @param occ_filtered the filtered occurrences (sorted by position)
 * @param occ_all all occurrences (sorted by position)
 * @param k the maximum number of errors
 * @return whether occ_filtered covers all occurrences in occ_all
 */
inline bool verify_edit_distance_coverage(
    const std::vector<aprx_occ_t<>>& occ_filtered,
    const std::vector<aprx_occ_t<>>& occ_all,
    uint64_t k
) {
    uint64_t max_dist = 2 * k + 1;
    uint64_t left = 0, right = 0;

    for (const auto& o : occ_all) {
        while (right < occ_filtered.size() && occ_filtered[right].pos() <= o.pos() + max_dist) right++;
        while (left < right && occ_filtered[left].pos() + max_dist < o.pos()) left++;
        if (left == right) return false;
        bool valid = false;

        for (uint64_t i = left; i < right; i++) {
            if (occ_filtered[i].err() <= o.err()) {
                valid = true;
                break;
            }
        }

        if (!valid) return false;
    }

    return true;
}

/**
 * @brief keeps a subset of the occurrences such that every occurrence has a kept one with an equal or smaller error
 *        count within max_dist positions: the error counts are processed in increasing order, and the leftmost
 *        occurrence of the current error count that is not covered yet is covered by keeping the rightmost occurrence
 *        at most max_dist positions after it whose error count is at most its own
 * @tparam mode whether the occurrences carry a CIGAR alignment (moved along with the kept occurrences)
 * @param Occ the occurrences (sorted by position); filtered in place
 * @param max_dist the maximum distance between an occurrence and a kept occurrence covering it
 */
template <cigar_mode_t mode = NO_CIGAR>
static void select_covering_occurrences(std::vector<aprx_occ_t<mode>>& Occ, uint64_t max_dist)
{
    enum : uint8_t { UNCOVERED, COVERED, KEPT };
    constexpr uint64_t num_errs = uint64_t(1) << aprx_occ_err_bits;
    const uint64_t n = Occ.size();
    std::vector<uint8_t> state(n, UNCOVERED);

    // by_err: the occurrences in increasing order of their error counts, then of their positions; the ones with error
    // count err are by_err[level_beg[err], level_beg[err + 1])
    std::array<uint64_t, num_errs + 1> level_beg {};
    for (const auto& occ : Occ) level_beg[occ.err() + 1]++;
    for (uint64_t err = 0; err < num_errs; err++) level_beg[err + 1] += level_beg[err];
    std::array<uint64_t, num_errs + 1> level_end = level_beg;
    std::vector<uint64_t> by_err(n);
    for (uint64_t i = 0; i < n; i++) by_err[level_end[Occ[i].err()]++] = i;

    for (uint64_t err = 0; err < num_errs; err++) {
        uint64_t hi = 0;

        for (uint64_t x = level_beg[err]; x < level_beg[err + 1]; x++) {
            const uint64_t i = by_err[x];
            if (state[i] != UNCOVERED) continue;
            hi = std::max(hi, i);
            while (hi + 1 < n && Occ[hi + 1].pos() <= Occ[i].pos() + max_dist) hi++;
            uint64_t j = hi;
            while (Occ[j].err() > err) j--;
            state[j] = KEPT;

            for (uint64_t t = j; t-- > 0 && Occ[j].pos() - Occ[t].pos() <= max_dist;) {
                if (state[t] == UNCOVERED && Occ[t].err() >= Occ[j].err()) state[t] = COVERED;
            }

            for (uint64_t t = j + 1; t < n && Occ[t].pos() - Occ[j].pos() <= max_dist; t++) {
                if (state[t] == UNCOVERED && Occ[t].err() >= Occ[j].err()) state[t] = COVERED;
            }
        }
    }

    uint64_t write_idx = 0;
    for (uint64_t i = 0; i < n; i++) if (state[i] == KEPT) Occ[write_idx++] = std::move(Occ[i]);
    Occ.resize(write_idx);
}

/** @brief the distance metric used for approximate pattern matching */
enum distance_metric_t : int8_t {
    HAMMING_DISTANCE = 0,
    EDIT_DISTANCE = 1
};

/**
 * @brief locates all occurrences of P in T with at most k errors w.r.t. the given distance metric
 * @tparam inp_t contiguous range type of the text and pattern
 * @param T the text
 * @param P the pattern
 * @param k the maximum number of errors
 * @param metric the distance metric
 * @return the occurrences of P in T with at most k errors
 */
template <std::ranges::contiguous_range inp_t>
static std::vector<aprx_occ_t<>> locate_naive(const inp_t& T, const inp_t& P, uint64_t k,
    distance_metric_t metric)
{
    assert(metric == HAMMING_DISTANCE || metric == EDIT_DISTANCE);
    if (metric == HAMMING_DISTANCE) return locate_hamming_dist(T, P, k);
    return locate_edit_dist(T, P, k);
}

// ############################# PARTITIONING #############################

/** @brief strategy for partitioning a pattern into the parts of a search scheme */
enum partitioning_t : uint8_t {
    _partition_uniform = 0, // p parts of equal length
    _partition_dynamic = 1  // part boundaries chosen per pattern, s.t. each part has few exact occurrences
};

class apm_costs_t;

/** @brief the options of an approximate-pattern-matching query that are not part of the search scheme */
struct apm_params_t {
    // strategy for partitioning the pattern
    partitioning_t partitioning = _partition_dynamic;

    /* true <=> a search starts with the search context of the exact match of its first part computed by the
       partitioning; false <=> it matches its first part symbol by symbol */
    bool start_at_part_contexts = true;

    /* the error model compiled for the pattern (see apm_costs_t); nullptr: every error costs 1 and the bound is the k
       of the search scheme. With a model, the k of the search scheme must be at least costs->max_errors(), and an
       occurrence is reported with its cost as its error count, in increasing order of the positions; for edit distance,
       every alignment within the bound that starts at a position i has a reported occurrence at most
       2 * costs->max_errors() + 1 positions away from i whose cost is at most its own */
    const apm_costs_t* costs = nullptr;
};

/**
 * @brief locates the occurrences of a pattern with at most best + x errors, where best is the minimum number
 *        of errors of an occurrence with at most max_k errors (all-best mapping)
 *
 * locate_k is called with increasing error bounds k, starting with k = 0, until it reports an occurrence;
 * k is increased by x + 2 for k < 5 and by x + 4 for k >= 5 (the steps of Columba, Renders et al.). If then
 * best + x > k, locate_k is called once more with k = min(best + x, max_k).
 *
 * @tparam mode whether the occurrences carry a CIGAR alignment
 * @tparam locate_fnc_t type of locate_k
 * @param max_k the maximum number of errors
 * @param x the number of error counts above best that are also reported
 * @param occ (output) the occurrences; the vector is cleared before each call of locate_k
 * @param locate_k function called as locate_k(k, occ), which appends all occurrences with at most k errors to occ
 * @return best, or max_k + 1 if there is no occurrence with at most max_k errors
 */
template <cigar_mode_t mode = NO_CIGAR, typename locate_fnc_t>
uint64_t locate_best(uint64_t max_k, uint64_t x, std::vector<aprx_occ_t<mode>>& occ, locate_fnc_t locate_k)
{
    for (uint64_t k = 0;;) {
        occ.clear();
        locate_k(k, occ);

        if (!occ.empty()) {
            uint64_t best = occ.front().err();
            for (const aprx_occ_t<mode>& o : occ) best = std::min<uint64_t>(best, o.err());

            // the occurrences with more than k and at most best + x errors have not been located yet
            if (best + x > k && k < max_k) { k = std::min(best + x, max_k); continue; }

            std::erase_if(occ, [&](const aprx_occ_t<mode>& o) { return o.err() > best + x; });
            return best;
        }

        if (k >= max_k) return max_k + 1;
        k = std::min(max_k, k + x + (k < 5 ? 2 : 4));
    }
}

#endif // MOVE_R_MISC_APM_HPP
