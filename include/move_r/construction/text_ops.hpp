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
 * @file text_ops.hpp
 * @brief functions that read or write all symbols of the input, implemented for each symbol width
 */

#ifndef MOVE_R_MOVE_R_CONSTRUCTION_TEXT_OPS_HPP
#define MOVE_R_MOVE_R_CONSTRUCTION_TEXT_OPS_HPP

#include <data_structures/elias_fano_bit_vector.hpp>
#include <data_structures/interleaved_bit_aligned_vectors.hpp>
#include <libsais40.hpp>
#include <misc/sort.hpp>
#include <move_r/construction/input_view.hpp>

#include <array>
#include <cstdint>
#include <vector>

/**
 * @brief replaces each effective symbol in T[0..n-2] with its input symbol
 * @tparam sym_t symbol type of the input
 * @param T the input, remapped to its effective alphabet
 * @param n the length of T (including the terminator)
 * @param p the number of threads to use
 * @param alphabet marks the distinct symbols of the input (see move_r's _alphabet)
 * @param order_inv order_inv[c - 1] = rank of the input symbol of the effective symbol c among the distinct
 *        symbols of the input; empty, if c is that rank (see move_r's _alphabet_order)
 */
template <typename sym_t>
void move_r_unmap_text_of(move_r_input_view& T, uint64_t n, uint64_t p,
    const elias_fano_bit_vector& alphabet, const interleaved_bit_aligned_vectors& order_inv)
{
    sym_t* text = reinterpret_cast<sym_t*>(T.data());
    const uint64_t num_ones = alphabet.num_ones();
    const bool reordered = !order_inv.empty();

    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < n - 1; i++) {
        uint64_t sym = uint64_t(text[i]);

        if (sym == 0) {
            text[i] = sym_t(0);
            continue;
        }

        uint64_t rank = reordered ? order_inv.template get<0, uint64_t>(sym - 1) : sym;
        text[i] = sym_t(rank <= num_ones
            ? alphabet.select_1(rank) : std::numeric_limits<uint64_t>::max());
    }
}

/**
 * @brief replaces each symbol c in T[0..n-2] with symbols[c - 1]
 * @tparam sym_t symbol type of the input
 * @param T the input, remapped to its effective alphabet
 * @param n the length of T (including the terminator)
 * @param p the number of threads to use
 * @param symbols [0..sigma-2] the distinct symbols of the input in increasing order
 */
template <typename sym_t>
void move_r_unmap_text_by_symbols_of(move_r_input_view& T, uint64_t n, uint64_t p,
    const std::vector<uint64_t>& symbols)
{
    sym_t* text = reinterpret_cast<sym_t*>(T.data());

    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < n - 1; i++) {
        text[i] = sym_t(symbols[uint64_t(text[i]) - 1]);
    }
}

// compares two symbols by decreasing frequency, symbols with equal frequency by increasing value
struct by_decreasing_frequency {
    const uint64_t* freq;

    inline bool operator()(uint64_t a, uint64_t b) const
    {
        return freq[a] != freq[b] ? freq[a] > freq[b] : a < b;
    }
};

#define MOVE_R_ORDER_ALPHABET_SORTS(prefix) \
    prefix template void sort_parallel<std::vector<uint64_t>::iterator, by_decreasing_frequency>( \
        std::vector<uint64_t>::iterator, std::vector<uint64_t>::iterator, by_decreasing_frequency, uint64_t); \
    prefix template void sort_sequential<std::vector<uint64_t>::iterator, by_decreasing_frequency>( \
        std::vector<uint64_t>::iterator, std::vector<uint64_t>::iterator, by_decreasing_frequency);

/**
 * @brief computes the effective symbols of the distinct symbols of an integer-alphabet input; the effective
 *        alphabet is ordered by decreasing frequency (symbols with equal frequency by increasing value)
 * @param freq [0..sigma-2] freq[i] = number of occurrences of the i-th smallest distinct symbol
 * @param reorder whether to reorder the alphabet (see move_r_params::reorder_alphabet)
 * @param effective [0..sigma-2] effective[i] stores the effective symbol of the i-th smallest distinct symbol;
 *        left empty, if the alphabet is not reordered (the effective symbol is then i + 1)
 * @param p the number of threads to use
 */
inline void move_r_order_alphabet(const std::vector<uint64_t>& freq, bool reorder,
    std::vector<uint64_t>& effective, uint64_t p)
{
    effective.clear();
    if (!reorder || freq.empty()) return;

    // [0..sigma-2] the indices of the distinct symbols, sorted by decreasing frequency (and by increasing index)
    no_init_vector<uint64_t> by_freq;
    resize_no_init(by_freq, freq.size());
    const uint64_t max_freq = *std::max_element(freq.begin(), freq.end());
    const uint8_t index_bits = std::bit_width(freq.size() - 1);

    if (std::bit_width(max_freq) + index_bits <= 64) {
        // the key (max_freq - freq[i]) * 2^index_bits + i orders by decreasing frequency, then by increasing index
        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < by_freq.size(); i++) by_freq[i] = ((max_freq - freq[i]) << index_bits) | i;

        radix_sort(by_freq.begin(), by_freq.end(), key_of_value {}, p);
        const uint64_t index_mask = (uint64_t { 1 } << index_bits) - 1;

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < by_freq.size(); i++) by_freq[i] &= index_mask;
    } else {
        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < by_freq.size(); i++) by_freq[i] = i;

        comparison_sort(by_freq.begin(), by_freq.end(), by_decreasing_frequency { freq.data() }, p);
    }

    resize_for_overwrite(effective, freq.size());

    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < by_freq.size(); i++) effective[by_freq[i]] = i + 1;
}

/**
 * @brief builds the suffix array of the input (see libsais40)
 * @tparam sym_t symbol type of the input
 * @param T the input
 * @param SA [0..n-1+fs] stores the suffix array
 * @param n the length of T
 * @param sigma the alphabet size for an integer alphabet, 0 for a byte alphabet
 * @param fs the extra space (in entries) after SA
 * @param p the number of threads to use
 */
template <typename sym_t>
void move_r_build_sa_of(move_r_input_view& T, uint40_t* SA, uint64_t n, uint64_t sigma,
    uint64_t fs, uint64_t p)
{
    libsais40(reinterpret_cast<const sym_t*>(T.data()), SA, n, sigma, fs, p);
}

/**
 * @brief builds the RLBWT from the suffix array
 * @tparam sym_t symbol type of the input
 * @param T the input, remapped to its effective alphabet
 * @param sa [0..n-1] the suffix array
 * @param n_p [0..p'] n_p[i_p] = starting position of thread i_p's section in L
 * @param r_p [0..p'] r_p[i_p] stores the number of runs starting in thread i_p's section
 * @param rlbwt [0..p'-1] rlbwt[i_p] stores the runs in thread i_p's section
 * @param C [0..p'-1][0..255] C[i_p][c] stores the number of occurrences of c in thread i_p's section (only
 *        for a byte alphabet)
 * @param byte_alphabet true <=> the input uses a byte alphabet
 * @param p_ p', the number of threads to use
 * @param prog progress of the phase, advanced by the number of processed positions
 */
template <typename sym_t>
void move_r_bwt_runs_from_sa_of(move_r_input_view& T, const no_init_vector<uint40_t>& sa,
    const std::vector<uint64_t>& n_p, std::vector<uint64_t>& r_p,
    std::vector<interleaved_bit_aligned_vectors>& rlbwt, std::vector<std::vector<uint64_t>>& C,
    bool byte_alphabet, uint64_t p_, phase_progress& prog)
{
    const sym_t* text = reinterpret_cast<const sym_t*>(T.data());

    auto add_run = [&](uint64_t i_p, uint64_t sym, uint64_t len) {
        std::array<uint64_t, interleaved_bit_aligned_vectors::num_vectors> values {};
        values[0] = sym;
        values[1] = len;
        rlbwt[i_p].emplace_back(values);
    };

    #pragma omp parallel num_threads(p_)
    {
        uint64_t i_p = omp_get_thread_num();
        // iteration range L[b..e-1] of thread i_p
        uint64_t b = n_p[i_p];
        uint64_t e = n_p[i_p + 1];
        uint64_t i_ = b; // start position of the last seen run in L
        uint64_t prev_sym = sa[b] == 0 ? 0 : uint64_t(text[uint64_t(sa[b]) - 1]);
        uint64_t cur_sym;

        // the progress is advanced after every grain positions
        const uint64_t grain = prog.grain();

        for (uint64_t base = b + 1; base < e; base += grain) {
            const uint64_t stop = std::min(e, base + grain);

            for (uint64_t i = base; i < stop; i++) {
                cur_sym = sa[i] == 0 ? 0 : uint64_t(text[uint64_t(sa[i]) - 1]);

                if (cur_sym != prev_sym) {
                    add_run(i_p, prev_sym, i - i_);
                    if (byte_alphabet) C[i_p][prev_sym] = C[i_p][prev_sym] + (i - i_);
                    prev_sym = cur_sym;
                    i_ = i;
                }
            }

            prog.advance(stop - base);
        }

        // add the run L[i'..e)
        add_run(i_p, prev_sym, e - i_);
        if (byte_alphabet) C[i_p][prev_sym] = C[i_p][prev_sym] + (e - i_);
        r_p[i_p] = rlbwt[i_p].size();
        rlbwt[i_p].shrink_to_fit();
    }
}

/**
 * @brief reverses T[0..n-2] in place
 * @tparam sym_t symbol type of the input
 * @param T the input
 * @param n the length of T (including the terminator)
 * @param p the number of threads to use
 */
template <typename sym_t>
void move_r_reverse_text_of(move_r_input_view& T, uint64_t n, uint64_t p)
{
    sym_t* text = reinterpret_cast<sym_t*>(T.data());
    uint64_t n_2 = (n - 1) / 2;

    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < n_2; i++) {
        std::swap(text[i], text[(n - i) - 2]);
    }
}

// defines the function name, which calls name_of<sym_t> for the symbol type sym_t of the symbol width of T
#define MOVE_R_TEXT_OP_DISPATCH(name, signature, arguments) \
    inline void name signature \
    { \
        move_r_with_symbol_type(T.symbol_width(), [&](auto type) { \
            using sym_t = typename decltype(type)::type; \
            name##_of<sym_t> arguments; \
        }); \
    }

MOVE_R_TEXT_OP_DISPATCH(move_r_unmap_text,
    (move_r_input_view& T, uint64_t n, uint64_t p, const elias_fano_bit_vector& alphabet,
        const interleaved_bit_aligned_vectors& order_inv),
    (T, n, p, alphabet, order_inv))

MOVE_R_TEXT_OP_DISPATCH(move_r_unmap_text_by_symbols,
    (move_r_input_view& T, uint64_t n, uint64_t p, const std::vector<uint64_t>& symbols),
    (T, n, p, symbols))

MOVE_R_TEXT_OP_DISPATCH(move_r_reverse_text,
    (move_r_input_view& T, uint64_t n, uint64_t p),
    (T, n, p))

MOVE_R_TEXT_OP_DISPATCH(move_r_build_sa,
    (move_r_input_view& T, uint40_t* SA, uint64_t n, uint64_t sigma, uint64_t fs, uint64_t p),
    (T, SA, n, sigma, fs, p))

MOVE_R_TEXT_OP_DISPATCH(move_r_bwt_runs_from_sa,
    (move_r_input_view& T, const no_init_vector<uint40_t>& sa, const std::vector<uint64_t>& n_p,
        std::vector<uint64_t>& r_p, std::vector<interleaved_bit_aligned_vectors>& rlbwt,
        std::vector<std::vector<uint64_t>>& C, bool byte_alphabet, uint64_t p_,
        phase_progress& prog),
    (T, sa, n_p, r_p, rlbwt, C, byte_alphabet, p_, prog))

#undef MOVE_R_TEXT_OP_DISPATCH

// the editor parser of VS Code (__INTELLISENSE__) does not resolve the extern declarations of constrained member
// templates; it does not need them
#if !defined(STRING_UTILS_NO_EXTERN_TEMPLATES) && !defined(__INTELLISENSE__)

#define MOVE_R_DECLARE_TEXT_OPS(sym_t) \
    extern template void move_r_unmap_text_of<sym_t>( \
        move_r_input_view&, uint64_t, uint64_t, const elias_fano_bit_vector&, \
        const interleaved_bit_aligned_vectors&); \
    extern template void move_r_unmap_text_by_symbols_of<sym_t>(move_r_input_view&, uint64_t, \
        uint64_t, const std::vector<uint64_t>&); \
    extern template void move_r_reverse_text_of<sym_t>(move_r_input_view&, uint64_t, uint64_t); \
    extern template void move_r_build_sa_of<sym_t>( \
        move_r_input_view&, uint40_t*, uint64_t, uint64_t, uint64_t, uint64_t); \
    extern template void move_r_bwt_runs_from_sa_of<sym_t>(move_r_input_view&, \
        const no_init_vector<uint40_t>&, const std::vector<uint64_t>&, std::vector<uint64_t>&, \
        std::vector<interleaved_bit_aligned_vectors>&, std::vector<std::vector<uint64_t>>&, \
        bool, uint64_t, phase_progress&);

MOVE_R_FOR_EACH_INPUT_WIDTH(MOVE_R_DECLARE_TEXT_OPS)
MOVE_R_ORDER_ALPHABET_SORTS(extern)

#undef MOVE_R_DECLARE_TEXT_OPS

#endif

#endif // MOVE_R_MOVE_R_CONSTRUCTION_TEXT_OPS_HPP
