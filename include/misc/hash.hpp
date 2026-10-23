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

#ifndef MOVE_R_MISC_HASH_HPP
#define MOVE_R_MISC_HASH_HPP

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <type_traits>
#include <misc/hashing.hpp>
#include <misc/utils.hpp>

/**
 * @brief hash functor for a value that converts to uint128_t (e.g. packed_uint_t, see move_r/construction/input_view.hpp);
 *        usable as a hash-map key hasher
 */
struct wide_hash {
    /**
     * @brief returns hash_64 of the XOR of the lower and the upper 64 bits of a value
     * @tparam val_t a type converting to uint128_t
     * @param value the value
     */
    template <typename val_t>
    uint64_t operator()(const val_t& value) const
    {
        uint128_t x = value;
        return hash_64(uint64_t(x) ^ uint64_t(x >> 64));
    }
};

/**
 * @brief combines the hash value hash with the value x (boost-style hash_combine)
 * @param hash the hash value to update
 * @param x the value to combine into hash
 */
inline static void hash_combine(uint64_t& hash, uint64_t x)
{
    hash ^= x + golden_ratio_64 + (hash << 6) + (hash >> 2);
}

/**
 * @brief hash functor for a std::span, hashing its elements (FNV-1a); usable as a hash-map key hasher
 */
struct span_hash {
    template <typename sym_t>
    uint64_t operator()(std::span<const sym_t> s) const
    {
        constexpr uint64_t fnv_offset_basis = 1469598103934665603ULL;
        constexpr uint64_t fnv_prime = 1099511628211ULL;
        uint64_t h = fnv_offset_basis;
        for (sym_t c : s) h = (h ^ (uint64_t) (std::make_unsigned_t<sym_t>) c) * fnv_prime;
        return h;
    }
};

/**
 * @brief equality functor comparing two std::spans element-wise; the matching key equality for span_hash
 */
struct span_eq {
    template <typename sym_t>
    bool operator()(std::span<const sym_t> a, std::span<const sym_t> b) const
    {
        return std::equal(a.begin(), a.end(), b.begin(), b.end());
    }
};

#endif // MOVE_R_MISC_HASH_HPP
