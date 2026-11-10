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
 * @file wide_word.hpp
 * @brief unsigned integers of 256, 512 or 1024 bits: the rows of an edit-distance matrix for more errors than are
 *        possible with uint128_t rows
 */

#ifndef MOVE_R_DATA_STRUCTURES_WIDE_WORD_HPP
#define MOVE_R_DATA_STRUCTURES_WIDE_WORD_HPP

#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <type_traits>

#if defined(__AVX2__)
#include <immintrin.h>
#endif

#define WIDE_WORD_INLINE __attribute__((always_inline)) inline

#if defined(__AVX2__)

/**
 * @brief unsigned integer of num_bits bits with the operators of a built-in unsigned integer, in AVX2 registers
 * @tparam num_bits number of bits, a multiple of 256
 */
template <uint16_t num_bits>
class wide_word {
    static_assert(num_bits % 256 == 0);

    static constexpr uint16_t num_vecs = num_bits / 256;
    static constexpr uint16_t num_limbs = num_bits / 64;

    __m256i vecs[num_vecs]; // the 256-bit parts, the least significant first

    // writes the 64-bit parts to limbs, the least significant first
    WIDE_WORD_INLINE void store(uint64_t* limbs) const
    {
        for (uint16_t i = 0; i < num_vecs; i++) _mm256_store_si256(reinterpret_cast<__m256i*>(limbs) + i, vecs[i]);
    }

    WIDE_WORD_INLINE void load(const uint64_t* limbs)
    {
        for (uint16_t i = 0; i < num_vecs; i++) vecs[i] = _mm256_load_si256(reinterpret_cast<const __m256i*>(limbs) + i);
    }

    // shifts by shift < 64 bits towards the most significant bit
    WIDE_WORD_INLINE void shift_left_bits(uint32_t shift)
    {
        const __m128i count = _mm_cvtsi32_si128(shift);
        const __m128i rest = _mm_cvtsi32_si128(64 - shift);
        __m256i lower = _mm256_setzero_si256(); // the rotated carries of the next less significant part

        for (uint16_t i = 0; i < num_vecs; i++) {
            const __m256i carries = _mm256_permute4x64_epi64(_mm256_srl_epi64(vecs[i], rest), 0b10010011);
            vecs[i] = _mm256_or_si256(_mm256_sll_epi64(vecs[i], count), _mm256_blend_epi32(carries, lower, 0b00000011));
            lower = carries;
        }
    }

    // shifts by shift < 64 bits towards the least significant bit
    WIDE_WORD_INLINE void shift_right_bits(uint32_t shift)
    {
        const __m128i count = _mm_cvtsi32_si128(shift);
        const __m128i rest = _mm_cvtsi32_si128(64 - shift);
        __m256i higher = _mm256_setzero_si256(); // the rotated carries of the next more significant part

        for (uint16_t i = num_vecs; i-- > 0;) {
            const __m256i carries = _mm256_permute4x64_epi64(_mm256_sll_epi64(vecs[i], rest), 0b00111001);
            vecs[i] = _mm256_or_si256(_mm256_srl_epi64(vecs[i], count), _mm256_blend_epi32(carries, higher, 0b11000000));
            higher = carries;
        }
    }

public:
    wide_word() = default;

    WIDE_WORD_INLINE wide_word(uint64_t value)
    {
        vecs[0] = _mm256_set_epi64x(0, 0, 0, static_cast<int64_t>(value));
        for (uint16_t i = 1; i < num_vecs; i++) vecs[i] = _mm256_setzero_si256();
    }

    // to bool: whether a bit is set; to an integer: the lowest bits
    template <std::integral int_t>
    WIDE_WORD_INLINE explicit operator int_t() const
    {
        if constexpr (std::is_same_v<int_t, bool>) {
            __m256i any = vecs[0];
            for (uint16_t i = 1; i < num_vecs; i++) any = _mm256_or_si256(any, vecs[i]);
            return !_mm256_testz_si256(any, any);
        } else {
            return static_cast<int_t>(_mm_cvtsi128_si64(_mm256_castsi256_si128(vecs[0])));
        }
    }

    WIDE_WORD_INLINE wide_word operator~() const
    {
        wide_word result;
        const __m256i ones = _mm256_set1_epi64x(-1);
        for (uint16_t i = 0; i < num_vecs; i++) result.vecs[i] = _mm256_xor_si256(vecs[i], ones);
        return result;
    }

    WIDE_WORD_INLINE wide_word& operator&=(const wide_word& other)
    {
        for (uint16_t i = 0; i < num_vecs; i++) vecs[i] = _mm256_and_si256(vecs[i], other.vecs[i]);
        return *this;
    }

    WIDE_WORD_INLINE wide_word& operator|=(const wide_word& other)
    {
        for (uint16_t i = 0; i < num_vecs; i++) vecs[i] = _mm256_or_si256(vecs[i], other.vecs[i]);
        return *this;
    }

    WIDE_WORD_INLINE wide_word& operator^=(const wide_word& other)
    {
        for (uint16_t i = 0; i < num_vecs; i++) vecs[i] = _mm256_xor_si256(vecs[i], other.vecs[i]);
        return *this;
    }

    // adds the 64-bit parts, then the carries: a part receives a carry if the part below overflows (generates), or is all
    // ones (propagates) and receives a carry itself; with one bit per part, these are the bits of ((generate << 1) + propagate)
    // that differ from propagate
    WIDE_WORD_INLINE wide_word& operator+=(const wide_word& other)
    {
        const __m256i sign = _mm256_set1_epi64x(INT64_MIN);
        const __m256i ones = _mm256_set1_epi64x(-1);
        uint32_t generate = 0;
        uint32_t propagate = 0;

        for (uint16_t i = 0; i < num_vecs; i++) {
            const __m256i sum = _mm256_add_epi64(vecs[i], other.vecs[i]);
            const __m256i overflow = _mm256_cmpgt_epi64(_mm256_xor_si256(vecs[i], sign), _mm256_xor_si256(sum, sign));
            generate |= uint32_t(_mm256_movemask_pd(_mm256_castsi256_pd(overflow))) << (4 * i);
            propagate |= uint32_t(_mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpeq_epi64(sum, ones)))) << (4 * i);
            vecs[i] = sum;
        }

        const uint32_t carries = ((generate << 1) + propagate) ^ propagate;
        const __m256i lanes = _mm256_set_epi64x(3, 2, 1, 0);
        const __m256i one = _mm256_set1_epi64x(1);

        for (uint16_t i = 0; i < num_vecs; i++) {
            const __m256i carry = _mm256_and_si256(_mm256_srlv_epi64(_mm256_set1_epi64x((carries >> (4 * i)) & 0b1111), lanes), one);
            vecs[i] = _mm256_add_epi64(vecs[i], carry);
        }

        return *this;
    }

    WIDE_WORD_INLINE wide_word& operator-=(const wide_word& other)
    {
        *this += ~other;
        return *this += wide_word(1);
    }

    WIDE_WORD_INLINE wide_word& operator<<=(uint32_t shift)
    {
        if (shift >= 64) {
            alignas(32) uint64_t limbs[num_limbs];
            store(limbs);
            const uint32_t limb_shift = shift / 64;
            for (int32_t i = num_limbs - 1; i >= 0; i--) limbs[i] = uint32_t(i) >= limb_shift ? limbs[i - limb_shift] : 0;
            load(limbs);
            shift %= 64;
        }

        shift_left_bits(shift);
        return *this;
    }

    WIDE_WORD_INLINE wide_word& operator>>=(uint32_t shift)
    {
        if (shift >= 64) {
            alignas(32) uint64_t limbs[num_limbs];
            store(limbs);
            const uint32_t limb_shift = shift / 64;
            for (uint32_t i = 0; i < num_limbs; i++) limbs[i] = i + limb_shift < num_limbs ? limbs[i + limb_shift] : 0;
            load(limbs);
            shift %= 64;
        }

        shift_right_bits(shift);
        return *this;
    }

    friend WIDE_WORD_INLINE wide_word operator&(wide_word a, const wide_word& b) { return a &= b; }
    friend WIDE_WORD_INLINE wide_word operator|(wide_word a, const wide_word& b) { return a |= b; }
    friend WIDE_WORD_INLINE wide_word operator^(wide_word a, const wide_word& b) { return a ^= b; }
    friend WIDE_WORD_INLINE wide_word operator+(wide_word a, const wide_word& b) { return a += b; }
    friend WIDE_WORD_INLINE wide_word operator-(wide_word a, const wide_word& b) { return a -= b; }
    friend WIDE_WORD_INLINE wide_word operator<<(wide_word a, uint32_t shift) { return a <<= shift; }
    friend WIDE_WORD_INLINE wide_word operator>>(wide_word a, uint32_t shift) { return a >>= shift; }

    friend WIDE_WORD_INLINE bool operator==(const wide_word& a, const wide_word& b)
    {
        __m256i diff = _mm256_xor_si256(a.vecs[0], b.vecs[0]);
        for (uint16_t i = 1; i < num_vecs; i++) diff = _mm256_or_si256(diff, _mm256_xor_si256(a.vecs[i], b.vecs[i]));
        return _mm256_testz_si256(diff, diff);
    }

    // returns whether bit b is set
    WIDE_WORD_INLINE bool bit(uint32_t b) const
    {
        alignas(32) uint64_t limbs[num_limbs];
        store(limbs);
        return (limbs[b / 64] >> (b % 64)) & 1;
    }

    // returns the number of set bits
    WIDE_WORD_INLINE int32_t count_ones() const
    {
        alignas(32) uint64_t limbs[num_limbs];
        store(limbs);
        int32_t count = 0;
        for (uint16_t i = 0; i < num_limbs; i++) count += std::popcount(limbs[i]);
        return count;
    }
};

#else

/**
 * @brief unsigned integer of num_bits bits with the operators of a built-in unsigned integer
 * @tparam num_bits number of bits, a multiple of 64
 */
template <uint16_t num_bits>
class wide_word {
    static_assert(num_bits % 64 == 0);

    static constexpr uint16_t num_limbs = num_bits / 64;

    uint64_t limbs[num_limbs]; // the 64-bit parts, the least significant first

public:
    wide_word() = default;

    WIDE_WORD_INLINE constexpr wide_word(uint64_t value) : limbs{}
    {
        limbs[0] = value;
    }

    // to bool: whether a bit is set; to an integer: the lowest bits
    template <std::integral int_t>
    WIDE_WORD_INLINE constexpr explicit operator int_t() const
    {
        if constexpr (std::is_same_v<int_t, bool>) {
            uint64_t any = 0;
            for (uint16_t i = 0; i < num_limbs; i++) any |= limbs[i];
            return any != 0;
        } else {
            return static_cast<int_t>(limbs[0]);
        }
    }

    WIDE_WORD_INLINE constexpr wide_word operator~() const
    {
        wide_word result;
        for (uint16_t i = 0; i < num_limbs; i++) result.limbs[i] = ~limbs[i];
        return result;
    }

    WIDE_WORD_INLINE constexpr wide_word& operator&=(const wide_word& other)
    {
        for (uint16_t i = 0; i < num_limbs; i++) limbs[i] &= other.limbs[i];
        return *this;
    }

    WIDE_WORD_INLINE constexpr wide_word& operator|=(const wide_word& other)
    {
        for (uint16_t i = 0; i < num_limbs; i++) limbs[i] |= other.limbs[i];
        return *this;
    }

    WIDE_WORD_INLINE constexpr wide_word& operator^=(const wide_word& other)
    {
        for (uint16_t i = 0; i < num_limbs; i++) limbs[i] ^= other.limbs[i];
        return *this;
    }

    WIDE_WORD_INLINE constexpr wide_word& operator+=(const wide_word& other)
    {
        __uint128_t carry = 0;

        for (uint16_t i = 0; i < num_limbs; i++) {
            carry += limbs[i];
            carry += other.limbs[i];
            limbs[i] = static_cast<uint64_t>(carry);
            carry >>= 64;
        }

        return *this;
    }

    WIDE_WORD_INLINE constexpr wide_word& operator-=(const wide_word& other)
    {
        return *this += ~other + wide_word(1);
    }

    WIDE_WORD_INLINE constexpr wide_word& operator<<=(uint32_t shift)
    {
        const uint32_t limb_shift = shift / 64;
        const uint32_t bit_shift = shift % 64;

        for (int32_t i = num_limbs - 1; i >= 0; i--) {
            uint64_t value = 0;

            if (uint32_t(i) >= limb_shift) {
                value = limbs[i - limb_shift] << bit_shift;
                if (bit_shift != 0 && uint32_t(i) > limb_shift) value |= limbs[i - limb_shift - 1] >> (64 - bit_shift);
            }

            limbs[i] = value;
        }

        return *this;
    }

    WIDE_WORD_INLINE constexpr wide_word& operator>>=(uint32_t shift)
    {
        const uint32_t limb_shift = shift / 64;
        const uint32_t bit_shift = shift % 64;

        for (uint32_t i = 0; i < num_limbs; i++) {
            uint64_t value = 0;

            if (i + limb_shift < num_limbs) {
                value = limbs[i + limb_shift] >> bit_shift;
                if (bit_shift != 0 && i + limb_shift + 1 < num_limbs) value |= limbs[i + limb_shift + 1] << (64 - bit_shift);
            }

            limbs[i] = value;
        }

        return *this;
    }

    friend WIDE_WORD_INLINE constexpr wide_word operator&(wide_word a, const wide_word& b) { return a &= b; }
    friend WIDE_WORD_INLINE constexpr wide_word operator|(wide_word a, const wide_word& b) { return a |= b; }
    friend WIDE_WORD_INLINE constexpr wide_word operator^(wide_word a, const wide_word& b) { return a ^= b; }
    friend WIDE_WORD_INLINE constexpr wide_word operator+(wide_word a, const wide_word& b) { return a += b; }
    friend WIDE_WORD_INLINE constexpr wide_word operator-(wide_word a, const wide_word& b) { return a -= b; }
    friend WIDE_WORD_INLINE constexpr wide_word operator<<(wide_word a, uint32_t shift) { return a <<= shift; }
    friend WIDE_WORD_INLINE constexpr wide_word operator>>(wide_word a, uint32_t shift) { return a >>= shift; }

    friend WIDE_WORD_INLINE constexpr bool operator==(const wide_word& a, const wide_word& b)
    {
        uint64_t diff = 0;
        for (uint16_t i = 0; i < num_limbs; i++) diff |= a.limbs[i] ^ b.limbs[i];
        return diff == 0;
    }

    // returns whether bit b is set
    WIDE_WORD_INLINE constexpr bool bit(uint32_t b) const
    {
        return (limbs[b / 64] >> (b % 64)) & 1;
    }

    // returns the number of set bits
    WIDE_WORD_INLINE constexpr int32_t count_ones() const
    {
        int32_t count = 0;
        for (uint16_t i = 0; i < num_limbs; i++) count += std::popcount(limbs[i]);
        return count;
    }
};

#endif

/**
 * @brief counts the number of set bits in x
 * @tparam num_bits number of bits of x
 * @param x a word
 * @return the number of set bits in x
 */
template <uint16_t num_bits>
WIDE_WORD_INLINE int32_t popcount(const wide_word<num_bits>& x)
{
    return x.count_ones();
}

/**
 * @brief returns whether bit b of x is set
 * @tparam word_t type of x, a built-in unsigned integer
 * @param x a word
 * @param b a bit position, less than the number of bits of word_t
 */
template <typename word_t>
WIDE_WORD_INLINE bool test_bit(const word_t& x, uint32_t b)
{
    return (x >> b) & 1;
}

template <uint16_t num_bits>
WIDE_WORD_INLINE bool test_bit(const wide_word<num_bits>& x, uint32_t b)
{
    return x.bit(b);
}

#endif
