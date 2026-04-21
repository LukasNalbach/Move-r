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

#pragma once

#include <algorithm>
#include <climits>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>
#include <bit>
#include <random>

#include <malloc_count.h>

__extension__ typedef unsigned __int128 uint128_t;

template <typename int_t, typename x_t, typename y_t>
inline static int_t abs_diff(x_t x, y_t y)
{
    return std::abs(int_t{x} - int_t{y});
}

inline static char uchar_to_char(uint8_t c)
{
    return *reinterpret_cast<char*>(&c);
}

inline static uint8_t char_to_uchar(char c)
{
    return *reinterpret_cast<uint8_t*>(&c);
}

template <typename T>
static bool contains(const std::vector<T>& vec, T val)
{
    return std::find(vec.begin(), vec.end(), val) != vec.end();
}

template <typename T>
static void remove(std::vector<T>& vec, T val)
{
    while (contains(vec, val)) {
        vec.erase(std::find(vec.begin(), vec.end(), val));
    }
}

template <typename T>
static bool is_subset_of(std::vector<T> first, std::vector<T> second)
{
    std::sort(first.begin(), first.end());
    std::sort(second.begin(), second.end());
    return std::includes(second.begin(), second.end(), first.begin(), first.end());
}

template <typename T>
class no_init {
    static_assert(std::is_fundamental<T>::value);

private:
    T v_;

public:
    no_init() noexcept { }

    constexpr no_init(T value) noexcept
        : v_ { value } { }

    constexpr operator T() const noexcept { return v_; }
};

template <typename T, typename Alloc = std::allocator<T>>
class default_init_allocator : public Alloc {
    using a_t = std::allocator_traits<Alloc>;
    using Alloc::Alloc;

public:
    template <typename U>
    struct rebind { };
    template <typename U>
    void construct(U* ptr) noexcept(std::is_nothrow_default_constructible<U>::value) { ::new (static_cast<void*>(ptr)) U; }
    template <typename U, typename... Args>
    void construct(U* ptr, Args&&... args) { }
};

static void no_init_resize(std::string& str, size_t size)
{
    (*reinterpret_cast<std::basic_string<char, std::char_traits<char>, default_init_allocator<char>>*>(&str)).resize(size);
}

template <typename T>
static void no_init_resize(std::vector<T>& vec, size_t size)
{
    (*reinterpret_cast<std::vector<no_init<T>>*>(&vec)).resize(size);
}

template <typename T1, typename T2>
static void no_init_resize(std::vector<std::pair<T1, T2>>& vec, size_t size)
{
    (*reinterpret_cast<std::vector<std::pair<no_init<T1>, no_init<T2>>>*>(&vec)).resize(size);
}

template <typename T>
static void no_init_resize(std::vector<std::tuple<T, T, T>>& vec, size_t size)
{
    (*reinterpret_cast<std::vector<std::tuple<no_init<T>, no_init<T>, no_init<T>>>*>(&vec)).resize(size);
}

enum direction_t : uint8_t {
    NO_DIR = 0,
    LEFT = 1,
    RIGHT = 2
};

template <direction_t dir>
inline static constexpr direction_t flip()
{
    if constexpr (dir == LEFT) return RIGHT;
    if constexpr (dir == RIGHT) return LEFT;
    return NO_DIR;
}

inline static direction_t flip(direction_t dir)
{
    if (dir == LEFT) return RIGHT;
    if (dir == RIGHT) return LEFT;
    return NO_DIR;
}

template <auto start, auto end, auto inc, class fnc_t>
inline static constexpr void for_constexpr(fnc_t fnc)
{
    if constexpr (start < end) {
        fnc(std::integral_constant<decltype(start), start>());
        for_constexpr<start + inc, end, inc>(fnc);
    }
}

template <auto... args, typename fnc_t>
constexpr void for_each_constexpr(fnc_t fnc)
{
    (fnc(std::integral_constant<decltype(args), args>{}), ...);
}

template <bool B, typename T>
struct constexpr_case {
    static constexpr bool value = B;
    using type = T;
};

template <bool B, typename TrueF, typename FalseF>
struct eval_if {
    using type = typename TrueF::type;
};

template <typename TrueF, typename FalseF>
struct eval_if<false, TrueF, FalseF> {
    using type = typename FalseF::type;
};

template <bool B, typename T, typename F>
using eval_if_t = typename eval_if<B, T, F>::type;

template <typename Head, typename... Tail>
struct constexpr_switch {
    using type = eval_if_t<Head::value, Head, constexpr_switch<Tail...>>;
};

template <typename T>
struct constexpr_switch<T> {
    using type = T;
};

template <bool B, typename T>
struct constexpr_switch<constexpr_case<B, T>> {
    static_assert(B, "!");
    using type = T;
};

template <typename Head, typename... Tail>
using constexpr_switch_t = typename constexpr_switch<Head, Tail...>::type;

template <typename uint_t>
inline static uint_t div_ceil(uint_t x, uint_t y)
{
    return x == 0 ? 0 : (1 + (x - 1) / y);
}

template <typename T>
inline static uint8_t byte_width(T val)
{
    return std::max<uint8_t>(1, div_ceil<uint8_t>(std::bit_width(uint64_t{val}), 8));
}

template <typename container_t>
static void log_contents(container_t container)
{
    if (container.empty()) return;

    for (uint64_t i = 0; i < container.size() - 1; i++) {
        std::cout << container[i] << ", ";
    }

    std::cout << container[container.size() - 1] << std::endl;
}

struct empty_t {};

template <typename T>
struct function_traits : public function_traits<decltype(&T::operator())> {};

template <typename class_t, typename return_t, typename... args>
struct function_traits<return_t(class_t::*)(args...) const> {
    static constexpr size_t arity = sizeof...(args);
    
    template <size_t N>
    using argument_type = std::tuple_element_t<N, std::tuple<args...>>;
};

template <typename word_t>
static int popcount(word_t x)
{
    static_assert(std::is_same_v<word_t, uint64_t> || std::is_same_v<word_t, __uint128_t>);
    
    if constexpr (std::is_same_v<word_t, __uint128_t>) {
        uint64_t high = x >> 64;
        uint64_t low = x;
        return __builtin_popcountll(high) + __builtin_popcountll(low);
    } else {
        return __builtin_popcountll(x);
    }
}