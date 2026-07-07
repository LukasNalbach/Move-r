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
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <bit>

#include <malloc_count.h>

#ifndef MOVE_R_USE_MALLOC_COUNT
#define MOVE_R_USE_MALLOC_COUNT 1
#endif

using uint128_t = __uint128_t;

/**
 * @brief returns the absolute difference of x and y
 * @tparam int_t signed integer type the operands are cast to before subtracting
 * @tparam x_t type of x
 * @tparam y_t type of y
 * @param x first operand
 * @param y second operand
 * @return the absolute difference |x - y|
 */
template <typename int_t, typename x_t, typename y_t>
inline static int_t abs_diff(x_t x, y_t y)
{
    return std::abs(int_t{x} - int_t{y});
}

/**
 * @brief reinterprets an unsigned byte as a (signed) char
 * @param c an unsigned byte
 * @return c reinterpreted as a char
 */
inline static char uchar_to_char(uint8_t c) { return *reinterpret_cast<char*>(&c); }

/**
 * @brief reinterprets a char as an unsigned byte
 * @param c a char
 * @return c reinterpreted as an unsigned byte
 */
inline static uint8_t char_to_uchar(char c) { return *reinterpret_cast<uint8_t*>(&c); }

/**
 * @brief reinterprets a one-byte symbol (char, uint8_t or int8_t) as an unsigned byte
 * @tparam sym_t a one-byte symbol type
 * @param c a symbol
 * @return c reinterpreted as an unsigned byte
 */
template <typename sym_t>
inline static uint8_t sym_to_uchar(sym_t c)
{
    static_assert(sizeof(sym_t) == 1);
    return *reinterpret_cast<uint8_t*>(&c);
}

/**
 * @brief returns whether vec contains val
 * @tparam T element type
 * @param vec a vector
 * @param val a value
 * @return whether vec contains val
 */
template <typename T>
static bool contains(const std::vector<T>& vec, T val) { return std::find(vec.begin(), vec.end(), val) != vec.end(); }

/**
 * @brief removes all occurrences of val from vec
 * @tparam T element type
 * @param vec a vector
 * @param val the value to remove
 */
template <typename T>
static void remove(std::vector<T>& vec, T val)
{
    while (contains(vec, val)) {
        vec.erase(std::find(vec.begin(), vec.end(), val));
    }
}

/**
 * @brief returns whether the values in first are a subset of the values in second
 * @tparam T element type
 * @param first a vector (passed by value, it is sorted internally)
 * @param second a vector (passed by value, it is sorted internally)
 * @return whether first is a subset of second
 */
template <typename T>
static bool is_subset_of(std::vector<T> first, std::vector<T> second)
{
    std::sort(first.begin(), first.end());
    std::sort(second.begin(), second.end());
    return std::includes(second.begin(), second.end(), first.begin(), first.end());
}

/**
 * @brief wrapper around a fundamental type T whose default constructor leaves the value uninitialized
 * @tparam T a fundamental type
 */
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

/**
 * @brief allocator that default-initializes (leaves uninitialized) the elements it constructs
 * @tparam T element type
 * @tparam Alloc underlying allocator
 */
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
    void construct(U*, Args&&...) { }
};

/**
 * @brief resizes str to size without initializing newly added characters
 * @param str the string to resize
 * @param size the new size
 */
inline void no_init_resize(std::string& str, size_t size)
{
    (*reinterpret_cast<std::basic_string<char, std::char_traits<char>, default_init_allocator<char>>*>(&str)).resize(size);
}

/**
 * @brief resizes vec to size without initializing newly added elements
 * @tparam T element type
 * @param vec the vector to resize
 * @param size the new size
 */
template <typename T>
static void no_init_resize(std::vector<T>& vec, size_t size) { (*reinterpret_cast<std::vector<no_init<T>>*>(&vec)).resize(size); }

/**
 * @brief resizes a vector of pairs to size without initializing newly added elements
 * @tparam T1 type of the first pair element
 * @tparam T2 type of the second pair element
 * @param vec the vector to resize
 * @param size the new size
 */
template <typename T1, typename T2>
static void no_init_resize(std::vector<std::pair<T1, T2>>& vec, size_t size)
{
    (*reinterpret_cast<std::vector<std::pair<no_init<T1>, no_init<T2>>>*>(&vec)).resize(size);
}

/**
 * @brief resizes a vector of 3-tuples to size without initializing newly added elements
 * @tparam T element type of the tuple
 * @param vec the vector to resize
 * @param size the new size
 */
template <typename T>
static void no_init_resize(std::vector<std::tuple<T, T, T>>& vec, size_t size)
{
    (*reinterpret_cast<std::vector<std::tuple<no_init<T>, no_init<T>, no_init<T>>>*>(&vec)).resize(size);
}

/**
 * @brief returns the file name component of a path (the part after the last '/' or '\\')
 * @param path the path
 * @return the file name (path itself if it contains no separator)
 */
inline std::string basename(const std::string& path) { return path.substr(path.find_last_of("/\\") + 1); }

/** @brief a search/reading direction */
enum direction_t : uint8_t {
    NO_DIR = 0,
    LEFT = 1,
    RIGHT = 2
};

/**
 * @brief flips a direction at compile time (LEFT <-> RIGHT, NO_DIR stays NO_DIR)
 * @tparam dir the direction to flip
 * @return the flipped direction
 */
template <direction_t dir>
inline static constexpr direction_t flip()
{
    if constexpr (dir == LEFT) return RIGHT;
    if constexpr (dir == RIGHT) return LEFT;
    return NO_DIR;
}

/**
 * @brief flips a direction (LEFT <-> RIGHT, NO_DIR stays NO_DIR)
 * @param dir the direction to flip
 * @return the flipped direction
 */
inline static direction_t flip(direction_t dir)
{
    if (dir == LEFT) return RIGHT;
    if (dir == RIGHT) return LEFT;
    return NO_DIR;
}

/**
 * @brief invokes fnc with the compile-time constant std::integral_constant for each value in [start, end) stepping by inc
 * @tparam start the first value
 * @tparam end the (exclusive) end value
 * @tparam inc the increment
 * @tparam fnc_t type of the function to invoke
 * @param fnc the function to invoke for each value
 */
template <auto start, auto end, auto inc, class fnc_t>
inline static constexpr void for_constexpr(fnc_t fnc)
{
    if constexpr (start < end) {
        fnc(std::integral_constant<decltype(start), start>());
        for_constexpr<start + inc, end, inc>(fnc);
    }
}

/**
 * @brief invokes fnc with the compile-time constant std::integral_constant for each value in the template argument pack
 * @tparam args the compile-time values to iterate over
 * @tparam fnc_t type of the function to invoke
 * @param fnc the function to invoke for each value
 */
template <auto... args, typename fnc_t>
constexpr void for_each_constexpr(fnc_t fnc)
{
    (fnc(std::integral_constant<decltype(args), args>{}), ...);
}

/**
 * @brief a single case of constexpr_switch: the type T is selected if the condition B holds
 * @tparam B the condition of the case
 * @tparam T the type selected if B holds
 */
template <bool B, typename T>
struct constexpr_case {
    static constexpr bool value = B;
    using type = T;
};

/**
 * @brief lazily evaluates ::type to TrueF::type if B holds, else to FalseF::type
 * @tparam B the condition
 * @tparam TrueF the case whose ::type is used if B holds
 * @tparam FalseF the case whose ::type is used otherwise
 */
template <bool B, typename TrueF, typename FalseF>
struct eval_if {
    using type = typename TrueF::type;
};

template <typename TrueF, typename FalseF>
struct eval_if<false, TrueF, FalseF> {
    using type = typename FalseF::type;
};

/** @brief alias for eval_if<B, T, F>::type */
template <bool B, typename T, typename F>
using eval_if_t = typename eval_if<B, T, F>::type;

/**
 * @brief selects, at compile time, the type of the first constexpr_case in the list whose condition holds
 * @tparam Head the first case
 * @tparam Tail the remaining cases
 */
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

/** @brief alias for constexpr_switch<Head, Tail...>::type */
template <typename Head, typename... Tail>
using constexpr_switch_t = typename constexpr_switch<Head, Tail...>::type;

/**
 * @brief computes the ceiling of x / y
 * @tparam uint_t unsigned integer type
 * @param x the dividend
 * @param y the divisor
 * @return ceil(x / y)
 */
template <typename uint_t>
inline static uint_t div_ceil(uint_t x, uint_t y) { return x == 0 ? 0 : (1 + (x - 1) / y); }

/**
 * @brief returns the number of bits needed to store the value val
 * @tparam T integer type of val
 * @param val a value
 * @return the number of bits needed to store val, as a uint8_t (suitable for interleaved-vector width lists)
 */
template <typename T>
inline static uint8_t bit_width(T val) { return std::bit_width(uint64_t(val)); }

/**
 * @brief returns the number of bytes needed to store the value val
 * @tparam T integer type of val
 * @param val a value
 * @return the number of bytes needed to store val (at least 1)
 */
template <typename T>
inline static uint8_t byte_width(T val)
{
    return std::max<uint8_t>(1, div_ceil<uint8_t>(std::bit_width(uint64_t{val}), 8));
}

/**
 * @brief prints the contents of a container to std::cout as a space-separated list
 * @tparam container_t container type (must support empty(), size() and operator[])
 * @param container the container to print
 */
template <typename container_t>
static void log_contents(container_t container)
{
    if (container.empty()) return;

    for (uint64_t i = 0; i < container.size() - 1; i++) {
        std::cout << container[i] << " ";
    }

    std::cout << container[container.size() - 1] << std::endl;
}

/**
 * @brief collects labeled rows of corresponding elements and prints them with the columns vertically
 *        aligned: each label is left-justified to the widest label and each column is right-justified
 *        to its widest element, with elements separated by single spaces (and no commas). Intended for
 *        logging corresponding arrays (e.g. the p, q and idx arrays of a move data structure) so that
 *        column i lines up across all rows.
 */
class aligned_log {
    std::vector<std::string> labels_;
    std::vector<std::vector<std::string>> rows_;

  public:
    /**
     * @brief adds a row labeled label whose i-th cell (for i in [0, count)) is the value cell(i)
     * @tparam cell_fnc_t type of the cell-value function
     * @param label the row label
     * @param count the number of cells in the row
     * @param cell function mapping a cell index to its (streamable) value
     */
    template <typename cell_fnc_t>
    void add_row(const std::string& label, uint64_t count, cell_fnc_t cell)
    {
        labels_.emplace_back(label);
        std::vector<std::string>& row = rows_.emplace_back();
        row.reserve(count);

        for (uint64_t i = 0; i < count; i++) {
            std::ostringstream s;
            s << cell(i);
            row.emplace_back(s.str());
        }
    }

    /**
     * @brief prints the collected rows with their columns vertically aligned to out, preceded by a
     *        header row showing the 0-based column indexes
     * @param out the output stream to print to
     */
    void print(std::ostream& out = std::cout) const
    {
        uint64_t label_width = 0;
        uint64_t num_columns = 0;

        for (const std::string& label : labels_)
            label_width = std::max<uint64_t>(label_width, label.size());
        for (const std::vector<std::string>& row : rows_)
            num_columns = std::max<uint64_t>(num_columns, row.size());

        // the 0-based column indexes participate in the column widths so that wide indexes still fit
        std::vector<uint64_t> column_width(num_columns, 0);
        for (uint64_t i = 0; i < num_columns; i++)
            column_width[i] = std::to_string(i).size();

        for (const std::vector<std::string>& row : rows_)
            for (uint64_t i = 0; i < row.size(); i++)
                column_width[i] = std::max<uint64_t>(column_width[i], row[i].size());

        // index header row (blank label area, then the 0-based column indexes)
        out << std::string(label_width, ' ');
        for (uint64_t i = 0; i < num_columns; i++) {
            std::string index = std::to_string(i);
            out << ' ' << std::string(column_width[i] - index.size(), ' ') << index;
        }
        out << std::endl;

        for (uint64_t r = 0; r < rows_.size(); r++) {
            out << labels_[r] << std::string(label_width - labels_[r].size(), ' ');

            for (uint64_t i = 0; i < rows_[r].size(); i++)
                out << ' ' << std::string(column_width[i] - rows_[r][i].size(), ' ') << rows_[r][i];

            out << std::endl;
        }
    }
};

/**
 * @brief logs a single array labeled label, vertically aligned and preceded by a 0-based index header
 * @tparam container_t container type (must support empty(), size() and operator[])
 * @param label the array's label
 * @param container the array to log
 */
template <typename container_t>
static void log_indexed(const std::string& label, const container_t& container)
{
    if (container.empty()) return;
    aligned_log log;
    log.add_row(label, container.size(), [&](uint64_t i) { return container[i]; });
    log.print();
}

/** @brief an empty placeholder type */
struct empty_t {};

/**
 * @brief extracts the arity and argument types of a callable (e.g. a lambda)
 * @tparam T the callable type
 */
template <typename T>
struct function_traits : public function_traits<decltype(&T::operator())> {};

template <typename class_t, typename return_t, typename... args>
struct function_traits<return_t(class_t::*)(args...) const> {
    static constexpr size_t arity = sizeof...(args);

    template <size_t N>
    using argument_type = std::tuple_element_t<N, std::tuple<args...>>;
};

/**
 * @brief counts the number of set bits in x
 * @tparam word_t word type (uint32_t, uint64_t or uint128_t)
 * @param x a word
 * @return the number of set bits in x
 */
template <typename word_t>
static int popcount(word_t x)
{
    static_assert(std::is_same_v<word_t, uint32_t> ||
                  std::is_same_v<word_t, uint64_t> ||
                  std::is_same_v<word_t, uint128_t>);
    
    if constexpr (std::is_same_v<word_t, uint128_t>) {
        uint64_t high = x >> 64;
        uint64_t low = x;
        return __builtin_popcountll(high) + __builtin_popcountll(low);
    } else if constexpr (std::is_same_v<word_t, uint64_t>) {
        return __builtin_popcountll(x);
    } else {
        return __builtin_popcount(x);
    }
}