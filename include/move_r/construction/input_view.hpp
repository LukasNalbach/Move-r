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

#ifndef MOVE_R_MOVE_R_CONSTRUCTION_INPUT_VIEW_HPP
#define MOVE_R_MOVE_R_CONSTRUCTION_INPUT_VIEW_HPP

#include <misc/files.hpp>
#include <misc/utils.hpp>
#include <misc/hash.hpp>
#include <misc/sort.hpp>
#include <omp.h>
#include <tsl/sparse_map.h>

#include <algorithm>
#include <bit>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <memory>
#include <string>
#include <stdexcept>
#include <type_traits>
#include <vector>

/**
 * @brief the unsigned symbol type with which an input with symbols of type sym_t is read (the unsigned
 *        integer type of the same width, or the packed type itself)
 * @tparam sym_t symbol type of the input: an integer type of at most 8 bytes or a packed unsigned integer type
 *         (uint24_t, uint40_t, uint48_t or uint56_t)
 */
template <typename sym_t>
struct move_r_symbol {
    static_assert(std::is_integral_v<sym_t> && sizeof(sym_t) <= sizeof(uint64_t),
        "move_r input symbols have to be integers of at most 8 bytes or packed unsigned integers");
    using type = std::make_unsigned_t<sym_t>;
};

template <uint8_t bytes> struct move_r_symbol<packed_unsigned<bytes>> { using type = packed_unsigned<bytes>; };

template <typename sym_t> using move_r_symbol_t = typename move_r_symbol<sym_t>::type;

// calls X(symbol type) for the unsigned symbol type of each width of an input symbol, 1 to 8 bytes
#define MOVE_R_FOR_EACH_INPUT_WIDTH(X) \
    X(uint8_t) \
    X(uint16_t) \
    X(uint24_t) \
    X(uint32_t) \
    X(uint40_t) \
    X(uint48_t) \
    X(uint56_t) \
    X(uint64_t)

/**
 * @brief calls fnc(std::type_identity<sym_t>()) with the symbol type sym_t of the given width
 * @param width width of an input symbol in bytes, in [1, 8]
 * @param fnc the function
 */
template <typename fnc_t>
inline decltype(auto) move_r_with_symbol_type(uint8_t width, fnc_t&& fnc)
{
    switch (width) {
        #define MOVE_R_WIDTH_CASE(sym_t) case sizeof(sym_t): return fnc(std::type_identity<sym_t>());
        MOVE_R_FOR_EACH_INPUT_WIDTH(MOVE_R_WIDTH_CASE)
        #undef MOVE_R_WIDTH_CASE
    }

    throw std::invalid_argument("move_r: the width of an input symbol must be 1 to " + std::to_string(sizeof(uint64_t))
        + " bytes");
}

/**
 * @brief type-erased view of a contiguous container storing the input; forwards the container operations
 *        that the construction needs and stores the width of a symbol
 */
class move_r_input_view {
private:
    void* object = nullptr;
    char* (*data_function)(void*) = nullptr;
    uint64_t (*size_function)(const void*) = nullptr;
    void (*resize_function)(void*, uint64_t) = nullptr;
    void (*reserve_function)(void*, uint64_t) = nullptr;
    void (*shrink_function)(void*) = nullptr;
    uint8_t sym_width = 1; // the width of one symbol in bytes

public:
    /**
     * @brief constructs a view of a container
     * @tparam container_t type of the container; its value type must be an integer type of at most 8 bytes or a
     *         packed unsigned integer type (see move_r_symbol)
     * @param input the container
     */
    template <typename container_t>
    explicit move_r_input_view(container_t& input)
        : object(&input)
        , data_function([](void* ptr) {
            auto& container = *static_cast<container_t*>(ptr);
            return reinterpret_cast<char*>(container.data());
        })
        , size_function([](const void* ptr) {
            return uint64_t(static_cast<const container_t*>(ptr)->size());
        })
        , resize_function([](void* ptr, uint64_t size) {
            // new symbols are not initialized, if the container supports it (see resize_for_overwrite)
            resize_for_overwrite(*static_cast<container_t*>(ptr), size);
        })
        , reserve_function([](void* ptr, uint64_t size) {
            if constexpr (requires(container_t& c) { c.reserve(uint64_t { 0 }); })
                static_cast<container_t*>(ptr)->reserve(size);
        })
        , shrink_function([](void* ptr) {
            if constexpr (requires(container_t& container) { container.shrink_to_fit(); })
                static_cast<container_t*>(ptr)->shrink_to_fit();
        })
    {
        using input_symbol_t = typename std::remove_cvref_t<container_t>::value_type;
        // checks the symbol type (see move_r_symbol)
        using checked_t [[maybe_unused]] = move_r_symbol_t<input_symbol_t>;
        sym_width = uint8_t(sizeof(input_symbol_t));
    }

    // returns the width of one symbol in bytes
    uint8_t symbol_width() const { return sym_width; }
    // returns whether the input uses a byte alphabet
    bool byte_alphabet() const { return sym_width == 1; }
    // returns a pointer to the symbols in the container
    char* data() { return data_function(object); }
    // returns a pointer to the symbols in the container
    const char* data() const { return data_function(object); }
    // returns the number of symbols in the container
    uint64_t size() const { return size_function(object); }
    // resizes the container to size symbols (new symbols are not initialized, if the container supports it)
    void resize(uint64_t size) { resize_function(object, size); }
    // reserves space for size symbols (if the container supports it)
    void reserve(uint64_t size) { reserve_function(object, size); }
    // removes all symbols
    void clear() { resize(0); }
    // frees unused capacity (if the container supports it)
    void shrink_to_fit() { shrink_function(object); }

    /**
     * @brief returns a reference to the i-th symbol
     * @tparam type the symbol type of the input
     * @param i [0..size-1] index of a symbol
     */
    template <typename type>
    inline type& at(uint64_t i)
    {
        return *reinterpret_cast<type*>(data() + i * uint64_t(sym_width));
    }

    // returns the i-th symbol
    inline uint64_t get(uint64_t i) const
    {
        const char* src = data() + i * uint64_t(sym_width);
        uint64_t value = 0;
        std::memcpy(&value, src, sym_width);
        return value;
    }

    // sets the i-th symbol to value
    inline void set(uint64_t i, uint64_t value)
    {
        char* dst = data() + i * uint64_t(sym_width);
        std::memcpy(dst, &value, sym_width);
    }

    // appends a symbol
    void push_back(uint64_t value)
    {
        uint64_t old_size = size();
        // reserves space for exactly one more symbol
        reserve(old_size + 1);
        resize(old_size + 1);
        set(old_size, value);
    }
};

/**
 * @brief view of a pattern stored in a contiguous container (a pointer and a length)
 * @tparam sym_t unsigned symbol type with the width of the symbols of the pattern (see move_r_symbol)
 */
template <typename sym_t>
class move_r_pattern_view {
private:
    const sym_t* ptr = nullptr;
    uint64_t len = 0;

public:
    using value_type = sym_t;

    move_r_pattern_view() = default;

    // constructs a view of the pattern symbols[0..size-1]
    move_r_pattern_view(const sym_t* symbols, uint64_t size) : ptr(symbols), len(size) { }

    /**
     * @brief constructs a view of the pattern stored in a contiguous container
     * @tparam pat_t type of the container; its symbols must have the width of sym_t
     * @param pattern the container
     */
    template <typename pat_t>
        requires(!std::is_same_v<std::remove_cvref_t<pat_t>, move_r_pattern_view>)
    move_r_pattern_view(const pat_t& pattern)
        : ptr(reinterpret_cast<const sym_t*>(std::data(pattern)))
        , len(uint64_t(std::size(pattern)))
    {
        // symbols of a signed type are read as the unsigned type of the same width (see move_r_symbol)
        static_assert(sizeof(typename std::remove_cvref_t<pat_t>::value_type) == sizeof(sym_t),
            "the symbols of the pattern must have the width of sym_t");
    }

    uint64_t size() const { return len; }
    bool empty() const { return len == 0; }
    const sym_t* data() const { return ptr; }
    const sym_t& operator[](uint64_t i) const { return ptr[i]; }
    const sym_t* begin() const { return ptr; }
    const sym_t* end() const { return ptr + len; }
};

// type of the pattern view for a pattern stored in a container of type pat_t
template <typename pat_t>
using move_r_pattern_view_t = move_r_pattern_view<
    move_r_symbol_t<typename std::remove_cvref_t<pat_t>::value_type>>;

// true <=> pat_t is a pattern view
template <typename pat_t> inline constexpr bool move_r_is_pattern_view = false;
template <typename sym_t> inline constexpr bool move_r_is_pattern_view<move_r_pattern_view<sym_t>> = true;

/**
 * @brief type of the container into which an input file is read for an input of type input_t: a vector of the
 *        symbol type of input_t, or a byte vector, if input_t is move_r_input_view
 */
template <typename input_t>
struct move_r_input_buffer {
    using type = no_init_vector<move_r_symbol_t<
        typename std::remove_cvref_t<input_t>::value_type>>;
};

template <> struct move_r_input_buffer<move_r_input_view> { using type = no_init_vector<uint8_t>; };

template <typename input_t> using move_r_input_buffer_t = typename move_r_input_buffer<input_t>::type;

/**
 * @brief calls f with a view of the input
 * @param input a container storing the input, or a move_r_input_view
 * @param f function to call with the view
 */
template <typename input_t, typename f_t>
inline void move_r_with_view(input_t& input, f_t&& f)
{
    if constexpr (std::is_same_v<std::remove_cvref_t<input_t>, move_r_input_view>) {
        f(input);
    } else {
        move_r_input_view view(input);
        f(view);
    }
}

/**
 * @brief maps the symbols of an integer-alphabet input to their effective symbols, the internal symbols
 *        [1..sigma-1] of the index (interface; implemented by move_r_int_alphabet_of for each symbol width)
 */
class move_r_int_alphabet {
public:
    virtual ~move_r_int_alphabet() = default;

    /** @brief collects the distinct symbols of T[0..n-2] and counts their occurrences */
    virtual void collect(move_r_input_view& T, uint64_t n) = 0;
    /**
     * @brief collects the distinct symbols of a file storing T[0..n-2] and counts their occurrences
     * @param path the file
     * @param n length of T (including the terminator, which the file does not store)
     * @param p the number of threads to use
     */
    virtual void collect_file(const std::string& path, uint64_t n, uint64_t p) = 0;
    /** @brief returns the number of distinct symbols */
    virtual uint64_t num_symbols() const = 0;
    /**
     * @brief stores the distinct symbols in increasing order
     * @param symbols [0..sigma-2] stores the distinct symbols
     * @param p the number of threads to use
     */
    virtual void sorted_symbols(std::vector<uint64_t>& symbols, uint64_t p) const = 0;
    /**
     * @brief stores the numbers of occurrences of the distinct symbols
     * @param symbols [0..sigma-2] the distinct symbols in increasing order
     * @param freq [0..sigma-2] freq[i] stores the number of occurrences of symbols[i]
     * @param p the number of threads to use
     */
    virtual void frequencies(const std::vector<uint64_t>& symbols,
        std::vector<uint64_t>& freq, uint64_t p) const = 0;
    /**
     * @brief maps each distinct symbol to its effective symbol
     * @param symbols [0..sigma-2] the distinct symbols in increasing order
     * @param effective [0..sigma-2] effective[i] is the effective symbol of symbols[i]; if it is empty,
     *        symbols[i] is mapped to i + 1 (see move_r's _alphabet_order)
     * @param p the number of threads to use
     */
    virtual void assign_effective_symbols(const std::vector<uint64_t>& symbols,
        const std::vector<uint64_t>& effective, uint64_t p) = 0;
    /** @brief replaces each symbol of T[0..n-2] with its effective symbol */
    virtual void map_text(move_r_input_view& T, uint64_t n, uint64_t p) const = 0;
    /**
     * @brief stores the effective symbols of symbols (after assign_effective_symbols)
     * @param src [0..len-1] the symbols, stored with the width of an input symbol
     * @param len number of symbols
     * @param dst [0..len-1] stores the effective symbols
     */
    virtual void map_symbols(const char* src, uint64_t len, uint64_t* dst) const = 0;

    /** @brief returns a map for an input with symbols of the given width (in bytes) */
    static std::unique_ptr<move_r_int_alphabet> make(uint8_t width);
};

/**
 * @brief an unsigned integer stored in width_bits / 8 bytes, implicitly convertible to and from uint128_t
 * @tparam width_bits the width in bits; a multiple of 8, at most 128
 */
template <uint8_t width_bits>
struct packed_uint_t {
    static_assert(std::endian::native == std::endian::little,
        "packed_uint_t stores the low bytes of a value, which assumes a little-endian layout");
    static_assert(width_bits % 8 == 0, "the width of a packed_uint_t is a whole number of bytes");
    static_assert(8 <= width_bits && width_bits <= 128, "a packed_uint_t is at most 128 bits wide");

    static constexpr uint8_t width = width_bits / 8; // bytes per value
    // the largest storable value
    static constexpr uint128_t max = width_bits == 128
        ? ~uint128_t(0) : (uint128_t { 1 } << width_bits) - 1;

    uint8_t bytes[width];

    // narrowest built-in unsigned integer type with at least width_bits bits; the bytes are copied from/to it
    using word_t = std::conditional_t<width_bits <= 32, uint32_t,
        std::conditional_t<width_bits <= 64, uint64_t, uint128_t>>;

    packed_uint_t() = default;

    packed_uint_t(uint128_t value)
    {
        assert(value <= max);
        const word_t word = word_t(value);
        std::memcpy(bytes, &word, width);
    }

    // returns the stored value
    operator uint128_t() const
    {
        word_t word = 0;
        std::memcpy(&word, bytes, width);
        return uint128_t(word);
    }

    friend bool operator==(const packed_uint_t& left, const packed_uint_t& right)
    {
        return std::memcmp(left.bytes, right.bytes, width) == 0;
    }
};

static_assert(sizeof(packed_uint_t<64>) == 8);
static_assert(alignof(packed_uint_t<64>) == 1);
static_assert(sizeof(packed_uint_t<40>) == uint40_t::width);

/**
 * @brief implementation of move_r_int_alphabet for one symbol width; the keys of the hash map have the
 *        width of sym_t
 */
template <typename sym_t>
class move_r_int_alphabet_of : public move_r_int_alphabet {
private:
    using key_t = packed_uint_t<8 * sizeof(sym_t)>;

    // symbol -> its number of occurrences (until assign_effective_symbols), then symbol -> its effective symbol
    tsl::sparse_map<key_t, uint48_t, wide_hash> map;

    // returns the key of a symbol (its bit pattern, without sign extension)
    static inline key_t key_of(sym_t sym) { return key_t(symbol_bit_pattern(sym)); }

public:
    void collect(move_r_input_view& T, uint64_t n) override
    {
        const sym_t* text = reinterpret_cast<const sym_t*>(T.data());

        // until assign_effective_symbols is called, the map stores the number of occurrences of each symbol
        for (uint64_t i = 0; i < n - 1; i++) {
            auto [it, inserted] = map.try_emplace(key_of(text[i]), 0);
            it.value() = uint48_t(uint64_t(it->second) + 1);
        }
    }

    void collect_file(const std::string& path, uint64_t n, uint64_t p) override
    {
        // each thread counts the symbols of a range of the file; the counts are added up afterwards
        std::vector<tsl::sparse_map<key_t, uint64_t, wide_hash>> counts(p);
        const uint64_t symbols_per_read = direct_reader::block / sizeof(sym_t);

        #pragma omp parallel num_threads(p)
        {
            const uint64_t i_p = omp_get_thread_num();
            const uint64_t beg = i_p * (n - 1) / p;
            const uint64_t end = (i_p + 1) * (n - 1) / p;
            direct_reader file(path);

            for (uint64_t pos = beg; pos < end; pos += symbols_per_read) {
                const uint64_t len = std::min(symbols_per_read, end - pos);
                const uint8_t* data = file.read(pos * sizeof(sym_t), (pos + len) * sizeof(sym_t));

                for (uint64_t i = 0; i < len; i++) {
                    sym_t sym;
                    std::memcpy(&sym, data + i * sizeof(sym_t), sizeof(sym_t));
                    counts[i_p][key_of(sym)]++;
                }
            }
        }

        for (const tsl::sparse_map<key_t, uint64_t, wide_hash>& counts_thr : counts) {
            for (const auto& [key, count] : counts_thr) {
                auto [it, inserted] = map.try_emplace(key, 0);
                it.value() = uint48_t(uint64_t(it->second) + count);
            }
        }
    }

    uint64_t num_symbols() const override { return map.size(); }

    void sorted_symbols(std::vector<uint64_t>& symbols, uint64_t p) const override
    {
        resize_for_overwrite(symbols, map.size());
        uint64_t sym_cur = 0;

        for (const auto& kv : map) {
            symbols[sym_cur] = uint64_t(kv.first);
            sym_cur++;
        }

        radix_sort(symbols.begin(), symbols.end(), key_of_value {}, p);
    }

    void frequencies(const std::vector<uint64_t>& symbols,
        std::vector<uint64_t>& freq, uint64_t p) const override
    {
        resize_for_overwrite(freq, symbols.size());

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < symbols.size(); i++) {
            freq[i] = uint64_t((*map.find(key_t(symbols[i]))).second);
        }
    }

    void assign_effective_symbols(const std::vector<uint64_t>& symbols,
        const std::vector<uint64_t>& effective, uint64_t p) override
    {
        uint64_t sigma = symbols.size() + 1;
        const bool identity = effective.empty();

        if (sigma > p) {
            #pragma omp parallel num_threads(p)
            {
                uint64_t i_p = omp_get_thread_num();

                uint64_t b = 1 + i_p * ((sigma - 1) / p);
                uint64_t e = i_p == p - 1 ? sigma - 1 : ((i_p + 1) * ((sigma - 1) / p));

                for (uint64_t i = b; i <= e; i++) {
                    map[key_t(symbols[i - 1])] = uint48_t(identity ? i : effective[i - 1]);
                }
            }
        } else {
            for (uint64_t i = 1; i < sigma; i++) {
                map[key_t(symbols[i - 1])] = uint48_t(identity ? i : effective[i - 1]);
            }
        }
    }

    void map_text(move_r_input_view& T, uint64_t n, uint64_t p) const override
    {
        sym_t* text = reinterpret_cast<sym_t*>(T.data());

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < n - 1; i++) {
            text[i] = sym_t(uint64_t((*map.find(key_of(text[i]))).second));
        }
    }

    void map_symbols(const char* src, uint64_t len, uint64_t* dst) const override
    {
        for (uint64_t i = 0; i < len; i++) {
            sym_t sym;
            std::memcpy(&sym, src + i * sizeof(sym_t), sizeof(sym_t));
            dst[i] = uint64_t((*map.find(key_of(sym))).second);
        }
    }
};

inline std::unique_ptr<move_r_int_alphabet> move_r_int_alphabet::make(uint8_t width)
{
    return move_r_with_symbol_type(width, [](auto type) -> std::unique_ptr<move_r_int_alphabet> {
        return std::make_unique<move_r_int_alphabet_of<typename decltype(type)::type>>();
    });
}

// the editor parser of VS Code (__INTELLISENSE__) does not resolve the extern declarations of constrained member
// templates; it does not need them
#if !defined(STRING_UTILS_NO_EXTERN_TEMPLATES) && !defined(__INTELLISENSE__)
#define MOVE_R_DECLARE_INT_ALPHABET(sym_t) extern template class move_r_int_alphabet_of<sym_t>;
MOVE_R_FOR_EACH_INPUT_WIDTH(MOVE_R_DECLARE_INT_ALPHABET)
#undef MOVE_R_DECLARE_INT_ALPHABET
#endif

#endif // MOVE_R_MOVE_R_CONSTRUCTION_INPUT_VIEW_HPP
