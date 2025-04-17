/**
 * rmq/rmq_bender_farach_colton.hpp
 * part of pdinklag/lzend
 * 
 * MIT License
 * 
 * Copyright (c) Patrick Dinklage
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

#ifndef _RMQ_RMQ_BENDER_FARACH_COLTON_HPP
#define _RMQ_RMQ_BENDER_FARACH_COLTON_HPP

#include <bit>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <memory>

#include "compare.hpp"

namespace rmq {

template<std::totally_ordered Value, std::unsigned_integral Index = std::make_unsigned_t<Value>, bool minimum_ = true>
class RMQBenderFarachColton {
private:
    static constexpr bool leq(Value const& a, Value const& b) { return compare<minimum_>(a, b); }

    Value const* data_;

    using Level = std::unique_ptr<Index[]>;
    std::unique_ptr<Level[]> rmq_;

public:
    RMQBenderFarachColton() : data_(nullptr) {
    }

    RMQBenderFarachColton(Value const* data, size_t const n) : data_(data) {
        size_t const num_levels = std::bit_width(n) - 1;
        rmq_ = std::make_unique<Level[]>(num_levels);

        // build first level
        if (n > 1) {
            rmq_[0] = std::make_unique<Index[]>(n-1);
            for(Index i = 0; i < n-1; i++) {
                rmq_[0][i] = leq(data[i], data[i+1]) ? i : i+1;
            }
        }

        // build higher levels
        for(size_t level = 1; level < num_levels; level++) {
            size_t const level_size = n - ((2ULL << level) - 1);
            rmq_[level] = std::make_unique<Index[]>(level_size);

            size_t const interval_size = 1ULL << level;
            for(Index i = 0; i < level_size; i++) {
                auto const a = rmq_[level-1][i];
                auto const b = rmq_[level-1][i+interval_size];
                rmq_[level][i] = leq(data[a], data[b]) ? a : b;
            }
        }
    }

    RMQBenderFarachColton(RMQBenderFarachColton&&) = default;
    RMQBenderFarachColton& operator=(RMQBenderFarachColton&&) = default;

    RMQBenderFarachColton(RMQBenderFarachColton const&) = delete;
    RMQBenderFarachColton& operator=(RMQBenderFarachColton const&) = delete;

    size_t rmq(size_t const i, size_t const j) const {
        assert(i <= j);
        if(i == j)[[unlikely]] return i;

        auto const d = j - i + 1;
        auto const level = std::bit_width(d) - 1;
        auto const interval_size = 1ULL << level;

        auto const a = rmq_[level-1][i];
        auto const b = rmq_[level-1][j + 1 - interval_size];
        return leq(data_[a], data_[b]) ? a : b;
    }

    size_t operator()(size_t const i, size_t const j) const {
        return rmq(i, j);
    }
};

}

#endif
