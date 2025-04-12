/**
 * rmq/rmq.hpp
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

#ifndef _RMQ_RMQ_HPP
#define _RMQ_RMQ_HPP

#include "rmq_bender_farach_colton.hpp"

#include <algorithm>

namespace rmq {

template<std::totally_ordered Value, size_t block_size_ = 64, std::unsigned_integral Index = uint32_t, bool minimum_ = true>
class RMQ {
private:
    static constexpr bool leq(Value const& a, Value const& b) { return compare<minimum_>(a, b); }
    static constexpr bool lt(Value const& a, Value const& b) { return compare_strict<minimum_>(a, b); }

    Value const* data_;
    std::unique_ptr<Index[]> block_min_pos_;
    std::unique_ptr<Value[]> block_min_;
    RMQBenderFarachColton<Value, Index, minimum_> block_rmq_;

    Index naive_rmq(Index beg, Index const end, Value& min) const {
        Index min_pos = beg;
        min = data_[beg];
        
        ++beg;
        while(beg < end) {
            if(lt(data_[beg], min)) {
                min = data_[beg];
                min_pos = beg;
            }
            ++beg;
        }
        return min_pos;
    }

public:
    RMQ() : data_(nullptr) {
    }

    RMQ(Value const* data, size_t const n) : data_(data) {
        size_t const num_blocks = (n-1) / block_size_ + 1;

        block_min_pos_ = std::make_unique<Index[]>(num_blocks);
        block_min_ = std::make_unique<Value[]>(num_blocks);

        for(size_t block = 0; block < num_blocks; block++) {
            Index const beg = block * block_size_;
            Index const end = std::min(beg + block_size_, n);

            Value min;
            auto const min_pos = naive_rmq(beg, end, min);

            block_min_[block] = min;
            block_min_pos_[block] = min_pos;
        }

        block_rmq_ = RMQBenderFarachColton<Value, Index, minimum_>(block_min_.get(), num_blocks);
    }

    RMQ(RMQ&&) = default;
    RMQ& operator=(RMQ&&) = default;

    RMQ(RMQ const&) = delete;
    RMQ& operator=(RMQ const&) = delete;

    Index rmq(Index const i, Index const j, Value& min) const {
        assert(i <= j);
        if(i == j)[[unlikely]] return i;

        // do naive scan for short intervals
        if(j - i <= 3 * block_size_) {
            return naive_rmq(i, j+1, min);
        }

        // determine minimum in block of i
        Index const left_end = (i / block_size_ + 1) * block_size_;
        Value left_min;
        Index const left_min_pos = naive_rmq(i, left_end, left_min);
        
        // determine minimum in block of j
        size_t const right_beg = (j / block_size_) * block_size_;
        Value right_min;
        Index const right_min_pos = naive_rmq(right_beg, j+1, right_min);

        // determine minimum in the between
        size_t const left_block = (i / block_size_) + 1;
        size_t const right_block = (j / block_size_) - 1;
        assert(left_block < right_block); // nb: we already ruled out short intervals
        Index const mid_min_pos = block_min_pos_[block_rmq_(left_block, right_block)];
        Value const mid_min = data_[mid_min_pos];

        // determine minimum of the three
        Index min_pos = left_min_pos;
        min = left_min;

        if(lt(mid_min, min)) {
            min = mid_min;
            min_pos = mid_min_pos;
        }

        if(lt(right_min, min)) {
            min = right_min;
            min_pos = right_min_pos;
        }
        
        return min_pos;
    }

    Index rmq(Index const i, Index const j) const {
        Value discard;
        return rmq(i, j, discard);
    }

    Index operator()(Index const i, Index const j) const {
        return rmq(i, j);
    }
};

}

#endif
