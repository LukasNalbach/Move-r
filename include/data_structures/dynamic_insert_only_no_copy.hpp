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

/**
 * @brief dynamically-growing insert-only no-copy data structure
 * @tparam T value type
 */
template <typename T>
class dynamic_insert_only_no_copy {
protected:
    std::vector<std::vector<T>> vectors; // vectors that store the elements

public:
    dynamic_insert_only_no_copy() = default;

    /**
     * @brief creates an empty data structure with a certain amount of elements reserved
     * @param size initially reserved number of elements
     */
    dynamic_insert_only_no_copy(uint64_t size)
    {
        vectors.resize(1);
        vectors.back().reserve(size);
    }

    /**
     * @brief clears all vectors
     */
    inline void clear()
    {
        vectors.clear();
        vectors.shrink_to_fit();
    }

    /**
     * @brief inserts an element into the data structure, doubles the number of elements reserved
     *        reserved by the data structure if it is full
     * @param v element
     * @return pointer to the element in the data structure
     */
    inline T* emplace_back(T&& v)
    {
        if (vectors.back().size() == vectors.back().capacity()) {
            size_t new_capacity = 2 * vectors.back().capacity();
            vectors.emplace_back(std::vector<T>());
            vectors.back().reserve(new_capacity);
        }

        vectors.back().emplace_back(v);
        return &(vectors.back().back());
    }

    /**
     * @brief inserts an element into the data structure, doubles the number of elements reserved
     *        reserved by the data structure if it is full
     * @param v element
     * @return pointer to the element in the data structure
     */
    inline T* emplace_back(T& v)
    {
        return emplace_back(std::move(v));
    }
};