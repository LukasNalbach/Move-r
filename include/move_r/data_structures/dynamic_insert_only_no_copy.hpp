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