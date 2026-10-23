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
 * @file value_writer.hpp
 * @brief writes fixed-width little-endian values to a file, several threads at different positions at the same time
 */

#ifndef MOVE_R_MISC_VALUE_WRITER_HPP
#define MOVE_R_MISC_VALUE_WRITER_HPP

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <utility>

#include <misc/files.hpp>
#include <misc/memory.hpp>
#include <misc/platform.hpp>

/**
 * @brief buffers values of bytes_per_value bytes (little endian) and writes them to consecutive positions of a file
 *        with O_DIRECT (see direct_positioned_file); objects writing to disjoint ranges of the file can be used by
 *        different threads (each object takes whole cache lines)
 * @tparam bytes_per_value the number of bytes of a value (at most 8)
 */
template <uint8_t bytes_per_value>
class alignas(cache_line_bytes) value_writer {
protected:
    static constexpr uint64_t sector = direct_positioned_file::sector;
    // default number of values per buffer: one huge page of values
    static constexpr uint64_t buffer_values = huge_page_bytes / bytes_per_value;

    direct_positioned_file* file; // the file
    uint64_t offset; // the offset in bytes of the value 0 in the file
    uint64_t next; // the position of the first buffered value
    uint64_t count = 0; // the number of buffered values
    uint64_t capacity; // the number of values per buffer
    /* the buffer, capacity * bytes_per_value + sector bytes, aligned to sector; the buffered values start at
       values = buffer + (the offset of the value next in the file) % sector, so the sectors of the file are aligned in
       it: [0..count-1] the value at position next + i is stored in values[i * bytes_per_value..(i + 1) *
       bytes_per_value - 1], little endian */
    char* buffer;
    char* values;

    // moves the buffered values to the position i (precondition: count = 0)
    void start_at(uint64_t i)
    {
        next = i;
        values = buffer + (offset + next * bytes_per_value) % sector;
    }

public:
    /**
     * @brief starts writing values at a position
     * @param file the file
     * @param offset the offset in bytes of the value 0 in the file
     * @param first the position of the first value written
     * @param max_values the maximum number of values per buffer
     */
    value_writer(direct_positioned_file& file, uint64_t offset, uint64_t first, uint64_t max_values = buffer_values)
        : file(&file), offset(offset), capacity(std::clamp<uint64_t>(max_values, 1, buffer_values)),
          buffer((char*) allocate_aligned(capacity * bytes_per_value + sector))
    {
        start_at(first);
    }

    value_writer(value_writer&& other) noexcept
        : file(other.file), offset(other.offset), next(other.next), count(std::exchange(other.count, 0)),
          capacity(other.capacity), buffer(std::exchange(other.buffer, nullptr)), values(other.values)
    {
    }

    ~value_writer()
    {
        flush();
        if (buffer != nullptr) platform::free_aligned(buffer);
    }

    value_writer(const value_writer&) = delete;
    value_writer& operator=(const value_writer&) = delete;
    value_writer& operator=(value_writer&&) = delete;

    // appends a value
    inline void push(uint64_t value)
    {
        std::memcpy(values + count * bytes_per_value, &value, bytes_per_value);
        if (++count == capacity) flush();
    }

    // writes a value to the position i (consecutive positions are buffered)
    inline void put(uint64_t i, uint64_t value)
    {
        if (i != next + count) {
            flush();
            start_at(i);
        }

        push(value);
    }

    // writes the buffered values to the file
    void flush()
    {
        if (count == 0) return;
        file->write(values, count * bytes_per_value, offset + next * bytes_per_value);
        next += count;
        count = 0;
        start_at(next);
    }
};

/**
 * @brief buffers values of bytes_per_value bytes (little endian) that are reported from right to left and writes them
 *        to consecutive positions of a file with O_DIRECT (see direct_positioned_file); objects writing to disjoint
 *        ranges of the file can be used by different threads (each object takes whole cache lines)
 * @tparam bytes_per_value the number of bytes of a value (at most 8)
 */
template <uint8_t bytes_per_value>
class alignas(cache_line_bytes) value_writer_backward {
protected:
    static constexpr uint64_t sector = direct_positioned_file::sector;
    // default number of values per buffer: one huge page of values
    static constexpr uint64_t buffer_values = huge_page_bytes / bytes_per_value;

    direct_positioned_file* file; // the file
    uint64_t offset; // the offset in bytes of the value 0 in the file
    uint64_t capacity; // the number of values per buffer
    // size in bytes of the part of the buffer the values end in: capacity values, rounded up to whole sectors
    uint64_t values_bytes;
    uint64_t end = 0; // the position after the last buffered value
    uint64_t count = 0; // the number of buffered values
    /* the buffer, values_bytes + sector bytes, aligned to sector; the buffered values end at values_end = buffer +
       values_bytes + (the offset of the value end in the file) % sector, so the sectors of the file are aligned in it:
       [0..count-1] the value at position end - count + i is stored in the bytes values_end[-(count - i) *
       bytes_per_value..-(count - i - 1) * bytes_per_value - 1] */
    char* buffer;
    char* values_end;

    // moves the buffered values to end at the position i (precondition: count = 0)
    void end_at(uint64_t i)
    {
        end = i;
        values_end = buffer + values_bytes + (offset + end * bytes_per_value) % sector;
    }

public:
    /**
     * @brief constructs a writer
     * @param file the file
     * @param offset the offset in bytes of the value 0 in the file
     * @param max_values the maximum number of values per buffer
     */
    value_writer_backward(direct_positioned_file& file, uint64_t offset, uint64_t max_values = buffer_values)
        : file(&file), offset(offset), capacity(std::clamp<uint64_t>(max_values, 1, buffer_values)),
          values_bytes((capacity * bytes_per_value + sector - 1) / sector * sector),
          buffer((char*) allocate_aligned(values_bytes + sector))
    {
        end_at(0);
    }

    value_writer_backward(value_writer_backward&& other) noexcept
        : file(other.file), offset(other.offset), capacity(other.capacity), values_bytes(other.values_bytes),
          end(other.end), count(std::exchange(other.count, 0)),
          buffer(std::exchange(other.buffer, nullptr)), values_end(other.values_end)
    {
    }

    value_writer_backward(const value_writer_backward&) = delete;
    value_writer_backward& operator=(const value_writer_backward&) = delete;
    value_writer_backward& operator=(value_writer_backward&&) = delete;

    ~value_writer_backward()
    {
        flush();
        if (buffer != nullptr) platform::free_aligned(buffer);
    }

    // writes a value to the position i (consecutive positions from right to left are buffered)
    inline void put(uint64_t i, uint64_t value)
    {
        if (i + 1 != end - count) {
            flush();
            end_at(i + 1);
        }

        count++;
        std::memcpy(values_end - count * bytes_per_value, &value, bytes_per_value);
        if (count == capacity) flush();
    }

    // writes the buffered values to the file
    void flush()
    {
        if (count == 0) return;
        file->write(values_end - count * bytes_per_value, count * bytes_per_value,
            offset + (end - count) * bytes_per_value);
        const uint64_t first = end - count;
        count = 0;
        end_at(first);
    }
};

#endif // MOVE_R_MISC_VALUE_WRITER_HPP
