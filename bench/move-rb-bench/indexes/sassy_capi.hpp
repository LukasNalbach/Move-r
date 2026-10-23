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

#ifndef MOVE_R_BENCH_MOVE_RB_BENCH_INDEXES_SASSY_CAPI_HPP
#define MOVE_R_BENCH_MOVE_RB_BENCH_INDEXES_SASSY_CAPI_HPP

#include <algorithm>
#include <cstdint>
#include <dlfcn.h>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>

#include <misc/files.hpp>
#include <misc/utils.hpp>

// ############################# the C ABI of bench/move-rb-bench/indexes/sassy-capi #############################

extern "C" {
    struct sassy_capi_handle;

    // creates a searcher for the profile "ascii", "dna" or "iupac" (returns null for any other profile);
    // trace: the matches are computed with an alignment (a CIGAR), all: all end positions with at most
    // k edits are reported (else only the rightmost position of each local minimum)
    using sassy_capi_new_t = sassy_capi_handle* (*)(const char* profile, bool trace, bool all);
    using sassy_capi_free_t = void (*)(sassy_capi_handle*);
    // searches pattern in text with at most k edits and returns the number of matches whose end
    // position is >= min_end; if out_bytes is not null, the size in bytes of sassy's match vector
    // is written to it
    using sassy_capi_search_t = uint64_t (*)(sassy_capi_handle*, const uint8_t* pattern, uint64_t m,
        const uint8_t* text, uint64_t len, uint64_t k, uint64_t min_end, uint64_t* out_bytes);
}

/**
 * @brief the shared library libsassy_capi.so (C interface to sassy, an approximate pattern matcher
 *        without an index, that supports only the edit distance)
 */
class sassy_library {
public:
    /**
     * @brief loads libsassy_capi.so
     * @param directory directory containing libsassy_capi.so
     * @return whether the library has been loaded
     */
    bool open(const std::string& directory)
    {
        if (lib != nullptr) return true;
        const std::string path = directory + "/libsassy_capi.so";
        lib = dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);

        if (lib == nullptr) {
            std::cerr << "warning: cannot load " << path << " (" << dlerror() << ")" << std::endl;
            return false;
        }

        new_searcher = reinterpret_cast<sassy_capi_new_t>(dlsym(lib, "sassy_capi_new"));
        free_searcher = reinterpret_cast<sassy_capi_free_t>(dlsym(lib, "sassy_capi_free"));
        search = reinterpret_cast<sassy_capi_search_t>(dlsym(lib, "sassy_capi_search"));

        if (new_searcher == nullptr || free_searcher == nullptr || search == nullptr) {
            std::cerr << "warning: " << path << " does not contain the sassy_capi symbols" << std::endl;
            dlclose(lib);
            lib = nullptr;
            return false;
        }

        return true;
    }

    // returns whether the library has been loaded
    bool is_open() const { return lib != nullptr; }

    sassy_capi_new_t new_searcher = nullptr; // sassy_capi_new
    sassy_capi_free_t free_searcher = nullptr; // sassy_capi_free
    sassy_capi_search_t search = nullptr; // sassy_capi_search

private:
    void* lib = nullptr; // handle of the library (it is never closed)
};

/**
 * @brief a sassy searcher for a text file; the text is read and searched in overlapping chunks
 */
class sassy_scanner {
public:
    // default size of a chunk in bytes, 256 MiB (the default of --sassy-chunk; measurement setup)
    static constexpr uint64_t default_chunk_size = uint64_t { 256 } << 20;

    /**
     * @brief creates a searcher for a text file
     * @param library the loaded sassy library
     * @param path path of the text file (plain symbols, as written by fasta-to-text)
     * @param profile "ascii", "dna" or "iupac"
     * @param trace whether the matches are computed with an alignment (a CIGAR)
     * @param all whether all end positions with at most k edits are reported (else only the rightmost
     *        position of each local minimum)
     * @param chunk_size size of a chunk in bytes
     */
    sassy_scanner(sassy_library& library, const std::string& path, const std::string& profile,
                  bool trace, bool all, uint64_t chunk_size = default_chunk_size)
        : library(library), path(path), chunk_size(chunk_size)
    {
        handle = library.new_searcher(profile.c_str(), trace, all);
        if (handle == nullptr) throw std::runtime_error("sassy: unknown profile " + profile);
        text_length = std::filesystem::file_size(path);
    }

    ~sassy_scanner() { if (handle != nullptr) library.free_searcher(handle); }

    sassy_scanner(const sassy_scanner&) = delete;
    sassy_scanner& operator=(const sassy_scanner&) = delete;

    // returns the length of the text in bytes
    uint64_t input_size() const { return text_length; }

    /**
     * @brief searches the whole text for the occurrences of a pattern with at most k edits. Consecutive
     *        chunks overlap by m + k - 1 bytes; while a chunk is searched, the next chunk is read by a
     *        second thread.
     * @param pattern the pattern
     * @param k maximum number of edits
     * @param time_search the time spent in sassy in ns (without reading the text) is added to it
     * @param peak_bytes is set to the maximum of its value and the sizes in bytes of the match vectors of
     *        the chunks
     * @return number of matches
     */
    uint64_t scan(const std::string& pattern, uint64_t k, uint64_t& time_search, uint64_t& peak_bytes)
    {
        const uint64_t m = pattern.size();
        const uint64_t overlap = m + k - 1; // (maximum length of a match) - 1

        // while the chunk in one buffer is searched, the next chunk is read into the other buffer
        chunk buffers[2];
        buffers[0].reader.open_file(path);
        buffers[1].reader.open_file(path);

        uint64_t matches = 0;
        uint64_t next_position = 0; // first position of the text that is not covered by a chunk yet
        uint64_t slot = 0; // index of the buffer of the current chunk

        // sets the range of the next chunk in buffers[buf]; returns false if the whole text is already covered
        auto next_chunk = [&](uint64_t buf) {
            if (next_position >= text_length) return false;
            chunk& c = buffers[buf];
            c.beg = next_position == 0 ? 0 : next_position - overlap;
            c.end = std::min(text_length, next_position + chunk_size);
            c.min_end = next_position - c.beg + 1;
            next_position = c.end;
            c.data = nullptr;
            return true;
        };

        auto read_chunk = [](chunk& c) { c.data = c.reader.read(c.beg, c.end); };

        if (!next_chunk(slot)) return 0;
        read_chunk(buffers[slot]);

        while (true) {
            chunk& current = buffers[slot];
            const uint64_t other = 1 - slot;

            // start reading the next chunk while this one is searched
            std::thread reader;
            const bool has_next = next_chunk(other);
            if (has_next) reader = std::thread([&]{ read_chunk(buffers[other]); });

            uint64_t bytes = 0;
            auto t1 = now();
            matches += library.search(handle,
                reinterpret_cast<const uint8_t*>(pattern.data()), m,
                current.data, current.end - current.beg, k, current.min_end, &bytes);
            time_search += time_diff_ns(t1);
            peak_bytes = std::max(peak_bytes, bytes);

            if (!has_next) break;
            reader.join();
            slot = other;
        }

        return matches;
    }

private:
    // a chunk T[beg..end-1] of the text
    struct chunk {
        direct_reader reader; // reads the chunk from the text file
        const uint8_t* data = nullptr; // the chunk
        // min_end: matches ending before position min_end of the chunk are not counted (they are
        // counted in the previous chunk)
        uint64_t beg = 0, end = 0, min_end = 0;
    };

    sassy_library& library;
    std::string path;
    sassy_capi_handle* handle = nullptr;
    uint64_t text_length = 0;
    uint64_t chunk_size = default_chunk_size;
};

#endif // MOVE_R_BENCH_MOVE_RB_BENCH_INDEXES_SASSY_CAPI_HPP
