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

#ifndef MOVE_R_BENCH_MOVE_RB_BENCH_INDEXES_COLUMBA_CAPI_HPP
#define MOVE_R_BENCH_MOVE_RB_BENCH_INDEXES_COLUMBA_CAPI_HPP

#include <cstddef>
#include <cstdint>

/* The C interface of a columba plugin. The two variants of columba (b-move and FM-index) are selected by a macro at
   compile time and define the same C++ symbols, so each variant is built into a shared library (a plugin) whose only
   exported symbol is its getter of columba_api_t (columba_rlc_api or columba_fm_api); the benchmarks load the plugins
   with dlopen (see with_columba_plugin). The times and the memory are measured in the benchmarks. */
extern "C" {

struct columba_api_t {
    const char* flavor; // name of the variant in the RESULT lines: columba_rlc (b-move) or columba (FM-index)
    int32_t max_k_hamming; // maximum number of errors of the search schemes for the hamming distance
    int32_t max_k_edit; // maximum number of errors of the search schemes for the edit distance

    // loads the index with base name base and the search scheme scheme; returns a handle, or nullptr if the index
    // cannot be loaded
    void* (*load)(const char* base, const char* scheme);
    // sets the search scheme
    void (*set_scheme)(void* handle, const char* scheme);
    // returns the length of the text without the terminator
    uint64_t (*input_size)(void* handle);

    // creates the read bundles of count patterns (patterns[i] of length lengths[i]); returns a handle of the bundles
    void* (*make_bundles)(void* handle, const char* const* patterns, const uint64_t* lengths, uint64_t count);
    // locates the pattern of bundle i with at most k errors (metric_edit: 1 = edit, 0 = hamming; with_cigar: 1 =
    // compute a CIGAR string for each occurrence) with the algorithm of columba; returns the number of occurrences
    // and writes the size in bytes of the occurrences to *occ_bytes and of the CIGAR strings to *cigar_bytes
    uint64_t (*locate)(void* handle, void* bundles, uint64_t i, int32_t metric_edit, int32_t k, int32_t with_cigar,
                       uint64_t* occ_bytes, uint64_t* cigar_bytes);

    // locates pat[0..m-1] with at most k errors with the algorithm of move_rb (via columba_apm_adapter) and removes
    // the duplicates as move-rb locate does; returns the number of occurrences and writes their size in bytes to
    // *occ_bytes
    uint64_t (*locate_apm)(void* handle, const char* pat, uint64_t m, int32_t metric_edit, int32_t k, uint64_t* occ_bytes);
    // returns the number of occurrences of pat[0..m-1] with at most k mismatches, computed with the algorithm of
    // move_rb (via columba_apm_adapter)
    uint64_t (*count_apm)(void* handle, const char* pat, uint64_t m, int32_t k);

    // creates a columba_apm_adapter for the exact bidirectional search of move-rb-bench-ext; returns its handle
    void* (*ext_make)(void* handle);
    // matches pat[0..m-1] starting at position start and adding the other symbols in the directions of
    // order[0..m-2] (0 = prepend, 1 = append), see ext_query.hpp; with_locate = 0: returns the number of occurrences;
    // with_locate = 1: also locates them and writes their size in bytes to *occ_bytes
    uint64_t (*ext_run)(void* ext, const char* pat, uint64_t m, uint64_t start, const uint8_t* order,
                        int32_t with_locate, uint64_t* occ_bytes);
    // frees the adapter of ext_make
    void (*ext_free)(void* ext);

    // frees the read bundles of make_bundles
    void (*free_bundles)(void* bundles);
    // frees the index of load
    void (*destroy)(void* handle);
};

// the getters of the plugins (each plugin defines one of them)
const columba_api_t* columba_rlc_api();
const columba_api_t* columba_fm_api();

} // extern "C"

#endif // MOVE_R_BENCH_MOVE_RB_BENCH_INDEXES_COLUMBA_CAPI_HPP
