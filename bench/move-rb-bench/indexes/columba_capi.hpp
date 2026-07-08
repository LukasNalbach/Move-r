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

#include <cstddef>
#include <cstdint>

// Thin C ABI used to drive one columba flavor from the (single) benchmark binary. columba's two flavors
// (run-length-compressed = b-move, and the FM-index) are selected by a compile-time macro and share global C++
// symbol names, so they cannot be linked into one executable directly. Instead each flavor is built as its own
// shared library that exposes ONLY this C ABI (its internal columba symbols are hidden), so the benchmark can
// link both libraries and measure both flavors in one process without a symbol clash. All heavy work (index
// loading and search) happens behind the ABI in the flavor library; timing, memory measurement (malloc_count)
// and output stay in the benchmark binary.
extern "C" {

struct columba_api_t {
    const char* flavor;      // RESULT algo suffix: "columba_rlc" (RLC / b-move) or "columba" (FM-index)
    const char* base_suffix; // appended to <index_dir>/<text_name> to form this flavor's index base name, so the
                             // two flavors' (flavor-specific) index files do not collide: "" for RLC, ".fm" for FM
    int max_k_hamming;       // maximum k the search scheme supports for hamming distance
    int max_k_edit;          // maximum k for edit distance

    // loads the index with base name @p base and search scheme @p scheme (k-mer table disabled); returns an
    // opaque handle, or nullptr if the index could not be loaded (missing files)
    void* (*load)(const char* base, const char* scheme);
    // switches the search scheme (e.g. minU for k<=7, columba for higher k); strategies are rebuilt lazily
    void (*set_scheme)(void* handle, const char* scheme);
    // length of the indexed text (excluding the sentinel)
    uint64_t (*input_size)(void* handle);

    // precomputes the (forward-strand) read bundles for @p count patterns (pointers + lengths), so bundle setup
    // is not counted in the per-pattern locate time. Returns an opaque bundle-set handle
    void* (*make_bundles)(void* handle, const char* const* patterns, const uint64_t* lengths, uint64_t count);
    // locates pattern @p i of the bundle set with at most @p k errors; metric_edit: 1 = edit, 0 = hamming;
    // want_cigar != 0 additionally builds each occurrence's CIGAR. Returns the occurrence count and writes the
    // occurrence-vector footprint to *occ_bytes and the CIGAR-data footprint to *cigar_bytes (both in bytes)
    uint64_t (*locate)(void* handle, void* bundles, uint64_t i, int metric_edit, int k, int want_cigar,
                       uint64_t* occ_bytes, uint64_t* cigar_bytes);
    // in-index match phase only (no locate) for edit distance; returns the undeduplicated FM-occurrence count
    uint64_t (*match_only_edit)(void* handle, void* bundles, uint64_t i, int k);

    // locates pattern @p pat (length m) with <=k errors by driving move_r's OWN apm algorithm on the columba
    // index (via columba_apm_adapter), then deduplicates exactly like move-rb. metric_edit: 1 = edit, 0 = hamming.
    // Returns the occurrence count and writes the occurrence-vector footprint to *occ_bytes. This measures the same
    // algorithm move_rb uses, but on columba's index instead of move_rb's.
    uint64_t (*locate_apm)(void* handle, const char* pat, uint64_t m, int metric_edit, int k, uint64_t* occ_bytes);

    // counts pattern @p pat (length m) with <=k HAMMING errors by driving move_r's OWN apm count algorithm on the
    // columba index (via columba_apm_adapter). columba has no cheap deduplicated count of its own, so this measures
    // the same count algorithm move_rb uses (SA-interval-width sum, no locate), but on columba's index. Returns the
    // (deduplicated) occurrence count. Hamming distance only (edit-distance count is unsupported).
    uint64_t (*count_apm)(void* handle, const char* pat, uint64_t m, int k);

    // ############################# bidirectional-extension benchmark (move-rb-bench-ext) #############################

    // builds a reusable bidirectional-extension adapter (columba_apm_adapter) over the loaded index; kept out of
    // the per-query timing. Returns an opaque handle to be passed to ext_run and freed with ext_free
    void* (*ext_make)(void* handle);
    // exact bidirectional search of pat[0..m) starting at position @p start and extending in the order given by
    // @p order (m-1 entries, 0 = prepend the next character on the left, 1 = append it on the right). want_locate:
    // 0 = count only (returns the occurrence count), 1 = also locate all occurrences (returns the count and writes
    // the occurrence-vector footprint to *occ_bytes)
    uint64_t (*ext_run)(void* ext, const char* pat, uint64_t m, uint64_t start, const uint8_t* order,
                        int want_locate, uint64_t* occ_bytes);
    void (*ext_free)(void* ext);

    void (*free_bundles)(void* bundles);
    void (*destroy)(void* handle);
};

// exported by libcolumba_rlc_plugin / libcolumba_fm_plugin respectively (each library defines only its own)
const columba_api_t* columba_rlc_api();
const columba_api_t* columba_fm_api();

} // extern "C"
