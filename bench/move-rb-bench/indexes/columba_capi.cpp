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

// One columba flavor's implementation of the columba_capi.hpp C ABI. Compiled twice -- once per flavor (the
// linked columba flavor library sets the RUN_LENGTH_COMPRESSION macro) -- into two shared libraries whose only
// exported symbol is the flavor's api getter (everything else, including all of columba, is hidden), so both
// can be linked into one benchmark executable. See columba_capi.hpp.

#include "columba_adapter.hpp"
#include "columba_apm_adapter.hpp"     // drive move_r's apm on the columba index
#include "columba_capi.hpp"
#include "../ext_query.hpp"            // ext_count / ext_locate (bidirectional-extension benchmark)

#include <algorithm>
#include <misc/search_schemes.hpp> // min_u_scheme
#include <cstdio>
#include <string>
#include <vector>

namespace {

// a precomputed set of forward-strand read bundles (one per pattern)
struct bundle_set_t {
    std::vector<ReadBundle> bundles;
};

void* cn_load(const char* base, const char* scheme)
{
    try {
        return new columba_native(base, scheme);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "note: columba load failed: %s\n", e.what());
        return nullptr; // missing / unreadable index files -> benchmark skips this flavor
    } catch (...) {
        return nullptr;
    }
}

void cn_set_scheme(void* handle, const char* scheme) { static_cast<columba_native*>(handle)->set_scheme(scheme); }

uint64_t cn_input_size(void* handle) { return static_cast<columba_native*>(handle)->input_size(); }

void* cn_make_bundles(void* handle, const char* const* patterns, const uint64_t* lengths, uint64_t count)
{
    columba_native* index = static_cast<columba_native*>(handle);
    bundle_set_t* bs = new bundle_set_t();
    bs->bundles.reserve(count);
    for (uint64_t i = 0; i < count; i++) {
        std::string pattern(patterns[i], patterns[i] + lengths[i]);
        bs->bundles.emplace_back(index->make_bundle(pattern)); // ReadBundle owns a copy of the sequence
    }
    return bs;
}

uint64_t cn_locate(void* handle, void* bundles, uint64_t i, int metric_edit, int k, int want_cigar,
                   uint64_t* occ_bytes, uint64_t* cigar_bytes)
{
    columba_native* index = static_cast<columba_native*>(handle);
    ReadBundle& bundle = static_cast<bundle_set_t*>(bundles)->bundles[i];
    const DistanceMetric metric = metric_edit ? EDIT : HAMMING;

    if (want_cigar) {
        std::vector<TextOcc> occs = index->locate_occ_cigar(metric, bundle, (uint64_t) k);
        *occ_bytes = occs.size() * sizeof(TextOcc);
        uint64_t cig = 0;
        for (TextOcc& o : occs) cig += o.getCigar().size();
        *cigar_bytes = cig;
        return occs.size();
    } else {
        std::vector<TextOcc> occs = index->locate_occ(metric, bundle, (uint64_t) k);
        *occ_bytes = occs.size() * sizeof(TextOcc);
        *cigar_bytes = 0;
        return occs.size();
    }
}

uint64_t cn_match_only_edit(void* handle, void* bundles, uint64_t i, int k)
{
    columba_native* index = static_cast<columba_native*>(handle);
    ReadBundle& bundle = static_cast<bundle_set_t*>(bundles)->bundles[i];
    return static_cast<uint64_t>(index->match_only(EDIT, bundle, (uint64_t) k));
}

uint64_t cn_locate_apm(void* handle, const char* pat, uint64_t m, int metric_edit, int k, uint64_t* occ_bytes)
{
    columba_native* native = static_cast<columba_native*>(handle);
    columba_apm_adapter adapter(native->get_index());
    std::string P(pat, pat + m);
    const search_scheme_t scheme = min_u_scheme((uint8_t) k); // same scheme move-rb uses in the benchmark

    std::vector<aprx_occ_t<uint64_t>> occ;
    auto collect = [&](aprx_occ_t<uint64_t> o) { occ.emplace_back(o); };
    if (metric_edit) {
        adapter.locate<EDIT_DISTANCE>(P, scheme, collect);
        // deduplicate exactly like move-rb (sort by position, then the coverage-preserving filter)
        std::sort(occ.begin(), occ.end(), [](const aprx_occ_t<uint64_t>& a, const aprx_occ_t<uint64_t>& b) { return a.pos < b.pos; });
        filter_edit_distance_occurrences<uint64_t>(occ, (uint64_t) k);
    } else {
        adapter.locate<HAMMING_DISTANCE>(P, scheme, collect);
    }
    *occ_bytes = occ.size() * sizeof(aprx_occ_t<uint64_t>);
    return occ.size();
}

// builds a columba_apm_adapter over the loaded index once, so the (per-query-timed) ext_run calls do not pay the
// adapter's construction cost. The adapter holds a reference to the index, which outlives it (freed via ext_free)
void* cn_ext_make(void* handle)
{
    columba_native* native = static_cast<columba_native*>(handle);
    return new columba_apm_adapter(native->get_index());
}

uint64_t cn_ext_run(void* ext, const char* pat, uint64_t m, uint64_t start, const uint8_t* order,
                    int want_locate, uint64_t* occ_bytes)
{
    columba_apm_adapter* adapter = static_cast<columba_apm_adapter*>(ext);
    std::string P(pat, pat + m);

    ext_plan_t plan;
    plan.start = start;
    plan.order.resize(m == 0 ? 0 : m - 1);
    for (uint64_t i = 0; i + 1 < m; i++) plan.order[i] = order[i] == 0 ? LEFT : RIGHT;

    if (!want_locate) return ext_count(*adapter, P, plan);

    std::vector<uint64_t> occ;
    ext_locate(*adapter, P, plan, [&](uint64_t o) { occ.push_back(o); });
    if (occ_bytes) *occ_bytes = occ.size() * sizeof(uint64_t);
    return occ.size();
}

void cn_ext_free(void* ext) { delete static_cast<columba_apm_adapter*>(ext); }

void cn_free_bundles(void* bundles) { delete static_cast<bundle_set_t*>(bundles); }
void cn_destroy(void* handle) { delete static_cast<columba_native*>(handle); }

const columba_api_t g_api = {
    columba_native::flavor_name,
#ifdef RUN_LENGTH_COMPRESSION
    "",    // RLC index uses the plain base name <index_dir>/<text_name>
#else
    ".fm", // FM index uses <index_dir>/<text_name>.fm so it does not collide with the RLC index
#endif
    MAX_K,      // hamming
    MAX_K_EDIT, // edit
    cn_load, cn_set_scheme, cn_input_size,
    cn_make_bundles, cn_locate, cn_match_only_edit, cn_locate_apm,
    cn_ext_make, cn_ext_run, cn_ext_free,
    cn_free_bundles, cn_destroy
};

} // namespace

// the flavor-unique, exported api getter (the only default-visibility symbol of the shared library)
#ifdef RUN_LENGTH_COMPRESSION
#define COLUMBA_API_GETTER columba_rlc_api
#else
#define COLUMBA_API_GETTER columba_fm_api
#endif

extern "C" __attribute__((visibility("default"))) const columba_api_t* COLUMBA_API_GETTER() { return &g_api; }
