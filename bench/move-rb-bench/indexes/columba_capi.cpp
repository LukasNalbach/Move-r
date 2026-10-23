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

// the C interface of columba_capi.hpp for one variant of columba (compiled once per variant and width of the positions)

#include <ips2ra.hpp>

#include "columba_adapter.hpp"
#include "columba_apm_adapter.hpp"
#include "columba_capi.hpp"
#include "../ext_query.hpp"

#include <algorithm>
#include <misc/search_schemes.hpp> // min_u_scheme
#include <cstdio>
#include <string>
#include <vector>

namespace {

// the read bundles of a pattern set: [0..num_patterns-1] the read bundle of the i-th pattern
struct bundle_set_t {
    std::vector<ReadBundle> bundles;
};

void* cn_load(const char* base, const char* scheme)
{
    try {
        return new columba_native(base, scheme);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "note: columba load failed: %s\n", e.what());
        return nullptr;
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
        bs->bundles.emplace_back(index->make_bundle(pattern));
    }
    return bs;
}

uint64_t cn_locate(void* handle, void* bundles, uint64_t i, int32_t metric_edit, int32_t k, int32_t with_cigar,
                   uint64_t* occ_bytes, uint64_t* cigar_bytes)
{
    columba_native* index = static_cast<columba_native*>(handle);
    ReadBundle& bundle = static_cast<bundle_set_t*>(bundles)->bundles[i];
    const DistanceMetric metric = metric_edit ? EDIT : HAMMING;

    std::vector<TextOcc> occs = index->locate_occ(metric, bundle, (uint64_t) k, with_cigar != 0);
    *occ_bytes = occs.size() * sizeof(TextOcc);
    uint64_t cig = 0;
    if (with_cigar) for (TextOcc& o : occs) cig += o.getCigar().size();
    *cigar_bytes = cig;
    return occs.size();
}

uint64_t cn_locate_apm(void* handle, const char* pat, uint64_t m, int32_t metric_edit, int32_t k, uint64_t* occ_bytes)
{
    columba_native* native = static_cast<columba_native*>(handle);
    columba_apm_adapter adapter(native->get_index());
    std::string P(pat, pat + m);
    const search_scheme_t scheme = min_u_scheme((uint8_t) k);

    std::vector<aprx_occ_t<>> occ;

    if (metric_edit) {
        adapter.locate(P, scheme, EDIT_DISTANCE, occ);
        // remove the duplicates as move-rb locate does
        ips2ra::sort(occ.begin(), occ.end(), [](const auto& occ) { return occ.pos(); });
        filter_edit_distance_occurrences(occ, (uint64_t) k);
    } else {
        adapter.locate(P, scheme, HAMMING_DISTANCE, occ);
    }
    *occ_bytes = occ.size() * sizeof(aprx_occ_t<>);
    return occ.size();
}

uint64_t cn_count_apm(void* handle, const char* pat, uint64_t m, int32_t k)
{
    columba_native* native = static_cast<columba_native*>(handle);
    columba_apm_adapter adapter(native->get_index());
    std::string P(pat, pat + m);
    const search_scheme_t scheme = min_u_scheme((uint8_t) k);
    return static_cast<uint64_t>(adapter.count_hamming_dist(P, scheme));
}

void* cn_ext_make(void* handle)
{
    columba_native* native = static_cast<columba_native*>(handle);
    return new columba_apm_adapter(native->get_index());
}

uint64_t cn_ext_run(void* ext, const char* pat, uint64_t m, uint64_t start, const uint8_t* order,
                    int32_t with_locate, uint64_t* occ_bytes)
{
    columba_apm_adapter* adapter = static_cast<columba_apm_adapter*>(ext);
    std::string P(pat, pat + m);

    ext_plan_t plan;
    plan.start = start;
    plan.order.resize(m == 0 ? 0 : m - 1);
    for (uint64_t i = 0; i + 1 < m; i++) plan.order[i] = order[i] == 0 ? LEFT : RIGHT;

    if (!with_locate) return ext_count(*adapter, P, plan);

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
    MAX_K,
    MAX_K_EDIT,
    cn_load, cn_set_scheme, cn_input_size,
    cn_make_bundles, cn_locate, cn_locate_apm, cn_count_apm,
    cn_ext_make, cn_ext_run, cn_ext_free,
    cn_free_bundles, cn_destroy
};

} // namespace

// the getter of the C interface (the only exported symbol of the plugin)
#ifdef RUN_LENGTH_COMPRESSION
#define COLUMBA_API_GETTER columba_rlc_api
#else
#define COLUMBA_API_GETTER columba_fm_api
#endif

extern "C" __attribute__((visibility("default"))) const columba_api_t* COLUMBA_API_GETTER() { return &g_api; }
