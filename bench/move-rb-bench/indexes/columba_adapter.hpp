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

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// columba is built in one of two flavors, selected at compile time by the
// RUN_LENGTH_COMPRESSION macro (set by the linked columba_rlc / columba_fm
// library): the run-length-compressed move structure (b-move) or the FM-index.
#ifdef RUN_LENGTH_COMPRESSION
#include <bmove/bmove.h>
#else
#include <fmindex/fmindex.h>
#endif
#include <definitions.h>
#include <indexhelpers.h>
#include <parameters/alignparameters.h>
#include <reads.h>
#include <searchstrategy.h>

/**
 * @brief Thin driver for columba's native lossless approximate-pattern-matching
 *        algorithm (Renders/Depuydt et al.), used to benchmark columba's own
 *        search against move_rb.
 *
 * columba is compiled in two flavors (chosen by the RUN_LENGTH_COMPRESSION
 * macro): the run-length-compressed bidirectional move structure (b-move, the
 * successor of the standalone b-move repository, now the RLC flavor of columba)
 * and the bidirectional FM-index. This driver wraps whichever flavor was linked
 * and exposes it through a small interface analogous to the (removed) b-move
 * driver.
 *
 * Like the previous b-move driver, it searches only the pattern itself and NOT
 * its reverse complement (setSearchReverseComplement(false)) and skips SAM-line
 * generation (setGenerateSAMLines(false)), so that only the matching is measured
 * -- comparable to move_rb, which is forward-strand only.
 */
class columba_native {
  public:
#ifdef RUN_LENGTH_COMPRESSION
    using columba_index_t = BMove; // run-length-compressed move structure (b-move)
    // name used in the RESULT lines to distinguish the flavor
    static constexpr const char* flavor_name = "columba_rlc";
#else
    using columba_index_t = FMIndex; // bidirectional FM-index
    static constexpr const char* flavor_name = "columba";
#endif

  protected:
    std::unique_ptr<columba_index_t> _index; // the wrapped columba index
    std::unique_ptr<SearchStrategy> _strategy_hamming; // created on first use
    std::unique_ptr<SearchStrategy> _strategy_edit; // created on first use
    std::string _scheme; // columba search scheme name (e.g. columba, minU, pigeon)
    std::string _id = ">p"; // a dummy sequence id reused for every pattern

    /**
     * @brief builds a search strategy for the given distance metric, configured
     *        for forward-strand-only all-mapping without SAM output
     */
    std::unique_ptr<SearchStrategy> make_strategy(DistanceMetric metric)
    {
        Parameters params;
        params.metric = metric;
        params.mMode = ALL;
        params.sMode = SINGLE_END;
        params.searchScheme = _scheme;
        params.pStrategy = UNIFORM; // match move-rb's uniform (equal-length) partitioning

        std::unique_ptr<SearchStrategy> strategy = params.createStrategy(*_index);
        strategy->setSearchReverseComplement(false);
        strategy->setGenerateSAMLines(false);
        return strategy;
    }

    /**
     * @brief returns the (lazily created) strategy for the given metric
     */
    SearchStrategy& strategy(DistanceMetric metric)
    {
        std::unique_ptr<SearchStrategy>& s =
            metric == EDIT ? _strategy_edit : _strategy_hamming;
        if (!s) s = make_strategy(metric);
        return *s;
    }

  public:
    /**
     * @brief loads the columba index from the given base path
     * @param base base name of the columba index files
     * @param scheme columba search scheme to use (columba, minU, pigeon, ...)
     *
     * The index is loaded with CIGAR generation enabled (noCIGAR = false) so
     * that the CIGAR-producing locate can be measured; plain locate never
     * requests the CIGAR (the edit filter does not compute one), so the plain
     * measurement is unaffected.
     *
     * The k-mer range table is DISABLED (wordSize = 0): populateTable then only
     * stores the full range for the empty string, so no exact-seed k-mer lookups
     * are used. This matches move_rb (which has no k-mer table) for a fair
     * comparison; results are identical, only the seeding shortcut is removed.
     */
    columba_native(const std::string& base, const std::string& scheme)
        : _scheme(scheme)
    {
        Parameters params; // for the index-construction numeric defaults
        const length_t word_size = 0; // 0 => k-mer range table disabled
#ifdef RUN_LENGTH_COMPRESSION
        _index = std::make_unique<BMove>(base, false /* verbose */,
                                         false /* noCIGAR */, word_size);
#else
        _index = std::make_unique<FMIndex>(
            base, params.inTextVerificationPoint, false /* noCIGAR */,
            params.sparsenessFactor, false /* verbose */, word_size);
#endif
    }

    /**
     * @brief switches the columba search scheme (e.g. minU for k <= 7, columba
     *        for higher k); the per-metric strategies are rebuilt lazily on next
     *        use. The index itself is not reloaded.
     */
    void set_scheme(const std::string& scheme)
    {
        if (scheme == _scheme) return;
        _scheme = scheme;
        _strategy_hamming.reset();
        _strategy_edit.reset();
    }

    /**
     * @brief returns the wrapped columba index (so move_r's apm can be driven on it via columba_apm_adapter)
     */
    inline columba_index_t& get_index() const { return *_index; }

    /**
     * @brief length of the indexed text (excluding the sentinel)
     */
    inline uint64_t input_size() const { return _index->getTextSize() - 1; }

    /**
     * @brief length of the indexed text including the sentinel
     */
    inline uint64_t text_size() const { return _index->getTextSize(); }

    /**
     * @brief wraps a pattern in a single-strand read bundle (its reverse
     *        complement is computed but never searched)
     */
    inline ReadBundle make_bundle(const std::string& P)
    {
        return ReadBundle(Read(_id, P, "*"));
    }

    /**
     * @brief locates the occurrences of the pattern in @p bundle with at most k
     *        errors w.r.t. the given distance metric, searching only the pattern
     *        (forward strand)
     * @param metric HAMMING or EDIT
     * @return the (native, deduplicated) in-text occurrences reported (no CIGARs)
     */
    inline std::vector<TextOcc> locate_occ(DistanceMetric metric, ReadBundle& bundle, uint64_t k)
    {
        Counters counters;
        std::vector<TextOcc> result;
        strategy(metric).matchApprox(bundle, (length_t) k, counters, result);
        return result;
    }

    /**
     * @brief locates the occurrences of the pattern in @p bundle with at most k
     *        errors w.r.t. the given distance metric, searching only the pattern
     *        (forward strand)
     * @param metric HAMMING or EDIT
     * @return the number of (native, deduplicated) occurrences reported
     */
    inline uint64_t locate(DistanceMetric metric, ReadBundle& bundle, uint64_t k)
    {
        return locate_occ(metric, bundle, k).size();
    }

    /**
     * @brief like locate_occ, but additionally computes each occurrence's CIGAR
     *        alignment (columba's generateCIGARS; for RLC from the matched
     *        string, for the FM-index from the located text) -- for measuring
     *        columba's CIGAR-producing locate against move_rb's
     * @param metric HAMMING or EDIT
     * @return the occurrences, each carrying a CIGAR string (occ.getCigar())
     */
    inline std::vector<TextOcc> locate_occ_cigar(DistanceMetric metric, ReadBundle& bundle, uint64_t k)
    {
        Counters counters;
        std::vector<TextOcc> result;
        SearchStrategy& s = strategy(metric);
        s.matchApprox(bundle, (length_t) k, counters, result);
        _index->generateCIGARS(result, counters, bundle);
        return result;
    }

    /**
     * @brief columba's MATCH phase only (matchApproxMatchOnly, no in-text
     *        locate): returns the number of in-index (FM) occurrences. For
     *        TIMING the match phase; the count is undeduplicated (do not use as
     *        a count). Forward strand only.
     */
    inline size_t match_only(DistanceMetric metric, ReadBundle& bundle, uint64_t k)
    {
        Counters counters;
        return strategy(metric).matchApproxMatchOnly(bundle, (length_t) k, counters);
    }

    // like match_only, but returns columba's NODE_COUNTER (number of index nodes visited during the match phase)
    inline uint64_t match_only_nodes(DistanceMetric metric, ReadBundle& bundle, uint64_t k)
    {
        Counters counters;
        strategy(metric).matchApproxMatchOnly(bundle, (length_t) k, counters);
        return counters.get(Counters::NODE_COUNTER);
    }
};
