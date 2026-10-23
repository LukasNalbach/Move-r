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

#ifndef MOVE_R_BENCH_MOVE_RB_BENCH_INDEXES_COLUMBA_ADAPTER_HPP
#define MOVE_R_BENCH_MOVE_RB_BENCH_INDEXES_COLUMBA_ADAPTER_HPP

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

// ALPHABET is set to 5 by the columba targets (see CMakeLists.txt); the fallback is for tools that parse this header
// on its own
#ifndef ALPHABET
#define ALPHABET 5
#endif

// the variant of columba is selected by RUN_LENGTH_COMPRESSION (set by the linked columba library): the run-length
// compressed bidirectional move data structure (b-move) or the bidirectional FM-index
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
 * @brief the approximate pattern matching of columba itself (b-move or FM-index, see RUN_LENGTH_COMPRESSION); the
 *        reverse complement of a pattern is not searched and no SAM lines are generated
 */
class columba_native {
  public:
#ifdef RUN_LENGTH_COMPRESSION
    using columba_index_t = BMove;
    static constexpr const char* flavor_name = "columba_rlc"; // name of the variant in the RESULT lines
#else
    using columba_index_t = FMIndex;
    static constexpr const char* flavor_name = "columba"; // name of the variant in the RESULT lines
#endif

  protected:
    std::unique_ptr<columba_index_t> _index; // the columba index
    std::unique_ptr<SearchStrategy> _strategy_hamming; // the search strategy for the hamming distance (see strategy)
    std::unique_ptr<SearchStrategy> _strategy_edit; // the search strategy for the edit distance (see strategy)
    std::string _scheme; // name of columba's search scheme (e.g. columba, minU, pigeon)
    std::string _id = ">p"; // the name of each pattern

    // returns a search strategy for a distance metric (all occurrences, single-end, uniform partitioning, only the
    // pattern itself, no SAM lines)
    std::unique_ptr<SearchStrategy> make_strategy(DistanceMetric metric)
    {
        Parameters params;
        params.metric = metric;
        params.mMode = ALL;
        params.sMode = SINGLE_END;
        params.searchScheme = _scheme;
        params.pStrategy = UNIFORM;

        std::unique_ptr<SearchStrategy> strategy = params.createStrategy(*_index);
        strategy->setSearchReverseComplement(false);
        strategy->setGenerateSAMLines(false);
        return strategy;
    }

    // returns the search strategy for a distance metric (created at its first use)
    SearchStrategy& strategy(DistanceMetric metric)
    {
        std::unique_ptr<SearchStrategy>& s = metric == EDIT ? _strategy_edit : _strategy_hamming;
        if (!s) s = make_strategy(metric);
        return *s;
    }

  public:
    /**
     * @brief loads a columba index (with CIGAR support and without the k-mer table: wordSize = 0)
     * @param base base name of the index files
     * @param scheme name of columba's search scheme (columba, minU, pigeon, ...)
     */
    columba_native(const std::string& base, const std::string& scheme)
        : _scheme(scheme)
    {
        Parameters params; // the default parameters of the index
        const length_t word_size = 0;
#ifdef RUN_LENGTH_COMPRESSION
        _index = std::make_unique<BMove>(base, false /* verbose */, false /* noCIGAR */, word_size);
#else
        _index = std::make_unique<FMIndex>(
            base, params.inTextVerificationPoint, false /* noCIGAR */,
            params.sparsenessFactor, false /* verbose */, word_size);
#endif
    }

    // sets the search scheme (the search strategies are created again at their next use)
    void set_scheme(const std::string& scheme)
    {
        if (scheme == _scheme) return;
        _scheme = scheme;
        _strategy_hamming.reset();
        _strategy_edit.reset();
    }

    // returns the columba index
    inline columba_index_t& get_index() const { return *_index; }

    // returns the length of the text without the terminator
    inline uint64_t input_size() const { return _index->getTextSize() - 1; }

    // returns the length of the text including the terminator
    inline uint64_t text_size() const { return _index->getTextSize(); }

    // returns a read bundle of a pattern
    inline ReadBundle make_bundle(const std::string& P) { return ReadBundle(Read(_id, P, "*")); }

    /**
     * @brief returns the occurrences in T of the pattern of a read bundle with at most k errors
     * @param metric HAMMING or EDIT
     * @param bundle the read bundle of the pattern
     * @param k the maximum number of errors
     * @param with_cigar whether a CIGAR string is computed for each occurrence (generateCIGARS)
     */
    inline std::vector<TextOcc> locate_occ(DistanceMetric metric, ReadBundle& bundle, uint64_t k, bool with_cigar)
    {
        Counters counters;
        std::vector<TextOcc> result;
        strategy(metric).matchApprox(bundle, (length_t) k, counters, result);
        if (with_cigar) _index->generateCIGARS(result, counters, bundle);
        return result;
    }
};

#endif // MOVE_R_BENCH_MOVE_RB_BENCH_INDEXES_COLUMBA_ADAPTER_HPP
