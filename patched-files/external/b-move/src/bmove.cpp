/******************************************************************************
 *  b-move: bidirectional move structure                                      *
 *  Copyright (C) 2020-2024 - Lore Depuydt <lore.depuydt@ugent.be> and        *
 *                            Luca Renders <luca.renders@ugent.be> and        *
 *                            Jan Fostier <jan.fostier@ugent.be>              *
 *                                                                            *
 *  This program is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU Affero General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Affero General Public License for more details.                       *
 *                                                                            *
 * You should have received a copy of the GNU Affero General Public License   *
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.     *
 ******************************************************************************/

#include "bmove.h"
#include "alphabet.h" // for Alphabet
#include "definitions.h"
#include "rindexhelpers.h"
#include "logger.h" // for Logger
#include "reads.h"

#include <algorithm>           // for max
#include <assert.h>            // for assert
#include <chrono>              // for duration, operator-, high_resolution_c...
#include <cstdint>             // for uint8_t
#include <ostream>             // for opera...
#include <sdsl/int_vector.hpp> // for int_v...
#include <stdexcept>           // for runti...
#include <string>              // for string
class MemoryMappedTextFile;
class Search;
class Substring;

using namespace std;

// ----------------------------------------------------------------------------
// PREPROCESSING ROUTINES
// ----------------------------------------------------------------------------

void BMove::fromFiles(const string& baseFile, bool verbose) {

    readTagCompAndCounts(baseFile, verbose);
    stringstream ss;
#ifndef LF_BENCHMARK_FUNCTIONALITY
    // Read BMove specific files
    if (verbose) {
        ss << "Reading " << baseFile << ".plcp...";
        logger.logInfo(ss);
    }
    if (!plcp.read(baseFile + ".plcp")) {
        throw runtime_error("Cannot open file: " + baseFile + ".plcp");
    }
#endif

#ifdef LF_MOVE_BIT_PACKED
    string LFextension = ".LFBP";
#else
    string LFextension = ".LFFull";
#endif

    if (verbose) {
        ss << "Reading " << baseFile << LFextension << "...";
        logger.logInfo(ss);
    }
    if (!move.load(baseFile)) {
        throw runtime_error("Error loading move file: " + baseFile +
                            LFextension);
    }

#ifndef LF_BENCHMARK_FUNCTIONALITY
    if (verbose) {
        ss << "Reading " << baseFile << ".smpf...";
        logger.logInfo(ss);
    }
    if (!readIntVector(baseFile + ".smpf", samplesFirst)) {
        throw runtime_error("Cannot open file: " + baseFile + ".smpf");
    }

    if (verbose) {
        ss << "Reading " << baseFile << ".smpl...";
        logger.logInfo(ss);
    }
    if (!readIntVector(baseFile + ".smpl", samplesLast)) {
        throw runtime_error("Cannot open file: " + baseFile + ".smpl");
    }

#ifdef PHI_MOVE
#ifdef PHI_MOVE_BALANCED
    string extension = ".balanced";
#else
    string extension = ".unbalanced";
#endif
    if (verbose) {
        ss << "Reading " << baseFile << ".move" << extension << ".prdf...";
        logger.logInfo(ss);
    }
    if (!predFirst.read(baseFile + ".move" + extension + ".prdf")) {
        throw runtime_error("Cannot open file: " + baseFile + ".move" +
                            extension + ".prdf");
    }

    if (verbose) {
        ss << "Reading " << baseFile << ".move" << extension << ".prdl...";
        logger.logInfo(ss);
    }
    if (!predLast.read(baseFile + ".move" + extension + ".prdl")) {
        throw runtime_error("Cannot open file: " + baseFile + ".move" +
                            extension + ".prdl");
    }
#else
    if (verbose) {
        ss << "Reading " << baseFile << ".og.prdf...";
        logger.logInfo(ss);
    }
    if (!predFirst.read(baseFile + ".og.prdf")) {
        throw runtime_error("Cannot open file: " + baseFile + ".og.prdf");
    }

    if (verbose) {
        ss << "Reading " << baseFile << ".og.prdl...";
        logger.logInfo(ss);
    }
    if (!predLast.read(baseFile + ".og.prdl")) {
        throw runtime_error("Cannot open file: " + baseFile + ".og.prdl");
    }
#endif

#ifndef PHI_MOVE
    if (verbose) {
        ss << "Reading " << baseFile << ".ftr...";
        logger.logInfo(ss);
    }
    if (!readIntVector(baseFile + ".ftr", firstToRun)) {
        throw runtime_error("Cannot open file: " + baseFile + ".ftr");
    }

    if (verbose) {
        ss << "Reading " << baseFile << ".ltr...";
        logger.logInfo(ss);
    }
    if (!readIntVector(baseFile + ".ltr", lastToRun)) {
        throw runtime_error("Cannot open file: " + baseFile + ".ltr");
    }
#else

#ifdef PHI_MOVE_BIT_PACKED
    extension = ".phiBP";
#else
    extension = ".phiFull";
#endif

#ifdef PHI_MOVE_BALANCED
    extension += ".balanced";
#else
    extension += ".unbalanced";
#endif

    if (verbose) {
        ss << "Reading " << baseFile << extension << "...";
        logger.logInfo(ss);
    }
    if (!phiMove.load(baseFile + extension)) {
        throw runtime_error("Cannot open file: " + baseFile + extension);
    }

    if (verbose) {
        ss << "Reading " << baseFile << extension << ".inv...";
        logger.logInfo(ss);
    }
    if (!phiInvMove.load(baseFile + extension + ".inv")) {
        throw runtime_error("Cannot open file: " + baseFile + extension +
                            ".inv");
    }
#endif
#endif

    if (verbose) {
        ss << "Reading " << baseFile << ".rev" << LFextension << "...";
        logger.logInfo(ss);
    }
    if (!moveR.load(baseFile + ".rev")) {
        throw runtime_error("Error loading reverse move file: " + baseFile +
                            ".rev" + LFextension);
    }

#ifndef LF_BENCHMARK_FUNCTIONALITY
    if (verbose) {
        ss << "Reading " << baseFile << ".rev.smpf...";
        logger.logInfo(ss);
    }
    if (!readIntVector(baseFile + ".rev.smpf", revSamplesFirst)) {
        throw runtime_error("Cannot open file: " + baseFile + ".rev.smpf");
    }

    if (verbose) {
        ss << "Reading " << baseFile << ".rev.smpl...";
        logger.logInfo(ss);
    }
    if (!readIntVector(baseFile + ".rev.smpl", revSamplesLast)) {
        throw runtime_error("Cannot open file: " + baseFile + ".rev.smpl");
    }
#endif

    readSequenceNamesAndPositions(baseFile, verbose);
}

// ----------------------------------------------------------------------------
// ROUTINES FOR ACCESSING DATA STRUCTURE
// ----------------------------------------------------------------------------

#ifndef PHI_MOVE
void BMove::phi(length_t& pos) const {
#ifdef PHI_BENCHMARK_FUNCTIONALITY
    // Start timing with rdtscp
    unsigned int start_aux;
    uint64_t start = __builtin_ia32_rdtscp(&start_aux);
#endif
    // Find the rank of the predecessor of pos in a circular manner.
    length_t predRank = predFirst.predecessorRankCircular(pos);
    // Select the predecessor position using its rank.
    length_t pred = predFirst.select(predRank);

    // Calculate the distance from the predecessor to the current position.
    length_t delta = pred < pos ? pos - pred : pos + 1;

    // Ensure that phi(SA[0]) is not called; the predecessor rank must be valid.
    assert(firstToRun[predRank] > 0);

    // Get the previous sample from the samplesLast array.
    length_t prev_sample = samplesLast[firstToRun[predRank] - 1];
    assert(prev_sample > 0);
    prev_sample--;

    // Calculate and return the new position, modulo the text length.
    pos = (prev_sample + delta) % textLength;

#ifdef PHI_BENCHMARK_FUNCTIONALITY
    // End timing with rdtscp
    unsigned int end_aux;
    uint64_t end = __builtin_ia32_rdtscp(&end_aux);

    // Check if the core has switched
    if (start_aux != end_aux) {
        cerr << "Core switch detected in phi after " << phi_call_count
             << " calls" << endl;

    } else {

        elapsed_phi += end - start;
        phi_call_count++;
    }
#endif
}

void BMove::phiInverse(length_t& pos) const {
#ifdef PHI_BENCHMARK_FUNCTIONALITY
    // Start timing with rdtscp
    unsigned int start_aux;
    uint64_t start = __builtin_ia32_rdtscp(&start_aux);
#endif
    // Find the rank of the predecessor of pos in a circular manner.
    length_t predRank = predLast.predecessorRankCircular(pos);
    // Select the predecessor position using its rank.
    length_t pred = predLast.select(predRank);

    // Calculate the distance from the predecessor to the current position.
    length_t delta = pred < pos ? pos - pred : pos + 1;

    // Ensure that phiInverse(SA[n-1]) is not called; the predecessor rank must
    // be valid.
    assert(lastToRun[predRank] < samplesFirst.size() - 1);

    // Get the next sample from the samplesFirst array.
    length_t prev_sample = samplesFirst[lastToRun[predRank] + 1];
    assert(prev_sample > 0);
    prev_sample--;

    // Calculate and return the new position, modulo the text length.
    pos = (prev_sample + delta) % textLength;

#ifdef PHI_BENCHMARK_FUNCTIONALITY
    // End timing with rdtscp
    unsigned int end_aux;
    uint64_t end = __builtin_ia32_rdtscp(&end_aux);

    // Check if the core has switched
    if (start_aux != end_aux) {
        cerr << "Core switch detected in phi after " << phi_call_count
             << " calls" << endl;

    } else {

        elapsed_phi += end - start;
        phi_call_count++;
    }
#endif
}
#endif

length_t BMove::computeToehold(const MoveRange& range, const length_t c) const {
#ifdef LF_BENCHMARK_FUNCTIONALITY
    return 0;
#endif

    length_t endRun = range.getEndRun();

    uint8_t endRunHead = move.getRunHead(endRun);

    if (endRunHead == c) {
        assert(samplesFirst[endRun] > 0);
        return samplesFirst[endRun] - 1;
    }

    length_t previousPos;
    length_t previousRun;

    move.walkToPreviousRun(range, previousPos, previousRun, c);

    assert(samplesLast[previousRun] > 0);
    return samplesLast[previousRun] - 1;
}

length_t BMove::computeToeholdRev(const MoveRange& range,
                                  const length_t c) const {
#ifdef LF_BENCHMARK_FUNCTIONALITY
    return 0;
#endif

    length_t endRun = range.getEndRun();

    length_t endRunHead = moveR.getRunHead(endRun);

    if (endRunHead == c) {
        assert(revSamplesFirst[endRun] > 0);
        return revSamplesFirst[endRun] - 1;
    }

    length_t previousPos;
    length_t previousRun;

    moveR.walkToPreviousRun(range, previousPos, previousRun, c);

    assert(revSamplesLast[previousRun] > 0);
    return revSamplesLast[previousRun] - 1;
}

// ----------------------------------------------------------------------------
//  APPROXIMATE PATTERN MATCHING
// ----------------------------------------------------------------------------

void BMove::updateRangeSARuns(SARangePair& ranges) const {
    move.computeRunIndices(ranges.getRangeSAMutable());
}

void BMove::findRangeWithExtraCharBackwardAuxiliary(
    length_t positionInAlphabet, SARange& parentBackwardRange,
    SARange& childBackwardRange) const {

    if (!parentBackwardRange.getRunIndicesValid()) {
        move.computeRunIndices(parentBackwardRange);
    }
    move.addChar(parentBackwardRange, childBackwardRange, positionInAlphabet);
}

bool BMove::findRangeWithExtraCharBackward(length_t posInAlpha,
                                           const SARangeBackwards& rangeOfP,
                                           SARangeBackwards& childRange) const {

    // first make the backward range by searching cP using B
    SARange range1, trivialRange = rangeOfP;
    findRangeWithExtraCharBackwardAuxiliary(posInAlpha, trivialRange, range1);

    if (range1.empty()) {
        childRange = SARangeBackwards(range1, 0, false, 0);
        return false;
    }

    // Check if range1 has the same width as the range of the parent
    if (trivialRange.width() == range1.width()) {
        childRange = SARangeBackwards(
            range1, rangeOfP.getToehold() - !rangeOfP.getToeholdRepresentsEnd(),
            rangeOfP.getToeholdRepresentsEnd(),
            rangeOfP.getOriginalDepth() + 1);
        return true;
    }

    // aP occurs for some a != c, the toehold must be updated
    length_t newToehold = computeToehold(trivialRange, posInAlpha);

    childRange = SARangeBackwards(range1, newToehold, false,
                                  rangeOfP.getOriginalDepth() + 1);
    return true;
}

bool BMove::findRangesWithExtraCharBackward(length_t positionInAlphabet,
                                            const SARangePair& rangesOfP,
                                            SARangePair& childRanges) const {

    // first make the backward range by searching cP using B
    SARange range1, trivialRange = rangesOfP.getRangeSA();
    findRangeWithExtraCharBackwardAuxiliary(positionInAlphabet, trivialRange,
                                           range1);

    // if the range is empty, return false
    if (range1.empty()) {
        childRanges = SARangePair(range1, range1, 0, false, 0);
        return false;
    }

    // Check if range1 has the same width as the range of the parent
    const SARange& otherRange = rangesOfP.getRangeSARev();
    if (trivialRange.width() == range1.width()) {
        // otherRange remains the same
        childRanges = SARangePair(range1, otherRange,
                                  rangesOfP.getToehold() -
                                      !rangesOfP.getToeholdRepresentsEnd(),
                                  rangesOfP.getToeholdRepresentsEnd(),
                                  rangesOfP.getOriginalDepth() + 1);
        return true;
    }

    // then make the less trivial range by counting the sizes of the ranges
    // of (dP) using B

    // first get the start of the range we are looking for of the parent
    length_t s = rangesOfP.getRangeSARev().getBegin();

    // get the start of the child within this range
    // find the number of occurrences of chars smaller than c in the parent
    // range
    length_t x = move.getCumulativeCounts(trivialRange, positionInAlphabet);

    // make the new range with width equal to that of the trivial range
    SARange range2 = SARange(s + x, s + x + range1.width(),
                             otherRange.getBeginRun(), otherRange.getEndRun());
    // set the run indices to be invalid
    range2.setRunIndicesValid(false);

    // aP occurs for some a != c, the toehold must be updated
    length_t newToehold = computeToehold(trivialRange, positionInAlphabet);

    // set the final SARangePair
    childRanges = SARangePair(range1, range2, newToehold, false,
                              rangesOfP.getOriginalDepth() + 1);
    return true;
}

bool BMove::findRangesWithExtraCharForward(length_t positionInAlphabet,
                                           const SARangePair& rangesOfP,
                                           SARangePair& childRanges) const {

    // first make the trivial range by searching (Pc)' using B'  if
    // searching forward we need to use B' so we need the reverse range
    SARange trivialRange = rangesOfP.getRangeSARev();

    // If the run indices of trivialRange are not valid (e.g., after a
    // directions switch), compute them
    if (!trivialRange.getRunIndicesValid()) {
        moveR.computeRunIndices(trivialRange);
    }

    // get the range of the child by adding ona character using move
    SARange range1;
    moveR.addChar(trivialRange, range1, positionInAlphabet);

    // if the range is empty, return false
    if (range1.empty()) {
        childRanges = SARangePair(range1, range1, 0, false, 0);
        return false;
    }

    // Check if range1 has the same width as the range of the parent
    const SARange& otherRange = rangesOfP.getRangeSA();
    if (trivialRange.width() == range1.width()) {
        // otherRange remains the same
        childRanges = SARangePair(otherRange, range1,
                                  rangesOfP.getToehold() +
                                      rangesOfP.getToeholdRepresentsEnd(),
                                  rangesOfP.getToeholdRepresentsEnd(),
                                  rangesOfP.getOriginalDepth() + 1);
        return true;
    }

    // then make the less trivial range by counting the size of the range of
    // (Pd)' using B' (forward)

    // first get the start of the range we are looking for of the parent
    length_t s = rangesOfP.getRangeSA().getBegin();

    // get the start of the child within this range
    // find the number of occurrences of chars smaller than c in the parent
    // range

    length_t x = moveR.getCumulativeCounts(trivialRange, positionInAlphabet);

    // make the new range with width equal to that of the trivial range
    SARange range2 = SARange(s + x, s + x + range1.width(),
                             otherRange.getBeginRun(), otherRange.getEndRun());
    // set the run indices to be invalid
    range2.setRunIndicesValid(false);

    // aP occurs for some a != c, the toehold must be updated
    length_t newToehold =
        textLength - 1 - computeToeholdRev(trivialRange, positionInAlphabet);

    // set the final SARangePair
    childRanges = SARangePair(range2, range1, newToehold, true,
                              rangesOfP.getOriginalDepth() + 1);
    return true;
}

bool BMove::findRangesWithExtraCharBackwardUniDirectional(
    length_t positionInAlphabet, const SARangePair& rangesOfP,
    SARangePair& rangesOfChild) const {

    // first make the trivial range by  searching cP using B
    SARange range1, trivialRange = rangesOfP.getRangeSA();
    findRangeWithExtraCharBackwardAuxiliary(positionInAlphabet, trivialRange,
                                           range1);

    // if the range is empty, return false
    if (range1.empty()) {
        rangesOfChild = SARangePair(range1, range1, 0, false, 0);
        return false;
    }

    // Check if range1 has the same width as the range of the parent
    if (trivialRange.width() == range1.width()) {
        // otherRange remains the same
        rangesOfChild = SARangePair(range1, SARange(),
                                    rangesOfP.getToehold() -
                                        !rangesOfP.getToeholdRepresentsEnd(),
                                    rangesOfP.getToeholdRepresentsEnd(),
                                    rangesOfP.getOriginalDepth() + 1);
        return true;
    }

    // aP occurs for some a != c, the toehold must be updated
    length_t newToehold = computeToehold(trivialRange, positionInAlphabet);

    // set the final SARangePair
    rangesOfChild = SARangePair(range1, SARange(), newToehold, false,
                                rangesOfP.getOriginalDepth() + 1);
    return true;
}

// ----------------------------------------------------------------------------
// ROUTINES FOR ACCESSING THE DATA STRUCTURE
// ----------------------------------------------------------------------------

SARangePair BMove::getRangeOfSingleChar(char c) const {
    int i = sigma.c2i(c);
    if (i < 0) {
        return SARangePair();
    }

    SARangePair pair = getCompleteRange();
    if ((unsigned int)i < sigma.size() - 1) {
        pair.getRangeSAMutable().setBegin(counts[i]);
        pair.getRangeSAMutable().setEnd(counts[i + 1]);
        pair.getRangeSARevMutable().setBegin(counts[i]);
        pair.getRangeSARevMutable().setEnd(counts[i + 1]);
        pair.getRangeSAMutable().setRunIndicesValid(false);
        pair.getRangeSARevMutable().setRunIndicesValid(false);
        return pair;
    }
    pair.getRangeSAMutable().setBegin(counts[i]);
    pair.getRangeSAMutable().setEnd(textLength);
    pair.getRangeSARevMutable().setBegin(counts[i]);
    pair.getRangeSARevMutable().setEnd(textLength);
    pair.getRangeSAMutable().setRunIndicesValid(false);
    pair.getRangeSARevMutable().setRunIndicesValid(false);
    return pair;
}

// ----------------------------------------------------------------------------
// POST-PROCESSING ROUTINES FOR APPROXIMATE MATCHING
// ----------------------------------------------------------------------------

void BMove::collectTextPositions(length_t firstPos, length_t originalDepth,
                                 std::vector<length_t>& positions) const {
#ifdef PHI_BENCHMARK_FUNCTIONALITY
    auto start = std::chrono::high_resolution_clock::now();
#endif

    length_t currentPos = firstPos;
#ifdef PHI_MOVE
    length_t textRun = predFirst.rank(firstPos);
#endif

    // Collect positions for the first while loop
    assert(currentPos < textLength);
    positions.push_back(currentPos);

    while (plcp[currentPos] >= originalDepth) {
#ifdef PHI_MOVE
        phiMove.phi(currentPos, textRun);
#else
        phi(currentPos);
#endif
        assert(currentPos < textLength);
        positions.push_back(currentPos);
    }

    // // Reverse the collected positions from the first while loop
    // std::reverse(positions.begin(), positions.end());

    currentPos = firstPos;
#ifdef PHI_MOVE
    textRun = predLast.rank(currentPos);
#endif

    // Collect positions for the second while loop
    while (currentPos != getInitialToehold()) {
#ifdef PHI_MOVE
        phiInvMove.phi(currentPos, textRun);
#else
        phiInverse(currentPos);
#endif
        if (plcp[currentPos] < originalDepth)
            break;
        assert(currentPos < textLength);
        positions.push_back(currentPos);
    }

#ifdef PHI_BENCHMARK_FUNCTIONALITY
    auto stop = std::chrono::high_resolution_clock::now();
    elapsed_locate += stop - start;
    locate_call_count++;
#endif
}

length_t BMove::locateCentered(const SARangePair& ranges, length_t& firstPos,
                               std::vector<length_t>& left,
                               std::vector<length_t>& right) const {

    length_t originalDepth = ranges.getOriginalDepth();

    assert(!ranges.getToeholdRepresentsEnd() ||
           ranges.getToehold() >= originalDepth - 1);

    firstPos = ranges.getToehold() -
               (ranges.getToeholdRepresentsEnd() ? (originalDepth - 1) : 0);

    length_t begin = ranges.getRangeSA().getBegin();

    // walk phi from the toehold occurrence (SA[c]) down to SA[begin], collecting
    // SA[c-1], SA[c-2], ..., SA[begin]; the number of steps equals c - begin.
    length_t currentPos = firstPos;
#ifdef PHI_MOVE
    length_t textRun = predFirst.rank(firstPos);
#endif
    while (plcp[currentPos] >= originalDepth) {
#ifdef PHI_MOVE
        phiMove.phi(currentPos, textRun);
#else
        phi(currentPos);
#endif
        assert(currentPos < textLength);
        left.push_back(currentPos);
    }

    length_t c = begin + left.size();

    // walk phi^{-1} from the toehold occurrence up to SA[end-1], collecting
    // SA[c+1], SA[c+2], ..., SA[end-1].
    currentPos = firstPos;
#ifdef PHI_MOVE
    textRun = predLast.rank(currentPos);
#endif
    while (currentPos != getInitialToehold()) {
#ifdef PHI_MOVE
        phiInvMove.phi(currentPos, textRun);
#else
        phiInverse(currentPos);
#endif
        if (plcp[currentPos] < originalDepth)
            break;
        assert(currentPos < textLength);
        right.push_back(currentPos);
    }

    return c;
}

void BMove::getTextPositionsFromSARange(
    const SARangePair& ranges, std::vector<length_t>& positions) const {

    assert(ranges.width() > 0);
    positions.reserve(ranges.width());

    assert(!ranges.getToeholdRepresentsEnd() ||
           ranges.getToehold() >= ranges.getOriginalDepth() - 1);

    length_t firstPos =
        ranges.getToehold() - (ranges.getToeholdRepresentsEnd()
                                   ? (ranges.getOriginalDepth() - 1)
                                   : 0);

    collectTextPositions(firstPos, ranges.getOriginalDepth(), positions);

    assert(ranges.width() == positions.size());
}

std::vector<length_t> BMove::getBeginPositions(const SARangeBackwards& rangeSA,
                                               length_t startDiff,
                                               length_t shift) const {

    assert(rangeSA.width() > 0);
    std::vector<length_t> positions;
    positions.reserve(rangeSA.width());

    assert(!rangeSA.getToeholdRepresentsEnd() ||
           rangeSA.getToehold() >= rangeSA.getOriginalDepth() - 1);

    length_t firstPos =
        rangeSA.getToehold() - (rangeSA.getToeholdRepresentsEnd()
                                    ? (rangeSA.getOriginalDepth() - 1)
                                    : 0);

    collectTextPositions(firstPos, rangeSA.getOriginalDepth(), positions);

    assert(rangeSA.width() == positions.size());
    return positions;
}