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
#ifndef INDEXINTERFACE_H
#define INDEXINTERFACE_H

#include "alphabet.h"
#include "bitparallelmatrix.h" // for BitParallelED64
#include "definitions.h"       // for length_t, PairStatus, Strand
#include "rindexhelpers.h"      // for SARangePair, Range, TextOcc
#include "reads.h"             // for ReadBundle
#include "substring.h"         // for Substring
#include "tkmer.h"             // for Kmer, KmerHash

#include <algorithm>                // for max
#include <parallel_hashmap/phmap.h> // phmap::parallel_flat_hash_map
#include <stddef.h>                 // for size_t
#include <stdint.h>                 // for uint16_t
#include <string>                   // for string, char_traits<>::pos_type
#include <utility>                  // for pair
#include <vector>                   // for vector

class MemoryMappedTextFile;
class Search;

class IndexInterface;
typedef bool (IndexInterface::*ExtraCharPtr)(length_t, const SARangePair&,
                                             SARangePair&) const;

/**
 * Class representing the bidirectional index and all operations on that
 * Index
 */
class IndexInterface {
  protected:
    // info about the text
    const std::string baseFile; // The base file of the reference text
    length_t textLength;        // the length of the text

    std::vector<length_t> counts; // the counts array of the reference genome

    Alphabet<ALPHABET> sigma; // the alphabet

    // direction variables
    thread_local static Direction dir; // the direction of the index
    thread_local static ExtraCharPtr
        extraChar; // pointer to extra char method (for direction)

    // stacks for search schemes
    thread_local static std::vector<std::vector<FMPosExt>>
        stacks; // stacks of nodes for the different searches
    thread_local static std::vector<BitParallelED64>
        matrices; // alignment matrices for the different parts

    thread_local static std::vector<BitParallelED128>
        matrices128; // big alignment matrices for the different parts

    // info about the current read/search
    thread_local static Strand
        strand; // The strand on which the current search is executed. If strand
                // is REVERSE_C_STRAND, then the read was reverse complemented
                // before starting the search and all found occurrences should
                // be labeled as reverse complemented.
    thread_local static PairStatus
        pairStatus; // The status of the read that is currently being searched
                    // in the pair.

    // sparse hash info
    const size_t wordSize = 10; // the size of the mers to be stored in a table
    phmap::parallel_flat_hash_map<Kmer, SARangePair, KmerHash>
        table; // hashtable that contains all wordSize-mers

    // assignment/SAM info
    std::vector<length_t> startPos;    // the start positions of the sequences
    std::vector<std::string> seqNames; // the names of the sequences

    // ----------------------------------------------------------------------------
    // PREPROCESSING ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Private helper function that reads in all the necessary files
     * @param baseFile the baseFile of the files that will be read in
     * @param verbose if true the steps will we written to cout
     */
    virtual void fromFiles(const std::string& baseFile, bool verbose) = 0;

    /**
     * Read the counts array, the build_tag info and the
     * compilation-bits info of the index.
     * @param baseFile The base name of the files that contain the info about
     * the index.
     * @param verbose If true, the steps will be written to cout.
     */
    void readTagCompAndCounts(const std::string& baseFile, bool verbose);

    /**
     * Read the sequence names and start positions of the sequences in the
     * index.
     * @param baseFile The base name of the files that contain the info about
     * the index.
     * @param verbose If true, the steps will be written to cout.
     */
    void readSequenceNamesAndPositions(const std::string& baseFile,
                                       bool verbose);

    /**
     * Read a binary file and stores content in array
     * @param filename File name
     * @param array Suffix array (contents will be overwritten)
     * @returns True if successful, false otherwise
     */
    static bool readArray(const std::string& filename,
                          std::vector<length_t>& array);

    /**
     * Populate the hash table
     * @param verbose if steps are written to cout
     */
    void populateTable(bool verbose);

    // ----------------------------------------------------------------------------
    // ROUTINES FOR ACCESSING DATA STRUCTURE
    // ----------------------------------------------------------------------------

    /**
     * Get the start positions of the sequences
     * @return The start positions vector
     */
    const std::vector<length_t>& getStartPositions() const {
        return startPos;
    }

    // ----------------------------------------------------------------------------
    //  APPROXIMATE PATTERN MATCHING
    // ----------------------------------------------------------------------------
    /**
     * Matches a search recursively with a depth first approach (each branch
     * of the tree is fully examined until the backtracking condition is
     * met) using the edit distance metric. This function uses all
     * optimizations for eliminating redundancy in the edit distance metric
     * @param search, the search to follow
     * @param startMatch, the approximate match found for all previous
     * partitions of the search
     * @param occ, a data structure with matches of the complete search, if
     * such a match is found it will be added to this data structure
     * @param parts, the parts of the pattern
     * @param counters the performance counters
     * @param idx, the index of the partition to match, defaults to 1 as an
     * exact search for the zeroth partition is assumed
     * @param descPrevDir, the descendants of the previous direction,
     * defaults to empty vector
     * @param initPrevDir, the initialization eds of the previous direction,
     * defaults to empty vector
     * @param descNotPrevDir, the descendants of the other direction,
     * defaults to empty vector
     * @param initNotPrevDir, the initialization eds of the other direction,
     * defaults to empty vector
     */
    void recApproxMatchEdit(
        const Search& search, const FMOcc& startMatch, Occurrences& occ,
        const std::vector<Substring>& parts, Counters& counters,
        const int& idx = 1,
        const std::vector<FMPosExt>& descPrevDir = std::vector<FMPosExt>(),
        const std::vector<uint16_t>& initPrevDir = std::vector<uint16_t>(),
        const std::vector<FMPosExt>& descNotPrevDir = std::vector<FMPosExt>(),
        const std::vector<uint16_t>& initNotPrevDir = std::vector<uint16_t>());

    /**
     * @brief Update the runs corresponding to the SA range
     *
     * @param ranges SARangePair for which the SA range runs must be updated
     */
    virtual void updateRangeSARuns(SARangePair& ranges) const = 0;

    /**
     * Finds the ranges of cP using the principle explained in the paper of
     * Lam.
     * @param positionInAlphabet the position in alphabet of the character
     * that is added in the front.
     * @param rangesOfP the ranges of pattern P.
     * @param childRanges the ranges corresponding to string cP, this
     * will be set during execution
     * @returns true if the range is not empty, false otherwise
     */
    virtual bool
    findRangesWithExtraCharBackward(length_t positionInAlphabet,
                                    const SARangePair& rangesOfP,
                                    SARangePair& childRanges) const = 0;

    /**
     * Finds the ranges of string Pc using the principle explained in the paper
     * of Lam
     * @param positionInAlphabet the position in alphabet of the character c
     * that is added in the back
     * @param rangesOfP the ranges of pattern P
     * @param childRanges the ranges corresponding to string Pc, this will be
     * set during execution
     * @returns true if the range is not empty, false otherwise
     */
    virtual bool
    findRangesWithExtraCharForward(length_t positionInAlphabet,
                                   const SARangePair& rangesOfP,
                                   SARangePair& childRanges) const = 0;

    /**
     * Finds the range of Pc using unidirectional backward matching. With
     * toeholds.
     * @param positionInAlphabet the position in alphabet of the character c
     * that is added in the back
     * @param rangeOfP the range over the suffix array of the text
     * corresponding to pattern P
     * @param childRange the range over the suffix array of text
     * corresponding to pattern Pc this will be set during execution
     * @returns true if the range is not empty, false otherwise
     */
    virtual bool
    findRangeWithExtraCharBackward(length_t positionInAlphabet,
                                   const SARangeBackwards& rangeOfP,
                                   SARangeBackwards& childRange

    ) const = 0;

    /**
     * Finds the range of cP and stores a 'dummy' range over the suffix array of
     * the reversed text
     * @param positionInAlphabet the position in alphabet of the character
     * that is added in the front.
     * @param rangesOfP the ranges of pattern P, only the backwards range is
     * used
     * @param childRanges the ranges corresponding to string cP, this
     * will be set during execution. The getRangeSARev() function will return a
     * dummy!
     * @returns true if the range is not empty, false otherwise
     */
    virtual bool findRangesWithExtraCharBackwardUniDirectional(
        length_t positionInAlphabet, const SARangePair& rangesOfP,
        SARangePair& rangesOfChild) const = 0;

    /**
     * Helper function for the approximate matching. This function fills in
     * the matrix for the current node at the current row and goes deeper
     * for the next part if that is necessary. The function returns true if
     * the search should backtrack.
     * @param bpED The pointer to the bit-parallel alignment matrix of the
     * current part.
     * @param cluster The cluster corresponding to the final column of the
     * matrix.
     * @param currentNode The node for which the matrix is filled in.
     * @param s The current search.
     * @param idx The index of the current part in the search.
     * @param parts The parts of the pattern.
     * @param occ Data structure with the in-index and in-text occurrences,
     * if an occurrence, either in-index or in-text, is found it will be
     * added to this data structure.
     * @param counters The performance counters.
     * @param initOther Edit distance values of known descendants in the
     * other direction.
     * @param descOther The known descendants in the other direction.
     * @param remainingDesc The remaining descendants on the current branch,
     * that are already created but aren't checked yet and might need to be
     * checked for the next part, defaults to an empty vector.
     * @return false if the search can continue along this branch for the
     * current part, true if the search should backtrack.
     */
    bool branchAndBound(IBitParallelED* bpED, Cluster& cluster,
                        const FMPosExt& currentNode, const Search& s,
                        const length_t& idx,
                        const std::vector<Substring>& parts, Occurrences& occ,
                        Counters& counters,
                        const std::vector<uint16_t>& initOther,
                        const std::vector<FMPosExt>& descOther,
                        const std::vector<FMPosExt>& remainingDesc = {});

    /**
     * Goes deeper in a search if a valid partial approximate match is found in
     * the cluster. This function is called by the recursive function's helper
     * branchAndBound function..
     * @param cluster The cluster corresponding to the final column of the
     * matrix.
     * @param nIdx The index of next part to in search s.
     * @param s The search.
     * @param parts The parts of the pattern.
     * @param occ Data structure with the in-index and in-text occurrences, if
     * an occurrence, either in-index or in-text, is found it will be added to
     * this data structure.
     * @param counters The performance counters.
     * @param descOtherD The known descendants in the other direction.
     * @param initOtherD Edit distance values of known descendants in the other
     * direction.
     * @param remDescThe remaining descendants on the current branch,
     * that are already created but aren't checked yet and might need to be
     * checked for the next part.
     */
    void goDeeper(Cluster& cluster, const length_t& nIdx, const Search& s,
                  const std::vector<Substring>& parts, Occurrences& occ,
                  Counters& counters, const std::vector<FMPosExt>& descOtherD,
                  const std::vector<uint16_t>& initOtherD,
                  const std::vector<FMPosExt>& remDesc);

    // ----------------------------------------------------------------------------
    // LOCATION ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Find the begin positions of the range that needs to be verified in-text
     * @param rangeSA The range over the SA of the partial in-index match that
     * needs to be verified in text.
     * @param startDiff The highest possible difference between the start of the
     * partial match and that of a valid full match in the text.
     * @param shift The shift of the node. Necessary in case the search in the
     * backwards direction has ended. (default = 0).
     * @returns The lowest possible begin position for each of the partial
     * matches in the range.
     */
    virtual std::vector<length_t>
    getBeginPositions(const SARangeBackwards& rangeSA, length_t startDiff,
                      length_t shift = 0) const = 0;

    // ----------------------------------------------------------------------------
    // EXTEND ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Creates the children nodes of the node with the given ranges and pushes
     * them to the stack.
     * @param ranges The ranges of the parent node to get the children of.
     * @param stack The stack to push the children on.
     * @param counters The performance counters.
     * @param row The row in the current matrix of the parentNode (defaults to
     * 0).
     */
    void extendFMPos(const SARangePair& ranges, std::vector<FMPosExt>& stack,
                     Counters& counters, length_t row = 0) const;

    /**
     * Creates the children nodes of the node with the given position and pushes
     * them to the stack.
     * @param pos The position of the parent node to get the children of.
     * @param stack The stack to push the children on.
     * @param counters The performance counters.
     */
    void extendFMPos(const FMPosExt& pos, std::vector<FMPosExt>& stack,
                     Counters& counters) const;

    // ----------------------------------------------------------------------------
    // SUBROUTINES FOR PATTERN MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Sets the index in the correct mode. All found occurrences will now be
     * be labeled as found along this strand and as an occurrence of this
     * read in the pair.
     * @param strand On which strand the current read lies
     * @param pairStatus The status of the current read in its pair. [default =
     * FIRST_IN_PAIR]
     */
    void setIndexInModeSubRoutine(Strand reverseComplement,
                                  PairStatus firstRead = FIRST_IN_PAIR) {
        strand = reverseComplement;
        pairStatus = firstRead;
    }

  public:
    // ----------------------------------------------------------------------------
    // INITIALIZATION ROUTINES
    // ----------------------------------------------------------------------------

    /**
     * Constructor
     * @param baseFile The base name of the files that contain the info about
     * the index.
     * @param verbose If true, the steps will be written to cout. [default =
     * true]
     * @param wordSize The size of the mers to be stored in the hashtable. Used
     * for quick look-ups of exact seeds. [default = 10]
     */
    IndexInterface(const std::string& baseFile, bool verbose = true,
                   length_t wordSize = 10)
        : baseFile(baseFile), wordSize(wordSize) {
    }

    /**
     * @brief Destructor
     *
     */
    virtual ~IndexInterface() = default;

    /**
     * Get the complete range of this index corresponding to an exact match of
     * the empty string.
     * @returns A SARangePair with both ranges the complete range of the
     * index.
     */
    virtual SARangePair getCompleteRange() const = 0;

    // ----------------------------------------------------------------------------
    // ROUTINES FOR ACCESSING THE DATA STRUCTURE
    // ----------------------------------------------------------------------------

 
    /**
     * Get the length of the text
     */
    const length_t getTextSize() const {
        return textLength;
    }

    /**
     * @returns the wordsize of the mers stored in the table
     */
    length_t getWordSize() const {
        return wordSize;
    }

    // ----------------------------------------------------------------------------
    // ADDITIONS FOR Move-r APPROXIMATE-PATTERN-MATCHING ADAPTER
    // ----------------------------------------------------------------------------

    /**
     * @returns the size of the alphabet (including the sentinel)
     */
    length_t alphabetSize() const {
        return sigma.size();
    }

    /**
     * Maps a character to its index in the alphabet.
     * @param c the character
     * @returns 0 for the sentinel, 1..alphabetSize()-1 for the other
     *          characters, or -1 if the character does not occur in the alphabet
     */
    int charToIdx(char c) const {
        return sigma.c2i(c);
    }

    /**
     * Maps an alphabet index back to its character.
     * @param i the alphabet index (in [0, alphabetSize()))
     * @returns the corresponding character
     */
    char idxToChar(int i) const {
        return sigma.i2c(i);
    }

    /**
     * Extends a match by one character in the given direction. This is a thin
     * public wrapper around the (protected) bidirectional extension routines so
     * that an external bidirectional searcher can drive the index.
     * @param positionInAlphabet index in the alphabet of the character to add
     * @param d FORWARD to append (Pc), BACKWARD to prepend (cP)
     * @param parent the ranges of the current match P
     * @param child output: the ranges of the extended match
     * @returns true if the extended match occurs in the text, false otherwise
     */
    bool extendRange(length_t positionInAlphabet, Direction d,
                     const SARangePair& parent, SARangePair& child) const {
        return (d == FORWARD)
                   ? findRangesWithExtraCharForward(positionInAlphabet, parent,
                                                    child)
                   : findRangesWithExtraCharBackward(positionInAlphabet, parent,
                                                     child);
    }

    /**
     * Find the range of an exact match of single character in this index.
     * Warning: this function assumes that the character is in the alphabet.
     * @returns the ranges of a single character in this index
     */
    virtual SARangePair getRangeOfSingleChar(char c) const = 0;

    /**
     * Looks up the SARangePair of the exact match of c p in the hashtable.
     * Warning: This function assumes p is of size wordSize.
     * @param p The substring to find the ranges of.
     * @returns the ranges corresponding to substring p, if no pair can be
     * found it returns empty ranges.
     */
    SARangePair lookUpInKmerTable(const Substring& p) const {
        auto it = table.find(Kmer(p.tostring()));
        return (it == table.end() || p.containsN()) ? SARangePair()
                                                    : it->second;
    }

    const std::vector<std::string>& getSeqNames() const {
        return seqNames;
    }

    /**
     * Finds the sequence name to which the given text-occurrence belongs. Also
     * changes the range of the text-occurrence to be relative to this sequence.
     * @param t The text-occurrence to find the sequence name of.
     * @param seqID the unique ID of the sequence to which the text-occurrence
     * belongs (if a sequence could be found)
     * @param counters The performance counters.
     * @param largestStratum the largest stratum of distance that is allowed.
     * @param metric the used distance metric. In case of hamming distance no
     * trimming will be performed.
     * @param pattern The pattern that was matched, needed for trimming
     * @returns SeqNameFound::FOUND if the name was found without
     * trimming, SeqNameFound::FOUND_WITH_TRIMMING if the name was found with
     * trimming and SeqNameFound::NOT_FOUND if no name could be found (i.e. t
     * spans two references and trimming cannot fix it).
     */
    SeqNameFound findSeqName(TextOcc& t, length_t& seqID, Counters& counters,
                             length_t largestStratum,
                             const DistanceMetric& metric,
                             const std::string& pattern) const;
    // ----------------------------------------------------------------------------
    // ROUTINES FOR EXACT MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Calculates the exact matches to the string in the index. Creates
     * corresponding text-occurrences and adds them to the output vector.
     * Make sure to set the index in the correct mode before calling this
     * function to correctly label the found occurrences.
     * @param s The string to match in the reference genome.
     * @param counters The performance counters.
     * @param tOcc The vector to which the found text occurrences will be added.
     */
    void exactMatchesOutput(const std::string& s, Counters& counters,

                            std::vector<TextOcc>& tOcc) const;

    /**
     * This function matches a string exactly and returns the ranges in the
     * suffix array of the text and of the reverse text. If the string cannot be
     * matched empty ranges will be returned.
     * Warning: this function assumes the direction of the index
     * corresponds to the direction of pattern!
     * @param pattern The string to match.
     * @param counters The performance counters.
     * @returns The pair of ranges corresponding to an exact match of this
     * pattern.
     */
    SARangePair matchStringBidirectionally(const Substring& pattern,
                                           Counters& counters) const {
        return matchStringBidirectionally(pattern, getCompleteRange(),
                                          counters);
    }

    /**
     * This function matches a string exactly starting form a previous start
     * range.
     * Warning: this function assumes the direction of the index corresponds to
     * the direction of pattern!
     * @param pattern The string to match.
     * @param startRange The range to start the match from.
     * @param counters The performance counters.
     * @returns The pair of ranges corresponding to an exact match of this
     * pattern.
     */
    SARangePair matchStringBidirectionally(const Substring& pattern,
                                           SARangePair startRange,
                                           Counters& counters) const;

    /**
     * Adds one character and updates the range. If the character cannot be
     * added the range will be set to an empty range.
     * Warning: this function assumes the direction of the index is set
     * correctly!
     * @param c The character to be added (in the current direction of the
     * index).
     * @param range The range to update.
     * @param counters The performance counters.
     */
    bool addChar(const char& c, SARangePair& range, Counters& counters) const;

    // ----------------------------------------------------------------------------
    // PREPROCESSING ROUTINES FOR APPROXIMATE MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Resizes the stacks vector to the required number of stacks and reserves
     * space on each stack, so that each stack can match the entire pattern
     * @param number The number of stacks required. This is equal to the number
     * of parts.
     * @param size The size of the pattern.
     */
    void reserveStacks(const length_t number, const length_t size) {
        stacks.resize(number);
        length_t stackSize = size * sigma.size();
        for (auto& stack : stacks) {
            stack.reserve(stackSize);
        }
    }
    /**
     * Reset the alignment matrices for the different partitions.
     * @param number The number of partitions.
     */
    void resetMatrices(const length_t number) {
        matrices.resize(2 * number);
        matrices128.resize(2 * number);

        for (auto& matrix : matrices) {
            matrix.reset();
        }
        for (auto& matrix : matrices128) {
            matrix.reset();
        }
    }

    // ----------------------------------------------------------------------------
    // ROUTINES FOR APPROXIMATE MATCHING
    // ----------------------------------------------------------------------------

    /**
     * Matches the pattern approximately via a very naive backtracking approach.
     * All found matches are at most a certain edit distance away from the
     * pattern and will be added to the return vector.
     * @param pattern the pattern to match
     * @param maxED the maximum edit distance
     * @param counters the performance counters
     * @param occ The vector to which the matches will be added.
     */
    void approxMatchesNaive(const std::string& pattern, length_t maxED,
                            Counters& counters, std::vector<TextOcc>& occ);

    /**
     * Sets the search direction of the fm-index
     * @param d the direction to search in, either FORWARD or BACKWARD
     * @param isUniDirectionalBackward if the search is unidirectional backward
     */
    void setDirection(Direction d, bool isUniDirectionalBackward) {
        dir = d;
        extraChar =
            (isUniDirectionalBackward)
                ? &IndexInterface::findRangesWithExtraCharBackwardUniDirectional
                : ((d == FORWARD)
                       ? &IndexInterface::findRangesWithExtraCharForward
                       : &IndexInterface::findRangesWithExtraCharBackward);
    }

    /**
     * Sets the index in the correct mode. All found occurrences will now be
     * be labeled as found along this strand and as an occurrence of this
     * read in the pair.
     * @param strand On which strand the current read lies
     * @param pairStatus The status of the current read in its pair. [default =
     * FIRST_IN_PAIR]
     */
    virtual void setIndexInMode(Strand reverseComplement,
                                PairStatus firstRead = FIRST_IN_PAIR) = 0;

    /**
     * Matches a search recursively with a depth first approach (each branch
     * of the tree is fully examined until the backtracking condition is
     * met) using the hamming distance metric.
     * @param s The search to follow.
     * @param startMatch The partial approximate match found for the previous
     * parts of the search
     * @param occ A data structure with approximate occurrences of the complete
     * pattern. Both in-index and in-text occurrences are stored. If a new
     * approximate occurrence is found, either in-index or in-text it will be
     * added to this data structure.
     * @param parts The parts of the pattern.
     * @param counters The performance counters.
     * @param idx The index of current part in the search, defaults to 1 as an
     * exact search for the zeroth part is assumed.
     */
    void recApproxMatchHamming(const Search& s, const FMOcc& startMatch,
                               Occurrences& occ,
                               const std::vector<Substring>& parts,
                               Counters& counters, const int& idx = 1);

    /**
     * Entry to the recursive approximate matching procedure for the edit
     * distance metric.
     * @param search The search to follow.
     * @param startMatch The match containing the SA ranges corresponding to
     * the match of the first part of the search
     * @param occ A data structure with approximate occurrences of the complete
     * pattern. Both in-index and in-text occurrences are stored. If a new
     * approximate occurrence is found, either in-index or in-text it will be
     * added to this data structure.
     * @param parts The parts of the pattern.
     * @param counters The performance counters.
     * @param idx The index of current part in the search, defaults to 1 as an
     * exact search for the zeroth part is assumed.
     */
    void recApproxMatchEditEntry(const Search& search, const FMOcc& startMatch,
                                 Occurrences& occ,
                                 const std::vector<Substring>& parts,
                                 Counters& counters, const int& idx = 1);

    // ----------------------------------------------------------------------------
    // POST-PROCESSING ROUTINES FOR APPROXIMATE MATCHING
    // ----------------------------------------------------------------------------

    /**
     * @brief Get the text positions corresponding to a suffix array range
     *
     * @param range The range in the suffix array
     * @param positions The vector to which the text positions will be added
     */
    virtual void
    getTextPositionsFromSARange(const SARangePair& ranges,
                                std::vector<length_t>& positions) const = 0;

    /**
     * Find unique text occurrences for the hamming distance. The in-index
     * occurrences are converted to in-text occurrences with CIGAR string.
     * Afterwards redundant occurrences are removed. Eventually all
     * non-redundant text-occurrences are returned.
     * @param occ The occurrences to find the unique text occurrences for.
     * @param counters The performance counters.
     * @returns the non-redundant text occurrences.
     */
    std::vector<TextOcc> getTextOccHamming(Occurrences& occ,
                                           Counters& counters) const;

    /**
     * Find unique text occurrences for the edit distance. The in-index
     * occurrences are converted to in-text occurrences.
     * Afterwards redundant occurrences are removed. For all non-redundant
     * in-text occurrences that do not have a CIGAR string a CIGAR string is
     * calculated. Eventually all non-redundant text-occurrences are returned
     * sorted.
     * @param occ The occurrences to find the unique text occurrences for.
     * @param maxED The maximal edit distance that is allowed for the search.
     * @param counters The performance counters.
     * @param bundle The read bundle with info about the read/reverse complement
     * @returns the non-redundant text occurrences.
     */
    std::vector<TextOcc>
    getUniqueTextOccurrences(Occurrences& occ, const length_t& maxED,
                             Counters& counters,
                             const ReadBundle& bundle) const;
};

#endif
