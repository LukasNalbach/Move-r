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
 * @file alignment_output_cli.hpp
 * @brief writes the approximate occurrences of reads as SAM and PAF records (move-rb locate -sam, -paf)
 */

#ifndef MOVE_R_MISC_ALIGNMENT_OUTPUT_CLI_HPP
#define MOVE_R_MISC_ALIGNMENT_OUTPUT_CLI_HPP

#include <algorithm>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include <data_structures/fasta_sequence_data.hpp>
#include <misc/apm.hpp>
#include <misc/cigar.hpp>
#include <misc/parallel_patterns.hpp>

// MAPQ of the primary SAM record of a read whose fewest errors only one occurrence has
inline constexpr uint16_t unique_mapq = 60;

// SAM flags (FLAG column): the record is a secondary alignment, the read is aligned as its reverse complement, the
// read is unmapped
inline constexpr uint16_t sam_flag_secondary = 0x100;
inline constexpr uint16_t sam_flag_reverse = 0x10;
inline constexpr uint16_t sam_flag_unmapped = 0x4;

// MAPQ of a PAF record (255: not available)
inline constexpr uint16_t paf_mapq = 255;

/**
 * @brief the occurrences of a read on one strand and their CIGARs
 * @tparam cm whether the occurrences carry a CIGAR
 */
template <cigar_mode_t cm>
struct strand_occurrences_t {
    std::vector<aprx_occ_t<cm>> occurrences; // the occurrences, sorted by position
    // the alignments of the occurrences: [0..c-1] CIGAR of the c distinct alignments; occurrences[i].cig_idx is an
    // index into it (only with cm == CIGAR)
    std::vector<cigar_t<uint8_t>> cigars;

    // removes the occurrences and the CIGARs
    void clear()
    {
        occurrences.clear();
        cigars.clear();
    }
};

/**
 * @brief writes the SAM record of an occurrence of a read
 * @tparam cm whether the occurrence carries a CIGAR (else the CIGAR and the tag MD are not written)
 * @param sam the output
 * @param sequences the sequences of the text
 * @param last_sequence index of the sequence of the previous occurrence (sequence_index searches from it)
 * @param qname the name of the read
 * @param occ the occurrence (it lies within one sequence: the search does not cross a separator)
 * @param cigars the CIGARs occ.cig_idx refers to
 * @param reverse whether the occurrence is one of the reverse complement of the read (sam_flag_reverse)
 * @param primary whether the record is the primary record of the read (else sam_flag_secondary)
 * @param mapq the mapping quality
 * @param seq the read as it was searched (forward strand or reverse complement); SEQ, only written with write_seq
 * @param qual the qualities of the read on the forward strand (empty: none); QUAL, reversed for the reverse strand
 * @param write_seq whether SEQ and QUAL are written (only for the primary record)
 */
template <cigar_mode_t cm>
void write_sam_record(string_appender& sam, const fasta_sequence_data& sequences, uint64_t& last_sequence,
    const std::string& qname, const aprx_occ_t<cm>& occ, const std::vector<cigar_t<uint8_t>>& cigars, bool reverse,
    bool primary, uint16_t mapq, const std::string& seq, const std::string& qual, bool write_seq)
{
    const uint64_t seq_idx = sequences.sequence_index(occ.pos(), last_sequence);
    const uint64_t pos_in_sequence = occ.pos() - sequences.sequence_start(seq_idx);
    const uint16_t flag = (primary ? 0 : sam_flag_secondary) | (reverse ? sam_flag_reverse : 0);
    sam << qname << '\t' << flag << '\t' << sequences.sequence_name(seq_idx) << '\t' << (pos_in_sequence + 1) << '\t'
        << mapq << '\t';
    if constexpr (cm == CIGAR) append_cigar(sam, cigars[occ.cig_idx]);
    else sam << '*';
    sam << "\t*\t0\t0\t";

    if (primary && write_seq) {
        sam << seq << '\t';
        if (qual.empty()) sam << '*';
        else if (reverse) for (auto it = qual.rbegin(); it != qual.rend(); ++it) sam << *it;
        else sam << qual;
    } else {
        sam << "*\t*";
    }

    // NM: the number of errors of the alignment, MD: the symbols of the text at the mismatches and deletions
    sam << "\tNM:i:" << occ.err();

    if constexpr (cm == CIGAR) {
        sam << "\tMD:Z:";
        append_md(sam, cigars[occ.cig_idx]);
    }

    sam << '\n';
}

/**
 * @brief writes the SAM records of a read: the first occurrence with the fewest errors (forward strand first) as
 *        primary record, the others as secondary records, or one unmapped record if it has no occurrence; MAPQ =
 *        unique_mapq if exactly one occurrence has the fewest errors, else 0
 * @tparam cm whether the occurrences carry a CIGAR
 * @param sam the output
 * @param sequences the sequences of the text
 * @param last_sequence index of the sequence of the previous occurrence (sequence_index searches from it)
 * @param qname the name of the read
 * @param seq the read
 * @param rc_seq the reverse complement of the read
 * @param qual the qualities of the read (empty: none)
 * @param fwd the occurrences of the read
 * @param rc the occurrences of its reverse complement
 * @param write_seq whether SEQ and QUAL are written in the primary record (-seq)
 * @param best_only whether only the primary record is written (-best)
 */
template <cigar_mode_t cm>
void write_sam_records(string_appender& sam, const fasta_sequence_data& sequences, uint64_t& last_sequence,
    const std::string& qname, const std::string& seq, const std::string& rc_seq, const std::string& qual,
    const strand_occurrences_t<cm>& fwd, const strand_occurrences_t<cm>& rc, bool write_seq, bool best_only)
{
    if (fwd.occurrences.empty() && rc.occurrences.empty()) {
        sam << qname << '\t' << sam_flag_unmapped << "\t*\t0\t0\t*\t*\t0\t0\t"
            << (write_seq ? seq : std::string("*")) << '\t'
            << (write_seq && !qual.empty() ? qual : std::string("*")) << '\n';
        return;
    }

    uint64_t min_errors = std::numeric_limits<uint64_t>::max();
    for (const auto& occ : fwd.occurrences) min_errors = std::min<uint64_t>(min_errors, occ.err());
    for (const auto& occ : rc.occurrences) min_errors = std::min<uint64_t>(min_errors, occ.err());
    uint64_t num_min_errors = 0;
    for (const auto& occ : fwd.occurrences) num_min_errors += occ.err() == min_errors;
    for (const auto& occ : rc.occurrences) num_min_errors += occ.err() == min_errors;
    const uint16_t primary_mapq = num_min_errors == 1 ? unique_mapq : 0;
    bool primary_written = false;

    // writes the records of the occurrences of one strand
    auto write_strand = [&](const strand_occurrences_t<cm>& strand, bool reverse, const std::string& strand_seq) {
        for (const auto& occ : strand.occurrences) {
            if (best_only && primary_written) return;
            const bool primary = !primary_written && occ.err() == min_errors;
            if (best_only && !primary) continue;
            primary_written |= primary;
            write_sam_record<cm>(sam, sequences, last_sequence, qname, occ, strand.cigars, reverse, primary,
                primary ? primary_mapq : 0, strand_seq, qual, write_seq);
        }
    };

    write_strand(fwd, false, seq);
    write_strand(rc, true, rc_seq);
}

/**
 * @brief writes one PAF record per occurrence of a read on one strand; columns: qname, qlen, 0, qlen, strand, tname,
 *        tlen, tstart, tend, matches, block length, paf_mapq, NM and (with a CIGAR) cg
 * @tparam cm whether the occurrences carry a CIGAR (else the matches are the length minus the errors)
 * @param paf the output
 * @param sequences the sequences of the text
 * @param last_sequence index of the sequence of the previous occurrence (sequence_index searches from it)
 * @param qname the name of the read
 * @param read_length the length of the read
 * @param strand the occurrences
 * @param reverse whether the occurrences are those of the reverse complement of the read
 */
template <cigar_mode_t cm>
void write_paf_records(string_appender& paf, const fasta_sequence_data& sequences, uint64_t& last_sequence,
    const std::string& qname, uint64_t read_length, const strand_occurrences_t<cm>& strand, bool reverse)
{
    for (const auto& occ : strand.occurrences) {
        const uint64_t seq_idx = sequences.sequence_index(occ.pos(), last_sequence);
        const uint64_t pos_in_sequence = occ.pos() - sequences.sequence_start(seq_idx);
        uint64_t matches = 0; // the number of matching symbols
        uint64_t block_length = 0; // the number of columns of the alignment

        if constexpr (cm == CIGAR) {
            for (const cigar_run_t<uint8_t>& run : strand.cigars[occ.cig_idx]) {
                block_length += run.len;
                if (run.op == cigar_op_t::MATCH) matches += run.len;
            }
        } else {
            block_length = occ.len();
            matches = occ.len() > occ.err() ? uint64_t(occ.len() - occ.err()) : 0;
        }

        paf << qname << '\t' << read_length << "\t0\t" << read_length << '\t' << (reverse ? '-' : '+') << '\t'
            << sequences.sequence_name(seq_idx) << '\t' << sequences.sequence_length(seq_idx) << '\t'
            << pos_in_sequence << '\t' << (pos_in_sequence + occ.len()) << '\t' << matches << '\t' << block_length
            << '\t' << paf_mapq << "\tNM:i:" << occ.err();

        if constexpr (cm == CIGAR) {
            paf << "\tcg:Z:";
            append_cigar(paf, strand.cigars[occ.cig_idx]);
        }

        paf << '\n';
    }
}

#endif // MOVE_R_MISC_ALIGNMENT_OUTPUT_CLI_HPP
