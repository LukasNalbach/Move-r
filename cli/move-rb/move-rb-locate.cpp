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

#include <filesystem>
#include <iostream>
#include <ips2ra.hpp>
#include <move_rb/move_rb.hpp>
#include <misc/fasta.hpp>
#include <misc/progress.hpp>

static constexpr int32_t min_args = 6;
int32_t arg_idx = 1;
uint64_t k = -1;
std::string scheme_str;
distance_metric_t dist_metr = NO_METRIC;
search_scheme_t search_scheme;
bool output_occurrences = false;
bool check_correctness = false;
bool sam_output = false;       // whether to write occurrences in SAM format (requires a FASTA-built index)
bool sam_include_seq = false;  // whether to include the read sequence (SEQ) in SAM records
bool sam_rc = true;            // whether to also align each read's reverse complement (SAM mode only)
bool fastq_input = false;      // whether the patterns file is FASTQ (variable-length reads with per-base qualities)
bool sam_best_only = false;    // whether to report only the best (minimum-error) alignment per read
bool sam_no_cigar = false;     // whether to omit CIGARs from SAM records (faster: skips the CIGAR computation)
bool paf_output = false;       // whether to write occurrences in PAF format (requires a FASTA-built index)
std::ofstream mf;
std::string path_index_file;
std::string path_patterns_file;
std::string path_output_file;
std::string path_sam_file;
std::string path_paf_file;
std::ifstream index_file;
std::ifstream patterns_file;
std::ifstream input_file;
std::ofstream output_file;
std::ofstream sam_file;
std::ofstream paf_file;
std::string name_text_file;
std::string path_input_file;
std::string command_line; // the full command line, for the SAM @PG header

/**
 * @brief prints the usage information and exits
 * @param msg an optional error message printed before the usage information
 */
void help(std::string msg)
{
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "move-rb-locate: locate all approximate occurrences of the input patterns." << std::endl << std::endl;
    std::cout << "usage: move-rb-locate [...] -d <metric> -s <scheme> [-k <mismatches>] <index_file> <patterns_file>" << std::endl;
    std::cout << "   -m <m_file> <text_name>    m_file is the file to write measurement data to," << std::endl;
    std::cout << "                              text_name should be the name of the original file" << std::endl;
    std::cout << "   -c <input_file>            checks correctness of each pattern occurrence against <input_file>" << std::endl;
    std::cout << "                              (input_file must be the file the index was built for)" << std::endl;
    std::cout << "   <metric>                   distance metric to use (hamming or edit)" << std::endl;
    std::cout << "   <scheme>                   search scheme to use (pigeon_hole, suffix_filter, min_u, 01 or path to a file)" << std::endl;
    std::cout << "   <mismatches>               maximum number of allowed mismatches (must be < 256); applies only" << std::endl;
    std::cout << "                              to pigeon_hole, suffix_filter, min_u and 01 search schemes" << std::endl;
    std::cout << "   -o <output_file>           write pattern occurrences to this file (in ASCII format; one line per pattern)" << std::endl;
    std::cout << "   -sam <sam_file>            write occurrences in SAM format to <sam_file> (requires an index built" << std::endl;
    std::cout << "                              with move-rb-build -f; uses the per-occurrence CIGAR)" << std::endl;
    std::cout << "   -seq                       include the read sequence (SEQ field) in SAM records (default: '*')" << std::endl;
    std::cout << "   -norc                      in SAM mode, do not align the reverse complement of each read" << std::endl;
    std::cout << "                              (by default both strands are aligned; reverse hits get SAM FLAG 0x10)" << std::endl;
    std::cout << "   -fastq                     read <patterns_file> as FASTQ (variable-length reads); the read name" << std::endl;
    std::cout << "                              becomes QNAME and the qualities are written to the SAM QUAL field" << std::endl;
    std::cout << "                              (requires -sam and -seq)" << std::endl;
    std::cout << "   -best                      in SAM mode, report only one best (minimum-error) alignment per read" << std::endl;
    std::cout << "   -nocigar                   in SAM mode, omit the CIGAR (write '*'); faster, as the per-occurrence" << std::endl;
    std::cout << "                              CIGAR is then not computed" << std::endl;
    std::cout << "   -paf <paf_file>            write occurrences in PAF format to <paf_file> (one line per occurrence;" << std::endl;
    std::cout << "                              requires an index built with move-rb-build -f). May be combined with -sam." << std::endl;
    std::cout << "   <index_file>               index file (with extension .move-r)" << std::endl;
    std::cout << "   <patterns_file>            file in move-rb-patterns format containing the patterns." << std::endl;
    exit(0);
}

/**
 * @brief parses the next command-line argument(s)
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments
 * @return whether there are further options to parse
 */
bool parse_args(char** argv, int argc)
{
    if (arg_idx >= argc - min_args) return false;
    std::string s = argv[arg_idx];
    if (s == "-d") return false;
    arg_idx++;

    if (s == "-c") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -c option.");
        check_correctness = true;
        path_input_file = argv[arg_idx++];
        input_file.open(path_input_file);
        if (!input_file.good()) help("error: cannot open <input_file>");
    } else if (s == "-m") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -m option.");
        std::string path_m_file = argv[arg_idx++];
        mf.open(path_m_file, std::filesystem::exists(path_m_file) ? std::ios::app : std::ios::out);
        if (!mf.good()) help("error: cannot open nor create <m_file>");
        name_text_file = argv[arg_idx++];
    }  else if (s == "-o") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -o option.");
        output_occurrences = true;
        path_output_file = argv[arg_idx++];
    } else if (s == "-sam") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -sam option.");
        sam_output = true;
        path_sam_file = argv[arg_idx++];
    } else if (s == "-seq") {
        sam_include_seq = true;
    } else if (s == "-norc") {
        sam_rc = false;
    } else if (s == "-fastq") {
        fastq_input = true;
    } else if (s == "-best") {
        sam_best_only = true;
    } else if (s == "-nocigar") {
        sam_no_cigar = true;
    } else if (s == "-paf") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -paf option.");
        paf_output = true;
        path_paf_file = argv[arg_idx++];
    } else {
        help("error: unrecognized '" + s + "' option");
    }

    return true;
}

/**
 * @brief returns the next command-line argument, exiting with an error if there is none
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments
 * @param what a description of the expected argument (used in the error message)
 * @return the next argument
 */
const char* next_arg(int argc, char** argv, const std::string& what)
{
    if (arg_idx >= argc) help("error: missing " + what);
    return argv[arg_idx++];
}

/**
 * @brief reads the next FASTQ record (the four lines: @name, sequence, +, quality) from in
 * @param in the FASTQ stream
 * @param name set to the read name (the header after '@', up to the first whitespace)
 * @param seq set to the read sequence
 * @param qual set to the per-base quality string (same length as seq)
 * @return whether a record was read (false at end of file or on a malformed record)
 */
static bool read_fastq(std::istream& in, std::string& name, std::string& seq, std::string& qual)
{
    std::string h, plus;
    if (!std::getline(in, h) || !std::getline(in, seq) || !std::getline(in, plus) || !std::getline(in, qual))
        return false;

    auto strip_cr = [](std::string& s){ if (!s.empty() && s.back() == '\r') s.pop_back(); }; // tolerate CRLF
    strip_cr(h); strip_cr(seq); strip_cr(qual);
    if (h.empty() || h[0] != '@') return false;

    uint64_t sp = h.find_first_of(" \t", 1);
    name = sp == std::string::npos ? h.substr(1) : h.substr(1, sp - 1);
    return true;
}

/**
 * @brief loads the index and benchmarks locating the input patterns
 * @tparam pos_t index integer type
 * @tparam support the move-r locate-support type
 */
template <typename pos_t, move_r_support support, cigar_mode_t cm = CIGAR>
void measure_locate()
{
    std::cout << std::setprecision(4);
    auto time = now();
    log_phase_start(true, time, "loading the index");
    using idx_t = move_rb<support, char, pos_t>;
    idx_t index;
    index.load(index_file);
    index_file.close();
    log_phase_end(true, time);
    index.log_data_structure_sizes();

    std::cout << std::endl << "searching patterns ... " << std::endl;

    // FASTQ reads are variable-length and self-delimiting, so they carry no patterns-format header; otherwise the
    // header gives the (fixed) pattern count and length
    uint64_t num_patterns = 0;
    uint64_t pattern_length = 0;

    if (!fastq_input) {
        std::string header;
        std::getline(patterns_file, header);
        num_patterns = number_of_patterns(header);
        pattern_length = patterns_length(header);

        if (pattern_length < search_scheme.p)
            help("error: pattern length < number of parts defined in the search scheme");
    }

    uint64_t reads_processed = 0; // number of reads actually read (FASTQ count is not known in advance)
    uint64_t num_occurrences = 0;
    uint64_t time_locate = 0;
    std::chrono::steady_clock::time_point t2, t3;
    std::string pattern;
    no_init_resize(pattern, pattern_length);
    std::vector<aprx_occ_t<pos_t>> occurrences;
    using sam_list_t = std::vector<aprx_occ_t<pos_t, cm>>; // occurrences (each holding a CIGAR index unless -nocigar)
    sam_list_t fwd_results; // forward-strand occurrences of the current read
    sam_list_t rc_results;  // reverse-strand occurrences (of the read's reverse complement)
    // one CIGAR per SA-interval per strand; an occurrence's occ.cig_idx indexes into its strand's vector (CIGAR mode)
    std::vector<cigar_t<char>> fwd_cigars, rc_cigars;
    uint64_t checksum = 0;
    uint64_t baseline_alloc = malloc_count_current();
    malloc_count_reset_peak();
    progress_meter meter(fastq_input ? 0 : num_patterns); // FASTQ count is unknown up front, so disable the meter

    const auto& seqs = index.seq_data();
    pos_t seq_cursor = 0; // sequence-index cursor, advanced incrementally across the reported occurrences

    // writes one SAM alignment record of the current read. seq/qual are the read's forward-strand sequence and
    // qualities; for a reverse-strand record (reverse) the caller passes the reverse-complemented sequence and this
    // reverses the qualities to match. SEQ/QUAL are only written for the primary record and only when -seq is set.
    auto emit_record = [&](const std::string& qname, const aprx_occ_t<pos_t, cm>& o, const std::vector<cigar_t<char>>& cigars,
        bool reverse, bool primary, uint16_t mapq, const std::string& seq, const std::string& qual)
    {
        pos_t seq_i = seqs.sequence_index(o.pos, seq_cursor);
        pos_t local = o.pos - seqs.sequence_start(seq_i);
        uint16_t flag = (primary ? 0 : 0x100) | (reverse ? 0x10 : 0);

        sam_file << qname << '\t' << flag << '\t' << seqs.sequence_name(seq_i)
                 << '\t' << (local + 1) << '\t' << mapq << '\t';
        if constexpr (cm == CIGAR) append_cigar(sam_file, cigars[o.cig_idx]);
        else sam_file << '*'; // -nocigar: CIGAR omitted (not computed)
        sam_file << "\t*\t0\t0\t";

        // SEQ and QUAL are only written for the primary record and only when -seq is set; a reverse-strand record's
        // qualities are reversed to match the reported (reverse-complemented) sequence
        if (primary && sam_include_seq) {
            sam_file << seq << '\t';
            if (qual.empty()) sam_file << '*';
            else if (reverse) for (auto it = qual.rbegin(); it != qual.rend(); ++it) sam_file << *it;
            else sam_file << qual;
        } else {
            sam_file << "*\t*";
        }

        sam_file << "\tNM:i:" << o.err; // NM = edit distance of this alignment
        if constexpr (cm == CIGAR) { sam_file << "\tMD:Z:"; append_md(sam_file, cigars[o.cig_idx]); } // reference-mismatch string
        sam_file << '\n';
    };

    // emits the SAM records of one read: one primary = a minimum-error alignment over both strands, the rest
    // secondary (the search itself never crosses a separator, so every occurrence lies within one sequence), or a
    // single unmapped record
    auto emit_sam = [&](const std::string& qname, const std::string& fwd_seq, const std::string& rc_seq,
        const std::string& qual, sam_list_t& fwd, sam_list_t& rc,
        const std::vector<cigar_t<char>>& fwd_cig, const std::vector<cigar_t<char>>& rc_cig)
    {
        if (fwd.empty() && rc.empty()) { // unmapped read
            sam_file << qname << "\t4\t*\t0\t0\t*\t*\t0\t0\t"
                     << (sam_include_seq ? fwd_seq : std::string("*")) << "\t"
                     << (sam_include_seq && !qual.empty() ? qual : std::string("*")) << "\n";
            return;
        }

        // the primary alignment is the first minimum-error occurrence over both strands; MAPQ reflects uniqueness
        pos_t min_err = std::numeric_limits<pos_t>::max();
        for (const auto& o : fwd) min_err = std::min<pos_t>(min_err, o.err);
        for (const auto& o : rc)  min_err = std::min<pos_t>(min_err, o.err);
        pos_t count_min = 0;
        for (const auto& o : fwd) count_min += (o.err == min_err);
        for (const auto& o : rc)  count_min += (o.err == min_err);
        uint16_t primary_mapq = count_min == 1 ? 60 : 0;

        if (sam_best_only) {
            // report only one best alignment: the first minimum-error occurrence (forward strand preferred)
            for (const auto& o : fwd)
                if (o.err == min_err) { emit_record(qname, o, fwd_cig, false, true, primary_mapq, fwd_seq, qual); return; }
            for (const auto& o : rc)
                if (o.err == min_err) { emit_record(qname, o, rc_cig, true, true, primary_mapq, rc_seq, qual); return; }
            return;
        }

        bool primary_assigned = false;

        auto emit_list = [&](sam_list_t& occs, const std::vector<cigar_t<char>>& cigars, bool reverse, const std::string& seq) {
            for (const auto& o : occs) {
                bool primary = !primary_assigned && o.err == min_err;
                uint16_t mapq = primary ? primary_mapq : 0;
                primary_assigned |= primary;
                emit_record(qname, o, cigars, reverse, primary, mapq, seq, qual);
            }
        };

        emit_list(fwd, fwd_cig, false, fwd_seq);
        emit_list(rc, rc_cig, true, rc_seq);
    };

    if (sam_output) {
        if (!index.seq_data().has_sequences())
            help("error: -sam requires an index built with move-rb-build -f (no sequence names/boundaries stored)");

        sam_file.open(path_sam_file);
        if (!sam_file.good()) help("error: could not create <sam_file>");

        sam_file << "@HD\tVN:1.6\tSO:unsorted\n";
        for (pos_t i = 0; i < index.seq_data().num_sequences(); i++)
            sam_file << "@SQ\tSN:" << index.seq_data().sequence_name(i)
                     << "\tLN:" << index.seq_data().sequence_length(i) << "\n";
        sam_file << "@PG\tID:move-rb\tPN:move-rb-locate\tCL:" << command_line << "\n";
    }

    if (paf_output) {
        if (!index.seq_data().has_sequences())
            help("error: -paf requires an index built with move-rb-build -f (no sequence names/boundaries stored)");
        paf_file.open(path_paf_file);
        if (!paf_file.good()) help("error: could not create <paf_file>");
    }

    pos_t paf_cursor = 0; // sequence-index cursor for the PAF reporting loop

    // writes one PAF line per occurrence: query name/length, query span (the whole read), strand, target
    // name/length, target span, number of matching bases, alignment block length, mapping quality, then the NM
    // (edit distance) and (in CIGAR mode) cg tags
    auto emit_paf = [&](const std::string& qname, uint64_t read_len, sam_list_t& occs,
        const std::vector<cigar_t<char>>& cigars, bool reverse) {
        for (const auto& o : occs) {
            pos_t seq_i = seqs.sequence_index(o.pos, paf_cursor);
            pos_t local = o.pos - seqs.sequence_start(seq_i);
            pos_t matches, aln;

            if constexpr (cm == CIGAR) {
                matches = 0; aln = 0;
                for (const cigar_run_t<char>& run : cigars[o.cig_idx]) { aln += run.len; if (run.op == cigar_op_t::MATCH) matches += run.len; }
            } else { // no CIGAR: approximate from the occurrence's length and error
                aln = o.len; matches = o.len > o.err ? pos_t(o.len - o.err) : 0;
            }

            paf_file << qname << '\t' << read_len << "\t0\t" << read_len << '\t' << (reverse ? '-' : '+') << '\t'
                     << seqs.sequence_name(seq_i) << '\t' << seqs.sequence_length(seq_i) << '\t'
                     << local << '\t' << (local + o.len) << '\t' << matches << '\t' << aln << "\t255\tNM:i:" << o.err;
            if constexpr (cm == CIGAR) { paf_file << "\tcg:Z:"; append_cigar(paf_file, cigars[o.cig_idx]); }
            paf_file << '\n';
        }
    };

    std::string qname, qual; // current read's name and (FASTQ only) quality string

    for (uint64_t i = 0; ; i++) {
        if (fastq_input) {
            if (!read_fastq(patterns_file, qname, pattern, qual)) break;
        } else {
            if (i >= num_patterns) break;
            patterns_file.read(pattern.data(), pattern_length);
            qname = "read" + std::to_string(i);
        }
        reads_processed++;

        if (sam_output || paf_output) {
            t2 = now();

            // locate one strand's occurrences (with CIGAR), then sort by position and remove redundancy. A read too
            // short for the scheme (fewer characters than parts) cannot be searched, so it yields no occurrences.
            auto locate_strand = [&](const std::string& read, sam_list_t& out, std::vector<cigar_t<char>>& cigars_out) {
                out.clear();
                cigars_out.clear();
                if (k > 0 && read.size() <= search_scheme.p) return;
                auto collect = [&](aprx_occ_t<pos_t, cm> occ){ out.emplace_back(std::move(occ)); };

                if (dist_metr == HAMMING_DISTANCE)
                    cigars_out = index.template locate<HAMMING_DISTANCE, cm>(read, search_scheme, collect);
                else
                    cigars_out = index.template locate<EDIT_DISTANCE, cm>(read, search_scheme, collect);

                ips2ra::sort(out.begin(), out.end(), [](const auto& o){ return o.pos; });
                filter_edit_distance_occurrences<pos_t, cm>(out, pos_t(k));
            };

            locate_strand(pattern, fwd_results, fwd_cigars);
            std::string rc_seq;

            if (sam_rc) {
                rc_seq = reverse_complement(pattern);
                locate_strand(rc_seq, rc_results, rc_cigars);
            } else {
                rc_results.clear();
                rc_cigars.clear();
            }

            t3 = now();
            time_locate += time_diff_ns(t2, t3);

            num_occurrences += fwd_results.size() + rc_results.size();
            checksum += fwd_results.size() + rc_results.size();
            if (sam_output) emit_sam(qname, pattern, rc_seq, qual, fwd_results, rc_results, fwd_cigars, rc_cigars);
            if (paf_output) {
                emit_paf(qname, pattern.size(), fwd_results, fwd_cigars, false);
                emit_paf(qname, pattern.size(), rc_results, rc_cigars, true);
            }
            meter.step();
            continue;
        }

        t2 = now();

        if (dist_metr == HAMMING_DISTANCE) {
            index.template locate<HAMMING_DISTANCE>(pattern, search_scheme, [&](aprx_occ_t<pos_t> occ){occurrences.emplace_back(occ);});
        } else {
            index.template locate<EDIT_DISTANCE>(pattern, search_scheme, [&](aprx_occ_t<pos_t> occ){occurrences.emplace_back(occ);});
            ips2ra::sort(occurrences.begin(), occurrences.end(), [](const auto& o){ return o.pos; });
            filter_edit_distance_occurrences<pos_t>(occurrences, k);
        }
        
        num_occurrences += occurrences.size();

        t3 = now();
        time_locate += time_diff_ns(t2, t3);

        if (check_correctness) {
            if (dist_metr == HAMMING_DISTANCE) {
                ips2ra::sort(occurrences.begin(), occurrences.end(), [](const auto& o){ return o.pos; });
            }
            
            std::string occ_str;

            for (auto occ : occurrences) {
                no_init_resize(occ_str, occ.len);
                input_file.seekg(occ.pos, std::ios::beg);
                input_file.read(occ_str.data(), occ_str.size());
                pos_t best_dist;

                if (dist_metr == HAMMING_DISTANCE) {
                    best_dist = hamming_dist_bounded<pos_t>(occ_str, pattern, k);

                    if (occ.err != best_dist || occ.err > k) {
                        std::cout << "error: wrong approximate occurrence " << occ.pos << ", '" << occ_str << "' with length " << occ.len <<
                            " of pattern '" << pattern << "' was reported with " << occ.err << " instead of " << best_dist << " errors" << std::endl;
                        exit(-1);
                    }
                } else {
                    best_dist = edit_dist_bounded<pos_t>(occ_str, pattern, k);

                    if (dist_metr == EDIT_DISTANCE && !(best_dist == occ.err && occ.err <= k)) {
                        std::cout << "error: correct approximate occurrence " << occ.pos << ", '" <<
                            occ_str << "' with length " << occ.len << " of pattern '" << pattern <<
                            "' was reported with " << occ.err;
                            
                        if (best_dist > k) {
                            std::cout << " > k = " << k << " errors" << std::endl;
                        } else {
                            std::cout << " instead of " << best_dist << " errors" << std::endl;
                        }

                        exit(-1);
                    }
                }
            }
        }

        if (output_occurrences) {
            if (dist_metr == HAMMING_DISTANCE) {
                ips2ra::sort(occurrences.begin(), occurrences.end(), [](const auto& o){ return o.pos; });
            }

            for (const auto& occ : occurrences) {
                output_file <<
                    "(pos=" << occ.pos <<
                    " len=" << occ.len <<
                    " err=" << occ.err << ") ";
            }

            output_file << std::endl;
        }

        checksum += occurrences.size();

        for (auto occ : occurrences) {
            checksum += occ.pos + occ.err;
        }

        occurrences.clear();
        meter.step();
    }

    meter.finish();
    patterns_file.close();
    if (fastq_input) num_patterns = reads_processed; // the actual read count is only known after reading the file
    std::cout << "checksum: " << checksum << std::endl;
    if (MOVE_R_USE_MALLOC_COUNT)
        std::cout << "additional memory consumption during the search phase: " << format_size(malloc_count_peak() - baseline_alloc) << std::endl;
    if (num_patterns != 0)
        std::cout << "average occurrences per pattern: " << (num_occurrences / num_patterns) << std::endl;
    std::cout << "number of patterns: " << num_patterns << std::endl;
    if (!fastq_input) std::cout << "pattern length: " << pattern_length << std::endl;
    std::cout << "maximum number of mismatches: " << k << std::endl;
    std::cout << "search scheme: " << scheme_str << std::endl;
    std::cout << "total number of occurrences: " << num_occurrences << std::endl;
    std::cout << "locate time: " << format_time(time_locate) << std::endl;
    if (num_patterns != 0)
        std::cout << "             " << format_time(time_locate / num_patterns) << "/pattern" << std::endl;
    if (num_patterns != 0 && pattern_length != 0)
        std::cout << "             " << format_time(time_locate / (num_patterns * pattern_length)) << "/character" << std::endl;
    if (num_occurrences != 0)
        std::cout << "             " << format_time(time_locate / num_occurrences) << "/occurrence" << std::endl;

    if (mf.is_open()) {
        mf << "RESULT";
        mf << " algo=locate_move_rb_" << move_r_support_suffix(support);
        mf << " dist_metric=" << (dist_metr == HAMMING_DISTANCE ? "hamming" : "edit");
        mf << " text=" << name_text_file;
        mf << " n=" << index.forward_index().input_size();
        mf << " pattern_length=" << pattern_length;
        mf << " num_patterns=" << num_patterns;
        mf << " max_mismatches=" << k;
        mf << " num_occurrences=" << num_occurrences;
        mf << " time_locate=" << time_locate;
        index.log_data_structure_sizes(mf);
        mf << std::endl;
        mf.close();
    }

    if (check_correctness) input_file.close();
    if (sam_output) {
        sam_file.close();
        std::cout << "wrote SAM output to " << path_sam_file << std::endl;
    }

    if (paf_output) {
        paf_file.close();
        std::cout << "wrote PAF output to " << path_paf_file << std::endl;
    }
}

/**
 * @brief dispatches measure_locate on the CIGAR mode: -nocigar (SAM only) skips the per-occurrence CIGAR computation
 * @tparam pos_t index integer type
 * @tparam support the move-r locate-support type
 */
template <typename pos_t, move_r_support support>
void measure_locate_dispatch()
{
    if ((sam_output || paf_output) && sam_no_cigar) measure_locate<pos_t, support, NO_CIGAR>();
    else                            measure_locate<pos_t, support, CIGAR>();
}

/**
 * @brief program entry point
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments
 * @return the exit code
 */
int main(int argc, char** argv)
{
    if (argc - 1 < min_args) help("");

    for (int32_t i = 0; i < argc; i++) command_line += (i == 0 ? "" : " ") + std::string(argv[i]); // for the SAM @PG header

    while (parse_args(argv, argc));

    std::string arg = next_arg(argc, argv, "'-d <metric>'");
    if (arg != "-d") help("error: expected '-d <metric>'");
    std::string dist_str = next_arg(argc, argv, "<metric>");
    if      (dist_str == "hamming") dist_metr = HAMMING_DISTANCE;
    else if (dist_str == "edit")    dist_metr = EDIT_DISTANCE;
    else help("error: invalid option after -d");

    arg = next_arg(argc, argv, "'-s <scheme>'");
    if (arg != "-s") help("error: expected '-s <scheme>'");
    scheme_str = next_arg(argc, argv, "<scheme>");
    bool is_default_scheme =
        scheme_str == "pigeon_hole" ||
        scheme_str == "suffix_filter" ||
        scheme_str == "min_u" ||
        scheme_str == "01";
    if (!is_default_scheme) {
        if (std::filesystem::exists(scheme_str)) {
            std::string file_content;
            uint64_t file_size = std::filesystem::file_size(scheme_str);
            no_init_resize(file_content, file_size);
            std::ifstream ifile(scheme_str);
            ifile.read(file_content.data(), file_size);
            search_scheme = parse_search_scheme(file_content);
        } else help("error: invalid option after -s");
    }
    
    if (is_default_scheme) {
        arg = next_arg(argc, argv, "'-k <mismatches>'");
        if (arg != "-k") help("error: expected '-k <mismatches>'");
        int32_t k_arg = atoi(next_arg(argc, argv, "<mismatches>"));
        k = k_arg;
        if (k_arg < 0) help("error: invalid k value");
        if (k_arg >= 256) help("error: k must be < 256");
        if      (scheme_str == "pigeon_hole")   search_scheme = pigeon_hole_scheme(k);
        else if (scheme_str == "suffix_filter") search_scheme = suffix_filter_scheme(k);
        else if (scheme_str == "min_u")         search_scheme = min_u_scheme(k);
        else if (scheme_str == "01")            search_scheme = zero_one_scheme(k);
    } else {
        k = search_scheme.k;
    }
    
    if (dist_metr == EDIT_DISTANCE && k > 20)
        help("error: k > 20 is not supported for edit distance");


    if (sam_best_only && !sam_output)
        help("error: -best requires -sam");

    if (fastq_input && !((sam_output && sam_include_seq) || paf_output))
        help("error: -fastq requires either -paf, or -sam with -seq (so the qualities go to the SAM QUAL field)");

    if (arg_idx + 2 > argc) help("error: missing <index_file> and/or <patterns_file>");
    path_index_file = argv[arg_idx];
    path_patterns_file = argv[arg_idx + 1];

    index_file.open(path_index_file);
    patterns_file.open(path_patterns_file);

    if (!index_file.good()) help("error: could not read <index_file>");
    if (!patterns_file.good()) help("error: could not read <patterns_file>");

    if (output_occurrences) {
        output_file.open(path_output_file);
        if (!output_file.good()) help("error: could not create <output_file>");
    }

    bool is_64_bit;
    index_file.read((char*) &is_64_bit, 1);
    move_r_support _support;
    index_file.read((char*) &_support, sizeof(move_r_support));
    index_file.seekg(0, std::ios::beg);

    if (_support == _count) {
        help("error: this index does not support locate");
    } else if (_support == _locate_move) {
        if (is_64_bit) measure_locate_dispatch<uint64_t, _locate_move>();
        else           measure_locate_dispatch<uint32_t, _locate_move>();
    } else if (_support == _locate_rlzsa) {
        if (is_64_bit) measure_locate_dispatch<uint64_t, _locate_rlzsa>();
        else           measure_locate_dispatch<uint32_t, _locate_rlzsa>();
    }
    
    return 0;
}