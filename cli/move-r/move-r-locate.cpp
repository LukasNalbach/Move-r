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
#include <move_r/move_r.hpp>
#include <misc/fasta.hpp>
#include <misc/progress.hpp>

int arg_idx = 1;
bool output_occurrences = false;
bool check_correctness = false;
bool sam_output = false; // whether to write occurrences in SAM format (requires a FASTA-built index)
bool bed_output = false; // whether to write occurrences in BED format (requires a FASTA-built index)
bool locate_rc = true;   // whether to also match the reverse complement (only in SAM/BED, i.e. DNA, mode)
std::string input;
std::ofstream mf;
std::string path_index_file;
std::string path_patterns_file;
std::string path_output_file;
std::string path_sam_file;
std::string path_bed_file;
std::ifstream index_file;
std::ifstream patterns_file;
std::ifstream input_file;
std::ofstream output_file;
std::ofstream sam_file;
std::ofstream bed_file;
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
    std::cout << "move-r-locate: locate all occurrences of the input patterns." << std::endl << std::endl;
    std::cout << "usage: move-r-locate [...] <index_file> <patterns>" << std::endl;
    std::cout << "   -m <m_file> <text_name>    m_file is the file to write measurement data to," << std::endl;
    std::cout << "                              text_name should be the name of the original file" << std::endl;
    std::cout << "   -c <input_file>            checks correctness of each pattern occurrence against <input_file>" << std::endl;
    std::cout << "                              (input_file must be the file the index was built for)" << std::endl;
    std::cout << "   -o <output_file>           write pattern occurrences to this file (in ASCII format; one line per pattern)" << std::endl;
    std::cout << "   -sam <sam_file>            write occurrences in SAM format to <sam_file> (requires an index built" << std::endl;
    std::cout << "                              from FASTA with move-r-build -f); the i-th pattern is named pat<i>," << std::endl;
    std::cout << "                              each occurrence is an exact forward alignment (CIGAR <m>=)" << std::endl;
    std::cout << "   -bed <bed_file>            write occurrences in BED format to <bed_file> (0-based half-open" << std::endl;
    std::cout << "                              intervals; requires a FASTA-built index)" << std::endl;
    std::cout << "   -norc                      in SAM/BED mode, do not also match the reverse complement of each" << std::endl;
    std::cout << "                              pattern (by default both strands are matched; reverse hits get strand '-')" << std::endl;
    std::cout << "   <index_file>               index file (with extension .move-r)" << std::endl;
    std::cout << "   <patterns_file>            file in pizza&chili format containing the patterns" << std::endl;
    exit(0);
}

/**
 * @brief parses the next command-line argument(s)
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments
 */
void parse_args(char** argv, int argc)
{
    std::string s = argv[arg_idx];
    arg_idx++;

    if (s == "-c") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -c option.");
        check_correctness = true;
        path_input_file = argv[arg_idx++];
        input_file.open(path_input_file);
        if (!input_file.good()) help("error: cannot open <input_file>");
    } else if (s == "-m") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -o option.");
        std::string path_m_file = argv[arg_idx++];
        mf.open(path_m_file, std::filesystem::exists(path_m_file) ? std::ios::app : std::ios::out);
        if (!mf.good()) help("error: cannot open nor create <m_file>");
        name_text_file = argv[arg_idx++];
    } else if (s == "-o") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -o option.");
        output_occurrences = true;
        path_output_file = argv[arg_idx++];
    } else if (s == "-sam") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -sam option.");
        sam_output = true;
        path_sam_file = argv[arg_idx++];
    } else if (s == "-bed") {
        if (arg_idx >= argc - 1) help("error: missing parameter after -bed option.");
        bed_output = true;
        path_bed_file = argv[arg_idx++];
    } else if (s == "-norc") {
        locate_rc = false;
    } else {
        help("error: unrecognized '" + s + "' option");
    }
}

/**
 * @brief loads the index and benchmarks locating the input patterns
 * @tparam pos_t index integer type
 * @tparam support the move-r locate-support type
 */
template <typename pos_t, move_r_support support>
void measure_locate()
{
    std::cout << std::setprecision(4);
    auto time = now();
    log_phase_start(true, time, "loading the index");
    using idx_t = move_r<support, char, pos_t>;
    idx_t index;
    index.load(index_file);
    log_phase_end(true, time);
    index_file.close();
    index.log_data_structure_sizes();

    if (check_correctness) {
        if (path_input_file == "") help("error: <input_file> not provided");
        log_phase_start(true, time, "\nloading input file");
        no_init_resize(input, index.input_size());
        read_from_file(input_file, input.data(), input.size());
        input_file.close();
        log_phase_end(true, time);
    }

    const auto& seqs = index.seq_data();
    if ((sam_output || bed_output) && !seqs.has_sequences())
        help("error: -sam/-bed require an index built from FASTA with move-r-build -f");

    std::cout << std::endl << "searching patterns ... " << std::endl;
    std::string header;
    std::getline(patterns_file, header);
    uint64_t num_patterns = number_of_patterns(header);
    uint64_t pattern_length = patterns_length(header);

    if (sam_output) {
        // SAM header: version, one @SQ line per sequence, and the program record
        sam_file << "@HD\tVN:1.6\tSO:unsorted\n";
        for (pos_t i = 0; i < seqs.num_sequences(); i++)
            sam_file << "@SQ\tSN:" << seqs.sequence_name(i) << "\tLN:" << seqs.sequence_length(i) << "\n";
        sam_file << "@PG\tID:move-r\tPN:move-r-locate\tCL:" << command_line << "\n";
    }

    pos_t seq_cursor = 0; // sequence-index cursor, advanced incrementally across the reported occurrences
    uint64_t num_occurrences = 0;
    uint64_t time_locate = 0;
    std::chrono::steady_clock::time_point t2, t3;
    std::string pattern;
    no_init_resize(pattern, pattern_length);
    std::vector<pos_t> occurrences;    // forward-strand occurrences of the current pattern
    std::vector<pos_t> rc_occurrences; // reverse-strand occurrences (of the pattern's reverse complement)
    // whether to also match the reverse complement: only meaningful for the biological (SAM/BED, i.e. DNA) output
    bool rc_mode = (sam_output || bed_output) && locate_rc;
    progress_meter meter(num_patterns);

    for (uint64_t i = 0; i < num_patterns; i++) {
        patterns_file.read(pattern.data(), pattern_length);
        t2 = now();
        occurrences = index.locate(pattern);
        if (rc_mode) rc_occurrences = index.locate(reverse_complement(pattern));
        t3 = now();
        time_locate += time_diff_ns(t2, t3);
        num_occurrences += occurrences.size() + rc_occurrences.size();

        if (check_correctness) {
            for (pos_t occ : occurrences) {
                if (input.substr(occ, pattern_length) != pattern) {
                    std::cout << "error: wrong occurrence: " << occ << ", '" <<
                        input.substr(occ, pattern_length) <<  "' of pattern '" <<
                        pattern << "'" << std::endl;
                    exit(-1);
                }
            }
        }

        if (output_occurrences || sam_output || bed_output)
            ips2ra::sort(occurrences.begin(), occurrences.end());
        if (rc_mode && !rc_occurrences.empty())
            ips2ra::sort(rc_occurrences.begin(), rc_occurrences.end());

        if (output_occurrences) {
            for (pos_t occ : occurrences) output_file << occ << " ";
            output_file << std::endl;
        }

        if (sam_output || bed_output) {
            // the i-th anonymous pizza&chili pattern is named pat<i>; each occurrence is an exact hit (CIGAR <m>=),
            // forward-strand ('+') for pattern hits and reverse-strand ('-') for reverse-complement hits
            std::string qname = "pat" + std::to_string(i);
            uint16_t mapq = (occurrences.size() + rc_occurrences.size()) == 1 ? 60 : 0; // unique => 60, multi => 0
            bool primary = true;

            auto emit = [&](const std::vector<pos_t>& occs, bool reverse) {
                for (pos_t occ : occs) {
                    pos_t seq_i = seqs.sequence_index(occ, seq_cursor);
                    pos_t local = occ - seqs.sequence_start(seq_i);

                    if (bed_output) // BED: 0-based, half-open [local, local + m)
                        bed_file << seqs.sequence_name(seq_i) << '\t' << local << '\t' << (local + pattern_length)
                                 << '\t' << qname << "\t0\t" << (reverse ? '-' : '+') << '\n';

                    if (sam_output) { // SAM: 1-based POS, exact CIGAR "<m>=", 0x10 for reverse, first hit primary
                        uint16_t flag = (primary ? 0 : 0x100) | (reverse ? 0x10 : 0);
                        sam_file << qname << '\t' << flag << '\t' << seqs.sequence_name(seq_i)
                                 << '\t' << (local + 1) << '\t' << (primary ? mapq : 0) << '\t'
                                 << pattern_length << "=\t*\t0\t0\t*\t*\n";
                    }

                    primary = false;
                }
            };

            emit(occurrences, false);
            if (rc_mode) emit(rc_occurrences, true);
        }

        occurrences.clear();
        rc_occurrences.clear();
        meter.step();
    }

    meter.finish();
    std::cout << "average occurrences per pattern: " << (num_occurrences / num_patterns) << std::endl;
    std::cout << "number of patterns: " << num_patterns << std::endl;
    std::cout << "pattern length: " << pattern_length << std::endl;
    std::cout << "total number of occurrences: " << num_occurrences << std::endl;
    std::cout << "locate time: " << format_time(time_locate) << std::endl;
    std::cout << "             " << format_time(time_locate / num_patterns) << "/pattern" << std::endl;
    if (num_occurrences != 0)
      std::cout << "             " << format_time(time_locate / num_occurrences) << "/occurrence" << std::endl;

    if (mf.is_open()) {
        mf << "RESULT";
        mf << " algo=locate_move_r_" << move_r_support_suffix(support);
        mf << " text=" << name_text_file;
        mf << " n=" << index.input_size();
        mf << " pattern_length=" << pattern_length;
        mf << " num_patterns=" << num_patterns;
        mf << " num_occurrences=" << num_occurrences;
        mf << " time_locate=" << time_locate;
        index.log_data_structure_sizes(mf);
        mf << std::endl;
        mf.close();
    }
}

/**
 * @brief program entry point
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments
 * @return the exit code
 */
int main(int argc, char** argv)
{
    if (argc < 3) help("");
    for (int i = 0; i < argc; i++) command_line += (i == 0 ? "" : " ") + std::string(argv[i]); // for the SAM @PG header
    while (arg_idx < argc - 2) parse_args(argv, argc);
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

    if (sam_output) {
        sam_file.open(path_sam_file);
        if (!sam_file.good()) help("error: could not create <sam_file>");
    }

    if (bed_output) {
        bed_file.open(path_bed_file);
        if (!bed_file.good()) help("error: could not create <bed_file>");
    }

    bool is_64_bit;
    index_file.read((char*) &is_64_bit, 1);
    move_r_support _support;
    index_file.read((char*) &_support, sizeof(move_r_support));
    index_file.seekg(0, std::ios::beg);

    if (_support == _count || _support == _locate_one) {
        help("error: this index does not support locate");
    } else if (_support == _locate_move) {
        if (is_64_bit) measure_locate<uint64_t, _locate_move>();
        else           measure_locate<uint32_t, _locate_move>();
    } else if (_support == _locate_rlzsa) {
        if (is_64_bit) measure_locate<uint64_t, _locate_rlzsa>();
        else           measure_locate<uint32_t, _locate_rlzsa>();
    }

    return 0;
}