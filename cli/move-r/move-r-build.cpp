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
#include <vector>
#include <misc/fasta.hpp>
#include <move_r/move_r.hpp>

int arg_idx = 1;
uint64_t n;
uint16_t a = 8;
uint16_t p = omp_get_max_threads();
std::string path_prefix_index_file;
move_r_construction_mode mode = _suffix_array;
move_r_support support = _locate_move;
std::ofstream mf_idx;
std::ofstream mf_mds;
std::ofstream index_file;
std::vector<std::string> path_input_files; // the input file(s); more than one is only allowed in FASTA mode
std::string name_text_file;
std::string path_index_file;
bool fasta_mode = false;              // whether the input is a (multi-sequence) FASTA file (DNA mode)
std::string fasta_allowed = "ACGT";   // the alphabet kept verbatim in FASTA mode; other bases are masked

/**
 * @brief prints the usage information and exits
 * @param msg an optional error message printed before the usage information
 */
void help(std::string msg)
{
    if (msg != "") std::cout << msg << std::endl;
    std::cout << "move-r-build: builds move-r." << std::endl << std::endl;
    std::cout << "usage: move-r-build [...] <input_file> [<input_file> ...]" << std::endl;
    std::cout << "   -c <mode>           construction mode: sa or bigbwt (default: sa)" << std::endl;
    std::cout << "   -o <base_name>      names the index file base_name.move-r(-rlzsa) (default: input_file)" << std::endl;
    std::cout << "   -s <support>        support: count, locate_move or locate_rlzsa" << std::endl;
    std::cout << "                       (default: locate_move)" << std::endl;
    std::cout << "   -p <integer>        number of threads to use during the construction of the index" << std::endl;
    std::cout << "                       (default: all threads)" << std::endl;
    std::cout << "   -a <integer>        balancing parameter; a must be an integer number and a >= 2 (default: 8)" << std::endl;
    std::cout << "   -m_idx <m_file_idx> m_file_idx is file to write measurement data of the index construction to" << std::endl;
    std::cout << "   -m_mds <m_file_mds> m_file_mds is file to write measurement data of the construction of the move" << std::endl;
    std::cout << "                       data structures to" << std::endl;
    std::cout << "   -f                  FASTA/DNA mode: read the <input_file>(s) as FASTA (strip headers, mask non-" << std::endl;
    std::cout << "                       allowed bases); multiple FASTA files may be given and are concatenated" << std::endl;
    std::cout << "   -A <alphabet>       in FASTA mode, the alphabet kept verbatim (default: ACGT); any other base is" << std::endl;
    std::cout << "                       replaced with '" << fasta_mask_symbol << "' and cannot be matched" << std::endl;
    std::cout << "   <input_file>        input file (multiple only in FASTA mode)" << std::endl;
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

    if (s == "-o") {
        if (arg_idx >= argc) help("error: missing parameter after -o option");
        path_prefix_index_file = argv[arg_idx++];
    } else if (s == "-p") {
        if (arg_idx >= argc) help("error: missing parameter after -p option");
        p = atoi(argv[arg_idx++]);
        if (p < 1) help("error: p < 1");
        if (p > omp_get_max_threads()) help("error: p > maximum number of threads");
    } else if (s == "-c") {
        if (arg_idx >= argc) help("error: missing parameter after -p option");
        std::string construction_mode_str = argv[arg_idx++];
        if (construction_mode_str == "sa") mode = _suffix_array;
        else if (construction_mode_str == "bigbwt") mode = _bigbwt;
        else help("error: invalid option for -c");
    } else if (s == "-s") {
        if (arg_idx >= argc) help("error: missing parameter after -s option");
        std::string support_str = argv[arg_idx++];
        if (support_str == "count") {
            support = _count;
        } else if (support_str == "locate_one") {
            support = _locate_one;
        } else if (support_str == "locate_move") {
            support = _locate_move;
        } else if (support_str == "locate_rlzsa") {
            support = _locate_rlzsa;
        } else help("error: unknown mode provided with -s option");
    } else if (s == "-a") {
        if (arg_idx >= argc) help("error: missing parameter after -a option");
        a = atoi(argv[arg_idx++]);
        if (a < 2) help("error: a < 2");
    } else if (s == "-m_idx") {
        if (arg_idx >= argc) help("error: missing parameter after -m_idx option");
        std::string path_mf_idx = argv[arg_idx++];
        mf_idx.open(path_mf_idx, std::filesystem::exists(path_mf_idx) ? std::ios::app : std::ios::out);
        if (!mf_idx.good()) help("error: cannot open nor create <m_file_idx>");
    } else if (s == "-m_mds") {
        if (arg_idx >= argc) help("error: missing parameter after -m_mds option");
        std::string path_mf_mds = argv[arg_idx++];
        mf_mds.open(path_mf_mds, std::filesystem::exists(path_mf_mds) ? std::ios::app : std::ios::out);
        if (!mf_mds.good()) help("error: cannot open nor create <m_file_mds>");
    } else if (s == "-f") {
        fasta_mode = true;
    } else if (s == "-A") {
        if (arg_idx >= argc) help("error: missing parameter after -A option");
        fasta_allowed = argv[arg_idx++];
    } else {
        help("error: unrecognized '" + s + "' option");
    }
}

/**
 * @brief builds the index and writes it to disk
 * @tparam pos_t index integer type
 * @tparam support the move-r locate-support type
 */
template <typename pos_t, move_r_support support>
void build()
{
    fasta_sequence_data<pos_t, char> fasta_data;       // the per-sequence metadata (FASTA mode only)
    std::string path_build_file = path_input_files[0]; // the file the index is built from (single non-FASTA input)
    std::string path_fasta_text;                       // the temporary preprocessed-text file (FASTA mode only)

    if (fasta_mode) {
        // stream the FASTA file(s) to a single temporary on-disk text (headers stripped, non-allowed bases masked);
        // multiple files are concatenated into one text, then the index is built from it as for a plain on-disk input
        path_fasta_text = path_index_file + ".tmp_text";

        auto time = now();
        log_phase_start(true, time, "preprocessing the FASTA file(s)");
        fasta_data = process_fasta<pos_t>(path_input_files, path_fasta_text, fasta_allowed);
        log_phase_end(true, time);
        std::cout << "read " << fasta_data.num_sequences() << " sequence(s) from "
                  << path_input_files.size() << " file(s)" << std::endl;

        path_build_file = path_fasta_text;
    }

    move_r<support, char, pos_t> index(path_build_file, {
        .file_input = true,
        .mode = mode,
        .num_threads = p,
        .a = a,
        .log = true,
        .mf_idx = mf_idx.is_open() ? &mf_idx : nullptr,
        .mf_mds = mf_mds.is_open() ? &mf_mds : nullptr,
        .name_text_file = name_text_file
    });

    if (fasta_mode) {
        // attach the per-sequence names/boundaries for sequence-relative SAM coordinates, then drop the temp file
        index.set_fasta_sequence_data(std::move(fasta_data));
        std::filesystem::remove(path_fasta_text);
    }

    auto time = now();
    log_phase_start(true, time, "serializing the index");
    index.serialize(index_file);
    log_phase_end(true, time);
}

/**
 * @brief program entry point
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments
 * @return the exit code
 */
int main(int argc, char** argv)
{
    if (argc < 2) help("");
    // parse the options, then collect all remaining (non-option) arguments as input file(s)
    while (arg_idx < argc && argv[arg_idx][0] == '-') parse_args(argv, argc);
    while (arg_idx < argc) path_input_files.emplace_back(argv[arg_idx++]);

    if (path_input_files.empty()) help("error: no <input_file> given");
    if (!fasta_mode && path_input_files.size() > 1)
        help("error: multiple input files are only supported in FASTA mode (-f)");
    for (const std::string& f : path_input_files)
        if (!std::filesystem::exists(f) || !std::filesystem::is_regular_file(f))
            help("error: <input_file> does not exist: " + f);
    if (path_prefix_index_file == "")
        path_prefix_index_file = path_input_files[0];

    std::cout << std::setprecision(4);
    name_text_file = path_input_files[0].substr(path_input_files[0].find_last_of("/\\") + 1);
    path_index_file = path_prefix_index_file.append(".move-r");
    if (support == _locate_rlzsa) path_index_file = path_index_file.append("-rlzsa");

    index_file.open(path_index_file);
    if (!index_file.good()) help("error: invalid input, could not create <index_file>");

    // an upper bound on the produced text length: the total size of the input file(s) (+1 for the sentinel); in FASTA
    // mode the actual text is smaller (headers become one byte each), but this suffices for the heuristics below
    n = 1;
    for (const std::string& f : path_input_files) n += std::filesystem::file_size(f);

    if (p > 1 && 1000 * p > n) {
        p = std::max<uint16_t>(1, n / 1000);
        std::cout << "n = " << n << ", warning: p > n/1000, setting p to n/1000 ~ " << std::to_string(p) << std::endl;
    } else {
        p = std::max<uint16_t>(1, std::min<uint64_t>({ (uint64_t)omp_get_max_threads(), n / 1000, p }));
    }

    std::cout << "building move-r of " << name_text_file << (path_input_files.size() > 1
        ? " (+" + std::to_string(path_input_files.size() - 1) + " more file(s))" : "");
    std::cout << " using " << format_threads(p) << " and a = " << a << std::endl;
    std::cout << "the index will be saved to " << path_index_file << std::endl << std::endl;

    if (mf_idx.is_open()) {
        mf_idx << "RESULT"
            << " algo=build_move_r_" << move_r_support_suffix(support)
            << " text=" << name_text_file
            << " num_threads=" << p
            << " a=" << a;
    }

    if (support == _count) {
        if (n < UINT_MAX) {
            build<uint32_t, _count>();
        } else {
            build<uint64_t, _count>();
        }
    } else if (support == _locate_move) {
        if (n < UINT_MAX) {
            build<uint32_t, _locate_move>();
        } else {
            build<uint64_t, _locate_move>();
        }
    } else if (support == _locate_rlzsa) {
        if (n < UINT_MAX) {
            build<uint32_t, _locate_rlzsa>();
        } else {
            build<uint64_t, _locate_rlzsa>();
        }
    }

    return 0;
}