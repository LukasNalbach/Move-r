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
 * @file build_cli.hpp
 * @brief the options, the input files, the FASTA preprocessing and the log of the build modes of move-r and
 *        move-rb
 */

#ifndef MOVE_R_MISC_BUILD_CLI_HPP
#define MOVE_R_MISC_BUILD_CLI_HPP

#include <omp.h>

#include <algorithm>
#include <array>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <data_structures/fasta_sequence_data.hpp>
#include <misc/cli_help.hpp>
#include <misc/fasta.hpp>
#include <misc/files.hpp>
#include <misc/log.hpp>
#include <misc/parse.hpp>
#include <misc/temp_files.hpp>
#include <move_r/move_r.hpp>

// ############################# OPTIONS #############################

// the construction mode of the build modes if -c is not given
inline constexpr move_r_construction_mode default_construction_mode = _pfp;

// the alphabet of -A if it is not given
inline constexpr std::string_view default_fasta_alphabet = "ACGT";

// the options shared by the build modes; a build mode derives its options from it
struct build_cli_t {
    uint64_t balancing_parameter = move_r_params().a; // -balance
    uint64_t num_threads = omp_get_max_threads(); // -p
    move_r_construction_mode mode = default_construction_mode; // -c
    bool construction_mode_given = false; // whether -c was given (see build_with_fallback)
    bool fasta_mode = false; // -f
    std::string fasta_alphabet = std::string(default_fasta_alphabet); // -A
    bool reorder_alphabet = true; // false <=> -keep-order (see move_r_params::reorder_alphabet)
    std::string index_base_name; // -o (default: the first input file)
    // the input files; more than one only in FASTA mode, then they are concatenated
    std::vector<std::string> path_input_files;
    std::string text_name; // the <name> of -m (default: the file name of the first input file)
    direct_ofstream measurement_file; // -m
    direct_ofstream measurement_file_mds; // -m-mds
    std::string path_index_file; // the file the index is written to
    direct_ofstream index_file; // the stream of path_index_file
    uint64_t input_bytes = 0; // the total size of the input files in bytes

    // returns the measurement file of -m (nullptr if -m is not given)
    std::ostream* measurement() { return measurement_file.is_open() ? &measurement_file : nullptr; }

    // returns the measurement file of -m-mds (nullptr if -m-mds is not given)
    std::ostream* measurement_mds() { return measurement_file_mds.is_open() ? &measurement_file_mds : nullptr; }
};

/**
 * @brief reads the value of an option as a decimal number >= 0 (see parse_double)
 * @param args the arguments, positioned after the option
 * @param option the option (for the error message)
 * @throws usage_error "<option> expects a number >= 0"
 */
inline double read_nonnegative_number(cli_arguments& args, const std::string& option)
{
    double value = 0;
    if (!parse_double(args.value(option), value) || !(value >= 0)) throw usage_error(option + " expects a number >= 0");
    return value;
}

/**
 * @brief reads one of the options -o, -p, -balance, -c, -m, -m-mds, -f, -A and -keep-order
 * @param args the arguments, positioned after the option
 * @param option the option
 * @param cli the options
 * @param with_mds_measurement whether the mode has the option -m-mds
 * @return whether option is one of these options
 */
inline bool read_build_option(cli_arguments& args, const std::string& option, build_cli_t& cli,
    bool with_mds_measurement)
{
    if (option == "-o") {
        cli.index_base_name = args.value(option);
    } else if (option == "-p") {
        cli.num_threads = args.integer(option, 1, omp_get_max_threads());
    } else if (option == "-balance") {
        cli.balancing_parameter = args.integer(option, min_balancing_parameter);
    } else if (option == "-c") {
        const std::string mode = args.value(option);
        if (mode == "pfp") cli.mode = _pfp;
        else if (mode == "sa") cli.mode = _suffix_array;
        else if (mode == "sa-space") cli.mode = _suffix_array_space;
        else throw usage_error("-c expects pfp, sa or sa-space");
        cli.construction_mode_given = true;
    } else if (option == "-m") {
        read_measurement_option(args, cli.measurement_file, cli.text_name);
    } else if (option == "-m-mds" && with_mds_measurement) {
        const std::string path = args.value(option);
        if (!open_measurement_file(cli.measurement_file_mds, path))
            throw std::runtime_error("cannot open or create " + path);
    } else if (option == "-f") {
        cli.fasta_mode = true;
    } else if (option == "-A") {
        cli.fasta_alphabet = args.value(option);
    } else if (option == "-keep-order") {
        cli.reorder_alphabet = false;
    } else {
        return false;
    }

    return true;
}

/**
 * @brief reads the input files (the remaining arguments) and sets index_base_name (if -o is not given), text_name
 *        (if -m is not given) and input_bytes
 * @param args the arguments, positioned at the first input file
 * @param cli the options
 * @throws usage_error if there is no input file, or several without -f; std::runtime_error "cannot read <path>"
 */
inline void read_input_files(cli_arguments& args, build_cli_t& cli)
{
    while (!args.empty()) cli.path_input_files.push_back(args.next());
    if (cli.path_input_files.empty()) throw usage_error("expected <input_file>");
    if (!cli.fasta_mode && cli.path_input_files.size() > 1) throw usage_error("several input files need -f");

    for (const std::string& file : cli.path_input_files) {
        if (!std::filesystem::is_regular_file(file)) throw std::runtime_error("cannot read " + file);
        cli.input_bytes += std::filesystem::file_size(file);
    }

    if (cli.index_base_name.empty()) cli.index_base_name = cli.path_input_files[0];
    if (cli.text_name.empty()) cli.text_name = std::filesystem::path(cli.path_input_files[0]).filename().string();
}

// ############################# HELP #############################

/**
 * @brief prints the options -f and -A
 * @param h the help printer
 * @param fasta_explanation the explanation of -f
 */
inline void print_fasta_options(help_printer& h, const std::string& fasta_explanation)
{
    h.option("-f", "", fasta_explanation);
    h.option("-A <alphabet>", "a string of symbols (default: " + std::string(default_fasta_alphabet) + ")",
        "with -f: the symbols kept verbatim");
}

/**
 * @brief prints the options -c, -p and -balance
 * @param h the help printer
 * @param pfp_explanation the explanation of the construction mode pfp
 */
inline void print_execution_options(help_printer& h, const std::string& pfp_explanation)
{
    // the limits of the modes in input symbols: 2^sa_bits - 1 for sa and sa-space, 2^48 - 2 for pfp
    const std::string sa_bits = std::to_string(std::bit_width(move_r_max_input_length_sa));
    h.option("-c <mode>", "pfp | sa | sa-space (default: pfp)", "pfp: " + pfp_explanation + "; inputs of at most 2^"
        + std::to_string(std::bit_width(move_r_max_input_length)) + " - 2 symbols\n"
        "sa: in memory with the suffix array; inputs of at most 2^" + sa_bits + " - 1 symbols\n"
        "sa-space: like sa, with idle data structures written to disk\n"
        "without -c, a failed pfp construction is repeated with sa-space if the input has at most 2^" + sa_bits
        + " - 1 symbols");
    h.option_threads();
    h.option("-balance <a>", "integer >= " + std::to_string(min_balancing_parameter) + " (default: "
        + std::to_string(move_r_params().a) + ")",
        "the balancing parameter of the move data structures");
}

/**
 * @brief prints the options -m and -m-mds
 * @param h the help printer
 * @param with_mds_measurement whether the mode has the option -m-mds
 */
inline void print_build_measurement_options(help_printer& h, bool with_mds_measurement)
{
    h.option_measurement();

    if (with_mds_measurement) {
        h.option("-m-mds <file>", "a path", "append measurement data of the construction of the move data "
            "structures to this file");
    }
}

// ############################# LOG #############################

/**
 * @brief limits the number of threads to max(1, n / move_r_min_symbols_per_thread), opens the output file and logs
 *        the start of the build: "building <what> of <file> (<size>) using <p> threads, <parameters>" and
 *        "output: <path>"; n is not known before the construction, its upper bound input_bytes + 1 is used
 * @param cli the options; path_index_file must be set
 * @param what what is built, e.g. "the move-r index"
 * @param parameters the parameters of the build, e.g. "support locate-move, a = 8"
 * @throws std::runtime_error "cannot create <path>"
 */
inline void begin_build(build_cli_t& cli, const std::string& what, const std::string& parameters)
{
    const uint64_t max_input_length = cli.input_bytes + 1;

    if (cli.num_threads > 1 && move_r_min_symbols_per_thread * cli.num_threads > max_input_length) {
        const uint64_t num_threads = std::max<uint64_t>(1, max_input_length / move_r_min_symbols_per_thread);
        log_out() << "warning: using " << format_threads(num_threads) << " instead of "
                  << format_number(cli.num_threads) << " (at most n/" << move_r_min_symbols_per_thread
                  << " threads, here n <= " << format_number(max_input_length) << ")" << std::endl;
        cli.num_threads = num_threads;
    }

    cli.index_file.open(cli.path_index_file);
    require_created(cli.index_file, cli.path_index_file);

    // real numbers are logged with 4 significant digits
    log_out() << std::setprecision(4);
    const std::string input_name = std::filesystem::path(cli.path_input_files[0]).filename().string();
    const uint64_t more_files = cli.path_input_files.size() - 1;
    log_out() << "building " << what << " of " << input_name;
    if (more_files > 0) log_out() << " and " << more_files << (more_files == 1 ? " more file" : " more files");
    log_out() << " (" << format_size(cli.input_bytes) << ") using " << format_threads(cli.num_threads) << ", "
              << parameters << std::endl
              << "output: " << cli.path_index_file << std::endl << std::endl;
}

/**
 * @brief writes the index to cli.path_index_file in a logged phase (its size is logged by the construction, before
 *        the sizes of its data structures)
 * @tparam index_t type of the index
 * @param index the index
 * @param cli the options
 */
template <typename index_t>
void write_index(const index_t& index, build_cli_t& cli)
{
    auto time = now();
    log_phase_start(true, time, "writing " + cli.path_index_file);
    index.serialize(cli.index_file);
    cli.index_file.close();
    log_phase_end(true, time, cli.measurement(), "time_serialize_index");
}

/**
 * @brief builds the index; if the construction with the prefix-free parsing fails, the construction mode has not
 *        been chosen with -c and the input is short enough for sa-space (see move_r_max_input_length_sa), the index
 *        is built with the construction mode sa-space
 * @tparam build_t type of build
 * @param mode the construction mode
 * @param construction_mode_given whether the construction mode has been chosen with -c
 * @param input_length the number of input symbols (see fasta_cli_t::input_length)
 * @param build build(mode) builds the index with the construction mode mode
 * @return the result of build
 */
template <typename build_t>
decltype(auto) build_with_fallback(move_r_construction_mode mode, bool construction_mode_given,
    uint64_t input_length, build_t build)
{
    try {
        return build(mode);
    } catch (const std::runtime_error& error) {
        if (mode != _pfp || construction_mode_given || input_length > move_r_max_input_length_sa) throw;
        end_phase_line();
        log_out() << std::endl << error.what() << std::endl
                  << "falling back to the construction mode sa-space" << std::endl;
        return build(_suffix_array_space);
    }
}

// ############################# FASTA INPUT #############################

/**
 * @brief the FASTA input of a build: its sequences, scanned once (see scan); a construction reads a single FASTA
 *        file itself, several FASTA files are concatenated into a temporary file that stores the preprocessed text
 *        (see construction_input_file)
 */
struct fasta_cli_t {
    build_cli_t& cli; // the options
    fasta_sequence_data data; // the names and the boundaries of the sequences
    // the frequencies of the bytes in the text: [0..255] freq[c] = number of occurrences of the byte c
    std::array<uint64_t, 256> freq {};
    uint64_t text_size = 0; // length of the text
    // sampled text positions with their offsets in the FASTA file (see pfp_rlbwt::fasta_sample)
    std::vector<pfp_rlbwt::fasta_sample> samples;
    std::string path_text; // path of the preprocessed text file (empty if none has been written)
    bool scanned = false; // false <=> freq and samples are not known (see scan)

    explicit fasta_cli_t(build_cli_t& cli) : cli(cli) {}

    /* scans the FASTA file(s) (does nothing without -f); a single FASTA file whose cache of the pfp output stores its
       sequences is not scanned (see pfp_cache_fasta_section) */
    void scan()
    {
        if (!cli.fasta_mode) return;

        if (cli.mode == _pfp && cli.path_input_files.size() == 1) {
            const std::string cache = move_r<>::pfp_cache_path(cli.path_input_files[0], cli.mode);

            if (read_fasta_from_pfp_cache(cache, text_size, data)) {
                log_out() << "read " << format_number(data.num_sequences()) << " sequence(s) from " << cache
                          << std::endl;
                return;
            }
        }

        auto time = now();
        log_phase_start(true, time, "scanning the FASTA file(s)");
        phase_progress progress(true, cli.input_bytes);
        fasta_scan_t scan = scan_fasta(cli.path_input_files, cli.fasta_alphabet, fasta_separator_symbol, true,
            cli.num_threads, &progress);
        progress.finish();
        log_phase_end(true, time, cli.measurement(), "time_scan_fasta");
        data = std::move(scan.data);
        freq = scan.freq;
        text_size = scan.text_size;
        samples = std::move(scan.samples);
        scanned = true;
        log_out() << "read " << format_number(data.num_sequences()) << " sequence(s) from "
                  << format_number(cli.path_input_files.size()) << " file(s)" << std::endl;
    }

    /**
     * @brief returns the path of the file to build the index from; writes the preprocessed text to a temporary file
     *        if the construction does not read the FASTA file itself
     * @param direct whether the construction reads the FASTA file itself
     */
    std::string construction_input_file(bool direct)
    {
        if (!cli.fasta_mode || direct) return cli.path_input_files[0];

        if (path_text.empty()) {
            path_text = temporary_file_path("move-r_fasta_", cli.path_input_files[0]);
            auto time = now();
            log_phase_start(true, time, "writing the preprocessed FASTA text");
            process_fasta(cli.path_input_files, path_text, cli.fasta_alphabet, fasta_separator_symbol, true,
                cli.num_threads);
            log_phase_end(true, time, cli.measurement(), "time_preprocess_fasta");
        }

        return path_text;
    }

    /**
     * @brief returns the description of the FASTA input for the construction
     * @param direct whether the construction reads the FASTA file itself (else the result is empty)
     */
    fasta_input_t fasta_input(bool direct) const
    {
        if (!cli.fasta_mode || !direct) return fasta_input_t {};
        return fasta_input_t {
            .text_len = text_size,
            .byte_freqs = scanned ? freq.data() : nullptr,
            .separator = (uint8_t) fasta_separator_symbol,
            .alphabet = cli.fasta_alphabet,
            .index_samples = scanned ? samples.data() : nullptr,
            .num_index_samples = samples.size(),
            .sequences = &data };
    }

    // returns the number of input symbols: the length of the FASTA text with -f (known after scan), else the size of
    // the input file
    uint64_t input_length() const { return cli.fasta_mode ? text_size : cli.input_bytes; }

    // deletes the preprocessed text file (if one has been written)
    void remove_preprocessed_copy()
    {
        if (!path_text.empty()) std::filesystem::remove(path_text);
    }
};

#endif // MOVE_R_MISC_BUILD_CLI_HPP
