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
 * @file apm_query_cli.hpp
 * @brief the approximate queries of the modes count and locate of move-rb: their options, their help
 *        and the queries of a patterns file
 */

#ifndef MOVE_R_MISC_APM_QUERY_CLI_HPP
#define MOVE_R_MISC_APM_QUERY_CLI_HPP

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <ips2ra.hpp>

#include <algorithms/move_rb/apm/model.hpp>
#include <misc/alignment_output_cli.hpp>
#include <misc/cli_help.hpp>
#include <misc/fasta.hpp>
#include <misc/files.hpp>
#include <misc/index_cli.hpp>
#include <misc/memory.hpp>
#include <misc/parallel_patterns.hpp>
#include <misc/pattern_reader.hpp>
#include <misc/query_cli.hpp>

// ############################# OPTIONS #############################

// the largest number of errors of a search scheme (-k)
inline constexpr uint64_t max_scheme_errors = std::numeric_limits<decltype(search_scheme_t::k)>::max();

// the largest number of errors of the edit distance (-k with -d edit)
inline constexpr uint64_t max_edit_errors = edit_distance_max_errors_128;

// the search scheme without -scheme and of the approximate queries of the interactive mode of move-rb
inline constexpr std::string_view default_scheme_name = "min-u";

// the error of the edit distance and of the error models for an index without a byte alphabet (more distinct bytes
// than internal symbols of one byte, see move_r::byte_alphabet)
inline const std::string edit_distance_needs_byte_alphabet = "has more distinct bytes than internal symbols of one "
    "byte; only the Hamming distance without an error model is supported";

/** @brief the options of the approximate queries of count and locate */
struct apm_query_cli_t {
    bool locate = true; // the mode locate (else count)
    bool given = false; // true <=> one of the options was given
    uint64_t max_errors = 0; // -k
    distance_metric_t metric = HAMMING_DISTANCE; // -d
    std::string scheme_name = std::string(default_scheme_name); // -scheme: the name of a built-in scheme or a file
    bool bestmap = false; // -bestmap: the number of errors is increased from 0 up to -k until an occurrence is found
    uint64_t bestmap_strata = 0; // -x: with -bestmap, the number of further error counts that are reported
    apm_params_t apm_params; // -partition, -no-seed
    apm_model_t apm_model; // -cost, -max, -exact, -class, -iupac, -qual, -tau, -cost-dp
    std::string path_hist_file; // -hist
    std::string path_sam_file; // -sam
    std::string path_paf_file; // -paf
    std::string path_bed_file; // -bed
    bool write_seq = false; // -seq: SEQ and QUAL are written in the primary SAM record of a read
    bool reverse_complements = true; // -rc: with -sam, -paf or -bed, the reverse complement of each pattern is searched too
    bool best_only = false; // -best: only the primary SAM record of a read is written
    bool no_cigar = false; // -no-cigar: no CIGAR is computed or written
    std::string tool; // the name of the tool, e.g. move-rb (for the program line of a SAM file)
    std::string command_line; // the command line (for the program line of a SAM file)
    // the search schemes of the queries of the patterns file: [0..e] schemes_by_errors[e'] = the scheme for e'
    // errors, with e = max_errors, if -bestmap or an error model is given; else only the scheme for max_errors errors
    std::vector<search_scheme_t> schemes_by_errors;
    direct_ofstream hist_file; // -hist
    direct_ofstream sam_file; // -sam
    direct_ofstream paf_file; // -paf
    direct_ofstream bed_file; // -bed
};

// the files and the names of the approximate queries of a patterns file
struct apm_query_run_t {
    query_options_t& options; // -p, -o, -m
    direct_ofstream& output_file; // -o
    const std::string& path_index_file;
    const std::string& path_patterns_file;
    std::string built_with; // how an index with sequences is built, e.g. "move-rb build in the FASTA mode"
    std::string measurement_name; // the name of the index in the measurement file (algo = <mode>_<measurement_name>)
};

// ############################# HELP #############################

/**
 * @brief prints the option -rc of locate
 * @param h the help printer
 */
inline void print_apm_input_options(help_printer& h)
{
    print_reverse_complements_option(h, true, "with -sam, -paf or -bed: also search the reverse complement of each "
        "pattern; its occurrences are written with the SAM flag 0x10 and the PAF and BED strand -");
}

/**
 * @brief prints the options -k, -d (locate), -scheme, -bestmap and -x
 * @param h the help printer
 * @param locate true <=> the help of locate
 */
inline void print_apm_search_options(help_printer& h, bool locate);

/**
 * @brief prints the group of the options of the error model
 * @param h the help printer
 * @param locate true <=> the help of locate
 */
inline void print_apm_model_options(help_printer& h, bool locate);

/**
 * @brief prints the options -sam, -paf, -bed, -seq, -best and -no-cigar (locate) or -hist (count)
 * @param h the help printer
 * @param locate true <=> the help of locate
 * @param built_with how an index with sequences is built, e.g. "move-rb build in the FASTA mode"
 */
inline void print_apm_output_options(help_printer& h, bool locate, const std::string& built_with);

/**
 * @brief prints the options -partition and -no-seed
 * @param h the help printer
 */
inline void print_apm_execution_options(help_printer& h);

/**
 * @brief reads an option of the approximate queries
 * @param args the arguments, positioned after the option
 * @param option the option
 * @param apm the options of the approximate queries
 * @return whether option is one of them
 * @throws usage_error if the value of the option is invalid
 */
inline bool read_apm_query_option(cli_arguments& args, const std::string& option, apm_query_cli_t& apm);

// ############################# SEARCH SCHEMES #############################

/**
 * @brief builds the search schemes of the queries of the patterns file (see apm_query_cli_t::schemes_by_errors) from
 *        -scheme and -k; a scheme file sets max_errors
 * @param apm the options of the approximate queries
 * @throws usage_error if -scheme names neither a built-in scheme nor a file, or with -bestmap and a scheme file
 */
inline void build_search_schemes(apm_query_cli_t& apm);

// returns the search scheme for at most the given number of errors (see apm_query_cli_t::schemes_by_errors)
inline const search_scheme_t& scheme_for_errors(const apm_query_cli_t& apm, uint64_t errors)
{
    return apm.schemes_by_errors[std::min<uint64_t>(errors, apm.schemes_by_errors.size() - 1)];
}

/**
 * @brief returns the search scheme of a query with at most the given number of errors, or of the error model compiled
 *        for its pattern
 * @param apm the options of the approximate queries
 * @param errors the maximum number of errors
 * @param costs the error model compiled for the pattern (nullptr: none)
 * @throws std::runtime_error if the scheme allows fewer errors than the error model
 */
inline const search_scheme_t& scheme_for_query(const apm_query_cli_t& apm, uint64_t errors, const apm_costs_t* costs);

/**
 * @brief returns the query parameters for a pattern: with an error model, the model compiled for the pattern (it is
 *        valid until the next call on the same thread)
 * @tparam index_t type of the index
 * @param apm the options of the approximate queries
 * @param index the index
 * @param pattern the pattern as it is searched
 * @param qual the qualities of the read the pattern was taken from
 * @param reverse true <=> the pattern is the reverse complement of the read
 * @param bound the maximum cost of an alignment
 * @param metric the distance metric
 */
template <typename index_t>
apm_params_t query_params(const apm_query_cli_t& apm, const index_t& index, const std::string& pattern,
    const std::string& qual, bool reverse, uint64_t bound, distance_metric_t metric);

// returns the details of the patterns line of the query summary, e.g. "search scheme min-u, at most 1 error (edit
// distance)"
inline std::string search_description(const apm_query_cli_t& apm)
{
    return "search scheme " + apm.scheme_name + ", at most " + format_number(apm.max_errors)
        + (apm.max_errors == 1 ? " error" : " errors")
        + (apm.metric == HAMMING_DISTANCE ? " (Hamming distance)" : " (edit distance)");
}

// returns the memory allocated during the queries in bytes (0 if malloc_count is not linked)
inline uint64_t additional_memory(uint64_t baseline)
{
    return is_malloc_count_linked ? malloc_count_peak() - baseline : 0;
}

/**
 * @brief checks the options that only apply to the queries of a patterns file and builds the search schemes
 * @param apm the options of the approximate queries
 * @throws usage_error if the options do not fit together
 */
inline void check_apm_query_options(apm_query_cli_t& apm);

// opens the files of -hist, -sam, -paf and -bed (see open_output_file)
inline void open_apm_query_files(apm_query_cli_t& apm)
{
    open_output_file(apm.hist_file, apm.path_hist_file);
    open_output_file(apm.sam_file, apm.path_sam_file);
    open_output_file(apm.paf_file, apm.path_paf_file);
    open_output_file(apm.bed_file, apm.path_bed_file);
}

// closes the files of -hist, -sam, -paf and -bed (see close_output_file)
inline void close_apm_query_files(apm_query_cli_t& apm)
{
    close_output_file(apm.hist_file, apm.path_hist_file);
    close_output_file(apm.sam_file, apm.path_sam_file);
    close_output_file(apm.paf_file, apm.path_paf_file);
    close_output_file(apm.bed_file, apm.path_bed_file);
}

// returns whether locate computes the CIGARs of the occurrences: with -sam or -paf, without -no-cigar
inline bool apm_locate_with_cigar(const apm_query_cli_t& apm)
{
    return (apm.sam_file.is_open() || apm.paf_file.is_open()) && !apm.no_cigar;
}

/**
 * @brief throws if an index cannot answer the approximate queries of a patterns file: the edit distance and the error
 *        models rank the symbols with the map of the bytes (see move_r::map_int)
 * @tparam index_t type of the index
 * @param apm the options of the approximate queries
 * @param index the index
 * @param path_index_file the index file
 * @throws std::runtime_error if the index has no byte alphabet but the queries need one
 */
template <typename index_t>
void require_apm_alphabet(const apm_query_cli_t& apm, const index_t& index, const std::string& path_index_file)
{
    if (!index.forward_index().byte_alphabet && (apm.metric == EDIT_DISTANCE || !apm.apm_model.is_default()))
        throw std::runtime_error(path_index_file + " " + edit_distance_needs_byte_alphabet);
}

// ############################# COUNT #############################

/**
 * @brief counts the approximate occurrences of the patterns of the patterns file and prints the summary
 * @tparam index_t type of the index (a move_rb that supports approximate pattern matching)
 * @tparam log_sizes_fnc_t type of log_sizes
 * @param index the index
 * @param apm the options of the approximate queries
 * @param run the files and the names of the queries
 * @param log_sizes called as log_sizes(measurement_file) to append the sizes of the data structures to the measurement
 *        line
 */
template <typename index_t, typename log_sizes_fnc_t>
void run_apm_count_queries(const index_t& index, apm_query_cli_t& apm, const apm_query_run_t& run,
    log_sizes_fnc_t log_sizes);

// ############################# LOCATE #############################

/**
 * @brief locates the occurrences of a read on one strand (with -sam, -paf or -bed), sorts them by position and removes the
 *        edit-distance occurrences that overlap a better one (see filter_edit_distance_occurrences); a read with at
 *        most p symbols (p: the number of parts of the search scheme) has no occurrences
 * @tparam cm whether the occurrences carry a CIGAR
 * @tparam index_t type of the index
 * @param apm the options of the approximate queries
 * @param index the index
 * @param read the read as it is searched
 * @param qual the qualities of the read
 * @param reverse whether read is the reverse complement of the read
 * @param strand is set to the occurrences and their CIGARs
 */
template <cigar_mode_t cm, typename index_t>
void locate_strand(const apm_query_cli_t& apm, const index_t& index, const std::string& read, const std::string& qual,
    bool reverse, strand_occurrences_t<cm>& strand);

/**
 * @brief locates the approximate occurrences of the patterns of the patterns file and prints the summary
 * @tparam cm whether the occurrences carry a CIGAR (with -sam or -paf without -no-cigar)
 * @tparam index_t type of the index (a move_rb that supports approximate pattern matching and locates all occurrences)
 * @tparam log_sizes_fnc_t type of log_sizes
 * @param index the index
 * @param apm the options of the approximate queries
 * @param run the files and the names of the queries
 * @param log_sizes called as log_sizes(measurement_file) to append the sizes of the data structures to the measurement
 *        line
 */
template <cigar_mode_t cm, typename index_t, typename log_sizes_fnc_t>
void run_apm_locate_queries(const index_t& index, apm_query_cli_t& apm, const apm_query_run_t& run,
    log_sizes_fnc_t log_sizes);

#include <misc/apm_query_cli.tpp>

#endif // MOVE_R_MISC_APM_QUERY_CLI_HPP
