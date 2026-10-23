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
 * @file extract.cpp
 * @brief the modes extract-sa, extract-bwt and extract-rlbwt of move-rb
 */

#include <stdexcept>
#include <string>
#include <type_traits>

#include <misc/extract_modes_cli.hpp>
#include <misc/move_rb_cli.hpp>
#include <misc/rlbwt_writer.hpp>
#include <move_rb/move_rb.hpp>

#include "modes.hpp"

namespace {

/**
 * @brief returns the description of an extract mode of move-rb
 * @param name the name of the mode
 * @param description what the mode writes
 * @param output the format of the output file and the indexes that provide it
 * @param reversed_explanation the explanation of -rev (empty: the mode has no option -rev)
 */
extract_mode_t move_rb_extract_mode(const std::string& name, const std::string& description, const std::string& output,
    const std::string& reversed_explanation)
{
    return { .tool = "move-rb", .name = name, .description = description, .output = output,
        .index_explanation = "an index built with move-rb build",
        .reversed_explanation = reversed_explanation };
}

// the explanation of -rev in the extract modes of move-rb
const std::string reversed_text_explanation = "write the output for the reversed text";

/**
 * @brief runs an extract mode of move-rb (see run_extract_mode and extract_from_index)
 * @tparam missing_fnc_t type of missing
 * @tparam write_fnc_t type of write
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments
 * @param mode the mode
 * @param missing returns the error message if an index type does not provide the output
 * @param write write(index, cli) writes the output for a move_rb index
 * @return the exit code
 */
template <typename missing_fnc_t, typename write_fnc_t>
int run_move_rb_extract(int argc, char** argv, const extract_mode_t& mode, missing_fnc_t missing, write_fnc_t write)
{
    return run_extract_mode(argc, argv, mode, [&](direct_ifstream& index_file, const extract_cli_t& cli) {
        const bool known = dispatch_move_rb_index(read_move_rb_support(index_file), [&]<move_r_support support>() {
            extract_from_index<move_rb<support>>(index_file, cli, "bidirectional index", missing, write);
        });

        if (!known) throw std::runtime_error(cli.path_index_file + " is not an index built with move-rb build");
    });
}

} // namespace

int move_rb_extract_sa(int argc, char** argv)
{
    return run_move_rb_extract(argc, argv, move_rb_extract_mode("extract-sa", move_rb_extract_sa_description,
        "The output file holds the suffix array of the text with the terminator: " + std::string(positions_format)
            + "; SA[0] = input size. It requires an index built with -s locate-move or locate-rlzsa. "
            "The index of the reversed text stores only the samples at the BWT run boundaries (see extract-rlbwt).",
        ""),
        []<typename index_t>(const extract_cli_t&) {
            using fwd_t = std::remove_cvref_t<decltype(std::declval<index_t>().forward_index())>;
            return std::string(fwd_t::supports_multiple_locate ? ""
                : "the index does not store the suffix array (see move-rb build -s)");
        },
        [](const auto& index, const extract_cli_t& cli) {
            using fwd_t = std::remove_cvref_t<decltype(index.forward_index())>;
            if constexpr (fwd_t::supports_multiple_locate) write_suffix_array(index.forward_index(), cli);
        });
}

int move_rb_extract_bwt(int argc, char** argv)
{
    return run_move_rb_extract(argc, argv, move_rb_extract_mode("extract-bwt", move_rb_extract_bwt_description,
        "The output file holds the BWT of the text (or of the reversed text, -rev) with the terminator: "
            + std::string(bwt_format) + ".", reversed_text_explanation),
        every_index_provides, [](const auto& index, const extract_cli_t& cli) {
            with_selected_index(index, cli, [&](const auto& selected) { write_bwt(selected, cli); });
        });
}

int move_rb_extract_rlbwt(int argc, char** argv)
{
    return run_move_rb_extract(argc, argv, move_rb_extract_mode("extract-rlbwt", move_rb_extract_rlbwt_description,
        std::string(rlbwt_file_format) + " Present fields: sa-beg and sa-end unless the index was built with -s "
            "count.", reversed_text_explanation),
        every_index_provides, [](const auto& index, const extract_cli_t& cli) {
            with_selected_index(index, cli, [&](const auto& selected) {
                write_rlbwt(selected, cli.path_output_file, cli.num_threads, true);
            });
        });
}
