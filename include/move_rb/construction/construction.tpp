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

#ifndef MOVE_R_MOVE_RB_CONSTRUCTION_CONSTRUCTION_TPP
#define MOVE_R_MOVE_RB_CONSTRUCTION_CONSTRUCTION_TPP

#include <move_rb/move_rb.hpp>

template <move_r_support support>
template <typename inp_t_>
move_rb<support>::move_rb(inp_t_& input, move_r_params params, bool log_summary)
{
    using input_t = std::remove_cvref_t<inp_t_>;
    using input_symbol_t = typename input_t::value_type;
    static_assert(std::is_integral_v<input_symbol_t>,
        "move_rb input symbols have to be integers");
    static_assert(sizeof(input_symbol_t) == 1 || sizeof(input_symbol_t) == 2
        || sizeof(input_symbol_t) == 4 || sizeof(input_symbol_t) == 8,
        "move_rb supports symbols of one, two, four or eight bytes");

    constexpr bool byte_input = sizeof(input_symbol_t) == 1;
    using view_sym_t = std::conditional_t<byte_input, uint8_t, input_symbol_t>;

    if (!params.input_file.empty()) {
        if constexpr (!byte_input) {
            throw std::invalid_argument(
                "move_rb: an alphabet wider than a byte has to be constructed from memory");
        }
        // for a FASTA file, the length of the text was computed by scan_fasta
        n = params.fasta.text_len != 0 ? params.fasta.text_len + 1
                                        : std::filesystem::file_size(params.input_file) + 1;

        no_init_vector<view_sym_t> file_buffer;
        move_r_input_view view(file_buffer);
        build(view, std::move(params), log_summary);
        return;
    }

    n = input.size() + 1;
    move_r_input_view view(input);
    build(view, std::move(params), log_summary);
}

#endif // MOVE_R_MOVE_RB_CONSTRUCTION_CONSTRUCTION_TPP
