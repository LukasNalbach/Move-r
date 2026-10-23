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
 * @file pfp_integer.hpp
 * @brief computes the RLBWT of a text with an integer alphabet with the prefix-free parsing (pfp-rlbwt)
 */

#ifndef MOVE_R_MOVE_R_CONSTRUCTION_PFP_INTEGER_HPP
#define MOVE_R_MOVE_R_CONSTRUCTION_PFP_INTEGER_HPP

#include <cstdint>
#include <string>

#include <move_r/construction/input_view.hpp>
#include <pfp_rlbwt/pfp_rlbwt.hpp>

/**
 * @brief a text whose symbols are read as their effective symbols: from memory, where they have been replaced by
 *        their effective symbols, or from a file, whose symbols are replaced while they are read
 */
struct pfp_integer_text_t {
    const char* text = nullptr; // [0..length-1] the text in memory (nullptr <=> the text is read from file)
    std::string file; // the file storing the text (if text is nullptr)
    uint8_t width = 0; // bytes per symbol
    uint64_t length = 0; // number of symbols
    // maps the symbols of the file to their effective symbols (nullptr <=> they are effective symbols)
    const move_r_int_alphabet* alphabet = nullptr;
    uint64_t max_symbol = 0; // largest effective symbol
};

/**
 * @brief writes the output of pfp_rlbwt::build for the text, with each effective symbol increased by offset
 * @param input the text
 * @param offset added to each effective symbol, so that all symbols are at least pfp_rlbwt::min_symbol (see
 *        construction::min_mapped_symbol)
 * @param reverse whether the output is computed for the reversed text
 * @param prefix prefix of the output files
 * @param prm parameters of pfp_rlbwt::build
 */
void build_rlbwt_pfp_integer(const pfp_integer_text_t& input, uint64_t offset, bool reverse,
    const std::string& prefix, const pfp_rlbwt::params& prm);

#endif // MOVE_R_MOVE_R_CONSTRUCTION_PFP_INTEGER_HPP
