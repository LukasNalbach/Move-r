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
 * @file modes.hpp
 * @brief the modes of the command-line tool move-rb and their descriptions (in the list of modes and at the start of
 *        the help of each mode); related modes share a translation unit, and a mode is called with argv[0] = its name
 */

#ifndef MOVE_R_CLI_MOVE_RB_MODES_HPP
#define MOVE_R_CLI_MOVE_RB_MODES_HPP

inline constexpr const char* move_rb_description = "builds and queries bidirectional move-r indexes";
inline constexpr const char* move_rb_build_description = "builds a bidirectional index of a text";
inline constexpr const char* move_rb_count_description
    = "counts the occurrences of patterns with at most k mismatches (Hamming distance)";
inline constexpr const char* move_rb_locate_description
    = "locates the occurrences of patterns with at most k errors (Hamming or edit distance)";
inline constexpr const char* move_rb_revert_description = "reconstructs the text from an index";
inline constexpr const char* move_rb_extract_sa_description = "writes the suffix array of the text";
inline constexpr const char* move_rb_extract_bwt_description = "writes the BWT of the text or of the reversed text";
inline constexpr const char* move_rb_extract_rlbwt_description
    = "writes the BWT runs of the text or of the reversed text and the values the index stores for them";
inline constexpr const char* move_rb_info_description
    = "prints the support type, the parameters, the sizes of the data structures and the modes of an index";

int move_rb_build(int argc, char** argv);
int move_rb_count(int argc, char** argv);
int move_rb_locate(int argc, char** argv);
int move_rb_revert(int argc, char** argv);
int move_rb_extract_sa(int argc, char** argv);
int move_rb_extract_bwt(int argc, char** argv);
int move_rb_extract_rlbwt(int argc, char** argv);
int move_rb_info(int argc, char** argv);

#endif // MOVE_R_CLI_MOVE_RB_MODES_HPP
