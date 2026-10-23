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

// explicit instantiations of the indexes. This file is compiled once per part: the build defines
// MOVE_R_INSTANTIATION_PART_<part> for each translation unit (see CMakeLists.txt), and each part instantiates
// one group of members for all index types.

#define MOVE_R_DECLARATIONS_ONLY

#if defined(MOVE_R_INSTANTIATION_PART_MOVE_R)

// explicit instantiation of the move_r classes; the construction phases are instantiated in
// move_r_construction_*.cpp (see move_r/construction/instantiations.hpp)

#include <move_r/construction/instantiations.hpp>

#define MOVE_R_DEFINE_INDEX(supp, enc) template class MOVE_R_INDEX_OF(supp, enc);

MOVE_R_CLASSES(MOVE_R_DEFINE_INDEX)

#elif defined(MOVE_R_INSTANTIATION_PART_MOVE_R_CONSTRUCTION_TEXT)

// explicit instantiation of the construction phases of move_r up to the RLBWT: the execution of the phases,
// the preprocessing of the input, the prefix-free parsing and the RLBWT (see move_r/construction/instantiations.hpp)

#include <move_r/construction/instantiations.hpp>
#include <move_r/construction/build_index.tpp>
#include <move_r/construction/preprocess.tpp>
#include <move_r/construction/pfp.tpp>
#include <move_r/construction/rlbwt.tpp>

#define MOVE_R_DEFINE_MEMBERS(supp, enc) \
    MOVE_R_CONSTRUCTION_DRIVER_MEMBERS(, supp, enc) \
    MOVE_R_CONSTRUCTION_TEXT_MEMBERS(, supp, enc)

MOVE_R_CLASSES(MOVE_R_DEFINE_MEMBERS)

#elif defined(MOVE_R_INSTANTIATION_PART_MOVE_R_CONSTRUCTION_MOVE)

// explicit instantiation of the construction phases of move_r that build the move data structures and the
// suffix array samples (see move_r/construction/instantiations.hpp)

#include <move_r/construction/instantiations.hpp>
#include <move_r/construction/move_lf.tpp>
#include <move_r/construction/move_phi.tpp>
#include <move_r/construction/load_store.tpp>

#define MOVE_R_DEFINE_MEMBERS(supp, enc) MOVE_R_CONSTRUCTION_MOVE_MEMBERS(, supp, enc)

MOVE_R_CLASSES(MOVE_R_DEFINE_MEMBERS)

#elif defined(MOVE_R_INSTANTIATION_PART_MOVE_R_CONSTRUCTION_SA_RLZ)

// explicit instantiation of the construction phase of move_r that builds the rlz of the differential suffix
// array (see move_r/construction/instantiations.hpp)

#include <move_r/construction/instantiations.hpp>
#include <rlz_ref/construction/construction.hpp>
#include <move_r/construction/sa_rlz.tpp>

#define MOVE_R_DEFINE_MEMBERS(supp, enc) MOVE_R_CONSTRUCTION_SA_RLZ_MEMBERS(, supp, enc)

MOVE_R_SA_RLZ_CLASSES(MOVE_R_DEFINE_MEMBERS)

#elif defined(MOVE_R_INSTANTIATION_PART_MOVE_RB_1)

// explicit instantiation of the first half of the move_rb classes, without the construction
// (see move_rb/instantiations.hpp)

#include <move_rb/instantiations.hpp>

#define MOVE_R_DEFINE_MOVE_RB(...) template class move_rb<__VA_ARGS__>;

MOVE_R_MOVE_RB_CLASSES_1(MOVE_R_DEFINE_MOVE_RB)

#elif defined(MOVE_R_INSTANTIATION_PART_MOVE_RB_2)

// explicit instantiation of the second half of the move_rb classes, without the construction
// (see move_rb/instantiations.hpp)

#include <move_rb/instantiations.hpp>

#define MOVE_R_DEFINE_MOVE_RB(...) template class move_rb<__VA_ARGS__>;

MOVE_R_MOVE_RB_CLASSES_2(MOVE_R_DEFINE_MOVE_RB)

#elif defined(MOVE_R_INSTANTIATION_PART_MOVE_RB_BUILD_1)

// explicit instantiation of build() and build_context_tables() of the first half of the move_rb classes
// (see move_rb/instantiations.hpp)

#include <move_rb/instantiations.hpp>
#include <move_rb/construction/build.tpp>
#include <move_rb/construction/context_tables.tpp>

#define MOVE_R_DEFINE_BUILD(...) MOVE_R_MOVE_RB_BUILD_MEMBERS(, __VA_ARGS__)

MOVE_R_MOVE_RB_CLASSES_1(MOVE_R_DEFINE_BUILD)

#elif defined(MOVE_R_INSTANTIATION_PART_MOVE_RB_BUILD_2)

// explicit instantiation of build() and build_context_tables() of the second half of the move_rb classes
// (see move_rb/instantiations.hpp)

#include <move_rb/instantiations.hpp>
#include <move_rb/construction/build.tpp>
#include <move_rb/construction/context_tables.tpp>

#define MOVE_R_DEFINE_BUILD(...) MOVE_R_MOVE_RB_BUILD_MEMBERS(, __VA_ARGS__)

MOVE_R_MOVE_RB_CLASSES_2(MOVE_R_DEFINE_BUILD)

#elif defined(MOVE_R_INSTANTIATION_PART_MOVE_RB_APM_1)

// explicit instantiation of the approximate pattern matching methods of the first half of the move_rb classes
// that support them (see move_rb/instantiations.hpp)

#include <move_rb/instantiations.hpp>

#define MOVE_R_DEFINE_APM(...) MOVE_R_MOVE_RB_APM_MEMBERS(, __VA_ARGS__)

MOVE_R_MOVE_RB_APM_CLASSES_1(MOVE_R_DEFINE_APM)

#elif defined(MOVE_R_INSTANTIATION_PART_MOVE_RB_APM_2)

// explicit instantiation of the approximate pattern matching methods of the second half of the move_rb classes
// that support them (see move_rb/instantiations.hpp)

#include <move_rb/instantiations.hpp>

#define MOVE_R_DEFINE_APM(...) MOVE_R_MOVE_RB_APM_MEMBERS(, __VA_ARGS__)

MOVE_R_MOVE_RB_APM_CLASSES_2(MOVE_R_DEFINE_APM)

#elif defined(MOVE_R_INSTANTIATION_PART_SORT)

// explicit instantiation of the sorting functions (see misc/sort.hpp) for the types of move_r

#include <move_r/construction/instantiations.hpp>

MOVE_R_ORDER_ALPHABET_SORTS()
MOVE_R_CONSTRUCTION_SORTS()

#elif defined(MOVE_R_INSTANTIATION_PART_TEXT_OPS)

// explicit instantiation of the text operations and of the alphabet map for each symbol width of the input
// (see move_r/construction/text_ops.hpp)

#include <move_r/move_r.hpp>

#define MOVE_R_DEFINE_TEXT_OPS(sym_t) \
    template class move_r_int_alphabet_of<sym_t>; \
    template void move_r_unmap_text_of<sym_t>( \
        move_r_input_view&, uint64_t, uint64_t, const elias_fano_bit_vector&, \
        const interleaved_bit_aligned_vectors&); \
    template void move_r_unmap_text_by_symbols_of<sym_t>(move_r_input_view&, uint64_t, \
        uint64_t, const std::vector<uint64_t>&); \
    template void move_r_reverse_text_of<sym_t>(move_r_input_view&, uint64_t, uint64_t); \
    template void move_r_build_sa_of<sym_t>( \
        move_r_input_view&, uint40_t*, uint64_t, uint64_t, uint64_t, uint64_t); \
    template void move_r_bwt_runs_from_sa_of<sym_t>(move_r_input_view&, \
        const no_init_vector<uint40_t>&, const std::vector<uint64_t>&, std::vector<uint64_t>&, \
        std::vector<interleaved_bit_aligned_vectors>&, std::vector<std::vector<uint64_t>>&, \
        bool, uint64_t, phase_progress&);

MOVE_R_FOR_EACH_INPUT_WIDTH(MOVE_R_DEFINE_TEXT_OPS)

#else
#error "no part of the instantiations selected (MOVE_R_INSTANTIATION_PART_<part>)"
#endif
