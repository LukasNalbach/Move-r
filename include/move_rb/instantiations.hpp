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
 * @file instantiations.hpp
 * @brief lists the instantiations of move_rb that are compiled by the instantiation library (src/instantiations.cpp)
 */

#ifndef MOVE_R_MOVE_RB_INSTANTIATIONS_HPP
#define MOVE_R_MOVE_RB_INSTANTIATIONS_HPP

#include <move_rb/move_rb.hpp>

// all types of move_rb, in two halves (see MOVE_R_CLASSES_1)
#define MOVE_R_MOVE_RB_CLASSES_1(X) \
    X(_count) \
    X(_locate_one)
#define MOVE_R_MOVE_RB_CLASSES_2(X) \
    X(_locate_move) \
    X(_locate_rlzsa)
#define MOVE_R_MOVE_RB_CLASSES(X) \
    MOVE_R_MOVE_RB_CLASSES_1(X) \
    MOVE_R_MOVE_RB_CLASSES_2(X)

// a step of build (see move_rb::build_state_t)
#define MOVE_R_MOVE_RB_BUILD_STEP(prefix, step, ...) \
    prefix template void move_rb<__VA_ARGS__>::step(move_rb<__VA_ARGS__>::build_state_t&);

// the non-template members of the construction (build.tpp and construction/context_tables.tpp)
#define MOVE_R_MOVE_RB_BUILD_MEMBERS(prefix, ...) \
    prefix template void move_rb<__VA_ARGS__>::build(move_r_input_view&, move_r_params, bool); \
    prefix template void move_rb<__VA_ARGS__>::reverse(move_r_input_view&, uint64_t, bool, std::ostream*, \
        const std::string&); \
    MOVE_R_MOVE_RB_BUILD_STEP(prefix, prepare_input, __VA_ARGS__) \
    MOVE_R_MOVE_RB_BUILD_STEP(prefix, build_byte_alphabet, __VA_ARGS__) \
    MOVE_R_MOVE_RB_BUILD_STEP(prefix, build_integer_alphabet, __VA_ARGS__) \
    MOVE_R_MOVE_RB_BUILD_STEP(prefix, build_sub_indexes, __VA_ARGS__) \
    MOVE_R_MOVE_RB_BUILD_STEP(prefix, load_forward_index, __VA_ARGS__) \
    MOVE_R_MOVE_RB_BUILD_STEP(prefix, build_SA_pos, __VA_ARGS__) \
    MOVE_R_MOVE_RB_BUILD_STEP(prefix, restore_input, __VA_ARGS__) \
    MOVE_R_MOVE_RB_BUILD_STEP(prefix, build_context_tables_phase, __VA_ARGS__) \
    prefix template void move_rb<__VA_ARGS__>::build_context_tables(double, uint64_t, bool, uint64_t); \
    prefix template void move_rb<__VA_ARGS__>::build_level_samples(::context_tables::level_t&, \
        const std::vector<move_rb<__VA_ARGS__>::sample_info_pack_t>&); \
    prefix template void move_rb<__VA_ARGS__>::build_reverse_order(::context_tables::level_t&, \
        std::vector<uint64_t>&, std::vector<uint64_t>&, const std::vector<uint64_t>&, const std::vector<uint64_t>&, \
        const no_init_vector<uint32_t>&, uint64_t, const std::vector<uint64_t>&, no_init_vector<uint32_t>&, \
        std::vector<uint64_t>&, uint64_t) const;

// the types that locate all occurrences with approximate pattern matching (see locate_apm), in two halves
#define MOVE_R_MOVE_RB_APM_CLASSES_1(X) \
    X(_locate_move)
#define MOVE_R_MOVE_RB_APM_CLASSES_2(X) \
    X(_locate_rlzsa)
#define MOVE_R_MOVE_RB_APM_CLASSES(X) \
    MOVE_R_MOVE_RB_APM_CLASSES_1(X) \
    MOVE_R_MOVE_RB_APM_CLASSES_2(X)

// the member templates for approximate pattern matching (locate_apm for each cigar mode)
#define MOVE_R_MOVE_RB_APM_MEMBERS(prefix, ...) \
    prefix template std::vector<cigar_t<uint8_t>> move_rb<__VA_ARGS__>::locate_apm<NO_CIGAR>( \
        const std::vector<uint8_t>&, const search_scheme_t&, distance_metric_t, \
        std::vector<aprx_occ_t<NO_CIGAR>>&, apm_params_t) const; \
    prefix template std::vector<cigar_t<uint8_t>> move_rb<__VA_ARGS__>::locate_apm<CIGAR>( \
        const std::vector<uint8_t>&, const search_scheme_t&, distance_metric_t, \
        std::vector<aprx_occ_t<CIGAR>>&, apm_params_t) const; \
    prefix template uint64_t move_rb<__VA_ARGS__>::count_hamming_dist<uint8_t>( \
        move_r_pattern_view<uint8_t>, const search_scheme_t&, apm_params_t) const; \
    prefix template std::vector<uint64_t> \
        move_rb<__VA_ARGS__>::count_hamming_dist_histogram<uint8_t>( \
        move_r_pattern_view<uint8_t>, const search_scheme_t&, apm_params_t) const;

// the editor parser of VS Code (__INTELLISENSE__) does not resolve the extern declarations of constrained member
// templates; it does not need them
#if !defined(STRING_UTILS_NO_EXTERN_TEMPLATES) && !defined(__INTELLISENSE__)

// extern template declarations of the class (includes all non-template members) and of the member templates
#define MOVE_R_DECLARE_MOVE_RB(...) extern template class move_rb<__VA_ARGS__>;
#define MOVE_R_DECLARE_MOVE_RB_APM(...) MOVE_R_MOVE_RB_APM_MEMBERS(extern, __VA_ARGS__)

MOVE_R_MOVE_RB_CLASSES(MOVE_R_DECLARE_MOVE_RB)
MOVE_R_MOVE_RB_APM_CLASSES(MOVE_R_DECLARE_MOVE_RB_APM)

#undef MOVE_R_DECLARE_MOVE_RB_APM
#undef MOVE_R_DECLARE_MOVE_RB

#endif

#endif // MOVE_R_MOVE_RB_INSTANTIATIONS_HPP
