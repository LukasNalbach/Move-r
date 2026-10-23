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

#ifndef MOVE_R_MOVE_R_CONSTRUCTION_INSTANTIATIONS_HPP
#define MOVE_R_MOVE_R_CONSTRUCTION_INSTANTIATIONS_HPP

/**
 * @file instantiations.hpp
 * @brief extern template declarations of the index types, which are compiled in the instantiation library
 *        (src/instantiations.cpp); define STRING_UTILS_NO_EXTERN_TEMPLATES if the library is not linked
 */

#include <move_r/move_r.hpp>

// the index type (support, mlf_enc) and its construction
#define MOVE_R_INDEX_OF(supp, enc) move_r<supp, enc>
#define MOVE_R_CONSTRUCTION_OF(supp, enc) MOVE_R_INDEX_OF(supp, enc)::construction

// all valid index types (support, mlf_enc), in two lists
#define MOVE_R_CLASSES_1(X) \
    X(_count, POS) \
    X(_locate_one, POS) \
    X(_locate_move, POS) \
    X(_locate_rlzsa, POS) \
    X(_count_bi, POS) \
    X(_locate_move_bi, POS)
#define MOVE_R_CLASSES_2(X) \
    X(_count, EF) \
    X(_locate_one, EF) \
    X(_locate_move, EF) \
    X(_locate_rlzsa, EF) \
    X(_locate_rlzsa_bi, POS) \
    X(_locate_one_bi, POS)
#define MOVE_R_CLASSES(X) \
    MOVE_R_CLASSES_1(X) \
    MOVE_R_CLASSES_2(X)

// the index types with an rlz of the differential suffix array (see has_SA_rlz)
#define MOVE_R_SA_RLZ_CLASSES(X) \
    X(_locate_rlzsa, POS) \
    X(_locate_rlzsa, EF) \
    X(_locate_rlzsa_bi, POS)

/* The phases of the construction (without member templates), grouped by the files defining them. Each phase
   is instantiated explicitly by one unit of the instantiation library, so a phase that is added to one of the
   files has to be added to the corresponding list. */

// the constructor, build_index and the log methods (build_index.tpp)
#define MOVE_R_CONSTRUCTION_DRIVER_MEMBERS(prefix, supp, enc) \
    prefix template MOVE_R_CONSTRUCTION_OF(supp, enc)::construction(MOVE_R_INDEX_OF(supp, enc)&, \
        move_r_input_view&, bool, move_r_params); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::choose_sa_samples_layout(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::start_construction(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::limit_num_threads(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::log_phase_end(const std::string&); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::log_phase_end_quiet(const std::string&); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::log_mds_phase_end(const std::string&, const std::string&, uint64_t); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::log_peak_mem_usage(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::build_index(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::log_finished(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::log_statistics();

// the input, the prefix-free parsing and the run-length encoded BWT (preprocess.tpp, pfp.tpp, rlbwt.tpp)
#define MOVE_R_CONSTRUCTION_TEXT_MEMBERS(prefix, supp, enc) \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::preprocess_t(const std::string&); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::read_t_from_file(const std::string&); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::unmap_t(); \
    prefix template std::array<uint8_t, 256> MOVE_R_CONSTRUCTION_OF(supp, enc)::pfp_alphabet_map() const; \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::release_pfp_output(); \
    prefix template bool MOVE_R_CONSTRUCTION_OF(supp, enc)::read_pfp_cache(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::write_pfp_cache(const std::string&); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::build_rlbwt_pfp(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::read_lcp_run_pfp(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::read_iphi_inv_pfp(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::build_sa(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::build_rlbwt_c(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::compute_C(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::store_rlbwt(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::load_rlbwt();

// the move data structures and the samples (move_lf.tpp, move_phi.tpp, load_store.tpp)
#define MOVE_R_CONSTRUCTION_MOVE_MEMBERS(prefix, supp, enc) \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::build_mlf(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::build_SA_run_from_iphi_inv(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::build_data_structures_over_l_(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::build_ilf(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::choose_revert_starts(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::build_sa_run_from_sa(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::build_iphi_inv_from_sa(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::build_mphi_inv(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::build_mphi(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::build_SA_Phi_inv(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::store_iphi_inv(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::load_iphi_inv(bool); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::build_iphi(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::sort_iphi_inv(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::sort_iphi(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::build_plcp_run(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::store_run_ranks(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::store_mlf(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::load_mlf(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::store_data_structures_over_l_(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::load_data_structures_over_l_();

// the rlz of the differential suffix array, for the index types in MOVE_R_SA_RLZ_CLASSES (sa_rlz.tpp)
#define MOVE_R_CONSTRUCTION_SA_RLZ_MEMBERS(prefix, supp, enc) \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::SA_rlz_reader::step_right(); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::SA_rlz_reader::seek(uint64_t); \
    prefix template uint64_t MOVE_R_CONSTRUCTION_OF(supp, enc)::SA_rlz_reader::SAd(uint64_t); \
    prefix template void MOVE_R_CONSTRUCTION_OF(supp, enc)::build_SA_rlz();

// the radix sorts used by the constructions of move_r and move_rb (see misc/sort.hpp; those of the move
// data structure are listed in move_data_structure/construction/construction.hpp)
using move_r_u64_pairs_it = std::vector<std::pair<uint64_t, uint64_t>>::iterator;
using move_r_u64_u32_pairs_it = std::vector<std::pair<uint64_t, uint32_t>>::iterator;

#define MOVE_R_CONSTRUCTION_SORTS(prefix) \
    MOVE_DATA_STRUCTURE_RADIX_SORT(prefix, move_r_u64_pairs_it, key_of_first) \
    MOVE_DATA_STRUCTURE_RADIX_SORT(prefix, move_r_u64_u32_pairs_it, key_of_first) \
    MOVE_DATA_STRUCTURE_RADIX_SORT(prefix, std::vector<uint64_t>::iterator, key_of_value)

// the editor parser of VS Code (__INTELLISENSE__) does not resolve the extern declarations of constrained member
// templates; it does not need them
#if !defined(STRING_UTILS_NO_EXTERN_TEMPLATES) && !defined(__INTELLISENSE__)

MOVE_R_CONSTRUCTION_SORTS(extern)

// declares the index type with all its members that are not templates (including its construction)
#define MOVE_R_DECLARE_INDEX(supp, enc) extern template class MOVE_R_INDEX_OF(supp, enc);

MOVE_R_CLASSES(MOVE_R_DECLARE_INDEX)

#undef MOVE_R_DECLARE_INDEX

#endif

#endif // MOVE_R_MOVE_R_CONSTRUCTION_INSTANTIATIONS_HPP
