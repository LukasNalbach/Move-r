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

#ifndef MOVE_R_MOVE_R_MOVE_R_HPP
#define MOVE_R_MOVE_R_MOVE_R_HPP

#include <iostream>
#include <string_view>
#include <data_structures/elias_fano_bit_vector.hpp>
#include <data_structures/interleaved_bit_aligned_vectors.hpp>
#include <data_structures/rmq_tree.hpp>
#include <move_data_structure/move_data_structure.hpp>
#include <data_structures/pred_succ_block.hpp>
#include <data_structures/pred_succ_occ.hpp>
#include <data_structures/run_sample_positions.hpp>
#include <data_structures/wavelet_tree_prefix_sum.hpp>
#include <rlz_ref/rlz.hpp>
#include <misc/utils.hpp>
#include <misc/portable_layout.hpp>
#include <misc/files.hpp>
#include <misc/strands.hpp>
#include <misc/value_writer.hpp>
#include <misc/search.hpp>
#include <data_structures/fasta_sequence_data.hpp>
#include <pfp_rlbwt/fasta.hpp>
#include <move_r/construction/input_view.hpp>
#include <move_r/construction/text_ops.hpp>
#include <misc/hash.hpp>
#include <misc/index_facts.hpp>
#include <misc/log.hpp>
#include <omp.h>
#include <tsl/sparse_map.h>
#include <concepts>
#include <functional>
#include <type_traits>
#include <variant>

// external symbol of the terminator: the bit pattern with only the highest bit set, which must not occur in the
// input (see move_r::map_symbol and move_r::unmap_symbol)
inline constexpr uint64_t TERMINATOR_SYMBOL = uint64_t(1) << (std::numeric_limits<uint64_t>::digits - 1);

// returns sym, or 0 if sym is TERMINATOR_SYMBOL (the BWT is exported with the terminator as 0, see move_r::BWT_range)
inline constexpr uint64_t terminator_as_zero(uint64_t sym) { return sym == TERMINATOR_SYMBOL ? 0 : sym; }

/**
 * @brief type of locate support
 */
enum move_r_support {
    // the values are stored in index files (see serialize) and must not be changed
    _count = 0, // only count support (no locate support)
    _locate_one = 1, // computes exactly one occurrence per pattern
    _locate_move = 2, // locate support with a move data structure for Phi^{-1} (M_Phi^{-1})
    _locate_rlzsa = 3, // locate support with a relative Lempel-Ziv encoding of the differential suffix array

    // ##################### SUB-INDEXES OF A BIDIRECTIONAL INDEX (DON'T USE DIRECTLY) #####################
    _count_bi = 4, // only count support
    _locate_move_bi = 5, // stores M_Phi and M_Phi^{-1}; M_Phi^{-1} stores the PLCP values as row 0
    _locate_rlzsa_bi = 6, // stores an rlz of the differential suffix array
    // stores only the suffix array samples at the BWT run boundaries (one occurrence per pattern)
    _locate_one_bi = 7
};

// the index files store the support in sizeof(move_r_support) bytes (the underlying type is int or unsigned int)
static_assert(sizeof(move_r_support) == sizeof(uint32_t));

/**
 * @brief returns the name of a locate-support type (used as the value of the option -s and in measurement files)
 * @param support a locate-support type
 * @return the name
 */
inline std::string_view move_r_support_name(move_r_support support)
{
    switch (support) {
        case _count:        return "count";
        case _locate_one:   return "locate_one";
        case _locate_rlzsa: return "locate_rlzsa";
        default:            return "locate_move";
    }
}

/**
 * @brief returns the file-name suffix associated with a locate-support type
 * @param support a locate-support type
 * @return the file-name suffix (empty for the sub-indexes of a bidirectional index)
 */
inline std::string_view move_r_support_file_suffix(move_r_support support)
{
    switch (support) {
        case _count:        return "-count";
        case _locate_one:   return "-locate-one";
        case _locate_move:  return "-locate";
        case _locate_rlzsa: return "-rlzsa";
        default:            return "";
    }
}

/**
 * @brief returns the name of the file a move-r index is stored in: <base>.move-r-<support>
 * @param base base name (usually the name of the input file)
 * @param support a locate-support type
 * @return the file name
 */
inline std::string move_r_index_file_name(const std::string& base, move_r_support support)
{
    return base + ".move-r" + std::string(move_r_support_file_suffix(support));
}

/**
 * @brief move-r construction mode
 */
enum move_r_construction_mode {
    // builds the RLBWT with the prefix-free parsing (pfp-rlbwt) and stores data structures on disk while they are not
    // needed
    _pfp,
    _suffix_array, // builds the suffix array in memory and stores no data structures on disk
    _suffix_array_space // builds the suffix array and stores some data structures on disk
};

// ############################# PARAMETERS OF THE DATA STRUCTURES OVER L' #############################

// maximum alphabet size for which the smaller symbols in a range of L' are counted with PSB_L', the byte alphabet; for
// larger alphabets, a wavelet tree over L' is built (see move_r::WTPS_L_)
inline constexpr uint64_t max_sigma_num_smaller_psb = uint64_t(1) << (8 * sizeof(uint8_t));

// minimum width (in runs) of a range of L' for which L_prev_min is queried (see L_prev_min_max_symbols)
inline constexpr uint64_t min_width_L_prev_min = 96;

// arity of L_prev_min
inline constexpr uint8_t L_prev_min_arity = 8;

/**
 * @brief returns the maximum width of a range of L' for which the extensions of a pattern are enumerated by
 *        scanning and sorting the range; for wider ranges, they are enumerated with PSB_L'
 * @param sigma alphabet size
 * @return sigma / ratio in runs, where ratio = (7 * floor(log_2(sigma)) - 29) / 2, limited to [1, 38]
 */
inline uint64_t max_width_scan_extensions(uint64_t sigma)
{
    const uint64_t lg = std::bit_width(sigma) - 1;
    const uint64_t ratio = std::clamp<uint64_t>(7 * lg >= 29 ? (7 * lg - 29) / 2 : 1, 1, 38);
    return std::max<uint64_t>(1, sigma / ratio);
}

/**
 * @brief returns the maximum number of distinct symbols that are reported with L_prev_min when enumerating
 *        the extensions of a pattern; if a range of L' contains more, the enumeration is repeated with one
 *        of the other methods (see max_width_scan_extensions)
 * @param sigma alphabet size
 * @param width width of the range of L' in runs
 * @return maximum number of symbols: width / 2, and at most sigma / 128 if the range is wider than
 *         max_width_scan_extensions(sigma); 0 <=> L_prev_min is not queried
 */
inline uint64_t L_prev_min_max_symbols(uint64_t sigma, uint64_t width)
{
    if (width < min_width_L_prev_min) return 0;
    if (width <= max_width_scan_extensions(sigma)) return width / 2;
    return std::min<uint64_t>(width / 2, sigma / 128);
}

/**
 * @brief parameters for building the index of the text described by a FASTA file; each header line and each
 *        run of symbols outside the allowed alphabet is replaced by one separator (see pfp_rlbwt::fasta_transform)
 */
struct fasta_input_t {
    uint64_t text_len = 0; // length of the text (0 <=> the input is not a FASTA file)
    // [0..255] number of occurrences of each byte in the text (nullptr: they are read from the pfp cache,
    // which then stores the sequences, see pfp_cache_fasta_section)
    const uint64_t* byte_freqs = nullptr;
    uint8_t separator = uint8_t(fasta_separator_symbol); // the separator symbol
    std::string alphabet = "ACGT"; // the symbols that are kept
    // sampled positions in the text with their offsets in the file (see pfp_rlbwt::fasta_sample); if nullptr,
    // pfp_rlbwt::fasta_source computes them
    const pfp_rlbwt::fasta_sample* index_samples = nullptr;
    uint64_t num_index_samples = 0; // number of samples in index_samples
    // the sequences of the FASTA file; stored in a pfp cache (see pfp_cache_fasta_section)
    const fasta_sequence_data* sequences = nullptr;
};

// log_data_structures of move_r and move_rb prints the arrays only for inputs up to this length (debugging)
inline constexpr uint64_t max_logged_input_length = 50;

// the largest number of input symbols of the construction modes sa and sa-space: they build the suffix array of
// the n = input symbols + 1 symbols (with the terminator) with libsais40
inline constexpr uint64_t move_r_max_input_length_sa = libsais40_max_n - 1;

// the largest number of input symbols: the positions in [0, n) are stored with at most 48 bits (see uint48_t)
inline constexpr uint64_t move_r_max_input_length = uint48_t::max - 1;

// default of move_r_params::revert_start_samples: revert_pieces_per_thread pieces for each of revert_max_threads
// threads
inline constexpr uint64_t revert_max_threads = 64;
inline constexpr uint64_t revert_pieces_per_thread = 64;
inline constexpr uint64_t default_revert_start_samples = revert_max_threads * revert_pieces_per_thread;

// the number of pieces per thread SA_range_in_pieces of move_r splits a range into (the same as for revert)
inline constexpr uint64_t sa_pieces_per_thread = revert_pieces_per_thread;

/**
 * @brief move-r construction parameters
 */
struct move_r_params {
    // file to build the index from; if empty, the index is built from the container passed to the constructor
    // (else only the symbol type of the container is used)
    std::string input_file = "";
    // file in whose directory the temporary files are stored; if empty, input_file is used (see tmp_files_beside)
    std::string tmp_files_beside = "";
    move_r_construction_mode mode = _suffix_array; // construction mode to use
    uint64_t num_threads = omp_get_max_threads(); // maximum number of threads to use during the construction
    uint64_t a = default_balancing_parameter; // balancing parameter, 2 <= a
    // (rlzsa only) maximum length of the reference of the rlz; 0 <=> chosen by the construction
    uint64_t rlz_max_reference_len = 0;
    /* alphabet size of the input; if set to 0, the symbols in the input are mapped to a compact alphabet; else
       (alphabet_size != 0), the input must already have a compact alphabet, and the symbols are not remapped */
    uint64_t alphabet_size = 0;
    /* true <=> the symbols are mapped in decreasing order of their frequencies (the most frequent symbol is
       mapped to the smallest internal symbol); false <=> the mapping preserves the order of the symbols */
    bool reorder_alphabet = true;
    // number of runs per block of PSB_L' in multiples of sigma (see pred_succ_block)
    uint64_t pred_succ_block_factor = 8;
    // number of occurrences of a symbol per block of PSO_L' (see pred_succ_occ)
    uint64_t pso_occs_per_block = 16;
    // minimum alphabet size for which L_prev_min is built (see move_r::L_prev_min)
    uint64_t min_sigma_L_prev_min = 4096;
    bool log = false; // controls, whether to print log messages
    std::ostream* measurement_file = nullptr; // measurement file for the index construction
    std::ostream* mf_mds = nullptr; // measurement file for the move data structure construction
    // suffix appended to every measurement_file key (e.g. "_fwd"/"_bwd" for the sub-indexes of a bidirectional index)
    std::string measurement_key_suffix = "";
    std::string text_name = ""; // name of the input file (used only for measurement output)
    fasta_input_t fasta = {}; // set if the input is a FASTA file (see fasta_input_t)

    /* false <=> with _pfp, no pfp cache of input_file is read or written (see move_r::pfp_cache_path); the output of
       the prefix-free parsing is written to temporary files (see tmp_files_beside), which are removed after use */
    bool pfp_cache = true;

    // true <=> an existing pfp cache of the input is not read, but overwritten
    bool ignore_pfp_cache = false;

    // true <=> the index of the reversed input is built (with _pfp, the input file is read from back to front)
    bool reversed = false;

    // false <=> with _pfp, the construction does not sample the suffix array at the ends of the BWT runs
    // (see move_r::has_run_end_samples)
    bool sample_run_ends = true;

    // the text is split into this many ranges, and the index stores a revert start in each range that contains a
    // BWT run start (see move_r::revert_starts); 0 <=> the index stores none and reverts with one thread
    uint64_t revert_start_samples = default_revert_start_samples;

    // (move_rb only) maximum size of the serialized context tables as a share of the size of the serialized index;
    // 0 <=> none are built
    // (see move_rb::build_context_tables)
    double context_table_share = 0.05;

    // (move_rb only) true <=> the context tables store only the patterns of the maximum length that
    // context_table_share permits, instead of all lengths up to it (see context_tables::single_len)
    bool context_tables_single_length = false;

    // ##################### SET BY A BIDIRECTIONAL INDEX FOR ITS SUB-INDEXES (DON'T USE DIRECTLY) #####################
    uint64_t* peak_memory_usage = nullptr; // variable to store the peak memory usage in
    // file storing the [0..255] alphabet map that is applied to the bytes of the input file while the prefix-free
    // parsing reads it (see construction::pfp_alphabet_map)
    std::string pfp_alphabet_map_file = "";
};

/**
 * @brief returns the file in whose directory the temporary files of the construction are stored
 * @param params construction parameters
 * @return params.tmp_files_beside, if set, else params.input_file; if empty, the files are stored in the
 *         temporary directory (see temporary_file_path)
 */
inline const std::string& tmp_files_beside(const move_r_params& params)
{
    return params.tmp_files_beside.empty() ? params.input_file : params.tmp_files_beside;
}

// empty type of a member that an index type does not have; [[no_unique_address]] members with distinct tags
// can share an address
template <uint8_t tag> struct absent_field_t { };

/**
 * @brief move-r index of a text T$ of length n, size O(r * a/(a-1)) words
 *
 * The index stores M_LF, the move data structure for LF over the BWT L of T$, with L' (the symbol of each of its
 * r' input intervals) as row 0, and PSO_L' for the backward search steps. Depending on the locate support, it also
 * stores the suffix array samples at the starts (and ends) of the r BWT runs, M_Phi^{-1} with r'' and M_Phi with r'''
 * input intervals, or an rlz of the differential suffix array. The symbols of T are mapped to the internal symbols [1..sigma-1]; 0 is the terminator $.
 *
 * Symbols are passed to and returned by the index as uint64_t. Whether the input uses a byte or an integer
 * alphabet is stored in the index (see byte_alphabet); the symbol type of the input is a template parameter
 * of the constructors only.
 *
 * @tparam support type of locate support
 * @tparam mlf_enc M_LF's interval-position encoding type (default: EF for _count to _locate_rlzsa, else POS); with
 *         EF, the backward search moves (interval, offset)-pairs (see move_data_structure::move_offset)
 */
template <move_r_support support = _locate_move, move_pos_encoding_t mlf_enc = (support <= _locate_rlzsa ? EF : POS)>
class move_r {
public:
    // ############################# TYPES, CONSTANTS AND PARAMETERS #############################

    // true <=> M_Phi and M_Phi^{-1} are stored; M_Phi^{-1} stores the PLCP value at the start of each input interval
    // as row 0
    static constexpr bool has_mphi = support == _locate_move_bi;
    static_assert(support < _count_bi || mlf_enc == POS,
        "a sub-index of a bidirectional index stores the positional M_LF encoding");

    // true <=> the suffix array samples are the only locate data structure the index stores (see _locate_one_bi)
    static constexpr bool stores_only_SA_samples = support == _locate_one_bi;

    // type of M_LF; it stores L' (the BWT run heads) as row 0
    using mlf_t = move_data_structure<mlf_enc, UNI>;
    // type of M_Phi and M_Phi^{-1}: the EF encoding (see move_pos_encoding_t)
    using mphi_t = move_data_structure<EF, UNI>;
    // number of additional rows: L' for M_LF, the PLCP values for M_Phi^{-1}
    static constexpr uint8_t mlf_num_rows = 1;
    static constexpr uint8_t mphi_inv_num_rows = has_mphi ? 1 : 0;

    // true <=> the index is a sub-index of a bidirectional index (see move_rb)
    static constexpr bool is_bidirectional = support >= _count_bi;

    // true <=> all extensions of a pattern can be enumerated at once (see move_rb's prepare_extend_all);
    // PSB_L' is then built
    static constexpr bool enumerates_extensions = is_bidirectional;

    // true <=> the index supports locate
    static constexpr bool supports_locate = support != _count && support != _count_bi;

    // true <=> the suffix array samples at the BWT run boundaries are built (see _SA_run)
    static constexpr bool builds_SA_samples = supports_locate;
    // true <=> the index supports locating multiple occurrences
    static constexpr bool supports_multiple_locate =
        supports_locate && support != _locate_one && !stores_only_SA_samples;
    // true <=> the index has an rlz of the differential suffix array
    static constexpr bool has_SA_rlz = support == _locate_rlzsa || support == _locate_rlzsa_bi;
    // true <=> the index supports random access to the suffix array (see SA)
    static constexpr bool has_SA_access = has_SA_rlz;

    // true <=> the suffix array is also sampled at the ends of the BWT runs (see move_r_params::sample_run_ends)
    static constexpr bool has_run_end_samples = is_bidirectional && builds_SA_samples;
    // true <=> the index contains move data structures for Phi/Phi^{-1}
    static constexpr bool has_locate_move = support == _locate_move || has_mphi;

protected:
    // reads L'[x] (for PSB_L' and PSO_L')
    struct L_reader {
        const move_r* idx = nullptr;
        inline uint64_t operator()(uint64_t x) const { return idx->L_(x); }
    };

    using psb_l_t = pred_succ_block<L_reader>;
    using pso_l_t = pred_succ_occ<L_reader>;
    using L_prev_min_t = rmq_tree<L_prev_min_arity>;

    // layout of the suffix array samples (also the value written by serialize_SA_samples): in SA_run and
    // SA_run_pos, or as the input intervals of M_Phi^{-1} in SA_Phi^{-1}
    enum sa_samples_layout_t : uint8_t { _at_run_boundaries = 1, _as_phi_inv_intervals = 2 };

    class construction;

    // ############################# INDEX VARIABLES #############################

    uint64_t n = 0; // the length of the input T$ (including the terminator)
    uint64_t sigma = 0; // the number of distinct symbols (including the terminator symbol 0) of T$
    uint64_t r = 0; // r, the number of runs in L
    uint64_t r_ = 0; // r', the number of input/output intervals in M_LF
    uint64_t r__ = 0; // r'', the number of input/output intervals in M_Phi^{-1}
    uint64_t r___ = 0; // r''', the number of input/output intervals in M_Phi
    uint64_t a = 0; // balancing parameter a >= 2; M_LF has at most r * a/(a-1) input intervals
    uint8_t omega_idx = 0; // bit width of an entry of SA_Phi^{-1}
    bool symbols_remapped = false; // true <=> the symbols of the input are mapped to internal symbols

public:
    // true <=> the input uses a byte alphabet; the alphabet is then mapped with _map_int and _map_ext, else
    // with _alphabet
    bool byte_alphabet = true;

protected:
    // ############################# INDEX DATA STRUCTURES #############################

    // the map of the bytes to the internal symbols (byte alphabet): [0..255] _map_int[c] = internal symbol of the byte
    // c (0 if c does not occur); empty for an integer alphabet
    std::vector<uint8_t> _map_int;
    // the map of the internal symbols to the bytes (byte alphabet): [0..255] _map_ext[s] = byte of the internal symbol
    // s (0 for the terminator and for values that are no internal symbol); empty for an integer alphabet
    std::vector<uint8_t> _map_ext;

    /* the distinct input symbols (integer alphabet): an Elias-Fano bit vector over [0..largest symbol] with a one at
       each distinct symbol except UINT64_MAX; the rank k of a symbol (its one is the k-th one) is mapped to its
       internal symbol with _alphabet_order */
    elias_fano_bit_vector _alphabet;
    // true <=> UINT64_MAX occurs in the input; it is not marked in _alphabet and has the largest rank
    bool _alphabet_has_max = false;

    // [0..sigma-2] _alphabet_order[k - 1] = internal symbol of the k-th smallest distinct symbol of the
    // input; empty <=> the internal symbol is k (see move_r_params::reorder_alphabet)
    interleaved_bit_aligned_vectors _alphabet_order;
    // [0..sigma-2] the inverse of _alphabet_order, i.e., _alphabet_order_inv[s - 1] = k <=> _alphabet_order[k - 1] = s
    interleaved_bit_aligned_vectors _alphabet_order_inv;

    // M_LF, the move data structure for LF; row 0: [0..r'-1] L'[x] = the symbol of the x-th input interval (see L_)
    mlf_t _M_LF;

    // data structures over L'; which of them are built depends on sigma and on the supported operations
    // (see build_data_structures_over_l_)
    psb_l_t _PSB_L_; // symbol predecessor/successor support for enumerating all symbols (see pred_succ_block)
    pso_l_t _PSO_L_; // symbol predecessor/successor support for backward search steps (see pred_succ_occ)
    wavelet_tree_prefix_sum _WTPS_L_; // wavelet tree for counting the smaller symbols in a range
    /* range-minimum tree over [0..r'-1] prev(x) = 1 + the largest y < x with L'[y] = L'[x] (0 if there is none); the
       distinct symbols of L'[b..e] are the symbols L'[x] with x in [b, e] and prev(x) <= b */
    L_prev_min_t _L_prev_min;

    // M_Phi^{-1}, the move data structure for Phi^{-1}; row 0 (has_mphi): [0..r''-1] the PLCP value at the start of
    // the x-th input interval
    mphi_t _M_Phi_inv;
    // M_Phi, the move data structure for Phi (has_mphi)
    mphi_t _M_Phi;
    /* [0..r'-1] if the x-th input interval of M_LF starts a BWT run, SA_Phi^{-1}[x] is the index of the input
       interval of M_Phi^{-1} that starts with SA[(M_LF.p[x] - 1) mod n] (its output interval starts with
       SA[M_LF.p[x]]); else SA_Phi^{-1}[x] = r''; omega_idx bits per entry */
    interleaved_bit_aligned_vectors _SA_Phi_inv;

    /* for the input intervals x in [0..r'-1] of M_LF: the indices in _SA_run of the sample at the start of the run
       starting at x (start_sample_index(x)) and of the sample at the end of the run ending at x (end_sample_index(x)); see
       run_sample_positions */
    run_sample_positions _SA_run_pos;
    /* the suffix array samples at the BWT run boundaries, run by run: SA at the first position of the run, then (if
       the run has length >= 2 and the ends of the runs are sampled) SA at its last position; bit_width(n) bits per
       sample; empty if the samples are stored in SA_Phi^{-1} */
    interleaved_bit_aligned_vectors _SA_run;

    // relative Lempel-Ziv encoding of the differential suffix array (see rlz_ref/rlz.hpp)
    rlz<_diff_access> _SA_rlz;

    // the starting positions, names and separator symbol of the sequences, empty unless the index was built from a
    // FASTA file (see data_structures/fasta_sequence_data.hpp)
    fasta_sequence_data _seq_data;

    /* the revert starts, sorted by text position: [0..m-1] field 0: a BWT run start i, field 1: SA[i], with
       1 <= SA[i] <= n-2; SA[i] is the smallest such value at a run start in its range, where [0, n) is split into
       move_r_params::revert_start_samples ranges of equal length; bit_width(n) bits per field */
    interleaved_bit_aligned_vectors _revert_starts;

    // true <=> the construction read the output of the prefix-free parsing from the pfp cache; not serialized
    bool _rlbwt_from_cache = false;

    // bit width of the PLCP row of M_Phi^{-1}; 0 <=> not known before it is built; not serialized
    uint8_t _plcp_row_width = 0;

    // ############################# INTERNAL METHODS #############################

    // member-wise copy/move of all members, then rebinds the data structures over L' that read it
    template <typename move_r_t>
    void assign_from(move_r_t&& other);

    /**
     * @brief builds the index from an input
     * @tparam input_t container type of the input
     * @param input the input
     * @param delete_input whether the input may be deleted during the construction
     * @param params construction parameters
     */
    template <typename input_t>
    void construct_from_input(input_t& input, bool delete_input, move_r_params params);

    /**
     * @brief builds _alphabet, _alphabet_order and _alphabet_order_inv (integer alphabet)
     * @param symbols the distinct symbols of the input in ascending order
     * @param order [0..sigma-2] order[i] = internal symbol of symbols[i]; if empty, the internal symbol of
     *        symbols[i] is i + 1
     */
    void assign_alphabet(const std::vector<uint64_t>& symbols, const std::vector<uint64_t>& order);

    /**
     * @brief returns the bit widths of the additional rows of M_Phi^{-1}
     * @param n length of the input
     * @return the width of the PLCP row (_plcp_row_width, or the bit width of n if it is 0); empty if there
     *         is no PLCP row
     */
    typename mphi_t::row_widths_t mphi_inv_row_widths(uint64_t n) const
    {
        if constexpr (has_mphi) {
            return { _plcp_row_width == 0 ? uint8_t(std::bit_width(n)) : _plcp_row_width };
        } else return {};
    }

    // returns the size of the data structures over L' in bytes
    inline uint64_t l_data_structures_bytes() const
    {
        return _PSB_L_.size_in_bytes() + _PSO_L_.size_in_bytes() + _WTPS_L_.size_in_bytes()
            + _L_prev_min.size_in_bytes();
    }

    // returns the size of the map from the internal symbols to the input symbols in bytes
    inline uint64_t map_ext_bytes() const
    {
        return byte_alphabet ? _map_ext.size() : _alphabet.size_in_bytes()
            + _alphabet_order.size_in_bytes() + _alphabet_order_inv.size_in_bytes();
    }

    /**
     * @brief sets s = SA[M_LF.p[x]] and s_ to the index of the input interval of M_Phi^{-1} containing s
     * @param x [0..r'-1] an input interval of M_LF that starts a BWT run
     * @param s variable to store SA[M_LF.p[x]] in
     * @param s_ variable to store the index of the input interval of M_Phi^{-1} containing s in
     */
    inline void init_phi_inv_at_run_start(uint64_t x, uint64_t& s, uint64_t& s_) const requires(has_locate_move);

    /**
     * @brief sets s = SA[b] and s_ to the index of the input interval of M_Phi^{-1} containing s, after a backward
     *        search (see backward_search_step)
     * @param b left limit of the suffix array interval
     * @param e right limit of the suffix array interval
     * @param s variable to store SA[b] in
     * @param s_ variable to store the index of the input interval of M_Phi^{-1} containing s in
     * @param hat_b_ap_y \hat{b}'_y
     * @param y y
     */
    inline void init_phi_inv(
        uint64_t b, uint64_t e,
        uint64_t& s, uint64_t& s_,
        uint64_t hat_b_ap_y, int64_t y) const requires(has_locate_move);

    // writes the layout of the suffix array samples (1 byte, see sa_samples_layout_t) and the samples to an output
    // stream
    void serialize_SA_samples(std::ostream& out) const;

    // reads the suffix array samples written by serialize_SA_samples
    void read_SA_samples(std::istream& in);

public:
    // ############################# CONSTRUCTORS #############################

    move_r() = default;

    move_r(const move_r& other) { assign_from(other); }
    move_r(move_r&& other) noexcept { assign_from(std::move(other)); }
    move_r& operator=(const move_r& other) { assign_from(other); return *this; }
    move_r& operator=(move_r&& other) noexcept { assign_from(std::move(other)); return *this; }

    /**
     * @brief constructs a move_r index of the input (consumed during the construction if passed as an rvalue)
     * @tparam input_t container type of the input: a contiguous container of an integral symbol type (e.g.,
     *         std::string, std::vector<uint8_t> or std::vector<int32_t>); a symbol type of one byte
     *         builds an index with a byte alphabet, a wider one an index with an integer alphabet (see
     *         byte_alphabet); signed symbols are reinterpreted (see symbol_idx)
     * @param input the input
     * @param params construction parameters
     */
    template <typename input_t>
    move_r(input_t& input, move_r_params params = {})
    {
        construct_from_input(input, false, std::move(params));
    }

    template <typename input_t>
    move_r(input_t&& input, move_r_params params = {})
    {
        construct_from_input(input, true, std::move(params));
    }

    // constructs a move_r index of a null-terminated string
    move_r(const char* input, move_r_params params = {}) : move_r(std::string(input), params) {}

    // ############################# MISC PUBLIC METHODS #############################

    /**
     * @brief returns the name of the file in which the output of the prefix-free parsing is cached (the pfp cache)
     * @param input_file file the index is built from
     * @param mode construction mode
     * @param reversed whether the index of the reversed input is built (see move_r_params)
     * @return input_file + ".bwt-runs" (or ".bwt-runs-rev" if reversed); empty if mode != _pfp or
     *         input_file is empty
     */
    static std::string pfp_cache_path(const std::string& input_file, move_r_construction_mode mode,
        bool reversed = false)
    {
        if (mode != _pfp || input_file.empty()) return "";
        return input_file + pfp_rlbwt::ext_runs + (reversed ? "-rev" : "");
    }

    /**
     * @brief returns the length n of the text T$ an index is built for
     * @param T the input (empty, if it is read from params.input_file)
     * @param params construction parameters
     */
    static uint64_t input_length(const move_r_input_view& T, const move_r_params& params);

    /**
     * @brief reinterprets a symbol of an integral type as an unsigned value of the same width
     * @tparam sym_t symbol type
     * @param sym symbol
     * @return sym reinterpreted as an unsigned value
     */
    template <typename sym_t>
    inline static uint64_t symbol_idx(sym_t sym)
    {
        static_assert(std::is_integral_v<sym_t> && sizeof(sym_t) <= sizeof(uint64_t),
            "move_r supports integral symbols of at most 64 bits");
        return symbol_bit_pattern(sym);
    }

    /**
     * @brief sets the internal and external alphabet maps (only for bidirectional indexes)
     * @param map_int [0..255] the map from the bytes to the internal symbols (see _map_int)
     * @param map_ext [0..255] the map from the internal symbols to the bytes (see _map_ext)
     */
    void set_alphabet_maps(const std::vector<uint8_t>& map_int, const std::vector<uint8_t>& map_ext)
        requires(is_bidirectional);

    /**
     * @brief sets the alphabet of an index with an integer alphabet (only for bidirectional indexes)
     * @param symbols the distinct symbols of the input in ascending order
     * @param order [0..sigma-2] order[i] = internal symbol of symbols[i]; if empty, the internal symbol
     *        of symbols[i] is i + 1 (see _alphabet_order)
     */
    void set_alphabet(const std::vector<uint64_t>& symbols,
        const std::vector<uint64_t>& order = {}) requires(is_bidirectional);

    /**
     * @brief builds _alphabet_order and _alphabet_order_inv
     * @param order [0..sigma-2] order[i] = internal symbol of the (i+1)-th smallest distinct symbol of the
     *        input; if empty, the alphabet is not reordered
     */
    void set_alphabet_order(const std::vector<uint64_t>& order);

    // returns whether the alphabet is reordered (see _alphabet_order)
    inline bool alphabet_reordered() const { return !_alphabet_order.empty(); }

    /**
     * @brief returns whether the internal symbols are in the order of the input symbols (see
     *        move_r_params::reorder_alphabet); else they are ordered by decreasing frequency
     */
    bool alphabet_in_input_order() const;

    // returns the length of the input (without the terminator symbol)
    inline uint64_t input_size() const { return n - 1; }
    // returns the number of distinct symbols in the input (without the terminator symbol)
    inline uint32_t alphabet_size() const { return sigma - 1; }
    // returns r, the number of runs in the BWT
    inline uint64_t num_bwt_runs() const { return r; }
    // returns the balancing parameter a
    inline uint64_t balancing_parameter() const { return a; }
    // returns the bit width of an entry of SA_Phi^{-1}
    inline uint8_t width_SA_Phi_inv() const requires(has_locate_move) { return omega_idx; }

    // returns the starting positions, names and separator symbol of the sequences of a FASTA input (see _seq_data)
    inline const fasta_sequence_data& seq_data() const { return _seq_data; }

    // returns whether the construction read the output of the prefix-free parsing from the pfp cache
    inline bool rlbwt_from_cache() const { return _rlbwt_from_cache; }

    // returns the number of revert starts (see _revert_starts)
    inline uint64_t num_revert_starts() const { return _revert_starts.size(); }
    // sets the sequences of a FASTA input (as computed by process_fasta); call after the construction
    inline void set_fasta_sequence_data(fasta_sequence_data seq_data) { _seq_data = std::move(seq_data); }

    // returns the rlz of the differential suffix array
    inline const rlz<_diff_access>& SA_rlz() const requires(has_SA_rlz) { return _SA_rlz; }
    // returns the number of phrases in the rlz
    inline uint64_t num_SA_rlz_phrases() const requires(has_SA_rlz) { return _SA_rlz.num_phrases(); }
    // returns the number of literal phrases in the rlz
    inline uint64_t num_SA_rlz_literal_phrases() const requires(has_SA_rlz) { return _SA_rlz.num_literal_phrases(); }
    // returns the number of copy phrases in the rlz
    inline uint64_t num_SA_rlz_copy_phrases() const requires(has_SA_rlz) { return _SA_rlz.num_copy_phrases(); }

    /**
     * @brief returns the size of the data structure in bytes
     * @return size of the data structure in bytes
     */
    uint64_t size_in_bytes() const;

    /**
     * @brief logs all data structures
     */
    void log_data_structures() const;

    /**
     * @brief logs the sizes of the data structures of the index with their shares of the index size, one per line
     *        (see log_component); unless a listing is in progress (e.g. of a move_rb), the list ends with the size of
     *        the rest of the index ("other")
     */
    void log_data_structure_sizes() const;

    /**
     * @brief adds the parameters the size of the index depends on to the facts of the info modes, one per line:
     *        n and sigma (with the terminator), r, n/r, k' of each move data structure, the number of phrases of the
     *        RLZ encoding and the number of FASTA sequences
     * @param facts the facts
     * @param suffix appended to the names of the quantities, e.g. "^R" for the index of the reversed text; n and sigma
     *        are only added if it is empty
     */
    void add_parameters(index_facts& facts, const std::string& suffix = "") const;

    /**
     * @brief logs the index data structure sizes to the output stream out
     * @param out an output stream
     * @param suffix suffix appended to each measurement key
     */
    void log_data_structure_sizes(std::ostream& out, std::string suffix = "") const;

    /**
     * @brief splits the input intervals of M_LF into parts of consecutive BWT runs
     * @param num_parts number of parts
     * @param p number of threads to use
     * @return [0..num_parts] (part_begin[t], part_run[t]): the first input interval of part t, which starts a BWT
     *         run, and the index of that run; a part can be empty; part_begin[num_parts] = r' and part_run[num_parts]
     *         = r
     */
    std::pair<std::vector<uint64_t>, std::vector<uint64_t>> split_into_run_parts(uint64_t num_parts, uint64_t p) const;

    // ############################# ACCESS METHODS #############################

    inline const mlf_t& M_LF() const { return _M_LF; }
    inline const mphi_t& M_Phi_inv() const requires(has_locate_move) { return _M_Phi_inv; }
    inline const mphi_t& M_Phi() const requires(has_mphi) { return _M_Phi; }

    // symbol predecessor/successor support over L' (see pred_succ_block)
    inline const psb_l_t& PSB_L_() const { return _PSB_L_; }
    // symbol predecessor/successor support over L' (see pred_succ_occ)
    inline const pso_l_t& PSO_L_() const { return _PSO_L_; }
    // wavelet tree over L' (see wavelet_tree_prefix_sum)
    inline const wavelet_tree_prefix_sum& WTPS_L_() const { return _WTPS_L_; }
    // range-minimum tree over the previous occurrences of the symbols in L' (see _L_prev_min)
    inline const L_prev_min_t& L_prev_min() const { return _L_prev_min; }

    // returns SA_Phi^{-1}[x] (see _SA_Phi_inv)
    inline uint64_t SA_Phi_inv(uint64_t x) const requires(support == _locate_move) { return _SA_Phi_inv[x]; }

    // returns the positions in SA_run of the samples of the input intervals of M_LF (see run_sample_positions)
    inline const run_sample_positions& SA_run_pos() const { return _SA_run_pos; }

    /**
     * @brief returns the suffix array sample at the starting position of the x-th input interval of M_LF
     * @param x [0..r'-1] an input interval of M_LF that starts a BWT run
     * @return SA[M_LF.p[x]]
     */
    inline uint64_t SA_run_start(uint64_t x) const
        requires(builds_SA_samples);

    /**
     * @brief returns the suffix array sample at the end position of the x-th input interval of M_LF
     * @param x [0..r'-1] an input interval of M_LF that ends a BWT run
     * @return SA[M_LF.p[x + 1] - 1]
     */
    inline uint64_t SA_run_end(uint64_t x) const
        requires(support != _count && support != _count_bi && support != _locate_rlzsa && support != _locate_one);

    /**
     * @brief returns L'[x]
     * @param x [0..r'-1]
     * @return L'[x]
     */
    inline uint64_t L_(uint64_t x) const { return _M_LF.template row<0>(x); }

    /**
     * @brief LF step with M_LF: i is a position (M_LF in POS) or the offset of a position within the x-th input
     *        interval of M_LF (M_LF in EF, see move_data_structure::move_offset)
     * @param i position or offset of j; is set to the position or offset of LF(j)
     * @param x [0..r'-1] index of the input interval of M_LF containing j; is set to the one containing LF(j)
     */
    inline void LF_step(uint64_t& i, uint64_t& x) const
    {
        if constexpr (mlf_enc == POS) _M_LF.move(i, x);
        else _M_LF.move_offset(i, x);
    }

    // returns whether the x-th input interval of M_LF starts a BWT run
    inline bool run_starts_at(uint64_t x) const { return x == 0 || L_(x - 1) != L_(x); }

    // returns whether the x-th input interval of M_LF ends a BWT run
    inline bool run_ends_at(uint64_t x) const { return x + 1 == r_ || L_(x + 1) != L_(x); }

    // returns the index of the first input interval of M_LF in the BWT run containing the x-th input interval
    inline uint64_t first_interval_of_run(uint64_t x) const
    {
        while (!run_starts_at(x)) x--;
        return x;
    }

    // returns the index of the last input interval of M_LF in the BWT run containing the x-th input interval
    inline uint64_t last_interval_of_run(uint64_t x) const
    {
        while (!run_ends_at(x)) x++;
        return x;
    }

    // returns the [0..255] map from the bytes to the internal symbols (byte alphabet only, see _map_int)
    inline const std::vector<uint8_t>& map_int() const { return _map_int; }

    /**
     * @brief maps a symbol to its internal symbol
     * @tparam sym_t symbol type
     * @param sym symbol
     * @return its internal symbol (0, the internal terminator, if sym is TERMINATOR_SYMBOL or does not occur in the
     *         input)
     */
    template <typename sym_t>
    inline uint64_t map_symbol(sym_t sym) const;

    /**
     * @brief maps an internal symbol to its input symbol
     * @param sym an internal symbol
     * @return its input symbol (TERMINATOR_SYMBOL, if sym is the terminator 0)
     */
    inline uint64_t unmap_symbol(uint64_t sym) const;

    // returns an internal symbol as the log_data_structures functions print it: the terminator (0) as '$', a symbol
    // of a byte alphabet as its character, else its input symbol
    inline std::string logged_symbol(uint64_t sym) const
    {
        if (sym == 0) return "$";
        if (byte_alphabet) return std::string(1, char(unmap_symbol(sym)));
        return std::to_string(unmap_symbol(sym));
    }

    // ############################# METHODS FOR MOVE_RB (SUB-INDEXES) #############################

    // returns the name of the data structure storing the suffix array samples (for log messages)
    const char* SA_samples_name() const
    {
        return _SA_run_pos.empty() && !_SA_Phi_inv.empty() ? "SA_Phi^{-1}" : "SA_run";
    }

    /**
     * @brief writes the suffix array samples to an output stream and deletes them from memory (see
     *        load_SA_samples)
     * @param out output stream to write the samples to
     */
    void store_SA_samples(std::ostream& out)
    {
        serialize_SA_samples(out);
        _SA_run_pos = run_sample_positions();
        _SA_run = interleaved_bit_aligned_vectors();
        if constexpr (support != _locate_move) _SA_Phi_inv = interleaved_bit_aligned_vectors();
    }

    /**
     * @brief reads the suffix array samples written by store_SA_samples
     * @param in input stream to read the samples from
     */
    void load_SA_samples(std::istream& in) { read_SA_samples(in); }

    /**
     * @brief writes M_Phi to an output stream and deletes it from memory (see load_M_Phi)
     * @param out output stream to write M_Phi to
     */
    void store_M_Phi(std::ostream& out) requires(has_mphi)
    {
        _M_Phi.serialize(out);
        _M_Phi = mphi_t();
    }

    /**
     * @brief reads M_Phi written by store_M_Phi
     * @param in input stream to read M_Phi from
     */
    void load_M_Phi(std::istream& in) requires(has_mphi) { _M_Phi.load(in, 0); }

    // ############################# QUERY METHODS #############################

    /**
     * @brief returns L[i], with the terminator as TERMINATOR_SYMBOL
     * @param i [0..input size]
     * @return L[i]
     */
    inline uint64_t BWT(uint64_t i) const;

    /**
     * @brief returns SA[i] (with random access to the suffix array, see has_SA_access)
     * @param i [0..input size]
     * @return SA[i]
     */
    uint64_t SA(uint64_t i) const
        requires(has_SA_access);

    /**
     * @brief reports SA[b..e] (with random access to the suffix array, see has_SA_access)
     * @tparam fnc_t type of report
     * @param b first position
     * @param e last position (b <= e)
     * @param report function called as report(SA[i]) for i = b, ..., e
     */
    template <typename fnc_t>
    void SA_range(uint64_t b, uint64_t e, fnc_t report) const
        requires(has_SA_access);

    /**
     * @brief the state of a backward search for a pattern P that is extended to the left (see prepend), and of the
     *        locate query of its occurrences; the suffix array interval of P is [b, e]
     */
    struct query_context_t {
        /**
         * @brief stores the variables needed to locate via Phi^{-1}-move-queries
         */
        struct phi_move_context_t {
            uint64_t i;  // current position in the suffix array interval
            uint64_t s;  // current suffix SA[i]
            uint64_t s_; // index of the input interval of M_Phi^{-1} containing s
        };

    protected:
        uint64_t l; // length of the currently matched pattern
        // b and e (M_LF in POS), or their offsets in the input intervals b_ and e_ of M_LF (M_LF in EF)
        uint64_t b, e;
        // the input intervals of M_LF containing b and e, and \hat{b'}_y, \hat{e'}_z (see backward_search_step)
        uint64_t b_, e_, hat_b_ap_y, hat_e_ap_z;
        int64_t y, z; // variables for backward search
        bool be_valid; // true <=> b and e are positions (not offsets in their input intervals), see ensure_be

        // decoding context for the rlz (present iff the index has an rlz); dec.pos() is the current
        // position in the suffix array interval and dec.value() the current suffix SA[dec.pos()]
        STRING_UTILS_NO_UNIQUE_ADDRESS std::conditional_t<has_SA_rlz,
            typename rlz<_diff_access>::decode_context_t, absent_field_t<0>> dec;

        // context for the Phi^{-1}-move-based locate variant (present iff has_locate_move)
        STRING_UTILS_NO_UNIQUE_ADDRESS std::conditional_t<has_locate_move,
            phi_move_context_t, absent_field_t<1>> phi;

        const move_r<support, mlf_enc>* idx; // index to query

        /**
         * @brief returns the current position in the suffix array interval
         * @return the current position
         */
        inline uint64_t pos() const requires(has_SA_rlz || has_locate_move)
        {
            if constexpr (has_SA_rlz) return dec.pos();
            else return phi.i; // has_locate_move
        }

        /**
         * @brief sets the current position in the suffix array interval; for the rlz, this also (re-)binds
         * the decoding context to the index's rlz
         * @param i the position to set
         */
        inline void set_pos(uint64_t i);

        // computes the positions b and e and sets the position of the next occurrence to b, if they are not
        // valid (see be_valid)
        inline void ensure_be();

    public:
        /**
         * @brief constructs a new query context for the index idx
         * @param idx an index
         */
        query_context_t(const move_r<support, mlf_enc>& idx)
        {
            this->idx = &idx;
            reset();
        }

        /**
         * @brief resets the query context to an empty pattern
         */
        inline void reset()
        {
            idx->init_backward_search(b, b_, e, e_, hat_b_ap_y, y, hat_e_ap_z, z);
            l = 0;
            be_valid = false; // b,e reconstructed lazily by ensure_be()
        }

        // returns the length of the currently matched pattern
        inline uint64_t length() const { return l; }

        /**
         * @brief returns the overall number of occurrences of the currently matched pattern
         * @return overall number of occurrences
         */
        inline uint64_t num_occ()
        {
            if constexpr (mlf_enc == EF) {
                if (!be_valid) return idx->M_LF().distance(b_, b, e_, e) + 1;
            }
            return e - b + 1;
        }

        /**
         * @brief returns the number of remaining (not yet reported) occurrences of the currently matched pattern
         * @return number of remaining occurrences
         */
        inline uint64_t num_occ_rem()
            requires(supports_multiple_locate)
        {
            ensure_be();
            return e >= pos() ? e - pos() + 1 : 0;
        }

        /**
         * @brief returns the suffix array interval of the currently matched pattern
         * @return suffix array interval
         */
        inline std::pair<uint64_t, uint64_t> sa_interval()
        {
            ensure_be();
            return std::make_pair(b, e);
        }

        /**
         * @brief prepends sym to the currently matched pattern P; if symP occurs in the input, true is
         * returned and the query context is adjusted to store the information for the pattern symP; else,
         * false is returned and the query context is not modified
         * @tparam sym_t symbol type
         * @param sym the symbol to prepend
         * @return whether symP occurs in the input
         */
        template <typename sym_t>
        bool prepend(sym_t sym);

        /**
         * @brief reports the next occurrence of the currently matched pattern
         * @return next occurrence
         */
        inline uint64_t next_occ() requires(supports_multiple_locate);

        /**
         * @brief reports one occurrence of the currently matched pattern
         * @return an occurrence
         */
        inline uint64_t one_occ() const requires(supports_locate);

        /**
         * @brief locates the remaining (not yet reported) occurrences of the currently matched pattern
         * @tparam report_fnc_t type of the report function
         * @param report function that is called with every remaining occurrence as a parameter
         */
        template <typename report_fnc_t>
        void locate(report_fnc_t report) requires(supports_multiple_locate);

        /**
         * @brief locates the remaining (not yet reported) occurrences of the currently matched pattern
         * @return vector containing the occurrences
         */
        std::vector<uint48_t> locate() requires(supports_multiple_locate);
    };

    // returns a query context for the index
    inline query_context_t query() const { return query_context_t(*this); }

    /**
     * @brief initializes the variables of a backward search for the empty pattern: the suffix array interval
     *        [b, e] = [0, n-1] is stored as b, e (M_LF in POS) or as their offsets in the input intervals b_ and
     *        e_ of M_LF (M_LF in EF)
     * @param b b, or its offset in the b_-th input interval of M_LF
     * @param b_ index of the input interval of M_LF containing b
     * @param e e, or its offset in the e_-th input interval of M_LF
     * @param e_ index of the input interval of M_LF containing e
     * @param hat_b_ap_y \hat{b}'_y, an input interval of M_LF that starts a BWT run, with
     *        SA[b] = SA[M_LF.p[\hat{b}'_y]] - (y + 1) (see one_occ)
     * @param y y
     * @param hat_e_ap_z \hat{e}'_z, an input interval of M_LF that ends a BWT run, with
     *        SA[e] = SA[M_LF.p[\hat{e}'_z + 1] - 1] - (z + 1)
     * @param z z
     */
    inline void init_backward_search(
        uint64_t& b, uint64_t& b_,
        uint64_t& e, uint64_t& e_,
        uint64_t& hat_b_ap_y, int64_t& y,
        uint64_t& hat_e_ap_z, int64_t& z) const;

    /**
     * @brief prepends sym to the currently matched pattern P, adjusts the variables to store
     *        the query context for the pattern symP and returns whether symP occurs in the input
     * @tparam sym_t symbol type
     * @param sym next symbol to match
     * @param b,b_,e,e_,hat_b_ap_y,y,hat_e_ap_z,z the variables of the backward search (see init_backward_search)
     * @return whether symP occurs in the input
     */
    template <typename sym_t>
    bool backward_search_step(
        sym_t sym,
        uint64_t& b, uint64_t& b_,
        uint64_t& e, uint64_t& e_,
        uint64_t& hat_b_ap_y, int64_t& y,
        uint64_t& hat_e_ap_z, int64_t& z) const
    {
        return backward_search_step_mapped(map_symbol(sym), b, b_, e, e_, hat_b_ap_y, y, hat_e_ap_z, z);
    }

    /**
     * @brief performs a backward search step (see backward_search_step) with an internal symbol
     * @param i_sym next symbol to match, as an internal symbol (see map_symbol)
     * @param b,b_,e,e_,hat_b_ap_y,y,hat_e_ap_z,z the variables of the backward search (see init_backward_search)
     * @return whether the extended pattern occurs in the input
     */
    bool backward_search_step_mapped(
        uint64_t i_sym,
        uint64_t& b, uint64_t& b_,
        uint64_t& e, uint64_t& e_,
        uint64_t& hat_b_ap_y, int64_t& y,
        uint64_t& hat_e_ap_z, int64_t& z) const;

    /**
     * @brief returns the number of occurrences of P in the input
     * @tparam sym_t symbol type of the pattern
     * @param P the pattern to count in the input
     * @return the number of occurrences of P in the input
     */
    template <typename sym_t>
    inline uint64_t count(move_r_pattern_view<sym_t> P) const;

    // counts a pattern stored in a contiguous container (see move_r_pattern_view)
    template <typename pat_t>
        requires(!move_r_is_pattern_view<std::remove_cvref_t<pat_t>>)
    inline uint64_t count(const pat_t& P) const { return count(move_r_pattern_view_t<pat_t>(P)); }

    // counts a null-terminated pattern
    inline uint64_t count(const char* P) const { return count(std::string_view(P)); }

    /**
     * @brief locates the pattern P in the input
     * @tparam report_fnc_t type of the report function
     * @tparam sym_t symbol type of the pattern
     * @param P the pattern to locate in the input
     * @param report function that is called with every occurrence of P in the input as a parameter
     */
    template <typename report_fnc_t, typename sym_t>
    inline void locate(move_r_pattern_view<sym_t> P, report_fnc_t report) const
        requires(supports_multiple_locate);

    // locates a pattern stored in a contiguous container (see move_r_pattern_view)
    template <typename report_fnc_t, typename pat_t>
        requires(supports_multiple_locate && !move_r_is_pattern_view<std::remove_cvref_t<pat_t>>)
    inline void locate(const pat_t& P, report_fnc_t report) const
    { locate<report_fnc_t>(move_r_pattern_view_t<pat_t>(P), report); }

    // locates a null-terminated pattern
    inline std::vector<uint48_t> locate(const char* P) const
        requires(supports_multiple_locate)
    { return locate(std::string_view(P)); }

    /**
     * @brief locates the pattern P in the input
     * @tparam pat_t type of the pattern
     * @param P the pattern to locate in the input
     * @return a vector containing the occurrences of P in the input
     */
    template <typename pat_t>
    inline std::vector<uint48_t> locate(const pat_t& P) const
        requires(supports_multiple_locate)
    {
        std::vector<uint48_t> Occ;
        locate(P, [&](uint64_t, uint64_t occ){Occ.emplace_back(occ);});
        return Occ;
    }

    // ############################# RETRIEVE-RANGE METHODS #############################

    // default number of pieces each thread of SA_range_in_pieces decodes alternately
    static constexpr uint64_t default_sa_strands = 1;

    // parameters for the retrieve-range methods (BWT_range, SA_range)
    struct retrieve_params {
        // [l, r] range of positions (0 <= l <= r <= n-1); l > r: all positions
        uint64_t l = 1;
        uint64_t r = 0;
        uint64_t num_threads = omp_get_max_threads(); // maximum number of threads to use
        /* (only for the methods that write to a file) maximum number of bytes of the output buffers; -1: about
           (r - l + 1) / retrieve_output_bytes_per_buffer_byte bytes per output value */
        int64_t max_bytes_alloc = -1;
        /* (internal, set by retrieve_range) number of positions of an output buffer block; the ranges of the threads
           start at multiples of it (relative to l); 0: no alignment */
        uint64_t buffer_align = 0;
        // (SA_range) the number of pieces each thread decodes alternately (see SA_range_in_pieces)
        uint64_t num_strands = default_sa_strands;
    };

    // default of retrieve_params::max_bytes_alloc: the output buffers take 1 / this of the size of the output
    static constexpr uint64_t retrieve_output_bytes_per_buffer_byte = 500;
    // minimum size of the output buffer of a thread in bytes: one sector (see retrieve_range)
    static constexpr uint64_t retrieve_min_buffer_bytes = direct_io_alignment;

protected:
    /**
     * @brief adjusts retrieve parameters, ensures (0 <= l <= r <= range_max); if l > r, it sets [l,r] <- [0,range_max]
     * @param params retrieve parameters to adjust
     * @param range_max maximum value for r
     */
    inline static void adjust_retrieve_params(retrieve_params& params, uint64_t range_max);

    /**
     * @brief calls fnc(i_p) for each i_p in [0..p-1], in an OpenMP parallel region of p threads if p > 1
     * @tparam fnc_t type of fnc
     * @param p the number of threads (p >= 1)
     * @param fnc the function, called with the thread index i_p
     */
    template <typename fnc_t>
    inline static void run_parallel_threads(uint64_t p, fnc_t&& fnc);

    /**
     * @brief returns the number of threads to use for a range retrieval: at least 1, and at most the requested number
     *        of threads, the number of available threads and max_useful
     * @param requested the number of threads requested via the retrieve parameters
     * @param max_useful an upper bound on the number of threads that is specific to the retrieve method
     * @return the number of threads to use
     */
    inline static uint64_t retrieve_num_threads(uint64_t requested, uint64_t max_useful)
    {
        return std::max<uint64_t>(1, std::min<uint64_t>(
            { (uint64_t) omp_get_max_threads(), (uint64_t) requested, max_useful }));
    }

    /**
     * @brief returns the start position of the range of thread i_p when [l, r] is split into p ranges; thread i_p
     *        processes [thread_range_start(i_p, ...), thread_range_start(i_p + 1, ...) - 1] (which may be empty); the
     *        starts of the ranges 1..p-1 are rounded up to multiples of align positions (relative to l)
     * @param i_p a thread index in [0..p]
     * @param p the number of threads
     * @param l the left limit of the range
     * @param r the right limit of the range
     * @param align the number of positions of an output buffer block (0 or 1: no alignment)
     * @return the start position of the range of thread i_p (l for i_p = 0, r + 1 for i_p >= p)
     */
    inline static uint64_t thread_range_start(uint64_t i_p, uint64_t p, uint64_t l, uint64_t r, uint64_t align);

    /**
     * @brief returns the number of threads revert_in_pieces uses: at least 1, and at most the requested number of
     *        threads, the number of available threads and the number of pieces
     * @param requested the requested number of threads
     */
    inline uint64_t revert_num_threads(uint64_t requested) const
    {
        return retrieve_num_threads(requested, _revert_starts.size() + 1);
    }

    /**
     * @brief returns the number of threads SA_range_in_pieces uses: at least 1, and at most the requested number of
     *        threads, the number of available threads and the number of BWT runs in the range (each piece starts at a
     *        BWT run start, so a thread skips about n/r positions)
     * @param params the adjusted parameters
     */
    inline uint64_t SA_range_num_threads(const retrieve_params& params) const
    {
        return retrieve_num_threads(params.num_threads, uint64_t(((params.r - params.l + 1) * double(r)) / double(n)));
    }

    /**
     * @brief calls retrieve(report, params), which reports the values in the range [params.l, params.r], and writes
     *        each reported value at index i - params.l of the file file_name through one output buffer per thread
     * @tparam output_t type of an output value (the file stores sizeof(output_t) bytes per value)
     * @tparam invoke_t type of retrieve
     * @param retrieve function called as retrieve(report, params); report(i, value) is called with i in
     *        [params.l, params.r] from a thread with a thread number below params.num_threads
     * @param file_name name of the file to write the output to (created if it does not exist; resized to the
     *        output)
     * @param params parameters (with params.l <= params.r)
     */
    template <typename output_t, typename invoke_t>
    void retrieve_range(invoke_t retrieve, std::string file_name, retrieve_params params) const;

public:
    /**
     * @brief returns the bwt in the range [l,r] (0 <= l <= r <= input size), else
     * if l > r, then the whole bwt is returned (default); the terminator is stored as 0 (see terminator_as_zero),
     * so if the input contained 0, the output is not equal to the real bwt
     * @param params parameters
     * @return the bwt range [l,r]
     */
    template <typename out_t = no_init_vector<uint8_t>>
    out_t BWT_range(retrieve_params params = {}) const;

    /**
     * @brief reports the characters in the bwt in the range [l,r] (0 <= l <= r <= input size), else if l > r, then all
     * characters of the bwt are reported (default); the terminator is reported as TERMINATOR_SYMBOL
     * @tparam report_fnc_t type of the report function
     * @param report function that is called with every tuple (i,c) as a parameter, where i in [l,r] and c = L[i]; if num_threads = 1,
     * then the values are reported from left to right, if num_threads > 1, the order may vary
     * @param params parameters
     */
    template <typename report_fnc_t>
    void BWT_range(report_fnc_t report, retrieve_params params = {}) const;

    /**
     * @brief writes the characters in the bwt in the range [l,r] to the file file_name (0 <= l <= r <= input size), else
     * if l > r, then the whole bwt is written (default); the terminator is written as 0 (see terminator_as_zero), so
     * if the input contained 0, the output is not equal to the real bwt
     * @param file_name name of the file to write the bwt to
     * @param params parameters
     */
    template <typename out_sym_t = uint8_t>
    void BWT_range(std::string file_name, retrieve_params params = {}) const
    {
        adjust_retrieve_params(params, n - 1);

        retrieve_range<out_sym_t>([this](auto&& report, retrieve_params p) {
            BWT_range([&](uint64_t i, uint64_t c) { report(i, terminator_as_zero(c)); }, p);
        }, file_name, params);
    }

    // default number of pieces each thread of revert_strands reverts alternately: 2 with M_LF in EF, else 1
    static constexpr uint64_t default_revert_strands = mlf_enc == EF ? 2 : 1;

    // a piece of the input that revert_strands reverts from right to left
    struct revert_strand {
        uint64_t j; // the text position of the next symbol to report
        uint64_t first; // the text position of the first symbol of the piece
        uint64_t x; // the input interval of M_LF containing ISA[j + 1]
        uint64_t d; // ISA[j + 1] (M_LF in POS) or its offset in the x-th input interval (EF)
        uint64_t w; // [0..num_threads * num_strands - 1] the strand slot, unique among all threads
    };

    /**
     * @brief reverts pieces of the input with up to num_threads threads; each thread reverts num_strands pieces
     *        alternately, one LF-step per piece and turn; a finished piece is replaced by the next piece not taken by
     *        any thread
     * @tparam piece_fnc_t type of piece
     * @tparam report_fnc_t type of report
     * @param num_pieces the number of pieces
     * @param num_threads the number of threads
     * @param num_strands the number of pieces each thread reverts alternately
     * @param piece piece(k, s) sets s.j, s.first, s.x and s.d for the k-th piece, k in [0..num_pieces-1]
     * @param report report(w, j, c) is called with c = input[j] for each text position j of each piece, from right
     *        to left within a piece; w is the strand slot of the piece
     */
    template <typename piece_fnc_t, typename report_fnc_t>
    void revert_strands(uint64_t num_pieces, uint64_t num_threads, uint64_t num_strands, piece_fnc_t piece,
        report_fnc_t report) const;

    /**
     * @brief returns the input, reverted with up to num_threads threads (see revert_in_pieces)
     * @param num_threads maximum number of threads to use
     * @param num_strands the number of pieces each thread reverts alternately
     * @return the input
     */
    template <typename out_t = no_init_vector<uint8_t>>
    out_t revert(uint64_t num_threads = omp_get_max_threads(), uint64_t num_strands = default_revert_strands) const
    {
        out_t input;
        resize_for_overwrite(input, n - 1);
        revert_in_pieces([&](uint64_t, uint64_t i, uint64_t c) { input[i] = c; }, num_threads, num_strands);
        return input;
    }

    /**
     * @brief reports the characters of the input from right to left (single-threaded, one LF-query per character)
     * @tparam report_fnc_t type of the report function
     * @param report function that is called with every tuple (i,c) as a parameter, where c = input[i]
     */
    template <typename report_fnc_t>
    void revert(report_fnc_t report) const requires(std::invocable<report_fnc_t&, uint64_t, uint64_t>);

    /**
     * @brief reports the characters of the input with up to num_threads threads: the input is split at the text
     *        positions of the revert starts (see _revert_starts), and each thread reverts num_strands pieces
     *        alternately from right to left (see revert_strands)
     * @tparam report_fnc_t type of the report function
     * @param report function called as report(w, i, c) with c = input[i] and the strand slot w in
     *        [0..num_threads * num_strands - 1], by several threads at the same time
     * @param num_threads maximum number of threads to use
     * @param num_strands the number of pieces each thread reverts alternately
     */
    template <typename report_fnc_t>
    void revert_in_pieces(report_fnc_t report, uint64_t num_threads, uint64_t num_strands) const;

    /**
     * @brief reverts the input with up to num_threads threads and writes it to a file, sizeof(out_sym_t) bytes per
     *        symbol (see revert_in_pieces); each strand slot writes through its own buffer, and the buffers of a
     *        thread take huge_page_bytes bytes
     * @param file_name name of the file to write the reverted input to
     * @param num_threads maximum number of threads to use
     * @param num_strands the number of pieces each thread reverts alternately
     */
    template <typename out_sym_t = uint8_t>
    void revert(std::string file_name, uint64_t num_threads = omp_get_max_threads(),
        uint64_t num_strands = default_revert_strands) const;

    /**
     * @brief rebuilds and returns the suffix array in the range [l,r] (0 <= l <= r <= input size),
     * else if l > r, then the whole suffix array is rebuilt (default)
     * @param params parameters
     * @return the suffix array range [l,r]
     */
    no_init_vector<uint48_t> SA_range(retrieve_params params = {}) const requires(supports_multiple_locate);

    /**
     * @brief reports the suffix array values in the range [l,r] (0 <= l <= r <= input size), else if l > r, then the
     * whole suffix array is reported (default), with up to params.num_threads threads in an unspecified order (see
     * SA_range_in_pieces)
     * @tparam report_fnc_t type of the report function
     * @param report function that is called with every tuple (i,s) as a parameter, where i in [l,r] and s = SA[i]
     * @param params parameters
     */
    template <typename report_fnc_t>
    void SA_range(report_fnc_t report, retrieve_params params = {}) const requires(supports_multiple_locate)
    {
        SA_range_in_pieces([&](uint64_t, uint64_t i, uint64_t s) { report(i, s); }, params);
    }

    /**
     * @brief returns the number of strand slots SA_range_in_pieces uses with the parameters params, i.e., the number
     *        of its threads times params.num_strands
     * @param params parameters
     */
    uint64_t SA_range_strand_slots(retrieve_params params) const requires(supports_multiple_locate)
    {
        adjust_retrieve_params(params, n - 1);
        return SA_range_num_threads(params) * std::max<uint64_t>(1, params.num_strands);
    }

    /**
     * @brief reports the suffix array values in the range [l,r] (0 <= l <= r <= input size), else if l > r, then the
     *        whole suffix array (default), with up to params.num_threads threads: the range is split at BWT run starts
     *        into pieces, and each thread decodes params.num_strands pieces alternately from left to right, one
     *        Phi^{-1}-step (rlz: one decoded value) per piece and turn (see run_strands)
     * @tparam report_fnc_t type of the report function
     * @param report function called as report(w, i, s) with s = SA[i] for each i in [l,r] and the strand slot w in
     *        [0..SA_range_strand_slots(params) - 1], by several threads at the same time; within a piece, i increases
     * @param params parameters
     */
    template <typename report_fnc_t>
    void SA_range_in_pieces(report_fnc_t report, retrieve_params params = {}) const
        requires(supports_multiple_locate);

    /**
     * @brief writes the values in the suffix array of the input in the range [l,r] to the file file_name (0 <= l <= r <= input size),
     * else if l > r, then the whole suffix array is written (default); each strand slot (see SA_range_in_pieces) writes
     * through its own buffer
     * @param file_name name of the file to write the suffix array to
     * @param params parameters
     */
    void SA_range(std::string file_name, retrieve_params params = {}) const requires(supports_multiple_locate);

    // ############################# SERIALIZATION METHODS #############################

    /**
     * @brief reads the locate support of a serialized index without changing the position of the stream
     * @param in an input stream positioned at the beginning of a serialized index
     * @return the locate support the index was built with
     */
    static move_r_support peek_support(std::istream& in);

    /**
     * @brief stores the index to an output stream
     * @param out output stream to store the index to
     */
    void serialize(std::ostream& out) const;

    /**
     * @brief reads a serialized index from an input stream
     * @param in an input stream storing a serialized index
     */
    void load(std::istream& in);

    /**
     * @brief stores the index to an output stream
     * @param os output stream
     * @return the output stream
     */
    std::ostream& operator>>(std::ostream& os) const
    {
        serialize(os);
        return os;
    }

    /**
     * @brief reads a serialized index from an input stream
     * @param is input stream
     * @return the input stream
     */
    std::istream& operator<<(std::istream& is)
    {
        load(is);
        return is;
    }
};

#include "construction/construction.hpp"
#include "move_r.tpp"
#include "queries.tpp"
// extern template declarations of the index types (has to be included last)
#include "construction/instantiations.hpp"

#endif // MOVE_R_MOVE_R_MOVE_R_HPP
