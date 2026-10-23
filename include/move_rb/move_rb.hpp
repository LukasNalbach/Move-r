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

#ifndef MOVE_R_MOVE_RB_MOVE_RB_HPP
#define MOVE_R_MOVE_RB_MOVE_RB_HPP

#include <misc/utils.hpp>
#include <misc/portable_layout.hpp>

#include <cctype>
#include <cstdlib>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <tsl/sparse_set.h>

#include <data_structures/context_tables.hpp>
#include <data_structures/elias_fano_bit_vector.hpp>
#include <move_r/move_r.hpp>
// search schemes, apm_params_t and query_support_t
#include <misc/apm.hpp>
#include <data_structures/fasta_sequence_data.hpp>

// approximate pattern matching with the Hamming distance on a bidirectional index idx_t
template <typename idx_t> class apm_hamming;

// approximate pattern matching with the edit distance on a bidirectional index idx_t
template <typename idx_t> class apm_edit;

/**
 * @brief returns the name of the file an index is stored in: <base>.move-rb<support suffix>
 * @param base base name (usually the name of the input file)
 * @param support type of locate support (see move_r_support_file_suffix)
 */
inline std::string move_rb_index_file_name(const std::string& base, move_r_support support)
{
    return base + ".move-rb" + std::string(move_r_support_file_suffix(support));
}

/**
 * @brief reads the type of locate support from the header of a serialized move_rb and resets the stream to the
 *        beginning of the index (see move_rb::serialize)
 * @param in input stream positioned at the beginning of a serialized move_rb
 * @return the type of locate support the index was built with
 */
inline move_r_support read_move_rb_support(std::istream& in)
{
    move_r_support support;
    std::streampos beg = in.tellg();
    in.read((char*) &support, sizeof(move_r_support));
    in.seekg(beg);
    return support;
}

// output of a locate query: the occurrences (_occ), or each occurrence with its suffix array position (_occ_and_sa_pos)
enum locate_output_t : uint8_t { _occ = 0, _occ_and_sa_pos = 1 };

/**
 * @brief bidirectional move-r index of size O(r * (a / (a - 1))) words: a move_r index over T (idx_fwd) and one over
 *        T^R (idx_bwd). Template arguments and members of a direction refer to T for LEFT and *_fwd, and to T^R for
 *        RIGHT and *_bwd. A search keeps the suffix array interval of the pattern in T and of the reversed pattern in
 *        T^R (see search_context_t).
 * @tparam support type of locate support (_count, _locate_one, _locate_move or _locate_rlzsa)
 */
template <move_r_support support = _locate_move>
class move_rb {
    public:

    template <typename> friend class apm_hamming;
    template <typename> friend class apm_edit;

    // ############################# TYPES, CONSTANTS AND PARAMETERS #############################

    // true <=> the positions of the SA-samples of the backward index in the forward suffix array are stored
    // (see _SA_pos); else, they are computed with Phi^{-1} (see compute_center)
    static constexpr bool stores_SA_pos = support == _locate_rlzsa;

    // type of idx_fwd; its M_LF uses the POS encoding (absolute input interval starting positions)
    using move_r_fwd_t = constexpr_switch_t<
        constexpr_case<support == _locate_move,  move_r<_locate_move_bi, POS>>,
        constexpr_case<support == _locate_rlzsa, move_r<_locate_rlzsa_bi, POS>>,
        constexpr_case<support == _locate_one,   move_r<_locate_one_bi, POS>>,
     /* constexpr_case<support == _count, */     move_r<_count_bi, POS>>;

    // type of idx_bwd; it stores M_LF and the SA-samples
    using move_r_bwd_t = std::conditional_t<support == _count, move_r<_count_bi, POS>, move_r<_locate_one_bi, POS>>;

    // true <=> the index supports locate
    static constexpr bool supports_locate = move_r_fwd_t::supports_locate;
    // true <=> the index supports locating all occurrences
    static constexpr bool supports_multiple_locate = move_r_fwd_t::supports_multiple_locate;
    static constexpr bool has_SA_rlz = move_r_fwd_t::has_SA_rlz; // true <=> the index has an rlz of SA^d
    static constexpr bool has_locate_move = move_r_fwd_t::has_locate_move; // true <=> Phi and Phi^{-1} use move data structures

    // true <=> M_Phi and M_Phi^{-1} of the forward index are stored in the index
    static constexpr bool has_mphi = move_r_fwd_t::has_mphi;

    // true <=> _SA_pos is stored
    static constexpr bool has_SA_pos = supports_locate && stores_SA_pos;

    using inp_t = std::vector<uint8_t>; // pattern type of the approximate pattern matching algorithms

    // every input_interval_sample_rate-th input interval starting position is stored in the sd-arrays S_*_p
    // (see sample_input_intervals)
    static constexpr uint64_t input_interval_sample_rate = 4;

    // maximum number of input intervals of M_LF that input_interval_near scans before it uses a predecessor query
    static constexpr uint64_t max_scanned_intervals = 16;

    // num_smaller_wavelet scans the range of input intervals instead of querying the wavelet tree over L', if
    // the range has at most this many input intervals per level of the wavelet tree
    static constexpr uint64_t max_scanned_intervals_per_wavelet_level = 96;

    // type of _SA_pos
    using sa_pos_t = std::conditional_t<has_SA_pos, interleaved_bit_aligned_vectors, std::monostate>;

    // output of a locate query (see ::locate_output_t)
    using locate_output_t = ::locate_output_t;
    static constexpr locate_output_t _occ = ::_occ;
    static constexpr locate_output_t _occ_and_sa_pos = ::_occ_and_sa_pos;

    // contexts of the queries (see search_context.hpp, locate_context.hpp and extend_context.hpp)
    template <locate_output_t output = _occ> struct locate_context_t;
    template <query_support_t query_support> struct search_context_t;
    template <query_support_t query_support> struct extend_context_t;

protected:
    // whether an SA-sample is at the start (BEG) or at the end (END) of its BWT run
    enum sample_t : uint8_t {BEG = 0, END = 1};

    /**
     * @brief the SA-sample of a search context, an occurrence of the pattern with its position in [b, e] (or a
     *        way to compute it):
     * - d == LEFT: the sample is SA_run_start(x) (t == BEG) resp. SA_run_end(x) (t == END) of idx_fwd minus i;
     *   its position is b + o
     * - d == RIGHT: the sample is n - v - (m - i + 1), where v is SA_run_start(x) resp. SA_run_end(x) of idx_bwd;
     *   its position is reached from SA_pos with m - i LF-steps (see compute_center)
     * - d == NO_DIR: the sample is x; its position is b + o (only kept if stores_SA_pos)
     * - d == BOTH_DIRS: the sample has not been computed yet; the interval contains an SA-sample of the index,
     *   which sample_at_run_boundary computes (the state of a context loaded from the context tables)
     *
     * o, x and i store values up to n (48 bits, see uint48_t; i is at most m). The fields are bit-fields in
     * the layout of a packed struct of GCC (see packed_bit_fields): d and t share the first byte, o, i and x take 6
     * bytes each.
     */
    struct sample_info_pack_t {
        using d_field = next_bit_field_t<void, std::bit_width(uint8_t(BOTH_DIRS))>;
        using t_field = next_bit_field_t<d_field, 8 * sizeof(direction_t) - d_field::width>;
        using o_field = next_bit_field_t<t_field, 8 * uint48_t::width>;
        using i_field = next_bit_field_t<o_field, 8 * uint48_t::width>;
        using x_field = next_bit_field_t<i_field, 8 * uint48_t::width>;

        packed_bit_fields<bytes_of_bit_fields<x_field>> bits;

        sample_info_pack_t() = default;

        sample_info_pack_t(direction_t d, sample_t t, uint64_t o, uint64_t i, uint64_t x)
        { bits.template assign<d_field, t_field, o_field, i_field, x_field>(d, t, o, i, x); }

        // the index the sample is computed with (LEFT: idx_fwd, RIGHT: idx_bwd), or NO_DIR, BOTH_DIRS
        direction_t d() const { return direction_t(bits.template get<d_field>()); }
        void set_d(direction_t value) { bits.template set<d_field>(value); }

        // BEG/END <=> the SA-sample of the index is at the start/end of its BWT run
        sample_t t() const { return sample_t(bits.template get<t_field>()); }
        void set_t(sample_t value) { bits.template set<t_field>(value); }

        // offset of the position of the sample from b
        uint64_t o() const { return bits.template get<o_field>(); }
        void set_o(uint64_t value) { bits.template set<o_field>(value); }

        // number of extensions in direction d since x was set (unused, if d == NO_DIR)
        uint64_t i() const { return bits.template get<i_field>(); }
        void set_i(uint64_t value) { bits.template set<i_field>(value); }

        // index of the input interval of M_LF that starts (t == BEG) resp. ends (t == END) the BWT run with the
        // SA-sample of the index (the value of the sample, if d == NO_DIR)
        uint64_t x() const { return bits.template get<x_field>(); }
        void set_x(uint64_t value) { bits.template set<x_field>(value); }
    };

    static_assert(sizeof(sample_info_pack_t) == sizeof(direction_t) + 3 * uint48_t::width);

    // a symbol c in a range [b_, e_] of input intervals of M_LF (see extend_context_t)
    struct symbol_runs_t {
        uint64_t sym; // the symbol c
        uint64_t first; // the first input interval x in [b_, e_] with L'[x] = c
        uint64_t last; // the last input interval x in [b_, e_] with L'[x] = c
    };

    /* the dir-interval [b, e] of a search context restricted to the input intervals b_..e_ of M_LF with the
       extension symbol; LF maps [b, e] to the dir-interval of the extended context */
    struct prime_tuple_t { uint64_t b, e, b_, e_; };

    // empty member type; the tag makes the types of the [[no_unique_address]] members distinct
    template <uint8_t tag> struct absent_t { };

    // indices of the fields of an SA-sample (see sample_info_pack_t) in the samples of the context tables; the
    // tables store the sample of a context only if its interval contains no SA-sample of the index (see
    // sample_at_run_boundary)
    static constexpr uint8_t _smpl_d = 0, _smpl_t = 1, _smpl_o = 2, _smpl_i = 3, _smpl_x = 4;

    // a pattern of a level of the context tables during their construction (see build_context_tables)
    struct context_table_row_t {
        uint64_t b, e; // suffix array interval of the pattern
        uint64_t b_R; // left limit of the suffix array interval of the reversed pattern
        uint32_t parent; // entry (in the previous level) of the pattern it was extended from
        bool stores_sample; // true <=> its interval contains no SA-sample of the index, so s is stored
        STRING_UTILS_NO_UNIQUE_ADDRESS std::conditional_t<supports_locate, sample_info_pack_t, absent_t<3>> s;
    };

    // marks a non-existing entry of the context tables (see context_tables::no_entry)
    static constexpr uint32_t no_ctx_entry = ::context_tables::no_entry;

    // direction, by which the context tables group their entries; they also group them by ctx_secondary_dir (see
    // context_tables)
    static constexpr direction_t ctx_primary_dir = LEFT;
    static constexpr direction_t ctx_secondary_dir = RIGHT;
    // keys of the context tables used for extensions in direction dir (see context_extension_t)
    template <direction_t dir> static constexpr context_extension_t ctx_key =
        dir == ctx_primary_dir ? _primary_extension : _secondary_extension;
    // order of the context tables, in which the entries of the dir-extensions of a pattern are contiguous
    // (see symbol_major_keys)
    template <direction_t dir> static constexpr context_order_t ctx_group_order =
        ctx_key<dir> == _primary_extension ? _reverse_order : _forward_order;

    // the state of the construction shared by its steps (see build)
    struct build_state_t;

    // ############################# DATA MEMBERS #############################

    uint64_t n = 0; // the length of T, including the terminator
    uint64_t sigma = 0; // the number of distinct symbols in T, including the terminator
    // number of texts (T, T^R), for which the construction read the pfp output from a cache (not serialized)
    uint8_t _pfp_cached_texts = 0;

    move_r_fwd_t idx_fwd; // move_r index of T
    move_r_bwd_t idx_bwd; // move_r index of T^R

    /* samples the input interval starting positions of M_LF of idx_fwd (resp. idx_bwd): an sd-array over
       [0..k], k = ceil(r' / input_interval_sample_rate), whose i-th value is M_LF.p(i * input_interval_sample_rate)
       for i < k and n for i = k; universe n + 1 (see sample_input_intervals and input_interval) */
    elias_fano_bit_vector _S_MLF_p_fwd;
    elias_fano_bit_vector _S_MLF_p_bwd;

    // samples the input interval starting positions of M_Phi resp. M_Phi^{-1} of idx_fwd (as _S_MLF_p_fwd); empty
    // unless has_mphi
    elias_fano_bit_vector _S_MPhi_p_fwd;
    elias_fano_bit_vector _S_MPhi_inv_p_fwd;

    /* the positions in the suffix array of T of the SA-samples of idx_bwd: [0..s-1], s = idx_bwd.SA_run_pos().num_samples();
       _SA_pos[k] = j with SA[j] = n - 1 - v_k, where v_k is the value of the k-th SA-sample of idx_bwd (numbered as
       in run_sample_positions: one at each start of a BWT run, one at each end of a BWT run of length > 1);
       width bit_width(n) */
    STRING_UTILS_NO_UNIQUE_ADDRESS sa_pos_t _SA_pos;

    // the start positions, names and separator symbol of the sequences of a multi-sequence (FASTA) input (see
    // fasta_sequence_data); empty unless the index was built from such an input
    fasta_sequence_data _seq_data;

    // the search contexts of all patterns up to some length (see context_tables); empty unless requested (see
    // move_r_params::context_table_share)
    ::context_tables _context_tables;

    // ############################# INTERNAL METHODS #############################

    // serializes a data structure, if the index has it (T != std::monostate)
    template <typename T> static void opt_serialize(const T& s, std::ostream& out)
    {
        if constexpr (!std::is_same_v<T, std::monostate>) s.serialize(out);
    }

    // loads a data structure, if the index has it (T != std::monostate)
    template <typename T> static void opt_load(T& s, std::istream& in)
    {
        if constexpr (!std::is_same_v<T, std::monostate>) s.load(in);
    }

    // returns the size of a data structure in bytes (0, if the index does not have it)
    template <typename T>
    static uint64_t opt_size(const T& s)
    {
        if constexpr (std::is_same_v<T, std::monostate>) return 0; else return s.size_in_bytes();
    }

    /**
     * @brief returns the position in the suffix array of T of an SA-sample of idx_bwd (see _SA_pos)
     * @param x index of an input interval of M_LF of idx_bwd starting (at_start = true) resp. ending
     *          (at_start = false) a BWT run
     * @param at_start true <=> the sample at the start of the run is used (else the sample at its end)
     */
    inline uint64_t SA_pos(uint64_t x, bool at_start) const requires(has_SA_pos)
    {
        const run_sample_positions& where = idx_bwd.SA_run_pos();
        return _SA_pos.template get<0, uint64_t>(at_start ? where.start_sample_index(x) : where.end_sample_index(x));
    }

    // returns idx_fwd (dir == LEFT) resp. idx_bwd (dir == RIGHT)
    template <direction_t dir>
    const auto& index() const
    {
        if constexpr (dir == LEFT) return idx_fwd; else return idx_bwd;
    }

    // returns the sd-array sampling the input interval starting positions of M_LF (see _S_MLF_p_fwd)
    template <direction_t dir>
    const elias_fano_bit_vector& S_MLF_p() const
    {
        if constexpr (dir == LEFT) return _S_MLF_p_fwd; else return _S_MLF_p_bwd;
    }

    // returns PLCP[Phi^{-1}(pos)], i.e., LCP[j + 1] for SA[j] = pos, where x is the index of the input interval of
    // M_Phi^{-1} containing pos
    inline uint64_t plcp_at_phi_inv(uint64_t pos, uint64_t x) const requires(has_mphi)
    {
        const auto& M = idx_fwd.M_Phi_inv();
        return M.template row<0>(x) - (pos - M.p(x));
    }

    /**
     * @brief returns the index of the input interval containing position i
     * @tparam fnc_t type of the function interval_start
     * @param i position in [0, n)
     * @param sd_arr sd-array sampling the input interval starting positions of a move data structure (see
     *        sample_input_intervals)
     * @param interval_start function returning the starting position of the x-th input interval of the same
     *        move data structure
     */
    template <typename fnc_t>
    inline static uint64_t input_interval(uint64_t i, const elias_fano_bit_vector& sd_arr, fnc_t interval_start);

    // returns the index of the input interval of M_LF containing i (see input_interval_near)
    template <direction_t dir>
    inline uint64_t mlf_interval(uint64_t i) const
    {
        const auto& M_LF = index<dir>().M_LF();
        return input_interval(i, S_MLF_p<dir>(), [&](uint64_t y) { return M_LF.p(y); });
    }

    // returns the index of the input interval of M_Phi containing pos
    inline uint64_t phi_interval(uint64_t pos) const requires(has_mphi)
    {
        const auto& M_Phi = idx_fwd.M_Phi();
        return input_interval(pos, _S_MPhi_p_fwd, [&](uint64_t y) { return M_Phi.p(y); });
    }

    // returns the index of the input interval of M_Phi^{-1} containing pos
    inline uint64_t phi_inv_interval(uint64_t pos) const requires(has_mphi)
    {
        const auto& M_Phi_inv = idx_fwd.M_Phi_inv();
        return input_interval(pos, _S_MPhi_inv_p_fwd, [&](uint64_t y) { return M_Phi_inv.p(y); });
    }

    // ############################# CONSTRUCTION #############################

    /**
     * @brief builds the index (see build_state_t for its steps)
     * @param input the input; it is mapped to the effective alphabet in place and restored afterwards
     * @param params construction parameters
     * @param log_summary true <=> log a summary of the construction (if params.log is set)
     */
    void build(move_r_input_view& input, move_r_params params = {}, bool log_summary = true);

    // the steps of build, in the order of their calls (see build_state_t)

    // reads the input file into memory or writes the input to a temporary file for the prefix-free parsing
    void prepare_input(build_state_t& st);
    // maps the bytes of the input to the effective alphabet (byte alphabet)
    void build_byte_alphabet(build_state_t& st);
    // maps the symbols of the input to the effective alphabet (integer alphabet)
    void build_integer_alphabet(build_state_t& st);
    // builds idx_fwd and idx_bwd; idx_fwd is stored on disk meanwhile in the modes _pfp and _suffix_array_space
    void build_sub_indexes(build_state_t& st);
    // builds the sd-arrays sampling the input interval starting positions of the move data structures of a direction
    template <direction_t dir> void build_samplings(build_state_t& st);
    // loads idx_fwd from disk, if it is stored there
    void load_forward_index(build_state_t& st);
    // builds _SA_pos
    void build_SA_pos(build_state_t& st);
    // maps the input back from the effective alphabet
    void restore_input(build_state_t& st);
    // builds the context tables; stores the SA-samples and M_Phi on disk meanwhile, if they are not read
    void build_context_tables_phase(build_state_t& st);
    /**
     * @brief logs and measures building an sd-array sampling the input interval starting positions of M
     * @param st the state of the construction
     * @param name name of the phase
     * @param key measurement key of the phase
     * @param M a move data structure
     * @param sampling output
     */
    template <typename mds_t>
    void build_sampling_phase(build_state_t& st, const std::string& name, const char* key, const mds_t& M,
        elias_fano_bit_vector& sampling);

    /**
     * @brief reverses T in memory
     * @param T the input
     * @param p number of threads to use
     * @param log true <=> print log messages
     * @param measurement stream to write the time of the phase to (or nullptr)
     * @param measurement_key measurement key of the phase
     */
    void reverse(move_r_input_view& T, uint64_t p, bool log = false,
        std::ostream* measurement = nullptr, const std::string& measurement_key = "");

    /**
     * @brief returns the sd-array sampling the input interval starting positions of a move data structure M:
     *        its i-th value is M.p(i * input_interval_sample_rate) for i < k = ceil(M.num_intervals() /
     *        input_interval_sample_rate) and n for i = k; universe n + 1
     * @param M a move data structure over [0, n)
     * @param p number of threads to use
     */
    template <typename mds_t>
    elias_fano_bit_vector sample_input_intervals(const mds_t& M, uint64_t p = 1) const
    {
        const uint64_t num_samples = div_ceil(M.num_intervals(), input_interval_sample_rate);

        return elias_fano_bit_vector(
            [&](uint64_t i) { return uint64_t(i < num_samples ? M.p(i * input_interval_sample_rate) : n); },
            num_samples + 1, n + 1, p);
    }

    // ############################# CONTEXT TABLES #############################

    /**
     * @brief sets a search context to an entry of the context tables; the values not stored in the entry
     *        are computed (see context_tables)
     * @param ctx the context to set
     * @param len length of the pattern (the level of the entry)
     * @param idx entry in the forward order
     * @param idx_R entry in the reverse order (unused, if the index has no backward index)
     */
    template <query_support_t query_support>
    void load_context(search_context_t<query_support>& ctx, uint64_t len, uint32_t idx, uint32_t idx_R) const;

    /**
     * @brief extends each pattern of level len - 1 of the context tables in direction dir with every symbol; the
     *        patterns are processed in parallel in rounds of round_parents patterns, each round in up to p chunks
     * @tparam dir extend direction
     * @param len the level of the extended patterns
     * @param num_parents the number of patterns of level len - 1
     * @param parent_perm [0..num_parents-1] the entry in reverse order of each entry in forward order of level
     *        len - 1
     * @param round_parents number of patterns per round
     * @param p number of threads to use
     * @param start_chunk function called as start_chunk(c) before chunk c of a round is processed
     * @param report function called as report(c, parent, sym, child) for each extension child of the pattern with
     *        entry parent in level len - 1 by the symbol sym in chunk c
     * @param end_round function called as end_round(num_chunks) after each round; if it returns false, no
     *        further round is processed
     */
    template <direction_t dir, typename start_t, typename report_t, typename end_t>
    void extend_level_parents(uint64_t len, uint64_t num_parents, const no_init_vector<uint32_t>& parent_perm,
        uint64_t round_parents, uint64_t p, start_t start_chunk, report_t report, end_t end_round) const;

    /**
     * @brief stores the SA-samples of the entries of a level of the context tables with the bit widths of their
     *        largest values (see _smpl_d)
     * @param lvl the level
     * @param samples the SA-samples in the order of the entries that store one
     */
    static void build_level_samples(::context_tables::level_t& lvl, const std::vector<sample_info_pack_t>& samples);

    /**
     * @brief builds the reverse order of a level of the context tables, i.e., its entries sorted by the left
     *        limits of their intervals in the suffix array of T^R
     * @param lvl the level; its forward order must be built
     * @param beg_R [0..num_entries-1] the left limit of the interval in the suffix array of T^R of each entry in forward
     *        order (deleted)
     * @param width [0..num_entries-1] e - b of each entry in forward order (deleted)
     * @param primary_keys [0..num_entries-1] the key sym * num_parents + parent of each entry in forward order
     * @param secondary_keys [0..num_entries-1] the key parent * sigma + sym of each entry in forward order, or empty
     * @param parent_perm see extend_level_parents
     * @param num_parents the number of entries of the previous level
     * @param packed_rev [0..num_entries-1] the packed reversed pattern of each entry in forward order, or empty
     * @param perm output, [0..num_entries-1] the entry in reverse order of each entry in forward order
     * @param packed_rev_sorted output, packed_rev in reverse order
     * @param p number of threads to use
     */
    void build_reverse_order(::context_tables::level_t& lvl, std::vector<uint64_t>& beg_R,
        std::vector<uint64_t>& width, const std::vector<uint64_t>& primary_keys, const std::vector<uint64_t>& secondary_keys,
        const no_init_vector<uint32_t>& parent_perm, uint64_t num_parents, const std::vector<uint64_t>& packed_rev,
        no_init_vector<uint32_t>& perm, std::vector<uint64_t>& packed_rev_sorted, uint64_t p) const;

    // ############################# CONSTRUCTORS #############################

public:
    move_rb() = default;

    /**
     * @brief constructs a bidirectional index of the input
     * @tparam inp_t_ container type of the input; any contiguous container of an integer symbol type of 1, 2,
     *         4 or 8 bytes (e.g., std::string, std::vector<uint8_t>, std::vector<int8_t>)
     * @param input the input
     * @param params construction parameters
     * @param log_summary true <=> log a summary of the construction (if params.log is set)
     */
    template <typename inp_t_>
    move_rb(inp_t_&& input, move_r_params params = {}, bool log_summary = true)
        : move_rb(input, std::move(params), log_summary) {}

    // constructs a bidirectional index of a NUL-terminated string
    move_rb(const char* input, move_r_params params = {}, bool log_summary = true)
        : move_rb(std::string(input), std::move(params), log_summary) {}

    /**
     * @brief constructs a bidirectional index of the input; if params.input_file is set, the input is read
     *        from that file (only for one-byte symbols)
     * @tparam inp_t_ container type of the input (see above)
     * @param input the input
     * @param params construction parameters
     * @param log_summary true <=> log a summary of the construction (if params.log is set)
     */
    template <typename inp_t_>
    move_rb(inp_t_& input, move_r_params params = {}, bool log_summary = true);

    // ############################# MISC PUBLIC METHODS #############################

    // returns the length of the input (without the terminator symbol)
    inline uint64_t input_size() const { return n - 1; }

    // returns the number of runs in the BWT of the input
    inline uint64_t num_bwt_runs() const { return idx_fwd.num_bwt_runs(); }

    // returns the number of distinct symbols in the input (without the terminator symbol)
    inline uint64_t alphabet_size() const { return sigma - 1; }

    // returns whether the index stores PSB_L' (see move_r::PSB_L_), which is used to enumerate all extensions
    inline bool has_PSB_L_() const { return !idx_fwd.PSB_L_().empty(); }

    // true <=> the index stores context tables (see context_tables)
    inline bool has_context_tables() const { return !_context_tables.empty(); }

    // returns the maximum pattern length, for which the contexts of all patterns are stored (0 <=> none are stored)
    inline uint64_t context_tables_max_length() const { return _context_tables.max_length(); }

    // returns the pattern length of the single-length table (0 if there is none, see context_tables)
    inline uint64_t context_tables_single_length() const { return _context_tables.single_length(); }

    // returns " (forward)" resp. " (backward)" for log messages
    template <direction_t dir>
    static constexpr const char* direction_suffix()
    { return dir == LEFT ? " (forward)" : " (backward)"; }

    /**
     * @brief attaches multi-sequence (FASTA) metadata to the index (as produced by process_fasta) for
     *        sequence-relative coordinates and names in the SAM output; call it after the construction
     * @param seq_data the multi-sequence metadata
     */
    inline void set_fasta_sequence_data(fasta_sequence_data seq_data) { _seq_data = std::move(seq_data); }

    // deletes the context tables
    inline void clear_context_tables() { _context_tables.clear(); }

    // returns the size of the index in bytes
    uint64_t size_in_bytes() const;

    /**
     * @brief logs the construction time, the peak memory usage and the construction throughput
     * @param time_construction construction time in nanoseconds
     * @param peak_mem_usage not used (the overall peak is logged)
     * @param mf measurement file to write the throughput to (or nullptr)
     */
    void log_construction_summary(uint64_t time_construction, uint64_t peak_mem_usage,
        std::ostream* mf = nullptr) const;

    /**
     * @brief logs the sizes of the data structures of one direction (the move_r index and the samplings of the
     *        input intervals) indented under the heading "forward index:" or "backward index:"
     * @tparam dir the direction
     */
    template <direction_t dir>
    void log_direction_sizes() const;

    // logs the size of the context tables, if the index has any (see log_component)
    void log_context_tables_size() const
    {
        if (!_context_tables.empty()) log_component("context tables", _context_tables.size_in_bytes());
    }

    // logs the lengths of the built context tables and the number of bits per context
    void log_built_context_tables() const;

    /**
     * @brief logs the sizes of the data structures with their shares of the index size, one per line (see
     *        log_component), indented under the headings "forward index:" and "backward index:", followed by the
     *        data structures of the whole index; the list ends with the size of the rest of the index ("other")
     */
    void log_data_structure_sizes() const;

    // logs the contents of all data structures (only for n <= 50)
    void log_data_structures() const;

    /**
     * @brief logs the context tables (see context_tables::log_data_structures), if the index has any
     * @param max_depth the patterns up to this length are logged in the tree
     * @param max_entries at most this many patterns are logged per tree or table
     */
    void log_context_tables(uint64_t max_depth = default_logged_tree_depth,
        uint64_t max_entries = default_logged_entries) const
    {
        if (_context_tables.empty()) {
            std::cout << "no context tables" << std::endl;
            return;
        }

        _context_tables.log_data_structures([&](uint64_t sym) { return idx_fwd.logged_symbol(sym); },
            ctx_primary_dir == LEFT, max_depth, max_entries);
    }

    /**
     * @brief adds the parameters the size of the index depends on to the facts of the info modes (see
     *        move_r::add_parameters; the quantities of the index of the reversed text get the suffix ^R)
     * @param facts the facts
     */
    void add_parameters(index_facts& facts) const;

    /**
     * @brief writes the sizes of the data structures as key=value pairs to the output stream out
     * @param out an output stream
     */
    void log_data_structure_sizes(std::ostream& out) const;

    // ############################# PUBLIC ACCESS METHODS #############################

    // returns idx_fwd, the move_r index of T
    inline const move_r_fwd_t& forward_index() const { return idx_fwd; }

    // returns idx_bwd, the move_r index of T^R
    inline const move_r_bwd_t& backward_index() const { return idx_bwd; }

    // maps a symbol to its symbol in the effective alphabet (see move_r::symbol_idx)
    template <typename sym_t>
    inline uint64_t map_symbol(sym_t sym) const { return idx_fwd.map_symbol(sym); }

    // maps a symbol of the effective alphabet to its symbol in the input
    inline uint64_t unmap_symbol(uint64_t sym) const { return idx_fwd.unmap_symbol(sym); }

    /**
     * @brief returns the start positions, names and separator symbol of the sequences of a multi-sequence
     *        (FASTA) input; empty (has_sequences() == false) unless the index was built from such an input
     */
    inline const fasta_sequence_data& seq_data() const { return _seq_data; }

    // returns the context tables
    inline const ::context_tables& context_tables() const { return _context_tables; }

    /**
     * @brief returns the index of the input interval of M_LF containing position i; searches at most
     *        max_scanned_intervals input intervals to the left and to the right of the x-th one, then uses a
     *        predecessor query
     * @tparam dir the direction
     * @param i [0..n-1] a position in the suffix array
     * @param x [0..r'-1] index of an input interval of M_LF near position i
     */
    template <direction_t dir>
    inline uint64_t input_interval_near(uint64_t i, uint64_t x) const;

    // ############################# QUERY METHODS #############################

    /**
     * @brief returns the search context of the empty pattern
     * @tparam query_support whether the context supports COUNT or LOCATE queries
     */
    template <query_support_t query_support>
    inline search_context_t<query_support> empty_context() const { return search_context_t<query_support>(*this); }

    /**
     * @brief searches a pattern from right to left; the search starts at a context stored in the context
     *        tables, if there is one (see context_of)
     * @tparam query_support whether the returned context supports COUNT or LOCATE queries
     * @param P the pattern
     * @param occurs is set to whether P occurs in the input
     * @return the search context of P (of the longest suffix of P that occurs in the input, if P does not occur)
     */
    template <query_support_t query_support = LOCATE, typename pat_t>
        requires(!move_r_is_pattern_view<std::remove_cvref_t<pat_t>>)
    search_context_t<query_support> search(const pat_t& P, bool& occurs) const
    { return search<query_support>(move_r_pattern_view_t<pat_t>(P), occurs); }

    template <query_support_t query_support = LOCATE, typename sym_t>
    search_context_t<query_support> search(move_r_pattern_view<sym_t> P, bool& occurs) const;

    /**
     * @brief returns the search context of the last single_length() symbols of P from the single-length
     *        context table (see context_tables::single_len); a search of P can start at this context
     * @tparam query_support whether the returned context supports COUNT or LOCATE queries
     * @param P the pattern
     * @param matched is set to the length of the pattern of the returned context; 0, if there is no
     *        single-length table, P is shorter than single_length() or the symbols do not occur in the input
     * @return the search context (the context of the empty pattern, if matched = 0)
     */
    template <query_support_t query_support, typename pat_t>
        requires(!move_r_is_pattern_view<std::remove_cvref_t<pat_t>>)
    search_context_t<query_support> context_of(const pat_t& P, uint64_t& matched) const
    { return context_of<query_support>(move_r_pattern_view_t<pat_t>(P), matched); }

    template <query_support_t query_support, typename sym_t>
    search_context_t<query_support> context_of(move_r_pattern_view<sym_t> P, uint64_t& matched) const;

    // counts the occurrences of a NUL-terminated pattern with at most k mismatches (Hamming distance)
    uint64_t count_hamming_dist(const char* P, const search_scheme_t& scheme) const
    { return count_hamming_dist(std::string_view(P), scheme); }

    /**
     * @brief counts the occurrences of a pattern with at most k mismatches (Hamming distance)
     * @tparam sym_t symbol type of the pattern
     * @param P the pattern to search
     * @param scheme the search scheme to use (provides k)
     * @param params parameters of the search (see apm_params_t)
     * @return number of occurrences of P in T with at most k mismatches
     */
    template <typename sym_t>
    uint64_t count_hamming_dist(move_r_pattern_view<sym_t> P, const search_scheme_t& scheme,
        apm_params_t params = {}) const;

    template <typename pat_t>
        requires(!move_r_is_pattern_view<std::remove_cvref_t<pat_t>>)
    uint64_t count_hamming_dist(const pat_t& P, const search_scheme_t& scheme, apm_params_t params = {}) const
    { return count_hamming_dist(move_r_pattern_view_t<pat_t>(P), scheme, params); }

    /**
     * @brief counts the occurrences of a pattern with at most k mismatches (Hamming distance), by number of
     *        mismatches
     * @tparam sym_t symbol type of the pattern
     * @param P the pattern to search
     * @param scheme the search scheme to use (provides k)
     * @param params parameters of the search (see apm_params_t)
     * @return [0..k] h[e] is the number of occurrences of P in T with exactly e mismatches
     */
    template <typename sym_t>
    std::vector<uint64_t> count_hamming_dist_histogram(move_r_pattern_view<sym_t> P,
        const search_scheme_t& scheme, apm_params_t params = {}) const;

    template <typename pat_t>
        requires(!move_r_is_pattern_view<std::remove_cvref_t<pat_t>>)
    std::vector<uint64_t> count_hamming_dist_histogram(const pat_t& P, const search_scheme_t& scheme,
        apm_params_t params = {}) const
    { return count_hamming_dist_histogram(move_r_pattern_view_t<pat_t>(P), scheme, params); }

    /**
     * @brief locates a pattern with at most k errors (w.r.t. dist_metr)
     * @tparam mode whether the occurrences store a CIGAR alignment
     * @param P the pattern to search (as bytes)
     * @param scheme the search scheme to use (provides k)
     * @param dist_metr distance metric (HAMMING_DISTANCE or EDIT_DISTANCE)
     * @param occ vector, to which the occurrences of P in T with at most k errors are appended (it is not
     *        cleared)
     * @param params parameters of the search (see apm_params_t)
     * @return in CIGAR mode, the CIGARs of P against its matched strings, one per suffix array interval
     *         (occ.cig_idx is an index in this vector); empty in NO_CIGAR mode
     */
    template <cigar_mode_t mode = NO_CIGAR>
    std::vector<cigar_t<uint8_t>> locate_apm(const std::vector<uint8_t>& P, const search_scheme_t& scheme,
        distance_metric_t dist_metr, std::vector<aprx_occ_t<mode>>& occ, apm_params_t params = {}) const
        requires(supports_multiple_locate);

    // converts the pattern to bytes and calls locate_apm
    template <cigar_mode_t mode = NO_CIGAR, typename sym_t>
    std::vector<cigar_t<uint8_t>> locate(move_r_pattern_view<sym_t> P, const search_scheme_t& scheme,
        distance_metric_t dist_metr, std::vector<aprx_occ_t<mode>>& occ, apm_params_t params = {}) const
        requires(supports_multiple_locate)
    {
        std::vector<uint8_t> P_(P.begin(), P.end());
        return locate_apm<mode>(P_, scheme, dist_metr, occ, params);
    }

    template <cigar_mode_t mode = NO_CIGAR, typename pat_t>
        requires(supports_multiple_locate
            && !move_r_is_pattern_view<std::remove_cvref_t<pat_t>>)
    std::vector<cigar_t<uint8_t>> locate(const pat_t& P, const search_scheme_t& scheme,
        distance_metric_t dist_metr, std::vector<aprx_occ_t<mode>>& occ, apm_params_t params = {}) const
    { return locate<mode>(move_r_pattern_view_t<pat_t>(P), scheme, dist_metr, occ, params); }

    /**
     * @brief locates a pattern with at most k errors (w.r.t. dist_metr)
     * @tparam sym_t symbol type of the pattern
     * @param P the pattern to search
     * @param scheme the search scheme to use (provides k)
     * @param dist_metr distance metric (HAMMING_DISTANCE or EDIT_DISTANCE)
     * @param params parameters of the search (see apm_params_t)
     * @return occurrences of P in T with at most k errors (w.r.t. dist_metr)
     */
    template <typename sym_t>
    std::vector<aprx_occ_t<>> locate(move_r_pattern_view<sym_t> P, const search_scheme_t& scheme,
        distance_metric_t dist_metr, apm_params_t params = {}) const
        requires(supports_multiple_locate)
    {
        std::vector<aprx_occ_t<>> occ;
        locate<NO_CIGAR>(P, scheme, dist_metr, occ, params);
        return occ;
    }

    template <typename pat_t>
        requires(supports_multiple_locate
            && !move_r_is_pattern_view<std::remove_cvref_t<pat_t>>)
    std::vector<aprx_occ_t<>> locate(const pat_t& P, const search_scheme_t& scheme,
        distance_metric_t dist_metr, apm_params_t params = {}) const
    { return locate(move_r_pattern_view_t<pat_t>(P), scheme, dist_metr, params); }

    // locates a NUL-terminated pattern with at most k errors (w.r.t. dist_metr)
    std::vector<aprx_occ_t<>> locate(const char* P, const search_scheme_t& scheme,
        distance_metric_t dist_metr, apm_params_t params = {}) const
        requires(supports_multiple_locate)
    { return locate(std::string_view(P), scheme, dist_metr, params); }

    // ############################# CONSTRUCTION METHODS #############################

    /**
     * @brief builds the context tables, i.e., stores the search contexts of all patterns up to the maximum
     *        length, for which the serialized tables fit into share * base bytes (a level is stored only
     *        completely); a search can then start at a stored context
     * @param share maximum size of the serialized tables as a fraction of base (0: the tables are deleted)
     * @param base size in bytes that share refers to (0: the size of the serialized index)
     * @param single_length true <=> store only the short patterns and the patterns of the maximum length
     *        that fits into the budget, the latter with the pattern as the key (see context_tables)
     * @param p number of threads to use
     */
    void build_context_tables(double share, uint64_t base = 0, bool single_length = false, uint64_t p = 1);

    // ############################# SERIALIZATION METHODS #############################

    // writes the index to an output stream
    void serialize(std::ostream& out) const;

    // reads a serialized index from an input stream; throws std::runtime_error, if it has another profile
    void load(std::istream& in);

    // writes the index to an output stream
    std::ostream& operator>>(std::ostream& os) const
    {
        serialize(os);
        return os;
    }

    // reads a serialized index from an input stream
    std::istream& operator<<(std::istream& is)
    {
        load(is);
        return is;
    }
};

#include <move_rb/search_context.hpp>
#include <move_rb/locate_context.hpp>
#include <move_rb/extend_context.hpp>
// build and build_context_tables are compiled only in the instantiation library (see move_rb/instantiations.hpp)
#ifndef MOVE_R_DECLARATIONS_ONLY
#include <move_rb/construction/build.tpp>
#include <move_rb/construction/context_tables.tpp>
#endif
#include <move_rb/construction/construction.tpp>
#include <move_rb/queries.tpp>
#include <move_rb/load_context.tpp>
#include <move_rb/move_rb.tpp>
#include <algorithms/move_rb/apm/hamming.hpp>
#include <algorithms/move_rb/apm/edit.hpp>
// the instantiations of move_rb in the instantiation library (must be included last)
#include <move_rb/instantiations.hpp>

#endif // MOVE_R_MOVE_RB_MOVE_RB_HPP
