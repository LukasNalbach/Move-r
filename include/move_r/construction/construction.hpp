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
#ifndef MOVE_R_MOVE_R_CONSTRUCTION_CONSTRUCTION_HPP
#define MOVE_R_MOVE_R_CONSTRUCTION_CONSTRUCTION_HPP

#include <optional>
#include <misc/utils.hpp>

#include <misc/sort.hpp>
#include <move_r/move_r.hpp>
#include <data_structures/interleaved_bit_aligned_vectors.hpp>
#include <move_data_structure/move_data_structure.hpp>
#include <move_r/construction/pfp_output.hpp>
#include <move_r/construction/text_alphabet.hpp>
#include <misc/files.hpp>
#include <misc/log.hpp>
#include <misc/temp_files.hpp>
#include <omp.h>

// minimum number of symbols of T$ per thread of a construction (see limit_num_threads)
inline constexpr uint64_t move_r_min_symbols_per_thread = 1000;

// storing, loading or deleting fewer bytes of a data structure is not logged (one sector of O_DIRECT)
inline constexpr uint64_t min_logged_spill_bytes = direct_io_alignment;

// returns the phase of a measurement key: the key without the prefix "time_" (see log_phase_end)
inline std::string measurement_phase(const std::string& mf_key)
{
    constexpr std::string_view time_prefix = "time_";
    return mf_key.starts_with(time_prefix) ? mf_key.substr(time_prefix.size()) : mf_key;
}

/**
 * @brief ends a phase that stores a data structure to disk, loads it from disk or deletes it, started with
 *        log_phase_start(log, time, ""): logs "<message>: <size>, in ~ <time>" (a stored data structure) or
 *        "<message>, in ~ <time>" if bytes >= min_logged_spill_bytes, else nothing; the measurement keys of
 *        log_phase_end are written in both cases
 * @param log whether logging is enabled
 * @param time the phase timer
 * @param message the message of the phase, e.g. "storing M_Phi"
 * @param bytes the number of bytes stored, loaded or deleted
 * @param stored whether the phase stores the data structure (only then its size is logged)
 * @param mf the measurement file (nullptr: none)
 * @param mf_key the measurement key of the phase (empty: none)
 */
inline void log_spill_phase_end(bool log, std::chrono::steady_clock::time_point& time, const std::string& message,
    uint64_t bytes, bool stored, std::ostream* mf, const std::string& mf_key)
{
    if (!log) return;

    if (bytes >= min_logged_spill_bytes) {
        log_out() << message;
        if (stored) log_out() << ": " << format_size(bytes);
        log_phase_end(log, time, mf, mf_key);
    } else {
        if (mf != nullptr && !mf_key.empty()) {
            const std::string phase = measurement_phase(mf_key);
            *mf << " " << mf_key << "=" << time_diff_ns(time, now()) << " peak_" << phase << "=" << phase_peak()
                << " held_" << phase << "=" << phase_allocated();
        }

        note_peak();
    }
}

/**
 * @brief construction of a move-r index; stores the input, the intermediate data structures and the
 *        index to build, and implements the construction phases
 *
 * The phases (see build_index):
 * 1. the RLBWT of T$ and the C-array, from the suffix array (modes _suffix_array and _suffix_array_space) or from
 *    the output of the prefix-free parsing (mode _pfp, see build_rlbwt_pfp);
 * 2. I_LF and M_LF;
 * 3. the suffix array samples at the run boundaries and I_Phi^{-1} (and I_Phi);
 * 4. M_Phi^{-1}, M_Phi or the rlz of the differential suffix array, depending on the support;
 * 5. the data structures over L'.
 * The sections of the threads in the phases 1-3 are ranges of consecutive BWT runs (see n_p and r_p).
 */
template <move_r_support support, move_pos_encoding_t mlf_enc>
class move_r<support, mlf_enc>::construction {
public:
    construction() = delete;
    construction(construction&&) = delete;
    construction(const construction&) = delete;
    construction& operator=(construction&&) = delete;
    construction& operator=(const construction&) = delete;

    // ############################# TYPES AND CONSTANTS #############################

    // type of the suffix array entries
    using sa_t = uint40_t;

    /* an open descriptor of the cache file; the cache is read through a path that names the opened file, as long as
       the descriptor is open (Linux: /proc/self/fd/<fd>; see platform::open_pinned and pfp_cache_file) */
    struct cache_descriptor_t {
        int fd = -1;
        std::string file; // the path the file was opened with

        cache_descriptor_t() = default;
        cache_descriptor_t(const cache_descriptor_t&) = delete;
        cache_descriptor_t& operator=(const cache_descriptor_t&) = delete;
        ~cache_descriptor_t() { if (fd >= 0) platform::close(fd); }

        // opens a file; the descriptor is kept only if the file can be reopened through it (see path)
        void open(const std::string& path_to_open)
        {
            if (fd >= 0) platform::close(fd);
            file = path_to_open;
            fd = platform::open_pinned(file);
            if (fd >= 0 && !std::filesystem::exists(path())) { platform::close(fd); fd = -1; }
        }

        // returns the path that reopens the file of the descriptor
        std::string path() const { return platform::reopen_path(fd, file); }
    };

    /**
     * @brief reader of the values of the differential suffix array SA^d (for the RLZ construction) that does not
     *        need the suffix array
     *
     * A reader starts at the SA sample at the start of a run and computes the following SA values with M_Phi^{-1}.
     * Reading consecutive values takes one Phi^{-1} step per value; other positions are reached with seek().
     * Each thread uses its own reader.
     */
    struct alignas(cache_line_bytes) SA_rlz_reader {
        construction* build = nullptr; // the construction
        // marks the starting positions of the input intervals of M_LF (and position n)
        const elias_fano_bit_vector* starts = nullptr;
        // marks the starting positions of the input intervals of M_Phi^{-1} (and position n)
        const elias_fano_bit_vector* phi_starts = nullptr;
        /* [0..(n-1) >> sampled_shift] sampled[j] = SA[j << sampled_shift]; only set for positions inside runs
           longer than 2^sampled_shift (see build_SA_rlz) */
        const uint64_t* sampled = nullptr;
        uint8_t sampled_shift = 0; // log_2 of the distance of the sampled positions

        bool has_pos = false; // true <=> pos and value are set
        bool has_prev_value = false; // true <=> prev_value = SA[pos - 1] is set
        uint64_t pos = 0; // the current position in the suffix array
        uint64_t value = 0; // SA[pos]
        uint64_t prev_value = 0; // SA[pos - 1]
        uint64_t phi_interval = 0; // the input interval of M_Phi^{-1} containing value

        SA_rlz_reader() = default;

        /**
         * @brief constructs a reader
         * @param build the construction
         * @param starts marks the starting positions of the input intervals of M_LF
         * @param phi_starts marks the starting positions of the input intervals of M_Phi^{-1}
         * @param sampled [0..(n-1) >> sampled_shift] the sampled SA values (see sampled), or null
         * @param sampled_shift log_2 of the sample distance
         */
        SA_rlz_reader(construction& build, const elias_fano_bit_vector& starts,
            const elias_fano_bit_vector& phi_starts, const uint64_t* sampled, uint8_t sampled_shift)
            : build(&build), starts(&starts), phi_starts(&phi_starts), sampled(sampled), sampled_shift(sampled_shift)
        { }

        // sets phi_interval to the input interval of M_Phi^{-1} containing value
        void locate_phi_interval()
        {
            auto& phi_inv = build->idx._M_Phi_inv;
            phi_interval = std::max<uint64_t>(1, phi_starts->rank_1(value)) - 1;
            while (value >= phi_inv.p(phi_interval + 1)) phi_interval++;
        }

        // advances to the next suffix array position (one Phi^{-1} step)
        void step_right();

        /**
         * @brief returns the first input interval of M_LF of the run containing a suffix array position
         * @param target [0..n-1] position in the suffix array
         */
        uint64_t first_interval_of_run_at(uint64_t target) const
        {
            const auto& mlf = build->idx._M_LF;
            uint64_t x = std::max<uint64_t>(1, starts->rank_1(target)) - 1;
            while (target >= mlf.p(x + 1)) x++;
            return build->idx.first_interval_of_run(x);
        }

        /**
         * @brief moves the reader to the start of a run (where its SA sample is stored)
         * @param x index of the input interval of M_LF at which the run starts
         */
        void enter_run(uint64_t x)
        {
            has_prev_value = false;
            pos = build->idx._M_LF.p(x);
            // with an rlz, the samples are always stored in _SA_run (see choose_sa_samples_layout)
            value = build->idx._SA_run.template get<0, uint64_t>(build->idx._SA_run_pos.start_sample_index(x));
        }

        /**
         * @brief moves the reader to a suffix array position; starts at the start of the run containing it or at
         *        the last sampled position before it in the run (see sampled)
         * @param target [0..n-1] position in the suffix array
         */
        void seek(uint64_t target);

        /**
         * @brief returns SA^d[i] = SA[i] - SA[i-1] + n (SA[0] for i = 0)
         * @param i [0..n-1] position in SA^d
         */
        uint64_t SAd(uint64_t i);
    };

    // ############################# DATA MEMBERS #############################

    uint8_t inp_width; // width of a symbol of T_inp in bytes
    // true <=> the input uses a byte alphabet; false for a byte input with more distinct bytes than internal symbols of
    // one byte (see preprocess_t)
    bool byte_alphabet;

    move_r_input_view T_inp; // view of the container storing T, or of T_wide
    // [0..|T|-1] T with two bytes per symbol, if its internal symbols do not fit into a byte (see widen_t)
    no_init_vector<uint16_t> T_wide;
    bool T_widened = false; // true <=> T_inp is a view of T_wide; then the input itself is not mapped
    move_r<support, mlf_enc>& idx; // the index to build
    // maps the symbols of an integer-alphabet input to the internal symbols (freed once T has been mapped)
    std::unique_ptr<move_r_int_alphabet> int_alphabet;
    uint64_t int_alphabet_bytes = 0; // size of int_alphabet in bytes

    move_r_params params; // construction parameters
    uint64_t p = 1; // the number of threads to use
    // p', the number of threads that build the RLBWT, C, I_LF and I_Phi^{-1} (1 for an integer alphabet)
    uint64_t p_ = 1;
    bool delete_T = false; // controls whether T should be deleted when not needed anymore
    std::string prefix_tmp_files = ""; // prefix of temporary files
    std::chrono::steady_clock::time_point time; // time of the start of the current construction phase
    std::chrono::steady_clock::time_point time_start; // time of the start of the construction
    uint64_t baseline_mem_usage = 0; // memory allocation at the start of the construction
    // measures the peak memory usage of this construction (only if logging is enabled, see peak_scope)
    std::optional<peak_scope> peak_measurement;
    /* with a pfp cache, the path of the cache file (see move_r::pfp_cache_path); else the prefix of the files the
       prefix-free parsing writes (prefix_tmp_files) */
    std::string pfp_prefix = "";
    // true <=> the output of the prefix-free parsing is stored in a pfp cache (see move_r::pfp_cache_path)
    bool pfp_cache = false;
    /* the table of contents of the runs file written by the prefix-free parsing: the cache file, or without a
       cache, pfp_prefix + pfp_rlbwt::ext_runs (see pfp_output_file) */
    pfp_rlbwt::run_data_file pfp_output;
    /* the number of construction steps that still read pfp_output; after the last one, the file is removed without a
       cache (see release_pfp_output) */
    uint8_t pfp_output_readers = 0;
    // true <=> the existing cache stores the SA samples (the LCP samples); a rebuilt cache stores them, too
    bool pfp_cache_sa_samples = false, pfp_cache_lcp_samples = false;
    cache_descriptor_t pfp_cache_descriptor; // the descriptor through which the cache file is read

    // true <=> the RLBWT and the SA samples are read from the output of the prefix-free parsing (mode _pfp)
    bool rlbwt_from_pfp = false;

    uint64_t n = 0; // the length of T$
    // true <=> the index is built for the reversed input (see move_r_params::reversed)
    bool reversed = false;
    uint64_t r = 0; // r, the number of runs in L
    uint64_t r_ = 0; // r', the number of input/output intervals in M_LF
    uint64_t r__ = 0; // r'', the number of input/output intervals in M_Phi^{-1}
    uint64_t r___ = 0; // r''', the number of input/output intervals in M_Phi

    // [0..n-1] SA, the suffix array of T$ (modes _suffix_array and _suffix_array_space)
    no_init_vector<sa_t> sa;
    /* the runs of L, section by section: [0..p'-1] RLBWT[i_p] stores the runs of L[n_p[i_p]..n_p[i_p+1]-1] in order;
       field 0: the symbol, field 1: the length */
    std::vector<interleaved_bit_aligned_vectors> RLBWT;
    // [0..p'] n_p[i_p] = start position of thread i_p's section of L, which is the start of a run; n_p[p'] = n
    std::vector<uint64_t> n_p;
    /* [0..p'] the runs of the sections: while the RLBWT is built (see build_rlbwt_c), r_p[i_p] = the number of runs of
       section i_p; afterwards r_p[i_p] = index of the first run of section i_p, and r_p[p'] = r */
    std::vector<uint64_t> r_p;

    /* [0..p'][0..s-1] with s = num_symbol_values(): until compute_C, C[i_p][c] = the number of occurrences of c in
       section i_p of L; afterwards C[i_p][c] = the number of occurrences of c in L[0..n_p[i_p]-1] for i_p < p' (advanced
       run by run by build_ilf), and C[p'][c] = the number of occurrences of the symbols smaller than c in L (the
       C-array) */
    std::vector<std::vector<uint64_t>> C;
    // [0..r-1] I_LF, the disjoint interval sequence for LF: pair x = (i, LF(i)), where i is the start of the x-th run
    mds_pairs I_LF;
    // [0..r-1] I_LF_sym[x] = the symbol of the x-th run (8 bits if sigma <= 256, else bit_width(sigma) bits)
    interleaved_bit_aligned_vectors I_LF_sym;
    /* [0..r-1] I_Phi, the disjoint interval sequence for Phi: pair = (SA[i], SA[i-1]), where i is the start of a run;
       while LCP_run is not empty, the second component of pair x is x (see build_iphi); sorted by the first
       component by sort_iphi */
    mds_pairs I_Phi;
    /* [0..r-1] I_Phi^{-1}, the disjoint interval sequence for Phi^{-1}: pair x = (SA[i-1], SA[i]), where i is the
       start of the x-th run (pair 0: (SA[n-1], SA[0])); sorted by the first component by sort_iphi_inv */
    mds_pairs I_Phi_inv;
    // [0..r-1] LCP_run[x] = LCP[i], where i is the start of the x-th run; bit_width(n) bits each
    interleaved_bit_aligned_vectors LCP_run;
    /* [0..r-1] the values of LCP_run, sorted by the text positions SA[i] of the run starts i (the order of I_Phi after
       sort_iphi), i.e., the PLCP values at the starts of the output intervals of I_Phi^{-1}; _plcp_row_width bits
       each */
    interleaved_bit_aligned_vectors plcp_pair;

    // layout of the suffix array samples at the run boundaries (see sa_samples_layout_t and choose_sa_samples_layout)
    sa_samples_layout_t samples_layout = _at_run_boundaries;
    /* [0..r-1] run_ranks[x] = the index of the pair of the x-th run in I_Phi^{-1} sorted by the first component;
       bytes_per_value(r) bytes each; stored in the file prefix_tmp_files + ".run_ranks" if run_ranks_on_disk */
    uint_bytes_array<1> run_ranks;
    bool run_ranks_on_disk = false; // true <=> the run ranks are stored in a temporary file (see store_run_ranks)
    /* marks the first input interval of M_Phi^{-1} of each pair of I_Phi^{-1} (sorted by the first component):
       select(i) = index of the first input interval of M_Phi^{-1} of the i-th pair */
    mds_piece_marks phi_inv_first_intervals;

    // ############################# INTERNAL METHODS #############################

    // returns whether the output of the prefix-free parsing is stored in a pfp cache
    bool uses_pfp_cache() const { return pfp_cache; }

    // returns the path through which the cache file is read (see pfp_cache_descriptor)
    std::string pfp_cache_file() const
    { return pfp_cache_descriptor.fd >= 0 ? pfp_cache_descriptor.path() : pfp_prefix; }

    // returns whether the input is read from a file
    bool file_input() const { return !params.input_file.empty(); }

    // returns whether T_inp is deleted once it is not needed anymore (see delete_T and T_widened)
    bool deletes_T() const { return delete_T || T_widened; }

    // returns whether the input is stored in memory (always, except if the prefix-free parsing reads it from the file)
    bool input_in_memory() const { return !move_r_builds_runs_file(params.mode) || !file_input(); }

    // returns whether the prefix-free parsing reads the FASTA file once and scans it (see fasta_input_t::read_once)
    bool reads_fasta_once() const { return params.fasta.read_once != nullptr && params.fasta.text_len == 0; }

    // returns whether the SA is also sampled at the ends of the runs (see move_r_params::sample_run_ends)
    bool samples_at_run_ends() const
    {
        if constexpr (!has_run_end_samples) return false;
        else return !rlbwt_from_pfp || params.sample_run_ends;
    }

    /* returns the smallest internal symbol a byte is mapped to before the BWT is built: move_r_min_mapped_symbol_pfp
       with _pfp (the symbols in the RLBWT read from the output of the prefix-free parsing are decreased to
       move_r_min_mapped_symbol and above, see build_rlbwt_c), else move_r_min_mapped_symbol */
    uint8_t min_mapped_symbol() const
    { return move_r_builds_runs_file(params.mode) ? move_r_min_mapped_symbol_pfp : move_r_min_mapped_symbol; }

    // returns the number of entries in a row of C (256 for a byte alphabet, else sigma)
    uint64_t num_symbol_values() const { return byte_alphabet ? 256 : idx.sigma; }

    /**
     * @brief returns T at index i interpreted as type
     * @tparam type the type to read the symbol at index i as; it must be the input's own symbol type
     * @param i an index into T
     * @return T at index i
     */
    template <typename type>
    inline type& T(uint64_t i)
    {
        return T_inp.template at<type>(i);
    }

    /**
     * @brief sets the length of the i-th run of thread i_p's section to len
     * @param i_p [0..p'-1] thread index
     * @param i [0..r_p[i_p]-1] index of the run in the section
     * @param len run length
     */
    inline void set_run_len(uint64_t i_p, uint64_t i, uint64_t len)
    {
        RLBWT[i_p].template set<1, uint64_t>(i, len);
    }

    /**
     * @brief returns the length of the i-th run of thread i_p's section
     * @param i_p [0..p'-1] thread index
     * @param i [0..r_p[i_p]-1] index of the run in the section
     * @return run length
     */
    inline uint64_t run_len(uint64_t i_p, uint64_t i)
    {
        return RLBWT[i_p].template get<1, uint64_t>(i);
    }

    /**
     * @brief returns the symbol of the i-th run of thread i_p's section
     * @param i_p [0..p'-1] thread index
     * @param i [0..r_p[i_p]-1] index of the run in the section
     * @return symbol
     */
    inline uint64_t run_sym(uint64_t i_p, uint64_t i)
    {
        return RLBWT[i_p].template get<0, uint64_t>(i);
    }

    /**
     * @brief appends a run of length len with the symbol sym to thread i_p's section of the RLBWT
     * @param i_p [0..p'-1] thread index
     * @param sym run symbol
     * @param len run length
     */
    inline void add_run(uint64_t i_p, uint64_t sym, uint64_t len)
    {
        std::array<uint64_t, interleaved_bit_aligned_vectors::num_vectors> values {};
        values[0] = uint64_t(sym);
        values[1] = len;
        RLBWT[i_p].emplace_back(values);
    }

    /**
     * @brief builds _SA_run_pos and allocates _SA_run (see run_sample_positions)
     * @param with_ends true <=> the SA is also sampled at the ends of the runs (see samples_at_run_ends)
     */
    void prepare_sa_run(bool with_ends)
    {
        idx._SA_run_pos = run_sample_positions(r_, [&](uint64_t x) { return idx.run_starts_at(x); },
            [&](uint64_t x) { return idx.M_LF().len(x); }, p, with_ends);
        idx._SA_run = interleaved_bit_aligned_vectors({ uint8_t(std::bit_width(n)) });
        idx._SA_run.resize_no_init(idx._SA_run_pos.num_samples(), p);
    }

    // ############################# LOGGING #############################

    /** @brief starts the construction: resets the timers and the memory baseline */
    void start_construction();

    /** @brief limits the number of threads to n / move_r_min_symbols_per_thread */
    void limit_num_threads();

    /**
     * @brief starts a construction phase: resets the phase timer and prints msg (if logging)
     * @param msg message describing the phase
     */
    void log_phase_start(const std::string& msg) { ::log_phase_start(params.log, time, msg); }

    /**
     * @brief ends a construction phase: writes the elapsed time to measurement_file under the key mf_key and the
     *        current size of the index under the key live_<mf_key> (if mf_key is non-empty), and logs the
     *        runtime (if logging)
     * @param mf_key measurement-file key for the elapsed time (empty => not written)
     */
    void log_phase_end(const std::string& mf_key = "");

    // ends a construction phase whose steps log themselves: like log_phase_end, but prints nothing
    void log_phase_end_quiet(const std::string& mf_key);

    /**
     * @brief ends a phase that stores a data structure to disk or loads it, started with log_phase_start(""): like
     *        log_phase_end, after logging "<message>: <size>" (stored) or "<message>" (loaded), if
     *        bytes >= min_logged_spill_bytes, else like log_phase_end_quiet
     * @param message the message of the phase, e.g. "storing M_Phi"
     * @param bytes the number of bytes stored or loaded
     * @param stored whether the phase stores the data structure (only then its size is logged)
     * @param mf_key measurement-file key of the phase
     */
    void end_spill_phase(const std::string& message, uint64_t bytes, bool stored, const std::string& mf_key)
    {
        if (params.log && bytes >= min_logged_spill_bytes) {
            std::cout << message;
            if (stored) std::cout << ": " << format_size(bytes);
            log_phase_end(mf_key);
        } else {
            log_phase_end_quiet(mf_key);
        }
    }

    // writes the start of the RESULT line of a move data structure construction to mf_mds (see log_mds_phase_end)
    void begin_mds_result(const std::string& algo)
    {
        if (params.log && params.mf_mds != nullptr) {
            *params.mf_mds << "RESULT"
                           << " algo=" << algo
                           << " text=" << params.text_name
                           << " num_threads=" << p
                           << " a=" << idx.a;
        }
    }

    /**
     * @brief ends a move-data-structure construction phase: ends the RESULT line in mf_mds (see begin_mds_result),
     *        writes the elapsed time, the peak memory usage and the number of intervals to measurement_file
     * @param time_key measurement-file key for the elapsed time
     * @param count_key measurement-file key for the number of intervals
     * @param count number of intervals in the built move data structure
     */
    void log_mds_phase_end(const std::string& time_key, const std::string& count_key, uint64_t count);

    /** @brief logs the peak memory usage until now, unless this index is part of a larger construction */
    void log_peak_mem_usage();

    /** @brief logs a construction summary */
    void log_finished();

    /** @brief logs statistics of T */
    void log_statistics();

    // ############################# STORING TO / LOADING FROM DISK #############################

    /**
     * @brief writes data structures to the temporary file prefix_tmp_files + ext in a phase of its own
     * @tparam fnc_t type of write
     * @param name name of the data structures (for the log message)
     * @param ext extension of the file
     * @param mf_key measurement-file key of the phase
     * @param write called as write(file) with a direct_ofstream; writes the data structures and frees them
     */
    template <typename fnc_t>
    void store_to_tmp_file(const std::string& name, const std::string& ext, const std::string& mf_key, fnc_t write)
    {
        log_phase_start("");
        const std::string file_name = prefix_tmp_files + ext;
        direct_ofstream file(file_name);
        write(file);
        file.close();
        drop_from_page_cache(file_name);
        end_spill_phase("storing " + name, std::filesystem::file_size(file_name), true, mf_key);
    }

    /**
     * @brief reads data structures from the temporary file prefix_tmp_files + ext in a phase of its own (see
     *        store_to_tmp_file)
     * @tparam fnc_t type of read
     * @param name name of the data structures (for the log message)
     * @param ext extension of the file
     * @param mf_key measurement-file key of the phase
     * @param remove_file true <=> the file is deleted afterwards
     * @param read called as read(file) with a direct_ifstream; reads the data structures
     */
    template <typename fnc_t>
    void load_from_tmp_file(const std::string& name, const std::string& ext, const std::string& mf_key,
        bool remove_file, fnc_t read)
    {
        log_phase_start("");
        const std::string file_name = prefix_tmp_files + ext;
        direct_ifstream file(file_name);
        read(file);
        file.close();
        const uint64_t bytes = std::filesystem::file_size(file_name);
        if (remove_file) std::filesystem::remove(file_name);
        end_spill_phase("reading " + name, bytes, false, mf_key);
    }

    /**
     * @brief stores a data structure to disk or loads it from disk in a separate phase, and logs the number of
     *        bytes written or read (see end_spill_phase)
     * @param verb "storing" or "loading"
     * @param name name of the data structure
     * @param stream the stream to write to or read from (see direct_ofstream)
     * @param fnc function that stores or loads the data structure
     */
    template <typename stream_t, typename fnc_t>
    void log_store_or_load(const char* verb, const char* name, stream_t& stream, fnc_t fnc);

    /** @brief writes the RLBWT to a temporary file and frees it */
    void store_rlbwt();

    /** @brief reads the RLBWT from its temporary file and deletes the file */
    void load_rlbwt();

    /** @brief writes I_Phi^{-1} to a temporary file and frees it */
    void store_iphi_inv();

    /**
     * @brief reads I_Phi^{-1} from its temporary file
     * @param remove_file true <=> the file is deleted afterwards
     */
    void load_iphi_inv(bool remove_file = true);

    /** @brief writes the run ranks to a temporary file and frees them (see run_ranks) */
    void store_run_ranks();

    /** @brief writes M_LF to a temporary file and frees it */
    void store_mlf();
    /** @brief reads M_LF from its temporary file and deletes the file */
    void load_mlf();
    /** @brief writes PSB_L', PSO_L' and L_prev_min to a temporary file and frees them */
    void store_data_structures_over_l_();
    /** @brief reads PSB_L', PSO_L' and L_prev_min from their temporary file and deletes the file */
    void load_data_structures_over_l_();

    // ############################# PREFIX-FREE PARSING #############################

    // returns the path through which the runs file written by the prefix-free parsing is read
    std::string pfp_output_file() const
    {
        return uses_pfp_cache() ? pfp_cache_file() : pfp_prefix + pfp_rlbwt::ext_runs;
    }

    /* returns the number of threads that read the SA or the LCP samples of the r runs from the runs file written by
       the prefix-free parsing: at most one per MiB (2^23 bits) of bit_width(n)-bit samples */
    uint64_t pfp_reader_threads() const { return std::clamp<uint64_t>((r * uint64_t(std::bit_width(n))) >> 23, 1, p); }

    // returns whether this construction reads the SA samples written by the prefix-free parsing
    bool reads_pfp_sa_samples() const { return chooses_revert_starts() || builds_SA_samples; }

    // returns whether this construction reads the LCP samples written by the prefix-free parsing
    bool reads_pfp_lcp_samples() const { return has_mphi; }

    /**
     * @brief called by each construction step that reads the runs file of the prefix-free parsing, when it has
     *        read it: after the last one, the file is removed (without a cache)
     */
    void release_pfp_output();

    /**
     * @brief reads the alphabet map that is applied to the bytes of the input file while the prefix-free parsing
     *        reads it (see move_r_params::pfp_alphabet_map_file)
     * @return [0..255] the symbol each byte is mapped to (all zero, if there is no map file)
     */
    std::array<uint8_t, 256> pfp_alphabet_map() const;

    /**
     * @brief computes the RLBWT, the SA samples at the run boundaries and the LCP values at the run starts with the
     *        prefix-free parsing (pfp-rlbwt), or uses the pfp cache
     */
    void build_rlbwt_pfp();

    /**
     * @brief computes the output of build_rlbwt_pfp from a FASTA file that is read once (see reads_fasta_once): the
     *        alphabet map is computed from the frequencies of the symbols once the text has been parsed (see
     *        preprocess_t), and n is set afterwards
     * @param output_prefix the prefix of the output files
     * @param prm the parameters of pfp_rlbwt::build
     * @param map (output) [0..255] the map that the parsing reads the bytes through
     */
    void build_rlbwt_pfp_fasta_once(const std::string& output_prefix, const pfp_rlbwt::params& prm,
        std::array<uint8_t, 256>& map);

    /** @brief builds I_Phi^{-1} from the SA samples written by the prefix-free parsing */
    void read_iphi_inv_pfp();

    /** @brief reads the LCP values at the run starts written by the prefix-free parsing into LCP_run */
    void read_lcp_run_pfp();

    /**
     * @brief reads the table of contents of the cache file and checks that the cache was built for the alphabet map
     *        of this construction
     * @return whether the cache can be used (it must store every section this construction reads)
     */
    bool read_pfp_cache();

    /**
     * @brief renames the runs file written by the prefix-free parsing to the cache file
     * @param output_file the runs file
     */
    void write_pfp_cache(const std::string& output_file);

    // ############################# CONSTRUCTOR #############################

    /**
     * @brief constructs a move_r index of the input
     * @param index the index to build
     * @param T view of the container storing T
     * @param delete_T controls whether T should be deleted once it is not needed anymore
     * @param params_ construction parameters
     */
    construction(move_r<support, mlf_enc>& index, move_r_input_view& T, bool delete_T,
        move_r_params params_);

    // ############################# CONSTRUCTION PHASES #############################

    /** @brief builds the index (runs the construction phases in the order of the construction mode) */
    void build_index();

    // sets samples_layout (see sa_samples_layout_t)
    void choose_sa_samples_layout();

    /**
     * @brief computes the alphabet of T and maps T to the internal symbols
     * @param T_file_name name of the file containing T (if T is not in memory)
     */
    void preprocess_t(const std::string& T_file_name);

    /**
     * @brief computes the alphabet map of the bytes from their frequencies (see preprocess_t): idx.sigma and
     *        idx._map_int
     * @param freq_thr [0..q-1][0..255] the frequencies of the bytes, counted by q threads
     */
    void map_byte_alphabet(std::vector<std::vector<uint64_t>>& freq_thr);

    // completes the alphabet maps of the bytes once T is read through idx._map_int (see preprocess_t)
    void finish_byte_alphabet();

    // replaces T_inp with a view of T_wide, which stores T with two bytes per symbol (see move_r_widen_bytes)
    void widen_t();

    /** @brief reads T from a file */
    void read_t_from_file(const std::string& T_file_name);

    /** @brief maps the internal symbols of T back to the input symbols */
    void unmap_t();

    /** @brief frees int_alphabet */
    void release_int_alphabet() { int_alphabet.reset(); }

    /** @brief builds the suffix array (uses libsais40) */
    void build_sa();

    /** @brief builds the RLBWT and C */
    void build_rlbwt_c();

    /** @brief turns the numbers of occurrences of the symbols in the sections in C into the ranks and the C-array (see C) */
    void compute_C();

    /** @brief frees the RLBWT */
    void release_rlbwt()
    {
        RLBWT.clear();
        RLBWT.shrink_to_fit();
    }

    /** @brief builds I_LF and I_LF_sym from the RLBWT */
    void build_ilf();

    // returns whether the index stores revert starts (see move_r::_revert_starts)
    bool chooses_revert_starts() const { return params.revert_start_samples != 0; }

    /**
     * @brief chooses the revert starts (see move_r::_revert_starts) from I_LF and the SA values at the run starts
     *        (from sa or from the output of the prefix-free parsing)
     */
    void choose_revert_starts();

    /** @brief builds M_LF */
    void build_mlf();

    /** @brief builds SA_run from the suffix array */
    void build_sa_run_from_sa();

    /** @brief builds I_Phi^{-1} (and SA_run) from the suffix array */
    void build_iphi_inv_from_sa();

    /** @brief builds SA_run from I_Phi^{-1} (if the RLBWT was read from a file) and frees the RLBWT */
    void build_SA_run_from_iphi_inv();

    /** @brief computes the LCP values at the run starts by comparing suffixes of T (see LCP_run) */
    void build_plcp_run();

    /**
     * @brief builds I_Phi from I_Phi^{-1}; if LCP_run is not empty, the second component of each pair is the
     *        index of its run (sort_iphi replaces it)
     */
    void build_iphi();

    /** @brief sorts I_Phi by the first components and builds plcp_pair from LCP_run (if LCP_run is not empty) */
    void sort_iphi();

    /**
     * @brief sorts I_Phi^{-1} by the first components, and computes the run ranks, if the samples are stored as
     *        input intervals of M_Phi^{-1} (see run_ranks)
     */
    void sort_iphi_inv();

    /** @brief builds M_Phi^{-1} */
    void build_mphi_inv();

    /** @brief builds M_Phi */
    void build_mphi();

    /** @brief builds SA_Phi^{-1} from the run ranks and phi_inv_first_intervals (see run_ranks) */
    void build_SA_Phi_inv();

    /**
     * @brief builds the rlz of the differential suffix array; unless SA is in memory, it builds a temporary
     *        M_Phi^{-1} first (see SA_rlz_reader)
     */
    void build_SA_rlz() requires(has_SA_rlz);

    /** @brief builds the data structures over L' (depending on the supported queries and on sigma) */
    void build_data_structures_over_l_();
};

// the definitions are omitted in translation units that only need the declarations (see instantiations.hpp)
#ifndef MOVE_R_DECLARATIONS_ONLY
#include "build_index.tpp"
#include "preprocess.tpp"
#include "pfp.tpp"
#include "rlbwt.tpp"
#include "move_lf.tpp"
#include "move_phi.tpp"
#include <rlz_ref/construction/construction.hpp>
#include "sa_rlz.tpp"
#include "load_store.tpp"
#endif

#endif // MOVE_R_MOVE_R_CONSTRUCTION_CONSTRUCTION_HPP
