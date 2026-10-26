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

#ifndef MOVE_R_MOVE_RB_CONSTRUCTION_BUILD_TPP
#define MOVE_R_MOVE_RB_CONSTRUCTION_BUILD_TPP

#include <cstdlib>
#include <exception>
#include <functional>
#include <misc/utils.hpp>

#include <hash_table6.hpp>
#include <misc/fasta.hpp>
#include <misc/files.hpp>
#include <move_r/construction/text_alphabet.hpp>
#include <pfp_rlbwt/pfp_rlbwt.hpp>

#include <move_rb/move_rb.hpp>

/**
 * @brief removes an input file written by the construction, together with the caches of the RLBWT and its samples
 *        and the files of the prefix-free parsing that are named after it
 * @param input_file path of the input file
 */
inline void remove_temporary_input(const std::string& input_file)
{
    std::error_code ec;
    if (input_file.empty()) return;

    // the input file and its caches (see move_r::pfp_cache_path)
    const std::string cache = input_file + pfp_rlbwt::ext_runs;
    for (const std::string& file : { input_file, cache, cache + "-rev" }) std::filesystem::remove(file, ec);
}

/**
 * @brief the state of the construction of a move_rb that its steps share. build calls the steps in this order:
 *        1. prepare_input: reads the input file into memory, or writes the input to a file for the pfp
 *        2. build_byte_alphabet or build_integer_alphabet: maps T to the effective alphabet
 *        3. build_sub_indexes: builds idx_fwd over T and idx_bwd over T^R
 *        4. build_samplings<RIGHT>: the data structures over T^R
 *        5. load_forward_index (modes _pfp and _suffix_array_space)
 *        6. build_samplings<LEFT>, build_SA_pos: the data structures over T
 *        7. restore_input
 *        8. build_context_tables_phase
 */
template <move_r_support support>
struct move_rb<support>::build_state_t {
    move_rb& idx; // the index that is built
    move_r_input_view& input; // the input
    move_r_params params; // construction parameters

    bool byte_input = false; // true <=> the symbols of the input are bytes
    bool log = false; // true <=> print log messages
    bool pfp = false; // true <=> mode _pfp
    bool sa_space = false; // true <=> mode _suffix_array_space
    // true <=> the sub-indexes read the input file through the alphabet map; an integer input is mapped in memory
    // and read from there (see move_r::construction::build_rlbwt_pfp)
    bool reads_input_file = false;
    bool file_input = false; // true <=> the input is given as a file (params.input_file)
    uint64_t p = 1; // number of threads
    std::chrono::steady_clock::time_point time; // start of the current phase
    std::string prefix_tmp_files; // prefix of the temporary files
    std::string input_file_name; // the input file, or the file the input is written to for the pfp

    std::string idx_fwd_file_name; // file idx_fwd is stored in while idx_bwd is built

    std::vector<uint8_t> map_int; // maps the bytes of the input to their effective symbols
    std::vector<uint8_t> map_ext; // maps the effective symbols to the bytes of the input
    std::vector<uint64_t> symbols; // the distinct symbols of an integer alphabet, in increasing order
    // [0..sigma-2] order[i] = effective symbol of symbols[i]; empty <=> it is i + 1 (see _alphabet_order in move_r)
    std::vector<uint64_t> order;
    std::string map_file_name; // the file the alphabet map is written to for the pfp

    // true <=> the pfp cache of T (T^R) is rebuilt (see move_r_params::ignore_pfp_cache)
    bool ignore_cache_fwd = false;
    bool ignore_cache_bwd = false;

    uint64_t* caller_peak = nullptr; // the peak memory usage variable of the caller (see move_r_params)
    uint64_t peak_memory_usage = 0; // peak memory usage of the sub-indexes and of build in bytes

    // ends the current phase (see log_phase_end) and writes live_<phase> (see measurement_phase), the size of the
    // index in bytes at the end of the phase, to the measurement file
    void finish_phase(const std::string& key)
    {
        ::log_phase_end(log, time, params.measurement_file, key);
        if (log && params.measurement_file != nullptr) *params.measurement_file << " live_" << measurement_phase(key) << "=" << idx.size_in_bytes();
    }

    /**
     * @brief ends a phase that stores a data structure to disk, loads it or deletes it, started with
     *        log_phase_start(log, time, "") (see ::log_spill_phase_end), and writes the size of the index under the
     *        key live_<phase> (see measurement_phase)
     * @param message the message of the phase, e.g. "storing M_Phi"
     * @param bytes the number of bytes stored, loaded or deleted
     * @param stored whether the phase stores the data structure (only then its size is logged)
     * @param key measurement-file key of the phase
     */
    void finish_spill_phase(const std::string& message, uint64_t bytes, bool stored, const std::string& key)
    {
        log_spill_phase_end(log, time, message, bytes, stored, params.measurement_file, key);
        if (log && params.measurement_file != nullptr) *params.measurement_file << " live_" << measurement_phase(key) << "=" << idx.size_in_bytes();
    }
};

template <move_r_support support>
void move_rb<support>::build(
    move_r_input_view& input, move_r_params params, bool log_summary)
{
    build_state_t st { .idx = *this, .input = input, .params = std::move(params) };
    st.byte_input = input.byte_alphabet();
    st.time = now();
    const auto time_start = st.time;
    st.log = st.params.log;
    st.pfp = move_r_builds_runs_file(st.params.mode);
    st.reads_input_file = st.pfp && st.byte_input;
    st.sa_space = !st.pfp && st.params.mode == _suffix_array_space;
    const uint64_t baseline_mem_usage = malloc_count_current();
    st.prefix_tmp_files = temporary_file_path("move-r_", tmp_files_beside(st.params));
    st.file_input = !st.params.input_file.empty();
    if (st.file_input) st.input_file_name = st.params.input_file;

    // files whose pages are removed from the page cache during the construction (see page_cache_dropper)
    std::vector<std::string> dropped_files; // the output files of the prefix-free parsing, every 100 ms
    std::vector<std::string> dropped_at_end; // the input file and the caches of the RLBWT, at the end
    std::vector<std::string> dropped_slowly; // the input file and the caches of the RLBWT, every 5 s

    // the pfp caches of T and T^R (see move_r::pfp_cache_path); empty without caches (see move_r_params::pfp_cache)
    const std::string cache_fwd = st.params.pfp_cache
        ? move_r_fwd_t::pfp_cache_path(st.input_file_name, st.params.mode) : std::string();
    const std::string cache_bwd = st.params.pfp_cache
        ? move_r_bwd_t::pfp_cache_path(st.input_file_name, st.params.mode, true) : std::string();

    for (const std::string& file : { st.input_file_name, cache_fwd, cache_bwd }) {
        if (file.empty()) continue;
        dropped_at_end.emplace_back(file);
        dropped_slowly.emplace_back(file);

        // the output file of the prefix-free parsing (see pfp_rlbwt::build)
        dropped_files.emplace_back(file + pfp_rlbwt::ext_runs);
    }

    page_cache_dropper dropper(dropped_files, dropped_at_end, std::chrono::milliseconds(100), dropped_slowly);

    // removes the input file written for the prefix-free parsing if the construction ends with an exception
    struct temporary_input_guard {
        const std::string& file;
        bool owned;
        ~temporary_input_guard() { if (owned && std::uncaught_exceptions() > 0) remove_temporary_input(file); }
    } input_guard { st.input_file_name, st.reads_input_file && !st.file_input };

    prepare_input(st);

    log_phase_start(st.log, st.time, "preprocessing T");
    st.p = st.params.num_threads;

    if (st.byte_input) build_byte_alphabet(st); else build_integer_alphabet(st);
    build_sub_indexes(st);

    // the data structures over T^R are built first (idx_fwd may be stored on disk meanwhile)
    build_samplings<RIGHT>(st);
    load_forward_index(st);

    build_samplings<LEFT>(st);
    if constexpr (has_SA_pos) build_SA_pos(st);

    if (st.reads_input_file && !st.file_input) remove_temporary_input(st.input_file_name);
    if (st.file_input) drop_from_page_cache(st.input_file_name, true);
    restore_input(st);

    // the context tables are built last, they are computed with searches in the index
    if (st.params.context_table_share > 0) build_context_tables_phase(st);

    st.peak_memory_usage = std::max(st.peak_memory_usage, overall_peak() - baseline_mem_usage);
    if (st.caller_peak != nullptr) *st.caller_peak = std::max(*st.caller_peak, st.peak_memory_usage);

    if (st.log && log_summary) {
        uint64_t time_construction = time_diff_ns(time_start, now());
        log_output_size("index", size_in_bytes(), n);
        log_data_structure_sizes();
        log_construction_summary(time_construction, st.peak_memory_usage, st.params.measurement_file);
    }
}

template <move_r_support support>
void move_rb<support>::prepare_input(build_state_t& st)
{
    move_r_input_view& input = st.input;

    if (!st.pfp && st.file_input) {
        log_phase_start(st.log, st.time, "reading T");

        // reserve space for the terminator that a sub-index appends (see move_r)
        input.reserve(n);
        input.resize(n - 1);
        phase_progress prog(st.log, n - 1);

        if (st.params.fasta.text_len != 0) {
            // read the text described by the FASTA file (see read_fasta_into)
            uint64_t num_read = read_fasta_into({ st.input_file_name },
                reinterpret_cast<uint8_t*>(input.data()), st.params.fasta.alphabet,
                char(st.params.fasta.separator), true, st.params.num_threads);

            if (num_read != n - 1)
                throw std::runtime_error("move_rb: the FASTA file changed while it was being read");
        } else if (!read_file_in_blocks(st.input_file_name, n - 1, input.data(),
                       ignore_block, 64 * 1024 * 1024, &prog)) {
            throw std::runtime_error("move_rb: cannot read the input file " + st.input_file_name);
        }

        prog.finish();
        st.finish_phase("time_read_input");
    } else if (st.reads_input_file && !st.file_input) {
        log_phase_start(st.log, st.time, "storing T");

        st.input_file_name = st.prefix_tmp_files + ".input";
        direct_ofstream input_ofile(st.input_file_name, std::ios::binary);
        phase_progress prog(st.log, n - 1);
        write_in_chunks(input_ofile, input.data(), n - 1, &prog);
        prog.finish();
        input_ofile.close();

        st.finish_phase("time_write_input");
    }
}

template <move_r_support support>
void move_rb<support>::build_byte_alphabet(build_state_t& st)
{
    // a FASTA file read once: the forward index computes the alphabet map after its parsing (see build_sub_indexes)
    if (st.params.fasta.read_once != nullptr && st.params.fasta.text_len == 0) {
        st.params.alphabet_size = 0;
        st.finish_phase("time_preprocess_t_outer");
        return;
    }

    move_r_input_view& input = st.input;
    const uint64_t p = st.p;
    // [0..p-1][0..255] the frequency of each byte in the part of the input counted by each thread
    std::vector<std::vector<uint64_t>> freq_thr(p, std::vector<uint64_t>(256, 0));
    bool alphabet_from_cache = false; // true <=> the symbol frequencies are read from the cache of T
    bool cache_bwd_read = false; // true <=> the alphabet map of the cache of T^R has been read
    std::array<uint8_t, 256> cache_map_fwd {}; // alphabet map the cache of T was built with
    std::array<uint8_t, 256> cache_map_bwd {}; // alphabet map the cache of T^R was built with

    // count the bytes of the input
    if (st.params.fasta.byte_freqs != nullptr) {
        // the symbol frequencies have been counted while scanning the FASTA file (see scan_fasta)
        for (uint16_t c = 0; c < 256; c++) freq_thr[0][c] = st.params.fasta.byte_freqs[c];
    } else if (st.pfp) {
        // read the symbol frequencies and the alphabet maps from the pfp caches, if they exist (see
        // read_alphabet_from_pfp_cache)
        const std::string cache_fwd = st.params.pfp_cache
            ? move_r_fwd_t::pfp_cache_path(st.input_file_name, st.params.mode) : std::string();
        const std::string cache_bwd = st.params.pfp_cache
            ? move_r_bwd_t::pfp_cache_path(st.input_file_name, st.params.mode, true) : std::string();
        std::vector<uint64_t> cache_freq;

        if (!cache_fwd.empty() && std::filesystem::exists(cache_fwd)) {
            alphabet_from_cache = read_alphabet_from_pfp_cache(cache_fwd, p, cache_freq, cache_map_fwd);
            st.ignore_cache_fwd = !alphabet_from_cache;
        }

        if (alphabet_from_cache && !cache_bwd.empty() && std::filesystem::exists(cache_bwd)) {
            std::vector<uint64_t> freq_bwd;
            cache_bwd_read = read_alphabet_from_pfp_cache(cache_bwd, p, freq_bwd, cache_map_bwd);
            st.ignore_cache_bwd = !cache_bwd_read;
        }

        if (alphabet_from_cache) {
            freq_thr[0] = cache_freq;
        } else if (st.params.fasta.text_len != 0) {
            throw std::invalid_argument("move_rb: the symbol frequencies of the FASTA input are missing");
        } else {
            // count the symbols of the input file in parallel (see count_bytes_of_file)
            phase_progress prog(st.log, n - 1);
            count_bytes_of_file(st.input_file_name, n - 1, p, freq_thr, prog);
            prog.finish();
        }
    } else {
        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < n - 1; i++) {
            freq_thr[omp_get_thread_num()][input.at<uint8_t>(i)]++;
        }
    }

    // smallest symbol of the effective alphabet; 0 is the terminator, 1 and 2 are reserved in the mode _pfp
    // (see min_mapped_symbol in move_r)
    const uint64_t min_symbol = move_r_builds_runs_file(st.params.mode) ? 3 : 1;
    // the alphabet size includes the terminator
    st.params.alphabet_size
        = 1 + compute_byte_alphabet_map(freq_thr, min_symbol, st.params.reorder_alphabet, st.map_int);
    sigma = st.params.alphabet_size;

    // a cache that was built with a different alphabet map is rebuilt
    if (alphabet_from_cache && !std::equal(cache_map_fwd.begin(), cache_map_fwd.end(), st.map_int.begin()))
        st.ignore_cache_fwd = true;
    if (cache_bwd_read && !std::equal(cache_map_bwd.begin(), cache_map_bwd.end(), st.map_int.begin()))
        st.ignore_cache_bwd = true;

    if (st.ignore_cache_fwd || st.ignore_cache_bwd) {
        log_note(st.log, std::string("the pfp cache of ")
            + (st.ignore_cache_fwd ? (st.ignore_cache_bwd ? "T and T^R" : "T") : "T^R")
            + " was built with a different alphabet map, rebuilding it");
    }

    st.finish_phase("time_preprocess_t_outer");

    // the pfp applies the alphabet map while reading the input file
    if (st.pfp) {
        log_phase_start(st.log, st.time, "storing the alphabet map");

        st.map_file_name = st.prefix_tmp_files + ".map";
        direct_ofstream map_ofile(st.map_file_name, std::ios::binary);
        write_in_chunks(map_ofile, (char*) st.map_int.data(), 256);
        map_ofile.close();
    } else {
        log_phase_start(st.log, st.time, "mapping T");

        #pragma omp parallel for num_threads(p)
        for (uint64_t i = 0; i < n - 1; i++) {
            input.at<uint8_t>(i) = st.map_int[input.at<uint8_t>(i)];
        }
    }

    // in the mode _pfp, the symbols of the BWT are 2 smaller than in the alphabet map (see min_mapped_symbol)
    if (move_r_builds_runs_file(st.params.mode)) {
        for (uint16_t c = 0; c < 256; c++) {
            if (st.map_int[c] != 0) st.map_int[c] -= 2;
        }
    }

    invert_byte_alphabet_map(st.map_int, st.map_ext);

    st.finish_phase(st.pfp ? "time_write_alphabet_map" : "time_map_t_outer");
}

template <move_r_support support>
void move_rb<support>::build_integer_alphabet(build_state_t& st)
{
    // the distinct symbols of the input and their frequencies; a symbol is identified by its bit pattern,
    // zero-extended to 64 bits (see symbol_bit_pattern)
    std::unique_ptr<move_r_int_alphabet> alphabet = move_r_int_alphabet::make(st.input.symbol_width());
    alphabet->collect(st.input, n);
    alphabet->sorted_symbols(st.symbols, st.p);

    if (std::binary_search(st.symbols.begin(), st.symbols.end(), TERMINATOR_SYMBOL))
        throw std::invalid_argument("move_rb: the input contains the symbol TERMINATOR_SYMBOL (bit pattern 2^63)");

    st.params.alphabet_size = st.symbols.size() + 1;
    sigma = st.params.alphabet_size;

    // compute order, which sorts the symbols by decreasing frequency (see move_r_order_alphabet)
    std::vector<uint64_t> freq; // [0..sigma-2] the frequency of symbols[i]
    alphabet->frequencies(st.symbols, freq, st.p);
    move_r_order_alphabet(freq, st.params.reorder_alphabet, st.order, st.p);
    freq.clear();
    freq.shrink_to_fit();

    // map each symbol of the input to its effective symbol
    alphabet->assign_effective_symbols(st.symbols, st.order, st.p);

    log_phase_start(st.log, st.time, "mapping T");

    alphabet->map_text(st.input, n, st.p);

    st.finish_phase("time_map_t_outer");
}

template <move_r_support support>
void move_rb<support>::take_fasta_read_once(build_state_t& st, move_r_params& sub_params)
{
    // the length of the text and the alphabet map, which the forward index has computed after its parsing
    const fasta_scan_t& scan = *sub_params.fasta.read_once;
    n = scan.text_size + 1;
    sigma = idx_fwd.alphabet_size() + 1;
    st.params.alphabet_size = sigma;
    st.map_int = idx_fwd.map_int();
    invert_byte_alphabet_map(st.map_int, st.map_ext);

    // the backward index reads T^R through the map of the bytes, whose symbols are 2 larger than in the BWT (see
    // build_byte_alphabet), with the samples of the forward parsing
    std::vector<uint8_t> map(st.map_int);
    for (uint8_t& c : map) if (c != 0) c += 2;
    st.map_file_name = st.prefix_tmp_files + ".map";
    direct_ofstream map_ofile(st.map_file_name, std::ios::binary);
    write_in_chunks(map_ofile, (char*) map.data(), map.size());
    map_ofile.close();
    sub_params.pfp_alphabet_map_file = st.map_file_name;
    sub_params.alphabet_size = sigma;
    sub_params.fasta.text_len = scan.text_size;
    sub_params.fasta.byte_freqs = scan.freq.data();
    sub_params.fasta.index_samples = scan.samples.data();
    sub_params.fasta.num_index_samples = scan.samples.size();
    sub_params.fasta.sequences = &scan.data;
}

template <move_r_support support>
void move_rb<support>::build_sub_indexes(build_state_t& st)
{
    move_r_input_view& input = st.input;

    // the sub-indexes write their peak memory usage to st.peak_memory_usage
    st.caller_peak = st.params.peak_memory_usage;
    st.params.peak_memory_usage = &st.peak_memory_usage;

    move_r_params sub_params = st.params;
    // the sub-indexes read the input file (see reads_input_file) or are built from the input in memory
    sub_params.input_file = st.reads_input_file ? st.input_file_name : "";
    sub_params.tmp_files_beside = tmp_files_beside(st.params);
    sub_params.pfp_alphabet_map_file = st.map_file_name;
    sub_params.ignore_pfp_cache = st.ignore_cache_fwd;
    sub_params.measurement_key_suffix = "_fwd";
    sub_params.sample_run_ends = supports_locate;

    try {
        log_message(st.log, "building the forward index:\n");
        const bool fasta_once = sub_params.fasta.read_once != nullptr && sub_params.fasta.text_len == 0;
        idx_fwd = move_r_fwd_t(input, sub_params);
        _pfp_cached_texts = idx_fwd.rlbwt_from_cache() ? 1 : 0;
        if (fasta_once) take_fasta_read_once(st, sub_params);
        if (st.byte_input) idx_fwd.set_alphabet_maps(st.map_int, st.map_ext); else idx_fwd.set_alphabet(st.symbols, st.order);

        // store the forward index on disk while the backward index is built
        if (st.pfp || st.sa_space) {
            log_phase_start(st.log, st.time, "storing the forward index");

            st.idx_fwd_file_name = st.prefix_tmp_files + ".move_r_fwd";
            direct_ofstream idx_fwd_ofile(st.idx_fwd_file_name);
            idx_fwd.serialize(idx_fwd_ofile);
            idx_fwd = move_r_fwd_t();
            idx_fwd_ofile.close();
            drop_from_page_cache(st.idx_fwd_file_name, true);
            log_spilled_size(st.log, st.idx_fwd_file_name);

            st.finish_phase("time_store_forward_index");
        }

        if (!st.pfp) reverse(input, st.p, st.log, st.params.measurement_file, "time_reverse_t_for_backward");

        log_message(st.log, "building the backward index:\n");
        sub_params.measurement_key_suffix = "_bwd";
        // the backward index is built for T^R (see move_r_params::reversed)
        sub_params.reversed = true;
        sub_params.ignore_pfp_cache = st.ignore_cache_bwd;
        idx_bwd = move_r_bwd_t(input, sub_params);
        if (idx_bwd.rlbwt_from_cache()) _pfp_cached_texts++;
        if (st.byte_input) idx_bwd.set_alphabet_maps(st.map_int, st.map_ext); else idx_bwd.set_alphabet(st.symbols, st.order);

        if (!st.pfp) reverse(input, st.p, st.log, st.params.measurement_file, "time_restore_t");

        if (!st.map_file_name.empty()) std::filesystem::remove(st.map_file_name);
    } catch (...) {
        std::error_code ec;
        if (st.reads_input_file && !st.file_input) remove_temporary_input(st.input_file_name);
        if (!st.idx_fwd_file_name.empty()) std::filesystem::remove(st.idx_fwd_file_name, ec);
        if (!st.map_file_name.empty()) std::filesystem::remove(st.map_file_name, ec);
        throw;
    }
}

template <move_r_support support>
template <typename mds_t>
void move_rb<support>::build_sampling_phase(build_state_t& st, const std::string& name,
    const char* key, const mds_t& M, elias_fano_bit_vector& sampling)
{
    log_phase_start(st.log, st.time, name);
    sampling = sample_input_intervals(M, st.p);
    st.finish_phase(key);
}

template <move_r_support support>
template <direction_t dir>
void move_rb<support>::build_samplings(build_state_t& st)
{
    if constexpr (dir == RIGHT) {
        build_sampling_phase(st, std::string("building S_MLF_p") + direction_suffix<RIGHT>(), "time_build_s_mlf_p_bwd",
            idx_bwd.M_LF(), _S_MLF_p_bwd);
    } else {
        if constexpr (has_mphi) {
            build_sampling_phase(st, std::string("building S_MPhi_p") + direction_suffix<LEFT>(),
                "time_build_s_mphi_p_fwd", idx_fwd.M_Phi(), _S_MPhi_p_fwd);
            build_sampling_phase(st, std::string("building S_MPhi^{-1}_p") + direction_suffix<LEFT>(),
                "time_build_s_mphi_inv_p_fwd", idx_fwd.M_Phi_inv(), _S_MPhi_inv_p_fwd);
        }

        build_sampling_phase(st, std::string("building S_MLF_p") + direction_suffix<LEFT>(), "time_build_s_mlf_p_fwd",
            idx_fwd.M_LF(), _S_MLF_p_fwd);
    }
}

template <move_r_support support>
void move_rb<support>::load_forward_index(build_state_t& st)
{
    if (st.pfp || st.sa_space) {
        log_phase_start(st.log, st.time, "reading the forward index");

        direct_ifstream idx_fwd_ifile(st.idx_fwd_file_name);
        idx_fwd.load(idx_fwd_ifile);
        idx_fwd_ifile.close();
        std::filesystem::remove(st.idx_fwd_file_name);

        st.finish_phase("time_load_forward_index");
    }
}

template <move_r_support support>
void move_rb<support>::build_SA_pos(build_state_t& st)
{
    if constexpr (has_SA_pos) {
        const uint64_t p = st.p;
        log_phase_start(st.log, st.time, "building H_SA");

        // H_SA maps the text position (in T) of each SA-sample of the backward index to the index of the
        // sample; a run of length 1 has one sample, every other run one at its start and one at its end
        // (see move_r::SA_run_pos)
        const uint64_t r_bwd = idx_bwd.M_LF().num_intervals();
        const run_sample_positions& sample_pos = idx_bwd.SA_run_pos();
        emhash6::HashMap<uint64_t, uint64_t, std::identity> H_SA;
        H_SA.reserve(sample_pos.num_samples());

        for (uint64_t i = 0; i < r_bwd; i++) {
            const bool starts = idx_bwd.run_starts_at(i);
            const bool ends = idx_bwd.run_ends_at(i);
            const bool single_pos = starts && ends && idx_bwd.M_LF().p(i + 1) - idx_bwd.M_LF().p(i) == 1;
            if (starts) H_SA.emplace((n - idx_bwd.SA_run_start(i) - 1), sample_pos.start_sample_index(i));
            if (ends && !single_pos) H_SA.emplace((n - idx_bwd.SA_run_end(i) - 1), sample_pos.end_sample_index(i));
        }

        st.finish_phase("time_build_h_sa");
        log_phase_start(st.log, st.time, "building SA_pos");

        _SA_pos = interleaved_bit_aligned_vectors({ uint8_t(std::bit_width(n)) });
        _SA_pos.resize_no_init(sample_pos.num_samples(), p);

        // each thread decodes the forward suffix array in a range of input intervals of M_LF that starts
        // at a BWT run start, and looks up each value in H_SA
        uint64_t r_fwd = idx_fwd.M_LF().num_intervals();
        phase_progress prog(st.log, n);
        // returns the first input interval of the range of thread i_p
        auto range_start = [&](uint64_t i_p) {
            return i_p == 0 ? 0 : i_p >= p ? r_fwd : idx_fwd.first_interval_of_run((i_p * r_fwd) / p);
        };

        #pragma omp parallel num_threads(p)
        {
            uint64_t i_p = omp_get_thread_num();
            uint64_t x_b = range_start(i_p);
            uint64_t x_e = range_start(i_p + 1);

            if (x_b < x_e) {
                uint64_t b = idx_fwd.M_LF().p(x_b);
                uint64_t e = idx_fwd.M_LF().p(x_e) - 1;

                auto report = [&](uint64_t i, uint64_t sa) {
                    auto it = H_SA.find(sa);
                    if (it != H_SA.end()) _SA_pos.template set_parallel<0, uint64_t>(it->second, i);
                };

                uint64_t sa = idx_fwd.SA_run_start(x_b); // = SA[M_LF.p(x_b)] = SA[b]
                report(b, sa);

                if (b < e) {
                    auto dec = idx_fwd.SA_rlz().decode_context();
                    dec.init_right(b + 1);
                    dec.set_value(sa);
                    dec.report_right(e, report);
                    prog.advance((e - b) + 1);
                }
            }
        }

        prog.finish();
        st.finish_phase("time_build_sa_sample_positions");
    }
}

template <move_r_support support>
void move_rb<support>::restore_input(build_state_t& st)
{
    move_r_input_view& input = st.input;

    // the input is restored only if it was mapped in memory
    if (st.reads_input_file || st.file_input) return;
    log_phase_start(st.log, st.time, "unmapping T");

    if (st.byte_input) {
        #pragma omp parallel for num_threads(st.p)
        for (uint64_t i = 0; i < n - 1; i++) {
            input.at<uint8_t>(i) = st.map_ext[input.at<uint8_t>(i)];
        }
    } else {
        // [0..sigma-2] by_effective[c - 1] = symbol of the input with the effective symbol c (see _alphabet_order
        // in move_r)
        std::vector<uint64_t> by_effective;

        if (!st.order.empty()) {
            resize_for_overwrite(by_effective, st.symbols.size());

            #pragma omp parallel for num_threads(st.p)
            for (uint64_t i = 0; i < st.symbols.size(); i++) by_effective[st.order[i] - 1] = st.symbols[i];
        }

        move_r_unmap_text_by_symbols(input, n, st.p, st.order.empty() ? st.symbols : by_effective);
    }

    st.finish_phase("time_unmap_t_outer");
}

template <move_r_support support>
void move_rb<support>::build_context_tables_phase(build_state_t& st)
{
    const uint64_t index_size = size_in_bytes();
    // the SA-samples and M_Phi are stored on disk while the context tables are built (in the modes _pfp and
    // _suffix_array_space)

    // a data structure stored on disk
    struct stored_t {
        std::string name; // name of the data structure
        std::string file; // file it is stored in
        std::function<void(std::istream&)> load; // loads it from the file
    };

    std::vector<stored_t> stored;

    // stores a data structure on disk; the measurement key is the same for all data structures
    auto store_to_disk = [&](std::string name, const char* suffix, auto&& store,
        std::function<void(std::istream&)> load) {
        log_phase_start(st.log, st.time, "");
        std::string file = st.prefix_tmp_files + suffix;
        direct_ofstream out(file);
        store(out);
        out.close();
        drop_from_page_cache(file, true);
        st.finish_spill_phase("storing " + name, std::filesystem::file_size(file), true, "time_store_sa_samples");
        stored.push_back({ std::move(name), std::move(file), std::move(load) });
    };

    if (st.pfp || st.sa_space) {
        if constexpr (move_r_fwd_t::builds_SA_samples) {
            store_to_disk(std::string(idx_fwd.SA_samples_name()) + direction_suffix<LEFT>(), ".sa_run_fwd",
                [&](std::ostream& out) { idx_fwd.store_SA_samples(out); },
                [&](std::istream& in) { idx_fwd.load_SA_samples(in); });
        }

        if constexpr (move_r_bwd_t::builds_SA_samples) {
            store_to_disk(std::string(idx_bwd.SA_samples_name()) + direction_suffix<RIGHT>(), ".sa_run_bwd",
                [&](std::ostream& out) { idx_bwd.store_SA_samples(out); },
                [&](std::istream& in) { idx_bwd.load_SA_samples(in); });
        }

        if constexpr (has_mphi) {
            store_to_disk(std::string("M_Phi") + direction_suffix<LEFT>(), ".mphi_fwd",
                [&](std::ostream& out) { idx_fwd.store_M_Phi(out); },
                [&](std::istream& in) { idx_fwd.load_M_Phi(in); });
        }
    }

    log_phase_start(st.log, st.time, "building the context tables");
    build_context_tables(st.params.context_table_share, index_size,
        st.params.context_tables_single_length, st.p);
    st.finish_phase("time_build_context_tables");
    if (st.log) log_built_context_tables();

    for (stored_t& entry : stored) {
        log_phase_start(st.log, st.time, "");
        direct_ifstream in(entry.file);
        entry.load(in);
        in.close();
        const uint64_t bytes = std::filesystem::file_size(entry.file);
        std::filesystem::remove(entry.file);
        st.finish_spill_phase("reading " + entry.name, bytes, false, "time_load_sa_samples");
    }
}

template <move_r_support support>
void move_rb<support>::reverse(move_r_input_view& T, uint64_t p, bool log,
    std::ostream* measurement, const std::string& measurement_key)
{
    auto time = now();
    log_phase_start(log, time, "reversing T");
    move_r_reverse_text(T, n, p);
    log_phase_end(log, time, measurement, measurement_key);
}

#endif // MOVE_R_MOVE_RB_CONSTRUCTION_BUILD_TPP
