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

#ifndef MOVE_R_MOVE_R_CONSTRUCTION_PFP_TPP
#define MOVE_R_MOVE_R_CONSTRUCTION_PFP_TPP

#include <cstring>
#include <filesystem>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>

#include <pfp_rlbwt/pfp_rlbwt.hpp>

#include <misc/fasta.hpp>
#include <pfp_rlbwt/fasta.hpp>
#include <move_r/construction/pfp_integer.hpp>
#include <move_r/move_r.hpp>

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::release_pfp_output()
{
    // the file is read with O_DIRECT, so a cache is not in the page cache
    if (--pfp_output_readers == 0 && !uses_pfp_cache()) std::filesystem::remove(pfp_output_file());
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
std::array<uint8_t, 256> move_r<support, mlf_enc>::construction::pfp_alphabet_map() const
{
    std::array<uint8_t, 256> map {};

    if (!params.pfp_alphabet_map_file.empty() && std::filesystem::exists(params.pfp_alphabet_map_file)) {
        direct_ifstream map_file(params.pfp_alphabet_map_file, std::ios::binary);
        read_in_chunks(map_file, (char*) map.data(), map.size());
    }

    return map;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
bool move_r<support, mlf_enc>::construction::read_pfp_cache()
{
    if (params.ignore_pfp_cache || reads_fasta_once()) return false;

    // the file is a cache <=> it has a table of contents, the runs and the alphabet map (see open_pfp_cache); it is
    // assumed to have been built for this text
    pfp_cache_alphabet alphabet {};
    if (!open_pfp_cache(pfp_cache_file(), pfp_output, alphabet)) return false;
    pfp_cache_sa_samples = pfp_output.has_sa_samples();
    pfp_cache_lcp_samples = pfp_output.has_lcp_samples();

    // the symbols of the cache must be the ones this construction reads the text with (see reorder_alphabet)
    if (alphabet.map != pfp_alphabet_map()) return false;

    // the cache must store the sections this index reads
    if (reads_pfp_sa_samples() && !pfp_cache_sa_samples) return false;
    if (reads_pfp_lcp_samples() && !pfp_cache_lcp_samples) return false;
    // a cache of a FASTA file stores its sequences, so that it can be used without scanning the FASTA file
    if (params.fasta.text_len != 0 && !pfp_output.contains(pfp_cache_fasta_section)) return false;
    return true;
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::write_pfp_cache(const std::string& output_file)
{
    pfp_cache_descriptor.open(output_file);
    std::error_code ec;
    std::filesystem::rename(output_file, pfp_prefix, ec);
    if (ec) throw std::runtime_error("move_r: cannot move " + output_file + " to " + pfp_prefix +
        " (" + ec.message() + ")");
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::build_rlbwt_pfp()
{
    // removes the map file, if it was written by this construction (see preprocess_t)
    auto remove_map_file = [&]() {
        if (!params.pfp_alphabet_map_file.empty() && params.pfp_alphabet_map_file == prefix_tmp_files + ".map") {
            std::error_code ec;
            std::filesystem::remove(params.pfp_alphabet_map_file, ec);
        }
    };

    // a cache written by another build after this build started is read through a descriptor, too
    if (uses_pfp_cache() && pfp_cache_descriptor.fd < 0 && std::filesystem::exists(pfp_prefix))
        pfp_cache_descriptor.open(pfp_prefix);

    // the steps that read the runs file: build_rlbwt_c, choose_revert_starts, read_lcp_run_pfp, read_iphi_inv_pfp
    pfp_output_readers = 1 + chooses_revert_starts() + reads_pfp_lcp_samples() + builds_SA_samples;

    if (uses_pfp_cache() && std::filesystem::exists(pfp_cache_file()) && read_pfp_cache()) {
        if (params.log) std::cout << "reading the PFP cache" << std::endl;
        idx._rlbwt_from_cache = true;
        remove_map_file();
        release_int_alphabet();
        return;
    }

    // the steps of pfp_rlbwt log themselves
    log_phase_start("");
    /* the prefix of the file the prefix-free parsing writes; with a cache, a temporary name in the directory of the
       cache file, which is renamed to the cache file at the end (see write_pfp_cache) */
    const std::filesystem::path cache_path(pfp_prefix);
    const std::string output_prefix = uses_pfp_cache()
        ? temporary_file_path(cache_path.filename().string() + ".",
            (cache_path.has_parent_path() ? cache_path : std::filesystem::path(".") / cache_path).string())
        : pfp_prefix;

    const std::string output_file = output_prefix + pfp_rlbwt::ext_runs;

    // removes the runs file written by the prefix-free parsing
    auto remove_output = [&] {
        std::error_code ec;
        std::filesystem::remove(output_file, ec);
    };

    /* a file input is read through the alphabet map (see preprocess_t; a FASTA file read once through the map computed
       after the parsing, see build_rlbwt_pfp_fasta_once); an input in memory has been remapped */
    std::array<uint8_t, 256> map = pfp_alphabet_map();
    const std::array<uint8_t, 256>* map_ptr = params.pfp_alphabet_map_file.empty() ? nullptr : &map;
    const bool fasta_once = reads_fasta_once();

    pfp_rlbwt::params prm;
    prm.num_threads = p;
    prm.format = pfp_rlbwt::file_format::plain;
    // the samples this index reads; a rebuilt cache keeps the samples of the existing one (see read_pfp_cache)
    prm.write_sa_samples = reads_pfp_sa_samples() || pfp_cache_sa_samples;
    prm.write_lcp_samples = reads_pfp_lcp_samples() || pfp_cache_lcp_samples;
    prm.log = params.log;

    // a cache stores the alphabet map and the sequences of a FASTA file, too (see pfp_cache_alphabet_section)
    if (uses_pfp_cache()) {
        prm.extra_sections.push_back({ pfp_cache_alphabet_section, [&map](std::ostream& out) {
            out.write((const char*) map.data(), map.size());
        } });

        if (fasta_once || (params.fasta.text_len != 0 && params.fasta.sequences != nullptr)) {
            prm.extra_sections.push_back({ pfp_cache_fasta_section, [&](std::ostream& out) {
                write_fasta_to_pfp_cache(out, params.fasta.text_len, *params.fasta.sequences);
            } });
        }
    }
    std::unique_ptr<pfp_rlbwt::text_source<uint8_t>> text;

    if (!byte_alphabet) {
        // an integer input in memory is replaced by its effective symbols, a file input is read through int_alphabet
        if (!file_input() && idx.symbols_remapped) {
            int_alphabet->map_text(T_inp, n, p);
            release_int_alphabet();
        }
    } else if (!file_input()) {
        text = std::make_unique<pfp_rlbwt::memory_source<uint8_t>>((const uint8_t*) T_inp.data(), n - 1, reversed);
    } else if (params.fasta.text_len != 0) {
        text = std::make_unique<pfp_rlbwt::fasta_source<uint8_t>>(params.input_file, params.fasta.text_len,
            pfp_rlbwt::fasta_transform(params.fasta.alphabet, params.fasta.separator), params.fasta.index_samples,
            params.fasta.num_index_samples, reversed, map_ptr, p);
    } else if (!fasta_once) {
        text = std::make_unique<pfp_rlbwt::file_source<uint8_t>>(params.input_file, reversed, map_ptr);
    }

    try {
        if (fasta_once) {
            build_rlbwt_pfp_fasta_once(output_prefix, prm, map);
        } else if (!byte_alphabet) {
            // the effective symbols are increased by min_mapped_symbol() - 1 (see build_rlbwt_c)
            build_rlbwt_pfp_integer({
                .text = file_input() ? nullptr : T_inp.data(),
                .file = params.input_file,
                .width = inp_width,
                .length = n - 1,
                .alphabet = file_input() ? int_alphabet.get() : nullptr,
                .max_symbol = idx.sigma - 1 }, min_mapped_symbol() - 1, reversed, output_prefix, prm);
            release_int_alphabet();
        } else {
            pfp_rlbwt::build(*text, output_prefix, prm);
        }
    } catch (...) {
        remove_output();
        remove_map_file();
        throw;
    }

    try {
        if (uses_pfp_cache()) write_pfp_cache(output_file);
        if (!pfp_output.open(pfp_output_file()))
            throw std::runtime_error("move_r: cannot read the table of contents of " + output_file);
    } catch (...) {
        remove_output();
        remove_map_file();
        throw;
    }

    remove_map_file();
    log_phase_end_quiet("time_pfp");
    log_peak_mem_usage();
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::build_rlbwt_pfp_fasta_once(const std::string& output_prefix,
    const pfp_rlbwt::params& prm, std::array<uint8_t, 256>& map)
{
    if (params.mode != _pfp)
        throw std::invalid_argument("move_r: only the construction mode pfp reads a FASTA file once");
    fasta_scan_t& scan = *params.fasta.read_once;
    std::vector<uint64_t> seq_starts;
    std::vector<std::string> seq_names;
    fasta_sequence_collector collect { 0, seq_starts, seq_names };

    pfp_rlbwt::fasta_ranges<fasta_chunk_t> text(params.input_file,
        pfp_rlbwt::fasta_transform(params.fasta.alphabet, params.fasta.separator),
        [&](const pfp_rlbwt::fasta_scan& text_scan) {
            // the text has been parsed: its length, frequencies, samples and sequences, then the alphabet map
            static_cast<pfp_rlbwt::fasta_scan&>(scan) = text_scan;
            // the index appends one sentinel to the text, so its length is text_size + 1 (see scan_fasta)
            scan.data.set_sequences(scan.text_size + 1, seq_starts, std::move(seq_names), char(params.fasta.separator));
            params.fasta.text_len = scan.text_size;
            params.fasta.byte_freqs = scan.freq.data();
            params.fasta.index_samples = scan.samples.data();
            params.fasta.num_index_samples = scan.samples.size();
            params.fasta.sequences = &scan.data;
            std::vector<std::vector<uint64_t>> freq_thr(1, std::vector<uint64_t>(scan.freq.begin(), scan.freq.end()));
            map_byte_alphabet(freq_thr);
            std::copy(idx._map_int.begin(), idx._map_int.end(), map.begin());
            finish_byte_alphabet();
            return map;
        },
        [&](fasta_chunk_t& chunk, uint64_t shift) { collect(chunk, shift); });

    pfp_rlbwt::build(text, output_prefix, prm);
    n = params.fasta.text_len + 1;
    idx.n = n;
    limit_num_threads();
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::read_lcp_run_pfp()
{
    log_phase_start("reading LCP_run");

    LCP_run = interleaved_bit_aligned_vectors({ uint8_t(std::bit_width(n)) });
    LCP_run.resize_no_init(r, p);
    // the maximum value in LCP_run is the maximum PLCP value, which determines the width of the PLCP rows
    uint64_t max_val = 0;
    const uint64_t threads = pfp_reader_threads();

    // each thread reads the LCP values of the runs [t * r / threads, (t + 1) * r / threads)
    #pragma omp parallel for num_threads(threads) schedule(static, 1) reduction(max : max_val)
    for (uint64_t t = 0; t < threads; t++) {
        const uint64_t beg = (t * r) / threads;
        const uint64_t end = ((t + 1) * r) / threads;
        pfp_rlbwt::run_stream runs(pfp_output, beg, end, false, true);
        pfp_rlbwt::run_record rec;

        for (uint64_t i = beg; i < end; i++) {
            runs.next(rec);
            LCP_run.template set_parallel<0, uint64_t>(i, rec.lcp);
            max_val = std::max(max_val, rec.lcp);
        }
    }

    idx._plcp_row_width = std::max<uint8_t>(1, std::bit_width(max_val));

    release_pfp_output();
    log_phase_end("time_read_run_lcp");
}

template <move_r_support support, move_pos_encoding_t mlf_enc>
void move_r<support, mlf_enc>::construction::read_iphi_inv_pfp()
{
    log_phase_start("reading I_Phi^{-1}");

    I_Phi_inv.reset(r, n);

    /* the prefix-free parsing has written the SA samples at the starts and at the ends of the runs (see
       pfp_rlbwt::run_record). The pair of run i is (the sample at the end of run i-1, the sample at the start of
       run i); the pair of run 0 contains the sample at the end of run r-1. Each thread reads the runs
       [t * r / threads, (t + 1) * r / threads). */
    const uint64_t threads = pfp_reader_threads();

    #pragma omp parallel for num_threads(threads) schedule(static, 1)
    for (uint64_t t = 0; t < threads; t++) {
        const uint64_t beg = (t * r) / threads;
        const uint64_t end = ((t + 1) * r) / threads;
        pfp_rlbwt::run_stream runs(pfp_output, beg, end, true, false);
        pfp_rlbwt::run_record rec;

        for (uint64_t i = beg; i < end; i++) {
            runs.next(rec);
            I_Phi_inv.set<1>(i, rec.sa_start);
            I_Phi_inv.set<0>(i + 1 == r ? 0 : i + 1, rec.sa_end);
        }
    }

    release_pfp_output();
    log_phase_end("time_read_iphi");
}

#endif // MOVE_R_MOVE_R_CONSTRUCTION_PFP_TPP
