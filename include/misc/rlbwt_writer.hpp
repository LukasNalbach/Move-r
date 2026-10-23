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
 * @file rlbwt_writer.hpp
 * @brief writes the BWT runs of a move_r or move_rb index and the values it stores for them (see the extract-rlbwt
 *        modes)
 */

#ifndef MOVE_R_MISC_RLBWT_WRITER_HPP
#define MOVE_R_MISC_RLBWT_WRITER_HPP

#include <algorithm>
#include <cstring>
#include <string>
#include <string_view>
#include <vector>

#include <omp.h>

#include <misc/log.hpp>
#include <misc/value_writer.hpp>
#include <move_r/move_r.hpp>

// the format of the files written by write_rlbwt_runs (printed by the help of the extract-rlbwt modes)
inline constexpr std::string_view rlbwt_file_format = "The output file starts with a header of 12 bytes: the number "
    "of bytes per value (5, one byte), a byte whose bits 0 and 1 state which of the fields sa-beg and sa-end are "
    "present, n = input size + 1 and the number r of runs (40 bits each). It is followed by one tuple "
    "(run-length, run-char, sa-beg, sa-end) of the present fields per run, 40 bits per value, little endian; for "
    "n >= 2^40, the first byte is 6, the header has 14 bytes, and every value has 48 bits. run-char is the symbol of "
    "the run (the terminator is written as 0), sa-beg and sa-end are the suffix array values at the first and the "
    "last position of the run.";

/**
 * @brief writes the BWT runs of a move_r index to a file (see rlbwt_file_format): a header, then for each run [i..j]
 *        of the BWT its length, its symbol and, as selected, SA[i] and SA[j]. The input intervals of M_LF are written in
 *        parts of consecutive runs in parallel (see move_r::split_into_run_parts); the progress over the input
 *        intervals is shown in the line of the running phase.
 * @tparam width bytes per value: 5, or 6 if the positions do not fit 40 bits (see wide_values)
 * @tparam sa_beg whether SA[i] is written for each run [i..j]
 * @tparam sa_end whether SA[j] is written for each run [i..j]
 * @param idx the move_r index
 * @param path the output file
 * @param p the number of threads
 * @param log whether the progress is shown
 */
template <uint8_t width, bool sa_beg, bool sa_end, typename index_t>
void write_rlbwt_runs(const index_t& idx, const std::string& path, uint64_t p, bool log)
{
    constexpr uint64_t fields = 2 + sa_beg + sa_end; // number of values per run
    // the header: width (1 byte), the flags (1 byte; bits 0 and 1: sa-beg, sa-end), n and r (width bytes each)
    constexpr uint64_t header_bytes = 2 + 2 * width;
    const uint64_t k = idx.M_LF().num_intervals(); // k, the number of input intervals of M_LF
    const uint64_t n = idx.input_size() + 1;

    /* the input intervals of M_LF are split into num_parts parts: at least one per thread, but each with at least one
       sector of output on average; each part starts at the first input interval of a run */
    const uint64_t num_parts = std::max<uint64_t>(1, std::min<uint64_t>(p,
        idx.num_bwt_runs() * fields * width / direct_io_alignment));
    /* the first input interval and the first run of each part: [0..num_parts] part_begin[t] = the first input interval
       of M_LF in part t (it starts a BWT run), part_begin[num_parts] = k; part_run[t] = the index of the first BWT
       run of part t, part_run[num_parts] = r */
    const auto parts = idx.split_into_run_parts(num_parts, p);
    const std::vector<uint64_t>& part_begin = parts.first;
    const std::vector<uint64_t>& part_run = parts.second;
    const uint64_t r = part_run[num_parts];

    direct_positioned_file file(path, header_bytes + r * fields * width);
    char header[header_bytes] = { char(width), char(sa_beg | (sa_end << 1)) };
    std::memcpy(header + 2, &n, width);
    std::memcpy(header + 2 + width, &r, width);
    file.write(header, header_bytes, 0);
    phase_progress progress(log, k);

    #pragma omp parallel for num_threads(std::min(p, num_parts)) schedule(dynamic, 1)
    for (uint64_t t = 0; t < num_parts; t++) {
        if (part_begin[t] == part_begin[t + 1]) continue;
        value_writer<width> out(file, header_bytes, part_run[t] * fields);

        for (uint64_t x_run = part_begin[t]; x_run < part_begin[t + 1];) {
            // the input intervals [x_run..x_end-1] form the run
            uint64_t x_end = x_run + 1;
            while (x_end < k && !idx.run_starts_at(x_end)) x_end++;
            uint64_t length = 0;

            for (uint64_t x = x_run; x < x_end; x++) {
                length += idx.M_LF().len(x);
                progress.at(x);
            }

            out.push(length);
            out.push(terminator_as_zero(idx.unmap_symbol(idx.L_(x_run))));
            if constexpr (sa_beg) out.push(idx.SA_run_start(x_run));
            if constexpr (sa_end) out.push(idx.SA_run_end(x_end - 1));
            x_run = x_end;
        }
    }

    if (!file.close_file()) throw std::runtime_error("cannot write " + path);
}

/**
 * @brief writes the BWT runs of a move_r index (or of a sub-index of a move_rb index) to a file with 40 or 48 bits per
 *        value (see write_rlbwt_runs); sa-beg is written if the index supports locate, sa-end if it stores the
 *        suffix array samples at the ends of the runs
 * @tparam index_t type of the move_r index
 * @param index the move_r index
 * @param path the output file
 * @param p the number of threads
 * @param log whether the progress is shown
 */
template <typename index_t>
void write_rlbwt(const index_t& index, const std::string& path, uint64_t p, bool log)
{
    constexpr bool sa_beg = index_t::supports_locate;
    constexpr bool sa_end = sa_beg && requires { index.SA_run_end(uint64_t { 0 }); };

    if (wide_values(index.input_size() + 1)) write_rlbwt_runs<uint48_t::width, sa_beg, sa_end>(index, path, p, log);
    else write_rlbwt_runs<uint40_t::width, sa_beg, sa_end>(index, path, p, log);
}

#endif // MOVE_R_MISC_RLBWT_WRITER_HPP
