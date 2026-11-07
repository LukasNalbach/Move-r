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
 * @file spill_log.hpp
 * @brief the log messages of the phases that store data structures to disk, load them or delete them
 */

#ifndef MOVE_R_MOVE_R_CONSTRUCTION_SPILL_LOG_HPP
#define MOVE_R_MOVE_R_CONSTRUCTION_SPILL_LOG_HPP

#include <chrono>
#include <ostream>
#include <string>
#include <misc/files.hpp>
#include <misc/log.hpp>
#include <misc/utils.hpp>

// storing, loading or deleting fewer bytes of a data structure is not logged (one sector of O_DIRECT)
inline constexpr uint64_t min_logged_spill_bytes = direct_io_alignment;

/**
 * @brief returns the log message of a phase that stores data structures to disk or reads them: "storing data
 *        structures to disk" resp. "reading data structures from disk" (other messages, e.g. "deleting M_Phi_inv",
 *        are returned unchanged)
 * @param message the action and the data structure, e.g. "storing M_Phi"
 */
inline std::string spill_message(const std::string& message)
{
    if (message.starts_with("storing ")) return "storing data structures to disk";
    if (message.starts_with("reading ")) return "reading data structures from disk";
    return message;
}

/**
 * @brief ends a phase that stores a data structure to disk, loads it from disk or deletes it, started with
 *        log_phase_start(log, time, ""): logs it with log_grouped as spill_message(message) (in one line with the
 *        phases of the same message logged directly before, e.g. "storing data structures to disk, in ~ <time>"), if
 *        bytes >= min_logged_spill_bytes; the measurement keys of log_phase_end are written in both cases
 * @param log whether logging is enabled
 * @param time the phase timer; set to now
 * @param message the action and the data structure, e.g. "storing M_Phi"
 * @param bytes the number of bytes stored, loaded or deleted
 * @param stored whether the phase stores the data structure (only then its size is logged)
 * @param mf the measurement file (nullptr: none)
 * @param mf_key the measurement key of the phase (empty: none)
 */
inline void log_spill_phase_end(bool log, std::chrono::steady_clock::time_point& time, const std::string& message,
    uint64_t bytes, bool stored, std::ostream* mf, const std::string& mf_key)
{
    log_grouped_phase_end(log, time, spill_message(message), 0, false, mf, mf_key, bytes >= min_logged_spill_bytes);
}

#endif // MOVE_R_MOVE_R_CONSTRUCTION_SPILL_LOG_HPP
