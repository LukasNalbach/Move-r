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

#pragma once

#include <cstdlib>
#include <string>

/**
 * @brief prepends the bundled Big-BWT driver's directory to PATH (once per process), so that any subsequent
 *        system("bigbwt ...") call finds the repo's own `bigbwt` -- and the sibling helpers it locates next to itself --
 *        regardless of the current working directory or where the calling executable is located.
 *
 * MOVE_R_BIGBWT_DIR is the driver's absolute directory, baked in at build time by CMake (it points into the source
 * tree, so it stays valid even if the tool executable is moved). If it is not defined (Big-BWT not bundled), this is a
 * no-op and system("bigbwt") falls back to a `bigbwt` on the existing PATH. Call this before invoking bigbwt.
 */
inline void ensure_bigbwt_on_path()
{
#ifdef MOVE_R_BIGBWT_DIR
    static const bool prepended = [] {
        const std::string dir = MOVE_R_BIGBWT_DIR;
        const char* current = std::getenv("PATH");
        const std::string path = (current != nullptr && *current != '\0') ? dir + ":" + current : dir;
        ::setenv("PATH", path.c_str(), 1);
        return true;
    }();
    (void) prepended;
#endif
}
