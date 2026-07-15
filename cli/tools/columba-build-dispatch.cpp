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

// Thin front-end for a columba index builder that automatically picks the 32- or 64-bit variant based on the input
// size: columba stores index positions in a compile-time length_t, so a 32-bit build (uint32) only supports inputs up
// to ~4 GiB but produces a smaller/faster index, while larger inputs require the 64-bit build. This wrapper sums the
// sizes of the -f/--fasta-files inputs and exec()s the matching variant (COLUMBA_BUILD_32 / COLUMBA_BUILD_64, given as
// paths relative to this executable's own directory), forwarding all arguments unchanged.

#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include <libgen.h>
#include <limits.h>
#include <sys/stat.h>
#include <unistd.h>

static uint64_t file_size(const char* path)
{
    struct stat st;
    return ::stat(path, &st) == 0 ? static_cast<uint64_t>(st.st_size) : 0;
}

int main(int argc, char** argv)
{
    // sum the sizes of the input FASTA file(s); the width choice only depends on the total text length, and the file
    // size (headers + newlines included) is an upper bound on it, so deciding on it is safe
    uint64_t total = 0;
    for (int32_t i = 1; i + 1 < argc; i++) {
        if (std::strcmp(argv[i], "-f") == 0 || std::strcmp(argv[i], "--fasta-files") == 0) {
            total += file_size(argv[i + 1]);
        }
    }

    // a uint32 length_t indexes inputs up to ~4 GiB; use the 32-bit build below that, the 64-bit build otherwise (and
    // whenever the size is unknown, so an unrecognised invocation still works for arbitrarily large inputs)
    const char* variant = (total > 0 && total < UINT_MAX) ? COLUMBA_BUILD_32 : COLUMBA_BUILD_64;

    // resolve the variant next to this wrapper (both live under build/cli/)
    char self[PATH_MAX];
    int64_t len = ::readlink("/proc/self/exe", self, sizeof(self) - 1);
    self[len > 0 ? len : 0] = '\0';
    const std::string dir = len > 0 ? std::string(::dirname(self)) : std::string(".");
    const std::string exe = dir + "/" + variant;

    std::vector<char*> args(argv, argv + argc);
    args[0] = const_cast<char*>(exe.c_str());
    args.push_back(nullptr);
    ::execv(exe.c_str(), args.data());

    std::fprintf(stderr, "columba build dispatcher: cannot exec %s: %s\n", exe.c_str(), std::strerror(errno));
    return 1;
}
