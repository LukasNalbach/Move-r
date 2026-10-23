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

// runs the columba index builder with 32-bit positions (COLUMBA_BUILD_32) if the input files (-f, --fasta-files) have
// less than 2^32 - 1 bytes in total, else the builder with 64-bit positions (COLUMBA_BUILD_64), with the same arguments;
// both builders are in the directory of this executable

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

// returns the size of a file in bytes (0 if it cannot be read)
static uint64_t file_size(const char* path)
{
    struct stat st;
    return ::stat(path, &st) == 0 ? static_cast<uint64_t>(st.st_size) : 0;
}

int main(int argc, char** argv)
{
    // total size of the input files in bytes (an upper bound for the length of the text)
    uint64_t total = 0;
    for (int32_t i = 1; i + 1 < argc; i++) {
        if (std::strcmp(argv[i], "-f") == 0 || std::strcmp(argv[i], "--fasta-files") == 0) {
            total += file_size(argv[i + 1]);
        }
    }

    // the 64-bit builder is also used if no input file is found
    const char* variant = (total > 0 && total < UINT_MAX) ? COLUMBA_BUILD_32 : COLUMBA_BUILD_64;

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
