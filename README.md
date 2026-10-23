# Move-r

This repository contains Move-r, an optimized and parallelized r-index [1, 2], and its bi-directional variant Move-rb
[7] with approximate pattern matching, together with command-line tools, tests and benchmarks.

- **Move-r** (`move_r`) counts and locates the occurrences of a pattern, gives access to the BWT and the suffix
  array and reconstructs the text. It locates with a move data structure for Phi^{-1} or with an RLZ of the
  differential suffix array (RLZSA [3, 4]). It works over byte alphabets and over integer alphabets.
- **Move-rb** (`move_rb`) extends a pattern on both sides. It counts and locates the occurrences of a pattern and finds
  the occurrences of a pattern with at most k errors (Hamming or edit distance, with search schemes [5, 6]), locating
  with a move data structure for Phi^{-1} or with the RLZSA.

The standalone RLZSA and LZEndSA indexes of [3, 4] are part of earlier releases of Move-r, e.g.
[ALENEX-2027](https://github.com/LukasNalbach/Move-r/tree/ALENEX-2027).

## Contents

- [Build](#build)
- [Command-line tools](#command-line-tools)
- [C++ usage](#c-usage)
- [Tests](#tests)
- [Benchmarks](#benchmarks)
- [References](#references)
- [License](#license)

## Build

Requirements: a C++20 compiler (GCC 13 or later, Clang 18 or later) on a 64-bit x86 CPU, CMake 3.21 or later,
OpenMP, Intel TBB (`libtbb-dev`) and zlib (`libz-dev`). Tested on Ubuntu 24.04 with GCC 13.3.

```shell
git clone https://github.com/LukasNalbach/Move-r.git
cd Move-r
mkdir build
cd build
cmake ..
make
```

`cmake` downloads the git submodules and applies the patches in [patches/](patches/) to them. The command-line tools
are written to `build/cli/`, the tests to `build/tests/` and the benchmarks to `build/bench/`.

The submodules used by the library are [string-utils](https://github.com/LukasNalbach/string-utils) (file I/O, logging,
sorting, bit vectors; it includes sdsl-lite, libsais, ips4o, ips2ra, gtl, emhash, sux, malloc_count and googletest),
[move-data-structure](https://github.com/LukasNalbach/move-data-structure) (the move data structure),
[rlz-ref](https://github.com/LukasNalbach/rlz-ref) (the RLZ of the differential suffix array),
[pfp-rlbwt](https://github.com/LukasNalbach/pfp-rlbwt) (the run-length encoded BWT with prefix-free parsing) and
[sparse-map](https://github.com/Tessil/sparse-map). The other submodules are used by the benchmarks only.

### CMake options

The `MOVE_R_BUILD_*_CLI` options and `MOVE_R_BUILD_TESTS` are `ON` if Move-r is the top-level project and `OFF` if it is
added by another project. The benchmark options are always `OFF` unless they are switched on; only then are the
submodules of the indexes they compare with downloaded and built.

| Option | Description |
| --- | --- |
| `MOVE_R_BUILD_MOVE_R_CLI`, `MOVE_R_BUILD_MOVE_RB_CLI` | build `move-r` and `move-rb` |
| `MOVE_R_BUILD_TOOLS_CLI` | build `gen-patterns`, `gen-repetitive`, `fasta-to-text`, `print-data-structures` |
| `MOVE_R_BUILD_MOVE_R_BENCH`, `MOVE_R_BUILD_MOVE_RB_BENCH` | build the benchmarks and the indexes they compare with (`OFF`) |
| `MOVE_R_BUILD_TESTS` | build the tests |
| `MOVE_R_USE_MALLOC_COUNT` | measure the peak memory usage with malloc_count |
| `MOVE_R_MARCH_NATIVE` | compile with `-march=native` (`ON`) |
| `MOVE_R_ENABLE_LTO` | link-time optimization (`OFF`) |
| `MOVE_R_USE_CCACHE` | use ccache, if it is installed |
| `MOVE_R_FAST_TESTS` | compile the tests with `-O0` (`OFF`) |
| `MOVE_R_EXTERN_TEMPLATES` | compile the index classes once in a library instead of in every translation unit (`ON`) |
| `MOVE_R_COLUMBA_BITS` | `both` (default), `32` or `64`: the index widths of columba built for the move-rb benchmarks |

### Windows

Native MSVC (`cl.exe`) is not supported. With LLVM/Clang (`winget install LLVM.LLVM`), Visual Studio Build Tools and
Git for Windows, from an *x64 Native Tools Command Prompt*:

```shell
cmake -S . -B build -G Ninja -DCMAKE_C_COMPILER=clang ^
    -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=Release
ninja -C build
```

With MinGW-w64 (POSIX threads, e.g. WinLibs), CMake and Ninja:

```shell
cmake -S . -B build -G Ninja -DCMAKE_C_COMPILER=gcc ^
    -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release
ninja -C build
```

The benchmarks against other indexes are not built on Windows. The executables of the Clang build need `libomp.dll`,
which is copied next to them. Blocks of at least 4 MiB are allocated in large pages if the user has the privilege
"Lock pages in memory" (`MOVE_R_HUGE_PAGES=0` turns this off).

## Command-line tools

`move-r` and `move-rb` have one mode per operation: `<tool>` lists the modes, `<tool> <mode> -h` prints the options of
a mode. A build mode writes the index to `<input_file>.<tool>[-<support>]` (`-o` sets another base
name); with `-f`, the input is read as FASTA (one or more files). Patterns are read in the Pizza&Chili format, as FASTA,
FASTQ or one per line, or from stdin. The query modes run on all available threads (`-p` sets their number).
Temporary files are written to the directory `MOVE_R_TMPDIR` (default: the temporary directory of the system).

### move-r

- **build**: builds a move-r index; `-s` selects the support (`count`, `locate-one`, `locate-move` or `locate-rlzsa`)
  and `-c` the construction mode (`sa`, `sa-space` or `pfp`). `move-r build -s locate-rlzsa text.txt`
- **count**, **locate**: count or locate the occurrences of patterns; for a FASTA index, `locate -sam` and `-bed`
  write SAM and BED. `move-r locate -o occ.txt text.txt.move-r-rlzsa patterns.txt`
- **revert**: writes the text, with `-p` threads. `move-r revert -p 16 text.txt.move-r-locate text.out`
- **extract-sa**, **extract-bwt**, **extract-rlbwt**: write the suffix array, the BWT or the BWT runs with the values
  the index stores for them. `move-r extract-sa text.txt.move-r-locate text.sa`
- **info**: prints the support, the parameters and the sizes of the data structures of an index.
  `move-r info text.txt.move-r-locate`

### move-rb

- **build**: builds a move-rb index; `-s` selects the support (`count`, `locate-one`, `locate-move` or
  `locate-rlzsa`), `-context-table` the share of the index size used for the search contexts of short patterns.
  `move-rb build -f -s locate-move genome.fa`
- **count**, **locate**: count or locate the occurrences of patterns with at most `-k` errors (`-d hamming` or
  `-d edit`, with a search scheme `-scheme`); `locate -sam` and `-paf` align reads against a FASTA index. Without a
  patterns file, the patterns are read interactively.
  `move-rb locate -k 2 -d edit -sam reads.sam genome.fa.move-rb-locate reads.fq`
- **revert**, **extract-sa**, **extract-bwt**, **extract-rlbwt**: as for move-r; `-rev` for the reversed text.
  `move-rb extract-bwt -rev genome.fa.move-rb-locate genome.bwt-rev`
- **info**: prints the support, the parameters and the sizes of the data structures of an index.

### Other tools

- **gen-patterns**: extracts random substrings of a file as a Pizza&Chili patterns file.
  `gen-patterns text.txt 16 10000 patterns.txt`
- **gen-repetitive**: generates a repetitive text. `gen-repetitive -kind versioned -n 1000000 -sigma 4 text.txt`
- **fasta-to-text**: writes the text that `build -f` indexes for a FASTA file. `fasta-to-text genome.fa genome.txt`
- **print-data-structures**: builds a data structure of a small input and prints it.
  `print-data-structures move-r -t abracadabra`

## C++ usage

Add Move-r to a CMake project and link the target `move_r`:

```cmake
add_subdirectory(Move-r/ EXCLUDE_FROM_ALL)
target_link_libraries(<target> PRIVATE move_r)
```

The use of each data structure, with complete example programs for all of its functions, is described in:

- [docs/move-r.md](docs/move-r.md): `move_r`, construction, storing, count, locate, incremental search, BWT and suffix
  array access, revert;
- [docs/move-rb.md](docs/move-rb.md): `move_rb`, search contexts (extension, locate), enumerating extensions, BWT and
  suffix array access, approximate pattern matching.

## Tests

```shell
cd build
ctest
```

The tests (`build/tests/test-*`) check random inputs: at most `MOVE_R_TEST_SECONDS` seconds (default 10) or
`MOVE_R_TEST_ITERATIONS` inputs (default 1000) per function. A failed test prints its seed; `MOVE_R_TEST_SEED=<seed>`
repeats its inputs.

## Benchmarks

All benchmark tools (in `build/bench/`) print a summary; `-m <file>` appends `RESULT key=value ...` lines to a file.

- **move-r-bench**: construction and queries of move-r and other r-indexes (run from the repository directory).
- **move-r-bench-ops**, **move-rb-bench-ops**: the operations of a move-r or move-rb index for increasing pattern
  lengths.
- **move-rb-bench-apm**, **move-rb-bench-ext**: approximate matching and bi-directional extension of move-rb compared
  with columba, b-move, the br-index and sassy; their patterns are generated by **move-rb-gen-apm-queries** and
  **move-rb-gen-ext-queries**. **bri-build**, **columba-build** and **columba-rlc-build** (in `build/cli/`) build the
  compared indexes.

The comparison of Move-r with other r-indexes is in [benchmarks/move-r.md](benchmarks/move-r.md); how to repeat the
measurements of [2] and [7] is described in [measurements/move-r/replicate.md](measurements/move-r/replicate.md) and
[measurements/move-rb/replicate.md](measurements/move-rb/replicate.md).

## References

[1] Takaaki Nishimoto and Yasuo Tabei. Optimal-time queries on bwt-runs compressed indexes.
In 48th International Colloquium on Automata, Languages, and Programming (ICALP), 2021. ([paper](https://arxiv.org/abs/2006.05104))

[2] Nico Bertram, Johannes Fischer and Lukas Nalbach. Move-r: Optimizing the r-index.
In Symposium on Experimental Algorithms (SEA), 2024. ([paper](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SEA.2024.1), [short slides](slides/move-r-slides-short.pdf), [long slides](slides/move-r-slides-long.pdf))

[3] Bella Zhukova. New space-time trade-offs for pattern matching with compressed indexes.
PhD thesis, University of Helsinki, Finland, 2024. ([thesis](http://hdl.handle.net/10138/570140))

[4] Patrick Dinklage, Johannes Fischer, Lukas Nalbach and Jan Zumbrink. RLZ-r and LZ-End-r: Enhancing Move-r.
In String Processing and Information Retrieval (SPIRE), 2025. ([paper](https://arxiv.org/abs/2507.17300))

[5] Gregory Kucherov, Kirill Salikhov and Dekel Tsur. Approximate string matching using a bidirectional index.
Theoretical Computer Science, 638:145–158, 2016. ([paper](https://arxiv.org/abs/1310.1440))

[6] Luca Renders, Kathleen Marchal and Jan Fostier. Dynamic partitioning of search patterns for approximate pattern matching using search schemes (Columba).
iScience, 24(7):102687, 2021. ([paper](https://www.cell.com/iscience/fulltext/S2589-0042%2821%2900697-1))

[7] Johannes Fischer and Lukas Nalbach. Move-rb: Faster Bi-Directional r-indexes and Approximate Pattern Matching.
Accepted at the SIAM Symposium on Algorithm Engineering and Experiments (ALENEX), 2027. ([paper](https://arxiv.org/abs/2609.30089))

## License

This project is licensed under the MIT License, see [LICENSE](LICENSE). The submodules in [external/](external/) have
their own licenses.
