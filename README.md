# Move-r

This repository contains a collection of uni- and bi-directional compressed text-index implementations and the data structures they are built from. It features move data structures [1], relative Lempel-Ziv encoded (differential) suffix arrays (RLZSA) [3, 4] and Lempel-Ziv-End compressed suffix arrays (LZEndSA) [4], together with command-line tools, examples, tests and benchmarks for all of them.

## Contents

- [Move-r](#move-r)
  - [Contents](#contents)
  - [Included Indexes \& Data Structures](#included-indexes--data-structures)
  - [Dependencies](#dependencies)
    - [External (install manually)](#external-install-manually)
    - [Included (git submodules)](#included-git-submodules)
  - [CLI Build Instructions](#cli-build-instructions)
  - [Supported Compilers and Systems](#supported-compilers-and-systems)
    - [oneTBB on Windows](#onetbb-on-windows)
    - [Building on Windows (Clang)](#building-on-windows-clang)
    - [Building on Windows (MinGW-w64)](#building-on-windows-mingw-w64)
  - [CMake Build Options](#cmake-build-options)
  - [Usage in C++](#usage-in-c)
    - [Cmake](#cmake)
    - [C++](#c)
      - [Move-r](#move-r-1)
      - [Move-rb (bi-directional)](#move-rb-bi-directional)
      - [Move Data Structure](#move-data-structure)
      - [Move-r Store \& Load](#move-r-store--load)
      - [Move-r Retrieval](#move-r-retrieval)
  - [CLI Tools](#cli-tools)
    - [move-r](#move-r-2)
    - [move-rb (bi-directional)](#move-rb-bi-directional-1)
    - [Biological (FASTA / DNA) support](#biological-fasta--dna-support)
    - [rlzsa](#rlzsa)
    - [lzendsa](#lzendsa)
    - [Tools](#tools)
    - [Competitor index builders](#competitor-index-builders)
    - [Benchmark tools](#benchmark-tools)
  - [Search Schemes](#search-schemes)
  - [Examples](#examples)
  - [Tests](#tests)
  - [Benchmarks](#benchmarks)
  - [References](#references)
  - [License](#license)

## Included Indexes & Data Structures

- **Move-r** — an optimized and parallelized implementation of the modified r-index OptBWTR described in [1] ([arxiv.org](https://arxiv.org/abs/2006.05104)) (see [benchmarks](benchmarks/move-r.md)). It supports `count`, `locate` and `revert` queries as well as random access to the suffix array (SA) and the Burrows-Wheeler-Transform (BWT). Its locate support uses either a move data structure (`locate_move`) or an optimized RLZSA (`locate_rlzsa`); there are also `count`-only and `locate_one` (one occurrence per pattern) modes. Move-r works over byte alphabets as well as over arbitrary integer alphabets.
- **Move-rb and Move-rb-rlzsa (bi-directional)** — the bi-directional variants of Move-r and Move-r-rlzsa. In addition to `count`, `locate` and `revert`, `move_rb` can extend the currently matched pattern to **both** the left and the right, and it supports **approximate pattern matching** under the hamming distance (`count`/`locate`) and the edit distance (`locate`) via configurable [search schemes](#search-schemes).
- **Move data structure** — an implementation of the move data structure [1] with the balancing algorithm described in [2]. A separate variant (`move_data_structure_l_`) additionally stores a string interleaved with the arrays needed for move queries (intended for storing the characters of the BWT (sub-)runs).
- **Optimized RLZSA** — an optimized RLZSA implementation [4] that can be constructed in O(r) additional space from an r-index and is smaller and faster than the original implementation from [3]. **This is the RLZSA used by the `locate_rlzsa` mode of Move-r and Move-rb** (i.e. the `.move-rb-rlzsa` indexes and the `move_rb_rlzsa` results in the benchmarks use this optimized RLZSA).
- **RLZSA and LZEndSA indexes** — an r-index combined with an LZEndSA [4] (`r-index-lzendsa`) and a faithful reimplementation of the RLZSA-based index described in [3] (`r-index-rlzsa`). Additionally, a plain RLZSA index and a plain LZEndSA index are included, in which pattern search is implemented using binary search over the compressed SA. **These standalone `rlzsa` / `r-index-rlzsa` indexes are the original, *unoptimized* RLZSA of [3]** — distinct from the optimized RLZSA above that Move-r/Move-rb's `locate_rlzsa` uses.
- **Internal data structures** — the indexes are built from a set of reusable, header-only data structures in [include/data_structures/](include/data_structures/): plain and hybrid bit vectors, byte- and bit-aligned interleaved vectors, an sd-array and a combined rank/select support. These are tested and benchmarked separately (see [Tests](#tests) and [Benchmarks](#benchmarks)).

## Dependencies

### External (install manually)

- [OpenMP](https://www.openmp.org/)
- [Intel TBB](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb.html)
- zlib (package `libz-dev`)
- Python 3.8+ with [psutil](https://pypi.org/project/psutil/) (package `python3-psutil`) — only needed at runtime for the `bigbwt` construction mode; Big-BWT itself is bundled (see below), so it does not have to be installed

### Included (git submodules)

Core dependencies, pulled in automatically and linked into the `move_r` library:

- [sdsl-lite](https://github.com/simongog/sdsl-lite) — succinct data structures (bundled with [libdivsufsort](https://github.com/simongog/libdivsufsort))
- [libsais](https://github.com/IlyaGrebnov/libsais) — suffix-array construction
- [ips4o](https://github.com/ips4o/ips4o) — for (parallel) sorting
- [ips2ra](https://github.com/ips4o/ips2ra) — in-place radix sort, used for the integer-key occurrence sorts
- [gtl](https://github.com/greg7mdp/gtl), [sparse-map](https://github.com/Tessil/sparse-map), [emhash](https://github.com/ktprime/emhash), [unordered_dense](https://github.com/martinus/unordered_dense) — hash maps
- [sux](https://github.com/vigna/sux) — for plain bit vectors
- [ordered](https://github.com/pdinklag/ordered) — for BTrees
- [rmq](https://github.com/pdinklag/rmq) — a practical range-minimum query data structure
- [malloc_count](https://github.com/ByteHamster-etc/malloc_count) — for measuring peak memory usage
- [googletest](https://github.com/google/googletest) — for unit tests
- [Big-BWT](https://gitlab.com/manzai/Big-BWT) — prefix-free-parsing based BWT and SA construction

Used only by the benchmark tools:

- [rcomp](https://github.com/kampersanda/rcomp), [r-index](https://github.com/alshai/r-index), [r-index-f](https://github.com/drnatebrown/r-index-f), [OnlineRlbwt](https://github.com/itomomoti/OnlineRlbwt), [block_RLBWT](https://github.com/saskeli/block_RLBWT) — uni-directional competitor r-indexes
- [columba](https://github.com/biointec/columba) — bi-directional competitor index, built in two flavors: the FM-index (*columba*) and run-length-compressed (*columba-rlc*).
- [br-index](https://github.com/U-Ar/br-index) — original bi-directional r-index [10]

## CLI Build Instructions
This implementation has been tested on Ubuntu 24.04 with GCC 13.3.0, `libtbb-dev`, `libomp-dev`, `python3-psutil` and `libz-dev` installed. 
```shell
git clone --recurse-submodules https://github.com/LukasNalbach/Move-r.git
cd Move-r
mkdir build
cd build
cmake ..
cp -rf ../patched-files/* ..
make
```
The `cp -rf ../patched-files/*` step copies the files in [patched-files/](patched-files/) over the corresponding submodule sources. These patches fix compilation errors in some dependencies (e.g. compatibility fixes needed under recent compilers) and apply the code patches required by the competitor indexes.

This creates the command-line tools (`move-r-build`, `move-r-count`, `move-r-locate`, ...) in the `build/cli/` folder, the example programs in `build/examples/`, and the test and benchmark binaries in `build/tests/` and `build/bench/`. The tools are described in the [CLI Tools](#cli-tools) section below.

## Supported Compilers and Systems

Move-r is developed on Linux and its own code is portable C++20; the platform differences are almost entirely in third-party dependencies. A 64-bit x86 CPU is required (the library uses 128-bit integers and the `cmpxchg16b` instruction).

| System | Compiler | Status |
| --- | --- | --- |
| Linux (x86-64) | GCC ≥ 13 | ✅ Primary development / test target |
| Linux (x86-64) | Clang ≥ 18 | ✅ Supported |
| Windows (x64) | LLVM/Clang (`clang++`, GNU driver) | ✅ **Recommended Windows toolchain** — builds and runs |
| Windows (x64) | MinGW-w64 (g++) | ✅ Supported — needs a MinGW-built oneTBB (see below) |

Native MSVC (`cl.exe`) is not supported: it lacks `__uint128_t` and the `__atomic_*` / `__builtin_*` intrinsics that move_r and its dependencies use throughout. Use `clang++` for MSVC-ABI-compatible binaries.

Everything Windows-specific is handled automatically by the `if(WIN32)` block in `CMakeLists.txt` (the `MSVC_COMPILER` / `_REENTRANT` defines, `-mcx16`, `libatomic` for MinGW, the OpenMP/libomp wiring, a `<sys/mman.h>`/`<sys/resource.h>` shim for `sux`/`sdsl`, and copying the required runtime DLLs next to the executables).

### oneTBB on Windows

ips4o's parallel sort links Intel oneTBB. Point `-DMOVE_R_WINDOWS_TBB_ROOT=<dir>` at an extracted oneTBB:

- **Clang** (`clang++`, MSVC ABI): use a prebuilt [`oneapi-tbb-*-win`](https://github.com/uxlfoundation/oneTBB/releases) release directly.
- **MinGW** (GNU ABI): the prebuilt release is MSVC-ABI and will not link. Use a MinGW-ABI oneTBB — e.g. MSYS2's `mingw-w64-x86_64-tbb` package, or build oneTBB from source with your MinGW (`cmake -G Ninja -DTBB_TEST=OFF -DCMAKE_INSTALL_PREFIX=<dir>`, then `cmake --build . --target install`) and point `MOVE_R_WINDOWS_TBB_ROOT` at the install prefix. Note: some oneTBB versions probe the assembler with a Unix `/dev/null` in `cmake/compilers/GNU.cmake`, which fails on Windows; if configuring oneTBB errors there, replace that `/dev/null` with `NUL`.

The build picks up whichever oneTBB you point at. If you build oneTBB **static** (add `-DBUILD_SHARED_LIBS=OFF`), TBB is linked *into* the executables and no `tbb`/`tbb12` DLL is needed beside them. A prebuilt/shared oneTBB instead ships its DLL next to the exe (copied automatically).

### Building on Windows (Clang)

Requirements: **LLVM/Clang** for Windows (`winget install LLVM.LLVM` — provides `clang`, `clang++`, `libomp`), **Visual Studio** / Build Tools (MSVC runtime + Windows SDK + bundled CMake and Ninja), and a prebuilt **oneTBB** (above). From an *x64 Native Tools Command Prompt* (`vcvars64`):

```shell
cmake -S . -B build -G Ninja ^
  -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ ^
  -DCMAKE_BUILD_TYPE=Release -DMOVE_R_MARCH_NATIVE=OFF ^
  -DMOVE_R_WINDOWS_TBB_ROOT="C:/path/to/oneapi-tbb-2023.0.0"
ninja -C build
```

### Building on Windows (MinGW-w64)

Requirements: **MinGW-w64** (e.g. WinLibs, POSIX-threads build — provides `gcc`, `g++`, `libgomp`), CMake + Ninja, and a **MinGW-built oneTBB** (above). The GCC/C++/OpenMP/pthreads/atomic runtimes are all statically linked (`-static`; libgomp's `dlopen`/`dlsym` are satisfied by MinGW's static `libdl`), so no MinGW `bin` on `PATH` is required at runtime. With a **static** oneTBB (`-DBUILD_SHARED_LIBS=OFF`, recommended) the executables are fully self-contained — they import only the Windows OS DLLs (`KERNEL32` and the UCRT `api-ms-win-crt-*`) and need nothing beside them. With a shared oneTBB, only its `tbb12.dll` sits beside the exe (copied there automatically).

```shell
cmake -S . -B build -G Ninja ^
  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ ^
  -DCMAKE_BUILD_TYPE=Release -DMOVE_R_MARCH_NATIVE=OFF ^
  -DMOVE_R_WINDOWS_TBB_ROOT="C:/path/to/tbb-mingw"
ninja -C build
```

**Windows limitations:**

- The `bigbwt` construction mode is unavailable — Big-BWT is POSIX-only (POSIX semaphores, System V shared memory, file `mmap`, a Python driver). Use the default (in-memory / libsais) construction, or run Big-BWT under WSL/MSYS2; all other construction modes work.
- Peak-memory reporting is disabled — `malloc_count` relies on glibc symbol interposition, so the no-op stub is used.
- The competitor-index benchmarks (grlBWT, columba, r-index, ...) are not ported to Windows. The CLIs, examples, unit tests and the internal data-structure benchmark all build; the CLIs, tests and examples run — the sole exception is anything using the `_bigbwt` construction mode (unavailable per the point above), including `example-move-r-queries`, which demonstrates that mode.

## CMake Build Options
The build can be tailored via the following options (all `ON` by default):

| Option | Description |
| --- | --- |
| `MOVE_R_BUILD_MOVE_R_CLI` | Build the move-r command-line tools |
| `MOVE_R_BUILD_MOVE_RB_CLI` | Build the move-rb command-line tools |
| `MOVE_R_BUILD_LZENDSA_CLI` | Build the LZEndSA command-line tools |
| `MOVE_R_BUILD_RLZSA_CLI` | Build the RLZSA command-line tools |
| `MOVE_R_BUILD_MOVE_R_BENCH` | Build `move-r-bench` (pulls in the comparison indexes) |
| `MOVE_R_BUILD_MOVE_RB_BENCH` | Build `move-rb-bench-apm`/`-ext` & `move-rb-gen-apm-queries`/`-ext-queries` (pulls in the comparison indexes) |
| `MOVE_R_BUILD_EXAMPLES` | Build the example programs |
| `MOVE_R_BUILD_TESTS` | Build the unit tests |
| `MOVE_R_BUILD_INTERNAL_BENCH` | Build the internal data-structure benchmarks |
| `MOVE_R_USE_MALLOC_COUNT` | Measure peak memory usage via malloc_count (overrides malloc/free) |
| `MOVE_R_MARCH_NATIVE` | Tune the build for the host CPU (`-march=native`); turn `OFF` for portable binaries |
| `MOVE_R_ENABLE_LTO` | Enable link-time optimization (LTO/IPO) for move-r targets |

The two benchmark options (`MOVE_R_BUILD_MOVE_R_BENCH` / `MOVE_R_BUILD_MOVE_RB_BENCH`) are forced `OFF` on Windows, since the comparison indexes are not ported there.

Two additional non-boolean cache variables are available:

| Variable | Default | Description |
| --- | --- | --- |
| `MOVE_R_COLUMBA_BITS` | `both` | columba index widths to build for the `move-rb-bench-*` tools: `both`, `32` or `64` (only used when `MOVE_R_BUILD_MOVE_RB_BENCH` is `ON`) |
| `MOVE_R_WINDOWS_TBB_ROOT` | *(empty)* | Root of an extracted oneAPI TBB (Windows) release; see [oneTBB on Windows](#onetbb-on-windows) |

## Usage in C++
### Cmake
```cmake
add_subdirectory(move_r/)
option(MOVE_R_BUILD_MOVE_R_CLI OFF)
option(MOVE_R_BUILD_MOVE_RB_CLI OFF)
option(MOVE_R_BUILD_LZENDSA_CLI OFF)
option(MOVE_R_BUILD_RLZSA_CLI OFF)
option(MOVE_R_BUILD_MOVE_R_BENCH OFF)
option(MOVE_R_BUILD_MOVE_RB_BENCH OFF)
option(MOVE_R_BUILD_EXAMPLES OFF)
option(MOVE_R_BUILD_TESTS OFF)
option(MOVE_R_BUILD_INTERNAL_BENCH OFF)
option(MOVE_R_USE_MALLOC_COUNT OFF)
```

### C++

#### Move-r
```c++
#include <move_r/move_r.hpp>

int main()
{
   // build an index
   move_r<> index("This is a test string");

   // build a 64-bit index (intended for large input strings > UINT_MAX
   // bytes ~ 4GB) with only count support, use the Big-BWT
   // construction algorithm, use at most 8 threads and set the
   // balancing parameter a to 4
   move_r<_count, char, uint64_t> index_2("a large string",
      { .mode = _bigbwt, .num_threads = 8, .a = 4 });

   // print the number of bwt runs in the input string
   std::cout << index.num_bwt_runs() << std::endl;

   // print the index size
   std::cout << format_size(index.size_in_bytes()) << std::endl;

   // print the number of occurences of a pattern
   std::cout << index.count("test") << std::endl;

   // store all occurences of a pattern in a vector
   auto Occ = index.locate("is");
   for (auto o : Occ) std::cout << o << ", ";
   std::cout << std::endl;

   // build an index for an integer vector using a relative
   // lempel-ziv encoded differential suffix array (rlzsa)
   move_r<_locate_rlzsa, int32_t> index_3({ 2, -1, 5, -1, 7, 2, -1 });

   // incrementally search the pattern [2,-1] in the input vector (from
   // right to left) and print the number of occurrences after each step
   auto query = index_3.query();
   query.prepend(-1);
   std::cout << query.num_occ() << std::endl;
   query.prepend(2);
   std::cout << query.num_occ() << std::endl;

   // print the suffix-array interval [b,e] of [2,-1]
   std::cout << "b = " << query.sa_interval().first
           << ", e = " << query.sa_interval().second << std::endl;

   // incrementally locate the occurrences of [2,-1] in the input vector
   while (query.num_occ_rem() > 0) {
      std::cout << query.next_occ() << ", " << std::flush;
   }

   // compute the longest suffix of [0,7,2] that occurs in the input vector
   std::vector<int32_t> pattern = { 0, 7, 2 };
   auto query_2 = index_3.query();
   uint32_t suffix = pattern.size();
   while (suffix > 0 && query_2.prepend(pattern[suffix - 1])) suffix--;
   std::cout << std::endl << suffix << std::flush;
}
```

#### Move-rb (bi-directional)
```c++
#include <move_rb/move_rb.hpp>

int main()
{
   // build a bi-directional move-r index (move_rb); in addition to count and
   // locate, it supports extending the currently matched pattern to both the
   // left and the right, as well as approximate pattern matching
   move_rb<> index("This is a test string");

   // build a bi-directional index that uses a relative lempel-ziv encoded
   // differential suffix array (rlzsa) for the locate support
   move_rb<_locate_rlzsa> index_rlzsa("This is a test string");

   // print the index size
   std::cout << format_size(index.size_in_bytes()) << std::endl;

   // get an empty (locate-) search context for the index
   auto ctx = index.empty_context<LOCATE>();

   // search the pattern "test" by extending it in both directions: start with
   // "e", extend to the left to "te", then to the right to "tes" and "test";
   // extend(...) returns the extended context and whether the extended pattern
   // occurs in the input
   auto [c1, ok1] = ctx.extend(index, 'e', LEFT);  // "e"
   auto [c2, ok2] = c1.extend(index, 't', LEFT);   // "te"
   auto [c3, ok3] = c2.extend(index, 's', RIGHT);  // "tes"
   auto [c4, ok4] = c3.extend(index, 't', RIGHT);  // "test"
   ctx = c4;

   // print the number of occurrences of "test"
   std::cout << ctx.num_occ() << std::endl;

   // print the suffix-array interval [b,e] of "test" in the forward text
   auto [b, e] = ctx.forward_sa_interval();
   std::cout << "b = " << b << ", e = " << e << std::endl;

   // enter the locate phase and store all occurrences of "test" in a vector
   auto loc_ctx = ctx.locate_phase();
   auto Occ = loc_ctx.locate(index, ctx);
   for (auto o : Occ) std::cout << o << ", ";
   std::cout << std::endl;

   // count the occurrences of "test" with at most 1 mismatch (w.r.t. the
   // hamming distance), using a pigeon-hole search scheme
   search_scheme_t scheme = pigeon_hole_scheme(1);
   std::cout << index.count_hamming_dist("test", scheme) << std::endl;

   // locate the occurrences of "best" with at most 1 error (w.r.t. the edit
   // distance); each occurrence stores its position, length and error count
   for (auto occ : index.locate<EDIT_DISTANCE>("best", min_u_scheme(1))) {
      std::cout << "(pos=" << occ.pos << " len=" << occ.len
              << " err=" << occ.err << ") ";
   }
   std::cout << std::endl;
}
```

#### Move Data Structure
```c++
#include <move_data_structure/move_data_structure.hpp>
#include <misc/utils.hpp>
#include <misc/log.hpp>

int main()
{
   // Build a move data structure from the disjoint interval
   // sequence I = (0,1),(1,0) with n = 2
   move_data_structure<uint32_t, POS> mds({ { 0, 1 }, { 1, 0 } }, 2);

   // create a pair to perform move queries with
   std::pair<uint32_t, uint32_t> ix { 0, 0 };

   // perform some move queries
   std::cout << to_string<>(ix = mds.move(ix)) << std::endl;
   std::cout << to_string<>(ix = mds.move(ix)) << std::endl;

   // build a move data structure that additionally interleaves an extra row (here a char column, intended
   // e.g. for the characters of the bwt (sub-)runs) via the row_ts... template parameter; the extra-row
   // widths (here { 8 } bits) are passed after the construction parameters

   // use at most 4 threads and set a = 2
   move_data_structure<uint32_t, POS, char> mds_str({
      { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 }, { 4, 0 } }, 8,
      { .num_threads = 4, .a = 2 }, { 8 });

   // this disjoint interval sequence is not 2-balanced, because the output
   // interval [0,3] contains 4 >= 2a = 4 input intervals

   // print the pairs of the resulting disjoint interval sequence
   for (uint32_t i = 0; i < mds_str.num_intervals(); i++) {
      std::cout << to_string<>({ mds_str.p(i), mds_str.q(i) });
   }

   // the balancing algorithm has added the pair (6, 2) to balance the sequence
}
```

#### Move-r Store & Load
```c++
#include <move_r/move_r.hpp>

int main()
{
   // build an index
   move_r<> index("This is a test string");

   // store an index in a file
   std::ofstream index_ofile("test_idx.move-r");
   index >> index_ofile;
   index_ofile.close();

   // load the same index into another move_r-object
   std::ifstream index_ifile("test_idx.move-r");
   move_r<> reloaded_index;
   reloaded_index << index_ifile;
   index_ifile.close();
}
```

#### Move-r Retrieval
```c++
#include <move_r/move_r.hpp>

int main()
{
   // build an index
   move_r<> index("This is a test string");

   // retrieve the range [8,17] of the original text and store
   // it in a string using at most 2 threads
   std::string reverted_range = index.revert_range(
      { .l = 8, .r = 17, .num_threads = 2 });
   for (auto c : reverted_range) std::cout << c;
   std::cout << std::endl;

   // print the original text from right to left without storing it
   // using 1 thread
   index.revert_range([](auto, auto c) { std::cout << c; }, { .num_threads = 1 });
   std::cout << std::endl;

   // retrieve the suffix array values in the range [2,6] using at
   // most 4 threads and store them in a vector
   std::vector<uint32_t> SA_range = index.SA_range(
      { .l = 2, .r = 6, .num_threads = 4 });
   for (auto s : SA_range) std::cout << s << ", ";
   std::cout << std::endl;

   // print SA[1]
   std::cout << index.SA(1) << std::endl;

   // retrieve the BWT in the range [7,14] from left to right
   // using 1 thread
   index.BWT_range([](auto, auto s) { std::cout << s << ", "; },
      { .l = 7, .r = 14, .num_threads = 1 });
   std::cout << std::endl;

   // print BWT[16]
   std::cout << index.BWT(16) << std::endl;
}
```

## CLI Tools
The command-line tools are built into the `build/cli/` folder. Each tool prints a detailed description of its parameters when run without arguments. Patterns files are expected in the [Pizza&Chili](http://pizzachili.dcc.uchile.cl/) format (which can be produced with `gen-patterns`).

### move-r

- **move-r-build** — builds a move-r index for an input file. `-s` selects the locate support (`count`, `locate_move` or `locate_rlzsa`), `-c` the construction mode (`sa` or `bigbwt`), `-o` the output base name, `-p` the number of threads and `-a` the balancing parameter; `-m_idx`/`-m_mds` write measurement data. `-f` reads the input(s) as FASTA (strip headers, mask non-`ACGT` bases; `-A` sets the kept alphabet); multiple FASTA files may be given and are concatenated into one index.
- **move-r-count** — counts the occurrences of each pattern in a patterns file (`-o` writes the counts as a `pat<i>\t<count>` TSV, `-m` writes measurement data).
- **move-r-locate** — locates the occurrences of each pattern; can verify them against the input (`-c`/`-i`) and write them to a file (`-o`). For a FASTA-built index it can also report the occurrences in the biological **SAM** (`-sam`) and **BED** (`-bed`) formats with sequence-relative coordinates; by default it matches each pattern on both strands (the reverse-complement hits get strand `-` / SAM flag `0x10`; `-norc` disables this).
- **move-r-revert** — reconstructs the original file from the index (`-im` reverts in memory, `-p` sets the number of threads).
- **move-r-query** — loads an index and answers `count`/`locate` queries typed interactively on stdin.

### move-rb (bi-directional)

- **move-rb-build** — builds a bi-directional move-r index (`-s`: `count`, `locate_move` or `locate_rlzsa`); other options as `move-r-build`. In FASTA mode (`-f`, optionally over several concatenated FASTA files) it also stores per-sequence names/boundaries for SAM output.
- **move-rb-count** — counts the approximate occurrences of each pattern w.r.t. the hamming distance using a search scheme (`-s pigeon_hole|suffix_filter|min_u|01|<file>`, `-k <mismatches>`). `-o` writes the counts as a TSV; `-hist` writes a per-mismatch-count histogram (`pat<i>\t<c_0>\t…\t<c_k>`).
- **move-rb-locate** — locates the approximate occurrences of each pattern (`-d hamming|edit`, `-s <scheme>`, `-k <errors>`); can verify them against the input (`-c`/`-i`). For a FASTA-built index it aligns reads (FASTQ via `-fastq`) on both strands and writes them in **SAM** (`-sam`, with `CIGAR`, `NM`/`MD` tags, `-seq` for `SEQ`/`QUAL`, `-nocigar` to skip the CIGAR, `-best` for one best hit, `-norc` to disable reverse-complement matching) and/or **PAF** (`-paf`) format.
- **move-rb-revert** — reconstructs the original file from the index.
- **move-rb-query** — loads an index and answers exact and approximate (`hamming-count`, `hamming-locate`, `edit-locate`) queries typed interactively on stdin.

### Biological (FASTA / DNA) support

Building either index with `-f` reads the input as (multi-sequence) FASTA — headers are stripped, non-`ACGT` bases masked to `N`, sequences delimited by a reserved separator that no search crosses — and stores each sequence's name and boundary so occurrences can be reported in sequence-relative coordinates. On such an index the locate tools emit the standard biological formats:

- **move-r-locate** (exact matches): `-sam` and `-bed`, both strands (`-norc` to disable).
- **move-rb-locate** (approximate matches): reads (incl. FASTQ) aligned on both strands, written as **SAM** (`CIGAR` with `=`/`X`/`I`/`D`, `NM:i` edit distance, `MD:Z` reference-mismatch string, optional `SEQ`/`QUAL`) and/or **PAF**; `-best` keeps only one best alignment per read.
- **move-rb-count**: `-hist` reports the occurrence count per number of mismatches.

### rlzsa

- **rlzsa-build** — builds a plain RLZSA index from a text file (`-d` sampling parameter, `--bigbwt`, `--f64`).
- **rlzsa-count** / **rlzsa-locate** — count / locate the occurrences of the patterns over a plain RLZSA index (binary search over the compressed SA); `rlzsa-locate --c` verifies against the text.
- **rlzsa-random-access** — measures the random-access time of an RLZSA (`-s` seed, `-q` queries, `-l` interval length).
- **rlzsa-convert** — converts an `r-index-rlzsa` index to a plain `rlzsa` index.
- **r-index-rlzsa-build** — builds an r-index combined with an RLZSA (`--r-index-samples` to use the SA-samples in the r-index instead of the literal phrases).
- **r-index-rlzsa-count** / **r-index-rlzsa-locate** — count / locate the occurrences of the patterns; `-c <input_file>` verifies the results.

### lzendsa

- **lzendsa-build** — builds a plain LZEndSA index from a text file (`-d` sampling parameter, `-h` maximum phrase length, `--bigbwt`).
- **lzendsa-count** / **lzendsa-locate** — count / locate the occurrences of the patterns over a plain LZEndSA index (binary search over the compressed SA); `lzendsa-locate --c` verifies against the text.
- **lzendsa-random-access** — measures the random-access time of an LZEndSA (`-s` seed, `-q` queries, `-l` interval length).
- **lzendsa-convert** — converts an `r-index-lzendsa` index to a plain `lzendsa` index.
- **r-index-lzendsa-build** — builds an r-index combined with an LZEndSA (`-h` maximum phrase length, `--lzend-samples` to use SA-samples at LZ-End phrase ends, `--f64`).
- **r-index-lzendsa-count** / **r-index-lzendsa-locate** — count / locate the occurrences of the patterns; `-c <input_file>` verifies the results.

### Tools

- **gen-patterns** — randomly extracts `<number>` substrings of length `<length>` from a file, producing a patterns file in Pizza&Chili format (an optional `<forbidden>` set of characters can be excluded).
- **print-data-structures** — builds one of `rlzsa`, `rlzsa_opt`, `lzendsa`, `mds`, `move_r` or `move_rb` from a small input (`-t`/`-f`/stdin) and prints its contents; useful for inspection and debugging.

### Competitor index builders

Index construction tools for the competitor indexes measured by the `move-rb-bench-*` tools (built into `build/cli/`):

- **bri-build** — builds the [br-index](https://github.com/U-Ar/br-index) (`.bri`) for an input file.
- **columba-build** — builds the columba FM-index. It automatically selects columba's 32- or 64-bit index-position type based on the total input size (a 32-bit index supports inputs up to ~4 GiB and is smaller/faster; larger inputs use the 64-bit build). Multiple FASTA files may be given with repeated `-f`. `-t` sets the number of OpenMP threads for suffix-array construction (libsais).
- **columba-rlc-build** — builds the run-length-compressed *columba-rlc* (b-move) index, with the same automatic 32-/64-bit selection and `-t` threading. `-p` builds the index via prefix-free parsing (using the bundled Big-BWT), which is recommended for large, repetitive inputs.

### Benchmark tools

- **move-r-bench** — benchmarks construction-, revert- and query performance of move-r against block-rlbwt (`-2`/`-v`/`-r`), r-index, r-index-f, rcomp-glfig and online-rlbwt; with `-a` it instead measures count/locate performance of move-r for a range of balancing parameters. Has to be run from the base folder.
- **move-rb-bench-apm** — benchmarks approximate count/locate of move_rb (move & rlzsa) against both columba flavors (`columba_rlc`/`columba`) and br-index, each passed by its own flag and measured only if given. Options select the metric (`--metric`), operation (`--only` count/locate/both), pattern subsets (`--k`/`--m`), scheme (`-s`), replay time (`--time`) and CIGAR-producing locate (`--cigar`); pattern sets come from `move-rb-gen-apm-queries`.
- **move-rb-gen-apm-queries** — generates the auto-calibrated pattern sets used by `move-rb-bench-apm`: for every (k, m) pair it writes pizza&chili-format files for count-hamming, locate-hamming and locate-edit by sampling text substrings and injecting up to k errors, with each file's pattern count calibrated (against a move_rb index) so its operation runs for about `--time` seconds. Options set the error counts/lengths (`-k`/`-m`), operations (`--only`), timing (`--time`/`--timeout`), minimum patterns (`--min`), forbidden characters (`-x`), scheme, seed and output dir.
- **move-rb-bench-ext** — benchmarks exact bidirectional extension (count) and locate on the same indexes, each query extending a random-start pattern left/right in random order. Count and locate use separate pattern files from `move-rb-gen-ext-queries`; options select the operation (`--only` count/locate/both), pattern lengths (`--m`) and replay time (`--time`). `RESULT` lines use the `move-rb-bench-apm` schema (`dist_metr=exact`).
- **move-rb-gen-ext-queries** — generates the pattern sets for `move-rb-bench-ext`: for every length m a `count` and a `locate` file of sampled text substrings, each auto-calibrated (against a move_rb index) to run for about `--time` seconds. Options set the lengths (`-m`), operations (`--only`), timing (`--time`/`--timeout`), minimum patterns (`--min`), forbidden characters (`-x`) and output dir.
- **bench-int-rank-select** — benchmarks the internal integer rank/select data structures (built into `build/bench/`).

## Search Schemes
Approximate pattern matching in `move_rb` is driven by *search schemes*. The built-in schemes `pigeon_hole`, `suffix_filter`, `min_u` and `01` are parameterized by the maximum number of errors `k`; alternatively a scheme can be read from a file. A collection of predefined scheme files (e.g. `kianfar`, `kuch_k+1`, `kuch_k+2`, `man_best`, each for several values of `k`) is provided in [search_schemes/](search_schemes/).

## Examples
Self-contained example programs corresponding to the C++ snippets above live in [examples/](examples/) and are built into `build/examples/`:

- `move_r-queries` — count / locate / incremental query usage of `move_r`
- `move-rb-queries` — bi-directional extension and approximate matching with `move_rb`
- `move_data_structure` — building and querying a move data structure
- `move_r-retrieval` — revert, SA and BWT retrieval
- `move_r-store-load` — storing an index to and loading it from disk

## Tests
GoogleTest-based unit tests live in [tests/](tests/) and are built into `build/tests/`:

- `test-move-r` — the uni-directional move-r index over all support and value types
- `test-move-rb` — the bi-directional index, including approximate pattern matching
- `test-move-data-structure` — the move data structure and its balancing algorithm
- `test-data-structures` — the internal bit vectors, interleaved vectors and rank/select support
- `test-sa-index` — the RLZSA / LZEndSA suffix-array indexes

## Benchmarks
A comparison of Move-r with other r-indexes (tested texts, query and construction performance) can be found in [benchmarks/move-r.md](benchmarks/move-r.md). The instructions on how to replicate the measurements presented in [2] can be found [here](measurements/move-r/replicate.md).

## References
[1] Takaaki Nishimoto and Yasuo Tabei. Optimal-time queries on bwt-runs compressed indexes.
In 48th International Colloquium on Automata, Languages, and Programming (ICALP), 2021. ([paper](https://arxiv.org/abs/2006.05104))

[2] Nico Bertram, Johannes Fischer and Lukas Nalbach. Move-r: Optimizing the r-index.
In Symposium on Experimental Algorithms (SEA), 2024. ([paper](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SEA.2024.1), [short slides](slides/move-r-slides-short.pdf), [long slides](slides/move-r-slides-long.pdf))

[3] Bella Zhukova. New space-time trade-offs for pattern matching with compressed indexes.
PhD thesis, University of Helsinki, Finland, 2024. ([thesis](http://hdl.handle.net/10138/570140))

[4] Patrick Dinklage, Johannes Fischer, Lukas Nalbach and Jan Zumbrink. RLZ-r and LZ-End-r: Enhancing Move-r.
In String Processing and Information Retrieval (SPIRE), 2024. ([paper](https://arxiv.org/abs/2507.17300))

[5] Gene Myers. A fast bit-vector algorithm for approximate string matching based on dynamic programming.
Journal of the ACM, 46(3):395–415, 1999. ([paper](https://doi.org/10.1145/316542.316550))

[6] Heikki Hyyrö. Explaining and extending the bit-parallel approximate string matching algorithm of Myers.
Technical Report A-2001-10, University of Tampere, 2001. ([report](https://researchportal.tuni.fi/en/publications/explaining-and-extending-the-bit-parallel-approximate-string-matc/))

[7] Gregory Kucherov, Kirill Salikhov and Dekel Tsur. Approximate string matching using a bidirectional index.
Theoretical Computer Science, 638:145–158, 2016. ([paper](https://arxiv.org/abs/1310.1440))

[8] Luca Renders, Kathleen Marchal and Jan Fostier. Dynamic partitioning of search patterns for approximate pattern matching using search schemes (Columba).
iScience, 24(7):102687, 2021. ([paper](https://www.cell.com/iscience/fulltext/S2589-0042%2821%2900697-1))

[9] Lore Depuydt, Luca Renders, Simon Van de Vyver, Lennart Veys, Travis Gagie and Jan Fostier. b-move: faster bidirectional character extensions in a run-length compressed index.
In 24th International Workshop on Algorithms in Bioinformatics (WABI), LIPIcs vol. 312, pages 10:1–10:18, 2024. ([paper](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.WABI.2024.10))

[10] Yuma Arakawa, Gonzalo Navarro and Kunihiko Sadakane. Bi-directional r-indexes.
In 33rd Annual Symposium on Combinatorial Pattern Matching (CPM), LIPIcs vol. 223, 2022. ([paper](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.CPM.2022.11))

## License
This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

The competitor indexes bundled as submodules under [external/](external/) retain their own licenses — see each submodule.
