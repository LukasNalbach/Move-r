# move-r
This is an optimized and parallelized implementation of the modified r-index described in [1] ([arxiv.org](https://arxiv.org/abs/2006.05104)).

## External Dependencies
- [OpenMP](https://www.openmp.org/)
- [intel TBB](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb.html)

## Included Dependencies
- [libsais](https://github.com/IlyaGrebnov/libsais)
- [abseil-cpp](https://github.com/abseil/abseil-cpp)
- [ips4o](https://github.com/ips4o/ips4o)
- [pasta::bit_vector](https://github.com/pasta-toolbox/bit_vector)
- [concurrentqueue](https://github.com/cameron314/concurrentqueue)
- [Big-BWT](https://github.com/alshai/Big-BWT)
- [sdsl-lite](https://github.com/simongog/sdsl-lite)

## CLI Build Instructions
This implementation has been tested on Ubuntu 20.04 with GCC 10.3.0, libtbb-dev, libomp-dev and libz-dev installed.
```
clone https://github.com/LukasNalbach/move-r.git
mkdir build
cd build
cmake ..
make
```
This creates five executeables in the build/cli/ folder:
- move-r-build
- move-r-count
- move-r-locate
- move-r-revert
- move-r-patterns

There is an explanation for each below.

## Usage in C++
### Cmake
```
add_subdirectory(move_r/)
set(MOVE_R_BUILD_CLI OFF)
set(MOVE_R_BUILD_BENCH_CLI OFF)
set(MOVE_R_BUILD_EXAMPLES OFF)
```

### C++
```
#include <move_r/move_r.hpp>

int main() {
    // build an index
    move_r<> index("This is a test string");

    // build a 64-bit index (intended for large input strings > UINT_MAX
    // bytes ~ 4GB) with revert support only, use the space-efficient
    // construction algorithm, use at most 8 threads and set the 
    // balancing parameter a to 4
    move_r<uint64_t> index_2("a large string",{revert},space,8,4);

    // print the number of bwt runs in the input string
    std::cout << index.num_bwt_runs() << std::endl;

    // print the index size
    std::cout << format_size(index.size_in_bytes()) << std::endl;

    // print the number of occurences of a pattern
    std::cout << index.count("test") << std::endl;

    // print all occurences of a pattern
    index.locate("is",[](auto o){std::cout << o << ", ";});

    // store all occurences of a pattern in a vector
    auto Occ = index.locate("test");

    std::cout << std::endl;
    for (auto o : Occ) {
        std::cout << o << ", ";
    }
}
```
There are more examples in the folder examples/, that show how to store, load and revert an index and access/retrieve the suffix array and the bwt. An example usage of the move data structure is also included.

## CLI-Usage
### move-r-build: builds move-r.
```
usage: move-r-build [options] <input_file>
   -c <mode>          construction mode: runtime or space (default: runtime)
   -o <base_name>     names the index file base_name.move-r (default: input_file)
   -s <op1> <op2> ... supported operations: revert, count and locate
                      (default: all operations)
   -p <integer>       number of threads to use during the construction of the index
                      (default: all threads)
   -a <integer>       balancing parameter; a must be an integer number and a >= 2 (default: 8)
   -m_idx <m_file>    m_file is file to write measurement data of the index construction to
   -m_mds <m_file>    m_file is file to write measurement data of the construction of the move
                      data structures to
   <input_file>       input file
```

### move-r-count: count all occurrences of the input patterns.
```
usage: move-r-count <index_file> <patterns_file>
   -m <m_file> <text_name>    m_file is the file to write measurement data to,
                              text_name should be the name of the original file
   <index_file>               index file (with extension .move-r)
   <patterns_file>            file in pizza&chili format containing the patterns.
```

### move-r-locate: locate all occurrences of the input patterns.
```
usage: move-r-locate [options] <index_file> <patterns>
   -c <input_file>            check correctness of each pattern occurrence on
                              this input file (must be the indexed input file)
   -m <m_file> <text_name>    m_file is the file to write measurement data to,
                              text_name should be the name of the original file
   -o <output_file>           write pattern occurrences to this file (ASCII)
   <index_file>               index file (with extension .move-r)
   <patterns_file>            file in pizza&chili format containing the patterns
```

### move-r-revert: reconstruct the original file from the index.
```
usage: move-r-revert [options] <index_file> <output_file>
   -im                        revert in memory; faster, but stores the whole
                              input in memory
   -p <integer>               number of threads to use while reverting
                              (default: greatest possible)
   -m <m_file> <text_name>    m_file is the file to write measurement data to,
                              text_name should be the name of the original file
   <index_file>               index file (with extension .move-r)
   <output_file>              output file
```

### move-r-patterns: generate patterns from a file.
```
usage: move-r-patterns <file> <length> <number> <patterns file> <forbidden>
       randomly extracts <number> substrings of length <length> from <file>,
       avoiding substrings containing characters in <forbidden>.
       The output file, <patterns file> has a first line of the form:
       # number=<number> length=<length> file=<file> forbidden=<forbidden>
       and then the <number> patterns come successively without any separator
```

### move-r-bench: benchmarks construction-, revert- and query-performance.
```
move-r-bench: benchmarks construction-, revert- and query-performance of move-r, r-index-f, rcomp, r-index
              (Prezza), r-index (Mun), OnlineRLBWT and DYNAMIC; has to be executed from the base folder.
usage: move-r-bench [options] <input_file> <patterns_file_1> <patterns_file_2> <num_threads>
   -c                 check for correctnes if possible; disables the -m option; will not print
                      runtime data if the runtime could be affected by checking for correctness
   -m <m_file>        writes measurement data to m_file
   <input_file>       input file
   <patterns_file_1>  file containing patterns (pattern length ~ number of occurrences) from <input_file>
                      to count and locate
   <patterns_file_2>  file containing patterns (pattern length << number of occurrences) from <input_file>
                      to locate
   <num_threads>      maximum number of threads to use
```
GitHub-repositories of the other indexes can be found in the folder `external/`.

## References
[1] Takaaki Nishimoto and Yasuo Tabei. Optimal-time queries on bwt-runs compressed indexes.
In 48th International Colloquium on Automata, Languages, and Programming (ICALP 2021),
volume 198, page 101. Schloss Dagstuhl–Leibniz-Zentrum für Informatik, 2021.
