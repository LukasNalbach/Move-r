# Move-r
This is an optimized and parallelized implementation of the modified r-index described in [1] ([arxiv.org](https://arxiv.org/abs/2006.05104)).

## External Dependencies
- [OpenMP](https://www.openmp.org/)
- [intel TBB](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onetbb.html)

## Included Dependencies
- [libsais](https://github.com/IlyaGrebnov/libsais)
- [abseil-cpp](https://github.com/abseil/abseil-cpp)
- [ips4o](https://github.com/ips4o/ips4o)
- [concurrentqueue](https://github.com/cameron314/concurrentqueue)
- [Big-BWT](https://gitlab.com/manzai/Big-BWT)
- [sdsl-lite](https://github.com/simongog/sdsl-lite)

## CLI Build Instructions
This implementation has been tested on Ubuntu 22.04 with GCC 11.4.0, libtbb-dev, libomp-dev, python3-psutil and libz-dev installed.
```shell
clone https://github.com/LukasNalbach/move-r.git
mkdir build
cd build
cmake ..
cp -rf ../patched-files/* ..
make
```
This creates six executeables in the build/cli/ folder:
- move-r-build
- move-r-count
- move-r-locate
- move-r-revert
- move-r-patterns
- move-r-bench

There is an explanation for each below.

## Usage in C++
### Cmake
```cmake
add_subdirectory(move_r/)
set(MOVE_R_BUILD_CLI OFF)
set(MOVE_R_BUILD_BENCH_CLI OFF)
set(MOVE_R_BUILD_EXAMPLES OFF)
set(MOVE_R_BUILD_TESTS OFF)
```

### C++

#### Move-r
```c++
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

#### Move Data Structure
```c++
#include <move_r/misc/utils.hpp>
#include <move_r/data_structures/move_data_structure/move_data_structure.hpp>
#include <move_r/data_structures/move_data_structure/move_data_structure_str.hpp>

int main() {
   // Build a move data structure from the disjoint interval
   // sequence I = (0,1),(1,0); => n = 2
   move_data_structure<> mds({{0,1},{1,0}},2);

   // create a pair to perform move queries with
   std::pair<uint32_t,uint32_t> ix{0,0};

   // perform some move queries
   std::cout << to_string<>(ix = mds.move(ix)) << std::endl;
   std::cout << to_string<>(ix = mds.move(ix)) << std::endl;

   // build a move_data_structure_str (intended for I_LF);
   // this move data structure additionally stores a string interleaved
   // with the arrays needed for performing move queries (intended for
   // storing the characters of the bwt (sub-)runs);

   // use at most 4 threads and set a := 2
   move_data_structure_str<> mds_str({{0,4},{1,5},{2,6},{3,7},{4,0}},8,4,2);

   // this disjoint interval sequence is not 2-balanced, because the output
   // interval [0,3] contains 4 >= 2a = 4 input intervals

   // print the pairs of the resulting disjoint interval sequence
   for (uint32_t i=0; i<mds_str.num_intervals(); i++) {
      std::cout << to_string<>(std::make_pair(mds_str.p(i),mds_str.q(i)));
   }
   
   // the balancing algorithm has added the pair (6,2) to balance the sequence
}
```

#### Move-r Store & Load
```c++
#include <move_r/move_r.hpp>

int main() {
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

   // load the same index into another move_r-object
   // but only with revert support
   index_ifile.open("test_idx.move-r");
   move_r<> reloaded_index_2;
   reloaded_index_2.load(index_ifile,{revert});
   index_ifile.close();
}
```

#### Move-r Retrieval
```c++
#include <move_r/move_r.hpp>

int main() {
   // build an index
   move_r<> index("This is a test string");

   // retrieve the range [8,17] of the original text and store 
   // it in a string using at most 2 threads
   std::string reverted_range = index.revert_range(8,17,2);
   
   std::cout << std::endl;
   for (auto c : reverted_range) {
      std::cout << c;
   }

   // print the original text from right to left without storing it
   // using 1 thread
   std::cout << std::endl;
   index.revert_range([](auto,auto c){std::cout << c;},0,20,1);

   // retrieve the suffix array values in the range [2,6] using at
   // most 4 threads and store them in a vector
   std::vector<uint32_t> sa_range = index.retrieve_sa_range(2,6,4);
   std::cout << std::endl;
   
   for (auto s : sa_range) {
      std::cout << s << ", ";
   }

   // print the suffix array values in the range [7,14] from right to left
   // using 1 thread
   std::cout << std::endl;
   index.retrieve_sa_range([](auto,auto s){std::cout << s << ", ";},7,14,1);
}
```

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
move-r-bench: benchmarks construction- and query performance of move-r, block-rlbwt-2, block-rlbwt-v,
              block-rlbwt-r, r-index, r-index-f, rcomp-glfig and online-rlbwt;
              has to be executed from the base folder.
usage 1: move-r-bench [options] <input_file> <patterns_file_1> <patterns_file_2>
   -c                 check for correctnes if possible; disables the -m option; will not print
                      runtime data if the runtime could be affected by checking for correctness
   -m <m_file>        writes measurement data to m_file
   <input_file>       input file
   <patterns_file_1>  file containing patterns (pattern length ~ number of occurrences) from <input_file>
                      to count and locate
   <patterns_file_2>  file containing patterns (pattern length << number of occurrences) from <input_file>
                      to locate
usage 2: move-r-bench -a [options] <input_file> <patterns_file_1> <patterns_file_2> <num_threads>
                      constructs move_r using <num_threads> threads and measures count- and locate
                      performance of move_r for a=2, a=4, ..., a=8192.
   -m <m_file>        writes measurement data to m_file
   <input_file>       input file
   <patterns_file_1>  file containing patterns (pattern length ~ number of occurrences) from <input_file>
                      to count and locate
   <patterns_file_2>  file containing patterns (pattern length << number of occurrences) from <input_file>
                      to locate
   <num_threads>   maximum number of threads to use
```

#### GitHub repositories of the other indexes
- [block-rlbwt](https://github.com/saskeli/block_RLBWT/tree/main)
- [r-index](https://github.com/alshai/r-index/tree/master)
- [r-index-f](https://github.com/drnatebrown/r-index-f/tree/master)
- [rcomp-glfig](https://github.com/kampersanda/rcomp/tree/main)
- [online-rlbwt](https://github.com/itomomoti/OnlineRlbwt/tree/master)

## References
[1] Takaaki Nishimoto and Yasuo Tabei. Optimal-time queries on bwt-runs compressed indexes.
In 48th International Colloquium on Automata, Languages, and Programming (ICALP 2021),
volume 198, page 101. Schloss Dagstuhl–Leibniz-Zentrum für Informatik, 2021.