# Move-r in C++

`move_r<support>` (header `move_r/move_r.hpp`) is the uni-directional index of a text T: it counts and locates the
occurrences of a pattern, gives access to the BWT and, depending on the support, to the suffix array, and reconstructs
T. Link your target with `move_r` (see the [main README](../README.md#c-usage)).

Every block below is a complete program. A comment `// -> x` at the end of a line or on the line above a statement
shows the line the statement prints, `// ~> x` an output that depends on the build (e.g. a size). The bi-directional
index is described in [move-rb.md](move-rb.md) and the move data structure in the
[README of move-data-structure](https://github.com/LukasNalbach/move-data-structure#usage-in-c).

## Contents

- [Construction](#construction)
- [Storing and loading](#storing-and-loading)
- [Counting and locating](#counting-and-locating)
- [Incremental search](#incremental-search)
- [Access to the BWT and the suffix array, reverting](#access-to-the-bwt-and-the-suffix-array-reverting)

## Construction

The template parameter `support` selects what the index stores besides M_LF (the move data structure for LF):

| support | queries | size |
| --- | --- | --- |
| `_count` | count | O(r) words |
| `_locate_one` | count, one occurrence per pattern | O(r) words |
| `_locate_move` (default) | count, locate (with M_Phi^{-1}) | O(r) words |
| `_locate_rlzsa` | count, locate, SA[i] (with an RLZ of the differential suffix array) | O(r) words + the RLZ |

The construction parameters (`move_r_params`) that matter most:

- `mode`: `_suffix_array` (default; builds the suffix array in memory, inputs of at most 2^38 - 1 symbols),
  `_suffix_array_space` (the same, stores some data structures on disk while they are not needed) or `_pfp` (builds
  the run-length encoded BWT with prefix-free parsing; less memory for repetitive inputs, up to 2^48 - 2 symbols);
- `num_threads`: the maximum number of threads (default: all);
- `a`: the balancing parameter (default 8); a smaller `a` makes M_LF larger and the queries faster;
- `reorder_alphabet`: `true` (default) maps the symbols by decreasing frequency, `false` keeps their order (then the
  suffix array is the one of the lexicographic order);
- `input_file`: builds the index from a file instead of the container passed to the constructor;
- `log`: prints the phases of the construction.

```cpp
#include <move_r/move_r.hpp>

int main()
{
    // from a string, with the default support _locate_move
    move_r<> index("abracadabra");
    std::cout << index.input_size() << std::endl; // -> 11
    std::cout << index.alphabet_size() << std::endl; // -> 5
    std::cout << index.num_bwt_runs() << std::endl; // -> 9
    std::cout << index.balancing_parameter() << std::endl; // -> 8
    std::cout << format_size(index.size_in_bytes()) << std::endl; // ~> 40 KB

    // count support only, built with the prefix-free parsing, at most 2
    // threads, balancing parameter 4
    std::string text;
    for (int i = 0; i < 100; i++)
        text += "number " + std::to_string(i % 10) + ", ";
    move_r<_count> index_count(text,
        { .mode = _pfp, .num_threads = 2, .a = 4 });
    std::cout << index_count.count("number 3") << std::endl; // -> 10

    // from a file (the string passed to the constructor only sets the symbol
    // type, one byte)
    std::ofstream("text.txt") << "abracadabra";
    std::string file = "text.txt";
    move_r<_locate_rlzsa> index_file(file, { .input_file = file });
    std::cout << index_file.count("abra") << std::endl; // -> 2

    // from a vector of integers (an integer alphabet; the symbols are mapped
    // to [1..sigma-1])
    std::vector<int32_t> values = { 7, -3, 1000000, 7, -3, 5 };
    move_r<> index_int(values);
    // -> 2
    std::cout << index_int.count(std::vector<int32_t>{ 7, -3 }) << std::endl;
    std::cout << index_int.byte_alphabet << std::endl; // -> 0
}
```

The second template parameter selects the encoding of the interval starting positions of M_LF: `EF` (default for
`_count` to `_locate_rlzsa`; smaller) or `POS` (faster), e.g. `move_r<_locate_move, POS>`.

## Storing and loading

```cpp
#include <filesystem>
#include <move_r/move_r.hpp>

int main()
{
    move_r<> index("abracadabra");

    // the file name the CLI uses: <base>.move-r-<support suffix>
    std::string file = move_r_index_file_name("abracadabra", _locate_move);
    std::cout << file << std::endl; // -> abracadabra.move-r-locate

    std::ofstream out(file);
    index.serialize(out); // or: index >> out;
    out.close();

    // read the support of a stored index without loading it
    std::ifstream in(file);
    // -> 1
    std::cout << (move_r<>::peek_support(in) == _locate_move) << std::endl;

    move_r<> loaded;
    loaded.load(in); // or: loaded << in;
    in.close();
    std::cout << loaded.count("a") << std::endl; // -> 5

    std::filesystem::remove(file);
}
```

## Counting and locating

`count` and `locate` take the pattern as a string, a `const char*` or any contiguous container of symbols of the
width of the input symbols. The occurrences are reported in the order of the suffix array.

```cpp
#include <move_r/move_r.hpp>

int main()
{
    move_r<> index("abracadabra");

    // number of occurrences
    std::cout << index.count("abra") << std::endl; // -> 2
    std::cout << index.count("abc") << std::endl; // -> 0

    // all occurrences as a vector
    for (auto occ : index.locate("abra")) std::cout << occ << " ";
    std::cout << std::endl; // -> 7 0

    // all occurrences reported as (i, SA[i]), where i is the position in the
    // suffix array
    index.locate(std::string("bra"), [](uint64_t i, uint64_t occ) {
        std::cout << occ << " ";
    });
    std::cout << std::endl; // -> 8 1

    // one occurrence per pattern (_locate_one stores only the suffix array
    // samples at the BWT run boundaries)
    move_r<_locate_one> index_one("abracadabra");
    auto query = index_one.query();
    std::string P = "cad";
    for (auto it = P.rbegin(); it != P.rend(); it++) query.prepend(*it);
    // -> 1 4
    std::cout << query.num_occ() << " " << query.one_occ() << std::endl;
}
```

## Incremental search

`query()` returns a search context that is extended to the left symbol by symbol (backward search).

```cpp
#include <move_r/move_r.hpp>

int main()
{
    move_r<_locate_rlzsa> index(std::vector<int32_t>{ 2, -1, 5, -1, 7, 2, -1 });
    auto query = index.query();

    // prepend returns whether the extended pattern occurs; if not, the context
    // is not changed
    // -> 1 3
    std::cout << query.prepend(-1) << " " << query.num_occ() << std::endl;
    // -> 1 2
    std::cout << query.prepend(2) << " " << query.num_occ() << std::endl;
    // -> 0 2
    std::cout << query.prepend(9) << " " << query.length() << std::endl;

    // suffix array interval [b, e] of the pattern (2, -1)
    auto [b, e] = query.sa_interval();
    std::cout << b << " " << e << std::endl; // -> 4 5

    // the occurrences one by one, then the remaining ones
    // -> 5 1
    std::cout << query.next_occ() << " " << query.num_occ_rem() << std::endl;
    for (auto occ : query.locate()) std::cout << occ << std::endl; // -> 0

    // the longest suffix of (0, 7, 2) that occurs: (7, 2)
    std::vector<int32_t> pattern = { 0, 7, 2 };
    query.reset();
    uint64_t start = pattern.size();
    while (start > 0 && query.prepend(pattern[start - 1])) start--;
    std::cout << start << std::endl; // -> 1
}
```

## Access to the BWT and the suffix array, reverting

The range methods take `retrieve_params`: the range `[l, r]` (all positions if `l > r`, the default) and the number of
threads. With more than one thread, the values are reported in no fixed order. `SA(i)` needs `_locate_rlzsa`; the
range methods for the suffix array need `_locate_move` or `_locate_rlzsa`. The BWT is returned with the terminator as 0
(`BWT(i)` returns it as `TERMINATOR_SYMBOL`). With `reorder_alphabet = false`, the suffixes are sorted
lexicographically. The index stores up to `move_r_params::revert_start_samples` (default 4096) BWT run starts with
their suffix array values; `revert` splits the text at these positions and reverts the pieces with up to the given
number of threads (default: all), each of which reverts `default_revert_strands` pieces alternately (2 with M_LF in
EF, else 1); `revert(report)` reports the symbols from right to left with one thread.
`BWT_range`, `SA_range` and `revert` return a `no_init_vector` (a `std::vector` that is resized without initializing
its elements, header `misc/memory.hpp` of string-utils); `BWT_range<out_t>` and `revert<out_t>` take another
container type.

```cpp
#include <filesystem>
#include <move_r/move_r.hpp>

int main()
{
    move_r<_locate_rlzsa> index("abracadabra", { .reorder_alphabet = false });

    // the BWT, as a vector and symbol by symbol
    for (uint8_t c : index.BWT_range()) std::cout << (c == 0 ? '$' : char(c));
    std::cout << std::endl; // -> ard$rcaaaabb
    index.BWT_range([](uint64_t i, uint64_t c) { std::cout << char(c); },
        { .l = 4, .r = 7, .num_threads = 1 });
    std::cout << std::endl; // -> rcaa
    // -> c 1
    std::cout << char(index.BWT(5)) << " "
              << (index.BWT(3) == TERMINATOR_SYMBOL) << std::endl;

    // the suffix array, as a vector and value by value
    for (auto s : index.SA_range()) std::cout << s << " ";
    std::cout << std::endl; // -> 11 10 7 0 3 5 8 1 4 6 9 2
    index.SA_range([](uint64_t i, uint64_t s) { std::cout << s << " "; },
        { .l = 2, .r = 4, .num_threads = 1 });
    std::cout << std::endl; // -> 7 0 3

    // random access to the suffix array (_locate_rlzsa)
    std::cout << index.SA(3) << std::endl; // -> 0
    index.SA_range(8, 10, [](uint64_t s) { std::cout << s << " "; });
    std::cout << std::endl; // -> 4 6 9

    // the text, as a vector (with 4 threads), from right to left symbol by
    // symbol, and into a file (with 2 threads)
    for (uint8_t c : index.revert(4)) std::cout << char(c);
    std::cout << std::endl; // -> abracadabra
    index.revert([](uint64_t i, uint64_t c) { std::cout << char(c); });
    std::cout << std::endl; // -> arbadacarba
    index.revert(std::string("reverted.txt"), 2);
    // -> 11
    std::cout << std::filesystem::file_size("reverted.txt") << std::endl;

    // the BWT (one byte per symbol) and the suffix array (8 bytes per value)
    // into files
    index.BWT_range(std::string("bwt.bin"));
    index.SA_range(std::string("sa.bin"));
    // -> 12 96
    std::cout << std::filesystem::file_size("bwt.bin") << " "
              << std::filesystem::file_size("sa.bin") << std::endl;
}
```
