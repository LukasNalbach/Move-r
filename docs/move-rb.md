# Move-rb in C++

`move_rb<support>` (header `move_rb/move_rb.hpp`) is the bi-directional index of a text T: a move-r index of T and one
of the reversed text T^R. A search context of a pattern P can be extended by a symbol on both sides. On top of this,
it matches patterns approximately with search schemes. Link your target with `move_r` (see the
[main README](../README.md#c-usage)).

Every block below is a complete program. A comment `// -> x` at the end of a line or on the line above a statement
shows the line the statement prints, `// ~> x` an output that depends on the build (e.g. a size).

## Contents

- [Construction](#construction)
- [Storing and loading](#storing-and-loading)
- [Search contexts: extension, locate](#search-contexts-extension-locate)
- [Enumerating extensions](#enumerating-extensions)
- [BWT, suffix array and text](#bwt-suffix-array-and-text)
- [Approximate pattern matching](#approximate-pattern-matching)

## Construction

The template parameter `support` is `_count`, `_locate_one` (one occurrence per pattern), `_locate_move` (default;
locate with M_Phi^{-1}) or `_locate_rlzsa` (locate with an RLZ of the differential suffix array).

The construction takes the parameters of move-r (`move_r_params`, see [move-r.md](move-r.md#construction)) and these:

- `context_table_share`: the share of the index size used for a table of the search contexts of all short patterns,
  at which a search starts (default 0.05, 0 builds none);
- `context_tables_single_length`: `true` stores the table for one pattern length only.

```cpp
#include <move_rb/move_rb.hpp>

int main()
{
    // from a string
    move_rb<> index("abracadabra");
    // -> 11 5
    std::cout << index.input_size() << " " << index.alphabet_size()
              << std::endl;
    std::cout << format_size(index.size_in_bytes()) << std::endl; // ~> 83 KB

    // with an RLZ of the differential suffix array, without context tables
    move_rb<_locate_rlzsa> index_rlzsa("abracadabra",
        { .context_table_share = 0 });
    std::cout << index_rlzsa.has_context_tables() << std::endl; // -> 0

    // from a file, with at most 2 threads
    std::ofstream("text.txt") << "abracadabra";
    std::string file = "text.txt";
    move_rb<_count> index_file(file, { .input_file = file, .num_threads = 2 });
    std::cout << index_file.num_bwt_runs() << std::endl; // -> 9

    // from a vector of integers
    move_rb<_locate_rlzsa> index_int(
        std::vector<uint32_t>{ 3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5 });
    bool occurs;
    auto ctx = index_int.search(std::vector<uint32_t>{ 3, 5 }, occurs);
    // -> 1 9
    std::cout << occurs << " " << ctx.locate_phase().locate()[0] << std::endl;
}
```

## Storing and loading

```cpp
#include <filesystem>
#include <move_rb/move_rb.hpp>

int main()
{
    move_rb<_locate_move> index("abracadabra");

    // the file name the CLI uses: <base>.move-rb<support suffix>
    std::string file = move_rb_index_file_name("abracadabra", _locate_move);
    // -> abracadabra.move-rb-locate
    std::cout << file << std::endl;

    std::ofstream out(file);
    index.serialize(out); // or: index >> out;
    out.close();

    // the support of a stored index
    std::ifstream in(file);
    std::cout << (read_move_rb_support(in) == _locate_move) << std::endl; // -> 1

    // load throws std::runtime_error if the index was built with another
    // support
    move_rb<_locate_move> loaded;
    loaded.load(in); // or: loaded << in;
    in.close();
    std::cout << loaded.input_size() << std::endl; // -> 11

    std::filesystem::remove(file);
}
```

## Search contexts: extension, locate

A search context stores the suffix array interval of P in T and of the reversed P in T^R. `empty_context<COUNT>()`
returns a context that only counts, `empty_context<LOCATE>()` one that can also locate. `extend(c, dir)` returns the
extended context and whether the extended pattern occurs (if not, the returned context is not valid).

```cpp
#include <move_rb/move_rb.hpp>

int main()
{
    using index_t = move_rb<_locate_move>;
    index_t index("abracadabra");

    // "c" -> "ca" -> "aca" -> "acad"
    auto [ctx_c, c_occurs] = index.empty_context<LOCATE>().extend('c', RIGHT);
    auto [ctx_ca, ca_occurs] = ctx_c.extend<RIGHT>('a');
    auto [ctx_aca, aca_occurs] = ctx_ca.extend('a', LEFT);
    auto [ctx_acad, occurs] = ctx_aca.extend('d', RIGHT);
    // -> 1 4 1
    std::cout << occurs << " " << ctx_acad.pattern_length() << " "
              << ctx_acad.num_occ() << std::endl;
    std::cout << std::get<1>(ctx_acad.extend('x', LEFT)) << std::endl; // -> 0

    // suffix array intervals of P in T and of the reversed P in T^R
    auto [b, e] = ctx_acad.forward_sa_interval();
    auto [b_R, e_R] = ctx_acad.backward_sa_interval();
    // ~> 6 6 5 5
    std::cout << b << " " << e << " " << b_R << " " << e_R << std::endl;

    // search a whole pattern from right to left (starting at the context
    // tables)
    bool found;
    auto ctx_a = index.search(std::string("a"), found);
    std::cout << found << " " << ctx_a.num_occ() << std::endl; // -> 1 5

    // locate: the occurrence of the sample of the context in O(1) time, then
    // all occurrences one by one
    auto loc = ctx_a.locate_phase();
    std::cout << (loc.sampled_occ() <= 10) << std::endl; // -> 1
    std::cout << loc.num_occ_rem() << std::endl; // -> 5
    uint64_t first = loc.next_occ();
    std::vector<uint48_t> rest = loc.locate();
    std::cout << first + 0 << " " << rest.size() << std::endl; // ~> 10 4

    // all occurrences reported with their suffix array positions
    ctx_a.locate_phase<index_t::_occ_and_sa_pos>().locate(
        [](uint64_t i, uint64_t occ) { std::cout << i << ":" << occ << " "; });
    std::cout << std::endl; // ~> 1:10 2:7 3:0 4:3 5:5

    // relations between contexts: "a" is a prefix of "acad" and a suffix of
    // "aca"
    // -> 1 1
    std::cout << ctx_a.is_prefix_of(ctx_acad) << " "
              << ctx_a.is_suffix_of(ctx_aca) << std::endl;

    // a counting context only counts; extend_one_sided maintains the interval
    // of one direction only
    auto [ctx_b, b_occurs] = index.empty_context<COUNT>().extend('b', RIGHT);
    auto [count_ctx, br_occurs] = ctx_b.extend('r', RIGHT);
    std::cout << count_ctx.num_occ() << std::endl; // -> 2
    auto [left_only, abr_occurs] = count_ctx.extend_one_sided<LEFT>('a');
    std::cout << left_only.num_occ_in<LEFT>() << std::endl; // -> 2
}
```

## Enumerating extensions

```cpp
#include <move_rb/move_rb.hpp>

int main()
{
    move_rb<> index("abracadabra");
    auto [ctx_a, occurs] = index.empty_context<COUNT>().extend('a', RIGHT);

    // all symbols c such that "ac" occurs (TERM: also the terminator, reported
    // as TERMINATOR_SYMBOL)
    ctx_a.for_each_extension<RIGHT, NO_TERM>([](uint64_t c, const auto& ext) {
        std::cout << char(c) << ext.num_occ() << " ";
        return true; // false stops the enumeration
    });
    std::cout << std::endl; // -> b2 c1 d1

    // the same with an extend context, here to the left: the symbols c such
    // that "ca" occurs
    auto ext = ctx_a.extend_phase();
    ext.prepare_extend_all<NO_TERM>(LEFT);

    while (ext.can_extend()) {
        char c = char(ext.next_symbol());
        std::cout << c << ext.extend_next().num_occ() << " ";
    }

    std::cout << std::endl; // -> r2 c1 d1
}
```

## BWT, suffix array and text

`forward_index()` returns the move-r index of T (see [move-r.md](move-r.md)) and `backward_index()` the one of T^R.

```cpp
#include <move_rb/move_rb.hpp>

int main()
{
    move_rb<_locate_move> index("abracadabra", { .reorder_alphabet = false });

    // the BWT and the suffix array of T and the text
    for (uint8_t c : index.forward_index().BWT_range())
        std::cout << (c == 0 ? '$' : char(c));
    std::cout << std::endl; // -> ard$rcaaaabb
    for (auto s : index.forward_index().SA_range()) std::cout << s << " ";
    std::cout << std::endl; // -> 11 10 7 0 3 5 8 1 4 6 9 2
    for (uint8_t c : index.forward_index().revert()) std::cout << char(c);
    std::cout << std::endl; // -> abracadabra
    // T^R, with 2 threads
    for (uint8_t c : index.backward_index().revert(2)) std::cout << char(c);
    std::cout << std::endl; // -> arbadacarba
}
```

## Approximate pattern matching

The index finds the occurrences of P with at most k errors under the Hamming distance (substitutions) or the edit
distance (substitutions, insertions and deletions) with a search scheme (header `misc/search_schemes.hpp`), which
splits P into parts and searches them in several orders. The built-in schemes for k errors are
`pigeon_hole_scheme(k)`, `suffix_filter_scheme(k)`, `zero_one_scheme(k)` and `min_u_scheme(k)` (usually the fastest);
`parse_search_scheme` reads one in the format of the files in [search_schemes/](../search_schemes/). An occurrence
(`aprx_occ_t`) stores its position, its length and its number of errors. Counting supports only the Hamming distance.

```cpp
#include <move_rb/move_rb.hpp>

int main()
{
    move_rb<> index("abracadabra");

    // count with at most 1 mismatch, and per number of mismatches
    // -> 2
    std::cout << index.count_hamming_dist("abrx", pigeon_hole_scheme(1))
              << std::endl;
    std::vector<uint64_t> hist = index.count_hamming_dist_histogram(
        std::string("abra"), min_u_scheme(1));
    std::cout << hist[0] << " " << hist[1] << std::endl; // -> 2 0

    // locate with at most 1 mismatch: (position, length, errors)
    for (auto occ : index.locate("abrx", min_u_scheme(1), HAMMING_DISTANCE))
        std::cout << occ << " ";
    std::cout << std::endl; // -> (0, 4, 1) (7, 4, 1)

    // locate with at most 1 error under the edit distance; a search scheme
    // from its text representation
    search_scheme_t scheme = parse_search_scheme(
        "p=2 k=1\n{0,1} {0,0} {0,1}\n{1,0} {0,0} {0,1}");
    for (auto occ : index.locate("acda", scheme, EDIT_DISTANCE))
        std::cout << occ << " ";
    std::cout << std::endl; // -> (3, 3, 1) (5, 3, 1)

    // uniform partitioning of P into the parts of the scheme instead of the
    // dynamic one
    apm_params_t uniform;
    uniform.partitioning = _partition_uniform;
    // -> 2
    std::cout << index.locate("abrx", zero_one_scheme(1), HAMMING_DISTANCE,
        uniform).size() << std::endl;

    // with a CIGAR per suffix array interval (= match, X mismatch,
    // I insertion, D deletion)
    std::vector<aprx_occ_t<CIGAR>> occs;
    auto cigars = index.locate<CIGAR>(std::string("acdab"),
        suffix_filter_scheme(1), EDIT_DISTANCE, occs);
    for (auto& occ : occs) {
        std::cout << occ.pos() << ":";
        for (auto& run : cigars[occ.cig_idx])
            std::cout << run.len << run.code();
        std::cout << " ";
    }
    std::cout << std::endl; // -> 3:2=1D3= 5:1=1I3=

    // the occurrences with the fewest errors (at most 3)
    std::vector<aprx_occ_t<>> best;
    uint64_t errors = locate_best(3, 0, best,
        [&](uint64_t k, std::vector<aprx_occ_t<>>& occ) {
            index.locate(std::string("dacab"), min_u_scheme(k),
                EDIT_DISTANCE, occ);
        });
    std::cout << errors << " " << best.size() << std::endl; // ~> 2 5
}
```

An error model (header `algorithms/move_rb/apm/model.hpp`) replaces the unit costs: costs of substitutions,
insertions, deletions and gap openings per range of P and per symbol, bounds on the number of operations in a range,
symbol classes and IUPAC codes. It is compiled for a pattern; the cost of an occurrence is reported as its error
count.

```cpp
#include <move_rb/move_rb.hpp>

int main()
{
    move_rb<> index("abracadabra");

    // an insertion or a deletion costs 2, the last two symbols of P have no
    // errors, an occurrence costs at most 2
    apm_model_t model;
    model.costs.push_back(parse_apm_cost_rule("id=2"));
    model.costs.push_back(parse_apm_cost_rule("sid=inf@-2:"));
    std::string P = "acxda";
    apm_costs_t costs(model, P, "", false, 2, EDIT_DISTANCE,
        index.forward_index().map_int());

    apm_params_t params;
    params.costs = &costs;
    auto occs = index.locate(P, min_u_scheme(costs.max_errors()), EDIT_DISTANCE,
        params);
    for (auto occ : occs) std::cout << occ << " ";
    std::cout << std::endl; // -> (3, 5, 1)
}
```
