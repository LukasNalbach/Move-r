# Reproducing the Move-rb measurements

This directory reproduces the construction, approximate-pattern-matching (APM) and raw
extension/enumeration measurements of the Move-rb paper for `move-rb`, `move-rb-rlzsa`,
`br-index`, `columba` and `columba-rlc` (b-move).

You provide the input texts in [`texts/`](texts/); the scripts build every index and write
all measurements to [`results/`](results/).

## 1. Prerequisites

**Build the project.** Compile Move-r once, following the top-level `README.md`
(tested on Ubuntu 24.04, GCC 13.3.0 and 14.2.0, with `libtbb-dev`, `libomp-dev`,
`python3-psutil`, `libz-dev`):

```shell
# from the repository root
mkdir build && cd build
cmake ..
cp -rf ../patched-files/* ..
make
```

This produces the CLI tools in `build/cli/` and the benchmark tools in `build/bench/`.
The reproduction scripts locate them automatically at `../../build` relative to this folder.
The `cp` step is only kept because the top-level `README.md` lists it; CMake applies the patched
submodule sources itself on every configure, so running it changes nothing.

**Compiler.** GCC 12 through 15 all work. Mixing versions does not: building with an older GCC
than the one your distribution's `libtbb` was compiled with fails to link with
`libtbb.so: undefined reference to __cxa_call_terminate@CXXABI_1.3.15`. Use the system compiler
unless you have a reason not to.

**Big-BWT.** `measure-text.sh` builds with `-c bigbwt`, and `columba-rlc` needs it too, so the
bundled Big-BWT in `external/Big-BWT` has to be there. Check that `external/Big-BWT/makefile`
exists after cloning and that `make` produced the `bigbwt` driver next to it.

**GNU time.** The competitor builds are wrapped in `/usr/bin/time -v` for wall time. Their
peak *heap* memory is measured with `malloc_count` (each build tool reports it directly,
consistent with how `move-rb` measures its own construction peak): `columba` /
`columba-rlc` and `bri-build` all link `malloc_count` and print a
`Peak memory usage during construction:` line, which is captured in the build log alongside
GNU time's resident-set size. `/usr/bin/time` is available by default on Ubuntu; the
scripts assume it exists.

**Memory / scale.** The paper uses 50 GB texts on a machine with 1 TB of RAM. `br-index`
and `columba` construct their index in RAM (`br-index` needs up to ~0.8 TB on the 50 GB
inputs), whereas `move-rb` streams intermediate data to disk. For a quick functional test,
use a small text; full-scale reproduction needs a large machine and free disk space.

## 2. Provide the texts and patterns

Place each input text in [`texts/`](texts/) under the name the scripts expect. All three are
archived on Zenodo together with this code, BSC-compressed
([10.5281/zenodo.22879117](https://doi.org/10.5281/zenodo.22879117)):

| Text       | File name         | Download | Notes                                        |
| ---------- | ----------------- | -------- | -------------------------------------------- |
| SARS-CoV-2 | `sars2.ACGT.50Gi` | [`sars2.ACGT.50Gi.bsc`](https://zenodo.org/records/22879117/files/sars2.ACGT.50Gi.bsc?download=1) (151 MB) | DNA, reduced to the `A`,`C`,`G`,`T` alphabet |
| chr19      | `chr19.ACGT.50Gi` | [`chr19.ACGT.50Gi.bsc`](https://zenodo.org/records/22879117/files/chr19.ACGT.50Gi.bsc?download=1) (10.3 GB) | DNA, reduced to `A`,`C`,`G`,`T`              |
| dewiki     | `dewiki.50Gi`     | [`dewiki.50Gi.bsc`](https://zenodo.org/records/22879117/files/dewiki.50Gi.bsc?download=1) (214 MB) | byte alphabet (German Wikipedia)             |

Each file decompresses to a 50 GiB text with [bsc](https://github.com/IlyaGrebnov/libbsc):

```shell
bsc d sars2.ACGT.50Gi.bsc texts/sars2.ACGT.50Gi
```

To reproduce with your own inputs, drop them into `texts/` and call `measure-text.sh`
directly (see below).

**Patterns.** The exact query sets used in the paper ship with this repository, packed in
[`patterns.7z`](patterns.7z) — the APM queries in `patterns/` and the extension queries in
`patterns_ext/`. Unpack them with:

```shell
7z x patterns.7z
```

These are machine-calibrated (their pattern counts were tuned on the paper's machine, see
§6), so unpacking them reproduces the paper's exact query sets. By default `measure-text.sh`
**regenerates** fresh patterns into `patterns/` / `patterns_ext/` on each run (calibrated to
*your* machine); unpack `patterns.7z` only if you want to inspect or reuse the original sets.

## 3. Run

**Everything, for all three texts:**

```shell
./measure-all.sh
```

and then, to turn the five per-index build logs into the single `results-build.txt`
that `charts/construction.tex` reads (see §6):

```shell
./make-results-build.sh
```

`measure-all.sh` clears the query-result files and calls `measure-text.sh` once per text
with the per-text parameters used in the paper:

| Text            | `-T` (calibration target, s) | `-M` (min patterns) | `-C` (columba) |
| --------------- | ---------------------------- | ------------------- | -------------- |
| sars2, chr19    | 5                            | 4                   | 1 (on)         |
| dewiki          | 1                            | 2                   | 0 (off)        |

`columba` and `columba-rlc` only support DNA, so they are skipped for `dewiki` (`-C 0`).

**A single text**, e.g. only chr19:

```shell
./measure-text.sh -t chr19.ACGT.50Gi -T 5 -M 4 -C 1
```

`measure-text.sh` options:

| Flag | Meaning                                                              | Default |
| ---- | ------------------------------------------------------------------- | ------- |
| `-t` | text file name inside `texts/` (required)                           | —       |
| `-T` | `move-rb-gen-*-queries --time`: per-set calibration target in seconds | 5       |
| `-M` | `move-rb-gen-*-queries --min`: minimum patterns per set             | 4       |
| `-C` | `1` = also build/measure columba & columba-rlc, `0` = skip them      | 1       |
| `-p` | build threads (the paper builds and queries single-threaded)         | 1       |

## 4. What each run does

For every text, `measure-text.sh`:

1. **Builds all indexes.**
   - `move-rb` and `move-rb-rlzsa` (`move-rb-build -c bigbwt`), which write their own
     construction metrics via `-m_idx`.
   - `br-index` (`bri-build -divsufsort`, as in the paper), `columba` (`columba-build`) and
     `columba-rlc` (four steps, as in upstream's `columba_build_pfp.sh`:
     `columba-rlc-build --preprocess`, `bigbwt` on the text, `bigbwt` on its reverse, then
     `columba-rlc-build --pfp`). Their wall time is captured
     with `/usr/bin/time -v`; their peak heap memory is reported by the build tools themselves
     via `malloc_count` (a `Peak memory usage during construction:` line in the build log),
     consistent with the `move-rb` construction peak.
2. **APM measurements.** `move-rb-gen-apm-queries` samples random substring patterns and
   auto-calibrates the pattern count `N` per `(k, m)` set so that `move-rb-rlzsa` runs for
   about `-T` seconds (but at least `-M` patterns). `move-rb-bench-apm` then benchmarks
   count and locate (`--cigar both --algo both`), replaying each set until at least 10 s
   elapse and reporting the average time per run.
3. **Extension measurements.** `move-rb-gen-ext-queries` and `move-rb-bench-ext` do the same
   for exact bidirectional extension and SA-interval enumeration.

Competitor **builds are non-fatal**: if one fails, the script prints a warning, skips that
index (the benchmarks load each index independently and simply omit missing ones), and
continues. The `move-rb`/`move-rb-rlzsa` builds are required and abort the run on failure.

## 5. Output

All results land in [`results/`](results/):

| File                              | Contents                                            |
| --------------------------------- | --------------------------------------------------- |
| `results-build-move-rb.txt`       | move-rb construction metrics (`-m_idx`)             |
| `results-build-move-rb-rlzsa.txt` | move-rb-rlzsa construction metrics (`-m_idx`)        |
| `results-build-br-index.txt`      | br-index build (GNU time + `malloc_count` peak)     |
| `results-build-columba.txt`       | columba build (GNU time + `malloc_count` peak)      |
| `results-build-columba-rlc.txt`   | columba-rlc (b-move) build (GNU time + `malloc_count` peak) |
| `results-apm.txt`                 | APM count/locate throughput (`RESULT` lines)        |
| `results-ext.txt`                 | raw extension + enumeration throughput (`RESULT`)   |
| `results-build.txt`               | the five build logs, assembled by `make-results-build.sh` |

`results-apm.txt`, `results-ext.txt` and `results-build.txt` are the inputs to the paper's
figures and tables (the `charts/*.tex` and `tables/*.tex` files import them via
`sqlplot-tools`) — see §6 below.

## 6. Regenerating the figures and tables

Every figure and table of the paper is generated from the `results-*.txt` files by
[`sqlplot-tools`](https://github.com/bingmann/sqlplot-tools): each file in
[`charts/`](charts/) and [`tables/`](tables/) carries its `% IMPORT-DATA` and
`%% SELECT` / `%% MULTIPLOT` / `%% TABULAR` directives at the top, and sqlplot-tools
rewrites the `\addplot` coordinates resp. the tabular rows below them in place.

| File                                  | Paper float     | Reads                                |
| ------------------------------------- | --------------- | ------------------------------------ |
| `tables/texts.tex`                    | Table 3.1       | `results-ext.txt`, `results-apm.txt` |
| `charts/construction.tex`             | Figure 3.1      | `results-build.txt`                  |
| `charts/apm.tex`                      | Figure 3.2      | `results-apm.txt`                    |
| `tables/mem_factor.tex`               | Table C.1       | `results-apm.txt`                    |
| `charts/raw_performance.tex`          | Figure D.1      | `results-ext.txt`                    |
| `charts/apm-same-alg.tex`             | Figure E.1      | `results-apm.txt`                    |
| `charts/apm_native_k{4,7,10,13}.tex`  | Figures F.1–F.4 | `results-apm.txt`                    |
| `charts/apm_samealg_k{4,7,10,13}.tex` | Figures F.5–F.8 | `results-apm.txt`                    |
| `tables/patterns_grid.tex`            | Tables G.1–G.5  | `results-ext.txt`, `results-apm.txt` |

The files are committed with the paper's numbers already substituted, so they compile as
they are. [`results-paper/`](results-paper/) holds the exact measurement data behind them.

**Build them.**

```shell
./make-plots.sh --paper     # rebuild the paper's figures from results-paper/
./make-plots.sh             # ... or from your own results/
```

`make-plots.sh` stages everything in `build-plots-paper/` resp. `build-plots/`, runs
sqlplot-tools over every chart and table there and builds the result into a single PDF
(`plots.pdf`) via [`plots.tex`](plots.tex) — one float per paper figure/table, captioned
with the number it carries in the paper. Nothing under `charts/` or `tables/` is modified,
so you can diff your numbers against the paper's:

```shell
diff -u charts/apm.tex build-plots/charts/apm.tex
```

It needs `sqlplot-tools` built with the SQLite backend (on your `PATH`, or pointed at by
`$SQLPLOT_TOOLS`) and a `pdflatex` with `pgfplots`, `siunitx`, `subcaption` and `caption`.
On Ubuntu:

```shell
sudo apt install cmake libboost-all-dev libsqlite3-dev libpq-dev \
                 texlive-latex-recommended texlive-pictures texlive-science
git clone https://github.com/bingmann/sqlplot-tools.git
cd sqlplot-tools && mkdir build && cd build && cmake .. && make
export SQLPLOT_TOOLS=$PWD/src/sqlplot-tools
```

(`libpq-dev` is needed even for the SQLite build — sqlplot-tools' CMake looks for
PostgreSQL unconditionally.)

[`plot-styles.tex`](plot-styles.tex) holds the marker and axis styles and the notation
macros the charts and tables use, extracted from the paper's preamble; include it if you
want to embed one of the charts in a document of your own.

`./make-plots.sh --paper` reproduces the committed `charts/` and `tables/` files byte for byte.

**Your own texts.** Every chart and table selects the paper's three texts by name
(`WHERE text = 'sars2.ACGT.50Gi'`, `'chr19.ACGT.50Gi'`, `'dewiki.50Gi'`). Measuring a text of
your own therefore leaves those queries with no rows, and sqlplot-tools reports
`MULTIPLOT() requires group column list`. `make-plots.sh` keeps the paper's version of such a
file, names it and carries on, so the rest of the PDF is still yours. To plot your own text,
replace the text names in the `%%` query block of the file you care about.

**sqlplot-tools and PostgreSQL.** A sqlplot-tools built with both backends tries PostgreSQL
first and falls back to SQLite3 on its own; the `Connection to PostgreSQL failed` line it
prints on the way is harmless.

**`results-build.txt`.** `charts/construction.tex` reads a single
`results/results-build.txt` holding one `RESULT` line per (index, text):

```
RESULT algo=build_move_rb_move text=sars2.ACGT.50Gi n=50000000001 time_build=22834000000000 peak_memory_usage=6099000000 size_index=3876616865
```

`algo` is one of `build_move_rb_move`, `build_move_rb_rlzsa`, `build_br_index`,
`build_columba` or `build_bmove` (= `columba-rlc`), `time_build` is in nanoseconds and
`peak_memory_usage` / `size_index` are in bytes.

`measure-text.sh` writes the five raw `results-build-*.txt` logs listed in §5;
[`make-results-build.sh`](make-results-build.sh) turns them into `results-build.txt`:

```shell
./make-results-build.sh
```

`move-rb` and `move-rb-rlzsa` report every field directly in their `-m_idx` record. For
`br-index`, `columba` and `columba-rlc` the fields are read out of the build log: the wall
time from `/usr/bin/time -v`, the peak heap from the `Peak memory usage during
construction:` line the build tools print via `malloc_count`, `n` from the input text and
`size_index` from the index files on disk. Builds that failed or are missing are skipped
with a warning, so a partial run still yields a usable chart.

The data behind the paper's Figure 3.1 is in
[`results-paper/results-build.txt`](results-paper/results-build.txt).

**`sigma`, `r` and `r_rev`.** `tables/texts.tex` takes the alphabet size and the two
compression rates $n/r$ and $n/\bwd{r}$ from the `build_move_rb_move` row of
`results-build.txt`, where `make-results-build.sh` carries them over from the `-m_idx` record.
They used to be literals in the query, which meant the column kept showing the paper's texts
whatever you had measured. `sigma` counts the sentinel, so the table prints `sigma - 1`.

The values in `results-paper/` were read back out of the index files the paper was measured
with. For sars2 and chr19 they give compression rates 0.03% above the ones printed in the
paper (1000.93 instead of 1000.64, and 1119.46 instead of 1119.13); dewiki matches exactly.
The two DNA indexes were rebuilt after the table had been made, which fits the pattern -- the
offset is the same in both columns of a text. The tables here follow the indexes.

## 7. Notes

- **Machine-dependent calibration.** `N` is calibrated to a target wall time on *your*
  machine, so the pattern counts and absolute throughputs will differ between systems.
- **Single-threaded.** All builds and queries run with one thread (`-p 1`).
- **columba input format.** `columba-build` / `columba-rlc-build` take their reference via
  `-f`, which validates the file *extension* (`.fasta`, `.fa`, `.FASTA`, `.FA`, `.fna`,
  `.FNA`), so a plain text is rejected whatever its content. `measure-text.sh` therefore
  writes a single-record FASTA copy next to the text (`texts/<text>.fa`) on the first run and
  reuses it afterwards; budget as much free disk as the text itself. They also do not accept
  a `-t` thread option.
- **Re-runs.** `measure-all.sh` truncates `results-apm.txt`, `results-ext.txt` and the
  `results-build-*.txt` logs before it starts; `measure-text.sh` appends, so running it
  directly several times accumulates rows.
