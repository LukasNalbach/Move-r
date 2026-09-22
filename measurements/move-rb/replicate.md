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
CMake copies the patched submodule sources into the submodules on every configure run. The `cp`
line above is the one from the top-level `README.md` and is not needed here.

**Compiler.** GCC 12, 13, 14 and 15 all work. Use the compiler of your distribution. If you
select an older GCC than the one your `libtbb` was built with, linking fails with
`libtbb.so: undefined reference to __cxa_call_terminate@CXXABI_1.3.15`.

**Big-BWT.** The scripts build the indexes with `-c bigbwt`, and `columba-rlc` needs Big-BWT as
well. After cloning, check that the file `external/Big-BWT/makefile` exists and that `make`
created the program `bigbwt` in that directory. If `external/Big-BWT` is empty, run
`git submodule update --init external/Big-BWT`.

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

Place each input text in [`texts/`](texts/) under the name the scripts expect. All three texts
are on Zenodo ([10.5281/zenodo.22879117](https://doi.org/10.5281/zenodo.22879117)), compressed
with bsc:

| Text       | File name         | Download | Notes                                        |
| ---------- | ----------------- | -------- | -------------------------------------------- |
| SARS-CoV-2 | `sars2.ACGT.50Gi` | [`sars2.ACGT.50Gi.bsc`](https://zenodo.org/records/22879117/files/sars2.ACGT.50Gi.bsc?download=1) (151 MB) | DNA, reduced to the `A`,`C`,`G`,`T` alphabet |
| chr19      | `chr19.ACGT.50Gi` | [`chr19.ACGT.50Gi.bsc`](https://zenodo.org/records/22879117/files/chr19.ACGT.50Gi.bsc?download=1) (10.3 GB) | DNA, reduced to `A`,`C`,`G`,`T`              |
| dewiki     | `dewiki.50Gi`     | [`dewiki.50Gi.bsc`](https://zenodo.org/records/22879117/files/dewiki.50Gi.bsc?download=1) (214 MB) | byte alphabet (German Wikipedia)             |

Each file decompresses to a text of 50 GiB. Decompress it with
[bsc](https://github.com/IlyaGrebnov/libbsc):

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
§7), so unpacking them reproduces the paper's exact query sets. By default `measure-text.sh`
**regenerates** fresh patterns into `patterns/` / `patterns_ext/` on each run (calibrated to
*your* machine); unpack `patterns.7z` only if you want to inspect or reuse the original sets.

## 3. Run

**Everything, for all three texts:**

```shell
./measure-all.sh
```

`measure-all.sh` writes one build log per index. Combine them into the single file
`results-build.txt`, which `charts/construction.tex` reads (see §6):

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
     `columba-rlc`. `columba-rlc` is built in the four steps of upstream's
     `columba_build_pfp.sh`: `columba-rlc-build --preprocess`, `bigbwt` on the text, `bigbwt`
     on the reversed text, and `columba-rlc-build --pfp`. The build time of these three is
     measured with `/usr/bin/time -v`. Their peak heap memory is reported by the build tools
     themselves via `malloc_count`, as a `Peak memory usage during construction:` line in the
     build log. This is the same quantity that `move-rb` reports for its own construction.
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

The files `charts/*.tex` and `tables/*.tex` read `results-apm.txt`, `results-ext.txt` and
`results-build.txt` and produce the figures and tables of the paper from them. Section 6
describes how.

## 6. Regenerating the figures and tables

The figures and tables of the paper are generated from the `results-*.txt` files with
[sqlplot-tools](https://github.com/bingmann/sqlplot-tools). Every file in [`charts/`](charts/)
and [`tables/`](tables/) begins with a block of `% IMPORT-DATA` and `%% SELECT`,
`%% MULTIPLOT` or `%% TABULAR` lines. sqlplot-tools loads the results file named there into an
SQLite database, runs the query, and writes the resulting `\addplot` coordinates or table rows
into the same file, below the query block.

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

Each of these files already contains the numbers of the paper, so it compiles without running
sqlplot-tools first. The measurement data behind these numbers is in
[`results-paper/`](results-paper/).

**Building the figures.**

```shell
./make-plots.sh --paper     # use results-paper/, the data of the paper
./make-plots.sh             # use results/, your own measurements
```

`make-plots.sh` copies the charts, the tables and the results into the directory
`build-plots-paper/` or `build-plots/`, runs sqlplot-tools on every chart and table in that
directory, and compiles all of them into one PDF `plots.pdf`. Each figure and table is on a
page of its own, with the number it has in the paper written in the caption.

The files in `charts/` and `tables/` are not modified, so you can compare your numbers with the
numbers of the paper:

```shell
diff -u charts/apm.tex build-plots/charts/apm.tex
```

Running `./make-plots.sh --paper` produces files that are byte for byte identical to the ones
in `charts/` and `tables/`.

**What you need.** sqlplot-tools with the SQLite backend, either in your `PATH` or in the
variable `SQLPLOT_TOOLS`, and pdflatex with the packages `pgfplots`, `siunitx`, `subcaption`
and `caption`. On Ubuntu:

```shell
sudo apt install cmake libboost-all-dev libsqlite3-dev libpq-dev \
                 texlive-latex-recommended texlive-pictures texlive-science
git clone https://github.com/bingmann/sqlplot-tools.git
cd sqlplot-tools && mkdir build && cd build && cmake .. && make
export SQLPLOT_TOOLS=$PWD/src/sqlplot-tools
```

`libpq-dev` is in this list because the CMake file of sqlplot-tools searches for PostgreSQL
even if you only build the SQLite backend. If your sqlplot-tools has both backends, it tries
PostgreSQL first and then uses SQLite. It prints the line `Connection to PostgreSQL failed`
when it does this. That line can be ignored.

[`plot-styles.tex`](plot-styles.tex) contains the colors, the markers, the axis settings and
the macros that the charts and tables use. Include this file if you want to use one of the
charts in a different document.

**If you measure your own text.** Every chart and table selects the three texts of the paper by
name, for example `WHERE text = 'sars2.ACGT.50Gi'`. If you measured a different text, these
queries find no rows, and sqlplot-tools reports `MULTIPLOT() requires group column list`. In
that case `make-plots.sh` prints the name of the file, keeps the version with the numbers of
the paper, and continues with the next file. To plot your own text, replace the text names in
the `%%` query block of that file.

**The file `results-build.txt`.** `charts/construction.tex` (Figure 3.1) reads one file with
one `RESULT` line per index and text:

```
RESULT algo=build_move_rb_move text=sars2.ACGT.50Gi n=50000000001 time_build=22834000000000 peak_memory_usage=6099000000 size_index=3876616865
```

`algo` is `build_move_rb_move`, `build_move_rb_rlzsa`, `build_br_index`, `build_columba` or
`build_bmove`, where `build_bmove` is `columba-rlc`. `time_build` is in nanoseconds,
`peak_memory_usage` and `size_index` are in bytes.

`measure-text.sh` writes one build log per index, the five `results-build-*.txt` files listed
in §5. [`make-results-build.sh`](make-results-build.sh) reads those five logs and writes
`results-build.txt`:

```shell
./make-results-build.sh
```

`move-rb` and `move-rb-rlzsa` write every field into their `-m_idx` record, so the script only
copies them. For `br-index`, `columba` and `columba-rlc` it reads the build time from the
output of `/usr/bin/time -v`, the peak memory from the line
`Peak memory usage during construction:` that these tools print, `n` from the size of the input
text, and `size_index` from the size of the index files. If a build failed or its log is
missing, the script prints a warning and omits that line.

The data of Figure 3.1 of the paper is in
[`results-paper/results-build.txt`](results-paper/results-build.txt).

**The columns `sigma`, n/r and n/r_rev of Table 3.1.** `tables/texts.tex` reads these three
values from the `build_move_rb_move` line of `results-build.txt`. `move-rb-build -m_idx` writes
them as `sigma`, `r` and `r_rev`, and `make-results-build.sh` copies them into
`results-build.txt`. `sigma` counts the sentinel character as well, so the table prints
`sigma - 1`.

The values in `results-paper/results-build.txt` were read out of the index files that the
measurements of the paper were made with. For sars2 and chr19 the resulting n/r and n/r_rev are
0.03% higher than the values printed in Table 3.1 of the paper: 1000.93 instead of 1000.64, and
1119.46 instead of 1119.13. For dewiki the values are the same. The reason for this difference
is not known. The tables in this repository show the values from the index files.

## 7. Notes

- **Machine-dependent calibration.** `N` is calibrated to a target wall time on *your*
  machine, so the pattern counts and absolute throughputs will differ between systems.
- **Single-threaded.** All builds and queries run with one thread (`-p 1`).
- **Input format of columba.** `columba-build` and `columba-rlc-build` read the text with `-f`
  and check the file extension (`.fasta`, `.fa`, `.FASTA`, `.FA`, `.fna`, `.FNA`). A file with
  a different extension is rejected, no matter what it contains. On the first run,
  `measure-text.sh` writes a FASTA copy of the text to `texts/<text>.fa` and uses that copy in
  later runs. The copy needs as much disk space as the text. Neither tool has an option for
  the number of threads.
- **Re-runs.** `measure-all.sh` truncates `results-apm.txt`, `results-ext.txt` and the
  `results-build-*.txt` logs before it starts; `measure-text.sh` appends, so running it
  directly several times accumulates rows.
