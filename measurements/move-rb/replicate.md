# Reproducing the Move-rb measurements

This directory reproduces the construction, approximate-pattern-matching (APM) and raw
extension/enumeration measurements of the Move-rb paper for `move-rb`, `move-rb-rlzsa`,
`br-index`, `columba` and `columba-rlc` (b-move).

You provide the input texts in [`texts/`](texts/); the scripts build every index and write
all measurements to [`results/`](results/).

## 1. Prerequisites

**Build the project.** Compile Move-r once, following the top-level `README.md`
(tested on Ubuntu 24.04, GCC 13.3.0, with `libtbb-dev`, `libomp-dev`, `python3-psutil`,
`libz-dev`):

```shell
# from the repository root
mkdir build && cd build
cmake ..
cp -rf ../patched-files/* ..
make
```

This produces the CLI tools in `build/cli/` and the benchmark tools in `build/bench/`.
The reproduction scripts locate them automatically at `../../build` relative to this folder.

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

Place each input text in [`texts/`](texts/) under the name the scripts expect:

| Text       | File name         | Download | Notes                                        |
| ---------- | ----------------- | -------- | -------------------------------------------- |
| SARS-CoV-2 | `sars2.ACGT.50Gi` | [Google Drive](https://drive.google.com/file/d/1PYwoLMdOjneMyVv4-B2rB0hRF9FiEAZU/view?usp=sharing) | DNA, reduced to the `A`,`C`,`G`,`T` alphabet |
| chr19      | `chr19.ACGT.50Gi` | [Google Drive](https://drive.google.com/file/d/1GzKUzfoh3i89M-1CuUiEvjwCWBKjHd6M/view?usp=sharing) | DNA, reduced to `A`,`C`,`G`,`T`              |
| dewiki     | `dewiki.50Gi`     | [Google Drive](https://drive.google.com/file/d/1y1ajO7JI3QVQAg_ud1rXkgRgihy0laQl/view?usp=sharing) | byte alphabet (German Wikipedia)             |

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
     `columba-rlc` (`columba-rlc-build -p`, prefix-free/Big-BWT). Their wall time is captured
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

`results-apm.txt` and `results-ext.txt` are the inputs to the paper's figures and tables
(the `charts/*.tex` and `tables/*.tex` files import them via `sqlplot-tools`).

## 6. Notes

- **Machine-dependent calibration.** `N` is calibrated to a target wall time on *your*
  machine, so the pattern counts and absolute throughputs will differ between systems.
- **Single-threaded.** All builds and queries run with one thread (`-p 1`).
- **columba input format.** `columba-build` / `columba-rlc-build` take their reference via
  `-f`, which expects a (multi-)FASTA file. If a plain text is rejected, wrap it as a
  single-record FASTA (one header line, then the sequence) before running, or adjust the
  `columba-*-build` invocation in `measure-text.sh`.
- **Re-runs.** `measure-all.sh` truncates `results-apm.txt`, `results-ext.txt` and the
  `results-build-*.txt` logs before it starts; `measure-text.sh` appends, so running it
  directly several times accumulates rows.
