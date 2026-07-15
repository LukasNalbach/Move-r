# Reproducing the Move-r measurements

This directory reproduces the construction and query (count / locate) measurements of the
Move-r paper, comparing `move-r` and `move-r-rlzsa` against the competitor r-indexes
(`block-rlbwt`, `r-index`, `r-index-f`, `rcomp-glfig`, `online-rlbwt`).

You provide the input texts in [`texts/`](texts/); the two pattern sets per text ship with
this repository, packed in [`patterns.7z`](patterns.7z). The scripts build the indexes and
write all measurements to [`results/`](results/).

## 1. Prerequisites

**Build the project** with the move-r benchmark enabled (it is `ON` by default), following
the top-level `README.md`:

```shell
# from the repository root
mkdir build && cd build
cmake .. -DMOVE_R_BUILD_MOVE_R_BENCH=ON
cp -rf ../patched-files/* ..
make
```

This produces the CLI tools (`move-r-build`, `move-r-locate`, `gen-patterns`) in
`build/cli/` and the benchmark tool `move-r-bench` in `build/bench/`. The reproduction
scripts locate them automatically at `../../build` relative to this folder.

## 2. Provide the texts and unpack the patterns

Download and decompress the texts, then place each into [`texts/`](texts/):

- [einstein.en.txt](https://pizzachili.dcc.uchile.cl/repcorpus/real/einstein.en.txt.7z)
- [english](https://pizzachili.dcc.uchile.cl/texts/nlang/english.gz)
- [chr19](https://drive.google.com/file/d/1GrCHHcc3zH56Q-c6WbI1N6qOh0sBD5DO/view?usp=sharing)
- [dewiki](https://drive.google.com/file/d/1GqvkN0FH6dkSxHZCXFPOr7I1iOBenUIZ/view?usp=sharing)
- [sars2](https://drive.google.com/file/d/134fLOpY1_3dFTdSSc_vW2qai4yKuyl2W/view?usp=sharing)

The two pattern sets per text ship with this repository, packed in
[`patterns.7z`](patterns.7z). Unpack them into `patterns/` before running:

```shell
7z x patterns.7z
```

This creates `patterns/<text>-patterns-bal` (pattern length ≈ number of occurrences) and
`patterns/<text>-patterns-phi` (pattern length ≪ number of occurrences). To measure other
texts, generate the two sets with `gen-patterns` (in `build/cli/`) using the same naming
scheme and drop them into `patterns/`.

## 3. Run

The scripts operate inside this directory (`measurements/move-r/`) regardless of where you
invoke them from.

**All texts:**

```shell
./measure-all-texts.sh -p <max_threads>
```

**A single text**, e.g. only einstein.en.txt:

```shell
./measure-text.sh -t einstein.en.txt -p <max_threads>
```

`<max_threads>` is the maximum number of threads to use (default 1; the paper builds and
queries single-threaded).

## 4. What each run does

For every text, `measure-text.sh`:

1. **Comparison benchmark** (`move-r-bench`): builds move-r, move-r-rlzsa and each
   competitor index and measures count/locate throughput on both pattern sets, writing
   `RESULT` lines to `results/results.txt`. `move-r-bench` builds some competitors via
   paths relative to the repository root, so the script runs it from there.
2. **Balancing-parameter sweep** (`move-r-bench -a`): count/locate performance of move-r
   for a range of thread counts and balancing parameters → `results/results-a.txt`.
3. **Construction phases** (`move-r-build`, a=8, threads 1..p): per-phase construction
   metrics via `-m_idx` / `-m_mds` → `results/results-move-r-build-phases.txt` and
   `results/results-mds-build-phases.txt`.
4. **Locate statistics** (`move-r-build -a 2` then `move-r-locate -m`) → `results/statistics.txt`.

## 5. Output

All results land in [`results/`](results/): `results.txt`, `results-a.txt`,
`results-move-r-build-phases.txt`, `results-mds-build-phases.txt` and `statistics.txt`.
They are the inputs to the paper's figures and tables; to import them into LaTeX, use
[sqlplot-tools](https://github.com/bingmann/sqlplot-tools).
