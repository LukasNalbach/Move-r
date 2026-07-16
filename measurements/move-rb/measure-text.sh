#!/usr/bin/env bash
#
# Reproduces all measurements for a single text.
#
# Layout (all relative to this script's directory):
#   texts/    <- reproducer drops the input text(s) here (e.g. texts/sars2.ACGT.50Gi)
#   results/  <- all measurement outputs are written here
#   indexes/, indexes_columba/, indexes_columba_rlc/, patterns/, patterns_ext/  <- working dirs
#
# The CLI/benchmark binaries are expected in the project's build/ folder (see README:
# "mkdir build && cd build && cmake .. && cp -rf ../patched-files/* .. && make").
#
# Usage: ./measure-text.sh -t <text> [-T <apm_time>] [-M <apm_min>] [-C <0|1>] [-p <threads>]
#   -t  text file name inside texts/
#   -T  gen-apm-queries --time  : per-set calibration target in seconds (default 5)
#   -M  gen-apm-queries --min   : minimum number of patterns per set    (default 4)
#   -C  1 = also build/measure columba & columba-rlc (DNA only), 0 = skip them (default 1)
#   -p  build threads (default 1; the paper builds and queries single-threaded)

set -euo pipefail

# ---------- paths ----------
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)
BUILD="$REPO_ROOT/build"
BRI_BUILD="$BUILD/cli/bri-build"   # Move-r builds its own bri-build into build/cli/ (see CMakeLists.txt)

TEXTS="$SCRIPT_DIR/texts"
RESULTS="$SCRIPT_DIR/results"
INDEXES="$SCRIPT_DIR/indexes"
INDEXES_COLUMBA="$SCRIPT_DIR/indexes_columba"
INDEXES_COLUMBA_RLC="$SCRIPT_DIR/indexes_columba_rlc"
PATTERNS="$SCRIPT_DIR/patterns"
PATTERNS_EXT="$SCRIPT_DIR/patterns_ext"

# always operate inside the bench directory, regardless of where the script is invoked from
cd "$SCRIPT_DIR"

# ---------- arguments ----------
apm_time=5
apm_min=4
run_columba=1
threads=1
t=""
while getopts t:T:M:C:p: flag; do
    case "$flag" in
        t) t=$OPTARG ;;
        T) apm_time=$OPTARG ;;
        M) apm_min=$OPTARG ;;
        C) run_columba=$OPTARG ;;
        p) threads=$OPTARG ;;
        *) echo "usage: $0 -t <text> [-T apm_time] [-M apm_min] [-C 0|1] [-p threads]" >&2; exit 1 ;;
    esac
done
: "${t:?usage: $0 -t <text> [-T apm_time] [-M apm_min] [-C 0|1] [-p threads]}"

mkdir -p "$RESULTS" "$INDEXES" "$PATTERNS" "$PATTERNS_EXT"
[ "$run_columba" = 1 ] && mkdir -p "$INDEXES_COLUMBA" "$INDEXES_COLUMBA_RLC"

if [ ! -f "$TEXTS/$t" ]; then
    echo "error: input text $TEXTS/$t not found" >&2
    exit 1
fi

# ---------- binary check ----------
required=( "$BUILD/cli/move-rb-build" \
           "$BUILD/bench/move-rb-gen-apm-queries" "$BUILD/bench/move-rb-bench-apm" \
           "$BUILD/bench/move-rb-gen-ext-queries" "$BUILD/bench/move-rb-bench-ext" "$BRI_BUILD" )
[ "$run_columba" = 1 ] && required+=( "$BUILD/cli/columba-build" "$BUILD/cli/columba-rlc-build" )
for bin in "${required[@]}"; do
    if [ ! -x "$bin" ]; then
        echo "error: missing binary '$bin' -- build the project first (see README)" >&2
        exit 1
    fi
done

# Build a competitor index under /usr/bin/time -v (wall time + peak RSS -> log).
# Never aborts the run: a failing competitor build is reported and skipped, and the
# benchmarks below simply omit that index (each index is loaded independently).
run_timed() {   # run_timed <label> <logfile> <cmd> [args...]
    local label=$1 log=$2; shift 2
    printf '# %s\n' "$*" >> "$log"
    if ! /usr/bin/time -v "$@" >> "$log" 2>&1; then
        echo "warning: building $label failed; continuing without it (see $log)" >&2
    fi
}

echo ">>> [$t] building indexes" >&2

# ============ 1. build all indexes ============
# move-rb and move-rb-rlzsa write their own construction metrics via -m_idx.
"$BUILD/cli/move-rb-build" -s locate_move  -c bigbwt -p "$threads" \
    -m_idx "$RESULTS/results-build-move-rb.txt"       -o "$INDEXES/$t" "$TEXTS/$t"
"$BUILD/cli/move-rb-build" -s locate_rlzsa -c bigbwt -p "$threads" \
    -m_idx "$RESULTS/results-build-move-rb-rlzsa.txt" -o "$INDEXES/$t" "$TEXTS/$t"

# br-index prints only its build time; /usr/bin/time -v adds the peak memory.
echo ">>> [$t] building br-index" >&2
run_timed "br-index" "$RESULTS/results-build-br-index.txt" \
    "$BRI_BUILD" -divsufsort -o "$INDEXES/$t" "$TEXTS/$t"

# columba / columba-rlc (b-move) only support DNA, so skip them for non-DNA texts (-C 0).
# NOTE: columba's -f expects a (multi-)FASTA reference; if the plain text is rejected,
#       wrap it as a single-record FASTA first.
if [ "$run_columba" = 1 ]; then
    echo ">>> [$t] building columba" >&2
    run_timed "columba" "$RESULTS/results-build-columba.txt" \
        "$BUILD/cli/columba-build"     -t "$threads"    -r "$INDEXES_COLUMBA/$t"     -f "$TEXTS/$t"
    echo ">>> [$t] building columba-rlc" >&2
    run_timed "columba-rlc" "$RESULTS/results-build-columba-rlc.txt" \
        "$BUILD/cli/columba-rlc-build" -t "$threads" -p -r "$INDEXES_COLUMBA_RLC/$t" -f "$TEXTS/$t"
fi

# assemble the competitor index flags shared by both benchmarks
index_flags=( --bri           "$INDEXES/$t.bri"
              --move-rb        "$INDEXES/$t.move-rb"
              --move-rb-rlzsa  "$INDEXES/$t.move-rb-rlzsa" )
if [ "$run_columba" = 1 ]; then
    index_flags+=( --columba-rlc "$INDEXES_COLUMBA_RLC/$t"
                   --columba     "$INDEXES_COLUMBA/$t" )
fi

# ============ 2. approximate pattern matching (count + locate) ============
echo ">>> [$t] APM: generating patterns (--time $apm_time, --min $apm_min) and measuring" >&2
"$BUILD/bench/move-rb-gen-apm-queries" \
    --min "$apm_min" --time "$apm_time" --timeout 10 --seed 42 \
    -n "$t" -o "$PATTERNS/" \
    "$TEXTS/$t" "$INDEXES/$t.move-rb-rlzsa"

"$BUILD/bench/move-rb-bench-apm" \
    --time 10 --cigar both --algo both \
    "${index_flags[@]}" \
    "$t" "$PATTERNS/" >> "$RESULTS/results-apm.txt"

# ============ 3. raw extension + SA-interval enumeration ============
echo ">>> [$t] extension: generating patterns and measuring" >&2
"$BUILD/bench/move-rb-gen-ext-queries" \
    --min "$apm_min" --time 1 \
    -n "$t" -o "$PATTERNS_EXT/" \
    "$TEXTS/$t" "$INDEXES/$t.move-rb-rlzsa"

"$BUILD/bench/move-rb-bench-ext" \
    --time 10 \
    "${index_flags[@]}" \
    "$t" "$PATTERNS_EXT/" >> "$RESULTS/results-ext.txt"

echo ">>> [$t] done" >&2
