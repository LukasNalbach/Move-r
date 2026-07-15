#!/usr/bin/env bash
#
# Reproduces the Move-r measurements for a single text.
#
# Layout (all relative to this script's directory, measurements/move-r/):
#   texts/     <- place the input text(s) here (e.g. texts/einstein.en.txt)
#   patterns/  <- two pattern sets per text: <text>-patterns-bal and <text>-patterns-phi
#   results/   <- all measurement outputs are written here
#   indexes/   <- working dir for the built indexes
#
# The CLI/benchmark binaries are expected in the project's build/ folder (see README:
# "mkdir build && cd build && cmake .. && make"). move-r-bench is built into build/bench/,
# the other tools into build/cli/.
#
# Usage: ./measure-text.sh -t <text> [-p <max_threads>]
#   -t  text file name inside texts/
#   -p  maximum number of threads (default 1; the paper builds/queries single-threaded)

set -euo pipefail

# ---------- paths ----------
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)
BUILD="$REPO_ROOT/build"

TEXTS="$SCRIPT_DIR/texts"
PATTERNS="$SCRIPT_DIR/patterns"
RESULTS="$SCRIPT_DIR/results"
INDEXES="$SCRIPT_DIR/indexes"

cd "$SCRIPT_DIR"

# ---------- arguments ----------
p=1
t=""
while getopts p:t: flag; do
    case "$flag" in
        p) p=$OPTARG ;;
        t) t=$OPTARG ;;
        *) echo "usage: $0 -t <text> [-p <max_threads>]" >&2; exit 1 ;;
    esac
done
: "${t:?usage: $0 -t <text> [-p <max_threads>]}"

mkdir -p "$RESULTS" "$INDEXES"

if [ ! -f "$TEXTS/$t" ]; then
    echo "error: input text $TEXTS/$t not found" >&2
    exit 1
fi

# ---------- binary check ----------
for bin in "$BUILD/bench/move-r-bench" "$BUILD/cli/move-r-build" "$BUILD/cli/move-r-locate"; do
    if [ ! -x "$bin" ]; then
        echo "error: missing binary '$bin' -- build the project first (see README)" >&2
        exit 1
    fi
done

# move-r-bench builds the competitor indexes (block_RLBWT, grlBWT, ...) via paths
# relative to the project root, so it must be run from there. Result/text/pattern paths
# are absolute, so the outputs still land in this directory's results/.
echo ">>> [$t] comparison benchmark (count / locate vs. competitor indexes)" >&2
( cd "$REPO_ROOT" && "$BUILD/bench/move-r-bench" -m "$RESULTS/results.txt" \
    "$TEXTS/$t" "$PATTERNS/$t-patterns-bal" "$PATTERNS/$t-patterns-phi" )

echo ">>> [$t] balancing-parameter sweep (count / locate for a range of thread counts)" >&2
( cd "$REPO_ROOT" && "$BUILD/bench/move-r-bench" -a -m "$RESULTS/results-a.txt" \
    "$TEXTS/$t" "$PATTERNS/$t-patterns-bal" "$PATTERNS/$t-patterns-phi" "$p" )

echo ">>> [$t] construction phases (a=8, threads 1..$p)" >&2
for ((i_p=1; i_p<=p; i_p*=2)); do
    "$BUILD/cli/move-r-build" -o "$INDEXES/$t-8" -p "$i_p" -a 8 \
        -m_idx "$RESULTS/results-move-r-build-phases.txt" \
        -m_mds "$RESULTS/results-mds-build-phases.txt" "$TEXTS/$t"
done

echo ">>> [$t] locate statistics (a=2 index)" >&2
"$BUILD/cli/move-r-build" -o "$INDEXES/$t-2" -p "$p" -a 2 "$TEXTS/$t"
"$BUILD/cli/move-r-locate" -m "$RESULTS/statistics.txt" "$t" "$INDEXES/$t-2.move-r" "$PATTERNS/$t-patterns-bal"
"$BUILD/cli/move-r-locate" -m "$RESULTS/statistics.txt" "$t" "$INDEXES/$t-2.move-r" "$PATTERNS/$t-patterns-phi"

echo ">>> [$t] done" >&2
