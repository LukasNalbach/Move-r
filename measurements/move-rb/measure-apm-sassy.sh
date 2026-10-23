#!/usr/bin/env bash
#
# Measures approximate pattern matching with move_rb and with sassy
# (https://github.com/RagnarGrootKoerkamp/sassy) on a FASTA file.
#
# Steps:
#   1. builds a move_rb index of the FASTA file in FASTA mode
#   2. writes the text of the index to a file (fasta-to-text); sassy searches this file and the
#      patterns are sampled from it
#   3. generates the pattern sets
#   4. runs move-rb-bench-apm
#
# sassy supports only locate queries under the edit distance; the count and hamming-distance
# results are those of move_rb only.
#
# Usage: ./measure-apm-sassy.sh -f <fasta> [options]
#   -f <fasta>    FASTA file (e.g. ~/texts/chr19)
#   -d <dir>      directory for the index and the text file (default: directory of the FASTA file)
#   -o <dir>      directory for the results (default: results/ in this script's directory)
#   -k <list>     numbers of errors             (default: 1,2,4,7,10,13)
#   -m <list>     pattern lengths               (default: 10,20,40,80,160,320,640,1280)
#   -S <count>    patterns per (k, m) for sassy (default: 2)
#   -P <list>     sassy profiles                (default: ascii; dna maps every byte to {A,C,G,T},
#                 i.e., also the separator)
#   -T <seconds>  gen-apm-queries --time: calibration target per pattern set (default: 1)
#   -R <seconds>  time each pattern set is replayed for move_rb (default: 10; sassy runs each
#                 pattern once)
#   -p <threads>  threads for the index construction (default: all)

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)
BUILD="$REPO_ROOT/build"

fasta=""
work=""
results="$SCRIPT_DIR/results"
k_list="1,2,4,7,10,13"
m_list="10,20,40,80,160,320,640,1280"
sassy_patterns=2
sassy_profiles="ascii"
calibration=1
replay=10
threads=""

while getopts f:d:o:k:m:S:P:T:R:p: flag; do
    case "$flag" in
        f) fasta=$OPTARG ;;
        d) work=$OPTARG ;;
        o) results=$OPTARG ;;
        k) k_list=$OPTARG ;;
        m) m_list=$OPTARG ;;
        S) sassy_patterns=$OPTARG ;;
        P) sassy_profiles=$OPTARG ;;
        T) calibration=$OPTARG ;;
        R) replay=$OPTARG ;;
        p) threads=$OPTARG ;;
        *) echo "usage: $0 -f <fasta> [-d dir] [-o dir] [-k list] [-m list] [-S count] [-P list] [-T s] [-R s] [-p threads]" >&2; exit 1 ;;
    esac
done
: "${fasta:?usage: $0 -f <fasta> [-d dir] [-o dir] [-k list] [-m list] [-S count] [-P list] [-T s] [-R s] [-p threads]}"

[ -f "$fasta" ] || { echo "error: $fasta does not exist" >&2; exit 1; }
name=$(basename "$fasta")
[ -n "$work" ] || work=$(cd "$(dirname "$fasta")" && pwd)
mkdir -p "$results" "$work"

index="$work/$name.move-rb-locate"
text="$work/$name.txt"
patterns="$work/$name.patterns-apm"
build_threads=()
[ -n "$threads" ] && build_threads=(-p "$threads")

# ---------- 1. index (FASTA mode) ----------
if [ ! -s "$index" ]; then
    echo ">>> [$name] building the move_rb index (FASTA mode)" >&2
    "$BUILD/cli/move-rb" build -f -c pfp -s locate-move "${build_threads[@]}" \
        -m "$results/results-build-$name.txt" "$name" "$fasta"
else
    echo ">>> [$name] the index already exists ($index)" >&2
fi

# ---------- 2. text of the index (searched by sassy, the patterns are sampled from it) ----------
if [ ! -s "$text" ]; then
    echo ">>> [$name] writing the text file for sassy" >&2
    "$BUILD/cli/fasta-to-text" "$fasta" "$text"
else
    echo ">>> [$name] the text file already exists ($text)" >&2
fi

# ---------- 3. pattern sets (-x: without the FASTA separator 0x7F) ----------
if [ ! -d "$patterns" ]; then
    echo ">>> [$name] generating the pattern sets (k = $k_list, m = $m_list)" >&2
    mkdir -p "$patterns"
    "$BUILD/bench/move-rb-gen-apm-queries" \
        -k "$k_list" -m "$m_list" -x $'\x7f' --time "$calibration" --min 2 --timeout 10 --seed 42 \
        -n "$name" -o "$patterns/" "$text" "$index"
else
    echo ">>> [$name] the pattern sets already exist ($patterns)" >&2
fi

# ---------- 4. measurement ----------
echo ">>> [$name] measuring move_rb and sassy" >&2
"$BUILD/bench/move-rb-bench-apm" \
    --time "$replay" --move-rb "$index" \
    --sassy "$text" --sassy-profile "$sassy_profiles" --sassy-patterns "$sassy_patterns" \
    -m "$results/results-apm-sassy.txt" \
    "$name" "$patterns/"

echo ">>> [$name] done; results in $results" >&2
