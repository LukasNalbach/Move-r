#!/usr/bin/env bash
#
# Runs the Move-r measurement suite for every text placed in texts/.
# Invoke it from anywhere; it operates inside its own directory (measurements/move-r/).
#
# Usage: ./measure-all-texts.sh [-p <max_threads>]
#   -p  maximum number of threads (default 1)

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

p=1
while getopts p: flag; do
    case "$flag" in
        p) p=$OPTARG ;;
        *) echo "usage: $0 [-p <max_threads>]" >&2; exit 1 ;;
    esac
done

mkdir -p results
# start from clean result files (build-phase logs are appended, so clear them too)
rm -f results/results.txt results/results-a.txt results/statistics.txt \
      results/results-move-r-build-phases.txt results/results-mds-build-phases.txt

for t in einstein.en.txt english dewiki sars2 chr19; do
    ./measure-text.sh -p "$p" -t "$t"
done
