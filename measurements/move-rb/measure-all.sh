#!/usr/bin/env bash
#
# Runs the full measurement suite for every text placed in texts/.
# Executed inside the bench directory (benchmarks/move-rb/); invoke it from anywhere.
#
# Per text we pass:
#   -T  gen-apm-queries --time (calibration target, s): 5 for the DNA texts, 1 for dewiki
#   -M  gen-apm-queries --min  (min patterns per set):  4 for the DNA texts, 2 for dewiki
#   -C  build/measure columba & columba-rlc:            1 for DNA, 0 for dewiki (DNA-only)
#
# Results land in results/:
#   results-build-move-rb.txt, results-build-move-rb-rlzsa.txt,
#   results-build-br-index.txt, results-build-columba.txt, results-build-columba-rlc.txt,
#   results-apm.txt, results-ext.txt

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

mkdir -p results
# start from clean query-result files (per-text build metrics are appended to their own logs)
rm -f results/results-apm.txt results/results-ext.txt results/results-build-*.txt

./measure-text.sh -t sars2.ACGT.50Gi -T 5 -M 4 -C 1
./measure-text.sh -t chr19.ACGT.50Gi -T 5 -M 4 -C 1
./measure-text.sh -t dewiki.50Gi     -T 1 -M 2 -C 0
