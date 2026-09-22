#!/usr/bin/env bash
#
# Regenerates every figure and table of the paper from a set of measurement results and
# builds them into a single PDF (plots.pdf).
#
# Usage:
#   ./make-plots.sh            use results/          (the output of ./measure-all.sh)
#   ./make-plots.sh --paper    use results-paper/    (the data behind the paper's figures)
#   ./make-plots.sh --no-pdf   only run sqlplot-tools, do not call pdflatex
#
# Nothing under charts/ or tables/ is modified: the run is staged in build-plots/ (or
# build-plots-paper/), so you can diff your numbers against the paper's, e.g.
#
#   diff -u charts/apm.tex build-plots/charts/apm.tex
#
# Requirements:
#   * sqlplot-tools (https://github.com/bingmann/sqlplot-tools), built with the SQLite
#     backend. Put it on your PATH or point $SQLPLOT_TOOLS at the binary.
#   * pdflatex with pgfplots, siunitx, subcaption and caption (TeX Live 2020+).

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

SQLPLOT=${SQLPLOT_TOOLS:-sqlplot-tools}
SOURCE=results
STAGE=build-plots
run_pdf=1

while [ $# -gt 0 ]; do
    case "$1" in
        -p|--paper)  SOURCE=results-paper; STAGE=build-plots-paper ;;
        -n|--no-pdf) run_pdf=0 ;;
        -h|--help)   sed -n '2,20p' "$0"; exit 0 ;;
        *) echo "error: unknown option '$1' (try --help)" >&2; exit 1 ;;
    esac
    shift
done

if ! command -v "$SQLPLOT" >/dev/null 2>&1; then
    echo "error: '$SQLPLOT' not found." >&2
    echo "       Build sqlplot-tools (https://github.com/bingmann/sqlplot-tools) with the" >&2
    echo "       SQLite backend and put it on your PATH, or set SQLPLOT_TOOLS=/path/to/it." >&2
    exit 1
fi

if [ ! -d "$SOURCE" ]; then
    echo "error: $SOURCE/ does not exist -- run ./measure-all.sh first, or pass --paper" >&2
    exit 1
fi

# ---------- stage ----------
rm -rf "$STAGE"
mkdir -p "$STAGE/results"
cp -r charts tables plot-styles.tex plots.tex "$STAGE/"
# the % IMPORT-DATA lines read results/<file>, so the chosen source is staged under that name
for f in results-apm.txt results-ext.txt results-build.txt; do
    [ -f "$SOURCE/$f" ] && cp "$SOURCE/$f" "$STAGE/results/$f"
done

# ---------- which chart/table needs which result file ----------
needs() { head -c 4096 "$1" | grep -o 'results/results-[a-z]*\.txt' | sort -u; }

cd "$STAGE"
stale=""
for f in charts/*.tex tables/*.tex; do
    skip=""
    for dep in $(needs "$f"); do
        [ -f "$dep" ] || skip="$skip $dep"
    done
    if [ -n "$skip" ]; then
        echo "skipping $f -- missing:$skip" >&2
        stale="$stale $f"
        continue
    fi
    if out=$("$SQLPLOT" "$f" 2>&1); then
        echo "regenerated $f"
    else
        # one chart failing must not cost you the others -- keep its committed version and go on
        echo "warning: sqlplot-tools failed on $f, keeping the paper's version" >&2
        sed 's/^/    /' <<< "$out" >&2
        stale="$stale $f"
    fi
done

if [ -n "$stale" ]; then
    echo >&2
    echo "note: these still hold the paper's numbers, not yours:" >&2
    for f in $stale; do echo "        $f" >&2; done
    echo "      Every chart and table filters on the paper's text names (WHERE text =" >&2
    echo "      'sars2.ACGT.50Gi' and so on), so with your own texts the queries return no rows" >&2
    echo "      and sqlplot-tools reports e.g. \"MULTIPLOT() requires group column list\"." >&2
    echo "      Edit the text names in the affected file, or measure the paper's texts. See §6." >&2
fi

# ---------- build the PDF ----------
if [ "$run_pdf" = 1 ]; then
    if ! command -v pdflatex >/dev/null 2>&1; then
        echo "warning: pdflatex not found -- the regenerated .tex files are in $STAGE/" >&2
        exit 0
    fi
    # twice, so the float numbering settles
    pdflatex -interaction=nonstopmode plots.tex > plots.build.log 2>&1 || true
    pdflatex -interaction=nonstopmode plots.tex > plots.build.log 2>&1 || true
    if [ -f plots.pdf ]; then
        echo
        echo "wrote $STAGE/plots.pdf"
    else
        echo "error: pdflatex failed -- see $STAGE/plots.build.log" >&2
        exit 1
    fi
fi
