#!/usr/bin/env bash
#
# Assembles results/results-build.txt -- the input of charts/construction.tex -- from the
# per-index build logs that measure-text.sh writes.
#
#   ./make-results-build.sh
#
# One RESULT line per (index, text):
#
#   RESULT algo=build_move_rb_move text=... n=... time_build=... peak_memory_usage=... size_index=...
#
# with time_build in nanoseconds and peak_memory_usage / size_index in bytes.
#
# move-rb and move-rb-rlzsa report these directly in their -m_idx record. For br-index,
# columba and columba-rlc the fields are read out of their build logs: the wall time from
# /usr/bin/time -v, the peak heap from the `Peak memory usage during construction:` line the
# build tools print via malloc_count, n from the input text and size_index from the index
# files on disk. Builds that failed or are missing are skipped with a warning.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

RESULTS=results
OUT="$RESULTS/results-build.txt"

if [ ! -d "$RESULTS" ]; then
    echo "error: $RESULTS/ does not exist -- run ./measure-all.sh first" >&2
    exit 1
fi

# ---------- helpers ----------

# size of a file, or 0 if it does not exist
fsize() { [ -f "$1" ] && stat -c %s "$1" || echo 0; }

# summed size of every file matching a prefix (columba writes several files per index)
prefix_size() {
    local total=0 f
    for f in "$1"*; do [ -f "$f" ] && total=$((total + $(stat -c %s "$f"))); done
    echo "$total"
}

# "Elapsed (wall clock) time (h:mm:ss or m:ss): 1:02.34" -> nanoseconds
elapsed_ns() {
    awk -F': ' '
        /Elapsed \(wall clock\) time/ {
            n = split($NF, p, ":")
            s = (n == 3) ? p[1]*3600 + p[2]*60 + p[3] : p[1]*60 + p[2]
            printf "%.0f", s * 1e9
            exit
        }'
}

# "Peak memory usage during construction: 12345 bytes"  (bri-build, plain byte count)
# "Peak memory usage during construction: 1.23 GiB"     (columba, human readable)
peak_bytes() {
    awk '
        /Peak memory usage during construction:/ {
            line = $0
            sub(/.*Peak memory usage during construction:[ \t]*/, "", line)
            split(line, p, /[ \t]+/)
            v = p[1] + 0
            u = p[2]
                 if (u ~ /^KiB/) v *= 1024
            else if (u ~ /^MiB/) v *= 1024 * 1024
            else if (u ~ /^GiB/) v *= 1024 * 1024 * 1024
            else if (u ~ /^TiB/) v *= 1024 * 1024 * 1024 * 1024
            if (v > max) max = v
        }
        END { printf "%.0f", max + 0 }'
}

# pull "key=value" out of a RESULT line
field() { sed -n "s/.*[[:space:]]$2=\([^[:space:]]*\).*/\1/p" <<< "$1"; }

# strip the quotes a logged command line may carry around a path
unquote() { local s=$1; s=${s#[\"\']}; s=${s%[\"\']}; printf '%s' "$s"; }

# emit <algo> <text> <n> <time_build> <peak> <size_index> [sigma] [r] [r_rev]
# The last three are only known for move-rb (the competitors do not report them); rows may
# therefore carry different fields, which sqlplot-tools handles -- the union of all keys becomes
# the table and a row that lacks one gets NULL there.
emit() {
    printf 'RESULT algo=%s text=%s n=%s time_build=%s peak_memory_usage=%s size_index=%s' \
        "$1" "$2" "$3" "$4" "$5" "$6" >> "$OUT"
    [ -n "${7:-}" ] && printf ' sigma=%s' "$7" >> "$OUT"
    [ -n "${8:-}" ] && printf ' r=%s' "$8" >> "$OUT"
    [ -n "${9:-}" ] && printf ' r_rev=%s' "$9" >> "$OUT"
    printf '\n' >> "$OUT"
}

: > "$OUT"
rows=0

# ---------- move-rb / move-rb-rlzsa: the -m_idx record already has every field ----------
for log in "$RESULTS/results-build-move-rb.txt" "$RESULTS/results-build-move-rb-rlzsa.txt"; do
    [ -f "$log" ] || { echo "warning: $log missing, skipping" >&2; continue; }
    while IFS= read -r line; do
        [ -n "$line" ] || continue
        algo=$(field "$line" algo)
        text=$(field "$line" text)
        n=$(field "$line" n)
        tb=$(field "$line" time_build)
        pk=$(field "$line" peak_memory_usage)
        sz=$(field "$line" size_index)
        # alphabet size and the two BWT run counts; tables/texts.tex derives sigma, n/r and
        # n/r_rev from them. Older -m_idx records do not have them, hence the empty fallback.
        sg=$(field "$line" sigma)
        r=$(field "$line" r)
        rrev=$(field "$line" r_rev)
        if [ -z "$n" ] || [ -z "$tb" ] || [ -z "$pk" ] || [ -z "$sz" ]; then
            echo "warning: $log has no overall metrics for '$algo $text' -- rebuild with a" >&2
            echo "         move-rb-build that writes n/time_build/peak_memory_usage/size_index" >&2
            continue
        fi
        emit "$algo" "$text" "$n" "$tb" "$pk" "$sz" "$sg" "$r" "$rrev"
        rows=$((rows + 1))
    done < <(grep '^RESULT' "$log" || true)
done

# ---------- br-index / columba / columba-rlc: read the fields out of the build log ----------
# each build appends "# <command line>" followed by the tool's output and GNU time's report
competitor() {  # competitor <algo> <logfile>
    local algo=$1 log=$2
    [ -f "$log" ] || { echo "warning: $log missing, skipping $algo" >&2; return; }

    local block cmd text_path index_path n tb pk sz
    while IFS= read -r block; do
        cmd=${block%%$'\x01'*}
        body=${block#*$'\x01'}
        body=${body//$'\x01'/$'\n'}

        # the input text is the last argument (bri-build) or the argument of -f (columba)
        if [[ $cmd == *" -f "* ]]; then
            text_path=$(unquote "$(sed -n 's/.* -f \([^ ]*\).*/\1/p' <<< "$cmd")")
            # columba is fed the FASTA copy measure-text.sh writes next to the text; report the
            # row under the text's own name and length so it joins the other indexes' rows
            for ext in .fasta .fa .FASTA .FA .fna .FNA; do
                if [ "${text_path%"$ext"}" != "$text_path" ] && [ -f "${text_path%"$ext"}" ]; then
                    text_path=${text_path%"$ext"}
                    break
                fi
            done
        else
            text_path=${cmd##* }
        fi
        # the index is the argument of -o (bri-build) or -r (columba)
        if [[ $cmd == *" -r "* ]]; then
            index_path=$(unquote "$(sed -n 's/.* -r \([^ ]*\).*/\1/p' <<< "$cmd")")
            sz=$(prefix_size "$index_path")
        else
            index_path=$(unquote "$(sed -n 's/.* -o \([^ ]*\).*/\1/p' <<< "$cmd")")
            sz=$(fsize "$index_path.bri")
        fi

        n=$(fsize "$text_path")
        tb=$(elapsed_ns <<< "$body")
        pk=$(peak_bytes <<< "$body")

        if [ "$n" = 0 ] || [ -z "$tb" ] || [ "$tb" = 0 ]; then
            echo "warning: no usable build record for $algo in $log (failed build?), skipping" >&2
            continue
        fi
        [ "$pk" = 0 ] && echo "warning: $algo on $(basename "$text_path"): no malloc_count peak in the log" >&2
        [ "$sz" = 0 ] && echo "warning: $algo on $(basename "$text_path"): index files not found, size_index=0" >&2

        emit "$algo" "$(basename "$text_path")" "$n" "$tb" "$pk" "$sz"
        rows=$((rows + 1))
    done < <(awk '/^# /{ if (n++) printf "\n"; printf "%s", substr($0,3); next } n { printf "\x01%s", $0 } END { if (n) printf "\n" }' "$log")
}

competitor build_br_index "$RESULTS/results-build-br-index.txt"
competitor build_columba  "$RESULTS/results-build-columba.txt"
competitor build_bmove    "$RESULTS/results-build-columba-rlc.txt"   # b-move in the paper

echo "wrote $OUT ($rows rows)"
