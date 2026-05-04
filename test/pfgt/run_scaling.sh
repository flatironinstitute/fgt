#!/usr/bin/env bash
# test/pfgt/run_scaling.sh
#
# Runs ./int2-pfgt-perf at a range of OMP_NUM_THREADS values and prints
# a Markdown-friendly scaling table to stdout. Saves a copy to
# ../../docs/perf/<host>-<utc-date>-pfgt-scaling.csv

set -euo pipefail

cd "$(dirname "$0")"

if [ ! -x ./int2-pfgt-perf ]; then
    echo "ERROR: build the perf binary first: cd ../../ && make test/pfgt-perf"
    exit 1
fi

# physical core count
case "$(uname -s)" in
    Darwin) NCORE=$(sysctl -n hw.physicalcpu) ;;
    Linux)  NCORE=$(getconf _NPROCESSORS_ONLN) ;;
    *)      NCORE=$(nproc 2>/dev/null || echo 4) ;;
esac

# build a thread-count list: 1, 2, 4, 8, ..., NCORE
THREAD_LIST=()
N=1
while [ "$N" -le "$NCORE" ]; do
    THREAD_LIST+=("$N")
    N=$((N * 2))
done
# also include physical core count if it isn't a power of two
LAST_IDX=$((${#THREAD_LIST[@]} - 1))
if [ "${THREAD_LIST[$LAST_IDX]}" -ne "$NCORE" ]; then
    THREAD_LIST+=("$NCORE")
fi

DIM=${FGT_DIM:-3}
NSRC=${FGT_NSRC:-1000000}
EPS=${FGT_EPS:-1d-6}
DELTA=${FGT_DELTA:-1d-4}

OUTDIR=../../docs/perf
mkdir -p "$OUTDIR"
HOST=$(hostname -s)
UTC=$(date -u +%Y%m%d-%H%M%S)
OUTFILE="$OUTDIR/${HOST}-${UTC}-pfgt-scaling.csv"

echo "nthreads,nsrc,dim,eps,total_s" > "$OUTFILE"
echo "# pfgt scaling on host=$HOST  dim=$DIM  nsrc=$NSRC  eps=$EPS  delta=$DELTA"
echo "| nthreads | total_s | speedup | efficiency |"
echo "|---------:|--------:|--------:|-----------:|"

T1=""
for T in "${THREAD_LIST[@]}"; do
    LINE=$(FGT_DIM=$DIM FGT_NSRC=$NSRC FGT_EPS=$EPS FGT_DELTA=$DELTA \
            OMP_NUM_THREADS=$T ./int2-pfgt-perf 2>/dev/null | tail -1)
    # strip any spaces inside the ES9.2 field
    LINE_CLEAN=$(echo "$LINE" | tr -d ' ')
    echo "$LINE_CLEAN" >> "$OUTFILE"
    TS=$(echo "$LINE_CLEAN" | awk -F, '{print $5}')
    if [ -z "$T1" ]; then T1=$TS; fi
    SPEEDUP=$(awk -v t1="$T1" -v ts="$TS" 'BEGIN { printf "%.2f", t1/ts }')
    EFF=$(awk -v t1="$T1" -v ts="$TS" -v t="$T" 'BEGIN { printf "%.0f%%", 100*t1/ts/t }')
    printf "| %8d | %7.3f | %7s | %10s |\n" "$T" "$TS" "$SPEEDUP" "$EFF"
done

echo
echo "Raw CSV: $OUTFILE"
