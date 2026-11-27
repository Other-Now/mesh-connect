#!/usr/bin/env bash
# Reproduce every number in README.md.
#   scripts/run_all.sh <mesh-connect binary> <work dir> [ccx]
# Needs ~300 MB RAM for the 1.2M-element runs and ~2.5 GB for the finest
# CalculiX lap model (LAP_SCALES=1,2,3,4).
set -euo pipefail
BIN=${1:?mesh-connect binary}
WORK=${2:?work dir}
CCX=${3:-ccx}
OUT=results
mkdir -p "$WORK" "$OUT"

"$BIN" proj --n 200000 > "$OUT/projection_accuracy.json"

"$BIN" gen "$WORK/big" --parts 8 --h 2 --welds 100000 --seed 1
"$BIN" check   "$WORK/big/mesh.bdf"                          > "$OUT/check_1m.json"
"$BIN" eval    "$WORK/big/mesh.bdf" "$WORK/big/labels.csv"   > "$OUT/eval_1m.json"
"$BIN" connect "$WORK/big/mesh.bdf" "$WORK/big/welds.csv" --out "$WORK/big/cwelds.bdf" > "$OUT/connect_1m.json"
for i in 1 2; do
  "$BIN" bench "$WORK/big/mesh.bdf" "$WORK/big/welds.csv" --threads 1,2,3,4 --reps 7 --brute 200 --ab-reps 11 \
    > "$OUT/bench_1m_run$i.json"
done

LAP_SCALES=${LAP_SCALES:-1,2,3,4} python3 validate/run_lap.py "$BIN" "$CCX" "$WORK/lap" "$OUT/lap_calculix.json"
echo "results written to $OUT/"
