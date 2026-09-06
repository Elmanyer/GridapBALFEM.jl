#!/bin/bash
# ==============================================================
#  run_all_th.sh — the Taylor-Hood nonlinear-instability campaign.
#  1-D ONLY, TAYLOR-HOOD ONLY (CLAUDE.md rule 2b). See run_th_*.sh headers.
#  Selection: pass a glob fragment, e.g.  bash run_all_th.sh sdirk
# ==============================================================
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
SEL="${1:-}"
JOBS="${JOBS:-6}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1
#  ⚠ PER-SELECTION summary, APPENDED. A shared file that is truncated at start lets a
#  second invocation (e.g. `run_all_th.sh rk4`) wipe the first's record and signal
#  "done" while the first batch is still running — which is exactly what happened
#  on 2026-09-05. Never share a completion file between concurrent batches.
LOG=output/local/logs; mkdir -p $LOG
TAG="${SEL:-all}"; SUMM="$LOG/batch_th_${TAG}_summary.txt"; : > "$SUMM"
for f in run/local/run_th_*${SEL}*.sh; do
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 5; done
  n=$(basename "$f" .sh)
  ( bash "$f" > "$LOG/${n}.log" 2>&1
    echo "$n exit=$?" >> "$SUMM" ) &
done
wait
echo "ALL TH DONE ($TAG)" >> "$SUMM"
