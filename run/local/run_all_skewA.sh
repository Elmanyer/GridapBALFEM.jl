#!/bin/bash
# ==============================================================
#  run_all_skewA.sh — Tier A of the skew-advection evaluation (6 cases).
#
#  Deliberately NOT globbed by run_all_1dnl.sh (which matches run_1dnl_*.sh):
#  these answer a different question and must not silently join that batch.
#
#  ⚠ dt IS NOT SCALED WITH dx. That is intentional: NONLINEAR_INSTABILITY.md F4
#  establishes the growth is dt-insensitive at fixed dx, so holding dt fixed keeps
#  the refinement triple a ONE-VARIABLE comparison. Tier C1 varies dt separately.
#
#  Cost (measured from the baselines): ~2 h per 80 s run at dx=0.25; the dx=0.125
#  legs are ~6-8 h. Six cases side by side fit this box (~1.9 GB RSS each), but
#  ONE PROCESS PER CASE, exiting at the end — rule 41: a long-lived Julia+Gridap
#  process degrades to 7-14 % CPU, and a shared worker would poison the results.
# ==============================================================
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
JOBS="${JOBS:-6}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1
LOG=output/local/logs; mkdir -p $LOG; : > $LOG/batch_skewA_summary.txt
for f in run/local/run_skewA_*.sh; do
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 5; done
  n=$(basename "$f" .sh)
  ( bash "$f" > "$LOG/${n}.log" 2>&1
    echo "$n exit=$?" >> $LOG/batch_skewA_summary.txt ) &
done
wait
echo "ALL SKEW TIER-A DONE" >> $LOG/batch_skewA_summary.txt
