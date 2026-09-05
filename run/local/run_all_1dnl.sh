#!/bin/bash
# ==============================================================
#  run_all_1dnl.sh — the 1-D NONLINEAR INSTABILITY study (4 tiers x 4 amplitudes).
#
#  Separate from run_all_1d.sh on purpose: that one globs `run_1d_*.sh` and runs
#  the six A=0.001 tier-agreement cases; this one globs `run_1dnl_*.sh`, which
#  does not match it. The two sets answer different questions and should not be
#  merged into one batch.
#
#  Each case is SEQUENTIAL (direct LU beats any MPI split of a 240x1 mesh) and the
#  cases run SIDE BY SIDE, so the cores go to case-level parallelism. With ny=1
#  these are ~3x smaller than the old ny=3 flume and far cheaper in wall time.
# ==============================================================
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
JOBS="${JOBS:-12}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1
LOG=output/local/logs; mkdir -p $LOG; : > $LOG/batch1dnl_summary.txt
for f in run/local/run_1dnl_*.sh; do
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 5; done
  n=$(basename "$f" .sh)
  ( bash "$f" > "$LOG/${n}.log" 2>&1
    echo "$n exit=$?" >> $LOG/batch1dnl_summary.txt ) &
done
wait
echo "ALL 1DNL DONE" >> $LOG/batch1dnl_summary.txt
