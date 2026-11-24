#!/bin/bash
# ==============================================================
#  run_1d_m3_M1_tol.sh — re-run of m3_M1_T160_n32 with a tight Newton tolerance (2026-11-10)
#
#  WHY. m3_M1_T160_n32 (run_1d_m3_M1.sh, nl_tol = 1e-5) converged in ONE Newton iteration per step
#  (r ≈ 8e-6, just under the tolerance; every other box run takes 2, to ≈ 1e-7). It lost box
#  momentum conservation: ∫Hū/L 1.9e-3 → 5.8e-2 in 98 s, linear in t, with a steady x-uniform shear
#  current (u_bed +0.29, u_surf −0.26 m/s) taking the energy (monitor E 1.00 → 2.90). Stopped at
#  98 s at the author's request. The M1 T = 3.0 s run (2 iterations) conserves momentum to 4 digits.
#
#  PREDICTION: with nl_tol = 1e-7 (≥ 2 iterations) the momentum stays constant and the n = 0 shear
#  current does not grow. If the shear current still grows with momentum conserved, it is a model
#  effect (a wave-driven mean shear), not the solver.
#
#  DESIGN: run_1d_m3_M1.sh verbatim for the T = 1.6 s case, plus BALFEM_NL_TOL = 1e-7.
#  Detach: setsid nohup run/local/run_1d_m3_M1_tol.sh > output/1D_M1_test/_logs/_batch_tol.out 2>&1 < /dev/null &
# ==============================================================
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
PROJ="$(pwd)"
PERIODS="${PERIODS:-200}"
MEM_PER_CASE="${MEM_PER_CASE:-5G}"
OUT=output/1D_M1_test
LOG=$OUT/_logs
mkdir -p "$LOG"
source "$PROJ/run/local/_slot_throttle.sh"

# case : name | M | T [s] | A [m] (k₀A = 0.1582) | NCELL        longest first
CASES=(
  "m3_M1_T160_n32_tol7  1  1.6  0.0799  32"
)

# ---- source snapshot ---------------------------------------------------------------------------
SNAP=$LOG/source_$(date '+%Y%m%d_%H%M%S')
mkdir -p "$SNAP"
git rev-parse HEAD > "$SNAP/commit.txt"
git diff HEAD > "$SNAP/tracked.diff"
git ls-files --others --exclude-standard -z | tar czf "$SNAP/untracked.tgz" --null -T -
echo "$(date '+%F %T')  M1 tol re-run source snapshot: $SNAP (commit $(cut -c1-10 "$SNAP/commit.txt") + diff)" | tee -a "$LOG/_batch.log"

run_case() {
  read -r name mvert twave awave ncell <<< "$1"
  [ -f "$LOG/$name.log" ] && { echo "skip $name (log exists)"; return; }
  (
    unset BALFEM_CIP_GU BALFEM_CIP_GE
    export BALFEM_REGIME=nonlinear BALFEM_NL_PRESSURE=0 BALFEM_FLAT_BED=1
    export BALFEM_BED_SLOPE=0 BALFEM_BED_NW=1
    export BALFEM_M=$mvert BALFEM_P_VERT=1
    export BALFEM_FE_ORDER=3 BALFEM_P_ETA=2 BALFEM_NCELL=$ncell BALFEM_NLAMBDA=1
    export BALFEM_AWAVE=$awave BALFEM_TWAVE=$twave BALFEM_D=3.5 BALFEM_DT=0.04 BALFEM_PERIODS=$PERIODS
    export BALFEM_SOLVER=theta BALFEM_USE_AD=0 BALFEM_NOISE=1e-8 BALFEM_NL_TOL=1e-7
    export BALFEM_CHECKPOINT_EVERY=250 BALFEM_CHECKPOINT_KEEP=2
    export BALFEM_SAVE_EVERY=5 BALFEM_DIAG_EVERY=5 BALFEM_PRINT_EVERY=25
    export BALFEM_OUTDIR=$PROJ/$OUT/$name
    export OPENBLAS_NUM_THREADS=2
    echo "case $name | m3: NL_PRESSURE=0 FLAT_BED=1, no CIP | M=$mvert T=$twave A=$awave NCELL=$ncell | Q3/Q2 | nl_tol 1e-7 | $PERIODS periods"
    echo "started $(date '+%F %T') | git $(git rev-parse --short HEAD)$(git diff --quiet -- src examples || echo ' + uncommitted edits')"
    systemd-run --user --scope --quiet --unit="m3M1tol-$name-$$" -p MemoryMax="$MEM_PER_CASE" \
      stdbuf -oL -eL julia --project="$PROJ" examples/local_1d/run_periodic_1d.jl
    echo "--- exit status $? at $(date '+%F %T') ---"
  ) > "$LOG/$name.log" 2>&1
}

for c in "${CASES[@]}"; do
  name=$(echo $c | cut -d' ' -f1)
  [ -f "$LOG/$name.log" ] && { echo "skip $name (log exists)"; continue; }
  with_slot run_case "$c"
  echo "$(date '+%F %T')  launch  $name  ($(slots_busy) busy, $(mem_free_gb) GB headroom)" | tee -a "$LOG/_batch.log"
done
wait
echo "$(date '+%F %T')  M1 tol re-run: finished" | tee -a "$LOG/_batch.log"
grep -L 'done:' "$LOG"/m3_*.log 2>/dev/null | sed 's/^/  NOT FINISHED: /' | tee -a "$LOG/_batch.log"
