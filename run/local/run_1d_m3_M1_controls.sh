#!/bin/bash
# ==============================================================
#  run_1d_m3_M1_controls.sh — the two one-axis controls of m3_M1_T300_n32 (2026-11-11)
#
#  WHY. m3_M1_T300_n32 (run_1d_m3_M1.sh: M = 1, T = 3.0 s, 32 cells/λ, CN, 600 s) stayed clean, against
#  the triad prediction (2ω₀/ω∞ = 0.88: a surface partner at 2ω₀ exists). But it differs from the
#  unstable campaign case (M = 2, T = 1.6 s, n32) in more than the triad channel: kd (1.69 vs 5.5),
#  dx in metres (0.41 vs 0.10), steps per period (75 vs 40), A (0.33 vs 0.10). Rule 14c: vary one axis.
#
#  CASES (each = m3_M1_T300_n32 with ONE change):
#    m3_M2_T300_n32   M = 2 instead of 1. Same T, kd, λ (13.09 vs 13.04 m), dx and steps per period, but
#                     the channel is CLOSED (2ω₀/ω∞ = 0.53). If the triad picture holds, neither grows; if
#                     M = 2 grows here, the growth needs the second layer, not the 2ω₀ partner.
#    m3_M1_T300_n64   64 cells/λ instead of 32: dx = 0.20 m, close to the 0.20 m of the n64 runs where the
#                     broadband grid-scale instability appeared (CH8_CAMPAIGN_TODO, 2026-11-09 interim).
#                     Tests whether the clean M1 T = 3.0 s result is a coarse-dx effect.
#  Steepness k₀A = 0.1582 kept; k₀ from each model's own dispersion relation (computed 2026-11-11):
#    P1LFE-1 T = 3.0: k₀ = 0.4817, λ = 13.043 m, A = 0.3284 m (as in run_1d_m3_M1.sh)
#    P1LFE-2 T = 3.0: k₀ = 0.4799, λ = 13.094 m, A = 0.3297 m
#  Both inside the 2 % range (M = 1: kd ≤ 3.28 C / 1.75 Cg; M = 2: 10.84 / 6.42).
#  nl_tol left at the driver default (1e-5), as in the reference run (which took 2 Newton iterations per
#  step); a run that drops to 1 iteration per step must be re-run at 1e-7 (see run_1d_m3_M1_tol.sh).
#
#  DESIGN otherwise run_1d_m3_M1.sh verbatim (m3: 𝓝 off, no CIP, flat bed; Q3/Q2; CN; dt = 0.04;
#  seed 1e-8; one wavelength; VTK every 0.2 s; 200 periods = 600 s).
#
#  SLOTS: run/local/_slot_throttle.sh (memory only). ANALYSIS: postprocessing (energy_history).
#  ⚠ Logs on persistent disk (rule 47); never edit this file while it runs (rule 38h). A case whose
#  log exists is skipped: delete the log to re-run it.
#
#  Detach: setsid nohup run/local/run_1d_m3_M1_controls.sh > output/1D_M1_test/_logs/_batch_controls.out 2>&1 < /dev/null &
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
  "m3_M1_T300_n64  1  3.0  0.3284  64"
  "m3_M2_T300_n32  2  3.0  0.3297  32"
)

# ---- source snapshot ---------------------------------------------------------------------------
SNAP=$LOG/source_$(date '+%Y%m%d_%H%M%S')
mkdir -p "$SNAP"
git rev-parse HEAD > "$SNAP/commit.txt"
git diff HEAD > "$SNAP/tracked.diff"
git ls-files --others --exclude-standard -z | tar czf "$SNAP/untracked.tgz" --null -T -
echo "$(date '+%F %T')  M1 controls source snapshot: $SNAP (commit $(cut -c1-10 "$SNAP/commit.txt") + diff)" | tee -a "$LOG/_batch.log"

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
    export BALFEM_SOLVER=theta BALFEM_USE_AD=0 BALFEM_NOISE=1e-8
    export BALFEM_CHECKPOINT_EVERY=250 BALFEM_CHECKPOINT_KEEP=2
    export BALFEM_SAVE_EVERY=5 BALFEM_DIAG_EVERY=5 BALFEM_PRINT_EVERY=25
    export BALFEM_OUTDIR=$PROJ/$OUT/$name
    export OPENBLAS_NUM_THREADS=2
    echo "case $name | m3: NL_PRESSURE=0 FLAT_BED=1, no CIP | M=$mvert T=$twave A=$awave NCELL=$ncell | Q3/Q2 | $PERIODS periods"
    echo "started $(date '+%F %T') | git $(git rev-parse --short HEAD)$(git diff --quiet -- src examples || echo ' + uncommitted edits')"
    systemd-run --user --scope --quiet --unit="m3M1c-$name-$$" -p MemoryMax="$MEM_PER_CASE" \
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
echo "$(date '+%F %T')  M1 controls: all cases finished" | tee -a "$LOG/_batch.log"
grep -L 'done:' "$LOG"/m3_M1_T300_n64.log "$LOG"/m3_M2_T300_n32.log 2>/dev/null | sed 's/^/  NOT FINISHED: /' | tee -a "$LOG/_batch.log"
