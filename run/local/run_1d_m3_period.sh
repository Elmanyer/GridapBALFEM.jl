#!/bin/bash
# ==============================================================
#  run_1d_m3_period.sh — Chapter 8: the wave-period test of the decay-triad mechanism (2026-11-08)
#
#  THE HYPOTHESIS (output/1D_stability_campaign/README.md, "What the m3 growth oscillates at"):
#  the m3 growth is a decay triad of the carrier's second harmonic (n = 2, 2ω₀) into a steady
#  vertical-shear mode (n_s, ω ≈ 0) and a backward surface wave (n_s − 2, ω ≈ 2ω₀). It needs a
#  surface wave at 2ω₀, and the P1LFE-2 surface branch saturates at ω∞ = 7.84 rad/s (d = 3.5).
#  At T = 1.6 s, 2ω₀ = 7.854 sits on that plateau.
#
#  PREDICTIONS (written before the runs; continuum detuning ω(k_n) − 2ω₀ from temporal_spectrum.jl):
#    T = 1.4 s   2ω₀ = 8.98 rad/s, 14 % above ω∞: channel CLOSED  → no steady short-wave growth
#    T = 1.5 s   2ω₀ = 8.38 rad/s,  7 % above ω∞: channel CLOSED  → no steady short-wave growth
#    T = 1.7 s   2ω₀ = 7.39 rad/s: OPEN at low n, fast n = 7 detuned −0.03 rad/s
#                → growth of the triad steady n = 9 + backward n = 7, at every mesh
#    T = 1.8 s   2ω₀ = 6.98 rad/s: open but detuned (n = 6 +0.10, n = 5 −0.16 rad/s) → weaker or none
#    control T = 1.6 s, 64 cells: must reproduce the October divergence (44.6 s)
#  Not covered by the prediction: the n = N_cells lattice channel (seen at n8) and the low-n
#  carrier channel (1, ω₀) → (−1, ω₀) + (2, 0), which are open at every T.
#
#  DESIGN: the campaign m3 case (P1LFE-2, Q3/Q2, CN, dt = 0.04, seed 1e-8, one wavelength, flat
#  bed, 𝓝 off, no CIP, VTK every 0.2 s) with T varied and, to keep the nonlinearity the same,
#  the steepness k₀A = 0.1582 fixed (A from the model wavenumber). Duration 320 s for every T
#  (growth rates are per second), so PERIODS = 320/T. Cells per wavelength 16, 32, 64.
#
#  SLOTS: run/local/_slot_throttle.sh (shared with every local 1-D launcher).
#  ANALYSIS: periodic_growth.jl (PG_NCELL, PG_PU=3, PG_NLAMBDA=1) for the bands;
#    output/1D_stability_campaign/_logs/temporal_spectrum.jl for the triad test.
#
#  ⚠ Logs on persistent disk (rule 47); never edit this file while it runs (rule 38h). A case whose
#  log exists is skipped: delete the log to re-run it.
#
#  Detach: MAXP=10 setsid nohup run/local/run_1d_m3_period.sh > output/1D_period_test/_logs/_batch.out 2>&1 &
# ==============================================================
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
PROJ="$(pwd)"
DURATION="${DURATION:-320}"
MEM_PER_CASE="${MEM_PER_CASE:-5G}"
OUT=output/1D_period_test
LOG=$OUT/_logs
mkdir -p "$LOG"
source "$PROJ/run/local/_slot_throttle.sh"

# case : name | T [s] | A [m] (k₀A = 0.1582) | NCELL        longest first (LPT); control first
CASES=(
  "m3_T160_n64  1.6  0.1000  64"
  "m3_T140_n64  1.4  0.0762  64"
  "m3_T150_n64  1.5  0.0877  64"
  "m3_T170_n64  1.7  0.1129  64"
  "m3_T180_n64  1.8  0.1266  64"
  "m3_T140_n32  1.4  0.0762  32"
  "m3_T150_n32  1.5  0.0877  32"
  "m3_T170_n32  1.7  0.1129  32"
  "m3_T180_n32  1.8  0.1266  32"
  "m3_T140_n16  1.4  0.0762  16"
  "m3_T150_n16  1.5  0.0877  16"
  "m3_T170_n16  1.7  0.1129  16"
  "m3_T180_n16  1.8  0.1266  16"
)

# ---- source snapshot ---------------------------------------------------------------------------
SNAP=$LOG/source_$(date '+%Y%m%d_%H%M%S')
mkdir -p "$SNAP"
git rev-parse HEAD > "$SNAP/commit.txt"
git diff HEAD > "$SNAP/tracked.diff"
git ls-files --others --exclude-standard -z | tar czf "$SNAP/untracked.tgz" --null -T -
echo "$(date '+%F %T')  period test source snapshot: $SNAP (commit $(cut -c1-10 "$SNAP/commit.txt") + diff)" | tee -a "$LOG/_batch.log"

run_case() {
  read -r name twave awave ncell <<< "$1"
  [ -f "$LOG/$name.log" ] && { echo "skip $name (log exists)"; return; }
  local periods
  periods=$(python3 -c "print($DURATION / $twave)")
  (
    unset BALFEM_CIP_GU BALFEM_CIP_GE
    export BALFEM_REGIME=nonlinear BALFEM_NL_PRESSURE=0 BALFEM_FLAT_BED=1
    export BALFEM_BED_SLOPE=0 BALFEM_BED_NW=1
    export BALFEM_FE_ORDER=3 BALFEM_P_ETA=2 BALFEM_NCELL=$ncell BALFEM_NLAMBDA=1
    export BALFEM_AWAVE=$awave BALFEM_TWAVE=$twave BALFEM_D=3.5 BALFEM_DT=0.04 BALFEM_PERIODS=$periods
    export BALFEM_SOLVER=theta BALFEM_USE_AD=0 BALFEM_NOISE=1e-8
    export BALFEM_CHECKPOINT_EVERY=250 BALFEM_CHECKPOINT_KEEP=2
    export BALFEM_SAVE_EVERY=5 BALFEM_DIAG_EVERY=5 BALFEM_PRINT_EVERY=25
    export BALFEM_OUTDIR=$PROJ/$OUT/$name
    export OPENBLAS_NUM_THREADS=2
    echo "case $name | m3: NL_PRESSURE=0 FLAT_BED=1, no CIP | T=$twave A=$awave NCELL=$ncell | Q3/Q2 | $periods periods = $DURATION s"
    echo "started $(date '+%F %T') | git $(git rev-parse --short HEAD)$(git diff --quiet -- src examples || echo ' + uncommitted edits')"
    systemd-run --user --scope --quiet --unit="m3T1d-$name-$$" -p MemoryMax="$MEM_PER_CASE" \
      stdbuf -oL -eL julia --project="$PROJ" examples/local_1d/run_periodic_1d.jl
    echo "--- exit status $? at $(date '+%F %T') ---"
  ) > "$LOG/$name.log" 2>&1
}

for c in "${CASES[@]}"; do
  name=$(echo $c | cut -d' ' -f1)
  [ -f "$LOG/$name.log" ] && { echo "skip $name (log exists)"; continue; }
  with_slot run_case "$c"
  echo "$(date '+%F %T')  launch  $name  ($(slots_busy) busy, $(mem_free_gb) GB free)" | tee -a "$LOG/_batch.log"
done
wait
echo "$(date '+%F %T')  period test: all cases finished" | tee -a "$LOG/_batch.log"
grep -L 'done:' "$LOG"/m3_*.log 2>/dev/null | sed 's/^/  NOT FINISHED: /' | tee -a "$LOG/_batch.log"
