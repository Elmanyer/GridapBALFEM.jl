#!/bin/bash
# ==============================================================
#  run_1d_m3_M1.sh — Chapter 8: the decay-triad test with ONE vertical element (P1LFE-1, 2026-11-10)
#
#  THE IDEA. The m3 growth is read as a decay triad: the carrier's second harmonic (n = 2, 2ω₀) feeds
#  a steady vertical-shear mode (ω ≈ 0) and a backward surface wave at ≈ 2ω₀. That needs a model
#  surface wave at 2ω₀. The surface branch saturates at ω∞ = √(g/d · DᵀB⁻¹D):
#    P1LFE-2 (M = 2): DᵀB⁻¹D = 22.0 → ω∞ = 7.86 rad/s → channel open for T ≥ 1.60 s (d = 3.5)
#    P1LFE-1 (M = 1): DᵀB⁻¹D = 8    → ω∞ = 4.74 rad/s → channel open for T ≥ 2.65 s
#  M = 1 still has a steady shear mode (u_bed = −u_surf: zero depth mean, flat surface), so the steady
#  member exists; what changes is whether a partner wave at 2ω₀ exists.
#
#  PREDICTIONS (written before the runs):
#    M1 T = 1.6 s   ω∞ / 2ω₀ = 0.60: channel CLOSED → no steady-shear growth. Growth here means the
#                   triad is not the whole story (other channels: n = N_cells lattice, the low-n
#                   carrier channel (1, ω₀) → (−1, ω₀) + (2, 0)).
#    M1 T = 3.0 s   ω∞ / 2ω₀ = 1.13: channel OPEN → growth of a steady shear mode + a backward wave
#                   near 2ω₀, as in M2 T = 1.6 s.
#  ⚠ M = 1 at T = 1.6 s is a poor wave model (λ 21 % short, kd 6.9 vs 5.5); the box is sized to the
#  model's own wavelength, so the test is self-consistent. At T = 3.0 s M = 1 is accurate (kd 1.69).
#
#  DESIGN: run_1d_m3_period.sh verbatim (m3: 𝓝 off, no CIP, flat bed; Q3/Q2; CN; dt = 0.04; seed 1e-8;
#  one wavelength; VTK every 0.2 s; 32 cells per wavelength) with M = 1 and T varied. Steepness
#  k₀A = 0.1582 fixed, k₀ from the P1LFE-1 dispersion relation:
#    T = 1.6 s: k₀ = 1.9795, λ = 3.174 m,  A = 0.0799 m
#    T = 3.0 s: k₀ = 0.4817, λ = 13.043 m, A = 0.3284 m (H/d = 0.19)
#  Duration 200 periods (320 s and 600 s): the growth at T = 3.0 s is expected to be slower.
#
#  SLOTS: run/local/_slot_throttle.sh (memory only). ANALYSIS: postprocessing/examples/periodic_viz.jl.
#  ⚠ Logs on persistent disk (rule 47); never edit this file while it runs (rule 38h). A case whose
#  log exists is skipped: delete the log to re-run it.
#
#  Detach: setsid nohup run/local/run_1d_m3_M1.sh > output/1D_M1_test/_logs/_batch.out 2>&1 < /dev/null &
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
  "m3_M1_T300_n32  1  3.0  0.3284  32"
  "m3_M1_T160_n32  1  1.6  0.0799  32"
)

# ---- source snapshot ---------------------------------------------------------------------------
SNAP=$LOG/source_$(date '+%Y%m%d_%H%M%S')
mkdir -p "$SNAP"
git rev-parse HEAD > "$SNAP/commit.txt"
git diff HEAD > "$SNAP/tracked.diff"
git ls-files --others --exclude-standard -z | tar czf "$SNAP/untracked.tgz" --null -T -
echo "$(date '+%F %T')  M1 test source snapshot: $SNAP (commit $(cut -c1-10 "$SNAP/commit.txt") + diff)" | tee -a "$LOG/_batch.log"

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
    systemd-run --user --scope --quiet --unit="m3M1-$name-$$" -p MemoryMax="$MEM_PER_CASE" \
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
echo "$(date '+%F %T')  M1 test: all cases finished" | tee -a "$LOG/_batch.log"
grep -L 'done:' "$LOG"/m3_*.log 2>/dev/null | sed 's/^/  NOT FINISHED: /' | tee -a "$LOG/_batch.log"
