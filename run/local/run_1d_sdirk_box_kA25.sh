#!/bin/bash
# ==============================================================
#  run_1d_sdirk_box_kA25.sh — the SDIRK box pair at a steeper carrier, k₀A = 0.25 (2026-11-11)
#
#  THE QUESTION. Under SDIRK_2_2, m3 at k₀A = 0.158 (run_1d_sdirk_box.sh) keeps the CN-unstable modes
#  n ≈ 14–18 at seed level to ≥ 170 s, while under CN they blow up at 246 s. Does the integrator's
#  dissipation still hold when the wave is more nonlinear? Each case = the run_1d_sdirk_box.sh case with
#  ONE change: A = 0.10 → 0.1581 m (k₀A 0.158 → 0.25; H/λ 0.050 → 0.080, 36 % → 56 % of the Miche
#  limit 0.142; second-order bound content ∝ (k₀A)² ×2.5). k₀ = 1.5818 (P1LFE-2, T = 1.6 s, d = 3.5).
#
#  PREDICTIONS (written before the runs). SDIRK_2_2 damps amplitude ≈ 0.36 (ωΔt)⁴ per step: the carrier
#  0.0054 s⁻¹, a 2ω₀ member 0.076 s⁻¹, a steady member nothing. The pump of the triad grows with the
#  carrier, ×1.6 (∝ k₀A) to ×2.5 (∝ (k₀A)²) at the start. At k₀A = 0.158 the measured early growth of n16
#  under SDIRK was 0.008 s⁻¹ (CN: 0.073), well below the simple triad estimate (0.046), so the margin is
#  unknown: growth here that was absent at 0.158 means the dissipation's margin is amplitude-limited.
#  The linear-eigenmode IC sheds more free harmonics at this steepness (expected larger "bound harmonics").
#  ⚠ Rule 15: no CN pair exists yet at k₀A = 0.25; these runs say what the dissipation does, not whether
#  the model is stable.
#
#  DESIGN: run_1d_sdirk_box.sh verbatim (Q3/Q2, T = 1.6 s, dt = 0.04, seed 1e-8, one wavelength, flat
#  bed, 200 periods, SDIRK_2_2, CIP = driver default: 0.3 with 𝓝, none without) except A.
#
#  SLOTS: run/local/_slot_throttle.sh (memory only). ⚠ Logs on persistent disk (rule 47); never edit
#  this file while it runs (rule 38h). A case whose log exists is skipped.
#  Detach: setsid nohup run/local/run_1d_sdirk_box_kA25.sh > output/1D_stability_campaign/_logs/_batch_sdirk_kA25.out 2>&1 < /dev/null &
# ==============================================================
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
PROJ="$(pwd)"
MEM_PER_CASE="${MEM_PER_CASE:-5G}"
PERIODS="${PERIODS:-200}"
CHECKPOINT_EVERY="${CHECKPOINT_EVERY:-250}"
OUT=output/1D_stability_campaign
LOG=$OUT/_logs
mkdir -p "$LOG"
source "$PROJ/run/local/_slot_throttle.sh"

SNAP=$LOG/source_$(date '+%Y%m%d_%H%M%S')
mkdir -p "$SNAP"
git rev-parse HEAD > "$SNAP/commit.txt"
git diff HEAD > "$SNAP/tracked.diff"
git ls-files --others --exclude-standard -z | tar czf "$SNAP/untracked.tgz" --null -T -
echo "$(date '+%F %T')  SDIRK kA25 box source snapshot: $SNAP (commit $(cut -c1-10 "$SNAP/commit.txt") + diff)" | tee -a "$LOG/_batch.log"

IMAGE="$PROJ/GridapBALFEM_sysimage_local.so"
src_hash() { ( cd "$PROJ/src" && find . -name '*.jl' -type f -print0 \
               | LC_ALL=C sort -z | xargs -0 -r sha256sum ) | sha256sum | awk '{print $1}'; }
JFLAGS=()
if [ -f "$IMAGE" ] && [ -f "$IMAGE.src.sha256" ] && [ "$(cat "$IMAGE.src.sha256")" = "$(src_hash)" ]; then
  JFLAGS=(-J "$IMAGE")
  echo "$(date '+%F %T')  sysimage: $IMAGE (src stamp matches)" | tee -a "$LOG/_batch.log"
else
  echo "$(date '+%F %T')  sysimage: NONE or STALE — every case JITs" | tee -a "$LOG/_batch.log"
fi

# case : name | NL_PRESSURE | NCELL | A [m]        longest first
CASES=(
  "m5_q32_n32_sdirk_kA25   1 32 0.1581"
  "m3_q32_n32_sdirk_kA25   0 32 0.1581"
)

run_case() {
  read -r name nlp ncell awave <<< "$1"
  [ -f "$LOG/$name.log" ] && { echo "skip $name (log exists)"; return; }
  (
    unset BALFEM_CIP_GU BALFEM_CIP_GE
    export BALFEM_REGIME=nonlinear BALFEM_NL_PRESSURE=$nlp BALFEM_FLAT_BED=1
    export BALFEM_BED_SLOPE=0 BALFEM_BED_NW=1
    export BALFEM_FE_ORDER=3 BALFEM_P_ETA=2 BALFEM_NCELL=$ncell BALFEM_NLAMBDA=1
    export BALFEM_AWAVE=$awave BALFEM_TWAVE=1.6 BALFEM_D=3.5 BALFEM_DT=0.04 BALFEM_PERIODS=$PERIODS
    export BALFEM_SOLVER=sdirk BALFEM_TABLEAU=SDIRK_2_2 BALFEM_USE_AD=0 BALFEM_NOISE=1e-8
    export BALFEM_CHECKPOINT_EVERY=$CHECKPOINT_EVERY BALFEM_CHECKPOINT_KEEP=2
    export BALFEM_SAVE_EVERY=5 BALFEM_DIAG_EVERY=5 BALFEM_PRINT_EVERY=25
    export BALFEM_OUTDIR=$PROJ/$OUT/$name
    export OPENBLAS_NUM_THREADS=2
    echo "case $name | NL_PRESSURE=$nlp FLAT_BED=1 NCELL=$ncell A=$awave | Q3/Q2 | SDIRK_2_2 | $PERIODS periods"
    echo "started $(date '+%F %T') | git $(git rev-parse --short HEAD)$(git diff --quiet -- src examples || echo ' + uncommitted edits')"
    systemd-run --user --scope --quiet --unit="sdirkA25-$name-$$" -p MemoryMax="$MEM_PER_CASE" \
      stdbuf -oL -eL julia "${JFLAGS[@]}" --project="$PROJ" examples/local_1d/run_periodic_1d.jl
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
echo "$(date '+%F %T')  SDIRK kA25 box: all cases finished" | tee -a "$LOG/_batch.log"
grep -L 'done:' "$LOG"/m[35]_q32_n32_sdirk_kA25.log 2>/dev/null | sed 's/^/  NOT FINISHED: /' | tee -a "$LOG/_batch.log"
