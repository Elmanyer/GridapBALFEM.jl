#!/bin/bash
# ==============================================================
#  run_1d_sdirk_box.sh — Chapter 8: does a dissipative integrator stabilise the box? (2026-11-10)
#
#  THE QUESTION. Under Crank–Nicolson, m3 (𝓝 off, no CIP) diverges at 246 s and m5 (𝓝 + CIP 0.3) bursts
#  at ≈ 245 s in the one-wavelength box at 32 cells/λ (output/1D_stability_campaign/README.md). Does the
#  numerical dissipation of SDIRK_2_2 (Gridap's DIRK22(1,0,1): A-stable, not L-stable) stop that?
#
#  PREDICTIONS (written before the runs). SDIRK_2_2 removes ≈ 0.75 (ωΔt)⁴ of a mode's amplitude per step:
#    carrier ω0 (ωΔt = 0.157): ≈ 0.011 s⁻¹ (measured 0.010–0.012 s⁻¹ in every v1 SDIRK box run)
#    triad fast member 2ω0 (0.31): ≈ 0.18 s⁻¹ (16× the carrier); steady member (ω = 0): none (R(0) = 1)
#  A decay triad with one damped member still grows, at −γ/2 + √(γ²/4 + σ0²) ≈ 0.03 s⁻¹ (from 0.075), and the
#  pump fades with the carrier (A ∝ e^{−0.011 t}: ≈ 3 % of the amplitude left after 320 s). Expected: m3
#  survives 200 periods, but mostly because the forcing wave dies, not because the instability is
#  controlled; read the carrier and the disturbance together (plot_energy_history, modes = :groups).
#  ⚠ Rule 15: no stability claim rests on SDIRK alone; these runs measure what the dissipation does, paired
#  with the Crank–Nicolson runs of the campaign (identical except BALFEM_SOLVER/TABLEAU).
#
#  DESIGN: run_1d_stability_campaign.sh verbatim (Q3/Q2, A = 0.10, T = 1.6 s, dt = 0.04, seed 1e-8, one
#  wavelength, flat bed, 200 periods, CIP = the driver default: 0.3 with 𝓝, none without), with
#  BALFEM_SOLVER=sdirk BALFEM_TABLEAU=SDIRK_2_2. Two implicit stages per step: ≈ 2× the CN cost
#  (CN: m3 n32 5 s/step, m5 n32 58 s/step during the throttled campaign).
#
#  SLOTS: run/local/_slot_throttle.sh (memory only). ⚠ Logs on persistent disk (rule 47); never edit this
#  file while it runs (rule 38h). A case whose log exists is skipped.
#  Detach: setsid nohup run/local/run_1d_sdirk_box.sh > output/1D_stability_campaign/_logs/_batch_sdirk.out 2>&1 < /dev/null &
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
echo "$(date '+%F %T')  SDIRK box source snapshot: $SNAP (commit $(cut -c1-10 "$SNAP/commit.txt") + diff)" | tee -a "$LOG/_batch.log"

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

# case : name | NL_PRESSURE | NCELL        longest first
CASES=(
  "m5_q32_n32_sdirk   1 32"
  "m3_q32_n32_sdirk   0 32"
)

run_case() {
  read -r name nlp ncell <<< "$1"
  [ -f "$LOG/$name.log" ] && { echo "skip $name (log exists)"; return; }
  (
    unset BALFEM_CIP_GU BALFEM_CIP_GE
    export BALFEM_REGIME=nonlinear BALFEM_NL_PRESSURE=$nlp BALFEM_FLAT_BED=1
    export BALFEM_BED_SLOPE=0 BALFEM_BED_NW=1
    export BALFEM_FE_ORDER=3 BALFEM_P_ETA=2 BALFEM_NCELL=$ncell BALFEM_NLAMBDA=1
    export BALFEM_AWAVE=0.10 BALFEM_TWAVE=1.6 BALFEM_D=3.5 BALFEM_DT=0.04 BALFEM_PERIODS=$PERIODS
    export BALFEM_SOLVER=sdirk BALFEM_TABLEAU=SDIRK_2_2 BALFEM_USE_AD=0 BALFEM_NOISE=1e-8
    export BALFEM_CHECKPOINT_EVERY=$CHECKPOINT_EVERY BALFEM_CHECKPOINT_KEEP=2
    export BALFEM_SAVE_EVERY=5 BALFEM_DIAG_EVERY=5 BALFEM_PRINT_EVERY=25
    export BALFEM_OUTDIR=$PROJ/$OUT/$name
    export OPENBLAS_NUM_THREADS=2
    echo "case $name | NL_PRESSURE=$nlp FLAT_BED=1 NCELL=$ncell | Q3/Q2 | SDIRK_2_2 | $PERIODS periods"
    echo "started $(date '+%F %T') | git $(git rev-parse --short HEAD)$(git diff --quiet -- src examples || echo ' + uncommitted edits')"
    systemd-run --user --scope --quiet --unit="sdirk1d-$name-$$" -p MemoryMax="$MEM_PER_CASE" \
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
echo "$(date '+%F %T')  SDIRK box: all cases finished" | tee -a "$LOG/_batch.log"
grep -L 'done:' "$LOG"/m[35]_q32_n32_sdirk.log 2>/dev/null | sed 's/^/  NOT FINISHED: /' | tee -a "$LOG/_batch.log"
