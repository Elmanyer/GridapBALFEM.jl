#!/bin/bash
# ==============================================================
#  run_1d_m3cip.sh — Chapter 8: model 3 WITH the first-order C⁰-CIP (campaign README, "Open questions")
#
#  THE QUESTION. m3 (campaign, no 𝓝, no CIP) and m5 (𝓝 + CIP 0.3) differ in TWO axes (rule 14c).
#  Each case here is the campaign m3 case with ONE change, γ_u = γ_η = 0.3 (DEFAULT_CIP_GAMMA), so
#  m3 ↔ m3_cip isolates the CIP and m3_cip ↔ m5 isolates 𝓝:
#    m3_q32_n8_cip    campaign m3 n8 grew in the high band (u, n 7–12) and rose E/E0 1.03 → 1.11
#                     from ~300 s; m5 n8 grew a long-wave vertical-structure mode from ~140 s
#  (n32 is the other arm named in the to-do; add it to CASES when wanted.)
#
#  DESIGN: identical to run_1d_stability_campaign.sh (P1LFE-2, Q3/Q2, A = 0.10, T = 1.6 s,
#  D = 3.5 m, dt = 0.04, CN, seed 1e-8, one wavelength, flat bed, 200 periods, VTK and diagnostics
#  every 5 steps) except BALFEM_CIP_GU = BALFEM_CIP_GE = 0.3 and the checkpoints: every 250 steps,
#  ALL kept (CHECKPOINT_KEEP=0), so a pre-onset restart is possible. Main tree (HEAD), which allows
#  the CIP without 𝓝 (setup_and_run attaches :jumpgrad whenever γ > 0).
#
#  SLOTS: run/local/_slot_throttle.sh (shared with the sponge sweep and the no-CIP controls).
#
#  ANALYSIS: PG_NCELL=<n> PG_PU=3 PG_NLAMBDA=1 postprocessing/examples/periodic_growth.jl on
#  output/1D_stability_campaign/<case>; compare with m3_q32_n<n> and m5_q32_n<n>.
#
#  ⚠ Logs on persistent disk (rule 47); never edit this file while it runs (rule 38h). A case whose
#  log exists is skipped: delete the log to re-run it.
#
#  Detach: nohup run/local/run_1d_m3cip.sh > output/1D_stability_campaign/_logs/_m3cip_batch.out 2>&1 &
# ==============================================================
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
PROJ="$(pwd)"
PERIODS="${PERIODS:-200}"
CIP_GAMMA="${CIP_GAMMA:-0.3}"
MEM_PER_CASE="${MEM_PER_CASE:-5G}"
OUT=output/1D_stability_campaign
LOG=$OUT/_logs
mkdir -p "$LOG"
source "$PROJ/run/local/_slot_throttle.sh"

# case : name | NCELL        (all: nonlinear, 𝓝 off, flat bed, γ = CIP_GAMMA, Q3/Q2, A = 0.10)
CASES=(
  "m3_q32_n8_cip    8"
)

# ---- source snapshot ---------------------------------------------------------------------------
SNAP=$LOG/source_m3cip_$(date '+%Y%m%d_%H%M%S')
mkdir -p "$SNAP"
git rev-parse HEAD > "$SNAP/commit.txt"
git diff HEAD > "$SNAP/tracked.diff"
git ls-files --others --exclude-standard -z | tar czf "$SNAP/untracked.tgz" --null -T -
echo "$(date '+%F %T')  m3+CIP source snapshot: $SNAP (commit $(cut -c1-10 "$SNAP/commit.txt") + diff)" | tee -a "$LOG/_batch.log"

run_case() {
  read -r name ncell <<< "$1"
  [ -f "$LOG/$name.log" ] && { echo "skip $name (log exists)"; return; }
  (
    export BALFEM_REGIME=nonlinear BALFEM_NL_PRESSURE=0 BALFEM_FLAT_BED=1
    export BALFEM_BED_SLOPE=0 BALFEM_BED_NW=1
    export BALFEM_CIP_GU=$CIP_GAMMA BALFEM_CIP_GE=$CIP_GAMMA
    export BALFEM_FE_ORDER=3 BALFEM_P_ETA=2 BALFEM_NCELL=$ncell BALFEM_NLAMBDA=1
    export BALFEM_AWAVE=0.10 BALFEM_TWAVE=1.6 BALFEM_D=3.5 BALFEM_DT=0.04 BALFEM_PERIODS=$PERIODS
    export BALFEM_SOLVER=theta BALFEM_USE_AD=0 BALFEM_NOISE=1e-8
    export BALFEM_CHECKPOINT_EVERY=250 BALFEM_CHECKPOINT_KEEP=0
    export BALFEM_SAVE_EVERY=5 BALFEM_DIAG_EVERY=5 BALFEM_PRINT_EVERY=25
    export BALFEM_OUTDIR=$PROJ/$OUT/$name
    export OPENBLAS_NUM_THREADS=2
    echo "case $name | NL_PRESSURE=0 FLAT_BED=1 NCELL=$ncell | Q3/Q2 | γ_u = γ_η = $CIP_GAMMA | $PERIODS periods"
    echo "started $(date '+%F %T') | git $(git rev-parse --short HEAD)$(git diff --quiet -- src examples || echo ' + uncommitted edits')"
    systemd-run --user --scope --quiet --unit="m3cip1d-$name-$$" -p MemoryMax="$MEM_PER_CASE" \
      stdbuf -oL -eL julia --project="$PROJ" examples/local_1d/run_periodic_1d.jl
    echo "--- exit status $? at $(date '+%F %T') ---"
  ) > "$LOG/$name.log" 2>&1
}

for c in "${CASES[@]}"; do
  name=$(echo $c | cut -d' ' -f1)
  [ -f "$LOG/$name.log" ] && { echo "skip $name (log exists)"; continue; }
  with_slot run_case "$c"
  echo "$(date '+%F %T')  launch  $name  (m3 + CIP; $(slots_busy) busy, $(mem_free_gb) GB free)" | tee -a "$LOG/_batch.log"
done
wait
echo "$(date '+%F %T')  m3+CIP cases finished" | tee -a "$LOG/_batch.log"
