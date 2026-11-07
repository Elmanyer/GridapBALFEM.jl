#!/bin/bash
# ==============================================================
#  run_1d_nocip_controls.sh — Chapter 8: the NO-CIP controls of the 1-D stability campaign
#
#  THE QUESTION. The campaign (run_1d_stability_campaign.sh) ran models 5/6 WITH the mandatory
#  first-order C⁰-CIP (γ_u = γ_η = 0.3). Is the CIP what keeps them stable? Each control here is
#  a campaign case with ONE change, γ_u = γ_η = 0 (rule 14c):
#    m6_q32_n{8,16,32,64}_nocip  no model-6 run without CIP exists (the periodic bed postdates
#                                the mandatory CIP)
#    m5_q32_n8_nocip             the October no-CIP m5 n8 run stopped at 160 s, bounded, while the
#                                campaign m5 n8 (CIP) grows from ~120 s: extend it to 320 s
#  The m5 n16/n32/n64 controls exist (output/local_1d/periodic/v2_p32_nlp1_n*_cn, Newton NaN at
#  81 / 27 / 9 s) and are not repeated.
#
#  WHY A WORKTREE. The current code refuses 𝓝 without the CIP (CLAUDE.md rule 1c: no back door
#  in the solver). These runs use a separate git worktree, ../GridapBALFEM_nocip.wt (branch
#  diag/nocip-control):
#    base 68919a3 (the parent of faa094f, which made the CIP mandatory), plus UNCOMMITTED
#      * the checkpoint feature (725b9ed, cherry-picked; driver conflicts resolved by keeping the
#        old γ default 0 and adding the `resume` name token),
#      * the periodic cosine bed of the box driver (df92a26).
#  At 68919a3, γ = 0 with 𝓝 keeps the bare skeleton geometry, so the kind-S facet terms are still
#  assembled. (Current HEAD with only the refusal removed would set skel = nothing and DROP them,
#  i.e. it would change two things at once.) The 𝓝 residual is otherwise identical to the
#  campaign's: ea78e64 → 1ae2695 touched only CIP enforcement, distributed facet bookkeeping and
#  the CIP-power diagnostic (README of the campaign, 2026-11-06).
#  The worktree has Manifest.toml copied from the main tree and symlinks Gridap.jl → the fork and
#  output → the main output/, so every result lands in the main output tree.
#
#  DESIGN: identical to run_1d_stability_campaign.sh (P1LFE-2, Q3/Q2, A = 0.10, T = 1.6 s,
#  D = 3.5 m, dt = 0.04, CN, seed 1e-8, one wavelength, 200 periods, BED_SLOPE 0.1 for m6, VTK and
#  diagnostics every 5 steps) except γ = 0 and the checkpoints: every 250 steps, ALL kept
#  (CHECKPOINT_KEEP=0), so a diagnostic restart from before an onset is always possible (the
#  campaign kept only 2 and lost every pre-onset state).
#  The divergence guard (|η| > 20 η_ref = 2 m) ends a failing run; that time is the measurement.
#
#  SLOTS: run/local/_slot_throttle.sh (shared with the sponge sweep): MAXP processes in total and
#  MemAvailable ≥ MIN_FREE_GB. Longest first (LPT, rule 43): the n8 cases may run the full 320 s
#  (~2–3 days); n16–n64 are expected to fail within minutes of simulated time.
#
#  ANALYSIS per run: PG_NCELL=<n> PG_PU=3 PG_NLAMBDA=1, postprocessing/examples/periodic_growth.jl
#  on output/1D_stability_campaign/<case>; compare with the campaign case without `_nocip`.
#
#  ⚠ Logs on persistent disk (rule 47); never edit this file while it runs (rule 38h). A case whose
#  log exists is skipped: delete the log to re-run it. Do not edit the worktree while cases run.
#
#  USAGE:  run/local/run_1d_nocip_controls.sh                     (MAXP=3 MIN_FREE_GB=6 defaults)
#  Detach: nohup run/local/run_1d_nocip_controls.sh > output/1D_stability_campaign/_logs/_nocip_batch.out 2>&1 &
#  Smoke:  SMOKE=1 run/local/run_1d_nocip_controls.sh   (0.25 periods = 10 steps of m6 n16 into
#          output/1D_stability_campaign/_smoke/, compare with the first rows of the campaign case)
# ==============================================================
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
PROJ="$(pwd)"
WT="$(cd "$PROJ/../GridapBALFEM_nocip.wt" 2>/dev/null && pwd)" ||
  { echo "worktree ../GridapBALFEM_nocip.wt not found"; exit 1; }
[ "$(git -C "$WT" rev-parse HEAD)" = "$(git rev-parse 68919a3)" ] ||
  { echo "worktree HEAD is not 68919a3"; exit 1; }
[ -f "$WT/src/checkpoint.jl" ] && [ -f "$WT/Manifest.toml" ] && [ -e "$WT/Gridap.jl" ] ||
  { echo "worktree is missing checkpoint.jl, Manifest.toml or the Gridap.jl link"; exit 1; }

PERIODS="${PERIODS:-200}"
BED_SLOPE="${BED_SLOPE:-0.1}"
MEM_PER_CASE="${MEM_PER_CASE:-5G}"
OUT=output/1D_stability_campaign
LOG=$OUT/_logs
mkdir -p "$LOG"
source "$PROJ/run/local/_slot_throttle.sh"

# case : name | FLAT_BED | NCELL        (all: nonlinear, 𝓝 on, γ = 0, Q3/Q2, A = 0.10)
CASES=(
  "m5_q32_n8_nocip    1 8"
  "m6_q32_n8_nocip    0 8"
  "m6_q32_n64_nocip   0 64"
  "m6_q32_n32_nocip   0 32"
  "m6_q32_n16_nocip   0 16"
)
if [ "${SMOKE:-0}" = 1 ]; then
  PERIODS=0.25
  OUT=$OUT/_smoke
  LOG=$OUT
  mkdir -p "$LOG"
  CASES=("m6_q32_n16_nocip_smoke   0 16")
fi

# ---- source snapshot of the worktree: base commit + uncommitted diff + untracked files ------
SNAP=$LOG/source_nocip_$(date '+%Y%m%d_%H%M%S')
mkdir -p "$SNAP"
git -C "$WT" rev-parse HEAD > "$SNAP/commit.txt"
git -C "$WT" diff HEAD > "$SNAP/tracked.diff"
( cd "$WT" && git ls-files --others --exclude-standard -z -- src examples test \
    | tar czf "$PROJ/$SNAP/untracked.tgz" --null -T - )
echo "$(date '+%F %T')  no-CIP worktree snapshot: $SNAP (68919a3 + diff + untracked)" | tee -a "$LOG/_batch.log"

run_case() {
  read -r name flat ncell <<< "$1"
  [ -f "$LOG/$name.log" ] && { echo "skip $name (log exists)"; return; }
  local slope=0
  [ "$flat" = 0 ] && slope=$BED_SLOPE
  (
    export BALFEM_REGIME=nonlinear BALFEM_NL_PRESSURE=1 BALFEM_FLAT_BED=$flat
    export BALFEM_BED_SLOPE=$slope BALFEM_BED_NW=1
    export BALFEM_CIP_GU=0 BALFEM_CIP_GE=0
    export BALFEM_FE_ORDER=3 BALFEM_P_ETA=2 BALFEM_NCELL=$ncell BALFEM_NLAMBDA=1
    export BALFEM_AWAVE=0.10 BALFEM_TWAVE=1.6 BALFEM_D=3.5 BALFEM_DT=0.04 BALFEM_PERIODS=$PERIODS
    export BALFEM_SOLVER=theta BALFEM_USE_AD=0 BALFEM_NOISE=1e-8
    export BALFEM_CHECKPOINT_EVERY=250 BALFEM_CHECKPOINT_KEEP=0
    export BALFEM_SAVE_EVERY=5 BALFEM_DIAG_EVERY=5 BALFEM_PRINT_EVERY=25
    export BALFEM_OUTDIR=$PROJ/$OUT/$name
    export OPENBLAS_NUM_THREADS=2
    echo "case $name | NL_PRESSURE=1 FLAT_BED=$flat BED_SLOPE=$slope NCELL=$ncell | Q3/Q2 | γ = 0 | $PERIODS periods"
    echo "started $(date '+%F %T') | worktree $WT @ $(git -C "$WT" rev-parse --short HEAD) + uncommitted (checkpoint, periodic bed)"
    cd "$WT" &&
    systemd-run --user --scope --quiet --unit="nocip1d-$name-$$" -p MemoryMax="$MEM_PER_CASE" \
      stdbuf -oL -eL julia --project="$WT" examples/local_1d/run_periodic_1d.jl
    echo "--- exit status $? at $(date '+%F %T') ---"
  ) > "$LOG/$name.log" 2>&1
}

for c in "${CASES[@]}"; do
  name=$(echo $c | cut -d' ' -f1)
  [ -f "$LOG/$name.log" ] && { echo "skip $name (log exists)"; continue; }
  with_slot run_case "$c"
  echo "$(date '+%F %T')  launch  $name  (no-CIP worktree; $(slots_busy) busy, $(mem_free_gb) GB free)" | tee -a "$LOG/_batch.log"
done
wait
echo "$(date '+%F %T')  no-CIP controls finished" | tee -a "$LOG/_batch.log"
