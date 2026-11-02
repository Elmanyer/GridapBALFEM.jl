#!/bin/bash
# ==============================================================
#  run_1d_stability_campaign.sh — Chapter 8: stability of the FULL nonlinear solver, closed 1-D box
#
#  THE QUESTION. Is the solver as defined after chapter 7 stable WITHOUT any further stabilisation?
#  The solver is 𝓝 (all eight components) + the Class-III kind-S facet terms + the MANDATORY
#  first-order C⁰-CIP (γ_u = γ_η = DEFAULT_CIP_GAMMA = 0.3, the driver default with 𝓝). Nothing
#  else: no filter, no extra penalty.
#
#  ARMS (Q3/Q2 only: 𝓝 refuses p_η < 2), at 8 / 16 / 32 / 64 cells per wavelength:
#    m5  model 5: nonlinear, 𝓝, flat bed
#    m6  model 6: nonlinear, 𝓝, VARIABLE bed — a periodic cosine bed d(x) = 3.5 + δ cos(k0 x),
#        one undulation per box (k_bed = k0: not Bragg-resonant), max slope BED_SLOPE (1:10,
#        δ ≈ 0.063 m). Mild by design: it exercises every ∇h and bed-Hessian term of 𝓝 without
#        making the bed the experiment.
#    m3  model 3: nonlinear, 𝓝 OFF, flat bed, NO skeleton terms — the null control in the same
#        batch (rule 38c), at 8 / 16 / 32 only. The 2026-10-06…09 ladder (100 periods, copied in
#        output/1D_stability_campaign/model3_reference/) left n32 ACCELERATING at t = 160 s
#        (mid/high/sub-element bands +0.10–0.12 s⁻¹ over the last 60 s) and n16 growing slowly
#        (+0.02–0.03); n64 diverged at 44.6 s and is not repeated.
#
#  DESIGN (rules 12c, 14c, 15, 38c, 38i): identical to run/local/run_1dper_v2_unstab_q32.sh, so
#  every arm compares with the model-3 reference: P1LFE-2, A = 0.10, T = 1.6 s, D = 3.5 m,
#  dt = 0.04, CRANK–NICOLSON, seed 1e-8 on every harmonic, one wavelength in the box, VTK and
#  diagnostics every 5 steps — except the duration: 200 PERIODS (320 s), to let a slow late
#  escape (v1 first-order CIP: +0.13 s⁻¹) show.
#
#  ANALYSIS per run (growth rate per band; read with x_at_max, rule 40):
#    PG_NCELL=<n> PG_PU=3 julia --project=postprocessing \
#        postprocessing/examples/periodic_growth.jl output/1D_stability_campaign/<case>
#
#  EXECUTION. Sequential LU, one process per case (a few hundred DOFs: rule 23), each in its own
#  systemd scope capped at MEM_PER_CASE so an OOM kill takes one case, never the session
#  (memory: 12 ranks × ~2.8 GB killed the terminal on 2026-10-30). Longest cases first (LPT,
#  rule 43). Uses the local sysimage when its src/ stamp matches the current src/ (rule 25),
#  otherwise plain JIT (and says so).
#  CHECKPOINTS every CHECKPOINT_EVERY steps (default 250 = 10 s simulated, 2 kept): a case lost to
#  a crash, OOM or reboot resumes with BALFEM_RESTART_FROM=<its dir> (u_n only, src/checkpoint.jl).
#  SOURCE SNAPSHOT: the commit, the tracked diff and the untracked files go to _logs/source_<date>/,
#  because the runs are launched from an uncommitted tree.
#  ⚠ Logs on persistent disk (rule 47); never edit this file while it runs (rule 38h). A case
#  whose log exists is skipped: delete the log to re-run it.
#
#  USAGE:  run/local/run_1d_stability_campaign.sh [MAX_PARALLEL]      (default 6)
#  Detach: nohup run/local/run_1d_stability_campaign.sh > output/1D_stability_campaign/_logs/_batch.out 2>&1 &
# ==============================================================
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
PROJ="$(pwd)"
MAXP="${1:-6}"
MEM_PER_CASE="${MEM_PER_CASE:-5G}"
PERIODS="${PERIODS:-200}"
BED_SLOPE="${BED_SLOPE:-0.1}"
CHECKPOINT_EVERY="${CHECKPOINT_EVERY:-250}"   # steps (= 10 s simulated); keep the 2 newest
OUT=output/1D_stability_campaign
LOG=$OUT/_logs
mkdir -p "$LOG"

# ---- source snapshot: the runs record a DIRTY tree, so keep exactly what they ran -----------
SNAP=$LOG/source_$(date '+%Y%m%d_%H%M%S')
mkdir -p "$SNAP"
git rev-parse HEAD > "$SNAP/commit.txt"
git diff HEAD > "$SNAP/tracked.diff"
git ls-files --others --exclude-standard -z | tar czf "$SNAP/untracked.tgz" --null -T -
echo "$(date '+%F %T')  source snapshot: $SNAP (commit $(cut -c1-10 "$SNAP/commit.txt") + diff + untracked)" | tee -a "$LOG/_batch.log"

# ---- sysimage: use it only if built from exactly this src/ ----------------------------------
IMAGE="$PROJ/GridapBALFEM_sysimage_local.so"
src_hash() { ( cd "$PROJ/src" && find . -name '*.jl' -type f -print0 \
               | LC_ALL=C sort -z | xargs -0 -r sha256sum ) | sha256sum | awk '{print $1}'; }
JFLAGS=()
if [ -f "$IMAGE" ] && [ -f "$IMAGE.src.sha256" ] && [ "$(cat "$IMAGE.src.sha256")" = "$(src_hash)" ]; then
  JFLAGS=(-J "$IMAGE")
  echo "$(date '+%F %T')  sysimage: $IMAGE (src stamp matches)" | tee -a "$LOG/_batch.log"
else
  echo "$(date '+%F %T')  sysimage: NONE or STALE — every case JITs (~25 min each)" | tee -a "$LOG/_batch.log"
fi

# case : name | NL_PRESSURE | FLAT_BED | NCELL        (Q3/Q2, A = 0.10 throughout)
#  Order: longest first (LPT) — 64 cells/λ costs the most per step; model 3 (~4 s/step) last.
CASES=(
  "m5_q32_n64   1 1 64"
  "m6_q32_n64   1 0 64"
  "m5_q32_n32   1 1 32"
  "m6_q32_n32   1 0 32"
  "m5_q32_n16   1 1 16"
  "m6_q32_n16   1 0 16"
  "m5_q32_n8    1 1 8"
  "m6_q32_n8    1 0 8"
  "m3_q32_n32   0 1 32"
  "m3_q32_n16   0 1 16"
  "m3_q32_n8    0 1 8"
)

run_case() {
  read -r name nlp flat ncell <<< "$1"
  [ -f "$LOG/$name.log" ] && { echo "skip $name (log exists)"; return; }
  local slope=0
  [ "$flat" = 0 ] && slope=$BED_SLOPE
  (
    export BALFEM_REGIME=nonlinear BALFEM_NL_PRESSURE=$nlp BALFEM_FLAT_BED=$flat
    export BALFEM_BED_SLOPE=$slope BALFEM_BED_NW=1
    export BALFEM_FE_ORDER=3 BALFEM_P_ETA=2 BALFEM_NCELL=$ncell BALFEM_NLAMBDA=1
    export BALFEM_AWAVE=0.10 BALFEM_TWAVE=1.6 BALFEM_D=3.5 BALFEM_DT=0.04 BALFEM_PERIODS=$PERIODS
    export BALFEM_SOLVER=theta BALFEM_USE_AD=0 BALFEM_NOISE=1e-8
    export BALFEM_CHECKPOINT_EVERY=$CHECKPOINT_EVERY BALFEM_CHECKPOINT_KEEP=2
    export BALFEM_SAVE_EVERY=5 BALFEM_DIAG_EVERY=5 BALFEM_PRINT_EVERY=25
    export BALFEM_OUTDIR=$OUT/$name
    export OPENBLAS_NUM_THREADS=2
    echo "case $name | NL_PRESSURE=$nlp FLAT_BED=$flat BED_SLOPE=$slope NCELL=$ncell | Q3/Q2 | $PERIODS periods"
    echo "started $(date '+%F %T') | git $(git rev-parse --short HEAD)$(git diff --quiet -- src examples || echo ' + uncommitted edits')"
    systemd-run --user --scope --quiet --unit="stab1d-$name-$$" -p MemoryMax="$MEM_PER_CASE" \
      stdbuf -oL -eL julia "${JFLAGS[@]}" --project="$PROJ" examples/local_1d/run_periodic_1d.jl
    echo "--- exit status $? at $(date '+%F %T') ---"
  ) > "$LOG/$name.log" 2>&1
}

for c in "${CASES[@]}"; do
  while [ "$(pgrep -fc '^[^ ]*julia .*run_periodic_1d[.]jl')" -ge "$MAXP" ]; do sleep 30; done
  echo "$(date '+%F %T')  launch  $(echo $c | cut -d' ' -f1)" | tee -a "$LOG/_batch.log"
  run_case "$c" &
  sleep 60          # stagger: let each process take the precompile lock in turn
done
wait
echo "$(date '+%F %T')  batch finished" | tee -a "$LOG/_batch.log"
