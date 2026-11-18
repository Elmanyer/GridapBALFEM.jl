#!/bin/bash
# ==============================================================
#  run_mixed_box_controls.sh — v3-mixed WITHOUT stabilisation in the closed 1-D box (branch mixed-vms-solver)
#
#  THE QUESTION. Is the mixed formulation (all eight 𝓝 components, 7 fields, component 4 ON, no CIP,
#  no VMS) unstable in the closed box, and does the onset move EARLIER with refinement, as v1's mixed
#  form did with component 4 off (Q2/Q1: 157 / 104 / 90 / 31 s at 8 / 16 / 32 / 64 cells/λ; Q3/Q2 n16:
#  26 s; archive/HISTORY_V1.md §5.2f)? These are also the unstabilised controls of every later VMS
#  run (rule 14c).
#
#  DESIGN (rule 38c): the environment of run/local/run_1d_stability_campaign.sh copied verbatim
#  (model 5, P1LFE-2, Q3/Q2, A = 0.10, T = 1.6 s, D = 3.5 m, dt = 0.04, CRANK–NICOLSON, seed 1e-8,
#  one wavelength, 200 periods, diagnostics/VTK every 5 steps, checkpoints every 250 steps); the
#  only change is BALFEM_FORMULATION=mixed (which drops the CIP). Reference arms: output/
#  1D_stability_campaign/m5_q32_n16, m5_q32_n32 (v2 primal + CIP; n32 burst at ≈ 245 s).
#
#  ARMS: m5 mixed at 16 and 32 cells/λ (the refinement pair).
#  ANALYSIS: PG_NCELL=<n> PG_PU=3 julia --project=postprocessing \
#                postprocessing/examples/periodic_growth.jl output/mixed_vms/box_controls/<case>
#
#  Memory-gated by _slot_throttle.sh (memory is the only cap, rule 47). A case whose log exists is
#  skipped. Never edit this file while it runs (rule 38h). Resume a lost case with
#  BALFEM_RESTART_FROM=<its dir>.
#  Detach: setsid nohup run/local/run_mixed_box_controls.sh > output/mixed_vms/box_controls/_logs/_batch.out 2>&1 < /dev/null &
# ==============================================================
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
PROJ="$(pwd)"
OUT=output/mixed_vms/box_controls
LOG=$OUT/_logs
mkdir -p "$LOG"
MEM_PER_CASE="${MEM_PER_CASE:-9G}"
PERIODS="${PERIODS:-200}"
export RUN_GB="${RUN_GB:-5}"
export SLOT_STAGGER="${SLOT_STAGGER:-60}"
source "$PROJ/run/local/_slot_throttle.sh"

#  Count every Julia driver of this branch as well as the shared ones.
PAT='^[^ ]*julia .*(run_periodic_1d|run_flume_1d|sponge_worker|compare_mixed_vms|vms_newton_probe)[.]jl'
slots_busy() { pgrep -fc "$PAT"; }
mem_free_gb() {
  local avail owed=0 rss age pid peak=$(( RUN_GB * 1048576 ))
  avail=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
  for pid in $(pgrep -f "$PAT"); do
    age=$(ps -o etimes= -p "$pid" 2>/dev/null) || continue
    [ -n "$age" ] && [ "$age" -lt "$YOUNG_S" ] || continue
    rss=$(awk '/^VmRSS:/ {print $2}' "/proc/$pid/status" 2>/dev/null) || continue
    [ -n "$rss" ] && [ "$rss" -lt "$peak" ] && owed=$(( owed + peak - rss ))
  done
  echo $(( (avail - owed) / 1048576 ))
}

# ---- source snapshot: the runs are launched from an uncommitted tree ------------------------
SNAP=$LOG/source_$(date '+%Y%m%d_%H%M%S')
mkdir -p "$SNAP"
git rev-parse HEAD > "$SNAP/commit.txt"
git diff HEAD > "$SNAP/tracked.diff"
git ls-files --others --exclude-standard -z | tar czf "$SNAP/untracked.tgz" --null -T -
echo "$(date '+%F %T')  source snapshot: $SNAP" | tee -a "$LOG/_batch.log"

# case : name | NCELL      (longest first)
CASES=(
  "m5_mixed_q32_n32  32"
  "m5_mixed_q32_n16  16"
)

run_case() {
  read -r name ncell <<< "$1"
  (
    export BALFEM_REGIME=nonlinear BALFEM_NL_PRESSURE=1 BALFEM_FLAT_BED=1
    export BALFEM_BED_SLOPE=0 BALFEM_BED_NW=1
    export BALFEM_FORMULATION=mixed BALFEM_VMS=0
    export BALFEM_FE_ORDER=3 BALFEM_P_ETA=2 BALFEM_NCELL=$ncell BALFEM_NLAMBDA=1
    export BALFEM_AWAVE=0.10 BALFEM_TWAVE=1.6 BALFEM_D=3.5 BALFEM_DT=0.04 BALFEM_PERIODS=$PERIODS
    export BALFEM_SOLVER=theta BALFEM_USE_AD=0 BALFEM_NOISE=1e-8
    export BALFEM_CHECKPOINT_EVERY=250 BALFEM_CHECKPOINT_KEEP=2
    export BALFEM_SAVE_EVERY=5 BALFEM_DIAG_EVERY=5 BALFEM_PRINT_EVERY=25
    export BALFEM_OUTDIR=$OUT/$name
    export OPENBLAS_NUM_THREADS=2
    echo "case $name | m5 MIXED (no CIP, no VMS) NCELL=$ncell | Q3/Q2 | A=0.10 | CN dt=0.04 | $PERIODS periods"
    echo "started $(date '+%F %T') | git $(git rev-parse --short HEAD)$(git diff --quiet -- src examples || echo ' + uncommitted edits')"
    systemd-run --user --scope --quiet --unit="mixbox-$name-$$" -p MemoryMax="$MEM_PER_CASE" \
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
echo "$(date '+%F %T')  batch finished" | tee -a "$LOG/_batch.log"
