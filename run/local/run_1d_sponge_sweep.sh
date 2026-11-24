#!/bin/bash
# ==============================================================
#  run_1d_sponge_sweep.sh — sponge-layer reflection sweep (CH8_CAMPAIGN_TODO.md §7)
#
#  THE QUESTION. The Snellius linear suite reflects 0.5 / 5 / 9 / 18 % of the amplitude from its
#  12 m, μ_max = 40 outflow sponge at T = 1.5 / 2.0 / 2.2 / 2.9 s (LINEAR_PRODUCTION_RESULTS.md).
#  Which (width W, strength μ_max) of the quadratic sponge brings Kr below 1 % over the Tp = 2 band,
#  at the smallest W?
#
#  STAGE A (84 cases): μ_max ∈ {4, 8, 16, 25, 40} × W ∈ {6, 12, 18, 24} m × T ∈ {1.5, 2.0, 2.5, 2.9} s,
#    plus the controls in the same batch (rule 38c):
#      floor    W = 48 m, μ_max = 4 at every T: what the relaxation zone and the generation
#               reflect on their own;
#      replica  W = 12, μ_max = 40 (inside the matrix): must reproduce the Snellius Kr(T), or this
#               1-D set-up does not represent the 2-D runs.
#  STAGE B (4 cases, after A): relaxation width 12 m at the best (W, μ_max) of stage A.
#    STAGE=B BEST_W=<W> BEST_MU=<μ> run/local/run_1d_sponge_sweep.sh
#
#  CASES: run/local/sponge_cases.py writes the table (geometry relax 6 m | interior 32 m | sponge W,
#  dx = 0.25, duration from the c_g transit, longest first). Settings: run/local/sponge_worker.jl
#  (the Snellius suite in 1-D: P1LFE-2, Q2/Q1, CN, dt = 0.02, linear, flat, BC generation + relaxation).
#  Code: the MAIN tree (models 1–2 are unaffected by every change since the Snellius runs).
#
#  EXECUTION. Workers, not one process per case: each worker (sponge_worker.jl) compiles once and
#  runs up to CASES_PER_WORKER cases pulled from the queue through claim directories (rules 41–44),
#  in its own systemd scope capped at MEM_PER_WORKER. Slots are shared with every other local 1-D
#  launcher through run/local/_slot_throttle.sh (MAXP processes, MemAvailable ≥ MIN_FREE_GB).
#  Estimated cost: 8.9e7 cell-steps in total; measured in the smoke run (2026-11-06) ≈ 3.3e-3 s
#  per cell-step + ≈ 6.5 min JIT per worker ⇒ ≈ 80 CPU-h, ≈ 14 h on 6 slots.
#  Measured in the real run (2026-11-07): 26 cases took ≈ 8 h on 2 workers. (The thermal throttle
#  above 3 busy cores seen then was a machine fault, cleared by a reboot; memory is the only slot
#  cap since 2026-11-10.) State, results and how to resume: markdown_files/SPONGE_SWEEP_RESULTS.md.
#
#  ANALYSIS:  julia --project=postprocessing postprocessing/examples/sponge_reflection.jl \
#                 output/1D_sponge_sweep [A|B]
#
#  ⚠ Logs on persistent disk (rule 47); never edit this file, sponge_worker.jl or sponge_cases.py
#  while the sweep runs (rule 38h). A case whose claim directory or log exists is skipped; to re-run
#  it delete output/1D_sponge_sweep/_logs/_claims/<case>/, its log and its output directory.
#
#  USAGE:  run/local/run_1d_sponge_sweep.sh                    (STAGE=A, MIN_FREE_GB=6, no process cap)
#  Detach: nohup run/local/run_1d_sponge_sweep.sh > output/1D_sponge_sweep/_logs/_batch_A.out 2>&1 &
#  Smoke:  SMOKE=1 run/local/run_1d_sponge_sweep.sh   (2 cases × 1 s into output/1D_sponge_sweep/_smoke)
# ==============================================================
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
PROJ="$(pwd)"
STAGE="${STAGE:-A}"
CASES_PER_WORKER="${CASES_PER_WORKER:-12}"
MEM_PER_WORKER="${MEM_PER_WORKER:-4G}"
OUT=output/1D_sponge_sweep
[ "${SMOKE:-0}" = 1 ] && OUT=$OUT/_smoke
LOG=$OUT/_logs
mkdir -p "$LOG"
source "$PROJ/run/local/_slot_throttle.sh"

# ---- the case table (kept with the logs: the analysis reads it) ------------------------------
TABLE=$PROJ/$LOG/_cases_$STAGE.tsv
if [ ! -f "$TABLE" ]; then
  case "$STAGE" in
    A) python3 -I run/local/sponge_cases.py A > "$TABLE" ;;
    B) [ -n "$BEST_W" ] && [ -n "$BEST_MU" ] || { echo "STAGE=B needs BEST_W and BEST_MU"; exit 1; }
       python3 -I run/local/sponge_cases.py B "$BEST_W" "$BEST_MU" > "$TABLE" ;;
    *) echo "STAGE must be A or B"; exit 1 ;;
  esac
  if [ "${SMOKE:-0}" = 1 ]; then          # the cheapest and the most expensive case, 1 s each
    { head -1 "$TABLE"; sed -n '2p' "$TABLE"; tail -1 "$TABLE"; } \
      | awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print; next} {$8 = 1.0; print}' > "$TABLE.tmp"
    mv "$TABLE.tmp" "$TABLE"
  fi
fi
NCASES=$(( $(wc -l < "$TABLE") - 1 ))

# ---- source snapshot -------------------------------------------------------------------------
SNAP=$LOG/source_$(date '+%Y%m%d_%H%M%S')
mkdir -p "$SNAP"
git rev-parse HEAD > "$SNAP/commit.txt"
git diff HEAD > "$SNAP/tracked.diff"
git ls-files --others --exclude-standard -z | tar czf "$SNAP/untracked.tgz" --null -T -
echo "$(date '+%F %T')  stage $STAGE: $NCASES cases ($TABLE); source $SNAP (commit $(cut -c1-10 "$SNAP/commit.txt") + diff)" | tee -a "$LOG/_batch.log"

# ---- sysimage: only if built from exactly this src/ (rule 25) --------------------------------
IMAGE="$PROJ/GridapBALFEM_sysimage_local.so"
src_hash() { ( cd "$PROJ/src" && find . -name '*.jl' -type f -print0 \
               | LC_ALL=C sort -z | xargs -0 -r sha256sum ) | sha256sum | awk '{print $1}'; }
JFLAGS=()
if [ -f "$IMAGE" ] && [ -f "$IMAGE.src.sha256" ] && [ "$(cat "$IMAGE.src.sha256")" = "$(src_hash)" ]; then
  JFLAGS=(-J "$IMAGE")
fi

claimed()   { find "$PROJ/$LOG/_claims" -mindepth 1 -maxdepth 1 -type d -name "sp${STAGE}_*" 2>/dev/null | wc -l; }   # this stage only
unclaimed() { echo $(( NCASES - $(claimed) )); }

start_worker() {
  local id=$1
  OPENBLAS_NUM_THREADS=2 systemd-run --user --scope --quiet --unit="sponge1d-$STAGE-$id-$$" \
    -p MemoryMax="$MEM_PER_WORKER" \
    stdbuf -oL -eL julia "${JFLAGS[@]}" --project="$PROJ" run/local/sponge_worker.jl \
      "$TABLE" "$PROJ/$OUT" "$PROJ/$LOG" "$CASES_PER_WORKER" >> "$LOG/_workers.out" 2>&1
}

w=0
while [ "$(unclaimed)" -gt 0 ]; do
  with_slot start_worker $((++w))
  echo "$(date '+%F %T')  worker $w started ($(unclaimed) of $NCASES unclaimed; $(slots_busy) busy, $(mem_free_gb) GB free)" | tee -a "$LOG/_batch.log"
done
wait
echo "$(date '+%F %T')  stage $STAGE: every case claimed and every worker finished" | tee -a "$LOG/_batch.log"
grep -L 'flume_1d \[.*\] done:' "$LOG"/sp*.log 2>/dev/null | sed 's/^/  NOT FINISHED: /' | tee -a "$LOG/_batch.log"
