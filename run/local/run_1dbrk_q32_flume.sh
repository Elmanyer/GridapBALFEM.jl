#!/bin/bash
# ==============================================================
#  run_1dbrk_q32_flume.sh — the broken formulation + C⁰-IP on the REAL 60 m flume (2026-10-01)
#
#  markdown_files/BROKEN_FORMULATION_PLAN.md §4.11. The closed periodic box (§4.5–4.10) showed
#  Q3/Q2 broken + C⁰-IP (γ_u = γ_η = 0.3) bounded under Crank–Nicolson; this is the flume test.
#  The flume adds the generation boundary, the relaxation zone and the sponge. Every earlier
#  :full failure pinned at the inflow, and the boundary ACCELERATES the interior instability
#  (CLAUDE.md §5.2f).
#
#  ENVIRONMENT copied VERBATIM from run_1dc3q_base_mixed.sh (rule 38c) — Q3/Q2, nx = 240
#  (dx = 0.25), dt = 0.04, A = 0.10, T = 1.6, relaxation inflow, flat bed, ny = 1 with walls,
#  c3_mask = gs — except for the three deliberate changes:
#    * mixed → broken (BALFEM_MIXED=0, BALFEM_BROKEN=1);
#    * the C⁰-IP penalty (BALFEM_CIP_GU/GE);
#    * SDIRK_2_2 → CRANK–NICOLSON (rule 15: SDIRK_2_2 masks this instability).
#  CONTROL in the same batch (rule 14c): broken, γ = 0, CN.
#  Records: Q3/Q2 projected died at 7.0 s, mixed at 16.0 s (both SDIRK).
#  ⚠ Rule 14: the flume fills in ≈ 36 s; nothing before that may be read as growth.
#
#  ⚠ Logs to persistent disk (rule 47); never edit this file while it runs (rule 38h).
#  USAGE: run/local/run_1dbrk_q32_flume.sh
# ==============================================================
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
LOG=output/local_1d/_logs
mkdir -p "$LOG"

# name | CIP_GU | CIP_GE
CASES=(
  "brk_q32_cipu0.3e0.3_cn  0.3 0.3"
  "brk_q32_cn              0   0"
)

run_case() {
  read -r name gu ge <<< "$1"
  [ -f "$LOG/$name.log" ] && { echo "skip $name (log exists)"; return; }
  (
    export BALFEM_WAVE_GEN=bc BALFEM_REGIME=nonlinear BALFEM_NL_PRESSURE=full BALFEM_FLAT_BED=1
    export BALFEM_RELAX=1 BALFEM_MPI=0 BALFEM_LX=60.0 BALFEM_NY=1 BALFEM_YBC=wall
    export BALFEM_AWAVE=0.1 BALFEM_TWAVE=1.6
    export BALFEM_SOLVER=theta
    export BALFEM_NX=240 BALFEM_LY=0.25 BALFEM_DT=0.04 BALFEM_USE_AD=0 BALFEM_PERIODS=100
    export BALFEM_FE_ORDER=3 BALFEM_P_ETA=2 BALFEM_NLP_INLOOP=0
    export BALFEM_SAVE_EVERY=5 BALFEM_DIAG_EVERY=5 BALFEM_PRINT_EVERY=25
    export BALFEM_C3_MASK=gs
    export BALFEM_MIXED=0 BALFEM_BROKEN=1 BALFEM_CIP_GU=$gu BALFEM_CIP_GE=$ge
    export BALFEM_OUTDIR=output/local_1d/$name
    source run/local/balfem_local.sh
    balfem_local_run examples/local_1d/run_flume_1d.jl
  ) > "$LOG/$name.log" 2>&1
}

for c in "${CASES[@]}"; do
  echo "$(date '+%F %T')  launch  $(echo $c | cut -d' ' -f1)"
  run_case "$c" &
  sleep 60
done
wait
echo "$(date '+%F %T')  batch finished"
