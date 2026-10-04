#!/bin/bash
# ==============================================================
#  run_1dper_batch_ghost.sh — box campaign of the DIRECT (volume) GHOST PENALTY (2026-10-02)
#
#  markdown_files/GHOST_PENALTY_PLAN.md §4. Closed x-periodic box, P1LFE-2, Q3/Q2, CRANK–NICOLSON
#  (rule 15), 100 periods (160 s), BALFEM_STAB=ghostvolume with γ_u = γ_η = γ.
#  γ* comes from the eigen-analysis (plan §5, T5) and is passed as $1.
#
#  ARMS (plan §4):
#    G1 :full broken  γ*      G2 :native  γ*      G3 :full broken  γ*/3      G4 :full broken  3γ*
#    G5 :full broken  γ*, 32 cells/λ        G6 :full broken  γ*, A = 0.15        G7 :native  γ*, 32 cells/λ
#  Controls already on record (same branch code): :full broken γ=0 died at 44 s; :native γ=0 grows
#  from ≈ 90 s; :jumpgrad order 1 (0.3, 0.3): :full mid-band growth from ≈ 110 s.
#
#  ⚠ CORE BUDGET: at most MAXTOT simulations on the machine IN TOTAL (default 10), counting every
#  running periodic AND flume solver, not only this batch's.
#  ⚠ Logs to persistent disk (rule 47); never edit this file while it runs (rule 38h).
#  USAGE:  run/local/run_1dper_batch_ghost.sh <gamma*> [MAXTOT]
# ==============================================================
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
G="${1:?usage: $0 <gamma*> [MAXTOT]}"
MAXTOT="${2:-10}"
OUT=output/local_1d/periodic
LOG=$OUT/_logs
mkdir -p "$LOG"
G3=$(python3 -c "print(f'{$G/3:.3g}')"); GX3=$(python3 -c "print(f'{3*$G:.3g}')")

# name | FE_ORDER | NCELL | AWAVE | NL_PRESSURE | BROKEN | GAMMA
CASES=(
  "p32_broken_ghost${G}_n16_cn_p100       3 16 0.10 full   1 $G"
  "p32_native_ghost${G}_n16_cn_p100       3 16 0.10 native 0 $G"
  "p32_broken_ghost${G3}_n16_cn_p100      3 16 0.10 full   1 $G3"
  "p32_broken_ghost${GX3}_n16_cn_p100     3 16 0.10 full   1 $GX3"
  "p32_broken_ghost${G}_n32_cn_p100       3 32 0.10 full   1 $G"
  "p32_broken_ghost${G}_n16_A15_cn_p100   3 16 0.15 full   1 $G"
  "p32_native_ghost${G}_n32_cn_p100       3 32 0.10 native 0 $G"
)

run_case() {
  read -r name fe ncell A nlp broken g <<< "$1"
  [ -f "$LOG/$name.log" ] && { echo "skip $name (log exists)"; return; }
  (
    export BALFEM_REGIME=nonlinear BALFEM_FLAT_BED=1 BALFEM_NL_PRESSURE=$nlp
    export BALFEM_MIXED=0 BALFEM_BROKEN=$broken
    export BALFEM_STAB=ghostvolume BALFEM_CIP_GU=$g BALFEM_CIP_GE=$g BALFEM_CIP_ORDER=1
    [ "$nlp" = "full" ] && export BALFEM_C3_MASK=gs
    export BALFEM_FE_ORDER=$fe BALFEM_P_ETA=$((fe - 1)) BALFEM_NCELL=$ncell BALFEM_NLAMBDA=1
    export BALFEM_AWAVE=$A BALFEM_TWAVE=1.6 BALFEM_DT=0.04 BALFEM_PERIODS=100
    export BALFEM_SOLVER=theta BALFEM_USE_AD=0 BALFEM_NOISE=1e-8
    export BALFEM_SAVE_EVERY=5 BALFEM_DIAG_EVERY=5 BALFEM_PRINT_EVERY=25
    export BALFEM_OUTDIR=$OUT/$name
    source run/local/balfem_local.sh
    balfem_local_run examples/local_1d/run_periodic_1d.jl
  ) > "$LOG/$name.log" 2>&1
}

running() { pgrep -fc '^[^ ]*/julia .*run_(periodic|flume)_1d[.]jl'; }
for c in "${CASES[@]}"; do
  while [ "$(running)" -ge "$MAXTOT" ]; do sleep 60; done
  echo "$(date '+%F %T')  launch  $(echo $c | cut -d' ' -f1)   (running before launch: $(running))"
  run_case "$c" &
  sleep 60
done
wait
echo "$(date '+%F %T')  batch finished"
