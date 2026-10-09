#!/bin/bash
# ==============================================================
#  run_1dper_v2_unstab_extra.sh — the A = 0.15 and 8-cell arms of run_1dper_v2_unstab.sh, launched at once (2026-10-06). The closed-box ladder on the v2 solver (2026-10-06)
#  (V2_SOLVER_PLAN.md §4 item 1: the v2 evidence for the stability chapter)
#
#  THE QUESTION. v1 showed, in the closed x-periodic box under Crank–Nicolson, that the
#  full-pressure discretisation grows grid-scale modes and diverges (mixed Q2/Q1: 157 / 104 /
#  90 / 31 s at 8 / 16 / 32 / 64 cells/λ; Q3/Q2 16 cells: 26 s; A = 0.15: 41 s; broken Q3/Q2
#  γ = 0: Newton failure at 44 s). Does the v2 solver reproduce it?
#
#  ⚠ NOT A ONE-VARIABLE REPEAT. Two things differ from every v1 arm, by design of v2:
#    * nl_pressure=1 carries ALL EIGHT 𝓝 components; every v1 :full box arm ran C3_MASK=gs,
#      i.e. without component 4 (and Class III by the mixed, projected or broken treatment);
#    * nl_pressure=0 carries NO 𝓝 component; the v1 twin was :native = {3,6,7,8}.
#  So "the same result" means the same qualitative signature (growth of the mid/high bands,
#  onset advancing with refinement, Q3/Q2 and A = 0.15 worse), not the same onset times.
#
#  DESIGN (rules 14c, 15, 38c, 38i): identical to v1 run_1dper_batch_cn{,2,3}.sh — P1LFE-2,
#  A = 0.10, T = 1.6 s, d = 3.5 m, dt = 0.04, CRANK–NICOLSON, seed 1e-8 on every harmonic,
#  100 periods, VTK/diagnostics every 5 steps, no stabiliser. Every nl_pressure=1 arm has its
#  nl_pressure=0 twin in the same batch.
#
#  Analysis per run:
#    PG_NCELL=<n> PG_PU=<p> julia --project=postprocessing \
#        postprocessing/examples/periodic_growth.jl output/local_1d/periodic/<case>
#
#  ⚠ Logs go to persistent disk (rule 47); never edit this file while it runs (rule 38h).
#  USAGE:  run/local/run_1dper_v2_unstab.sh [MAX_PARALLEL]  (default 5)
# ==============================================================
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
MAXP="${1:-9}"
OUT=output/local_1d/periodic
LOG=$OUT/_logs
mkdir -p "$LOG"

# case : name | FE_ORDER | NCELL | AWAVE | NL_PRESSURE      (v1 CN record in the comment)
#  Order: the discriminating cells first.
CASES=(
  "v2_p21_nlp1_n16_A15_cn  2 16 0.15 1"   # mixed 41 s
  "v2_p21_nlp0_n16_A15_cn  2 16 0.15 0"   # :native flat
  "v2_p21_nlp1_n8_cn       2 8  0.10 1"   # mixed 157 s
  "v2_p21_nlp0_n8_cn       2 8  0.10 0"   # :native flat
)

run_case() {
  read -r name fe ncell A nlp <<< "$1"
  [ -f "$LOG/$name.log" ] && { echo "skip $name (log exists)"; return; }
  (
    export BALFEM_REGIME=nonlinear BALFEM_FLAT_BED=1 BALFEM_NL_PRESSURE=$nlp
    export BALFEM_FE_ORDER=$fe BALFEM_P_ETA=$((fe - 1)) BALFEM_NCELL=$ncell BALFEM_NLAMBDA=1
    export BALFEM_AWAVE=$A BALFEM_TWAVE=1.6 BALFEM_DT=0.04 BALFEM_PERIODS=100
    export BALFEM_SOLVER=theta BALFEM_USE_AD=0 BALFEM_NOISE=1e-8
    export BALFEM_SAVE_EVERY=5 BALFEM_DIAG_EVERY=5 BALFEM_PRINT_EVERY=25
    export BALFEM_OUTDIR=$OUT/$name
    export OPENBLAS_NUM_THREADS=2
    source run/local/balfem_local.sh
    balfem_local_run examples/local_1d/run_periodic_1d.jl
  ) > "$LOG/$name.log" 2>&1
}

for c in "${CASES[@]}"; do
  while [ "$(pgrep -fc '^[^ ]*/julia .*run_periodic_1d[.]jl')" -ge "$MAXP" ]; do sleep 30; done
  echo "$(date '+%F %T')  launch  $(echo $c | cut -d' ' -f1)"
  run_case "$c" &
  sleep 60          # stagger: let each process take the precompile lock in turn
done
wait
echo "$(date '+%F %T')  batch finished"
