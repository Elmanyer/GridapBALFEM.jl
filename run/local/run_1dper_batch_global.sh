#!/bin/bash
# ==============================================================
#  run_1dper_batch.sh — the x-PERIODIC interior-stability campaign (2026-09-24)
#
#  THE QUESTION. Every :full failure on the 60 m flume pins at the inflow or the relaxation
#  zone. On a periodic box of ONE model wavelength there is no inflow, no relaxation, no
#  sponge and no source, and the continuum answer is known: nothing grows (no Benjamin–Feir
#  sidebands fit, and κa ≈ 0.16 is superharmonically stable). So:
#     :full grows here          ⇒ INTERIOR discretisation instability (→ grid-scale damping)
#     :full stays bounded here  ⇒ the flume failure is BOUNDARY-driven (→ fix the inflow/relaxation)
#  and the dx ladder separates a grid-scale mode (σ rising as dx falls) from a physical one.
#  LaTeX: chapter "Stability Analysis", §7.7 (what is not yet established).
#
#  DESIGN (rules 14c, 38c): every :full arm has its :native control at the SAME settings in
#  the SAME batch; the mixed and projected treatments of ∇𝖲 are both run (C3_MASK=gs, the
#  arm that carries the flume failure). All: P1LFE-2, flat, A=0.10 unless stated, T=1.6,
#  SDIRK_2_2, dt=0.04, 100 periods, seed 1e-8 on every harmonic.
#
#  Analysis per run:
#    PG_NCELL=<n> PG_PU=<p> julia --project=postprocessing \
#        postprocessing/examples/periodic_growth.jl output/local_1d/periodic/<case>
#
#  ⚠ Logs go to persistent disk (rule 47); never edit this file while it runs (rule 38h).
#  USAGE:  run/local/run_1dper_batch_global.sh [MAX_PARALLEL]  (default 4)
#  ⚠ Successor of run_1dper_batch.sh (2026-09-24): the slot limit counts ALL running
#  periodic solvers (not just this shell's children), and a case is skipped once its
#  LOG exists — so it can take over a batch that is already running without duplicating it.
# ==============================================================
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
MAXP="${1:-4}"
OUT=output/local_1d/periodic
LOG=$OUT/_logs
mkdir -p "$LOG"

# case : name | FE_ORDER | NCELL | AWAVE | NL_PRESSURE | MIXED | extra env
CASES=(
  # ---- core triad, Q2/Q1, dx = λ/16 ≈ 0.25 m (the flume's resolution) --------------
  "p21_native_n16     2 16 0.10 native 0"
  "p21_mixed_n16      2 16 0.10 full   1 BALFEM_C3_MASK=gs"
  "p21_proj_n16       2 16 0.10 full   0 BALFEM_C3_MASK=gs"
  # ---- dx ladder, Q2/Q1 ----------------------------------------------------------------
  "p21_native_n32     2 32 0.10 native 0"
  "p21_mixed_n32      2 32 0.10 full   1 BALFEM_C3_MASK=gs"
  "p21_mixed_n8       2 8  0.10 full   1 BALFEM_C3_MASK=gs"
  "p21_native_n8      2 8  0.10 native 0"
  # ---- pairing: Q3/Q2 ------------------------------------------------------------------
  "p32_native_n16     3 16 0.10 native 0"
  "p32_mixed_n16      3 16 0.10 full   1 BALFEM_C3_MASK=gs"
  # ---- amplitude -----------------------------------------------------------------------
  "p21_mixed_n16_A15  2 16 0.15 full   1 BALFEM_C3_MASK=gs"
  "p21_native_n16_A15 2 16 0.15 native 0"
  # ---- finest Q2/Q1 level last ---------------------------------------------------------
  "p21_mixed_n64      2 64 0.10 full   1 BALFEM_C3_MASK=gs"
)

run_case() {
  read -r name fe ncell A nlp mixed extra <<< "$1"
  [ -f "$LOG/$name.log" ] && { echo "skip $name (log exists)"; return; }
  (
    export BALFEM_REGIME=nonlinear BALFEM_FLAT_BED=1 BALFEM_NL_PRESSURE=$nlp BALFEM_MIXED=$mixed
    export BALFEM_FE_ORDER=$fe BALFEM_P_ETA=$((fe - 1)) BALFEM_NCELL=$ncell BALFEM_NLAMBDA=1
    export BALFEM_AWAVE=$A BALFEM_TWAVE=1.6 BALFEM_DT=0.04 BALFEM_PERIODS=100
    export BALFEM_SOLVER=sdirk BALFEM_TABLEAU=SDIRK_2_2 BALFEM_USE_AD=0 BALFEM_NOISE=1e-8
    export BALFEM_SAVE_EVERY=5 BALFEM_DIAG_EVERY=5 BALFEM_PRINT_EVERY=25
    export BALFEM_OUTDIR=$OUT/$name
    [ -n "$extra" ] && export $extra
    source run/local/balfem_local.sh
    balfem_local_run examples/local_1d/run_periodic_1d.jl
  ) > "$LOG/$name.log" 2>&1
}

for c in "${CASES[@]}"; do
  #  anchored on the executable (argv[0] ends in /julia): an unanchored pattern also matched
  #  any shell whose command line merely MENTIONED it, and stalled the queue (2026-09-24)
  while [ "$(pgrep -fc '^[^ ]*/julia .*run_periodic_1d[.]jl')" -ge "$MAXP" ]; do sleep 30; done
  echo "$(date '+%F %T')  launch  $(echo $c | cut -d' ' -f1)"
  run_case "$c" &
  sleep 60          # stagger: let each process take the precompile lock in turn
done
wait
echo "$(date '+%F %T')  batch finished"
