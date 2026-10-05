#!/bin/bash
# ==============================================================
#  run_1dper_batch_broken.sh — SCREENING batch S1 of the broken formulation + C⁰-IP (2026-09-30)
#
#  THE QUESTIONS (markdown_files/BROKEN_FORMULATION_PLAN.md §T7; LaTeX §6.4):
#    (a) does the BROKEN (distributional) Class-III operator — cellwise Hessians + skeleton
#        layer, no projection, no auxiliary unknowns — change the interior grid-scale
#        instability that the mixed formulation showed under Crank–Nicolson?
#          prediction of the audit: NO (the layer terms are consistent, not sign-definite);
#    (b) does the C⁰ interior penalty J_h remove it, and at which γ_u?
#          estimate: σ_grid ≈ γ√(gd)/h ≈ 23γ s⁻¹ at h = 0.25 m against the measured +0.48 s⁻¹
#          ⇒ γ ≳ 0.02, carrier damping σ_c ≲ 1e-3 s⁻¹;
#    (c) does the same penalty cure the :native instability (Class-III-free, so not a
#        regularity problem at all)?
#
#  DESIGN. Closed x-periodic box, ONE model wavelength, P1LFE-2, A = 0.10, dt = 0.04,
#  CRANK–NICOLSON (rule 15 — no stability claim may rest on SDIRK_2_2), c3_mask = gs (the
#  mixed campaign's setting). 40 periods (64 s): the unstabilised controls on record diverge at
#  26 s (Q3/Q2 mixed) and ≈ 90 s (Q2/Q1 32 cells mixed). The γ=0 broken arms are the
#  same-batch controls for the CIP arms (rule 14c); the :native twin is in the batch too.
#
#  Analysis per run:
#    PG_NCELL=<n> PG_PU=<p> julia --project=postprocessing \
#        postprocessing/examples/periodic_growth.jl output/local_1d/periodic/<case>
#
#  ⚠ Logs go to persistent disk (rule 47); never edit this file while it runs (rule 38h).
#  USAGE:  run/local/run_1dper_batch_broken.sh [MAX_PARALLEL]   (default 4)
# ==============================================================
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
MAXP="${1:-4}"
OUT=output/local_1d/periodic
LOG=$OUT/_logs
mkdir -p "$LOG"

# case : name | FE_ORDER | NCELL | AWAVE | NL_PRESSURE | BROKEN | CIP_GU | CIP_GE | PERIODS
#  ⚠ REVISED BEFORE LAUNCH (2026-09-30), from the LINEAR EIGEN-ANALYSIS of the box operator:
#  the C⁰-IP penalty does NOT damp oscillatory modes wavenumber-selectively — it OVERDAMPS the
#  non-C¹ part of the space (the penalty limit is a C¹ constraint) and leaves the oscillatory
#  modes that fit in the C¹ subspace almost undamped (≤ 0.016 s⁻¹ for γ_u, falling like 1/γ).
#  The γ_η penalty removes the upper half of the η band (above the CELL Nyquist). So the arms
#  test the two regimes: TRANSITION (γ = 0.03) and CONSTRAINT (γ = 0.3), BOTH fields.
CASES=(
  "p32_broken_n16_cn                3 16 0.10 full   1 0    0    40"
  "p32_broken_cipu0.3e0.3_n16_cn    3 16 0.10 full   1 0.3  0.3  40"
  "p32_native_cipu0.3e0.3_n16_cn    3 16 0.10 native 0 0.3  0.3  40"
  "p32_native_n16_cn_s1             3 16 0.10 native 0 0    0    40"
  "p32_broken_cipu0.03e0.03_n16_cn  3 16 0.10 full   1 0.03 0.03 40"
  "p21_broken_n32_cn                2 32 0.10 full   1 0    0    40"
  "p21_broken_cipu0.3e0.3_n32_cn    2 32 0.10 full   1 0.3  0.3  40"
)

run_case() {
  read -r name fe ncell A nlp broken gu ge periods <<< "$1"
  [ -f "$LOG/$name.log" ] && { echo "skip $name (log exists)"; return; }
  (
    export BALFEM_REGIME=nonlinear BALFEM_FLAT_BED=1 BALFEM_NL_PRESSURE=$nlp
    export BALFEM_MIXED=0 BALFEM_BROKEN=$broken BALFEM_CIP_GU=$gu BALFEM_CIP_GE=$ge
    [ "$nlp" = "full" ] && export BALFEM_C3_MASK=gs
    export BALFEM_FE_ORDER=$fe BALFEM_P_ETA=$((fe - 1)) BALFEM_NCELL=$ncell BALFEM_NLAMBDA=1
    export BALFEM_AWAVE=$A BALFEM_TWAVE=1.6 BALFEM_DT=0.04 BALFEM_PERIODS=$periods
    export BALFEM_SOLVER=theta BALFEM_USE_AD=0 BALFEM_NOISE=1e-8
    export BALFEM_SAVE_EVERY=5 BALFEM_DIAG_EVERY=5 BALFEM_PRINT_EVERY=25
    export BALFEM_OUTDIR=$OUT/$name
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
