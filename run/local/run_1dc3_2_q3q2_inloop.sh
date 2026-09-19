#!/bin/bash
# ==============================================================
#  CLASS-III TREATMENT FACTORIAL — run #2: Q3/Q2 IN-LOOP — isolates the LAG. THE MEASUREMENT THAT MATTERS
#
#  Branch `new-classIII-treatment`. Design: markdown_files/NEW_TREATMENT.md Part C.4.
#
#  THE QUESTION. `nl_pressure=:full` dies on a FLAT bed at t = 12.60 s at Q2/Q1 with
#  the projections frozen one step behind. Two things changed on this branch, and this
#  four-run factorial separates them from the pairing:
#      #0  Q2/Q1  lagged   <- CONTROL, must reproduce 12.60 s
#      #1  Q3/Q2  lagged   <- the PAIRING alone
#      #2  Q3/Q2  in-loop  <- the LAG removed (static condensation)
#      #3  Q3/Q2  in-loop  <- repeat of #2
#
#  ⚠ RUN #0 FIRST AND CHECK IT GIVES 12.60 s. A factorial whose control does not
#  reproduce the known number is measuring something else (rule 38c, and rule 12b's
#  root-cause-of-the-root-cause: the failing configuration was never the verified one).
#
#  ⚠ 100 PERIODS, NOT 16. The flume fills in 36 s at c_g = 1.25 m/s, so NOTHING before
#  t ~ 36 s may be read as growth (rule 14 — it has caught four measurements). And a run
#  that ends early proves nothing: the 1:2 bar looked stable at 80 s and died at 137 s
#  (rule 12c). Survival is only meaningful at >= 100 periods.
#
#  ⚠ The algebraic reduction is EXACT (test_class3_reduction.jl 19/19,
#  test_class3_residual_parity.jl 4/4 at 3e-16), so it cannot be what moves any onset
#  here. It is present in every arm, including the control.
# ==============================================================
source "$(dirname "${BASH_SOURCE[0]}")/balfem_local.sh"

# ---- physics: IDENTICAL in all four runs -----------------------------------
export BALFEM_WAVE_GEN=bc
export BALFEM_REGIME=nonlinear
export BALFEM_NL_PRESSURE=full
export BALFEM_FLAT_BED=1
export BALFEM_RELAX=1
export BALFEM_MPI=0
export BALFEM_LX=60.0
export BALFEM_NY=1
export BALFEM_YBC=wall
export BALFEM_AWAVE=0.1
export BALFEM_TWAVE=1.6
export BALFEM_SOLVER=sdirk
export BALFEM_TABLEAU=SDIRK_2_2
export BALFEM_NX=240            # dx = 0.25, as the 12.60 s baseline
export BALFEM_LY=0.25           # Ly = dx keeps the cross-cell isotropic (rule 12/22)
export BALFEM_DT=0.04
export BALFEM_USE_AD=0          # hand Jacobian — rule 17b: it does not set the onset
export BALFEM_PERIODS=100

# ---- the two varied knobs --------------------------------------------------
export BALFEM_FE_ORDER=3
export BALFEM_P_ETA=2        # Taylor-Hood (rule 2b)
export BALFEM_NLP_INLOOP=1  # 1 = projections from the CURRENT Newton iterate

# common output time grid (0.2 s) so onsets compare on the same sampling
export BALFEM_SAVE_EVERY=5
export BALFEM_DIAG_EVERY=5
export BALFEM_PRINT_EVERY=25

export BALFEM_OUTDIR=output/local_1d/c3_2_Q3Q2_inloop

balfem_local_run examples/local_1d/run_flume_1d.jl
