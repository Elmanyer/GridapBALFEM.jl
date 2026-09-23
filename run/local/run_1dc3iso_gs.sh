#!/bin/bash
# ==============================================================
#  CLASS-III TREATMENT FACTORIAL — CLASS-III ISOLATION, arm 'gs': the GRAD-S FAMILY ONLY (components {1,2,5} collapsed + N2's admissible remainder)
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
export BALFEM_PERIODS=25

# ---- the two varied knobs --------------------------------------------------
export BALFEM_FE_ORDER=2
export BALFEM_P_ETA=1        # Taylor-Hood (rule 2b)
export BALFEM_NLP_INLOOP=0  # 1 = projections from the CURRENT Newton iterate

# common output time grid (0.2 s) so onsets compare on the same sampling
export BALFEM_SAVE_EVERY=5
export BALFEM_DIAG_EVERY=5
export BALFEM_PRINT_EVERY=25

export BALFEM_OUTDIR=output/local_1d/c3iso_gs

# ---- THE ISOLATION KNOB ----------------------------------------------------
#  ⚠ THIS MUST SIT BEFORE `balfem_local_run` (rule 38h). The first version of these
#  launchers appended it AFTER the run line, where it is a silent NO-OP: all three arms
#  then produced the CONTROL's trajectory, bit-identical to every printed digit
#  (t = 12.600000, eta = 1.25518087e-01, u_max = 2.50462190, Newton 51). That
#  digit-for-digit agreement IS the tell rule 38h describes.
#
#  Measured statically on the assembled residual (test_class3_split.jl, flat bed, Q3/Q2):
#      |grad-S arm| = 1.2837e-02   |grad-b arm| = 1.1271e-05   ratio ~1140x
#      partition error |gs + gb - both| / |both| = 8.65e-15  (the split is exact)
#  So Class III is NUMERICALLY almost entirely the grad-S contraction. This batch asks
#  whether it is also DYNAMICALLY the carrier of the instability.
#  CONTROL: c3_0_Q2Q1_lagged died at t = 12.600000 with u_max = 2.5046.
export BALFEM_C3_MASK=gs

balfem_local_run examples/local_1d/run_flume_1d.jl
