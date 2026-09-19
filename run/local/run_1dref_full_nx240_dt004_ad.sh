#!/bin/bash
# ==============================================================
#  :full FLAT — dx x dt REFINEMENT, EXACT AD JACOBIAN   [ad]
#
#  THE QUESTION. `nl_pressure=:full` on a FLAT bed dies at t = 12.6 s (8 wave
#  periods) with a VELOCITY-LED mode at the generation boundary: u_max goes
#  0.45 -> 2.50 in ~1 s while eta stays bounded at 0.12, against a :native
#  control that holds u_max = 0.42 for 100 periods (OPEN_ISSUES.md 0c).
#
#  RULE 38b IS THE DISCRIMINATOR, AND IT IS THE WHOLE POINT OF THIS BATCH:
#      refinement DELAYS onset  => under-resolution; the mode is a resolution
#                                  artefact and a finer mesh cures it;
#      refinement ADVANCES onset => a GRID-SCALE discretisation-stability
#                                  problem -- refinement makes it worse, which
#                                  is never under-resolution.
#  The equal-order instability of rule 12b was identified exactly this way.
#
#  THE FACTORIAL. 2 (dx) x 2 (dt) x 2 (Jacobian) = 8 runs, and the BASELINE IS
#  IN THE BATCH rather than remembered (rule 38c). dx and dt are separated on
#  purpose: an L-stable integrator's damping grows with dt, so a dt-only change
#  can move onset for reasons that have nothing to do with the spatial mode
#  (rule 12c/15). Reading the 2x2 tells you which axis carries it.
#
#  WHY BOTH JACOBIANS. Measured on the baseline: hand and exact-AD produce
#  IDENTICAL eta to five digits and identical r0, and die at the SAME step --
#  the trajectory belongs to the residual, not the Jacobian (rule 17b). The AD
#  arm is therefore a CONTROL: it rules out the refinement result being an
#  artefact of the incomplete jacobian_u, which omits the {1,2,4,5} blocks
#  that :full adds to the residual. It is NOT expected to change the answer.
#  ⚠ AD costs ~10x per step (measured 120-156 s vs 12-14 s at nx=240).
#
#  DURATION 16 periods = 25.6 s, i.e. 2x the baseline onset. Long enough that
#  survival is evidence of delay, short enough to be affordable in the AD arm.
#  ⚠ A run that COMPLETES 25.6 s has not been shown stable -- only that onset
#  moved past 2x. Rule 12b's standing lesson: a threshold may just be a run that
#  ended too early (the 1:2 bar looked stable at 80 s and died at 137 s).
# ==============================================================
source "$(dirname "${BASH_SOURCE[0]}")/balfem_local.sh"

# ---- physics: IDENTICAL in all 8 runs -------------------------------------
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
export BALFEM_FE_ORDER=2
export BALFEM_P_ETA=1            # Q2/Q1 Taylor-Hood (rule 2b)
export BALFEM_SOLVER=sdirk
export BALFEM_TABLEAU=SDIRK_2_2
export BALFEM_PERIODS=16

# ---- the three varied knobs ------------------------------------------------
export BALFEM_NX=240            # dx = 0.25
#  Ly = dx keeps the single cross-cell ISOTROPIC (rule 12). Holding Ly fixed
#  while halving dx would give 2:1 cells, and rule 22 measured that anisotropy
#  costs ~760 GMRES iterations against ~480. With ny=1 + walls the solution is
#  y-invariant and Uy = 0 exactly, so Ly changes the CONDITIONING, not the physics.
export BALFEM_LY=0.25
export BALFEM_DT=0.04
export BALFEM_USE_AD=1        # EXACT AD Jacobian

#  VTK and diagnostics on a COMMON TIME GRID across the batch (0.2 s), so onset
#  times are compared on the same sampling rather than on the same step index.
export BALFEM_SAVE_EVERY=5
export BALFEM_DIAG_EVERY=5
export BALFEM_PRINT_EVERY=10

export BALFEM_OUTDIR=output/local_1d/P1LFE-2_1d_bcplane_nl_full_flat_Q2Q1-nx240_A0.1_T1.6_dt0.04_ad

balfem_local_run examples/local_1d/run_flume_1d.jl
