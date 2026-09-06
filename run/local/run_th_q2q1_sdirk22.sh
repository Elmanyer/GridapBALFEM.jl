#!/bin/bash
#  TH1 - Q2/Q1, SDIRK_2_2. Minimal change from the failing run: same order, TH pairing.
# ==============================================================
#  TAYLOR-HOOD NONLINEAR-INSTABILITY CAMPAIGN  (1-D ONLY, TH ONLY)
#
#  THE FINDING THAT PROMPTED IT. Every nonlinear-instability run before 2026-09-05
#  used EQUAL-ORDER Q2/Q2 spaces, while the ENTIRE verified scope of this solver
#  (30/30 spatial MMS studies, five vertical bases, eight models) was measured on
#  Q3/Q2. eta enters momentum undifferentiated via div(v) after IBP, so it plays the
#  pressure role of a Stokes system and equal-order continuous spaces are inf-sup
#  deficient. Their classic failure mode is a checkerboard at lambda ~ 2*dx --
#  exactly where the measured growth spectrum peaks (2-3*dx). CLAUDE.md rule 2b.
#
#  So the blown-up configuration and the verified configuration were NEVER THE SAME
#  DISCRETISATION. This campaign re-runs the failing case on Taylor-Hood pairs only.
#
#  REFERENCE (equal order, Q2/Q2, SDIRK_2_2, the same A and dx):
#      DIED at t = 26.4 s, eta -> 5.06   (delivered amplitude ~0.102)
#  A Taylor-Hood run that completes 80 s near 0.10 is the signal.
#
#  ⚠ ALL RUNS HERE ARE 1-D (ny=1, solid walls) AND TAYLOR-HOOD. No exceptions.
#  ⚠ Do NOT compare any of these against an equal-order result and call the
#    difference physics -- that is a two-variable comparison (rule 2b corollary).
source "$(dirname "${BASH_SOURCE[0]}")/balfem_local.sh"
export BALFEM_WAVE_GEN=bc
export BALFEM_REGIME=nonlinear
export BALFEM_NL_PRESSURE=none
export BALFEM_FLAT_BED=1
export BALFEM_RELAX=1
export BALFEM_MPI=0
export BALFEM_LX=60.0
export BALFEM_NX=240
export BALFEM_NY=1
export BALFEM_YBC=wall
export BALFEM_LY=0.25
export BALFEM_AWAVE=0.1
export BALFEM_PERIODS=50
export BALFEM_SAVE_EVERY=25
export BALFEM_DIAG_EVERY=5
export BALFEM_FE_ORDER=2
export BALFEM_P_ETA=1
export BALFEM_SOLVER=sdirk
export BALFEM_TABLEAU=SDIRK_2_2
export BALFEM_OUTDIR=output/local_1d/th_q2q1_sdirk22
balfem_local_run examples/local_1d/run_flume_1d.jl
