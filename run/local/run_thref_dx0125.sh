#!/bin/bash
#  Q2/Q1 Taylor-Hood, SDIRK_2_2, nx=480 — FINE leg — THE DECISIVE ONE (equal order died here at t~10 s)
# ==============================================================
#  TAYLOR-HOOD dx-REFINEMENT LEG  (1-D ONLY, TH ONLY — CLAUDE.md rule 2b)
#
#  THE DEFINING SIGNATURE OF THE INSTABILITY IS THAT REFINING dx MAKES IT WORSE:
#  on EQUAL-ORDER Q2/Q2 the run died at t~24 s for dx=0.25 and t~10 s for dx=0.125,
#  at an identical delivered amplitude. That inversion is what excluded ordinary
#  under-resolution and forced a discretisation-stability explanation.
#
#  ⚠ A CURE MUST THEREFORE INVERT IT. It is NOT enough that a Taylor-Hood run
#  completes at one mesh: dx=0.125 must be NO WORSE than dx=0.25. If the fine leg
#  still degrades, the pairing has delayed the mode, not removed it.
#
#  Reference points, EQUAL ORDER, same A and integrator:
#      dx=0.50   essentially suppressed
#      dx=0.25   DIED t=26.4 s, eta -> 5.06
#      dx=0.125  DIED t~10 s   (control run skewA4: eta=0.79 at t=8.4)
#  Taylor-Hood dx=0.25 (TH1) is flat at 0.105 through t=14.8 so far.
#
#  dt IS HELD FIXED at 0.04 across the ladder ON PURPOSE: NONLINEAR_INSTABILITY.md
#  F4 establishes the growth is dt-insensitive at fixed dx, so holding dt keeps this
#  a ONE-VARIABLE comparison. The dt axis is a separate test.
source "$(dirname "${BASH_SOURCE[0]}")/balfem_local.sh"
export BALFEM_WAVE_GEN=bc
export BALFEM_REGIME=nonlinear
export BALFEM_NL_PRESSURE=none
export BALFEM_FLAT_BED=1
export BALFEM_RELAX=1
export BALFEM_MPI=0
export BALFEM_LX=60.0
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
export BALFEM_NX=480
export BALFEM_OUTDIR=output/local_1d/thref_dx0125
balfem_local_run examples/local_1d/run_flume_1d.jl
