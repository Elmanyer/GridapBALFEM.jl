#!/bin/bash
#  A3 — A1 at dx=0.125 (fine leg), SKEW ON. THE decisive run: refinement must no longer make it worse.
# ==============================================================
#  TIER-A of the SKEW-ADVECTION evaluation
#  (building_files/SKEW_SYMMETRIC_ADVECTION_PLAN.md §5.2)
#
#  THE QUESTION. The plain Galerkin advection operator has an exact energy
#  production −½∫∇·(Hū)(ΣMᵢⱼuᵢ·uⱼ) that the discrete continuity equation cannot
#  cancel, because the cancellation needs continuity tested against a QUARTIC that
#  is not in the η space. BALFEM_SKEW=1 adds the consistent Temam-type correction
#  that removes it. Does that clear the grid-scale instability?
#
#  ⚠ READ A1/A2/A3 AS A SET, NOT INDIVIDUALLY. The baseline signature is that
#  REFINING dx MAKES IT WORSE (blow-up t≈24 s at dx=0.25, t≈10 s at dx=0.125).
#  Success is not "A1 completes" — it is the dx-ORDERING DISAPPEARING: all three
#  complete AND their eta_max(t) curves converge toward each other as dx falls.
#  If A1 completes but A3 still dies, the leak is reduced, not removed.
#
#  ⚠ RULE 38c. These scripts are the reference 1dnl launchers copied VERBATIM with
#  exactly one variable added. Do not "tidy" them by rebuilding the variable list:
#  a partial list previously flipped generation from :bc to the interior source,
#  which delivers 2.8x the requested amplitude, and produced a confident result
#  pointing the wrong way.
#
#  Baseline (skew off) for reference — see the plan §5.1 for the full table:
#     nl_none A=0.05  survives 80 s but DRIFTS 0.0512 → 0.0821
#     nl_none A=0.10  DIES at t=26.4 s (eta → 5.06)
#     linear  A=0.10  flat to four digits for 80 s
# ==============================================================
source "$(dirname "${BASH_SOURCE[0]}")/balfem_local.sh"
export BALFEM_WAVE_GEN=bc
export BALFEM_FLAT_BED=1
export BALFEM_RELAX=1
export BALFEM_MPI=0
export BALFEM_LX=60.0
export BALFEM_NY=1
export BALFEM_YBC=wall
export BALFEM_LY=0.25
export BALFEM_PERIODS=50
export BALFEM_SAVE_EVERY=25
export BALFEM_DIAG_EVERY=5
export BALFEM_REGIME=nonlinear
export BALFEM_NL_PRESSURE=none
export BALFEM_SKEW=1
export BALFEM_NX=480
export BALFEM_AWAVE=0.1
export BALFEM_OUTDIR=output/local_1d/skewA3_none_A010_dx0125
balfem_local_run examples/local_1d/run_flume_1d.jl
