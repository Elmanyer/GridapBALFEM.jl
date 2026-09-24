#!/bin/bash
#  ==== 2026-09-23 RK4 INTEGRATOR TEST (projected :full) — Q3/Q2, SDIRK_2_2 — the same-dt CONTROL for the Q3/Q2 RK4 arm ====
#  Derived from run_1dref_full_nx240_dt0005_hand_T100.sh (SDIRK_2_2, dt=0.005, Q2/Q1,
#  onset t = 40.40 s). Only the knobs named above differ.
# ==============================================================
#  :full FLAT — dx x dt REFINEMENT, HAND (pre-derived) JACOBIAN   [hand]
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
export BALFEM_FE_ORDER=3
export BALFEM_P_ETA=2            # Taylor-Hood (rule 2b)
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
export BALFEM_USE_AD=0        # HAND (pre-derived) Jacobian

#  VTK and diagnostics on a COMMON TIME GRID across the batch (0.2 s), so onset
#  times are compared on the same sampling rather than on the same step index.
export BALFEM_SAVE_EVERY=5
export BALFEM_DIAG_EVERY=5
export BALFEM_PRINT_EVERY=10

export BALFEM_OUTDIR=output/local_1d/P1LFE-2_1d_bcplane_nl_full_flat_Q2Q1-nx240_A0.1_T1.6_dt0.04_hand


# ==============================================================
#  EXTENDED dt LADDER — dt = 0.005s, hand Jacobian, dx = 0.25 UNCHANGED
#
#  LADDER POSITION. This fills the gap between the factorial's dt = 0.02 and the
#  dt = 0.001 run, giving four points at FIXED dx = 0.25:
#      dt = 0.04  -> onset 12.6 s
#      dt = 0.02  -> onset 16.8 s
#      dt = 0.005 -> THIS RUN
#      dt = 0.001 -> running
#  Four points is the minimum that can distinguish a trend from a pair of numbers,
#  and this one is the only member of the ladder that both completes in a working
#  day AND sits far enough from 0.02 to move the answer.
#
#  WHAT IT TESTS. The Class-III blocks {1,2,4,5} are carried by FROZEN L2
#  projections refreshed once per step, so their lag error falls with dt; that
#  predicts onset keeps moving out. Against it, SDIRK_2_2 is L-stable and its
#  numerical damping ALSO falls with dt (rule 12c), which exposes a growing mode
#  sooner. The two effects are separable only because dx is held fixed, and the
#  ladder is where they separate: a projection-driven mode gives a monotone trend,
#  a damping-driven one turns over.
#
#  ⚠ It REPLACES the dt = 1e-4 run, killed deliberately: at 36 h of wall clock per
#  simulated second it would have needed ~19 days to reach the baseline onset and
#  could not have produced a number. Its partial output remains in
#  ..._dt0.0001_hand/ and should be read as an aborted run, not a result.
#
#  COST: 5,120 steps to T_final = 25.6 s at ~13 s/step = ~18 h complete;
#  it passes the dt = 0.02 onset at 16.8 s in roughly 12 h.
#  Output cadence pinned to 0.2 s of SIMULATED time (40 steps), as elsewhere.
# ==============================================================
export BALFEM_DT=0.005
export BALFEM_SAVE_EVERY=40
export BALFEM_DIAG_EVERY=40
export BALFEM_PRINT_EVERY=80
export BALFEM_USE_AD=0
export BALFEM_OUTDIR=output/local_1d/P1LFE-2_1d_bcplane_nl_full_flat_Q2Q1-nx240_A0.1_T1.6_dt0.005_hand


# ==============================================================
#  LONG RE-RUN OF THE dt = 5e-3 CASE — 100 PERIODS
#
#  WHY. The 16-period run COMPLETED 5121/5121 steps to t = 25.6 s, but it was NOT
#  stable -- it was still growing when the window ended:
#      window        mean eta    mean u_max   Newton
#      t =  0-6       0.0742       0.2612      3.19
#      t =  6-12      0.1124       0.4140      4.00
#      t = 12-18      0.1179       0.4420      4.00
#      t = 18-25.6    0.1247       0.4755      4.00     <- still climbing
#  against the :native control's 0.393 over the same times. So 25.6 s is a BOUND
#  on the onset, not a stability result (rule 12b: a threshold may just be a run
#  that ended too early -- the 1:2 bar looked clean at 80 s and died at 137 s).
#
#  WHAT THIS RUN DECIDES. Onset at fixed dx = 0.25 goes 12.60 s (dt = 0.04) ->
#  21.20 s (dt = 0.02), a factor 1.68 per halving of dt. Two more halvings
#  extrapolates to ~60 s, which this window covers with margin:
#      * crashes near 60 s   => the dt trend is geometric and continues, and the
#        frozen-projection lag is behaving exactly as the hypothesis predicts;
#      * completes 160 s flat => dt = 5e-3 is genuinely stable at this dx, which
#        is a much stronger statement and directly comparable with the :native
#        flat 100-period control (eta 0.10494 -> 0.10300, Newton 5.12);
#      * crashes much earlier or later => the geometric reading is wrong and the
#        dt dependence needs its own explanation.
#
#  ⚠ 100 periods is chosen over 50 for COMPARABILITY: it is the same duration as
#  the :native control, so "stable" would mean the same thing for both. The
#  answer is expected around t ~ 60 s, i.e. ~2 days in, well before the ~105 h
#  the full window would take.
#
#  ⚠ SEPARATE OUTPUT DIRECTORY. BALFEM_OUTDIR bypasses unique_output_dir, so
#  reusing the 16-period path would silently overwrite a finished run (rule 2c --
#  that destroyed six runs on 2026-09-06).
export BALFEM_PERIODS=100
export BALFEM_OUTDIR=output/local_1d/c3r_q3q2_proj_sdirk

balfem_local_run examples/local_1d/run_flume_1d.jl
