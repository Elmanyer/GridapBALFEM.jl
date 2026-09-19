#!/bin/bash
# ==============================================================
#  1-D NONLINEAR PRODUCTION RUN — nl_pressure=native, bar bed
#
#  PURPOSE: assess whether the FULLY NONLINEAR model stays stable over a long
#  integration, at production amplitude, with VTK saved densely enough to watch
#  the solution evolve. This is the run the suite does not have: §5.7 item 4 of
#  CLAUDE.md — "a long-duration nonlinear regression test. The mode needed 50-80 s
#  to emerge and no test runs that long; this is the gap that let it reach
#  production."
#
#  ⚠ FOUR CASES, ONE FACTORIAL — nl_pressure {native, full} x bed {flat, bar}.
#  Read them together: a divergence in one case is not evidence until the run
#  that differs in exactly ONE axis is beside it (rule 14c).
#
#  DISCRETISATION: Q2/Q1 Taylor-Hood (rules 2b / 12b). Chosen over Q3/Q2 because
#  this is a STABILITY question and Q2/Q1 is (a) the lowest Taylor-Hood pair, so
#  the most demanding of the inf-sup fix, and (b) the exact discretisation of the
#  80 s reference that resolved the instability
#  (output/local_1d/nonlinear_stability_taylor_hood/integrators/th_q2q1_sdirk22:
#  2000/2000 steps, settled eta 0.1059, Newton 3.99 it/step). This run extends
#  that baseline to 160 s and adds the pressure tier and the bathymetry.
#  ⚠ Q3/Q2 remains the CONVERGENCE-optimal production pairing (§5.3); a stable
#  Q2/Q1 result is not automatically a Q3/Q2 result.
#
#  AMPLITUDE: A = 0.10 m (kd = 5.5, ka ~ 0.16) — the amplitude rule 13 was lifted
#  to, and 100x the old cap. It is ~37 % of the Miche limit at the bar crest.
#
#  DURATION: 100 wave periods = 160 s = 4000 steps. The flume fills in 36 s
#  (45 m at c_g = 1.25 m/s, rule 14) and the ramp takes 2 T, so ~120 s = 75 T of
#  steady state is left — enough to fit an envelope growth rate over a window
#  that the 80 s reference could only half-cover.
#
#  SEQUENTIAL BY MEASUREMENT, not oversight: a direct LU beats every MPI split at
#  this size (6.45 s/step at 1 rank against 13-19 s/step at 2-12), and sequential
#  keeps the point gauges. The four cases run SIDE BY SIDE instead — see
#  run_all_1dprod.sh.
# ==============================================================
source "$(dirname "${BASH_SOURCE[0]}")/balfem_local.sh"

export BALFEM_WAVE_GEN=bc          # Dirichlet boundary generation, regular plane wave
export BALFEM_REGIME=nonlinear
export BALFEM_NL_PRESSURE=native
export BALFEM_FLAT_BED=0
export BALFEM_RELAX=1              # bc inflow: sponge_wL=0, the relaxation zone absorbs back-radiation

# ---- geometry: 60 m x 0.25 m, dx = 0.25 m (16 cells/lambda at kd=5.5) -----
export BALFEM_MPI=0
export BALFEM_LX=60.0
export BALFEM_NX=240
export BALFEM_NY=1                 # 1-D default: one cell across with SOLID WALLS (rule 12)
export BALFEM_YBC=wall
export BALFEM_LY=0.25              # = dx, so the single cell stays isotropic

# ---- BATHYMETRY: a big, square-cross-section submerged bar ----------------
#  h = 2.0 m on a still depth of 3.5 m, so the bar BLOCKS 57 % OF THE WATER
#  COLUMN and the crest depth is 1.5 m. Deliberately large: the point is that
#  the bathymetric effect be unmistakable in the VTK series, not marginal.
#    * span   x in [26, 34] m — 8 m = 2 carrier wavelengths of flat crest, so
#      the crest supports a standing-wave pattern rather than acting as a single
#      scatterer, and both shoulders are far from the sponge (45 m) and from the
#      inflow relaxation zone (one wavelength = 4 m).
#    * shoulder 0.5 m — the tanh crosses ~90 % of the height in ~1.5 m = 6 cells,
#      i.e. SQUARE-SHOULDERED to within what dx = 0.25 can carry. It is NOT a
#      vertical step and must not be made one: the residual carries grad(h)
#      explicitly, so a discontinuous bed is not representable by the model.
#  Local wave conditions over the crest: kd goes 5.5 -> 2.4, still far inside
#  P1LFE-2's applicable kd ~ 10.9, and Miche gives a_max ~ 0.27 m against the
#  0.10 m run here. So the case is demanding but neither out of the model's
#  validity range nor breaking.
export BALFEM_HBAR=2.0
export BALFEM_XBAR=30.0
export BALFEM_WBAR=4.0             # HALF-width: the bar spans xbar +/- wbar
export BALFEM_SBAR=0.5             # shoulder length; << wbar => box cross-section

# ---- wave, duration, discretisation ---------------------------------------
export BALFEM_AWAVE=0.1            # rule 13 is LIFTED: 0.10 m is the validated Taylor-Hood amplitude
export BALFEM_TWAVE=1.6            # kd = 5.5, lambda = 4.0 m
export BALFEM_PERIODS=100          # 160 s = 4000 steps at dt = 0.04
export BALFEM_DT=0.04
export BALFEM_FE_ORDER=2           # Q2 velocity
export BALFEM_P_ETA=1              # Q1 surface  => Taylor-Hood (rule 2b)
export BALFEM_SOLVER=sdirk
export BALFEM_TABLEAU=SDIRK_2_2
#  ⚠ SDIRK_2_2 is L-stable, i.e. DISSIPATIVE BY CONSTRUCTION (rule 15/12c): a run
#  that completes may be one whose numerical dissipation exceeded the growth rate.
#  A clean result here therefore wants the :theta or explicit-RK4 confirmation the
#  instability campaign used before "stable" is claimed unqualified.

# ---- output: VTK dense enough to ANIMATE ----------------------------------
#  save_every=5 at dt=0.04 => one snapshot every 0.2 s = 8 frames per wave
#  period, 800 snapshots over the run. The .pvd index is rewritten after every
#  snapshot (write_pvd_index), so the series opens in ParaView while the run is
#  still going and survives a kill.
export BALFEM_SAVE_EVERY=5
export BALFEM_WRITE_W=1            # w_s(sigma) reconstruction
export BALFEM_WRITE_PRESSURE=1     # p_s(sigma) reconstruction
export BALFEM_DIAG_EVERY=5         # diagnostics.csv row every 5 steps
export BALFEM_PRINT_EVERY=10


# ==============================================================
#  DIAGNOSTIC OVERRIDE — ONE VARIABLE: the bar's SHOULDER STEEPNESS.
#
#  WHY. run_1dprod_nl_native_bar.sh (shoulder 0.5 m) died at t = 37.2 s with
#  Newton stalled at a residual of 0.25 -- a HARD stall, not the :full pair's
#  near-miss at 2e-05. The growth was pinned at x = 33.5-34.25 m, i.e. the bar's
#  DOWNWAVE SHOULDER, with eta climbing 0.193 -> 0.237 -> 0.274 -> 0.372 there
#  while the background stayed at 0.104, and x_at_max did NOT migrate with the
#  crest (rule 40: never read growth without x_at_max -- a pinned location is a
#  stationary mode, not a travelling wave). The flat-bed run of the SAME physics
#  tier was clean at t = 54 s, 4 Newton it/step. So the bathymetry is the variable.
#
#  ⚠ THE SQUARE BAR'S SLOPE IS 2:1. max|grad h| = hbar/(2*sramp) = 2.0/1.0 = 2.0,
#  a 63-degree face. That is what "square cross-section" COSTS, and it is far
#  outside the mild-slope regime depth-integrated models are normally exercised on
#  (a Beji-Battjes bar is 1:20). The failure may therefore be the bathymetry I
#  specified rather than the nonlinear model.
#
#  THE TEST. sramp 0.5 -> 2.0 makes max|grad h| = 0.5 (1:2) -- still steep, 4x
#  gentler, same height, same span, same everything else:
#     * survives past t = 37 s  => the 2:1 face was the driver, and the model is
#       stable on this bathymetry class;
#     * fails at the same place => the shoulder steepness is not the variable and
#       the grad-h treatment itself is implicated.
#  ⚠ This does NOT replace the square-bar case, which is the one that was asked
#  for. It identifies what killed it.
export BALFEM_SBAR=2.0
#  The output grammar (rule 2c) carries no shoulder field, so without an explicit
#  path this would land beside the square-bar run as a bare _v2 and the two would
#  be indistinguishable after the fact -- exactly the failure the <discr> token
#  exists to prevent. Named explicitly instead, keeping the grammar and appending
#  the one thing that differs.
export BALFEM_OUTDIR=output/local_1d/P1LFE-2_1d_bcplane_nl_native_bar_Q2Q1-nx240_A0.1_T1.6_shoulder2.0

balfem_local_run examples/local_1d/run_flume_1d.jl
