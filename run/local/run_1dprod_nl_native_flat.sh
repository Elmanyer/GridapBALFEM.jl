#!/bin/bash
# ==============================================================
#  1-D NONLINEAR PRODUCTION RUN — nl_pressure=native, flat bed
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
export BALFEM_FLAT_BED=1
export BALFEM_RELAX=1              # bc inflow: sponge_wL=0, the relaxation zone absorbs back-radiation

# ---- geometry: 60 m x 0.25 m, dx = 0.25 m (16 cells/lambda at kd=5.5) -----
export BALFEM_MPI=0
export BALFEM_LX=60.0
export BALFEM_NX=240
export BALFEM_NY=1                 # 1-D default: one cell across with SOLID WALLS (rule 12)
export BALFEM_YBC=wall
export BALFEM_LY=0.25              # = dx, so the single cell stays isotropic

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

balfem_local_run examples/local_1d/run_flume_1d.jl
