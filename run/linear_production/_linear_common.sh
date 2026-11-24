# ============================================================================
#  _linear_common.sh — settings shared by every run/linear_production launcher
#
#  Sourced AFTER run/balfem_env.sh and BEFORE balfem_run (CLAUDE.md rule 38h:
#  an export after the balfem_run line is a no-op). Change a production-wide
#  choice here, once, rather than in seven launchers.
#
#  PHYSICS: the linear model (models 1 flat / 2 variable bed). Its residual and
#  both Jacobians are identical to v1 entry by entry (regression snapshot
#  lin_flat / lin_bed / lin_bed_q21: rel diff 0.00e+00), and it is verified at
#  optimal order in 2-D, sequential and distributed (archive/HISTORY_V1.md §5.1).
# ============================================================================

export BALFEM_REGIME=linear
export BALFEM_NL_PRESSURE=0           # v2 is Boolean 0/1; the v1 names none/native/full are refused

# --- Discretisation: the v1 small-suite pair, P1LFE-2 + Q2/Q1 Taylor-Hood -----
export BALFEM_M=2
export BALFEM_P_VERT=1
export BALFEM_FE_ORDER=2              # p_eta = FE_ORDER-1 is derived and gated by balfem_env.sh
export BALFEM_DT=0.02

# --- Integrator: Crank-Nicolson ----------------------------------------------
#  The linear model has no instability for a dissipative integrator to mask, so
#  the default SDIRK_2_2 only costs here. It loses ~0.75(w dt)^4 per step (rule
#  15): ~1e-5/step at the peak (Tp=2 s) but ~2e-4/step at the kd~6 tail, i.e.
#  ~20 % of the short waves over one domain transit -- exactly the part of a
#  directional/irregular sea worth looking at. CN is non-dissipative and has ONE
#  implicit stage per step against SDIRK_2_2's two, so it is also ~2x cheaper.
#  The sponge and relaxation zone peak at mu_max*dt = 0.8, far from the stiff
#  limit where CN's R(inf) = -1 would ring.
#  To revert to the v1 integrator: BALFEM_SOLVER=sdirk.
export BALFEM_SOLVER=theta

# --- Production safety --------------------------------------------------------
#  Refuse to start on a sysimage that no longer matches src/ (balfem_env.sh):
#  a stale image silently runs the OLD solver, found out only after the run.
export BALFEM_STRICT_SYSIMAGE=1
export BALFEM_PRINT_EVERY=50
