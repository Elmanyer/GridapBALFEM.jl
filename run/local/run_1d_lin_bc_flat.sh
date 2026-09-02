#!/bin/bash
# LOCAL quasi-1D flume - lin_bc_flat
#   Dirichlet boundary generation in isolation
# SEQUENTIAL, and that is a MEASURED choice, not an oversight. Strong scaling on
# this exact 240x3 mesh (20.2k DOFs, 40 steps):
#     ranks    1(LU)     2       4       6      12
#     s/step    6.45   14.44   13.10   19.53   18.27
#     gmres      --    758     758     758     715-774
# The GMRES iteration count is RANK-INDEPENDENT (~760), so the cost is the weak
# Jacobi preconditioner, not the decomposition: a direct LU factorisation beats
# it 2-3x at this problem size. The 12 cores are used by running several CASES
# side by side instead (see run_all_1d.sh), which is strictly faster than
# 12-way-decomposing one case. Sequential also KEEPS the point gauges.
source "$(dirname "${BASH_SOURCE[0]}")/balfem_local.sh"

export BALFEM_WAVE_GEN=bc
export BALFEM_REGIME=linear
export BALFEM_NL_PRESSURE=none
export BALFEM_FLAT_BED=1
export BALFEM_RELAX=1

# 12-rank sizing: 60 x 3 m, dx=0.25 (16 cells/lambda at kd=5.5), 16 periods
export BALFEM_MPI=0           # direct LU: measured 2-3x faster than any MPI split here
export BALFEM_LX=60.0
export BALFEM_NX=240
# 1-D HORIZONTAL DEFAULT: ny=1 with SOLID WALLS, not 3 cells with periodicity.
# For a normal-incidence wave the exact solution has Uy == 0, and the wall
# condition Uy = 0 is EXACTLY consistent with it -- it approximates nothing.
# :periodic merely PERMITS Uy == 0 while also admitting a family of y-periodic
# modes a true 1-D model does not have. ny=1 + :wall is also 2.8x cheaper:
# 7215 free DOFs against 20202 for ny=3 periodic on this same 240-cell mesh,
# and a direct LU costs more than linearly in DOFs. Ly=0.25 keeps the single
# cell isotropic (dx=dy=0.25). Use :periodic only for genuinely oblique or
# short-crested content, where a solid wall would reflect.
export BALFEM_NY=1
export BALFEM_YBC=wall
export BALFEM_LY=0.25
# DURATION. The driver's transit-aware default for boundary generation is 26 T,
# which counts the 22.5 T fill (45 m to the right sponge at c_g = 1.25 m/s) plus
# ~3.5 T. It does NOT separately count the 2 T Hann ramp, so by this project's
# own convention (transit + 3 T settle, after the ramp) the earliest readable
# steady state is 27.5 T -- and 26 T leaves only ~1.5 T of it. Pinned to 30 T
# Pinned to 50 T (80 s), which leaves ~40 s = 25 T of steady state after the
# 36 s fill and the 2 T ramp -- enough to fit an envelope growth rate, which the
# ~5 T that 30 T would leave is not. Pinned rather than defaulted so all six
# 1-D cases share one duration and stay comparable.
export BALFEM_PERIODS=50
export BALFEM_SAVE_EVERY=10
export BALFEM_DIAG_EVERY=5

balfem_local_run examples/local_1d/run_flume_1d.jl
