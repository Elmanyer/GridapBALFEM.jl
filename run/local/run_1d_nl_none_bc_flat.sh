#!/bin/bash
# LOCAL quasi-1D flume - nl_none_bc_flat
#   Isolates nonlinear advection from the nonlinear pressure
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
export BALFEM_REGIME=nonlinear
export BALFEM_NL_PRESSURE=none
export BALFEM_FLAT_BED=1
export BALFEM_RELAX=1        # bc inflow: sponge_wL=0, so the relaxation zone is what absorbs the back-radiation

# sizing: 60 x 3 m, dx=0.25 (16 cells/lambda at kd=5.5)
# DURATION IS TRANSIT-BASED. Boundary generation injects at x=0, not at the
# interior source's x=18 m, so the wave has 45 m to cross before it reaches the
# right sponge instead of 27 m. At kd=5.5 the group speed is only c_g=1.25 m/s
# (CLAUDE.md rule 14, "the c_g transit trap"), so filling takes 36 s = 22.5 T,
# and +2 T ramp +3 T settle puts the earliest readable steady state at 27.5 T.
# 30 periods is that plus a small margin. The 16 periods these cases used with
# the INTERIOR source were correct for it (18.5 T needed) and would leave a
# boundary-generated run reading a still-filling flume.
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
export BALFEM_PERIODS=50
export BALFEM_SAVE_EVERY=10
export BALFEM_DIAG_EVERY=5

balfem_local_run examples/local_1d/run_flume_1d.jl
