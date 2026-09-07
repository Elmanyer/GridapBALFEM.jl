#!/bin/bash
# SMALL run — bcplane nonlinear :full, flat bed, A=0.10 — AMPLITUDE LADDER point.
#  THE AMPLITUDE LADDER. Benjamin-Feir growth goes as (ka)^2, so the envelope
#  rate must scale with the SQUARE of the delivered amplitude. Predicted, anchored
#  on the measured +0.087 1/s at A=0.1:
#
#        A_nom    delivered eta   predicted rate   e-folding
#        0.05         ~0.088        0.022 1/s        32 s
#        0.10         ~0.175        0.087 1/s         8 s     <- measured
#        0.15         ~0.263        0.196 1/s       3.5 s
#
#  A factor 9 in rate across a factor 3 in amplitude is a signature very little
#  else reproduces. Fit the rate the same way for all three: log-linear on the
#  per-period envelope peak of eta_max_int over the last 10 s.
#
#  >> WHAT THIS LADDER CANNOT DO. The quasi-Newton gap in jacobian_u ALSO scales
#  >> as A^2, so an A^2 law is necessary but NOT sufficient evidence for BF. The
#  >> tolerance runs (run_nl_full_bc_plane_{tight,floor}tol_small.sh) are what
#  >> separate those two; this ladder and those runs are a PAIR.
#
#  Expect the A=0.15 case to diverge before 60 s -- that is fine and expected.
#  The rate is measured during growth, not at saturation.
#  SNELLIUS ROME SIZING: 28 ranks x 4 GiB = 112 GiB = exactly 4/8 of a node.
#  Was 32 ranks, which billed 5/8 (memory fraction 4.57/8 rounds up). The rome node is
#  128 cores / 224 GiB divided into EIGHTHS, and the LARGER of the core/memory fraction
#  sets the bill -- at 4 GiB/rank memory always wins, so the cost-optimal count is 7k
#  ranks for tier k/8. See run/SNELLIUS_ROME_LAUNCH_CONFIGS.md.
#SBATCH --job-name="BALFEM_nl_full_bcplane_A0.10"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
# Memory: 4 GB/core is required and MEASURED (this exact case peaked at 2882 MB
# on 32 ranks / 250 cells per rank). 32 x 4 GB = 128 GB/node, well inside the
# ~224 GB a rome node can allocate. See building_files/OPEN_ITEMS.md section 1.
#SBATCH --mem-per-cpu=4G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source $HOME/GridapBALFEM.jl/run/balfem_env.sh
# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
export BALFEM_PX=7
export BALFEM_PY=4            # 7*4 = 28 ranks
# --- the reference configuration, held fixed across the whole ladder ---
export BALFEM_REGIME=nonlinear
export BALFEM_NL_PRESSURE=full
export BALFEM_FLAT_BED=1
export BALFEM_PERIODS=30      # 60 s; the original run stopped at 32 s
export BALFEM_AWAVE=0.10
#  OUTDIR set explicitly for every rung so the three sit together and the A=0.10
#  rung cannot overwrite the existing 16-period production run, whose auto-tag
#  (which carries A but NOT the duration) is otherwise identical.
export BALFEM_OUTDIR=$HOME/GridapBALFEM.jl/output/small_bcplane_nlfull_ampladder_A0.10_P1LFE-2
balfem_run 32 examples/distributed_small/run_bc_plane_small.jl
