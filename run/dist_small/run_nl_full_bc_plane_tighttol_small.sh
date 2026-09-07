#!/bin/bash
# SMALL run — bcplane nonlinear :full A=0.1, TIGHT Newton tolerance (1e-7).
#  WHY THIS TEST. jacobian_u OMITS EVERY nonlinear-pressure (N) derivative by
#  design -- problem.jl states it: "the N blocks add to the residual but not
#  here". Newton is therefore QUASI-Newton, and the size of the omitted block
#  scales as A^2. That matters because the envelope growth under investigation
#  ALSO scales as A^2 (Benjamin-Feir goes as (ka)^2), so an amplitude sweep alone
#  CANNOT separate the two. This test can: a quasi-Newton artefact is a property
#  of the SOLVE and must move when the tolerance moves; genuine modulational
#  instability is a property of the OPERATOR and must not.
#
#  READING IT. Compare the envelope growth rate against the production run
#  (nl_tol=1e-5, rate +0.087 1/s):
#    rate essentially unchanged  -> the growth is in the operator, i.e. physical;
#    rate changes materially     -> Newton is converging to the fixed point of
#                                   the wrong map (CLAUDE.md rule 5).
#
#  REACHABILITY. The production run achieved a minimum Newton residual of
#  4.6e-08 (median 7.1e-07) while merely stopping at its 1e-5 threshold, so the
#  quasi-Newton floor is at most 4.6e-08.
#  1e-7 is 100x tighter than production and sits ABOVE the observed 4.6e-08
#  floor, so it should converge. Its companion run at 1e-9 probes below the
#  floor deliberately.
#SBATCH --job-name="BALFEM_nl_full_bcplane_tol1e7"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
# 4 GiB/rank is REQUIRED (measured peak 3633 MB/rank; 2 GB/core was OOM-killed).
# Sizing: 28 ranks/node x 4 GiB = 112 GiB/node = 4/8 of a rome node (28 total).
# Rules and tiers: run/SNELLIUS_ROME_LAUNCH_CONFIGS.md
#SBATCH --mem-per-cpu=4G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source $HOME/GridapBALFEM.jl/run/balfem_env.sh
# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
export BALFEM_PX=7
export BALFEM_PY=4            # 7*4 = 28 ranks
# --- the reference configuration, held fixed across this whole family ---
export BALFEM_REGIME=nonlinear
export BALFEM_NL_PRESSURE=full
export BALFEM_FLAT_BED=1
#  30 periods = 60 s. The original run stopped at 32 s, just as the envelope
#  growth became visible (rate +0.087 1/s, e-folding 7.9 s). 60 s leaves ~18 s of
#  steady state after transit+ramp+settle, i.e. ~2.3 further e-foldings.
export BALFEM_PERIODS=30
export BALFEM_AWAVE=0.1
export BALFEM_NL_TOL=1e-7     # production default is 1e-5
export BALFEM_NL_ITER=100     # was 50: a tighter tolerance needs more iterations
export BALFEM_OUTDIR=$HOME/GridapBALFEM.jl/output/small_bcplane_nlfull_A0.1_nltol1e-7_P1LFE-2
balfem_run 28 examples/distributed_small/run_bc_plane_small.jl
