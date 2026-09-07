#!/bin/bash
# SMALL run — BC-generated LINEAR plane wave, flat bed, A=0.1, EXTENDED to 30 periods.
#
# PURPOSE. The missing cell of the modulational-instability test. The nonlinear
# :full run at this exact amplitude and generation grows exponentially after
# t~25 s (envelope rate +0.087 1/s, e-folding 7.9 s); the linear runs available
# so far are all at A=0.001, so "linear is stable" is confounded with amplitude.
# This run is linear at the SAME delivered amplitude, on the SAME sponge, so a
# flat envelope here isolates the growth as a NONLINEAR effect and simultaneously
# rules out sponge reflection (a linear boundary process would show here too).
#  SNELLIUS ROME SIZING: 28 ranks x 4 GiB = 112 GiB = exactly 4/8 of a node.
#  Was 32 ranks, which billed 5/8 (memory fraction 4.57/8 rounds up). The rome node is
#  128 cores / 224 GiB divided into EIGHTHS, and the LARGER of the core/memory fraction
#  sets the bill -- at 4 GiB/rank memory always wins, so the cost-optimal count is 7k
#  ranks for tier k/8. See run/SNELLIUS_ROME_LAUNCH_CONFIGS.md.
#SBATCH --job-name="BALFEM_lin_bcplane_bigamp"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
# Memory: the rome node default (2 GB/core) was TESTED (2026-08) and OOM-killed.
# 4 GB/core is required and is MEASURED: the small-domain runs peaked at 3633
# MB/rank, so 3 GB/core would OOM too. RSS rose over the first ~100 steps then
# PLATEAUED -- no per-step leak. See building_files/OPEN_ITEMS.md section 1.
# A rome node advertises 256 GB but SLURM can allocate only ~224 GB, so
# 64 ranks/node x 4 GB = 256 GB/node is REFUSED at submit time. 32 ranks x 4 GB
# = 128 GB/node fits and is the proven small-domain configuration.
#SBATCH --mem-per-cpu=4G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
source $HOME/GridapBALFEM.jl/run/balfem_env.sh
# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
export BALFEM_PX=7
export BALFEM_PY=4            # 7*4 = 28 ranks
# --- Case-specific overrides ---
export BALFEM_REGIME=linear
export BALFEM_NL_PRESSURE=none
export BALFEM_FLAT_BED=1
export BALFEM_AWAVE=0.1
#  30 periods = 60 s. Transit at c_g=1.58 is 31.7 s, +4 s ramp +3T settle = 42 s,
#  leaving ~18 s of steady state = 2.3 e-foldings of the nonlinear growth rate.
#  The original 16-period run stopped just as the instability became visible.
export BALFEM_PERIODS=30
balfem_run 32 examples/distributed_small/run_bc_plane_small.jl
