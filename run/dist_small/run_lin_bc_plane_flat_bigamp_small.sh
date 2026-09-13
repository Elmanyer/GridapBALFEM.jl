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
#SBATCH --job-name="BALFEM_lin_bcplane_bigamp"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
# 4 GiB/rank is REQUIRED (measured peak 3633 MB/rank; 2 GB/core was OOM-killed).
# Sizing: 28 ranks/node x 4 GiB = 112 GiB/node = 4/8 of a rome node (28 total).
# Rules and tiers: run/SNELLIUS_ROME_LAUNCH_CONFIGS.md
# ⚠ --mem, NOT --mem-per-cpu (measured 2026-09-13). Slurm counts --mem-per-cpu
# against ALL ALLOCATED CPUs and rounds the allocation up to 16 cores, so
# --ntasks-per-node=28 --mem-per-cpu=4G actually asked for 32x4 = 128 GiB
# and was billed 5/8 (80 CPUs) to run 28 ranks. --mem is a flat per-node total:
#   112 GiB / 28 ranks = 4.00 GiB per rank, billed 4/8 (64 CPUs).
# See run/SNELLIUS_ROME_LAUNCH_CONFIGS.md §4.
#SBATCH --mem=112G
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
balfem_run 28 examples/distributed_small/run_bc_plane_small.jl
