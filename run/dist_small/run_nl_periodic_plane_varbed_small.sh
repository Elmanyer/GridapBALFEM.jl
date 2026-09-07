#!/bin/bash
# SMALL run — nonlinear full plane wave over a bar, A=0.1
#SBATCH --job-name="BALFEM_nl_plane_varbed"
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

# --- Case-specific overrides (only what differs from the script's base) ---
export BALFEM_FLAT_BED=0

balfem_run 28 examples/distributed_small/run_periodic_plane_small.jl
