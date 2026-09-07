#!/bin/bash

#SBATCH --job-name="GridapBALFEM"
#SBATCH --partition=rome
#SBATCH --time=72:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=84
#SBATCH --ntasks-per-node=42
#SBATCH --cpus-per-task=1
# 4 GiB/rank is REQUIRED (measured peak 3633 MB/rank; 2 GB/core was OOM-killed).
# Sizing: 42 ranks/node x 4 GiB = 168 GiB/node = 6/8 of a rome node (84 total).
# Rules and tiers: run/SNELLIUS_ROME_LAUNCH_CONFIGS.md
#SBATCH --mem-per-cpu=4G
#SBATCH --output=GridapBALFEM.%j.out
#SBATCH --error=GridapBALFEM.%j.err

source $HOME/GridapBALFEM.jl/run/balfem_env.sh

# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
#  Set explicitly: the driver defaults to 32*4 = 128, which no longer matches
#  the rank count after the memory-driven cut to 96.
export BALFEM_PX=14
export BALFEM_PY=6            # 14*6 = 84 ranks (2 nodes)

balfem_run 84 examples/distributed/run_plane_wave_dist.jl
