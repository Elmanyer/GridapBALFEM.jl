#!/bin/bash
# SMALL run — nonlinear full directional sea, flat, Hs=0.2 (BASE)
#SBATCH --job-name="BALFEM_nl_directional"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --ntasks=42
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=42
#SBATCH --cpus-per-task=1
# 4 GiB/rank is REQUIRED (measured peak 3633 MB/rank; 2 GB/core was OOM-killed).
# Sizing: 42 ranks/node x 4 GiB = 168 GiB/node = 6/8 of a rome node (42 total).
# Rules and tiers: run/SNELLIUS_ROME_LAUNCH_CONFIGS.md
#SBATCH --mem-per-cpu=4G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source $HOME/GridapBALFEM.jl/run/balfem_env.sh

# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
export BALFEM_PX=7
export BALFEM_PY=6            # 7*6 = 42 ranks

# --- Case-specific overrides (only what differs from the script's base) ---
# (base config — no overrides needed)

balfem_run 42 examples/distributed_small/run_directional_sea_small.jl
