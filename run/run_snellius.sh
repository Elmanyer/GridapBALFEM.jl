#!/bin/bash

#SBATCH --job-name="GridapBALFEM"
#SBATCH --partition=rome
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=42
#SBATCH --ntasks-per-node=42
#SBATCH --cpus-per-task=1
# 4 GiB/rank is REQUIRED (measured peak 3633 MB/rank; 2 GB/core was OOM-killed).
# Sizing: 42 ranks x 4 GiB = 168 GiB = 6/8 of ONE rome node.
# Rules and tiers: run/SNELLIUS_ROME_LAUNCH_CONFIGS.md
# ⚠ --mem, NOT --mem-per-cpu (measured 2026-09-13). Slurm counts --mem-per-cpu
# against ALL ALLOCATED CPUs and rounds the allocation up to 16 cores, so
# --ntasks-per-node=42 --mem-per-cpu=4G actually asked for 48x4 = 192 GiB
# and was billed 7/8 (112 CPUs) to run 42 ranks. --mem is a flat per-node total:
#   168 GiB / 42 ranks = 4.00 GiB per rank, billed 6/8 (96 CPUs).
# See run/SNELLIUS_ROME_LAUNCH_CONFIGS.md §4.
#SBATCH --mem=168G
#SBATCH --output=GridapBALFEM.%j.out
#SBATCH --error=GridapBALFEM.%j.err

source $HOME/GridapBALFEM.jl/run/balfem_env.sh

# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
#  Set explicitly: the driver defaults to 32*4 = 128, which no longer matches
#  the rank count after the memory-driven cut to 96.
export BALFEM_PX=7
export BALFEM_PY=6            # 7*6 = 42 ranks

balfem_run 42 examples/distributed/run_plane_wave_dist.jl
