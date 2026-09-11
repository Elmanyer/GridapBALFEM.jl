#!/bin/bash
# SMALL run — linear directional sea, flat, Hs=0.2
#SBATCH --job-name="BALFEM_lin_directional"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --ntasks=56
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --cpus-per-task=1
# 4 GiB/rank is REQUIRED (measured peak 3633 MB/rank; 2 GB/core was OOM-killed).
# Sizing: 56 ranks/node x 4 GiB = 224 GiB/node = 8/8 of a rome node (56 total).
# Rules and tiers: run/SNELLIUS_ROME_LAUNCH_CONFIGS.md
#
# >> WHY 56 AND NOT 42 (2026-09-11). Ly 20 -> 40 m doubles the mesh to
# >> 200x160 = 32000 cells. Memory, not cores, sets the tier: the measured cost
# >> is ~3.8 MB per cell per rank above a ~1450 MB Julia+Gridap baseline, so
# >>     42 ranks -> 762 cells/rank -> ~4350 MB  => OVER the 4 GiB cap
# >>     56 ranks -> 571 cells/rank -> ~3620 MB  => 88 % of the cap
# >> 56 is also the cost-optimal full-node count (7k ranks for tier k/8).
# >> IF IT IS OOM-KILLED (exit 137, or 'oom-kill event' in .err): go to 2 nodes,
# >> --ntasks=112 --nodes=2 --ntasks-per-node=56 with PX=14 PY=8, which halves
# >> cells/rank to 286. Do NOT coarsen the mesh: dx=0.25 m is 14.7 cells per the
# >> shortest wavelength (lambda=3.67 m at kd=5.99) and that short-wave tail is
# >> exactly what a directional run exists to carry.
#SBATCH --mem-per-cpu=4G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source $HOME/GridapBALFEM.jl/run/balfem_env.sh

# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
export BALFEM_PX=8
export BALFEM_PY=7            # 8*7 = 56 ranks

# --- Case-specific overrides (only what differs from the script's base) ---
export BALFEM_REGIME=linear
export BALFEM_NL_PRESSURE=none

balfem_run 56 examples/distributed_small/run_directional_sea_small.jl
