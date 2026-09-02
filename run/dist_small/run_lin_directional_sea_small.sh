#!/bin/bash
# SMALL run — linear directional sea, flat, Hs=0.2
#SBATCH --job-name="BALFEM_lin_directional"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --ntasks=48
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
# Memory: the rome node default (2 GB/core) was TESTED (2026-08) and the job was
# still OOM-killed. 4 GB/core is required, and is now MEASURED rather than
# argued: the 2026-08-30 directional run peaked at 3633 MB/rank, so 3 GB/core
# would OOM as well. RSS rose over the first ~100 steps and then PLATEAUED --
# no per-step leak. See building_files/OPEN_ITEMS.md section 1.
#
# RANK COUNT HERE IS MEMORY-BOUND, NOT PHYSICS-BOUND. A rome node advertises
# 256 GB but SLURM can only allocate ~224 GB of it, so the old request of
# 64 ranks x 4 GB = 256 GB/node was REFUSED at submit time:
#   "Requested node configuration is not available"
# 48 ranks x 4 GB = 192 GB/node fits, and is the decomposition the nonlinear
# twin actually ran at (cpu_grid=(8,6), see output/small_directional_*). Keep
# the linear and nonlinear directional runs on the SAME grid so the linear one
# stays a valid control. To recover all 64 ranks instead, spread them over two
# nodes: set nodes to 2 and ntasks-per-node to 32 (= 128 GB/node).
# Do NOT drop mem-per-cpu on the argument that the sysimage removes the
# per-rank compile: that argument was tested and refuted.
#SBATCH --mem-per-cpu=4G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source $HOME/GridapBALFEM.jl/run/balfem_env.sh

# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
export BALFEM_PX=8
export BALFEM_PY=6            # 8*6 = 48 ranks

# --- Case-specific overrides (only what differs from the script's base) ---
export BALFEM_REGIME=linear
export BALFEM_NL_PRESSURE=none

balfem_run 48 examples/distributed_small/run_directional_sea_small.jl
