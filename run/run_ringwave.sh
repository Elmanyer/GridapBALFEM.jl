#!/bin/bash

#SBATCH --job-name="BALFEM_ringwave"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
# Memory: the rome node default (2 GB/core) was TESTED (2026-08) and the job was
# still OOM-killed. 4 GB/core is required, and is now MEASURED rather than
# argued: the 2026-08-30 small-domain runs peaked at 3633 MB/rank, so 3 GB/core
# would OOM as well. RSS rose over the first ~100 steps and then PLATEAUED --
# no per-step leak. See building_files/OPEN_ITEMS.md section 1.
#
# RANK COUNT HERE IS MEMORY-BOUND, NOT PHYSICS-BOUND. A rome node advertises
# 256 GB but SLURM can only allocate ~224 GB of it, so the old request of
# 64 ranks/node x 4 GB = 256 GB/node was REFUSED at submit time:
#   "Requested node configuration is not available"
# Cut to 48 ranks/node x 4 GB = 192 GB/node (64 -> 48 ranks over 1 node(s)),
# rather than spreading the old rank count over more nodes.
# Do NOT drop mem-per-cpu on the argument that the sysimage removes the
# per-rank compile: that argument was tested and refuted.
#SBATCH --mem-per-cpu=4G
#SBATCH --output=GridapBALFEM.%j.out
#SBATCH --error=GridapBALFEM.%j.err

source $HOME/GridapBALFEM.jl/run/balfem_env.sh

# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
export BALFEM_PX=8
export BALFEM_PY=6              # 8*6 = 48 ranks; mesh 200x200 -> 25x33 cells/rank

# --- Case-specific knobs (defaults shown; uncomment to override) ------------
# export BALFEM_L=200           # basin side length [m]
# export BALFEM_D=3.5           # still-water depth [m]
# export BALFEM_TWAVE=1.6       # wave period [s]
# export BALFEM_AWAVE=0.001     # amplitude [m]
# export BALFEM_SPONGE=20       # sponge width, all 4 sides [m]
# export BALFEM_MUMAX=10        # sponge strength
# export BALFEM_PERIODS=20      # run length in wave periods
# export BALFEM_NX=200; export BALFEM_NY=200

# --- Solver knobs (RungeKutta :SDIRK_2_2 defaults; bump if Newton stalls) ----
# export BALFEM_LS_MAXITER=4000 # raise GMRES cap if convergence warnings appear
# export BALFEM_NL_TOL=1e-6

balfem_run 48 examples/distributed/run_ring_wave_dist.jl
