#!/bin/bash

#SBATCH --job-name="BALFEM_periodicplanewave"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
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
export BALFEM_PX=14
export BALFEM_PY=6            # 14*6 = 84 ranks (2 nodes)

# --- Case-specific knobs (periodic-width flume; defaults shown) --------------
#     The y-edges are periodic (set in the example) -> no lateral sponges.
# export BALFEM_LX=400; export BALFEM_LY=20   # domain size [m]
# export BALFEM_D=3.5           # still-water depth [m]
# export BALFEM_TWAVE=1.6       # wave period [s]
# export BALFEM_AWAVE=0.001     # amplitude [m]
# export BALFEM_XWM=40          # wavemaker x [m]
# export BALFEM_SPONGE=40       # x-end sponge width [m]
# export BALFEM_MUMAX=5         # sponge strength
# export BALFEM_PERIODS=50      # run length in wave periods
# export BALFEM_NX=2000; export BALFEM_NY=100

# --- Solver knobs (RungeKutta :SDIRK_2_2 defaults; bump if Newton stalls) ----
# export BALFEM_LS_MAXITER=4000
# export BALFEM_NL_TOL=1e-6

balfem_run 84 examples/distributed/run_periodic_plane_wave_dist.jl
