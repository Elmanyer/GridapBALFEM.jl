#!/bin/bash

#SBATCH --job-name="BALFEM_bathymetry"
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
export BALFEM_PX=14
export BALFEM_PY=6            # 14*6 = 84 ranks (2 nodes)

# --- Case-specific knobs (plane wave shoaling over a tanh submerged bar) -----
# export BALFEM_LX=300; export BALFEM_LY=20   # domain size [m]
# export BALFEM_D0=3.5          # offshore depth [m]
# export BALFEM_HBAR=2.0        # bar height [m]
# export BALFEM_XBAR=150        # bar centre x [m]
# export BALFEM_WBAR=25         # bar half-width [m]
# export BALFEM_TWAVE=2.5       # wave period [s]
# export BALFEM_AWAVE=0.001     # amplitude [m]
# export BALFEM_XWM=40          # wavemaker x [m]
# export BALFEM_SPONGE=35       # sponge width [m]
# export BALFEM_PERIODS=40      # run length in wave periods
# export BALFEM_NX=1500; export BALFEM_NY=100

# --- Sea-bed geometry (variable bed → ∇h terms ON) --------------------------
# export BALFEM_FLAT_BED=0         # 0 = variable bathymetry (default for the bar); 1 = flat bed (∇h≡0)
# export BALFEM_NL_PRESSURE=native # nonlinear pressure: none | native | full

# --- Solver knobs (RungeKutta :SDIRK_2_2 defaults; bump if Newton stalls) ----
# export BALFEM_LS_MAXITER=4000
# export BALFEM_NL_TOL=1e-6

balfem_run 84 examples/distributed/run_bathymetry_dist.jl
