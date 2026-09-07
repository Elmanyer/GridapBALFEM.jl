#!/bin/bash

#SBATCH --job-name="BALFEM_ichump"
#SBATCH --partition=rome
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=14
#SBATCH --ntasks-per-node=14
#SBATCH --cpus-per-task=1
# 4 GiB/rank is REQUIRED (measured peak 3633 MB/rank; 2 GB/core was OOM-killed).
# Sizing: 14 ranks/node x 4 GiB = 56 GiB/node = 2/8 of a rome node (14 total).
# Rules and tiers: run/SNELLIUS_ROME_LAUNCH_CONFIGS.md
#SBATCH --mem-per-cpu=4G
#SBATCH --output=GridapBALFEM.%j.out
#SBATCH --error=GridapBALFEM.%j.err

source $HOME/GridapBALFEM.jl/run/balfem_env.sh

# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
export BALFEM_PX=7
export BALFEM_PY=2            # 7*2 = 14 ranks

# --- Case-specific knobs (Gaussian hump released from rest, CLOSED basin) ----
#     (x_wall_bc=true is forced inside the example — do not disable it)
# export BALFEM_L=100           # basin side [m]
# export BALFEM_D=3.5           # depth [m]
# export BALFEM_A0=0.01         # hump amplitude [m]
# export BALFEM_SIGMA0=4.0      # hump half-width [m]
# export BALFEM_TFINAL=40       # final time [s]
# export BALFEM_NX=100; export BALFEM_NY=100

# --- Solver knobs (RungeKutta :SDIRK_2_2 defaults; bump if Newton stalls) ----
# export BALFEM_LS_MAXITER=4000
# export BALFEM_NL_TOL=1e-6

balfem_run 14 examples/distributed/run_ic_hump_dist.jl
