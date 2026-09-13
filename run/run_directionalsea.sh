#!/bin/bash

#SBATCH --job-name="BALFEM_directionalsea"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
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
export BALFEM_PX=7
export BALFEM_PY=6            # 7*6 = 42 ranks

# --- Sea state (WaveSpec JONSWAP x cosine-power spreading — short-crested) ---
#     directional seas force y_wall_bc=:open inside the example (lateral sponges)
# export BALFEM_HS=0.002        # significant wave height [m]
# export BALFEM_TP=1.6          # peak period [s]
# export BALFEM_GAMMA=3.3       # JONSWAP peakedness
# export BALFEM_NFREQ=21        # frequency samples (bins = nf-1)
# export BALFEM_NTHETA=7        # angle samples (bins = n-1)
# export BALFEM_SPREAD_STD=20   # spreading sigma_theta [deg]
# export BALFEM_THETA_MAX=60    # angular truncation +/- [deg]
# export BALFEM_SEED=20260723   # phase seed (reproducible, rank-deterministic)

# --- Inflow relaxation zone (generation/absorption) -------------------------
#     ON for directional runs: the oblique :model boundary data is not an exact
#     discrete 2D eigenmode, so a bare clamped Dirichlet inflow reflects the
#     residual mismatch and it resonates into an exponential instability. The
#     relaxation zone absorbs that mismatch at the inflow (relaxes toward the
#     incident field). Width 0 -> one peak wavelength.
export BALFEM_RELAX=1
export BALFEM_RELAX_W=8       # relaxation-zone width [m] (~2 peak wavelengths); 0 -> 1 peak λ

# --- Domain / run length ----------------------------------------------------
# export BALFEM_LX=400; export BALFEM_LY=200  # domain size [m] (Ly wide: oblique paths)
# export BALFEM_D=3.5           # still-water depth [m]
# export BALFEM_SPONGE=40       # right + lateral sponge width [m]
# export BALFEM_MUMAX=5         # sponge strength
# export BALFEM_PERIODS=100     # run length in Tp units
# export BALFEM_WRITE_W=0; export BALFEM_WRITE_PRESSURE=0   # field output OFF (large runs)

# --- Solver knobs (RungeKutta :SDIRK_2_2 defaults; bump if Newton stalls) ----
# export BALFEM_LS_MAXITER=4000
# export BALFEM_NL_TOL=1e-6

balfem_run 42 examples/distributed/run_directional_sea_dist.jl
