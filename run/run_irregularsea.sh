#!/bin/bash

#SBATCH --job-name="BALFEM_irregularsea"
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

# --- Sea state (WaveSpec JONSWAP, seeded phases — long-crested) --------------
# export BALFEM_HS=0.002        # significant wave height [m] (linear regime)
# export BALFEM_TP=1.6          # peak period [s]
# export BALFEM_GAMMA=3.3       # JONSWAP peakedness
# export BALFEM_NFREQ=21        # frequency samples (bins = nf-1)
# export BALFEM_SEED=20260723   # phase seed (reproducible, rank-deterministic)
# export BALFEM_FMIN_FAC=2.5    # fmin = 1/(FMIN_FAC*Tp)
# export BALFEM_FMAX_FAC=0.75   # fmax = 1/(FMAX_FAC*Tp); keep kd(fmax) <= kd_app
# export BALFEM_TRAMP=          # Hann ramp [s] (unset -> 2*Tp)

# --- Domain / run length ----------------------------------------------------
# export BALFEM_LX=400; export BALFEM_LY=20   # domain size [m]
# export BALFEM_D=3.5           # still-water depth [m]  (kd_p~3.6 at Tp=1.6)
# export BALFEM_SPONGE=40       # right-sponge width [m] (>= 4 peak wavelengths)
# export BALFEM_MUMAX=5         # sponge strength
# export BALFEM_PERIODS=200     # run length in Tp units
# export BALFEM_WRITE_W=0; export BALFEM_WRITE_PRESSURE=0   # field output OFF (large runs)

# --- Solver knobs (RungeKutta :SDIRK_2_2 defaults; bump if Newton stalls) ----
# export BALFEM_LS_MAXITER=4000
# export BALFEM_NL_TOL=1e-6

balfem_run 84 examples/distributed/run_irregular_sea_dist.jl
