#!/bin/bash

#  SNELLIUS ROME SIZING: 42 ranks x 4 GiB = 168 GiB = exactly 6/8 of a node.
#  Was 48 ranks, which billed 7/8 (memory fraction 6.86/8 rounds up). The rome node is
#  128 cores / 224 GiB divided into EIGHTHS, and the LARGER of the core/memory fraction
#  sets the bill -- at 4 GiB/rank memory always wins, so the cost-optimal count is 7k
#  ranks for tier k/8. See run/SNELLIUS_ROME_LAUNCH_CONFIGS.md.
#SBATCH --job-name="BALFEM_irregularsea"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --nodes=2
#SBATCH --ntasks=84
#SBATCH --ntasks-per-node=42
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
# Cut to 48 ranks/node x 4 GB = 192 GB/node (128 -> 96 ranks over 2 node(s)),
# rather than spreading the old rank count over more nodes.
# Do NOT drop mem-per-cpu on the argument that the sysimage removes the
# per-rank compile: that argument was tested and refuted.
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

balfem_run 96 examples/distributed/run_irregular_sea_dist.jl
