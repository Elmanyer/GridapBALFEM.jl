#!/bin/bash
# PRODUCTION — linear long-crested irregular sea (JONSWAP, BC-generated), FLAT bed, Hs=0.2, Tp=2.0
#   domain 50 x 20 m (y-walls, normal incidence), mesh 200 x 40 (dx 0.25, dy 0.5), Q2/Q1, CN
#   60 Tp = 120 s = 6000 steps (settles at ~22 Tp; ~38 Tp of statistically steady record)
#SBATCH --job-name="BALFEM_lin_irregular_flat"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
# MEMORY: --mem only, never --mem-per-cpu (run/SNELLIUS_ROME_LAUNCH_CONFIGS.md §2, §4).
#   112 GiB / 28 ranks = 4.00 GiB per rank, billed 4/8 of a rome node (64 CPUs).
#   Same mesh and rank grid as the bc-plane cases (v1 peak 2.6 GB/rank, 24-27 s/step
#   under SDIRK_2_2) -> ~45 h at the worst case, ~25 h expected under CN.
#SBATCH --mem=112G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source $HOME/GridapBALFEM.jl/run/balfem_env.sh
source "$BALFEM_PROJ/run/linear_production/_linear_common.sh"

# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
export BALFEM_PX=7
export BALFEM_PY=4            # 7*4 = 28 ranks

# --- Case ---------------------------------------------------------------------
export BALFEM_FLAT_BED=1
export BALFEM_HS=0.2 BALFEM_TP=2.0
export BALFEM_NFREQ=15
export BALFEM_SEED=20260723   # fixed sea-state seed: the flat and bar runs see the same sea
export BALFEM_LX=50.0 BALFEM_LY=20.0
export BALFEM_NX=200  BALFEM_NY=40
export BALFEM_PERIODS=60
export BALFEM_SAVE_EVERY=10   # 0.2 s snapshots

balfem_run 28 examples/distributed_small/run_irregular_sea_small.jl
