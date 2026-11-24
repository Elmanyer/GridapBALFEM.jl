#!/bin/bash
# PRODUCTION — linear BC-generated plane wave, VARIABLE bed (submerged bar), A=0.1, T=2.0
#   domain 50 x 20 m (y-periodic), mesh 200 x 40 (dx 0.25, dy 0.5), Q2/Q1, CN
#   bar: d=3.5, height 1.5 (crest depth 2.0), centre x=26, half-width 6, tanh ramp 2 m
#   40 periods = 80 s = 4000 steps (settles at 15 T; 25 T of steady record)
#SBATCH --job-name="BALFEM_lin_bcplane_bar"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
# MEMORY: --mem only, never --mem-per-cpu (run/SNELLIUS_ROME_LAUNCH_CONFIGS.md §2, §4).
#   112 GiB / 28 ranks = 4.00 GiB per rank, billed 4/8 of a rome node (64 CPUs).
#   v1 measured peak for this exact case: 2.5 GB/rank. Wall: 27 s/step under
#   SDIRK_2_2 -> ~30 h at the worst case, ~17 h expected under CN.
#SBATCH --mem=112G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source $HOME/GridapBALFEM.jl/run/balfem_env.sh
source "$BALFEM_PROJ/run/linear_production/_linear_common.sh"

# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
export BALFEM_PX=7
export BALFEM_PY=4            # 7*4 = 28 ranks

# --- Case ---------------------------------------------------------------------
export BALFEM_FLAT_BED=0      # variable bathymetry => the driver builds the y-invariant bar
export BALFEM_HBAR=1.5 BALFEM_XBAR=26.0 BALFEM_WBAR=6.0
export BALFEM_AWAVE=0.1       # linear: amplitude only scales the solution; 0.1 matches the nonlinear twins
export BALFEM_TWAVE=2.0
export BALFEM_LX=50.0 BALFEM_LY=20.0
export BALFEM_NX=200  BALFEM_NY=40
export BALFEM_PERIODS=40
export BALFEM_SAVE_EVERY=10   # 0.2 s snapshots = 10 per period

balfem_run 28 examples/distributed_small/run_bc_plane_small.jl
