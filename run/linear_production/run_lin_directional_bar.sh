#!/bin/bash
# PRODUCTION — linear short-crested directional sea (JONSWAP x cos-type spreading,
#   BC-generated), VARIABLE bed (shore-parallel submerged bar), Hs=0.2, Tp=2.0,
#   7 directions, std 15 deg, |theta| <= 40 deg
#   domain wide_small: 50 x 40 m, lateral sponges 2 m, right sponge 12 m (y open)
#   bar: d=3.5, height 1.5 (crest depth 2.0), centre x=26, half-width 6, tanh ramp 2 m
#   mesh 200 x 160 (square 0.25 m cells), Q2/Q1, CN
#   50 Tp = 100 s = 5000 steps
#   Geometry and diagnostics notes: see run_lin_directional_flat.sh (its twin; the
#   bed is the single axis between them).
#SBATCH --job-name="BALFEM_lin_directional_bar"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=42
#SBATCH --ntasks-per-node=42
#SBATCH --cpus-per-task=1
# MEMORY: --mem only, never --mem-per-cpu (run/SNELLIUS_ROME_LAUNCH_CONFIGS.md §2, §4).
#   168 GiB / 42 ranks = 4.00 GiB per rank, billed 6/8 of ONE rome node (96 CPUs).
#   ONE NODE, NOT v1's 2 x 28: same per-rank budget, 6/8 billed instead of 8/8, no
#   inter-node traffic, and it matches the layout the numbers below were measured on.
#   v1's measured wide run was 1 node x 40 ranks (8x5) on this exact 200 x 160 mesh:
#   800 cells/rank, peak 2.84 GB (max over ranks). Here 42 ranks give 762 cells/rank,
#   i.e. ~2.8 GB = 70 % of the 4 GiB budget.
#   Flat twin, v1 wall at 40 ranks: 73.5 s/step under SDIRK_2_2; the bar adds ~10 %
#   GMRES work -> ~112 h at the worst case (tight against the 119:59 wall), ~55-65 h
#   expected under CN.
#   If it is OOM-killed (exit 137): 1 node x 56 ranks, --exclusive --mem=0 (8/8, the
#   whole node), PX=8 PY=7.
#SBATCH --mem=168G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source $HOME/GridapBALFEM.jl/run/balfem_env.sh
source "$BALFEM_PROJ/run/linear_production/_linear_common.sh"

# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
export BALFEM_PX=7
export BALFEM_PY=6            # 7*6 = 42 ranks; 7.14 x 6.67 m subdomains (rule 22)

# --- Case ---------------------------------------------------------------------
export BALFEM_FLAT_BED=0      # variable bathymetry => the driver builds the shore-parallel bar
export BALFEM_HBAR=1.5 BALFEM_XBAR=26.0 BALFEM_WBAR=6.0
export BALFEM_HS=0.2 BALFEM_TP=2.0
export BALFEM_NFREQ=15 BALFEM_NTHETA=7
export BALFEM_SPREAD_STD=15 BALFEM_THETA_MAX=40
export BALFEM_SEED=20260723   # fixed sea-state seed: the flat and bar runs see the same sea
# wide_small geometry
export BALFEM_LX=50.0 BALFEM_LY=40.0
export BALFEM_NX=200  BALFEM_NY=160
export BALFEM_SPONGE_Y=2.0
export BALFEM_PERIODS=50
export BALFEM_SAVE_EVERY=20   # 0.4 s snapshots; halves VTK I/O on the largest mesh

balfem_run 42 examples/distributed_small/run_directional_sea_small.jl
