#!/bin/bash
# PRODUCTION — linear short-crested directional sea (JONSWAP x cos-type spreading,
#   BC-generated), FLAT bed, Hs=0.2, Tp=2.0, 7 directions, std 15 deg, |theta| <= 40 deg
#   domain wide_small: 50 x 40 m, lateral sponges 2 m, right sponge 12 m (y open)
#   mesh 200 x 160 (square 0.25 m cells), Q2/Q1, CN
#   50 Tp = 100 s = 5000 steps (oblique paths are longer than the 22 Tp normal transit)
#
#   WIDE_SMALL IS THE DEFAULT HERE. On the 20 m-wide small domain the 8 m lateral
#   sponges covered 16 of 20 m of the Dirichlet inflow edge, and the undamped strip
#   (4 m = 0.45 lambda_eff) was too narrow for the transverse structure of the sea to
#   form. Ly=40 with 2 m laterals leaves 36 m = 4.0 lambda_eff. The geometry is spelled
#   out below although the driver defaults to it, so the launcher cannot inherit a
#   different one silently.
#   WATCH eta_max_damped/eta_max_int in diagnostics.csv: in the v1 linear wide run it
#   sat at 0.8-1.35, and the >1 peaks were at x ~ 2-3 m, i.e. inside the inflow
#   RELAXATION zone (counted as "damped"), not lateral reflection. If the maximum
#   moves into the lateral sponges, try BALFEM_MUMAX=16 (gentler), never stronger.
#SBATCH --job-name="BALFEM_lin_directional_flat"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --nodes=2
#SBATCH --ntasks=56
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
# MEMORY: --mem only, never --mem-per-cpu (run/SNELLIUS_ROME_LAUNCH_CONFIGS.md §2, §4).
#   --mem is PER NODE: 112 GiB / 28 ranks = 4.00 GiB per rank on each of 2 nodes,
#   billed 4/8 + 4/8 = 8/8 (128 CPUs). v1 linear wide run of this exact mesh and
#   grid: peak 2.84 GB/rank, 73.5 s/step under SDIRK_2_2 -> ~102 h at the worst
#   case, ~55-60 h expected under CN.
#   If it is OOM-killed (exit 137): --ntasks=84 --ntasks-per-node=42 --mem=168G,
#   PX=14 PY=6. Do NOT coarsen: dx 0.25 is 14.7 cells per the shortest wavelength.
#SBATCH --mem=112G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source $HOME/GridapBALFEM.jl/run/balfem_env.sh
source "$BALFEM_PROJ/run/linear_production/_linear_common.sh"

# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
export BALFEM_PX=8
export BALFEM_PY=7            # 8*7 = 56 ranks; 6.25 x 5.71 m subdomains (rule 22)

# --- Case ---------------------------------------------------------------------
export BALFEM_FLAT_BED=1
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

balfem_run 56 examples/distributed_small/run_directional_sea_small.jl
