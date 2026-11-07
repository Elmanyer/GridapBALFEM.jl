#!/bin/bash
# PRODUCTION RERUN — linear short-crested directional sea (JONSWAP x cos-type spreading,
#   BC-generated), FLAT bed, Hs=0.2, Tp=2.0, 7 directions, std 15 deg, |theta| <= 40 deg
#   domain wide_small: 50 x 40 m, lateral sponges 2 m, right sponge 12 m (y open)
#   mesh 200 x 160 (square 0.25 m cells), Q2/Q1, CN
#   50 Tp = 100 s = 5000 steps
#   Identical to run_lin_directional_flat.sh except the preconditioner (:schwarz) and the
#   checkpoints. Why the rerun: markdown_files/LINEAR_PRODUCTION_RESULTS.md (the Jacobi run,
#   job 28614763, was killed by the wall at t = 69.2 s).
#   This flat run is also the reference for run_lin_directional_bar_schwarz.sh and
#   run_lin_directional_bigbar.sh (same sea: BALFEM_SEED).
#   PILOT FIRST (the :schwarz cost at 42 ranks is unmeasured; ~1 h, 100 steps):
#       BALFEM_TFINAL=2.0 sbatch --time=02:00:00 run/linear_production/run_lin_directional_flat_schwarz.sh
#   then read s/step and gmres= in the .out before submitting the full run.
#SBATCH --job-name="BALFEM_lin_directional_flat_schwarz"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=42
#SBATCH --ntasks-per-node=42
#SBATCH --cpus-per-task=1
# MEMORY: --mem only, never --mem-per-cpu (run/SNELLIUS_ROME_LAUNCH_CONFIGS.md §2, §4).
#   168 GiB / 42 ranks = 4.00 GiB per rank, billed 6/8 of ONE rome node. The Jacobi run peaked
#   at 2.77 GB/rank; the own-block LU adds one sparse factorisation of ~762 cells' DOFs per rank
#   (tens of MB). If it is OOM-killed (exit 137): 1 node x 56 ranks, --exclusive --mem=0, PX=8 PY=7
#   (a checkpoint cannot be resumed on a different rank grid).
#SBATCH --mem=168G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source $HOME/GridapBALFEM.jl/run/balfem_env.sh
source "$BALFEM_PROJ/run/linear_production/_linear_common.sh"

# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
export BALFEM_PX=7
export BALFEM_PY=6            # 7*6 = 42 ranks; 7.14 x 6.67 m subdomains (rule 22)

# --- Linear solver: own-block LU preconditioner (src/timeloop_dist.jl, OwnBlockLUSolver) ---
#  The Jacobi run needed ~930 GMRES iterations per solve (3 solves/step, 125-137 s/step): Jacobi
#  cannot remove the R_P grad-div conditioning, kappa ~ (H p / h_min)^2, and this mesh has the
#  smallest h_min of the suite. :schwarz is an exact sparse LU of each rank's own block,
#  refreshed every Newton iteration (15-19 GMRES iterations on the 4-rank model-5 flume; no
#  coarse space, so expect more at 42 ranks). It changes nothing above the ls_rtol level
#  (rule 24). Override at submit time: BALFEM_PRECOND=jacobi sbatch <this file>
export BALFEM_PRECOND=${BALFEM_PRECOND:-schwarz}

# --- Checkpoints (src/checkpoint.jl) -------------------------------------------
#  u_n every 500 steps (10 s simulated), the 2 newest kept, in <output_dir>/checkpoints/.
#  If the wall kills the job, resume on the SAME rank grid (T_final stays absolute):
#      BALFEM_RESTART_FROM=<output_dir of the killed run> sbatch <this file>
#  The resumed run writes to a new directory carrying a `resume` token.
export BALFEM_CHECKPOINT_EVERY=${BALFEM_CHECKPOINT_EVERY:-500}
export BALFEM_CHECKPOINT_KEEP=${BALFEM_CHECKPOINT_KEEP:-2}

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
export BALFEM_NAME_EXTRA=schwarz

balfem_run 42 examples/distributed_small/run_directional_sea_small.jl
