#!/bin/bash
# PRODUCTION — BIG BAR — linear long-crested irregular sea (JONSWAP, BC-generated), VARIABLE bed, Hs=0.2, Tp=2.0
#   domain 50 x 20 m (y-walls, normal incidence), mesh 200 x 40 (dx 0.25, dy 0.5), Q2/Q1, CN
#   bar: d(x) = 3.5 - 0.5*2.85*[tanh((x-20)/3) - tanh((x-32)/3)]: true crest depth 0.753 m
#   (HBAR is nominal: the tanh does not saturate at SRAMP = 3), max slope 0.47, plateau
#   x = 24.5-27.5, bed within 1.5 % of 3.5 m outside x = 14-38 (clean stretch 6-14 m upwave).
#   WHY A BIGGER BAR: the 2.0 m crest is kd 2.1 at T = 2 s, nearly deep water: the first runs
#   showed only -3.5 % amplitude over the crest and 0.23 rad of phase (markdown_files/
#   LINEAR_PRODUCTION_RESULTS.md). Linear theory at the 0.75 m crest (kd 0.99 at T = 2 s):
#   wavelength -31 %, +3.2 rad of phase over 10 m, shoaling Ks 0.92. Ks bottoms out near
#   kd ~ 1, so the AMPLITUDE stays a weak signal; phase and wavelength are the strong ones.
#   Reference: run_lin_irregular_flat.sh (completed 2026-10; same forcing).
#   60 Tp = 120 s = 6000 steps (settles at ~22 Tp; ~38 Tp of statistically steady record)
#SBATCH --job-name="BALFEM_lin_irregular_bigbar"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
# MEMORY: --mem only, never --mem-per-cpu (run/SNELLIUS_ROME_LAUNCH_CONFIGS.md §2, §4).
#   112 GiB / 28 ranks = 4.00 GiB per rank, billed 4/8 of a rome node (64 CPUs).
#   Same mesh and rank grid as the bc-plane cases (v1 peak 2.5 GB/rank, 27 s/step
#   under SDIRK_2_2 on the bar) -> ~45 h at the worst case, ~27 h expected under CN.
#SBATCH --mem=112G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source $HOME/GridapBALFEM.jl/run/balfem_env.sh
source "$BALFEM_PROJ/run/linear_production/_linear_common.sh"

# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
export BALFEM_PX=7
export BALFEM_PY=4            # 7*4 = 28 ranks

# --- Linear solver: own-block LU preconditioner (src/timeloop_dist.jl, OwnBlockLUSolver) ---
#  GMRES + Jacobi took 515-597 iterations per solve on this mesh (3 solves/step). :schwarz is an
#  exact sparse LU of each rank's own block, refreshed every Newton iteration; it changes nothing
#  above the ls_rtol level (rule 24), so the flat reference run with Jacobi stays comparable.
#  Override at submit time: BALFEM_PRECOND=jacobi sbatch <this file>
export BALFEM_PRECOND=${BALFEM_PRECOND:-schwarz}

# --- Checkpoints (src/checkpoint.jl) -------------------------------------------
#  u_n every 500 steps (10 s simulated), the 2 newest kept, in <output_dir>/checkpoints/.
#  If the wall kills the job, resume on the SAME rank grid (T_final stays absolute):
#      BALFEM_RESTART_FROM=<output_dir of the killed run> sbatch <this file>
export BALFEM_CHECKPOINT_EVERY=${BALFEM_CHECKPOINT_EVERY:-500}
export BALFEM_CHECKPOINT_KEEP=${BALFEM_CHECKPOINT_KEEP:-2}

# --- Case ---------------------------------------------------------------------
export BALFEM_FLAT_BED=0      # variable bathymetry => the driver builds the y-invariant bar
export BALFEM_HBAR=2.85 BALFEM_XBAR=26.0 BALFEM_WBAR=6.0
export BALFEM_SRAMP=3.0      # gentler ramps than the default WBAR/3 = 2: max slope 0.47 instead of 0.71
export BALFEM_HS=0.2 BALFEM_TP=2.0
export BALFEM_NFREQ=15
export BALFEM_SEED=20260723   # fixed sea-state seed: the flat and bar runs see the same sea
export BALFEM_LX=50.0 BALFEM_LY=20.0
export BALFEM_NX=200  BALFEM_NY=40
export BALFEM_PERIODS=60
export BALFEM_SAVE_EVERY=10   # 0.2 s snapshots
export BALFEM_NAME_EXTRA=bigbar

balfem_run 28 examples/distributed_small/run_irregular_sea_small.jl
