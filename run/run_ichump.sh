#!/bin/bash

#  SNELLIUS ROME SIZING: 14 ranks x 4 GiB = 56 GiB = exactly 2/8 of a node.
#  Was 16 ranks, which billed 3/8 (memory fraction 2.29/8 rounds up). The rome node is
#  128 cores / 224 GiB divided into EIGHTHS, and the LARGER of the core/memory fraction
#  sets the bill -- at 4 GiB/rank memory always wins, so the cost-optimal count is 7k
#  ranks for tier k/8. See run/SNELLIUS_ROME_LAUNCH_CONFIGS.md.
#SBATCH --job-name="BALFEM_ichump"
#SBATCH --partition=rome
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=14
#SBATCH --ntasks-per-node=14
#SBATCH --cpus-per-task=1
# Memory: the rome node default (2 GB/core) was TESTED (2026-08) and the job was
# still OOM-killed. This request is therefore NOT provisional headroom for a
# compile spike -- it is required until the consumption is attributed; see
# building_files/OPEN_ITEMS.md section 1. Sized so it FITS a rome
# node (256 GB / 128 cores): 16 ranks x 4 GB over 1 node(s) = 64 GB/node.
# Do NOT drop this on the argument that the sysimage removes the per-rank
# compile: that argument was tested and refuted.
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

balfem_run 16 examples/distributed/run_ic_hump_dist.jl
