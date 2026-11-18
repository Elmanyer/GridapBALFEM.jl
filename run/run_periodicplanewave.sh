#!/bin/bash

#SBATCH --job-name="BABALFEM_periodicplanewave"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --nodes=2
#SBATCH --ntasks=128
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
# Memory: the rome node default (2 GB/core) was TESTED (2026-08) and the job was
# still OOM-killed. This request is therefore NOT provisional headroom for a
# compile spike -- it is required until the consumption is attributed; see
# building_files/LOCAL_VALIDATION_PLAN.md section 2.2. Sized so it FITS a rome
# node (256 GB / 128 cores): 128 ranks x 4 GB over 2 node(s) = 256 GB/node.
# Do NOT drop this on the argument that the sysimage removes the per-rank
# compile: that argument was tested and refuted.
#SBATCH --mem-per-cpu=4G
#SBATCH --output=GridapBABALFEM.%j.out
#SBATCH --error=GridapBABALFEM.%j.err

source $HOME/GridapBABALFEM.jl/run/babalfem_env.sh

# --- MPI process grid (PX*PY MUST equal the rank count given to babalfem_run) ---
export BABALFEM_PX=32
export BABALFEM_PY=4            # 32*4 = 128 ranks; mesh 2000x100

# --- Case-specific knobs (periodic-width flume; defaults shown) --------------
#     The y-edges are periodic (set in the example) -> no lateral sponges.
# export BABALFEM_LX=400; export BABALFEM_LY=20   # domain size [m]
# export BABALFEM_D=3.5           # still-water depth [m]
# export BABALFEM_TWAVE=1.6       # wave period [s]
# export BABALFEM_AWAVE=0.001     # amplitude [m]
# export BABALFEM_XWM=40          # wavemaker x [m]
# export BABALFEM_SPONGE=40       # x-end sponge width [m]
# export BABALFEM_MUMAX=5         # sponge strength
# export BABALFEM_PERIODS=50      # run length in wave periods
# export BABALFEM_NX=2000; export BABALFEM_NY=100

# --- Solver knobs (RungeKutta :SDIRK_2_2 defaults; bump if Newton stalls) ----
# export BABALFEM_LS_MAXITER=4000
# export BABALFEM_NL_TOL=1e-6

babalfem_run 128 examples/distributed/run_periodic_plane_wave_dist.jl
