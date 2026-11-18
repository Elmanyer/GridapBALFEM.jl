#!/bin/bash
# LOCAL small 2-D - nl_plane_flat
#   Nonlinear advection only.
#   Cluster sibling: run/dist_small/run_nl_none_periodic_plane_small.sh
# MPI on a 12-rank partition (4x3). NO point gauges - judge it from diagnostics.csv.
source "$(dirname "${BASH_SOURCE[0]}")/babalfem_local.sh"

export BABALFEM_WAVE_GEN=line
export BABALFEM_REGIME=nonlinear
export BABALFEM_NL_PRESSURE=none
export BABALFEM_FLAT_BED=1

export BABALFEM_PX=4
export BABALFEM_PY=3            # 4*3 = 12 ranks
export BABALFEM_LX=40.0
export BABALFEM_LY=15.0
export BABALFEM_NX=96           # 96/4 = 24 cells per rank in x
export BABALFEM_NY=36           # 36/3 = 12 cells per rank in y
export BABALFEM_PERIODS=10
export BABALFEM_SAVE_EVERY=10
export BABALFEM_DIAG_EVERY=5

babalfem_local_mpi 12 examples/local_2d/run_small_2d.jl
