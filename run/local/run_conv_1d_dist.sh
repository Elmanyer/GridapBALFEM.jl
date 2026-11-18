#!/bin/bash
# Suite B — quasi-1D, DISTRIBUTED (4 ranks). GMRES+Jacobi at ls_rtol=1e-13 so the
# linear solve cannot cap the spatial error; verify with run_conv_tolerance_gate.sh.
set -uo pipefail; cd "$(dirname "${BASH_SOURCE[0]}")/../.."
source run/local/babalfem_local.sh; mkdir -p output/local/logs
export BABALFEM_CONV_PU="${PU:-2,3,4}" BABALFEM_CONV_DOMAIN=d1 BABALFEM_CONV_MODE="${MODE:-both}"
export BABALFEM_CONV_LEVELS="${LEVELS:-4}" BABALFEM_CONV_NX0="${NX0:-16}" BABALFEM_CONV_DIST=1
export BABALFEM_CONV_PX="${PX:-4}" BABALFEM_CONV_PY="${PY:-1}"
export BABALFEM_CONV_LSRTOL="${LSRTOL:-1e-13}" BABALFEM_CONV_NLTOL="${NLTOL:-1e-13}"
export BABALFEM_CONV_DT="${DT:-1e-5}" BABALFEM_CONV_NSTEPS="${NSTEPS:-100}"
export BABALFEM_CONV_OUT=output/local/mms_conv
babalfem_local_mpi "${RANKS:-4}" examples/local_mms/run_mms_matrix.jl 2>&1 | tee output/local/logs/conv_1d_dist.log
