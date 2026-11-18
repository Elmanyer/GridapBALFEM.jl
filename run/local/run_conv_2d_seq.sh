#!/bin/bash
# Suite C — 2-D, SEQUENTIAL. Q4 fine level dominates the cost; keep NX0 modest.
set -uo pipefail; cd "$(dirname "${BASH_SOURCE[0]}")/../.."
source run/local/babalfem_local.sh; mkdir -p output/local/logs
export BABALFEM_CONV_PU="${PU:-2,3,4}" BABALFEM_CONV_DOMAIN=d2 BABALFEM_CONV_MODE="${MODE:-static}"
export BABALFEM_CONV_LEVELS="${LEVELS:-3}" BABALFEM_CONV_NX0="${NX0:-6}" BABALFEM_CONV_NY0="${NY0:-4}"
export BABALFEM_CONV_DIST=0 BABALFEM_CONV_DT="${DT:-1e-5}" BABALFEM_CONV_NSTEPS="${NSTEPS:-100}"
export BABALFEM_CONV_OUT=output/local/mms_conv
babalfem_local_run examples/local_mms/run_mms_matrix.jl 2>&1 | tee output/local/logs/conv_2d_seq.log
