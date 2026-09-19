#!/bin/bash
# ==============================================================
#  run_all_1dprod.sh — the four 1-D nonlinear production runs, side by side
#
#  nl_pressure {native, full} x bed {flat, bar}, Q2/Q1 Taylor-Hood, A = 0.10 m,
#  100 wave periods (160 s, 4000 steps), VTK every 0.2 s.
#
#  WHY FOUR AT ONCE AND NOT 12 RANKS EACH. A direct LU beats every MPI split on
#  this mesh (6.45 s/step at 1 rank against 13-19 s/step at 2-12 ranks: the GMRES
#  count is rank-independent, so the cost is the weak Jacobi preconditioner, not
#  the decomposition). Four sequential cases therefore finish the SET sooner than
#  four 12-rank cases in series -- and sequential keeps the point gauges, which
#  the distributed driver does not evaluate at all.
#
#  THREADS. OPENBLAS_NUM_THREADS=3 per process, 4 processes = 12 of this box's 16
#  cores, leaving 4 for the OS and analysis. Left at the default every process
#  would claim all 16 and they would fight.
#
#  MEMORY. ~1.7 GB per process at the start. ⚠ Rule 41: a long-lived Julia+Gridap
#  process GROWS -- the vertical-basis campaign measured 1.5 -> 3.9 GB over 14 h
#  and had to kill workers to stay off swap. Four processes at 4 GB is 16 GB of
#  this machine's 30 GB, so watch `free -g` rather than assuming.
#
#  A KILLED RUN IS STILL USABLE: diagnostics.csv is appended as it goes and the
#  .pvd index is rewritten after every snapshot, so ParaView opens the series
#  mid-run.
# ==============================================================
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
logdir="$root/output/local_1d/_logs"
mkdir -p "$logdir"
stamp="$(date '+%Y%m%d_%H%M')"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-3}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-3}"

for case in nl_native_flat nl_native_bar nl_full_flat nl_full_bar; do
    log="$logdir/prod_${case}_${stamp}.log"
    echo "launching $case -> $log"
    nohup stdbuf -oL -eL bash "$here/run_1dprod_${case}.sh" > "$log" 2>&1 &
    echo "  pid $!"
    #  Stagger: four simultaneous cold starts contend for the same precompile
    #  cache and the same 16 cores while none of them is solving yet.
    sleep 45
done
wait
