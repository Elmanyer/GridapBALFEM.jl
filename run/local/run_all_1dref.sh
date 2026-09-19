#!/bin/bash
# ==============================================================
#  run_all_1dref.sh — the :full flat dx x dt x Jacobian refinement factorial
#
#  8 runs: nx {240,480} x dt {0.04,0.02} x Jacobian {hand, exact AD}, 16 periods.
#  Rule 38b is the discriminator: refinement DELAYING onset means
#  under-resolution, refinement ADVANCING it means a grid-scale problem.
#
#  LAUNCH ORDER IS BY COST, CHEAPEST FIRST (SJF), because the deliverable here is
#  COVERAGE of the 2x2, not makespan (rule 43). The four hand-Jacobian runs are
#  ~10x cheaper per step and answer the question on their own; the AD arm is the
#  control that rules out the incomplete jacobian_u, and it is expected to agree.
#
#  ⚠ THREADS AND MEMORY. 2 BLAS threads x 8 = 16, this box's core count. Measured
#  RSS: ~2 GB for a hand run, ~3.9 GB for an AD run at nx=240 -- and rule 41 says
#  they GROW. Eight at once could reach ~24 GB of 30. The runs are therefore
#  started in two waves: the four hand runs first, then the AD arm as slots free.
#  Do NOT raise this to 3 threads without watching `free -g`.
# ==============================================================
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
logdir="$root/output/local_1d/_logs"; mkdir -p "$logdir"
stamp="$(date '+%Y%m%d_%H%M')"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-2}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"

launch () {
    local case="$1"
    local log="$logdir/ref_${case}_${stamp}.log"
    nohup stdbuf -oL -eL bash "$here/run_1dref_${case}.sh" > "$log" 2>&1 &
    echo "  launched $case  pid $!  -> $(basename "$log")"
    sleep 30          # stagger: eight simultaneous cold starts contend for the
                      # precompile cache while none of them is solving yet
}

echo "WAVE 1 — hand (pre-derived) Jacobian, cheapest first"
for c in full_nx240_dt004_hand full_nx240_dt002_hand full_nx480_dt004_hand full_nx480_dt002_hand; do
    launch "$c"
done

echo "WAVE 2 — exact AD Jacobian (~10x per step; the control arm)"
for c in full_nx240_dt004_ad full_nx240_dt002_ad full_nx480_dt004_ad full_nx480_dt002_ad; do
    launch "$c"
done
wait
