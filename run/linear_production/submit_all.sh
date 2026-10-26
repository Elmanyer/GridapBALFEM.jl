#!/bin/bash
# ============================================================================
#  submit_all.sh — submit the linear production campaign on Snellius
#
#    bash run/linear_production/submit_all.sh            # all seven cases
#    bash run/linear_production/submit_all.sh bcplane    # only launchers matching a pattern
#
#  The SLURM .out/.err files land in $BALFEM_PROJ/output/slurm_logs/linear_production
#  (the job's working directory is the submit directory, and %x.%j.out is relative).
#  Build the sysimage FIRST: every launcher sets BALFEM_STRICT_SYSIMAGE=1 and
#  refuses a stale one.
# ============================================================================
set -euo pipefail

: "${BALFEM_PROJ:=$HOME/GridapBALFEM.jl}"
here="$BALFEM_PROJ/run/linear_production"
logs="$BALFEM_PROJ/output/slurm_logs/linear_production"
mkdir -p "$logs"

pattern="${1:-}"
cd "$logs"
for f in "$here"/run_lin_*.sh; do
    [[ -z "$pattern" || "$(basename "$f")" == *"$pattern"* ]] || continue
    echo -n "$(basename "$f"): "
    sbatch "$f"
done
