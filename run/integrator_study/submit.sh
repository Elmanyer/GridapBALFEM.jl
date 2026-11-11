#!/bin/bash
# ============================================================================
#  submit.sh — submit tiers of the time-integrator stability study on Snellius (rome)
#
#    bash run/integrator_study/submit.sh P            # the pilot pack (do this first)
#    bash run/integrator_study/submit.sh A            # a tier: P A B C D E (see cases.py)
#    bash run/integrator_study/submit.sh A C E        # several tiers
#    bash run/integrator_study/submit.sh status       # per-case state from the output tree
#
#  Each tier becomes ONE array job whose indices are its pack ids in packs.tsv; at most
#  BALFEM_MAX_CONCURRENT packs (default 8 = one full rome node's worth) run at once.
#  Output: $BALFEM_STUDY_OUT/<case>/seg<NN>/ (default $BALFEM_PROJ/output/integrator_study);
#  SLURM logs: $BALFEM_PROJ/output/slurm_logs/integrator_study.
#  Before the first submission: build the sysimage (compile/README.md; this study uses only
#  P1LFE-2, Nσ = 3, the default BALFEM_WARMUP_VERT), and check the disk quota (≈ 160 MB per case).
# ============================================================================
set -euo pipefail
: "${BALFEM_PROJ:=$HOME/GridapBALFEM.jl}"
: "${BALFEM_STUDY_OUT:=$BALFEM_PROJ/output/integrator_study}"
: "${BALFEM_MAX_CONCURRENT:=8}"
here="$BALFEM_PROJ/run/integrator_study"
logs="$BALFEM_PROJ/output/slurm_logs/integrator_study"
mkdir -p "$logs" "$BALFEM_STUDY_OUT"

[ $# -ge 1 ] || { sed -n '4,9p' "$0"; exit 1; }

if [ "$1" = "status" ]; then
    printf "%-40s %-5s %s\n" case segs state
    tail -n +2 "$here/cases.tsv" | cut -f1,2 | while IFS=$'\t' read -r name tier; do
        d="$BALFEM_STUDY_OUT/$name"
        nseg=$( (ls -d "$d"/seg[0-9][0-9] 2>/dev/null || true) | wc -l)
        state="not started"
        [ "$nseg" -gt 0 ] && state="running or cut"
        [ -f "$d/FINISHED" ] && state="finished"
        [ -f "$d/FAILED" ] && state="ERROR: $(head -1 "$d/FAILED")"
        printf "%-40s %-5s %s\n" "$name" "$nseg" "$state"
    done
    exit 0
fi

cd "$logs"
for tier in "$@"; do
    ids=$(awk -F'\t' -v t="$tier" 'NR > 1 && $2 == t {print $1}' "$here/packs.tsv" | paste -sd, -)
    [ -n "$ids" ] || { echo "tier $tier: no packs (tiers: P A B C D E)" >&2; continue; }
    echo -n "tier $tier, packs $ids: "
    sbatch --array="${ids}%${BALFEM_MAX_CONCURRENT}" --job-name="BALFEM_integ_$tier" \
           --export=ALL,BALFEM_STUDY_OUT="$BALFEM_STUDY_OUT" "$here/run_pack.sh"
done
