#!/bin/bash
# ============================================================================
#  submit.sh — the ONE submission script of the time-integrator stability study (Snellius rome)
#
#    bash run/integrator_study/submit.sh P              # the pilot pack (do this first)
#    bash run/integrator_study/submit.sh A              # a tier: P A B C D E F G H (see cases.py)
#    bash run/integrator_study/submit.sh A C E          # several tiers
#    bash run/integrator_study/submit.sh status         # per-case state from the output tree
#    bash run/integrator_study/submit.sh requeue-quota  # make the quota-killed cases resumable
#    bash run/integrator_study/submit.sh resume [T…]    # every pack (of tiers T…, default all) with
#                                                       #   a case neither finished nor failed
#    BALFEM_CASE_FILTER='_T160_.*_(crz34|sdirk43)' bash run/integrator_study/submit.sh resume
#                                                       # only the cases whose name matches (regex)
#
#  One table (cases.tsv, packs.tsv; 2026-11-23 the main matrix and the TR-BDF2 / generalised-α
#  "extra" set were joined: submit_extra.sh is gone). Each submission is ONE array job whose indices
#  are pack ids of packs.tsv; at most BALFEM_MAX_CONCURRENT packs (default 8) run at once.
#
#  OUTPUT on the project storage, never in $HOME (2026-11-21: the home quota killed 161 cases):
#    $BALFEM_STUDY_OUT/<case>/seg<NN>/   default $BALFEM_OUTPUT_ROOT/integrator_study
#    SLURM logs: $BALFEM_OUTPUT_ROOT/slurm_logs/integrator_study
#  with BALFEM_OUTPUT_ROOT = $STORAGE/GridapBALFEM_output (run/storage.sh; STORAGE exported in
#  ~/.bashrc, Snellius: /gpfs/admin/projectlinks/prjs1960; output/ of the repo is its link).
#  Volume: ≈ 200 MB per case (VTK with sub-cell samples, 8 per period) → ≈ 55 GB for the study.
#  Before the first submission: build the sysimage (compile/README.md; this study uses only
#  P1LFE-2, Nσ = 3, the default BALFEM_WARMUP_VERT).
#  CLUSTER: BALFEM_CLUSTER = snellius (default) | blue (DelftBlue); sbatch flags in cluster.sh.
#
#  RESUMING. A case directory is resumed by run_pack.sh unless it has FINISHED or FAILED. The quota
#  failures of 2026-11-21 left FAILED markers that are not results; `requeue-quota` renames exactly
#  those (message LightXML XMLWriteError or "Disk quota exceeded") to FAILED.quota, so `resume`
#  picks the cases up from their newest checkpoint (new segment seg02). Blow-ups (Newton) keep
#  their FAILED marker. Generalised-α resumes O(Δt²) off an uninterrupted run (cases.py).
# ============================================================================
set -euo pipefail
: "${BALFEM_PROJ:=$HOME/GridapBALFEM.jl}"
here="$BALFEM_PROJ/run/integrator_study"
source "$BALFEM_PROJ/run/storage.sh"
balfem_storage_check
source "$here/cluster.sh"
: "${BALFEM_STUDY_OUT:=$BALFEM_OUTPUT_ROOT/integrator_study}"
: "${BALFEM_MAX_CONCURRENT:=8}"
logs="$BALFEM_OUTPUT_ROOT/slurm_logs/integrator_study"
mkdir -p "$logs" "$BALFEM_STUDY_OUT"
cases_tsv="$here/cases.tsv"; packs_tsv="$here/packs.tsv"

[ $# -ge 1 ] || { sed -n '4,13p' "$0"; exit 1; }

#  state of one case: finished | failed | quota | partial (has segments) | new
case_state() {
    local d="$BALFEM_STUDY_OUT/$1"
    if   [ -f "$d/FINISHED" ]; then echo finished
    elif [ -f "$d/FAILED" ]; then
        if grep -q -E 'XMLWriteError|Disk quota exceeded' "$d/FAILED"; then echo quota; else echo failed; fi
    elif ls -d "$d"/seg[0-9][0-9] >/dev/null 2>&1; then echo partial
    else echo new; fi
}

submit_ids() {   # <label> <comma-separated pack ids>
    local label="$1" ids="$2"
    [ -n "$ids" ] || { echo "$label: nothing to submit"; return; }
    echo -n "$label, packs $ids: "
    ( cd "$logs" && sbatch --array="${ids}%${BALFEM_MAX_CONCURRENT}" --job-name="BALFEM_integ_$label" \
           ${STUDY_SBATCH_FLAGS[@]+"${STUDY_SBATCH_FLAGS[@]}"} \
           --export=ALL,BALFEM_STUDY_OUT="$BALFEM_STUDY_OUT",BALFEM_CLUSTER="${BALFEM_CLUSTER:-snellius}" \
           "$here/run_pack.sh" )
}

case "$1" in
status)
    printf "%-40s %-5s %-9s %s\n" case segs state detail
    tail -n +2 "$cases_tsv" | cut -f1 | while read -r name; do
        d="$BALFEM_STUDY_OUT/$name"
        nseg=$( (ls -d "$d"/seg[0-9][0-9] 2>/dev/null || true) | wc -l)
        st=$(case_state "$name"); det=""
        if [ "$st" = failed ] || [ "$st" = quota ]; then
            det=$( (grep -m1 '^ERROR' "$d/FAILED" || true) | cut -c1-70)
        fi
        printf "%-40s %-5s %-9s %s\n" "$name" "$nseg" "$st" "$det"
    done
    echo "--- totals:"; tail -n +2 "$cases_tsv" | cut -f1 | while read -r n; do case_state "$n"; done | sort | uniq -c
    ;;
requeue-quota)
    tail -n +2 "$cases_tsv" | cut -f1 | while read -r name; do
        [ "$(case_state "$name")" = quota ] || continue
        mv "$BALFEM_STUDY_OUT/$name/FAILED" "$BALFEM_STUDY_OUT/$name/FAILED.quota"
        echo "  $name: FAILED → FAILED.quota (resumable)"
    done
    echo "done; now: bash $0 resume"
    ;;
resume)
    shift; tiers="${*:-PABCDEFGH}"; tiers="${tiers// /}"
    ids=""; nq=0
    while IFS=$'\t' read -r pid tier names; do
        [[ "$tiers" == *"$tier"* ]] || continue
        IFS=',' read -r -a arr <<< "$names"
        for n in "${arr[@]}"; do
            [ -n "${BALFEM_CASE_FILTER:-}" ] && ! grep -qE "$BALFEM_CASE_FILTER" <<< "$n" && continue
            s=$(case_state "$n")
            if [ "$s" = partial ] || [ "$s" = new ]; then ids="${ids:+$ids,}$pid"; break; fi
            [ "$s" = quota ] && nq=$((nq + 1))
        done
    done < <(tail -n +2 "$packs_tsv")
    [ "$nq" -eq 0 ] || { echo "quota-killed cases still carry their FAILED marker: run '$0 requeue-quota' first" >&2; exit 1; }
    submit_ids resume "$ids"
    ;;
*)
    for tier in "$@"; do
        ids=$(awk -F'\t' -v t="$tier" 'NR > 1 && $2 == t {print $1}' "$packs_tsv" | paste -sd, -)
        [ -n "$ids" ] || { echo "tier $tier: no packs in packs.tsv" >&2; continue; }
        submit_ids "$tier" "$ids"
    done
    ;;
esac
