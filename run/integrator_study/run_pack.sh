#!/bin/bash
# ============================================================================
#  run_pack.sh — one PACK of the time-integrator stability study on a rome 1/8 node
#
#  One SLURM array task = one pack = up to 8 sequential box cases run side by side
#  (examples/local_1d/run_periodic_1d.jl, one Julia process each, the sysimage loaded with -J).
#  The case matrix and the packs: run/integrator_study/{cases,packs}.tsv (cases.py).
#  Submit with submit.sh, never by hand (it maps tiers to array indices).
#
#  SEGMENTS. A case may need more than the 120 h wall. Each case therefore runs in segments
#  <case>/seg01, seg02, …: a segment is stopped by `timeout` before the wall, and the next one
#  resumes from the newest checkpoint (BALFEM_RESTART_FROM; checkpoints every 250 steps, so a
#  segment loses at most 250 steps). When a pack ends with a case cut by the timeout, the pack
#  resubmits itself (at most BALFEM_MAX_SEG segments). Per segment the launcher records
#  `exit_status`; a case is
#    finished  when its log has the driver's `done:` line (completed OR stopped by the divergence
#              guard: both are results),
#    failed    when it exited with an error (no retry: failures reproduce, rule 45),
#    cut       when the timeout stopped it (exit 124) — resumed by the next segment.
#  The analysis stitches the segments (each has its own diagnostics.csv and VTK series).
#
#  MEMORY: --mem only (run/SNELLIUS_ROME_LAUNCH_CONFIGS.md §2): 28 GiB = 1/8 node for 8 cases
#  (≈ 3.5 GiB each; a 𝓝 box case peaked at 2.8–3.0 GB locally). One thread per case.
# ============================================================================
#SBATCH --job-name=BALFEM_integ
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=28G
#SBATCH --output=%x.%A_%a.out
#SBATCH --error=%x.%A_%a.err

source "${BALFEM_PROJ:-$HOME/GridapBALFEM.jl}/run/balfem_env.sh"
set +e                                   # balfem_env.sh sets -e; a failing case must not end the pack
here="$BALFEM_PROJ/run/integrator_study"
: "${BALFEM_STUDY_OUT:=$BALFEM_PROJ/output/integrator_study}"
: "${BALFEM_SEG_HOURS:=119}"             # timeout of one segment (the wall is 119:59)
: "${BALFEM_MAX_SEG:=6}"                 # resubmissions stop after this many segments
pack="${SLURM_ARRAY_TASK_ID:?run through submit.sh (an array task id is needed)}"

line=$(awk -F'\t' -v p="$pack" 'NR > 1 && $1 == p' "$here/packs.tsv")
[ -n "$line" ] || { echo "no pack $pack in packs.tsv" >&2; exit 1; }
IFS=',' read -r -a names <<< "$(cut -f3 <<< "$line")"
echo "$(date '+%F %T') pack $pack (tier $(cut -f2 <<< "$line")): ${#names[@]} cases on $(hostname)"

#  the sysimage: present and fresh, or BALFEM_NO_SYSIMAGE=1 (then every case compiles)
jflag=()
if [ "${BALFEM_NO_SYSIMAGE:-0}" != "1" ]; then
    [ -f "$BALFEM_SYSIMAGE" ] || { echo "sysimage not found: $BALFEM_SYSIMAGE" >&2; exit 1; }
    BALFEM_STRICT_SYSIMAGE=1 balfem_check_sysimage_freshness || exit 1
    jflag=(-J"$BALFEM_SYSIMAGE")
fi

#  one segment of one case (runs in the background)
run_segment() {
    local name="$1" row case_dir last seg prev restart
    row=$(awk -F'\t' -v n="$name" 'NR > 1 && $1 == n' "$here/cases.tsv")
    [ -n "$row" ] || { echo "case $name not in cases.tsv" >&2; return; }
    IFS=$'\t' read -r _ tier model nlp flat slope T kA ncell dt periods solver tableau theta rho cost <<< "$row"
    case_dir="$BALFEM_STUDY_OUT/$name"; mkdir -p "$case_dir"
    [ -f "$case_dir/FINISHED" ] && { echo "  $name: finished"; return; }
    [ -f "$case_dir/FAILED" ]   && { echo "  $name: failed earlier (see FAILED)"; return; }
    last=$(ls -d "$case_dir"/seg[0-9][0-9] 2>/dev/null | sort | tail -1)
    seg=$(printf "%02d" $(( ${last:+10#${last##*seg}} + 1 )))
    #  resume from the newest segment that holds a checkpoint
    restart=""
    for prev in $(ls -d "$case_dir"/seg[0-9][0-9] 2>/dev/null | sort -r); do
        ls "$prev"/checkpoints/*.toml >/dev/null 2>&1 && { restart="$prev"; break; }
    done
    local segdir="$case_dir/seg$seg"; mkdir -p "$segdir"
    echo "  $name: segment $seg${restart:+ (resume from $(basename "$restart"))}"
    (
        export BALFEM_REGIME=nonlinear BALFEM_NL_PRESSURE="$nlp" BALFEM_FLAT_BED="$flat"
        export BALFEM_BED_SLOPE="$slope" BALFEM_BED_NW=1
        export BALFEM_M=2 BALFEM_P_VERT=1 BALFEM_FE_ORDER=3 BALFEM_P_ETA=2
        export BALFEM_NCELL="$ncell" BALFEM_NLAMBDA=1 BALFEM_D=3.5
        export BALFEM_TWAVE="$T" BALFEM_KA="$kA" BALFEM_DT="$dt" BALFEM_PERIODS="$periods"
        export BALFEM_SOLVER="$solver" BALFEM_TABLEAU="$tableau" BALFEM_THETA_T="$theta" BALFEM_RHO_INF="$rho"
        export BALFEM_NOISE=1e-8 BALFEM_USE_AD=0
        export BALFEM_SAVE_EVERY=5 BALFEM_DIAG_EVERY=5 BALFEM_PRINT_EVERY=25
        export BALFEM_CHECKPOINT_EVERY=250 BALFEM_CHECKPOINT_KEEP=2
        export BALFEM_OUTDIR="$segdir"
        unset BALFEM_CIP_GU BALFEM_CIP_GE                    # the driver default: 0.3 with 𝓝, none without
        [ -n "$restart" ] && export BALFEM_RESTART_FROM="$restart"
        export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
        cd "$BALFEM_PROJ"
        timeout --signal=TERM "${BALFEM_SEG_HOURS}h" \
            stdbuf -oL -eL julia ${jflag[@]+"${jflag[@]}"} --project="$BALFEM_PROJ" \
            examples/local_1d/run_periodic_1d.jl > "$segdir/run.log" 2>&1
        st=$?
        echo "$st" > "$segdir/exit_status"
        if grep -q "done:" "$segdir/run.log"; then
            touch "$case_dir/FINISHED"
        elif [ "$st" -ne 124 ]; then
            { echo "segment $seg exit $st at $(date '+%F %T')"; grep -m3 '^ERROR' "$segdir/run.log"; } > "$case_dir/FAILED"
        fi
    ) &
}

for n in "${names[@]}"; do
    run_segment "$n"
    sleep 20                                          # stagger the starts (sysimage mmap, disk)
done
wait
echo "$(date '+%F %T') pack $pack: all cases exited"

#  resubmit the pack if a case was cut by the timeout and has segments left
again=0
for n in "${names[@]}"; do
    d="$BALFEM_STUDY_OUT/$n"
    [ -f "$d/FINISHED" ] || [ -f "$d/FAILED" ] && continue
    nseg=$(ls -d "$d"/seg[0-9][0-9] 2>/dev/null | wc -l)
    if [ "$nseg" -lt "$BALFEM_MAX_SEG" ]; then again=1
    else echo "  $n: $nseg segments without finishing — not resubmitted (raise BALFEM_MAX_SEG)"; fi
done
if [ "$again" = 1 ]; then
    cd "$BALFEM_PROJ/output/slurm_logs/integrator_study" 2>/dev/null || true
    sbatch --array="$pack" --export=ALL "$here/run_pack.sh"
    echo "pack $pack resubmitted for its next segment"
fi
