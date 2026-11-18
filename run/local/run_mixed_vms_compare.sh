#!/bin/bash
# ==============================================================
#  run_mixed_vms_compare.sh — the mixed + VMS comparison arms in parallel (branch mixed-vms-solver)
#
#  One Julia process per arm of examples/local_1d/compare_mixed_vms.jl (BALFEM_CMP_ARMS), plus the
#  Newton probe examples/local_1d/vms_newton_probe.jl, admitted by the shared memory gate of
#  _slot_throttle.sh (memory is the only cap, rule 47). When every process has ended, the summary
#  mode rebuilds output/mixed_vms/comparison/summary.txt from the arms on disk.
#
#  An arm whose log exists is skipped (delete the log to rerun it). Never edit this file while it
#  runs (rule 38h).
#  Detach: setsid nohup run/local/run_mixed_vms_compare.sh > output/mixed_vms/_logs/_batch.out 2>&1 < /dev/null &
# ==============================================================
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
PROJ="$(pwd)"
LOG=output/mixed_vms/_logs
mkdir -p "$LOG"
MEM_PER_CASE="${MEM_PER_CASE:-9G}"
export RUN_GB="${RUN_GB:-5}"            # a fresh process with the mixed 7-field operator peaks higher than an 𝓝 box case
export SLOT_STAGGER="${SLOT_STAGGER:-20}"
source "$PROJ/run/local/_slot_throttle.sh"

#  The shared gate counts only the box/flume/sponge drivers; count these scripts as well.
PAT='^[^ ]*julia .*(run_periodic_1d|run_flume_1d|sponge_worker|compare_mixed_vms|vms_newton_probe)[.]jl'
slots_busy() { pgrep -fc "$PAT"; }
mem_free_gb() {
  local avail owed=0 rss age pid peak=$(( RUN_GB * 1048576 ))
  avail=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
  for pid in $(pgrep -f "$PAT"); do
    age=$(ps -o etimes= -p "$pid" 2>/dev/null) || continue
    [ -n "$age" ] && [ "$age" -lt "$YOUNG_S" ] || continue
    rss=$(awk '/^VmRSS:/ {print $2}' "/proc/$pid/status" 2>/dev/null) || continue
    [ -n "$rss" ] && [ "$rss" -lt "$peak" ] && owed=$(( owed + peak - rss ))
  done
  echo $(( (avail - owed) / 1048576 ))
}

#  Arms, longest first (the mixed arms carry 7 fields and the longest JIT).
ARMS="${ARMS:-m6:mixed_vms m5:mixed_vms m6:mixed m5:mixed m6:primal_cip m5:primal_cip m4:primal_vms}"

run_job() {   # name script [env assignments…]
  local name=$1 script=$2; shift 2
  (
    export "$@"
    export OPENBLAS_NUM_THREADS=1
    echo "job $name | $script | $* | started $(date '+%F %T') | git $(git rev-parse --short HEAD)$(git diff --quiet -- src examples || echo ' + uncommitted edits')"
    systemd-run --user --scope --quiet --unit="mixvms-$name-$$" -p MemoryMax="$MEM_PER_CASE" \
      stdbuf -oL -eL julia --project="$PROJ" "$script"
    echo "--- exit status $? at $(date '+%F %T') ---"
  ) > "$LOG/$name.log" 2>&1
}

launch() {    # name script [env…]
  local name=$1
  [ -f "$LOG/$name.log" ] && { echo "skip $name (log exists)"; return; }
  with_slot run_job "$@"
  echo "$(date '+%F %T')  launch  $name  ($(slots_busy) busy, $(mem_free_gb) GB headroom)" | tee -a "$LOG/_batch.log"
}

launch probe_m3 examples/local_1d/vms_newton_probe.jl BALFEM_PROBE_MODEL=3 BALFEM_PROBE_STEPS=3
for a in $ARMS; do
  launch "cmp_${a/:/_}" examples/local_1d/compare_mixed_vms.jl BALFEM_CMP_ARMS="$a"
done
wait
echo "$(date '+%F %T')  all jobs ended; building the summary" | tee -a "$LOG/_batch.log"
BALFEM_CMP_MODE=summary stdbuf -oL julia --project="$PROJ" examples/local_1d/compare_mixed_vms.jl > "$LOG/summary_build.log" 2>&1
echo "$(date '+%F %T')  SUMMARY-DONE" | tee -a "$LOG/_batch.log"
