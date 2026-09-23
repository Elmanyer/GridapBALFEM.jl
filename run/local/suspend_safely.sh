#!/bin/bash
# ==============================================================
#  SUSPEND THIS WORKSTATION WITHOUT LOSING RUN EVIDENCE
#
#  Suspend (unlike reboot) preserves BOTH the running solvers and /tmp.
#  This script only adds the pre-flight that makes a bad resume survivable:
#    1. copy every solver log from the volatile /tmp scratchpad onto disk
#       (rule 47: this machine's tmpfiles.d carries `D /tmp` -> wiped on boot)
#    2. sync(1), so diagnostics/VTK/docs are on the platter, not in page cache
#    3. then suspend
#
#  It KILLS NOTHING. The solvers keep running across the suspend.
# ==============================================================
set -u
PROJ="/home/pmanyerfuertes/Documents/TU-Delft/WP2/BALFEM_model/GridapBALFEM.jl"
ARCHIVE="$PROJ/output/local_1d/_logs_2026-09-23"
cd "$PROJ" || exit 1

echo "--- pre-flight $(date '+%F %T') ---"

# 1. locate the live scratchpad from a running solver's own stdout fd, so this
#    keeps working even if the session directory name changes.
SCRATCH=""
for pid in $(pgrep -f "bin/julia.*run_flume_1d" 2>/dev/null); do
    t=$(readlink "/proc/$pid/fd/1" 2>/dev/null)
    [ -n "$t" ] && { SCRATCH=$(dirname "$t"); break; }
done
[ -z "$SCRATCH" ] && SCRATCH=$(ls -td /tmp/claude-*/*/*/scratchpad 2>/dev/null | head -1)

mkdir -p "$ARCHIVE"
if [ -n "$SCRATCH" ] && [ -d "$SCRATCH" ]; then
    echo "  scratchpad: $SCRATCH"
    cp -p "$SCRATCH"/*.log "$ARCHIVE"/ 2>/dev/null
    n=0
    for d in output/local_1d/c3[vq]_*_mixed output/local_1d/c3[vq]_*_proj; do
        [ -d "$d" ] || continue
        b=$(basename "$d")
        [ -f "$SCRATCH/$b.log" ] && { cp -p "$SCRATCH/$b.log" "$d/run.log"; n=$((n+1)); }
    done
    echo "  archived $(ls "$ARCHIVE"/*.log 2>/dev/null | wc -l) logs; refreshed $n arm run.log files"
else
    echo "  WARNING: scratchpad not found — logs NOT refreshed"
fi

# 2. record what was alive, so the state is readable after the resume
{
  echo "=== suspend at $(date '+%F %T') ==="
  for d in output/local_1d/c3[vq]_*_mixed output/local_1d/c3[vq]_*_proj; do
      [ -f "$d/diagnostics.csv" ] || continue
      b=$(basename "$d")
      pgrep -f "run_1d${b}\.sh" >/dev/null && st=LIVE || st=dead
      echo "  $b $st last=$(tail -1 "$d/diagnostics.csv" | cut -d, -f2,3,7,12)"
  done
} >> "$ARCHIVE/suspend_state.log"
echo "  state written to $ARCHIVE/suspend_state.log"

# 3. flush to disk
echo "  syncing filesystem..."
sync
echo "  synced."

echo "  solvers still running: $(pgrep -cf 'bin/julia.*run_flume_1d')"
echo "--- suspending now ---"
sleep 2
systemctl suspend
