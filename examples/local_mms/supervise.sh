#!/bin/bash
#  supervise.sh — keep ONE shard slot busy with fresh processes.
#
#  The shard runner exits voluntarily once it exceeds VBC_MAX_LIFE_S (a long-lived
#  Julia+Gridap process degrades to GC-bound: measured 3.6-4.5 GB and 7-14 % CPU
#  after days, one level taking 6.2 days that takes seconds when fresh). This loop
#  restarts it. Resume logic means a fresh process picks up whatever is still
#  missing, so restarting is free apart from one JIT.
#
#  Exits when a run completes without doing any work (nothing left for this slot).
#  usage: supervise.sh <slot> <nslots> <logdir> [extra env already exported]
slot=$1; nslots=$2; logdir=$3
cd "$(dirname "$0")/../.." || exit 1     # repo root, wherever it is checked out
for gen in $(seq 1 100); do
  log=$logdir/slot${slot}_gen${gen}.log
  VBC_SHARD=$slot VBC_NSHARD=$nslots julia --project=. \
      examples/local_mms/run_vbasis_shard.jl > $log 2>&1
  rc=$?
  n=$(grep -c "DONE \|ERROR " $log 2>/dev/null)
  echo "[sup $slot] gen$gen rc=$rc completed=$n $(date '+%m-%d %H:%M')" >> $logdir/sup_$slot.log
  # nothing assigned, or nothing accomplished and no lifetime cut -> queue is empty
  if grep -q "0 cases assigned" $log 2>/dev/null; then
    echo "[sup $slot] queue empty — stopping $(date '+%m-%d %H:%M')" >> $logdir/sup_$slot.log; break
  fi
  if [ "$n" -eq 0 ] && ! grep -q "LIFETIME" $log 2>/dev/null; then
    echo "[sup $slot] no progress and no lifetime cut (rc=$rc) — stopping" >> $logdir/sup_$slot.log; break
  fi
done
