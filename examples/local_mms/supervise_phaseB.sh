#!/bin/bash
#  Supervised Phase-B batch. Workers are RECYCLED after PHASEB_MAX_STUDIES so
#  RSS cannot drift (measured: 1.5 -> 4.6 GB over ~20 h, CPU falling to 7-14 %).
#  Resume in the runner makes a restart cost one JIT rather than the shard.
#  Distinct shard ids never share jobs (round-robin partition), and the same id
#  is never run twice concurrently, so no CSV is written by two processes.
cd "$(dirname "$0")/../.." || exit 1
NW=${NW:-5}; NSHARD=12
export PHASEB_MAX_STUDIES=${PHASEB_MAX_STUDIES:-2}
OUT=output/local/mms/phaseB
rows() { local t=0; for f in $OUT/shard_*.csv; do [ -f "$f" ] && t=$((t + $(wc -l < "$f") - 1)); done; echo $t; }
round=0
while :; do
  before=$(rows); round=$((round+1))
  echo "[sup] round $round  rows=$before/40  $(date +%H:%M:%S)"
  for i in $(seq 0 $((NW-1))); do
    s=$(( (round*NW + i) % NSHARD ))
    nohup julia --project=. examples/local_mms/run_phaseB_shard.jl $s $NSHARD \
      >> $OUT/shard_$(printf %02d $s).log 2>&1 &
  done
  wait
  after=$(rows)
  echo "[sup] round $round done  rows=$after/40"
  [ "$after" -ge 40 ] && { echo "[sup] COMPLETE"; break; }
  if [ "$after" -le "$before" ]; then
    stall=$((stall+1)); [ "$stall" -ge 3 ] && { echo "[sup] STALLED (3 rounds, no new rows)"; break; }
  else stall=0; fi
done
