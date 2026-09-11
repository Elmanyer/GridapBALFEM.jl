#!/bin/bash
#  Memory-AWARE Phase-B supervisor.
#
#  The first attempt at this batch ran 12 workers on a 31 GB machine; RSS drifted
#  1.5 -> 4.6 GB per worker over ~20 h and swap hit 7850/8191 MB, at which point
#  every worker was GC-bound and throughput was WORSE than a third as many.
#  This version cannot repeat that:
#
#    * it starts a worker ONLY when free RAM exceeds MIN_FREE_MB;
#    * a WATCHDOG kills the single largest worker whenever swap passes
#      SWAP_MAX_MB -- recoverable, because the runner resumes from the CSVs;
#    * PHASEB_MAX_STUDIES=1 recycles a worker after every study, so RSS drift is
#      bounded by ONE study instead of by the whole shard.
#
#  Distinct shard ids never share jobs (round-robin partition of a cost-sorted
#  queue) and this script never starts an id that is already running, so no CSV
#  is ever written by two processes.
cd "$(dirname "$0")/../.." || exit 1
NW=${NW:-8}
NSHARD=12
MIN_FREE_MB=${MIN_FREE_MB:-6000}
SWAP_MAX_MB=${SWAP_MAX_MB:-4500}
export PHASEB_MAX_STUDIES=${PHASEB_MAX_STUDIES:-1}
OUT=output/local/mms/phaseB

rows() { local t=0; for f in $OUT/shard_*.csv; do [ -f "$f" ] && t=$((t + $(wc -l < "$f") - 1)); done; echo $t; }
ids()  { pgrep -af 'run_phaseB_shard[.]jl [0-9]' | sed 's/.*shard[.]jl //' | awk '{print $1}'; }
freemb() { free -m | awk 'NR==2{print $7}'; }
swapmb() { free -m | awk 'NR==3{print $3}'; }

( while :; do
    if [ "$(swapmb)" -gt "$SWAP_MAX_MB" ]; then
      v=$(for p in $(pgrep -f 'run_phaseB_shard[.]jl [0-9]'); do
            echo "$(awk '/VmRSS/{print $2}' /proc/$p/status 2>/dev/null) $p"; done | sort -rn | head -1 | awk '{print $2}')
      if [ -n "$v" ]; then echo "[wd] swap $(swapmb)MB > $SWAP_MAX_MB -- killing fattest worker $v"; kill -9 "$v"; fi
      sleep 120
    fi
    sleep 30
  done ) & WD=$!
trap 'kill $WD 2>/dev/null' EXIT

while :; do
  r=$(rows)
  if [ "$r" -ge 42 ]; then echo "[sup] COMPLETE rows=$r/42"; break; fi
  cur=$(ids | wc -l)
  echo "[sup] $(date +%H:%M:%S) rows=$r/42 workers=$cur free=$(freemb)MB swap=$(swapmb)MB"
  if [ "$cur" -lt "$NW" ] && [ "$(freemb)" -gt "$MIN_FREE_MB" ]; then
    busy=" $(ids | tr '\n' ' ')"
    for s in $(seq 0 $((NSHARD-1))); do
      case "$busy" in *" $s "*) continue;; esac
      #  skip shards whose queue is empty -- see run_phaseB_shard.jl. Picking the
      #  lowest free id without this spins on an exhausted shard forever and
      #  starves the ids that still have work.
      [ -f "$OUT/$(printf 'exhausted_%02d' $s)" ] && continue
      nohup julia --project=. examples/local_mms/run_phaseB_shard.jl $s $NSHARD \
        >> $OUT/shard_$(printf %02d $s).log 2>&1 &
      echo "[sup] launched shard $s"
      break
    done
  fi
  sleep 90
done
