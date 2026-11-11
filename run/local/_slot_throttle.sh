# ==============================================================
#  _slot_throttle.sh — shared slot gate for the local 1-D launchers (sourced, not run)
#
#  Every local 1-D launcher counts the SAME processes, so launchers started side by side share
#  the machine instead of each assuming it owns it:
#    * the box driver   examples/local_1d/run_periodic_1d.jl  (main tree AND the no-CIP worktree)
#    * the flume driver examples/local_1d/run_flume_1d.jl
#    * the sponge worker run/local/sponge_worker.jl
#  A slot is free when the memory left after every running case reaches its expected peak is
#  ≥ MIN_FREE_GB (the 2026-10-30 OOM: memory is the binding limit on this workstation; one 𝓝
#  case ≈ 2.5–2.8 GB, a linear flume worker ≈ 1.5–2 GB). A case starts small and grows during
#  JIT, so each running case younger than its peak is charged RUN_GB − RSS on top of what it
#  already holds; MemAvailable alone would admit several cases before any of them has grown.
#  The check and the launch happen under one flock, followed by a stagger, so two launchers
#  can never both see the same free slot.
#
#  USAGE (in a launcher):
#    source "$PROJ/run/local/_slot_throttle.sh"
#    with_slot <command> [args…]     # blocks until a slot is free, then runs the command in the
#                                    # background; the PID is in $SLOT_PID
# ==============================================================
#  MEMORY IS THE ONLY CAP (author, 2026-11-10). The 2026-11-07 thermal throttle (> 3 busy cores →
#  1.1 GHz) was a machine fault that a reboot cleared; the old MAXP = 3 default is withdrawn.
#  MAXP stays as an optional process cap for one batch (unset or 0 = none).
MAXP="${MAXP:-0}"
MIN_FREE_GB="${MIN_FREE_GB:-6}"
RUN_GB="${RUN_GB:-3}"                   # GB: expected peak RSS of one case (𝓝 ≈ 2.8)
SLOT_STAGGER="${SLOT_STAGGER:-60}"      # s: let a new process take the precompile lock first
SLOT_LOCK="${SLOT_LOCK:-$PROJ/output/.local_slots.lock}"

slots_busy() {
  pgrep -fc '^[^ ]*julia .*(run_periodic_1d|run_flume_1d|sponge_worker)[.]jl'
}
mem_free_gb() {   # MemAvailable minus the growth still owed by the running cases, in GB
  local avail owed=0 rss pid peak=$(( RUN_GB * 1048576 ))
  avail=$(awk '/^MemAvailable:/ {print $2}' /proc/meminfo)
  for pid in $(pgrep -f '^[^ ]*julia .*(run_periodic_1d|run_flume_1d|sponge_worker)[.]jl'); do
    rss=$(awk '/^VmRSS:/ {print $2}' "/proc/$pid/status" 2>/dev/null) || continue
    [ -n "$rss" ] && [ "$rss" -lt "$peak" ] && owed=$(( owed + peak - rss ))
  done
  echo $(( (avail - owed) / 1048576 ))
}
slot_free() {
  { [ "$MAXP" -le 0 ] || [ "$(slots_busy)" -lt "$MAXP" ]; } && [ "$(mem_free_gb)" -ge "$MIN_FREE_GB" ]
}

with_slot() {
  exec 9>"$SLOT_LOCK"
  while true; do
    flock 9
    if slot_free; then
      "$@" 9>&- &          # the case must not inherit the lock descriptor
      SLOT_PID=$!
      sleep "$SLOT_STAGGER"
      flock -u 9
      return 0
    fi
    flock -u 9
    sleep 30
  done
}
