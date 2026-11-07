# ==============================================================
#  _slot_throttle.sh — shared slot gate for the local 1-D launchers (sourced, not run)
#
#  Every local 1-D launcher counts the SAME processes, so launchers started side by side share
#  the machine instead of each assuming it owns it:
#    * the box driver   examples/local_1d/run_periodic_1d.jl  (main tree AND the no-CIP worktree)
#    * the flume driver examples/local_1d/run_flume_1d.jl
#    * the sponge worker run/local/sponge_worker.jl
#  A slot is free when fewer than MAXP of them run AND MemAvailable ≥ MIN_FREE_GB (the
#  2026-10-30 OOM: memory, not cores, is the binding limit on this workstation; one 𝓝 case
#  ≈ 2.5–2.8 GB, a linear flume worker ≈ 1.5–2 GB).
#  The check and the launch happen under one flock, followed by a stagger, so two launchers
#  can never both see the same free slot.
#
#  USAGE (in a launcher):
#    source "$PROJ/run/local/_slot_throttle.sh"
#    with_slot <command> [args…]     # blocks until a slot is free, then runs the command in the
#                                    # background; the PID is in $SLOT_PID
# ==============================================================
#  MAXP default 3 (2026-11-07): with more than 3 busy cores this workstation reaches 100 °C within
#  minutes and the firmware drops every core to ~1.1 GHz (9–10 W) for ~5 min at a time, so 8
#  processes run at roughly half speed each and the campaign cases at a third (CH8_CAMPAIGN_TODO §10).
MAXP="${MAXP:-3}"
MIN_FREE_GB="${MIN_FREE_GB:-6}"
SLOT_STAGGER="${SLOT_STAGGER:-60}"      # s: let a new process take the precompile lock first
SLOT_LOCK="${SLOT_LOCK:-$PROJ/output/.local_slots.lock}"

slots_busy() {
  pgrep -fc '^[^ ]*julia .*(run_periodic_1d|run_flume_1d|sponge_worker)[.]jl'
}
mem_free_gb() {
  awk '/^MemAvailable:/ {printf "%d", $2/1048576}' /proc/meminfo
}
slot_free() {
  [ "$(slots_busy)" -lt "$MAXP" ] && [ "$(mem_free_gb)" -ge "$MIN_FREE_GB" ]
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
