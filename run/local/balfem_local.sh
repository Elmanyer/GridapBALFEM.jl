#!/bin/bash
# ==============================================================
#  balfem_local.sh — shared helper for LOCAL runs (this workstation, 12-rank partitions)
#
#  The local counterpart of run/balfem_env.sh. That helper is cluster-only: it
#  loads Snellius/DelftBlue modules and resolves the prebuilt system image.
#  Locally there is no sysimage and no scheduler, so this helper only has to
#  resolve the project, cap the rank count, and launch.
#
#  SOURCE IT, then call one of:
#      balfem_local_run  <script.jl>            # sequential (KEEPS point gauges)
#      balfem_local_mpi  <nranks> <script.jl>   # MPI, nranks ≤ BALFEM_MAX_RANKS
#
#  Sequential is the default for anything that measures: the distributed driver
#  evaluates no point gauges (timeloop_dist.jl:21-22), so an MPI local run can
#  only be judged from the diagnostics CSV. Use balfem_local_mpi when exercising
#  the DISTRIBUTED path is the point of the run.
#
#  PARTITION SIZE. This workstation has 16 cores; a simulation is given a
#  12-rank partition, leaving 4 cores for the OS and for the analysis tooling.
#  The consequence is that only ONE simulation runs at a time — 12 ranks per
#  case is a deliberate trade of case throughput for per-case size/resolution,
#  not a free speed-up. (Two 6-rank partitions side by side would finish a
#  multi-case set sooner; a single 12-rank partition runs a bigger case.)
#
#  Knobs
#    BALFEM_PROJ        project directory      (default: this file's ../..)
#    BALFEM_MAX_RANKS   hard cap on ranks      (default 12 of this box's 16 cores)
#    JULIA            julia binary           (default: julia on PATH)
#    MPIEXECJL        MPI launcher           (default ~/.julia/bin/mpiexecjl)
#
#  Notes
#   * ALWAYS use mpiexecjl, never the system mpiexec: the latter fails on this
#     machine with a PMIx version mismatch.
#   * MPI_Finalize prints a benign OFI error from the WiFi NIC and exits 143.
#     balfem_local_mpi treats 143 as success — do not "fix" that by hiding real
#     failures: any other non-zero status is still reported.
#   * Julia buffers stdout hard when redirected, so stdbuf is used throughout;
#     otherwise a running job's log file looks empty until it exits.
# ==============================================================

_balfem_local_here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export BALFEM_PROJ="${BALFEM_PROJ:-$(cd "$_balfem_local_here/../.." && pwd)}"
BALFEM_MAX_RANKS="${BALFEM_MAX_RANKS:-12}"
JULIA="${JULIA:-julia}"
MPIEXECJL="${MPIEXECJL:-$HOME/.julia/bin/mpiexecjl}"

# --- Taylor-Hood pairing guard (CLAUDE.md rule 2b) -------------------------
#  Resolve and VALIDATE the horizontal element pair before anything expensive
#  starts. The solver gates this too (check_taylor_hood), but on the cluster that
#  is a queue wait plus a JIT away, and on a workstation it is ~30 min of
#  compilation — a bad pairing should cost seconds, not an allocation.
#
#  WHY IT IS A HARD REQUIREMENT: eta enters momentum undifferentiated, via div(v)
#  after the integration by parts, so it plays the pressure role of a Stokes system
#  and equal-order continuous spaces are inf-sup deficient. Equal order is what
#  produced the year-long "nonlinear instability" -- an unbounded grid-scale mode at
#  lambda ~ 2*dx that got WORSE under refinement. CLAUDE.md rules 2b and 12b.
#
#  Every launcher gets this by sourcing the helper, so the pairing is enforced and
#  RECORDED IN THE LOG for local, cluster, sequential and distributed runs alike --
#  which is what makes a run's discretisation auditable after the fact.
balfem_require_taylor_hood() {
    export BALFEM_FE_ORDER="${BALFEM_FE_ORDER:-2}"
    export BALFEM_P_ETA="${BALFEM_P_ETA:-$((BALFEM_FE_ORDER - 1))}"
    if [ "$BALFEM_P_ETA" -lt 1 ] 2>/dev/null; then
        echo "FATAL: BALFEM_P_ETA=$BALFEM_P_ETA — the surface order must be >= 1." >&2
        return 2
    fi
    if [ "$BALFEM_FE_ORDER" -ne $((BALFEM_P_ETA + 1)) ] 2>/dev/null; then
        echo "FATAL: NON-TAYLOR-HOOD pairing Q$BALFEM_FE_ORDER/Q$BALFEM_P_ETA." >&2
        echo "  BALFE-M requires BALFEM_FE_ORDER = BALFEM_P_ETA + 1 (velocity one order" >&2
        echo "  above the surface). Equal order is inf-sup deficient here and caused the" >&2
        echo "  2026-09 grid-scale instability — CLAUDE.md rules 2b and 12b." >&2
        echo "  Use FE_ORDER=$((BALFEM_P_ETA + 1)) P_ETA=$BALFEM_P_ETA, or FE_ORDER=$BALFEM_FE_ORDER P_ETA=$((BALFEM_FE_ORDER - 1))." >&2
        return 2
    fi
    echo "  elements : Q${BALFEM_FE_ORDER}/Q${BALFEM_P_ETA} (Taylor-Hood, verified pairing)"
    return 0
}

balfem_local_banner() {
    echo "--------------------------------------------------------------"
    echo " GridapBALFEM local run"
    echo "   project : $BALFEM_PROJ"
    echo "   script  : $1"
    echo "   mode    : $2"
    echo "   started : $(date '+%Y-%m-%d %H:%M:%S')"
    echo "--------------------------------------------------------------"
}

# --- sequential ------------------------------------------------------------
balfem_local_run() {
    local script="$1"
    [ -n "$script" ] || { echo "balfem_local_run: no script given" >&2; return 2; }
    [ -f "$BALFEM_PROJ/$script" ] || [ -f "$script" ] || {
        echo "balfem_local_run: script not found: $script" >&2; return 2; }
    [ -f "$script" ] || script="$BALFEM_PROJ/$script"

    balfem_require_taylor_hood || return 2
    balfem_local_banner "$script" "sequential (gauges available)"
    stdbuf -oL -eL "$JULIA" --project="$BALFEM_PROJ" "$script"
    local rc=$?
    echo "--- exit status $rc ---"
    return $rc
}

# --- MPI -------------------------------------------------------------------
balfem_local_mpi() {
    local n="$1"; local script="$2"
    [ -n "$n" ] && [ -n "$script" ] || {
        echo "balfem_local_mpi: usage balfem_local_mpi <nranks> <script.jl>" >&2; return 2; }
    if [ "$n" -gt "$BALFEM_MAX_RANKS" ]; then
        echo "balfem_local_mpi: refusing $n ranks — this machine is capped at" \
             "$BALFEM_MAX_RANKS (raise BALFEM_MAX_RANKS deliberately if you mean it)" >&2
        return 2
    fi
    [ -f "$script" ] || script="$BALFEM_PROJ/$script"
    [ -f "$script" ] || { echo "balfem_local_mpi: script not found: $2" >&2; return 2; }
    [ -x "$MPIEXECJL" ] || {
        echo "balfem_local_mpi: $MPIEXECJL not found. Install it with" >&2
        echo "    julia -e 'using MPI; MPI.install_mpiexecjl()'" >&2
        return 2; }

    balfem_require_taylor_hood || return 2
    balfem_local_banner "$script" "MPI, $n ranks (NO point gauges — read diagnostics.csv)"
    export BALFEM_MPI=1
    stdbuf -oL -eL "$MPIEXECJL" --project="$BALFEM_PROJ" -n "$n" \
        "$JULIA" --project="$BALFEM_PROJ" "$script"
    local rc=$?
    # 143 = SIGTERM during MPI_Finalize: the documented benign OFI cleanup error
    # on this machine. The computation has already completed at that point.
    if [ "$rc" -eq 143 ]; then
        echo "--- exit status 143 (benign MPI_Finalize/OFI cleanup — treated as success) ---"
        return 0
    fi
    echo "--- exit status $rc ---"
    return $rc
}
