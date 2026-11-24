#!/bin/bash
# ============================================================================
#  run_precluster_checks.sh — the two local checks before the integrator study goes to Snellius
#  (markdown_files/CH8_CAMPAIGN_TODO.md, "Integrator study: hypotheses", 2026-11-12)
#
#    bash run/local/run_precluster_checks.sh            # both, side by side, detached
#
#  CHECK 1 (≈ 10–30 min, postprocessing env): what grew in the m3 k0A = 0.25 SDIRK_2_2 blow-up
#    (output/1D_stability_campaign/_logs/kA25_mode.jl): the growing harmonics, their frequencies,
#    the decay-triad test, and the pump coupling Γ that tier F (Crouzeix / RK4 Δt ladder) faces.
#  CHECK 2 (hours, main env): the frozen-state spectrum of m5 n32 against γ and kA, with the
#    verdict of each (integrator, Δt) (examples/local_1d/frozen_gamma_sweep.jl). Sets GAMMAS_G in
#    run/integrator_study/cases.py around the predicted boundary.
#  Logs: output/precluster/check{1,2}_*.log; markers .done / .failed next to them.
# ============================================================================
set -u
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
out="$ROOT/output/precluster"; mkdir -p "$out"

run_check() {    # name, log, command...
    local name="$1" log="$2"; shift 2
    (
        echo "$(date '+%F %T') start $name" >> "$out/_checks.log"
        "$@" > "$log" 2>&1
        st=$?
        echo "$(date '+%F %T') end   $name (exit $st)" >> "$out/_checks.log"
        [ "$st" -eq 0 ] && touch "$log.done" || touch "$log.failed"
    ) &
}

export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
run_check check1 "$out/check1_kA25_mode.log" \
    julia --project=postprocessing output/1D_stability_campaign/_logs/kA25_mode.jl
run_check check2 "$out/check2_frozen_gamma_n32.log" \
    env FGS_NCELL=32 FGS_OUT="$out/frozen_gamma_sweep" julia --project=. examples/local_1d/frozen_gamma_sweep.jl
wait
