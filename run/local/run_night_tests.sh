#!/bin/bash
# ==============================================================
#  run_night_tests.sh — re-baseline the whole sequential suite on TAYLOR-HOOD.
#
#  WHY. Until 2026-09-06 the `p_eta = 0` sentinel meant EQUAL ORDER, so every test
#  that called setup_and_run without an explicit p_eta ran Q2/Q2. The sentinel now
#  resolves to Taylor-Hood (rule 2b), so those tests have changed discretisation and
#  their pinned constants -- measured on the deficient pairing -- must be
#  RE-MEASURED, not re-thresholded. This runs them all and records what each does.
#
#  ⚠ The verification tier (test_mms_*, test_jacobians_ad, test_linear_newton_gate)
#  always passed p_eta explicitly and was ALREADY Taylor-Hood, so any failure there
#  is a real regression, not a re-baselining artefact. The report separates them.
#
#  ⚠ SCRIPT HARDENING: the whole body is wrapped in { ... } and the file ends with
#  `exit`. Bash reads a script incrementally BY BYTE OFFSET, so replacing the file
#  while it runs (a git checkout, an edit) makes it resume at that offset in the new
#  content and re-execute. That happened on 2026-09-06 and relaunched six finished
#  runs over their own output. A brace-wrapped body is parsed in full before any of
#  it executes, which makes that impossible.
# ==============================================================
{
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
JOBS="${JOBS:-6}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1
LOG=output/local/logs/night; mkdir -p "$LOG"
SUMM="$LOG/summary.txt"; : > "$SUMM"

#  Ordered cheapest-first so the fast signals land early (SJF: the deliverable here
#  is COVERAGE, not makespan -- CLAUDE.md rule 43).
TESTS=(
  test_taylor_hood.jl test_primitives.jl test_vertical.jl test_dispersion_curve.jl
  test_waveinput.jl test_mms_forcing.jl test_mms_forcing_nonlinear.jl
  test_basic.jl test_conservation.jl test_sloshing.jl test_shallow_water.jl
  test_dispersion.jl test_energy.jl test_convergence.jl test_vertical_profile.jl
  test_selfconsistency.jl test_nlpressure.jl test_dispersion_nonlinear.jl
  test_bc_generation.jl test_bc_spectrum.jl test_linear_newton_gate.jl
  test_mms_convergence.jl test_mms_convergence_nonlinear.jl test_jacobians_ad.jl
)
for t in "${TESTS[@]}"; do
  [ -f "test/$t" ] || { echo "$t MISSING" >> "$SUMM"; continue; }
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 5; done
  ( out="$LOG/${t%.jl}.log"
    timeout 10800 stdbuf -oL julia --project=. "test/$t" > "$out" 2>&1
    rc=$?
    #  Verdict from GATE OUTPUT, never from the exit code (rule 35): a clean exit is
    #  not evidence a test ran, and these print their own PASS/FAIL tallies.
    res=$(grep -oE '[0-9]+ PASS,? +[0-9]+ FAIL' "$out" | tail -1)
    [ -z "$res" ] && res=$(grep -oE 'Results: +[0-9]+ PASS, +[0-9]+ FAIL' "$out" | tail -1)
    [ $rc -eq 124 ] && res="TIMEOUT(3h)"
    printf '%-38s rc=%-3s %s\n' "$t" "$rc" "${res:-<no gate line>}" >> "$SUMM"
  ) &
done
wait
echo "ALL NIGHT TESTS DONE" >> "$SUMM"
}
exit
