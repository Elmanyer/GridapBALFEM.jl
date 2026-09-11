#!/bin/bash
# ==============================================================
#  run_distributed_tests.sh — the 4-rank MPI test tier, SEQUENTIALLY.
#
#  One at a time on purpose: each needs 4 ranks at ~2 GB, so two concurrently is
#  ~16 GB and would race the recovery runs for memory. These are the tests that
#  cover the CLUSTER path, which no sequential test reaches.
#
#  ⚠ test_mms_distributed_parity is the important one: it gates the distributed
#  MMS branch against the sequential one, so it is what would catch the two paths
#  disagreeing about the Taylor-Hood pairing after the p_u/p_eta rename.
#
#  ⚠ Body wrapped in { } and ending in `exit` so bash parses it whole — editing or
#  git-checkout'ing a running script otherwise makes it resume at a byte offset and
#  re-execute (that destroyed six finished runs on 2026-09-06).
#
#  MPI_Finalize prints a benign OFI error and exits 143 on this machine; treated as
#  success, exactly as balfem_local_mpi does.
# ==============================================================
{
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
MPIEXECJL="${MPIEXECJL:-$HOME/.julia/bin/mpiexecjl}"
LOG=output/local/logs/dist; mkdir -p "$LOG"; SUMM="$LOG/summary.txt"; : > "$SUMM"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 JULIA_NUM_THREADS=1
for t in test_basic_distributed.jl test_nlpressure_distributed.jl \
         test_bc_generation_distributed.jl test_mms_distributed_parity.jl; do
  [ -f "test/$t" ] || { echo "$t MISSING" >> "$SUMM"; continue; }
  out="$LOG/${t%.jl}.log"
  timeout 7200 stdbuf -oL "$MPIEXECJL" --project=. -n 4 julia --project=. "test/$t" > "$out" 2>&1
  rc=$?
  [ $rc -eq 143 ] && rc=0                    # benign MPI_Finalize/OFI cleanup
  #  Verdict from GATE OUTPUT, never the exit code (rule 35).
  res=$(grep -oE '[0-9]+ PASS,? +[0-9]+ FAIL' "$out" | tail -1)
  [ -z "$res" ] && res=$(grep -oiE '(PASS|FAIL)[^ ]*$' "$out" | tail -1)
  [ $rc -eq 124 ] && res="TIMEOUT(2h)"
  printf '%-38s rc=%-3s %s\n' "$t" "$rc" "${res:-<no gate line>}" >> "$SUMM"
done
echo "ALL DISTRIBUTED TESTS DONE" >> "$SUMM"
}
exit
