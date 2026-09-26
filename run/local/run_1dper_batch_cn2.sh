#!/bin/bash
# ==============================================================
#  run_1dper_batch_cn2.sh — Crank–Nicolson RESOLUTION LADDER, Q2/Q1 (2026-09-26)
#
#  THE QUESTION. Every periodic run so far used SDIRK_2_2, which is L-stable and damps an
#  oscillatory mode by ≈ γ⁴(ωΔt)⁴ per step. Two readings depend on that damping:
#    (1) the carrier loses energy at ≈ 0.010 /s in EVERY run — SDIRK at the carrier explains
#        only 1.1e-4 /s, so either energy is transferred to fast modes the integrator kills,
#        or the spatial discretisation itself dissipates;
#    (2) the element-scale growth measured in :full (σ_E ≈ 0.10 /s at Q2/Q1 λ/h=32; divergence
#        at Q3/Q2 λ/h=16) is a NET rate — the fastest :full modes reach ρ ≈ 20–40 rad/s, where
#        SDIRK damping (≈ 0.08 /s at ρ = 20.6) is comparable with the growth.
#  Crank–Nicolson (θ = 1/2) has |R(iω Δt)| ≡ 1: it damps nothing, so it exposes the spatial
#  rate unmasked and tells where the carrier energy goes.
#     carrier decay persists under CN   ⇒ the SPATIAL discretisation dissipates
#     carrier decay vanishes under CN   ⇒ it was SDIRK acting on transferred energy
#     :full growth faster under CN      ⇒ SDIRK was masking part of the interior instability
#
#  DESIGN (rules 14c, 38c): identical to run_1dper_batch_global.sh except BALFEM_SOLVER=theta
#  (θ = 0.5 is the default); every :full arm has its :native twin in the same batch. The two
#  cells are the ones that discriminate: Q2/Q1 at λ/h = 32 (slow growth under SDIRK) and
#  Q3/Q2 at λ/h = 16 (diverged at 35 s under SDIRK).
#  ⚠ On the mixed layout the auxiliary rows carry no time derivative; the θ-method enforces
#  them at the midpoint state. That is a bounded (sign-alternating) inconsistency, not a drift,
#  but a Nyquist-in-time oscillation of 𝖦 is its signature if one appears.
#
#  Analysis per run:
#    PG_NCELL=<n> PG_PU=<p> julia --project=postprocessing \
#        postprocessing/examples/periodic_growth.jl output/local_1d/periodic/<case>
#
#  ⚠ Logs go to persistent disk (rule 47); never edit this file while it runs (rule 38h).
#  USAGE:  run/local/run_1dper_batch_cn.sh [MAX_PARALLEL]  (default 4; the limit counts ALL
#  running periodic solvers, so it shares the machine with run_1dper_batch_global.sh)
# ==============================================================
cd "$(dirname "${BASH_SOURCE[0]}")/../.." || exit 1
MAXP="${1:-4}"
OUT=output/local_1d/periodic
LOG=$OUT/_logs
mkdir -p "$LOG"

# case : name | FE_ORDER | NCELL | AWAVE | NL_PRESSURE | MIXED | extra env
#  Second CN batch. The first (run_1dper_batch_cn.sh) showed that SDIRK_2_2 masked growth:
#  Q2/Q1 mixed at λ/h = 32 diverged at ≈ 90 s under CN (bounded under SDIRK), and the SDIRK run at
#  λ/h = 64 showed NO growth — consistent with SDIRK damping its fastest modes at ≈ 10 /s. This
#  batch completes the Q2/Q1 ladder (16, 32, 64) without integrator damping, each with its twin.
CASES=(
  "p21_mixed_n16_cn   2 16 0.10 full   1 BALFEM_C3_MASK=gs"
  "p21_native_n16_cn  2 16 0.10 native 0"
  "p21_mixed_n64_cn   2 64 0.10 full   1 BALFEM_C3_MASK=gs"
  "p21_native_n64_cn  2 64 0.10 native 0"
)

run_case() {
  read -r name fe ncell A nlp mixed extra <<< "$1"
  [ -f "$LOG/$name.log" ] && { echo "skip $name (log exists)"; return; }
  (
    export BALFEM_REGIME=nonlinear BALFEM_FLAT_BED=1 BALFEM_NL_PRESSURE=$nlp BALFEM_MIXED=$mixed
    export BALFEM_FE_ORDER=$fe BALFEM_P_ETA=$((fe - 1)) BALFEM_NCELL=$ncell BALFEM_NLAMBDA=1
    export BALFEM_AWAVE=$A BALFEM_TWAVE=1.6 BALFEM_DT=0.04 BALFEM_PERIODS=100
    export BALFEM_SOLVER=theta BALFEM_USE_AD=0 BALFEM_NOISE=1e-8
    export BALFEM_SAVE_EVERY=5 BALFEM_DIAG_EVERY=5 BALFEM_PRINT_EVERY=25
    export BALFEM_OUTDIR=$OUT/$name
    [ -n "$extra" ] && export $extra
    source run/local/balfem_local.sh
    balfem_local_run examples/local_1d/run_periodic_1d.jl
  ) > "$LOG/$name.log" 2>&1
}

for c in "${CASES[@]}"; do
  while [ "$(pgrep -fc '^[^ ]*/julia .*run_periodic_1d[.]jl')" -ge "$MAXP" ]; do sleep 30; done
  echo "$(date '+%F %T')  launch  $(echo $c | cut -d' ' -f1)"
  run_case "$c" &
  sleep 60          # stagger: let each process take the precompile lock in turn
done
wait
echo "$(date '+%F %T')  batch finished"
