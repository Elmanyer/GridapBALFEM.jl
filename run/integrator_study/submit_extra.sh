#!/bin/bash
# ============================================================================
#  submit_extra.sh — submit tiers of the EXTRA integrator set: TR-BDF2 and generalised-α
#  (ρ∞ = 0 and 0.5), added 2026-11-18 after the main matrix had been submitted (cases.py, EXTRA SET)
#
#    bash run/integrator_study/submit_extra.sh P          # the extra pilot (3 cases, s/step per scheme)
#    bash run/integrator_study/submit_extra.sh E A G      # tiers: P A B C D E G
#    bash run/integrator_study/submit_extra.sh status     # per-case state of the extra set
#
#  Same as submit.sh with BALFEM_STUDY_SET=extra: tables cases_extra.tsv / packs_extra.tsv, job names
#  BALFEM_integ_extra_<tier>, output in the same $BALFEM_STUDY_OUT (the case names differ by the
#  integrator tag: trbdf2, galpha0, galpha50). Same cluster rules (cluster.sh, BALFEM_CLUSTER).
# ============================================================================
set -euo pipefail
: "${BALFEM_PROJ:=$HOME/GridapBALFEM.jl}"
[ $# -ge 1 ] || { sed -n '6,8p' "$0"; exit 1; }
BALFEM_STUDY_SET=extra exec bash "$BALFEM_PROJ/run/integrator_study/submit.sh" "$@"
