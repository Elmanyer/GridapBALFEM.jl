#!/bin/bash
# 1-D NONLINEAR INSTABILITY STUDY -- nl_full, A=0.05
#   + N components {1,2,4,5}. Component 4 is :full-ONLY, and {4,7} are what carry AMPLITUDE DISPERSION -- the third-order effect that drives Benjamin-Feir. If BF is the mechanism, growth should be strongest here
#
#  THE QUESTION. A bcplane 2-D run at nonlinear/:full grew exponentially
#  (envelope rate +0.087 1/s, e-folding 7.9 s) and would have diverged. It is NOT
#  sponge reflection: the linear run on the identical sponge settles to an
#  envelope max/min of 1.07, and the sponge-region amplitude stays flat while the
#  interior doubles. The leading hypothesis is Benjamin-Feir modulational
#  instability. This flume tests it in 1-D, where the phenomenon actually lives
#  and a case costs minutes rather than the 10-20 h a 2-D run takes.
#
#  WHY AN AMPLITUDE LADDER. BF growth goes as (ka)^2, so the rate must scale with
#  the SQUARE of the delivered amplitude -- a signature little else reproduces.
#  At kd=5.50 (well above the 1.363 threshold), lambda=4.00 m, Miche H/lambda=0.14:
#      A=0.05  17.9% of Miche   BF e-folding 57.1 s
#      A=0.10  35.7%                        14.3 s
#      A=0.15  53.6%                         6.3 s
#      A=0.20  71.5%                         3.6 s
#  This case is A=0.05 (17.9% of Miche, e-folding 57.1 s).
#  >> The existing run_1d_*.sh cases run at A=0.001, where the e-folding is
#  >> 142845 s (40 h). They are a tier-AGREEMENT check and cannot show this at all.
#
#  >> WHAT THE LADDER CANNOT DO ALONE. The quasi-Newton gap in jacobian_u (which
#  >> omits every N derivative by design) ALSO scales as A^2, so an A^2 law is
#  >> necessary but not sufficient. The tier axis is what separates them: a solver
#  >> artefact should not care which N components are present, whereas BF should
#  >> be strongest in :full, the only tier carrying a complete amplitude-dispersion
#  >> term.
#
#  Sequential, ny=1 + solid walls, boundary generation, 50 periods (~40 s of
#  steady state after fill and ramp -- enough to fit an envelope growth rate).
#  At A=0.05 that is under one e-folding, so treat that rung as a stability
#  reference rather than a rate measurement.
source "$(dirname "${BASH_SOURCE[0]}")/balfem_local.sh"
export BALFEM_WAVE_GEN=bc
export BALFEM_REGIME=nonlinear
export BALFEM_NL_PRESSURE=full
export BALFEM_FLAT_BED=1
export BALFEM_RELAX=1
export BALFEM_MPI=0
export BALFEM_LX=60.0
export BALFEM_NX=240
export BALFEM_NY=1
export BALFEM_YBC=wall
export BALFEM_LY=0.25
export BALFEM_AWAVE=0.05
export BALFEM_PERIODS=50
export BALFEM_SAVE_EVERY=25
export BALFEM_DIAG_EVERY=5
balfem_local_run examples/local_1d/run_flume_1d.jl
