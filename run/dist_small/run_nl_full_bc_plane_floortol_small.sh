#!/bin/bash
# SMALL run — bcplane nonlinear :full A=0.1, Newton tolerance 1e-9 (FLOOR PROBE).
#  WHY THIS TEST. jacobian_u OMITS EVERY nonlinear-pressure (N) derivative by
#  design -- problem.jl states it: "the N blocks add to the residual but not
#  here". Newton is therefore QUASI-Newton, and the size of the omitted block
#  scales as A^2. That matters because the envelope growth under investigation
#  ALSO scales as A^2 (Benjamin-Feir goes as (ka)^2), so an amplitude sweep alone
#  CANNOT separate the two. This test can: a quasi-Newton artefact is a property
#  of the SOLVE and must move when the tolerance moves; genuine modulational
#  instability is a property of the OPERATOR and must not.
#
#  READING IT. Compare the envelope growth rate against the production run
#  (nl_tol=1e-5, rate +0.087 1/s):
#    rate essentially unchanged  -> the growth is in the operator, i.e. physical;
#    rate changes materially     -> Newton is converging to the fixed point of
#                                   the wrong map (CLAUDE.md rule 5).
#
#  REACHABILITY. The production run achieved a minimum Newton residual of
#  4.6e-08 (median 7.1e-07) while merely stopping at its 1e-5 threshold, so the
#  quasi-Newton floor is at most 4.6e-08.
#  1e-9 sits BELOW the observed 4.6e-08 floor ON PURPOSE. Two informative
#  outcomes, and a failure here is a RESULT, not a wasted run:
#    it converges       -> the quasi-Newton gap is smaller than the production
#                          logs suggested, and 2a's comparison is the clean one;
#    Newton stalls      -> the gap is pinned between 1e-9 and 4.6e-08 at this
#                          amplitude, quantifying for a PRODUCTION run what the
#                          MMS work measured only for a manufactured solution
#                          (where the stall floor jumped ~500000x between
#                          manufactured amplitudes 0.4 and 0.8).
#  Read the stall value out of the .out warnings, not just the exit status.
#SBATCH --job-name="BALFEM_nl_full_bcplane_tol1e9"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
# Memory: 4 GB/core is required and MEASURED (this exact case peaked at 2882 MB
# on 32 ranks / 250 cells per rank). 32 x 4 GB = 128 GB/node, well inside the
# ~224 GB a rome node can allocate. See building_files/OPEN_ITEMS.md section 1.
#SBATCH --mem-per-cpu=4G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source $HOME/GridapBALFEM.jl/run/balfem_env.sh
# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
export BALFEM_PX=8
export BALFEM_PY=4            # 8*4 = 32 ranks
# --- the reference configuration, held fixed across this whole family ---
export BALFEM_REGIME=nonlinear
export BALFEM_NL_PRESSURE=full
export BALFEM_FLAT_BED=1
#  30 periods = 60 s. The original run stopped at 32 s, just as the envelope
#  growth became visible (rate +0.087 1/s, e-folding 7.9 s). 60 s leaves ~18 s of
#  steady state after transit+ramp+settle, i.e. ~2.3 further e-foldings.
export BALFEM_PERIODS=30
export BALFEM_AWAVE=0.1
export BALFEM_NL_TOL=1e-9
export BALFEM_NL_ITER=100
export BALFEM_OUTDIR=$HOME/GridapBALFEM.jl/output/small_bcplane_nlfull_A0.1_nltol1e-9_P1LFE-2
balfem_run 32 examples/distributed_small/run_bc_plane_small.jl
