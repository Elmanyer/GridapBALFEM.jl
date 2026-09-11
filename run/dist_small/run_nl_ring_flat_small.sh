#!/bin/bash
# SMALL run — nonlinear full ring wave, flat, A=0.1 (BASE)
#SBATCH --job-name="BALFEM_nl_ring"
#
# NOTE ON THE PRE-2026-09-11 OUTPUT. Unlike the irregular and directional cases,
# this one is NOT affected by the per-rank sea-state bug: the ring is driven by
# an INTERIOR Gaussian source (wave_gen=:inner_res) and never calls
# build_airy_state, so no WaveSpec randomness enters it. Its archived run
# remains a valid record of the :full source-region failure (died t=6.76,
# max pinned on the source at x=20.0, u_max 0.69 -> 1.90, Newton 6 -> 32).
# This relaunch is a dt/dx refinement test against that record, not a repair.
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --ntasks=56
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --cpus-per-task=1
# 4 GiB/rank is REQUIRED (measured peak 3633 MB/rank; 2 GB/core was OOM-killed).
# Sizing: 56 ranks/node x 4 GiB = 224 GiB/node = 8/8 of a rome node (56 total).
# Rules and tiers: run/SNELLIUS_ROME_LAUNCH_CONFIGS.md
#
# >> COST ARITHMETIC FOR THE REFINEMENT BELOW (2026-09-11). dx 0.25 -> 0.20
# >> raises the mesh 25600 -> 40000 cells (1.56x) and dt 0.02 -> 0.01 doubles the
# >> step count to 2400, so the job is ~2.9x the old one. 42 ranks would not fit
# >> the wall; 56 gives 714 cells/rank -> ~3090 MB (75 % of cap, ring measured
# >> 2.3 MB/cell/rank) and ~152 s/step -> ~101 h against the 119:59 limit.
# >> A FINER MESH NEEDS A SECOND NODE: dx=0.167 (nx=ny=240) is 57600 cells,
# >> ~3820 MB/rank at 56 ranks (93 % of cap, too close) and ~147 h. Run that as
# >> --ntasks=112 --nodes=2 --ntasks-per-node=56, PX=14 PY=8, ~73 h.
#SBATCH --mem-per-cpu=4G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source $HOME/GridapBALFEM.jl/run/balfem_env.sh

# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
export BALFEM_PX=8
export BALFEM_PY=7            # 8*7 = 56 ranks

# --- Case-specific overrides (only what differs from the script's base) ---
#  ⚠ WHY dt AND dx BOTH MOVE. The 2026-09 run died at t=6.76 with the maximum
#  pinned on the source point (x=20.0) and u_max going 0.69 -> 1.90 while Newton
#  went 6 -> 32 iterations. Halving dt tests the time discretisation; refining dx
#  tests whether the source region was under-resolved.
export BALFEM_DT=0.01
export BALFEM_NX=200          # was 160 -> dx 0.25 -> 0.20 m
export BALFEM_NY=200          # square cells kept
#  dt/dx are not fields of the output-name grammar (rule 2c), so name them here
#  or this run silently becomes the old one's _v2.
export BALFEM_NAME_EXTRA=dt0.01,dx0.20

#  >> WHAT NEITHER KNOB ADDRESSES, STATED SO IT IS NOT MISREAD AS FIXED. The
#  >> interior source over-delivers (rule 14b): A_wave=0.1 put eta at 0.35-0.45 m
#  >> on the source by t=3.8 and 0.95 m at death, against a Miche breaking
#  >> amplitude of ~0.44 m at kd=3.53 (H_max = 0.142*lambda*tanh(kd) = 0.88 m).
#  >> The run was AT the breaking limit from t~3.8 and 2.2x past it when it died,
#  >> and the model has no breaking closure. If this refinement still fails near
#  >> the source, the next variable is amplitude, not resolution:
#  >> BALFEM_AWAVE=0.025 brings the delivered crest to ~0.11 m, a quarter of Miche.

balfem_run 56 examples/distributed_small/run_ring_small.jl
