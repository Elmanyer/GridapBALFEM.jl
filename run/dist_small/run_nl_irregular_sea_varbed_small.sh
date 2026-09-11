#!/bin/bash
# SMALL run — nonlinear FULL irregular sea over a bar, Hs=0.2 (BASELINE)
#SBATCH --job-name="BALFEM_nl_irregular_varbed"
#
# ⚠ THE PRE-2026-09-11 OUTPUT OF THIS CASE IS INVALID, NOT A BASELINE. It ran
# with the per-rank sea-state bug (build_airy_state seeded the phases but not the
# angular spreading), so its Dirichlet inflow carried a different sea on every
# rank. `unique_output_dir` will suffix this run _v2 rather than overwrite it;
# that _v2 is the FIRST valid run of this configuration, not a second attempt.
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
# >> WHY ALL FOUR IRREGULAR JOBS MOVED 28 -> 56 RANKS (2026-09-11). They are a
# >> FACTORIAL, and the whole point is that they differ in exactly one axis each
# >> (CLAUDE.md rules 14c / 38c). The dt/2 members need 5200 steps, which at 28
# >> ranks is ~170 h against a 119:59 wall; putting only those on 56 ranks would
# >> have made the partition a second variable. So every member runs 56/(14,4).
# >> Memory is ample here: 16000/56 = 286 cells/rank -> ~2210 MB, 54 % of cap.
#SBATCH --mem-per-cpu=4G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source $HOME/GridapBALFEM.jl/run/balfem_env.sh

# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
#  (14,4) on 200x80: 14.3 x 20 cells/rank = 3.57 x 5.0 m subdomains, aspect 0.71.
#  The alternative (8,7) gives 6.25 x 2.86 m, aspect 2.19 — more halo, and rule 22
#  says anisotropy is paid for in the linear solve.
export BALFEM_PX=14
export BALFEM_PY=4            # 14*4 = 56 ranks

# --- Case-specific overrides (only what differs from the script's base) ---
#  The bar toe at x=20 FOCUSES the transverse mode: its flat twin and this run
#  agreed to 3 digits for 30 s, then this one died at t=31.96 while the flat one
#  ran to t=52. flat_bed is the single axis between them.
export BALFEM_FLAT_BED=0
balfem_run 56 examples/distributed_small/run_irregular_sea_small.jl
