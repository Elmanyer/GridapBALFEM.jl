#!/bin/bash
# SMALL run — nonlinear FULL irregular sea, flat, Hs=0.2 (BASELINE / control for the twins)
#SBATCH --job-name="BALFEM_nl_irregular"
#
# ⚠ THE PRE-2026-09-11 OUTPUT OF THIS CASE IS INVALID, NOT A BASELINE. It ran
# with the per-rank sea-state bug (build_airy_state seeded the phases but not the
# angular spreading), so its Dirichlet inflow carried a different sea on every
# rank. `unique_output_dir` will suffix this run _v2 rather than overwrite it;
# that _v2 is the FIRST valid run of this configuration, not a second attempt.
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --ntasks=56
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
# 4 GiB/rank is REQUIRED (measured peak 3633 MB/rank; 2 GB/core was OOM-killed).
# Sizing: 28 ranks/node x 4 GiB = 112 GiB/node = 4/8 of a rome node, x2 nodes
#         = 8/8 billed for 56 ranks total.
# ⚠ 2 NODES x 28, NOT 1 x 56. Same 56 ranks and the SAME BILL (4/8 + 4/8 = 8/8),
# but 112 GiB per node instead of 224. The 8/8 row of SNELLIUS_ROME_LAUNCH_CONFIGS.md
# §2 is arithmetically correct and PRACTICALLY UNSUBMITTABLE: a rome node allocates
# ~224 GiB, so asking for exactly 224 leaves zero headroom and Slurm refuses it at
# submit time (same failure mode as the 64 x 4G = 256 GiB trap documented there).
# 4/8 is the tier 23 other launchers in this repo already use.
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
#  WHAT TO READ WHEN IT FINISHES. Not eta — eta stayed bounded at 0.13-0.14 m
#  through the entire 2026-09 failure. Read the TRANSVERSE symmetry of u3x: this
#  case is exactly y-invariant (long-crested, y_wall_bc=:wall, y-uniform bed), so
#  std_y(u3x)/|mean_y(u3x)| must stay at round-off. It previously sat at 1e-3
#  from step 1 (that was the per-rank sea-state bug, now fixed) and then ran
#  0.027 -> 0.057 -> 0.114 over t = 50, 51, 52 while Newton still converged in 6
#  iterations. The run ENDED at t=52; it did not survive.
# (base config — no overrides needed)
balfem_run 56 examples/distributed_small/run_irregular_sea_small.jl
