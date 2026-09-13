#!/bin/bash
# SMALL run — nonlinear NATIVE irregular sea, flat, Hs=0.2 (physics-tier control)
#SBATCH --job-name="BALFEM_nl_native_irregular"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --ntasks=56
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
# 4 GiB/rank is REQUIRED (measured peak 3633 MB/rank; 2 GB/core was OOM-killed).
# Delivered as --mem=112G on each of 2 nodes / 28 ranks = 4 GiB per rank.
# Sizing: 112 GiB/node = 4/8 of a rome node, x2 nodes = 8/8 billed, 56 ranks.
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
# ⚠ --mem, NOT --mem-per-cpu. Slurm on rome counts --mem-per-cpu against ALL
# ALLOCATED CPUs, and the allocation is rounded up to the 16-core minimum — so
# --ntasks-per-node=28 is given 32 cores and 4G/cpu becomes 32*4 = 128 GiB, not
# the 112 intended, and the request lands on a node configuration that does not
# exist:
#   sbatch: error: All allocated CPUs are counted when using --mem-per-cpu,
#           which is a multiple of the minimum allocation size.
#   sbatch: error: Batch job submission failed: Requested node configuration
#           is not available
# --mem is a flat per-node total and does not depend on how many cores Slurm
# decides to hand over. 112 GiB/node = 4/8 of a rome node; x2 nodes = 8/8 billed
# for 56 ranks, i.e. 4 GiB per rank exactly.
#SBATCH --mem=112G
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
#  ⚠ THE POINT OF THIS RUN. `:full` is the ONE physics tier the analytic MMS does
#  not verify (models 7-8): its {{1,2,4,5}} blocks are assembled in the residual
#  but ABSENT from jacobian_u, and its frozen projections are refreshed from the
#  PREVIOUS step (update_nlp_state! runs after the step), i.e. an explicitly
#  lagged term inside an otherwise implicit scheme. Every run in the 2026-09
#  small-domain batch was `:full`, so the batch had NO control and could not say
#  whether the transverse velocity instability belongs to the model or to that
#  tier. `:native` is the production tier and IS verified. This is the control.
#  Read: y-symmetry of u3x (see the note in the flat launcher).
export BALFEM_NL_PRESSURE=native
balfem_run 56 examples/distributed_small/run_irregular_sea_small.jl
