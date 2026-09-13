#!/bin/bash
# SMALL run — nonlinear NATIVE irregular sea over a bar, Hs=0.2 (physics-tier control)
#SBATCH --job-name="BALFEM_nl_native_irregular_varbed"
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
export BALFEM_FLAT_BED=0
balfem_run 56 examples/distributed_small/run_irregular_sea_small.jl
