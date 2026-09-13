#!/bin/bash
# SMALL run — linear directional sea, flat, Hs=0.2
#SBATCH --job-name="BALFEM_lin_directional"
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
# >> WHY 56 AND NOT 42 (2026-09-11). Ly 20 -> 40 m doubles the mesh to
# >> 200x160 = 32000 cells. Memory, not cores, sets the tier: the measured cost
# >> is ~3.8 MB per cell per rank above a ~1450 MB Julia+Gridap baseline, so
# >>     42 ranks -> 762 cells/rank -> ~4350 MB  => OVER the 4 GiB cap
# >>     56 ranks -> 571 cells/rank -> ~3620 MB  => 88 % of the cap
# >> IF IT IS OOM-KILLED (exit 137, or 'oom-kill event' in .err): go to 2 nodes,
# >> --ntasks=84 --nodes=2 --ntasks-per-node=42 (6/8 per node) with PX=14 PY=6,
# >> which cuts cells/rank to 381 (~2900 MB). ⛔ NOT 56/node — that is 224 GiB on a
# >> node and is refused at submit time (see the SBATCH block). Do NOT coarsen the mesh: dx=0.25 m is 14.7 cells per the
# >> shortest wavelength (lambda=3.67 m at kd=5.99) and that short-wave tail is
# >> exactly what a directional run exists to carry.
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
export BALFEM_PX=8
export BALFEM_PY=7            # 8*7 = 56 ranks

# --- Case-specific overrides (only what differs from the script's base) ---
export BALFEM_REGIME=linear
export BALFEM_NL_PRESSURE=none

balfem_run 56 examples/distributed_small/run_directional_sea_small.jl
