#!/bin/bash
# SMALL run — nonlinear full directional sea, flat, Hs=0.2 (BASE)
#SBATCH --job-name="BALFEM_nl_directional"
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
# >> WHY 56 AND NOT 42 (2026-09-11). Ly 20 -> 40 m doubles the mesh to
# >> 200x160 = 32000 cells. Memory, not cores, sets the tier: the measured cost
# >> is ~3.8 MB per cell per rank above a ~1450 MB Julia+Gridap baseline, so
# >>     42 ranks -> 762 cells/rank -> ~4350 MB  => OVER the 4 GiB cap
# >>     56 ranks -> 571 cells/rank -> ~3620 MB  => 88 % of the cap
# >> 56 is also the cost-optimal full-node count (7k ranks for tier k/8).
# >> IF IT IS OOM-KILLED (exit 137, or 'oom-kill event' in .err): go to 2 nodes,
# >> --ntasks=112 --nodes=2 --ntasks-per-node=56 with PX=14 PY=8, which halves
# >> cells/rank to 286. Do NOT coarsen the mesh: dx=0.25 m is 14.7 cells per the
# >> shortest wavelength (lambda=3.67 m at kd=5.99) and that short-wave tail is
# >> exactly what a directional run exists to carry.
#SBATCH --mem-per-cpu=4G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source $HOME/GridapBALFEM.jl/run/balfem_env.sh

# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
export BALFEM_PX=8
export BALFEM_PY=7            # 8*7 = 56 ranks

# --- Case-specific overrides (only what differs from the script's base) ---
# (base config — no overrides needed)

# Snapshot buffers add to peak RSS and this mesh runs at ~88 % of the cap, so
# halve the VTK cadence. 130 snapshots over 52 s is still ample for animation.
export BALFEM_SAVE_EVERY=20

balfem_run 56 examples/distributed_small/run_directional_sea_small.jl
