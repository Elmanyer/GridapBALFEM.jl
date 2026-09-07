#!/bin/bash
# SMALL run — BC-generated NONLINEAR :native plane wave, flat bed, A=0.1, 30 periods.
#
# PURPOSE. Amplitude dispersion -- the mechanism that drives Benjamin-Feir --
# is a THIRD-ORDER Stokes effect, and per StokesWaveFourierAnalysis.tex it is
# carried by components {4,7} of the nonlinear pressure. Component 7 is in
# :native; component 4 is :full-ONLY. So :full is the only tier with a complete
# amplitude-dispersion term, which would explain why only :full has shown the
# instability. PREDICTION: this run should show a substantially weaker envelope
# growth than its :full twin, or none. If instead it grows the same way, the
# mechanism is advection rather than the {1,2,4,5} pressure block.
#  SNELLIUS ROME SIZING: 28 ranks x 4 GiB = 112 GiB = exactly 4/8 of a node.
#  Was 32 ranks, which billed 5/8 (memory fraction 4.57/8 rounds up). The rome node is
#  128 cores / 224 GiB divided into EIGHTHS, and the LARGER of the core/memory fraction
#  sets the bill -- at 4 GiB/rank memory always wins, so the cost-optimal count is 7k
#  ranks for tier k/8. See run/SNELLIUS_ROME_LAUNCH_CONFIGS.md.
#SBATCH --job-name="BALFEM_nl_native_bcplane_bigamp"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
# Memory: the rome node default (2 GB/core) was TESTED (2026-08) and OOM-killed.
# 4 GB/core is required and is MEASURED: the small-domain runs peaked at 3633
# MB/rank, so 3 GB/core would OOM too. RSS rose over the first ~100 steps then
# PLATEAUED -- no per-step leak. See building_files/OPEN_ITEMS.md section 1.
# A rome node advertises 256 GB but SLURM can allocate only ~224 GB, so
# 64 ranks/node x 4 GB = 256 GB/node is REFUSED at submit time. 32 ranks x 4 GB
# = 128 GB/node fits and is the proven small-domain configuration.
#SBATCH --mem-per-cpu=4G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
source $HOME/GridapBALFEM.jl/run/balfem_env.sh
export BALFEM_PX=7
export BALFEM_PY=4            # 7*4 = 28 ranks
export BALFEM_REGIME=nonlinear
export BALFEM_NL_PRESSURE=native
export BALFEM_FLAT_BED=1
export BALFEM_AWAVE=0.1
export BALFEM_PERIODS=30
balfem_run 32 examples/distributed_small/run_bc_plane_small.jl
