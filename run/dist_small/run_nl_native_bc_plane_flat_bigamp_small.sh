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
#SBATCH --job-name="BALFEM_nl_native_bcplane_bigamp"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
# 4 GiB/rank is REQUIRED (measured peak 3633 MB/rank; 2 GB/core was OOM-killed).
# Sizing: 28 ranks/node x 4 GiB = 112 GiB/node = 4/8 of a rome node (28 total).
# Rules and tiers: run/SNELLIUS_ROME_LAUNCH_CONFIGS.md
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
balfem_run 28 examples/distributed_small/run_bc_plane_small.jl
