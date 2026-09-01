#!/bin/bash
# SMALL run — BC-generated NONLINEAR :native plane wave, SHALLOW (d=0.85), A=0.1, 30 periods.
#  SHALLOW WATER: d=0.85 m at T=2.0 s gives kd=1.079, which is 21 % BELOW the
#  Benjamin-Feir threshold kd=1.363. Below that threshold the envelope equation
#  is DEFOCUSING and a uniform train is modulationally STABLE, so if the growth
#  seen at kd=3.53 is genuine modulational instability it must DISAPPEAR here.
#  That is the experiment.
#
#  Why depth and not period. Lowering kd by raising T would give lambda=22 m on a
#  50 m domain -- under 3 wavelengths, and a 12 m sponge would be only 0.55
#  lambda, violating the rule that the sponge must cover the longest component.
#  Lowering d instead SHORTENS lambda to 4.95 m: 10 wavelengths in the domain and
#  a sponge of 2.4 lambda. Resolution improves rather than degrades.
#
#  Safety margins at A=0.1, d=0.85: H/d = 0.235 (depth-limited breaking ~0.78),
#  H/lambda = 0.040 against a Miche limit of 0.111 (36 %), Ursell = 8.0 (below 26,
#  so still Stokes-like rather than cnoidal). Nothing here is near breaking.
#
#  OUTDIR IS SET EXPLICITLY. The driver builds its output tag from
#  regime/nl_pressure/bed/A/T and does NOT include the depth, so a shallow run
#  would otherwise overwrite its deep-water twin of the same tier.
#SBATCH --job-name="BALFEM_nl_native_bcplane_shallow"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
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
export BALFEM_PX=8
export BALFEM_PY=4            # 8*4 = 32 ranks
export BALFEM_REGIME=nonlinear
export BALFEM_NL_PRESSURE=native
export BALFEM_FLAT_BED=1
export BALFEM_AWAVE=0.1
export BALFEM_D=0.85          # kd = 1.079  (deep twin: d=3.5, kd=3.527)
#  30 periods = 60 s. Transit at c_g=1.86 is 26.8 s, +4 s ramp +3T settle = 37 s,
#  leaving ~23 s of steady state -- more margin than the deep cases get.
export BALFEM_PERIODS=30
export BALFEM_OUTDIR=$HOME/GridapBALFEM.jl/output/small_bcplane_nonlinear_native_shallow_d0.85_A0.1_T2.0_P1LFE-2
balfem_run 32 examples/distributed_small/run_bc_plane_small.jl
