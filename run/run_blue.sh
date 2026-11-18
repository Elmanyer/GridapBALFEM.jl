#!/bin/bash

#SBATCH --job-name="GridapBABALFEM"
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3900M
#SBATCH --output=GridapBABALFEM.%j.out
#SBATCH --error=GridapBABALFEM.%j.err

export BABALFEM_CLUSTER=blue
source $HOME/GridapBABALFEM.jl/run/babalfem_env.sh

babalfem_run 16 examples/distributed/run_plane_wave_dist.jl
