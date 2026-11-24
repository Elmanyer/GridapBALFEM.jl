#!/bin/bash
# PRODUCTION — linear ring wave (interior Gaussian point source), FLAT bed, A=0.1, T=2.0
#   domain 40 x 40 m, source at the centre, 4-side sponges 7 m (1.12 lambda),
#   mesh 160 x 160 (square 0.25 m cells — isotropy matters for a radial pattern), Q2/Q1, CN
#   30 periods = 60 s = 3000 steps (settles at 8.1 T; ~22 T of steady record)
#   NOTE: the interior source over-delivers (rule 14b, eta/A ~ 2-3); in the linear
#   model that is a pure scale factor, not a defect.
#SBATCH --job-name="BALFEM_lin_ring_flat"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=42
#SBATCH --ntasks-per-node=42
#SBATCH --cpus-per-task=1
# MEMORY: --mem only, never --mem-per-cpu (run/SNELLIUS_ROME_LAUNCH_CONFIGS.md §2, §4).
#   168 GiB / 42 ranks = 4.00 GiB per rank, billed 6/8 of a rome node (96 CPUs).
#   25600 cells / 42 = 610 cells/rank -> ~2.9 GB/rank (v1 linear: ~2.4 MB/cell
#   above a ~1.45 GB baseline), 72 % of the cap.
#   WHY ONE NODE x 42 AND NOT THE DRIVER'S 8x8=64: 64 x 4 GiB = 256 GiB does not
#   fit a node, and one node beats two when the problem fits (§4).
#   Wall time is the least certain of the suite: the ring carries walls and the
#   largest cells/rank. If it runs short, go to 56 = 2 nodes x 28 at --mem=112G
#   with PX=8 PY=7 (billed 8/8).
#SBATCH --mem=168G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source $HOME/GridapBALFEM.jl/run/balfem_env.sh
source "$BALFEM_PROJ/run/linear_production/_linear_common.sh"

# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
export BALFEM_PX=7
export BALFEM_PY=6            # 7*6 = 42 ranks; 22.9 x 26.7 cells per rank

# --- Case ---------------------------------------------------------------------
export BALFEM_FLAT_BED=1      # ring: flat bed only (the driver has no bar variant)
export BALFEM_AWAVE=0.1
export BALFEM_TWAVE=2.0
export BALFEM_LX=40.0 BALFEM_LY=40.0
export BALFEM_NX=160  BALFEM_NY=160
export BALFEM_SPONGE_X=7.0 BALFEM_SPONGE_Y=7.0
export BALFEM_PERIODS=30
export BALFEM_SAVE_EVERY=10   # 0.2 s snapshots = 10 per period

balfem_run 42 examples/distributed_small/run_ring_small.jl
