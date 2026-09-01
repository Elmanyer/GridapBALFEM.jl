#!/bin/bash
# WIDE run — LINEAR directional sea, flat bed, Hs=0.2, on a WIDENED domain with
# MUCH THINNER lateral sponges. 48 ranks (the cluster's ceiling at 4 GB/core).
#
#  WHY. On the 20 m-wide small domain the lateral sponges were 8 m each, so
#  16 of 20 m of the domain -- and 16 of 20 m of the Dirichlet inflow boundary --
#  sat inside a sponge. The undamped strip was only 4 m, which is 0.45 of the
#  energy-carrying wavelength (lambda_eff ~ 8.9 m at T ~ 1.5 Tp). A short-crested
#  sea with +/-40 deg spreading cannot develop its transverse structure in less
#  than one oblique wavelength, so the interior never carried the intended Airy
#  sea state: it was damped before it could form.
#
#  THE CHANGE. Ly 20 -> 40 m and the lateral sponges 8 -> 4 m. The undamped strip
#  goes 4 -> 32 m, i.e. 0.45 -> 3.6 lambda_eff. There is now ample room for the
#  directional structure to exist.
#
#  >> THE TRADE-OFF, STATED PLAINLY. A 4 m lateral sponge is 0.45 lambda_eff --
#  >> half of the 8 m it replaces, and still under one wavelength.
#  >> Sponge strength SATURATES past mu_max ~ 5*omega (= 15.7 here, against the
#  >> default 40), so beyond that point WIDTH is the only lever -- and this sponge
#  >> is now very narrow and very strong, which is the configuration closest to a
#  >> partially reflecting wall. Expect imperfect lateral absorption.
#  >>
#  >> WATCH eta_max_damped / eta_max_int in diagnostics.csv. If it climbs above
#  >> ~1, the lateral sponges are reflecting rather than absorbing. The fix is
#  >> then EITHER a wider Ly (keeping the strip generous) OR a gentler mu_max --
#  >> a lower mu spread over a thin layer reflects less than a violent one, and
#  >> mu_max=40 buys nothing over ~16 anyway. Set BALFEM_MUMAX=16 to try that.
#  >> Do NOT respond by making the sponge stronger; that is the saturated axis.
#
#  MESH AND COST. dx = dy = 0.25 m is KEPT: coarsening to hold the cell count
#  down would under-resolve the short-wave tail (kd up to 5.99, lambda 3.67 m),
#  which is exactly the part of the spectrum this run exists to preserve.
#  200 x 160 = 32000 cells, twice the small-domain case, on the same 48 ranks.
#
#  >> MEMORY. The linear directional case peaked at 2642 MB/rank at 333
#  >> cells/rank. Subtracting the ~1450 MB Julia+Gridap baseline gives ~3.6 MB
#  >> per cell per rank, so at 800 cells/rank (32000/40) this run is predicted at
#  >> ~4310 MB -- 84 % of the 5 GB cap. Raising the request from 4 to 5 GB/core
#  >> is what buys that headroom: at 4 GB/core this mesh needs 105 % of the cap
#  >> on 40 ranks and 94 % on 48, neither of which leaves any margin.
#  >> If it is OOM-killed (exit 137), the cheapest fix that preserves the physics
#  >> is to coarsen ~12 %: BALFEM_NX=176 BALFEM_NY=140 (dx=0.284, dy=0.286)
#  >> gives 513 cells/rank and ~3290 MB, at 12.9 cells per shortest wavelength
#  >> instead of 14.7. Reducing Ly would defeat the purpose of the run.
#
#  LINEAR ON PURPOSE. Its nonlinear twin diverged at t=3.14 of 52 s, inside the
#  ramp, with the maximum pinned at the inflow and dmp/int of 8-12. Establishing
#  that the WIDE geometry carries a clean linear sea state is the prerequisite
#  for re-testing the nonlinear case; running the nonlinear one first would
#  confound the geometry change with the instability under investigation.
#
#  OUTDIR IS SET EXPLICITLY: the driver's tag carries regime/tier/bed/Hs but not
#  the geometry, so this would otherwise overwrite the 20 m-wide linear run.
#SBATCH --job-name="BALFEM_lin_directional_wide"
#SBATCH --partition=rome
#SBATCH --time=119:59:00
#SBATCH --ntasks=40
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --cpus-per-task=1
# Memory: 5 GB/core requested (up from the 4 GB the rest of the suite uses),
# because this run carries twice the cells of the small-domain cases on fewer
# ranks. THE RANK COUNT IS SET BY THE PER-NODE CEILING, NOT BY THE CORE LIMIT:
#   48 x 5 GB = 240 GB/node  -> REFUSED. A rome node advertises 256 GB but SLURM
#                               can allocate only ~224 GB; this is the same wall
#                               that rejected 64 x 4 GB = 256 GB.
#   40 x 5 GB = 200 GB/node  -> the largest that fits, and what is used here.
# FALLBACK if 200 GB/node is still refused: 36 ranks on a (6,6) grid = 180 GB,
# which is below the 192 GB already proven to submit. That gives 889 cells/rank
# (~4630 MB, 90 % of the 5 GB cap) instead of 800 (~4310 MB, 84 %).
#SBATCH --mem-per-cpu=5G
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

source $HOME/GridapBALFEM.jl/run/balfem_env.sh

# --- MPI process grid (PX*PY MUST equal the rank count given to balfem_run) ---
#  (8,5) on 200x160 divides EXACTLY both ways: 25 x 32 cells per rank, i.e.
#  6.25 x 8.00 m subdomains (aspect 1.28). The alternative (10,4) gives
#  5 x 10 m subdomains (aspect 2.00) and more halo for the same rank count.
export BALFEM_PX=8
export BALFEM_PY=5            # 8*5 = 40 ranks; 200/8 = 25, 160/5 = 32 cells each (exact)

# --- Case-specific overrides ---
export BALFEM_REGIME=linear
export BALFEM_NL_PRESSURE=none
export BALFEM_FLAT_BED=1
export BALFEM_LY=40.0         # was 20.0
export BALFEM_NY=160          # keeps dy = 0.25 m
export BALFEM_SPONGE_Y=4.0    # was 8.0 -> undamped strip 4 m -> 32 m (3.6 lambda_eff)
#  Halve the VTK write frequency: snapshot buffers add to peak RSS and this run
#  has only ~6 % headroom. 130 snapshots over 52 s is still ample for animation.
export BALFEM_SAVE_EVERY=20
export BALFEM_OUTDIR=$HOME/GridapBALFEM.jl/output/small_directional_linear_none_flat_Hs0.2_WIDE_P1LFE-2

balfem_run 40 examples/distributed_small/run_directional_sea_small.jl
