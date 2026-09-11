#!/bin/bash
# WIDE run — LINEAR directional sea, flat bed, Hs=0.2.
#
# ⚠ 2026-09-11: THE "WIDE" GEOMETRY IS NOW THE BASE SCRIPT'S DEFAULT (Ly=40,
# lateral sponges 2 m, ny=160). The three overrides below are therefore no
# longer changes — they are kept, spelled out, because a launcher that states
# its geometry cannot be mistaken for one that inherited a different one. That
# ambiguity is what let equal-order and Taylor-Hood runs sit side by side for
# months (CLAUDE.md rule 12b). This file's remaining purpose is the LINEAR
# control for the nonlinear directional pair.
#
#  WHY. On the 20 m-wide small domain the lateral sponges were 8 m each, so
#  16 of 20 m of the domain -- and 16 of 20 m of the Dirichlet inflow boundary --
#  sat inside a sponge. The undamped strip was only 4 m, which is 0.45 of the
#  energy-carrying wavelength (lambda_eff ~ 8.9 m at T ~ 1.5 Tp). A short-crested
#  sea with +/-40 deg spreading cannot develop its transverse structure in less
#  than one oblique wavelength, so the interior never carried the intended Airy
#  sea state: it was damped before it could form.
#
#  THE CHANGE. Ly 20 -> 40 m and the lateral sponges 8 -> 2 m. The undamped strip
#  goes 4 -> 36 m, i.e. 0.45 -> 4.0 lambda_eff. There is now ample room for the
#  directional structure to exist.
#
#  >> THE TRADE-OFF, STATED PLAINLY. A 2 m lateral sponge is 0.22 lambda_eff --
#  >> a quarter of the 8 m it replaces, and only 8 cells at dy = 0.25 m.
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
#  200 x 160 = 32000 cells, on 56 ranks. Memory arithmetic: see the SBATCH block.
#
#  LINEAR ON PURPOSE. Its nonlinear twin diverged at t=5.52 of 52 s, just after
#  the ramp, with every |u3x|>0.5 node at x = 0.25-0.50 and dmp/int of 3-5. Establishing
#  that the WIDE geometry carries a clean linear sea state is the prerequisite
#  for re-testing the nonlinear case; running the nonlinear one first would
#  confound the geometry change with the instability under investigation.
#
#  OUTDIR IS SET EXPLICITLY: the driver's tag carries regime/tier/bed/Hs but not
#  the geometry, so this would otherwise overwrite the 20 m-wide linear run.
#SBATCH --job-name="BALFEM_lin_directional_wide"
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
#  (8,7) on 200x160: 25 x 22.9 cells per rank = 6.25 x 5.71 m subdomains,
#  aspect 1.09 — the most isotropic split of 56 (rule 22). (14,4) would give
#  3.57 x 10 m, aspect 0.36, and far more halo for the same rank count.
export BALFEM_PX=8
export BALFEM_PY=7            # 8*7 = 56 ranks

# --- Case-specific overrides ---
export BALFEM_REGIME=linear
export BALFEM_NL_PRESSURE=none
export BALFEM_FLAT_BED=1
export BALFEM_LY=40.0         # was 20.0
export BALFEM_NY=160          # keeps dy = 0.25 m
export BALFEM_SPONGE_Y=2.0    # was 8.0 -> undamped strip 4 m -> 36 m (4.0 lambda_eff)
#  Halve the VTK write frequency: snapshot buffers add to peak RSS and this run
#  has only ~6 % headroom. 130 snapshots over 52 s is still ample for animation.
export BALFEM_SAVE_EVERY=20
export BALFEM_OUTDIR=$HOME/GridapBALFEM.jl/output/small_directional_linear_none_flat_Hs0.2_WIDE_P1LFE-2

balfem_run 56 examples/distributed_small/run_directional_sea_small.jl
