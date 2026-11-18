# Distributed algebraic LFE-M example runs (cluster-ready)

MPI-parallel wave-propagation cases for the **algebraic (stacked-layout) 2D LFE-M solver**
(`GridapBABALFEM.jl/src/GridapBABALFEM.jl`), all writing the **full field set** to VTK: surface elevation
`eta`, per-σ-node horizontal velocity components `u<j>x`/`u<j>y`, and the reconstructed
**vertical velocity** `w_s<σ>` and **total pressure** `p_s<σ>` at every vertical σ-node
(`write_w`/`write_pressure` on by default).

Unlike the old solver's distributed examples, these run the **full nonlinear physics through the
single Gridap path** (`setup_and_run_distributed`, hand Jacobians + GMRES+Jacobi+Newton) —
no owned V⊗H loop, no linear-core restriction. New-physics flags are exposed as env vars:
`BABALFEM_PFULL` (full `P¹L¹+P²L²+P³L³` leading pressure) and `BABALFEM_NLP68` (nonlinear-pressure
components 6–8).

## Cases

| Script | Physics |
|--------|---------|
| `run_plane_wave_dist.jl` | long-crested plane wave, line-source wavemaker, long flume |
| `run_periodic_plane_wave_dist.jl` | plane wave in a **periodic-width flume** (`y_wall_bc=:periodic` — top/bottom edges identified, no lateral sponges) |
| `run_ring_wave_dist.jl`  | radial ring wave, point-source wavemaker, square basin, 4-side sponge |
| `run_ic_hump_dist.jl`    | Gaussian hump released from rest, **closed basin** (`x_wall_bc=true` — mandatory for IC problems) |
| `run_bathymetry_dist.jl` | shoaling of a plane wave over a smooth submerged bar (variable `d(x)`, slope-pressure package on) |
| `run_irregular_sea_dist.jl` | **long-crested JONSWAP sea via Dirichlet boundary generation** (WaveSpec.jl, no wavemaker; seeded phases → rank-identical component table) |
| `run_directional_sea_dist.jl` | **short-crested JONSWAP × cosine-power spreading via Dirichlet BCs** (η, 𝖴x AND 𝖴y prescribed; `y_wall_bc=:open` + lateral sponges) |

Each runs at **any** vertical resolution (`BABALFEM_M`) and any core count. All knobs are
environment variables with sensible defaults — see `_dist_common.jl` (shared) and each
script's header (case-specific).

## Launch

Use Julia's own MPI launcher (system `mpiexec` fails here with a PMIx mismatch). **The MPI rank
count `-n` MUST equal `BABALFEM_PX · BABALFEM_PY`**, and keep `BABALFEM_NX` divisible by `BABALFEM_PX`,
`BABALFEM_NY` by `BABALFEM_PY`.

```bash
# From the project root. Plane wave, M=2, 8 ranks (8×1 grid):
BABALFEM_M=2 BABALFEM_PX=8 BABALFEM_PY=1 \
  ~/.julia/bin/mpiexecjl --project=. -n 8 julia --project=. \
  GridapBABALFEM.jl/examples/distributed/run_plane_wave_dist.jl

# Ring wave at M=3, 16 ranks (4×4), bigger mesh:
BABALFEM_M=3 BABALFEM_PX=4 BABALFEM_PY=4 BABALFEM_NX=400 BABALFEM_NY=400 \
  ~/.julia/bin/mpiexecjl --project=. -n 16 julia --project=. \
  GridapBABALFEM.jl/examples/distributed/run_ring_wave_dist.jl

# Quick local smoke of the plumbing (2 ranks, tiny mesh, 1 period):
BABALFEM_PX=2 BABALFEM_PY=1 BABALFEM_NX=60 BABALFEM_NY=4 BABALFEM_LX=60 BABALFEM_LY=4 \
BABALFEM_PERIODS=1 BABALFEM_SAVE_EVERY=40 BABALFEM_XWM=12 BABALFEM_SPONGE=12 \
  ~/.julia/bin/mpiexecjl --project=. -n 2 julia --project=. \
  GridapBABALFEM.jl/examples/distributed/run_plane_wave_dist.jl

# Irregular JONSWAP sea via Dirichlet generation, 128 ranks:
BABALFEM_M=2 BABALFEM_PX=32 BABALFEM_PY=4 BABALFEM_HS=0.002 BABALFEM_TP=1.6 BABALFEM_NFREQ=21 \
  ~/.julia/bin/mpiexecjl --project=. -n 128 julia --project=. \
  GridapBABALFEM.jl/examples/distributed/run_irregular_sea_dist.jl

# Directional (short-crested) sea, 128 ranks:
BABALFEM_M=2 BABALFEM_PX=16 BABALFEM_PY=8 BABALFEM_NTHETA=7 BABALFEM_SPREAD_STD=20 \
  ~/.julia/bin/mpiexecjl --project=. -n 128 julia --project=. \
  GridapBABALFEM.jl/examples/distributed/run_directional_sea_dist.jl
```

### Dirichlet-generation guidance (the two `*_sea_dist` scripts)

- **Spectral band inside the model:** choose `BABALFEM_FMAX_FAC` so the shortest bin stays inside
  the vertical model's applicable band (`kd(fmax) ≲ kd_app(M)`: 10.9/39.2/127.9 for
  LFE-2/3/4) AND is resolved by ≥ 6 cells per wavelength. Out-of-band components fall back
  to the Airy wavenumber with a warning and are simply not propagated correctly.
- **Bin spacing:** `BABALFEM_SAMPLING=uniform` (default) keeps bins well separated, so offline
  per-component DFT analysis is leakage-free over long windows; `energy` (equal-energy bins)
  concentrates resolution at the peak but packs bins tightly — prefer it only for spectrum-
  shape studies analysed with Welch PSDs (postprocessing `psd_welch`).
- **Reproducibility:** the phase seed `BABALFEM_SEED` makes the realisation deterministic —
  identical on every rank (WaveInput snapshots the seeded phases into plain arrays) and
  across reruns; change it for ensemble statistics.
- **Amplitude regime:** keep component amplitudes within the linear-stability rule
  (`Hs/2 ≲ 0.001` scaled) for long fully nonlinear runs.
- **Absorption:** never place a plain sponge on the generation side (it damps the incident
  wave — the driver warns); use `BABALFEM_RELAX=1` for a generation/absorption relaxation zone
  instead when re-reflection at the inflow matters.

> `MPI_Finalize` prints a benign OFI/WiFi-NIC error and exits 143 on the development machine —
> the computation completes correctly beforehand. Julia buffers stdout when redirected to a
> file; the drivers `flush` at every print interval so logs can be followed live.

## On a SLURM cluster

Ready-made launchers for these scripts live in **`run/`** (one per case). They all run against the
prebuilt system image via the shared helper `run/babalfem_env.sh` — without it every rank JIT-compiles
the FEM stack (~30–45 min and a ~4–8 GB/rank memory spike). Build the image once
(`cd compile && sbatch compile_snellius.sh`), then submit:

```bash
sbatch run/run_planewave.sh      # or run_ringwave.sh, run_irregularsea.sh, …
```

To write a new launcher, copy the nearest one and change only the `#SBATCH` header, the `BABALFEM_*`
overrides, and the `babalfem_run` target:

```bash
#!/bin/bash
#SBATCH --job-name=babalfem_plane
#SBATCH --partition=rome        # the sysimage removes the compile spike that needed fat_rome
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=08:00:00
#SBATCH --output=babalfem_plane_M%a.out
#SBATCH --array=2,3             # job array: one run per M (SLURM_ARRAY_TASK_ID)
                                # (no --mem-per-cpu: take the node default, 2 GB/core on rome)

source $HOME/GridapBABALFEM.jl/run/babalfem_env.sh

export BABALFEM_M=$SLURM_ARRAY_TASK_ID
export BABALFEM_PX=8; export BABALFEM_PY=4     # 8*4 = 32 = --ntasks
export BABALFEM_NX=4000; export BABALFEM_NY=200
export BABALFEM_SAVE_EVERY=25
export BABALFEM_OUTDIR=$HOME/GridapBABALFEM.jl/output/plane_M${BABALFEM_M}

babalfem_run 32 examples/distributed/run_plane_wave_dist.jl
```

`babalfem_run <nranks> <script.jl>` (path relative to the project root) loads the cluster modules the
image was built with, verifies the image, and expands to the full
`mpiexecjl --project=… -n … julia --project=… -J<sysimage> …` invocation. Set `BABALFEM_NO_SYSIMAGE=1`
to fall back to the JIT path while the image is stale or rebuilding. Details: `compile/README.md`.

## Validation against the old solver

The sequential algebraic solver is oracle-equivalent to `LFE-M_2D_solver` at machine precision
(`GridapBABALFEM.jl/test/test_equivalence.jl`); the distributed path is validated by
`GridapBABALFEM.jl/test/test_basic_distributed.jl` (4 ranks, linear + fully nonlinear, sequential
agreement within GMRES tolerance). For a cross-solver check at scale, run the old
`run_plane_wave_dist.jl` and this `run_plane_wave_dist.jl` with identical env settings and
compare `eta` snapshots/gauge sections in ParaView (identical field naming).
