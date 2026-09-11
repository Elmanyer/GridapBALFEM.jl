# RUNNING.md — how to run the solver, locally and on a cluster

Related: [`CONFIGURATION.md`](CONFIGURATION.md) (what the knobs mean and cost),
[`TEST_SUITE.md`](TEST_SUITE.md) (running the tests).

---

## 1. Minimal usage

```julia
using GridapBALFEM          # a real package — never include()
using LinearAlgebra         # dot/norm are NOT re-exported

diags, vert, prob = setup_and_run(
    M = 2, d_val = 3.5, T_wave = 1.6, A_wave = 0.001,
    regime = :nonlinear, nl_pressure = :native, flat_bed = true,
)
```

Distributed, from an `mpiexecjl -n 4` context:

```julia
diags, vert, prob = setup_and_run_distributed(cpu_grid = (2,2), M = 2, d_val = 3.5, ...)
```

**MPI must be launched with `~/.julia/bin/mpiexecjl`** — the system `mpiexec` fails with a PMIx
version mismatch. `MPI_Finalize` prints a benign OFI error on this machine and exits 143.

---

## 2. Examples

| directory | contents |
|---|---|
| `examples/` | sequential: `plane_wave`, `ring_wave`, `periodic_plane_wave`, `bc_plane_wave`, `bc_irregular_sea` (JONSWAP + gauge CSV), `bc_directional_sea` |
| `examples/distributed/` | 6 env-configurable cluster scripts (plane/ring wave, IC hump, bathymetry, irregular sea, directional sea) + `_dist_common.jl` with `build_airy_state()` |
| `examples/distributed_small/` | **5 parametric** small-domain (50×20 m) scripts grouped by generation type: `run_periodic_plane_small` (line source), `run_ring_small` (point source), `run_bc_plane_small` (`:bc_gen`), `run_directional_sea_small`, `run_irregular_sea_small` |
| `examples/validation/` | physical benchmarks: `stokes_harmonics`, `submerged_bar`, `solitary_wave`, `ring_spreading`, `bichromatic_sideband`, plus `dispersion_sweep.jl` and `spectral_fidelity.jl` |
| `examples/local_1d/run_flume_1d.jl` | the parametric **quasi-1D flume** — `BALFEM_WAVE_GEN ∈ {inner,bc,sea}` |
| `examples/local_2d/run_small_2d.jl` | the parametric small 2-D case (25×10 m, 6 ranks as 3×2), `BALFEM_WAVE_GEN ∈ {line,point,bc,sea}` |
| `examples/local_mms/` | MMS convergence drivers |
| `examples/inspect_run.jl` | stdlib-only reader turning a run's `diagnostics.csv` into a health verdict (works on cluster output) |

> **The solver is structurally 2-D.** The "quasi-1D flume" is a narrow domain with `ny ≥ 3` and
> y-periodic sides — that is how a 1-D horizontal problem is posed here. It is **not** a 1-D model.

Each parametric script bakes in the common geometry/mesh/partition/solver defaults and reads the
physics from the environment with `get!` base defaults, so the banner and the solver can never
disagree. The output directory is auto-tagged by configuration so cases do not clobber each other.

---

## 3. Local launchers (`run/local/`, ≤ 6 cores)

`balfem_local.sh` is the helper — project resolution, a hard rank cap, `balfem_local_run`
(sequential) / `balfem_local_mpi` (MPI), and exit-143 handling.

* **7 × `run_1d_*.sh`** — **sequential**, because that was measured faster than any MPI split at
  20 k DOFs (`CONFIGURATION.md` §6) and it restores point gauges.
* **8 × `run_2d_*.sh`** — 12-rank MPI.
* `run_all_1d.sh` runs the 1-D set side by side; `bench_solver_config.sh` compares
  preconditioner × `ls_rtol` configurations on one fixed case, reporting iterations, wall time
  **and** `max|η|`.
* `run/local/README.md` documents how to read `diagnostics.csv`.

Deliberately shares no code with `run/balfem_env.sh`, which is cluster-only. The point is a
minutes-long feedback loop instead of a multi-hour one.

---

## 4. Cluster launchers (`run/`, `run/dist_small/`)

> ### ⚠ JOB SIZING AND COST — read [`run/SNELLIUS_ROME_LAUNCH_CONFIGS.md`](../run/SNELLIUS_ROME_LAUNCH_CONFIGS.md) first
>
> A Snellius `rome` node is **128 cores / 224 GiB, divided into EIGHTHS** (1/8 = 16 cores / 28 GiB).
> Slurm bills `max(cores/128, memory/224)` **rounded up to the next eighth**, evaluating cores and
> memory separately — **the larger fraction wins**. BALFE-M needs ≈4 GiB/rank (measured peak
> 3633 MB), which is well above the 1.75 GiB/core a node provides pro-rata, so **memory always sets
> the tier** and the cost-optimal rank count is **`7k` for tier `k/8`**:
>
> | tier | ranks @4G | | tier | ranks @4G |
> |---|---|---|---|---|
> | 1/8 | 7 | | 5/8 | 35 |
> | 2/8 | 14 | | 6/8 | 42 |
> | 3/8 | 21 | | 7/8 | 49 |
> | 4/8 | 28 | | 8/8 | 56 |
>
> All 38 launchers were resized to these bins on 2026-09-07 (32→28, 48→42, 40→35, 16→14), each
> saving one eighth. **`--mem-per-cpu=5G` is never valid on rome** — it pushes almost every rank
> count over a boundary. The full rules, worked examples and traps are in that file.


**All 29 launchers run against the prebuilt sysimage**, through the shared helper
**`run/balfem_env.sh`**: it loads the cluster modules the image was built with, resolves and
verifies `GridapBALFEM_sysimage.so`, and exposes

```bash
balfem_run <nranks> <script.jl>     # wraps: mpiexecjl --project=… -n … julia --project=… -J<sysimage> …
```

A launcher is therefore **only** its `#SBATCH` header, its `BALFEM_*` overrides, and one
`balfem_run` call.

| helper knob | effect |
|---|---|
| `BALFEM_PROJ` | project path (default `$HOME/GridapBALFEM.jl`) |
| `BALFEM_CLUSTER` | `snellius` / `blue` |
| `BALFEM_SYSIMAGE` | explicit image path |
| `BALFEM_NO_SYSIMAGE=1` | escape hatch — drop `-J` and JIT-compile (for a stale image or a mid-rebuild) |
| `BALFEM_STRICT_SYSIMAGE=1` | upgrade the staleness warning to a hard abort (use for long production jobs) |

**A missing image fails the job immediately** with the build command, rather than silently falling
back to a ~45 min/rank compile.

**Contents:** 9 production cases in `run/` (plane / ring / periodic-plane wave, IC hump, bathymetry,
irregular + directional sea, plus the DelftBlue `run_blue.sh`), and **20 small-domain
observation/comparison cases** in `run/dist_small/`, each overriding only what changes (regime /
`nl_pressure` / `flat_bed`→bar / amplitude / period). Small-domain geometry: 50×20 m, `d=3.5`,
`T=2.0` (`kd≈3.5`), `save_every=10`, partition **8×4 = 32 ranks** (plane/irregular, 200×40) or
**8×8 = 64** (directional 200×80; ring 160×160 on a 40×40 m domain — square cells throughout,
because GMRES iteration count grows with mesh anisotropy).

**Resources:** every Snellius launcher runs on `rome` with `--mem-per-cpu=4G`, sized to fit a node
(256 GB / 128 cores), so ranks are spread `nodes × ntasks-per-node` and no launcher requests more
than 256 GB/node. The 4 GB/core is **required, not compile headroom** — see `OPEN_ISSUES.md` §1.

---

## 5. The sysimage (`compile/`)

Builds a precompiled system image so ranks **load** the solver instead of JIT-compiling it,
removing the ~30–45 min/rank compile.

* `set_preferences.jl` pins MPI.jl to the **system** OpenMPI via `use_system_binary()`.
* `compile.jl` bakes the deps and traces `warmup.jl`, with an MPI preflight and
  **`include_transitive_dependencies=false`** — load-bearing: it stops PackageCompiler force-loading
  `OpenMPI_jll`, whose baked initializer would `dlopen` the JLL artifact and crash the image at
  launch (`undefined symbol: opal_single_threaded`).
* `compile_snellius.sh` / `compile_blue.sh` run the chain; walkthrough in `compile/README.md`.

**Verify a built image:**

```bash
strings GridapBALFEM_sysimage.so | grep -i OpenMPI_jll   # MUST be empty
strings GridapBALFEM_sysimage.so | grep libmpi           # only the system libmpi.so
```

**Staleness detection.** The image bakes a *compiled copy* of the solver, so editing `src/*.jl`
without rebuilding runs old code silently. `compile_*.sh` calls `balfem_write_sysimage_stamp` after
a successful build, writing `GridapBALFEM_sysimage.so.src.sha256` (a hash of all `src/*.jl`);
`balfem_check_sysimage_freshness` recomputes and compares before every launch. Edits, additions and
deletions trip it; a bare `touch`, a re-clone or a content-restoring `git checkout` do not. Images
predating the stamp fall back to an mtime comparison (coarser, can cry wolf).
**Rebuilding remains manual — staleness is detected, not prevented.**

---

## 6. Environment variables

| variable | selects |
|---|---|
| `BALFEM_P_VERT` | vertical polynomial order (default 1); each parametric script derives `model_name = "P$(p_vert)LFE-$(M)"` for its banner and output directory |
| `BALFEM_P_ETA` | surface FE order. **Default `BALFEM_FE_ORDER − 1` (Taylor-Hood).** `check_taylor_hood` RAISES on any other pairing — equal order is inf-sup deficient and caused the 2026-09 instability (`CLAUDE.md` rules 2b, 12b) |
| `BALFEM_NX/NY`, `BALFEM_PX/PY` | mesh and partition — `NX % PX == 0`, `NY % PY == 0`, `-n == PX·PY` |
| `BALFEM_QUAD_EXTRA` | extra horizontal quadrature degree on top of the default `2·max(p_h,p_η)+2`. **0 is the default and is one degree SHORT of the degree-7 nonlinear advection integrand at `Q2/Q2`** — see `CLAUDE.md` rule 12b |
| `BALFEM_SOLVER`, `BALFEM_TABLEAU` | integrator |
| `BALFEM_NL_ITER`, `BALFEM_NL_TOL` | Newton |
| `BALFEM_LS_MAXITER`, `BALFEM_LS_RTOL`, `BALFEM_KRYLOV_M` | linear solve |
| `BALFEM_DIAG_EVERY`, `BALFEM_DIV_FACTOR` | diagnostics and the divergence guard |
| `BALFEM_WAVE_GEN` | generation mechanism in the parametric local scripts |
| `BALFEM_TESTS` | test tier (`fast` / `default` / `all` / explicit files) |
| `MMS_NL_LEVELS` | shorten the nonlinear MMS studies |

---

## 7. Postprocessing

```julia
# separate environment: postprocessing/
using GridapBALFEMPost
sim = load_simulation("output/<case>")     # solution.pvd + CSV → WaveSimulation
```

See `ARCHITECTURE.md` §7 for the module inventory. `postprocessing/README.md` is the user guide.
