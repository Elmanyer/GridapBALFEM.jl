# GridapBALFEM.jl — the algebraic BALFE-M wave solver (v2)

A loop-free, sequential-and-distributed Gridap implementation of the **BALFE-M**
(*Basis-Agnostic Layer-integrated Finite Element*, `M` vertical elements) family of depth-integrated,
non-hydrostatic wave models. It generalises the **LFE-M** model of Yang & Liu (2024, *J. Fluid Mech.*
999, A32) from a fixed piecewise-linear vertical basis to an arbitrary one. A concrete member is named
**P*p*LFE-*M*** after the basis it uses, so P1LFE-2 is the two-layer piecewise-linear model (the
counterpart of Yang & Liu's LFE-2); *LFE-M* is reserved for their published models.

Repository: `GridapBALFEM.jl` (branch **`main`**, the v2 solver since 2026-11-24; v1 is frozen on `v1-solver`,
tag `v1_final_solver`). Module: `src/GridapBALFEM.jl` → `module GridapBALFEM`.

> **Sister solver.** This repository is the primal formulation stabilised by the C⁰ interior penalty
> (CIP). The mixed formulation of the same models (conforming without CIP) with residual-based (VMS)
> stabilisation is a separate solver, [`GridapBALFEMixed.jl`](https://github.com/Elmanyer/GridapBALFEMixed.jl),
> split from this repository on 2026-11-23 (history: tag `archive/mixed-vms-solver-pre-split`).

**Where things are written down.**

| question | file |
|---|---|
| the map, physics switches and standing rules | [`CLAUDE.md`](CLAUDE.md) |
| where the solver stands, what comes next | [`markdown_files/STATUS.md`](markdown_files/STATUS.md) |
| the open stability problem | [`markdown_files/STABILITY.md`](markdown_files/STABILITY.md), [`INTEGRATOR_STUDY_RESULTS.md`](markdown_files/INTEGRATOR_STUDY_RESULTS.md) |
| everything else that is open | [`markdown_files/OPEN_ISSUES.md`](markdown_files/OPEN_ISSUES.md) |
| every other document | [`markdown_files/INDEX.md`](markdown_files/INDEX.md) |
| the mathematics | the LaTeX project `latex_docs/GridapBALFEM_documentation/` (Overleaf, gitignored here; `CLAUDE.md` §2), condensed in [`MODEL.md`](markdown_files/MODEL.md) |

---

## The core idea

The model represents the horizontal velocity as a sum of vertical basis functions,
`u(x, σ, t) = Σⱼ uⱼ(x, t) φⱼ(σ)`, so every term of the depth-integrated equations carries a sum over
the vertical (layer) index. `w` and the non-hydrostatic pressure are eliminated analytically into
precomputed σ-tensors; there is no pressure Poisson solve. Instead of one scalar FE field per layer
and component (`1 + 2Nσ` fields) with Julia loops in the residual, this solver **stacks the layer
index into the FE value type**:

```
MultiField = [η, 𝖴x, 𝖴y]        (exactly 3 fields, not 1+2Nσ)
𝖴x, 𝖴y  ∈  VectorValue{Nσ}       (all Nσ layers' x-/y-velocity in one field value)
```

Every vertical sum becomes one algebraic operation on `VectorValue{Nσ}` / `TensorValue{Nσ,Nσ}` /
`ThirdOrderTensorValue{Nσ,Nσ,Nσ}` constants:

| model sum | stacked Gridap expression |
|---|---|
| `Σⱼ Mᵢⱼ u̇ⱼ` | `𝗠 ⋅ Ut` (matvec) |
| `Σⱼ Φⱼ uⱼ` (depth average) | `𝚽 ⋅ Ux` |
| `Σₖⱼ 𝓜ᵢₖⱼ (uₖ·∇uⱼ)` | `double_contraction(𝗠3, Ux⊗∂x(Ux) + Uy⊗∂y(Ux))` |

The residual has no vertical-index loops, is well-typed by construction, and assembles **~3.6× faster
with ~13× fewer allocations** than the per-layer form. Because `Operation` is forwarded for
`DistributedCellField`, the same residual serves sequential and MPI runs
([`ARCHITECTURE.md`](markdown_files/ARCHITECTURE.md) §1; LaTeX chapter 7).

---

## Physics implemented

The complete global residual of the LaTeX weak form (chapters 6–7):

- mass continuity, with an optional interior Gaussian wavemaker source;
- `H`-weighted acceleration, nonlinear horizontal and vertical advection (`𝓜`, `𝓖`);
- gravity in the integrated-by-parts energy form (rest-state exact on any bed);
- the **leading pressure `R_P`**, which is the model's entire frequency dispersion;
- the linear bed-slope / surface-slope non-hydrostatic package (`𝓐/𝓚`);
- the **nonlinear pressure `𝓝`, all eight components** when `nl_pressure = true`:
  - components {3,6,7,8} assembled directly;
  - {1,2,4,5}: the ∇h half by exact integration by parts, the Class-III half by cellwise Hessians
    completed by the kind-S facet terms on the skeleton, with an exact Jacobian (`src/broken.jl`);
  - **the first-order C⁰-CIP on η and every 𝖴 is part of the model whenever `𝓝` is on** (default
    γ = 0.3, `p_eta ≥ 2`); the code refuses `𝓝` without it.
- boundaries: solid walls, open, x- and y-periodic, quadratic sponges (damping velocity **and** η),
  generation/absorption relaxation zones;
- wave generation, `wave_gen`: `:inner_res` (interior source: line ⇒ plane wave, point ⇒ ring wave) or
  `:bc_gen` (time-dependent Dirichlet data: a regular wave, a `WaveInput`, or a **WaveSpec `AiryState`**
  sea state), with the discrete-eigenmode (`:model`) or Airy (`:airy`) vertical polarisation and a Hann
  start-up ramp ([`WAVE_GENERATION.md`](markdown_files/WAVE_GENERATION.md));
- time integration (`build_time_integrator`): `:theta` (θ = 0.5 is Crank–Nicolson), `:sdirk` / `:rk`
  (any Gridap tableau; default `SDIRK_2_2`) and `:gen_alpha` (ρ∞). Unusable tableaus are refused.

### Selecting the model — three orthogonal switches

| control | values | meaning |
|---|---|---|
| `regime` | `:linear` \| `:nonlinear` | linearised core / full nonlinear core with advection |
| `nl_pressure` | `false` \| `true` | `𝓝` off / all eight components (needs `:nonlinear` and `p_eta ≥ 2`; carries the mandatory CIP) |
| `flat_bed` | `Bool` | `∇h ≡ 0`: every ∇h term dropped at one control point |

They give six models:

| model | regime | bed | `𝓝` | status |
|---|---|---|---|---|
| 1 | linear | flat | — | MMS optimal order, 1-D and 2-D, sequential and distributed |
| 2 | linear | variable | — | MMS optimal order |
| 3 | nonlinear | flat | off | MMS to `nx = 32` |
| 4 | nonlinear | variable | off | MMS to `nx = 32` |
| 5 | nonlinear | flat | on | outside the verified MMS scope; grid-scale unstable in long box runs |
| 6 | nonlinear | variable | on | as model 5 |

The verification record and its limits: [`VERIFIED_SCOPE.md`](markdown_files/VERIFIED_SCOPE.md) and
[`OPEN_ISSUES.md`](markdown_files/OPEN_ISSUES.md) §1, §4. The stability of models 3–6 in the closed
periodic box is the open problem ([`STABILITY.md`](markdown_files/STABILITY.md)).

---

## Repository layout

```
GridapBALFEM.jl/
├── src/                     the package (20 files): vertical.jl, tensors.jl, horizontal.jl,
│                            problem.jl (residual + hand Jacobians), nlpressure.jl, broken.jl
│                            (Class III + facet terms + CIP), timeloop{,_dist}.jl,
│                            utilities{,_dist}.jl (setup_and_run*), configuration.jl, checkpoint.jl,
│                            monitor.jl, reconstruct.jl, waveinput.jl, vopt.jl, mms*.jl, errors.jl
├── test/                    runtests.jl (tiers fast | default | all; verdict from gate output),
│                            test_*.jl, local/ (quasi-1-D machinery gates, smoke_integrators.jl),
│                            cluster/, v2_migration/regression_snapshot.jl, yl_collapse_wip/
├── examples/                plane / ring / periodic / bc_* sequential cases, inspect_run.jl,
│                            local_1d/ (flume, closed periodic box, eigen-analyses), local_2d/,
│                            local_mms/, distributed/, distributed_small/, validation/
├── run/                     balfem_env.sh (cluster helper), storage.sh + setup_storage.sh,
│                            run_*.sh (production cases), linear_production/, integrator_study/,
│                            local/ (workstation launchers, balfem_local.sh)
├── compile/                 cluster sysimage build (Snellius, DelftBlue, local)
├── postprocessing/          GridapBALFEMPost: its own environment; VTK/CSV analysis, growth spectra,
│                            3-D reconstruction for ParaView
├── markdown_files/          the written record (tracked); archive/ holds v1 and completed plans
├── Project.toml / Manifest.toml
└── CLAUDE.md
```

Present locally but gitignored: `output/` (v2 runs; on a cluster a link to the project storage),
`latex_docs/` (a link to the LaTeX project), `Gridap.jl/` (the Gridap fork, below) and
`WaveSpec.jl/` (the vendored sea-state package, GitHub version).
The v1 output is no longer here: since 2026-11-26 it is at the workspace root, `../output_v1/`,
indexed by `../output_v1/CATALOGUE.md` (every v1 run labelled with a reconstructed `CONFIGURATION.txt`).

---

## Environment

`GridapBALFEM.jl/` is a Julia package (`name = "GridapBALFEM"`) **and** the working environment. Load
it by name, never with `include`:

```julia
# julia --project=/path/to/GridapBALFEM.jl
using GridapBALFEM
```

Being a package is what lets PackageCompiler bake the solver and its Gridap specialisations into the
cluster system image, and gives sequential runs a cached precompile.

* **Versions:** Gridap 0.20.8 (fork), GridapDistributed 0.4.17, GridapSolvers 0.7.1,
  PartitionedArrays 0.3.5, MPI 0.20.26, WaveSpec 0.3.2; Julia 1.10.
* **⚠ Gridap is a fork** (`Elmanyer/Gridap.jl`, branch `fix-transient-multifield-ad`, based on
  v0.20.8), loaded from `Gridap.jl/` with `Pkg.develop`. It carries two AD patches: transient
  multifield (`fa860899c`) and skeleton integrals on transient fields (`3928b98b9`, not yet pushed).
  Gate A0 of `test/test_jacobians_ad.jl` checks the first. AD is an oracle for the Jacobians only,
  never for the residual. Details: [`CONFIGURATION.md`](markdown_files/CONFIGURATION.md) §2.
* **WaveSpec.jl** must track its GitHub version: the tagged release's `change_seed!` reads a
  non-existent field and breaks every sea-state run.

---

## Quick start

```julia
using GridapBALFEM

# Sequential: nonlinear plane wave from an interior source over a flat bed
diags, vert, prob = setup_and_run(
    M           = 2,                            # vertical elements (P1LFE-2)
    domain      = ((0.0, 60.0), (0.0, 10.0)),
    partition   = (120, 20),
    p_u         = 2,                            # p_eta = p_u − 1 (Taylor–Hood) is enforced
    h_val       = 3.5, T_wave = 1.6, A_wave = 0.01,
    x_wm        = 12.0,                         # line source (y_wm = nothing)
    sponge_wL   = 12.0, sponge_wR = 12.0, mu_max = 5.0,
    T_final     = 12.8, dt = 0.02,
    regime      = :nonlinear,
    nl_pressure = false,                        # true: all eight 𝓝 components (+ mandatory CIP, p_u ≥ 3)
    flat_bed    = true,
    solver_type = :theta,                       # CN for diagnosis; default :sdirk (SDIRK_2_2) is dissipative
    gauges      = [(30.0, 5.0)],
    save_every  = 20,
)
```

Distributed runs use `setup_and_run_distributed(cpu_grid = (px, py), …)` with the same physics
keywords, launched with Julia's own MPI launcher:

```bash
~/.julia/bin/mpiexecjl --project=. -n 4 julia --project=. my_script.jl
```

Every run writes, into a directory named by `output_dir_name` (it never overwrites), a
`run_manifest.toml` (machine record), a `CONFIGURATION.txt` (readable identity card, including the
integrator's carrier damping), `diagnostics.csv` and VTK snapshots. Long runs checkpoint
(`checkpoint_every`) and restart (`restart_from`). The closed x-periodic box used for all stability
work is `examples/local_1d/run_periodic_1d.jl`.

---

## Running

* **Workstation:** `run/local/` launchers built on `run/local/balfem_local.sh`; parallel cases are
  admitted by free memory (`_slot_throttle.sh`). See [`RUNNING.md`](markdown_files/RUNNING.md).
* **Cluster (Snellius, DelftBlue):** every SLURM launcher goes through `run/balfem_env.sh`, which loads
  the modules, uses the prebuilt system image (`compile/`, rebuild after any `src/` edit: the helper
  warns on a stale image) and **refuses to run unless output goes to the project storage**:
  `STORAGE` exported in `~/.bashrc`, and the repository's `output/` a link to
  `$STORAGE/GridapBALFEM_output` (one-time `run/setup_storage.sh`). Launcher sets: `run/run_*.sh`
  (production cases), `run/linear_production/`, `run/integrator_study/` (each with a README). Memory:
  `--mem-per-cpu=4G` is required. Details: [`RUNNING.md`](markdown_files/RUNNING.md) §4–5,
  [`compile/README.md`](compile/README.md).

---

## Tests and verification

```bash
julia --project=. test/runtests.jl                    # default tier (fast + medium)
BALFEM_TESTS=fast julia --project=. test/runtests.jl  # seconds
test/local/run_local_tests.sh                         # quasi-1-D machinery gates
```

`runtests.jl` judges each file by its **gate output**, never by its exit code. The distributed tests
(e.g. `test_nlp_distributed_parity.jl`, models 5–6 sequential ↔ distributed) run by hand under
`mpiexecjl` or on PartitionedArrays' debug backend.

Only the **analytic MMS** (`src/mms.jl`, `test/test_mms_*.jl`; its forcing never touches
`problem.jl`, enforced by a grep gate) verifies the residual. Everything else is self-consistency or
physics agreement. The one external oracle is the collapse of BALFE-M at `p = 1` onto Yang & Liu's
LFE-M (`test_yl_collapse.jl`). The test inventory and what each gate cannot detect:
[`TEST_SUITE.md`](markdown_files/TEST_SUITE.md) (scores there are v1's; the v2 rebaseline is
`V2_SOLVER_PLAN.md` step 8). Known failures: [`OPEN_ISSUES.md`](markdown_files/OPEN_ISSUES.md) §2.

---

## Design notes for anyone extending this solver

The full list is the standing rules in [`CLAUDE.md`](CLAUDE.md) §7. The ones that bite first:

- **Taylor–Hood is mandatory** (`p_eta = p_u − 1`, enforced): equal order is inf-sup deficient and
  produced the v1 "nonlinear instability".
- **`fe_order ≥ 2`**: Q1 zeroes `R_P` and all non-hydrostatic physics.
- **`B_stored = −B̃`** is negative definite; the explicit `(−1)` factors are load-bearing.
- **`∂R/∂u̇` is exact; the nonlinear `∂R/∂u` is quasi-Newton by choice.** An omission is benign only if
  it is higher order in amplitude; measure with `test/test_jacobians_ad.jl`.
- **Never apply `∇` to an `Operation`-composed expression containing a test basis**; expand by hand.
- **Dissipation is a legitimate stabiliser only when named and measured.** `SDIRK_2_2` (the default) is
  A-stable but not L-stable and damps the carrier; diagnose under Crank–Nicolson and report every
  dissipative remedy with its carrier cost (rule 15).
- **Any CellField a skeleton term reads on both sides must include ghosts** in distributed runs.
- **Distributed:** `NewtonSolver(GMRESSolver(Pr = Jacobi))`, never `lu()`; initial conditions with
  `interpolate_everywhere`; keep `ConsecutiveMultiFieldStyle`; solid-wall Dirichlet BCs include the
  corner tags.
- **The system image bakes a compiled copy of `src/`**: rebuild it after any edit, or launch with
  `BALFEM_NO_SYSIMAGE=1`.
