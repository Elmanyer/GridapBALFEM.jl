# CONFIGURATION.md — the current solver configuration and what it costs

**What this file is.** The settings the solver runs at, the evidence behind each, and the measured
performance envelope. Everything here is a measurement, not a preference.

Related: [`ARCHITECTURE.md`](ARCHITECTURE.md) (what the settings configure),
[`RUNNING.md`](RUNNING.md) (how to launch), [`OPEN_ISSUES.md`](OPEN_ISSUES.md).

---

## 1. Dependency stack

| package | version |
|---|---|
| Gridap | **0.20.x** (a fork — §2) |
| GridapDistributed | 0.4.17 |
| GridapSolvers | 0.7.1 |
| PartitionedArrays | 0.3.5 |
| MPI | 0.20.26 |

`[compat]` admits **both** minors (`Gridap = "0.19, 0.20"`, `GridapSolvers = "0.6, 0.7"`) **on
evidence, not as a hedge**: `test_basic` gives *identical* results under `0.19.11 + 0.6.2` and
`0.20.8 + 0.7.1` — max η 0.00410 m, gauge 0.00212 m, same Newton count — and the package precompiles
cleanly on both. The Gridap minor has no effect on this solver's results.

This environment and the cluster both run 0.20.x / 0.7.1. The parent `../` is a **separate**
environment for the legacy solvers and stays on 0.19.11 / 0.6.2; it does not need to match.

**A sysimage is only valid for the versions it was built against — build it in the environment you
run in.**

> **Identifying a cluster's versions without shell access:** Julia package directory slugs are
> content-addressed and stable across machines, so a cluster traceback names its versions
> (`GridapSolvers/WuCdi` = 0.7.1, `PartitionedArrays/MVmxR` = 0.3.5, `MPI/pvbg6` = 0.20.26).

**`WaveSpec.jl` is vendored and tracks the GitHub repository version, not a tagged release** — the
release's `change_seed!` reads a non-existent `state.spec` field and breaks every sea-state run.
It is `Pkg.develop`ed in both environments and re-exported by `GridapBALFEM`.

Calling scripts must `using LinearAlgebra` explicitly — `dot`, `norm` etc. are not re-exported.

---

## 2. The Gridap fork (transient multifield AD)

Gridap is pinned to **`Elmanyer/Gridap.jl` @ `fix-transient-multifield-ad`**, one commit on the
`v0.20.8` tag, `Manifest.toml` pins it.

**Root cause of the bug it fixes** — it was never `ForwardDiff.Dual`.
`time_derivative(::TransientMultiFieldCellField)` (`src/ODEs/TransientCellFields.jl:213`) builds its
third constructor argument with `map(cellfield, derivatives...)`. **`map` returns a `Tuple` when fed
a Tuple and a `Vector` when fed an array-like MultiField container**, while the struct field is
declared `transient_single_fields::Vector{<:TransientCellField}`. The hand path passes a
`MultiFieldFEFunction` (array-like ⇒ Vector ⇒ works); the AD path passes a plain Tuple ⇒
`MethodError: no applicable constructor`.

**The fix is one additive outer constructor** — purely additive, no existing call site changes
behaviour, and a wider dispatch surface survives upgrades better than editing a function body:

```julia
function TransientMultiFieldCellField(
  cellfield, derivatives::Tuple, transient_single_fields::Tuple
)
  TransientMultiFieldCellField(
    cellfield, derivatives, collect(TransientCellField, transient_single_fields)
  )
end
```

**Consequences.** `use_ad=true` works, so AD differentiates the same assembled residual and is a
usable **oracle for the hand Jacobians** — but for the Jacobians only (`VERIFIED_SCOPE.md` §5). The
hand Jacobians remain the fast production path.

**The standing guard that the pin has not rotted is gate A0 of `test/test_jacobians_ad.jl`** — an
instant `hasmethod` check for the fork's Tuple-argument constructor, which stock Gridap does not
have. Before it, nothing in the repo exercised the fork at all, so the pin could have rotted
silently and the loss would only have surfaced mid-debugging.

---

## 3. Time integrator

**Default `RungeKutta(:SDIRK_2_2)`** — fully implicit, L-stable, 2nd order. Robust in the stiff
deep-water regime, and **dissipative by construction**: see `TEST_SUITE.md` §4 for the four tests
that must pin `solver_type=:theta` and the measured difference.

Driver kwargs (both drivers): `solver_type=:sdirk`, `tableau=:SDIRK_2_2`, `nl_iter=50`,
`nl_tol=1e-5`, distributed `ls_maxiter=1000` / `krylov_m=100` / `ls_rtol=1e-5`, `print_every`,
`check_every=50`, `check_tol=1e-8`.

---

## 4. The tolerance ladder

**Defaults: `ls_rtol = 1e-5`, `nl_tol = 1e-5`.** No run script overrides them. The two behave
completely differently, and conflating them is the mistake this section exists to prevent.

### `ls_rtol` is free

One fixed 2-D case (96×36, `:full`, 12 ranks, 30 steps):

| `ls_rtol` | GMRES/solve | s/step | Newton/step | `max\|η\|` |
|---|---|---|---|---|
| 1e-9 (old default) | 491–515 | 227.5 | 3.17 | 3.065280e-03 |
| 1e-7 | 358–363 | 206.1 | 3.17 | 3.065280e-03 |
| **1e-6** | **294–314** | **191.9** | **3.17** | **3.065280e-03** |
| **1e-5** | **238–253** | — | unchanged | unchanged |

A 40 % then a further 19 % cut in linear iterations for an answer identical to seven significant
figures. **The control that matters is `Newton/step`, and it did not move.** Relaxing a linear
tolerance frequently just relocates work — Newton gets worse steps and needs more iterations. Here
it demonstrably does not, which means the extra linear accuracy was being *discarded*.

### `nl_tol` is NOT free — it is a step function

Newton runs an integer number of iterations, so only threshold crossings matter. At 1e-6 it needs
**3** iterations/step; at 1e-5 **and** 1e-4 it needs **2**, and those two are **bit-identical to 12
digits**. Dropping that third iteration moves `max|η|` by **1.9e-5 (2-D) / 3.8e-5 (`test_basic`)**
relative. **Do not describe this as "no effect".**

`nl_tol = 1e-5` is adopted on an **error budget, not a null result**: the dropped iteration polishes
a `~3e-6` residual, while the `O(Δt²)` time-discretisation error already in the answer is
`‖R‖∞ ≈ 1.8e-3` (measured by the run's own residual checker) — **~600× larger**. Tests needing a
sharper answer pin `nl_tol=1e-8`.

Note the adopted pair leaves **no separation** between `ls_rtol` and `nl_tol`, deliberately. The
older "keep the linear solve one order tighter" rule is sound general guidance but is not what
governs accuracy here, and `ls_rtol` was measured not to limit Newton anywhere in 1e-9…1e-5.

> **Reference values that move with tolerance:** `test_basic`'s Newton count is **240** at
> `nl_tol=1e-5` (it was 408 at 1e-6). max η and gauge amplitude are unchanged.

---

## 5. The linear solve is the entire cost — and it is a preconditioner problem

Newton is at its floor: **1–2 iterations per RK stage everywhere, zero non-convergences in 13 runs**.
GMRES meanwhile needs **392–774 iterations per solve**.

| case | DOFs | ranks | GMRES/solve |
|---|---|---|---|
| 2-D BC plane | 98.6 k | 12 | 392–457 |
| 2-D plane (lin & nl) | 98.6 k | 12 | 451–517 |
| 2-D `:full` over bar | 98.6 k | 12 | 524–556 |
| 2-D ring (solid lateral walls) | 98.6 k | 12 | **654–695** |
| quasi-1-D flume | 20.2 k | 2 / 4 / 6 / 12 | **758 / 758 / 758 / 715–774** |

Three independently diagnostic facts:

1. **Rank-independent** — 758 iterations at 2, 4, 6 *and* 12 ranks. Over-decomposition would grow
   with rank count; this does not.
2. **The smaller problem needs MORE iterations** — 20 k DOFs → ~760, 99 k → ~480. The difference is
   **mesh anisotropy**: the flume has 4:1 cells (`dx=0.25, dy=1.0`) where the 2-D mesh is isotropic.
   Point-Jacobi degrades badly on stretched elements.
3. **Boundary conditions matter** — solid lateral walls cost ~40 % more than periodic.

> **Design rule: keep horizontal cells near-isotropic, or pay for it in the linear solve.**

Headroom is thinner than it looks: `ls_maxiter = 1000` against a measured 774 is **1.3×**.

**Both cheap preconditioner drop-ins are ruled out by measurement:**

* **Additive Schwarz** dies with `SingularException` — `SchwarzLinearSolver(LUSolver())` factorises
  each rank's local block, but `partition(::PSparseMatrix)` returns own+**ghost** rows with the
  ghost rows unassembled, i.e. structurally zero.
* **Symmetric Gauss–Seidel** completed **zero steps in 12 h**, against 1.9 h for a *complete* Jacobi
  run.

What remains, in order of effort: **field-split/Schur** (medium) → **geometric multigrid** (high) →
**restricted Schwarz** (needs library support GridapSolvers 0.7.1 does not expose). The
**ring-vs-plane gap on an identical mesh is an isolated, cheap, reproducible test bed** — a
preconditioner that closes it likely helps everywhere.

`run/local/bench_solver_config.sh` compares configurations on one fixed case and reports iterations,
wall time **and** `max|η|`, so a cheaper setting is rejected if it changes the answer. It brackets
the adopted `ls_rtol=1e-5` rather than searching for it, and the two measured-dead preconditioners
are opt-in behind `BENCH_FAILED_PRECOND=1`.

---

## 6. Partitioning must follow the problem, not the core count

Strong scaling on an identical 240×3 mesh, 20 202 DOFs, 40 steps/point:

| ranks | 1 (direct LU) | 2 | 4 | 6 | 12 |
|---|---|---|---|---|---|
| s/step | **6.45** | 14.44 | 13.10 | 19.53 | 18.27 |
| GMRES | — | 758 | 758 | 758 | 715–774 |

**Sequential beats every MPI configuration by 2–3× at that size**, and more ranks do not help. For
the 98.6 k-DOF 2-D problem the ordering reverses: 12 ranks (55 s/step) beat 6 (132), 4 (151), 2 (143).

> **Use MPI when the problem is big enough to give each rank a meaningful share, and a direct LU
> otherwise. Spend spare cores on running more *cases*, not on decomposing one small case further.**

Sequential also restores point gauges, which the distributed driver lacks. The 1-D launchers are
sequential for exactly this reason.

---

## 7. Measured costs (this workstation, 16 cores)

| quantity | measured |
|---|---|
| `using GridapBALFEM` | **2.8 s** (native precompile) |
| first `setup_and_run` in a session | ≈ **207 s** of JIT |
| warm, 4 669 DOFs, sequential | 0.349 s/step |
| warm, 20 202 DOFs, sequential LU | 6.45 s/step |
| warm, 98 600 DOFs, 12 ranks | 55 s/step (linear) |
| 2-D `:full` over a bar | 385 s/step |
| diagnostics overhead | −0.8 % … +2.2 % — **within noise** |

**Concurrency is memory-bandwidth bound, not core bound.** Six concurrent Julia processes buy
≈**2×** throughput, not 4–6×; three at once each run at ~⅓ solo speed. **Budget `JOBS = 2–4`.**

**Memory: the baseline dominates at this problem size.** RSS across 13 local cases: **1.4 GB at
first sample → 2.5 GB typical → 3.75 GB peak**. A Julia + Gridap + GridapBALFEM process costs
**≈1.4 GB before solving anything**, and within a run memory is essentially flat (+2 % over 400
steps) — so **no per-step leak is supported at this scale**. Cluster consequence: a 2 GB/core budget
would leave ~0.6 GB/rank for the computation. `--mem-per-cpu=4G` is in every launcher and a
node-default 2 GB/core job was OOM-killed; see `OPEN_ISSUES.md` §1.

---

## 8. Physics settings that are not free choices

* ~~**`A_wave ≤ 0.001 m`** keeps fully nonlinear runs stable over long integrations.~~
  ✅ **WITHDRAWN 2026-09-06.** The cap was the equal-order instability in disguise, not a property
  of the model. On Taylor-Hood (now enforced) the 1-D flume runs 50 wave periods at
  **`A_wave = 0.10 m`** (`κa ≈ 0.16`, `kd = 5.5`) with the mean `η` flat to the fourth decimal,
  sitting +5.2 % above the linear control — the second-order Stokes crest elevation. That is 100×
  the old cap. ⚠ Physical limits (Miche, breaking) still apply; the model has no breaking closure.
  See `CLAUDE.md` rules 2b and 12b.
* **`nl_pressure=:native` is the production tier.** The whole `{1,2,4,5}` hierarchy contributes
  **0.013 % (1-D) / 0.094 % (2-D)** on top of advection's 0.77 % / 1.84 % at `A=1e-3` — the most
  expensive part of the residual changing the answer in the fourth significant figure. `:full` also
  carries a mesh-independent velocity-error floor (`VERIFIED_SCOPE.md` §4).
* **Sponge strength saturates: past `μ_max ≈ 5ω`, WIDTH is the lever, not strength.**

  | `μ_max` | 5 | 20 | 40 |
  |---|---|---|---|
  | `R²` of the damping law | 0.9997 | 0.9996 | 0.9991 |
  | `rate/μ_max` | 0.537 | 0.495 | 0.479 |

  The envelope follows `ln a ∝ −μ_max(x−x_R)³/(3w²c_g)` with excellent fidelity, but the constant
  falls 11 % across `μ_max/ω = 1.27 → 10.2`. The law is derived from `ω → ω + iμ`, valid only for
  `μ ≪ ω`; beyond that the sponge becomes an evanescent barrier and extra stiffness buys
  progressively less absorption. Reflection is small and does **not** grow with strength (1.0–1.4 %).
  **The constant cannot be validated here:** a weak-damping probe (`μ_max = 1`, `μ/ω = 0.25`)
  *diverges*, because a sponge that weak no longer holds the η-dominated boundary mode —
  **the WKB regime and open-boundary stability are mutually exclusive at a free outflow.**

  > **Reading rule this cost a measurement:** a decay-ratio sequence that accelerates and then
  > abruptly flattens (`1.04 1.5 2.8 6.1 17.9 2.1`) has hit a **noise floor**; it has not changed
  > physics. Exclude saturated stations before fitting.

* **Sponge width must cover the LONGEST component** (`kd_min` ⇒ `λ_max`), not the peak.
* **The `c_g` transit trap is the standing hazard of the deep-water regime.** At `kd = 5.5`,
  `c_g = 1.25 m/s`, so filling a 45 m flume takes 36 s = 22.5 periods. **Budget
  `t_settle ≈ (x_sponge − x_source)/c_g + 3T` before reading any steady state.** It produced
  misleading measurements three separate times; the 1-D driver's default duration is now
  transit-aware (16 periods interior, **26** boundary/sea).

---

## 8b. Horizontal quadrature degree — a default that is short for the nonlinear terms

`src/utilities.jl:685` sets the horizontal measure to `Measure(trian, 2*max(p_h, p_η) + 2)`, with a
comment claiming this "integrates the nonlinear (product) terms exactly enough". **For the linear
terms it does; for the nonlinear ones it does not.** At the common `Q2/Q2` pairing the rule is
degree 6, while the advection integrand `φᵢ·(u_k·∇u_j)·H` counts

| `φᵢ` | `u_k` | `∇u_j` | `H = d+η` | total |
|---|---|---|---|---|
| 2 | 2 | 1 | 2 | **7** |

Confirmed directly on a one-cell integral of `x⁸` (exact `1/9`):

| degree | value | error |
|---|---|---|
| 6 | 0.111088435374 | 2.27e-05 |
| 10 | 0.111111111111 | 1.39e-17 |
| 14 | 0.111111111111 | 2.78e-17 |

`quad_extra` (kwarg on `setup_and_run`, `BALFEM_QUAD_EXTRA` in the 1-D driver) adds to this degree.

⚠ **It was tested as a cure for the nonlinear instability and REFUTED** (2026-09-05,
`CLAUDE.md` rule 12b): against an exact control, degrees 6, 10
and 14 give growth curves agreeing to **four decimal places**. The under-integration is real but
**dynamically irrelevant** — do not describe raising the degree as a stability measure.

The deficit is still worth correcting on **correctness** grounds, and it is cheap to justify: because
raising the degree changes only how exactly the residual is integrated and not the residual itself,
**no recorded MMS rate is invalidated** by the change. The default is nonetheless left at `0` for now
— the cost is real (quadrature points per cell rise steeply with degree) and the benefit is a
fourth-decimal correction, so it should be changed deliberately with a cost measurement, not as a
side effect of this investigation.

---

## 9. Workflow hazards worth knowing

* **Launchers re-export their own environment.** Exporting `BALFEM_PERIODS=1` *before* calling a
  launcher does nothing — the launcher's own `export` wins. A short probe must bypass the launcher
  and call the script directly.
* **Long single evaluations lose their output.** A 30-minute-plus REPL/MCP call is aborted
  client-side while the Julia process continues, and the PASS/FAIL lines are lost. Run suites
  through their runners, which give each test its own process and log. `diagnostics.csv` is written
  continuously and survives such aborts.
* **Revise does not hot-swap signature changes.** Adding a keyword argument silently produces
  pre-edit results until the session is restarted. After any signature or struct change, restart
  before trusting a number.
