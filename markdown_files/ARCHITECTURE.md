# ARCHITECTURE.md — solver structure and workflow

**What this file is.** How the code is organised and how a run flows through it: the stacked-value
design decision that shapes everything, the file map, the FE spaces and boundary conditions, the
time loops, the distributed path, and the runtime instrumentation.

Related: [`MODEL.md`](MODEL.md) (the maths), [`CONFIGURATION.md`](CONFIGURATION.md) (settings and
measured performance), [`WAVE_GENERATION.md`](WAVE_GENERATION.md) (sources and boundaries).

---

## 1. The core design decision: stack the layer index into the value type

Everything turns on one decision: **the vertical (layer) index lives in the FE value type, not in a
Julia array.**

* Velocity is **two vector-valued fields** `𝖴x, 𝖴y ∈ VectorValue{Nσ}` (all layers' x- resp.
  y-velocity). The MultiField is `[η, 𝖴x, 𝖴y]` — **3 fields**, not `1+2Nσ`.
* The static vertical arrays become **constant tensors**: `M^V → TensorValue{Nσ,Nσ}`,
  `𝓜^V, 𝓖^V → ThirdOrderTensorValue{Nσ,Nσ,Nσ}`; the component-indexed `A^V/K^V/𝓐^V/𝓚^V` are split
  per component into 3 (resp. 8) such constants.
* Consequently **every layer sum `Σ_j`, `Σ_{kj}` is a matvec or tensor double-contraction** — the
  residual contains no vertical-index loops.
* The spatial index is only 2-D and is written explicitly: `e_x=(1,0)`, `e_y=(0,1)`,
  `∂_x f ≡ e_x⋅∇f`.
* The MultiField is touched only for `η=U[1]`, `𝖴x=U[2]`, `𝖴y=U[3]` — **no `Nσ`-decomposition**.

**Why.** Expressing the layer sums as native tensor contractions makes the residual well-typed by
construction: no `VectorValue{3}`-of-vectors rank mismatches, no illegal `U[2:end]` slices of a
`TransientMultiFieldCellField`, and no gradient-of-a-divergence or `∂²η` on the trial space. It is
also fast — the stacked integrand measured **~3.6× faster with ~13× fewer allocations** than the
fused per-layer form — and the same CellField algebra is forwarded transparently to
`DistributedCellField`, so **one residual and one pair of Jacobians serve both sequential and MPI
execution**. That is what removed the old solver's need for a separate hand-written owned-loop
distributed path.

### Verified Gridap facts (rely on them)

1. `∇` of a `VectorValue{Nσ}` field is `TensorValue{2,Nσ}` — **spatial index first**,
   `(∇f)[d,j] = ∂f_j/∂x_d`. So `∂_x f = e_x⋅∇f`. **Use `e⋅∇f`, not `∇f⋅e`.**
2. `double_contraction(𝓣::ThirdOrderTensorValue, S::TensorValue{Nσ,Nσ})` contracts the **trailing
   two** indices → `VectorValue{Nσ}` = `Σ_{k,j} 𝓣_{ikj} S_{kj}`. Native.
3. `W::VectorValue{Nσ} ⋅ 𝓣::ThirdOrderTensorValue` contracts the **first** index →
   `TensorValue{Nσ,Nσ}`. Native.
4. `⊗`, `⊙` (Frobenius), and `Operation(VectorValue)(a,b)` (build a `VectorValue{2}` from two scalar
   CellFields) are all native.

**Consequence: the residual needs no custom contraction primitives.** The only hand-written helpers
are the build-time constructors turning assembled Float arrays into constant tensors (`src/tensors.jl`).

### Variable / operator dictionary

| symbol | type | meaning | Gridap |
|---|---|---|---|
| `η, H` | scalar CellField | free surface, total depth `H = d+η` | `U[1]`; `H = d_cf + η` |
| `𝖴x, 𝖴y` | `VectorValue{Nσ}` | stacked layer velocities | `U[2]`, `U[3]` |
| `𝖶x, 𝖶y` | `VectorValue{Nσ}` | test functions | `V[2]`, `V[3]` |
| `d_cf` | scalar | still-water depth `d(x)`; `∇h = ∇d`, `∇²h` analytic | `CellField(d_func, Ωₕ)` |
| `𝚽` | `VectorValue{Nσ}` const | depth-average weights `Φ_j` | `alg_to_vec(Φvec)` |
| `𝗠` | `TensorValue{Nσ,Nσ}` const | vertical mass `M^V` | `alg_to_tensor2(M2)` |
| `𝗠3, 𝗚3` | `ThirdOrderTensorValue` const | advection, index `[i,k,j]` | `alg_to_tensor3(…)` |
| `A[c], K[c]` | 3 × `TensorValue{Nσ,Nσ}` | linear pressure per component | component split |
| `P[c]` | 3 × `TensorValue{Nσ,Nσ}` | leading pressure (`P[3] = −𝗕`) | `alg_to_tensor2(vert.P[:,:,c])` |
| `𝗔3[c], 𝗞3[c]` | 8 × `ThirdOrderTensorValue` | nonlinear pressure per component | component split |
| `DU`, `DW` | `VectorValue{Nσ}` | per-layer divergence of trial / test | `∂x(𝖴x)+∂y(𝖴y)` |
| `S` | `VectorValue{Nσ}` | `∇·(H u_j) = H·DU + u_j·∇H` | product rule |
| `ū, W̄` | `VectorValue{2}` | depth-averaged velocity / test | `Operation(VectorValue)(…)` |

**Coding rule.** Never apply `∇` to an `Operation`-composed expression containing a test basis
(block-array `copyto!` is not implemented). Expand by hand using linearity of the vertical
contraction in the test: `∂_a(W⋅𝓣) = (∂_aW)⋅𝓣`, and the collapse `Σ_a Ψ⋅∂_aU_a = Ψ⋅DU`. Same policy
for state expressions — obtain the bed Hessian once via `∇∇(d_cf)` and compose afterwards.
See `src/nlpressure.jl`.

---

## 2. `src/` file map

| file | contents |
|---|---|
| `GridapBALFEM.jl` | module entry: deps, includes, exports; re-exports `WaveSpec` |
| `vertical.jl` | σ-mesh pre-compute — all vertical tensors incl. `P`, `Pcal` (`MODEL.md` §2). Plus `assemble_dispersion_tensors` — the **`(Φ, Mmat, B)` subset only**, for σ-mesh optimisation: the full path assembles `3·8·N⁴` extra integrals that dispersion never touches, and an optimiser calls its objective thousands of times |
| `tensors.jl` | constant-tensor constructors + `Operation` helpers |
| `horizontal.jl` | mesh + stacked FE spaces (distributed-safe dispatch; transient-Dirichlet inflow variants) |
| `problem.jl` | `BALFEMProblem`, `resolve_physics`, `global_residual`, `jacobian_u`, `jacobian_u_t` |
| `nlpressure.jl` | the eight 𝓝 components: native set, exact-IBP ∇h half, frozen-projection half |
| `timeloop.jl` | sequential ODE operator, solver factory, time loop, VTK |
| `timeloop_dist.jl` | distributed mesh builder, GMRES+Jacobi+Newton, distributed time loop |
| `utilities.jl` | `setup_and_run` (the sequential driver), sponge, sources, dispersion helpers. `resolve_cbdy` is now the ONE place σ-element boundaries are chosen (shared with the distributed driver and the MMS drivers). Linear wave properties added 2026-08-21: `model_R` (R, R′, R″), `airy_R`, `wave_properties`, `property_errors`, `applicable_range` — `C_g` and `γ`, which the codebase previously lacked entirely |
| `utilities_dist.jl` | `setup_and_run_distributed` (`with_mpi` wrapper) |
| `waveinput.jl` | Dirichlet boundary wave generation + WaveSpec coupling |
| `reconstruct.jl` | `w_s<σ>` / `p_s<σ>` VTK reconstruction (`MODEL.md` §9) |
| `monitor.jl` | `SolverMonitor`, `RunDiagnostics`, `ResidualChecker` (§6) |
| `mms.jl` | the analytic MMS forcing — **must never reference `problem.jl`** (grep-gated) |
| `mms_driver.jl` | `run_mms_case`, `run_conv_study`, `run_mms_case_distributed`. All three take the vertical basis `(M, p_vert, c_bdy)` as parameters and the distributed one takes the model switches — both were hard-wired until 2026-08-21 (`COMPLETED_VBASIS_STUDY.md` §1 prereq A) |
| `errors.jl` | error norms and convergence-rate fitting |

**The solver is a real Julia package** (`name = "GridapBALFEM"`,
`uuid = 43e94d05-4d7d-4679-96a4-d46e2615da34`) and is loaded with **`using GridapBALFEM`**, never
`include()`. This is load-bearing for the cluster: an `include()`d module lives in a throw-away
`Main.GridapBALFEM`, so neither the solver nor — far more expensive — the Gridap FEM specialisations
keyed on its types survive into a sysimage, and every rank recompiles them. As a package it also
precompiles natively (~3.5 s to load). If you add a `using X` to any `src/*.jl`, add it to
`[deps]` too.

This directory is **both the package and the working environment** — tests, examples and the compile
tooling run directly against it — so `Test`, `BlockArrays`, `MPIPreferences` and `Preferences` are
kept in `[deps]` rather than `[extras]`.

---

## 3. FE spaces and boundary conditions (`horizontal.jl`)

```julia
reffe_H = ReferenceFE(lagrangian, Float64, orderH)                  # η
reffe_U = ReferenceFE(lagrangian, VectorValue{Nσ,Float64}, orderU)  # stacked layer velocity
X = MultiFieldFESpace([V_H, V_Ux, V_Uy])
```

* **`p_u ≥ 2` is mandatory** — `Q1` zeroes the `R_P` dispersion term. (`fe_order`/`p_horizontal`
  were renamed to `p_u` on 2026-09-06; the surface order is `p_eta`.)
* **The pairing must be Taylor-Hood, `p_u = p_eta + 1`**, enforced by `check_taylor_hood` in
  `src/horizontal.jl` — called from `build_fe_spaces` and from both `setup_and_run` drivers
  before any expensive setup. Equal order is inf-sup deficient here (`CLAUDE.md` rule 2b).
* **`p_eta`** selects a lower surface order (Taylor–Hood-style pairing). It **defaults to equal
  order**, so nothing changes until opted into; exposed through both drivers and the run scripts as
  `BALFEM_P_ETA`. See `VERIFIED_SCOPE.md` §3 before switching production.
* Keep **`ConsecutiveMultiFieldStyle`** — `BlockMultiFieldStyle` breaks Jacobi's `diag`.

**Boundary conditions:**

| control | values | effect |
|---|---|---|
| `y_wall_bc` | `:wall` / `:open` / `:periodic` | solid wall `𝖴y=0` / natural / y-edges identified |
| `x_wall_bc` | `Bool` | solid x-walls |

Periodicity is a property of the **mesh topology**, not of the FE-space Dirichlet data, so the
symbol splits into two orthogonal actions: the mesh builder takes `isperiodic=(false, y_periodic)`,
and at FE-space level `:periodic` behaves exactly like `:open` (no y-Dirichlet — the mesh already
identifies the DOFs). The residual, Jacobians and time loop are untouched.

**Rules that are not optional:**

* **Solid-wall Dirichlet BCs must include the corner tags.** Omitting them leaves corner DOFs
  unconstrained and the run diverges exponentially.
* **IC-release problems need `x_wall_bc=true`.** A free x-wall together with the dispersion term
  forms a spurious-forcing mode that an initial perturbation excites directly, so closed-basin
  initial-condition cases run with solid x-walls.
* **`ny ≥ 3` is mandatory for a y-periodic mesh** (Gridap `CartesianGrids.jl:39`); `ny=1,2` are
  rejected at mesh construction.
* **Distributed initial conditions** use `interpolate_everywhere([0.0 for _ in 1:n_fields], U)`.
  `FEFunction(U, zeros(…))` creates a local array and fails distributed.

---

## 4. Time integration

`build_ode_solver` / `build_ode_solver_distributed`:

* **Default: `RungeKutta(nls, ls, dt, :SDIRK_2_2)`** — fully implicit, L-stable, 2nd order,
  diagonally implicit; robust in the stiff deep-water regime.
* `solver_type=:theta` selects Crank–Nicolson.

> ⚠ **The default integrator is DISSIPATIVE by construction.** Any test measuring a
> non-dissipative property (energy conservation, amplitude transfer, temporal order) **must pin
> `solver_type=:theta`**. See `TEST_SUITE.md` §4 for the four tests that do and the measured
> difference.

**Transient API** (Gridap 0.20): `TransientFEOperator(res, jac, jac_t, U, V)`, `res(t,u,v)` with
`∂t(u)`, `TransientCellField` in `Gridap.ODEs`, `solve(solver, op, t0, tF, u0)` yielding `(t, uh)`.

---

## 5. The distributed path

**One code path.** The residual and both Jacobians are byte-identical between sequential and
distributed execution — `Operation` is forwarded for `DistributedCellField`, FE bases and
`TransientCellField` are wrapped per part, and `CellField(f, trian)` / `createvtk` / `createpvd`
have distributed methods.

What differs:

| aspect | sequential | distributed |
|---|---|---|
| linear solve | `LUSolver` (direct) | `GMRESSolver(krylov_m; restart=true, maxiter, Pr=JacobiLinearSolver())` |
| nonlinear solve | `NLSolver(...; method=:newton)` | `NewtonSolver(...)` (GridapSolvers) |
| `max|η|` | direct reduction | `own_values(PVector)` + `reduce(max, …; init=0.0)` |
| frozen-projection mass solve | `lu()` | `CGSolver(JacobiLinearSolver())` |
| gauges | point evaluation available | not available (inter-rank point search) |

**Conventions that must not be rediscovered:**

1. A direct LU factorisation does **not** scale to partitioned matrices at cluster size.
2. `norm(PVector, Inf)` is broken in PartitionedArrays 0.3.5 — use the own-values reduction.
3. **`krylov_m` and `ls_maxiter` are different bounds.** `GMRESSolver`'s first *positional* argument
   is `m`, the Krylov **basis size** — a memory bound, allocating `m+1` distributed vectors plus a
   dense `(m+1)×m` Hessenberg per rank, up front. `maxiter` is a *keyword* iteration budget whose
   library default is **100**. `restart=true` is load-bearing: the default `restart=false` lets the
   basis grow past `m`, i.e. unbounded memory. Symptom of getting this wrong: **`gmres=` pinned at
   exactly the same number every step** with Newton needing 8–24 iterations instead of 3–5.
4. Frozen-projection RHS/solution vectors must be allocated **from the matrix**
   (`allocate_in_range`/`allocate_in_domain` + in-place `assemble_vector!`): an independently
   assembled vector is only isomorphic to the matrix's `PRange`, not identical, and `solve!`'s
   internal `mul!` asserts exact equality.
5. `BALFEM_NX` divisible by `BALFEM_PX`, `BALFEM_NY` by `BALFEM_PY`, and `-n == PX·PY`.
6. All rank-0-only printing behind `i_am_main(ranks)`; `mkpath` on rank 0 then `MPI.Barrier`.
7. Launch with `~/.julia/bin/mpiexecjl` — the system `mpiexec` fails with a PMIx version mismatch.
   `MPI_Finalize` prints a benign OFI error on this machine and exits 143.
8. Julia buffers stdout when redirected to a file; use explicit `flush` or poll the file.

---

## 6. Runtime instrumentation (`src/monitor.jl`)

Overhead measured **within noise** (−0.8 % sampling every step, +2.2 % every 10th).

* **`SolverMonitor`** — a transparent `NonlinearSolver` wrapper (pass via `monitor=`). Per step:
  Newton iterations, initial→final residual, convergence flag, GMRES counts, nonlinear-solve wall
  time. Under Runge–Kutta it accumulates over stages (`nl_iters` is the sum, `ncalls` the stage
  count). Tracks **linear-solver saturation** (`lin_min`/`lin_max`, `lin_sat` flag → WARN). The
  banner is built from values read back **out of the constructed solver objects** (`ls.m`,
  `ls.log.tols.maxiter`, …), not from the caller's kwargs, so a mis-passed argument shows on line 1.
* **Field diagnostics** (`RunDiagnostics`, `field_diagnostics`) — sampled every `diag_every` steps
  and written to `<output_dir>/diagnostics.csv`:
  * **where** `max|η|` sits (`x_at_max`) and its **interior vs damped-zone** split;
  * `max|u|` and `|u|/|η|`;
  * mass `∫η` and energy, with drift against the **t=0** state (baselines seeded from `u0`, so a
    forced run's gain is visible);
  * per-rank RSS, current and peak, reduced over ranks;
  * a **relative divergence guard**: abort at `div_factor · eta_ref` (default 20×), `eta_ref`
    inferred from the forcing via `resolve_eta_ref`.
* **`ResidualChecker`** — every `check_every` steps the governing equations are reassembled through
  a separate code path: (a) the θ-scheme discrete residual (meaningful only for a single-stage
  scheme, so it runs **only under `solver_type=:theta`**; `res_theta = NaN` under SDIRK), and (b) the
  instantaneous PDE residual, i.e. the `O(Δt)` time-discretisation error, for any integrator.

Driver kwargs: `diag_every`, `diag_csv`, `eta_ref`, `div_factor`; env `BALFEM_DIAG_EVERY`,
`BALFEM_DIV_FACTOR`. Read a run with `julia --project=. examples/inspect_run.jl <output_dir>`
(stdlib-only, works on cluster output).

**Three diagnostic-interpretation rules** (each learned from a measurement that misled first):

1. **Never read `growth` without `x_at_max`.** Rising `eta_max` with a *moving* `x_at_max` is a
   filling domain; with a *pinned* `x_at_max` at a boundary it is the spurious boundary mode.
2. **`dmp/int` means "≫1 is trouble", never "<1 is fine".** Healthy band measured **0.29–0.96**; the
   boundary mode gives ≈65.
3. **Mass drift is an invariant only in a closed, unforced basin** (1.2e-17 there). In a forced run
   it measures *injected* mass and should scale with `A` — confirmed, 10× amplitude gave exactly 10×
   the drift.

---

## 7. Output and postprocessing

**VTK** per component (`eta, u1x, u1y, …`) plus reconstructed `w_s<σ>` / `p_s<σ>` fields, indexed by
`solution.pvd`. VTK filenames replace the decimal point with an underscore
(`sol_t_1_6000.vtu`) to avoid a WriteVTK extension warning.

**`postprocessing/GridapBALFEMPost`** is a self-contained library with **its own environment**
(ReadVTK, Plots+GR, FFTW, Interpolations — pinned separately from the solver, and with **no
dependency on the solver**). It reads `solution.pvd` / `sol_t_*.vtu` + CSV into a `WaveSimulation`,
auto-`regularize!`ing the duplicated Q2 node cloud onto a Cartesian grid. Modules: `io, probes,
spectral, diagnostics, reconstruct, plotting, seastate` — gauges/DFT/celerity/harmonics/radial/
conservation; heatmap/animation/Hovmöller/dispersion/profile plots; Welch PSD, JONSWAP overlay,
spectral moments/Hs, zero-upcrossing heights, Rayleigh exceedance. `reconstruct.jl` rebuilds
`w(σ)`/`p_nh(σ)` from the stored velocity modes at any σ using the analytic σ-basis and Gauss
quadrature (no Gridap), matching the solver's own `w_s` to 4–8 %.
