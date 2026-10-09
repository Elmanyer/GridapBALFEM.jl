# CLAUDE.md — `GridapBALFEM.jl/` · the 2D BALFE-M algebraic wave solver (v2)

> ## ⇨ START HERE
>
> 1. This file: the map, the physics switches and the standing rules. It stays at the repository root,
>    because Claude Code auto-loads `CLAUDE.md` from the root only.
> 2. **[`markdown_files/STATUS.md`](markdown_files/STATUS.md): where the solver stands and what comes
>    next.**
> 3. [`markdown_files/STABILITY.md`](markdown_files/STABILITY.md) (the open problem) and
>    [`markdown_files/OPEN_ISSUES.md`](markdown_files/OPEN_ISSUES.md) (everything else).
>
> Reference documents: [`INDEX.md`](markdown_files/INDEX.md). Completed plans and full v1 records:
> [`markdown_files/archive/`](markdown_files/archive/README.md). ⚠ `markdown_files/` is **tracked**;
> `latex_docs/` is **gitignored**, so anything load-bearing goes in `markdown_files/`.
>
> **v1 and v2 are separated (2026-10-05).**
>
> | | v1 (frozen) | v2 (active) |
> |---|---|---|
> | code | branches `main`, `v1-solver`; tag `v1_final_solver` | branch **`v2-solver`** |
> | LaTeX | `latex_docs/BALFEM_models_v1/` | `latex_docs/BALFEM_models_v2/` (Overleaf) |
> | outputs | `output_v1/` | `output/` |
> | record | [`archive/HISTORY_V1.md`](markdown_files/archive/HISTORY_V1.md) (v1 `CLAUDE.md`, verbatim) | here + `STATUS.md` |

---

## 0. What this is

**BALFE-M** (*Basis-Agnostic Layer-integrated Finite Element*, `M` vertical elements) generalises Yang &
Liu (2024, *JFM* 999 A32) LFE-M to an **arbitrary vertical FE basis**. It is a depth-integrated,
**non-hydrostatic** free-surface wave model.
* The water column `σ ∈ [0,1]` carries `u_h = Σ_j u_j(x,t) φ_j(σ)`, with `Nσ = M·p + 1` vertical
  nodes and modes `j = 1..Nσ`.
* `w` and `p_nh` are eliminated **analytically** into precomputed σ-tensors. There is no pressure
  Poisson solve.
* What remains is `Nσ+2` coupled 2-D PDEs in `(H, u_1…u_Nσ)`, `H = d + η`, solved by Gridap on a
  horizontal FE mesh: **small dense vertical algebra ⊗ large sparse horizontal FE**.

**Three names, never mixed** (LaTeX `par: nomenclature`):
* **BALFE-`M`**: our family.
* **P`p`LFE-`M`**: a concrete member we run; always `p = 1` so far.
* **LFE-`M`**: Yang & Liu's published models, used only when citing them.

## 1. Repository map

| path | what it is |
|---|---|
| `Project.toml` / `Manifest.toml` | the package `GridapBALFEM`, loaded with **`using GridapBALFEM`, never `include()`**; also the working environment |
| `src/` (18 files) | `problem.jl`: `BALFEMProblem`, `resolve_physics`, the loop-free residual and hand Jacobians · `nlpressure.jl`: `𝓝` {3,6,7,8} and the ∇h half · `broken.jl`: broken Class III + exact Jacobian, the skeleton stabilisers `:jumpgrad`/`:ghostvolume` · `utilities.jl`: `setup_and_run`, `output_dir_name`, `check_v1_env`, `write_run_manifest` · `vopt.jl`, `mms.jl`, `vertical.jl`, … (`ARCHITECTURE.md` §2) |
| `test/` | `runtests.jl` judges by gate output (never exit codes); `test/local/`, `test/cluster/`; `test_skeleton_ad.jl`; `test/v2_migration/regression_snapshot.jl` (v1→v2 gate, recording in `output/v2_baseline/`) |
| `examples/` | `local_1d/` (`run_flume_1d.jl`, `run_periodic_1d.jl` the closed box, `stability_eig.jl`, `cip_eigen_analysis.jl`, `cip_eigen_order_rule.jl`), `local_2d/`, `local_mms/`, `distributed/`, `validation/` |
| `run/` | `local/balfem_local.sh`, `balfem_env.sh` (cluster), SLURM launchers, `SNELLIUS_ROME_LAUNCH_CONFIGS.md`; v2 box launchers `run/local/run_1dper_v2_unstab{,_extra,_q32}.sh` |
| `compile/` | cluster sysimage chain (`RUNNING.md` §5); `set_preferences.jl` adds the forks by URL |
| `output/` · `output_v1/` | gitignored. `output/` = v2 only; every run writes `run_manifest.toml` |
| `postprocessing/` | `GridapBALFEMPost`, its own environment; `examples/periodic_growth.jl` gives the band growth rates |
| `WaveSpec.jl/` | vendored sea-state synthesis (GitHub version) |
| `Gridap.jl/` | the Gridap **fork** (`Elmanyer/Gridap.jl` @ `fix-transient-multifield-ad`, v0.20.8) with two AD patches: transient multifield (`fa860899c`) and skeleton integrals on transient fields (`3928b98b9`, **not pushed**; loaded locally via `Pkg.develop`) |
| `latex_docs/` | `BALFEM_models_v1/` (frozen), `BALFEM_models_v2/` (Overleaf), `CFC2027_abstract/`, `doc_figures/` |

## 2. The LaTeX project

* **`latex_docs/BALFEM_models_v2/` is the single source of mathematical truth.** It is its own git
  repository (remote `overleaf`). **The author edits it in Overleaf**, so always run
  `git fetch overleaf` and fast-forward before editing; the author may stash local edits when
  pulling.
* **Chapter map** (`LATEX_STRUCTURE.md`):
  * model 1–3;
  * checks without a solver 4–5;
  * weak form on the broken domain 6;
  * multi-field implementation 7;
  * stability 8;
  * stabilisation 9;
  * solver functionalities 10;
  * validation 11;
  * appendices A–C, which are v1 records.
* Chapter 6 rests on **one** identity, `eq: broken integration by parts identity`.
* Build in a scratch copy: `pdflatex`, `biber`, `pdflatex` ×2 → 272 pp, 0 errors. One pre-existing
  undefined reference remains.
* ⚠ Never global-replace `H^2`/`L^2`: it hits the depth `H²`, `h²` and `\partial^2` (33 lines to
  restore on 2026-10-07).
* Labels are `\label{chap:|sec:|subsec: name}`; cross-references are `\S\ref{…}`.
* If you change the maths, change `MODEL.md` too.
* Notation bridge (`MODEL.md` §2): LaTeX `M^V, 𝓜^V, 𝓖^V, A^V, K^V, 𝓐^V, 𝓚^V, Φ` ↔ code
  `Mmat, Mcal, Gcal, A, K, Acal, Kcal, Phi/D/C`; the unit basis is stored as `w_j = −φ_j_int`.

## 3. Features

* **Solver core.**
  * The stacked `[η,𝖴x,𝖴y]` loop-free residual with hand Jacobians (exact `∂R/∂u̇`; quasi-Newton
    `∂R/∂u`, rule 5).
  * Integrators: `:sdirk` (`SDIRK_2_2`, default — **dissipative, rule 15**), `:theta` (CN),
    generalised-α, explicit RK.
  * Sequential LU+Newton; distributed GMRES+Jacobi+Newton.
* **Physics.** Advection, the full leading pressure `R_P`, and `𝓝` with all eight components when
  `nl_pressure=true`:
  * `{3,6,7,8}`: direct assembly;
  * `{1,2,4,5}`: the ∇h half by exact integration by parts, and the **Class III** half by the
    **broken formulation** (`broken_class3_residual`), i.e. cellwise `∇∇` plus one skeleton integral
    of the jumps, with an exact Jacobian;
  * the `{1,2,5}` reduction is exact to 4.4e-16.
  * ⚠ The leading-pressure skeleton term `⟦𝒫⟧` is excluded on purpose (`h⁻⁴`-conditioned mass).
* **Skeleton stabilisers** (`attach_skeleton!`), one per run:
  * `:jumpgrad`: normal-derivative jumps of order ≤ `cip_order` ≤ 2;
  * `:ghostvolume`: direct ghost penalty, every order; uniform Cartesian meshes, sequential only.

  Knobs `cip_gamma_u/_eta` (scaled by `d√(gd)`, `√(gd)`, `h^(s−2)`). **Being superseded by
  filtering** (`STATUS.md` §4).
* **Boundaries.** Wavemaker source, Dirichlet generation with WaveSpec `AiryState` sea states and
  the `:model` eigenmode polarisation; sponge (it must damp η); relaxation; walls; x/y periodic.
* **Diagnostics.**
  * Runtime monitor and an independent governing-equation residual check.
  * `w_s`/`p_s` reconstruction.
  * Closed x-periodic box with band-energy growth rates.
  * Frozen-state and penalised-operator eigen-analyses.
* **Verification tools.** Analytic MMS (`src/mms.jl`, grep-gated independence from `problem.jl`;
  models 5–6 opt-in), `test_jacobians_ad`, the Yang & Liu collapse (`test_yl_collapse.jl`, the only
  external oracle), `test_broken_formulation`, `test_skeleton_ad`, the regression snapshot.
* **Linear wave properties** (`src/utilities.jl`): `model_R`, `wave_properties` (`C`, `C_g`, `γ`),
  `applicable_range`. All nine published applicable ranges reproduced to < 1 %.
* **Vertical-grid optimisation** (`src/vopt.jl`): Yang & Liu's total error functional, eq. (3.10),
  exactly (rule 46). `M = 2` is reproduced to 0.0036.
  * `Ω` is an input, never derived; at `p ≥ 2` it is an open choice, so pass `c_bdy`.
  * `resolve_cbdy` gives the published set at `p = 1` and a uniform split at `p ≥ 2`.
* **Postprocessing** (`GridapBALFEMPost`): VTK/CSV analysis, from-modes reconstruction, sea-state
  statistics.

## 4. Physics selection

| control | values | meaning |
|---|---|---|
| `regime` | `:linear` \| `:nonlinear` | linearised core / full nonlinear core with advection |
| `nl_pressure` | `false` \| `true` | `𝓝` off / **all eight components** (requires `:nonlinear`; refuses a `Symbol`) |
| `flat_bed` | `Bool` | `∇h ≡ 0`: every ∇h term dropped |

* Six models:
  * 1–2: linear (flat / variable bed);
  * 3–4: nonlinear, `𝓝` off;
  * 5–6: nonlinear, `𝓝` on.
* `𝓝` is all-or-nothing because no amplitude or `kd` ordering separates its components
  (`V2_SOLVER_PLAN.md` §0).
* Numerical options:
  * `stabilization`, `cip_*`;
  * `solver_type` (**stability claims under `:theta` only**);
  * `use_ad` (diagnosis only).
* Environment variables: `BALFEM_NL_PRESSURE`, `BALFEM_STAB`, `BALFEM_CIP_*`, `BALFEM_SOLVER`,
  `BALFEM_OUTDIR`.
* `check_v1_env()` refuses the v1 knobs.
* Output names carry `nlp0`/`nlp1` (`output_dir_name`, rule 2c).

## 5. Current Implementation Stage

*2026-10-09. Full picture: [`markdown_files/STATUS.md`](markdown_files/STATUS.md).*

**Working and verified.**
* The v2 code base: broken Class III only; Boolean `nl_pressure`; v1 code deleted.
* v1→v2 equivalence entry by entry (12/12).
* Linear models 1–2 at optimal order: 1-D and 2-D, sequential and distributed, five vertical bases.
* Nonlinear models 3–4 MMS to `nx = 32`.
* The Yang & Liu collapse.
* `test_jacobians_ad` 14/14; fast suite 9/9; `test_broken_formulation` 35/35.

**The open problem: the nonlinear discretisation is unstable at the grid scale.**
* The v2 closed-box ladder (CN, A = 0.10) reproduces it:
  * with `𝓝`: Q2/Q1 at 64 cells/λ diverges at 81 s; Q3/Q2 diverges at 82 / 27 / 9 s for
    16 / 32 / 64 cells;
  * **without `𝓝`, too**: Q3/Q2 at 64 cells diverges at 44.6 s.
* **The standard C⁰ interior penalty is not enough.** It removes the fastest growth, but nonlinear
  mid-band modes keep growing (v1: +0.13 s⁻¹ from ≈ 110 s).
* Higher-order penalties held the v1 box. But their order must grow with the element order, which
  Gridap cannot provide beyond 2 (Q4/Q3 has no window), and the ghost penalty is
  Cartesian/sequential/costly.
* **Decision 2026-10-09: filtering, as Yang & Liu do.** Options F1–F4 are in `STATUS.md` §4;
  evidence is in `STABILITY.md`.

**Other open items** (`OPEN_ISSUES.md`):
* nonlinear `p_η` order reduction at Q3/Q2;
* `test_skeleton_ad` ghost arm (`_shift`);
* `test_selfconsistency` step 12;
* the variable-bed lee-shoulder mode;
* models 5–6 outside the verified scope;
* distributed broken path (step 10);
* Gridap fork commit unpushed;
* stale sysimage.

**Committed 2026-10-09 (not pushed):** the documentation restructure, `cip_eigen_order_rule.jl` and
the v2 box launchers. LaTeX chapter 8 still needs the Q3/Q2 rows.

## 6. The design decision everything rests on

**The vertical (layer) index lives in the FE value type, not in a Julia array.**
* Velocity is `𝖴x, 𝖴y ∈ VectorValue{Nσ}`, so the MultiField is `[η, 𝖴x, 𝖴y]`: three fields, not
  `1+2Nσ`.
* The vertical arrays are constant `TensorValue` / `ThirdOrderTensorValue` objects, and every layer
  sum is a contraction. **The residual has no vertical-index loops.**
* The result is well-typed by construction, ~3.6× faster and with ~13× fewer allocations than the
  per-layer form.
* Because `Operation` is forwarded for `DistributedCellField`, **one residual serves sequential and
  MPI runs.** (`ARCHITECTURE.md` §1)

---

## 7. Standing rules

*Condensed 2026-10-09. The full text and evidence of every rule are in
[`archive/STANDING_RULES_FULL.md`](markdown_files/archive/STANDING_RULES_FULL.md); rule numbers are
unchanged, because code and documents cite them. Rules earned on v1 still apply: read `:full` as
`nl_pressure=true` and `:native` as the v1 partial `𝓝`.*

### Model and residual

1. **`R_P` is the entire frequency dispersion.** Only the boundary part of its integration by parts
   vanishes; the volume part must be assembled.

1b. **The `C⁰` objection applies to the trial space, not the test space.**
    * The model carries `∂³u`; the weak form carries `∂²u`.
    * The distributional `∂²` of a `C⁰` trial field is the cellwise Hessian plus a skeleton Dirac
      layer, so typing `∇∇(u)` drops the layer. v2 assembles it (`src/broken.jl`).
    * **Raising `fe_order` does not help:** Lagrange `H1` is `C⁰` for every `p`.
    * A test function's `∂²v` is computable cellwise.
    * The task is a *consistent* broken formulation, not a smoother space.

2. **`fe_order ≥ 2`.** Q1 zeroes `R_P` and all non-hydrostatic physics.

2b. **Always use Taylor–Hood, `p_eta = p_u − 1`; never equal order.**
    * η plays the Stokes-pressure role, so equal order is inf-sup deficient and produces a `2dx`
      checkerboard.
    * Enforced by `check_taylor_hood`.
    * A comparison across pairings is not a one-variable comparison.

2c. **Every output directory is named by `output_dir_name`**:
    `<model>_<domain>_<wave>_<regime>_<nlp>_<bed>_<discr>_<amplitude>_<period>[_extra]`.
    * Fields are never omitted.
    * Impossible combinations are refused.
    * `unique_output_dir` never overwrites.
    * `BALFEM_OUTDIR` overrides, for scratch work.

2d. **Raising the polynomial order changes two things at once:** what the discretisation can
    represent, and the effective resolution.
    * The residual is unchanged, so it is never "more physics".
    * A stable low-order run may be solving a partially masked operator.
    * To separate the two effects, match both DOFs and `dx`.

3. **`B_stored = −B̃` is negative definite.**
    * The explicit `(−1)` factors are load-bearing.
    * The elementwise sign `B ≤ 0` holds only at `p = 1`; never take `abs.(B)`.

4. **Assembly invariant: every classification row has exactly one consumer**, guarded by the
   conjunction of its activation conditions; the `𝓐/𝓚` packages are alternatives, never addends.
   A flat-bed regression can never test ∇h code.

5. **`∂R/∂u̇` is exact; nonlinear `∂R/∂u` is quasi-Newton by choice** (the gap vanishes at order
   ≈ 1.1 in amplitude).
    * An omission is benign only if it is higher order in amplitude.
    * Measure with `test_jacobians_ad`; never assume.

6. **Never apply `∇` to an `Operation`-composed expression containing a test basis.** Expand by hand.

6b. **When inserting a helper above a function, anchor above its docstring**, not on the `function`
    line. Otherwise the package fails to precompile.

7. **Nested closures must declare `local`.** An assignment otherwise overwrites the enclosing
   variable. The symptom can be flattering (a model that looks perfect). Re-run the audit after
   adding any nested helper.

### Boundaries and stability

8. **Solid-wall Dirichlet BCs must include the corner tags.**

9. **IC-release problems need `x_wall_bc=true`.**

10. **The sponge must damp η** (`+∫ μ q η`): the open-boundary spurious mode is η-dominated.

11. **Sponge strength saturates past `μ_max ≈ 5ω`; width is the lever**, sized for the longest
    component.

12. **1-D cases use `ny = 1`, `y_wall_bc=:wall`, `Ly = dx`.**
    * The wall is exact for normal incidence, and costs 2.8× fewer DOFs than `ny=3` periodic.
    * `ny ≥ 3` only for y-periodic meshes.
    * **A 1-D case is always normal-incidence.** Directional settings error in the 1-D driver.

12b. **The 2026-09 "nonlinear instability" was an equal-order (`Q2/Q2`) artefact.** Taylor–Hood
     removed it.
     * Measurements on equal order are **void**, not superseded.
     * Lesson: diff the failing configuration against the verified one, parameter by parameter,
       before diagnosing.

12c. **"Stable" means nothing without `dx`, `dt` and a duration.** Dissipation can mask growth. A
     masking objection is retired by integrators of different character agreeing, not by a `dt`
     ladder.

13. ~~`A_wave ≤ 0.001`~~ **Lifted:** that cap was rule 12b's artefact. Physical limits (Miche,
    breaking) remain.

14. **The `c_g` transit trap.** Wait `(x_sponge − x_source)/c_g + 3T` before reading a steady state.

### Solver and execution

14b. **The interior source does not deliver `A_wave`.** `η/A` = 3.32 or 2.13 depending on geometry;
     rescale every `:inner_res` result.

14c. **A crash is not evidence until its control runs.** Pair every failure with the run that
     differs in exactly one axis.

15. **`SDIRK_2_2` is Gridap's `DIRK22(1,0,1)`.**
    * It is A-stable but not L-stable, and loses ≈ `0.75(ωΔt)⁴` per step: 100× the textbook SDIRK2.
    * **No stability claim may rest on it.** Pair every stability run with Crank–Nicolson
      (`BALFEM_SOLVER=theta`).
    * Never fix a non-dissipative test by moving a threshold.
    * Read the tableau before quoting its damping.

16. **Distributed linear solve: `NewtonSolver(GMRESSolver(Pr=Jacobi))`.**

17. **`krylov_m` (memory) and `ls_maxiter` (time) are different bounds.** `restart=true` is
    load-bearing. The symptom of getting it wrong: the same `gmres=` count every step.

17b. **The incomplete `jacobian_u` does not affect stability (proven).** Hand and AD Jacobians crash
     at the same time; only the Newton count differs. Never diagnose an instability by completing
     the Jacobian, and never read a Newton stall as the cause.

18. `norm(PVector, Inf)` is broken; reduce over `own_values`.

19. Distributed ICs use `interpolate_everywhere`.

20. Keep `ConsecutiveMultiFieldStyle`: block style breaks Jacobi's `diag`.

21. Launch MPI with `~/.julia/bin/mpiexecjl`. Exit 143 and the OFI error at finalise are benign.

22. Keep cells near-isotropic: a 4:1 mesh costs ~1.6× the GMRES iterations.

23. Use MPI only when each rank gets a meaningful share; spend spare cores on more cases.

24. **`nl_tol` changes the answer** (integer Newton counts); `ls_rtol` does not.

25. A sysimage is valid only for what it was built against; rebuild after any `src/` edit.

26. Julia buffers stdout when it is redirected; flush or poll.

27. Revise does not hot-swap signature or struct changes; restart.

### Long campaigns

41. **Long-lived Julia+Gridap processes degrade** (RSS 1.5 → 3.9 GB, CPU down to 7 %). Bound the
    worker lifetime and supervise restarts.

42. **Use independent processes, not `pmap`.** It is slow, and its stdout is invisible.

43. **Never slice a cost-sorted queue contiguously.** Use LPT for makespan, SJF for coverage.

44. **Make in-flight work visible** with PID-keyed claim files in a shared directory.

45. **Switch `RETRY_ERRORS` off once a failure is established.** Failures reproduce bit-identically.

46. **Reproducing a published optimisation (vopt):**
    * the median is over the design scan;
    * `E_cg` is first power;
    * `E_shoal` is signed, with no `|·|`; its sign cancels only with the matching median;
    * κ was a fit read off its own fit and was deleted;
    * `M ≥ 3` cannot be recovered from the publication, and that is settled.

    Use `Δc` = 1e-3 / 5e-3 / 3e-2 for `M` = 2 / 3 / 4.

47. **This workstation loses run evidence.**
    * `/tmp` is wiped on boot, so redirect logs into `output/`.
    * Ending the GNOME session kills detached runs unless linger is enabled (it is).
    * Suspend is safe; reboot is not.
    * A traceback ending in `createpvd` is not a VTK failure: `grep -m1 '^ERROR'`.

### Testing and measurement

28. **"The suite passes" is not "the model is verified."** Only the analytic MMS can detect a
    self-consistently wrong residual.

29. **AD is an oracle for the Jacobian, never for the residual.**

30. **A bounds check is not a value check.**

31. **A conjunction guard needs a case that satisfies the conjunction and asserts a value.**

32. **A refinement study measures the dominant error.** Verify isolation in both directions;
    guards come in pairs.

33. **Read the pairwise rate sequence, not the fitted slope.** A rising sequence with the error
    still falling is pre-asymptotic.

34. **Resolution principle:** a test validates a term only if it can resolve that term's
    contribution.

35. **Batch verdicts come from gate output, never from exit codes.**

36. **Do not collapse per-field measurements into one pass/fail.**

37. **A reference-vs-reference comparison needs every discretisation parameter matched.**

38. **An error naming a library type is where the bug surfaced, not where it lives.** Print the
    types first.

38b. **Refinement that makes things worse is a discretisation-stability problem, never
     under-resolution.** Measure this first.

38c. **Copy the reference environment verbatim and vary one variable**, with the null control in the
     same batch.

38d. **Verify a new knob is live before spending the run on it.**

38e. **A reversal too large to be the effect under test is a bug signal.**

38f. **Identify an instability by what grows, never by what is largest.** An extremum at the window
     edge is a finding about the window.

38g. **A null result needs its effect size measured against the effect being explained.** Print
     precision is part of the experiment.

38h. **The launcher is part of the experiment.**
     * An override after the `balfem_local_run` line is a no-op.
     * Never edit a shell script while bash runs it.
     * `run_flume_1d.jl` exits 1 on success (undefined `tag`).

38i. **To separate an interior instability from a boundary feed, remove the boundary:** use the
     closed x-periodic box, under Crank–Nicolson.

39. **Test a diagnosis against a case it cannot explain.**

39b. **A construction that cannot commit the suspected error is the strongest discriminator.** When a
     structural fix fails to move a signature, retire the hypothesis and everything queued on it.

40. **Diagnostics.**
    * Never read `growth` without `x_at_max`.
    * `dmp/int ≫ 1` is trouble; `< 1` is not proof of health.
    * Mass drift is an invariant only in a closed basin.

---

## 8. Conventions for editing here

* Run Julia via the **`julia-mcp`** tool. Runs longer than 30 min go detached, logging into
  `output/`.
* If you add `using X` to `src/`, add it to `[deps]`.
* New design records go in `markdown_files/`; update `STATUS.md` when the state changes; move
  finished plans to `archive/`.
