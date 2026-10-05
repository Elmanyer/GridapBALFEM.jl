# V2_SOLVER_PLAN.md — from the v1 solver to v2: one Class-III formulation, one switch

*Written 2026-10-05; updated the same day after the v1/v2 split.*

**Where everything lives:**

| | v1 (frozen) | v2 (active) |
|---|---|---|
| solver code | GitHub `GridapBALFEM.jl`: branches `main` and `v1-solver` (`d7a4df9`), tag `v1_final_solver` | branch `v2-solver`, this plan |
| LaTeX | `latex_docs/BALFEM_models_v1/` (GitHub `BALFEM_models`, tag `v1_final_solver`) | `latex_docs/BALFEM_models_v2/`, its own Overleaf project; structure plan in [`LATEX_STRUCTURE.md`](LATEX_STRUCTURE.md) |
| run outputs | `output_v1/` (22 GB, gitignored) | `output/` (started empty) |

**The main changes, at a glance:**

1. **One Class-III formulation: broken.**
   * The projected treatment (frozen `L²` projections) is deleted.
   * The mixed treatment (auxiliary `𝖦`, `𝖥` unknowns, 5–7 fields) is deleted.
   * The broken formulation (cellwise Hessians plus a skeleton layer, 3 fields, exact Jacobian) is
     no longer opt-in. It *is* the Class-III assembly.
2. **`nl_pressure` becomes a Boolean.**
   * `:none` / `:native` / `:full` are dropped.
   * `true` assembles **all eight** `𝓝` components; `false` assembles none.
   * Passing a `Symbol` is an error.
3. **No component mask.** `c3_mask` and `BALFEM_C3_MASK` are removed, so `∇𝖲` and `∇𝖻` are always
   assembled together. Every v1 stability run had `∇𝖻` (component 4) off.
4. **Six models instead of eight.** `regime × flat_bed × nl_pressure`, with the linear regime
   forcing `false`. The two `:native` models disappear, and the two `nl_pressure=true` models get
   their first MMS verification.
5. **Removed options error loudly.** The keywords `broken`, `mixed`, `p_aux` and `c3_mask` are
   removed, as are the environment variables `BALFEM_MIXED`, `_P_AUX`, `_C3_MASK`, `_NLP_INLOOP`
   and `_BROKEN`. Setting any of them is an error, never silently ignored.
6. **Distributed path, interim.** `nl_pressure=false` only, until the broken formulation is ported
   to MPI (step 10). v1's distributed `:full` relied on the projection.
7. **Stabilisation unchanged in this refactor.** `:jumpgrad` and `:ghostvolume` stay as they are.
   Choosing the default, assembling the penalty matrix once, and porting it to general meshes and
   MPI come after (§4).
8. **Output naming and provenance.** The `<nlp>` token becomes `nlp0`/`nlp1`. Every run writes a
   `run_manifest.toml` with its commit and full configuration (step 9).
9. **v1 launchers removed from the v2 tree** (still in the tag). The docs are restructured, with
   the v1 history moved to `HISTORY_V1.md`.

---

## 0. Why v2

v1 accumulated options faster than they could be validated:
* three Class-III treatments: projected, mixed and broken;
* three nonlinear-pressure tiers: `:none`, `:native`, `:full`;
* a component mask (`c3_mask`) that silently changed what "`:full`" meant. Every stability and
  stabilisation run used `C3_MASK=gs`, i.e. component 4 (`∇𝖻`) was omitted.

Results from all these configurations sit side by side, and in practice they can no longer be told
apart. v2 removes the options that have no physical meaning or are superseded, and keeps one exact
formulation.

| v1 | v2 | reason |
|---|---|---|
| `nl_pressure ∈ {:none, :native, :full}` | `nl_pressure::Bool` | `:native` is defined by a *numerical* criterion (the `𝓝` components first-order on `C⁰`), not a physical one. No ordering in amplitude or `kd` separates components {1,2,4,5} from {3,6,7,8}. Component 3 is kept and 1/2/4/5 dropped although all five are horizontal advection of parts of `w`. Component 2 is split in half. |
| `c3_mask` (∇𝖲 arm, ∇𝖻 arm) | removed | All eight components are active or none is. A partial operator is not a model. |
| projected Class III (frozen `L²` projections) | removed | It is an approximation (one-step lag, project-then-differentiate), it is the least stable treatment (closed box: dies at 8.8 s), and its cost advantage is gone: the broken formulation is also 3-field. |
| mixed Class III (`𝖦`, `𝖥` auxiliary unknowns) | removed | Exact, but 5–7 fields, a quasi-Newton `∂R/∂u`, unusable with explicit RK, and sequential only. The broken formulation is exact with 3 fields and an exact Jacobian. |
| broken Class III, opt-in (`broken=true`) | **the** Class-III treatment, implicit when `nl_pressure=true` | Consistent on `C⁰` (LaTeX §6.4); exact Jacobian `broken_class3_jacobian`; skeleton integrals are distributable. |
| stabilisation `:jumpgrad` / `:ghostvolume` | **unchanged in this refactor** | It is orthogonal to the formulation. Choosing the default, precomputing the penalty matrix, and porting it to general meshes and MPI are follow-ups (§4). |

**Unchanged:**
* the linear models;
* the `𝓛` operator;
* advection, gravity and the mass equation;
* the native components {3,6,7,8}, which are always assembled when `nl_pressure=true`;
* the exact integration by parts of the `∇h` (`𝓐`) half of {1,2,4,5} (`nlp_gradh_contrib`);
* boundary generation, sponge, relaxation, periodic BCs;
* the vertical module and `vopt.jl`;
* postprocessing.

Everything removed remains reproducible from tag `v1_final_solver`.

---

## 1. The target, precisely

**Physics selection**

| control | values | meaning |
|---|---|---|
| `regime` | `:linear` \| `:nonlinear` | unchanged |
| `nl_pressure` | `false` \| `true` | `𝓝` off / **all eight components** on. Requires `regime=:nonlinear`. |
| `flat_bed` | `Bool` | unchanged |

When `nl_pressure=true`, the eight `𝓝` components are assembled as follows:

| block | {3,6,7,8} | {1,2,4,5} |
|---|---|---|
| bed-slope `𝓐` (prefactor `H∂h`) | direct (`nlp_direct_contrib`) | exact integration by parts onto the test (`nlp_gradh_contrib`) |
| surface-slope `𝓚` (prefactor `H∂H`) | direct | **broken**: cellwise `∇∇` + skeleton layer, **both arms** `∇𝖲` and `∇𝖻` |
| leading `𝓟` (prefactor `H²`, tested by `∇·V`) | direct | **broken**, both arms |

* **Jacobians:** `∂R/∂u̇` stays exact. `∂R/∂u` includes `broken_class3_jacobian`, the exact
  linearisation of the broken blocks.
* **Skeleton context:** `nl_pressure=true` *requires* a skeleton context on the problem. The residual
  fails loudly if it is missing, never silently.

**Models.** The MMS / AD numbering becomes 6 models instead of 8: `regime × flat_bed × nl_pressure`,
with the linear regime forcing `nl_pressure=false`.

| v2 model | regime | bed | `nl_pressure` | v1 equivalent |
|---|---|---|---|---|
| 1 | linear | flat | — | 1 |
| 2 | linear | variable | — | 2 |
| 3 | nonlinear | flat | false | 3 (`:none`) |
| 4 | nonlinear | variable | false | 4 (`:none`) |
| 5 | nonlinear | flat | **true** | 7 (`:full`), now broken and with all components |
| 6 | nonlinear | variable | **true** | 8 (`:full`), same |

v1 models 5–6 (`:native`) are dropped.

---

## 2. Step-by-step implementation

Each step ends with:
* the package precompiling;
* the gates listed for that step passing;
* a **commit** whose message names the step (`v2 step N: …`).

No step is merged with another, so any regression can be bisected to one step.

### Step 0 — Baseline from v1, before touching code

Measure the reference numbers that must survive the refactor, on the current tree (identical to
`v1_final_solver`). Write them to `output/v2_baseline/` and record them in this file (§3):
* `test_jacobians_ad` — models 1–4 (`:none`) and 7–8 (`:full`);
* the MMS orders and error magnitudes, models 1–4, Q3/Q2, one ladder each;
* `test_broken_formulation` — all 32 gate values, in particular G1, G4, G9b, G10;
* `test_linear_newton_gate`, `test_conservation`, `test_energy` values;
* one short closed-box run (Q3/Q2, Crank–Nicolson, 5 periods): broken `:full` with **both** arms
  (`C3_MASK=both`), with the ghost penalty at γ = 0.01. Its `diagnostics.csv` becomes the bitwise
  reference for step 2.

### Step 1 — The `nl_pressure` switch becomes a Boolean

* **Files:** `src/problem.jl` (`resolve_physics`, `build_problem`, `build_problem_raw`,
  `BALFEMProblem`), `src/utilities.jl` (`setup_and_run`), `src/utilities_dist.jl`,
  `src/mms_driver.jl` (`run_mms_case`, `run_model_case`), `src/mms.jl`, `src/timeloop*.jl`
  (banners).
* `resolve_physics(; regime, nl_pressure::Bool=false, flat_bed)`:
  * returns one flag, `nl_pressure`, replacing `nl_pressure68` and `nl_pressure_full`;
  * errors on `regime=:linear, nl_pressure=true`;
  * **errors on a `Symbol`** with a message pointing to this file. Old call sites must fail, not
    silently convert.
* Every `prob.nl_pressure68` / `prob.nl_pressure_full` test is replaced by `prob.nl_pressure`.
* **Default:** `nl_pressure=false`, the same as v1's `:none`. Decision point D-a (§5).
* **Gates:** package loads; models 1–4 reproduce the step-0 baseline **bitwise**, since nothing in
  their path changes.

### Step 2 — One Class-III assembly: broken, all components

* **Files:** `src/problem.jl` (`global_residual`, `jacobian_u`), `src/broken.jl`,
  `src/nlpressure.jl`.
* The Class-III branch in `global_residual` becomes a single path: if `prob.nl_pressure` is set,
  assemble `nlp_direct_contrib`, then `nlp_gradh_contrib` (unless `flat_bed`), then the broken
  `𝓚`/`𝓟` blocks. The `mixed → broken → projected` cascade is deleted.
* Remove `c3_mask` from `BALFEMProblem` and from every signature. `_c3_sum`, `_c3_layer` and
  `broken_class3_cell_fields` assemble both arms unconditionally.
* Remove the `broken` keyword: brokenness is implied by `nl_pressure=true`.
* `jacobian_u` always adds `broken_class3_jacobian` when `nl_pressure=true`. Its `∇𝖻` arm must be
  checked, because all v1 runs had it off (gate below).
* **Skeleton attachment.** `build_problem` gains a `model` keyword. When `nl_pressure=true` it builds
  the skeleton context itself (`build_skeleton_ctx`), and `global_residual` asserts the context
  exists.
  * `attach_skeleton!` remains for the stabilisation settings only: γ, order, method.
  * Decision point D-b: alternatively, `setup_and_run` attaches it.
* **Gates:**
  * G1 / G10 of `test_broken_formulation` with **both arms on**;
  * a new gate: the `∇𝖻`-arm Jacobian against finite differences on a non-flat state with
    `∇H ≠ 0`;
  * the step-0 closed-box run reproduces its `diagnostics.csv` to round-off;
  * models 5–6 (`nl_pressure=true`) pass the AD Jacobian comparison.

### Step 3 — Delete the mixed formulation

* **Files:**
  * delete `src/mixed.jl`;
  * `src/horizontal.jl`: the auxiliary spaces `Vaux` and `p_aux`, and the 5/7-field layouts;
  * `src/utilities.jl`: the `mixed`, `p_aux` and `n_aux` branches, `mixed_consistent_ic`, the
    mixed ODE operator;
  * `src/timeloop.jl`: the mixed operator and solver paths;
  * `src/monitor.jl`, `src/errors.jl`: mixed-aware diagnostics (`_n_multifields`);
  * `src/mms_driver.jl`;
  * `src/GridapBALFEM.jl`: the include and the exports (`global_residual_mixed`,
    `build_ode_operator_mixed`, `make_initial_conditions_mixed`, `mixed_n_aux`,
    `mixed_coupling_jacobian`, `mixed_delta_S`, `mixed_consistent_ic`).
* `_n_multifields`-style helpers collapse to the fixed 3-field layout.
* **Tests:**
  * delete `test_mixed_formulation.jl`, `test_mixed_jacobian.jl` and `test_diagnostics_mixed.jl`;
  * remove `test_class3_split.jl`, or reduce it to a broken-only statement of the `{1,2,5}`
    reduction (it measured the two arms through the mask);
  * update `runtests.jl`.
* **Gate:** `grep -rn "mixed\|p_aux\|n_aux" src/` returns only unrelated hits (e.g. `vopt.jl`); the
  full `:fast` suite passes.

### Step 4 — Delete the projected treatment

* **Files:**
  * `src/nlpressure.jl`: delete `build_nlp_ctx`, `update_nlp_state!`, `refresh_nlp_state!`,
    `nlp_enable_inloop!`, `nlp_plain_iterate`, `NLP_REFRESH_COUNT`, `nlp_frozen_N`,
    `nlp_gradH_frozen_contrib` and `nlp_P_frozen_contrib`;
  * keep `nlp_direct_contrib`, `nlp_gradh_contrib`, `alg_bed_hessian`, `nlp_class3_reduced_fields`,
    `_c3_sum`, `nlp_gradH_reduced_contrib` and `nlp_P_reduced_contrib`, which the broken path uses;
    rename where "reduced" no longer contrasts with anything;
  * `src/timeloop.jl`, `src/timeloop_dist.jl`: the nlp-context priming and the per-step refresh;
  * `src/utilities.jl`, `src/utilities_dist.jl`, `src/mms_driver.jl`: nlp-context construction;
  * the distributed CG + Jacobi mass solves for the projections;
  * the `NLP_INLOOP` environment variable.
* **Tests:** delete `test_nlp_inloop.jl`. Review `test_class3_residual_parity.jl`, `test_equivalence.jl`
  and `test_nlpressure.jl`: keep what tests the operator, drop what tests the projection. Update
  `test_nlpressure_distributed.jl` (see step 5).
* **Gate:** `grep -rn "nlp_ctx\|update_nlp_state\|frozen\|inloop" src/` is empty; the suite passes.

### Step 5 — The distributed path

Removing the projection leaves the distributed path with **no** Class-III implementation.
* **Interim (this step):** `setup_and_run_distributed` accepts `nl_pressure=false` only, and errors
  on `true` with a pointer to step 10. Models 1–4 stay fully distributed.
* **Tests:** `test_basic_distributed`, `test_mms_distributed_parity` (models 1–4) and the cluster
  tests; `test_nlpressure_distributed` is disabled with a reason until step 10.
* **Gate:** the distributed suite on 4 ranks (`mpiexecjl`), models 1–4, matches sequential to the
  v1 parity tolerance.

### Step 6 — MMS and the model numbering

* **Files:** `src/mms.jl` (forcing per model), `src/mms_driver.jl`,
  `test/test_mms_convergence*.jl`, `test/test_mms_forcing_nonlinear.jl`, `test_jacobians_ad.jl`,
  `examples/local_mms/*`.
* The forcing for `nl_pressure=true` is v1's `:full` forcing, which already contains all eight
  components. The `:native` forcing is removed.
* Renumber the models to 1–6 (§1).
* **Gates:**
  * models 1–4 reproduce the step-0 orders and magnitudes;
  * `test_jacobians_ad` covers all 6 models;
  * the **first MMS of models 5–6** (`nl_pressure=true`, broken, exact Jacobian), Q3/Q2,
    unstabilised and at γ\* — the first verification of the full nonlinear pressure.

### Step 7 — Drivers, launchers and the output naming

* **Environment variables:**
  * `BALFEM_NL_PRESSURE` takes `0`/`1`;
  * remove `BALFEM_MIXED`, `BALFEM_P_AUX`, `BALFEM_C3_MASK`, `BALFEM_NLP_INLOOP` and
    `BALFEM_BROKEN`. Each removed variable **errors if set** (rule 12's lesson: an ignored knob is
    worse than a refused one).
* **Drivers:** `examples/local_1d/run_periodic_1d.jl`, `run_flume_1d.jl`, `stability_eig.jl`,
  `cip_eigen_analysis.jl`, `examples/local_2d/*`, `examples/distributed*/*`, `examples/*.jl`.
* **`output_dir_name`:** the `<nlp>` token becomes `nlp0` / `nlp1`; the `broken` / `c3gs` extras
  disappear; the stabilisation token stays (`ghost0.01`, `jg2e-3o2`). Fields are still never
  omitted (rule 2c).
* **Launchers:** `run/local/` has 122 v1 launchers, almost all for configurations that no longer
  exist. Remove them from the v2 tree; they are preserved in the tag. Keep `balfem_local.sh`,
  `run/balfem_env.sh` and `run/SNELLIUS_ROME_LAUNCH_CONFIGS.md`. Write v2 launchers as campaigns
  need them. Decision point D-c.
* **Gate:** every driver runs a 2-step smoke test in each of its supported configurations.

### Step 8 — The test suite, rebaselined

* Run the full `:fast` + `:slow` suite, the `test/local/` suite and the distributed suite.
* Compare every reference constant against the step-0 baseline:
  * models 1–4 must agree to round-off;
  * any change in models 5–6 is expected and must be **re-measured, not copied**.
* Update `TEST_SUITE.md` with the v2 scores.
* **Gate:** every registered test passes; the known v1 failures are either fixed or re-listed with
  their v2 status.

### Step 9 — Run provenance (proposed)

Every run writes `run_manifest.toml` into its output directory, containing:
* the git commit, a dirty-tree flag and the branch;
* every keyword argument of `setup_and_run` / `run_mms_case`, and every `BALFEM_*` environment
  variable;
* Julia and Gridap versions, host and start time.

A run started from a dirty tree prints a warning in its banner. This is the measure that prevents a
repeat of the v1 mixing of results; it costs one small function in `src/utilities.jl`.

### Step 10 — Distributed broken formulation (later; independent of steps 0–9)

Port the skeleton context and the broken Class-III blocks to `DistributedDiscreteModel`.
GridapDistributed supports `SkeletonTriangulation` and `jump`/`.plus`/`.minus` on distributed
fields; facets on partition boundaries see both cells through the ghost layer. Then lift step 5's
refusal.
* **Gate:** a sequential-vs-distributed parity test of the residual and Jacobian, models 5–6.
* The ghost-penalty distributed port and the GMRES + Jacobi conditioning with the penalty belong to
  the stabilisation follow-ups (§4).

### Step 11 — Documentation

* `CLAUDE.md`: move the v1 history (§5.2b–f, the void and superseded findings, the `:native`
  production-tier text) into a new `markdown_files/HISTORY_V1.md`. Keep the map, the standing rules
  that remain valid, and a v2 status section.
* `MODEL.md` §5–6: the switch table and the component assembly table of §1.
* `ARCHITECTURE.md`, `VERIFIED_SCOPE.md` (6 models; scope re-stated after step 8),
  `CONFIGURATION.md`, `RUNNING.md` (new environment variables), `INDEX.md`.
* `NEW_TREATMENT.md`, `BROKEN_FORMULATION_PLAN.md`, `GHOST_PENALTY_PLAN.md`: mark them as v1
  records.
* Paths: v1 results are now under `output_v1/`. Every `output/...` reference in the markdown that
  points to a v1 result is updated.
* The LaTeX restructure of `latex_docs/BALFEM_models_v2/` (Overleaf) follows `LATEX_STRUCTURE.md`
  revision 4, edited by the author. The broken
  formulation becomes the only implemented treatment, and the projected treatment the appendix.

---

## 3. Baseline values (filled in at step 0)

*To be recorded here: model, quantity, v1 value, tolerance for v2.*

---

## 4. After the refactor — the first v2 campaign (not part of this plan's code changes)

In order; each depends on the previous one:
1. **The unstabilised closed-box ladder, all eight components** (Q3/Q2 and Q2/Q1, Crank–Nicolson):
   the v2 evidence for chapter 8. It replaces the v1 `gs`-mask results.
2. **The stabilised closed box, all eight components:** the ghost penalty at γ ∈ {0.0033, 0.01,
   0.03}, `:jumpgrad` order 2 as a control, 16 and 32 cells/λ, A = 0.10 and 0.15.
3. **MMS of models 5–6 at γ\*** (with step 6 if not already done): does the penalty preserve the
   orders?
4. **The flume with the ghost penalty**, then **bathymetry (bar)**, then the **amplitude ladder**.
5. **Stabilisation engineering:**
   * choose the default method;
   * assemble the constant penalty matrix once;
   * build the ghost-volume extension as a per-cell physical affine map (general meshes);
   * the distributed ghost penalty.

---

## 5. Decision points for the author

* **D-a.** Default of `nl_pressure`. Proposed: `false`, i.e. explicit opt-in. The full model needs
  Q3/Q2 and a stabiliser, so a silent default would mislead.
* **D-b.** Where the skeleton context is built: in `build_problem` via a `model` keyword (proposed),
  or by the caller.
* **D-c.** The v1 launchers: remove them from the v2 tree (proposed; they stay reachable via the
  tag), or move them to `run/v1/`.
* **D-d.** Whether `nl_pressure=true` without a stabiliser is allowed. Proposed: allowed with a
  banner warning, because chapter 8 needs unstabilised runs.

---

## 6. Execution record (2026-10-05)

**Decisions adopted** (the proposals of §5, applied as written):
* **D-a.** `nl_pressure` defaults to `false`.
* **D-b.** `build_problem(…; model, quad_degree)` builds the skeleton when `nl_pressure=true`. A
  full-pressure problem without it is refused.
* **D-c.** The v1 launchers were removed: 120 in `run/local/` and the 34 in `run/dist_small/`. So
  were the v1 campaign scripts in `examples/local_mms/`: the vbasis campaign/shard/report, the
  Phase-B shard and supervisors, the quadrature probe and the broken-MMS check. All remain in the
  tag. Generic drivers were ported.
* **D-d.** `nl_pressure=true` without a stabiliser runs, with a banner warning.

**Deviation from the step order.** Steps 1–4 could not each compile on their own: `mixed.jl` and the
projection code read the very struct fields steps 1–2 remove (`c3_mask`, `nl_pressure_full`,
`nlp_state`, `nlp_ctx`). They were therefore done as ONE source change, gated by the stronger
criterion: every configuration of the regression snapshot reproduces v1 entry by entry.
Steps 5–9 and 11 followed in the same pass; step 10 is open.

**What changed, by file.**
* **`src/problem.jl`.** `BALFEMProblem` now has `nl_pressure::Bool` in place of
  `nl_pressure68`/`_full`, and no `nlp_state`, `c3_mask` or `nlp_ctx`. `resolve_physics` refuses a
  `Symbol`. `build_problem` takes `model`/`quad_degree`. The residual has a single Class-III path
  (`broken_class3_residual`), and `jacobian_u` always adds `broken_class3_jacobian` when
  `nl_pressure` is set.
* **`src/broken.jl`.** One skeleton record (`skeleton_nt`). `attach_skeleton!` has no `broken`
  keyword, and both arms are always on. New: `broken_class3_residual`.
* **`src/nlpressure.jl`.** The projection machinery is deleted.
* **`src/mixed.jl`.** Deleted.
* **`src/horizontal.jl`, `src/monitor.jl`.** No auxiliary fields.
* **Time loops.** No projection hooks.
* **`src/utilities.jl`.**
  * `setup_and_run` takes the Boolean switch and loses the mixed/broken/mask/in-loop options; it
    warns on unstabilised full-pressure runs.
  * New `check_v1_env` and `write_run_manifest`.
  * The `nlp0`/`nlp1` naming token.
* **`src/utilities_dist.jl`.** Refuses `nl_pressure=true`; writes the manifest on rank 0.
* **`src/mms.jl`, `src/mms_driver.jl`.** Boolean switch; all eight components when on; no projection
  context.
* **Tests.**
  * Deleted: mixed formulation / mixed Jacobian / mixed diagnostics, in-loop projections, the mask
    split, projected parity.
  * Ported: Jacobians vs AD (six models), the nonlinear MMS (models 3–6, with 5–6 rate-gated at
    `a_eta = 0.4`), the nonlinear MMS forcing, nlpressure, selfconsistency, broken formulation,
    distributed parity, and the local and cluster tests.
  * `test_nlpressure_distributed` disabled until step 10.
* **Drivers.** `check_v1_env` on load; `nl_pressure_flag()`; the validation examples moved to
  Q3/Q2 plus the ghost penalty; `stability_eig.jl` ported to models `:off`/`:on` with a
  skeleton-aware FD sparsity pattern.
* **Docs.** `CLAUDE.md` rebuilt for v2, with the v1 version preserved as `HISTORY_V1.md`. MODEL,
  ARCHITECTURE, RUNNING, CONFIGURATION and INDEX updated. The v1-record banners added to
  VERIFIED_SCOPE, TEST_SUITE, OPEN_ISSUES, PLANNED_CAMPAIGNS, NEW_TREATMENT,
  BROKEN_FORMULATION_PLAN and GHOST_PENALTY_PLAN.
