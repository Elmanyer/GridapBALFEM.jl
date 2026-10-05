# CLAUDE.md — `GridapBALFEM.jl/` · the 2D BALFE-M algebraic wave solver (v2)

> ## ⇨ START HERE
>
> This file is the map and the standing rules. **It stays at the repository root**, because Claude
> Code auto-loads `CLAUDE.md` from the root and parent directories only. The rest of the tracked
> documentation lives in **`markdown_files/`**; start from [`INDEX.md`](markdown_files/INDEX.md).
>
> **v1 AND v2 ARE SEPARATED IN CODE, DOCUMENTATION AND OUTPUT (2026-10-05):**
>
> | | v1 (frozen) | v2 (active) |
> |---|---|---|
> | code | branches `main` and `v1-solver` (`d7a4df9`), tag `v1_final_solver` | branch **`v2-solver`** |
> | LaTeX | `latex_docs/BALFEM_models_v1/` (GitHub `BALFEM_models`, tag `v1_final_solver`) | `latex_docs/BALFEM_models_v2/` (own Overleaf project) |
> | outputs | `output_v1/` (gitignored) | `output/` |
> | this file | [`HISTORY_V1.md`](markdown_files/HISTORY_V1.md) — the v1 `CLAUDE.md`, verbatim | here |
>
> **What v2 is.** One exact Class-III formulation (the **broken** one: cellwise Hessians plus the
> skeleton layer, 3 fields, exact Jacobian), and **`nl_pressure::Bool`** — all eight `𝓝`
> components on, or none. v1's projected and mixed treatments, its `:native` tier and its component
> mask are gone; every v1 knob is refused, not ignored. Plan and status of the transition:
> [`V2_SOLVER_PLAN.md`](markdown_files/V2_SOLVER_PLAN.md).
>
> | document | what it answers |
> |---|---|
> | [`V2_SOLVER_PLAN.md`](markdown_files/V2_SOLVER_PLAN.md) | **v2: what changed and why, step by step**, the baseline values and the first v2 campaign |
> | [`LATEX_STRUCTURE.md`](markdown_files/LATEX_STRUCTURE.md) | **the v2 document's structure** (`BALFEM_models_v2`): chapter order and the reasoning |
> | [`INDEX.md`](markdown_files/INDEX.md) | which file to open; each recurring topic's one authoritative home |
> | [`MODEL.md`](markdown_files/MODEL.md) | the maths: σ-tensors, the global residual term by term, `𝓝`, Jacobians |
> | [`ARCHITECTURE.md`](markdown_files/ARCHITECTURE.md) | code structure: stacked layout, `src/` map, FE spaces, time loops, distributed path |
> | [`VERIFIED_SCOPE.md`](markdown_files/VERIFIED_SCOPE.md) | is it correct, and how far does the claim reach? |
> | [`TEST_SUITE.md`](markdown_files/TEST_SUITE.md) | what is checked, and what would slip through |
> | [`OPEN_ISSUES.md`](markdown_files/OPEN_ISSUES.md) | what is known wrong or unexplained |
> | [`CONFIGURATION.md`](markdown_files/CONFIGURATION.md) · [`RUNNING.md`](markdown_files/RUNNING.md) | settings and their evidence · how to launch (local, cluster, env vars) |
> | [`WAVE_GENERATION.md`](markdown_files/WAVE_GENERATION.md) | sources, Dirichlet generation, sponge, WaveSpec coupling |
> | [`HISTORY_V1.md`](markdown_files/HISTORY_V1.md) | **the whole v1 record**: every campaign, measurement and lesson up to the freeze |
> | v1 design records | [`NEW_TREATMENT.md`](markdown_files/NEW_TREATMENT.md) (projected/mixed), [`BROKEN_FORMULATION_PLAN.md`](markdown_files/BROKEN_FORMULATION_PLAN.md), [`GHOST_PENALTY_PLAN.md`](markdown_files/GHOST_PENALTY_PLAN.md) (stabilisation), [`PLANNED_CAMPAIGNS.md`](markdown_files/PLANNED_CAMPAIGNS.md), [`COMPLETED_VBASIS_STUDY.md`](markdown_files/COMPLETED_VBASIS_STUDY.md), [`MMS_VBASIS_CAMPAIGN.md`](markdown_files/MMS_VBASIS_CAMPAIGN.md), [`HORIZONTAL_CONVERGENCE.md`](markdown_files/HORIZONTAL_CONVERGENCE.md), [`CAMPAIGN_COST.md`](markdown_files/CAMPAIGN_COST.md), [`OUTPUT_NAMING_PROPOSAL.md`](markdown_files/OUTPUT_NAMING_PROPOSAL.md) |
>
> ⚠ **`markdown_files/` is TRACKED; `latex_docs/` is GITIGNORED** (each LaTeX folder is its own
> git repository). Anything load-bearing belongs in `markdown_files/`.

---

## 0. What this folder is

The production home of the 2D BALFE-M solver: the self-contained serial + distributed package
`src/GridapBALFEM.jl` (stacked `[η,𝖴x,𝖴y]` layout, loop-free residual), its test suite (`test/`),
sequential and cluster examples (`examples/`), the standalone postprocessing library
(`postprocessing/GridapBALFEMPost`), the vendored sea-state package (`WaveSpec.jl/`), and the design
and paper material in `markdown_files/`.

**BALFE-M** — *Basis-Agnostic Layer-integrated Finite Element*, `M` vertical elements — generalises
Yang & Liu (2024, *JFM* 999 A32) LFE-M to an **arbitrary vertical FE basis**. It is a
depth-integrated, **non-hydrostatic** free-surface wave model. The water column `σ ∈ [0,1]` is
discretised with `M` vertical elements of order `p`, so `u_h(x,σ,t) = Σ_j u_j(x,t) φ_j(σ)`. The
vertical velocity `w` and the non-hydrostatic pressure `p_nh` are eliminated **analytically** — no
pressure Poisson solve — into small, precomputable vertical σ-tensors. What remains is `Nσ+2`
coupled 2-D PDEs in `(H, u_1…u_Nσ)`, `H = d + η`, solved on a horizontal FE mesh by Gridap. The
structure is **(small dense vertical algebra) ⊗ (large sparse horizontal FE)**.

`Nσ = num_free_dofs(U_phi) = M·p + 1` = the number of vertical nodes; velocity modes are 1-based
`j = 1..Nσ`, one per node.

**Three names, kept strictly distinct in all prose** (defined in the LaTeX at
`\label{par: nomenclature}`):

* **BALFE-`M`** — the family this project derives, arbitrary vertical basis.
* **P`p`LFE-`M`** — a concrete member implemented/run/tabulated. Every instantiation here is `p=1`.
* **LFE-`M`** — Yang & Liu's published piecewise-linear models, used only when citing or comparing.

Tables comparing our numbers against theirs must not label both sides the same way.

---

## 1. Repository map

| path | what it is |
|---|---|
| `Project.toml` / `Manifest.toml` | the Julia package — `name = "GridapBALFEM"`. Loaded with **`using GridapBALFEM`, never `include()`**. This directory is both the package and the working environment. `[compat]` admits two Gridap minors on measured evidence — `CONFIGURATION.md` §1 |
| `src/` | the solver package, **18 files** (v2 deleted `mixed.jl` and the projection machinery of `nlpressure.jl`). `problem.jl` — `BALFEMProblem`, `resolve_physics`, the loop-free residual and hand Jacobians; `nlpressure.jl` — `𝓝` {3,6,7,8} and the exact-IBP ∇h half; `broken.jl` — Class III (broken), its exact Jacobian, and the skeleton stabilisers `:jumpgrad` / `:ghostvolume`; `utilities.jl` — `setup_and_run`, output naming, `check_v1_env`, `write_run_manifest`. Map: `ARCHITECTURE.md` §2 |
| `test/` | the suite (`runtests.jl` judges on gate output, never exit codes) + `test/local/` + `test/cluster/` + **`test/v2_migration/regression_snapshot.jl`** (the v1→v2 entry-by-entry gate). Inventory: `TEST_SUITE.md` |
| `examples/` | sequential examples, `distributed/` (cluster drivers + `_dist_common.jl`, which runs `check_v1_env` on load), `distributed_small/`, `validation/`, `local_1d/` (`run_flume_1d.jl`, `run_periodic_1d.jl` — the closed box, `stability_eig.jl` + `run_stability_eig.jl` — frozen eigen-analysis, `cip_eigen_analysis.jl` — penalty γ design), `local_2d/`, `local_mms/` (generic MMS study drivers). v1 campaign scripts removed (tag keeps them) |
| `run/` | `balfem_env.sh` (cluster), `local/balfem_local.sh` (workstation), the top-level SLURM launchers, and [`SNELLIUS_ROME_LAUNCH_CONFIGS.md`](run/SNELLIUS_ROME_LAUNCH_CONFIGS.md) (job sizing: memory sets the tier, 7k ranks for tier k/8). **v1's 156 campaign launchers were removed** (`run/local/README.md`) |
| `compile/` | the cluster sysimage build chain — `RUNNING.md` §5. ⚠ `Manifest.toml` is gitignored; `set_preferences.jl` `Pkg.add`s the two forks by URL |
| `output/` · `output_v1/` | ⚠ **gitignored**. `output/` holds **v2 results only**; every v2 run writes `run_manifest.toml` (commit, dirty flag, full config, `BALFEM_*` env). `output_v1/` is the archived v1 output (22 GB) |
| `postprocessing/` | `GridapBALFEMPost` — own environment, no dependency on the solver. `examples/periodic_growth.jl` turns a periodic-box run into band energies and growth rates |
| `WaveSpec.jl/` · `Gridap.jl/` | vendored sea-state synthesis (GitHub version — the release's `change_seed!` is broken) · the Gridap **fork** (`fix-transient-multifield-ad`) that makes transient-multifield AD work |
| `latex_docs/BALFEM_models_v1/` · `_v2/` | the v1 document (frozen) and the v2 document (own Overleaf project, restructured per `LATEX_STRUCTURE.md`). `_v1_old/` is an older v1 clone. Also `CFC2027_abstract/`, `doc_figures/` (`generate_doc_plots.ipynb`) |
| `GridapSWE.jl/` · `GridapEmbedded.jl/` | untracked / gitignored references, not dependencies |

---

## 2. The LaTeX project

Two documents, two repositories: **`latex_docs/BALFEM_models_v1/`** (frozen with the v1 solver,
GitHub `BALFEM_models`, tag `v1_final_solver`) and **`latex_docs/BALFEM_models_v2/`** (its own
Overleaf project, started as a copy of v1). v2 is restructured following
[`LATEX_STRUCTURE.md`](markdown_files/LATEX_STRUCTURE.md): model (1–3) → checks without a solver
(4–5) → horizontal weak form and regularity (6) → multi-field implementation (7) → stability of
the unstabilised discretisation (8) → stabilisation (9) → solver functionalities (10) → validation
(11); the projected treatment is Appendix A. **The author edits it in Overleaf.**

**Working rules (both documents):**



* ⚠ **COMPARE MTIMES BEFORE EDITING.** The author edits in Overleaf and exports the `.zip`, so the
  ZIP is newer after any round-trip. Correct procedure: unzip to scratch, `diff`, sync the folder
  **from the newer side** (back it up first), edit the folder, re-zip.
* ⚠ **After re-zipping, extract and compile FROM THE EXTRACT.** A `zip -x '*.pdf'` intended to skip
  a compiled `main.pdf` also matched `Figures/*.pdf` and silently shipped a zip missing four figures
  the document uses. Three clean compiles missed it because each ran against a *copy of the folder*,
  never against the zip — **the artefact that was verified was not the artefact being shipped.**
  The folder holds no build artefacts, so the zip needs no exclusion list at all.
* Conventions: structural labels `\label{chap:|sec:|subsec: name}`; cross-references written
  `\S\ref{…}`.
* Bibliography is **biblatex + biber**, 16 entries. ⚠ `biblatex`/`biber` are **not installed on this
  workstation** (only `bibtex`), so the production bibliography cannot be compile-tested here —
  validate via a bibtex shim in a scratch copy.
* Compilation: `pdflatex` ×3 → 0 errors, 0 undefined refs, 147 pp, with `pstricks` (not installed,
  unused) commented out **in the scratch copy only**; `multirow` required.

Table 4.1's Yang & Liu column matches ours: all three `C` entries agree to three significant figures
(our error at their `kd` is 2.04 / 1.99 / 2.00 %), and the `C_g`/`γ` offsets are **tolerance
conventions** — theirs implies ≈2.5 % on `C_g` and ≈0.09 absolute on `γ` against our 2 % and 0.02.

---


## 3. Feature summary

**Solver core.** The stacked `[η,𝖴x,𝖴y]` loop-free residual + hand Jacobians; time integrators
`RungeKutta(:SDIRK_2_2)` (default), `:theta` (Crank–Nicolson), generalised-α, explicit RK —
sequential LU+Newton, distributed GMRES+Jacobi+Newton. The full nonlinear physics: advection, the
full leading pressure `R_P`, and the nonlinear pressure `𝓝` with **all eight components** when
`nl_pressure=true`. Wavemaker / sponge / relaxation / wall / periodic BCs (x- and y-periodic);
runtime monitoring plus an independent governing-equation residual checker; `w_s`/`p_s` VTK
reconstruction; a `run_manifest.toml` per run.

**The nonlinear pressure `𝓝` (v2).** All eight components or none:
* `{3,6,7,8}` — first order everywhere (`c=3`'s second derivative is the analytic bed Hessian);
  assembled directly in the 𝓐, 𝓚 and 𝓟 blocks (`nlp_direct_contrib`);
* `{1,2,4,5}`, bed-slope 𝓐 half — exact integration by parts onto the test (`nlp_gradh_contrib`);
* `{1,2,4,5}`, surface-slope 𝓚 and leading 𝓟 halves (**Class III**) — the **broken formulation**
  (`broken_class3_residual`): the distributional gradients of `𝖲 = ∇·(H𝗎)` and `𝖻 = 𝗎·∇H`,
  i.e. cellwise `∇∇` Hessians plus one skeleton integral of their jumps; both arms (∇𝖲, ∇𝖻)
  always. Consistent on `C⁰` (LaTeX broken audit); the `{1,2,5}` reduction collapses three
  components onto one contraction (exact, 4.4e-16). **Exact Jacobian** (`broken_class3_jacobian`);
  the other `𝓝` blocks stay quasi-Newton (O(A²), rule 5). The skeleton record is built by
  `build_problem(…; model, quad_degree)`; a full-pressure problem cannot exist without it.
* ⚠ The leading-pressure skeleton term `⟦𝒫⟧` is excluded on purpose (h⁻⁴-conditioned effective mass).

**Skeleton stabilisation** (`attach_skeleton!`; any model). `stabilization = :jumpgrad` (normal-
derivative jumps, orders `cip_order ≤ 2` — Gridap's derivative limit) or `:ghostvolume` (the direct
ghost penalty: each cell's polynomial extended into its neighbour; every order `0…p`; uniform
Cartesian meshes, sequential). Coefficients `cip_gamma_u`, `cip_gamma_eta` scaled by `τ_u = d√(gd)`,
`τ_η = √(gd)` and `h^(s−2)`; linear, so the Jacobian is exact. γ is designed with the linear
eigen-analysis (`cip_eigen_analysis.jl`). v1 evidence (component 4 omitted — to be repeated):
the closed box at Q3/Q2 under Crank–Nicolson holds 100 periods with ghost γ ∈ [0.0033, 0.03] or
jump order ≤ 2 at γ = 2e-3; a first-order penalty only delays the failure; a Q1 surface cannot be
penalised (it locks). `HISTORY_V1.md` §5.2g, `GHOST_PENALTY_PLAN.md` §5.

**Dirichlet boundary wave generation + WaveSpec coupling.** Regular, multichromatic, or WaveSpec
`AiryState` stochastic sea states (seeded phases ⇒ rank-deterministic). The `:model`
discrete-eigenmode polarization prescribes an exact discrete transport, so a generated wave is a
solution of the discrete equations at the boundary and radiates cleanly.

**Verification.** The analytic MMS (`src/mms.jl`; independence from `problem.jl` enforced by a grep
gate) over the **six v2 models** (`regime × flat_bed × nl_pressure`); `test_jacobians_ad.jl` (hand
vs AD Jacobians, all six models); the linear one-Newton-iteration gate; the Yang & Liu collapse
(`test_yl_collapse.jl` — the only external oracle); `test_broken_formulation.jl` (layer identity,
both Jacobians against FD, stabiliser algebra); and the v1→v2 regression snapshot.

**Vertical grid optimisation** (`src/vopt.jl`) — the Yang & Liu **total relative-error functional**,
`E_total = Ē_c + Ē_cg + Ē_shoal + Ē_u + Ē_w`, each term integrated over `kd` against
`W = exp[(2^−kd − 2^−π) log 5]`. **Corrected 2026-09-08 to match eq. (3.10) of the paper exactly**
(rule 46); results in `output/local/vopt/`.
* `SigmaBasis` is a **self-contained analytic** piecewise-Lagrange σ-basis. Required because `E_u`
  and `E_w` put an absolute value *inside* the σ-integral, so no Gram tensor can absorb them and
  `φ_j(σ)` must be evaluated pointwise in the optimiser's inner loop — and building it here keeps
  Gridap's DOF numbering off the load path. Validated against `assemble_dispersion_tensors` through
  the permutation-**invariant** `R(μ)`: agreement **2e-15** across `p=1,2`, `M=1..4`.
* **The five terms are NOT of a common form, and two were wrong here for months.** `E_c` is in the
  **squared** celerity with `|·|`; `E_cg` is **first power** with `|·|`; `E_shoal` is
  `exp[∫(γ_e−γ_m)W/kd] − 1`, **signed, no `|·|`**; `E_u`, `E_w` carry `|·|` inside `dσ`. Verified at
  400 dpi against the PDF. See rule 46 for why the signed form is self-consistent.
* **The median is over the DESIGN SCAN, not an auxiliary sample.** Yang & Liu scan `c₂` at increment
  `0.001` and read the optimum off that same sweep; `vopt_medians`/`optimise_cbdy` now do likewise,
  one scan serving both roles. This replaced a Kronecker low-discrepancy sweep of the *width*
  simplex — a different population, hence a different argmin (rule 46).
* **`Ω` IS AN INPUT AND IS NEVER DERIVED.** `W` saturates at 0.834 rather than decaying, so the
  `kd`-integral does not converge and `Ω` is a genuine design band. Yang & Liu publish it for `p=1`
  (Table 1: `Ω_kd` = 8, 24, 80). ⚠ **`VOPT_KAPPA`, `vopt_Omega`, `optimised_cbdy` and
  `DEFAULT_CBDY_P` were DELETED 2026-09-08** — κ was fitted, then read off its own fit (rule 46).
  For `p ≥ 2` no published band exists, so **choosing `Ω` is an open question**: state it explicitly
  and pass `c_bdy`. `resolve_cbdy` returns the published set at `p=1` and a **uniform** split at
  `p ≥ 2` — visibly undesigned, rather than a mesh that looks authoritative but is not.
* **Validation — `M=2` reproduced to 0.0036** at the published band on the paper's own increment.
  `M=3`/`M=4` land 0.0189/0.0260 away and **that is settled, not open**: the optimum is converged in
  the scan increment (rule 46), and the reference specifies its population only for `M=2`.

**Linear wave properties** (added 2026-08-29, `src/utilities.jl`). `model_R`
(`R, R′, R″` from one factorisation), `airy_R`, `wave_properties` (`C`, `C_g`, `γ` for model and
Airy), `property_errors` (each already in the units of its own tolerance) and `applicable_range`.
`C_g` and `γ` did not previously exist in code — only `C`, via `dispersion_ratio`/`applicable_kd`.
Calibrated against `StokesWaveFourierAnalysis.tex` Table 4.1: **all nine published applicable ranges
reproduced to <1 %**. `assemble_dispersion_tensors` (`src/vertical.jl`) is the cheap `(Φ, Mmat, B)`
path these use — the full assembly adds `3·8·N⁴` integrals that dispersion never touches.

**Postprocessing.** VTK/CSV → analysis and plots, from-modes `w(σ)`/`p_nh(σ)` reconstruction, and a
sea-state module (Welch PSD, JONSWAP overlay, Hs, Rayleigh exceedance), validated against solver
output.

**Stability tooling.** The **closed x-periodic box** (`x_periodic=true`;
`examples/local_1d/run_periodic_1d.jl`): one model wavelength, no inflow, relaxation, sponge or
source, the discrete eigenmode plus a deterministic 1e-8 seed, sub-cell VTK, band-energy growth
rates (`postprocessing/examples/periodic_growth.jl`). The **frozen-state eigen-analysis**
(`stability_eig.jl`: coloured-FD `J*`, `M*` with a skeleton-aware sparsity pattern) and the
**penalty design analysis** (`cip_eigen_analysis.jl`).

---

## 4. Physics selection — three orthogonal controls, plus numerical options

| control | values | meaning |
|---|---|---|
| `regime` | `:linear` \| `:nonlinear` | linearised core, no advection / full nonlinear core + advection |
| `nl_pressure` | `false` \| `true` | `𝓝` off / **all eight components** (requires `:nonlinear`) |
| `flat_bed` | `Bool` | `true` ⇔ `∇h ≡ 0` (every ∇h-term dropped, ∇η/dispersion kept) |

`resolve_physics` maps these onto the internal flags (`linearised, advection, lin_pressure, P_full,
nl_pressure, flat_bed`) and **refuses a `Symbol` for `nl_pressure`** (the v1 tiers). `flat_bed` acts
at one control point (`dhx,dhy = flat_bed ? 0 : ∂h`). Full term table: `MODEL.md` §6. The six
models: 1–2 linear (flat / variable bed), 3–4 nonlinear `nl_pressure=false`, 5–6 nonlinear
`nl_pressure=true`.

**Why `𝓝` is all-or-nothing.** v1's `:native` = `{3,6,7,8}` was defined by a numerical criterion
(first order on `C⁰`): it kept component 3 and dropped 1/2/4/5 although all five are horizontal
advection of parts of `w`, and it split component 2. No ordering in amplitude or `kd` separates
them, so a partial `𝓝` is not a model. `V2_SOLVER_PLAN.md` §0.

| numerical option | values | meaning |
|---|---|---|
| `stabilization`, `cip_gamma_u/_eta`, `cip_order`, `cip_hexp` | `:jumpgrad` \| `:ghostvolume`; γ ≥ 0 | skeleton stabiliser; active iff a γ > 0 |
| `solver_type` | `:sdirk` \| `:theta` \| … | ⚠ stability claims only under `:theta` (Crank–Nicolson), rule 15 |
| `use_ad` | `Bool` | AD Jacobians (diagnosis only) |

Environment (drivers): `BALFEM_NL_PRESSURE=0|1`, `BALFEM_STAB`, `BALFEM_CIP_GU/_GE/_HEXP/_ORDER`.
`check_v1_env()` refuses `BALFEM_MIXED`, `_P_AUX`, `_C3_MASK`, `_NLP_INLOOP`, `_BROKEN` and a v1 tier
name. Output names carry `nlp0` / `nlp1` (`output_dir_name`).

---

## 5. Current Implementation Stage

*v2, 2026-10-05.* The transition is recorded step by step in `V2_SOLVER_PLAN.md`; the v1 status and
its history are in `HISTORY_V1.md` §5.

### 5.0 Status at a glance

**Done (code).** Steps 1–5, 7, 9 of the plan: the Boolean switch; one Class-III assembly (broken,
both arms, exact Jacobian); mixed and projected code deleted; distributed path restricted to
`nl_pressure=false` (refused otherwise); MMS and tests renumbered to six models; drivers, launchers
and output naming; `check_v1_env`; run manifests. See §5.1 for the gates.

**Inherited from v1 and still to be re-confirmed on v2** (the refactor is designed not to touch
them; the regression snapshot checks it entry by entry):
* linear models 1–2: optimal order in both fields, 1-D and 2-D, sequential and distributed;
* nonlinear models 3–4 (`nl_pressure=false`): theoretical order;
* the Yang & Liu collapse at `p=1` to round-off; the vertical-basis independence of the orders.

**Open, in order** (`V2_SOLVER_PLAN.md` §4):
1. the **first MMS of models 5–6** (`nl_pressure=true`, broken, exact) — never verified in v1;
2. the **unstabilised closed-box ladder with all eight components** (chapter 8's v2 evidence;
   v1's runs omitted component 4);
3. the **stabilised box**, ghost γ ∈ {0.0033, 0.01, 0.03} and `:jumpgrad` order 2, 16 and 32 cells/λ,
   A = 0.10 and 0.15; the MMS orders at γ\*;
4. the flume, the bar bathymetry, the amplitude ladder;
5. step 10 — the broken formulation (and stabiliser) distributed; step 11 remainder — the docs
   listed in the plan;
6. inherited open items: the nonlinear `p_η` order reduction at Q3/Q2 (`OPEN_ISSUES.md` §0b), the
   Q2/Q1 velocity shortfall (`HISTORY_V1.md` §5.3), the cluster production suite.

### 5.1 Migration gates

*Filled in as they run; values in `V2_SOLVER_PLAN.md` §3.*

---

## 6. The design decision everything rests on

**The vertical (layer) index lives in the FE value type, not in a Julia array.**

* Velocity is two vector-valued fields `𝖴x, 𝖴y ∈ VectorValue{Nσ}`. The MultiField is
  `[η, 𝖴x, 𝖴y]` — **3 fields, not `1+2Nσ`**.
* The static vertical arrays become constant `TensorValue` / `ThirdOrderTensorValue` objects.
* Every layer sum `Σ_j`, `Σ_{kj}` is therefore a matvec or tensor double-contraction — **the
  residual contains no vertical-index loops**.
* The MultiField is touched only for `η=U[1]`, `𝖴x=U[2]`, `𝖴y=U[3]`.

This makes the residual well-typed by construction, ~3.6× faster with ~13× fewer allocations than
the fused per-layer form, and — because `Operation` is forwarded for `DistributedCellField` —
**one residual and one pair of Jacobians serve both sequential and MPI execution**. Details,
including the verified Gridap contraction facts and the variable dictionary: `ARCHITECTURE.md` §1.

---

## 7. Standing rules

> **Reading them in v2.** These rules were earned on v1 and are carried over unchanged, because they
> are lessons about the model, the method and the workflow, not about v1's options. Where a rule
> cites `:full` read "`nl_pressure=true`"; where it cites `:native`, the projected or the mixed
> treatment, or `c3_mask`, it describes a v1 configuration (tag `v1_final_solver`,
> `HISTORY_V1.md`). Paths `output/…` in a rule's evidence are now under `output_v1/`.

These are the rules that cost something to learn. Each is stated where it is enforced; the
supporting measurement is in the linked document.

### Model and residual

1. **`R_P` is the entire frequency dispersion of the model.** Only the *boundary* part of its
   integration by parts vanishes; the volume part must be assembled. Without it the model degenerates
   to non-dispersive shallow water.
1b. **THE `C⁰` DISTRIBUTIONAL OBJECTION APPLIES TO THE *TRIAL* SPACE, NOT TO THE TEST SPACE.**
    ⚠ THE MODEL CARRIES `∂³u`; THE WEAK FORM CARRIES `∂²u`. Measured 2026-09-23 by a
    `ε(x−x₀)^k` probe: `p_nh` depends on `∂²u` and not `∂³u`; the strong-form momentum residual —
    ours and Yang & Liu's ⌊2.30⌋ alike — depends on `∂³u`. `R_P`'s integration by parts moves one
    derivative onto the test function, which is what leaves the trial field needing only a Hessian.
    Do not "correct" either count into the other; they describe different objects.
    The distributional `∂²` of a `C⁰` **trial** field is `{∂²u}` cell-wise plus a **Dirac layer on the
    skeleton** weighted by the jump `[∂ₙu]`; cell quadrature sees only the first half, so typing
    `∇∇(u)` into an integrand silently drops the layer. v1 carried the Class-III components
    `{1,2,4,5}` by frozen `L²` projections (`OPEN_ISSUES.md` §0c) or by mixed unknowns; **v2 assembles
    the layer explicitly — the broken formulation, `src/broken.jl`.**
    ⚠ **RAISING `fe_order` DOES NOT FIX THIS.** `ReferenceFE(lagrangian, …, p)` with
    `conformity=:H1` is exactly `C⁰` for **every** `p` — nodal DOFs match values across a face, never
    normal derivatives. `Q2`, `Q3`, `Q4` are all `C⁰`; the "higher order ⇒ smoother" intuition comes
    from splines (degree `p` ⇒ `C^{p−1}`), which is IGA, not Lagrangian FEM. `p ≥ 2` is *necessary*
    (a `Q1` cell-wise Hessian keeps only the mixed `∂ₓ∂ᵧ` entry and vanishes **entirely in 1-D**, so
    the term would be absent rather than approximated — the same shape as rule 2) but nowhere near
    sufficient.
    But a **test** function is a chosen, known object: `∂²v` is a computable cellwise polynomial,
    non-zero from `Q2` up and discontinuous only across faces, which is exactly what `C⁰` interior
    penalty exists to handle. ⚠ `GlobalResidual.tex` §`sec: pressure operator
    implementation` currently rules out integrating the leading-pressure block by parts "for exactly
    the same distributional reason as for the trial fields" — **that equivalence is too strong and is
    the one thing in that section that does not hold.** ⚠ Nor was this ever a missing-API problem:
    Gridap exposes `∇∇` (already used for the analytic bed Hessian, `src/nlpressure.jl:35`) and the
    full skeleton machinery (`SkeletonCellFieldPair`, jump/mean). **The task is to choose a
    formulation that is CONSISTENT with a broken `∂²` — volume term plus skeleton jump terms plus a
    penalty — not to find a library or a polynomial order that makes `∂²` classical.**
2. **`fe_order ≥ 2`.** `Q1` elements zero `R_P` and disable all non-hydrostatic physics.

2b. **ALWAYS USE TAYLOR-HOOD HORIZONTAL PAIRS: `p_eta = fe_order − 1`. NEVER EQUAL ORDER.**
    `η` enters the momentum equation **undifferentiated**, through `∇·v` after the integration by
    parts, so it plays exactly the role pressure plays in a Stokes system. Equal-order continuous
    spaces are therefore **inf-sup (LBB) deficient**, and their classic failure mode is a spurious
    **checkerboard mode at `λ ≈ 2·dx`** that refinement does not remove.

    **The evidence was already in hand and was not acted on.** The analytic MMS measures order `p`
    rather than `p+1` on equal-order spaces, and **the entire verified scope of this solver —
    all 30/30 spatial studies, five vertical bases, eight models — was measured on `Q3/Q2`.**

    ⚠ **Every nonlinear-instability run before 2026-09-05 used `Q2/Q2`** (`BALFEM_P_ETA` unset, and
    the sentinel then meant equal order). So the configuration that blew up and the configuration
    that was verified **were never the same discretisation**, and the measured instability sits at
    `λ ≈ 2–3·dx` — exactly where the equal-order spurious mode lives. That is not proof of cause,
    but it means no equal-order result can be quoted as a property of the model.

    **Enforced, not merely documented (2026-09-05):** the `p_eta = 0` sentinel in `setup_and_run`,
    `setup_and_run_distributed` and `run_flume_1d.jl` now resolves to `p_u − 1`. ⚠ **The sentinel
    was subsequently REMOVED entirely** (`0cf8cb5`): `p_eta` defaults to `p_u − 1` and
    `check_taylor_hood` REJECTS `p_eta < 1`, so equal order can only be requested explicitly and is
    then refused. `p_horizontal` was renamed `p_u` in the same commit.
    ⚠ Suite reference constants measured on equal order may shift; that is intended.

    **Corollary: a comparison across the pairing is not a one-variable comparison.** Never set an
    equal-order result beside a Taylor-Hood one and attribute the difference to physics.
3. **`B_stored = −B̃ ≤ 0`**, and the explicit `(−1)` factors in the `R_P` and slope-pressure terms
   are load-bearing. ⚠ **The invariant is NEGATIVE DEFINITENESS, not the elementwise sign.**
   `B̃ = ∫φᵢ_int φⱼ_int` is a Gram matrix, so `B_stored = −B̃ ≺ 0` for every `(M,p)` — but the
   *elementwise* `B ≤ 0` holds only because `φ_int ≥ 0` for `p = 1`. At **`p = 2` some entries are
   positive** (measured: P2LFE-2 and P2LFE-3), and that is correct, not a defect. No code depends on
   the elementwise sign — `dispersion_ratio`, `applicable_kd` and `model_celerity` all form
   `M − kd²·B_stored`, which is sign-correct for any `p`. Do not "fix" a p≥2 basis by taking `abs.(B)`.
4. **THE ASSEMBLY INVARIANT: every classification row must have exactly ONE consumer, guarded by the
   CONJUNCTION of its three activation conditions.** A guard testing only the bed condition or only
   the amplitude condition is a defect whenever the physics has a second representation in the other
   regime — which, for the `𝓐/𝓚` packages, it always does. They are **alternatives, never addends**.
   **Corollary: a flat-bed regression can never test `∇h` code.** (`MODEL.md` §7)
5. **`∂R/∂u̇` is EXACT in all eight models; nonlinear `∂R/∂u` is QUASI-NEWTON by choice**, its gap
   vanishing at order 1.11–1.16 in amplitude. **An omission is benign only if it is HIGHER ORDER IN
   AMPLITUDE** — a block whose prefactor does not scale with the solution is an `O(1)` error in the
   effective mass matrix, and Newton then converges to the fixed point of the *wrong map*. **Never
   assume; measure** with `test_jacobians_ad.jl`'s amplitude-scaling gate. Do not "complete" the
   omissions without re-measuring every nonlinear reference value.
6. **Never apply `∇` to an `Operation`-composed expression containing a test basis** — expand by
   hand via `∂_a(W⋅𝓣) = (∂_aW)⋅𝓣`.
6b. **WHEN INSERTING A HELPER ABOVE AN EXISTING FUNCTION, ANCHOR ON ITS DOCSTRING, NOT ON ITS
   `function` LINE.** Julia refuses a triple-quoted docstring immediately followed by another
   string literal (`cannot document the following expression`), so a new docstring+function
   inserted between an existing docstring and its `function` leaves the OUTER docstring
   documenting a *string* — and the whole package fails to precompile. ⚠ **This cost three
   separate edit-compile cycles in ONE session** (`NLP_REFRESH_COUNT`, `_n_multifields`,
   `mixed_aux_jacobian`): the natural anchor for a textual insert is the `function` line, and
   the docstring above it is invisible to that pattern. Anchor above the opening quotes, or
   append at end of file.
7. **NESTED CLOSURES MUST DECLARE `local`.** In Julia, a nested function assigning a name already
   local to an enclosing function **assigns the enclosing variable**. This codebase is full of long
   functions with nested helpers, so the hazard is **structural**. One instance cost two days.
   **Re-run the mechanical audit after adding any nested helper.** (`VERIFIED_SCOPE.md` §7)
   *Third instance, 2026-08-21, `wave_properties`:* two closures each called their phase speed `C`,
   which was also the enclosing function's return value, so the last γ evaluation — the **Airy** one
   — overwrote it. **Note the shape of the symptom: `|C/Ce − 1|` collapsed to ~1e-12 for every mesh
   at every `kd`, i.e. the model looked PERFECT exactly where it is worst, and `C_g` and `γ`
   reproduced their published values throughout.** A capture bug can be invisible in every channel
   but one, and the channel it corrupts can fail in the flattering direction.

2c. **EVERY RUN'S OUTPUT DIRECTORY IS NAMED BY `output_dir_name`, AND THE NAME CARRIES THE
    DISCRETISATION.** Grammar (`markdown_files/OUTPUT_NAMING_PROPOSAL.md`):

    ```
    <model>_<domain>_<wave>_<regime>_<nlp>_<bed>_<discr>_<amplitude>_<period>[_<extra>…]
    P1LFE-2_1d_bcplane_nl_full_flat_Q2Q1_A0.1_T1.6
    ```

    One generator in `src/utilities.jl`, used by all 14 drivers; the old per-driver prefixes
    (`flume_`, `small2d_`, `small_`) named the *script*, not the case. **Fields are never omitted and
    an absent field never means a default** — that convention is what let equal-order runs sit beside
    Taylor-Hood ones for months without anyone seeing it (rule 12b), and the `<discr>` token exists
    for exactly that reason.

    * `output_dir_name` **refuses impossible combinations** rather than labelling them: directional
      content on a 1-D domain (rule 12), `:linear` with `nl_pressure ≠ :none`, and any
      non-Taylor-Hood pairing (it calls `check_taylor_hood`).
    * `unique_output_dir` suffixes `_v2`, `_v3`, … and **never overwrites**. Not cosmetic: a
      re-executed batch wrote six finished runs over their own output on 2026-09-06 and destroyed
      them.
    * `BALFEM_OUTDIR` still overrides, for scratch work. The standard is the default, not a cage.

2d. **RAISING THE POLYNOMIAL ORDER IS NOT A ONE-VARIABLE CHANGE, AND IT IS NOT "MORE PHYSICS".**
    Going `Q2/Q1 → Q3/Q2` moves **two** things at once: (i) what the discretisation can *represent*
    — with `η ∈ Q2`, `∇η` is piecewise linear instead of piecewise constant, so operator content
    that was structurally absent becomes live (at `Q1`, `∂²η ≡ 0` **identically** in 1-D); and
    (ii) the **effective resolution**, since `Q3` on the same cells is ~2.1× the DOFs, and by rule
    38b more resolution *advances* onset in a grid-scale failure. ⚠ **The residual is unchanged** —
    the same `𝓝` components, the same terms — so never describe an order change as "activating more
    nonlinear effects". **Corollary, and it is the one that costs:** a run that looks *stable* at low
    order may simply be solving a **partially masked operator**, so a low-order pass is weaker
    evidence than it appears. To separate (i) from (ii), match DOF count *and* `dx` — matching only
    DOFs leaves the grid scale free (measured 2026-09-23: at ~30k DOFs, Q3/Q2 at `dx`=0.25 died
    2.3× earlier than Q2/Q1 at `dx`=0.125, which bounds but does not isolate the order effect).
    This is rule 2b's corollary in the other direction, and §5.2e is the worked instance.

### Boundaries and stability

8. **Solid-wall Dirichlet BCs must include the corner tags** — otherwise 4 corner DOFs are
   unconstrained and the run diverges exponentially.
9. **IC-release problems need `x_wall_bc=true`** — a free x-wall plus the dispersion term is a
   spurious-forcing mode that an initial perturbation excites directly.
10. **The open-boundary spurious mode is η-dominated, so the sponge MUST damp η** (`+∫ μ q η`), not
    just velocity. No value of `mu_max` absorbs it otherwise. (`WAVE_GENERATION.md` §3)
11. **Sponge strength saturates past `μ_max ≈ 5ω` — WIDTH is the lever.** And the width must cover
    the **longest** component (`kd_min` ⇒ `λ_max`), not the peak.
12. **`ny ≥ 3` is mandatory for a y-PERIODIC mesh** (Gridap `CartesianGrids.jl:39`) — and for a
    periodic mesh only. **`:wall` and `:open` accept `ny = 1`.**

    **THE DEFAULT FOR ANY 1-D HORIZONTAL CASE IS `ny = 1` WITH `y_wall_bc=:wall`**, not three cells
    with periodicity. The solver is structurally 2-D, so a 1-D problem is a narrow flume; one cell
    across with solid walls is both the cheapest and the more correct way to pose it:
    * for a normal-incidence wave the exact solution has `𝖴y ≡ 0`, and the wall condition `𝖴y = 0`
      is **exactly consistent** with it — it approximates nothing. `:periodic` merely *permits*
      `𝖴y ≡ 0` while admitting a family of y-periodic modes a true 1-D model does not have, whose
      shortest member has wavelength `Ly` and can therefore sit inside the physical band;
    * measured on the 240-cell flume: **7215 free DOFs at `ny=1`/`:wall` against 20202 at
      `ny=3`/`:periodic`** — 2.8×, and a direct LU costs more than linearly in DOFs. The wall pins
      the bottom and top `𝖴y` node layers (2886 constrained = 2 levels × 481 x-nodes × `Nσ`);
    * set `Ly = dx` so the single cell stays isotropic.

    Use `:periodic` only where the case genuinely carries oblique or short-crested content, which a
    solid wall would reflect. `examples/local_1d/run_flume_1d.jl` defaults to `ny=1`, `Ly=0.25`,
    `BALFEM_YBC=wall`, and all six `run/local/run_1d_*.sh` use it.

    ⚠ **AND A 1-D CASE IS ALWAYS NORMAL-INCIDENCE — never oblique, never short-crested.** A 1-D
    domain carries one propagation direction; oblique content has a transverse wavenumber
    `k_y = k sin θ`, and a flume one cell across cannot represent it. The request is not refused by
    the mathematics — it is **silently aliased** onto a normal-incidence wave at the wrong
    wavenumber, with the transverse component dropped, producing a wrong answer that runs to
    completion and looks plausible. This is why `build_airy_state` is called **without**
    `directional=true` in the 1-D driver, and why that driver now **errors** if any of
    `BALFEM_WAVE_DIR` (non-zero), `BALFEM_NTHETA`, `BALFEM_SPREAD_STD`, `BALFEM_THETA_MAX` or
    `BALFEM_DIRECTIONAL` is set — those variables were previously *ignored*, which is the worse
    failure. Directional content belongs in the 2-D driver
    (`run_directional_sea_small.jl`: `y_wall_bc=:open` plus lateral sponges).

    **Corollary: the `:wall` default above is not a compromise.** It is exact precisely *because*
    1-D cases are normal-incidence — the two rules support each other.
13. ~~**`A_wave ≤ 0.001 m`** for stable long fully-nonlinear integrations.~~ ✅ **LIFTED 2026-09-06 —
    THIS WAS THE INSTABILITY WEARING A DISGUISE.** The cap was never a property of the model; it was
    the equal-order pairing's growth rate being slow enough at tiny amplitude to finish a run
    (rule 12b). On **Taylor-Hood** the 1-D flume completes 50 wave periods at `A_wave = 0.10 m`
    (`κa ≈ 0.16`, `kd = 5.5`) with the mean `η` flat to the fourth decimal — **100× the old cap**.
    ⚠ This does not license unlimited amplitude: physical limits (Miche, breaking) still apply, and
    the model has no breaking closure. It removes a *numerical* restriction, not a physical one.

12b. ✅ **THE "NONLINEAR INSTABILITY" WAS AN EQUAL-ORDER ARTEFACT — RESOLVED 2026-09-06.**
    *(This is the complete record; `NONLINEAR_INSTABILITY.md` was folded in here and deleted.)*

    **The phenomenon, as it appeared for a year.** Fully nonlinear runs grew an unbounded
    free-surface mode while the *same case* linear was flat for 80 s at 4× the amplitude. Growth
    appeared in the interior, not at a boundary. **Refining `dx` made it worse** — blow-up at
    t≈24 s for `dx=0.25`, t≈10 s for `dx=0.125`, barely at all at `dx=0.50` — which is why it read
    as a discretisation-stability problem rather than under-resolution (rule 38b). The measured
    growth spectrum peaked at **`λ ≈ 2–3·dx`**, gaining 10³–10⁴× while the carrier gained ~6×.

    **The cause: equal-order `Q2/Q2` horizontal spaces.** `η` reaches the test function only through
    `∇·v`, so it plays the pressure role of a Stokes system and the pairing is subject to the
    inf-sup (LBB) condition. Equal order is deficient and admits a spurious checkerboard at
    `λ ≈ 2·dx` — exactly the measured peak. On **Taylor-Hood** the mode does not exist:

    | `dx` | equal order `Q2/Q2` | Taylor-Hood `Q2/Q1` (peak / settled t≥30) |
    |---|---|---|
    | 0.50 | suppressed | ✅ 80 s — 0.1169 / 0.11055 |
    | 0.25 | **died t=26.4, η→5.06** | ✅ 80 s — 0.1115 / 0.10593 |
    | 0.125 | **died t≈10** | ✅ 80 s — 0.1112 / **0.10481** |

    **The refinement signature INVERTED** — settled amplitude now *decreases* monotonically with
    `dx`, converging on the delivered 0.102. Newton holds at ~4 iterations/step against 30–57 as the
    equal-order runs came apart, and `x_at_max` migrates with the crest instead of pinning at one
    station. Confirmed on `Q3/Q2` as well, and on **three integrators** — `SDIRK_2_2` (dissipative, rule 15),
    `SDIRK_3_3`, and explicit `EXRK_RungeKutta_4_4` (near dissipation-free, which is what retires the
    `dt`-masking objection; RK4 tracked SDIRK to within 1–5 % with the same envelope, the offset
    being exactly the SDIRK_2_2 damping).

    **Ten hypotheses were refuted before the pairing was questioned** — worth keeping, because each
    is a real property of the solver:

    | hypothesis | verdict, on evidence |
    |---|---|
    | sponge reflection; domain length | onset identical at `Lx` = 40 and 90 |
    | boundary generation artefact | the interior source fails too, at matched *delivered* amplitude |
    | inflow relaxation zone | widening it delays onset, never prevents it |
    | CFL / time step | 4× `dt` invisible at fixed `dx` |
    | incomplete hand Jacobian | exact AD agrees to 4 dp *through the divergence* (rule 17b) |
    | Benjamin–Feir | a physical rate cannot depend on `dx`; measured σ ∝ A⁴, BF is A² |
    | under-resolution | refinement makes it worse |
    | quadrature aliasing | degrees 6/10/14 agree to 4 dp against an exact control |
    | Galerkin advection lacks an energy sink | the operator's production is *exactly* the continuity defect, and that defect is 5 orders too small — the correction moved the solution <5e-5 while the mode grew 10³–10⁴× |

    ⚠ **THE ROOT CAUSE OF THE ROOT CAUSE: the configuration that blew up was never the configuration
    that was verified.** The MMS campaign ran `Q3/Q2`; every instability run ran `Q2/Q2`, because the
    launchers left `p_eta` unset and the sentinel then meant equal order. The two were compared for
    months as though they were the same solver. **Before diagnosing a numerical failure, diff the
    failing configuration against the verified one — parameter by parameter, including the ones
    nobody thought to set.** This is now enforced by `check_taylor_hood` (rule 2b), so it cannot
    recur.

    ⚠ **Everything measured about the mode is void, not superseded** — the `σ ∝ A⁴` rates, the
    "no amplitude threshold" claim, the tier ordering, the `dt`-masking numbers, the `A_wave ≤ 0.001`
    cap (old rule 13). They characterise a discretisation no longer in use. **Discard them; do not
    re-explain them.**

    **Method lessons that cost time and generalise:**
    * **Copy the reference environment verbatim and vary one variable explicitly.** `env -i` plus a
      partial list once flipped generation from `:bc` to the interior source (2.8× the amplitude) and
      produced a confident result pointing the wrong way.
    * **Put the null-treatment control in the same batch**, never against a remembered baseline.
    * **Verify a knob is LIVE — and then that it is BIG ENOUGH TO MATTER.** A dead parameter gives
      three identical curves; a live but negligible one gives a clean negative that means nothing.
      Measuring the treatment against the effect size costs seconds.
    * **A reversal too large to be the effect under test is a bug signal, not a finding.**
    * **A "threshold" may just be a run that ended too early.** Stability claims need a duration.
    * **Separate the phenomenon from its symptoms.** Newton stalling at `‖r‖≈1.2` was real and
      reproducible but was the quasi-Newton Jacobian failing to track an already-diverging solution.
    * **Identify an instability by what GROWS, never by what is LARGEST** — ranking by amplitude
      found the carrier's harmonics and read as evidence *against* a grid mode; ranking by gain found
      it. An extremum on the edge of the search window is a finding about the window.
      (`postprocessing/examples/growth_spectrum.jl` implements both guards.)

12c. **"STABLE" IS MEANINGLESS WITHOUT `dx`, `dt` AND DURATION.** A dissipative integrator's numerical
    dissipation grows with `dt` and can *mask* a growing mode, so **a run that completes may simply
    be one whose dissipation exceeded the growth rate.** This is rule 15's trap in a new guise, and
    the rule stands on its own merits — but note how it was DISCHARGED for 12b, because that is the
    template: the same case was run on `SDIRK_2_2` (2nd order, dissipative), `SDIRK_3_3` (3rd order) and
    **explicit `EXRK_RungeKutta_4_4`** (essentially non-dissipative). All three stayed flat, and RK4
    tracked SDIRK_2_2 to within 1–5 % with the same envelope shape — the offset being exactly the
    SDIRK_2_2 damping. **Three integrators of different order and stability character agreeing is what
    retires a dissipation-masking objection; a `dt` ladder on one scheme is not.**
    ⚠ Gridap DOES provide explicit tableaux (`Gridap.jl/src/ODEs/ODESolvers/TableausEX.jl`, 18 of
    them incl. classical RK4). They cost a mass solve per stage here, since `∂R/∂u̇` carries the `B`
    dispersion term and is not the identity — but they run, and they are the right tool for exactly
    this question.
14. **The `c_g` transit trap.** At `kd = 5.5`, `c_g = 1.25 m/s` — filling a 45 m flume takes 22.5
    periods. Budget `t_settle ≈ (x_sponge − x_source)/c_g + 3T` before reading any steady state. It
    caught three separate measurements.

### Solver and execution

14b. **THE INTERIOR SOURCE DOES NOT DELIVER `A_wave`, AND THE FACTOR IS GEOMETRY-DEPENDENT.**
    Measured `η/A`: **3.32 on the 50 m small domain, 2.13 on the 60 m 1-D flume**; Dirichlet
    generation delivers ~1.0–1.05. The ratio is a *linear* property of the source calibration —
    identical to three digits between a `linear/:none` run at `A=0.1` and a `nonlinear/:full` run at
    `A=0.001` — so it is not a nonlinear artefact and cannot be tuned away. **Every `:inner_res`
    result must be rescaled before it is read**, and "3.3×" must not be quoted as a constant.

14c. **A CRASH IS NOT EVIDENCE UNTIL ITS CONTROL RUNS.** The small-domain suite is a factorial for
    a reason: a `nonlinear/:full` run diverging at `A=0.1` looks amplitude-driven, but the **linear**
    run at the same *delivered* 0.332 m — 76 % of the Miche limit — completes. Amplitude alone is
    survivable; it takes amplitude *and* nonlinearity. Conversely the directional case died at
    **5.7 %** of Miche, lower than three cases that completed, so it is a different failure
    (generation-region, velocity-led) wearing the same symptom. Pair every failure with the run
    that differs in exactly one axis.

15. **THE DEFAULT INTEGRATOR `SDIRK_2_2` IS STRONGLY DISSIPATIVE — AND IT IS NOT THE SCHEME ITS NAME
    SUGGESTS (corrected 2026-09-26).** Gridap's `:SDIRK_2_2` is `DIRK22(1,0,1)`
    (`Gridap.jl/src/ODEs/ODESolvers/TableausDIM.jl`): `A = [1 0; −1 1]`, `b = (½, ½)`,
    `R(z) = (1 − z − z²/2)/(1 − z)²`. It is A-stable but **not** L-stable (`R(∞) = −½`), and it removes
    `(¾)y⁴/(1+y²)² ≈ 0.75 y⁴` of a neutral mode's energy per step, `y = ωΔt` — **100× the textbook
    L-stable `γ = 1 − 1/√2` SDIRK2** (0.0074 y⁴), which every earlier note here had assumed. At
    `Δt` = 0.04: carrier 0.011 s⁻¹ (exactly the carrier decay of every periodic run), `ω∞` 0.15 s⁻¹,
    fastest `:full` modes 2–10 s⁻¹.
    * Any test measuring a non-dissipative property must pin `solver_type=:theta`. **Do not remove those
      pins, and never fix such a failure by moving a threshold.** Recognise it by: amplitude damped
      while **phase is correct**, error **growing with frequency**, and **refining the mesh does not
      help**. (`TEST_SUITE.md` §4)
    * ⛔ **NO STABILITY CLAIM MAY REST ON `SDIRK_2_2` ALONE.** It masked the `:full` interior instability:
      Q2/Q1 at 32 cells/λ ran 100 periods bounded under SDIRK and **diverged at ≈ 90 s under
      Crank–Nicolson**; the 64-cell SDIRK run showed no growth at all and **diverged at 31 s under
      Crank–Nicolson**; and SDIRK hid a `:native` instability altogether (Q2/Q1 64 cells/λ diverges at
      146 s under CN). Pair every stability run with a
      Crank–Nicolson repeat (`BALFEM_SOLVER=theta`), which has `|R(iy)| ≡ 1`.
    * Before quoting any tableau's damping, read the tableau in `TableausDIM.jl` — never infer it from
      the name.
16. **Distributed linear solve = `NewtonSolver(GMRESSolver(Pr=Jacobi))`.** A direct LU does not scale
    to partitioned matrices at cluster size.
17. **`krylov_m` (basis size, memory) and `ls_maxiter` (iteration budget, time) are different
    bounds**, and `restart=true` is load-bearing. Symptom of getting this wrong: **`gmres=` pinned at
    exactly the same number every step.** (`ARCHITECTURE.md` §5)
17b. **⛔ THE INCOMPLETE `jacobian_u` DOES NOT AFFECT STABILITY. PROVEN, NOT ARGUED (2026-09-16).**
    A converged Newton step is a property of the RESIDUAL, not the Jacobian: the Jacobian sets the
    *path* to the root and the *cost*, the root is defined by `r(u)=0`. So an incomplete but
    convergent Jacobian **cannot change the answer** — only the iteration count.

    **THE DECISIVE MEASUREMENT.** Every cell of the `:full` flat `dx`×`dt` factorial (§5.2c) was run
    TWICE — once with the hand Jacobian, which omits the `{1,2,4,5}` blocks, and once with
    `BALFEM_USE_AD=1`, the exact AD Jacobian of the same residual:

    | `nx` | `dt` | hand onset | AD onset | Newton, last 5 steps |
    |---|---|---|---|---|
    | 240 | 0.04 | **12.60 s** | **12.60 s** | hand 10,12,15,19,51 · AD 4,4,4,4,6 |
    | 480 | 0.04 | **5.00 s**  | **5.00 s**  | hand 6,6,8,12,58 · AD 4,4,4,6,10 |
    | 480 | 0.02 | **8.20 s**  | **8.40 s**  | hand 6,8,8,10,16 · AD 4,4,4,6,12 |

    **A CRASHING SIMULATION CRASHES IN THE SAME PLACE WITH EITHER JACOBIAN.** Onset is identical in
    two cells and one output interval (0.2 s) apart in the third; `u_max` at failure agrees to three
    significant figures (2.5046 vs 2.4971; 1.6764 vs 1.6702); `η` agrees to five digits and `r0` to
    three at every sample up to failure. **The only difference is the iteration count** — the hand
    Jacobian degrades 6 → 12 → 19 → cap while AD holds ~4 to the last step, so the quasi-Newton gap
    costs TIME and nothing else.

    ⚠ **CONSEQUENCES, both directions.**
    * **Never diagnose an instability, a divergence or a wrong answer by completing the Jacobian.**
      Completing `jacobian_u` would buy iteration count, not stability. §5.7 item 5 previously said
      the opposite and was corrected on this evidence.
    * **Never read Newton stalling as the primary fault.** A stall at 50 iterations is the
      quasi-Newton Jacobian failing to track a solution that is *already* diverging underneath it —
      the symptom, not the cause.

18. **`norm(PVector, Inf)` is broken** (PartitionedArrays 0.3.5) — reduce over `own_values`.
19. **Distributed ICs use `interpolate_everywhere`**, never `FEFunction(U, zeros(…))`.
20. **Keep `ConsecutiveMultiFieldStyle`** — `BlockMultiFieldStyle` breaks Jacobi's `diag`.
21. **Launch MPI with `~/.julia/bin/mpiexecjl`**; the system `mpiexec` fails on a PMIx mismatch.
    `MPI_Finalize` prints a benign OFI error and exits 143 on this machine.
22. **Keep horizontal cells near-isotropic, or pay for it in the linear solve** — a 4:1-celled mesh
    needs ~760 GMRES iterations against ~480 for an isotropic one.
23. **Use MPI when the problem is big enough to give each rank a meaningful share, and a direct LU
    otherwise.** Spend spare cores on more *cases*, not on decomposing one small case further.
24. **`nl_tol` is a step function and DOES change the answer** (integer Newton counts); `ls_rtol`
    does not, anywhere in 1e-9…1e-5. Do not conflate them. (`CONFIGURATION.md` §4)
25. **A sysimage is only valid for the versions it was built against** — build it in the environment
    you run in, and rebuild after **any** `src/*.jl` edit. Staleness is *detected*, not prevented.
26. **Julia buffers stdout when redirected to a file**; use `flush` or poll.
27. **Revise does not hot-swap signature changes** — restart before trusting a number after any
    signature or struct change.

### Long campaigns (earned 2026-08-21…30, the vertical-basis campaign)

41. **A LONG-LIVED Julia+Gridap PROCESS DEGRADES TO USELESSNESS.** Measured: workers grew
    1.5 → 2.5 → 3.9 GB over ~14 h and, left for days, fell to **7-14 % CPU** — GC-bound, not
    compute-bound. One spent **6.2 days on an `nx=8` level that takes seconds when fresh**, which is
    why a campaign sat at 25/36 for a week while "still running". **Bound worker lifetime and RSS and
    have a supervisor restart them** (`run_vbasis_shard.jl` + `supervise.sh`); resume logic makes a
    restart cost one JIT. ⚠ The cap is checked BETWEEN studies, so a single long study can still get
    there. This also **unseats the "flat memory ⇒ H3 not supported" argument** in `OPEN_ISSUES.md` §1,
    which was measured over 400 steps on one tiny case.
42. **`Distributed`/`pmap` was NOT usable for this workload; independent processes were.** Three
    `pmap` launches each completed ONE study in >4 h while the identical `run_mms_case` calls ran at
    full speed in a plain process. Sharded ordinary processes writing their own CSVs fixed it — and
    fixed observability too: under `pmap` the workers' prints are relayed through the master, whose
    stdout buffer is never flushed, so there is **no way to tell slow from hung**.
43. **CONTIGUOUS slicing of a COST-SORTED queue is the worst partition for makespan** (measured 12×
    imbalance: 1.7 h vs 20.8 h). Use LPT to minimise makespan for a fixed set, **SJF when the
    deliverable is coverage** — and beware that LPT schedules the cheap high-value studies LAST.
44. **In-flight work is invisible unless you make it visible.** A case is only "done" when it writes
    rows, so restarting slots re-run what others are mid-way through. Use **PID-keyed claim files in
    a SHARED directory** — per-output-dir claims let batches collide (one study ran on three slots).
    A claim honoured only while its PID is alive self-heals; a plain lock file would not.
45. **`RETRY_ERRORS` must be switched OFF once a failure is established.** A deterministic failure
    re-queued forever burns slots on work already understood — and these reproduce **bit-identically**
    (`‖r‖ = 0.27484014031572` twice), so a retry is not a new sample.

46. ✅ **REPRODUCING A PUBLISHED OPTIMISATION: THE THREE THINGS THAT WENT WRONG (2026-09-08).**
    `src/vopt.jl` reproduced Yang & Liu Table 1 to only 0.0085 at `M=2`. Three defects, none of them
    numerical, and the order in which they were found is itself the lesson.

    **(a) The median was over the wrong population — this was the big one.** Each of the five terms
    is *already* `kd`-integrated, so it is one number per candidate mesh and its median can only be
    over *meshes*. The code took it over a Kronecker low-discrepancy sweep of the **width simplex**;
    the paper takes it over **the scan it optimises on** ("c₂ … from 0 to 1 with an increment of
    0.001"). The medians are the **trade-off weights**, so a different population is a different
    optimisation problem, not a different approximation of the same one. Fixing it: 0.7365 → 0.7316
    against a published 0.728.
    ⚠ The LaTeX had a *third* reading — the half-mass point in `kd`, `∫₀^m E_X W = ½∫_Ω E_X W`.
    It matches at `M=2` and is **2.4–2.6× worse** at `M=3,4`. Do not reinstate it.

    **(b) `E_cg` was squared.** (3.10) takes `C_g` to the **first** power; `E_c` alone is squared.
    The asymmetry is the paper's own and is easy to flatten by copy-paste — it had been flattened in
    both the code and the LaTeX.

    **(c) `E_shoal` had a spurious `abs`.** (3.10) carries **no** `|·|` there (verified at 400 dpi;
    the next line shows `|u_m−u_e|` with unmistakable bars). ⚠ **The signed form is correct AND
    self-consistent, for a reason that is easy to "fix" wrongly:** `exp(I)−1` is single-signed across
    the design population (`γ_m > γ_e` throughout), so its median is negative too and **the sign
    cancels in the quotient** — the normalised term is positive and identical to putting a modulus
    around the whole term. This holds **only** if the median is over the `exp(I)−1` population.
    Normalising by `median(exp(I)) ≈ +0.81` breaks the cancellation, inverts the term's role, and
    drives the optimum to the scan edge (`c₂→0`). *That was a live defect I introduced while fixing
    (c), and it produced a confident, entirely wrong conclusion that the paper's equation was
    pathological.*

    **`M ≥ 3` DOES NOT REPRODUCE THE PUBLISHED NODES, AND THAT IS SETTLED.** Measured, at fixed band:
    `M=3` over `Δc` = 1e-2 → 1e-3 (4 851 → 498 501 designs) drifts **3e-4**, the last two levels
    identical to 4 dp; `M=4` over 3e-2 → 5e-3 (4 960 → 1 293 699 designs) is **unchanged to 4 dp**.
    The gaps (1.9e-2, 2.6e-2) are 60–100× the drift. The reference specifies its population only for
    `M=2`, so **the `M≥3` interfaces are not recoverable from what was published**. Ours are a
    defensible member of that family from a stated, reproducible population.

    **κ WAS A FIT READ OFF ITS OWN FIT — deleted, do not reintroduce.** `VOPT_KAPPA = 2.17` closed
    `Ω` for `p ≥ 2` via `Ω·Δσ_top/p = κ`. But `Ω` had first been *tuned* to reproduce Table 1
    (`Ω* = 7.60, 29.0, 100.0`) and the near-constancy of the product noticed afterwards. It is also
    35 % from the analytic `Δσ_top ≈ 2.94/kd_max` the LaTeX derives from the same boundary-layer
    argument. **A constant obtained by fitting cannot then validate the thing it was fitted to.**
    ⚠ The "`M=3,4` are best fitted by exactly 1.25× the published band" coincidence was an artefact
    of the same circularity — discard it.

    **Cost:** the scan is `C(1/Δc − 1, M−1)` designs, so the paper's `Δc = 1e-3` is exact at `M=2`
    (999) and impossible at `M=4` (1.66e8, ~8 days). Use `Δc` = 1e-3 / 5e-3 / 3e-2 for `M` = 2/3/4:
    converged optima to 1e-4 in minutes. **One notebook run cost 461 min for numbers identical to
    four decimals.**

47. ⚠ **THIS WORKSTATION LOSES RUN EVIDENCE IN THREE WAYS, ALL LEARNED 2026-09-23.**
    * **`/tmp` IS EMPTIED ON EVERY BOOT** (`tmpfiles.d` carries `D /tmp`). Long solver runs launched
      from an agent session write stdout into a scratchpad under `/tmp`, so **a reboot destroys
      every run log** — and the log is the *only* place holding the config banner, the per-step
      Newton trace and **the reason a run stopped**. `diagnostics.csv` records none of that; it
      simply ends. **Redirect batch stdout into `output/…` directly, never to a scratchpad.**
    * **A GNOME/Wayland session cannot restart its shell in place.** `org.gnome.Shell@wayland.service`
      *is* the compositor; restarting it is a logout to GDM, not a reload. ⚠ **And with `Linger=no`
      and a single session, ending that session stops `user@<uid>.service` and kills `app.slice` —
      i.e. every detached solver dies with the GUI**, `nohup` notwithstanding, because they live in
      `vte-spawn-*.scope` units under that slice. **`loginctl enable-linger <user>` decouples them**
      and is the precondition for touching the desktop at all (enabled here 2026-09-23).
    * **SUSPEND IS SAFE; REBOOT IS NOT.** Suspend preserves both the running solvers and `/tmp`.
      A reboot loses both. When a GUI freeze forces the choice, **suspend**, and persist logs first
      either way.
    ⚠ **A traceback ending in `paraview_collection`/`createpvd` is NOT a VTK failure** — that is the
    enclosing `do`-block in `run_time_loop`. The real cause is higher up; use
    `grep -m1 '^ERROR' run.log`. This misled the first reading of four crashed arms.

### Testing and measurement

28. **"The suite passes" is not "the model is verified."** Most of the suite is
    **self-consistency**: writing `R = R_true + E`, the error appears on both sides and cancels
    identically, so such a test passes for **any** residual. Only the analytic MMS — whose forcing
    never touches `problem.jl` — can detect a self-consistently wrong residual. (`VERIFIED_SCOPE.md` §0)
29. **AD is an oracle for the JACOBIAN, never for the RESIDUAL.** "AD converges where hand fails,
    therefore the residual is correct" is an invalid inference — AD differentiates the *same*
    assembled residual.
30. **A BOUNDS CHECK ON THE RIGHT CONFIGURATION IS NOT A VALUE CHECK.** A gate that only asserts
    nothing exploded will pass while the quantity it computes moves 58 %.
31. **For a guard that is a CONJUNCTION (`nonlinear ∧ ∇h ≠ 0`), the suite needs a case satisfying the
    conjunction AND asserting a value.** Conjuncts satisfied separately, or together but only
    bounded, both look green.
32. **A refinement study measures the rate of whichever error DOMINATES — verify isolation in BOTH
    directions before interpreting any slope.** A saturated slope and a wrong coefficient produce the
    *same* observable. **Guards come in pairs.**
33. **Read the pairwise rate SEQUENCE, not the fitted slope**, and check error *magnitude* before
    trusting a fine-level high-order rate. **Three independent instances in one campaign** where the
    fit would have given the wrong conclusion and the sequence gave the right one (fit 3.515/3.729/
    3.603 against a true 3.94/3.88). The diagnostic that generalises: **a rising sequence with the
    error still dropping is PRE-ASYMPTOTIC; a flat sequence with a stalled error is not.** And a
    *fitted* slope through points that were never on an asymptotic curve (P2LFE-2 M4: `e_u` starts
    LARGER than the solution amplitude) is **meaningless, not low** — never report it as a rate.
34. **THE RESOLUTION PRINCIPLE: a test validates a term only if it can RESOLVE that term's
    contribution.** Always ask: *if this term were wrong, would this test have noticed?*
35. **A batch runner must take its verdict from GATE OUTPUT, never from exit codes.** A clean exit
    code is not evidence a test ran.
36. **A mirrored guard that collapses two per-field measurements into one pass/fail can hide the
    very asymmetry it exists to detect.**
37. **Two tests sharing a bathymetry function and a physics tier can still be different problems.**
    A reference-vs-reference comparison is meaningful only when *every* discretisation parameter
    matches.
38. **When an error message names a library type, that is where the bug SURFACED, not where it
    lives.** Print the type of *every* input to the failing expression before searching that library.
38b. **REFINEMENT THAT MAKES THINGS WORSE IS A DISCRETISATION-STABILITY PROBLEM, NEVER
    UNDER-RESOLUTION.** Under-resolution improves with `h`; a grid-scale mode does not. This single
    check separated the two whole classes of explanation faster than anything else in the
    investigation, and it is the first thing to measure whenever a nonlinear run diverges.

38c. **COPY THE REFERENCE ENVIRONMENT VERBATIM AND VARY ONE VARIABLE EXPLICITLY.** `env -i` plus a
    partial list of variables silently changed wave generation from `:bc` to the interior source,
    which delivers **2.8×** the requested amplitude — so the treatment ran at triple the intended
    forcing and produced a confident result pointing the wrong way. **And put the null-treatment
    control in the SAME batch**: without it, a difference from a remembered baseline cannot be
    attributed to the treatment.

38d. **VERIFY A NEW KNOB ACTUALLY DOES SOMETHING BEFORE SPENDING THE RUN ON IT.** A dead parameter
    yields identical curves — a clean, confident, entirely wrong negative result. A seconds-long
    unit check (here, integrating `x⁸` at each quadrature degree) buys the whole experiment.

38f. **IDENTIFY AN INSTABILITY BY WHAT GROWS, NEVER BY WHAT IS LARGEST.** Ranking a spectrum by
    *amplitude* found the carrier's own harmonics at `λ≈1.7 m` — a wavelength independent of `dx`,
    which reads as evidence *against* a grid mode. Ranking the same data by *gain* found the real
    mode at `λ≈2–3·dx`. An unstable mode is small for most of its life; that is why it goes unnoticed
    until it dominates. **Corollary: when a reported extremum sits at the edge of the search window,
    the window is the finding** — a `0.95·k_Nyq` cut put one case's "peak" exactly on the boundary.

38e. **A REVERSAL TOO LARGE TO BE THE EFFECT UNDER TEST IS A BUG SIGNAL, NOT A FINDING.** Extra
    quadrature turning a 0.11 plateau into 1.60 is not a plausible quadrature effect; that
    implausibility is what prompted the check that found the real cause.

38g. **A NULL RESULT IS NOT A FINDING UNTIL THE EFFECT SIZE IS MEASURED AGAINST THE EFFECT BEING
    EXPLAINED.** Rule 38d says verify the knob is LIVE; this is the other half. The 2026-09-12
    quadrature probe returned `e_η` columns that were **bitwise identical** between treatment and
    control at every level — which reads as "dead knob, result void" (I called it that) and equally
    as "refuted". Both readings were wrong: the knob was live, and its effect simply fell below the
    CSV's 7-significant-digit print precision. **Identical printed output distinguishes nothing**;
    only assembling the quantity under test at both settings does. When it was assembled, the crime
    was real (3.5e-08 in the residual) and six orders too small to explain a 3.1e-07 rate deficit —
    which IS the finding, and is the same shape as the skew-advection refutation in rule 12b.
    **Corollary: print precision is part of the experiment.** A comparison whose resolution is set by
    `%.6e` cannot detect anything below 1e-7 relative, however carefully the runs were controlled.

38h. **THE LAUNCHER IS PART OF THE EXPERIMENT. Three ways a run silently became a different run,
    all on 2026-09-15, all costing hours.**
    * **An override appended AFTER the `balfem_local_run` line is a NO-OP.** A diagnostic built that
      way ran the production configuration instead, and reported it under the diagnostic's name.
      The tell was rule 38d's: η, Newton counts, `r0` and the final residual **identical to every
      printed digit** against the control. An exact Jacobian cannot reproduce a quasi-Newton
      iteration count to the digit — *check that the knob moved something before believing the run*.
    * **NEVER EDIT A SHELL SCRIPT WHILE BASH IS EXECUTING IT.** Bash reads a script incrementally by
      byte offset, so rewriting the file under a running interpreter shifts the content beneath it:
      the fixed launcher's relocated run-line was re-executed by the *old* process, silently starting
      a second solver that competed for RAM for an hour before anyone noticed. Write a NEW file.
    * ⚠ **`exit status 1` DOES NOT MEAN THE RUN FAILED.** `examples/local_1d/run_flume_1d.jl:350`
      ends with `@printf(..., tag, ...)` and **`tag` is defined nowhere** — a leftover from the
      `output_dir_name` rename (rule 2c). Every *successful* run therefore throws `UndefVarError`
      after writing every solve result, VTK file and diagnostics row. A batch runner that reads exit
      codes will score a completed 100-period run as a failure (rule 35). One-line fix, not yet made.

38i. **TO SEPARATE AN INTERIOR INSTABILITY FROM A BOUNDARY FEED, REMOVE THE BOUNDARY.** A failure that
    pins at the inflow cannot be attributed by any amount of flume refinement, because every flume run
    carries the inflow, the relaxation zone and the sponge. The **closed x-periodic box** (one model
    wavelength, no inflow, relaxation, sponge or source, deterministic 1e-8 seed, sub-cell sampling,
    band energies relative to the carrier) has a known null — nothing grows — and it measures the growth
    over full wave cycles (a Floquet-type rate, which frozen eigenvalues cannot give). Run it with its
    `:native` twin **and** under Crank–Nicolson (rule 15). It settled in three days what two weeks of
    flume factorials could not (§5.2f).
39. **Test a diagnosis against a case it cannot explain, rather than looking harder where it points.**
39b. **A CONSTRUCTION THAT *CANNOT* COMMIT THE SUSPECTED ERROR IS THE STRONGEST DISCRIMINATOR
    AVAILABLE — AND ITS NULL RESULT IS A REFUTATION, NOT A DISAPPOINTMENT.** Rule 39 says test a
    diagnosis where it cannot reach; this is its constructive form. For a year the `dx` refinement
    signature of `:full` was attributed to **recovering `∂²` of a `C⁰` field** (the skeleton Dirac
    layer that cell quadrature never sees). The mixed formulation was built to remove the *lag*, but
    it also **structurally cannot form that second derivative at all** — `𝖦 ≈ ∇𝖲` is a genuine
    unknown defined by an integrated-by-parts weak equation. **The `dx` signature survived it
    unchanged.** That single observation refutes an explanation that no amount of further refinement
    studies on the projected path could have touched, because every such run shares the error being
    blamed. ⚠ **Read the consequence honestly and immediately:** the refutation also invalidates the
    *remedy* queued behind it — `C⁰` interior penalty exists to make a broken `∂²` consistent, and a
    signature measured where there is no `∂²` to break cannot justify it. **When a structural fix
    fails to move a signature, retire the hypothesis AND everything scheduled on its authority**, in
    the same edit, before the queued work is launched on a dead rationale.
40. **Diagnostics interpretation:** never read `growth` without `x_at_max`; `dmp/int` means "≫1 is
    trouble", never "<1 is fine"; mass drift is an invariant **only** in a closed, unforced basin.
    (`ARCHITECTURE.md` §6)

---

## 8. Conventions for editing here

* **`latex_docs/BALFEM_models_v2/` is the single source of mathematical truth** (the v1 document is
  frozen). `MODEL.md` mirrors its notation.
  If you change the maths, change both.
* **Notation bridge, the most common source of confusion.** LaTeX `M^V, 𝓜^V, 𝓖^V, A^V, K^V, 𝓐^V,
  𝓚^V, Φ, φ_j` ↔ code `Mmat, Mcal, Gcal, A, K, Acal, Kcal, Phi/D/C`, with the unit basis stored as
  `w_j = −φ_j_int`. Always cross-reference `MODEL.md` §2 when quoting tensor names.
* **If you add a `using X` to any `src/*.jl`, add it to `[deps]`.**
* **`markdown_files/` IS TRACKED** (16 files, since 2026-09-11) — put load-bearing design records
  there and cross-link them from this file. ⚠ **`latex_docs/` is the gitignored tree**, so nothing
  that matters may live only in the LaTeX. (This bullet previously said the opposite about
  `markdown_files/`; it was stale and is corrected.)
* **Run Julia via the `julia-mcp` tool.**
