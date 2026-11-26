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
> `latex_docs/` is **gitignored** (it only holds a link to the LaTeX project, §2), so anything
> load-bearing goes in `markdown_files/`.
>
> **Sister solver (split 2026-11-23).** The mixed formulation of the BALFE-M models (auxiliary
> unknowns for the 𝓝 gradients, conforming without CIP) with residual-based (VMS) stabilisation is
> its own repository, [`GridapBALFEMixed.jl`](https://github.com/Elmanyer/GridapBALFEMixed.jl). Its
> plans, results and conclusions live there, not here. Its development history in this repository is
> kept by the tag `archive/mixed-vms-solver-pre-split`. This repository is the primal formulation with
> the C⁰ interior penalty (CIP) only.
>
> **v1 and v2 are separated (2026-10-05).**
>
> | | v1 (frozen) | v2 (active) |
> |---|---|---|
> | code | branch `v1-solver`; tag `v1_final_solver` | branch **`main`** (v2 since 2026-11-24; `v2-solver` merged in full 2026-11-26) |
> | LaTeX | `../LATEX_DOCS/BALFEM_models_v1/` | `latex_docs/GridapBALFEM_documentation/` (Overleaf, §2) |
> | outputs | `../output_v1/` (workspace root since 2026-11-26; index `CATALOGUE.md`) | `output/` |
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
| `src/` (20 files) | `problem.jl`: `BALFEMProblem`, `resolve_physics`, the loop-free residual and hand Jacobians · `nlpressure.jl`: `𝓝` {3,6,7,8} and the ∇h half · `broken.jl`: Class III cellwise Hessians + kind-S facet terms + exact Jacobian, the C⁰-CIP `:jumpgrad` (**mandatory with 𝓝**) and `:ghostvolume`; sequential **and distributed** (ghost sequential only) · `utilities.jl`: `setup_and_run`, `output_dir_name`, `check_v1_env`, `write_run_manifest` · `configuration.jl`: `write_run_configuration` (`CONFIGURATION.txt`) · `vopt.jl`, `mms.jl`, `vertical.jl`, … (`ARCHITECTURE.md` §2) |
| `test/` | `runtests.jl` judges by gate output (never exit codes); `test/local/`, `test/cluster/`; `test_skeleton_ad.jl`; `test_nlp_distributed_parity.jl` (models 5–6, sequential ↔ distributed, debug or MPI backend); `test/v2_migration/regression_snapshot.jl` (v1→v2 gate, recording in `output/v2_baseline/`) |
| `examples/` | `local_1d/` (`run_flume_1d.jl`, `run_periodic_1d.jl` the closed box, `stability_eig.jl`, `cip_eigen_analysis.jl`, `cip_eigen_order_rule.jl`), `local_2d/`, `local_mms/`, `distributed/`, `validation/` |
| `run/` | `local/balfem_local.sh`, `balfem_env.sh` (cluster), SLURM launchers, `SNELLIUS_ROME_LAUNCH_CONFIGS.md`; `linear_production/` (v2 linear Snellius suite, single-node 1×42 directional); v2 box launchers `run/local/run_1dper_v2_unstab{,_extra,_q32}.sh` (⚠ their nlp1 rows are historical: 𝓝 now carries the mandatory CIP); **`integrator_study/`** (the time-integrator box study: `cases.py` → ONE `cases.tsv`/`packs.tsv` (275 cases), ONE `submit.sh` (tiers, `status`, `requeue-quota`, `resume`, `BALFEM_CASE_FILTER`), `run_pack.sh` with checkpoint segments, `cluster.sh` = per-cluster sbatch flags incl. the memory request, `BALFEM_CLUSTER=snellius\|blue`; `analysis/` = the 2026-11-23 analysis scripts); **`storage.sh`** + **`setup_storage.sh`**: cluster output on the project storage (`STORAGE` in `~/.bashrc`; `output/` → `$STORAGE/GridapBALFEM_output`; checked by `balfem_env.sh` and every submit script) |
| `compile/` | cluster sysimage chain (`RUNNING.md` §5): `compile_snellius.sh` and `compile_blue.sh` (DelftBlue compute-p2, `load_modules_blue.sh`: OpenMPI 5.0.10 + juliaup Julia 1.10); `set_preferences.jl` adds the forks by URL and binds MPI to the system library |
| `output/` · `../output_v1/` | gitignored. `../output_v1/` (moved to the workspace root 2026-11-26) = v1 only, every run labelled with a reconstructed `CONFIGURATION.txt`, index `../output_v1/CATALOGUE.md`. `output/` = v2 only, **catalogued 2026-11-26: index `output/CATALOGUE.md`, every run labelled with `RUN_LABEL.txt` (stage, code state, outcome, validity), tools `output/_catalogue/`**; every run writes `run_manifest.toml` (machine record) and `CONFIGURATION.txt` (readable identity card: a one-line `ID`, model, physics with the derived kd/kA/λ/dispersion errors, discretisation, integrator with its carrier damping, output; same sections in every run; 2026-11-11) |
| `postprocessing/` | `GridapBALFEMPost`, its own environment; `examples/periodic_growth.jl` gives the band growth rates; `examples/reconstruct_3d_vtk.jl` writes 3-D `.vts` series (u, v, w, p, p_nh, p_dyn) from a run's stored modes for ParaView |
| `WaveSpec.jl/` | vendored sea-state synthesis (GitHub version) |
| `Gridap.jl/` | the Gridap **fork** (`Elmanyer/Gridap.jl` @ `fix-transient-multifield-ad`, v0.20.8) with two AD patches: transient multifield (`fa860899c`) and skeleton integrals on transient fields (`3928b98b9`, **not pushed**; loaded locally via `Pkg.develop`) |
| `latex_docs/` | gitignored; only the link `GridapBALFEM_documentation` → `../LATEX_DOCS/GridapBALFEM_documentation` (§2). The other LaTeX projects (v1, the superseded combined `BALFEM_models_v2`, the derivation, the sister solver's documentation, `CFC2027_abstract`, `doc_figures`) are siblings in `../LATEX_DOCS/` |

## 2. The LaTeX project

* **`latex_docs/GridapBALFEM_documentation/` is the single source of mathematical truth** (since
  2026-11-23; a link to `../LATEX_DOCS/GridapBALFEM_documentation`). It is its own git repository
  (remote **`origin`** = its Overleaf project). **The author edits it in Overleaf**, so always run
  `git fetch origin` and fast-forward before editing; the author may stash local edits when pulling.
  * **Chapters 1–5 (the model derivation) are not written here.** They are the folder
    `BALFEM_derivation/`, which comes from the separate Overleaf project
    `../LATEX_DOCS/BALFEM_models_derivation` (shared with the sister solver's
    `GridapBALFEMixed_documentation`). Edit those chapters in that project; `main.tex` here reads them
    with `\input{BALFEM_derivation/…}` and also takes its commands (`newcommands_config`),
    `hyperref_config` and bibliography.
  * `../LATEX_DOCS/BALFEM_models_v2/` (remote `overleaf`, last commit 2026-11-09) is the superseded
    combined document; do not edit it.
* **Chapter map** (`LATEX_STRUCTURE.md`):
  * model 1–3;
  * checks without a solver 4–5;
  * horizontal weak form 6: H² conformity, the 𝓛/𝓝 pressure blocks, C⁰-CIP (§6.4, incl. §6.4.5 accuracy
    of the cellwise evaluation), the global residual `eq: multifield general residual`;
  * multi-field implementation 7 (the complete model only);
  * multi-fidelity implementation 8 (2026-11-09: the switches `regime`/`nl_pressure`/`flat_bed`, the six
    models, the term classification with activation conditions, their code realisation; moved out of 6–7);
  * stability 9;
  * stabilisation 10;
  * solver functionalities 11;
  * validation 12;
  * appendices A–C, which are v1 records.
* **Chapter 6 (rev. 2026-11-06):** Class III of blocks C and L is evaluated **cellwise** and completed
  by the kind-S **distributed-Hessian facet terms** `ℱ^C`, `ℱ^L` (§6.5 `sec: cellwise accuracy`; stacked
  form `eq: multifield general residual facet terms`, §6.7). They are **in the global residual
  `eq: multifield general residual`, marked optional**: dropping them stays consistent but costs one
  order (and at `p_η = 1` the η curvature). The **first-order C⁰-CIP on η and every 𝖴 is part of the weak
  form**. The reduced weights `𝒲^K`, `𝒲^P` carry component 5 **transposed** (`T⁽⁵⁾_{ijk}`, as in the code
  and `eq: class III reduction`). Chapter 7's implemented residual = chapter 6 with the facet terms always
  on (`eq: implemented global residual`, `tab: implemented vs weak form`). The broken identity
  `eq: broken integration by parts identity` is the tool that classifies facet terms (kind 0/S/T).
* Build in a scratch copy: `pdflatex`, `biber`, `pdflatex` ×2 → 282 pp, 0 errors. Chapters 6, 7, 10 have
  no undefined references; chapters 8–9, `ModelVerification`, `StokesWaveFourierAnalysis` and the
  appendices still cite removed broken-audit labels (they need rewriting, not redirecting).
* ⚠ Never global-replace `H^2`/`L^2`: it hits the depth `H²`, `h²` and `\partial^2` (33 lines to
  restore on 2026-10-07).
* Labels are `\label{chap:|sec:|subsec:|subsubsec: name}`; cross-references are `\S\ref{…}`.
* **Hyperlinked (2026-10-30).** `hyperref` (+ `bookmark`, `hypertexnames=false`) is the LAST package
  in `main.tex`; biblatex has `backref=true`. Every `\ref`/`\eqref`/`\cite` and TOC entry is a link,
  bookmarks are numbered. Keep it that way:
  * **every** chapter/section/subsection/subsubsection carries a `\label` (all 70 missing ones added);
  * a heading with math or `\cite` writes its bookmark text with
    `\texorpdfstring{$\opNkj$}{N\_kj}` (literal UTF-8 Greek, `\texttwosuperior`, `\textzerosuperior`;
    `\textsigma` etc. are NOT available);
  * refer to a chapter by name with `\hyperref[chap: …]{text}`; never hand-type a number for an
    internal target (the hand-typed ones left are all Yang & Liu's).
* If you change the maths, change `MODEL.md` too.
* Notation bridge (`MODEL.md` §2): LaTeX `M^V, 𝓜^V, 𝓖^V, A^V, K^V, 𝓐^V, 𝓚^V, Φ` ↔ code
  `Mmat, Mcal, Gcal, A, K, Acal, Kcal, Phi/D/C`; the unit basis is stored as `w_j = −φ_j_int`.

## 3. Features

* **Solver core.**
  * The stacked `[η,𝖴x,𝖴y]` loop-free residual with hand Jacobians: `∂R/∂u̇` exact; `∂R/∂u` exact for
    models 1–4 (η-derivatives of the pressure terms completed, `a0282bc`), quasi-Newton only in the `𝓝`
    blocks {3,6,7,8} and the ∇h half of models 5–6 (0.5–0.65 % gap, rule 5).
  * Integrators (`build_time_integrator`, 2026-11-11): `:theta` (θ; 0.5 = CN), `:sdirk` = `:rk` (any
    Gridap tableau, diagonally implicit or explicit `EXRK_*`; default `SDIRK_2_2`, **dissipative, rule
    15**), `:gen_alpha` (ρ∞). Refused before the build: IMEX tableaus, unknown names, θ ∉ (0,1],
    ρ∞ ∉ [0,1]; warned: defective tableaus and those that amplify neutral oscillations.
    ⚠ **Explicit tableaus are unusable with 𝓝** (2026-11-12): every one failed on m5 at T/40 in the
    smoke test. Explicit RK only for models without CIP. (The −9.5e3 s⁻¹ / Δt ≈ 2.6e-4 s limits
    once quoted are void: `stability_eig.jl`'s FD pattern missed CIP entries until 2026-11-13,
    `STABILITY.md` §6b.)
  * Sequential LU+Newton; distributed GMRES+Jacobi+Newton, **all six models** (𝓝 distributed since
    2026-10-30; facet size and count built per rank).
* **Physics.** Advection, the full leading pressure `R_P`, and `𝓝` with all eight components when
  `nl_pressure=true`:
  * `{3,6,7,8}`: direct assembly;
  * `{1,2,4,5}`: the ∇h half (block B) by exact integration by parts, and the **Class III** half
    (blocks C, L) by cellwise `∇∇` plus the kind-S facet term (`broken_class3_residual`), with an
    exact Jacobian;
  * **the first-order C⁰-CIP is mandatory** (2026-10-30): `build_problem` attaches it on both fields at
    `DEFAULT_CIP_GAMMA = 0.3`; `attach_skeleton!`, `global_residual` and `jacobian_u` refuse 𝓝 without
    it, and `p_eta < 2` is refused;
  * the `{1,2,5}` reduction is exact to 4.4e-16.
  * ⚠ The leading-pressure skeleton term `⟦𝒫⟧` is excluded on purpose (`h⁻⁴`-conditioned mass).
* **Skeleton penalties** (`attach_skeleton!`), one per run:
  * `:jumpgrad`: normal-derivative jumps of order ≤ `cip_order` ≤ 2 (first order = the C⁰-CIP of
    LaTeX §6.4, `γ τ h_F² ⟦∂ₙ·⟧⟦∂ₙ·⟧`, `τ_u = d√(gd)`, `τ_η = √(gd)`, `M^V`-weighted); sequential and
    distributed;
  * `:ghostvolume`: direct ghost penalty, every order; uniform Cartesian meshes, sequential only.

  Knobs `cip_gamma_u/_eta` (`nothing` ⇒ 0.3 with 𝓝, 0 otherwise; the manifest records the resolved
  value), `cip_hexp`, `cip_order`. Env: `BALFEM_CIP_GU/_GE` (unset ⇒ the default).
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
  * **3-D reconstruction for ParaView** (`examples/reconstruct_3d_vtk.jl <run> [nσ] [every]`,
    2026-11-14): w and p (p_nh, p_dyn) rebuilt from the stored velocity modes and η of a 2-D run, on
    z = −d + σ(d+η). Default: the σ-nodes of the vertical mesh only (`<run>/recon3d/`); `nσ ≥ 2` adds
    uniform interpolated levels (`<run>/recon3d_s<nσ>/`). ∂ₓ by 4th-order FD on the vertex grid, ∂ₜ by
    central FD across snapshots (coarse output cadence loses high frequencies); p_nh is the flat-bed
    formula.

## 4. Physics selection

| control | values | meaning |
|---|---|---|
| `regime` | `:linear` \| `:nonlinear` | linearised core / full nonlinear core with advection |
| `nl_pressure` | `false` \| `true` | `𝓝` off / **all eight components** (requires `:nonlinear` and `p_eta ≥ 2`; carries the mandatory C⁰-CIP; refuses a `Symbol`) |
| `flat_bed` | `Bool` | `∇h ≡ 0`: every ∇h term dropped |

* Six models:
  * 1–2: linear (flat / variable bed);
  * 3–4: nonlinear, `𝓝` off;
  * 5–6: nonlinear, `𝓝` on.
* `𝓝` is all-or-nothing because no amplitude or `kd` ordering separates its components
  (`V2_SOLVER_PLAN.md` §0).
* Numerical options:
  * `stabilization`, `cip_*`;
  * `solver_type` (**diagnosis under `:theta`; dissipative remedies state their carrier cost**, rule 15);
  * `use_ad` (diagnosis only).
* Environment variables: `BALFEM_NL_PRESSURE`, `BALFEM_STAB`, `BALFEM_CIP_*`, `BALFEM_SOLVER`
  (`theta|sdirk|rk|gen_alpha`), `BALFEM_TABLEAU`, `BALFEM_THETA_T`, `BALFEM_RHO_INF`, `BALFEM_DT`,
  `BALFEM_OUTDIR`. The output name carries the integrator (`theta`, `theta0.60`, `galpha0.50`, a
  tableau name); `test/local/smoke_integrators.jl` runs every choice through the box driver.
* `check_v1_env()` refuses the v1 knobs.
* Output names carry `nlp0`/`nlp1` (`output_dir_name`, rule 2c).

## 5. Current Implementation Stage

*Updated 2026-11-26. Full picture: [`markdown_files/STATUS.md`](markdown_files/STATUS.md) (its 2026-11-11
update is the current route); actions: [`CH8_CAMPAIGN_TODO.md`](markdown_files/CH8_CAMPAIGN_TODO.md)
"RESUME HERE".*

### Repository and records (2026-11-26)
* **`main` is the v2 solver** (merged 2026-11-24, in full 2026-11-26); `v2-solver` is deleted, v1 is frozen
  on `v1-solver`. ⚠ The Snellius clone still tracks `v2-solver`: `git fetch --prune && git switch main &&
  git pull`, then rebuild the sysimage (`src/problem.jl` changed with `a0282bc`).
* **All output is catalogued.** `output/CATALOGUE.md` (637 v2 runs: stage v2.0–v2.d, code state, outcome,
  validity; `RUN_LABEL.txt` in every run, reconstructed `CONFIGURATION.txt` for runs before 2026-11-11) and
  `../output_v1/CATALOGUE.md` (376 v1 runs). Rebuild with `output/_catalogue/` after new runs land.
* **`jacobian_u` complete in η** for the nonlinear models (`a0282bc`, 2026-11-18, on `main` since 11-23):
  models 1–4 hand ≡ AD, `test_jacobians_ad` 14/14. Every nonlinear run in `output/` predates it (rule 17b:
  stability unaffected).

### Working and verified
* **The v2 code base:** Class III cellwise + kind-S facet terms; Boolean `nl_pressure`; v1 code deleted.
  The implemented residual matches LaTeX chapters 6–7 term by term (re-read 2026-10-29, signs included);
  v1→v2 equivalence entry by entry (12/12).
* **Linear models 1–2** at optimal order: 1-D and 2-D, sequential and distributed, five vertical bases.
  **Nonlinear models 3–4** MMS to `nx = 32`. The Yang & Liu collapse.
* **Mandatory first-order C⁰-CIP with 𝓝** (2026-10-30, author): default γ_u = γ_η = 0.3; refused
  without it or with `p_eta < 2`; every driver records the resolved γ in the manifest.
* **The full nonlinear model distributed** (skeleton facet terms and `:jumpgrad` on
  `DistributedDiscreteModel`; `:ghostvolume` stays sequential). The depth is built with the ghost layer
  (`depth_cf`): **any CellField a skeleton term reads on both sides must include ghosts.**
  * Gates: `test_nlp_distributed_parity` 14/14 (debug backend and `mpiexecjl -n 4`, ≤ 2e-12);
    `test_broken_formulation` 40/40; `test_jacobians_ad` 14/14; regression snapshot 12/12;
    `test_skeleton_ad` S0, S1, S2 `:jumpgrad` pass.
* **Time integrators** (`build_time_integrator`, 2026-11-11): θ, any Gridap RK tableau, generalised-α;
  the carrier damping of the chosen scheme is printed in `CONFIGURATION.txt`.
* **Run records:** every run writes `run_manifest.toml` and `CONFIGURATION.txt` (2026-11-11).

### Chapter 8 — stability of the CIP-mandatory solver (1-D closed box)
**CN campaign done 2026-11-07** (`STABILITY.md` §3b; record `output/1D_stability_campaign/README.md`).
* Models 5 and 6 run 200 periods clean at 16 and 64 cells/λ. At 32 cells/λ a late burst (245 s flat
  bed, ~295 s cosine bed); at 8 cells/λ a slow long-wave growth drains the carrier. dx dependence is
  not monotone. Model 3 diverged at n32 (246 s) and grew at n16 and n8.
* **Every growing content is a vertical-structure velocity mode** (layer velocities grow, ū and η
  barely move): a stabiliser must act on each u_j.
* **No-CIP controls:** m6 without CIP failed (Newton NaN) at 8.4 / 25.8 / 74.8 s (n64 / n32 / n16),
  as m5 did; at n8 the no-CIP m5 and m6 reached 320 s bounded (E/E₀ 1.04 / 1.01), while m5 n8 *with*
  CIP shows the late growth (not yet identified per harmonic).

### Chapter 9 — stabilisation: integrator dissipation first (author, 2026-11-11)
* **Local SDIRK_2_2 box runs, m5 n32, Δt = 0.04, 200 periods** (finished 2026-11-17;
  `output/1D_stability_campaign/m5_q32_n32_sdirk{,_kA25}`, analysis
  `_logs/harmonic_breakdown_m5_sdirk_2026-11-17.log`):
  * no burst at kA 0.158 or 0.25, but the carrier kept only 2.8 % / 2.3 % of its energy (pure-sine
    prediction 3.2 %): "stabilised at a cost of ≈ 97 % of the carrier energy" (rule 15);
  * the CN burst's −2ω₀ component is gone; a **steady (ω ≈ 0) vertical-structure velocity growth**
    continues underneath at +0.005–0.014 s⁻¹ (bands; single harmonics up to 0.04 s⁻¹), the same pace as
    under CN before its burst. **The integrator cannot damp it** (ω = 0);
  * open: "SDIRK broke the mechanism" vs "SDIRK starved it of carrier" are not separated by these runs
    (rule 12c); the integrator study's Crouzeix 3/4 arm (damps 2ω₀, keeps the carrier) separates them.
* m3 n32 kA 0.25 under SDIRK_2_2 diverged at 172 s (pre-cluster check 1).
* **Cluster integrator study** (`run/integrator_study/`, 275 cases: tiers P–H × CN, SDIRK_2_2, SDIRK_4_3,
  Crouzeix 3/4, TR-BDF2, generalised-α ρ∞ 0 / 0.5; RK3/RK4 on m3) **ran on Snellius 2026-11-17 → 22,
  analysed 2026-11-23: [`INTEGRATOR_STUDY_RESULTS.md`](markdown_files/INTEGRATOR_STUDY_RESULTS.md).**
  * ⚠ **161 cases were killed by the Snellius HOME quota** (VTK output) at a median 49 % of their time.
    Their data is used, and their checkpoints are intact. The launchers now write to the project storage
    (rule 47c; set up on Snellius 2026-11-24); the resubmission plan is drafted (§6 of the record),
    **not submitted** (author).
  * **No integrator holds m5/m6 at kd 5.5, n32, kA ≥ 0.25 while keeping the carrier.** Survivors remove
    90–100 % of its energy. The gentle schemes (73–94 % kept over 200 T) delay the burst 1.2–1.9× at kA
    0.32; at kA 0.158 / 0.25 they were cut before CN's onset (open).
  * **The burst is a steady (ω ≈ 0) vertical-structure velocity mode at λ ≈ 2 cells** (n32 and n64),
    absent at n16 and faster at n64 (rule 38b). **A larger CIP γ brings it earlier** (CN onset 101 / 76
    / 60 / 51 periods at γ 0.1 / 0.3 / 0.5 / 1, kA 0.25): H2 falsified. H3: Δt, not γ. H1 (m3)
    confirmed at kA 0.158 only.
  * Under CN, m5/m6 n32 **saturate** after the burst (≈ 10 % spurious grid-scale velocity, 200 periods
    "completed"). kd 2.9 / 1.7 are clean under CN up to kA 0.25.
  * ⚠ Unexplained: on m3 the multi-stage schemes (SDIRK_4_3, RK3, RK4) die earlier than CN at kd 1.7
    kA 0.25.
* **Consequence:** the integrator route is closed for the kd 5.5 regime (`STABILITY.md` §6b result). The
  remedy must act on the steady grid-scale vertical-structure velocity at any frequency: a weak-form
  term on each u_j (open for this solver). The mixed formulation with residual-based (VMS)
  stabilisation is a separate solver, `GridapBALFEMixed.jl` (see the note at the top).

### Other open items (`OPEN_ISSUES.md`)
* Models 5–6 outside the verified MMS scope (`p_u` 3.94 → 1.97 on fine meshes, cause unknown);
  nonlinear `p_η` order reduction at Q3/Q2.
* Not yet run: a distributed *transient* run of models 5–6. Fast suite 9/9 not re-run since 2026-10-30.
* `test_skeleton_ad` ghost arm (`_shift`) and S3; `test_selfconsistency` step 12; the variable-bed
  lee-shoulder mode.
* LaTeX chapters 8–9 and appendices still cite removed broken-audit labels.
* Sponge sweep: stage A complete (84/84, 2026-11-08) and stage B analysed (`SPONGE_SWEEP_RESULTS.md`).
* `src/configuration.jl` reports any `WaveInput` (incl. a one-component regular wave) as a sea state in
  `CONFIGURATION.txt` (found 2026-11-26 while cataloguing; not patched, author's call).
* Gridap fork skeleton-AD commit unpushed; rebuild the sysimage after any `src/` edit (rule 25).

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
      layer, so typing `∇∇(u)` drops the layer. Dropping it is **consistent** (kind S) but costs one
      order: those terms converge as `O(h^(p−1))` instead of `O(h^p)`, and at `p_η = 1` the η-curvature
      part is absent at every resolution (LaTeX §6.4.5). v2 assembles it (`src/broken.jl`).
    * **Raising `fe_order` does not help:** Lagrange `H1` is `C⁰` for every `p`.
    * A test function's `∂²v` is computable cellwise.
    * The task is a *consistent* broken formulation, not a smoother space.

1c. **With `𝓝`, the first-order C⁰-CIP is part of the model, not an option** (author, 2026-10-30;
    LaTeX §6.4). Both fields, `γ > 0` (default 0.3), `p_η ≥ 2`. The code refuses anything else;
    never add a back door. Tests that need the operator without it subtract `stab_contrib` explicitly.
    The kind-S facet terms are a *different* object (linear in the jump, part of the operator); CIP
    (quadratic, symmetric) cannot replace them.

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

5. **`∂R/∂u̇` is exact; nonlinear `∂R/∂u` is quasi-Newton by choice only in the `𝓝` blocks**
   {3,6,7,8} and the ∇h half (the gap vanishes at order ≈ 1.1 in amplitude). The η-derivatives of the
   pressure terms are complete (`a0282bc`, 2026-11-18).
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

15. **Dissipation is a legitimate stabiliser; name it and measure it** (rewritten 2026-11-11, author).
    `SDIRK_2_2` is Gridap's `DIRK22(1,0,1)` (A-stable, not L-stable, `R(∞) = −½`; not the textbook
    SDIRK2). It removes ≈ `0.36(ωΔt)⁴` of a mode's amplitude per step (`0.75(ωΔt)⁴` of its energy):
    at Δt = 0.04, T = 1.6 s the carrier loses 0.0054 s⁻¹ in amplitude, ≈ 82 % over 200 periods.
    * **Diagnosis** (what the discretisation grows, and how fast) uses Crank–Nicolson
      (`BALFEM_SOLVER=theta`): it neither adds nor removes energy, so nothing hides the growth.
    * **Remedy** runs may use any dissipation (integrator, CIP, filter). Each reports the carrier's
      amplitude loss over the target duration (the accuracy cost), the damping at the growing modes'
      frequencies or wavenumbers, and the CN run it is compared with. A dissipative run is
      "stabilised at cost X", never "stable" without that number.
    * Integrator damping acts by **frequency**; it cannot reach steady (ω = 0) grid modes directly.
    * Never fix a non-dissipative test by moving a threshold. Read the tableau before quoting its
      damping. ⚠ Gridap's `SDIRK_3_3` is only first order, and `SDIRK_LobattoIIIC_2_2` is the explicit
      Heun scheme (2026-11-11): do not use them.

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
    * **Memory is the only cap on local parallel runs** (author, 2026-11-10). `_slot_throttle.sh`
      admits a case while MemAvailable, minus the growth still owed by cases younger than `YOUNG_S` = 40 min (`RUN_GB` = 3 each),
      stays ≥ `MIN_FREE_GB` = 6. The 2026-11-07 thermal throttle (> 3 busy cores → 1.1 GHz) was a
      machine fault cleared by a reboot; the old `MAXP = 3` rule is withdrawn.

47c. **Cluster output goes to the project storage, never to `$HOME`** (2026-11-23). The home quota
     killed 161 of 275 integrator-study cases on 2026-11-21. `STORAGE` is exported in `~/.bashrc`
     (Snellius: `/gpfs/admin/projectlinks/prjs1960`); `output/` on the cluster is a link to
     `$STORAGE/GridapBALFEM_output` (`run/storage.sh`, one-time `run/setup_storage.sh`).
     `balfem_env.sh` and every submit script refuse otherwise. Set up on Snellius 2026-11-24. Size a campaign's output before
     submitting it (≈ 200 MB per 1-D box case with sub-cell VTK).

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
  * The MCP call aborts after **30 min without output** (the Julia process keeps running): redirect a
    long `include` to a log file (`redirect_stdio(stdout=…)`) so the result survives.
  * Sessions are keyed by the **resolved** project path: a trailing slash or a symlink reaches the same
    (possibly busy) session. Revise is loaded, so `src/` edits apply at the next evaluation.
  * Distributed code paths can be tested without MPI on PartitionedArrays' `with_debug` backend.
* If you add `using X` to `src/`, add it to `[deps]`.
* New design records go in `markdown_files/`; update `STATUS.md` when the state changes; move
  finished plans to `archive/`.
