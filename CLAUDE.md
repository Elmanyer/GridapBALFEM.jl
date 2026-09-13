# CLAUDE.md — `GridapBALFEM.jl/` · the 2D BALFE-M algebraic wave solver

> ## ⇨ START HERE
>
> This file is the map and the standing rules. **It stays at the repository root**, because Claude
> Code auto-loads `CLAUDE.md` from the root and parent directories only — anywhere else and the
> rules silently stop reaching a session started here. The rest of the tracked documentation lives
> in **`markdown_files/`**; start from [`INDEX.md`](markdown_files/INDEX.md), which says which file
> answers which question.
>
> | document | what it answers |
> |---|---|
> | [`INDEX.md`](markdown_files/INDEX.md) | **which file to open** — the taxonomy, and the recurring topics with their one authoritative home |
> | [`MODEL.md`](markdown_files/MODEL.md) | the maths: σ-tensors, the global residual term by term, nonlinear pressure, Jacobians |
> | [`ARCHITECTURE.md`](markdown_files/ARCHITECTURE.md) | code structure and workflow: stacked layout, `src/` map, FE spaces, time loops, distributed path |
> | [`VERIFIED_SCOPE.md`](markdown_files/VERIFIED_SCOPE.md) | **is it correct, and how far does the claim reach?** — the MMS campaign, the Jacobian oracle, scope boundaries |
> | [`TEST_SUITE.md`](markdown_files/TEST_SUITE.md) | **what is checked, and what would slip through?** — every gate and its blind spots |
> | [`OPEN_ISSUES.md`](markdown_files/OPEN_ISSUES.md) | **what is known wrong or unexplained?** — defects in work already done |
> | [`PLANNED_CAMPAIGNS.md`](markdown_files/PLANNED_CAMPAIGNS.md) | **what runs next?** — the MMS campaigns, their ladders, and the code they need |
> | [`COMPLETED_VBASIS_STUDY.md`](markdown_files/COMPLETED_VBASIS_STUDY.md) | the design record of the finished vertical-basis study |
> | [`CONFIGURATION.md`](markdown_files/CONFIGURATION.md) | the settings, the evidence behind each, measured performance |
> | [`WAVE_GENERATION.md`](markdown_files/WAVE_GENERATION.md) | sources, Dirichlet generation, sponge, WaveSpec coupling |
> | [`RUNNING.md`](markdown_files/RUNNING.md) | how to launch: local, cluster, sysimage, env vars |
> | [`HORIZONTAL_CONVERGENCE.md`](markdown_files/HORIZONTAL_CONVERGENCE.md) | **the horizontal pairing × physics convergence matrix** — 1-D complete (18 studies), 2-D partial; why no pairing is deficient in both fields |
> | [`MMS_VBASIS_CAMPAIGN.md`](markdown_files/MMS_VBASIS_CAMPAIGN.md) | results of the `(M,p)` × 8-model convergence matrix |
> | [`CAMPAIGN_COST.md`](markdown_files/CAMPAIGN_COST.md) | measured run costs, memory bands, queue limits, scheduling lessons |
> | [`OUTPUT_NAMING_PROPOSAL.md`](markdown_files/OUTPUT_NAMING_PROPOSAL.md) | the output-directory grammar — implemented in `output_dir_name` |
>
> ⚠ **`markdown_files/` is TRACKED; `latex_docs/` is GITIGNORED.** Anything load-bearing belongs
> here, not in the LaTeX. (Until 2026-09-11 all of this sat in a gitignored `markdown_files/`.)
>
> **✅ THE GRID-SCALE NONLINEAR INSTABILITY IS RESOLVED (2026-09-06). IT WAS THE EQUAL-ORDER
> ELEMENT PAIRING, NOT THE MODEL.** `η` enters the momentum equation *undifferentiated*, through
> `∇·v` after the integration by parts, so it plays exactly the role pressure plays in a Stokes
> system. The runs that blew up used **equal-order `Q2/Q2`** spaces, which are **inf-sup (LBB)
> deficient**; their classic failure mode is a spurious checkerboard at `λ ≈ 2·dx`, and the measured
> growth spectrum peaked at `λ ≈ 2–3·dx`. Switching to a **Taylor-Hood** pairing removes it:
>
> | `dx` | equal-order `Q2/Q2` | Taylor-Hood `Q2/Q1` |
> |---|---|---|
> | 0.50 | suppressed | ✅ 80 s, peak 0.1169 |
> | 0.25 | **died t=26.4, η→5.06** | ✅ 80 s, peak 0.1115 |
> | 0.125 | **died t≈10** | ✅ peak 0.1112 |
>
> **The refinement signature INVERTED** — peak amplitude now *decreases* monotonically with `dx`
> (0.1169 → 0.1115 → 0.1112, converging on the delivered 0.102) where it used to get worse.
> Confirmed on **three integrators** (`SDIRK_2_2`, `SDIRK_3_3`, explicit `EXRK_RungeKutta_4_4`), so
> it is not `dt`-dissipation masking, and on **both** Taylor-Hood pairings. The verified `Q3/Q2`
> reproduces it. `Q2/Q1` at `A=0.10` holds a mean `η` of 0.1058–0.1059 over t=30–80 (fourth decimal,
> 50 wave periods) with Newton at **3.99 iterations/step**, against 30–57 as the equal-order run tore
> itself apart. It sits **+5.2 % above the linear control** — the second-order Stokes crest elevation
> the model exists to produce, a regime previously unreachable because the run died first.
>
> ⚠ **THE FIX COSTS NOTHING.** No residual term changes: `Q2/Q1` and `Q3/Q2` are the pairings the MMS
> campaign already certifies (30/30 spatial, five bases, eight models). Enforced by **rule 2b** — the
> `p_eta = 0` sentinel now resolves to Taylor-Hood everywhere.
>
> ⚠ **THE OLD PHENOMENOLOGY IS VOID, NOT SUPERSEDED.** The `σ ∝ A⁴` growth rates, the "no amplitude
> threshold" finding, the `:none`/`:native`/`:full` tier ordering and the `dt`-masking observations
> were all measured on a discretisation that is no longer used. Do not re-explain them; discard them.
> ⚠ **The three refuted-cure findings remain valid and are worth keeping**: quadrature aliasing is not
> the driver (degrees 6/10/14 agree to 4 dp), the hand Jacobian is not at fault (exact AD agrees to
> 4 dp), and the **advection operator's energy production is exactly the continuity defect** —
> derived, and measured to be five
> orders too small to drive the mode.
> Full account, including the ten refuted hypotheses and the method lessons: **rule 12b** below.
> (`NONLINEAR_INSTABILITY.md` was folded into it and deleted — one record, in the tracked file.)
>
> **One-line status (2026-09-02).** The solver is feature-complete in serial and distributed. The
> analytic MMS verifies **six of the eight models** — all four `:none` and both `:native` — at
> theoretical order, on **five vertical bases**. Suite: sequential 20/20 files, distributed 13/13,
> Jacobian-vs-AD 17/17 over 8 models, nonlinear MMS 8/8, `test/local/` 50/50 — one failure, a
> gate-window specification defect, not a solver defect.
>
> **✅ THE `:full` MMS DEFECT IS FIXED (2026-09-01).** `run_mms_case` never built an `nlp` context,
> so the frozen `{1,2,4,5}` blocks were **absent, not lagged**, and on a flat bed `:full`
> degenerated silently to `:native`. Both MMS drivers now build the context, and both time loops
> **prime it from the initial condition** (`update_nlp_state!` runs only *after* a step, so step 1
> assembled as if from rest — exact for a rest start, wrong for the MMS `u0 = u*(t0)`). Measured on
> P1LFE-2 Model 7: `e_u` fell **1.19e-03 → 1.37e-06 (~870×)** and `p_u` went from a flat −0.00 to
> **1.948**. So the recorded floor measured **omission**, and `:full` is *not* "MMS-unverifiable by
> construction" — it converges, at the projection error. ⚠ Consequence: with the blocks now in the
> residual but still absent from `jacobian_u`, the quasi-Newton gap has a **cliff in amplitude** —
> Newton stalls at 9.2e-04 at the campaign's `a_eta=0.8` and only falls to ~2e-09 at `a_eta ≤ 0.4`.
> Tier-3 studies must drop the manufactured amplitude, and their floors are **not** comparable with
> the previously recorded ones.
>
> **✅ THE VERTICAL-BASIS CAMPAIGN IS COMPLETE (2026-08-30).** 83 studies over five vertical bases
> (P1LFE-2/3/4, P2LFE-1/2; `Nσ` = 3, 4, 5) × the eight model configurations. ⚠ Its Phase-1 σ-meshes
> were derived by multi-property **minimax and are SUPERSEDED TWICE OVER** (rule 46) — `src/vopt.jl` now implements the
> Yang & Liu total-relative-error functional instead (§3), a different objective, so the `p ≥ 2`
> node sets have moved (P2LFE-2: `c₁` 0.8794 → 0.8298). **Headline: `η` reaches 2.999-3.000 and `u` reaches 4.000 on
> EVERY basis and EVERY model in the spatial matrix (30/30) — the order of accuracy is independent
> of the vertical basis, which is the direct evidence for the basis-agnosticism the model is named
> for.** Full results, method and every correction: `markdown_files/MMS_VBASIS_CAMPAIGN.md`; single
> data file `output/local/mms_campaign/campaign_results.csv`.

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
| `Project.toml` / `Manifest.toml` | the Julia package manifest — `name = "GridapBALFEM"`, `uuid = 43e94d05-4d7d-4679-96a4-d46e2615da34`. Loaded with **`using GridapBALFEM`, never `include()`** (§5). This directory is **both the package and the working environment**, so `Test`, `BlockArrays`, `MPIPreferences`, `Preferences` are in `[deps]`, not `[extras]`. `[compat]` admits two Gridap minors **on measured evidence** — see `CONFIGURATION.md` §1 |
| `src/` | the solver package — **17 files** (`vopt.jl` added 2026-09-01), mapped in `ARCHITECTURE.md` §2. `horizontal.jl` carries `check_taylor_hood` (rule 2b); `utilities.jl` carries `output_dir_name`/`unique_output_dir` (rule 2c) |
| `test/` | **30** test files + `runtests.jl` + `test/cluster/` + `test/local/` — inventory and scores in `TEST_SUITE.md`. `test_mms_distributed_parity.jl` (4 ranks) gates the distributed MMS path against the sequential one; `test_taylor_hood.jl` (`:fast`, 13/13) gates the horizontal pairing — Taylor-Hood accepted, equal order REJECTED |
| `examples/` | 7 sequential + `distributed/` (7 cluster scripts + `_dist_common.jl`), `distributed_small/` (5 parametric), `validation/` (7), `local_1d/`, `local_2d/`, `local_mms/` (**11** — the parametric MMS studies, the vertical-basis campaign (`run_vbasis_campaign.jl`, `run_vbasis_shard.jl`, `report_vbasis_campaign.jl`), and the Phase-B batch: `run_phaseB_shard.jl` + `supervise_phaseB2.sh`, a resume-capable runner with a memory-aware supervisor), `inspect_run.jl` — `RUNNING.md` §2 |
| `run/` | 10 production SLURM launchers + `run/dist_small/` (**30** small-domain 2-D cases — the 7 superseded 1-D twins were deleted 2026-09-02) + `run/local/` (**58** scripts: the original 30 incl. **6** 1-D cases + `run_all_1d.sh`, plus the **16** `run_1dnl_*.sh` nonlinear-instability ladder — 4 physics tiers × A∈{0.05,0.10,0.15,0.20} — + `run_all_1dnl.sh`; ⚠ **all 16 of those ran EQUAL ORDER and their results are void — see rule 12b**; plus the **6** `run_skewA_*.sh` (energy-consistent advection, refuted) and the **9** Taylor-Hood campaign scripts `run_th_*.sh` / `run_thref_*.sh` + `run_all_th.sh` that resolved the instability), all through `run/balfem_env.sh` (cluster) or `balfem_local.sh` (workstation) — `RUNNING.md` §3–4. **Both helpers now gate the Taylor-Hood pairing before launch and print it into the banner** (rule 2b). ⚠ **Job sizing is governed by [`run/SNELLIUS_ROME_LAUNCH_CONFIGS.md`](run/SNELLIUS_ROME_LAUNCH_CONFIGS.md)** — TRACKED, because it decides what every job COSTS. A `rome` node is 128 cores / 224 GiB divided into **eighths** (1/8 = 16 cores / 28 GiB), and the bill is `max(cores/128, mem/224)` rounded UP to the next eighth. At the 4 GiB/rank this solver needs, **memory always sets the tier**, so the cost-optimal count is **7k ranks for tier k/8** — 28 ranks bills 4/8 where 32 billed 5/8, and 42 bills 6/8 where 48 billed 7/8. All 38 launchers were resized on 2026-09-07; `--mem-per-cpu=5G` is never valid here |
| `compile/` | the cluster sysimage build chain — `RUNNING.md` §5. ⚠ **`Manifest.toml` is GITIGNORED**, so a fresh cluster checkout has only `Project.toml`; `set_preferences.jl` therefore `Pkg.add`s both forks by URL (`Elmanyer/Gridap.jl` @ `fix-transient-multifield-ad` for transient-multifield AD, `Elmanyer/WaveSpec.jl` for the `change_seed!` fix) — `Pkg.instantiate()` would resolve the BROKEN upstream versions. Consequence: builds track the branch tips and are not reproducible across time |
| `output/` | ⚠ **GITIGNORED** — nothing here survives a fresh clone, so anything load-bearing must be copied into a tracked file. Result records worth knowing about: `output/local/vopt/` (the corrected mesh-optimisation results + `README.md`, rule 46), `output/local/mms_vbasis/` (vertical-basis convergence, **10/10 PASS** on the published `p=1` meshes — accepted, not to be re-run), `output/local/mms_campaign/` (the 2026-08-30 five-basis campaign, 30/30 spatial), `output/local/mms_phaseB/` (live batch + `INVALID_P2LFE-2.md`) |
| `postprocessing/` | `GridapBALFEMPost` — self-contained, own environment, **no dependency on the solver** |
| `WaveSpec.jl/` | vendored stochastic sea-state synthesis (CMOE-TUDelft). Tracks the **GitHub repository version, not a tagged release** — the release's `change_seed!` is broken |
| `Gridap.jl/` | the vendored **fork** (`Elmanyer/Gridap.jl` @ `fix-transient-multifield-ad`, one commit on `v0.20.8`) making transient-multifield AD work — `CONFIGURATION.md` §2 |
| `latex_docs/BALFEM_models/` | **THE CURRENT LaTeX project** — the authoritative derivation (§2). ⚠ the `.zip` no longer exists in the checkout; re-export from Overleaf before assuming a round-trip |
| `latex_docs/LFEM_discretisation.zip` | **SUPERSEDED** pre-rename LaTeX, kept for provenance as a zip only. Do not edit |
| `latex_docs/CFC2027_abstract/` | conference abstract (CFC 2027); its class needs the `newtx` fonts. **Renamed from `CFC2027_LFEMultilayer_abstract/`** |
| `latex_docs/doc_figures/` | `generate_doc_plots.ipynb` — one cell per figure, activating the repo project by walking up to `Project.toml`. Reproduces Yang & Liu figs 2 and 3 from `assemble_dispersion_tensors` + `model_R`; recovers their published `kd_app` (10.84 / 39.23 / 127.91) |
| `markdown_files/` | the **eleven** documents above, plus the LaTeX projects, `doc_figures/` and the abstract. **Gitignored** |
| `GridapSWE.jl/` | vendored reference implementation, **untracked** — not a dependency of this package |
| `../LFE-M_2D_solver/` | a genuinely external per-layer implementation of the same weak form. Legacy, deliberately not renamed |

---

## 2. The LaTeX project (`BALFEM_models/`)

*Derivation, Study and Implementation of Basis-Agnostic Layer-integrated Finite Element (BALFE-$M$)
Models.* `main.tex` inputs `SigmaEulerModel.tex`; the chapter *Vertical Multilayer Discretisation*
over five sections in `VerticalFESemiDiscretisation/` (`VerticalFEapprox` — carrying the nomenclature
paragraph — `wDerivation`, `pDerivation`, `VerticalProjection`, `VerticalSemiDiscreteSystem` incl.
the flat-bed reduction); then `LinearModel.tex`, `StokesWaveFourierAnalysis.tex` (the analytical
core: the Stokes–Fourier hierarchy, the dispersion functional `R(μ) = Φᵀ(M+μ|B|)⁻¹Φ` and its four
basis-independent properties, group velocity and shoaling gradient, the second-order bound harmonic
reproducing the published `kd=6.0` for P1LFE-2, third-order solvability, and the vertical grid
optimisation `Δσ_top ≈ 2.94/kd_max`); then `NumericalImplementation/` — `GlobalResidual`,
`GridapImplementation.tex` and `ValidationTests.tex`.

> **`GridapImplementation.tex` §`subsec: term classification` IS THE RESIDUAL'S SPECIFICATION.**
> Every term of the full model tagged by amplitude order × bed-slope class × activation condition on
> each of the three switches. The table **factorises** — `regime` is a truncation in amplitude
> order, `flat_bed` a projection onto `∇h ≡ 0`, `nl_pressure` a component filter on `𝓝` — which is
> the formal statement that the three switches are orthogonal. The `O(ε)` rows are exactly eight,
> and they *are* the linearised model.

**Working rules:**

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

**Solver core.** The stacked `[η,𝖴x,𝖴y]` loop-free residual + hand Jacobians; time loops defaulting
to the fully-implicit `RungeKutta(:SDIRK_2_2)` (`:theta` Crank–Nicolson selectable) — sequential
LU+Newton, distributed GMRES+Jacobi+Newton; the full nonlinear physics (advection, the full leading
pressure `R_P`, all eight `𝓝` components) in **both** regimes. Wavemaker / sponge / wall / periodic
BCs; runtime monitoring plus an independent governing-equation residual checker; `w_s`/`p_s` VTK
reconstruction.

**Dirichlet boundary wave generation + WaveSpec coupling.** Regular, multichromatic, or WaveSpec
`AiryState` stochastic sea states (seeded phases ⇒ rank-deterministic). The `:model`
discrete-eigenmode polarization prescribes an exact discrete transport, so a generated wave is a
solution of the discrete equations at the boundary and radiates cleanly.

**Verification.** The analytic MMS (`src/mms.jl`), whose independence from `problem.jl` is enforced
by a grep gate; `test_jacobians_ad.jl`, comparing the hand Jacobians against AD of the same residual
matrix-by-matrix and gating the nonlinear branch on *amplitude scaling*; the linear
one-Newton-iteration gate on a sloping bed; and the **vertical-basis campaign** (§5), which
reproduces the whole verified scope on five vertical bases rather than one.

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

---

## 4. Physics selection — three orthogonal controls

| control | values | meaning |
|---|---|---|
| `regime` | `:linear` \| `:nonlinear` | linearised core, no advection / full nonlinear core + advection |
| `nl_pressure` | `:none` \| `:native` \| `:full` | 𝓝 off / components `{3,6,7,8}` / `+{1,2,4,5}` |
| `flat_bed` | `Bool` | `true` ⇔ **`∇h ≡ 0`** (every ∇h-term dropped, ∇η/dispersion kept); `false` = variable bathymetry |

`resolve_physics` maps these onto the seven internal booleans; `build_problem_raw` is the low-level
escape hatch. Pressure content is intrinsic to the model: `P_full = advection`,
`lin_pressure = advection ∨ ¬flat_bed`. `flat_bed` acts at a **single control point** —
`dhx,dhy = flat_bed ? 0 : ∂h` in `global_residual`/`jacobian_*`. Full term table: `MODEL.md` §6.

**`nl_pressure=:native` is the production tier.** The whole `{1,2,4,5}` hierarchy contributes
**0.013 % (1-D) / 0.094 % (2-D)** on top of advection's 0.77 % / 1.84 % at `A=1e-3`, and `:full`
carries a mesh-independent velocity-error floor (`VERIFIED_SCOPE.md` §4).

---

## 5. Current Implementation Stage

*Rewritten 2026-09-13, after Campaign C, the 2-D Q4/Q3 tier and the quadrature probe.*

### 5.0 Status in one paragraph

The solver is **feature-complete in serial and distributed form** and has been so since 2026-09-02:
stacked loop-free residual, hand Jacobians, the full nonlinear physics, SDIRK/θ integrators, every
boundary treatment, Dirichlet generation with WaveSpec coupling. **The LINEAR models are in excellent
shape** — optimal order in both fields, on five vertical bases, in 1-D and 2-D, sequential and
distributed. **The NONLINEAR models are correct but carry one unexplained order reduction** in the
surface elevation at fine mesh, isolated as of 2026-09-13 to the advection block and *not* explained
by any of the six candidate causes tested. Production wave runs are stable at realistic amplitude
since the Taylor-Hood fix, with the exception of boundary-generation failures that are a separate,
configuration-level problem.

### 5.1 LINEAR models (1, 2) — verified, no open defects

| property | evidence |
|---|---|
| `p_η` optimal | **3.000** at Q3/Q2, **4.000** at Q4/Q3, **2.000** at Q2/Q1 — all six 1-D studies |
| `p_u` optimal | **3.99** at Q3/Q2 (rising 3.50→3.86→3.96→3.99) |
| 2-D agreement | `p_η` 3.0000 (Q3/Q2), 3.9990 (Q4/Q3) — matches 1-D to 3–4 digits |
| transient ≡ static | agree to 3–4 digits, so the measured rates are clean *spatial* rates |
| seq ≡ distributed | agree to 3–4 digits (the parity `test_mms_distributed_parity` gates) |
| one Newton iteration | `test_linear_newton_gate` 10/10 — exactly 1 iteration/stage on a **sloping** bed, residual 3.9e-15 |
| five vertical bases | `MMS_VBASIS_CAMPAIGN.md`: 30/30 spatial at optimal order, `Nσ` = 3, 4, 5 |

⚠ **The one blemish is shared with the nonlinear models and is not a linear defect:** at **Q2/Q1**
`p_u` converges to **2.003** against an optimal 3 — a full order short, on every model, linear
included. See §5.3.

### 5.2 NONLINEAR models (3–6) — correct, stable, one order-reduction defect

**What is established.**

| property | evidence |
|---|---|
| residual correctness | analytic MMS verifies models 3–6 (`:none` and `:native`) at theoretical order on the measured ladder |
| Jacobians | `test_jacobians_ad` 17/17 over all 8 models; `∂R/∂u̇` exact, nonlinear `∂R/∂u` quasi-Newton by choice with the gap vanishing at order 1.11–1.16 in amplitude |
| `:native` is free | `:none` vs `:native` agree to **~1e-4 in `p_η`** at every pairing, 1-D and 2-D — the `{1,2,4,5}` hierarchy is dynamically negligible |
| bathymetry is free | flat vs variable bed agree to ~4e-4 in rate |
| long-run stability | 50 wave periods at `A = 0.10` m (`κa ≈ 0.16`, `kd = 5.5`), mean `η` flat to the 4th decimal over t = 30–80, Newton **3.99 it/step**, **+5.2 %** above the linear control — the 2nd-order Stokes crest |
| Q4/Q3 clean | `p_η` **4.000–4.002** in 1-D and **3.9997–4.0004** in 2-D, all six models |

**⛔ THE OPEN DEFECT: `p_η` loses an order at Q3/Q2 on fine meshes.**

```
nonlinear  pairwise p_eta   2.971 -> 2.922 -> 2.757 -> 2.450     (nx = 4..64, optimal 3)
linear     same ladder      2.990 -> 2.997 -> 2.999 -> 3.000
```

The error still falls; the *rate* decays, monotonically, **worsening with refinement**. Reproduced
in 2-D (2.5565). `p_u` on the same runs moves the opposite way (2.77 → 3.69, rising). Six causes
eliminated by measurement — algebra, `𝓝`, `∇h`, the linear core, 1-D posing, quadrature — leaving
the advection block. **Full account and next step: `OPEN_ISSUES.md` §0b.**

### 5.3 The Q2/Q1 velocity shortfall — universal, cause unattributed

`p_u = 2.003–2.012` against an optimal 3 on **all six models**, while `p_η` is exactly 2.000. The
sequence *descends* onto 2.00 and `e_u ≈ 1.6–3.3e-05` sits five orders above the algebraic floor, so
it is a converged rate, not saturation. Identical in 2-D (**1.9998**, and to four digits for linear
and nonlinear alike).

⚠ **It is not established that this is a defect.** The "optimum" is the Taylor-Hood L² pattern
(`u → p_u+1`, `η → p_u`) inherited by analogy because `η` plays the pressure role; **those estimates
have never been derived for this dispersive depth-integrated system.** `HORIZONTAL_CONVERGENCE.md`
§4 sets out why the matrix as a whole does not localise a failure — every pairing is optimal in at
least one field, none is deficient in both, and the deficiency *moves* between fields as `p_u`
increments. **Practical consequence: use Q3/Q2, the cheapest pairing optimal in both fields.**

### 5.4 Verification campaigns — everything run, with results

| campaign | scope | result |
|---|---|---|
| **Vertical-basis** (2026-08-30) | 83 studies, 5 bases × 8 models | **30/30 spatial** at optimal order; order independent of the vertical basis — the evidence for basis-agnosticism |
| **Phase B** (to 2026-09-11) | 38 studies (T7–T10) | extended ladders, `:full` floor, `p≥2` bases, the horizontal pairings |
| **Campaign C** (2026-09-12) | **18 studies / 84 runs, 1-D**, 3 pairings × 6 models | **18 OK, 0 errors, 44.1 core-h** — full matrix in `HORIZONTAL_CONVERGENCE.md` §2 |
| **C3, 2-D Q4/Q3** (2026-09-12) | 6 studies / 24 runs | **6/6 OK, 77.6 core-h**, all at optimal `p_η`; closed the largest gap in the pairing study |
| **Quadrature probe** (2026-09-13) | dose-response `q = 0,2,4` × 2 beds, + residual assembly | **⛔ NOT quadrature** — knob live (1.4e-08) but ~7 orders too small |

### 5.5 Test suite — last full run 2026-08-18/19, re-baselined on Taylor-Hood 2026-09-07

**Sequential 20/20 files · distributed 13/13 gates on 4 ranks · `test_jacobians_ad` 17/17 over 8
models · nonlinear MMS 8/8 · `test/local/` 50/50.** Per-file scores: `TEST_SUITE.md` §2. Notable
gates: `test_taylor_hood` 13/13 (rejects equal order), `test_linear_newton_gate` 10/10,
`test_conservation` mass drift **9.8e-16**, `test_energy` **6.0e-14**.

**Three known failures, all understood, none a solver defect:** `test_mms_convergence` G7 (a
gate-window specification defect, rules 32/33); `test_nlpressure` G1/G3 (test-side — G1 compared an
analytic value against an FE-interpolated one, fixed by raising that test to `Q3/Q2`; G3 is a stale
constant needing re-measurement); `test_bc_generation` case C (10.1 % against a 10 % gate).

⚠ **The suite could not distinguish equal-order from Taylor-Hood at all** — all 17 physics/smoke
tests passed unchanged across that discretisation change, because they run at tiny amplitude over
short durations. **That blind spot is why the instability reached production, and it is still open.**

### 5.6 Cluster production runs (Snellius, 2026-09, `nl_pressure=:full`)

| case | outcome |
|---|---|
| irregular sea, flat | ✅ completed 2600/2600 steps to t=52 |
| irregular sea, bar | ⛔ NaN at t=31.96 |
| directional sea, flat / bar | ⛔ NaN at t=5.52 / 6.26 |
| ring wave | ⛔ NaN at t=6.76 |

**Diagnosed, all four:** a **transverse, grid-scale, velocity-led** mode — η stayed bounded at
0.13–0.14 m throughout while `u` exploded, and the y-symmetry error of a y-invariant problem grew
1e-3 → 0.887 with an e-folding of ~0.5 s. **Its seed was a real bug, now fixed**
(`build_airy_state` seeded the phases but not the angular spreading, so every MPI rank generated a
different sea; the directional inflow carried six uncorrelated seas with O(1) jumps at the rank
cuts). Two further causes were separate: the directional runs also had the Dirichlet inflow
overlapping 80 % of the lateral sponge, and the ring exceeded the Miche breaking limit by 2.2×.
⚠ **The flat irregular run did not "pass" — it ended at T_final while its own asymmetry was at 0.114
and climbing.** All launchers are re-specified and none has been re-run.

### 5.7 Still to check — ordered by what would change a conclusion

1. ⛔ **The nonlinear `p_η` order reduction** (`OPEN_ISSUES.md` §0b). Next step is **derivation, not
   another run**: check the advection block term-by-term against `BALFEM_models/`. Cheap
   discriminators if wanted: amplitude sweep (`a_eta = 0.8 → 0.4 → 0.2`), and Q4/Q3 extended to
   `nx = 64`.
2. 🔴 **Derive the expected orders for this system.** Until that exists, "sub-optimal" in every
   convergence table is measured against an assumed optimum, and §5.3 cannot be resolved by runs.
3. 🔴 **Re-run the cluster suite** with the seed fix, the new directional geometry and the corrected
   sponges. Nothing on record post-dates those fixes.
4. 🔴 **A long-duration nonlinear regression test.** The mode needed 50–80 s to emerge and no test
   runs that long; this is the gap that let it reach production.
5. 🟠 `:full` (models 7–8) is **not MMS-verifiable as built** — the `{1,2,4,5}` blocks are in the
   residual but absent from `jacobian_u`. Closing it means completing `jacobian_u`, not refining a
   mesh. **Every diverged cluster run above was `:full`, and the batch had no `:native` control.**
6. 🟠 **The `:sdirk` linear-model temporal deficit** — `pw_u ≈ 1.69` against `:theta`'s 1.99, not
   spatial-floor contamination, unexplained.
7. 🟠 Per-basis `dt` ladders for temporal studies; no MPI tier in `runtests.jl`; preconditioner
   replacement (the largest performance item); `Nσ=5` nonlinear temporal is a solve limit.
8. 🟠 **Cluster checkout and sysimage rebuild** — `src/mms_driver.jl` changed 2026-09-12, so the
   content stamp is stale (functionally irrelevant: no cluster driver touches the MMS path).

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

These are the rules that cost something to learn. Each is stated where it is enforced; the
supporting measurement is in the linked document.

### Model and residual

1. **`R_P` is the entire frequency dispersion of the model.** Only the *boundary* part of its
   integration by parts vanishes; the volume part must be assembled. Without it the model degenerates
   to non-dispersive shallow water.
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
    station. Confirmed on `Q3/Q2` as well, and on **three integrators** — `SDIRK_2_2` (L-stable),
    `SDIRK_3_3`, and explicit `EXRK_RungeKutta_4_4` (near dissipation-free, which is what retires the
    `dt`-masking objection; RK4 tracked SDIRK to within 1–5 % with the same envelope, the offset
    being exactly the L-stable damping).

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

12c. **"STABLE" IS MEANINGLESS WITHOUT `dx`, `dt` AND DURATION.** An L-stable integrator's numerical
    dissipation grows with `dt` and can *mask* a growing mode, so **a run that completes may simply
    be one whose dissipation exceeded the growth rate.** This is rule 15's trap in a new guise, and
    the rule stands on its own merits — but note how it was DISCHARGED for 12b, because that is the
    template: the same case was run on `SDIRK_2_2` (L-stable, 2nd order), `SDIRK_3_3` (3rd order) and
    **explicit `EXRK_RungeKutta_4_4`** (essentially non-dissipative). All three stayed flat, and RK4
    tracked SDIRK_2_2 to within 1–5 % with the same envelope shape — the offset being exactly the
    L-stable damping. **Three integrators of different order and stability character agreeing is what
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

15. **The default integrator (`SDIRK_2_2`) is L-stable, i.e. DISSIPATIVE BY CONSTRUCTION.** Any test
    measuring a non-dissipative property must pin `solver_type=:theta`. **Do not remove those pins,
    and never fix such a failure by moving a threshold.** Recognise it by: amplitude damped while
    **phase is correct**, error **growing with frequency**, and **refining the mesh does not help** —
    that last check separates it from genuine under-resolution. (`TEST_SUITE.md` §4)
16. **Distributed linear solve = `NewtonSolver(GMRESSolver(Pr=Jacobi))`.** A direct LU does not scale
    to partitioned matrices at cluster size.
17. **`krylov_m` (basis size, memory) and `ls_maxiter` (iteration budget, time) are different
    bounds**, and `restart=true` is load-bearing. Symptom of getting this wrong: **`gmres=` pinned at
    exactly the same number every step.** (`ARCHITECTURE.md` §5)
17b. **A CONVERGED NEWTON STEP IS A PROPERTY OF THE RESIDUAL, NOT THE JACOBIAN.** The Jacobian sets
    the *path* to the root and the *cost*; the root is defined by `r(u)=0`. So an incomplete but
    convergent Jacobian **cannot change the answer** — only the iteration count, or whether it
    converges at all. Measured: exact-AD and quasi-Newton hand Jacobians agree to **four decimal
    places** through an entire divergence (0.8100 vs 0.8102), while costing 6 vs 49 iterations.
    **Corollary: never diagnose a wrong ANSWER by replacing the Jacobian** — and do not read Newton
    stalling as the primary fault when the solution underneath is already diverging.

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

39. **Test a diagnosis against a case it cannot explain, rather than looking harder where it points.**
40. **Diagnostics interpretation:** never read `growth` without `x_at_max`; `dmp/int` means "≫1 is
    trouble", never "<1 is fine"; mass drift is an invariant **only** in a closed, unforced basin.
    (`ARCHITECTURE.md` §6)

---

## 8. Conventions for editing here

* **`BALFEM_models/` is the single source of mathematical truth.** `MODEL.md` mirrors its notation.
  If you change the maths, change both.
* **Notation bridge, the most common source of confusion.** LaTeX `M^V, 𝓜^V, 𝓖^V, A^V, K^V, 𝓐^V,
  𝓚^V, Φ, φ_j` ↔ code `Mmat, Mcal, Gcal, A, K, Acal, Kcal, Phi/D/C`, with the unit basis stored as
  `w_j = −φ_j_int`. Always cross-reference `MODEL.md` §2 when quoting tensor names.
* **If you add a `using X` to any `src/*.jl`, add it to `[deps]`.**
* **`markdown_files/` is gitignored** — nothing in it is version-controlled, so anything load-bearing
  belongs in the tracked `CLAUDE.md`, `README.md` or test files instead.
* **Run Julia via the `julia-mcp` tool.**
