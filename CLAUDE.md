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
> | [`HORIZONTAL_CONVERGENCE_1D.md`](markdown_files/HORIZONTAL_CONVERGENCE_1D.md) | **the 1-D horizontal pairing × physics convergence matrix** — 18 studies, and why no pairing is deficient in both fields |
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

**Working.** Feature-complete in both serial and distributed forms: the stacked loop-free residual +
hand Jacobians, the full nonlinear physics, the SDIRK/θ integrators, all boundary treatments, and
Dirichlet boundary wave generation with WaveSpec coupling.

**✅ NEW (2026-09-06) — long fully-nonlinear runs are now stable, and that is the headline change.**
The grid-scale instability that made every long nonlinear run unusable was the **equal-order element
pairing** (rules 2b, 12b). On Taylor-Hood the 1-D flume runs **50 wave periods at `A=0.10`
(`κa ≈ 0.16`, `kd = 5.5`)** with a mean `η` flat to the fourth decimal over t=30–80 and Newton at
3.99 iterations/step — and sitting **+5.2 % above the linear control**, i.e. reproducing the
second-order Stokes crest elevation. This is the regime the model was built for and it was
unreachable until now. It cost **no residual change and no re-verification**.

**✅ THE SUITE IS RE-BASELINED ON TAYLOR-HOOD (2026-09-07).** The whole sequential suite plus the
local tier was re-run after the pairing change. ⚠ **The predicted mass re-baselining did NOT
happen** — of the 17 physics/smoke tests that silently moved `Q2/Q2 → Q2/Q1`, **every one passed
without touching a threshold.**

| result | |
|---|---|
| 19 of 21 completed sequential files | ✅ pass |
| `test_linear_newton_gate` | 10/10 — exactly 1 Newton iter/stage, residual 3.9e-15 |
| `test_conservation` | mass drift **9.8e-16** |
| `test_energy` | energy drift **6.0e-14**, envelope −0.01 % |
| `test_taylor_hood` (new) | 13/13 |
| local tier: `reststate` / `sponge` / `relaxation` | 8/0, 18/0, 9/0 |
| MMS spatial, linear (`G6`) | `p_η` 2.995, `p_u` 3.770 → 4 |
| MMS nonlinear, Models 3 and 4 | `p_η` 2.996, `p_u` 3.995 / 3.997 |

⚠ **Why the thresholds absorbed a discretisation change is itself the finding**: these tests run at
tiny amplitude over short durations, the one regime where the equal-order deficiency never
expressed itself. **The suite could not distinguish the two discretisations at all** — the same
blind spot that let the instability reach production. The coverage gap is real and still open even
though the re-baselining was not needed.

**Three failures, all understood, none a solver defect:**
* `test_mms_convergence` **G7** — the *pre-existing, documented* gate-window defect. The pairwise
  rates read 2.008 then 0.757: the temporal ladder's finest `dt` has sunk to the spatial floor, and
  `G10` confirms it (a 1.5× mesh refinement moves `e_eta` by 0.5 %). Rules 32/33 — the fitted 1.382
  averages a good rate with a saturated one. `G6`, `G8`, `G9`, `G10` all pass.
* `test_nlpressure` **G1/G3** — test-side. G1 compares an analytic evaluation against an
  FE-*interpolated* one, which is an identity only if the space represents the manufactured field
  exactly; `eta_f` is quadratic, so under `Q2/Q1` the Q1 surface cannot hold it and the 2.05e-03 gap
  **was the interpolation error, not a broken identity**. Fixed by raising that test to `Q3/Q2`.
  G3 is a stale constant (0.00286 measured at equal order vs 0.00316) — **re-measure, not
  re-threshold**.
* `test_bc_generation` **case C** — delivered amplitude 10.1 % against a 10 % gate, linear regime.
  Threshold untouched pending a refinement check; `test_relaxation_1d` passing 9/0 argues it is
  marginal rather than structural.

**Verified scope — six of eight models** (`Q3/Q2`, 1-D static unless noted). ⚠ **The table below is
the ORIGINAL single-basis (P1LFE-2) campaign.** Every one of these rows has since been reproduced on
four further vertical bases — see `MMS_VBASIS_CAMPAIGN.md` — so "verified" now means verified across
`Nσ` = 3, 4 and 5, not at one basis:

| model | `regime` / `flat_bed` / `nl_pressure` | `p_η` (opt 3) | `p_u` (opt 4) | |
|---|---|---|---|---|
| 1 | `:linear` / flat / `:none` | `p_e+1` exactly | `p_u+1` exactly | ✅ |
| 2 | `:linear` / **variable** / `:none` | 3.000 | 4.000 | ✅ |
| 3 | `:nonlinear` / flat / `:none` | 2.996 | 3.995 | ✅ |
| 4 | `:nonlinear` / **variable** / `:none` | 2.996 | 3.997 | ✅ |
| 5 | `:nonlinear` / flat / **`:native`** | 2.996 | 3.997 | ✅ |
| 6 | `:nonlinear` / **variable** / **`:native`** | 2.996 | 3.998 | ✅ |
| 7–8 | `:nonlinear` / any / **`:full`** | 2.99 → 2.59 | **−0.00** | ⛔ not MMS-verifiable **by construction** |

Model 2 additionally confirmed transient (2.999/3.998) and 2-D (3.000/3.963).

> **Say "the `:none` and `:native` models are verified."** Never bare *"the residual is verified"*
> (which would wrongly include `:full`), and never *"the `𝓝` tiers are verified"* (same error).

**Vertical-basis campaign — COMPLETE 2026-08-30** (`MMS_VBASIS_CAMPAIGN.md`; data in
`output/local/mms_campaign/campaign_results.csv`). 83 studies over five bases × the eight models.

*Phase 1 — optimised σ-meshes*, by the multi-property minimax of `StokesWaveFourierAnalysis.tex`
§sec: vertical grid optimisation (inner Chebyshev problem + outer bisection, over `C`, `C_g`, `γ`).
Calibrated on two independent standards: Table 4.1's nine ranges (<1 %) and the published
band-dependent optima `c₁ = 0.702/0.802/0.860` at `K = 5/10/20` (to 2e-4).

| basis | `Nσ` | optimised `c_bdy` | `kd_app` (multi-property; **γ binds throughout**) |
|---|---|---|---|
| P1LFE-2 | 3 | `[0, 0.8064, 1]` | 8.60 |
| P2LFE-1 | 3 | `[0, 1]` (no free parameter) | 3.01 |
| P1LFE-3 | 4 | `[0, 0.7597, 0.9339, 1]` | 24.95 |
| P1LFE-4 | 5 | `[0, 0.7809, 0.9335, 0.9820, 1]` | 92.22 |
| P2LFE-2 | 5 | `[0, 0.8794, 1]` | 21.84 |

*Phase 2 — the convergence matrix.* **SPATIAL: 30/30, `η` → 2.999-3.000 and `u` → 4.000 on EVERY
basis and EVERY model.** The order of accuracy is **independent of the vertical basis** — the direct
quantitative evidence for the property the model family is named for, previously resting on P1LFE-2
alone. TEMPORAL: clean second order at `Nσ` = 3 and 4; `Nσ = 5` nonlinear is a scope limit (below).
Tier 3 (`:full` floor) 10/10.

**Three results worth quoting, and one scope limit:**
1. **At fixed `Nσ`, grading beats raising the order** — 8.60 vs 3.01 at `Nσ=3`, 92.22 vs 21.84 at
   `Nσ=5` (2.9× and 4.2×). Independently confirms `main.tex`'s statement, with the margin *widening*.
2. **P1LFE-4's apparent `u`-shortfall is PRE-ASYMPTOTIC**, not suboptimal: `pw_u` runs
   3.33 → 3.45 → 3.78 → **3.94** to `nx=128` on three models. Answers `OPEN_ISSUES.md` §6 for this case.
3. **The `:full` floor GROWS with `Nσ`** (1.2 → 2.0 → 2.9e-03 in the `p=1` family) **and depends on
   basis shape at high `Nσ`** (1.6× gap at `Nσ=5`). ⚠ It is an **omission** floor, not a
   frozen-projection one — see the `:full` open item.
4. ⛔ **`Nσ = 5` nonlinear TEMPORAL is not measurable on this machine.** The degradation is monotone
   in model complexity: linear ✅ → nonlinear/flat ⚠ (rate falls to a floor) → nonlinear/∇h ⛔
   (saturates or never converges) → `:native` ⛔ (NaN). **That ordering is the evidence it is a
   nonlinear-SOLVE limit, not an operator defect** — the linear models on the identical basis and
   mesh stay textbook.

**Suite** (all measured 2026-08-18/19, not carried over): sequential **20/20 files**, distributed
**13/13 gates** on 4 ranks, `test_jacobians_ad` **17/17** over 8 models, `test_mms_convergence_nonlinear`
**8/8**, `test/local/` **50/50**. Per-file scores: `TEST_SUITE.md` §2.

**Production runs — the small-domain suite is READY TO LAUNCH** (audited 2026-08-19). The 5
parametric scripts and 20 two-dimensional launchers were checked mechanically against the driver's
keyword set and against each case's own geometry. Fixed in that pass: transit-aware durations (two
case families would otherwise have ended while the domain was still filling), the inflow relaxation
zone (documented as on, coded as off, on all three `:bc_gen` scripts, which also run
`sponge_wL=0`), sponge widths sized against wavelength rather than domain fraction (ring
0.64 → 1.12 λ; directional lateral 0.45 → 0.90 λ_eff; long-period case `Lx` 50 → 70 m), and a guard
refusing oblique generation on a y-periodic domain. Boundary conditions verified correct on all
five cases.

**1-D cases run LOCALLY AND SEQUENTIALLY — this is a decision, not a fallback** (2026-08-19). The
solver is structurally 2-D, so a "1-D" problem is a narrow y-periodic strip whose cross-section is
pinned at Gridap's periodic minimum (`ny=3`). Direct LU cost scales with the *front width*, which
that cross-section fixes, so **the DOF count cancels from the serial-vs-distributed ratio** and no
domain length or `dx` puts a genuine 1-D case on the distributed side. Refining `ny` would move it
there, but the solution is exactly y-invariant, so that buys parallel efficiency with work that
produces no physics. Consequence: `run/local/run_1d_*.sh` (sequential, with point gauges) are the
supported path. The seven cluster twins that used to sit in `run/dist_small/` were **deleted
2026-09-02**: they carried no warning in their own headers, so the only record that they were dead
was this paragraph — a file that looks launchable, in the launcher directory, next to 30 live ones.
Their physics is preserved case-for-case in `run/local/`; git history holds their cluster sizing
(200 m / `nx=800` / 60 periods at 32 ranks) if it is ever wanted.

**✅ CAMPAIGN C IS COMPLETE (2026-09-12) — 18/18 1-D studies, 84/84 solver runs, 0 errors, 44.1
core-hours** (`PLANNED_CAMPAIGNS.md` §2; data `output/local/mms/phaseB/shard_NN_v2.csv`). Two
results settle long-standing questions and one opens a defect:

* **The Q2/Q1 one-order velocity shortfall is UNIVERSAL** — `p_u = 2.003–2.012` against an optimal 3
  on all six models, with `p_η` exactly optimal. Tier 2 had it on flat-bed `:none` only, which by
  rule 4 could not exercise `∇h`; it now holds on **variable bathymetry** and on **`:native`, the
  production tier**. **Q3/Q2 is the cheapest pairing optimal in both fields.**
* **`:native` is free** — `:none` vs `:native` differ by ~1e-4 in `p_η` on every matched pair.
* ⛔ **NONLINEAR `p_η` LOSES AN ORDER on fine meshes** — 2.971 → 2.922 → 2.757 → **2.450** at Q3/Q2
  while linear models on the identical ladder hold 3.000. Not algebraic, not `𝓝`, not `∇h`, not
  quadrature (each eliminated by measurement), so it sits in the **advection block**. It is a
  fine-mesh effect beyond `nx=32`, which is why the 30/30 verified spatial studies never saw it.
  Full account and next step: **`OPEN_ISSUES.md` §0b**.

**Open work** — nothing is half-built in the solver; the open items are verification gaps,
performance, and follow-through. Full list with decisive next steps: `OPEN_ISSUES.md`; studies
designed but not begun: `COMPLETED_VBASIS_STUDY.md`.

**▶ PHASE-B MMS BATCH RUNNING — 35 valid rows as of 2026-09-09** (`output/local/mms_phaseB/`, 12
shard slots, `run_phaseB_shard.jl` + `supervise_phaseB2.sh`). A shard is a **worker slot, not a task
type**: one 34-job queue is cost-sorted (SJF) and dealt round-robin, shard `s` owning indices with
`(i−1) mod 12 == s`, each an independent process (rule 42). `PHASEB_MAX_STUDIES=1` recycles a worker
after every study, which is why shard indices rotate — that is rule 41 being *enforced*, not drift.
Tasks: T7 P1LFE-4 extended ladder, T8 the `:full` pair at `a_eta=0.4`, T9/T9b tier 2 (the horizontal
pairings; T9b re-specifies `Q4/Q3` on a **shorter** ladder because at `nx=128` its errors sit on the
double-precision floor), T10 the `p ≥ 2` bases.
⚠ **P2LFE-2 IS CANCELLED from T8 and T10 (2026-09-08) and is NOT to be regenerated.** Its mesh came
from `DEFAULT_CBDY_P[(2,2)] = [0, 0.8298, 1]`, produced by the superseded objective (rule 46), so
those studies refine on a mesh that no longer exists. The **6 rows already written are void** —
`output/local/mms_phaseB/INVALID_P2LFE-2.md` carries the filter. Everything else is unaffected:
`p=1` resolves to the published `DEFAULT_CBDY`, and P2LFE-1 is `M=1` with no free interface.
⚠ Timings are legitimately long — a single mesh level is a 100-step run at 120–150 s/step, so one
CSV row can cost most of a day. The `solve time 00:00 (0.0%)` and `Newton iters: 0` in those logs are
**instrumentation gaps, not stalls**: `mms_driver.jl` installs no `SolverMonitor`, so those counters
are never incremented. Wall time is real.
**Provisional, and three results already matter:**
* **`Q2/Q1` has a genuine one-order velocity shortfall** — `p_u = 2.001` against an optimal 3,
  reproduced on two models over a five-level ladder to `nx=128`, with `e_u ≈ 8e-06`, five orders
  above round-off. `p_η` is optimal on the same runs. `Q3/Q2` is optimal (3.000/3.998).
* ⚠ **`Q4/Q3` in 1-D is UNMEASURABLE, not defective.** Its `p_u = 2.505` sits on `e_u = 1.4e-10`,
  the double-precision floor. A saturated error and a genuine low rate give the same slope — it
  needs a *shorter* ladder (`nx = 8…32`). Do not quote that rate.
* **Nonlinear `p_η` degrades to ≈2.17–2.27 on fine meshes** where linear models on the same basis
  and ladder give 2.9998. It appears both on `p=2` bases at four levels and on `p=1` at five, so it
  is a **fine-mesh** effect, not a property of the vertical order — a second-order component takes
  over once the third-order part has decayed. It is **not** caused by the new σ-meshes: P2LFE-1 has
  no free interface at all and degrades identically.

* ✅ **THE NONLINEAR GRID-SCALE INSTABILITY — RESOLVED 2026-09-06.** It was the **equal-order
  `Q2/Q2`** pairing (inf-sup deficient; spurious checkerboard at `λ ≈ 2·dx`, matching the measured
  `λ ≈ 2–3·dx` growth peak), not the model, the residual or the integrator. Taylor-Hood removes it
  and the **refinement signature inverts**. Costs nothing — `Q2/Q1` and `Q3/Q2` are already MMS-
  verified. Enforced by rule 2b; the complete account is rule 12b.
  * ✅ **DONE 2026-09-07 — the 17 physics/smoke tests were re-run on Taylor-Hood and ALL PASSED
    without a threshold change.** The re-baselining that looked owed was not needed; what the
    exercise actually exposed is that the suite cannot tell the two discretisations apart (see §5
    above). Original note follows for provenance:
  * ⚠ ~~**FOLLOW-UP OWED — the 17 physics/smoke tests re-baselined.**~~ `test_basic`, `test_dispersion`,
    `test_dispersion_nonlinear`, `test_nlpressure`, `test_sloshing`, `test_conservation`,
    `test_convergence`, `test_shallow_water`, `test_bc_generation{,_distributed}`, `test_bc_spectrum`,
    `test_basic_distributed`, `test_nlpressure_distributed`, and `test/local/` `test_2d_reduces_to_1d`,
    `test_boundary_modes_1d`, `test_relaxation_1d`, `test_reststate_1d`, `test_sponge_1d` all call
    `setup_and_run` **without** an explicit `p_eta`, so they moved from `Q2/Q2` to `Q2/Q1` when the
    sentinel flipped. Their pinned constants were measured on equal order and must be **re-measured,
    not re-thresholded**. ✅ The verification tier is UNAFFECTED: `test_mms_convergence*`,
    `test_jacobians_ad`, `test_linear_newton_gate`, `test_mms_distributed_parity` and
    `test_taylor_hood` all pass `p_eta` explicitly (or test it directly) and were already Taylor-Hood.
  * ⚠ **A long-duration nonlinear regression is still missing.** No test runs long enough to have
    caught this (the mode needed 50–80 s to emerge); that gap is why it reached production. The
    Taylor-Hood 80 s flume run is the natural basis for one.
* ✅ **the energy-consistent (skew-symmetric) advection correction — DERIVED, VERIFIED, REFUTED as
  the cure** (branch `fix-nonlinear-instabilities`, commit `25be653`). Kept because the *derivation*
  is a permanent asset even though the cure failed: the advection block's energy production is
  **exactly** the continuity defect,
  `n(U;U,U) = −½∫∇·(Hū)(Σ M_ij u_i·u_j) + ½∮flux`, resting on the basis-agnostic σ-tensor identity
  `½(𝓖_ikj + 𝓖_jki) = ½𝓜_ikj − ½Φ_k M_ij` (true because `ψ_k = σΦ_k − varphi_k` vanishes at **both**
  ends of the water column). Measured at the time to 8 significant figures against the assembled
  operator (closed basin: `pi_adv` = −3.04e-03 vs `pi_res` = −3.83e-11) and ★ to 3.9e-15 on five
  vertical bases.

  ⚠ **THE FEATURE ITSELF WAS REVERTED (`0e556d8`) AND IS NOT IN THE CODE.** `skew_advection`, its
  residual/Jacobian blocks, the `pi_adv`/`pi_res`/`div_flux` diagnostics and
  `test_skew_advection.jl` are all gone. Two reasons, both decisive: it is **not the cure** (treated
  and control blew up identically to 5–6 significant figures — the defect it removes is 0.4–17 % of
  a production that is itself `O(1e-3)`, five orders too small to matter), and **the LaTeX derives
  no such term** — `BALFEM_models/` is the single source of mathematical truth, and a residual term
  it does not contain is one nobody can check. The derivation is recorded here because it is a
  permanent result about the operator; the code is not.
  **Method lesson (now standing): verifying a knob is LIVE is not enough — it must be BIG ENOUGH TO
  MATTER.** Measuring `‖𝒞‖` against the energy growth rate on a stored baseline would have cost
  seconds and pre-empted 12 h of runs.
* ✅ **THE vopt DISCREPANCY IS RESOLVED AND κ IS DELETED (2026-09-08).** The functional had three
  departures from eq. (3.10); correcting them takes `M=2` from 0.0085 to **0.0036** of the published
  `c₂`. The `M≥3` residual is **not** a numerical artefact — it is converged in the scan increment —
  but a consequence of the reference specifying its median population only for `M=2`. `VOPT_KAPPA`
  and its machinery are gone: κ was fitted to Table 1 and then read off its own fit. Complete
  account: **rule 46**. Data: `output/local/vopt/`. ⚠ Consequence: `DEFAULT_CBDY_P` no longer exists,
  so **`p ≥ 2` has no default mesh** — that is deliberate, and a `p ≥ 2` design needs a justified `Ω`
  first.
* 🟢 **`quad_extra`** (`setup_and_run` kwarg, `BALFEM_QUAD_EXTRA` in the 1-D driver; **threaded
  through `run_conv_study → run_mms_case → Measure` on 2026-09-12** — the MMS path had hard-coded the
  degree, so the claim below had never been testable there). The default `2·max(p_u,p_η)+2` is
  genuinely **one short** of the nonlinear integrands.
  ⚠ **"Raising it changes no MMS rate" is now MEASURED, and the reason is effect size, not deadness.**
  Assembling the residual at degrees 8/10/12/16 (`output/local/mms/quadrature_test/knob_liveness.log`):
  the LINEAR core moves 4e-14 across all degrees (exact at the default), the NONLINEAR core moves
  **3.49e-08** from degree 8→10 and 4e-14 thereafter — so the crime is real and the knob is live. Its
  effect on converged `e_η` is ≤4.7e-07 relative (~4.6e-13 absolute), against the **3.1e-07 absolute**
  change needed to restore third order: **six orders too small**. Default stays 0.
  ⚠ **Gridap's `Measure(trian, d)` is exact to degree `d+1`, not `d`** (measured: 8→x⁹, 10→x¹¹,
  12→x¹³). Read exactness from a probe, never off the argument — a liveness check built on `x⁹`
  cannot fail at either degree and silently certifies nothing.
* ✅ **the two `src/mms_driver.jl` defects are FIXED** (found 2026-08-19, fixed 2026-08-21).
  **A1** — `run_conv_study` hard-coded `assemble_vertical_tensors(M, 1, [0,0.728,1])`, so `p_vert`
  was not a parameter and any `M≠2` threw. It now takes `p_vert` and `c_bdy`, resolving through the
  single shared `resolve_cbdy` that `setup_and_run` and its distributed twin also use; the tag
  carries `P{p}LFE-{M}` so an `(M,p)` sweep cannot produce two studies under one label. Verified:
  `M=3` runs and reaches optimal order (P1LFE-3, Q2/Q1: `p_η` 1.994/2, `p_u` 2.997/3).
  **A2** — `run_mms_case_distributed` hard-coded **Model 1** (`regime=:linear, nl_pressure=:none,
  flat_bed=true` + `mms_forcing_stage1`) while `run_conv_study` built its tag from the *requested*
  switches, so a distributed 8-model campaign returned 8 identical Model-1 studies under 8 different
  labels, **all passing**. The switches and `hfun` are now parameters feeding the general
  `mms_forcing`, from the same variables that feed the solver. `test/test_mms_distributed_parity.jl`
  (4 ranks) gates the two branches against each other on two non-Model-1 configurations — with a
  **separation negative control first**, because a bare parity check would pass with the defect
  present. This unblocks the vertical-basis convergence study (`COMPLETED_VBASIS_STUDY.md` §1 prerequisite A);
  the study driver is `examples/local_mms/run_vertical_basis_study.jl`.
* ✅ **the MMS `:full` omission is FIXED** (found 2026-08-21, fixed 2026-09-01). Both MMS drivers now
  build the `nlp` context, and both time loops **prime it from the IC** — `update_nlp_state!` runs
  only *after* an accepted step, so step 1 assembled as if from rest (exact for a rest start, wrong
  for `u0 = u*(t0)`). `e_u` fell **1.19e-03 → 1.37e-06** and `p_u` went **−0.00 → 1.948**, so the
  old floor measured **omission**, not the frozen-projection lag `VERIFIED_SCOPE.md` §4 describes.
  ⚠ **A second, deeper limit is now exposed**: with the blocks in the residual but absent from
  `jacobian_u`, the quasi-Newton gap has a **cliff in amplitude** — Newton stalls at 9.2e-04 at the
  campaign's `a_eta = 0.8` and reaches ~2e-09 only at `a_eta ≤ 0.4`. `run_conv_study` gained an
  `a_eta` kwarg; **tier-3 floors are not comparable across amplitudes.**
* ✅ **cluster memory attribution — H3 (per-step leak) is REFUTED at production scale.** The
  2026-08/09 small-domain runs show RSS rising over the first ~100 steps and then **plateauing**
  (plane 1894→~3230 MB, directional 1546→~2800 MB), with a peak of **3633 MB/rank**. That is H4:
  the baseline footprint simply exceeds 2 GB/core, so **4 GB/core is required and permanent**, and
  3 GB would OOM. ⚠ Distinct from the *long-lived-process* drift of rule 41, which is real and
  unbounded — these are per-run measurements, not multi-hour worker lifetimes.
* 🔴 `test_mms_convergence` G7 — a gate-window specification decision, not a fix
* 🟠 **the `:sdirk` LINEAR-model temporal deficit is unexplained.** Across the completed matrix
  `:sdirk` reaches `pw_u ≈ 1.99` on the nonlinear models but only ≈1.69 on the LINEAR ones at
  `dt0=0.15`, while `:theta` gives ≈1.99 for both. **It is not spatial-floor contamination** (the
  temporal errors sit four orders above the floor at `nx=36`), and it is **not a dissipation penalty
  in general** — a matched pair at `dt=0.05` has `:sdirk` and `:theta` agreeing to 0.002. The one
  genuinely open question the vertical-basis campaign leaves behind.
* 🟠 **a per-basis `dt` ladder is required for temporal studies, and is not yet in the test suite.**
  A single ladder diverges outright on the richest basis; `run_vbasis_shard.jl` scales `dt0` with
  `Nσ`, but `run_conv_study` and the gates still use one fixed ladder.
* 🟠 no MPI tier in `runtests.jl` (cost three stale reference constants)
* 🟠 preconditioner replacement — the single biggest performance item
* 🟠 four run-output gaps
* ✅ *closed by the vertical-basis campaign:* the `(M,p)` convergence study itself
  (`COMPLETED_VBASIS_STUDY.md` §1), tier 3's `:full`-floor-vs-`Nσ` question, and `OPEN_ISSUES.md` §6's
  velocity-shortfall question for the vertical-basis case (pre-asymptotic — the horizontal `Q2/Q1`
  and `Q4/Q3` pairings are still untested on an extended ladder, and the recipe is now cheap)
* naming follow-through outside this checkout: ✅ the GitHub repo is renamed and the local remote
  URL was updated to `git@github.com:Elmanyer/GridapBALFEM.jl.git` (2026-09-06); **still outstanding:**
  the cluster checkout and the sysimage rebuild

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
