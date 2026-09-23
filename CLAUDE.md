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
> | [`NEW_TREATMENT.md`](markdown_files/NEW_TREATMENT.md) | ⚙ **the Class-III replacement** (branch `new-classIII-treatment`): the exact `{1,2,5}` algebraic reduction, the in-loop static-condensation projections, the derivation, the implementation plan and the test programme |
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
> **⚙ THE `:full` CLASS-III INSTABILITY IS HALF-EXPLAINED (2026-09-23, branch
> `mixed-formulation-solver`).** Replacing the frozen `L²` Class-III projections with a genuine
> mixed unknown `𝖦 ≈ ∇𝖲` — a formulation that **never differentiates a `C⁰` field twice** — moves
> Q2/Q1 flat `:full` at `A`=0.10 from a **12.6 s** death to **>100 s (62.9 periods), flat in the
> fourth decimal with Newton constant at 6**, and beats the projected control in **all eight**
> matched pairs. **So the projection was a large real contributor.**
> ⛔ **BUT IT IS A DELAY, NOT A CURE.** The `A`=0.15 twin of that same run **diverged at ~87 s** with
> the identical velocity-led signature. A pass at one amplitude is not a pass; never quote the
> `A`=0.10 trace without it.
> ⛔ **It is not the whole cause.** Every **Q3/Q2** arm still dies, *earlier* than its Q2/Q1 twin and
> earlier even at matched DOF count; and the **`dx` refinement signature survives**, which refutes
> the broken-Hessian recovery explanation of §5.2c **as stated** (rule 39b) and removes the stated
> rationale for `C⁰`-IP.
> ⚠ **DO NOT ATTRIBUTE THE Q3/Q2 FAILURE TO CLASS III YET** — every Q3/Q2 arm run is `:full` and
> **no `:native` Q3/Q2 control exists** (rule 14c). That one cheap run gates every further
> conclusion. Full account, with all four hypotheses and their verdicts: **§5.2e**; design and
> per-arm detail: `markdown_files/NEW_TREATMENT.md` Part H.
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
| `test/` | ⚠ **+`test_yl_collapse.jl` 2026-09-16** — stage 1 of the Yang & Liu collapse (40/40, 0.5 s), the **only gate that checks our derivation against an independently published one**; everything else is self-consistency (rule 28). Plus `test/yl_collapse_wip/` (stages 2–3, self-contained, not a gate). **30** test files + `runtests.jl` + `test/cluster/` + `test/local/` — inventory and scores in `TEST_SUITE.md`. `test_mms_distributed_parity.jl` (4 ranks) gates the distributed MMS path against the sequential one; `test_taylor_hood.jl` (`:fast`, 13/13) gates the horizontal pairing — Taylor-Hood accepted, equal order REJECTED |
| `examples/` | 7 sequential + `distributed/` (7 cluster scripts + `_dist_common.jl`), `distributed_small/` (5 parametric), `validation/` (7), `local_1d/`, `local_2d/`, `local_mms/` (**11** — the parametric MMS studies, the vertical-basis campaign (`run_vbasis_campaign.jl`, `run_vbasis_shard.jl`, `report_vbasis_campaign.jl`), and the Phase-B batch: `run_phaseB_shard.jl` + `supervise_phaseB2.sh`, a resume-capable runner with a memory-aware supervisor), `inspect_run.jl` — `RUNNING.md` §2 |
| `run/` | ⚠ **+9 local launchers 2026-09-15**: `run_1dprod_*.sh` (the §5.2b campaign — `nl_{native,full}_{flat,bar}`, the shoulder ladder `_bar_s10/_s15/_gentle`, the `_ad` diagnostic) + `run_all_1dprod.sh`. 10 production SLURM launchers + `run/dist_small/` (**30** small-domain 2-D cases — the 7 superseded 1-D twins were deleted 2026-09-02) + `run/local/` (**58** scripts: the original 30 incl. **6** 1-D cases + `run_all_1d.sh`, plus the **16** `run_1dnl_*.sh` nonlinear-instability ladder — 4 physics tiers × A∈{0.05,0.10,0.15,0.20} — + `run_all_1dnl.sh`; ⚠ **all 16 of those ran EQUAL ORDER and their results are void — see rule 12b**; plus the **6** `run_skewA_*.sh` (energy-consistent advection, refuted) and the **9** Taylor-Hood campaign scripts `run_th_*.sh` / `run_thref_*.sh` + `run_all_th.sh` that resolved the instability), all through `run/balfem_env.sh` (cluster) or `balfem_local.sh` (workstation) — `RUNNING.md` §3–4. **Both helpers now gate the Taylor-Hood pairing before launch and print it into the banner** (rule 2b). ⚠ **Job sizing is governed by [`run/SNELLIUS_ROME_LAUNCH_CONFIGS.md`](run/SNELLIUS_ROME_LAUNCH_CONFIGS.md)** — TRACKED, because it decides what every job COSTS. A `rome` node is 128 cores / 224 GiB divided into **eighths** (1/8 = 16 cores / 28 GiB), and the bill is `max(cores/128, mem/224)` rounded UP to the next eighth. At the 4 GiB/rank this solver needs, **memory always sets the tier**, so the cost-optimal count is **7k ranks for tier k/8** — 28 ranks bills 4/8 where 32 billed 5/8, and 42 bills 6/8 where 48 billed 7/8. All 38 launchers were resized on 2026-09-07; `--mem-per-cpu=5G` is never valid here |
| `compile/` | the cluster sysimage build chain — `RUNNING.md` §5. ⚠ **`Manifest.toml` is GITIGNORED**, so a fresh cluster checkout has only `Project.toml`; `set_preferences.jl` therefore `Pkg.add`s both forks by URL (`Elmanyer/Gridap.jl` @ `fix-transient-multifield-ad` for transient-multifield AD, `Elmanyer/WaveSpec.jl` for the `change_seed!` fix) — `Pkg.instantiate()` would resolve the BROKEN upstream versions. Consequence: builds track the branch tips and are not reproducible across time |
| `output/` | ⚠ **GITIGNORED** — nothing here survives a fresh clone, so anything load-bearing must be copied into a tracked file. Result records worth knowing about: `output/local/vopt/` (the corrected mesh-optimisation results + `README.md`, rule 46), `output/local/mms_vbasis/` (vertical-basis convergence, **10/10 PASS** on the published `p=1` meshes — accepted, not to be re-run), `output/local/mms_campaign/` (the 2026-08-30 five-basis campaign, 30/30 spatial), `output/local/mms_phaseB/` (live batch + `INVALID_P2LFE-2.md`) |
| `postprocessing/` | `GridapBALFEMPost` — self-contained, own environment, **no dependency on the solver** |
| `WaveSpec.jl/` | vendored stochastic sea-state synthesis (CMOE-TUDelft). Tracks the **GitHub repository version, not a tagged release** — the release's `change_seed!` is broken |
| `Gridap.jl/` | the vendored **fork** (`Elmanyer/Gridap.jl` @ `fix-transient-multifield-ad`, one commit on `v0.20.8`) making transient-multifield AD work — `CONFIGURATION.md` §2 |
| `latex_docs/BALFEM_models/` | **THE CURRENT LaTeX project** — the authoritative derivation (§2). ⚠ the `.zip` no longer exists in the checkout; re-export from Overleaf before assuming a round-trip |
| `latex_docs/LFEM_discretisation.zip` | **SUPERSEDED** pre-rename LaTeX, kept for provenance as a zip only. Do not edit |
| `latex_docs/CFC2027_abstract/` | conference abstract (CFC 2027); its class needs the `newtx` fonts. **Renamed from `CFC2027_LFEMultilayer_abstract/`** |
| `latex_docs/doc_figures/` | `generate_doc_plots.ipynb` — one cell per figure, activating the repo project by walking up to `Project.toml`. Reproduces Yang & Liu figs 2 and 3 from `assemble_dispersion_tensors` + `model_R`; recovers their published `kd_app` (10.84 / 39.23 / 127.91) |
| `markdown_files/` | the documents listed in the START HERE table — **16 files, TRACKED** (`git ls-files markdown_files/`). ⚠ An earlier revision of this row and of §8 said *gitignored*; that was true until 2026-09-11 and is now wrong in the direction that loses work — load-bearing content belongs **here**, and `latex_docs/` is the gitignored one |
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
BCs (the 1-D bar now carries `BALFEM_SBAR`, its shoulder length: `max|∇h| = hbar/(2·sramp)`,
small ⇒ square cross-section — the ladder parameter of §5.2b); runtime monitoring plus an independent governing-equation residual checker; `w_s`/`p_s` VTK
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
*Amended 2026-09-15 with the 1-D nonlinear production campaign — §5.2b, §5.6b.*
*Amended 2026-09-23 with the mixed-formulation campaign — §5.2e, and §5.7 item 0 rewritten.*

### 5.0 Status in one paragraph

The solver is **feature-complete in serial and distributed form** and has been so since 2026-09-02:
stacked loop-free residual, hand Jacobians, the full nonlinear physics, SDIRK/θ integrators, every
boundary treatment, Dirichlet generation with WaveSpec coupling. **The LINEAR models are in excellent
shape** — optimal order in both fields, on five vertical bases, in 1-D and 2-D, sequential and
distributed. **The NONLINEAR models are correct but carry one unexplained order reduction** in the
surface elevation at fine mesh, isolated as of 2026-09-13 to the advection block and *not* explained
by any of the six candidate causes tested. Production wave runs are stable at realistic amplitude
since the Taylor-Hood fix, with the exception of boundary-generation failures that are a separate,
configuration-level problem. ⚠ **AMENDED 2026-09-15:** that last sentence is now known to hold only
for `nl_pressure ∈ {:none, :native}` ON A FLAT BED. The 1-D production campaign (§5.2b) found two
distinct, reproducible instabilities — `:full` fails on a flat bed in 8 wave periods, and ANY
variable bed grows a lee-shoulder mode whose onset time is set by `|∇h|` — while the `:native`
flat-bed case completed **100 wave periods at `A = 0.10` m**.

⚠ **AMENDED 2026-09-23 — `:full` IS NO LONGER UNCONDITIONALLY UNSTABLE, AND THE CAUSE IS ONLY
PARTLY IDENTIFIED (§5.2e).** Replacing the frozen `L²` Class-III projections with a genuine mixed
unknown (`𝖦 ≈ ∇𝖲`) turns the Q2/Q1 flat `:full` case from a **12.6 s** death into **≥61 wave
periods flat on the `:native` trace**, and outlives the projected control in **all eight** matched
pairs. So the projection was a large real contributor. **But it was not the whole cause:** every
**Q3/Q2** arm still dies — *earlier* than its Q2/Q1 twin, and earlier even at matched DOF count —
and the **`dx` refinement signature survives a formulation that never differentiates a `C⁰` field
twice**, which refutes the broken-Hessian recovery explanation of §5.2c as stated. ⚠ **The
attribution to Class III at Q3/Q2 is NOT established**: every Q3/Q2 arm run is `:full` and no
`:native` Q3/Q2 control exists (rule 14c). That one cheap run gates every further conclusion.

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
| long-run stability | **100 wave periods** at `A = 0.10` m (`κa ≈ 0.16`, `kd = 5.5`) — extended 2026-09-15 from the original 50. `η` settles 0.10494 → 0.10300 across four windows, Newton flat at **5.12 it/step** start to finish. ⚠ FLAT BED AND `:native` ONLY (§5.2b) |
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

### 5.2b ⛔ THE 1-D NONLINEAR PRODUCTION CAMPAIGN (2026-09-15) — two instabilities

Eight sequential runs, `P1LFE-2`, Q2/Q1 Taylor-Hood, SDIRK_2_2, boundary-generated regular wave,
`A = 0.10` m, `kd = 5.5`, 60 m flume at `dx = 0.25`, `ny = 1` + walls, 100 periods (4000 steps),
VTK every 0.2 s. Launchers `run/local/run_1dprod_*.sh`; output under `output/local_1d/`.

| case | outcome |
|---|---|
| `:native`, flat | ✅ **4000/4000 steps to t = 160 s**, η 0.10494 → 0.10300, Newton 5.12, 800 VTK frames |
| `:native`, bar | ⛔ lee-shoulder mode; onset set by bed slope (ladder below) |
| `:full`, flat | ⛔ t = 12.6 s |
| `:full`, bar | ⛔ t = 12.6 s — **dies before the wave reaches the bar**, so it is the flat failure |

**⛔ DEFECT 1 — `:full` IS AN OPERATOR FAILURE, NOT A QUASI-NEWTON ONE. THIS IS THE IMPORTANT ONE.**
The hand-Jacobian run stalls at t = 12.6 s (Newton 6 → 12 → 19 → 50-cap at `‖r‖` = 1.9e-05), which
reads as the known incomplete `jacobian_u`. **It is not.** Re-run with `BALFEM_USE_AD=1` — the
*complete* Jacobian — Newton holds **4 it/step to residual 1e-09 all the way to the same step**, and
then the state goes to **NaN**. `r0` grows identically in both (1.09 → 1.47 → 2.33), i.e. the
trajectory belongs to the residual, not the Jacobian (rule 17b).

The mechanism is **velocity-led at the generation boundary**, with the `:native` control flat:

```
                    :full (both Jacobians)      :native (control, same t)
t = 11.2   u_max =  0.4465                       0.4194
t = 12.4            0.9266                       0.4196
t = 12.6            2.4971  <- 6x in ~1 s        0.4273
   eta_max          0.124 (BOUNDED throughout)   0.109
   location         x = 0.50 m, the inflow       --
```

⚠ **CONSEQUENCE: completing `jacobian_u` will NOT make `:full` stable** — AD *is* the completed
Jacobian. It may still be needed for MMS verification; it is not the fix for this. §5.7 item 5 was
written on the opposite assumption and is corrected there.
⚠ **This is a cheap 1-D reproduction of the cluster signature** (§5.6: "velocity-led — η stayed
bounded while `u` exploded"), in 13 s of simulated time, and **with** the `:native` control that
batch lacked.

**⛔ DEFECT 2 — ANY VARIABLE BED GROWS A LEE-SHOULDER MODE; `|∇h|` SETS THE RATE, NOT A THRESHOLD.**
Five runs, identical but for the bar's shoulder length (height 2.0 m on `d` = 3.5 m, span 26–34 m,
crest depth 1.5 m; `max|∇h| = hbar/(2·sramp)`):

| shoulder | `max\|∇h\|` | face | onset | periods |
|---|---|---|---|---|
| 0.5 m | 2.0 | 63° | t = 37.2 s | 23 |
| 1.0 m | 1.0 | 45° | t = 57.8 s | 36 |
| 1.5 m | 0.67 | 34° | t = 93.8 s | 59 |
| 2.0 m | 0.5 | 27° | t = 137.2 s | 86 |
| flat | 0 | — | **none** | ✅ 100 |

One mode, four growth rates: growth **pinned at the downwave shoulder** (x ≈ 32–34.5 in all four),
`u_max` 2.2–3.3 against the flat control's steady 0.42, η to 0.8–3.0 m from a 0.10 background,
Newton degrading 12 → 15 → 20 → cap. **Onset falls monotonically with slope and the flat bed is the
limit point of the same family** — a rate, not a threshold, which is the signature of an instability
in the formulation rather than of a badly-posed bathymetry. The `t = 93.8` point was a PREDICTION
(it had to fall between 58 and 137) and was met.

⚠ **NOT YET MEASURED, AND IT IS THE DECIDING ONE:** halve `dx` on a bar case. Rule 38b — refinement
DELAYING onset means under-resolution of the bed; refinement ADVANCING it means a grid-scale problem
in the ∇h terms, which the pinned location and velocity involvement already suggest.

⚠ **Caveats on the positive result.** SDIRK_2_2 is L-stable, so the flat-bed run's slow −1.9e-05 m/s
decline is consistent with numerical damping and "stable" is not unqualified until `:theta` or
explicit RK4 confirms it (rule 12c). And it is Q2/Q1, not the convergence-optimal Q3/Q2.

### 5.2c THE `:full` FLAT REFINEMENT FACTORIAL (2026-09-16) — `dx` ADVANCES ONSET, `dt` DELAYS IT

2 (`dx`) × 2 (`dt`) × 2 (Jacobian) on the case that dies at t = 12.6 s, plus an extended `dt` ladder.
`P1LFE-2`, Q2/Q1, SDIRK_2_2, `A` = 0.10 m, `kd` = 5.5, 60 m flume, 16 periods. Launchers
`run/local/run_1dref_full_*.sh`.

**Onset time (s), hand Jacobian; AD in parentheses where finished:**

| | `dt` = 0.04 | `dt` = 0.02 |
|---|---|---|
| `dx` = 0.25  | **12.60** (12.60) | **21.20** (running) |
| `dx` = 0.125 | **5.00** (5.00)   | **8.20** (8.40) |

* **Halving `dx` ADVANCES onset ×2.5–2.6.** Rule 38b: refinement making things worse is a
  **grid-scale discretisation-stability problem, never under-resolution**.
* **Halving `dt` DELAYS onset ×1.6–1.7.**
* **Not a CFL effect**: (0.25, 0.04) and (0.125, 0.02) hold `dx/dt` fixed and give 12.60 vs 8.20.
* **The two Jacobians agree** — see rule 17b, which this factorial proves.

**The extended `dt` ladder at fixed `dx` = 0.25** (four points, hand Jacobian):

| `dt` | 0.04 | 0.02 | 0.005 | 0.001 |
|---|---|---|---|---|
| onset (s) | 12.60 | **21.20** | **40.40** | > 25.6, uninformative |

⚠ **THE FIRST THREE ARE ONSETS; THE FOURTH IS NOT.** The `dt` = 5e-3 case was re-run to 100 periods
and **crashed at t = 40.40 s** (25.2 wave periods), with `x_at_max` migrating at `c_g` up to x ≈ 46.8
and then **snapping back to x ≈ 41 and pinning** while η exploded — post-fill growth rate
**σ = +0.34/s** against the **+0.022/s** transit artefact, a factor of fifteen. The `dt` = 1e-3 run
completed its 16-period window (25 600 steps, Newton flat at 2) but **25.6 s is inside the 36 s
fill**, so it carries no information about stability at all and its ">25.6" is not an onset.
⚠ Note also that the 0.04 and 0.02 onsets both fall **inside** the fill, and only the 5e-3 one is
post-fill — which is the likeliest reason the mode pins at the inflow in the first two and near the
front in the third.

⚠ **The `dt` = 5e-3 entry is a BOUND, AND ITS 16-PERIOD WINDOW WAS TOO SHORT TO READ AT ALL.**
It completed 5121/5121 steps to t = 25.6 s. ⚠ **A first reading called it "still growing" (η 0.112 →
0.125 across windows); THAT WAS RULE 14's TRANSIT TRAP AND IS WITHDRAWN.** The flume fills in
**36 s** (45 m to the right sponge at `c_g` = 1.25 m/s), so a 25.6 s window ends *mid-fill*:
`max|η|` creeps up because the front is still advancing, and `x_at_max` migrates at `c_g` rather than
sitting still. The 100-period re-run shows exactly that — `x_at_max` 15 → 22 → 29 → 35 → 42 m at
~1.1 m/s — and a naive exponential fit over the fill returns **+0.022/s**, against the **+0.028/s**
the driver already documents as the spurious transit rate. ⚠ **NOTHING before t ≈ 36 s in ANY
boundary-generated run on this flume may be read as growth** (rule 14; this is the fourth
measurement it has caught).

**READING.** Onset moves out as `dt` falls and in as `dx` falls. The Class-III blocks `{1,2,4,5}`
are assembled through **frozen `L²` projections refreshed once per step**, so their lag error falls
with `dt` — while the recovered **second** derivative of a `C⁰` field is not an `L²` object at all
(it carries a Dirac layer on the element skeleton), so its recovery does not improve with `h`.
⚠ **The competing explanation is dead**: the `𝓝` operator is now verified against the governing
equations (§5.2d), so this is about the **assembly**, not the derivation. It remains a hypothesis;
what is established is the two refinement signs and the Jacobian irrelevance.

⚙ **AMENDED 2026-09-23 — THE `dx` HALF OF THIS READING IS REFUTED (§5.2e H4).** The `dx` signature
**survives the mixed formulation**, which never differentiates a `C⁰` field twice and therefore
cannot commit the recovery error blamed below. The **`dt`/lag** half stands and was confirmed (both
`dt` arms improved once the lag was removed by construction). **Read the next paragraph as the
historical hypothesis it now is, not as the explanation.**

⚠ **THE TWO SIGNS ARE TWO DIFFERENT ERRORS AND MUST NOT BE CONFLATED.** The **recovery** error (the
recovered `∇s` is not the true `∇s`; the exact `∂²` of a `C⁰` field is `{∂²u}` cell-wise **plus a
Dirac layer on the skeleton**, which no cell quadrature sees) carries the `dx` signature. The **lag** error (`update_nlp_state!` runs only *after* a step,
so the residual is evaluated against the previous step's state) carries the `dt` signature — and it
is the stronger, cleaner trend. **The lag is intrinsic to *freezing*, not to projecting**, and can be
removed without changing the function spaces at all.
⚠ **Every run in this factorial is Q2/Q1**, where `η` is piecewise linear so `∂²η ≡ 0` *identically*
in 1-D — so component 4's free-surface half, and the `∇H` half of `∇s`, contribute nothing there —
and `H·u` is cubic per element, so the `∂²(Hu)` being recovered is only piecewise **linear**. That is
the crudest discretisation that can carry Class-III content at all, and the Q3/Q2 repeat has not been
run. Full account: `OPEN_ISSUES.md` §0c.

### 5.2d ✅ BALFE-M AT `p=1` COLLAPSES ONTO YANG & LIU'S LFE-M — VERIFIED TO ROUND-OFF (2026-09-16)

The supplementary material of Yang \& Liu — free with the open-access article, now kept beside the
PDF — publishes the **non-linear** coefficients: §A the vertical velocity, §B the pressure, §C the
residual. That makes the collapse checkable coefficient by coefficient, and it was.

| object | source | worst relative difference |
|---|---|---|
| continuity weights `Φ` | (2.31) | **0** |
| linear coefficients `A`, `B`, `D` | Appendix A | **7.8e-15** |
| vertical velocity `w` (non-linear) | §A (2.19)+(A.1)–(A.2) | **1.0e-14** |
| non-hydrostatic pressure `p_nh` | §B, via `∂p/∂σ = −ρH·Dw/Dt` | **1.7e-15** |
| weighted momentum residual | (2.24)+(2.26) | **1.2e-13** |

**⚠ TWO OF THESE ARE PHYSICS CHECKS, NOT CROSS-CHECKS.** The pressure and residual stages test our
package against the vertical momentum equation and the projection of (2.12) — so an error that our
derivation and theirs happened to *share* would still have been caught. The last stage uses the
**solver's own assembled tensors** (`Mmat, Mcal, Gcal, A, K, P, Acal, Kcal, Pcal`), so it verifies
`src/vertical.jl` and not only the LaTeX.

**This is the only external oracle the project has.** Everything else in the suite is
self-consistency (rule 28) and would pass for a consistently-wrong residual.
**The `𝓝` package, including the Class-III set `{1,2,4,5}`, is verified.**

⚠ **SCOPE, and it must be quoted with the claim.** The comparison is **1-DH only** (their
supplementary is "in 1DH version for brevity"), so the 2-D structure — the tensor character of
`∇uⱼ`, cross terms, both components of the `𝓜/𝓖` contractions — is NOT exercised. It is a **`p=1`**
statement by construction, so it says nothing about the basis-agnostic claim. And it compares
**operators at a prescribed state**: the Class-III *assembly* is untouched, which is precisely what
§5.2c is about.

⚠ **The state must satisfy depth-integrated continuity.** `pDerivation.tex` eliminates `∂H/∂t` via
`∂H/∂t = −Σⱼ∇·(Huⱼ)Φⱼ`, both directly (Θ components 1,2) and inside `ω` (6,7,8), so our package is
an identity ONLY on that constraint; Yang & Liu keep `H^(0,1)` explicit and are state-independent.
A test with `H` and `uⱼ` prescribed independently makes OUR side fail by `O(ε²)` **as an artefact of
the test**. This cost a day and a retracted defect claim; it is the first thing to check if the
comparison is ever re-run. Full account: `latex_docs/BALFEM_models/ModelVerification.tex` (chapter 4)
and `test/yl_collapse_wip/README.md`.

### 5.2e ⛔ THE MIXED-FORMULATION CAMPAIGN (2026-09-22/23) — THE PROJECTION IS A LARGE PART OF THE CAUSE, BUT NOT ALL OF IT

*Branch `mixed-formulation-solver`. Design: `NEW_TREATMENT.md` Parts F–H. 15 arms on the real
60 m flume, `nl_pressure=:full`, flat bed, SDIRK_2_2, boundary-generated regular wave, 100 periods
requested. Every arm paired with its own projected control in the same batch (rules 14c, 38c).*

**THE CONSTRUCTION UNDER TEST.** The Class-III blocks `{1,2,5}` collapse algebraically onto a single
contraction in `∇𝖲`, `𝖲 = ∇·(Hu)` (§A of `NEW_TREATMENT.md`, exact to 4.4e-16). The **projected**
path recovers `∇𝖲` by a frozen `L²` projection refreshed once per step. The **mixed** path carries
`𝖦 ≈ ∇𝖲` as a genuine FE unknown in a 5-field layout `[η,𝖴x,𝖴y,𝖦x,𝖦y]`, defined by its own weak
equation integrated by parts, `∫𝖦·Ψ = −∫𝖲 ∇·Ψ + ∮𝖲Ψ·n`. **The mixed path never forms a second
derivative of a `C⁰` field at all** — that is the entire point of it, and it is what makes the
comparison a clean test of the projection rather than of the operator.

#### The four hypotheses, and what each measurement returned

**H1 — "the frozen projection is what destabilises `:full`."**
Test: replace it with the mixed unknown, everything else held fixed, against controls in the same
batch. **PARTIALLY CONFIRMED — the effect is large and real in all eight matched pairs, but it is a
DELAY, not a cure** (the `A`=0.15 arm went on to diverge at ~87 s; see the final-state block below):

| arm | pairing | varied | mixed onset | projected onset | gain |
|---|---|---|---|---|---|
| `base` | Q2/Q1 | — | ✅ **running >97.4 s** | 12.60 (documented) | **>7.7×** |
| `amp` | Q2/Q1 | `A`=0.15 | ✅ running >93.2 s ⚠ | 5.8 | >16× |
| `dt` | Q2/Q1 | `dt`=0.02 | ✅ running >57.3 s | 21.2 | >2.7× |
| `dx` | Q2/Q1 | `nx`=480 | ⛔ 36.6 | 5.0 | 7.3× |
| `base` | Q3/Q2 | — | ⛔ **16.0** | 7.0 | 2.3× |
| `amp` | Q3/Q2 | `A`=0.15 | ⛔ 9.6 | 4.6 | 2.1× |
| `dt` | Q3/Q2 | `dt`=0.02 | running 19.3 | 10.4 | >1.9× |
| `dx` | Q3/Q2 | `nx`=480 | ⛔ 6.4 (NaN) | 3.4 | 1.9× |

The Q2/Q1 base arm holds `η` at **0.1112–0.1115** and `u_max` at **0.402–0.412** flat across
t = 40–100 (six consecutive 10 s windows) with Newton at 5–6 — i.e. it is not merely surviving, it
is on the `:native` reference trace, at the tier that previously died in 12.6 s. **This is the
strongest single result on the branch.**

**H2 — "Q3/Q2 will improve matters, because Q2/Q1 under-represents the Class-III content"**
(`OPEN_ISSUES.md` §0c option 0; §5.7 item 0 step 0; the reasoning was that `η ∈ Q1` makes `∂²η ≡ 0`
identically in 1-D, so part of the operator never contributes).
**⛔ REFUTED, AND WITH THE SIGN REVERSED. Every Q3/Q2 arm dies earlier than its Q2/Q1 twin** —
16.0 vs >97.4, 9.6 vs >93.2, 6.4 vs 36.6, and the `dt` pair is tracking the same way. Raising the
polynomial order makes the instability **worse, monotonically, in all four cells.**

⚠ **The mechanism is not "more nonlinear effects."** The residual is the same operator at both
pairings — same `:full`, same eight `𝓝` components, no physics is switched on. What Q3/Q2 changes is
what the discretisation can **represent**: with `η ∈ Q2`, `∇η` is piecewise linear instead of
piecewise constant, so the `∇η·u` half of `𝖲` — and hence the Class-III content `𝖦` is built from —
becomes a real object for the first time. **The Q2/Q1 runs are therefore not "the stable case"; they
are solving a partially masked operator, and the >97.4 s result must be quoted with that caveat.**

**H3 — "the Q3/Q2 penalty is just a resolution effect"** (Q3 on the same cells is 2.1× the DOFs, and
by rule 38b more resolution advances onset in this failure).
Test: the batch already contains a near-matched-DOF pair, which was not designed but is decisive:

| run | pairing | `nx` | `dx` | free DOFs | onset |
|---|---|---|---|---|---|
| `c3v_dx_mixed` | Q2/Q1 | 480 | 0.125 | **29781** | **36.6 s** |
| `c3q_base_mixed` | Q3/Q2 | 240 | 0.25 | **31710** | **16.0 s** |

**At matched problem size, Q3/Q2 still dies 2.3× earlier.** So the order penalty is **not** merely
DOF count or effective resolution — the polynomial order carries something of its own, which points
back at H2's representable-content reading. ⚠ **Partial isolation only:** the two rows have the same
DOF count but different `dx` (0.125 vs 0.25), so the grid scale — where a `λ ≈ 2–3·dx` mode would
live — is *not* matched. A proper isolation needs Q3/Q2 at `nx` = 120.

**H4 — "the `dx` signature is the recovery error, so the mixed formulation should remove it."**
This was the standing explanation of §5.2c: the exact `∂²` of a `C⁰` field is `{∂²u}` cell-wise plus
a Dirac layer on the skeleton that cell quadrature never sees, and that error does not improve with
`h`. **⛔ REFUTED. The `dx` signature SURVIVES the mixed formulation:** `c3v_dx_mixed` (nx=480) died
at 36.6 s while `c3v_base_mixed` (nx=240) is past 97.4 s — halving `dx` still advances onset, by
≥2.7×, in a formulation that **never differentiates a `C⁰` field twice.**

⚠ **THIS IS THE MOST CONSEQUENTIAL DEDUCTION OF THE CAMPAIGN, AND IT IS A DEDUCTION, NOT A
MEASUREMENT.** If the mixed path structurally cannot commit the recovery error and the `dx`
signature persists anyway, then at least one of the following holds, and they are not yet separated:
1. the `dx` signature was **never** the broken-Hessian recovery error;
2. `𝖦`'s own FE approximation of `∇𝖲` carries an `h`-dependent error of the same sign;
3. there is a genuine **grid-scale instability in the coupled `(η,u,𝖦)` system** — a pairing
   problem in the new `𝖦↔u` block rather than in `η↔u`.
**Reading (3) has direct support in the code and none of it reassuring:** `Vaux` is built from
`reffe_U` (`src/horizontal.jl:59`), so **`𝖦` sits in the velocity space at equal order with `u`** at
both pairings. `𝖦`'s own block is a Gram matrix and therefore coercive, so this is not the classic
inf-sup failure of rule 2b — but the coupled pairing has never been analysed, and the one thing that
would test it, `𝖦` one order lower, has still never been run (`NEW_TREATMENT.md` F.3.2, open since
it was written).

Conversely, the **`dt` half did improve** — both `dt` = 0.02 arms outlive their controls (>57.3 vs
21.2; 19.3 vs 10.4) — which is consistent with the lag error being genuinely removed by construction,
since the mixed unknown is solved for at the current iterate rather than frozen from the last step.

#### ⚙ FINAL STATE AT 2026-09-23 15:29 — AMPLITUDE SPLITS THE Q2/Q1 RESULT IN TWO

*Captured before a forced machine suspend/restart. These are the last numbers taken; the four LIVE
arms had not finished, and `output/local_1d/<arm>/diagnostics.csv` + `run.log` are on persistent
disk.*

| arm | state | t | η_max | u_max | Newton |
|---|---|---|---|---|---|
| `c3v_base_mixed` | LIVE | **100.6 s (62.9 T)** | 0.1116 | 0.403 | 6 |
| `c3v_amp_mixed` | LIVE | 95.2 s | **0.2300** | **1.133** | **15** |
| `c3v_dt_mixed` | LIVE | 58.8 s | 0.1290 | 0.417 | 6 |
| `c3q_dt_mixed` | LIVE | 19.7 s | 0.1312 | 0.501 | 8 |

**✅ `A = 0.10`, Q2/Q1, mixed — PASSES 100 s FLAT, AND THE WINDOWED TRACE IS THE EVIDENCE:**

```
 t= 40- 50  eta 0.1112  u 0.4023  NL 6        t= 80- 90  eta 0.1114  u 0.4112  NL 6
 t= 50- 60  eta 0.1112  u 0.4043  NL 6        t= 90-100  eta 0.1115  u 0.4127  NL 6
 t= 60- 70  eta 0.1112  u 0.4076  NL 6        t=100-110  eta 0.1116  u 0.4082  NL 6
 t= 70- 80  eta 0.1113  u 0.4097  NL 6
```

η flat in the **fourth decimal** over seven consecutive windows, `u` flat in the third, Newton
**constant at 6**, all well past the 36 s fill (rule 14). **This is a genuine post-fill stability
result, not a run that merely had not failed yet** — the distinction rule 12c exists to enforce.

**⛔ BUT `A = 0.15` DIVERGES, AND IT IS THE SAME MODE.** `c3v_amp_mixed`, identical in every other
knob, broke away in its last window:

```
 t= 60- 70  eta 0.1865  u 0.6949  NL  8
 t= 70- 80  eta 0.2113  u 0.7379  NL 10
 t= 80- 90  eta 0.2072  u 0.8514  NL  8
 t= 90-100  eta 0.2780  u 1.1443  NL 15   <- eta +34%, u +34%, Newton ~2x in one window
```

Onset ≈ **t = 85–90 s (53–56 periods)**, `η` rising *with* `u` and Newton degrading — the §5.2b
velocity-led signature, arriving late instead of not at all.

⚠ **THE CONCLUSION THIS FORCES, AND IT IS SHARPER THAN "THE PROJECTION WAS A CONTRIBUTOR": THE
MIXED FORMULATION DELAYS THE INSTABILITY, IT DOES NOT REMOVE IT.** At `A` = 0.10 the delay exceeds
the 100 s observation window and the run looks cured; at `A` = 0.15 the same construction buys
**≈ 16×** over its projected control (5.8 s → ~87 s) and then fails anyway. **A pass at one
amplitude is therefore not a pass** — and reading `c3v_base_mixed` alone, which is what a single-arm
test would have done, would have produced a confident and wrong "fixed" verdict. ⚠ Do not restate
the `base` result without the `amp` result beside it.

⚠ **This is an amplitude dependence at the `:full` tier and it is NOT old rule 13** (the
`A_wave ≤ 0.001` cap, lifted 2026-09-06 as an equal-order artefact). That was Taylor-Hood curing an
inf-sup failure; this is a distinct amplitude sensitivity in the Class-III path, measured **on**
Taylor-Hood, at 100× that cap. Both `amp` arms (Q2/Q1 and Q3/Q2) die, and both die earlier than
their `base` twins — so amplitude is a **second independent axis** alongside pairing and `dx`.

⚠ **`c3v_dt_mixed` and `c3q_dt_mixed` are both showing early warning** — Newton maxima creeping
4 → 6 and 4 → 12 respectively, with `c3q_dt_mixed`'s `u` at 0.434 → 0.707 in one window. Neither had
failed at capture; neither may be quoted as a pass.

⚠ **THESE FOUR ARMS WERE TERMINATED BY AN OPERATOR REBOOT AT 2026-09-23 15:38, NOT BY DIVERGENCE.**
A frozen GNOME session forced a restart; suspend was tried first and failed. Final states:
`c3v_base_mixed` **t = 101.4 s (63.4 T), η 0.1054, u 0.411, Newton 6 — still flat and healthy**;
`c3v_amp_mixed` t = 95.6 s (diverging, η 0.211, u 0.849, Newton 14); `c3v_dt_mixed` t = 59.2 s;
`c3q_dt_mixed` t = 19.8 s. **Do NOT read `c3v_base_mixed` stopping at 101.4 s as a failure** — its
diagnostics simply end where the machine went down. It never reached the 160 s / 100-period target,
so its result is a **lower bound**: ≥63.4 periods stable, not a completed 100-period run. Re-running
it is the cheapest way to convert that bound into the regression-gate reference trace §5.7 item 4
needs. All artifacts survived (`output/local_1d/<arm>/` + `_logs_2026-09-23/`).

#### What is established, and what is not

✅ **The projection was a large, real contributor.** Removing it buys 1.9×–>16× in every matched
pair, and turns the Q2/Q1 `A`=0.10 flat case from a 12.6 s death into **>100 s (62.9 periods) flat
in the fourth decimal with Newton constant at 6**, post-fill.
✅ **There is now one configuration in which the full eight-component `𝓝` runs long at production
amplitude.** That did not exist before this branch, and it is the only candidate reference trace for
a `:full` regression gate (§5.7 item 4).
⛔ **It is a DELAY, NOT A CURE.** The `A`=0.15 twin of that very run diverged at ~87 s with the same
velocity-led signature. **`:full` remains unstable; the mixed formulation moves the onset.**
⛔ **The projection was not the whole cause.** Q3/Q2 dies at every setting, `dx` still advances onset
in a projection-free formulation, and **amplitude is a third independent axis** — both `amp` arms
die, each earlier than its `base` twin.
⛔ **The `dx` explanation of §5.2c is refuted as stated** and has no replacement yet (H4).

⚠ **THE DECIDING CONTROL HAS NOT BEEN RUN — THIS IS A RULE 14c GAP AND IT INVALIDATES ANY
CLASS-III ATTRIBUTION AT Q3/Q2.** Every Q3/Q2 arm in this batch is `BALFEM_NL_PRESSURE=full`; the
only `:native` runs on record are Q2/Q1. So "Q3/Q2 `:full` dies at 16 s" **cannot presently
distinguish**:
* Class III is the carrier and Q3/Q2 is merely where it is finally represented → a `:native` Q3/Q2
  arm survives 100 periods; from
* Q3/Q2 is unstable here for a reason unrelated to Class III → `:native` dies too, and §5.7 item 0,
  §0c and this entire branch are aimed at the wrong object.
**A `:native` Q3/Q2 arm is the single highest-value run available and it is cheap. Run it before
drawing any further Class-III conclusion.**

⚠ **`c3v_amp_mixed` IS UNDER WATCH AND SHOULD NOT BE QUOTED AS A PASS.** Its `u_max` window maxima
run **0.695 → 0.738 → 0.851 → 0.951** over t = 60–100 with `η` bounded at 0.205–0.211 — the
velocity-led signature of §5.2b beginning. It may yet join the failures; at the time of writing it
has not.

⚠ **Cost note.** The Q3/Q2 mixed arms quoted ETAs of **186–242 h** for 100 periods at ~31.7k DOFs
with an AD-coupled 5-field Jacobian. None could have finished regardless of stability; any future
Q3/Q2 stability claim needs either hand Jacobians for the auxiliary rows or a shorter target.

⚠ **All eleven failures share one signature** — `x_at_max` snapping back and **pinning at the
inflow** (x = 0.50, 0.25, 0.12, 0.56) while `η` stays bounded and `u` runs away — i.e. the same
velocity-led generation-boundary mode as §5.2b and §5.6. The mixed formulation changed **when** it
happens, never **what** happens. (This is readable at all only because `build_run_diagnostics` was
made mixed-aware on this branch; before that, mixed failures had a time but no place —
`NEW_TREATMENT.md` F.3.1.)

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

### 5.6b The 1-D reproduction of that signature (2026-09-15)

⚠ **The velocity-led mode above is NOT exclusive to the cluster, to 2-D, or to a sea state.** §5.2b
reproduces it in a 1-D flume with a regular wave in 13 s of simulated time, and — unlike the cluster
batch — **with a `:native` control that survives 100 periods**. Every diverged run in §5.6 was
`:full`; the 1-D campaign shows `:full` failing on a FLAT bed with η bounded and `u` exploding at the
generation boundary, and shows an exact-AD Jacobian failing identically. **Re-running the cluster
suite should therefore carry a `:native` arm**, or it cannot distinguish the seed bug it was
re-specified for from the `:full` operator defect.

### 5.7 Still to check — ordered by what would change a conclusion

0. ⛔ **THE CLASS-III ASSEMBLY IS NO LONGER A SUFFICIENT EXPLANATION — AND THE DECIDING CONTROL IS
   ONE CHEAP RUN.** *(Rewritten 2026-09-23 on the mixed-formulation campaign, §5.2e.)*
   The operator is verified (§5.2d) and the Jacobian is exonerated (rule 17b), so the assembly was
   the last suspect standing. The mixed formulation (`𝖦 ≈ ∇𝖲` as a genuine unknown, which **never
   forms `∂²` of a `C⁰` field**) tested it directly and returned a **split verdict**:
   * ✅ the projection is a **large real contributor** — mixed outlives projected in **all eight**
     matched pairs, and Q2/Q1 flat `:full` goes from 12.6 s to **≥61 wave periods, flat and on the
     `:native` trace**. `:full` is therefore *not* unconditionally unstable;
   * ⛔ but it is **not the whole cause** — every **Q3/Q2** arm still dies, *earlier* than its Q2/Q1
     twin, and the **`dx` signature survives a formulation that cannot commit the recovery error**.
   ⚠ **Two previously-recorded next steps are now spent, one of them refuted.** Step (0) "run `:full`
   at Q3/Q2 — free" **was run and REFUTED with the sign reversed**: Q3/Q2 is uniformly worse, so the
   `∂²η ≡ 0` masking argument explains why Q2/Q1 looked better, not why Q3/Q2 fails. Step (1)
   **de-lag** is discharged by the mixed construction and did help (both `dt` arms improved). Step
   (2) the **algebraic reduction** is done and exact. ⚠ **Step (3) `C⁰`-IP was conditioned on "if the
   `dx` signature survives 0–2" — it has survived, but H4 of §5.2e removes its rationale**: `C⁰`-IP
   is a consistency fix for a broken `∂²`, and the surviving `dx` sign was measured in a formulation
   that has no `∂²` to break. **Do not launch `C⁰`-IP on the old justification.**
   ⚠ **And "Rejected: the mixed formulation (system size)" is superseded** — it was built anyway and
   is the source of every result above (`OPEN_ISSUES.md` §0c option 4 is re-scored accordingly).
   **THE ONE RUN TO DO NEXT, BEFORE ANY OTHER CLASS-III WORK:** a **`:native` arm at Q3/Q2**. Every
   Q3/Q2 arm so far is `:full` and the only `:native` runs on record are Q2/Q1, so the batch cannot
   tell "Class III is the carrier" from "Q3/Q2 is unstable here for an unrelated reason" (rule 14c).
   Second: **`𝖦` one order below `u`** — the aux space is currently `reffe_U`, i.e. equal order with
   the velocity (`src/horizontal.jl:59`), and that pairing has never been analysed.
1. ⛔ **The nonlinear `p_η` order reduction** (`OPEN_ISSUES.md` §0b). Next step is **derivation, not
   another run**: check the advection block term-by-term against `BALFEM_models/`. Cheap
   discriminators if wanted: amplitude sweep (`a_eta = 0.8 → 0.4 → 0.2`), and Q4/Q3 extended to
   `nx = 64`.
2. 🔴 **Derive the expected orders for this system.** Until that exists, "sub-optimal" in every
   convergence table is measured against an assumed optimum, and §5.3 cannot be resolved by runs.
3. 🔴 **Re-run the cluster suite** with the seed fix, the new directional geometry and the corrected
   sponges. Nothing on record post-dates those fixes.
4. 🔴 **A long-duration nonlinear regression test.** The mode needed 50–80 s to emerge and no test
   runs that long; this is the gap that let it reach production. ⚠ **PARTLY DISCHARGED 2026-09-15:**
   §5.2b provides the reference trace — `:native` flat, 100 periods, η 0.10494 → 0.10300, Newton
   5.12 — but it is a *run*, not a *gate*. It is also the duration the bar ladder needed: the 1:2
   bar looked stable at t = 80 s and died at t = 137 s, so a regression built at 80 s would have
   passed it. **Any gate written from this must run to at least 100 periods.**
5. ⛔ `:full` (models 7–8) is **not MMS-verifiable as built** — the `{1,2,4,5}` blocks are in the
   residual but absent from `jacobian_u`. ⚠ **THE SECOND SENTENCE OF THIS ITEM WAS WRONG AND IS
   WITHDRAWN.** It read "closing it means completing `jacobian_u`, not refining a mesh", which
   assumed the incomplete Jacobian was also what makes `:full` diverge. **Rule 17b now proves the
   opposite across the whole factorial**: hand and exact-AD Jacobians crash in the SAME PLACE and
   differ only in iteration count. Completing `jacobian_u` is needed for MMS verification and buys
   **nothing** for stability. ⚠ And the operator is not at fault either — §5.2d verifies `𝓝`,
   Class III included, against the governing equations. **What is left is the Class-III ASSEMBLY**
   (frozen `L²` projections), which is what §5.2c's `dx`/`dt` signs point at.
   Promoted 🟠 → ⛔: it is the tier every diverged cluster run used, and it now has a 13-second
   reproduction with a passing control.
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
1b. **THE `C⁰` DISTRIBUTIONAL OBJECTION APPLIES TO THE *TRIAL* SPACE, NOT TO THE TEST SPACE.**
    ⚠ THE MODEL CARRIES `∂³u`; THE WEAK FORM CARRIES `∂²u`. Measured 2026-09-23 by a
    `ε(x−x₀)^k` probe: `p_nh` depends on `∂²u` and not `∂³u`; the strong-form momentum residual —
    ours and Yang & Liu's ⌊2.30⌋ alike — depends on `∂³u`. `R_P`'s integration by parts moves one
    derivative onto the test function, which is what leaves the trial field needing only a Hessian.
    Do not "correct" either count into the other; they describe different objects.
    The distributional `∂²` of a `C⁰` **trial** field is `{∂²u}` cell-wise plus a **Dirac layer on the
    skeleton** weighted by the jump `[∂ₙu]`; cell quadrature sees only the first half, so typing
    `∇∇(u)` into an integrand silently drops the layer. That is why the Class-III components
    `{1,2,4,5}` are carried by frozen `L²` projections (`OPEN_ISSUES.md` §0c).
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

* **`BALFEM_models/` is the single source of mathematical truth.** `MODEL.md` mirrors its notation.
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
