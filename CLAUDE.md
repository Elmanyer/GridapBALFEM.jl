# CLAUDE.md — `GridapBALFEM.jl/` · the 2D BALFE-M algebraic wave solver

> ## ⇨ START HERE
>
> This file is the map and the standing rules. The detail lives in ten focused documents in
> `building_files/` (**gitignored — not under version control**):
>
> | document | what it answers |
> |---|---|
> | [`MODEL.md`](building_files/MODEL.md) | the maths: σ-tensors, the global residual term by term, nonlinear pressure, Jacobians |
> | [`ARCHITECTURE.md`](building_files/ARCHITECTURE.md) | code structure and workflow: the stacked layout, `src/` map, FE spaces, time loops, distributed path, diagnostics |
> | [`VERIFICATION.md`](building_files/VERIFICATION.md) | what is *proven*: the MMS campaign, verified scope, the Jacobian oracle, scope boundaries |
> | [`TEST_SUITE.md`](building_files/TEST_SUITE.md) | every test, its score, and **what it cannot detect** |
> | [`CONFIGURATION.md`](building_files/CONFIGURATION.md) | the settings, the evidence behind each, measured performance |
> | [`WAVE_GENERATION.md`](building_files/WAVE_GENERATION.md) | sources, Dirichlet generation, sponge, WaveSpec coupling |
> | [`RUNNING.md`](building_files/RUNNING.md) | how to launch: local, cluster, sysimage, env vars |
> | [`OPEN_ITEMS.md`](building_files/OPEN_ITEMS.md) | open work, each with its decisive next step |
> | [`PENDING_TASKS.md`](building_files/PENDING_TASKS.md) | studies designed but not begun, with their blocking prerequisites |
> | [`MMS_VBASIS_CAMPAIGN.md`](building_files/MMS_VBASIS_CAMPAIGN.md) | the vertical-basis campaign: the `(M,p)` × 8-model convergence matrix (its Phase-1 minimax σ-meshes are superseded by `src/vopt.jl`) |
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
> were derived by multi-property **minimax and are SUPERSEDED** — `src/vopt.jl` now implements the
> Yang & Liu total-relative-error functional instead (§3), a different objective, so the `p ≥ 2`
> node sets have moved (P2LFE-2: `c₁` 0.8794 → 0.8298). **Headline: `η` reaches 2.999-3.000 and `u` reaches 4.000 on
> EVERY basis and EVERY model in the spatial matrix (30/30) — the order of accuracy is independent
> of the vertical basis, which is the direct evidence for the basis-agnosticism the model is named
> for.** Full results, method and every correction: `building_files/MMS_VBASIS_CAMPAIGN.md`; single
> data file `output/local/mms_campaign/campaign_results.csv`.

---

## 0. What this folder is

The production home of the 2D BALFE-M solver: the self-contained serial + distributed package
`src/GridapBALFEM.jl` (stacked `[η,𝖴x,𝖴y]` layout, loop-free residual), its test suite (`test/`),
sequential and cluster examples (`examples/`), the standalone postprocessing library
(`postprocessing/GridapBALFEMPost`), the vendored sea-state package (`WaveSpec.jl/`), and the design
and paper material in `building_files/`.

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
| `src/` | the solver package — **17 files** (`vopt.jl` added 2026-09-01), mapped in `ARCHITECTURE.md` §2 |
| `test/` | 28 test files + `runtests.jl` + `test/cluster/` + `test/local/` — inventory and scores in `TEST_SUITE.md`. `test_mms_distributed_parity.jl` (4 ranks) gates the distributed MMS path against the sequential one |
| `examples/` | 7 sequential + `distributed/` (7 cluster scripts + `_dist_common.jl`), `distributed_small/` (5 parametric), `validation/` (7), `local_1d/`, `local_2d/`, `local_mms/` (**11** — the parametric MMS studies, the vertical-basis campaign (`run_vbasis_campaign.jl`, `run_vbasis_shard.jl`, `report_vbasis_campaign.jl`), and the Phase-B batch: `run_phaseB_shard.jl` + `supervise_phaseB2.sh`, a resume-capable runner with a memory-aware supervisor), `inspect_run.jl` — `RUNNING.md` §2 |
| `run/` | 10 production SLURM launchers + `run/dist_small/` (**30** small-domain 2-D cases — the 7 superseded 1-D twins were deleted 2026-09-02) + `run/local/` (30 scripts incl. **6** 1-D cases + `run_all_1d.sh`), all through `run/balfem_env.sh` (cluster) or `balfem_local.sh` (workstation) — `RUNNING.md` §3–4. ⚠ Every launcher is sized to **≤48 cores/node and ≤224 GB/node**: a rome node advertises 256 GB but SLURM allocates only ~224, so `64 × 4 GB` is refused at submit time |
| `compile/` | the cluster sysimage build chain — `RUNNING.md` §5 |
| `postprocessing/` | `GridapBALFEMPost` — self-contained, own environment, **no dependency on the solver** |
| `WaveSpec.jl/` | vendored stochastic sea-state synthesis (CMOE-TUDelft). Tracks the **GitHub repository version, not a tagged release** — the release's `change_seed!` is broken |
| `Gridap.jl/` | the vendored **fork** (`Elmanyer/Gridap.jl` @ `fix-transient-multifield-ad`, one commit on `v0.20.8`) making transient-multifield AD work — `CONFIGURATION.md` §2 |
| `BALFEM_models/` + `.zip` | **THE CURRENT LaTeX project** — the authoritative derivation (§2) |
| `LFEM_discretisation/` + `.zip` | **SUPERSEDED** pre-rename LaTeX, kept for provenance. Do not edit |
| `CFC2027_LFEMultilayer_abstract/` | conference abstract (CFC 2027); its class needs the `newtx` fonts |
| `building_files/` | the seven documents above. **Gitignored** |
| `algebraic_balfem2D.jl`, `test_algebraic_balfem2D.jl`, `test_2HDmodel.ipynb` | the single-file prototype and an early notebook, kept for reference; `src/` is the maintained form |
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

**Vertical grid optimisation** (`src/vopt.jl`, added 2026-09-01) — the Yang & Liu **total
relative-error functional**, replacing the earlier minimax objective:
`E_total = Ē_c + Ē_cg + Ē_shoal + Ē_u + Ē_w`, each term integrated over `kd` against
`W = exp[(2^−kd − 2^−π) log 5]` and normalised by its **median over the sampled design population**
(the medians are the trade-off weights, so the argmin depends on them).
* `SigmaBasis` is a **self-contained analytic** piecewise-Lagrange σ-basis. Required because `E_u`
  and `E_w` put an absolute value *inside* the σ-integral, so no Gram tensor can absorb them and
  `φ_j(σ)` must be evaluated pointwise in the optimiser's inner loop — and building it here keeps
  Gridap's DOF numbering off the load path. Validated against `assemble_dispersion_tensors` through
  the permutation-**invariant** `R(μ)`: agreement **2e-15** across `p=1,2`, `M=1..4`.
* **The design band is a calibrated fixed point, not a knob.** `W` saturates at 0.834 rather than
  decaying, so the `kd`-integral does not converge and `Ω` is a real design band; a single `Ω`
  reproduces `M=2` but misses `M=3,4` by 0.20/0.33. What is invariant is
  `κ = Ω·Δσ_top/p = 2.06, 2.24, 2.20` → **`VOPT_KAPPA = 2.17`**, solved by bisection (`vopt_Omega`).
* **Validation — reproduces Yang & Liu Table 1 to ≤0.014** (`M=2,3,4`); fitting `Ω` directly gives
  **0.002**, which bounds the band closure's share of the error. `DEFAULT_CBDY` (`p=1`) is therefore
  **kept at the published values**; `DEFAULT_CBDY_P` holds the new `p ≥ 2` optima (P2LFE-2 =
  `[0, 0.8298, 1]`), and `resolve_cbdy(M, c_bdy, p=1)` is `p`-aware via a trailing positional.

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
carries a mesh-independent velocity-error floor (`VERIFICATION.md` §4).

---

## 5. Current Implementation Stage

**Working.** Feature-complete in both serial and distributed forms: the stacked loop-free residual +
hand Jacobians, the full nonlinear physics, the SDIRK/θ integrators, all boundary treatments, and
Dirichlet boundary wave generation with WaveSpec coupling.

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
   3.33 → 3.45 → 3.78 → **3.94** to `nx=128` on three models. Answers `OPEN_ITEMS.md` §6 for this case.
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

**Open work** — nothing is half-built in the solver; the open items are verification gaps,
performance, and follow-through. Full list with decisive next steps: `OPEN_ITEMS.md`; studies
designed but not begun: `PENDING_TASKS.md`.

**🔄 PHASE-B MMS BATCH IN FLIGHT** (`output/local/mms_phaseB/`, 25/40 at 2026-09-02). Four tasks:
T7 P1LFE-4 on an extended ladder, T8 the `:full` pair at `a_eta=0.4` with the projections now
assembled, T9 tier 2 (the horizontal pairings), T10 the `p ≥ 2` bases on the new σ-meshes.
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
  present. This unblocks the vertical-basis convergence study (`PENDING_TASKS.md` §1 prerequisite A);
  the study driver is `examples/local_mms/run_vertical_basis_study.jl`.
* ✅ **the MMS `:full` omission is FIXED** (found 2026-08-21, fixed 2026-09-01). Both MMS drivers now
  build the `nlp` context, and both time loops **prime it from the IC** — `update_nlp_state!` runs
  only *after* an accepted step, so step 1 assembled as if from rest (exact for a rest start, wrong
  for `u0 = u*(t0)`). `e_u` fell **1.19e-03 → 1.37e-06** and `p_u` went **−0.00 → 1.948**, so the
  old floor measured **omission**, not the frozen-projection lag `VERIFICATION.md` §4 describes.
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
  (`PENDING_TASKS.md` §1), tier 3's `:full`-floor-vs-`Nσ` question, and `OPEN_ITEMS.md` §6's
  velocity-shortfall question for the vertical-basis case (pre-asymptotic — the horizontal `Q2/Q1`
  and `Q4/Q3` pairings are still untested on an extended ladder, and the recipe is now cheap)
* naming follow-through outside this checkout: GitHub repo (still `GridapLFEM.jl`; the remote URL
  is stale but redirects), cluster checkout, sysimage rebuild

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
   **Re-run the mechanical audit after adding any nested helper.** (`VERIFICATION.md` §7)
   *Third instance, 2026-08-21, `wave_properties`:* two closures each called their phase speed `C`,
   which was also the enclosing function's return value, so the last γ evaluation — the **Airy** one
   — overwrote it. **Note the shape of the symptom: `|C/Ce − 1|` collapsed to ~1e-12 for every mesh
   at every `kd`, i.e. the model looked PERFECT exactly where it is worst, and `C_g` and `γ`
   reproduced their published values throughout.** A capture bug can be invisible in every channel
   but one, and the channel it corrupts can fail in the flattering direction.

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
13. **`A_wave ≤ 0.001 m`** for stable long fully-nonlinear integrations.
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
    there. This also **unseats the "flat memory ⇒ H3 not supported" argument** in `OPEN_ITEMS.md` §1,
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

### Testing and measurement

28. **"The suite passes" is not "the model is verified."** Most of the suite is
    **self-consistency**: writing `R = R_true + E`, the error appears on both sides and cancels
    identically, so such a test passes for **any** residual. Only the analytic MMS — whose forcing
    never touches `problem.jl` — can detect a self-consistently wrong residual. (`VERIFICATION.md` §0)
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
* **`building_files/` is gitignored** — nothing in it is version-controlled, so anything load-bearing
  belongs in the tracked `CLAUDE.md`, `README.md` or test files instead.
* **Run Julia via the `julia-mcp` tool.**
