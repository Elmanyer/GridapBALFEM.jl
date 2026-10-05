# BROKEN_FORMULATION_PLAN.md — the broken (skeleton) Class-III formulation and the C⁰-IP stabiliser

> **⚙ v1 RECORD (2026-10-05).** This file documents the broken formulation and first-order C⁰-IP on v1 — the broken formulation is now v2's ONLY Class-III treatment, but every campaign here used v1 options (`broken=true` opt-in, `C3_MASK=gs`). It is kept as the record of how v1 reached its results;
> the v2 solver (branch `v2-solver`) differs — one Class-III formulation (broken), `nl_pressure::Bool`
> with all eight `𝓝` components, no component mask (`V2_SOLVER_PLAN.md`). Paths `output/…` here are
> now under `output_v1/`.

*Branch `broken-formulation-solver`, from `main` @ `58ef9bb` (2026-09-30). The mathematics is
LaTeX §6.4 "Broken Weak Formulation: Term-by-Term Audit"
(`latex_docs/BALFEM_models/NumericalImplementation/BrokenAudit.tex`); this file is the
implementation and campaign record. Companion records: `NEW_TREATMENT.md` (projected and mixed
treatments), `OPEN_ISSUES.md` §0e (the instability), `PLANNED_CAMPAIGNS.md` §6c (stabilisation
campaign).*

---

## 0. What is being built, and why

The audit's conclusions fix the scope:

1. **The Galerkin residual is consistent.** No test-derivative-jump (kind T) facet term is missing from it — jumps of the C⁰ test functions themselves are identically zero (kind 0).
   Every block except two is complete in the broken sense.
2. **The two incomplete blocks are Class III in the surface-slope (𝓚) and leading-pressure (𝓟)
   blocks.** Their `∇𝖲`, `∇𝖻` carry a skeleton layer `−Σ_F ⟦·⟧ δ_F`. The projected path smooths it
   away; the mixed path contains it only after an L² projection. **The broken formulation assembles
   it exactly:** cellwise Hessians plus one skeleton integral, with no auxiliary unknowns.
3. **The leading-pressure skeleton term `⟦𝒫ᵢ⟧·vᵢ` is excluded.** It is consistent, but it degrades the
   effective mass matrix from `h⁻²` to `h⁻⁴` conditioning and destroys convergence (Table 6.1 of
   the LaTeX).
4. **Stability is sought from the C⁰-IP penalty** `J_h`, the only sign-definite skeleton term. It is
   independent of the Class-III treatment and also applies to `:native`.

It is delivered as two **orthogonal** options, so that each can be compared against what exists:

| option | replaces / adds | knobs | combines with |
|---|---|---|---|
| `broken=true` | the Class-III 𝓚/𝓟 treatment (instead of `projected` or `mixed=true`) | — | `nl_pressure=:full` only |
| C⁰-IP penalty | adds `J_h` to the residual and the exact block to `∂R/∂u` | `cip_gamma_u`, `cip_gamma_eta`, `cip_hexp` | any tier (`:none`, `:native`, `:full`); any Class-III treatment |

The `(broken, γ)` pair spans the comparison the campaign needs: broken alone (does distributional
completeness change stability?), CIP alone on `:native` (does the penalty remove the
Class-III-free instability?), and broken + CIP (the culminating C⁰-IP weak form).

## 1. Design decisions (recorded before coding)

* **D1 — Layout stays `[η, 𝖴x, 𝖴y]`.** No auxiliary fields, so the hand Jacobians, the
  distributed layout and every driver index are unchanged.
* **D2 — Form (A), not form (B).** Form (A) needs only cellwise *trial* Hessians and first
  *test* derivatives. Form (B) would need cellwise test Hessians plus a mandatory kind-T facet
  term. They are the same operator (LaTeX eq. "form A" ≡ "form B" − ∮).
* **D3 — Cellwise ∂ₐ𝖲 and ∂ₐ𝖻 are hand-expanded from `∇∇`** of the raw fields (`∇∇(𝖴)` is a
  `ThirdOrderTensorValue{2,2,Nσ}`, `∇∇(η)` a `TensorValue{2,2}`). `∇` of an `Operation`-composed
  field is not implemented (rule 6; re-confirmed 2026-09-30). The bed Hessian is analytic
  (`alg_bed_hessian`).
* **D4 — The volume part reuses the reduced contributors** `nlp_gradH_reduced_contrib` and
  `nlp_P_reduced_contrib` verbatim, with `GU → GU_𝒯` and `N4 → N4_𝒯`. The skeleton part uses the
  same weights `WK3/WP3/K3[4]/P3[4]` and the same `c3_mask` logic, so the three treatments differ
  **only** in how `∇𝖲` and `∇𝖻` are obtained.
* **D5 — The facet test factor is averaged:** `{∂ₐH}`, `{D_W}` and `{Wₐ}`. `Wₐ` is continuous, so
  its mean is its value. `H` and `𝖴_n` are continuous and taken from the `+` side.
  `n_F = nΛ.plus` and `⟦f⟧_F = f⁺ − f⁻`.
* **D6 — The broken facet terms are quasi-Newton in `jacobian_u`,** exactly like every other
  Class-III block (rule 17b: the Jacobian sets the path, not the root). `use_ad=true` gives the
  exact Jacobian; AD through skeleton integrals was verified on 2026-09-30 (AD vs FD 7e-11).
* **D7 — The penalty uses still-water scales** `τ_u = d√(g d)`, `τ_η = √(g d)` with `d = h(x)`, and
  weights the layer contraction with `Mv`. It is therefore **linear**, with an **exact** Jacobian
  block and nothing in `∂R/∂u̇`. `h_F` is the facet measure, and the exponent `cip_hexp` defaults to 2.
* **D8 — Skeleton objects are built once** (`attach_skeleton!`) and stored on
  `prob.skel::RefValue{Any}`, following the `nlp_state`/`nlp_ctx` pattern. With `nothing`, every
  path is bit-identical to `main`.
* **D9 — Sequential only in this branch,** like the mixed path. On an x-periodic mesh the
  identified edge is an interior facet: the skeleton has 16 facets on 16 cells, total length Lx
  (checked 2026-09-30).
* **D10 — `broken` is refused with `mixed=true` and with `nlp_inloop=true`.** It is a different
  treatment of the same blocks, not an addend.

## 2. Tasks (sequential)

### T1 — Skeleton infrastructure (`src/broken.jl`)
* `alg_hess(u, a, b)`: the `(a,b)` Hessian component of a scalar or stacked field, via `∇∇`.
* `build_skeleton_ctx(model; degree)`: returns `(Λ, dΛ, nΛ, nx, ny, hF)`, with `hF` from
  `get_cell_measure(Λ)`.
* `attach_skeleton!(prob, model; broken, cip_gamma_u, cip_gamma_eta, cip_hexp, degree)`: stores
  the NamedTuple on `prob.skel[]`.
* `BALFEMProblem` gets one field, `skel::Base.RefValue{Any}` (`nothing` by default); its single
  constructor call is in `build_problem_raw`.

### T2 — Cellwise Class-III fields
`broken_class3_cell_fields(prob, d_cf, η, H, dHx, dHy, Ux, Uy, DU, S)` returns `(GU_𝒯, SD, N4_𝒯)`
from

    ∂ₐ𝖲 = ∂ₐH·D + H·∂ₐD + ∂ₐ𝖻,      ∂ₐD = ∂ₐ∂ₓ𝖴x + ∂ₐ∂ᵧ𝖴y
    ∂ₐ𝖻 = (∂ₐ∂ₓH)𝖴x + ∂ₓH·∂ₐ𝖴x + (∂ₐ∂ᵧH)𝖴y + ∂ᵧH·∂ₐ𝖴y,   ∂ₐ∂_bH = ∂ₐ∂_b h + ∂ₐ∂_b η

### T3 — Skeleton layer of the 𝓚 and 𝓟 blocks
`broken_class3_skeleton_contrib(prob, sk, H, dHx, dHy, Ux, Uy, S, bf, Wx, Wy, DW)` assembles the
two facet integrals of LaTeX eq. (6.35), with the `c3_mask` bits applied identically to
`_c3_sum`.

### T4 — C⁰-IP penalty
* `cip_contrib(prob, sk, η, Ux, Uy, q, Wx, Wy)`: the residual term.
* `cip_jacobian(prob, sk, dη, dUx, dUy, q, Wx, Wy)`: the exact Jacobian block. It is the same
  function, because `J_h` is linear.
* Wired into `global_residual` whenever `prob.skel[]` has γ > 0, and into `jacobian_u`.

### T5 — Wiring
* **`global_residual`, Class-III branch:** `aux !== nothing` → mixed; `sk.broken` → broken
  (T2 + T3); otherwise projected.
* **`setup_and_run`:**
  * new kwargs `broken`, `cip_gamma_u`, `cip_gamma_eta`, `cip_hexp`;
  * the guards of D10;
  * no projection context when `broken`;
  * banner lines naming the Class-III treatment and the penalty;
  * the skeleton attached after `build_problem`.
* **Periodic driver:** `BALFEM_BROKEN`, `BALFEM_CIP_GU`, `BALFEM_CIP_GE`, `BALFEM_CIP_HEXP`, and
  output-name extras `broken`, `cipu<γ>` and `cipe<γ>`.
* **Exports** in `src/GridapBALFEM.jl`.

### T6 — Unit gates (`test/test_broken_formulation.jl`, registered `:medium`)
| gate | what it proves | oracle |
|---|---|---|
| G1 layer identity | cellwise ∂ₐ𝖲 (hand Hessian expansion) − facet layer ≡ the mixed weak gradient `−∫𝖲∂ₐψ + ∮𝖲ψnₐ`, both directions, for 𝖲 and 𝖻, on a sloping bed with a walled y-boundary | exact integration-by-parts identity, to round-off |
| G2 layer live | the broken Class-III facet contribution is non-zero and of the same order as the volume part on a rough state (rule 38d) | — |
| G3 no leakage | `broken=true` with the penalty off leaves `:none`/`:native` residuals **bitwise** unchanged, and `skel = nothing` reproduces `main` | bitwise |
| G4 broken ≈ mixed | at a smooth state the broken Class-III residual agrees with the mixed one to the projection error, and both differ from projected by the lag/layer | assembled vectors |
| G5 CIP algebra | symmetric PSD; zero on constants and on a global polynomial of degree ≤ p (non-periodic mesh); the η-penalty conserves mass (`q ≡ 1`) | exact |
| G6 CIP selective | the Rayleigh quotient on a λ = 2dx mode exceeds the resolved-carrier value by > 10⁴ | — |
| G7 Jacobian | the hand `jacobian_u` with CIP equals the FD oracle on the CIP block (rel < 1e-6), and the gap on the broken facet block shrinks with amplitude (the D6 omission) | central FD of the residual |
| G8 smoke | a 20-step periodic box run with `broken + CIP` converges with Newton ≤ 10 | — |

### T7 — Stability campaign in the closed x-periodic box (Crank–Nicolson)
Template: `run/local/run_1dper_batch_cn*.sh`. The launcher is new,
`run/local/run_1dper_batch_broken.sh`. Stdout goes to `output/local_1d/periodic/_logs/`
(rule 47). Each `:full` arm has a `:native` twin in the same batch (rule 14c).
`P1LFE-2`, `A` = 0.10, `dt` = 0.04, `c3_mask = gs` (the mixed campaign's setting), analysis with
`periodic_growth.jl`.

**Screening batch (S1), 40 periods (64 s).** The two cells where the unstabilised solver fails
fastest under CN:

| case | pairing, cells/λ | tier / treatment | γ_u | control on record (CN) |
|---|---|---|---|---|
| `p32_broken_n16_cn` | Q3/Q2, 16 | `:full` broken | 0 | mixed diverged 26 s |
| `p32_broken_cip{1e-2,3e-2,1e-1}_n16_cn` | Q3/Q2, 16 | `:full` broken | ladder | — |
| `p32_native_cip3e-2_n16_cn` | Q3/Q2, 16 | `:native` | 3e-2 | `:native` grew +0.10 s⁻¹ |
| `p21_broken_n32_cn` | Q2/Q1, 32 | `:full` broken | 0 | mixed diverged ≈ 90 s |
| `p21_broken_cip3e-2_n32_cn` | Q2/Q1, 32 | `:full` broken | 3e-2 | — |

> ⚠ **SUPERSEDED BEFORE LAUNCH (2026-09-30), by the linear eigen-analysis in §4.2.** The damping
> estimate below treats `J_h` as a smooth, wavenumber-graded damper. It is not one: it acts as a
> **C¹ constraint** and leaves the oscillatory modes that fit in the C¹ subspace almost undamped.
> The arms actually launched test the transition (γ = 0.03) and constraint (γ = 0.3) regimes
> with **both** fields penalised; see `run/local/run_1dper_batch_broken.sh`.

The γ range comes from the grid-scale damping estimate `σ ≈ γ√(gd)/h ≈ 23γ s⁻¹` at `h` = 0.25 m.
It must beat the measured +0.48 s⁻¹ (Q3/Q2 `:full`, CN), so γ ≳ 0.02 up to an O(1) constant.
The carrier estimate `σ_c ≈ γ√(gd)k(kh)^(2p+1)` is below 1e-3 s⁻¹ throughout.

**Acceptance** (`OPEN_ISSUES.md` §0e):
1. the mid/high-band growth rate is ≤ 0 over the run;
2. the carrier energy decay is ≤ 1e-3 s⁻¹ (the CN `:native` reference: |σ_E| < 4e-5 s⁻¹);
3. MMS orders are preserved (T8);
4. the amplitude ceiling is measured (Stage 2, `A` = 0.15).

**Extension batch (S2):** the surviving γ at 100 periods, plus the Q2/Q1 64-cell `:native` case
(diverged at 146 s under CN).

### T8 — Consistency check of the penalty (MMS)
Plumb `cip_gamma_u` into `run_mms_case`, and run the linear and `:native` 1-D Q3/Q2 spatial ladder
at the chosen γ against γ = 0: the orders must be unchanged. Deferred until S1 selects γ.

### T9 — Documentation
* `CLAUDE.md` §5 status line;
* `OPEN_ISSUES.md` §0e;
* this file's results section;
* LaTeX chapter 8, once results exist.

## 3. Risks and how each is caught
* **Sign or normal convention of the layer:** G1 is an exact identity against the independently
  verified mixed definition.
* **Hessian component order** (`[a,b,j]`): G1 again. The identity fails at O(1) if `a`/`b` or the
  slot order is swapped. An asymmetric state is used (rule: asymmetric test states).
* **A dead knob** (rule 38d): G2 and G3, and every campaign arm prints its knobs in the banner and
  the output name.
* **The penalty masking the carrier:** acceptance criterion 2 is measured, not assumed.
* **Cost:** a skeleton integral adds one facet loop. The 3-field system is the projected path's
  size, so broken arms should be cheaper than mixed per step.

## 4. Results

### 4.1 Unit gates (T6), measured 2026-09-30

| gate | result |
|---|---|
| G1 layer identity | cellwise Hessians + skeleton layer vs the mixed weak gradient: **≤ 6.7e-12** (smooth) and **≤ 3.8e-13** (rough), for 𝖲 and 𝖻 in x and y. The cellwise part alone is off by 9e-4 to 1.4e-2 (smooth) and **0.50–0.84 (rough)**: the layer is load-bearing exactly where grid-scale content lives |
| G2 layer live | on a rough state the skeleton integral is **65 %** of the broken Class-III 𝓚+𝓟 contribution |
| G3 no leakage | an all-off attach leaves `skel = nothing`; `:native` is bitwise unchanged; broken without `:full` is refused |
| G4 broken vs projected | at a smooth state (projection at the same state, no lag), the relative difference of the Class-III functional is 1.39e-5 → 1.51e-6 → 1.71e-7 at nx = 6, 12, 24: **order ≈ 3.1**. Two consistent discretisations converging to each other |
| G5 C⁰-IP algebra | asymmetry 8e-17; λ_min −3e-13 (PSD); **1.9e-16** on a global cubic/quadratic polynomial; the η-penalty conserves mass (1.7e-13) |
| G6 selectivity | Rayleigh quotient λ = 2dx / carrier = **2.6e6** (Q3/Q2, 16 cells/λ) |
| G7 Jacobian | linear + C⁰-IP: hand vs FD **3.1e-11** (exact). `:full` broken + C⁰-IP: the quasi-Newton gap is 0.21 / 0.11 / 0.054 at a = 1, ½, ¼, i.e. **order 1 in amplitude** (the deliberate omission, D6) |
| smoke | Q3/Q2 box, broken + C⁰-IP (0.3, 0.3), CN, 40 steps: Newton 4.7/step, mass drift ~1e-17, 1.8 s/step |

### 4.2 ⚠ What the penalty actually does: a linear eigen-analysis of the box operator

The linearised flat-bed operator `−M_eff⁻¹(K + J_h)` was diagonalised in the Q3/Q2 box with
16 cells/λ (960 unknowns, 95 oscillatory wave pairs, carrier ω = 3.913 against 2π/1.6 = 3.927).

| penalty | oscillatory modes | their damping (s⁻¹) | overdamped (non-oscillatory) modes |
|---|---|---|---|
| none | 190 | 0 (round-off) | — |
| γ_u = 0.01 | 190 | ≤ **0.016** | 288, rates 0.036–137 |
| γ_u = 0.1 | 190 | ≤ 0.003 | 288, rates 0.36–1.4e3 |
| γ_u = 1 | 190 | ≤ 0.0003 | 288, rates 3.6–1.4e4 |
| γ_η = 0.01 | 112 | up to 6.9 (18 modes > 0.5) | transition |
| γ_η = 0.1 or 1 | **94** | ≤ 5e-5 | half of the wave modes |
| γ_u = 0.1, γ_η = 1 | 94 | ≤ 1e-5 | — |

**Reading.** As a quadratic form, `J_h` is grid-scale selective (G6). As a *dynamical* damper it
is not:
* the oscillatory modes re-arrange into the subspace where the penalised jumps vanish (C¹ in the
  penalised field) and stay undamped, with damping falling like 1/γ;
* the velocity penalty overdamps the non-oscillatory velocity modes outside C¹;
* the surface penalty removes the wave modes above the **cell** Nyquist (190 → 94). These are the
  ones only C⁰ quadratics can represent.

**The penalty limit is therefore a C¹ constraint, not a band-limited damping.** The ≥ 0.5 s⁻¹
acceptance criterion cannot be met by damping oscillatory modes in the C¹ subspace. What the
penalty *can* do is remove the element-scale half of the band. That is where the spectral-reach
analysis (LaTeX ch. 8) puts the C⁰ excess over Yang & Liu's FD. The carrier is untouched
(damping < 2e-7 s⁻¹ in every case).

Consequence for T7: the screening tests whether the growing modes live in the removed subspace,
in the constraint regime (γ = 0.3) and the transition regime (γ = 0.03), with both fields
penalised. The smooth-damping estimate `σ ≈ 23γ` was wrong and is withdrawn. The LaTeX property
(iv) was corrected in the same edit.

### 4.3 Screening S1: closed x-periodic box, Crank–Nicolson (INTERIM, 2026-09-30 16:50)

Q3/Q2, 16 cells/λ, P1LFE-2, `A` = 0.10, `dt` = 0.04, seed 1e-8, `c3_mask = gs`. Launcher
`run/local/run_1dper_batch_broken.sh`; logs `output/local_1d/periodic/_logs/`. Rates are energy
growth rates σ_E (s⁻¹) from `periodic_growth.jl`, with E(end)/E(start) over the fit window.
**Read the ratios with the rates:** a positive fitted slope with a ratio < 1 is carrier-phase
breathing, not growth.

**`:full`, t = 10–24 s** (mixed on record, CN: diverged at 26 s, +0.56/+0.57/+0.59 over 2–9 s):

| arm | η mid / high / sub | ratio | u mid / high / sub | ratio | Newton | verdict |
|---|---|---|---|---|---|---|
| broken, γ = 0 | +0.24 / +0.19 / +0.19 | ×16–23 | +0.21 / +0.20 / +0.22 | ×5–15 | rising, 10 at 24 s | **unstable** (as the audit predicted: completeness ≠ stability) |
| broken + C⁰-IP (0.3, 0.3) | +0.03 / +0.01 / +0.01 | ×0.07–0.52 | +0.05 / +0.04 / +0.04 | ×0.11–0.97 | flat, 5 | **bounded past 26 s**; high bands breathe at 1e-9 |
| broken + C⁰-IP (0.03, 0.03) | (to 14.6 s) | ×0.36–1.08 | | ×0.45–0.63 | 6 | no net growth yet |

**`:native`, t = 10–24 s:** with C⁰-IP (0.3, 0.3) the η bands stay at the seed level
(3–5e-10) to t = 45 s, against the control's rise to 1e-7. The surface-velocity high band still
rises (σ fit +0.099), but **sub-linearly and decelerating**: its increments per 5 s shrink from
4.0e-8 to 2.1e-8 and √E/t falls monotonically. ⚠ The unpenalised control's exponential phase
only starts at t ≈ 90 s (mid band ×3–4 per 10 s over 90–150 s), so this arm runs 100 periods
(`p32_native_cipu0.3e0.3_n16_cn_p100`).

**Operational notes.**
* The same-batch `:native` control reproduced the main-branch run to every printed digit
  (η band-4 energy 3.58e-08 at 10.2 s, 9.04e-08 at 20.2 s) — a live check of G3. It was stopped
  as a duplicate.
* Cost: broken Q3/Q2 box ≈ 15–24 s per step (`:native` ≈ 6–8 s). Gridap assembly dominates, as
  before.
* The 40-period window is too short for `:native`, whose instability emerges after ≈ 90 s.

**MMS (T8), Q3/Q2 spatial ladder:** linear p_η/p_u = 2.995/3.770 without the penalty and
**3.089/3.771 with γ_u = γ_η = 0.3**, so the orders are preserved. `:native` pending.

**Still running at the time of writing:**
* broken γ = 0 (to divergence or 64 s);
* broken + C⁰-IP (0.3) and (0.03), to 64 s;
* `:native` + C⁰-IP, to 160 s;
* the queued Q2/Q1 32-cell pair;
* the `:native` MMS pair.

### 4.4 Screening S1 update (2026-09-30 19:45)

| arm | fit window | η mid/high/sub σ_E (ratio) | u mid/high/sub σ_E (ratio) | state |
|---|---|---|---|---|
| broken, γ = 0 | 20–40 s | +0.29/+0.16/+0.18 (×18–183) | +0.27/+0.24/+0.25 (×141–161) | **unstable**, η_max 0.131, Newton 10, not yet diverged (mixed died at 26 s) |
| broken + C⁰-IP (0.3, 0.3) | 20–50 s | −0.005/−0.002/−0.001 (×0.28–0.87) | −0.006/+0.008/+0.008 (×0.33–1.6) | **bounded**, Newton 5 |
| broken + C⁰-IP (0.03, 0.03) | 20–38 s | −0.042/−0.005/−0.004 (×0.001–0.28) | −0.033/−0.001/0.000 | **bounded**, Newton 6 |
| `:native` control (main) | 90–111 s | +0.099/+0.030/+0.030 (×1.4–5.5) | +0.121/+0.032/+0.048 (×2.1–15) | exponential phase |
| `:native` + C⁰-IP (0.3, 0.3) | 90–111 s | +0.005/+0.003/+0.005 (×1.2–1.3) | **−0.040/−0.047/−0.047 (×0.38–0.51)** | **no growth**; the early velocity rise has turned into decay |

`:native` MMS (T8): p_η/p_u = 2.888/3.057 at γ = 0; the γ = 0.3 arm is running.

### 4.5 Screening S1, Q3/Q2 arms — final (2026-09-30 22:40)

| arm | end state | late-window rates | envelope (per-8 s max of band energies) |
|---|---|---|---|
| broken, γ = 0 | ⛔ **Newton failure at t ≈ 44 s** (η 0.239, Newton 19 → cap at 50) | +0.16…+0.29 over 20–40 s | exponential |
| broken + C⁰-IP (0.3, 0.3) | ✅ **64 s completed**, Newton 4.87/step | ≤ +0.007 over 30–64 s | **flat** in every band; u band 4 rises 3.2e-9 → 6.5e-9 by 24 s, then plateaus |
| broken + C⁰-IP (0.03, 0.03) | running, t ≈ 59 s, Newton 5 | ≤ 0 over 20–38 s | flat to 40 s |
| `:native` control (main) | 160 s, growing | **+0.07…+0.13** over 110–160 s (mid-band energy ×920) | exponential from ≈ 90 s |
| `:native` + C⁰-IP (0.3, 0.3) | ✅ **160 s completed**, Newton 2.00/step | **≤ +0.005** over 110–160 s (ratios 0.56–1.6) | no growth; carrier decay ≤ 4e-5 s⁻¹ |

⚠ **Read the mid-band start/end ratio with care:** that band breathes by ×30 with the carrier
phase, so an end-point ratio of ×48 coexists with a flat envelope. Judge by the envelope.

**Q2/Q1 32-cell pair** (launched by the batch at 20:40 and 21:29; 40 periods). ⚠ An operator
relaunch at 22:25 started duplicates writing to the same log and output paths. They were
killed within ~2 min. `diagnostics.csv` and the VTK series are intact (checked: rows from step 1,
one per 5 steps), but the two `.log` files were truncated, so their banner and early Newton trace
are lost.

### 4.6 Q2/Q1, 32 cells/λ (2026-10-01 00:50) — ⚠ THE η-PENALTY IS INADMISSIBLE ON A Q1 SURFACE

Band rates over 10–60 s (`PG_NCELL=32 PG_PU=2`):

| arm | carrier σ_E | η mid/high | u mid/high | reading |
|---|---|---|---|---|
| mixed (earlier CN control) | 0 | +0.046 / +0.005 (×14) | +0.038 / +0.017 | slow growth; diverged ≈ 90 s |
| broken, γ = 0 (40 periods) | 0 | **+0.058** / +0.016 (×17) | +0.042 / +0.013 | growing like mixed |
| broken + C⁰-IP (0.3, 0.3) | **−0.021** (energy ×0.35) | −0.037 / −0.022 | +0.019 / −0.025 | **everything decays, the carrier included: fails acceptance 2** |

**Cause, confirmed on the linearised box operator.** Q2/Q1 at 32 cells, carrier damping
(amplitude rate, s⁻¹):

| penalty | carrier damping | oscillatory modes left |
|---|---|---|
| none | 0 | 126 |
| γ_u = 0.3 | 1e-8 | 126 |
| γ_η = 0.3 | 0.0105 (energy 0.021 — **matches the measurement**) | 42 |
| γ_η = 0.03 | 0.00105 | 74 |

The C¹ subspace of a piecewise-**linear** field is the global linear functions, so the η-penalty
on a Q1 surface has no room for any wave: in the penalty limit it forbids the carrier itself. On a
Q2 surface (Q3/Q2) the C¹ quadratic splines contain the carrier, which is why it survived there.
**Rule: never apply γ_η at `p_η = 1`.** The velocity penalty leaves the carrier alone at both
pairings.

Follow-up launched 00:50, 100 periods, same batch:
* `p21_broken_cipu0.3_n32_cn_p100` (γ_u = 0.3, γ_η = 0);
* `p21_broken_n32_cn_p100` (γ = 0 control).

**MMS (T8), complete** — `output/local/broken_mms/mms_cip.csv`, Q3/Q2 spatial ladder:

| model | γ | p_η | p_u | e_η (finest) | e_u (finest) |
|---|---|---|---|---|---|
| linear | 0 | 2.995 | 3.770 | 2.48e-5 | 8.52e-6 |
| linear | 0.3 | 3.089 | 3.771 | 2.51e-5 | 8.19e-6 |
| `:native` | 0 | 2.888 | 3.057 | 2.90e-5 | 3.03e-5 |
| `:native` | 0.3 | 3.089 | 3.591 | 2.51e-5 | 1.76e-5 |

The penalty preserves the orders and does not degrade the errors (acceptance 3 met at Q3/Q2).

### 4.7 Velocity-only penalty at Q3/Q2 (launched 2026-10-01 00:45)

Does γ_u alone suffice, or was γ_η doing the work in §4.5? Q3/Q2, 16 cells/λ, CN, γ_u = 0.3,
γ_η = 0; banners verified:
* `p32_broken_cipu0.3_n16_cn` (`:full` broken, 40 periods). Controls: broken γ = 0, failed at
  44 s; broken (0.3, 0.3), bounded to 64 s.
* `p32_native_cipu0.3_n16_cn_p100` (`:native`, 100 periods). Controls: `:native` γ = 0,
  exponential from ≈ 90 s; `:native` (0.3, 0.3), no growth to 160 s.
* Transition-regime twins, γ_u = 0.03, γ_η = 0 (launched 2026-10-01 ~00:50, banners verified):
  * `p32_broken_cipu0.03_n16_cn` (40 periods);
  * `p32_native_cipu0.03_n16_cn_p100` (100 periods);
  * `p21_broken_cipu0.03_n32_cn_p100` (Q2/Q1, 32 cells/λ, 100 periods).
  Together with the γ_u = 0.3 arms they give a two-point γ_u ladder at each pairing.

### 4.8 Velocity-only ladder — interim (2026-10-01 08:45)

| pairing, tier | γ_u | γ_η | reached | growth | verdict so far |
|---|---|---|---|---|---|
| Q2/Q1 n32 `:full` broken | 0 | 0 | 64 s (Newton 37) | mid +0.18, envelope ×20 over 50–70 s | ⛔ unstable (onset ≈ 50 s) |
| Q2/Q1 n32 `:full` broken | 0.3 | 0 | 97 s | ≤ 1e-3; **envelope flat to 100 s**; carrier −6e-5 s⁻¹ | ✅ bounded, carrier kept |
| Q2/Q1 n32 `:full` broken | 0.03 | 0 | 86 s | ≤ +0.024 (η mid fit); envelope flat to 80 s, **first uptick 80–90 s** (η b4 6.1→6.7e-9, u b4 1.9→2.7e-9) | ⚠ watch |
| Q3/Q2 n16 `:full` broken | 0.3 | 0 | 53 s | −0.002…−0.003 | ✅ bounded |
| Q3/Q2 n16 `:full` broken | 0.03 | 0 | 38 s (Newton 9, r0 0.22) | **+0.14…+0.15** all bands, ×8–25 | ⛔ unstable, like γ = 0 |
| Q3/Q2 n16 `:native` | 0.3 | 0 | 136 s | **mid +0.12 (×100–300)**, high +0.03…+0.04 | ⛔ **same as the unpenalised control** |
| Q3/Q2 n16 `:native` | 0.03 | 0 | 131 s | mid +0.12 (×110–180) | ⛔ same as control |

**Reading.**
* At Q3/Q2 the velocity penalty alone does **not** stabilise `:native`, at either γ. The
  `:native` success of §4.5 came from the **η-penalty**, which on a Q2 surface keeps the carrier.
* For `:full` the velocity penalty suppresses the growth at γ_u = 0.3, on both pairings so far,
  but not at 0.03 (Q3/Q2).
* So the stabilising content is pairing- and tier-dependent:
  * **Q3/Q2:** needs γ_η (plus γ_u = 0.3 for `:full`, unless γ_η alone suffices — untested);
  * **Q2/Q1:** cannot use γ_η (it kills the carrier, §4.6), and γ_u = 0.3 holds `:full` so far.
* The Q2/Q1 `:native` case at 64 cells/λ (diverged at 146 s unpenalised) has not been run with
  any penalty.
* The 40- and 100-period Q2/Q1 controls agree to every printed digit at step 1600: the runs are
  deterministic.

### 4.9 Velocity-only ladder — update (2026-10-01 12:15)

Envelope (per-20 s maxima of band energies) is the criterion:

| pairing, tier | γ_u | state | envelope | verdict |
|---|---|---|---|---|
| Q2/Q1 n32 `:full` broken | 0 | 64 s | takes off ≈ 50 s | ⛔ |
| Q2/Q1 n32 `:full` broken | 0.03 | 123 s, η_max 0.135, Newton 14 | flat to 80 s, then ×10³ by 100–120 s | ⛔ onset delayed ≈ 30 s |
| Q2/Q1 n32 `:full` broken | **0.3** | 154 s, Newton 5 | flat to ≈ 110 s, then **u band 4 1.8e-10 → 1.6e-8 (×80) over 110–160 s**, η band 3 ×4 | ⛔ **onset delayed ≈ 60 s, not removed** |
| Q3/Q2 n16 `:full` broken | **0.3** | ✅ 64 s complete, Newton 5.4 | flat to 64 s | bounded within the window (but the window ends before the Q2/Q1 onset time) |
| Q3/Q2 n16 `:full` broken | 0.03 | 58 s, Newton 11 | +0.15 s⁻¹ | ⛔ |
| Q3/Q2 n16 `:native` | 0.3 / 0.03 | 160 s complete, η_max 0.130 / 0.139 | as the unpenalised control | ⛔ no effect |

**Conclusion of the velocity-only ladder.**
* The ⟦∂ₙ𝖴⟧ penalty **delays** the `:full` instability, and the delay grows with γ_u: onset
  ≈ 50 → 80 → 110 s at γ_u = 0, 0.03, 0.3 on Q2/Q1. It does not remove it.
* On `:native` it has no effect.
* This is what the linear eigen-analysis (§4.2) predicted: velocity oscillatory modes in the C¹
  subspace are undamped, so the growing modes can live there.
* The only configuration bounded over 160 s so far is **Q3/Q2 `:native` with γ_η = 0.3**
  (§4.5), where the η-penalty removes the surface's above-cell-Nyquist half of the wave band and
  the Q2 surface keeps the carrier.
* ⚠ Every `:full` "bounded to 64 s" result (§4.5, §4.7) is now **provisional**: 64 s is shorter
  than the delayed onset seen here.

### 4.10 The Q3/Q2 `:full` confirmations (launched 2026-10-01 ~12:20, banners verified)

100 periods, CN, n16, broken:
* `p32_broken_cipu0.3e0.3_n16_cn_p100` (γ_u = γ_η = 0.3) — extends §4.5 past the delayed onset
  times of §4.9;
* `p32_broken_cipe0.3_n16_cn_p100` (γ_u = 0, γ_η = 0.3) — is the surface penalty alone
  sufficient?

**On Q2/Q1 (the author's position, 2026-10-01): a Q1 surface is not a viable discretisation of
this fourth-order problem.** Evidence:
* `η ∈ Q1` makes ∂²η vanish on every cell in 1-D, so the Class-III free-surface curvature lives
  only in the skeleton layer ⟦∂ₙη⟧;
* the only stabiliser that worked (γ_η) is inadmissible on Q1 (§4.6);
* the velocity penalty only delays the `:full` onset (§4.9);
* (caveat) Q2/Q1 `:native` was stable under CN up to 32 cells/λ for 100 periods (CLAUDE.md
  §5.2f), so the claim is about the nonlinear Class-III tier and refinement, not about every
  configuration.

### 4.11 Flume test, Q3/Q2 `:full` broken + C⁰-IP (launched 2026-10-01 13:10)

`run/local/run_1dbrk_q32_flume.sh`. The environment is copied from `run_1dc3q_base_mixed.sh`
(60 m flume, nx = 240, dx = 0.25, dt = 0.04, A = 0.10, relaxation inflow, flat bed, ny = 1,
c3_mask = gs). The deliberate changes: broken instead of mixed, the penalty, and **Crank–Nicolson**
instead of SDIRK_2_2. 100 periods requested.
* `brk_q32_cipu0.3e0.3_cn` (γ_u = γ_η = 0.3);
* `brk_q32_cn` (γ = 0, same-batch control).

Banners verified: 14 406 free DOFs, 239 interior facets, θ = 0.5. Logs in
`output/local_1d/_logs/`.

**Cost:** ≈ 27–31 s per step with eight simulations sharing the machine (projected 4 000 steps
≈ 30–35 h; the printed ETA includes JIT and overestimates early). Milestones:
* 16 s, where mixed died under SDIRK: ≈ 3.5 h;
* 36 s, the flume fill time (rule 14): ≈ 8 h;
* 64 s: ≈ 13 h.

The `run_flume_1d.jl` knobs `BALFEM_BROKEN` / `BALFEM_CIP_GU` / `BALFEM_CIP_GE` /
`BALFEM_CIP_HEXP` are passed on the **sequential** path only. The MPI path errors if they are set.

### 4.12 Status 2026-10-01 22:50

* **Flume** (§4.11): both arms are past t = 16 s, where mixed died under SDIRK.
  * C⁰-IP: t = 19.0 s, η_max 0.125, Newton 6;
  * control: t = 18.0 s, η_max 0.109, Newton 9.
  * Now ≈ 45–54 s/step, so 160 s is ≈ 2 days.
  * Still inside the 36 s fill: nothing readable yet (rule 14).
* **Box, Q3/Q2 broken (0.3, 0.3), 100 periods:** envelope flat to 64 s, identical to the
  40-period run (deterministic). The window that decides it — beyond the Q2/Q1 delayed onset at
  ≈ 110 s — is still ahead.
* **Box, Q3/Q2 broken, surface-only (γ_η = 0.3):** bounded to 63 s; η bands at seed level. The u
  high band sits ~20× above the (0.3, 0.3) arm (7e-8 to 1.2e-7) but is not growing; Newton 8
  (vs 5).
* **Velocity-only ladder, final:**
  * Q2/Q1 γ_u = 0.03: ⛔ Newton failure at t ≈ 129 s;
  * Q2/Q1 γ_u = 0.3: completed 160 s but **grew from ≈ 120 s** (u band 4 ×200, η band 3 ×10 by
    140–160 s);
  * Q3/Q2 γ_u = 0.03: completed 64 s with η_max 0.153 and Newton 19, ⛔ growing.

### 4.13 ⛔ The `:full` box runs grow late, in the MID band (2026-10-02 09:45)

Box, Q3/Q2, n16, CN, fit 115–160 s (σ_E in s⁻¹, η / u):

| arm | mid (5–12) | high (13–24) | sub-element (25–48) | onset |
|---|---|---|---|---|
| `:full` broken (0.3, 0.3) | **+0.13 / +0.14** (×128 / ×254) | +0.007 / +0.014 | +0.07 / +0.09 | ≈ 110 s |
| `:full` broken, γ_η = 0.3 only | **+0.13 / +0.09** (×32 / ×25) | +0.009 / +0.016 | +0.07 / +0.05 | ≈ 90–110 s |
| `:native` (0.3, 0.3) | −0.017 / +0.006 (×1.1 / ×2.0) | 0 / −0.006 | 0 / −0.008 | none — ✅ |
| `:native` control | +0.12 / +0.12 | +0.07 / +0.06 | +0.07 / +0.08 | ≈ 90 s |

**Reading.**
* For `:full` the penalty removes the GRID-scale growth and delays the onset (γ = 0 died at 44 s;
  (0.3, 0.3) grows from ≈ 110 s). A mid-band mode at **λ ≈ 0.33–0.8 m (4–10 nodes per
  wavelength)** then grows at the same rate as the unpenalised `:native` instability.
* These wavelengths are representable by C¹ splines, so — as the eigen-analysis (§4.2) said —
  the first-derivative jump penalty cannot reach them.
* The `:full` mid band also carries a forced level 1000× above `:native`'s (1e-8 vs 1e-11): the
  Class-III terms feed it directly.

**Flume (§4.11):**
* control (γ = 0): ⛔ **Newton failure at t ≈ 38 s**, η_max 0.155 at x = 43.6 m, Newton 16 by
  36–40 s;
* C⁰-IP: t = 51 s, η_max settling 0.155 → 0.114 after the front passes, Newton ≤ 9. Its
  readable window (> 36 s) has just opened; the box onset suggests watching to ≥ 120 s.

**Status of the stabiliser:**
* ✅ the only clean 160 s pass is **`:native` + (γ_u, γ_η) = (0.3, 0.3) at Q3/Q2**;
* for `:full`, the C⁰-IP on first-derivative jumps **delays but does not remove** the instability.

**Next candidates:**
* higher-order jump penalties (⟦∂ₙ²𝖴⟧, ⟦∂ₙ³𝖴⟧ — Burman–Ern), which also act on the C¹ subspace;
* a spectral or low-pass filter on the mid/high band;
* a resolution ladder of the late `:full` mid-band rate, to tell a grid-scale mechanism (rate
  changes with h) from a continuum one (rate converges).

---

## 5. Where we stand and how to proceed (2026-10-02)

### 5.1 What is established

* **Broken Class-III operator** (form A: cellwise Hessians + skeleton layer): implemented,
  verified to round-off (G1), consistent with the projected treatment (G4, order 3). On its own it
  is **not** a stabiliser — the box at Q3/Q2 died at 44 s, the flume at 38 s — as the audit
  predicted.
* **C⁰-IP on first-derivative jumps** (`γ_u`, `γ_η`):
  * removes the GRID-scale growth everywhere and delays the `:full` onset (box: 44 s → ≈ 110 s;
    flume: still alive at 51 s, where the control died at 38 s);
  * does **not** remove the late `:full` growth, which sits in the **mid band** (λ ≈ 0.33–0.8 m,
    4–10 nodes per wavelength) at the same rate as the unpenalised `:native` instability
    (+0.13 s⁻¹);
  * the velocity penalty alone only delays the instability (§4.9); on `:native` it does nothing;
  * ✅ the only clean 160 s pass is `:native` + (γ_u, γ_η) = (0.3, 0.3) at Q3/Q2.
* **The penalty acts through its KERNEL** (§4.2). Oscillatory modes that fit in
  K = {J = 0} are not damped; at large γ the solution is forced into K. A first-derivative penalty
  leaves K = C¹ splines, which carry the mid band — hence the escape.
* **Q1 surfaces are excluded from the surface penalty**:
  * K(⟦∂ₙη⟧ on Q1) = global linear functions, i.e. constants in the box;
  * the carrier is damped (measured σ_E = −0.021 s⁻¹ at γ_η = 0.3, and reproduced exactly by
    the linear operator);
  * together with §4.9, this supports the position that a Q1 surface is not a viable
    discretisation of this fourth-order problem in the `:full` tier.

### 5.2 Lead 1 — high-order jump penalty (Burman–Ern hp-CIP)

    J_h = Σ_F Σ_{j=1..p_field} γ τ h_F^{2j} ⟦∂ₙʲ·⟧⟦∂ₙʲ·⟧      (both η and 𝖴)

* **Kernel:** a single polynomial of degree p across each two-cell patch, i.e. a global polynomial;
  in the periodic box, the constants. **No escape subspace.**
* **Weight per order:** h^{2j} makes every order contribute
  * O(γ√(gd)/h) of damping on grid-scale modes,
  * ∝ (kh)^{2p+2} on resolved modes,

  so all orders share one consistency order and the MMS rates should survive (to be measured).
* **Selectivity at Q3/Q2, 16 cells/λ:** carrier kh ≈ 0.39, mid band kh ≈ 2–5. With (kh)⁸ the
  damping ratio is ~10⁷–10⁹, so a wide γ window should exist in principle.
* ⚠ **LOCKING** (the price of a trivial kernel): with γ too large the solution is forced into K —
  a constant — and the carrier is suppressed, which is the Q1-surface failure on every field. The
  penalty must stay in the damping regime; the eigen-analysis must show the window:
  * mid-band damping ≳ 0.13 s⁻¹;
  * carrier damping ≪ 1/T_run ≈ 0.006 s⁻¹.

  Burman–Hansbo–Larson ("locking-free ghost penalty", Numer. Math. 2025/26) show that γ → ∞ is
  safe only when K is a deliberately rich space; a global polynomial is not one.

### 5.3 The derivative limit in Gridap (checked 2026-10-02)

* `Fields/FieldsInterfaces.jl:85–86` defines only `gradient(f, Val{1})` and `gradient(f, Val{2})`.
* `Polynomials/PolynomialInterfaces.jl` evaluates `FieldGradientArray{N,<:PolynomialBasis}` for
  N = 1, 2 only.
* Measured: on a `MonomialBasis` order 3 throws `no method matching gradient(::Monomial, ::Val{3})`,
  and the same holds for analytic `Function` fields.
* So **⟦∂ₙ³𝖴⟧ (needed for Q3 velocity) is not available**. Orders 1–2 are complete for η at
  Q3/Q2 and for both fields at Q2/Q1.

### 5.4 GridapEmbedded.jl (cloned at the repo root, v0.9.12) — no help for the derivatives

* `GhostSkeleton(cut)` (`src/Interfaces/EmbeddedDiscretizations.jl:506–547`) is only a facet mask
  → `SkeletonTriangulation(model, mask)`, which equals our skeleton on an uncut mesh.
* There is no penalty operator in `src/`. Every ghost penalty in its examples and tests is
  user-written and **first-order on P1 elements**:
  * `(γ h) jump(n⋅∇u) jump(n⋅∇v)`;
  * `(β h³) jump(n⋅∇p) jump(n⋅∇q)` for the pressure.
* It builds on Gridap's derivatives, so it inherits the order-2 limit. AgFEM (cell aggregation)
  is a different mechanism.

### 5.5 Lead 2 — the direct (volume) ghost penalty

    J_h = Σ_F γ τ h⁻¹ ∫_{T₁∪T₂} (u₁ − E u₂)·Mv (v₁ − E v₂)

E is the canonical polynomial extension of the neighbour's local polynomial.

* **Origin and use:** Preuß (2018); analysed in Lehrenfeld & Olshanskii (M2AN 2019). Used for
  moving-domain Stokes, FSI contact, and Maxwell ("direct extension stabilization", 2024). The
  reference implementation is believed to be ngsxfem (not verified). It is in neither Gridap nor
  GridapEmbedded.
* **Equivalence:** as a norm, it is equivalent to Σ_{j=0..p} h^{2j+1}‖⟦∂ₙʲu⟧‖²_F. It has the
  same kernel as the full hp-CIP and covers **every order without any derivative**. The h⁻¹
  scaling reproduces our h^{2j} family.
* **Weakly consistent** (O(h^{p+1}) on the interpolant of the smooth solution) rather than
  strongly consistent. The abstract theory (Gürkan & Massing 2019) gives optimal rates for
  elliptic, parabolic, Stokes and advection-reaction problems. For our dispersive system this is
  **to be measured with the MMS ladder**.
* **Implementation:** on the uniform Cartesian mesh, E is a ±1 reference-coordinate shift across
  the facet. That is custom code, but no library change.
* **Same locking caveat as Lead 1.**

### 5.6 Other leads, not yet pursued

* Extend the Gridap fork's polynomial bases to third derivatives. This modifies a library, so it
  needs the author's approval.
* A low-pass or spectral filter on the mid/high band.
* **A resolution ladder of the late `:full` mid-band rate:**
  * a rate that converges with h ⇒ a property of the continuum model, and any penalty is then a
    regularisation;
  * a rate that changes with h ⇒ a discretisation effect.

### 5.7 Plan (in order)

1. **Implement hp-CIP orders j ≤ min(2, p_field)** for η and 𝖴 (knob `cip_order`, weights
   h^{s+2(j−1)}, s = 2).
2. **Linear eigen-analysis** of the box operator (Q3/Q2 n16 and Q2/Q1 n32): damping of each
   oscillatory mode against its dominant wavenumber band, over a γ ladder, for orders 1 and 1–2.
   Read off the window: mid band ≥ 0.13 s⁻¹ with carrier ≪ 0.006 s⁻¹.
3. If the window exists, run Q3/Q2 box runs (`:full` broken and `:native`, CN, 160 s) at the
   chosen γ, plus the MMS order check.
4. If the missing ⟦∂ₙ³𝖴⟧ matters (mid band undamped with orders ≤ 2 at Q3/Q2), implement
   Lead 2.

### 5.8 hp-CIP orders ≤ 2 — implemented, and the eigen-analysis (2026-10-02)

**Code:**
* `attach_skeleton!(…; cip_order, p_u, p_eta)`: orders j ≤ min(cip_order, p_field), weight
  h_F^(s+2(j−1)); `CIP_MAX_ORDER = 2`, order 3 refused.
* `_jdn2`: ⟦∂ₙ²f⟧ = n·(∇∇f⁺ − ∇∇f⁻)·n from `alg_hess`.
* Knob `cip_order` in `setup_and_run`, `run_mms_case`, and `BALFEM_CIP_ORDER` in both 1-D
  drivers.
* Gate G8 in `test_broken_formulation.jl`: symmetric PSD; zero on a global polynomial of the trial
  degree; the order-2 terms raise the rank; order 3 refused.

**Analysis:** `examples/local_1d/cip_eigen_analysis.jl` → `output/local_1d/cip_eigen/cip_eigen.csv`.
* The linearised flat-bed box operator −M_eff⁻¹(K + J_h) is diagonalised.
* Oscillatory modes are classified by the dominant x-harmonic into the `periodic_growth.jl` bands;
  numbers are ENERGY damping (s⁻¹, comparable with σ_E).
* The carrier is tracked as the x-wave nearest ω = 3.913.
* ⚠ **x-uniform modes are excluded.** The one-cell-wide flume has a transverse standing mode
  (ω ≈ 7.55, η varying in y only) that the x-facet penalty cannot see and a y-invariant run never
  excites. Left in, it read as an "undamped mid-band mode" with a random band label — that cost one
  wrong reading of the first ladder.

**Q3/Q2, 16 cells/λ, both fields, orders ≤ 2 — A WINDOW EXISTS:**

| γ | carrier damping (energy lost in 160 s) | mid band min (median) | high band min | oscillatory pairs |
|---|---|---|---|---|
| 1e-4 | 8.2e-6 | ~0 (1.3) | 1.7 | 95 (pure damping) |
| 3e-4 | 2.4e-5 | ~0 | 5.2 | 95 |
| **1e-3** | **7.6e-5 (1.2 %)** | **0.27** (1.6) | 33 | 73 |
| **2e-3** | **1.5e-4 (2.4 %)** | **0.55** | 67 | 54 |
| **3e-3** | **2.25e-4 (3.5 %)** | **0.82** | 101 | 50 |
| 1e-2 | 7.5e-4 (11 %) | ~0 (some mid modes escape) | 337 | 48 — locking sets in |

* Carrier damping ≈ 0.075·γ (perturbative, linear in γ).
* Every x-varying mid/high-band mode is damped above the +0.13 s⁻¹ late `:full` growth rate for
  γ ∈ [1e-3, 3e-3], with margin ×2…×6.

**Without the η-penalty, or with order 1 only, there is no window at Q3/Q2:**
* velocity-only, orders ≤ 2: the high band reaches 0.2–2 s⁻¹, but mid-band modes stay
  undamped;
* order 1: the C¹ escape — mid-band minima ~0 at every γ.

**Q2/Q1, 32 cells/λ: NO WINDOW at any γ or order.** Mid-band modes stay undamped until the
carrier itself is damped at the 1e-3 s⁻¹ level (γ ≈ 0.03). This is consistent with §4.6 and with
the position that Q1 surfaces are not viable for `:full`.

**Caveat:** the order-3 velocity jump (⟦∂ₙ³𝖴⟧) is still missing at Q3/Q2. The window exists
**without** it because the η-penalty at orders 1–2 has a trivial kernel and every η-carrying mode
is damped. Velocity-dominated modes with 𝖴 ∈ C² cubic splines could still escape; in the linear
operator none of the x-modes does in the window.

**Next (§5.7 step 3):** Q3/Q2 box runs, `:full` broken and `:native`, CN, 160 s, with
(γ_u, γ_η) = (2e-3, 2e-3) and `cip_order = 2` (plus 1e-3 / 3e-3), and the MMS order check at the
chosen γ.

### 5.9 Box runs with the hp-CIP window value (launched 2026-10-02 11:58)

Q3/Q2, n16, CN, 100 periods, (γ_u, γ_η) = (2e-3, 2e-3), `cip_order = 2`; banners verified
(`j ≤ 2` on both fields, 16 interior facets):
* `p32_broken_cip2e-3o2_n16_cn_p100` (`:full` broken). Controls on the same branch: γ = 0
  (dead at 44 s); first-order (0.3, 0.3) (mid-band growth from ≈ 110 s, §4.13).
* `p32_native_cip2e-3o2_n16_cn_p100` (`:native`). Controls: γ = 0 (growth from ≈ 90 s);
  first-order (0.3, 0.3) (no growth to 160 s).

The `:native` arm is stepping at Newton 2 with mass drift ~1e-17: the order-2 jumps run inside the
time loop. Expected: `:native` ≈ 8 h, `:full` ≈ 20 h.
