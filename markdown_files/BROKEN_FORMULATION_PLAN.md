# BROKEN_FORMULATION_PLAN.md — the broken (skeleton) Class-III formulation and the C⁰-IP stabiliser

*Branch `broken-formulation-solver`, from `main` @ `58ef9bb` (2026-09-30). The mathematics is
LaTeX §6.4 "Broken Weak Formulation: Term-by-Term Audit"
(`latex_docs/BALFEM_models/NumericalImplementation/BrokenAudit.tex`); this file is the
implementation and campaign record. Companion records: `NEW_TREATMENT.md` (projected and mixed
treatments), `OPEN_ISSUES.md` §0e (the instability), `PLANNED_CAMPAIGNS.md` §6c (stabilisation
campaign).*

---

## 0. What is being built, and why

The audit's conclusions fix the scope:

1. **The Galerkin residual is consistent.** No test-jump (kind T) facet term is missing from it.
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
