# `yl_collapse_wip/` — stage 2 of the Yang & Liu operator comparison, IN PROGRESS

Stage 1 (vertical velocity `w`) is **complete and passing** — it lives in the real
suite as `test/test_yl_collapse.jl`, 40/40 at ~1e-14.

`stage2_pressure_wip.jl` compares the non-hydrostatic pressure against
(2.22)–(2.23) + supplementary (B.1)–(B.4). **It is not a gate yet and must not be
added to `runtests.jl`.** State as of 2026-09-15:

| finding | status |
|---|---|
| Sign convention: their polynomial part is `−` ours | **established** — `Δ/\|p\|` was exactly 1.994 before the flip |
| Linear pressure collapses exactly | **established** — after the flip the residual scales as `ε²`, not `ε` |
| Quadratic (`𝓝·Λ`) part differs | **open** — clean `O(ε²)` residual, coefficient 0.87, constant to 2 % over `ε = 1 … 0.01` |

⚠ **The open item is NOT yet attributable.** The probability mass is on transcription
in this file — (B.1)–(B.4) are four long formulas, and `Θ_kj` components 6 and 7
carry `φ'_k` — not on the derivation. Do not quote it as a discrepancy with
\cite{LFEM}.

**The decisive next step** avoids re-reading either transcription: stage 1 proved
our `w` ≡ their `w` to 1e-14, so compute `Dw/Dt` *directly from that common,
verified `w` field* and test which pressure satisfies `∂p_nh/∂σ = −ρH Dw/Dt`.
Whichever side fails owns the error, and the test needs neither §B nor `Θ`.

Also note (2.23) as printed adds `Σₙ p_{M,n}`, which does not satisfy the dynamic
boundary condition `p|_{σ=1} = 0`; this file fixes the constant from (2.5) plus
interface continuity instead, which is what the surrounding text says (2.23) encodes.

---

## UPDATE — the `Dw/Dt` discriminator has run (2026-09-15)

`stage2b_dwdt_discriminator.jl`. Builds `Dw/Dt` from the **common, stage-1-verified** `w` field and
asks which pressure satisfies the vertical momentum equation `∂p_nh/∂σ = −ρH·Dw/Dt`.

| | max residual | scale |
|---|---|---|
| **Yang & Liu's `p`** | **2e-14 … 2e-13** | 2–7 |
| **our `p_nh`** | **3.16** (ε=1), 3.10e-02 (ε=0.1) | 6.44, 0.42 |

Theirs is round-off. Ours is ~50 % of scale and falls 102× when ε falls 10× — **O(ε²)**, the
quadratic package `𝓝·Λ`. ⚠ **The error is on our side.**

### Attribution — how far it is pinned down

* `Θ_kj`: the transcription used in this test was checked against the solver's own tensors,
  `Pcal[i,k,j,c] = ∫ Θ_kj[c] φᵢ_int dσ` from `src/vertical.jl`. **Identical, all 8 components,
  max diff 4.5e-16.** So the Θ profiles here ARE the code's, and the failure is in the *shared*
  derivation, not in this script's reading of Θ.
* `𝓝_kj`: transcribed here from the LaTeX; **not yet cross-checked against `problem.jl`**. That is
  the one remaining gap before the defect can be localised to Θ vs 𝓝.

### Why this matters beyond the chapter

`𝓝` is the non-linear pressure package, and `{1,2,4,5}` of it is Class III. A derivation defect there
is a **candidate explanation for the `:full` instability that has nothing to do with the frozen
projections** — the operator itself would be inconsistent with the `w` it is built from.

⚠ **Not patched.** Per the project's standing policy this is reported, not fixed: the defect appears
to originate in previously validated code/derivation (`src/vertical.jl`, `pDerivation.tex`), and the
right response is a decision, not a workaround.

### Next step to localise it

Zero the eight `𝓝`/`Θ` components one at a time and re-measure the residual; the identity
`∂p_nh/∂σ = −ρH·Dw/Dt` is self-contained (it needs only the verified `w`), so it can be used as an
oracle to say which component is wrong.

---

## UPDATE 2 — the transcript is CONFIRMED CORRECT by the author, and the defect is OURS (2026-09-15)

`stage2c_localise_N.jl`.

### The chain of evidence, each link measured

| link | evidence |
|---|---|
| our `w` ≡ their `w` | stage 1, **1e-14**, M = 2,3,4 |
| `Dw/Dt` built from that common `w` is correct | their `p` satisfies `∂p_nh/∂σ = −ρH·Dw/Dt` to **2e-14** using it |
| §A + §B transcription correct | confirmed by the author, **and** independently: their `p` built from it satisfies the identity |
| my `Θ_kj` ≡ the solver's | `Pcal` rebuilt from it matches `src/vertical.jl`, all 8 components, **4.5e-16** |
| my `𝓝_kj` ≡ the LaTeX | checked line by line against `eq: def Nkj nh pressure derivative`, 8/8 |

⛔ **Conclusion: our `Σ_j 𝓛_j·θ_j + Σ_kj 𝓝_kj·Θ_kj` does not equal `Dw/Dt`.** The linear half is
exact; the quadratic half is wrong by ~50 % of scale.

### Localisation

* **NOT the `∇h` terms.** Flat bed gives `|E| = 1.181`, variable bed `1.135` — essentially
  identical, so components 3 and 6 are exonerated.
* **NOT a per-component coefficient slip.** Least squares of `E ≈ Σ_c α_c C_c` over 200 samples,
  two meshes, four stations: relative residual **0.43** (all 8) and **0.47** (flat bed, 6 active).
  No single component does better than 0.76. **The eight σ-profiles do not span the error**, so a
  term type is missing or a profile is wrong — a structural defect, not a constant.

### Why nothing caught this before — and it is not a contradiction

`:native` is MMS-verified at optimal order, and that remains true. The MMS verifies that the solver
solves the equations **the specification encodes**; `mms_forcing` is derived from the same `𝓝`/`Θ`,
so a consistently-wrong specification passes it (rule 28, in its strongest form). This test verifies
the **specification against the physics** — the vertical momentum equation — which is a different
question and the one no existing gate asks.

Magnitude is why it stayed invisible: `𝓝`'s whole contribution is 0.013 % (1-D) at `A = 1e-3`, the
amplitude most of the verified scope was measured at. It is `O(ε²)`, so it grows quadratically with
amplitude.

⚠ **Candidate explanation for the `:full` instability that is INDEPENDENT of the frozen projections.**
Not established — but it now competes with the projection hypothesis and should be resolved first,
because if the operator is wrong the projection question is secondary.

⚠ **NOT PATCHED.** Standing policy: this originates in previously validated derivation and code
(`VerticalFESemiDiscretisation/pDerivation.tex`, `src/vertical.jl`, `src/problem.jl`). Reported for a
decision. The next step is a derivation review of `pDerivation.tex` — specifically the expansions
feeding `Θ_kj` components {1,2,4,5,7,8} — not another numerical experiment.


---

# ⚠⚠ RETRACTION — UPDATE 3 (2026-09-15). THE DERIVATION IS CORRECT.

**Everything in "UPDATE 2" claiming a defect in `pDerivation.tex` / `src/vertical.jl` is WITHDRAWN.**
The review of the expansions found the flaw, and it was in the TEST, twice over.

## Final result

`Σ_j 𝓛_j·θ_j + Σ_kj 𝓝_kj·Θ_kj  ==  Dw/Dt`, **worst relative difference 1.7e-15**, over
M = 2,3,4 (4 meshes), 3 stations, 2 amplitudes. The expansions are exact.

## The two errors, both mine

1. **The test state violated the depth-integrated continuity constraint.**
   `pDerivation.tex` eliminates `∂H/∂t` using `∂H/∂t = -Σ_j ∇·(Hu_j)Φ_j`, both directly (Θ
   components 1, 2) and inside `ω` (6, 7, 8). Our package is therefore an identity **only on that
   constraint**. Yang & Liu keep `H^(0,1)` explicit, so theirs holds for any state — which is
   exactly why theirs passed and ours did not. Fixed by building `H(x,t) = H0(x) + (t-t0)G(x)`
   with `G = -Σ_j ∂_x[H0 u_j]Φ_j`.
2. **A closure bug that silently dropped half a product rule.** The helper
   `d1(f) = ForwardDiff.derivative(z->f(z,t), x)` returns a NUMBER, so
   `(z,τ) -> d1(Hu(k))*u_j(z,τ)` captured `∂_x(Hu_k)` as a constant and computed
   `∂_x(Hu_k)·∂_x u_j` instead of `∂_x[(∂_x(Hu_k))u_j]` — losing the `u_j ∂_x²(Hu_k)` half of
   `𝓝_2`. This is CLAUDE.md rule 7's hazard in a new dress: not a shadowed name, but a *value*
   captured where a *function* was meant. The original `yl_stage2b.jl` had it right; the rewritten
   helper did not.

## How each was localised (the method is the reusable part)

* Expansion-by-expansion: `u·∇w` exact at 2e-16, `(ω/H)∂w/∂σ` exact at 6e-17, `∂w/∂t` wrong at
  7.8e-01 — isolating one of three.
* Then unsubstituted-vs-substituted within `∂w/∂t`: the raw product rule exact at 2e-16, the
  continuity-substituted form wrong — isolating the substitution.
* Then the constraint itself (exact, diff 0) and term 1 (exact, 5e-17) — leaving term 2, which is
  where the closure bug was.

## Stage 2 is closed by transitivity

No further §B run is needed: our package ≡ `Dw/Dt` (1.7e-15) and their `p` satisfies
`∂p_nh/∂σ = -ρH·Dw/Dt` (2e-14), and both pressures vanish at `σ=1`. Therefore
**our `p_nh` ≡ their `p_nh`**. §A and §B both collapse.

## Consequence for the `:full` instability

The competing "the operator is wrong" hypothesis is **dead**. The `𝓝` package is verified against
the vertical momentum equation. **The frozen-projection hypothesis is again the leading
explanation**, and the `dx`/`dt` evidence (onset advancing under `dx`-refinement, retreating under
`dt`-refinement) stands unchallenged by this work.


---

# ✅ STAGE 3 PASSES — THE COLLAPSE IS COMPLETE (2026-09-15)

`stage3_weighted_residual.jl`.

**§C was not transcribed, deliberately.** §C is (2.24) expanded with `w̃` and `p̃`, and stages 1–2
verified both of those fields against ours. So §C carries no information beyond (2.24) itself, while
transcribing five more long formulas would add risk without adding evidence. The test runs against
(2.24)+(2.26) directly.

| | |
|---|---|
| **ours** | the tensor-form semi-discrete momentum residual, assembled from the **solver's own** `Mmat, Mcal, Gcal, A, K, P, Acal, Kcal, Pcal` |
| **theirs** | `H · ∫₀¹ φᵢ R(σ) dσ`, with `R` the pointwise residual of (2.12) built from the verified `u`, `w`, `p_nh` |
| **result** | **worst relative difference 1.2e-13**, M = 2 and 3, 2 stations, 2 amplitudes, every DOF |

⚠ This is the first stage that tests the **assembled tensors** rather than a hand-rolled Θ. The
`𝓜/𝓖` advection tensors and the `A/K/𝓐/𝓚/P/𝓟` pressure tensors are all exercised.

## The collapse, complete

| stage | object | worst relative difference |
|---|---|---|
| 1 | vertical velocity `w` (§A) | 1.1e-14 |
| 2 | non-hydrostatic pressure `p_nh` (§B) | 1.7e-15 |
| 3 | weighted momentum residual (§C / 2.24+2.26) | 1.2e-13 |

Together with the linear coefficient check (`A`, `B`, `D` to 8e-15), **BALFE-M at p = 1 reproduces
LFE-M in full — linear and non-linear, coefficients and operators.**

## What this settles

* The `𝓝` package, **including Class III `{1,2,4,5}`**, is verified against an independently
  published derivation. The "wrong operator" hypothesis for the `:full` instability is closed.
* The **frozen-projection hypothesis is the remaining explanation**, and the `dx`/`dt` ladder
  evidence stands.
* ⚠ Note what is NOT tested here: this compares OPERATORS at a prescribed state. It says nothing
  about how the Class-III blocks are *assembled* in the solver (frozen projections), which is
  precisely the open question.


---

# STAGE 3 REDONE AGAINST THE PUBLISHED ⌊2.30⌋ (2026-09-17)

The earlier stage 3 compared our tensor-form residual against `H·∫φᵢ R dσ`, i.e. against the
**projection of the governing equation** ⌊2.12⌋ built from our own verified fields. That verified our
tensor assembly, but it was **not** a comparison against Yang & Liu's published model — and it left
the §C transcript decorative.

`yl_collapse_stages23.jl` now compares against **⌊2.30⌋ itself**:

| | |
|---|---|
| theirs | `Σₙ (R_{k−1,n}F₁(k,n) + R_{k,n}F₂(k,n))`, with `R_{k,n}` from supplementary §C (C.1)–(C.5) and `F₁,F₂` from ⌊2.28⌋–⌊2.29⌋ — nothing of ours enters beyond the prescribed state |
| ours | the `p=1` semi-discrete momentum residual from the **solver's own** vertical tensors |
| result | **worst relative difference 2.1e-12**, 42 comparisons (2 meshes × 3 stations × 2 amplitudes × every DOF) |

**Free check that gates the weighting integrals:** `F₁(k,0)+F₂(k,0)` must equal `Φₖ`, since both are
`∫φₖ dσ` computed two ways. At `c₂ = 0.728`: 0.364000000000 / 0.500000000000 / 0.136000000000 against
`Φ = (0.364, 0.5, 0.136)` — exact.

## Two conventions had to be resolved, and both are findings about the paper

**(a) Their published `p_{k,n}` are the coefficients of `−p_nh/(ρH)`, not `+p_nh`.** This follows from
their own ⌊3.4⌋ (`p_{k,1} = (w_{k,0})_t`), whereas integrating ⌊2.13⌋ down from the surface gives a
σ¹ coefficient of `−w_{0,t}`. §C, however, is written in terms of the **true** pressure polynomial:
the σᵐ coefficient must be `(1−m)H_x P_m + H P_m^(1,0) + (m+1)h' P_{m+1}`, which reproduces all five
printed tails of (C.1)–(C.5) exactly. **§B and §C therefore carry opposite sign conventions**, and §C
must be fed with `−p_{k,n}`.

**(b) The surface condition already supplies the hydrostatic share of `P_0`.** Expanding `g(1−σ)` puts
`−g` into `P_1`; ⌊2.5⌋ then fixes `P_0 = −Σ_{n≥1}P_n`, which already contains the compensating `+g`.
Writing `P_0 = g + C` double-counts gravity and appears as a constant offset of exactly `g·H_x` in
`R_{k,0}` — which is how it was found (offset − `g·η_x`, divided by `h'`, came out 9.82 / — / 9.81 at
three stations).

⚠ **(b) was my bug, the fourth in this exercise.** All four were localised the same way: reduce the
discrepancy to a *structure* — a scaling in ε, a constant in σ, a per-coefficient sign — and let the
structure name the term. That method has now worked four times where inspection did not.
