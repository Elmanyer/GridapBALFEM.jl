# LATEX_STRUCTURE.md — the chapter order of `latex_docs/BALFEM_models/` and why

> **Revision 4 (2026-10-05) — PROPOSED, not yet applied to the LaTeX.**
>
> **What changes.**
> * The two implementation chapters are reorganised around the development *story* of the
>   horizontal discretisation, not around the code.
> * The old chapter 6 (*Gridap Multi-field Residual*) is split in two:
>   * the **mathematics** of the horizontal Galerkin weak form and its regularity;
>   * the **multi-field implementation** of that weak form.
> * Stabilisation gets a chapter of its own, after the stability analysis that proves it is needed.
> * The Gridap solver functionalities move to the end, after the core discretisation has been
>   shown to work.
>
> Revisions 1–3 are summarised under "History" at the end; their full text is in git.
> The open decisions this plan needs from the author are listed in §5.

---

## 1. The principle

**The document follows the development of the solver. At each step something is built, then
analysed; the analysis finds an issue, and the next step fixes it.** A chapter uses only what the
chapters before it have built. A result is placed at the earliest point where the machinery it
needs exists.

The core of the document, chapters 6–9, is one chain of four questions. Each chapter answers one
question, and its answer raises the next:

```
 6  Is the horizontal weak form of the exact model well defined on C⁰ elements, term by term?
      → classify every term (𝓛, 𝓝 classes I/II/III) by the regularity it demands; only on the
        broken (element-wise) domain is it consistent; two consistent formulations result
 7  How is that weak form written in Gridap and advanced in time?
      → the stacked multi-field layout, every term written in it, time integration and Newton,
        and the unit vertical basis that supplies the tensors
 8  Is the resulting discretisation stable?
      → NO: consistent, but unstable at the grid scale; it needs a sink
 9  What sink, and does it work?
      → a scale-selective penalty, chosen by a linear eigen-analysis; validated in the closed box
```

**The key distinction the order protects: consistency is not stability.** Chapter 6 makes the
discretisation *consistent*: every distributional term is accounted for. Chapter 8 then shows that
a consistent discretisation can still be unstable.

If stabilisation were discussed alongside the regularity fix, the two would be conflated. That is
exactly the confusion the project lived through: the `C⁰`-Hessian problem was held responsible
for the instability for a year (rule 39b). Keeping them in separate chapters, in this order, is
what lets the reader see that they are two different problems.

The machinery grows in six layers:

| layer | what exists at that point | chapters |
|---|---|---|
| **I. The model** | the continuous equations, the vertical FE basis and its tensors | 1–3 |
| **II. Checking the model** | the above plus published equations and analytic wave theory; **no horizontal discretisation, solver or time integration** | 4–5 |
| **III. Discretising the model horizontally** | the Galerkin weak form on `C⁰` elements, term by term with its regularity (6); its stacked multi-field realisation and time integration (7) | 6–7 |
| **IV. Assessing the discretisation** | a running solver, in a closed periodic box (**no boundary machinery needed**) | 8 |
| **V. Fixing the discretisation** | the stability diagnosis, and the stabilisers measured against it | 9 |
| **VI. Making it a usable solver** | wave generation, absorption, boundary conditions, workflow, distributed path, validation | 10–11 |

Layers II and IV still answer the two questions that must be kept apart: *is the model right?* and
*is the discretisation right?* What is new is that layer IV no longer depends on the boundary
machinery. The closed periodic box needs no wave generation, sponge or relaxation zone. So the
stability analysis can come *before* the chapter that introduces them, and must. The flume
experiments could never separate an interior instability from a boundary feed (rule 38i). Placing
the boundary machinery after the stability analysis encodes that lesson in the order of the
document.

---

## 2. The proposed order

| # | chapter title (working) | file(s) | label | status |
|---|---|---|---|---|
| 1 | σ-Euler Model | `SigmaEulerModel.tex` | `chap: Sigma Euler Model` | unchanged |
| 2 | Vertical Multilayer Discretisation | `VerticalFESemiDiscretisation/*.tex` | `chap: Vertical Discretisation` | unchanged |
| 3 | Linearised Model | `LinearModel.tex` | `chap: linearised model` | unchanged |
| 4 | Model Verification | `ModelVerification.tex` | `chap: model verification` | unchanged |
| 5 | Stokes-Wave-Type Fourier Analysis | `StokesWaveFourierAnalysis.tex` | `chap: stokes fourier analysis` | unchanged |
| **6** | **Horizontal Weak Formulation and Regularity** | new `HorizontalFormulation/WeakForm.tex` (+ `BrokenAudit.tex`) | `chap: horizontal weak form` | **new: the mathematics of old ch. 6, incl. the 𝓛/𝓝 blocks and the term classification** |
| **7** | **Multi-field Residual Implementation** | `NumericalImplementation/GlobalResidual.tex` (reduced) | `chap: Gridap Global Residual` (kept) | **split from old ch. 6: layout, unit basis, time integration** |
| **8** | **Stability of the Unstabilised Discretisation** | `SolverValidation/StabilityAnalysis.tex` (revised) | `chap: stability analysis` (kept) | **revised, ends on "needs stabilisation"** |
| **9** | **Stabilisation** | new `SolverValidation/Stabilisation.tex` | `chap: stabilisation` | **new** |
| **10** | **Gridap Solver Implementation** | `NumericalImplementation/GridapImplementation.tex` (revised) | `chap: Gridap Implementation` (kept) | **moved to the end** |
| 11 | Validation of the BALFE-$M$ Solver | `SolverValidation/ValidationTests.tex` | `chap: Validation Tests` | placement proposed (§5, D1) |
| App. A | Projected Treatment of the Class-III Terms | `ClassIIIterms.tex` (after `\appendix`) | `chap: class III projection` (kept) | **decided: appendix** |

**Labels:** every existing label is kept when its content moves. `\ref`s therefore survive the
move, and only the chapter wording ("in the previous chapter", "chapter 7") needs editing. The
folder names (`NumericalImplementation/`, `SolverValidation/`) may be renamed later. That is
cosmetic and best done in a separate commit, so that the move and the rename stay distinguishable
in the history.

---

## 3. Why each chapter sits where it does — and what justifies the step

### Layer I–II — chapters 1–5 (unchanged)

The model is derived (1–3), checked against Yang & Liu coefficient by coefficient (4: worst
relative difference 1.2e-13 on the weighted momentum residual), and analysed as a model (5:
dispersion functional `R(μ)`, applicable ranges reproduced to < 1 %). None of it needs a
horizontal discretisation. After chapter 5 the model's accuracy is known independently of any
solver. Revisions 1–3 established this part; nothing in revision 4 touches it.

### Chapter 6 — Horizontal Weak Formulation and Regularity

**Question.** The model of chapter 2 is now discretised in the horizontal with continuous (`C⁰`)
Lagrange elements. Which of its terms are admissible on such spaces, what does each one demand,
and what must the weak form contain to be consistent?

**Why here.** It is the first step after the model is fully understood. Every later chapter
discretises *this* weak form, so its correctness must be settled first, and on mathematical
grounds only: no code, no runs.

**Why the 𝓛/𝓝 blocks and the term classification belong here, not in the implementation
chapter.** The regularity problem is not a property of the weak form as a whole. It sits in
*specific terms*: the components of the nonlinear pressure operator `𝓝` that carry second
derivatives of the trial fields. A reader can only see why the broken formulation is needed, and
which terms it changes, once every term has been written and classified by the derivatives it
puts on the trial and test functions. That classification is the core of the regularity argument.
Placing it in an implementation chapter would leave the argument without its evidence.

**Contents, in order.**
1. **The Galerkin residual on `Ω_h`: continuity, acceleration, gravity, advection.**
   * Horizontal FE spaces, test functions, and integration by parts of each block.
   * Source: old ch. 6 §"Weak Form system of Residuals", and the unused draft
     `HorizontalDiscretisation.tex` (§§ horizontal FE approximation, discretisation of the
     momentum equation, fully discrete weak form, tensor-product structure). That draft is
     commented out in `main.tex` and should be mined, then deleted or archived.
2. **The horizontal spaces: the Taylor–Hood pairing.**
   * `η` enters the momentum equation undifferentiated, through `∇·v`, so it plays the pressure
     role and the pairing must satisfy the inf-sup condition.
   * Source: old ch. 8 `subsec: horizontal FE spaces`, moved here because it is a property of the
     weak form, not of a stability experiment.
   * Justifying result: equal-order `Q2/Q2` grows a checkerboard at `λ ≈ 2–3·dx` and dies at
     t = 26.4 s (`dx` = 0.25); Taylor–Hood `Q2/Q1` runs 80 s flat (rule 12b).
3. **The linear pressure operator `𝓛`.**
   * The weak form of each block, and its admissibility: the derivatives it puts on trial and on
     test functions.
   * `R_P` integrated by parts is the whole frequency dispersion of the model (rule 1). Its volume
     part must be assembled; only its boundary part vanishes.
   * `𝓛` is first order in the trial fields once integrated by parts, so it is admissible on `C⁰`.
   * Source: old ch. 6 §"Implementation of `𝓛` pressure blocks" → "Admissibility of weak
     formulation terms" and the operator definitions. The Gridap realisation of these blocks goes
     to chapter 7.
4. **The nonlinear pressure operator `𝓝`: classification of its components.**
   * The eight components, grouped into Class I, II and III by their prefactor (`∇h`-type,
     surface-slope-type, leading-pressure-type) and by their order in the trial derivatives.
   * **Class III `{1,2,4,5}` carries Hessians of the trial fields.** The model carries `∂³u` and the
     weak form `∂²u` (rule 1b). These are the terms the rest of the chapter is about.
   * Source: old ch. 6 §"Implementation of `𝓝` pressure blocks" → "Classification of components",
     "Treatment of the Class-III terms", and the per-class treatments (integration by parts for the
     bed-slope prefactor; the surface-slope and leading-pressure terms).
5. **The term-by-term classification and the model switches.**
   * Every term of the full model, tagged by amplitude order × bed-slope class × activation
     condition. The table factorises: `regime` is a truncation in amplitude order, `flat_bed` a
     projection onto `∇h ≡ 0`, `nl_pressure` a component filter on `𝓝`. That is why the three
     switches are orthogonal, and the `O(ε)` rows *are* the linearised model.
   * It belongs here because it is a statement about the weak form (which terms exist in which
     model), not about code. Chapter 8 needs it to compare `:native` against `:full`.
   * Source: old ch. 7 §"Model setups" and §"Term-by-term anatomy". The *code* realisation of the
     switches (`resolve_physics`, the internal booleans) goes to chapter 7.
6. **The regularity problem.**
   * On `C⁰` elements, the distributional `∂²` of a trial field is the cellwise Hessian **plus a
     Dirac layer on the skeleton**. Cell quadrature silently drops the layer.
   * Raising the polynomial order does not help: Lagrange `Qp` is `C⁰` for every `p`.
   * Source: old ch. 6 §"Model Regularity: Broken Weak Formulation" (the author's section) and its
     subsections.
7. **The broken weak form and the term-by-term audit.**
   * Integration by parts element by element, going through the blocks of §§1–5 again, and the
     classification of facet terms:
     * kind 0 — identically zero on `C⁰`, e.g. `⟦v⟧`;
     * kind S — trial jumps, zero at the exact solution;
     * kind T — test-derivative jumps, mandatory.
   * Source: `BrokenAudit.tex` **minus its C⁰-IP subsection**, which moves to chapter 9.
   * Justifying results:
     * the cellwise-plus-layer form and the integrated-by-parts form are identical (gate G1,
       7e-12);
     * the layer is not negligible (G2);
     * the leading-pressure skeleton term `⟦𝒫⟧` must **not** be added: in the 1-D toy it makes
       `∂R/∂u̇` `h⁻⁴`-conditioned and non-convergent.
8. **Two consistent treatments of the Class-III gradients.**
   * **Broken:** the distributional gradient, assembled directly (cellwise Hessians plus one
     skeleton integral). 3 fields.
   * **Mixed:** `𝖦 ≈ ∇𝖲`, `𝖥 ≈ ∇𝖻` as auxiliary unknowns defined by integrated-by-parts weak
     equations. 5 or 7 fields. The mixed `𝖦` *is* the `L²` projection of the distributional
     gradient: it contains the skeleton layer, verified to 2e-16.
   * The reduction of components `{1,2,5}` to a single contraction in `∇𝖲`, exact to 4.4e-16.
     This is what makes both treatments affordable.
   * Source: old ch. 6 §"Mixed Formulation" (its weak form) and the Class-III treatment
     subsections.
   * A short paragraph on the projected treatment, as a non-consistent approximation, pointing to
     the appendix (D2).

**Conclusion of the chapter.** Every term of the model is classified by the regularity it demands.
Only the Class-III components exceed what `C⁰` spaces carry. The horizontal weak form of the exact
model is consistent on `C⁰` elements if and only if those components are written on the broken
domain. Two exact formulations exist (broken and mixed); the frozen projection is an approximation
of them. **Nothing is said about stability.**

⚠ **One sentence in the author's section must change for this chapter to be true.** The section
currently ends: *"the discrete system will only yield a stable and convergent scheme provided that
the broken weak form properly incorporates the required interior facet contributions."*
* Including the facet terms is **necessary for consistency**, but it is **not sufficient for
  stability**. The broken formulation includes them exactly, and `:full` broken still dies at 44 s
  in the closed box without a penalty (chapter 8).
* Suggested rewording: "...will only yield a **consistent** scheme provided that...; its stability
  is a separate question, examined in chapter 8."
* This sentence is the pivot of the whole chain. Leaving it as written would make chapter 8
  contradict chapter 6.

### Chapter 7 — Multi-field Residual Implementation

**Question.** How is the weak form of chapter 6 written in Gridap for an arbitrary vertical basis,
with no loops over vertical indices, and how is it advanced in time?

**Why here.** It turns chapter 6 into something that runs. It must precede chapter 8, which needs
a running solver, but it needs nothing from the boundary machinery: it contains only what the
closed periodic box uses. The *mathematics* of each term was settled in chapter 6; this chapter
only realises it.

**Order of the chapter.** First the layout is justified (§1). Then, immediately, every term is
written in that layout (§2), so the reader sees the layout pay off on the terms they already know
from chapter 6. Time integration follows (§3). The unit vertical basis closes the chapter (§4): it
is the computation that supplies the constant tensors every expression of §2 contracts against.
Until then the tensors can be taken as given, as chapter 2 defines them.

**Contents.**
1. **Multi-field block compact form; the stacked `[η, 𝖴x, 𝖴y]` layout.**
   * The vertical index lives in the FE value type: 3 fields, not `1 + 2Nσ`. The static vertical
     arrays become constant `TensorValue` / `ThirdOrderTensorValue` objects, and every vertical sum
     becomes a tensor contraction.
   * The three reasons for the layout: the layer index must be a value-type axis; fusing both
     indices would blow up the tensor rank; boundary conditions.
   * Source: old ch. 6 §"Multi-field block compact form", §"Field layout".
   * Justifying result: ~3.6× faster with ~13× fewer allocations than the per-layer form, and one
     residual serves both the sequential and MPI paths (CLAUDE.md §6).
2. **The terms in the stacked layout.** Every block classified in chapter 6, written as a stacked,
   loop-free Gridap residual:
   * the first-order building blocks;
   * the left-hand-side residuals (continuity, acceleration, gravity, advection);
   * the `𝓛` blocks;
   * the `𝓝` blocks, including the Class-III terms in their broken form (cellwise Hessians via
     `∇∇` plus the skeleton integral) and their mixed form (the 5/7-field layout and the Gram
     block);
   * the code realisation of the model switches (`resolve_physics`, one control point for
     `flat_bed`).
   * Source: old ch. 6 §"Implementation of Left-Hand-Side Residuals", §"Implementation of `𝓛`
     Residuals", §"Implementation of `𝓝` Residuals" and the mixed-formulation assembly.
   * The Gridap constraints that shape the code are stated here: never apply `∇` to an `Operation`
     containing a test basis (rule 6); the remarks on rest-state-safe gravity and the linearised
     reductions.
3. **Time integration.**
   * The semi-discrete system `M(u)u̇ + K(u) = 0` as a Gridap transient operator.
   * The integrators: Crank–Nicolson, `SDIRK_2_2`, explicit RK. Explicit RK cannot run the mixed
     formulation, because its `∂R/∂u̇` is singular on the auxiliary rows.
   * The Newton loop and its Jacobians:
     * `∂R/∂u̇` is exact in every model;
     * `∂R/∂u` is quasi-Newton for the mixed treatment, and exact for the broken one
       (`broken_class3_jacobian`).
   * Justifying results: `test_jacobians_ad` 17/17 over 8 models; broken Class-III Jacobian
     against finite differences 3.5e-11 (G10).
   * Rule 17b, stated here once: the Jacobian changes the cost of a step, never the converged
     answer.
   * Moved here from old ch. 7 §"Solver algorithm workflow" because chapter 8 depends on it: the
     stability verdict is integrator-dependent (rule 15). The *analysis* of the integrators
     (amplification factors, `SDIRK_2_2` as `DIRK22(1,0,1)`) stays in chapter 8.
4. **The unit vertical basis `{φⱼ}` and the vertical tensors** (end of chapter).
   * The unit basis BVP on the σ-mesh, then the Gram, `𝓜`, `𝓖`, `𝓐`, `𝓚` and `𝓟` tensors that §2
     contracts against.
   * Source: old ch. 7 §"Unit vertical basis `{φⱼ}` calculation".

**Conclusion of the chapter.** A loop-free, verified realisation of the weak form of chapter 6 for
any vertical basis, integrated in time and runnable in a closed periodic domain.

### Chapter 8 — Stability of the Unstabilised Discretisation

**Question.** Does the discretisation of chapters 6–7 keep perturbations bounded? If not, where
does the growth come from?

**Why here.** Only now does a running, consistent discretisation exist to test. Stability must
precede validation (Lax–Richtmyer: a consistent scheme converges if and only if it is stable), and
it must precede any change to the formulation: an approximation can only be judged against the
stability of the exact discretisation (revision 2's lesson).

**Contents (the existing chapter, revised).**
1. **What is being asked**
   * the three kinds of growth: physical, ill-posedness, numerical;
   * semi-discrete versus fully discrete stability;
   * the toolbox table.
   * Unchanged.
2. **The reference discretisation.**
   * Why Yang & Liu's five-point finite differences with RK4 are stable: a five-point stencil
     capped at `1.37/Δx` with a null at `2Δx`, plus a Shapiro filter.
   * Unchanged. It sets up the contrast.
3. **The FE discretisation spectrally.**
   * Spectral reach: `C⁰` elements represent 3.5–4.4/Δx, beyond the node Nyquist, with no null.
   * Dispersion saturation.
   * The equal-order episode, kept as a short worked example: the method works.
   * The Taylor–Hood subsection moves to chapter 6 and leaves a pointer.
4. **The integrator.**
   * `SDIRK_2_2` removes ≈ `0.75 (ωΔt)⁴` of a neutral mode's energy per step. That explains every
     carrier decay of 0.010–0.012 s⁻¹ and masks the instability.
   * **Every stability statement from here on is made under Crank–Nicolson.**
   * Justifying result: Q2/Q1 `:full` at 32 cells/λ is bounded for 100 periods under `SDIRK_2_2`,
     but diverges at ≈ 90 s under Crank–Nicolson.
5. **The closed periodic box: the main experiment of the chapter.**
   * The test: one wavelength, no boundaries; the continuum answer is "nothing grows".
   * Results, under Crank–Nicolson:
     * mixed `:full`, Q2/Q1: diverges at 157 / 104 / 90 / 31 s at 8 / 16 / 32 / 64 cells/λ. Onset
       advances with every refinement: the grid-scale signature (rule 38b).
     * projected `:full`: dies at 8.8 s;
     * **broken `:full`, Q3/Q2: dies at 44 s.**
     * `:native`: diverges at 64 cells/λ (Q2/Q1, 146 s) and grows at Q3/Q2 with 16 cells/λ
       (+0.10 s⁻¹).
   * The broken result is new and must be added. Together with the mixed one, it shows that
     **both consistent formulations of chapter 6 are unstable**. The instability is therefore not
     an artefact of one Class-III treatment.
6. **Linearised stability about a frozen state.**
   * Frozen growth rates are not predictive. The spectral radius is: the Class-III Doppler branch
     `ω ≈ kU + ω∞`, with `ρ ∝ A/h_e`.
   * Unchanged.
7. **Synthesis: the mechanism and the verdict.**
   * The nonlinear terms, Class III above all, transfer carrier energy up the wavenumber ladder at a
     rate proportional to `kU`.
   * `C⁰` elements represent that band with no null and no dissipation, so nothing removes what
     arrives there.
   * Integrator damping acts on frequency, not wavenumber, and vanishes as `Δt → 0`, so it is a
     mask, not a remedy.
   * **The chapter ends with the verdict: the discrete operator of the full nonlinear model is
     consistent but unstable, and `:native` shares the instability more weakly. A
     wavenumber-selective energy sink is required.**
   * Today's §"Remedies the analysis admits" becomes that closing statement plus a forward
     reference. The remedies themselves move to chapter 9.

**What leaves the chapter.**
* **The flume phenomenology** (§"The nonlinear instability of the full model": the experiment and
  the flume phenomenology) needs wave generation and absorption, which are now introduced in
  chapter 10. See decision D3.
* **The projected-treatment stability record** is already in `ClassIIIterms.tex`.

### Chapter 9 — Stabilisation

**Question.** Which wavenumber-selective sink, with what strength, and does it remove the
instability without degrading the carrier?

**Why here.** Chapter 8 proves a sink is needed and states the rate it has to beat. A remedy can
only be designed against a measured requirement: here, mid-band damping greater than the transfer
rate of +0.13 s⁻¹ (A = 0.10), with carrier loss ≪ 1/T_run.

**Contents.**
1. **The design requirement**, taken from chapter 8: scale selectivity, i.e. strong damping in the
   mid and high bands and negligible damping of the carrier. Plus the four acceptance criteria:
   * damping margin;
   * carrier preserved over 100 periods;
   * convergence orders preserved;
   * the amplitude ceiling stated.
2. **The candidates, with pros and cons:**
   * **Explicit filters** (a Shapiro analogue): what Yang & Liu do; not tried here.
   * **Integrator damping** (generalised-α with `ρ∞`): a mask, not a remedy (chapter 8 §4).
   * **First-order C⁰-IP / CIP edge stabilisation** (Burman–Hansbo).
     * Consistent, linear, exact Jacobian.
     * **But its null space is the C¹ subspace.** Linear eigen-analysis: oscillatory modes retreat
       into it undamped, with damping ≤ 0.016 s⁻¹ that *falls* like 1/γ. It is a constraint, not a
       damper.
     * Box results: removes the grid-scale growth, but only **delays** `:full`, with late mid-band
       growth at +0.13 s⁻¹ from ≈ 110 s.
     * Velocity-only penalties fail.
     * **A Q1 surface cannot be penalised**: the null space is empty, so the carrier locks
       (−0.021 s⁻¹). Q2/Q1 has no γ window at all, so `:full` requires Q3/Q2.
     * Source: `BrokenAudit.tex` §"The C⁰ interior penalty weak form", plus
       `BROKEN_FORMULATION_PLAN.md` §4.
   * **High-order jump penalty** (`:jumpgrad`, orders ≤ 2; Burman–Ern hp scaling).
     * Shrinks the null space to C² functions.
     * Limited by Gridap to FE derivatives of order ≤ 2, so at Q3 the third order is missing.
   * **Direct (volume) ghost penalty** (`:ghostvolume`; Preuß 2018, Lehrenfeld–Olshanskii 2019).
     * The jump of each cell's polynomial extension into its neighbour.
     * Equals the full jump family to every order, with cross terms: identity (★), verified to
       3e-16. So it is the complete high-order penalty, written without derivatives.
     * Implementation via shifted cell polynomials (`AffineField`). Its present limits are a
       uniform Cartesian mesh and a sequential run only.
3. **Choosing γ: the linear eigen-analysis as a design tool.**
   * The γ window per method and pairing:
     * jump penalty, order ≤ 2: [1e-3, 3e-3], with mid-band damping 0.27–0.82 s⁻¹;
     * ghost penalty: [0.01, 0.03], with mid-band damping 0.85–2.55 s⁻¹.
   * Q2/Q1 has no window.
   * The measured carrier damping matches the prediction to ±8 %, so the tool is reliable.
   * γ is dimensionless, scaled with `τ` and `h^(s−2)`. It is calibrated per discretisation
     family, not per run.
4. **Results in the closed box** (Q3/Q2, Crank–Nicolson, 100 periods):

   | arm | result | carrier energy lost over 160 s |
   |---|---|---|
   | jump penalty, order ≤ 2, γ = 2e-3: `:full` and `:native` | ✅ 160 s | ≈ 2 % |
   | ghost penalty, γ = 0.0033 / 0.01 / 0.03, `:full` | ✅ 160 s each, decaying | 0.5 / 1.5 / 4.3 % |
   | ghost penalty, γ = 0.01, 32 cells/λ: `:full` and `:native` | ✅ 160 s | — |
   | ghost penalty, `A` = 0.15 | decaying at 137 s; final result pending | — |
   | controls in the same batch | first-order penalty: ×100–250 growth; `:native` with γ = 0: ×920 | — |

5. **The choice of method, and what remains.**
   * The choice and its reasons.
   * Still to do:
     * the MMS orders at γ\*;
     * the flume;
     * bathymetry;
     * the amplitude ladder;
     * cost (assemble the constant penalty matrix once);
     * general meshes and the distributed path for the ghost penalty.
   * These are written as they are completed.

**Conclusion of the chapter.** The discretisation of the exact model is stable in the closed box
once a scale-selective penalty is added; the remaining items extend the claim beyond the box.

### Chapter 10 — Gridap Solver Implementation

**Question.** What else does a usable solver need, beyond the stabilised core?

**Why here.** Everything in it is a *functionality* added around a core that is now known to be
consistent (chapter 6), efficiently assembled (chapter 7) and stable (chapters 8–9). Introducing
it earlier would put boundary mechanisms into the stability analysis. That is exactly what made
two weeks of flume factorials inconclusive.

**Contents.**
* Wave generation: the internal source, and Dirichlet boundary generation with WaveSpec coupling.
* Absorption: sponge and relaxation zone.
* Boundary conditions.
* The augmented residual.
* The summary of the implemented global residual, including the stabilisation term.
* The solver algorithm workflow, minus what moved to chapter 7.
* Software architecture and the distributed path, which currently implements only the projected
  treatment.
* **The flume demonstrations of the stabilised solver**: the first-order-penalty flume held
  η_max = 0.110 flat from 60 s to 160 s, where the unpenalised flume died at 38 s. This is the
  natural home for the flume work (D3).

### Chapter 11 — Validation (placement proposed)

The validation campaign as it exists:
* semi-analytical dynamics;
* MMS models 1–6;
* grid and time convergence;
* vertical reconstruction;
* physical benchmarks against Yang & Liu.

It needs both the stabilised solver (chapter 9) and the boundary machinery (chapter 10: the
generation gates, the physical benchmarks), so it comes last.

⚠ The existing MMS results were measured **without** a stabiliser. Chapter 11 must state that, and
add the MMS at γ\* once it is run. Until then, the stabiliser's effect on the convergence orders
is the open acceptance criterion of chapter 9.

---

## 4. Migration map — where every current section goes

| current location | section | goes to |
|---|---|---|
| `GlobalResidual.tex` | Multi-field block compact form | **7** |
| `GlobalResidual.tex` | Weak Form system of Residuals (continuity, momentum) | **6** §1 |
| `GlobalResidual.tex` | Model Regularity: Broken Weak Formulation (+ subsections) | **6** §6–7 (reword the stability sentence) |
| `BrokenAudit.tex` | Setting … Class III, Summary | **6** §7–8 |
| `BrokenAudit.tex` | The C⁰ interior penalty weak form | **9** §2 |
| `GlobalResidual.tex` | Field layout `[H, 𝖴x, 𝖴y]` | **7** §1 |
| `GlobalResidual.tex` | Implementation of LHS residuals; building blocks | **7** §2 |
| `GlobalResidual.tex` | `𝓛` blocks: operator, admissibility | **6** §3; the Gridap residual subsections go to **7** §2 |
| `GlobalResidual.tex` | `𝓝` blocks: classification of components, Class I/II/III treatment | **6** §4 (classification and treatment) / **7** §2 (the `𝓝` residual implementation) |
| `GlobalResidual.tex` | Mixed Formulation | **6** §8 (weak form) / **7** §2 (assembly) and §3 (Jacobian blocks) |
| `HorizontalDiscretisation.tex` (unused draft) | horizontal FE approximation, fully discrete weak form, tensor-product structure | mine into **6** §1, then retire |
| `GridapImplementation.tex` | Wave generation; Absorption; Boundary conditions; Augmented residual; Summary | **10** |
| `GridapImplementation.tex` | Solver algorithm workflow | integrators and Newton → **7** §3; the rest → **10** |
| `GridapImplementation.tex` | Software architecture | **10** |
| `GridapImplementation.tex` | Model setups; term-by-term anatomy | **6** §5 (classification, switches as a statement about the weak form) / **7** §2 (code realisation) |
| `GridapImplementation.tex` | Unit vertical basis calculation | **7** §4 (end of chapter) |
| `StabilityAnalysis.tex` | §8.1 what is asked; reference discretisation | **8** |
| `StabilityAnalysis.tex` | FE horizontal discretisation: Taylor–Hood subsection | **6** §2 (pointer left behind) |
| `StabilityAnalysis.tex` | FE horizontal discretisation: the other subsections | **8** |
| `StabilityAnalysis.tex` | The nonlinear instability of the full model (flume) | **10** or appendix (D3), with a motivating paragraph left in 8 |
| `StabilityAnalysis.tex` | Closed periodic flume; frozen linearisation | **8** (+ the broken `:full` 44 s result) |
| `StabilityAnalysis.tex` | Synthesis: mechanism, what is established | **8** §7 |
| `StabilityAnalysis.tex` | Remedies the analysis admits | **9** §2 (8 keeps the verdict) |
| `ValidationTests.tex` | everything | **11** (D1) |
| `ClassIIIterms.tex` | everything | **Appendix A** |
| `BROKEN_FORMULATION_PLAN.md` §4–5, `GHOST_PENALTY_PLAN.md` §5 | the stabilisation campaign and results | **9** §2–5 (new LaTeX text) |

---

## 5. Decisions needed from the author

* **D1 — Validation placement.**
  * **Recommended: last (chapter 11).** It needs both the stabilised core and the boundary
    machinery.
  * Alternative: between 9 and 10, splitting off the boundary-dependent gates. That cuts a
    coherent chapter in two for little gain.
* ~~**D2 — The projected treatment.**~~ **Decided (2026-10-05): Appendix A.** Its old
  justification for existing, a cheaper 3-field system, is gone now that the broken formulation is
  an exact 3-field system. It is kept because the distributed path still uses it, and because its
  stability record is part of the history (projected box: 8.8 s).
  * In `main.tex`: `\appendix` before `\input{ClassIIIterms}`, so its label survives and its
    numbering becomes "A".
  * Its chapter introduction is reworded from "a change of formulation, judged last" to "a record
    of the first treatment implemented, and the one the distributed path still uses".
  * Chapter 6 §8 mentions it in one paragraph; chapter 10 points to it from the distributed path.
* **D3 — The flume instability record.**
  * The flume campaigns (§5.2b–e of CLAUDE.md) are where the instability was first seen, but they
    cannot attribute it.
  * **Recommended:** keep one motivating paragraph and one figure in chapter 8 §5 ("first seen on
    the flume, always pinned at the inflow; the box removes the boundary to test the interior"),
    with a forward reference. Move the detailed record to chapter 10, beside the stabilised flume
    results, or to the appendix with the projected record, since most of those runs were
    projected or mixed.
* ~~**D4 — Unit vertical basis computation.**~~ **Decided (2026-10-05): the end of chapter 7
  (§4)**, after the terms that consume its tensors.
* **D6 — Where the model switches live.** Decided with D4: the term classification and the switches,
  as a statement about which terms exist in which model, go to **chapter 6 §5**. Their code
  realisation (`resolve_physics`, the internal booleans) goes to **chapter 7 §2**. Revisit this if
  splitting the table across two chapters reads badly.
* **D5 — Chapter titles.** The working titles above are descriptive. They can be shortened once
  the content settles.

---

* **D7 — The `:native` / `:full` tiers (open, 2026-10-05).**
  * `:native` is defined by a *numerical* criterion: the `𝓝` components that are first order on
    `C⁰`. It is not a physical one, and no ordering in amplitude or `kd` separates it from `:full`.
  * Proposal: replace `nl_pressure ∈ {:none, :native, :full}` with a Boolean (all eight components
    on or off). Keep a diagnostic mask for the stability chapter's "Class III removed" control.
  * If adopted, chapters 6 §5 and 8 must present `:native` as a diagnostic decomposition, not as a
    model.
  * ⚠ Related: every `:full` box and flume run so far used `C3_MASK=gs`, i.e. **component 4 (`∇𝖻`)
    omitted**. Before chapter 9 can claim "the full model is stabilised", the penalised arms must be
    repeated with all eight components.

## 6. Loose ends carried forward

* Chapters 7 and 9 (old numbering) say "`Q1` zeroes the B-dispersion term". This was never
  re-derived, and is not obviously true in 2-D for the mixed derivatives.
* `ModelVerification.tex:113` references `subsec: comparison point derived`, which is undefined.
  It is the only undefined reference in the build.
* `\Gridap` is undefined at line 239 of the author's regularity section. This is a compile error
  that predates the broken audit.
* Chapter 8's periodic-box table needs the broken `:full` result (Q3/Q2, Crank–Nicolson, died at
  44 s). Chapter 10 / the appendix needs the RK4 paragraph transcribed (Q2/Q1 40.2 s, Q3/Q2 15.8 s,
  against `SDIRK_2_2` 40.4 / 15.6 s).
* `CLAUDE.md` §2 describes the old chapter order. Update it when revision 4 is applied, not before.

---

## History

* **Revision 1 (2026-09-25).** The Overleaf restructure grouped the document into layers:
  * model (1–3);
  * checking the model with no solver (4–5);
  * implementation (6–7);
  * assessment (8–9).

  Its principle: each chapter uses only what precedes it. Stability comes before validation
  (Lax–Richtmyer).
* **Revision 2 (2026-09-25).** The implementation chapters presented the **mixed** formulation as
  *the* treatment, and the projected one moved, intact, to a final chapter 10. The reason: the
  projection had been implemented before any stability analysis of the exact model, so a year of
  stability results described the approximation rather than the model. Lesson: an approximation of
  the model belongs *after* the analysis of the exact discretisation.
* **Revision 3 (2026-09-25/26).** A consistency pass:
  * a Taylor–Hood section was added;
  * 11 dangling references were retargeted;
  * stale statements were corrected;
  * chapter 8 was expanded into a teaching chapter with "method at a glance" boxes, the
    `SDIRK_2_2` tableau correction and nine figures.
* **Revision 4 (2026-10-05, this file, proposed; amended the same day).** The split described above. The
  amendment moved the `𝓛`/`𝓝` blocks and the term classification into chapter 6, leaving chapter 7 the
  stacked-layout realisation, the unit vertical basis and time integration. It follows from:
  * the broken audit, which made the regularity question a chapter of its own;
  * the finding that both consistent formulations are unstable, so stability is a separate question
    from regularity;
  * the stabilisation campaign, which made the remedy a chapter of its own.
