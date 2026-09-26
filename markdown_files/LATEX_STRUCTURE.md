# LATEX_STRUCTURE.md — the chapter order of `latex_docs/BALFEM_models/` and why

*Recorded 2026-09-25, after the restructure made in Overleaf (`BALFEM_models` commit `4e1cef9`,
"Update on Overleaf"), which moved `StabilityAnalysis.tex` and `ValidationTests.tex` into a new
`SolverValidation/` folder and placed both after the two implementation chapters.*

> ## ⚙ REVISION 2 (2026-09-25): THE SOLVER IS IMPLEMENTED WITH THE MIXED FORMULATION; THE PROJECTION MOVES TO THE END
>
> **What changed.** In the first structure, the implementation chapters (6–7) presented the
> Class-III terms `{1,2,4,5}` as treated by **frozen `L²` projections**: `𝖲 = ∇·(H𝗎)` and
> `𝖻 = 𝗎·∇H` projected onto a `C⁰` space once per accepted step and differentiated during the next,
> with the mixed formulation mentioned only as a rejected, "too expensive" alternative. From
> revision 2 on:
> * chapters 6–7 present the **mixed formulation** as *the* implementation: `𝖦ₐ ≈ ∂ₐ𝖲` and
>   `𝖥ₐ ≈ ∂ₐ𝖻` are genuine FE unknowns defined by their own weak equations (integration by parts
>   built in, boundary term assembled), solved monolithically with `[η, 𝖴x, 𝖴y]` at every Newton
>   iteration;
> * every description, equation, workflow step, table row, test and stability result that belongs
>   to the **projected** treatment is **moved, not deleted**, to a new final chapter,
>   `ClassIIIterms.tex` (chapter 10, after the validation chapter), and replaced in its original
>   place by the mixed-formulation counterpart or by a pointer to chapter 10.
>
> **Why — the error this corrects.** The projection was implemented *first*, directly, as the way to
> make the inadmissible second derivatives computable — **without the step that has to precede any
> approximation of that kind: a stability analysis of the discretisation of the exact model.** The
> projection is not the exact model; it adds two approximations (a one-step lag, and the
> project-then-differentiate recovery of a derivative of a discontinuous field), and its stability
> was never established before it became the production path. The stability study of chapter 8 —
> the current work — has since been carried out on the **mixed** formulation precisely because it
> solves the exact derived BALFE-M system (after discretisation) with no projection and no added
> approximation. It has shown the projection to be far less stable than the mixed formulation in a
> closed periodic box (projected `:full` dies at t = 8.4 s, where the mixed one runs 100 periods at the
> same resolution) — and, as of 2026-09-26, that **the mixed formulation is itself unstable** at finer
> meshes, higher order and larger amplitude once the integrator's damping is removed: the discretisation
> of the exact model needs stabilisation (`OPEN_ISSUES.md` §0e). That finding is only attributable to
> the discretisation *because* the analysis was done on the exact model.
>
> **The natural progression the document now follows** is therefore:
> 1. **implement the exact model** in the solver — the mixed formulation (chapters 6–7);
> 2. **analyse its stability** and the issues that appear (chapter 8), and **validate** it (chapter 9);
> 3. only then, **potentially change the formulation** — e.g. project the Class-III terms because it is
>    cheaper (no auxiliary unknowns, a 3-field instead of a 7-field system), or for any other reason —
>    and assess what that change costs in accuracy and stability (chapter 10).
>
> An approximation of the model belongs *after* the analysis of the model's discretisation, because
> only that analysis can tell whether the approximation is benign. Placing the projection first made
> every later result about `:full` a statement about the approximation, not about the model — which is
> exactly how a year of stability work ended up diagnosing the projection.
>
> **What still uses the projection, stated where it applies.** The distributed (MPI) path implements
> only the projected treatment — the mixed formulation is sequential only — and the `:full` MMS
> campaign and several `:full` tests were run with projections. Those passages move to chapter 10;
> where the chapters 6–9 must mention them, they point there.

## The principle

**The document follows the development of the project, step by step: each chapter uses only what
the chapters before it have built.** A result is placed at the earliest point where the machinery it
needs exists, and no earlier. The machinery grows in five layers, and the chapters are grouped
accordingly:

| layer | what exists at that point | chapters |
|---|---|---|
| **I. The model** | continuous equations, the vertical FE basis and its tensors | 1–3 |
| **II. Checking the model** | the same, plus published equations and analytic wave theory — **no horizontal discretisation, no solver, no time integration** | 4–5 |
| **III. Implementing the model** | the horizontal FE discretisation in Gridap, the multi-field residual, the solver | 6–7 |
| **IV. Assessing the implementation** | a running solver: time integration, runs, convergence studies | 8–9 |
| **V. Changing the formulation** | everything above, including the stability analysis of the exact model | 10 |

The two questions a reader must be able to keep apart — *is the model right?* and *is the solver
right?* — are answered in separate blocks (II and IV), with the implementation between them. A
defect found in block IV can therefore be attributed to the implementation, because block II has
already established the model it implements.

## The order

| # | chapter title | file | label |
|---|---|---|---|
| 1 | σ-Euler Model | `SigmaEulerModel.tex` | `chap: Sigma Euler Model` |
| 2 | Vertical Multilayer Discretisation | `main.tex` + `VerticalFESemiDiscretisation/{VerticalFEapprox, wDerivation, pDerivation, VerticalProjection, VerticalSemiDiscreteSystem}.tex` | `chap: Vertical Discretisation` |
| 3 | Linearised Model | `LinearModel.tex` | `chap: linearised model` |
| 4 | Model Verification | `ModelVerification.tex` | `chap: model verification` |
| 5 | Stokes-Wave-Type Fourier Analysis | `StokesWaveFourierAnalysis.tex` | `chap: stokes fourier analysis` |
| 6 | Gridap Multi-field Residual | `NumericalImplementation/GlobalResidual.tex` | `chap: Gridap Global Residual` |
| 7 | Gridap solver Implementation | `NumericalImplementation/GridapImplementation.tex` | `chap: Gridap Implementation` |
| 8 | Stability Analysis | `SolverValidation/StabilityAnalysis.tex` | `chap: stability analysis` |
| 9 | Validation of the BALFE-$M$ Solver | `SolverValidation/ValidationTests.tex` | `chap: Validation Tests` |
| 10 | Projected Treatment of the Class-III Terms | `ClassIIIterms.tex` | `chap: class III projection` |

The folder names carry the grouping: `VerticalFESemiDiscretisation/` is layer I, `NumericalImplementation/`
is layer III and `SolverValidation/` is layer IV.

## Why each chapter sits where it does

### Layer I — the model (chapters 1–3)

**1. σ-Euler Model.** The starting point: the governing equations mapped onto the σ-coordinate. Nothing
precedes it because nothing else is needed.

**2. Vertical Multilayer Discretisation.** The vertical FE approximation `u_h = Σ_j u_j(x,t) φ_j(σ)`,
the analytical elimination of `w` and `p_nh`, the vertical projection that produces the structural
tensors, and the resulting semi-discrete system in `(H, u_1 … u_Nσ)`. This is the **full BALFE-M model**:
fully nonlinear, variable bathymetry, arbitrary vertical basis. It is discrete in the vertical only; the
horizontal direction is still continuous.

**3. Linearised Model.** The linearisation of chapter 2. It comes directly after the full model because
it is derived *from* it, and before any check because the first checks (chapter 4, and the dispersion
analysis of chapter 5) are stated on the linearised system.

### Layer II — checking the model, before any solver exists (chapters 4–5)

**4. Model Verification.** The empirical verification that the derivation of chapters 2–3 is correct:
at `p = 1`, BALFE-M collapses onto Yang & Liu's published LFE-M coefficient by coefficient — the linear
coefficients, the nonlinear vertical velocity, the non-hydrostatic pressure and the weighted momentum
residual, to round-off. It compares **operators at prescribed states**, so it needs neither the
horizontal discretisation nor a solver nor any time integration. It therefore belongs immediately after
the derivation it verifies, and before anything is built on that derivation.

**5. Stokes-Wave-Type Fourier Analysis.** Almost fully analytical: the Stokes–Fourier hierarchy, the
dispersion functional `R(μ) = Φᵀ(M + μ|B|)⁻¹Φ` and its basis-independent properties, group velocity,
shoaling, the second-order bound harmonic and the vertical grid optimisation. It depends **exclusively
on the vertical tensors** — there is still no horizontal discretisation. It closes the study of the
model as a model: after it, the model's accuracy is known independently of any solver.

### Layer III — implementing the verified model (chapters 6–7)

Only after the model has been derived, verified and analysed is a solver needed. These two chapters
describe how the model of chapters 2–3, verified in chapter 4, becomes a new solver in Gridap.

**6. Gridap Multi-field Residual.** The horizontal discretisation enters here for the first time: the
multi-field FE spaces, the single scalar global residual Gridap consumes, the stacked `[η, 𝖴x, 𝖴y]`
layout, the separation of the pressure operators into components, and the treatment of the Class-III
terms (projection or mixed unknown).

**7. Gridap solver Implementation.** Everything else the solver needs to run: wave generation (internal
source, Dirichlet boundary generation), absorption (sponge, relaxation zone), boundary conditions, the
augmented residual, the software architecture and workflow, the model switches and the term-by-term
classification, and the unit vertical basis computation. It follows chapter 6 because it builds on the
residual defined there.

### Layer IV — assessing the implementation (chapters 8–9)

With a solver in hand, the natural next step is to assess it. Both chapters need a running solver —
time integration and actual runs — which is why they come after layer III and not before.

**8. Stability Analysis.** Whether the discretisation of chapters 6–7 keeps perturbations bounded, for
the fully nonlinear model with variable bathymetry: why Yang & Liu's finite-difference/RK4 scheme is
stable, what a `C⁰` finite-element discretisation changes spectrally, the observed nonlinear
instability, the projected-versus-mixed treatment of the Class-III terms, and the linearised stability
analysis about a frozen state. It assesses the Class-III treatments of chapter 6 and the boundary
mechanisms of chapter 7, so it must follow both.

**9. Validation of the BALFE-M Solver.** The validation campaign: semi-analytical dynamics against
closed-form results, the manufactured-solution convergence studies, grid and time convergence,
vertical-profile reconstruction, and physical benchmarks against Yang & Liu.

**Why stability precedes validation.** By the Lax–Richtmyer equivalence theorem, a consistent
discretisation converges if and only if it is stable. The validation campaign measures convergence and
agreement; those measurements are meaningful only for a configuration known to be stable, and the
stability chapter is what states which configurations are (and which, such as `:full` with the
projected Class-III treatment, are not). Reading stability first also tells the reader why the
validation campaign runs the tiers it runs.

### Layer V — changing the formulation (chapter 10)

**10. Projected Treatment of the Class-III Terms.** The frozen-`L²`-projection alternative to the mixed
formulation: its definition and implementation, its place in the solver workflow and in the
distributed path, its verification record, and its stability record (the projected-versus-mixed
factorial, the refinement factorial, the periodic-box failure). It comes last because it is a
*change* of the formulation of chapters 6–7, and a change of formulation can only be judged against
the stability analysis and validation of the exact one (chapters 8–9). See REVISION 2 above.

## Loose ends

*Resolved in revision 2 (2026-09-25):* the stale comment above `\input{SolverValidation/StabilityAnalysis}`
in `main.tex` and the header comment of `StabilityAnalysis.tex` (both said "directly after the
multi-field residual") were rewritten; the abstract now mentions the mixed formulation and chapter 10.

*Resolved in revision 3 (2026-09-25, global consistency pass, uncommitted in `BALFEM_models`):*
* **The Taylor–Hood pairing was used from chapter 8 on but never introduced.** New section
  `sec: horizontal FE spaces` in chapter 6, placed between the layout choice and the LHS residuals.
  It covers the pairing, why `η` sits one order lower (the inf-sup condition), and the order of the
  auxiliary spaces. Chapters 7–9 now point back to it rather than restating it.
* The user's restructure of chapter 6 had removed `sec: pressure operator implementation`, leaving
  11 dangling refs. Each was retargeted by context, and `\subsubsec` (a compile error) is now `\paragraph`.
* Stale or wrong statements corrected:
  * the chapter 7 sponge remarks said η is not damped, while the equation damps it;
  * `GridapLFEM` / `p_horizontal` → `GridapBALFEM` / `p_u`;
  * the Level-3 MMS rates in chapter 9 (η is `h^p`, not `h^{p+1}`);
  * the chapter 9 `A ≤ 0.001` gotcha (an equal-order artefact);
  * chapter 4's claim that MMS establishes the Class-III assembly (it is not yet run on the mixed `:full`);
  * the claim in chapter 8 that SDIRK explains the carrier's energy loss (it is 100× too small).
* Added back-links from the implementation to the derivation:
  * momentum terms → the vertical projection;
  * building blocks → the `p_nh` derivation;
  * the unit basis → its BVP;
  * the Doppler branch → chapter 5;
  * the Newton gate → the MMS model-2 section.

*Revision 3b (2026-09-26): chapter 8 expanded for readability.* The chapter was too dense, so
it now teaches its methods as well as reporting them.
* A new §8.1 comes first. It covers:
  * the three kinds of growth: physical, ill-posedness and numerical;
  * semi-discrete versus fully discrete stability;
  * a toolbox table with, for each method, the question it answers, what it assumes and what it
    cannot see.
* Each method section opens with a "method at a glance" box: question, hypotheses, procedure, how to
  read the outcome, blind spots. The boxes cover the Fourier analysis, the FE spectrum, the flume
  experiments, the periodic box and the frozen linearisation.
* Added content:
  * a von Neumann derivation;
  * a table of amplification factors for CN, SDIRK_2_2 and RK4;
  * the Bloch-branch and min–max explanation of the FE spectral reach;
  * a table predicting how onset moves for each kind of growth;
  * the periodic box split into why / construction / null / measurement / outcome table / results;
  * Floquet versus frozen eigenvalues, with the two textbook examples;
  * the energy method and non-normality;
  * a synthesis table of what each method contributed.
* Section order is unchanged: observation → isolation (periodic box) → explanation (linearisation).
  Every old label is kept.

*Still open:*
* Chapters 7 and 9 say "`Q1` zeroes the B-dispersion term". This was not re-derived in revision 3
  and was left unchanged; it is not obviously true in 2-D for the mixed derivatives.
* **Pre-existing, unrelated to the restructure:** `ModelVerification.tex:113` references
  `subsec: comparison point derived`, a label defined nowhere — the only undefined reference in the
  build (unchanged since `BALFEM_models` commit `36ae631`).
* Chapter 8's periodic-box table now carries the SDIRK_2_2 campaign and the first Crank–Nicolson
  repeats (2026-09-26); the remaining Crank–Nicolson cases (`run_1dper_batch_cn{2,3}.sh`) must be added
  when they finish, together with the figure cells. Chapter 10's RK4 paragraph: both RK4 runs have
  finished (Q2/Q1 40.2 s, Q3/Q2 15.8 s, against SDIRK 40.4 / 15.6 s) — to be transcribed.

*Resolved 2026-09-26:* `CLAUDE.md` §2 now describes `SolverValidation/`, `StabilityAnalysis.tex` and
`ClassIIIterms.tex`; `INDEX.md` lists this file; chapter 8 corrected for Gridap's `SDIRK_2_2` tableau
(A-stable, not L-stable, ≈ `0.75(ωΔt)⁴` per step) and gained nine figures (plus two in chapter 5).
