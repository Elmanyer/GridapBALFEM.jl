# LATEX_STRUCTURE.md — the v2 LaTeX document (`latex_docs/GridapBALFEM_documentation/`): structure and state

*Condensed 2026-10-09; locations and chapter numbers updated 2026-11-24. The full reasoning behind the
order (revision 4, with the migration map and history) is
[`archive/LATEX_STRUCTURE_rev4.md`](archive/LATEX_STRUCTURE_rev4.md). The structure below is
**applied**.*

**Where it lives (since 2026-11-23).** `latex_docs/GridapBALFEM_documentation` is a link to
`../LATEX_DOCS/GridapBALFEM_documentation`, its own git repository whose remote `origin` is the
Overleaf project. The author edits it in Overleaf: always `git fetch origin` and fast-forward before
editing (`CLAUDE.md` §2). Chapters 1–5 are the folder `BALFEM_derivation/`, taken from the separate
Overleaf project `../LATEX_DOCS/BALFEM_models_derivation` (shared with the sister solver's
`GridapBALFEMixed_documentation`); edit them there. `../LATEX_DOCS/BALFEM_models_v2/` is the
superseded combined document.

## 1. The principle

The document follows the development of the solver: each step builds something, analyses it, and
the next step fixes what the analysis finds. A chapter uses only what precedes it. Chapters 6–9 form
one chain of four questions, and **keep consistency apart from stability** (the confusion that cost a
year, rule 39b):

```
 6  Is the weak form well defined on C⁰ elements?    → only on the broken domain; the broken formulation
 7  How is it written in Gridap and advanced in time? → stacked layout, Newton, unit vertical basis
 8  Which reduced models does the solver run?         → three switches, six models, term activation
 9  Is the discretisation stable?                     → NO: consistent but unstable at the grid scale
10  What sink, and does it work?                      → open (STATUS.md, STABILITY.md §6b)
```

The closed periodic box needs no boundary machinery, so the stability chapters come **before** the
wave-generation chapter (11).

## 2. Chapter map (`main.tex`)

| # | chapter | file | state |
|---|---|---|---|
| 1 | σ-Euler model | `BALFEM_derivation/SigmaEulerModel.tex` | stable (derivation project) |
| 2 | Vertical multilayer discretisation | `BALFEM_derivation/VerticalFESemiDiscretisation/*.tex` | stable (derivation project) |
| 3 | Linearised model | `BALFEM_derivation/LinearModel.tex` | stable (derivation project) |
| 4 | Model verification (Yang & Liu collapse) | `BALFEM_derivation/ModelVerification.tex` | stable (derivation project) |
| 5 | Stokes-wave-type Fourier analysis | `BALFEM_derivation/StokesWaveFourierAnalysis.tex` | stable (derivation project) |
| 6 | Horizontal weak formulation and regularity | `HorizontalFormulation/WeakForm.tex` | rev. 2026-11-06 (cellwise Class III + kind-S facet terms, C⁰-CIP §6.4) |
| 7 | Multi-field residual implementation | `NumericalImplementation/GlobalResidual.tex` | the complete model only |
| 8 | Multi-fidelity implementation | `NumericalImplementation/MultiFidelity.tex` | new 2026-11-09 (switches, six models, term classification) |
| 9 | Stability of the unstabilised discretisation | `SolverValidation/StabilityAnalysis.tex` | v2 Q2/Q1 rows in; **Q3/Q2 rows and the 2026-11 box / integrator results to add** (`STABILITY.md` §3b, §6b) |
| 10 | Stabilisation | `SolverValidation/Stabilisation.tex` | describes the ghost penalty; **to be rewritten** (the integrator route is closed, `INTEGRATOR_STUDY_RESULTS.md` §5) |
| 11 | Gridap solver implementation (generation, absorption, BCs, distributed) | `NumericalImplementation/GridapImplementation.tex` | at the end on purpose |
| 12 | Validation | `SolverValidation/ValidationTests.tex` | v1 MMS results |
| A | Projected Class-III treatment (v1 record) | `ClassIIIterms.tex` | appendix |
| B | Mixed formulation (v1 record) | `Appendices/MixedFormulation.tex` | appendix; the v2 mixed solver is `GridapBALFEMixed.jl` |
| C | Jump-penalty stabilisation (first-order C⁰-IP, hp `:jumpgrad`) | `Appendices/JumpPenalty.tex` | appendix; its γ window predates the classifier correction (`STABILITY.md` §5) |

**Chapter 6, as the author designed it:**
* **6.1** The broken framework on `C⁰` Lagrange spaces: tessellation, traces, the **single** identity
  `eq: broken integration by parts identity`, regularity and admissibility, facet kinds.
* **6.2** Stacked multi-field notation and the Gridap global residual.
* **6.3** One term-by-term broken derivation (`sec: broken weak form term by term derivation`). It
  ends with the broken global residual `eq: multifield general residual`. The leading-pressure term
  is the only integration by parts, and it refers to the 6.1 identity.
* **6.4** The `𝓛` operator.
* **6.5** The `𝓝` operator: classification, Classes I–III. The Class-III forms (A)/(B) are
  corollaries of the identity.
* **6.6** Model setups.

The Taylor–Hood section is in chapter 9. LaTeX commit `5381f4f` (in `BALFEM_models_v2`) holds the identity-based rewrite.

## 3. Conventions and known loose ends

* Every label is kept when its content moves; v1 results carry "(v1)" in their captions.
* Undefined reference `subsec: comparison point derived` (`ModelVerification.tex`), pre-existing.
* "`Q1` zeroes the B-dispersion term" was never re-derived for the 2-D mixed derivatives.
* Never global-replace `H^2`/`L^2`: it hits the water depth `H²`, `h²` and `\partial^2` (`CLAUDE.md` §2).
