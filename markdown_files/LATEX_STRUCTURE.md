# LATEX_STRUCTURE.md — the v2 LaTeX document (`latex_docs/BALFEM_models_v2/`): structure and state

*Condensed 2026-10-09. The full reasoning behind the order (revision 4, with the migration map and
history) is [`archive/LATEX_STRUCTURE_rev4.md`](archive/LATEX_STRUCTURE_rev4.md). The structure below
is **applied**. The document is its own git repository (remote `overleaf`) and the author edits it in
Overleaf: always `git fetch overleaf` and fast-forward before editing (see `CLAUDE.md` §2).*

## 1. The principle

The document follows the development of the solver: each step builds something, analyses it, and
the next step fixes what the analysis finds. A chapter uses only what precedes it. Chapters 6–9 form
one chain of four questions, and **keep consistency apart from stability** (the confusion that cost a
year, rule 39b):

```
 6  Is the weak form well defined on C⁰ elements?  → only on the broken domain; the broken formulation
 7  How is it written in Gridap and advanced in time? → stacked layout, Newton, unit vertical basis
 8  Is the discretisation stable?                     → NO: consistent but unstable at the grid scale
 9  What sink, and does it work?                      → (being rewritten: filtering, STATUS.md §4)
```

The closed periodic box needs no boundary machinery, so the stability chapters come **before** the
wave-generation chapter (10).

## 2. Chapter map (`main.tex`)

| # | chapter | file | state |
|---|---|---|---|
| 1 | σ-Euler model | `SigmaEulerModel.tex` | stable |
| 2 | Vertical multilayer discretisation | `VerticalFESemiDiscretisation/*.tex` | stable |
| 3 | Linearised model | `LinearModel.tex` | stable |
| 4 | Model verification (Yang & Liu collapse) | `ModelVerification.tex` | stable |
| 5 | Stokes-wave-type Fourier analysis | `StokesWaveFourierAnalysis.tex` | stable |
| 6 | Horizontal weak formulation and regularity | `HorizontalFormulation/WeakForm.tex` | **restructured 2026-10-06/08** (below) |
| 7 | Multi-field residual implementation | `NumericalImplementation/GlobalResidual.tex` | restructured |
| 8 | Stability of the unstabilised discretisation | `SolverValidation/StabilityAnalysis.tex` | v2 Q2/Q1 rows in; **Q3/Q2 rows to add** (`STABILITY.md` §3) |
| 9 | Stabilisation | `SolverValidation/Stabilisation.tex` | describes the ghost penalty as adopted; **to be rewritten around filtering** |
| 10 | Gridap solver implementation (generation, absorption, BCs, distributed) | `NumericalImplementation/GridapImplementation.tex` | moved to the end |
| 11 | Validation | `SolverValidation/ValidationTests.tex` | v1 MMS results |
| A | Projected Class-III treatment (v1 record) | `ClassIIIterms.tex` | appendix |
| B | Mixed formulation (v1 record) | `Appendices/MixedFormulation.tex` | appendix |
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

The Taylor–Hood section is in chapter 8. LaTeX commit `5381f4f` holds the identity-based rewrite.

## 3. Conventions and known loose ends

* Every label is kept when its content moves; v1 results carry "(v1)" in their captions.
* Undefined reference `subsec: comparison point derived` (`ModelVerification.tex`), pre-existing.
* "`Q1` zeroes the B-dispersion term" was never re-derived for the 2-D mixed derivatives.
* Never global-replace `H^2`/`L^2`: it hits the water depth `H²`, `h²` and `\partial^2` (`CLAUDE.md` §2).
