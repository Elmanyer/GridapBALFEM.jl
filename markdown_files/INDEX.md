# `markdown_files/` — the project's written record

Sixteen documents and this index, plus [`CLAUDE.md`](../CLAUDE.md) at the repository root, which is the map and the standing rules. This
index exists so the right file is obvious without opening four of them.

⚠ **This directory is tracked; `latex_docs/` is gitignored.** Anything load-bearing therefore
belongs here, not in the LaTeX.

## The four that are easy to confuse

They are separated by **the question each answers**, and each carries a scope header saying what it
must not contain:

| file | the question | example of what belongs |
|---|---|---|
| [`VERIFIED_SCOPE.md`](VERIFIED_SCOPE.md) | **Is it correct, and how far does the claim reach?** | "six of eight models verified at theoretical order on five vertical bases" |
| [`TEST_SUITE.md`](TEST_SUITE.md) | **What is checked, and what would slip through?** | "this gate is a bounds check, so it passes while the quantity moves 58 %" |
| [`OPEN_ISSUES.md`](OPEN_ISSUES.md) | **What is known wrong, missing, or unexplained?** | "the `:sdirk` linear-model temporal deficit is unexplained" |
| [`PLANNED_CAMPAIGNS.md`](PLANNED_CAMPAIGNS.md) | **What runs next, and what must be true first?** | "the 2-D pairing tier, its short ladder, and the code it needs" |

**The boundary that was being crossed:** a run nobody has launched is a **plan**, not a gap. It
belongs in `PLANNED_CAMPAIGNS.md`. On 2026-09-11 the "Runs not yet made" list moved out of
`OPEN_ISSUES.md` for exactly this reason.

### Renames, 2026-09-11

| was | now | why |
|---|---|---|
| `NEXT_TESTS.md` | `PLANNED_CAMPAIGNS.md` | it plans campaigns, not only tests |
| `VERIFICATION.md` | `VERIFIED_SCOPE.md` | the content is the *scope* of what is proven, not the act |
| `OPEN_ITEMS.md` | `OPEN_ISSUES.md` | "items" read as a to-do list; these are defects |
| `PENDING_TASKS.md` | `COMPLETED_VBASIS_STUDY.md` | ⚠ **the name had become false** — its only study finished 2026-08-30 |

## Everything else

| file | what it holds |
|---|---|
| [`CLAUDE.md`](../CLAUDE.md) | **start here** — repository map, physics switches, and the standing rules |
| [`MODEL.md`](MODEL.md) | the mathematics: σ-tensors, the global residual term by term, Jacobians |
| [`ARCHITECTURE.md`](ARCHITECTURE.md) | code structure: the stacked layout, `src/` map, FE spaces, time loops, distributed path |
| [`CONFIGURATION.md`](CONFIGURATION.md) | the settings, the evidence behind each, measured performance |
| [`WAVE_GENERATION.md`](WAVE_GENERATION.md) | sources, Dirichlet generation, sponge, WaveSpec coupling |
| [`RUNNING.md`](RUNNING.md) | how to launch: local, cluster, sysimage, env vars |
| [`HORIZONTAL_CONVERGENCE.md`](HORIZONTAL_CONVERGENCE.md) | the horizontal pairing × physics convergence matrix, 1-D (18 studies, complete) and 2-D (partial), and why it does not localise a failure |
| [`MMS_VBASIS_CAMPAIGN.md`](MMS_VBASIS_CAMPAIGN.md) | results of the `(M,p)` × 8-model convergence matrix |
| [`CAMPAIGN_COST.md`](CAMPAIGN_COST.md) | measured run costs, memory bands, queue limits, scheduling lessons |
| [`OUTPUT_NAMING_PROPOSAL.md`](OUTPUT_NAMING_PROPOSAL.md) | the output-directory grammar — accepted and implemented |
| [`NEW_TREATMENT.md`](NEW_TREATMENT.md) | the Class-III treatments and every campaign run on them: reduction, in-loop projections, the **mixed formulation**, the flume factorial, the closed periodic box (Parts A–I) |
| [`V2_SOLVER_PLAN.md`](V2_SOLVER_PLAN.md) | **v2 solver plan**: the main changes (broken formulation only, `nl_pressure::Bool`, all eight `𝓝` components) and the step-by-step implementation |
| [`LATEX_STRUCTURE.md`](LATEX_STRUCTURE.md) | **v2 document structure** (`latex_docs/BALFEM_models_v2/`): the chapter order and why — weak form and regularity, multi-field implementation, stability, stabilisation, solver functionalities, validation |
| [`COMPLETED_VBASIS_STUDY.md`](COMPLETED_VBASIS_STUDY.md) | design record of the finished vertical-basis study |

## Recurring themes, and where each is authoritative

Several topics appear in many files. To avoid contradictory copies, each has **one** home:

| topic | authoritative in | others should link, not restate |
|---|---|---|
| the Taylor-Hood requirement and the equal-order instability | `CLAUDE.md` rules 2b / 12b | |
| what the MMS proves, and its boundaries | `VERIFIED_SCOPE.md` | |
| the `Q2/Q1` one-order velocity shortfall | `PLANNED_CAMPAIGNS.md` §2 — **measured on all six models 2026-09-12** | |
| the nonlinear `p_η` order reduction | `OPEN_ISSUES.md` §0b | |
| ⛔ **the `:full` discretisation is unstable (interior, grid-scale) — the open stabilisation item**: evidence, mechanism, candidates, acceptance criteria | `OPEN_ISSUES.md` §0e (evidence `CLAUDE.md` §5.2f; analysis LaTeX chapter 8) | |
| the stabilisation campaign (filter, CIP, the box ladder under Crank–Nicolson, the gate) | `PLANNED_CAMPAIGNS.md` §6c | |
| ⚠ Gridap's `SDIRK_2_2` is `DIRK22(1,0,1)`: A-stable, **not** L-stable, ≈ `0.75(ωΔt)⁴` per step — stability claims need Crank–Nicolson | `CLAUDE.md` rule 15 | |
| the closed x-periodic box as the discriminator between interior and boundary instabilities | `CLAUDE.md` rule 38i, §5.2f | |
| the `:full` Class-III instability — how the operator, the Jacobian and the assembly were ruled out, what "Class III" means, why projection was chosen | `OPEN_ISSUES.md` §0c (history; superseded as the active item by §0e) | |
| why the Class-III terms cannot simply be integrated by parts, and the trial-vs-test correction | `CLAUDE.md` rule 1b, `OPEN_ISSUES.md` §0c §5/§7 | |
| the earlier ordered fixes for `:full` (Q3/Q2 → de-lag → algebraic reduction → `C⁰`-IP) — spent or superseded | `PLANNED_CAMPAIGNS.md` §6b items 0a–0c | |
| ⚙ **the Class-III treatments** — the exact `{1,2,5}` reduction, the in-loop projections (broken, off by default), the **mixed formulation** (implemented, `src/mixed.jl`) and their campaigns | `NEW_TREATMENT.md` (branches `new-classIII-treatment` → `mixed-formulation-solver`) | |
| ⛔ the incomplete Jacobian does NOT affect stability (proven) | `CLAUDE.md` rule 17b | |
| the `dx`×`dt`×Jacobian refinement factorial | `CLAUDE.md` §5.2c | |
| the collapse of `p=1` BALFE-M onto Yang & Liu's LFE-M | `CLAUDE.md` §5.2d, `VERIFIED_SCOPE.md` §0b | |
| the variable-bed lee-shoulder mode, and `\|∇h\|` as its rate | `OPEN_ISSUES.md` §0d | |
| the 100-period flat-bed stability result, and its caveats | `CLAUDE.md` §5.2b | |
| stability follow-ups still owed | `PLANNED_CAMPAIGNS.md` §6c (the §6b list is spent) | |
| the pairing × physics convergence matrix (1-D and 2-D), and the "optimal in one field or the other" pattern | `HORIZONTAL_CONVERGENCE.md` | |
| measured study cost and memory bands per pairing | `CAMPAIGN_COST.md` §2 | |
| why `:full` cannot reach optimal order | `PLANNED_CAMPAIGNS.md` §3 | |
| the vopt correction and the deleted κ | `CLAUDE.md` rule 46 | |
| the algebraic error floor (~1e-10) | `PLANNED_CAMPAIGNS.md` §0 | |
| run cost and memory bands | `CAMPAIGN_COST.md` | |
