# `markdown_files/` — the project's written record

Thirteen documents plus [`CLAUDE.md`](CLAUDE.md), which is the map and the standing rules. This
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
| [`PLANNED_CAMPAIGNS.md`](PLANNED_CAMPAIGNS.md) | **What runs next, and what must be true first?** | "18 studies, 84 solver runs, and the driver changes they need" |

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
| [`CLAUDE.md`](CLAUDE.md) | **start here** — repository map, physics switches, and the standing rules |
| [`MODEL.md`](MODEL.md) | the mathematics: σ-tensors, the global residual term by term, Jacobians |
| [`ARCHITECTURE.md`](ARCHITECTURE.md) | code structure: the stacked layout, `src/` map, FE spaces, time loops, distributed path |
| [`CONFIGURATION.md`](CONFIGURATION.md) | the settings, the evidence behind each, measured performance |
| [`WAVE_GENERATION.md`](WAVE_GENERATION.md) | sources, Dirichlet generation, sponge, WaveSpec coupling |
| [`RUNNING.md`](RUNNING.md) | how to launch: local, cluster, sysimage, env vars |
| [`MMS_VBASIS_CAMPAIGN.md`](MMS_VBASIS_CAMPAIGN.md) | results of the `(M,p)` × 8-model convergence matrix |
| [`CAMPAIGN_COST.md`](CAMPAIGN_COST.md) | measured run costs, memory bands, queue limits, scheduling lessons |
| [`OUTPUT_NAMING_PROPOSAL.md`](OUTPUT_NAMING_PROPOSAL.md) | the output-directory grammar — accepted and implemented |

## Recurring themes, and where each is authoritative

Several topics appear in many files. To avoid contradictory copies, each has **one** home:

| topic | authoritative in | others should link, not restate |
|---|---|---|
| the Taylor-Hood requirement and the equal-order instability | `CLAUDE.md` rules 2b / 12b | |
| what the MMS proves, and its boundaries | `VERIFIED_SCOPE.md` | |
| the `Q2/Q1` one-order velocity shortfall | `PLANNED_CAMPAIGNS.md` §2 (it is what the next campaign tests) | |
| why `:full` cannot reach optimal order | `PLANNED_CAMPAIGNS.md` §3 | |
| the vopt correction and the deleted κ | `CLAUDE.md` rule 46 | |
| the algebraic error floor (~1e-10) | `PLANNED_CAMPAIGNS.md` §0 | |
| run cost and memory bands | `CAMPAIGN_COST.md` | |
