# `markdown_files/` — the project's written record

*Restructured 2026-10-09: completed plans and full campaign records moved to [`archive/`](archive/README.md);
the current state lives in a few short files.* ⚠ This directory is **tracked**; `latex_docs/` is
gitignored, so anything load-bearing belongs here.

## Reading order for a new agent

1. [`../CLAUDE.md`](../CLAUDE.md): the repository map, physics switches and standing rules.
2. [`STATUS.md`](STATUS.md): **where the solver stands, the open problem, the next step.**
3. [`STABILITY.md`](STABILITY.md): the instability and every stabilisation attempt, condensed.
4. [`OPEN_ISSUES.md`](OPEN_ISSUES.md): every other open item, with its next step.

## Reference (current)

| file | the question it answers |
|---|---|
| [`MODEL.md`](MODEL.md) | the maths: σ-tensors, the residual term by term in Gridap, `𝓝`, Jacobians |
| [`ARCHITECTURE.md`](ARCHITECTURE.md) | the code: stacked layout, `src/` map, FE spaces, time loops, distributed path, monitoring |
| [`CONFIGURATION.md`](CONFIGURATION.md) | settings and the evidence behind each; the Gridap fork; tolerances; measured cost |
| [`RUNNING.md`](RUNNING.md) | how to launch: local, cluster, sysimage, environment variables |
| [`WAVE_GENERATION.md`](WAVE_GENERATION.md) | sources, Dirichlet generation, sponge, relaxation, WaveSpec |
| [`VERIFIED_SCOPE.md`](VERIFIED_SCOPE.md) | is it correct, and how far does the claim reach? (MMS, AD oracle, Yang & Liu collapse; v1 measurements, models 1–4 unchanged in v2) |
| [`TEST_SUITE.md`](TEST_SUITE.md) | what each test checks, and what would slip through (v1 scores; v2 status in `STATUS.md`) |
| [`V2_SOLVER_PLAN.md`](V2_SOLVER_PLAN.md) | what v2 changed from v1 and why; the execution record; open step 10 (distributed) |
| [`LATEX_STRUCTURE.md`](LATEX_STRUCTURE.md) | the v2 LaTeX document: chapter map, state of each chapter |
| [`OUTPUT_NAMING_PROPOSAL.md`](OUTPUT_NAMING_PROPOSAL.md) | the output-directory grammar (implemented; rule 2c) |

## Where each recurring topic lives

| topic | home |
|---|---|
| current state, the filter decision and its options | `STATUS.md` |
| the grid-scale instability: evidence, mechanism, penalties tried, acceptance criteria | `STABILITY.md` |
| Taylor–Hood requirement / the equal-order artefact | `CLAUDE.md` rules 2b, 12b |
| `SDIRK_2_2` is dissipative; stability needs Crank–Nicolson | `CLAUDE.md` rule 15 |
| the incomplete Jacobian does not affect stability | `CLAUDE.md` rule 17b |
| `C⁰` trial vs test regularity; why the broken layer | `CLAUDE.md` rule 1b, LaTeX chapter 6 |
| nonlinear `p_η` order reduction; lee-shoulder mode; test failures | `OPEN_ISSUES.md` §1–3 |
| the `p = 1` collapse onto Yang & Liu | `VERIFIED_SCOPE.md` §0b |
| vertical-grid optimisation (vopt) | `CLAUDE.md` §3 and rule 46 |
| any v1 measurement in full | `archive/` (start at `archive/README.md`) |
