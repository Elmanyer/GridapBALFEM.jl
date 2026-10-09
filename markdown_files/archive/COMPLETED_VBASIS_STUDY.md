# COMPLETED_VBASIS_STUDY.md — planned studies not yet started

> **Scope — what belongs here.** The *design record of one completed study*: the vertical-basis
> convergence study, its prerequisites and its reasoning, kept so the design can be reused.
>
> ⚠ **This file was called `PENDING_TASKS.md`, which had become wrong** — its only study completed
> on 2026-08-30. Results: [`MMS_VBASIS_CAMPAIGN.md`](MMS_VBASIS_CAMPAIGN.md) and
> `output/local/mms/vbasis_campaign_2026-08-30/`. New work goes in
> [`PLANNED_CAMPAIGNS.md`](PLANNED_CAMPAIGNS.md).


Work that is **designed but not begun**, kept separate from [`OPEN_ISSUES.md`](OPEN_ISSUES.md)
(which tracks gaps in work already done). Each entry states the goal, the prerequisites that block
it, and the decisive design — so it can be picked up without re-deriving the reasoning.

---

## 1. Vertical-basis convergence study — rates vs `(M, p)` ✅ COMPLETE (2026-08-30)

> **RUN AND FINISHED.** 83 studies; results, method, corrections and scope limits in
> [`MMS_VBASIS_CAMPAIGN.md`](MMS_VBASIS_CAMPAIGN.md). Headline: the spatial order of accuracy is
> **independent of the vertical basis** (30/30 at `η`→3, `u`→4 across `Nσ` = 3, 4, 5). Tier 3
> answered (the `:full` floor grows with `Nσ`). One scope limit stands: `Nσ = 5` nonlinear TEMPORAL
> is not measurable on this machine. The section below is kept as the design record.

> **Status change 2026-08-21.** Prerequisite A is **CLEARED** — both `mms_driver.jl` defects are
> fixed and the study driver exists: `examples/local_mms/run_vertical_basis_study.jl` (tiers 1–3,
> `VB_BASES` / `VB_MODELS`). Prerequisite B (optimised σ-meshes for `p ≥ 2`) remains open but is
> **NOT blocking**: node positions change the error CONSTANT, never the ORDER, so tier 1 can run on
> `resolve_cbdy`'s fallback today. Two data points already exist from the fix verification, both
> Q2/Q1 1-D static Model 1: **P1LFE-3 → `p_η` 1.994, `p_u` 2.997** and **P2LFE-2 → `p_η` 1.994,
> `p_u` 2.996**, i.e. theoretical order at `M ≠ 2` and at `p = 2`. That is encouraging, not the
> study: it is one model, one pairing, two levels.

### The goal

The convergence campaign to date varies **only the horizontal** FE order: `Q2/Q1`, `Q3/Q2`, `Q4/Q3`
(see [`VERIFIED_SCOPE.md`](VERIFIED_SCOPE.md) §3). Every study, every model, every gate in the entire
repository runs a single vertical configuration — **P1LFE-2 (`M=2`, `p=1`, `Nσ=3`)**.

That is a gap in exactly the claim the project is named for. **BALFE-M asserts basis-agnosticism —
that the model family works for an *arbitrary* vertical FE basis** — and the evidence for it is
currently one member. The study measures convergence rates across `(M, p)` and answers whether the
discretisation reaches theoretical order independently of the vertical basis.

### ⚠ This study belongs LOCAL, not on the cluster

Decided 2026-08-19, and worth recording because the opposite is tempting.

**MMS measures the order of accuracy under *mesh refinement*. The domain size does not enter the
rate at all** — `h = Lx/nx`, so a bigger domain at fixed `nx` is simply a coarser mesh. Cluster
resources buy nothing a laptop cannot deliver for a rate study; they only buy *more refinement
levels* and *2-D* variants, both of which are cost, not new information about the operator.

Two further ceilings make "throw cluster at it" actively wrong:

* **1-D `Q4/Q3` already sits on the round-off floor** (`e_u`: 1.79e-11 → 1.14e-11 over the last
  refinement). Double precision, not compute, is the binding constraint there.
* The forcing cost scales as **`Nσ²`** (the `Ψ` block carries an `Nσ²×8` component sum that the
  outer gradient then differentiates), so the expensive axis is the *vertical* resolution, which is
  a serial cost per quadrature point — it does not parallelise away.

**Keep every MMS convergence study sequential and local.**

### Cost of the vertical axis

| configuration | `Nσ = M·p+1` | forcing cost vs P1LFE-2 |
|---|---|---|
| P1LFE-2 | 3 | 1.0× |
| P1LFE-3 | 4 | 1.8× |
| P1LFE-4 | 5 | 2.8× |
| P2LFE-2 | 5 | 2.8× |
| P2LFE-3 | 7 | 5.4× |
| P2LFE-4 | 9 | 9.0× |

Multiply by the tier factor already measured: `:native` forcing costs ~4–6× `:none`. A P2LFE-4
`:native` study is therefore ~50× a P1LFE-2 `:none` one. **Prune the grid before running it.**

### ✅ Prerequisite A — the two `src/mms_driver.jl` defects (FIXED 2026-08-21)

Found 2026-08-19 during the design of this study; both were interface changes, not call-site fixes.
Recorded in full because the *shape* of each is worth recognising again.

**A1 — `run_conv_study` could not vary `M` or `p`.** It hard-coded

```julia
vert = assemble_vertical_tensors(M, 1, [0.0, 0.728, 1.0])
```

so `p_vert` was not a parameter at all and `c_bdy` was pinned to the **M=2** node set, against
`assemble_vertical_tensors`'s `length(c_bdy) == M+1` assertion (`src/vertical.jl:90`) — the
signature accepted `M`, advertised it as a degree of freedom, and ignored it in the vertical build,
so `run_conv_study(M=3)` threw immediately.

*Fixed:* `p_vert` and `c_bdy` are kwargs; `c_bdy === nothing` resolves through the new shared
`resolve_cbdy(M, c_bdy)` (`src/utilities.jl`) that `setup_and_run` and `setup_and_run_distributed`
now also use — one definition, so a default valid for one `M` cannot be hard-wired anywhere again.
The study tag carries `P{p}LFE-{M}` and the return value carries `M, p_vert, c_bdy, Nsigma`.

**A2 — every DISTRIBUTED MMS study silently ran Model 1.** `run_mms_case_distributed` contained

```julia
prob = build_problem(vert; g=g, h_bathy=(x -> d), regime=:linear,
                     nl_pressure=:none, flat_bed=true,
                     mms_src=mms_forcing_stage1(f, vert, d, g))
```

as **literals**, accepted no `regime`/`flat_bed`/`nl_pressure`/`hfun`, and used the Model-1 closed
form rather than the general `mms_forcing` — while `run_conv_study` still built its printed tag from
the *requested* switches.

> **Consequence: a campaign over 8 models returned 8 identical Model-1 studies under 8 different
> labels, and all of them PASSED** — Model 1 converges optimally. It failed silently, confidently,
> and inside the one tool whose purpose is to catch self-consistently wrong results.

*Fixed:* the switches and `hfun` are parameters, the forcing comes from `mms_forcing`, and both
branches of `run_conv_study` are called from ONE shared keyword bundle so they cannot diverge again.
**Gate:** `test/test_mms_distributed_parity.jl` (4 ranks, in `runtests.jl`'s `MPI_TESTS`) compares
the two branches on Model 2 (moves `flat_bed`) and Model 5 (moves `regime` and `nl_pressure`) —
**preceded by a SEPARATION negative control**, because a bare parity check would have passed with
the defect present on any configuration whose answer coincides with Model 1's. Its reference is
recomputed in-process, never pinned, so it cannot go stale.

> **The lesson worth carrying: a label is not evidence.** The tag was built from the requested
> switches and read correctly through an entire campaign that computed one model eight times.
> Whenever a runner prints a configuration it did not itself hand to the thing under test, that gap
> is where a silent-wrong-answer lives.

### 🟠 Prerequisite B — optimised σ-meshes for each `(M, p)`  (open, but NOT blocking tier 1)

**Status: the author is reviewing the vertical mesh optimisation section of
`BALFEM_models/StokesWaveFourierAnalysis.tex` before this is specified.**

`DEFAULT_CBDY` (`src/utilities.jl`) tabulates ELEMENT BOUNDARIES for `M = 1…4` **at `p = 1` only** —
these are Yang & Liu (2024) Table 1, derived for piecewise-linear elements. For any `p ≥ 2` there is
no published optimum. `resolve_cbdy` returns that table's entry for the requested `M` (valid as a
mesh at any `p`, since for `p ≥ 2` the extra nodes are interior to each element) and falls back to a
**uniform** `LinRange(0,1,M+1)` for `M > 4`.

The LaTeX derivation supplies the machinery to compute them: the dispersion functional
`R(μ) = Φᵀ(M + μ|B|)⁻¹Φ`, its four basis-independent properties, and the grid-optimisation result
`Δσ_top ≈ 2.94/kd_max`. Generalising that to arbitrary `p` and tabulating the resulting node sets is
a piece of work in its own right, and it feeds two things: this study, and the applicable-`kd`
table that characterises every new family member.

> **⚠ Do not conflate the two questions.** The optimised node positions change the error
> **constant**, never the **order**. So:
> * the convergence study can proceed on *any* reasonable node set (uniform is fine) — it tests the
>   generality of the **discretisation**;
> * the node optimisation is a **dispersion-accuracy** question, measured by applicable `kd`
>   (`test_dispersion_curve`'s territory), not by a rate.
>
> Reporting a rate study as evidence that the optimised meshes are right — or vice versa — would be
> a category error. They are independent axes and should be presented as such.

### Design — runnable now (A cleared; B affects the constant, not the order)

Implemented as `examples/local_mms/run_vertical_basis_study.jl`. `VB_BASES="M:p,…"` and
`VB_MODELS="1,…,8"` select the grid; the three tiers below are the intended pruning.

Three tiers, so the expensive parts are only paid for if the cheap ones justify them.

| tier | grid | studies | purpose |
|---|---|---|---|
| **1 — the paper claim** | 6 verifiable models × `(M,p) ∈ {(2,1),(3,1),(4,1),(2,2),(3,2)}`, `Q3/Q2`, 1-D static, 4 levels | 30 | does theoretical order hold independently of the vertical basis? |
| **2 — the open question** | 3 pairings × 2 models, `(M,p)=(2,1)`, **2-D**, 5 levels | 6 | is the `Q2/Q1` and `Q4/Q3` velocity shortfall genuinely suboptimal or merely pre-asymptotic? |
| **3 — the `:full` floor** | 2 models × 2 `(M,p)` | 4 | ✅ **ANSWERED 2026-08-29: YES, it grows** — 1.2 → 2.0 → 2.9 e-03 for `Nσ` = 3 → 4 → 5 (`MMS_VBASIS_CAMPAIGN.md` §2.6). ⚠ It is an OMISSION floor, not a frozen-projection one — see §2.1 |

**Constraints that must be carried into any run script:**

* **Tier 3 is NOT a rate study.** `:full` pins `p_u` at −0.00 by construction
  ([`VERIFIED_SCOPE.md`](VERIFIED_SCOPE.md) §4). Report the **floor value**, labelled as such. Enabling
  it as a rate gate would manufacture a defect that does not exist.
* **Tolerances must be ~8 orders tighter than production**: `ls_rtol ≈ 1e-13`, `nl_tol ≈ 1e-12…1e-14`.
  The driver defaults are already correct; the production defaults (`1e-5`) would make the study
  measure the *solver*, not the discretisation. Note `nl_tol=1e-14` is unreachable for
  `mode=:transient` — use `1e-12` there, or keep the study static.
* **Re-run the tolerance-independence gate per tier**, not once for the campaign. Without it every
  rate in the tier is void.
* **Read the pairwise rate sequence, not the fitted slope**, and check error *magnitude* before
  trusting any fine-level high-order rate ([`VERIFIED_SCOPE.md`](VERIFIED_SCOPE.md) §6).

### What it would close

* "Only one vertical resolution tested" — the standing scope limit in
  [`TEST_SUITE.md`](TEST_SUITE.md) §8. (The *capability* gap is closed: every MMS test and example
  now takes `MMS_M`/`MMS_PVERT` or `VB_BASES`. The *measurement* gap closes only when the sweep runs.)
* The `Q2/Q1` / `Q4/Q3` velocity shortfall ([`OPEN_ISSUES.md`](OPEN_ISSUES.md) §6).
* Direct quantitative support for the basis-agnosticism claim in the LaTeX and the CFC 2027 abstract.
