# OPEN_ISSUES.md — open work, each with its decisive next step

> **Scope — what belongs here.** *Known defects and gaps in work already done*, each with its
> decisive next step. Answers **"what is known to be wrong, missing, or unexplained?"**
>
> **What does NOT belong here:** work that is simply not started → [`PLANNED_CAMPAIGNS.md`](PLANNED_CAMPAIGNS.md)
> (a run nobody has launched is a plan, not a gap); what is proven → [`VERIFIED_SCOPE.md`](VERIFIED_SCOPE.md);
> test coverage → [`TEST_SUITE.md`](TEST_SUITE.md).


Ordered by whether the next step is *decided* or still needs a judgement call. Nothing in the solver
is half-built; every item below is a gap in verification, performance, or follow-through.

Related: [`VERIFIED_SCOPE.md`](VERIFIED_SCOPE.md), [`TEST_SUITE.md`](TEST_SUITE.md),
[`CONFIGURATION.md`](CONFIGURATION.md).

---

## 0. ✅ Grid-scale instability in the fully nonlinear model — RESOLVED 2026-09-06

**It was the equal-order `Q2/Q2` horizontal pairing, not the model.** `η` enters momentum
undifferentiated via `∇·v`, so it plays the pressure role of a Stokes system and equal-order
continuous spaces are inf-sup deficient, admitting a spurious checkerboard at `λ ≈ 2·dx` — exactly
where the measured growth spectrum peaked. On Taylor-Hood the mode does not exist and the refinement
signature **inverts** (settled `η` 0.11055 → 0.10593 → 0.10481 as `dx` halves, converging on the
delivered 0.102, where refinement used to move the blow-up from t≈24 s to t≈10 s).

The fix cost **no residual change and no re-verification** — `Q2/Q1` and `Q3/Q2` are the pairings the
MMS campaign already certifies. It is now enforced by `check_taylor_hood` (`CLAUDE.md` rule 2b), so
the configuration cannot be selected again. The complete account — ten refuted hypotheses, the
campaign results across three meshes and three integrators, and the method lessons — is **`CLAUDE.md`
rule 12b**, which absorbed the former `NONLINEAR_INSTABILITY.md`.

⚠ **The quantitative phenomenology is void, not superseded**: the `σ ∝ A⁴` rates, the "no amplitude
threshold" claim, the tier ordering and the `dt`-masking numbers all describe a discretisation no
longer in use.

**What this leaves open, and it is real work:**
* **The 17 physics/smoke tests need re-baselining onto Taylor-Hood.** They call `setup_and_run`
  without an explicit `p_eta`, so they moved from `Q2/Q2` to `Q2/Q1` when the default changed. Their
  pinned constants were measured on equal order and must be **re-measured, not re-thresholded**.
  ✅ The verification tier is unaffected — it always passed `p_eta` explicitly.
* **A long-duration nonlinear regression is still missing.** No test runs long enough to have caught
  this (the mode needed 50–80 s to emerge); that gap is why it reached production. The Taylor-Hood
  80 s flume run is the natural basis for one.
* **Production runs made on equal order should be re-examined**, since they were made on a
  discretisation now known to be deficient.

---

## 0b. 🔴 THE NONLINEAR MODELS LOSE AN ORDER IN `p_η` — narrowed to the advection block

> Found 2026-09-11 by Campaign C, cause narrowed 2026-09-12. Letter-suffixed so the ten
> sections below keep their numbers (the §4 renumbering hazard, in miniature).

**THE OBSERVATION.** Campaign C, 1-D, Q3/Q2, `a_eta = 0.8`, five-level ladder `nx = 4…64`,
`nl_tol = 1e-14`. Same basis, same ladder, same mesh; only `regime` differs:

| | pairwise `p_η` | `e_η` ratios (8.00 = optimal 3rd order) |
|---|---|---|
| **linear** (models 1, 2) | 2.990, 2.997, 2.999, **3.000** | 7.95  7.99  8.00  8.00 |
| **nonlinear** (models 3, 4, 5, 6) | 2.971, 2.922, 2.757, **2.450** | 7.84  7.58  6.76  5.46 |

⚠ **This is an ORDER REDUCTION, not a floor.** `e_η` keeps falling — it is the *rate* that decays,
from 3 toward 2, and it decays MONOTONICALLY WITH REFINEMENT. That is the signature of a fixed
lower-order error component overtaking the third-order one (rule 33), and it means the problem gets
worse on better meshes, not better.

⚠ **`p_u` moves the OTHER WAY on the same runs** — 2.77 → 3.19 → 3.58 → 3.69, rising toward its
optimal 4 while `p_η` falls. Any explanation has to account for both signs at once; several
otherwise-plausible ones do not.

**WHERE IT DOES AND DOES NOT SHOW, and why that is consistent.**

| pairing | optimal `p_η` | measured | reading |
|---|---|---|---|
| Q2/Q1 | 2 | 2.000 exactly | a lower-order component at order 2 CANNOT be seen — it is degenerate with the optimal rate |
| Q3/Q2 | 3 | 2.450 and falling | visible, because the third-order part decays past it |
| Q4/Q3 | 4 | 4.000–4.002 | ladder stops at `nx=32`; discretisation error still dominates |

So the absence at Q2/Q1 and Q4/Q3 is expected under the same explanation, not evidence against it.

**WHAT IT IS NOT — five eliminations, each with its own measurement.**

| ruled out | evidence |
|---|---|
| algebraic contamination (rule 32) | Newton reaches **~1e-15 at every level of every model** (`nl_tol=1e-14`); from the per-level `diagnostics.csv` the §4b work added |
| the `𝓝` pressure blocks | `:none` vs `:native` agree to **4 digits** — m3 2.4498 vs m5 2.4499, m4 2.4494 vs m6 2.4496 |
| bed-slope (`∇h`) code | flat vs variable bed agree to 4 digits; reproduced on an independent process (m4 control, 2.4494) |
| the linear core | linear models are textbook 3.000 on the identical basis, ladder and mesh |
| **quadrature under-integration** | refuted quantitatively — see below |
| the test harness generally | the `quad_extra` control reproduced Campaign C to **5 significant figures**, so the driver edit perturbed nothing |

**THE QUADRATURE HYPOTHESIS, AND WHY IT FAILED — worth keeping, because the shape of the refutation
is the lesson.** The default MMS quadrature is `2·max(p_u,p_η)+2` = degree 8 at Q3/Q2, which in
Gridap is a 5-point Gauss rule **exact to degree 9** (measured; do not read exactness off the
`Measure` argument). Counting from `problem.jl`, per direction: the linear core is ≤ 8 and exact,
while the nonlinear advection `H·(𝓜₃:(U⊗∇U))·W` (line 149) and the nonlinear pressure `H²·(sP·DW)`
(line 110) are **degree 10** — genuinely under-integrated. `src/mms_driver.jl` hard-coded the degree,
so no MMS study had ever been able to vary it; `quad_extra` is now threaded through
`run_conv_study → run_mms_case → Measure`.

Assembling the residual directly at degrees 8/10/12/16 (`output/local/mms/quadrature_test/
knob_liveness.log`):

```
LINEAR    core  ‖R‖ changes 4e-14 across ALL degrees        -> exact at the default
NONLINEAR core  ‖R‖ changes 3.49e-08 from degree 8 -> 10,
                then 4e-14 for 10 -> 12 -> 16               -> under-integrated, exact from 10
```

So the knob is **live** and the crime is **real**. It is also **six orders of magnitude too small**:

| | |
|---|---|
| quadrature error in the residual | 3.5e-08 relative |
| its effect on converged `e_η` | ≤ 4.7e-07 relative ≈ **4.6e-13 absolute** at `nx=64` |
| change needed to restore 3rd order at `nx=64` | `e_η` 9.77e-07 → 6.67e-07 = **3.1e-07 absolute** |

And directly: model 4 at `quad_extra=4`, integrand integrated exactly, returns `p_η = 2.4494` and
ratios `7.840 / 7.577 / 6.756 / 5.462` — **identical to its own control**. Model 3 at `quad_extra=2`
likewise.

⚠ **Two method traps fired during this probe and both are worth remembering.** (a) The first
liveness check integrated `x⁹`, which a degree-8 rule already does exactly — a check that *could not
fail*, and it passed a knob whose reach was still unproven. (b) The `e_η` columns came back BITWISE
identical, which reads as "dead knob, result void" — I called it that — when it was really a live
knob whose effect fell below the CSV's 7-significant-digit print precision. **Neither a liveness
check nor a null result means anything until the effect size is measured against the effect being
explained** (rule 12b's standing lesson, and the same shape as the skew-advection refutation).

**WHAT REMAINS, AND THE DECISIVE NEXT STEP.** By elimination the error sits in the **nonlinear
advection block** (`problem.jl` §"nonlinear advection (𝓜/𝓖 block)", lines 139–149) — a second-order
consistency error in `𝓜₃/𝓖` as discretised. The next step is **derivation, not another run**: check
that block term-by-term against `BALFEM_models/` (`GlobalResidual`, and the `𝓖` identity recorded in
`CLAUDE.md` under the reverted skew-advection work), since the LaTeX is the single source of
mathematical truth and a residual term it does not contain is one nobody can check. Two cheap
discriminators exist if a run is wanted later:

* **amplitude scaling** — a consistency error in the advection block scales with the manufactured
  amplitude. Re-run model 3 at `a_eta = 0.4` and `0.2`: if the crossover moves to finer mesh, the
  offending term is quadratic in the solution; if the rate degrades identically, it is not.
* **extend Q4/Q3 to `nx = 64`** — it should then show the same decay once its fourth-order part has
  decayed past the second-order component. ⚠ Check `e_u` against the ~1e-10 algebraic floor first
  (`PLANNED_CAMPAIGNS.md` §0); at Q4/Q3 that level may be unmeasurable.

**Scope of the damage.** This affects *measured convergence rates of the nonlinear models*, not the
verified scope as stated: `VERIFIED_SCOPE.md` records 30/30 spatial studies at optimal order, and
those were run at Q3/Q2 with **four levels to `nx=32`**, where `p_η` still reads 2.76–3.00. The
degradation is a **fine-mesh** effect that the shorter ladder never reached. The claim "the `:none`
and `:native` models are verified" therefore stands *as far as it was measured*, but it should now
carry the qualifier that nonlinear `p_η` degrades beyond `nx = 32` at Q3/Q2.

---

## 1. 🔴 Cluster memory attribution

`--mem-per-cpu=4G` is in every launcher and is **required**: a `rome` job at the node-default
2 GB/core was OOM-killed with the same error as the earlier crashes. So the compile-spike
explanation is not sufficient on its own. ***Why* is the open question, not *whether*.**

Four hypotheses, one decisive measurement each:

| # | hypothesis | decisive measurement | verdict if true |
|---|---|---|---|
| H1 | per-rank JIT (image absent / stale / incomplete) | **when** the RSS peak occurs — a compile peak is **before step 1**; confirm `-J` was on the command line and the freshness check did not warn | rebuild the image, re-test 2 GB/core |
| H2 | GMRES cache over-allocation | was the job's commit after the `krylov_m` fix? predicted ≈5.6 MB/rank at `krylov_m=100` | if pre-fix, re-run post-fix |
| H3 | per-step leak (Gridap caches, `:full` frozen projections, VTK buffers) | RSS rising with step index | escalate; bisect by disabling `write_w`/`write_pressure`, then `nl_pressure` |
| **H4** | **baseline footprint simply exceeds 2 GB/core** | RSS flat and already > 2 GB after step 1 | **4 GB/core is correct and permanent** |

**H4 leads on local evidence.** On a *tiny* sequential case (3366–4669 DOFs) the first diagnostics
sample reads **1390 MB**, and a 400-step run went 2046 → 2089 MB (**+2 %, flat**). A Julia + Gridap +
GridapBALFEM process costs ≈1.4 GB before solving anything, so 2 GB/core leaves ~0.6 GB/rank for the
computation. Flat memory also means **H3 is not supported at this scale** — though a 32-rank
`:full` run exercises paths this case does not, so H3 is not eliminated.

⚠ Those are **local** numbers and must not be reported as the cluster answer: the cluster runs
against a ~1 GB sysimage (code mapped and shared) and MPI adds per-rank buffers.

> **NEW EVIDENCE FOR H3, 2026-08-22 — the flat-memory result does NOT hold at length.** The
> vertical-basis campaign ran 8-10 independent sequential processes for ~16 h. Every one grew
> steadily and none plateaued:
>
> | elapsed | typical process RSS |
> |---|---|
> | start | 1.5 GB |
> | +4 h | 2.0-2.5 GB |
> | +10 h | 2.5-3.1 GB |
> | +14 h | up to **3.9 GB** |
>
> Machine swap went 0 → 630 → 996 → 1820 → 2171 MB over the same period and **three processes had to
> be killed** to keep the run off swap. This is a *sequential* Julia+Gridap process — no MPI, no
> distributed assembly — so it isolates the growth to the per-process solve path rather than to
> anything rank-related.
>
> That does not confirm H3 by itself (a long-lived process accumulating Gridap caches and compiled
> specialisations is not the same as a per-step leak), but it **removes the "flat memory ⇒ H3 not
> supported" argument above**, which was measured over 400 steps on one tiny case. The decisive
> measurement is unchanged — RSS against step index within a single run — and it can now be taken
> locally, cheaply, from any of these shard logs, without cluster access.

**Evidence still to collect (needs cluster access):**

```bash
sacct -j <jobid> --format=JobID,JobName,State,ExitCode,MaxRSS,MaxRSSTask,AveRSS,Elapsed
seff <jobid>
grep -c '^\[dist\] step' <job>.out        # how far it got before the kill
grep -i 'sysimage\|stale\|-J ' <job>.out  # was the image used, was it fresh
```

Exit code **137** = OOM-kill; an `oom-kill event` line in `.err` confirms it.

**Decision rule.** Drop to 2 GB/core only if global-max RSS stays below ~1.6 GB/rank for a
full-length run of the largest small-domain case **and** the peak is not a pre-step-1 compile spike.
Any monotonic growth (H3) escalates instead — it means run length is bounded by memory, which no
per-core request fixes. **Cheapest decisive job:** `run/dist_small/run_lin_periodic_plane_small.sh`.

---

## 2. 🔴 The MMS path never assembles the `:full` frozen projections

Found 2026-08-21 while designing the vertical-basis campaign. **Reported, deliberately not patched**
— it is a design-level question about an existing interface, not a usage error.

`run_time_loop` defaults `nlp = nothing` (`src/timeloop.jl:130`) and **`run_mms_case` never passes
one** (`src/mms_driver.jl:99–105`); only `setup_and_run` builds a context (`src/utilities.jl:863`).
`src/problem.jl:396` gates the frozen `{1,2,4,5}` contribution on `st !== nothing`, so through the
MMS driver the `nl_pressure_full` branch adds **nothing at all on a flat bed**, and only the `𝓐` ∇h
IBP half on a sloping one.

The studies stay valid — `mms_forcing` is selected by the same symbol and *does* compute `{1,2,4,5}`
exactly, so forcing and solver really do encode different operators — but **what the floor measures
is OMISSION, not LAG**. That invalidates two statements in `VERIFIED_SCOPE.md` §4: that the lag is
"negligible, now *measured* rather than argued" (it was never exercised), and that this floor
characterises a *production* `:full` run (where the context exists and the projections genuinely are
frozen-and-lagged).

**Decisive next step:** add an `nlp` kwarg to `run_mms_case` — or build the context internally when
`nl_pressure=:full` — re-run the two `:full` studies and compare. **Floor drops sharply ⇒ it was
omission; floor barely moves ⇒ the original wording was accidentally right.** Either way the
`:full`-is-not-MMS-verifiable conclusion probably survives; only its reason is in question.

---

## 3. 🔴 `test_mms_convergence` G7 — needs a specification decision, not a fix

Diagnosed in full in `TEST_SUITE.md` §6: the temporal window is contaminated at both ends, `η`'s
asymptotic window is one refinement wide, and the cause is structural (the `Q3/Q2` pairing puts `η`
in the lower-order space). **Deliberately not re-specified, because the fix is a runtime-cost
decision:** re-specify at ~3× runtime, gate the two fields differently, or refine to `nx ≥ 48`.

**Do not widen the ±0.3 tolerance.** Also fix **G10**, which already computes the per-field
contamination (3.4 % for `u` vs 0.5 % for `η`) but collapses both into one 10 % threshold.

---

## 4. 🟠 No MPI tier in `runtests.jl`

The distributed set — now **four** files, `test_mms_distributed_parity.jl` added 2026-08-21 — needs
`mpiexecjl` and is only ever run by hand. That cost **three stale reference constants**, all found
the first time they were re-run in a session.

**Fix:** an MPI tier that runs them when `mpiexecjl` resolves and **skips LOUDLY** when it does not.
A silent skip recreates the same blind spot. `MPI_TESTS` in `runtests.jl` already lists all four with
their rank counts, so the tier has its inventory; what is missing is execution and verdict capture
(from **gate output**, never exit codes).

---

## 5. 🟠 Preconditioner replacement — the single biggest performance item

GMRES needs ~300–770 iterations where a well-preconditioned solve of this size should need tens.
Both cheap drop-ins are measured dead (`CONFIGURATION.md` §5). Remaining, in order of effort:
**field-split/Schur** → **geometric multigrid** → **restricted Schwarz** (needs library support
GridapSolvers 0.7.1 does not expose).

**Cheap reproducible test bed:** the ring-vs-plane gap on an identical mesh (654–695 vs 451–517).
A preconditioner that closes it likely helps everywhere.

---

## 6. 🟠 Four run-output gaps

In priority order:

1. **The discrete-equation check does not run under the default integrator.** It is `is_theta`-gated,
   so `res_theta = NaN` under SDIRK. Highest value: it is the only *in-run* verification that the
   accepted state satisfies the discretised governing equations.
2. **No per-Newton-iteration linear counts.** `SolverMonitor` samples the `ConvergenceLog` once per
   `solve!`, i.e. per RK stage, so `gmres 451/517` is min/max over *stages*. A Newton step whose
   first linear solve is cheap and second expensive is invisible.
3. **No wall-time breakdown.** `t_solve` lumps Jacobian assembly and the linear solve together, so a
   55 s/step 2-D cost cannot be attributed from the log alone.
4. **The linear solve's achieved tolerance is not recorded** — only its iteration count. A solve that
   stopped on `atol` rather than `rtol` is indistinguishable from one that converged well.

---

## 7. Verification gaps

* **`Q2/Q1` and `Q4/Q3` velocity shortfall** — ✅ **answered for the vertical-basis case
  (2026-08-29): PRE-ASYMPTOTIC.** P1LFE-4's `pw_u` ran 3.33 → 3.45 → 3.78 over `nx = 8…64` and
  reaches **3.94 at `nx = 128`**, error still falling ~15× per refinement
  (`MMS_VBASIS_CAMPAIGN.md` §2.5). The method is the transferable part: **extend the ladder and read
  the pairwise SEQUENCE**, never a slope fitted over a window that never reached the asymptotic
  regime — the fit read 3.515/3.729 against a true 3.94.
  ⚠ Still open for the **horizontal FE pairings** (`Q2/Q1`, `Q4/Q3`), which were not re-run on an
  extended ladder. The same diagnostic applies and is now cheap to run.
* **Error-vs-DOF study at production resolution** — **must precede any change to the default FE
  pairing.** `Q3/Q3` was 40× more accurate than `Q3/Q2` at `nx=24` despite the worse rate.
* **`:full` dynamic sensitivity has never been measured.** The block's relative effect scales as
  `A²`: 0.013 % at `A=1e-3` but **~1.3 % at `A=1e-2`**. The 2-D bar case already ran at `A=0.01`
  with `:full`; **its `:none` counterpart was never run**, and that single pair would both resolve
  the difference and confirm the `A³` scaling dynamically.
* **Only one vertical resolution and one `kd`** exercised locally (P1LFE-2, `Nσ=3`, `kd=5.5`).
  The convergence half of this is designed AND now **unblocked** — see
  [`COMPLETED_VBASIS_STUDY.md`](COMPLETED_VBASIS_STUDY.md) §1: both `src/mms_driver.jl` defects were fixed 2026-08-21,
  every MMS driver/test/example takes `(M, p_vert)`, and the sweep driver is
  `examples/local_mms/run_vertical_basis_study.jl`. **Decisive next step: run tier 1.** Two spot
  checks already reach theoretical order off the default basis (P1LFE-3 and P2LFE-2, Q2/Q1 1-D
  static Model 1: `p_η` 1.994, `p_u` 2.997/2.996).
  A `T` sweep measuring sponge reflection at `kd = 1, 3, 5.5` would directly test the requirement
  that the sponge cover the longest component — that half is still untouched.
* A run-and-reconstruct **pressure** profile test; a **distributed-gauge** utility.

---

## 8. Runs not yet made → moved

The list of runs that are designed and unblocked but not yet executed now lives in
[`PLANNED_CAMPAIGNS.md`](PLANNED_CAMPAIGNS.md) §7, which is the single place forward work is
tracked. **This file is for gaps in work already done**; a run that has simply not been launched is
not a gap, it is a plan.

---

## 9. Naming follow-through (outside this checkout)

The `GridapLFEM → GridapBALFEM` rename is complete, committed and pushed. Remaining:

* the **GitHub repository** is still named `GridapLFEM.jl` — rename it, then
  `git remote set-url origin …GridapBALFEM.jl.git`;
* the **cluster checkout** must be renamed to match (`run/balfem_env.sh` defaults
  `BALFEM_PROJ=$HOME/GridapBALFEM.jl`);
* the **sysimage must be rebuilt** as `GridapBALFEM_sysimage.so`.

**Deliberately NOT renamed:** `../LFE-M_2D_solver/` and its `LFEModel2D` module (a real legacy
directory), and `CFC2027_LFEMultilayer_abstract.pdf` (an existing artefact).

---

## 10. Rules for anyone doing a repository-wide rename here

Three earned the hard way, and all three are cheap to respect:

1. **Use a single regex alternation, longest-first, with protected tokens mapping to themselves.**
   An *ordered cascade* of find/replace rules re-matches its own output
   (`GridapLFEM → GridapBALFEM → GridapBABALFEM`).
2. **Renaming file CONTENTS and file NAMES are two jobs.** 39 launchers ended up sourcing a helper
   under its new name while the file still had the old one. Verify by **resolving every sourced
   path**, not by grepping.
3. **A blanket rename corrupts semantics exactly where the distinction matters most** — it rewrote
   the block *defining* the three model names. Guard such passages with an `assert count == 1`.
