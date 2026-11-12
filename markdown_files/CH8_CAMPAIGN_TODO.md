# Chapter 8 campaign: the working to-do list (opened 2026-10-31)

The ordered path from the current state to a running 1-D stability campaign. Work top to bottom.
When a task finishes, tick it, add the date and a one-line result, and continue with the next open
task. A task marked **ASK** waits for the author. Do not skip it, and while it waits, keep working
on later tasks that do not depend on it.

Standing constraints:
* No commits or pushes unless the author asks in that message.
* Stop and report on any design error in existing code.
* Runs longer than 30 min go detached, logging into `output/`.
* **Memory is the only cap on parallel runs** (2026-11-10; the 3-process thermal rule of §10 is
  withdrawn). Give each one its own systemd scope (`workstation-ram-budget`).
* Any `src/` edit makes the sysimage stale (rule 25).

---

## ⇨ RESUME HERE (written 2026-11-07, before a planned machine shut-down)

> **2026-11-11: NEXT = the time-integrator dissipation study** (author's decision; plan, candidates and
> predictions in [`STABILITY.md`](STABILITY.md) §6b). CN has shown the discrete operator unstable;
> before any new weak-form term, measure what dissipative integration buys and costs. Rule 15 is
> rewritten accordingly (remedy runs report the carrier cost). Actions, in order:
> 1. **In flight** (launched 2026-11-10/11): `SDIRK_2_2` box at Δt = 0.04 for m3 and m5 n32
>    (`run_1d_sdirk_box.sh`) and at k₀A = 0.25 (`run_1d_sdirk_box_kA25.sh`); the M = 1 triad controls
>    (`run_1d_m3_M1_controls.sh`: M = 2 at T = 3.0 s, M = 1 at n64). Analyse with `energy_history`
>    (groups + modes 14–18, 32) against the CN twins.
> 2. **No runs:** frequencies of the growing modes of every CN-unstable run (temporal spectrum), the
>    fastest discrete frequency (`stability_eig.jl`, explicit limits), and the damping/carrier-cost table
>    per candidate tableau and Δt (STABILITY §6b).
> 3. **ASK: launch the cluster study** `run/integrator_study/` (2026-11-11, replaces the local
>    m3 ladder): 170 box cases in tiers P (pilot) A–H, 23 packs (F, G, H and γ added 2026-11-12:
>    see "Integrator study: hypotheses" below). m5/m6 run CN, SDIRK_2_2, SDIRK_4_3,
>    SDIRK_Crouzeix_3_4 (re-added 2026-11-12 as RK4's stand-in); m3 (tier E) runs CN, SDIRK_2_2,
>    SDIRK_4_3, RK4, RK3 (TR-BDF2 and generalised-α dropped 2026-11-11). × kA up to
>    0.32 × kd 5.5/3.0/1.7 × n16/32/64 × dt. Packs of 8 on rome 1/8 nodes, automatic checkpoint segments.
>    **2026-11-12: explicit RK is unaffordable with the CIP** (damping modes to −9.5e3 s⁻¹ at n32, so
>    Δt ≲ 2.6e-4 s, 150× below T/40; STABILITY §6b): the explicit pair left m5/m6 (author).
>    ⚠ m3 k₀A = 0.25 under SDIRK_2_2 diverged at 172 s (growth from ≈ 144 s while the carrier had
>    decayed to E = 0.24): at kA ≥ 0.25 the less dissipative schemes are expected to fail sooner on
>    m3; m5/m6 (CIP) are not predicted by it (2026-11-12).
>    Order: sysimage build → pilot (P) → E → A → C → B → D. README in that folder.
>    (The 2-step smoke of every tableau: `test/local/smoke_integrators.jl`.)
> 4. **ASK:** the hard case m3 T = 1.4/1.5 s n64 (CN 8.7/12.9 s) under the best candidates.
> 5. **ASK:** m5/m6 n32 with the best one or two (cost: m5 ≈ 1–2 weeks per run under SDIRK; the
>    Gridap performance mode, the speed-up kept, is expected to cut ≈ 20 %; the depth-field memo was
   reverted 2026-11-12 after its benchmark showed no gain).
> 6. Write the outcome in STABILITY §6b and LaTeX chapter 9. Weak-form stabilisers (STATUS §4, SUPG/VMS)
>    only if no integrator passes §6 within the accuracy budget.
>
> ⚠ Excluded tableaus: Gridap's `SDIRK_3_3` is first order and `SDIRK_LobattoIIIC_2_2` is explicit
> Heun (library defects found 2026-11-11, not patched).

### Integrator study: hypotheses, tests and expected results (2026-11-12)

> **Design (author, 2026-11-12).** The cluster study (`run/integrator_study/`, 170 cases, tiers P–H)
> now varies three things: the integrator, Δt and the C⁰-CIP penalty γ. Models 3–4 run **without CIP**
> (their weak form has no non-conforming term to stabilise); models 5–6 keep the mandatory CIP and vary
> γ. The aim is to find the **stability region** of the solver in (scheme, Δt, γ) at each steepness.
> Before submitting, two local checks (`run/local/run_precluster_checks.sh`) test the premises.

**What the design rests on (measured, or derived from the tableaus).**
* *The m3 mechanism is a decay triad* (2026-11-08, below): the carrier's second harmonic (n = 2,
  2ω′) pumps a **steady** vertical-shear member n_s (ω ≈ 0; n_s = 15–16 at n32, λ ≈ 2 cells) and a
  **fast** member n_s − 2 at ≈ 2ω′ travelling −x. Under CN both grow at Γ ≈ 0.075 /s (|û|) at
  k₀A = 0.158, n32.
* *Integrator damping is frequency-selective.* The amplitude lost per step is ≈ c (ωΔt)^q: q = 4 for
  SDIRK_2_2, SDIRK_4_3 and RK3 (3rd order or less), q ≈ 6 for Crouzeix 3/4 and RK4 (4th order).
  At ω = 0 it is zero for every scheme, so **no integrator can damp the steady member.**
* *Damping one member of a triad only slows it.* With the steady member undamped and the fast member
  damped at δ, the pair grows at σ_eff = −δ/2 + √(δ²/4 + Γ²) ≈ Γ²/δ. That is positive for every δ:
  the integrator buys time, it cannot stop the triad.
  Check against the runs: SDIRK_2_2 at T/40 damps 2ω₀ at δ = 0.076 /s, so σ_eff = 0.046 /s at
  k₀A = 0.158. The 16 e-folds from the 1e-8 seed to O(0.1) then take ≈ 350 s, and the run survived
  320 s ✓. At k₀A = 0.25, Γ ≈ 0.19 /s if the pump scales with A² (an estimate), giving ≈ 100 s;
  the run died at 172 s. The order is right; the carrier's own decay weakens the pump and delays it.
* *Models 5–6 are different.* The CIP damps by **wavenumber**, so it can reach the steady member at
  the grid scale. With both members damped (δ_s by the CIP; δ_f by the CIP and the integrator), a
  pumped triad has a **true threshold**: it is stable when δ_s · δ_f ≥ Γ².
  This is the reason to vary γ and the integrator together.

**H1 (tier F, m3, no CIP): aim 4th-order damping at the 2ω₀ member.**
A 4th-order scheme damps 2ω₀ about 50–60× more than the carrier; a 3rd-order one only ≈ 15×.
Δt sets the level. From the tableaus (amplitude, T = 1.6 s; Γ = 0.075 measured, 0.19 and 0.31 for
k₀A = 0.25 and 0.32 scaled ∝ A²):

| scheme, Δt | carrier left after 200 T | δ at 2ω₀ (1/s) | 2ω₀ : carrier damping | σ_eff (1/s) at k₀A 0.158 / 0.25 / 0.32 | predicted blow-up (s) |
|---|---|---|---|---|---|
| CN | 100 % | 0 | — | 0.075 / 0.19 / 0.31 | 213 / 84 / 52 |
| SDIRK_2_2, T/40 (E) | 18 % | 0.076 | 14 | 0.046 / 0.156 / 0.274 | 347 / 103 / 58 |
| Crouzeix, T/40 | 95 % | 0.008 | 51 | ≈ CN | 225 / 86 / 52 |
| Crouzeix, T/20 | 29 % | 0.114 | 29 | 0.037 / 0.141 / 0.258 | 430 / 113 / 62 |
| Crouzeix, T/10 | 0 % | 0.546 | 10 | 0.010 / 0.060 / 0.140 | 1583 / 269 / 114 |
| RK4, T/20 | 97 % | 0.005 | 62 | ≈ CN | 221 / 85 / 52 |
| RK4, T/10 | 44 % | 0.140 | 55 | 0.033 / 0.132 / 0.248 | 492 / 121 / 65 |

*Expected.*
* (i) At k₀A = 0.158, Crouzeix T/20 and RK4 T/10 survive the 200 periods. They keep more carrier
  than SDIRK_2_2 (29 % and 44 %, against 18 %), and **RK4 T/10 is the best trade-off.** Crouzeix
  T/40 and RK4 T/20 behave like CN.
* (ii) At k₀A ≥ 0.25 every scheme fails before 320 s unless it destroys the carrier. That would be
  the measured limit of the integrator route on m3.
* (iii) Selectivity falls as Δt grows (Crouzeix 51 → 10 from T/40 to T/10), so there is a best Δt;
  larger is not always better.
* *Falsified if* a scheme survives k₀A = 0.25 with > 30 % of the carrier. That would mean the triad
  is not the limit there, or Γ is far below the A² estimate. Check 1 decides which.
* Blow-up times are rough: the onset delay is ignored, and the pump weakens as the carrier decays.

**H2 (tier G, m5, γ map): the CIP reaches what the integrators cannot.**
Raising γ damps the grid-scale content at any frequency, the steady member included. Together with
an integrator damping the fast member, a threshold δ_s(γ) · δ_f ≥ Γ² can be crossed.
*Expected.*
* Under CN, the growth slows as γ rises. Above some γ*, the k₀A = 0.158 n32 burst (245 s at γ = 0.3)
  disappears.
* If the growth **escapes** to the mid band instead (wavelengths of several cells, where the CIP is
  weak), it comes later and slower. v1 saw +0.13 /s; the v2 frozen analysis predicted the same.
  The selective integrators then decide.
* In each k₀A plane, the stable region lies at large γ and strong damping, with a boundary roughly
  γ × (damping at 2ω₀) ≈ const. It moves to larger γ as k₀A rises (Γ ∝ A²).
* γ costs the carrier almost nothing. The logged CIP drain is ≈ 1e-8 /s at γ = 0.3 and scales with
  γ. Any accuracy cost shows on the resolved short waves (n = 2, 3), not on the carrier: compare
  them with tier A.
* γ = 0.1 should be worse than 0.3 (earlier growth). It brackets the boundary from below.
* A stiffer penalty may slow Newton; the pilot's CN γ = 3 case measures s/step and the Newton count.

**H3 (tier H, m5 k₀A = 0.32, γ × Δt): once γ has done its part, Δt sets the integrator's share.**
Damping per unit time scales as Δt^(q−1): 8× from T/80 to T/40 for a 3rd-order scheme.
*Expected.*
* At γ = 0.3 (tier D) every scheme fails at k₀A = 0.32.
* At γ = 1 and 3, SDIRK_2_2 and Crouzeix survive longer at T/20 than at T/80. Crouzeix keeps more
  carrier at the same survival.
* If T/80 survives as long as T/20, the integrator is not the active ingredient: γ is.

**Controls.**
* CN is the diagnosis pair in every cell (rule 15).
* Tier A (γ = 0.3, T/40) is the reference for G and H; tier E is the reference for F.
* Every remedy run reports its carrier loss next to its survival time (`CONFIGURATION.txt`,
  `energy_history`).

**Pre-cluster checks** (`run/local/run_precluster_checks.sh`; logs `output/precluster/`):
1. *What grew in the m3 k₀A = 0.25 SDIRK_2_2 blow-up* (`_logs/kA25_mode.jl`, temporal spectrum):
   the fastest-growing harmonics, their frequencies, and the triad test. From the measured growth,
   Γ = √(σ(σ + δ)).
   - Expected: the triad (n_s = 15–16), with Γ ≈ 0.15–0.2 /s.
   - If it is not the triad, H1's table does not apply and tier F is re-planned before submission.
   **Result (2026-11-12, `output/precluster/check1_kA25_mode.log` + `check1_followup.log`).**
   The run shows **two mechanisms, and the triad is not the one that killed it.**
   * *13–140 s: decay triads, held by SDIRK_2_2.* Steady members at n = 12, 15, 20 and fast members
     n − 2 at −2.00 ω₀ are phase-locked (spread ≤ 0.01). They grow at 0.06–0.08 /s early, falling to
     0.01–0.02 /s as the damped carrier weakens the pump (|η̂₂| 2.5e-3 → 1.2e-3). By 150 s they are at
     only ~1e-6. Measured Γ = √(σ(σ+δ)) ≈ 0.1 /s, not the 0.19 of the A² estimate in H1's table.
   * *≈ 150–166 s: an explosive burst.* It grows at σ ≈ 0.7–1.0 /s at n = 13–18 (≈ 2 cells) in η, ū
     and u₃ alike, from ~1e-6 to ~1e-2 in ≈ 10 s. Its frequencies are spread over 0.1–1.5 ω₀ and the
     triad lock breaks. It is **oscillatory at or below the carrier frequency**, so no integrator can
     damp it selectively: damping it costs the carrier as much.
   * The trigger is not identified. The triad band was still at ~1e-6. One candidate is the
     sub-element band (n = 60 and 92, images of n = 4 or 28), which grew at ≈ 0.03 /s to ~3e-5 just
     before the burst.
   * Consequences for H1: at k₀A = 0.158 the triad remains the mechanism, and tier F tests it as
     designed. At k₀A ≥ 0.25 the expected failure is this burst, not the triad outrunning the damping;
     tier F there measures how soon the burst comes when less carrier is drained. The burst sits at
     the grid scale (n ≈ 16 of 32 cells), which is where the CIP acts. That supports the γ tiers on
     m5/m6. **No m3 + CIP runs (author, 2026-11-12):** if m3/m4 need grid-scale stabilisation, it
     will come from another method (SUPG, VMS), not the CIP, which stays tied to the 𝓝 models.
   * The script's automatic triad pick (n = 49 → Γ = 0.60 /s) is wrong: n = 49 is a sub-element image
     taken in the burst window, and its partner is at noise level. Ignore that line of the log.
2. *Frozen-state γ sweep of m5 n32* (`examples/local_1d/frozen_gamma_sweep.jl`, k₀A {0.158, 0.25,
   0.32} × γ {0.1 … 10}): σ_max(γ), where the growth moves (n, ω), and the verdict of each
   (scheme, Δt). It sets `GAMMAS_G` to bracket the predicted γ*.
   ⚠ A frozen state misses the time-dependent triad. It screens the grid-scale growth and the
   damping budgets; it gives no verdict on the run.


> **2026-11-08: the m3 mechanism and the period test.** The temporal-spectrum analysis of every m3 run
> found a **decay triad**: the carrier's second harmonic (n = 2, 2ω₀) feeds a steady vertical-shear
> pattern (n_s, ω ≈ 0) and a backward short surface wave (n_s − 2, ω ≈ 2ω₀ ≈ ω∞); equal growth rates,
> frequency closure, locked triad phase. Record: `output/1D_stability_campaign/README.md` ("What the m3
> growth oscillates at"); LaTeX chapter 9 §9.2 (`sec: stabilisation triad`, 2026-11-08, not pushed to
> Overleaf). **Running since 19:02: the period test** (`run/local/run_1d_m3_period.sh`, 13 cases,
> MAXP = 10, output `output/1D_period_test/`). When it ends: band growth with periodic_growth.jl, the
> triad test with `temporal_spectrum.jl` (`ts_triad`, pair 9/−7 at T = 1.7 s), then fill the "result"
> column of `tab: triad period test`. Sponge stage B finished 2026-11-08 11:49 (not analysed yet).
>
> **2026-11-09 11:30, interim.** Sponge stage B analysed: relaxation 12 m gives the same Kr as 6 m to
> < 0.1 % (`SPONGE_SWEEP_RESULTS.md` §0b). Period test: the T = 1.6 n64 control reproduced the
> October divergence exactly (44.56 s). **T = 1.4 and 1.5 at n64 diverged FASTER (8.72 s, 12.92 s),
> against the triad prediction** ("closed" periods). It is not the triad: the same grid-scale band
> n ≈ 40–50 (≈ 1.3–1.6 cells per wavelength) grows in all three n64 runs, but at T = 1.4/1.5 every
> harmonic of the band grows together at σ ≈ 1.2–1.7 /s, with mixed frequencies (single-frequency fit
> residual 1.5–2.3) near ω₀, and no steady/fast pair; at T = 1.6 the band shows the clean triad
> (n 47 steady, n 45/−48 at 2ω₀, σ = 0.28, residual < 0.4). So at n64 a second, faster, broadband
> grid-scale instability opens as T shortens (onset 3.6 / 4.4 / ≈ 14 s at T = 1.4 / 1.5 / 1.6).
> Still running at 11:23: T1.7/T1.8 n64 (≈ 175 s, alive), all n32 (211–231 s, alive), T1.4/T1.5 n16
> (181–216 s, alive); T1.7/T1.8 n16 queued.

**State.** Nothing is running. The stability campaign (§5) finished, and every case was analysed:
results in `output/1D_stability_campaign/README.md` ("FINAL RESULTS") and `STABILITY.md` §3b. The
sponge sweep (§7, §9.2) is **paused** at 26 of 84 cases: results and resume procedure in
`SPONGE_SWEEP_RESULTS.md`. The main tree is committed (2026-11-07, not pushed).

> **LAUNCHED 2026-11-07 12:48–12:56 at the author's request, with MAXP = 8 (author: "fill at least
> the 8 cores"; this overrides the MAXP = 3 throttle rule of §10 for this batch; MIN_FREE_GB = 4).**
> * 3(b), n8 only: `m3_q32_n8_cip` (new launcher `run/local/run_1d_m3cip.sh`; γ = 0.3, all
>   checkpoints kept), log `output/1D_stability_campaign/_logs/m3_q32_n8_cip.log`.
> * 5: the sponge sweep was resumed (RESERVED-STOP claims removed; 7 workers), `_batch_A.out`.
> * 3(a): `run_1d_nocip_controls.sh` was started after the 8 slots filled. It waits in the slot gate, so the
>   5 no-CIP cases start as sponge workers exit (`_nocip_batch.out`).
> * Check the clock (§10). At 12:56, during JIT, it read 3.4 GHz / 22 W.

**Next actions, in order.** Each needs the author's go-ahead to launch. Slots are limited by memory only (2026-11-10, §10).
1. ~~**ASK: commit** the pending tree.~~ **Done 2026-11-07** (8aa27f7, 58071da, 2348286 on `v2-solver`, not pushed; the no-CIP worktree `../GridapBALFEM_nocip.wt` stays uncommitted by design). That covers the checkpoint/bed/solver work, the launchers,
   `_slot_throttle.sh` (MAXP default then 3; memory-only since 2026-11-10), the sponge files, the no-CIP launcher and these docs.
   The run manifests otherwise keep recording a dirty tree.
2. **ASK: the assembly speed-up (§1)** before any further 𝓝 run. Models 5 and 6 cost 52–69 s/step
   against 4–5 for model 3. The measured `_get_cache` fix gives −29 %, and Gridap performance mode
   should remove another ~20 %. A 320 s model-5 run takes ~3.5 days even at full clock.
   **Outcome (2026-11-12, author):** option (b), the depth-field memo (`7c42a87`), was reverted: the
   A/B benchmark (`output/verify_depth_memo/bench_debug.log`, debug mode, bit-identical states) gave
   no gain, m5 0.99×, m6 0.77×, m3 0.97×. Building the depth fields once does not hit the cache-scan
   cost the 2026-10-31 experiment removed. **The speed-up kept is Gridap performance mode** (option
   (a), `compile/set_preferences.jl`, the cluster sysimage default). Option (c) is not pursued.
3. **Chapter-8 discriminating runs** (questions in the campaign README, "Open questions"):
   * (a) the no-CIP controls, §9.1 (prepared, never launched). They give the missing CIP axis for
     model 6 and for m5 n8 past 160 s.
   * (b) model 3 *with* the CIP at n32 and n8: separates 𝓝 from the CIP (rule 14c).
   * (c) resume m6 n32 and m3 n8 from their t = 320 s checkpoints to ~480 s: does the onset saturate
     or diverge?
   * (d) the n32 anomaly: an n24/n48 ladder, and n32 at dt = 0.02.
4. **Analysis without runs:** the dispersion of the vertical-structure (ω∞) branches about the
   carrier. Is the vertical-structure growth a continuum (model-PDE) resonance or a discrete one?
   The per-layer tool is `output/1D_stability_campaign/_logs/vertical_structure.jl`. All 11 cases are analysed.
5. **Resume the sponge sweep** (≈ 12–15 h at 2–3 workers). First the replica check (W = 12, μ = 40
   must reproduce Snellius), then stage B (§9.3), then the §7 "Apply" decision before Snellius round 2 (§8).
6. **LaTeX chapter 8:** add the Q3/Q2 rows (STATUS §6) and the campaign results (§3b).
   Chapter 9: the stabiliser must act on every u_j, and the rates to beat are 0.02–0.24 /s with
   onsets at 140–300 s.

## 0. Already in flight (check these first)

- [x] (2026-10-31; P1 measured with the profiler on; P2 crashed because the driver had gained the
  checkpoint kwargs mid-run, which is harmless) **Profile of model 5 in the 8-cell box.**
  * Scope `prof-m5-box2`, log `output/profile_m5_box/profile.log`.
  * Outputs: `profile_flat.txt`, `profile_tree.txt`, and the `>>> SUMMARY` line (P1 with VTK and
    diagnostics vs P2 without).
- [x] (2026-10-31) PASS: 2 steps converge (2 Newton iterations each), bed d ∈ [3.437, 3.563] m, slope
  0.100, 480 free DOFs, about 29 s/step (the same assembly overhead as model 5). **Smoke test of model 6 with the cosine bed, 2 steps.**
  * Scope `stab1d-smoke`, log `output/1D_stability_campaign/_smoke/m6_n8.log`.
  * Pass: both steps converge, the `BED` banner line is printed, and there is no error.

## 1. The assembly overhead (from the profile)

Evidence: about 44 s of assembly per Newton iteration for 14 k DOFs (`krylov_conditioning/mpi_check`),
and 31 s/step in the 8-cell box. The step time is nearly independent of the cell count, so the cost
is fixed, not proportional to the problem size.

- [x] (2026-10-31) **Finding.** About 54 % of all samples are spent BUILDING Gridap's lazy-array
  caches (`array_cache`), not evaluating anything. That cost is per assembled term, not per cell,
  which is why the step time barely depends on the mesh. It has two parts:
  * ~30 %: `Gridap.Arrays._get_cache` (LazyArrays.jl:155). When a node's objectid is not in the
    cache dict, it scans the WHOLE dict with a structural `same_branch` comparison. The cost is
    O(nodes²) per term, and the node that keeps missing is `∇(d_cf)` of the depth CellField:
    `alg_dx(d_cf)` creates a new node at every call site, even on a flat bed (`0.0*alg_dx(d_cf)`).
  * ~20 %: Gridap's debug-mode `@check` (Helpers/Macros.jl:61). This is the default execution mode;
    `Gridap.Helpers.set_performance_mode()` removes it.
  Experiment `output/profile_m5_box/cache_experiment.jl` (no repository change): the same run with
  `_get_cache` reduced to the objectid lookup. **RESULT (2026-10-31):** 35.9 → 25.5 s/step
  (−29 %; 18.0 → 12.7 s per Newton iteration), and the answer is **bit-identical** (final max η
  equal to all 16 digits). The machine was loaded by the campaign's compilation, so the absolute
  numbers are noisy; the ratio was measured in one process under the same load. Next, while the
  campaign runs: Gridap performance mode in a separate environment copy (the @check share), and
  the number of DAG nodes per term.
  Read the profile and name where the time goes: residual, `jacobian_u` or `jacobian_u_t`;
  volume terms or skeleton terms; repeated compilation or construction; VTK and diagnostics.
- [ ] If the cause is fixable without changing validated maths:
  * fix it;
  * re-run the profile to measure the speed-up;
  * re-run `test_jacobians_ad`, `test_broken_formulation` and the regression snapshot. The answer must
    not change.
- [ ] If the fix touches validated code or its design: **ASK** before changing anything.
  The options for the author:
  * (a) Gridap performance mode, a preference with no code change;
  * (b) build the depth CellFields (d, ∇h, ∇∇h, and constant zeros on a flat bed) once per problem,
    in `problem.jl`/`broken.jl`/`nlpressure.jl`, then verify with the regression snapshot and
    `test_jacobians_ad` — **tried (`7c42a87`) and reverted 2026-11-12: no gain (item 2 above)**;
  * (c) make `_get_cache` objectid-only in the Gridap fork, a library change.

## 2. Checkpoint and resume (u_n only; the author chose this on 2026-10-31)

Design (agreed):
* A checkpoint holds t, the step number, the free DOF values of u_n, the git commit and a
  configuration fingerprint.
* Resuming rebuilds the problem, refuses a fingerprint mismatch, and continues from t0 to `T_final`.
* It does **not** store Gridap's Newton-guess buffers (u_{n−1} for CN, the stage slopes for SDIRK).
  A resumed run is therefore a valid continuation but not bit-identical: it differs at the Newton
  tolerance level. Measure that difference; never assume it.

- [x] **Format** (2026-10-31): raw little-endian Float64 plus a hand-written TOML header, the
  same no-new-dependency style as `run_manifest.toml`. Files: `checkpoints/ck_<step>.{toml,f64}`,
  or `.part<k>.f64` per rank. The header is written last and atomically, and it stores Σv and Σv².
- [x] (2026-10-31, `src/checkpoint.jl`, both time loops) **Writing (sequential).** `checkpoint_every` kwarg in `setup_and_run` and `run_time_loop`:
  * write `checkpoint_<step>.jld2` every N steps, at the end of the run, and at divergence;
  * keep only the latest K, to bound disk use.
- [x] (2026-10-31: `CHECKPOINT_HARD_KEYS`/`_SOFT_KEYS`, the free DOF count and a 7×3 bathymetry sample) **Fingerprint.** Every quantity that fixes the discrete problem:
  * M, p_vertical, c_bdy, p_u, p_eta, domain, partition, periodicity and BCs;
  * regime, nl_pressure, flat_bed and a bathymetry sample;
  * the CIP γ/order/hexp and the stabilization;
  * dt, solver_type, tableau and θ;
  * the free DOF count.
  `dt`, the integrator, `nl_tol` and the output cadence are allowed to change only with
  `restart_allow_change = true`, and the manifest then records that change.
- [x] (2026-10-31; mass0/energy0 are carried over, so diagnostics stay relative to the original t = 0) **Resuming (sequential).** `restart_from = <checkpoint>`:
  * build everything as usual, verify the fingerprint, and set u0 from the stored free values;
  * t0 = t_checkpoint, and the step counter continues from there;
  * forcing (source, BC generation, ramp) is evaluated at absolute t, so check that the time loop
    receives t0 everywhere.
- [ ] **Output continuity.**
  * Write to a new output directory (`unique_output_dir`, never overwrite).
  * The manifest gets a `[restart]` section: parent directory, checkpoint, step, t, and the commit
    of both the parent and this run.
  * `diagnostics.csv` and `band_energy.csv` continue in absolute t.
  * Teach `GridapBALFEMPost` / `periodic_growth.jl` to follow a chain of run directories. **(open)**
- [x] (2026-10-31; the part count is checked before MPI starts) **Distributed.** Each rank writes its own values plus the partition description. Resuming
  requires the same `cpu_grid` and refuses anything else. Re-partitioning is out of scope; write
  it down in OPEN_ISSUES.
- [x] (2026-10-31: `checkpoint_kw()` in `_dist_common.jl`, plus a `resume` name token) **Drivers.** `BALFEM_CHECKPOINT_EVERY` and `BALFEM_RESTART_FROM` in `_dist_common.jl`,
  `run_flume_1d.jl` and `run_periodic_1d.jl`.
- [x] **Gates:** `output/restart_gate/test_restart.log` (sequential, R1–R8): **16/16 PASS** (2026-10-31).
  Resumed vs straight at step 10: m3 CN 1.4e-8, m3 SDIRK 3.4e-6, m5 CN 5.5e-8. Every reloaded
  state is bit-identical, and all four refusals fire. The remaining open items: the runtests tier entry
  and the postprocessing chain.
  `test_restart_distributed.log` (2 MPI ranks, D1–D4): **4/4 PASS** (2026-10-31). The resumed run
  differs from the straight run by a relative 4.6e-6 at step 10 (model 3, forced flume, CN, `:schwarz`,
  nl_tol 1e-8 and ls_rtol 1e-5, so the GMRES tolerance sets the level). Then add `test_restart.jl` to the
  `runtests.jl` tier, and the MPI one to its list of MPI tests.
  **Gate `test/test_restart.jl`.** Cases: model 3 and model 5, in the box, under CN and SDIRK;
  sequential, and distributed on the debug backend.
  * One run of N steps straight through, against N/2 steps + checkpoint + resume for N/2 steps.
  * Report max|Δu| and |Δη|/|η| at the end; the expected size is the Newton-tolerance level.
  * Confirm that a resumed run with an unchanged state reproduces the time and step counters.
  * The refusals fire on a mismatched fingerprint and on a mismatched `cpu_grid`.
  * Verdict from the gate output, never from the exit code (rule 35).
- [ ] **Optional.** `restart_from_vtu` for old runs. Their sub-cell VTK holds every Lagrange
  node in Float64. Only if the gate shows that u_n-only resumption is acceptable.
- [ ] Docs: `README.md` (kwargs), `RUNNING.md` (how to resume), `CLAUDE.md` §3 features and §5.

## 3. Housekeeping before the rebuild

- [ ] Fast suite `test/runtests.jl`. It has not been run since the 2026-10-30 changes; read the
  verdict from the gate output.
- [ ] Re-run `kc4_ownblocklu.jl` if `timeloop_dist.jl` changed in task 2.
- [ ] **Docs for the solver changes of 2026-10-31:**
  * `CONFIGURATION.md` §5: `:schwarz` now works as `OwnBlockLUSolver`, with 15–19 GMRES iterations on
    the production case against Jacobi's cap. `:gs` is dropped.
  * `OPEN_ISSUES.md`: the preconditioner item.
  * `CLAUDE.md`: §3, rule 16, and the integrator list (`:gen_alpha` and `:rk3` removed).
  * `STATUS.md`.
  * `STABILITY.md` §3: model-3 Q3/Q2 n32 was **accelerating** at 160 s (+0.10–0.12 s⁻¹ over the last
    60 s), not just "growing +0.07". n16 was growing slowly (+0.02–0.03).
  * Rule 24 caveat: at Q3/Q2, Jacobi at `ls_rtol = 1e-5` leaves a 4–5 % solution error
    (`kc4_ownblocklu.log`).
- [ ] **ASK: distributed default preconditioner.** Keep `:jacobi` as the default, or switch to
  `:schwarz`? `:schwarz` is measured better on 1-D; on 2-D at more ranks it is unmeasured.
- [ ] **ASK: commit before the campaign?** The run manifests otherwise record a dirty tree (rule 2c
  reproducibility). The pending set:
  * the integrator removal;
  * `OwnBlockLUSolver` and the removal of `:gs`;
  * the periodic cosine bed;
  * checkpoint and resume;
  * `compile/compile_local.{jl,sh}` and the campaign launcher.

## 4. Local sysimage (optional, but it saves ~25 min of compilation per campaign case)

**Skipped for the campaign (2026-10-31, author: "prioritise launching").** `src/` changed that day,
and a rebuild costs about 4 h, more than the ~25 min of JIT each campaign slot pays. Rebuild after
the campaign, when `src/` has settled.

- [ ] Only once `src/` is final (tasks 1–3), and with nothing else heavy running.
- [ ] `compile/warmup.jl` dry run: no `gen_alpha`/`rk3`/`gs` cases; the `:schwarz` case is still there.
  Consider adding a model-6 x-periodic box case with the cosine bed, so the campaign path is traced.
- [ ] Build:
  `systemd-run --user --scope --unit=balfem-sysimage -p MemoryMax=24G nohup compile/compile_local.sh > output/sysimage_build/build.log 2>&1 &`
  * The 2026-10-30 build was OOM-killed at 16 G in the final compile, after 3 h 40 min of warmup.
  * If 24 G is not enough, reduce the warmup to the campaign path (sequential models 3, 5 and 6, CN, box).
- [ ] Pass: the warmup summary shows **0 FAILED**, `GridapBALFEM_sysimage_local.so` exists, and
  the `.src.sha256` stamp matches.
- [ ] Smoke: one 2-step box case with `-J`, checking that the first step takes seconds, not about
  25 min.

## 5. The stability campaign

- [x] **LAUNCHED 2026-10-31 14:29.**
  * Scope `stab1d-campaign`, logs in `output/1D_stability_campaign/_logs/`.
  * Source snapshot in `_logs/source_20261031_142905/` (commit 1ae2695 + tracked diff + untracked
    tarball), because the tree is uncommitted.
  * No sysimage; checkpoints every 250 steps, the 2 newest kept.
  * Resume a lost case with `BALFEM_RESTART_FROM=output/1D_stability_campaign/<case>`, deleting
    its log first (the launcher skips existing logs).
- [ ] Launch:
  `nohup run/local/run_1d_stability_campaign.sh 6 > output/1D_stability_campaign/_logs/_batch.out 2>&1 &`
  * 11 cases: m5 and m6 at n64/32/16/8, and m3 at n32/16/8.
  * Q3/Q2, CN, A = 0.10, 200 periods, seed 1e-8.
  * m6 has a cosine bed with slope 1:10.
  * Check `_batch.log`, which says whether the sysimage was used.
- [x] Checkpoints are on (launcher `CHECKPOINT_EVERY=250`).
- [x] (2026-10-31 → 11-07) Monitor: one status line per case (t, max η, Newton iterations, s/step,
  RSS); divergences, with t and x_at_max; memory against the 5 G per-case cap. Interim checks are in
  the campaign `README.md`.
- [x] (2026-11-07) Per run, when it ends: `periodic_growth.jl` (`_logs/final_growth_<case>.log`;
  driver `_logs/final_growth_all.jl`), band tables over time (`_logs/band_table.py`,
  `final_band_table_2026-11-07.txt`), and per-layer energies (`_logs/vertical_structure.jl`,
  `vertical_structure_2026-11-07.log`).
- [x] (2026-11-07) Compared against `model3_reference/` and the m3 reruns, and written into
  `STABILITY.md` §3b and the README ("FINAL RESULTS"). Reported to the author. **Next: LaTeX chapter 8
  (author).**

## 6. Parked (waiting for a decision)

- [ ] Re-launching the 1-D flume CIP-monitor runs (m5/m6, A = 0.02/0.05). Sequential LU, or
  4 ranks with `BALFEM_PRECOND=schwarz`, now that it works.
- [ ] 2-D with 𝓝: `:schwarz` at more ranks, which has no coarse space; check its iteration growth
  with the rank count.
- [ ] γ calibration of the CIP (chapter 8 route, step 2).

## 7. Sponge-layer reflection campaign (opened 2026-11-02; outside Chapter 8)

> **2026-11-07: stage A paused at 26 of 84 cases.** Results, deductions (Kr ≈ 0.21 μ_max/(ω k² W²)
> in the weak-damping regime; W = 24, μ_max = 8 passes at every T) and the resume procedure:
> [`SPONGE_SWEEP_RESULTS.md`](SPONGE_SWEEP_RESULTS.md).

**Why.** The Snellius linear production runs (`LINEAR_PRODUCTION_RESULTS.md`) reflect from the
12 m outflow sponge (μ_max = 40, quadratic profile). The reflection is 0.5 % of the amplitude at
T = 1.5 s, 5 % at 2.0 s, 9 % at 2.2 s and 18 % at 2.9 s. It shows up as a ±5 % standing-wave
ripple over the whole interior at T = 2 s. The ring's 7 m sponges return about 10 %.

**Hypothesis.** The far wall is not the source: e^{−2 μ_max W/(3 c_g)} ≈ e⁻²⁰⁰ at T = 2 s. The
reflection comes from where μ rises at the front of the sponge. It grows with μ_max/ω and with
λ/W, which matches the measured rise with period. μ_max = 40 is 2.5–3.6× past the saturation
μ_max ≈ 5ω (rule 11), so the extra damping buys nothing and adds reflection. Expected optimum:
μ_max ≈ 2–5ω of the longest component that matters, and W ≈ 1.5–2 λ of that component.

**Design** (linear model, cheap; run it when the stability campaign frees cores, or on Snellius):
- [x] (2026-11-06, prepared: §9.2) Driver: `examples/local_1d/run_flume_1d.jl`, sequential
  (`BALFEM_MPI=0`); run through workers (`run/local/sponge_worker.jl`, several cases per process,
  one 4 G scope per worker) instead of one process per case.
  * `BALFEM_WAVE_GEN=bc`, `BALFEM_REGIME=linear`, flat bed.
  * P1LFE-2, Q2/Q1 (`BALFEM_FE_ORDER=2`), `BALFEM_D=3.5`, dt = 0.02, CN (`BALFEM_SOLVER=theta`).
  * ny = 1 with walls (the driver default), dx = 0.25 m (NX = LX/0.25), `BALFEM_RELAX=1`.
  * These match the Snellius suite.
- [ ] Geometry: relaxation zone `BALFEM_RELAX_W` (6 m), then a FIXED 32 m interior, then
  `BALFEM_SPONGE_R` = W. So LX = 6 + 32 + W, and the measuring stretch (7–37 m) is the same in
  every case.
- [ ] Matrix, stage A: μ_max (`BALFEM_MUMAX`) ∈ {4, 8, 16, 25, 40} × W ∈ {6, 12, 18, 24} m ×
  T ∈ {1.5, 2.0, 2.5, 2.9} s, i.e. 80 runs. Each run lasts until a steady state, then 10 periods
  (rule 14: (x_sponge − x_inflow)/c_g + round trip + 10 T). 6 in parallel, about 9 h at an
  estimated 0.5 s/step. Write a launcher `run/local/run_1d_sponge_sweep.sh` in the style of the
  stability-campaign launcher (systemd scopes, LPT order, skip existing logs).
- [ ] Controls in the same batch (rule 38c):
  * **floor**: a very long, gentle sponge (W = 48 m, μ_max = 4). It measures what the inflow
    relaxation zone and the generation reflect on their own.
  * **round-1 replica**: W = 12, μ_max = 40. It must reproduce the Snellius Kr(T) (5 % at 2 s,
    18 % at 2.9 s), or the 1-D set-up is not representative.
- [ ] Stage B (after A): around the best (W, μ_max), run relaxation width
  `BALFEM_RELAX_W` ∈ {6, 12} m, since the inflow zone is also under one long wavelength.
- [x] (2026-11-06, written; first used on the smoke cases) Analysis script `postprocessing/examples/sponge_reflection.jl`:
  * complex amplitude at ω over the last 10 periods, from vertex η (SAVE_EVERY 5);
  * least-squares fit of incident + reflected over 7–37 m (the fit in `LINEAR_PRODUCTION_RESULTS.md`),
    giving Kr;
  * the residual of the fit, and the amplitude left at the wall end;
  * table and heat map of Kr(μ_max, W) per T.
- [ ] Pass criterion: Kr ≤ 1 % over T = 1.5–2.9 s (Tp = 2 band) at the smallest W that achieves
  it. If no quadratic sponge reaches it, report the best one and propose a smoother profile (cubic
  or exponential). That is a `src/utilities.jl` `make_sponge` change and needs the author's
  approval.
- [ ] Apply. All the distributed_small drivers already read `BALFEM_SPONGE_R` and `BALFEM_MUMAX`
  (and the directional one `BALFEM_SPONGE_Y`). If the outflow sponge changes:
  * the round-1 flat references must be rerun with it, or the big-bar comparisons lose their
    reference;
  * so decide this BEFORE submitting the big-bar launchers (§8).
  * The ring's 7 m sponges and the directional 2 m lateral sponges also need sizing.
  * Oblique incidence on the lateral sponges is not covered by a 1-D sweep. Follow-up: a 2-D
    plane-wave case at θ = 20–40°.

## 8. Linear production, round 2 (prepared 2026-11-02; `run/linear_production/README.md`)

- [ ] Decide first whether §7 changes the sponge. If it does, edit the launchers and add flat
  references.
- [ ] Push `v2-solver` (needs the author's go-ahead), pull on Snellius, rebuild the sysimage
  (`src/` changed; `BALFEM_STRICT_SYSIMAGE=1`).
- [ ] Pilot `:schwarz` at 42 ranks:
  `BALFEM_TFINAL=2.0 sbatch --time=02:00:00 run/linear_production/run_lin_directional_flat_schwarz.sh`.
  Read `gmres=` (it must stay well below the restart length krylov_m = 100) and s/step. If it is
  worse than Jacobi, fall back with `BALFEM_PRECOND=jacobi` plus checkpoint resubmission.
- [ ] Submit `submit_all.sh schwarz`, then `submit_all.sh bigbar`. Never run it without a pattern:
  that would resubmit all seven round-1 cases.
- [ ] Analyse as in round 1 (the merged-piece reader, bar/flat ratios, incident/reflected fits).
  For the big bar, add the phase and wavelength over the crest and, in the directional run, the
  refraction of the oblique components.

## 9. Queue for the freed cores (prepared 2026-11-06, NOT launched: the author launches)

> **State 2026-11-07:**
> * §9.1, the no-CIP controls: **never launched**.
> * §9.2, the sponge sweep: launched 11-06 16:37 and stopped twice for CPU throttling (§10).
>   Relaunched 11-07 02:25 at MAXP = 3, then paused at 26 of 84 cases at ~11:00 for the machine
>   shut-down.
> * §9.3: not started.

State on 2026-11-06 ~10:00: the campaign's 11 cases have all started. Five are still running
(m6 n64, m5 n8, m6 n8, m3 n16, m3 n8). At 11:45, m3 n16 and m3 n8 had reached about 200 s of
320 s, so they finish at about 15:00 if they do not diverge first. The three batches below are
prepared and smoke-tested. They all go through one shared slot gate (`run/local/_slot_throttle.sh`:
at most `MAXP` = 6 local 1-D Julia processes in total, and `MemAvailable` ≥ `MIN_FREE_GB` = 6). So
each launcher can be started at any time: it waits for free slots and never oversubscribes the
machine, even with the campaign and the other launcher running beside it.

### 9.1 The no-CIP controls (Chapter 8; launch first)

Five cases with γ_u = γ_η = 0, otherwise identical to the campaign: `m6_q32_n{8,16,32,64}_nocip`
(no model 6 without CIP exists) and `m5_q32_n8_nocip` (the October no-CIP m5 n8 stopped at 160 s).
* Code: the git worktree `../GridapBALFEM_nocip.wt` (branch `diag/nocip-control`) = 68919a3, the
  parent of the mandatory-CIP commit faa094f, plus the checkpoint feature and the periodic bed,
  applied **uncommitted**. Why this base and not HEAD minus the refusal: the launcher header.
  `Manifest.toml` is copied from the main tree, and `Gridap.jl` and `output` are symlinks into it.
  Do not edit the worktree while cases run.
* Checkpoints every 250 steps, all kept (`CHECKPOINT_KEEP=0`), so that a pre-onset restart is
  possible.
* Launch:
  `nohup run/local/run_1d_nocip_controls.sh > output/1D_stability_campaign/_logs/_nocip_batch.out 2>&1 &`
* Expected cost: n16–n64 should fail within minutes of simulated time, if they behave like m5
  without CIP (9 / 27 / 81 s). The n8 pair may run the full 320 s, about 2–3 days each.
* Afterwards: run `periodic_growth.jl` on each `output/1D_stability_campaign/*_nocip` case
  (`PG_NCELL=<n> PG_PU=3 PG_NLAMBDA=1`). Then tabulate the failure time and the band rates against
  the CIP case of the same name, and write the with/without table into `STABILITY.md` (model 5
  table: campaign `README.md`, 2026-11-06).
* Gate for the worktree: if the 10-step smoke test (`_smoke/m6_q32_n16_nocip_smoke.log`) matches
  the campaign case it should match (see the smoke record below), the worktree reproduces the
  campaign operator.

### 9.2 Sponge sweep, stage A (§7; launch second, or together)

84 cases (80 matrix + 4 floor), linear, cheap.
* Files:
  * `run/local/run_1d_sponge_sweep.sh` (launcher);
  * `run/local/sponge_cases.py` (case table; writes `output/1D_sponge_sweep/_logs/_cases_A.tsv`);
  * `run/local/sponge_worker.jl` (worker: one JIT per worker, up to 12 cases each, cases claimed
    through `_logs/_claims/`);
  * `postprocessing/examples/sponge_reflection.jl` (analysis).
* Launch:
  `mkdir -p output/1D_sponge_sweep/_logs && nohup run/local/run_1d_sponge_sweep.sh > output/1D_sponge_sweep/_logs/_batch_A.out 2>&1 &`
* Monitor `_logs/_batch.log`, `_logs/_workers.out` (one start/done line per case) and the per-case
  logs.
* Analyse:
  `julia --project=postprocessing postprocessing/examples/sponge_reflection.jl output/1D_sponge_sweep A`
  It writes `kr_A.csv` and the heat maps, and prints Kr(μ, W) per T, the floor and the replica.
* Check, in this order:
  1. the replica reproduces the Snellius Kr(T) (≈ 5 % at 2.0 s, 18 % at 2.9 s);
  2. the floor is well below 1 %;
  3. the fit residuals are below Kr/2 (`?` marks);
  4. then the pass criterion of §7.

### 9.3 Sponge sweep, stage B (after 9.2)

`STAGE=B BEST_W=<W> BEST_MU=<μ> nohup run/local/run_1d_sponge_sweep.sh > output/1D_sponge_sweep/_logs/_batch_B.out 2>&1 &`
* 4 cases: relaxation width 12 m at the best stage-A pair. Their own fit stretch is [13, 43] m.
* Analyse with argument `B`, and compare with the stage-A row at the same (W, μ), which has a
  relaxation width of 6 m.

### 9.4 Not queued (each needs a decision or other resources)

* The m3 n8 no-CIP 200-period run is the campaign case `m3_q32_n8`. It is already running (started
  2026-11-06 06:15), so nothing new is needed.
* §6 parked:
  * the 1-D flume CIP-monitor reruns (m5/m6, A = 0.02/0.05);
  * the 2-D `:schwarz` rank-count scaling with 𝓝;
  * the γ calibration of the CIP.
* §8 Snellius round 2. It depends on §9.2/9.3 for the sponge, and on pushing `v2-solver`.
* Diagnostic rerun of m5 n32 from just before its onset. It needs `restart_from_vtu`, an optional
  item of §2 not built yet, because the campaign kept only 2 checkpoints.
* Per-node u_j analysis of the m5 n8 VTK, to test the vertical-mode hypothesis. This is an
  analysis, not a run.

### 9.5 Smoke record (2026-11-06, run beside the 5 campaign cases with `MAXP` raised to 7/8)

* **No-CIP worktree.** `m6_q32_n16_nocip_smoke`, 10 steps.
  * Exit 0. The `periodic_1d [...] done: 10 steps` line is printed, along with the cosine-bed
    banner and the "nl_pressure=true WITHOUT a skeleton stabilisation" warning.
  * The manifest records `cip_gamma_u = cip_gamma_eta = 0.0`, branch `diag/nocip-control`.
  * Against the campaign `m6_q32_n16` (CIP), at step 5:
    * the energy ratio is identical (0.999653721);
    * the mass is identical;
    * `res0` is 2.8256e-2 vs 2.8216e-2;
    * η_max is 0.1005631 vs 0.1005222, a 4e-4 relative difference. At t = 0.2 s that is the size
      of the CIP's contribution (`cip_rate` ≈ 2e-5).
  * Conclusion: the worktree reproduces the campaign operator minus the CIP.
  * Cost: the first step (JIT) takes 31 min, about 46 s per step after that.
* **Sponge sweep.** Two workers, `spA_T1.5_W48_mu4` (nx 344) and `spA_T2.9_W6_mu40` (nx 176),
  50 steps each.
  * Both print `flume_1d [...] done:`. The claims, `_workers.out` and the verdict line all behave
    as designed.
  * Cost: JIT ≈ 6.5 min per worker, then 1.07 s/step and 0.61 s/step, i.e. ≈ 3.3e-3 s per
    cell-step (with 7–8 processes running).
  * Estimate for stage A (8.9e7 cell-steps): ≈ 80 CPU-h, ≈ 14 h on 6 slots. This replaces the
    0.5 s/step guess.
  * `sponge_reflection.jl` is not exercised by the smoke: 1 s holds no 10 analysis periods. Its
    first real test is the first finished stage-A case. Run it as soon as a few cases are done,
    before the sweep ends.
* The smoke outputs stay in `output/1D_stability_campaign/_smoke/` and `output/1D_sponge_sweep/_smoke/`.
  They do not interfere with the real runs (separate directories and claim directories).


## 10. The workstation throttles above 3 busy cores (measured 2026-11-06/07)

HP ZBook Power 16 G11 A, Ryzen 7 8845HS: 8 cores / 16 threads, 30 GB.
* **3 busy single-thread Julia runs:** 27–30 W, 4.6–4.8 GHz, 94–100 °C. Held for 20–25 min in two
  tests, no slow phase.
* **8 busy runs:** full clock for 1–3 min until the chip has sat at 100 °C, then the firmware drops
  every core to ~1.1 GHz and 9–10 W (54 °C) for ~5 min, and the cycle repeats.
  * Each run goes at roughly half speed.
  * The 𝓝 campaign cases went from ~38 to 110–160 s/step.
* **Not the cause** (all checked): the charger (the HP adapter changed nothing), the OS power
  profile (performance; toggling it changed nothing), lid or dock.
* Logs: `output/1D_sponge_sweep/_logs/_clock_*.txt`. Read the clock with `grep MHz /proc/cpuinfo` and
  the package power with `cat /sys/class/hwmon/hwmon4/power1_average` (amdgpu, µW).
* ~~**Rule: MAXP = 3**~~ **Withdrawn 2026-11-10** (author): the throttling was a machine fault that a
  reboot cleared and no longer occurs. Memory is the only cap; `_slot_throttle.sh` defaults to no process cap.
* Better cooling (raised rear, cooling pad) may allow more; measure before relying on it.
* The user has sudo. `ryzenadj` (not installed) could read the firmware limits. That was not pursued;
  the author called the throttle a firmware bug.
