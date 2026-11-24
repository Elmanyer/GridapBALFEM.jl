# INTEGRATOR_STUDY_RESULTS.md — the Snellius time-integrator study: results, conclusions, how to resume

*Written 2026-11-23 from the downloaded campaign (`output/snellius/integrator_study/`, gitignored). The
study's design, cases and hypotheses: [`run/integrator_study/README.md`](../run/integrator_study/README.md),
[`CH8_CAMPAIGN_TODO.md`](CH8_CAMPAIGN_TODO.md) "Integrator study: hypotheses", [`STABILITY.md`](STABILITY.md)
§6b. Analysis scripts: `run/integrator_study/analysis/` (tables in
`output/snellius/integrator_study_analysis/`).*

**The question** (author, 2026-11-11): Crank–Nicolson (CN) shows the discrete nonlinear operator unstable in
every closed-box run. How much of that growth does time-integrator dissipation hold, and at what cost to
the carrier wave? 275 box cases: models 3, 5, 6; kd 5.5 / 2.9 / 1.7 (T 1.6 / 2.2 / 3.0 s); kA 0.158 /
0.25 / 0.32; 16 / 32 / 64 cells per wavelength; Δt = T/20, T/40, T/80 (T/10 on m3); CIP γ 0.1–1; nine
integrators (CN, SDIRK_2_2, SDIRK_4_3, Crouzeix 3/4, TR-BDF2, generalised-α ρ∞ = 0 and 0.5, and RK3 / RK4
on m3 only). 200 periods each, P1LFE-2, Q3/Q2, seed 1e-8.

---

## 0. The campaign was cut by the home quota

| outcome | cases | meaning |
|---|---|---|
| DONE | 59 | reached 200 periods (not necessarily clean: §2) |
| GUARD | 33 | stopped by the divergence guard (|η| > 20 η_ref): a blow-up result |
| NEWTON | 22 | Newton failed: a blow-up result (NaN), except m3 Crouzeix T/10 kA 0.32 (a stall at step 3) |
| **QUOTA** | **161** | killed 2026-11-21 05:00 → 11-22: the VTK output filled the Snellius **home** quota (LightXML write error on `solution.pvd`, `Disk quota exceeded`; 7 without a marker) |

* The cut cases reached a median **49 %** of their target time. Their data up to the cut is valid and is
  used below. **All 161 newest checkpoints are intact** (length and checksums verified on the download).
* No case needed a second segment.
* Remaining cost to finish them, at the measured s/step: **≈ 1.8 × 10⁴ case-hours** (§6).
* Cause: ≈ 200 MB per case (VTK with sub-cell samples, 8 per period; 256 k files in all) written under
  `$HOME`. Since 2026-11-23 every cluster launcher writes to the project storage (`run/storage.sh`, §6).

**Cost per step** (m5 n32 T/40, Snellius rome, one thread): CN 33–40 s; generalised-α ρ∞ 0 / 0.5 31–45 s
(one implicit stage, ≈ CN); SDIRK_2_2 48–67 s; TR-BDF2 71–105 s; Crouzeix 3/4 99–179 s; SDIRK_4_3 140–179 s.
The Newton count per step of the multi-stage schemes grows with kA (SDIRK_4_3: 8 → 15).

## 1. Method

* **Box harmonics** û_n(t) of η, ū, u₁, u₂, u₃ from every VTK (exact Fourier transform on the periodic box).
* **The surface / vertical-structure split.** At harmonic n, the linear P1LFE-2 operator
  (A + k²d²B̃) u̇ = −ikg D η has one surface-wave velocity profile s = M⁻¹D (M = A + k²d²B̃) and two steady
  **vertical-structure (VS)** profiles with Dᵀu = 0 (no depth-average). They are M-orthogonal, so
  u = α s + v uniquely; |v|_M (m/s) is the VS amplitude. Bands: low n = 1–4, **grid** n = 5 … node Nyquist
  (3·ncell/2), sub-element above it.
* **Onset** of a burst: the first period whose period-rms of a band exceeds 10× its median over periods
  10–40. Carrier = |η̂₁|. Rates σ are amplitude rates (1/s).
* Just before each onset: the fastest-growing harmonics, their observed frequency (phase fit) and their
  η/v ratio (`burst_modes.log`).

## 2. What grows

1. **The burst is a steady, grid-scale, vertical-structure mode.** Before every m5/m6 onset the fastest
   harmonics are n ≈ 15–22 at 32 cells/λ and n ≈ 32–40 at 64 cells/λ: **a wavelength of ≈ 2 cells at both
   resolutions**. Their sub-element images n ± 32k (48, 50, 80, 82) grow with them. The observed frequency
   is ω ≈ 0 and |η̂|/|v| ≈ 0.03–0.2: velocity, not surface. Growth σ ≈ 0.10–0.16 s⁻¹ at n32
   (kA 0.25–0.32), 0.2–0.47 s⁻¹ at n64. Model 3 bursts the same way at kA ≤ 0.25 (n = 16–18 and images);
   at kA 0.32 at n = 11–14 with ω 0.1–0.7 ω₀.
2. **After the burst, m5/m6 under CN saturate rather than blow up.** All n32 T 1.6 CN runs reached 200
   periods ("DONE"), but they carry a grid-band VS velocity of 2–10 cm/s (≈ 10 % of the carrier's orbital
   velocity) and lost 12–28 % of the carrier energy (kA 0.25: E/E₀ 0.75), drained by the CIP.
   **Completed is not clean**: read the band content, never the survival alone (rule 38f).
3. **The steady VS field outlives the carrier, and nothing damps it.** SDIRK_2_2 and generalised-α 0 at
   T/20, kA 0.32: the carrier falls from 9e-2 to 3e-4 m, while the low-band VS velocity stays at
   ≈ 9e-3 m/s for all 200 periods. The grid band grows 3.5e-4 → 2e-3 m/s while the carrier is there, then
   stays; the sub-element band keeps growing slowly (2e-5 → 2.6e-4 m/s) after the carrier is gone. At T/40
   under the strong schemes the grid-band VS grows +0.006 s⁻¹ in absolute terms (1e-4 → 4.5e-4 m/s over
   150 periods): the slow steady growth found 2026-11-17 under SDIRK_2_2, now also under generalised-α 0.

## 3. Results by tier (m5/m6 at CIP γ = 0.3 unless stated)

### 3.1 The carrier cost (rule 15), measured

Carrier amplitude decay rate before any burst, m5 kd 5.5 n32 T/40 (CN: 0 to −1e-4 s⁻¹, CIP + nonlinear):

| scheme | SDIRK_2_2 | gen-α 0 | SDIRK_4_3 | Crouzeix 3/4 | gen-α 0.5 | TR-BDF2 |
|---|---|---|---|---|---|---|
| σ_carrier (s⁻¹) | 0.0056–0.0060 | 0.0038–0.0040 | 0.0003–0.0005 | 0.0002–0.0005 | 0.0002–0.0003 | 0.0001 |
| carrier energy kept over 200 T | 2–3 % | 8–9 % | 73–83 % | 73–88 % | 83–88 % | 94 % |

These match the tableau predictions of `run/integrator_study/README.md` (whose "carrier kept" column is the
amplitude). At kd 2.9 and 1.7 the rates scale with ω (T/40 fixes ωΔt). At T/20, SDIRK_2_2 and
generalised-α 0 keep 0.05 % of the carrier energy after 200 T.

### 3.2 Tier A: steepness and depth (m5, n32, T/40)

Burst onset period / periods reached (`–` = no onset within the record; q = cut by the quota):

| kd 5.5 (T 1.6) | CN | SDIRK_2_2 | SDIRK_4_3 | Crouzeix | TR-BDF2 | gen-α 0 | gen-α 0.5 |
|---|---|---|---|---|---|---|---|
| kA 0.158 | 167/200 | –/160 q | –/64 q | –/84 q | –/85 q | –/200 | –/175 q |
| kA 0.25 | 76/200 | –/146 q | –/58 q | –/73 q | –/82 q | –/185 q | 124/174 q |
| kA 0.32 | 39/200 | –/135 q | 48/55 q | **–/71 q** | 48/72 q | –/168 q | 74/150 q |

* **kA 0.158 and 0.25: undecided.** Every gentle scheme (≥ 73 % carrier kept) was cut **before** CN's own
  onset. The 2026-11-17 question ("SDIRK_2_2 broke the mechanism" or "starved it of carrier") is still
  open there.
* **kA 0.32:** SDIRK_4_3 and TR-BDF2 delay the burst 1.2× and slow it (σ 0.11–0.13 against CN 0.15 s⁻¹).
  Generalised-α 0.5 delays it 1.9× (σ 0.08). Crouzeix 3/4 shows no onset up to 71 periods (≥ 1.8×) with
  96 % of the carrier kept, the one favourable data point (cut there). The strong schemes show no burst,
  at a cost of 90–97 % of the carrier energy.
* **Depth.** The burst belongs to kd 5.5, where 2ω₀ ≈ ω∞.
  * kd 2.9: CN is clean for 200 periods at kA ≤ 0.25; at kA 0.32 it bursts at period 107.
  * kd 1.7: CN is clean at kA ≤ 0.25.
  * kd 1.7, kA 0.32 (H/d 0.38, Ursell 5.3) fails early and differently: under CN u_max goes 1.5 → 9.5 m/s
    and Newton returns NaN at 42 s; under Crouzeix it fails at 2 s (Newton count 9 → 29 per step from the
    first step); generalised-α 0 and 0.5 were alive at 168 / 142 periods when cut.

### 3.3 Tier B: resolution (rule 38b)

| m5 kd 5.5 | n16 | n32 (tier A) | n64 |
|---|---|---|---|
| kA 0.25 | no burst under any scheme; CN clean to 200 T | CN burst at 76 | burst at 32–58 and Newton blow-up under CN, SDIRK_4_3, Crouzeix, TR-BDF2, gen-α 0.5 |
| kA 0.32 | no burst; CN clean to 200 T | CN burst at 39 | burst at 21–31, blow-up under the same five |

Only generalised-α 0 survives n64 (200 periods, 7 % of the carrier energy); SDIRK_2_2 was cut at 138
periods without an onset. **Refinement makes it worse, and the mode sits at 2 cells at every resolution:
a discretisation-stability problem, not under-resolution.** (n16 being clean is not a remedy: there the
same physical wave has half the grid wavenumbers, rule 2d.)

### 3.4 Tier C: variable bed (m6)

Model 6 behaves as m5: CN onset 75 (kA 0.25) and 39 (kA 0.32); the same γ trend as tier G.

### 3.5 Tier D and the CN Δt control (m5 kd 5.5 kA 0.32)

* **CN:** onset 80 at T/20, 39 at T/40, 40 at T/80. The T/40 result is converged in time. The coarse CN
  step, which has no dissipation, *delays* the burst (σ 0.066 against 0.147 s⁻¹): CN's phase error detunes
  the growth.
* **At T/80 every scheme bursts at 33–55 periods**: damping per unit time falls as Δt^(q−1).
* **At T/20:**
  * SDIRK_2_2 and generalised-α 0 reach 200 periods with the carrier gone, and burst anyway at 113–120;
  * SDIRK_4_3: no onset to 108 periods, carrier energy 34 %;
  * Crouzeix: none to 140, carrier 10 %;
  * TR-BDF2: onset 55, carrier 56 %;
  * generalised-α 0.5: onset 85, reaches 200 periods with 38 %.

### 3.6 Tier E: the m3 null control (no 𝓝, no CIP)

* **kd 5.5:** every scheme dies at every kA, except SDIRK_2_2 and generalised-α 0 at kA 0.158 (3 % / 11 %
  carrier energy left). At kA 0.32 everything dies within 22–37 periods (generalised-α 0: 75).
* **kd 1.7:**
  * kA 0.158: clean under every scheme;
  * kA 0.25: clean under CN, TR-BDF2, generalised-α 0 and 0.5, and SDIRK_2_2, **but SDIRK_4_3 (103 periods),
    RK3 (60) and RK4 (47) die where CN lives**;
  * kA 0.32: only SDIRK_2_2 and generalised-α 0 survive.
* ⚠ **Unexplained (rule 38e).** The multi-stage schemes (SDIRK_4_3, Crouzeix T/40, RK3, RK4) die *earlier*
  than CN on m3 at kd 1.7 kA 0.25, and reach onset earlier at kd 5.5 kA 0.158, although they add damping.

## 4. The hypotheses

* **H1 (m3, tier F: aim 4th-order damping at the 2ω₀ member): confirmed.**
  * kA 0.158: Crouzeix T/20 (9 % carrier energy kept), Crouzeix T/10 (0 %) and RK4 T/10 (22 %, **the best
    trade-off, as predicted**) survive 200 periods. Crouzeix T/40 (124 periods) and RK4 T/20 (135) die
    like CN (156). The predicted carrier levels match: 29 % / 44 % in amplitude = 8 % / 19 % in energy.
  * kA 0.25: only Crouzeix T/10 survives, with the carrier gone; RK4 T/10 dies at 112.
  * kA 0.32: everything dies (Crouzeix T/10 stalls in Newton at step 3).
  * The falsification clause (surviving kA 0.25 with > 30 % carrier) is not triggered. Not seen: a best
    Δt (prediction iii); more damping always lasted longer.
* **H2 (m5, tier G: the CIP reaches what the integrators cannot): falsified.** The CIP does not remove the
  burst, and **a larger γ makes it come sooner and grow faster** at kA 0.25:

  | γ | 0.1 | 0.3 | 0.5 | 1 |
  |---|---|---|---|---|
  | CN onset (period) | 101 | 76 | 60 | 51 |
  | CN burst growth σ (s⁻¹) | 0.069 | 0.100 | 0.120 | 0.159 |
  | m6 CN onset | — | 75 | 61 | 51 |

  * Generalised-α 0.5 and TR-BDF2 follow the same order. At kA 0.158 and 0.32, γ changes nothing.
  * The mode stays at the 2-cell scale and does not escape to the mid band.
  * The frozen check 2 (2026-11-13) had already shown σ_max rising slightly with γ above 0.3; the runs
    confirm the direction.
* **H3 (tier H: once γ has done its part, Δt sets the integrator's share): Δt is the active ingredient,
  γ is not.**
  * SDIRK_2_2 at T/20 survives 200 periods at γ 0.3, 0.5 and 1 (carrier gone); at T/80 it bursts at 46–55.
  * Crouzeix at T/20 shows no onset up to 140–145 periods (carrier 9–10 %); at T/80 it bursts at 33.
  * At a given Δt the onsets agree within 3 periods across γ. Every survival was paid with the carrier.

## 5. Conclusions

1. **No integrator, in what was run, holds the full model (m5/m6) at kd 5.5, 32 cells/λ, kA ≥ 0.25 for
   200 periods while keeping the carrier.** The schemes that last are those that remove 90–100 % of its
   energy (SDIRK_2_2, generalised-α ρ∞ = 0, any scheme at T/20). By STABILITY.md §6 that is "stabilised at a
   cost of ≥ 90 % of the carrier", outside any accuracy budget.
2. **The growing content is steady (ω ≈ 0).** Frequency-selective damping cannot reach it. The integrators
   only slow its growth (gentle schemes 1.2–1.9× later at kA 0.32) or take away its pump (the carrier).
3. **It is a grid mode.** Its wavelength is 2 cells at every resolution, it is absent at 16 cells/λ and
   faster at 64 (rule 38b). That points at the spatial discretisation, not at the time integration.
4. **The C⁰-CIP does not control it.** At the γ studied (0.1–1) it makes it come earlier at kA 0.25, and
   it does nothing at kA 0.158 and 0.32. H2 is falsified.
5. **The steady VS field persists once excited**, with or without the carrier. No term in the model damps
   steady layer-velocity shear at the grid scale.
6. **Depth matters.** The burst lives at kd 5.5 (2ω₀ ≈ ω∞). At kd 2.9 and 1.7, CN is clean for 200 periods
   up to kA 0.25. The kd 1.7, kA 0.32 failure is a separate, early, steep-shallow event.
7. **Model 3 confirms the m3 picture** (H1): the triad is held by 4th-order damping aimed at 2ω₀ only at
   kA 0.158, and only with most of the carrier gone at kA ≥ 0.25.
8. **CN "DONE" runs of m5/m6 are not clean.** They saturate at ≈ 10 % spurious grid-scale velocity.
9. **Open, because of the quota:**
   * whether the gentle schemes hold or only delay the burst at kA 0.158 and 0.25 (they were cut before
     CN's onset);
   * the Crouzeix kA 0.32 case beyond 71 periods.
10. **Unexplained:** the multi-stage schemes failing earlier than CN on m3 (§3.6).

**Direction implied (for the author's decision; STABILITY.md §6b's own rule):** "weak-form stabilisers
only if no integrator passes §6 within the accuracy budget". On this evidence none does at kd 5.5. The
remedy has to act on the steady grid-scale vertical-structure velocity, at any frequency: a term on each
layer velocity u_j at the 2Δx scale (VMS/SUPG-type, a filter F1–F4 of STATUS §4, or a VS-targeted
dissipation). Which one this solver adopts is open. (The mixed formulation with VMS is studied in the
separate solver `GridapBALFEMixed.jl`.)

## 6. Resubmission plan (drafted 2026-11-23; **not submitted**)

### 6.1 What changed in the launchers (2026-11-23)

* **One case set, one submission script.** `cases.py` writes one `cases.tsv` (275 cases) and one
  `packs.tsv` (36 packs); `submit.sh` is the only entry point. `cases_extra.tsv`, `packs_extra.tsv` and
  `submit_extra.sh` are gone. Case names are unchanged, so the existing output tree resumes as it is; the
  pack ids changed.
* **All cluster output on the project storage.**
  * `STORAGE` is exported in `~/.bashrc` (Snellius: `/gpfs/admin/projectlinks/prjs1960`).
  * The output root is `$STORAGE/GridapBALFEM_output`, and the repository's `output/` on the cluster is a
    link to it (`run/storage.sh`).
  * `balfem_env.sh` (every SLURM launcher), the integrator `submit.sh`, `linear_production/submit_all.sh`
    and `test/cluster/slurm_validation.sh` **refuse to run** if `STORAGE` is unset or `output/` is a real
    directory in home.
  * The sysimage build is exempt (`BALFEM_SKIP_STORAGE=1`).
  * The integrator study writes `$BALFEM_OUTPUT_ROOT/integrator_study/<case>/seg<NN>/`, and its SLURM
    logs go to `$BALFEM_OUTPUT_ROOT/slurm_logs/integrator_study/`.
* **New `submit.sh` commands:**
  * `status`: per case finished / failed / quota / partial / new, with totals;
  * `requeue-quota`: renames the quota-caused `FAILED` markers to `FAILED.quota`. It deletes nothing, and
    only markers whose message is the LightXML write error or `Disk quota exceeded` are touched;
  * `resume [tiers]`: submits every pack holding a case that is neither finished nor failed, and refuses
    while quota markers remain;
  * `BALFEM_CASE_FILTER=<regex>`: restricts `resume` and `run_pack.sh`, including its resubmissions, to
    the matching cases.
* Tested in a sandbox (fake project, fake `STORAGE`, stub `sbatch`/`julia`, the real markers):
  * the storage refusals work, and `setup_storage.sh` copies, keeps a backup and links;
  * `status` reproduces 92 finished / 22 failed / 154 quota / 7 partial;
  * `requeue-quota` turns them into 161 resumable cases;
  * the filter selects cases;
  * a resumed case starts `seg02` from `seg01`'s checkpoint on the storage.

### 6.2 Steps on Snellius

```bash
# 0. nothing of the old campaign still queued (pack 18 resubmitted itself on 2026-11-22: check it)
squeue -u $USER                       # cancel stale BALFEM_integ_* jobs; their pack ids are the OLD ones
# 1. storage (once)
echo 'export STORAGE=/gpfs/admin/projectlinks/prjs1960' >> ~/.bashrc && source ~/.bashrc
cd $HOME/GridapBALFEM.jl && git pull  # the launchers of 2026-11-23
bash run/setup_storage.sh             # output/ → $STORAGE/GridapBALFEM_output, link, home backup kept
#   check the copy, then free the home quota: rm -rf $HOME/GridapBALFEM.jl/output.home_backup_<date>
# 2. sysimage: rebuild if src/ changed since the campaign's build (rule 25; strict check in run_pack.sh)
# 3. make the quota-killed cases resumable, look, resume
bash run/integrator_study/submit.sh requeue-quota
bash run/integrator_study/submit.sh status | tail -6     # expect 92 finished, 22 failed, 161 partial
BALFEM_CASE_FILTER='…' bash run/integrator_study/submit.sh resume     # by priority, below
```

⚠ The campaign's runs were made from a dirty tree at `b1affeb308` (`CONFIGURATION.txt`). A resumed
segment must run the same solver. If `src/` changed since then, the checkpoint fingerprint check refuses
a different discrete problem, but a changed residual with the same fingerprint would not be caught:
resume from the same commit, or accept the change explicitly.

### 6.3 Priorities (remaining case-hours at the measured s/step; ≈ 1.8 × 10⁴ in all)

| # | group | cases | case-h | longest case | what it decides | filter |
|---|---|---|---|---|---|---|
| 1 | gentle schemes (SDIRK_4_3, Crouzeix, TR-BDF2, gen-α 0.5) at kd 5.5 n32 T/40, tiers A, G, C | 64 | 8 900 | 363 h | conclusion 9: hold or delay past CN's onset at kA 0.158 / 0.25; the Crouzeix kA 0.32 case | `'^m[56]_T160_kA[0-9]+_n32_(sdirk43\|crz34\|trbdf2\|galpha50)(_g[0-9]+)?$'` |
| 2 | CN cases cut (m6 tier C, m5 T 2.2 kA 0.32, T/80) | 7 | 190 | 82 h | the CN references of the above | `'_cn(_dt80\|_g[0-9]+)?$'` |
| 3 | strong schemes (SDIRK_2_2, gen-α 0) at kd 5.5 n32 | 28 | 1 000 | 99 h | the late steady growth with the carrier gone (low value: cost already known) | `'^m[56]_T160_kA[0-9]+_n32_(sdirk22\|galpha0)(_g[0-9]+)?$'` |
| 4 | resolution (tier B) | 12 | 1 140 | 254 h | n16 under the dissipative schemes; SDIRK_2_2 at n64 | `'_n(16\|64)_'` |
| 5 | Δt tiers D, H | 15 | 3 040 | 554 h | T/20 vs T/80 for the gentle schemes | `'_dt(20\|80)'` |
| 6 | kd 2.9 / 1.7 (tier A) | 35 | 3 670 | 311 h | confirmation where CN is already clean (lowest value) | `'^m5_T(220\|300)_'` |

(The `\|` in the filters is the Markdown table escape: type a plain `|` in the shell. Each filter was
checked against the cut cases: 64 / 7 / 28 / 12 / 16 / 36 matches; groups 5 and 6 also match one CN case
of group 2, which is harmless: a finished case is skipped.)

Group 1 alone answers the study's open question. At 8 cases per 1/8 rome node and ≈ 2 SBU per case-hour,
it costs ≈ 1.8 × 10⁴ SBU and needs ≈ 4 segments of 119 h for its longest case (363 h).
Expected volume: ≈ 200 MB per case.
