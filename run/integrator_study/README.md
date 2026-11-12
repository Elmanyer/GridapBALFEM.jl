# `run/integrator_study/` — the time-integrator stability study (Snellius `rome`)

**Question** (STABILITY.md §6b, author 2026-11-11): Crank–Nicolson has shown the discrete nonlinear
operator unstable in every box run. How much of that growth does numerical dissipation from the time
integration hold, and at what cost to the carrier? The study comes before any new weak-form term.
It leans on the regime that challenges the full nonlinear solver: steep waves (kA up to 0.32).

## The set-up (every case)

* **Box:** the closed x-periodic box (`examples/local_1d/run_periodic_1d.jl`), one model wavelength,
  flat bed (m6: cosine bed, max slope 0.1), seed 1e-8 on every harmonic.
* **Discretisation:** P1LFE-2 (Nσ = 3), Q3/Q2, `ncell` cells per wavelength.
* **CIP penalty γ** (γ_u = γ_η, a column of `cases.tsv`, passed as `BALFEM_CIP_GU/GE`): 0.3 with 𝓝
  (m5/m6, mandatory) unless the case varies it (tiers C, G, H); **0 for m3/m4**, whose weak form has no
  non-conforming term (author, 2026-11-12).
* **Time:** 200 periods, dt = T/40 unless the case varies it (tiers D, F, H). At T/40, ωΔt = 0.157 at
  every T, so each integrator removes the same share of the carrier per period whatever the period.
* **Steepness:** `BALFEM_KA` fixes k₀A with k₀ the model's own wavenumber. The amplitude follows:

| T [s] | kd (model) | 2ω/ω∞ | A at kA = 0.158 / 0.25 / 0.32 [m] | H/λ vs Miche | H/d at kA 0.32 |
|---|---|---|---|---|---|
| 1.6 | 5.54 | 1.00 | 0.100 / 0.158 / 0.202 | 35 / 56 / 72 % | 0.12 |
| 2.2 | 2.94 | 0.73 | 0.188 / 0.298 / 0.381 | 36 / 56 / 72 % | 0.22 |
| 3.0 | 1.68 | 0.53 | 0.329 / 0.521 / 0.667 | 38 / 60 / 77 % | 0.38 (Ursell 5.3) |

  All within the P1LFE-2 ranges (kd ≤ 10.84 for C, ≤ 6.42 for Cg). The kd 5.5 cases sit at the
  triad threshold (2ω = ω∞); kd 3.0 and 1.7 are below it.

## Integrators (`cases.py` INTEGRATORS)

Five in the matrix (author, 2026-11-11: TR-BDF2, Crouzeix 3/4 and generalised-α dropped, explicit
RK3 added; 2026-11-12: the explicit pair runs on m3 only, see below, and SDIRK_Crouzeix_3_4 joins on
m5/m6 as RK4's implicit stand-in):

| tag | scheme | order / stages | loss per step: carrier (ωΔt 0.157) / 2ω₀ / ωΔt = 1 | role |
|---|---|---|---|---|
| `cn` | θ-method θ = 0.5 | 2 / 1 | 0 / 0 / 0 | **diagnosis pair** (rule 15): what grows with no damping |
| `sdirk22` | `SDIRK_2_2` (DIRK22(1,0,1)) | 2 / 2, R(∞) = −½ | 2.2e-4 / 3.0e-3 / 9.9e-2 | strong damping at every frequency; the carrier pays (58 % per 100 periods) |
| `sdirk43` | `SDIRK_4_3` | 3 / 4, **L-stable** | 1.3e-5 / 2.0e-4 / 1.7e-2 | the only L-stable scheme: removes stiff modes outright, carrier cost 5 % per 100 periods |
| `crz34` | `SDIRK_Crouzeix_3_4` (= Gridap's `SDIRK_Norsett_3_4`) | 4 / 3, R(∞) = 0.63 | 6.2e-6 / 3.1e-4 / 4.7e-2 | **m5/m6 only**: 4th order like RK4, damping rising ≈ y⁶ (2ω₀ 50× the carrier, against 16× for the 3rd-order schemes); 2.4 % of the carrier per 100 periods; CIP modes damped only ≈ 40 % per step |
| `rk4` | `EXRK_RungeKutta_4_4` | 4 / 4, explicit | 1.0e-7 / 6.6e-6 / 6.1e-3 | Yang & Liu's scheme: nearly non-dissipative; unstable if a discrete mode has ωΔt > 2.83. **m3 only** |
| `rk3` | `EXRK_SSP_3_3` | 3 / 3, explicit | 2.5e-5 / 3.9e-4 / 2.8e-2 | the dissipative explicit scheme: SDIRK_2_2's (ωΔt)⁴ shape at 1/9 of its level (10 % of the carrier per 100 periods); unstable if a discrete mode has ωΔt > √3. **m3 only** |

The two explicit schemes are a pair. Their stability limits differ (√3 against 2.83), so a mode
between the two limits blows up RK3 and leaves RK4 standing. Each explicit stage still solves a
system here, because the mass operator ∂R/∂u̇ carries the dispersion term B; the stages are not
cheap. To run a scheme in the pilot only, put it in `PILOT_EXTRA` in `cases.py`.

**Explicit schemes run on m3 only (author, 2026-11-12).** The mandatory C⁰-CIP of m5/m6 is a
damping term: it adds modes with large negative *real* eigenvalues. Frozen-state spectra
(`stability_eig.jl`, P1LFE-2 Q3/Q2, A = 0.1, T = 1.6; m3 with γ = 0.3, the CIP on its own):

| cells/λ | fastest damping rate | stable Δt RK3 (\|z\| < 2.51) | stable Δt RK4 (\|z\| < 2.79) | vs T/40 = 0.04 s |
|---|---|---|---|---|
| 8 | 3.9e3 s⁻¹ | 6.5e-4 s | 7.2e-4 s | 60× smaller |
| 16 | 4.7e3 s⁻¹ | 5.3e-4 s | 5.9e-4 s | 75× smaller |
| 32 | 9.5e3 s⁻¹ | 2.6e-4 s | 2.9e-4 s | 150× smaller |

Without the CIP (m3) the fastest discrete mode is a wave at ≈ 9 s⁻¹ (imaginary), so T/40 sits well
inside both limits. With it, a 200-period n32 case needs ≈ 1.6 million explicit steps: years of wall
time. The smoke test confirms the limit: at Δt = 0.04, m5 grows ×10⁶ in energy in one RK4 step. So
RK3/RK4 are in tier E (m3) only; on the full model the comparison is between implicit schemes.
An IMEX split (CIP implicit, the rest explicit) would lift the limit; it is not implemented
(the drivers refuse IMEX tableaus).

Not in the matrix:
* θ > 0.5: first order, 2.4e-3 per step at θ = 0.6, so the carrier is gone within a few hundred steps;
* the schemes that amplify oscillations: explicit first/second order, the DIRK Radau and Lobatto IIIC variants;
* Gridap's defective `SDIRK_3_3` and `SDIRK_LobattoIIIC_2_2` (rule 15).

## The matrix (`cases.tsv`, 170 cases in 23 packs)

| tier | what | cases |
|---|---|---|
| **P** pilot | m5 n32 T 1.6 kA 0.25, the 4 implicit integrators + CN at γ = 3, 100 steps: s/step per scheme on the node, and whether a stiff penalty slows Newton | 5 |
| **A** core | m5 n32 × T {1.6, 2.2, 3.0} × kA {0.158, 0.25, 0.32} × 4 implicit integrators | 36 |
| **B** resolution | m5 T 1.6 × kA {0.25, 0.32} × n {16, 64} × 4 (rule 38b: does refinement help or hurt?) | 16 |
| **C** variable bed | m6 n32 T 1.6 × kA {0.25, 0.32} × 4; kA 0.25 also at γ {1, 3} × 4 | 16 |
| **D** time step | m5 n32 T 1.6 kA 0.32 × 4 implicit integrators × dt {T/20, T/80} | 8 |
| **E** null control | m3 (𝓝 off, no CIP) n32 × T {1.6, 3.0} × kA {0.158, 0.25, 0.32} × cn, sdirk22, sdirk43, rk4, rk3 (RK3/RK4 at T/40) | 30 |
| **F** selective damping (H1) | m3 n32 T 1.6 × kA {0.158, 0.25, 0.32} × {Crouzeix at T/40, T/20, T/10; RK4 at T/20, T/10}: aim the 4th-order damping at the 2ω₀ triad member | 15 |
| **G** γ map (H2) | m5 n32 T 1.6 × kA {0.158, 0.25, 0.32} × γ {0.1, 1, 3} (`GAMMAS_G`, provisional) × 4 implicit integrators | 36 |
| **H** γ × dt (H3) | m5 n32 T 1.6 kA 0.32 × γ {1, 3} × {SDIRK_2_2, Crouzeix} × dt {T/20, T/80} | 8 |

Tiers P–D, G, H (m5/m6) run `cn`, `sdirk22`, `sdirk43`, `crz34`; tier E (m3) runs `cn`, `sdirk22`,
`sdirk43`, `rk4`, `rk3`. **The hypotheses H1–H3, why each tier tests them and what each should show:
`markdown_files/CH8_CAMPAIGN_TODO.md`, "Integrator study: hypotheses".** `GAMMAS_G` is to be reset
after pre-cluster check 2 (`run/local/run_precluster_checks.sh`) to bracket its predicted boundary.

Regenerate after any change: `python3 run/integrator_study/cases.py` (it writes both tables and
prints the per-tier relative cost). Names: `<model>_T<100T>_kA<1000kA>_n<ncell>_<integrator>[_dt<spp>][_g<100γ>]`
(`_dt<spp>` when dt ≠ T/40, `_g<100γ>` when γ is not the model's default).

## How it runs

* **Packs.** One array task runs one pack of up to 8 sequential cases side by side on a **rome 1/8
  node** (each tier is cut into the fewest packs, as equal in size as possible: 23 packs of 5–8) (16 cores, `--mem=28G`, ≈ 3.5 GiB per case; a 𝓝 box case peaked at 2.8–3.0 GB locally),
  one thread each, Julia loaded from the **sysimage** (`-J`, strict freshness check). Cases of one
  tier and model are cost-sorted and grouped, so a pack's cases end at similar times.
* **Segments.** `run_pack.sh` stops each case with `timeout` at `BALFEM_SEG_HOURS` (119 h). The next
  segment `<case>/seg<NN+1>` resumes from the newest checkpoint (every period: `ckpt_every` = spp
  steps in `cases.tsv`, so a cut loses at most one period, ≈ 1.5–3 h for the costliest cases). A
  pack with a case cut this way **resubmits itself** (up to `BALFEM_MAX_SEG` = 10 segments).
  Estimated segments per case (local timings, performance mode): 80 cases need 1, 31 need 2,
  26 need 3, 29 need 4, 3 need 5, and m5 kA 0.32 SDIRK_4_3 at T/80 about 7. Every segment queues
  again, so submit the long tiers (A, G, C, D, H) early. Markers in
  `<case>/`:
  * `FINISHED`: the driver printed `done:` (completed, or stopped by the divergence guard; both are
    results);
  * `FAILED`: an error. A Newton `NaN` here is a **blow-up result**, not a fault; read the message.
  The stub test of this logic (finish, cut + resume, fail, resubmission limit) passed on 2026-11-11.
* **Output.** `$BALFEM_STUDY_OUT/<case>/seg<NN>/` (default `$BALFEM_PROJ/output/integrator_study`):
  `CONFIGURATION.txt`, `run_manifest.toml`, `diagnostics.csv`, VTK and diagnostics every
  `save_every` = max(1, spp ÷ 8) steps, i.e. 8–10 samples per period whatever dt is (the
  temporal-spectrum analysis needs them: at 2 samples per period 2ω₀ would alias to 0). About 160 MB
  per case, about 30 GB in all: check the quota. `run.log`, `exit_status`.

## Launch

```bash
# on Snellius, after pulling v2-solver
cd $HOME/GridapBALFEM.jl/compile && sbatch compile_snellius.sh      # sysimage: performance mode, Nσ = 3
bash $HOME/GridapBALFEM.jl/run/integrator_study/submit.sh P         # pilot: one pack, ~100 steps per case
bash $HOME/GridapBALFEM.jl/run/integrator_study/submit.sh status    # per-case state
# read the pilot's s/step per integrator (run.log), re-estimate the cost below, then:
bash $HOME/GridapBALFEM.jl/run/integrator_study/submit.sh E F     # cheap: m3 null control + the H1 ladder
bash $HOME/GridapBALFEM.jl/run/integrator_study/submit.sh A         # the core
# reset GAMMAS_G in cases.py from pre-cluster check 2, regenerate, commit, pull, then:
bash $HOME/GridapBALFEM.jl/run/integrator_study/submit.sh G C
bash $HOME/GridapBALFEM.jl/run/integrator_study/submit.sh B D H
```

## Cost (to be replaced by the pilot's numbers)

The only measurements are local: m5 n32 under CN at ≈ 58 s/step on a loaded laptop, in Gridap debug
mode. The cluster sysimage runs in performance mode (≈ −20 %), the only speed-up kept: the
depth-field memo was reverted on 2026-11-12, its benchmark showed no gain. m3 runs at ≈ 5 s/step. An SDIRK_2_2 step costs ≈ 2.5 CN steps.
In "m3 CN case" units (8000 steps, ≈ 11 h), `cases.py` weights the tiers A 981, B 436, C 480, D 273,
E 88, F 26, G 981, H 275, about 3540 units in all.

* **Order of magnitude:** ≈ 3.9 × 10⁴ case-hours. Each case is billed 2 cores (1/8 node per 8
  cases), so ≈ 8 × 10⁴ SBU before any speed-up, ≈ 6.4 × 10⁴ with performance mode alone.
* **Tiers A and G:** each about 28 % of that; the m3 tiers E and F together 3 %.
* **SDIRK_4_3 and Crouzeix 3/4** (4 and 3 stages) are about two thirds of the bill.
* **Upper bound:** every case that diverges stops early, and at kA 0.32 most CN cases are expected
  to.

**Decide the scope after the pilot.**

## Analysis

Per case, stitch the segments in time order. Each has its own `diagnostics.csv` and VTK series, and
`load_periodic` reads one segment directory. Then:
* `energy_history(r; modes = :groups)`: the carrier against the disturbance;
* the growing modes' rates against the CN pair;
* the carrier loss reported in `CONFIGURATION.txt`.

Judge by STABILITY.md §6 (σ ≤ 0 over the run, carrier loss within the author's budget).
