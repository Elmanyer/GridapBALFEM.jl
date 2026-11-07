# Sponge-layer reflection sweep, stage A: partial results (2026-11-07)

*The design, the question and the pass criterion are in [`CH8_CAMPAIGN_TODO.md`](CH8_CAMPAIGN_TODO.md)
§7 and §9.2. This file holds what was measured. The sweep is **paused, not finished**: 26 of 84 cases
are done. How to resume is in §4.*

## 1. Set-up (as run)

* **Launcher and code.**
  * Launcher: `run/local/run_1d_sponge_sweep.sh`.
  * Workers: `run/local/sponge_worker.jl`, one JIT per worker and up to 12 cases each.
  * Case table: `output/1D_sponge_sweep/_logs/_cases_A.tsv`, written by `run/local/sponge_cases.py`.
  * Analysis: `postprocessing/examples/sponge_reflection.jl`.
* **Physics:** the Snellius linear suite in 1-D.
  * P1LFE-2, Q2/Q1, CN, dt = 0.02, d = 3.5 m, flat bed, linear regime.
  * Dirichlet (BC) generation of a regular wave, A = 0.1 m, `:model` profile, 6 m inflow relaxation zone.
  * ny = 1 with walls, dx = 0.25 m.
* **Geometry:** relaxation 6 m | fixed 32 m interior | quadratic outflow sponge of width W, ending at a
  wall. Kr is fitted over x ∈ [7, 37] m in every case.
* **Analysis:** complex amplitude at ω over the last 10 periods; least-squares fit of the incident and
  reflected waves with k scanned ±5 % around Airy. The outputs are:
  * Kr = |a_R|/|a_I|;
  * the relative rms fit residual;
  * the generated amplitude |a_I|/A;
  * |a| at the wall relative to |a_I| (`a_end`).
* **Code state:** commit 2557ab8 plus the uncommitted diff. A source snapshot is in
  `_logs/source_2026110{6,7}_*`. Models 1–2 are untouched by every change since the Snellius runs.

## 2. Results (26 cases; `output/1D_sponge_sweep/kr_A.csv`, `kr_A_2026-11-07.log`, heat maps `kr_A_T*.png`)

Kr in % of the incident amplitude. "floor" = W = 48 m, μ_max = 4: what generation and relaxation
reflect on their own.

| T [s] (λ [m], ω [1/s]) | floor | W = 18: μ = 4 / 8 / 16 / 25 / 40 | W = 24: μ = 4 / 8 / 16 / 25 / 40 |
|---|---|---|---|
| 1.5 (3.5, 4.19) | 0.003 (residual-limited) | 0.020 / 0.039 / 0.079 / 0.12 / 0.21 | 0.011 / 0.022 / 0.044 / 0.070 / 0.11 |
| 2.0 (6.2, 3.14) | 0.011 | not run | 0.046 / 0.093 / 0.19 / 0.30 / 0.32 |
| 2.5 (9.6, 2.51) | 0.031 | not run | 0.12 / 0.26 / 0.31 / **1.05** / **3.03** |
| 2.9 (12.4, 2.17) | 0.054 | not run | 0.23 / 0.26 / not run / not run / not run |

Quality of the measurement:
* **Fit residual.**
  * At T = 1.5 s the residual is a constant ≈ 2.4e-4 in every case, the floor included. So a Kr below
    about 5e-4 is not resolved at that period, and the script marks it `?`.
  * At T ≥ 2 s the residual is 2.5e-5 to 9e-4, below Kr/2 in every case.
* **Discrete wavenumber:** k_fit/k_Airy − 1 = +0.89 / +0.46 / +0.23 / +0.12 % at T = 1.5 / 2.0 / 2.5 /
  2.9 s. This is the P1LFE-2 dispersion error, which grows with kd (kd = 6.3 at T = 1.5 s). It is the same
  in every case at a given T.
* **Generated amplitude:** |a_I|/A = 1.001–1.018. BC generation delivers A to within 2 %.
* **Absorption:** `a_end` is 1e-9 to 1e-4 except where the sponge is too weak for the wavelength. At
  T = 2.9 s, W = 24 m, μ = 4 it is **4.4 %**, and at T = 2.5 s, μ = 4 it is 0.6 %.

**Not yet run:**
* the round-1 **replica** (W = 12, μ = 40), which is the check that this 1-D set-up reproduces the
  Snellius Kr (5 % at 2.0 s, 18 % at 2.9 s);
* every W = 6 and W = 12 case;
* W = 18 at T ≥ 2 s;
* T = 2.9 s, W = 24, μ ≥ 16.

## 3. Observations and deductions

1. **The floor is far below 1 %** (0.003–0.054 %, rising with T). Generation and the 6 m relaxation
   zone are not a limitation at this level.
2. **In the weak-damping regime, Kr is linear in μ_max and falls as 1/W².**
   * At T = 1.5 s, Kr/μ_max is constant to ±5 % from μ = 4 to 40 at both widths.
   * The W = 18 / W = 24 ratio is 1.79–1.85, against (24/18)² = 1.78.
3. **Across periods the data follow Kr ≈ C · μ_max / (ω k² W²) with C ≈ 0.21.**
   * Fitted C: 0.212 (T = 1.5, W = 18, μ = 4), 0.215 (T = 1.5, W = 24, μ = 40), 0.211 (T = 2.0, W = 24,
     μ = 4).
   * Predicted Kr at μ = 4, W = 24: 0.13 % at T = 2.5 s (measured 0.12) and 0.26 % at 2.9 s (measured 0.23).
   * A first-order (WKB) estimate of the reflection from the jump in μ'' = 2μ_max/W² at the sponge front
     gives the same scaling, with C = 1/8 to 1/4 depending on whether the damping enters one equation or
     both. The model damps both η and u (rule 10). The *reason* for the scaling is an estimate, not a
     derivation checked in detail, but the scaling itself is measured.
   * **Reading:** the quadratic sponge reflects at its front, where μ'' jumps, as the §7 hypothesis said.
     The wall is not the source, except in the under-damped case of item 5.
4. **The linear law breaks down once μ_max/ω exceeds about 3–6.**
   * At T = 2.5 s, Kr is 0.31 % at μ = 16 (μ/ω = 6.4, below the law's 0.54 %), then 1.05 % at μ = 25
     (μ/ω = 10) and 3.0 % at μ = 40 (μ/ω = 16, 2.2× the law).
   * At T = 2.0 s, μ = 40 (μ/ω = 13) gives 0.32 %, below the law's 0.46 %.
   * So beyond the saturation μ_max ≈ 5ω (rule 11), extra strength buys no absorption, and from
     μ/ω ≈ 10 it buys steeply more reflection. This is the regime of the Snellius sponge: W = 12, μ = 40
     gives μ/ω = 13–18 over T = 2–2.9 s.
5. **The lower bound on μ_max is absorption.** At T = 2.9 s, W = 24 m (1.9 λ), μ = 4 (μ/ω = 1.8), 4.4 %
   of the amplitude still reaches the wall. Kr (0.23 %) is still small, but the margin is gone. With
   μ = 8 (μ/ω = 3.7), `a_end` is 0.28 %.
6. **The best measured setting is W = 24 m, μ_max = 8.** Kr = 0.022 / 0.093 / 0.26 / 0.26 % at
   T = 1.5 / 2.0 / 2.5 / 2.9 s, so it **passes the § 7 criterion (Kr ≤ 1 %) at every T, with a 4× margin**.
   The Snellius sponge (12 m, μ = 40) reflects 5–18 %.
7. **Predictions for the cases not yet run (from item 3; to be confirmed, not used as results):**
   * W = 18, μ = 8 at T = 2.9 s: about 0.9 %, marginal;
   * W = 12, μ = 8 at T = 2.9 s: about 2 %, failing.
   * The smallest W that passes at T = 2.9 s is therefore probably 18–24 m (1.5–2 λ), in line with the
     §7 expectation W ≈ 1.5–2 λ and μ_max ≈ 2–5 ω.

## 4. State of the sweep and how to resume

* **Stopped 2026-11-07 ~11:00** for the machine shut-down (author).
  * 26 cases are finished (`flume_1d [...] done:` in their logs).
  * No case is partial in the live directories.
  * Interrupted attempts from the earlier stops were *moved*, not deleted, to
    `_logs/_interrupted_20261106_{1830,1900}/` and `_logs/_interrupted_20261107_0240/`. They are
    debris: the cases they hold were re-run, or are still in the queue.
* **Queue lock.** The 58 unstarted cases hold a claim directory whose `pid` file reads
  `RESERVED-STOP …`. That is what stopped the workers. **To resume:**
  ```bash
  cd output/1D_sponge_sweep/_logs/_claims && grep -l RESERVED-STOP */pid | xargs -n1 dirname | xargs rm -r
  cd - && MAXP=3 MIN_FREE_GB=4 setsid nohup run/local/run_1d_sponge_sweep.sh >> output/1D_sponge_sweep/_logs/_batch_A.out 2>&1 < /dev/null &
  ```
  The launcher skips every case with a log, so the 26 finished cases are not re-run.
* **Use MAXP = 3** (counting every local 1-D Julia process) on this workstation. Above 3 busy cores the
  CPU overheats and the firmware throttles it to 1.1 GHz for minutes at a time. See
  [`CH8_CAMPAIGN_TODO.md`](CH8_CAMPAIGN_TODO.md) §10.
* **Cost at full speed** (3 busy cores, 4.6 GHz): 0.3–2 s per step depending on nx. A worker's 12 cases
  take about 6 h. The remaining 58 cases are the cheaper end of the queue (the queue is ordered longest
  first): about 12–15 h on 2–3 workers.
* **Re-analyse at any time** (it skips unfinished cases):
  `julia --project=postprocessing postprocessing/examples/sponge_reflection.jl output/1D_sponge_sweep A`.
* **Check order once the sweep is complete:**
  1. the replica reproduces Snellius (rule 38c), or this set-up does not represent the 2-D runs;
  2. floor ≪ 1 % (done);
  3. residuals;
  4. the pass criterion, at the smallest W.
  Then stage B (relaxation width 12 m at the chosen pair), and then the §7 "Apply" step before the
  Snellius round 2 (§8).
