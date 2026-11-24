# Linear production suite (Snellius) — results check, 2026-11-02

Runs: `run/linear_production/` (models 1–2, P1LFE-2, Q2/Q1, CN, dt = 0.02, GMRES + Jacobi), launched at
commit 148371c + uncommitted edits. Data: `output/snellius/linear_production/` (19 GB, downloaded
2026-11-02). Analysis: vertex values of η from the per-rank `.vtu` pieces merged onto the node grid
(Δx = 0.25 m; Δy = 0.5 m for the 50 × 20 m cases).

## Completion

| case | status | wall | s/step | GMRES/solve | Newton/step |
|---|---|---|---|---|---|
| bcplane flat / bar | complete, 80 s | 45.0 / 47.7 h | 40 / 43 | 515–558 | 2.95 |
| irregular flat / bar | complete, 120 s | 70.5 / 73.9 h | 42 / 44 | 550–597 | 2.98 |
| ring | complete, 60 s | 22.2 h | 27 | ≈ 400 | 1.00 |
| **directional flat / bar** | **killed by the time limit at t = 69.2 / 66.8 s of 100 s** | 120 h | 22 → 125–137 | 344 → 900–980 | 1 → 2.9 |

The directional cost estimate (55–65 h) was wrong by about 2×. Two reasons:
* GMRES + Jacobi needs about 930 iterations per solve once the domain is full. That mesh is
  isotropic at 0.25 m, while the other cases use 0.25 × 0.5 m, and the Jacobi conditioning goes as
  (H·p/h_min)².
* Newton takes 3 iterations per step.

No solve hit the 1000 cap. Rerun options: `precond = :schwarz` (OwnBlockLUSolver, 15–19 iterations
on model 5), or restart from a checkpoint (none was written: the runs predate the feature).

## Physics checks: all pass

* **Stability.** All bounded. Energy stays steady once the domain is full (bcplane 28.44–28.46 over
  the last 40 s; ring 39.5). The response is periodic to 2e-4 (bcplane) and 2e-3 (ring) between
  consecutive 10 s windows.
* **Linearity and 1-D invariance.** |a(2ω)| ≤ 6e-5 against a(ω) ≈ 0.1. The y-variation of |a| in
  the y-periodic plane-wave cases is 3.5e-11.
* **Dispersion (bcplane).** The measured k is within 0.45 % of Airy theory on the flat bed (kd = 3.5)
  and within 0.15 % over the bar crest (kd = 2.1).
* **Ring.** √r·|a| = 0.32–0.33 for r = 4–12 m (1/√r spreading). The pattern is exactly 4-fold and
  mirror symmetric (to 1e-5). The amplitude at the source is 0.41 m for A = 0.1 (interior source,
  rule 14b).

## Bar vs flat (bar: d = 3.5 − 0.75[tanh((x−20)/2) − tanh((x−32)/2)], crest 2.0 m at x ≈ 22–30)

At T = 2 s the waves are intermediate-to-deep over the whole bar (kd 3.5 → 2.1), so linear theory
predicts a small effect: shoaling Ks = 0.960, a 3 % shorter wavelength, and negligible reflection.
The runs reproduce this:

| quantity | bcplane (T = 2) | irregular (Hs) | directional (Hs, y 10–30) | linear theory |
|---|---|---|---|---|
| over the crest, bar/flat | 0.965 | 0.964 | 0.964 | 0.960 |
| downstream (34–38 m) | 1.002 | 0.998–0.999 | 1.000 | 1 |
| phase gained across the bar | 0.23 rad | — | — | ≈ 0.3 rad |

Per frequency (irregular, incident/reflected fitted on 6–17 m and over the crest):
* f ≥ 0.45 Hz (most of the energy): crest transmission matches Ks(f) to ≤ 1 %, and downstream
  agrees to ≤ 0.5 %.
* f = 0.35–0.42 Hz (T 2.4–2.9 s, a ≤ 0.017 m): crest 0.91–0.93 against Ks 0.95–0.97, downstream
  0.94–0.97. These fit windows are shorter than λ/2, so the numbers are uncertain.
* The bar adds no measurable reflection. The counter-propagating fraction is identical with and
  without the bar at every frequency.

## Unusual behaviour (none is a solver defect)

1. **The outflow sponge (12 m, μ = 40) reflects 5 % of the amplitude at T = 2 s.** This gives a
   ±5 % standing-wave ripple (period L/2) over the whole interior, identical in the flat and bar
   cases. The reflection grows with period: 0.5 % at 1.5 s, 9 % at 2.2 s, 18 % at 2.9 s. 12 m is
   about one wavelength of the long tail (rule 11: width is the lever).
2. **The "irregular" sea is exactly periodic with period 60 s** (η(t + 30) = −η(t) to 4e-3; all the
   variance sits at odd multiples of 1/60 Hz). This comes from the uniform sampling of 15
   frequencies, Δf = 1/30 Hz, at mid-bins. The 120 s irregular run therefore holds two copies of
   one 60 s record, and statistics beyond 60 s add no information. Fix: more frequencies, or random
   or jittered frequency sampling.
3. **Ring: about 10 % of the amplitude returns from the 7 m sponges** (≈ 1.1 L, walls behind them).
   This gives an azimuthal modulation of up to ±25–30 % at fixed r (axes vs diagonals). Its 4-fold
   symmetry and its change with r show it is interference with boundary reflections, not grid
   anisotropy.
4. **Directional: Hs varies by 55–77 % across y.** The variation is already present in the imposed
   relaxation zone (x = 1–5 m): same-frequency components from different directions form
   time-invariant interference patterns. Shadow zones lie within 4–8 m of the open lateral
   boundaries, so the usable width is y ≈ 6–32 m.
5. **Equation-check warnings** (ring 59 ×, directional 3 ×; θ-scheme ‖R‖∞ ≈ 6e-8 > check_tol 1e-8).
   They occur on steps where Newton stops after one iteration with ‖r‖₂ just under nl_tol = 1e-5.
   This is a tolerance mismatch between check_tol and nl_tol, not an error in the solution.
6. **Linear problem, 3 Newton iterations per step** in the BC-generated cases: each GMRES solve
   reduces the residual by only ls_rtol = 1e-5 from r0 ≈ 15. This is a cost matter (3 assemblies + 3
   Krylov solves per step), not an accuracy one.
