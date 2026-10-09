## 7. Standing rules

> **Reading them in v2.** These rules were earned on v1 and are carried over unchanged, because they
> are lessons about the model, the method and the workflow, not about v1's options. Where a rule
> cites `:full` read "`nl_pressure=true`"; where it cites `:native`, the projected or the mixed
> treatment, or `c3_mask`, it describes a v1 configuration (tag `v1_final_solver`,
> `HISTORY_V1.md`). Paths `output/…` in a rule's evidence are now under `output_v1/`.

These are the rules that cost something to learn. Each is stated where it is enforced; the
supporting measurement is in the linked document.

### Model and residual

1. **`R_P` is the entire frequency dispersion of the model.** Only the *boundary* part of its
   integration by parts vanishes; the volume part must be assembled. Without it the model degenerates
   to non-dispersive shallow water.
1b. **THE `C⁰` DISTRIBUTIONAL OBJECTION APPLIES TO THE *TRIAL* SPACE, NOT TO THE TEST SPACE.**
    ⚠ THE MODEL CARRIES `∂³u`; THE WEAK FORM CARRIES `∂²u`. Measured 2026-09-23 by a
    `ε(x−x₀)^k` probe: `p_nh` depends on `∂²u` and not `∂³u`; the strong-form momentum residual —
    ours and Yang & Liu's ⌊2.30⌋ alike — depends on `∂³u`. `R_P`'s integration by parts moves one
    derivative onto the test function, which is what leaves the trial field needing only a Hessian.
    Do not "correct" either count into the other; they describe different objects.
    The distributional `∂²` of a `C⁰` **trial** field is `{∂²u}` cell-wise plus a **Dirac layer on the
    skeleton** weighted by the jump `[∂ₙu]`; cell quadrature sees only the first half, so typing
    `∇∇(u)` into an integrand silently drops the layer. v1 carried the Class-III components
    `{1,2,4,5}` by frozen `L²` projections (`OPEN_ISSUES.md` §0c) or by mixed unknowns; **v2 assembles
    the layer explicitly — the broken formulation, `src/broken.jl`.**
    ⚠ **RAISING `fe_order` DOES NOT FIX THIS.** `ReferenceFE(lagrangian, …, p)` with
    `conformity=:H1` is exactly `C⁰` for **every** `p` — nodal DOFs match values across a face, never
    normal derivatives. `Q2`, `Q3`, `Q4` are all `C⁰`; the "higher order ⇒ smoother" intuition comes
    from splines (degree `p` ⇒ `C^{p−1}`), which is IGA, not Lagrangian FEM. `p ≥ 2` is *necessary*
    (a `Q1` cell-wise Hessian keeps only the mixed `∂ₓ∂ᵧ` entry and vanishes **entirely in 1-D**, so
    the term would be absent rather than approximated — the same shape as rule 2) but nowhere near
    sufficient.
    But a **test** function is a chosen, known object: `∂²v` is a computable cellwise polynomial,
    non-zero from `Q2` up and discontinuous only across faces, which is exactly what `C⁰` interior
    penalty exists to handle. ⚠ `GlobalResidual.tex` §`sec: pressure operator
    implementation` currently rules out integrating the leading-pressure block by parts "for exactly
    the same distributional reason as for the trial fields" — **that equivalence is too strong and is
    the one thing in that section that does not hold.** ⚠ Nor was this ever a missing-API problem:
    Gridap exposes `∇∇` (already used for the analytic bed Hessian, `src/nlpressure.jl:35`) and the
    full skeleton machinery (`SkeletonCellFieldPair`, jump/mean). **The task is to choose a
    formulation that is CONSISTENT with a broken `∂²` — volume term plus skeleton jump terms plus a
    penalty — not to find a library or a polynomial order that makes `∂²` classical.**
2. **`fe_order ≥ 2`.** `Q1` elements zero `R_P` and disable all non-hydrostatic physics.

2b. **ALWAYS USE TAYLOR-HOOD HORIZONTAL PAIRS: `p_eta = fe_order − 1`. NEVER EQUAL ORDER.**
    `η` enters the momentum equation **undifferentiated**, through `∇·v` after the integration by
    parts, so it plays exactly the role pressure plays in a Stokes system. Equal-order continuous
    spaces are therefore **inf-sup (LBB) deficient**, and their classic failure mode is a spurious
    **checkerboard mode at `λ ≈ 2·dx`** that refinement does not remove.

    **The evidence was already in hand and was not acted on.** The analytic MMS measures order `p`
    rather than `p+1` on equal-order spaces, and **the entire verified scope of this solver —
    all 30/30 spatial studies, five vertical bases, eight models — was measured on `Q3/Q2`.**

    ⚠ **Every nonlinear-instability run before 2026-09-05 used `Q2/Q2`** (`BALFEM_P_ETA` unset, and
    the sentinel then meant equal order). So the configuration that blew up and the configuration
    that was verified **were never the same discretisation**, and the measured instability sits at
    `λ ≈ 2–3·dx` — exactly where the equal-order spurious mode lives. That is not proof of cause,
    but it means no equal-order result can be quoted as a property of the model.

    **Enforced, not merely documented (2026-09-05):** the `p_eta = 0` sentinel in `setup_and_run`,
    `setup_and_run_distributed` and `run_flume_1d.jl` now resolves to `p_u − 1`. ⚠ **The sentinel
    was subsequently REMOVED entirely** (`0cf8cb5`): `p_eta` defaults to `p_u − 1` and
    `check_taylor_hood` REJECTS `p_eta < 1`, so equal order can only be requested explicitly and is
    then refused. `p_horizontal` was renamed `p_u` in the same commit.
    ⚠ Suite reference constants measured on equal order may shift; that is intended.

    **Corollary: a comparison across the pairing is not a one-variable comparison.** Never set an
    equal-order result beside a Taylor-Hood one and attribute the difference to physics.
3. **`B_stored = −B̃ ≤ 0`**, and the explicit `(−1)` factors in the `R_P` and slope-pressure terms
   are load-bearing. ⚠ **The invariant is NEGATIVE DEFINITENESS, not the elementwise sign.**
   `B̃ = ∫φᵢ_int φⱼ_int` is a Gram matrix, so `B_stored = −B̃ ≺ 0` for every `(M,p)` — but the
   *elementwise* `B ≤ 0` holds only because `φ_int ≥ 0` for `p = 1`. At **`p = 2` some entries are
   positive** (measured: P2LFE-2 and P2LFE-3), and that is correct, not a defect. No code depends on
   the elementwise sign — `dispersion_ratio`, `applicable_kd` and `model_celerity` all form
   `M − kd²·B_stored`, which is sign-correct for any `p`. Do not "fix" a p≥2 basis by taking `abs.(B)`.
4. **THE ASSEMBLY INVARIANT: every classification row must have exactly ONE consumer, guarded by the
   CONJUNCTION of its three activation conditions.** A guard testing only the bed condition or only
   the amplitude condition is a defect whenever the physics has a second representation in the other
   regime — which, for the `𝓐/𝓚` packages, it always does. They are **alternatives, never addends**.
   **Corollary: a flat-bed regression can never test `∇h` code.** (`MODEL.md` §7)
5. **`∂R/∂u̇` is EXACT in all eight models; nonlinear `∂R/∂u` is QUASI-NEWTON by choice**, its gap
   vanishing at order 1.11–1.16 in amplitude. **An omission is benign only if it is HIGHER ORDER IN
   AMPLITUDE** — a block whose prefactor does not scale with the solution is an `O(1)` error in the
   effective mass matrix, and Newton then converges to the fixed point of the *wrong map*. **Never
   assume; measure** with `test_jacobians_ad.jl`'s amplitude-scaling gate. Do not "complete" the
   omissions without re-measuring every nonlinear reference value.
6. **Never apply `∇` to an `Operation`-composed expression containing a test basis** — expand by
   hand via `∂_a(W⋅𝓣) = (∂_aW)⋅𝓣`.
6b. **WHEN INSERTING A HELPER ABOVE AN EXISTING FUNCTION, ANCHOR ON ITS DOCSTRING, NOT ON ITS
   `function` LINE.** Julia refuses a triple-quoted docstring immediately followed by another
   string literal (`cannot document the following expression`), so a new docstring+function
   inserted between an existing docstring and its `function` leaves the OUTER docstring
   documenting a *string* — and the whole package fails to precompile. ⚠ **This cost three
   separate edit-compile cycles in ONE session** (`NLP_REFRESH_COUNT`, `_n_multifields`,
   `mixed_aux_jacobian`): the natural anchor for a textual insert is the `function` line, and
   the docstring above it is invisible to that pattern. Anchor above the opening quotes, or
   append at end of file.
7. **NESTED CLOSURES MUST DECLARE `local`.** In Julia, a nested function assigning a name already
   local to an enclosing function **assigns the enclosing variable**. This codebase is full of long
   functions with nested helpers, so the hazard is **structural**. One instance cost two days.
   **Re-run the mechanical audit after adding any nested helper.** (`VERIFIED_SCOPE.md` §7)
   *Third instance, 2026-08-21, `wave_properties`:* two closures each called their phase speed `C`,
   which was also the enclosing function's return value, so the last γ evaluation — the **Airy** one
   — overwrote it. **Note the shape of the symptom: `|C/Ce − 1|` collapsed to ~1e-12 for every mesh
   at every `kd`, i.e. the model looked PERFECT exactly where it is worst, and `C_g` and `γ`
   reproduced their published values throughout.** A capture bug can be invisible in every channel
   but one, and the channel it corrupts can fail in the flattering direction.

2c. **EVERY RUN'S OUTPUT DIRECTORY IS NAMED BY `output_dir_name`, AND THE NAME CARRIES THE
    DISCRETISATION.** Grammar (`markdown_files/OUTPUT_NAMING_PROPOSAL.md`):

    ```
    <model>_<domain>_<wave>_<regime>_<nlp>_<bed>_<discr>_<amplitude>_<period>[_<extra>…]
    P1LFE-2_1d_bcplane_nl_full_flat_Q2Q1_A0.1_T1.6
    ```

    One generator in `src/utilities.jl`, used by all 14 drivers; the old per-driver prefixes
    (`flume_`, `small2d_`, `small_`) named the *script*, not the case. **Fields are never omitted and
    an absent field never means a default** — that convention is what let equal-order runs sit beside
    Taylor-Hood ones for months without anyone seeing it (rule 12b), and the `<discr>` token exists
    for exactly that reason.

    * `output_dir_name` **refuses impossible combinations** rather than labelling them: directional
      content on a 1-D domain (rule 12), `:linear` with `nl_pressure ≠ :none`, and any
      non-Taylor-Hood pairing (it calls `check_taylor_hood`).
    * `unique_output_dir` suffixes `_v2`, `_v3`, … and **never overwrites**. Not cosmetic: a
      re-executed batch wrote six finished runs over their own output on 2026-09-06 and destroyed
      them.
    * `BALFEM_OUTDIR` still overrides, for scratch work. The standard is the default, not a cage.

2d. **RAISING THE POLYNOMIAL ORDER IS NOT A ONE-VARIABLE CHANGE, AND IT IS NOT "MORE PHYSICS".**
    Going `Q2/Q1 → Q3/Q2` moves **two** things at once: (i) what the discretisation can *represent*
    — with `η ∈ Q2`, `∇η` is piecewise linear instead of piecewise constant, so operator content
    that was structurally absent becomes live (at `Q1`, `∂²η ≡ 0` **identically** in 1-D); and
    (ii) the **effective resolution**, since `Q3` on the same cells is ~2.1× the DOFs, and by rule
    38b more resolution *advances* onset in a grid-scale failure. ⚠ **The residual is unchanged** —
    the same `𝓝` components, the same terms — so never describe an order change as "activating more
    nonlinear effects". **Corollary, and it is the one that costs:** a run that looks *stable* at low
    order may simply be solving a **partially masked operator**, so a low-order pass is weaker
    evidence than it appears. To separate (i) from (ii), match DOF count *and* `dx` — matching only
    DOFs leaves the grid scale free (measured 2026-09-23: at ~30k DOFs, Q3/Q2 at `dx`=0.25 died
    2.3× earlier than Q2/Q1 at `dx`=0.125, which bounds but does not isolate the order effect).
    This is rule 2b's corollary in the other direction, and §5.2e is the worked instance.

### Boundaries and stability

8. **Solid-wall Dirichlet BCs must include the corner tags** — otherwise 4 corner DOFs are
   unconstrained and the run diverges exponentially.
9. **IC-release problems need `x_wall_bc=true`** — a free x-wall plus the dispersion term is a
   spurious-forcing mode that an initial perturbation excites directly.
10. **The open-boundary spurious mode is η-dominated, so the sponge MUST damp η** (`+∫ μ q η`), not
    just velocity. No value of `mu_max` absorbs it otherwise. (`WAVE_GENERATION.md` §3)
11. **Sponge strength saturates past `μ_max ≈ 5ω` — WIDTH is the lever.** And the width must cover
    the **longest** component (`kd_min` ⇒ `λ_max`), not the peak.
12. **`ny ≥ 3` is mandatory for a y-PERIODIC mesh** (Gridap `CartesianGrids.jl:39`) — and for a
    periodic mesh only. **`:wall` and `:open` accept `ny = 1`.**

    **THE DEFAULT FOR ANY 1-D HORIZONTAL CASE IS `ny = 1` WITH `y_wall_bc=:wall`**, not three cells
    with periodicity. The solver is structurally 2-D, so a 1-D problem is a narrow flume; one cell
    across with solid walls is both the cheapest and the more correct way to pose it:
    * for a normal-incidence wave the exact solution has `𝖴y ≡ 0`, and the wall condition `𝖴y = 0`
      is **exactly consistent** with it — it approximates nothing. `:periodic` merely *permits*
      `𝖴y ≡ 0` while admitting a family of y-periodic modes a true 1-D model does not have, whose
      shortest member has wavelength `Ly` and can therefore sit inside the physical band;
    * measured on the 240-cell flume: **7215 free DOFs at `ny=1`/`:wall` against 20202 at
      `ny=3`/`:periodic`** — 2.8×, and a direct LU costs more than linearly in DOFs. The wall pins
      the bottom and top `𝖴y` node layers (2886 constrained = 2 levels × 481 x-nodes × `Nσ`);
    * set `Ly = dx` so the single cell stays isotropic.

    Use `:periodic` only where the case genuinely carries oblique or short-crested content, which a
    solid wall would reflect. `examples/local_1d/run_flume_1d.jl` defaults to `ny=1`, `Ly=0.25`,
    `BALFEM_YBC=wall`, and all six `run/local/run_1d_*.sh` use it.

    ⚠ **AND A 1-D CASE IS ALWAYS NORMAL-INCIDENCE — never oblique, never short-crested.** A 1-D
    domain carries one propagation direction; oblique content has a transverse wavenumber
    `k_y = k sin θ`, and a flume one cell across cannot represent it. The request is not refused by
    the mathematics — it is **silently aliased** onto a normal-incidence wave at the wrong
    wavenumber, with the transverse component dropped, producing a wrong answer that runs to
    completion and looks plausible. This is why `build_airy_state` is called **without**
    `directional=true` in the 1-D driver, and why that driver now **errors** if any of
    `BALFEM_WAVE_DIR` (non-zero), `BALFEM_NTHETA`, `BALFEM_SPREAD_STD`, `BALFEM_THETA_MAX` or
    `BALFEM_DIRECTIONAL` is set — those variables were previously *ignored*, which is the worse
    failure. Directional content belongs in the 2-D driver
    (`run_directional_sea_small.jl`: `y_wall_bc=:open` plus lateral sponges).

    **Corollary: the `:wall` default above is not a compromise.** It is exact precisely *because*
    1-D cases are normal-incidence — the two rules support each other.
13. ~~**`A_wave ≤ 0.001 m`** for stable long fully-nonlinear integrations.~~ ✅ **LIFTED 2026-09-06 —
    THIS WAS THE INSTABILITY WEARING A DISGUISE.** The cap was never a property of the model; it was
    the equal-order pairing's growth rate being slow enough at tiny amplitude to finish a run
    (rule 12b). On **Taylor-Hood** the 1-D flume completes 50 wave periods at `A_wave = 0.10 m`
    (`κa ≈ 0.16`, `kd = 5.5`) with the mean `η` flat to the fourth decimal — **100× the old cap**.
    ⚠ This does not license unlimited amplitude: physical limits (Miche, breaking) still apply, and
    the model has no breaking closure. It removes a *numerical* restriction, not a physical one.

12b. ✅ **THE "NONLINEAR INSTABILITY" WAS AN EQUAL-ORDER ARTEFACT — RESOLVED 2026-09-06.**
    *(This is the complete record; `NONLINEAR_INSTABILITY.md` was folded in here and deleted.)*

    **The phenomenon, as it appeared for a year.** Fully nonlinear runs grew an unbounded
    free-surface mode while the *same case* linear was flat for 80 s at 4× the amplitude. Growth
    appeared in the interior, not at a boundary. **Refining `dx` made it worse** — blow-up at
    t≈24 s for `dx=0.25`, t≈10 s for `dx=0.125`, barely at all at `dx=0.50` — which is why it read
    as a discretisation-stability problem rather than under-resolution (rule 38b). The measured
    growth spectrum peaked at **`λ ≈ 2–3·dx`**, gaining 10³–10⁴× while the carrier gained ~6×.

    **The cause: equal-order `Q2/Q2` horizontal spaces.** `η` reaches the test function only through
    `∇·v`, so it plays the pressure role of a Stokes system and the pairing is subject to the
    inf-sup (LBB) condition. Equal order is deficient and admits a spurious checkerboard at
    `λ ≈ 2·dx` — exactly the measured peak. On **Taylor-Hood** the mode does not exist:

    | `dx` | equal order `Q2/Q2` | Taylor-Hood `Q2/Q1` (peak / settled t≥30) |
    |---|---|---|
    | 0.50 | suppressed | ✅ 80 s — 0.1169 / 0.11055 |
    | 0.25 | **died t=26.4, η→5.06** | ✅ 80 s — 0.1115 / 0.10593 |
    | 0.125 | **died t≈10** | ✅ 80 s — 0.1112 / **0.10481** |

    **The refinement signature INVERTED** — settled amplitude now *decreases* monotonically with
    `dx`, converging on the delivered 0.102. Newton holds at ~4 iterations/step against 30–57 as the
    equal-order runs came apart, and `x_at_max` migrates with the crest instead of pinning at one
    station. Confirmed on `Q3/Q2` as well, and on **three integrators** — `SDIRK_2_2` (dissipative, rule 15),
    `SDIRK_3_3`, and explicit `EXRK_RungeKutta_4_4` (near dissipation-free, which is what retires the
    `dt`-masking objection; RK4 tracked SDIRK to within 1–5 % with the same envelope, the offset
    being exactly the SDIRK_2_2 damping).

    **Ten hypotheses were refuted before the pairing was questioned** — worth keeping, because each
    is a real property of the solver:

    | hypothesis | verdict, on evidence |
    |---|---|
    | sponge reflection; domain length | onset identical at `Lx` = 40 and 90 |
    | boundary generation artefact | the interior source fails too, at matched *delivered* amplitude |
    | inflow relaxation zone | widening it delays onset, never prevents it |
    | CFL / time step | 4× `dt` invisible at fixed `dx` |
    | incomplete hand Jacobian | exact AD agrees to 4 dp *through the divergence* (rule 17b) |
    | Benjamin–Feir | a physical rate cannot depend on `dx`; measured σ ∝ A⁴, BF is A² |
    | under-resolution | refinement makes it worse |
    | quadrature aliasing | degrees 6/10/14 agree to 4 dp against an exact control |
    | Galerkin advection lacks an energy sink | the operator's production is *exactly* the continuity defect, and that defect is 5 orders too small — the correction moved the solution <5e-5 while the mode grew 10³–10⁴× |

    ⚠ **THE ROOT CAUSE OF THE ROOT CAUSE: the configuration that blew up was never the configuration
    that was verified.** The MMS campaign ran `Q3/Q2`; every instability run ran `Q2/Q2`, because the
    launchers left `p_eta` unset and the sentinel then meant equal order. The two were compared for
    months as though they were the same solver. **Before diagnosing a numerical failure, diff the
    failing configuration against the verified one — parameter by parameter, including the ones
    nobody thought to set.** This is now enforced by `check_taylor_hood` (rule 2b), so it cannot
    recur.

    ⚠ **Everything measured about the mode is void, not superseded** — the `σ ∝ A⁴` rates, the
    "no amplitude threshold" claim, the tier ordering, the `dt`-masking numbers, the `A_wave ≤ 0.001`
    cap (old rule 13). They characterise a discretisation no longer in use. **Discard them; do not
    re-explain them.**

    **Method lessons that cost time and generalise:**
    * **Copy the reference environment verbatim and vary one variable explicitly.** `env -i` plus a
      partial list once flipped generation from `:bc` to the interior source (2.8× the amplitude) and
      produced a confident result pointing the wrong way.
    * **Put the null-treatment control in the same batch**, never against a remembered baseline.
    * **Verify a knob is LIVE — and then that it is BIG ENOUGH TO MATTER.** A dead parameter gives
      three identical curves; a live but negligible one gives a clean negative that means nothing.
      Measuring the treatment against the effect size costs seconds.
    * **A reversal too large to be the effect under test is a bug signal, not a finding.**
    * **A "threshold" may just be a run that ended too early.** Stability claims need a duration.
    * **Separate the phenomenon from its symptoms.** Newton stalling at `‖r‖≈1.2` was real and
      reproducible but was the quasi-Newton Jacobian failing to track an already-diverging solution.
    * **Identify an instability by what GROWS, never by what is LARGEST** — ranking by amplitude
      found the carrier's harmonics and read as evidence *against* a grid mode; ranking by gain found
      it. An extremum on the edge of the search window is a finding about the window.
      (`postprocessing/examples/growth_spectrum.jl` implements both guards.)

12c. **"STABLE" IS MEANINGLESS WITHOUT `dx`, `dt` AND DURATION.** A dissipative integrator's numerical
    dissipation grows with `dt` and can *mask* a growing mode, so **a run that completes may simply
    be one whose dissipation exceeded the growth rate.** This is rule 15's trap in a new guise, and
    the rule stands on its own merits — but note how it was DISCHARGED for 12b, because that is the
    template: the same case was run on `SDIRK_2_2` (2nd order, dissipative), `SDIRK_3_3` (3rd order) and
    **explicit `EXRK_RungeKutta_4_4`** (essentially non-dissipative). All three stayed flat, and RK4
    tracked SDIRK_2_2 to within 1–5 % with the same envelope shape — the offset being exactly the
    SDIRK_2_2 damping. **Three integrators of different order and stability character agreeing is what
    retires a dissipation-masking objection; a `dt` ladder on one scheme is not.**
    ⚠ Gridap DOES provide explicit tableaux (`Gridap.jl/src/ODEs/ODESolvers/TableausEX.jl`, 18 of
    them incl. classical RK4). They cost a mass solve per stage here, since `∂R/∂u̇` carries the `B`
    dispersion term and is not the identity — but they run, and they are the right tool for exactly
    this question.
14. **The `c_g` transit trap.** At `kd = 5.5`, `c_g = 1.25 m/s` — filling a 45 m flume takes 22.5
    periods. Budget `t_settle ≈ (x_sponge − x_source)/c_g + 3T` before reading any steady state. It
    caught three separate measurements.

### Solver and execution

14b. **THE INTERIOR SOURCE DOES NOT DELIVER `A_wave`, AND THE FACTOR IS GEOMETRY-DEPENDENT.**
    Measured `η/A`: **3.32 on the 50 m small domain, 2.13 on the 60 m 1-D flume**; Dirichlet
    generation delivers ~1.0–1.05. The ratio is a *linear* property of the source calibration —
    identical to three digits between a `linear/:none` run at `A=0.1` and a `nonlinear/:full` run at
    `A=0.001` — so it is not a nonlinear artefact and cannot be tuned away. **Every `:inner_res`
    result must be rescaled before it is read**, and "3.3×" must not be quoted as a constant.

14c. **A CRASH IS NOT EVIDENCE UNTIL ITS CONTROL RUNS.** The small-domain suite is a factorial for
    a reason: a `nonlinear/:full` run diverging at `A=0.1` looks amplitude-driven, but the **linear**
    run at the same *delivered* 0.332 m — 76 % of the Miche limit — completes. Amplitude alone is
    survivable; it takes amplitude *and* nonlinearity. Conversely the directional case died at
    **5.7 %** of Miche, lower than three cases that completed, so it is a different failure
    (generation-region, velocity-led) wearing the same symptom. Pair every failure with the run
    that differs in exactly one axis.

15. **THE DEFAULT INTEGRATOR `SDIRK_2_2` IS STRONGLY DISSIPATIVE — AND IT IS NOT THE SCHEME ITS NAME
    SUGGESTS (corrected 2026-09-26).** Gridap's `:SDIRK_2_2` is `DIRK22(1,0,1)`
    (`Gridap.jl/src/ODEs/ODESolvers/TableausDIM.jl`): `A = [1 0; −1 1]`, `b = (½, ½)`,
    `R(z) = (1 − z − z²/2)/(1 − z)²`. It is A-stable but **not** L-stable (`R(∞) = −½`), and it removes
    `(¾)y⁴/(1+y²)² ≈ 0.75 y⁴` of a neutral mode's energy per step, `y = ωΔt` — **100× the textbook
    L-stable `γ = 1 − 1/√2` SDIRK2** (0.0074 y⁴), which every earlier note here had assumed. At
    `Δt` = 0.04: carrier 0.011 s⁻¹ (exactly the carrier decay of every periodic run), `ω∞` 0.15 s⁻¹,
    fastest `:full` modes 2–10 s⁻¹.
    * Any test measuring a non-dissipative property must pin `solver_type=:theta`. **Do not remove those
      pins, and never fix such a failure by moving a threshold.** Recognise it by: amplitude damped
      while **phase is correct**, error **growing with frequency**, and **refining the mesh does not
      help**. (`TEST_SUITE.md` §4)
    * ⛔ **NO STABILITY CLAIM MAY REST ON `SDIRK_2_2` ALONE.** It masked the `:full` interior instability:
      Q2/Q1 at 32 cells/λ ran 100 periods bounded under SDIRK and **diverged at ≈ 90 s under
      Crank–Nicolson**; the 64-cell SDIRK run showed no growth at all and **diverged at 31 s under
      Crank–Nicolson**; and SDIRK hid a `:native` instability altogether (Q2/Q1 64 cells/λ diverges at
      146 s under CN). Pair every stability run with a
      Crank–Nicolson repeat (`BALFEM_SOLVER=theta`), which has `|R(iy)| ≡ 1`.
    * Before quoting any tableau's damping, read the tableau in `TableausDIM.jl` — never infer it from
      the name.
16. **Distributed linear solve = `NewtonSolver(GMRESSolver(Pr=Jacobi))`.** A direct LU does not scale
    to partitioned matrices at cluster size.
17. **`krylov_m` (basis size, memory) and `ls_maxiter` (iteration budget, time) are different
    bounds**, and `restart=true` is load-bearing. Symptom of getting this wrong: **`gmres=` pinned at
    exactly the same number every step.** (`ARCHITECTURE.md` §5)
17b. **⛔ THE INCOMPLETE `jacobian_u` DOES NOT AFFECT STABILITY. PROVEN, NOT ARGUED (2026-09-16).**
    A converged Newton step is a property of the RESIDUAL, not the Jacobian: the Jacobian sets the
    *path* to the root and the *cost*, the root is defined by `r(u)=0`. So an incomplete but
    convergent Jacobian **cannot change the answer** — only the iteration count.

    **THE DECISIVE MEASUREMENT.** Every cell of the `:full` flat `dx`×`dt` factorial (§5.2c) was run
    TWICE — once with the hand Jacobian, which omits the `{1,2,4,5}` blocks, and once with
    `BALFEM_USE_AD=1`, the exact AD Jacobian of the same residual:

    | `nx` | `dt` | hand onset | AD onset | Newton, last 5 steps |
    |---|---|---|---|---|
    | 240 | 0.04 | **12.60 s** | **12.60 s** | hand 10,12,15,19,51 · AD 4,4,4,4,6 |
    | 480 | 0.04 | **5.00 s**  | **5.00 s**  | hand 6,6,8,12,58 · AD 4,4,4,6,10 |
    | 480 | 0.02 | **8.20 s**  | **8.40 s**  | hand 6,8,8,10,16 · AD 4,4,4,6,12 |

    **A CRASHING SIMULATION CRASHES IN THE SAME PLACE WITH EITHER JACOBIAN.** Onset is identical in
    two cells and one output interval (0.2 s) apart in the third; `u_max` at failure agrees to three
    significant figures (2.5046 vs 2.4971; 1.6764 vs 1.6702); `η` agrees to five digits and `r0` to
    three at every sample up to failure. **The only difference is the iteration count** — the hand
    Jacobian degrades 6 → 12 → 19 → cap while AD holds ~4 to the last step, so the quasi-Newton gap
    costs TIME and nothing else.

    ⚠ **CONSEQUENCES, both directions.**
    * **Never diagnose an instability, a divergence or a wrong answer by completing the Jacobian.**
      Completing `jacobian_u` would buy iteration count, not stability. §5.7 item 5 previously said
      the opposite and was corrected on this evidence.
    * **Never read Newton stalling as the primary fault.** A stall at 50 iterations is the
      quasi-Newton Jacobian failing to track a solution that is *already* diverging underneath it —
      the symptom, not the cause.

18. **`norm(PVector, Inf)` is broken** (PartitionedArrays 0.3.5) — reduce over `own_values`.
19. **Distributed ICs use `interpolate_everywhere`**, never `FEFunction(U, zeros(…))`.
20. **Keep `ConsecutiveMultiFieldStyle`** — `BlockMultiFieldStyle` breaks Jacobi's `diag`.
21. **Launch MPI with `~/.julia/bin/mpiexecjl`**; the system `mpiexec` fails on a PMIx mismatch.
    `MPI_Finalize` prints a benign OFI error and exits 143 on this machine.
22. **Keep horizontal cells near-isotropic, or pay for it in the linear solve** — a 4:1-celled mesh
    needs ~760 GMRES iterations against ~480 for an isotropic one.
23. **Use MPI when the problem is big enough to give each rank a meaningful share, and a direct LU
    otherwise.** Spend spare cores on more *cases*, not on decomposing one small case further.
24. **`nl_tol` is a step function and DOES change the answer** (integer Newton counts); `ls_rtol`
    does not, anywhere in 1e-9…1e-5. Do not conflate them. (`CONFIGURATION.md` §4)
25. **A sysimage is only valid for the versions it was built against** — build it in the environment
    you run in, and rebuild after **any** `src/*.jl` edit. Staleness is *detected*, not prevented.
26. **Julia buffers stdout when redirected to a file**; use `flush` or poll.
27. **Revise does not hot-swap signature changes** — restart before trusting a number after any
    signature or struct change.

### Long campaigns (earned 2026-08-21…30, the vertical-basis campaign)

41. **A LONG-LIVED Julia+Gridap PROCESS DEGRADES TO USELESSNESS.** Measured: workers grew
    1.5 → 2.5 → 3.9 GB over ~14 h and, left for days, fell to **7-14 % CPU** — GC-bound, not
    compute-bound. One spent **6.2 days on an `nx=8` level that takes seconds when fresh**, which is
    why a campaign sat at 25/36 for a week while "still running". **Bound worker lifetime and RSS and
    have a supervisor restart them** (`run_vbasis_shard.jl` + `supervise.sh`); resume logic makes a
    restart cost one JIT. ⚠ The cap is checked BETWEEN studies, so a single long study can still get
    there. This also **unseats the "flat memory ⇒ H3 not supported" argument** in `OPEN_ISSUES.md` §1,
    which was measured over 400 steps on one tiny case.
42. **`Distributed`/`pmap` was NOT usable for this workload; independent processes were.** Three
    `pmap` launches each completed ONE study in >4 h while the identical `run_mms_case` calls ran at
    full speed in a plain process. Sharded ordinary processes writing their own CSVs fixed it — and
    fixed observability too: under `pmap` the workers' prints are relayed through the master, whose
    stdout buffer is never flushed, so there is **no way to tell slow from hung**.
43. **CONTIGUOUS slicing of a COST-SORTED queue is the worst partition for makespan** (measured 12×
    imbalance: 1.7 h vs 20.8 h). Use LPT to minimise makespan for a fixed set, **SJF when the
    deliverable is coverage** — and beware that LPT schedules the cheap high-value studies LAST.
44. **In-flight work is invisible unless you make it visible.** A case is only "done" when it writes
    rows, so restarting slots re-run what others are mid-way through. Use **PID-keyed claim files in
    a SHARED directory** — per-output-dir claims let batches collide (one study ran on three slots).
    A claim honoured only while its PID is alive self-heals; a plain lock file would not.
45. **`RETRY_ERRORS` must be switched OFF once a failure is established.** A deterministic failure
    re-queued forever burns slots on work already understood — and these reproduce **bit-identically**
    (`‖r‖ = 0.27484014031572` twice), so a retry is not a new sample.

46. ✅ **REPRODUCING A PUBLISHED OPTIMISATION: THE THREE THINGS THAT WENT WRONG (2026-09-08).**
    `src/vopt.jl` reproduced Yang & Liu Table 1 to only 0.0085 at `M=2`. Three defects, none of them
    numerical, and the order in which they were found is itself the lesson.

    **(a) The median was over the wrong population — this was the big one.** Each of the five terms
    is *already* `kd`-integrated, so it is one number per candidate mesh and its median can only be
    over *meshes*. The code took it over a Kronecker low-discrepancy sweep of the **width simplex**;
    the paper takes it over **the scan it optimises on** ("c₂ … from 0 to 1 with an increment of
    0.001"). The medians are the **trade-off weights**, so a different population is a different
    optimisation problem, not a different approximation of the same one. Fixing it: 0.7365 → 0.7316
    against a published 0.728.
    ⚠ The LaTeX had a *third* reading — the half-mass point in `kd`, `∫₀^m E_X W = ½∫_Ω E_X W`.
    It matches at `M=2` and is **2.4–2.6× worse** at `M=3,4`. Do not reinstate it.

    **(b) `E_cg` was squared.** (3.10) takes `C_g` to the **first** power; `E_c` alone is squared.
    The asymmetry is the paper's own and is easy to flatten by copy-paste — it had been flattened in
    both the code and the LaTeX.

    **(c) `E_shoal` had a spurious `abs`.** (3.10) carries **no** `|·|` there (verified at 400 dpi;
    the next line shows `|u_m−u_e|` with unmistakable bars). ⚠ **The signed form is correct AND
    self-consistent, for a reason that is easy to "fix" wrongly:** `exp(I)−1` is single-signed across
    the design population (`γ_m > γ_e` throughout), so its median is negative too and **the sign
    cancels in the quotient** — the normalised term is positive and identical to putting a modulus
    around the whole term. This holds **only** if the median is over the `exp(I)−1` population.
    Normalising by `median(exp(I)) ≈ +0.81` breaks the cancellation, inverts the term's role, and
    drives the optimum to the scan edge (`c₂→0`). *That was a live defect I introduced while fixing
    (c), and it produced a confident, entirely wrong conclusion that the paper's equation was
    pathological.*

    **`M ≥ 3` DOES NOT REPRODUCE THE PUBLISHED NODES, AND THAT IS SETTLED.** Measured, at fixed band:
    `M=3` over `Δc` = 1e-2 → 1e-3 (4 851 → 498 501 designs) drifts **3e-4**, the last two levels
    identical to 4 dp; `M=4` over 3e-2 → 5e-3 (4 960 → 1 293 699 designs) is **unchanged to 4 dp**.
    The gaps (1.9e-2, 2.6e-2) are 60–100× the drift. The reference specifies its population only for
    `M=2`, so **the `M≥3` interfaces are not recoverable from what was published**. Ours are a
    defensible member of that family from a stated, reproducible population.

    **κ WAS A FIT READ OFF ITS OWN FIT — deleted, do not reintroduce.** `VOPT_KAPPA = 2.17` closed
    `Ω` for `p ≥ 2` via `Ω·Δσ_top/p = κ`. But `Ω` had first been *tuned* to reproduce Table 1
    (`Ω* = 7.60, 29.0, 100.0`) and the near-constancy of the product noticed afterwards. It is also
    35 % from the analytic `Δσ_top ≈ 2.94/kd_max` the LaTeX derives from the same boundary-layer
    argument. **A constant obtained by fitting cannot then validate the thing it was fitted to.**
    ⚠ The "`M=3,4` are best fitted by exactly 1.25× the published band" coincidence was an artefact
    of the same circularity — discard it.

    **Cost:** the scan is `C(1/Δc − 1, M−1)` designs, so the paper's `Δc = 1e-3` is exact at `M=2`
    (999) and impossible at `M=4` (1.66e8, ~8 days). Use `Δc` = 1e-3 / 5e-3 / 3e-2 for `M` = 2/3/4:
    converged optima to 1e-4 in minutes. **One notebook run cost 461 min for numbers identical to
    four decimals.**

47. ⚠ **THIS WORKSTATION LOSES RUN EVIDENCE IN THREE WAYS, ALL LEARNED 2026-09-23.**
    * **`/tmp` IS EMPTIED ON EVERY BOOT** (`tmpfiles.d` carries `D /tmp`). Long solver runs launched
      from an agent session write stdout into a scratchpad under `/tmp`, so **a reboot destroys
      every run log** — and the log is the *only* place holding the config banner, the per-step
      Newton trace and **the reason a run stopped**. `diagnostics.csv` records none of that; it
      simply ends. **Redirect batch stdout into `output/…` directly, never to a scratchpad.**
    * **A GNOME/Wayland session cannot restart its shell in place.** `org.gnome.Shell@wayland.service`
      *is* the compositor; restarting it is a logout to GDM, not a reload. ⚠ **And with `Linger=no`
      and a single session, ending that session stops `user@<uid>.service` and kills `app.slice` —
      i.e. every detached solver dies with the GUI**, `nohup` notwithstanding, because they live in
      `vte-spawn-*.scope` units under that slice. **`loginctl enable-linger <user>` decouples them**
      and is the precondition for touching the desktop at all (enabled here 2026-09-23).
    * **SUSPEND IS SAFE; REBOOT IS NOT.** Suspend preserves both the running solvers and `/tmp`.
      A reboot loses both. When a GUI freeze forces the choice, **suspend**, and persist logs first
      either way.
    ⚠ **A traceback ending in `paraview_collection`/`createpvd` is NOT a VTK failure** — that is the
    enclosing `do`-block in `run_time_loop`. The real cause is higher up; use
    `grep -m1 '^ERROR' run.log`. This misled the first reading of four crashed arms.

### Testing and measurement

28. **"The suite passes" is not "the model is verified."** Most of the suite is
    **self-consistency**: writing `R = R_true + E`, the error appears on both sides and cancels
    identically, so such a test passes for **any** residual. Only the analytic MMS — whose forcing
    never touches `problem.jl` — can detect a self-consistently wrong residual. (`VERIFIED_SCOPE.md` §0)
29. **AD is an oracle for the JACOBIAN, never for the RESIDUAL.** "AD converges where hand fails,
    therefore the residual is correct" is an invalid inference — AD differentiates the *same*
    assembled residual.
30. **A BOUNDS CHECK ON THE RIGHT CONFIGURATION IS NOT A VALUE CHECK.** A gate that only asserts
    nothing exploded will pass while the quantity it computes moves 58 %.
31. **For a guard that is a CONJUNCTION (`nonlinear ∧ ∇h ≠ 0`), the suite needs a case satisfying the
    conjunction AND asserting a value.** Conjuncts satisfied separately, or together but only
    bounded, both look green.
32. **A refinement study measures the rate of whichever error DOMINATES — verify isolation in BOTH
    directions before interpreting any slope.** A saturated slope and a wrong coefficient produce the
    *same* observable. **Guards come in pairs.**
33. **Read the pairwise rate SEQUENCE, not the fitted slope**, and check error *magnitude* before
    trusting a fine-level high-order rate. **Three independent instances in one campaign** where the
    fit would have given the wrong conclusion and the sequence gave the right one (fit 3.515/3.729/
    3.603 against a true 3.94/3.88). The diagnostic that generalises: **a rising sequence with the
    error still dropping is PRE-ASYMPTOTIC; a flat sequence with a stalled error is not.** And a
    *fitted* slope through points that were never on an asymptotic curve (P2LFE-2 M4: `e_u` starts
    LARGER than the solution amplitude) is **meaningless, not low** — never report it as a rate.
34. **THE RESOLUTION PRINCIPLE: a test validates a term only if it can RESOLVE that term's
    contribution.** Always ask: *if this term were wrong, would this test have noticed?*
35. **A batch runner must take its verdict from GATE OUTPUT, never from exit codes.** A clean exit
    code is not evidence a test ran.
36. **A mirrored guard that collapses two per-field measurements into one pass/fail can hide the
    very asymmetry it exists to detect.**
37. **Two tests sharing a bathymetry function and a physics tier can still be different problems.**
    A reference-vs-reference comparison is meaningful only when *every* discretisation parameter
    matches.
38. **When an error message names a library type, that is where the bug SURFACED, not where it
    lives.** Print the type of *every* input to the failing expression before searching that library.
38b. **REFINEMENT THAT MAKES THINGS WORSE IS A DISCRETISATION-STABILITY PROBLEM, NEVER
    UNDER-RESOLUTION.** Under-resolution improves with `h`; a grid-scale mode does not. This single
    check separated the two whole classes of explanation faster than anything else in the
    investigation, and it is the first thing to measure whenever a nonlinear run diverges.

38c. **COPY THE REFERENCE ENVIRONMENT VERBATIM AND VARY ONE VARIABLE EXPLICITLY.** `env -i` plus a
    partial list of variables silently changed wave generation from `:bc` to the interior source,
    which delivers **2.8×** the requested amplitude — so the treatment ran at triple the intended
    forcing and produced a confident result pointing the wrong way. **And put the null-treatment
    control in the SAME batch**: without it, a difference from a remembered baseline cannot be
    attributed to the treatment.

38d. **VERIFY A NEW KNOB ACTUALLY DOES SOMETHING BEFORE SPENDING THE RUN ON IT.** A dead parameter
    yields identical curves — a clean, confident, entirely wrong negative result. A seconds-long
    unit check (here, integrating `x⁸` at each quadrature degree) buys the whole experiment.

38f. **IDENTIFY AN INSTABILITY BY WHAT GROWS, NEVER BY WHAT IS LARGEST.** Ranking a spectrum by
    *amplitude* found the carrier's own harmonics at `λ≈1.7 m` — a wavelength independent of `dx`,
    which reads as evidence *against* a grid mode. Ranking the same data by *gain* found the real
    mode at `λ≈2–3·dx`. An unstable mode is small for most of its life; that is why it goes unnoticed
    until it dominates. **Corollary: when a reported extremum sits at the edge of the search window,
    the window is the finding** — a `0.95·k_Nyq` cut put one case's "peak" exactly on the boundary.

38e. **A REVERSAL TOO LARGE TO BE THE EFFECT UNDER TEST IS A BUG SIGNAL, NOT A FINDING.** Extra
    quadrature turning a 0.11 plateau into 1.60 is not a plausible quadrature effect; that
    implausibility is what prompted the check that found the real cause.

38g. **A NULL RESULT IS NOT A FINDING UNTIL THE EFFECT SIZE IS MEASURED AGAINST THE EFFECT BEING
    EXPLAINED.** Rule 38d says verify the knob is LIVE; this is the other half. The 2026-09-12
    quadrature probe returned `e_η` columns that were **bitwise identical** between treatment and
    control at every level — which reads as "dead knob, result void" (I called it that) and equally
    as "refuted". Both readings were wrong: the knob was live, and its effect simply fell below the
    CSV's 7-significant-digit print precision. **Identical printed output distinguishes nothing**;
    only assembling the quantity under test at both settings does. When it was assembled, the crime
    was real (3.5e-08 in the residual) and six orders too small to explain a 3.1e-07 rate deficit —
    which IS the finding, and is the same shape as the skew-advection refutation in rule 12b.
    **Corollary: print precision is part of the experiment.** A comparison whose resolution is set by
    `%.6e` cannot detect anything below 1e-7 relative, however carefully the runs were controlled.

38h. **THE LAUNCHER IS PART OF THE EXPERIMENT. Three ways a run silently became a different run,
    all on 2026-09-15, all costing hours.**
    * **An override appended AFTER the `balfem_local_run` line is a NO-OP.** A diagnostic built that
      way ran the production configuration instead, and reported it under the diagnostic's name.
      The tell was rule 38d's: η, Newton counts, `r0` and the final residual **identical to every
      printed digit** against the control. An exact Jacobian cannot reproduce a quasi-Newton
      iteration count to the digit — *check that the knob moved something before believing the run*.
    * **NEVER EDIT A SHELL SCRIPT WHILE BASH IS EXECUTING IT.** Bash reads a script incrementally by
      byte offset, so rewriting the file under a running interpreter shifts the content beneath it:
      the fixed launcher's relocated run-line was re-executed by the *old* process, silently starting
      a second solver that competed for RAM for an hour before anyone noticed. Write a NEW file.
    * ⚠ **`exit status 1` DOES NOT MEAN THE RUN FAILED.** `examples/local_1d/run_flume_1d.jl:350`
      ends with `@printf(..., tag, ...)` and **`tag` is defined nowhere** — a leftover from the
      `output_dir_name` rename (rule 2c). Every *successful* run therefore throws `UndefVarError`
      after writing every solve result, VTK file and diagnostics row. A batch runner that reads exit
      codes will score a completed 100-period run as a failure (rule 35). One-line fix, not yet made.

38i. **TO SEPARATE AN INTERIOR INSTABILITY FROM A BOUNDARY FEED, REMOVE THE BOUNDARY.** A failure that
    pins at the inflow cannot be attributed by any amount of flume refinement, because every flume run
    carries the inflow, the relaxation zone and the sponge. The **closed x-periodic box** (one model
    wavelength, no inflow, relaxation, sponge or source, deterministic 1e-8 seed, sub-cell sampling,
    band energies relative to the carrier) has a known null — nothing grows — and it measures the growth
    over full wave cycles (a Floquet-type rate, which frozen eigenvalues cannot give). Run it with its
    `:native` twin **and** under Crank–Nicolson (rule 15). It settled in three days what two weeks of
    flume factorials could not (§5.2f).
39. **Test a diagnosis against a case it cannot explain, rather than looking harder where it points.**
39b. **A CONSTRUCTION THAT *CANNOT* COMMIT THE SUSPECTED ERROR IS THE STRONGEST DISCRIMINATOR
    AVAILABLE — AND ITS NULL RESULT IS A REFUTATION, NOT A DISAPPOINTMENT.** Rule 39 says test a
    diagnosis where it cannot reach; this is its constructive form. For a year the `dx` refinement
    signature of `:full` was attributed to **recovering `∂²` of a `C⁰` field** (the skeleton Dirac
    layer that cell quadrature never sees). The mixed formulation was built to remove the *lag*, but
    it also **structurally cannot form that second derivative at all** — `𝖦 ≈ ∇𝖲` is a genuine
    unknown defined by an integrated-by-parts weak equation. **The `dx` signature survived it
    unchanged.** That single observation refutes an explanation that no amount of further refinement
    studies on the projected path could have touched, because every such run shares the error being
    blamed. ⚠ **Read the consequence honestly and immediately:** the refutation also invalidates the
    *remedy* queued behind it — `C⁰` interior penalty exists to make a broken `∂²` consistent, and a
    signature measured where there is no `∂²` to break cannot justify it. **When a structural fix
    fails to move a signature, retire the hypothesis AND everything scheduled on its authority**, in
    the same edit, before the queued work is launched on a dead rationale.
40. **Diagnostics interpretation:** never read `growth` without `x_at_max`; `dmp/int` means "≫1 is
    trouble", never "<1 is fine"; mass drift is an invariant **only** in a closed, unforced basin.
    (`ARCHITECTURE.md` §6)

---

