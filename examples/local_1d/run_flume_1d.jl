# ==============================================================
#  run_flume_1d.jl — PARAMETRIC quasi-1D flume (local + cluster)
#
#  ONE script for every 1-D-horizontal case. The physics, the wave-generation
#  mechanism and the geometry are all selected by environment variables (set in
#  the launcher); nothing is hard-coded per case. See
#  markdown_files/TEST_SUITE.md (the local-validation plan itself is gone) §5.
#
#  WHAT "1-D horizontal" MEANS HERE. The solver is structurally 2-D
#  (CartesianDiscreteModel on a rectangle, Ex/Ey throughout the residual), so a
#  1-D horizontal problem is posed as a NARROW FLUME: a long, thin domain with
#  y-periodic (or y-wall) lateral boundaries and the minimum number of cells
#  across. This needs no solver change and is the repository's established way
#  of doing it (test_bc_spectrum.jl does the same).
#
#  It is NOT a 1-D model, and the difference matters when interpreting results:
#    * the 𝖴y DOFs exist and are solved. For a y-invariant state they are ~0,
#      but they still enter the Jacobian, the GMRES spectrum and the cost;
#    * y_wall_bc=:periodic admits y-periodic modes that a true 1-D model has
#      none of — good for isolating streamwise behaviour from lateral-wall
#      artefacts, but it is an extra mode family, not fewer;
#    * ny CANNOT go below 3 with :periodic — Gridap asserts "a minimum of 3
#      elements is required in any periodic direction"
#      (Gridap/src/Geometry/CartesianGrids.jl:39). ny=1 and ny=2 are rejected
#      at mesh construction. :wall accepts ny=2, but ny=3 is used for both so
#      the wall/periodic comparison is like-for-like.
#
#  CONFIG (env; defaults = the linear flat-bed reference case)
#    BALFEM_WAVE_GEN     inner | bc | sea            inner
#                        inner = interior Gaussian line source (plane wave)
#                        bc    = Dirichlet boundary generation, regular wave
#                        sea   = Dirichlet boundary generation, WaveSpec JONSWAP
#    BALFEM_REGIME       linear | nonlinear          linear
#    BALFEM_NL_PRESSURE  none | native | full        none
#    BALFEM_FLAT_BED     1 flat | 0 submerged bar    1
#    BALFEM_MPI          0 sequential | 1 MPI        1   (sequential keeps GAUGES;
#                                                       the distributed driver
#                                                       has none — see the plan §4.2)
#    BALFEM_PX           MPI ranks in x (BALFEM_MPI=1) 12  (px·1 must equal mpiexec -n)
#    BALFEM_LX/LY        domain [m]                  60 / 3
#    BALFEM_NX/NY        cells                       240 / 3  (dx=0.25 = 16 cells/λ)
#    BALFEM_D            still-water depth [m]       3.5
#    BALFEM_TWAVE        period [s]                  1.6   (⇒ kd=5.5, λ=4.0 m)
#    BALFEM_AWAVE        amplitude [m]               0.001
#    BALFEM_PERIODS      duration in wave periods    16 inner / 26 bc,sea (transit-based)
#    BALFEM_RELAX        inflow relaxation zone      1 for bc/sea, 0 for inner
#    BALFEM_XWM          interior source position    sponge_wL + 6 m (must clear the sponge)
#    BALFEM_HBAR/XBAR/WBAR   bar shape (FLAT_BED=0)  1.0 / 30 / 5
#                        height / centre / HALF-width [m] (bar spans xbar±wbar)
#    BALFEM_SBAR         bar shoulder length [m]     wbar/3
#                        small ⇒ square-shouldered (box) bar; large ⇒ trapezoid
#  plus every knob of examples/distributed/_dist_common.jl (solver, tolerances,
#  sea state, output).
#
#  RUN
#    local (12-rank)   : run/local/run_1d_<case>.sh
#    by hand           : BALFEM_MPI=1 BALFEM_PX=12 ~/.julia/bin/mpiexecjl --project=. -n 12 \
#                          julia --project=. examples/local_1d/run_flume_1d.jl
#    sequential+gauges : BALFEM_MPI=0 julia --project=. examples/local_1d/run_flume_1d.jl
#    cluster           : sbatch run/dist_small/run_1d_<case>.sh
# ==============================================================

include(joinpath(@__DIR__, "..", "distributed", "_dist_common.jl"))

# base-case physics (launchers override; get! ⇒ the banner and the solver can
# never disagree, because both read the same resolved ENV)
get!(ENV, "BALFEM_REGIME",      "linear")
get!(ENV, "BALFEM_NL_PRESSURE", "none")
get!(ENV, "BALFEM_FLAT_BED",    "1")

wave_gen_kind = lowercase(genv("BALFEM_WAVE_GEN", "inner"))
wave_gen_kind in ("inner", "bc", "sea") ||
    error("BALFEM_WAVE_GEN must be inner, bc or sea (got $wave_gen_kind)")
use_mpi = genv_b("BALFEM_MPI", 1)

# ---- geometry / numerics -------------------------------------------------
M       = genv_i("BALFEM_M", 2)
#  Vertical BASIS ORDER. The model is named P{p_vert}LFE-{M}: `Pp` is the
#  vertical Lagrange order and `M` the number of vertical elements, so the run
#  says which member of the BALFE-M family it actually exercises. Default p=1
#  reproduces the piecewise-linear models of Yang & Liu.
p_vert  = genv_i("BALFEM_P_VERT", 1)
model_name = "P$(p_vert)LFE-$(M)"
Lx, Ly  = genv_f("BALFEM_LX", 60.0), genv_f("BALFEM_LY", 0.25)
nx, ny  = genv_i("BALFEM_NX", 240), genv_i("BALFEM_NY", 1)
feord   = genv_i("BALFEM_FE_ORDER", 2)
#  ⚠ BALFEM_P_ETA defaults to BALFEM_FE_ORDER-1 (Taylor-Hood). There is no longer a
#  sentinel: check_taylor_hood REJECTS p_eta<1 and any non-Taylor-Hood pairing. eta enters momentum undifferentiated (via div(v) after IBP), so it
#  plays the pressure role of a Stokes system: equal-order continuous spaces are inf-sup
#  deficient, the analytic MMS measures order p rather than p+1 on them, and the ENTIRE
#  verified scope of this solver was measured on Q3/Q2. Every nonlinear-instability run
#  before this date used equal order -- i.e. a discretisation the campaign never verified.
#  Equal order is still reachable with an EXPLICIT BALFEM_P_ETA = BALFEM_FE_ORDER, but it
#  is not a supported production configuration. CLAUDE.md rule 2b.
p_eta   = genv_i("BALFEM_P_ETA", feord - 1)
d       = genv_f("BALFEM_D", 3.5)
Twave   = genv_f("BALFEM_TWAVE", 1.6)
Awave   = genv_f("BALFEM_AWAVE", 0.001)
dt      = genv_f("BALFEM_DT", 0.04)
#  DEFAULT DURATION IS SET FROM THE TRANSIT TIME, and it differs by generation
#  type. At kd=5.5 the group velocity is only c_g≈1.25 m/s, so filling the flume
#  from the source to the far sponge takes:
#      interior source at x≈18 → sponge at 45 : 27 m / 1.25 = 21.6 s = 13.5 T
#      boundary source at x=0  → sponge at 45 : 45 m / 1.25 = 36.0 s = 22.5 T
#  Running a BC case for 16 T therefore leaves the far half of the domain EMPTY,
#  and `max|η|` keeps creeping up as the front advances — which reads as a
#  spurious "growth rate" (+0.028/s measured) even though the amplitude is
#  rock-steady at A. Measured that way once; defaults now cover the transit.
periods = genv_f("BALFEM_PERIODS", wave_gen_kind == "inner" ? 16.0 : 26.0)
Tfinal  = haskey(ENV, "BALFEM_TFINAL") ? genv_f("BALFEM_TFINAL", 0.0) : periods*Twave
save_ev = genv_i("BALFEM_SAVE_EVERY", 10)
mumax   = genv_f("BALFEM_MUMAX", 40.0)

#  ny ≥ 3 is required ONLY for a PERIODIC y-direction — Gridap asserts "a minimum
#  of 3 elements is required in any periodic direction" (CartesianGrids.jl:39).
#  With :wall or :open there is no such constraint and ny=1 is legal, which is
#  the DEFAULT for a genuinely 1-D horizontal case: see the note below.
ybc_sym = Symbol(genv("BALFEM_YBC", "wall"))
ybc_sym in (:wall, :open, :periodic) ||
    error("BALFEM_YBC must be wall, open or periodic (got $ybc_sym)")
(ybc_sym !== :periodic || ny >= 3) && (ny >= 1) ||
    error("BALFEM_NY=$ny is invalid: :periodic needs ny ≥ 3 (Gridap's " *
          "periodic-direction minimum); :wall and :open accept ny ≥ 1.")

#  ---- 1-D CASES ARE NORMAL-INCIDENCE, BY PHYSICS AND BY CONSTRUCTION --------
#
#  A 1-D horizontal domain has ONE propagation direction. An obliquely incident
#  wave has a transverse wavenumber k_y = k sin(theta), i.e. structure ACROSS the
#  flume — and a flume one cell wide cannot represent it. What such a request
#  actually produces is not an oblique wave but an aliased normal-incidence one
#  at the wrong wavenumber, with the transverse component silently dropped: a
#  wrong answer that still runs to completion and looks plausible.
#
#  The same holds for a SHORT-CRESTED sea. `build_airy_state(d)` is called here
#  WITHOUT `directional=true` on purpose; a spread spectrum would put energy at
#  k_y != 0 that this geometry cannot carry.
#
#  This driver exposes no direction knob at all, so oblique content cannot be
#  requested through the intended interface. The guard below exists because a
#  user CAN still export the directional variables the 2-D scripts read, and
#  today they would be SILENTLY IGNORED — which is the worse failure. Refuse
#  rather than warn: a warning is not read until the run has been paid for.
#
for v in ("BALFEM_WAVE_DIR", "BALFEM_NTHETA", "BALFEM_SPREAD_STD",
          "BALFEM_THETA_MAX", "BALFEM_DIRECTIONAL")
    haskey(ENV, v) || continue
    #  WAVE_DIR = 0 is normal incidence and therefore harmless.
    v == "BALFEM_WAVE_DIR" && abs(genv_f(v, 0.0)) < 1e-12 && continue
    error("""
    $v=$(ENV[v]) is set, but this is a 1-D horizontal flume.

    A 1-D domain carries ONE propagation direction. Oblique or short-crested
    content has a transverse wavenumber k_y = k*sin(theta), which a flume one
    cell across cannot represent — the request would be silently aliased onto a
    normal-incidence wave rather than refused, which is why this errors.

    For directional content use the 2-D directional-sea driver instead:
    examples/distributed_small/run_directional_sea_small.jl (y_wall_bc=:open
    with lateral sponges). Unset $v to run this flume.""")
end

#  ---- THE DEFAULT FOR A 1-D HORIZONTAL CASE: ny = 1 with y_wall_bc = :wall ----
#
#  The solver is structurally 2-D, so a 1-D problem is posed as a narrow flume.
#  The cheapest CORRECT way to do that is ONE cell across with SOLID WALLS, not
#  three cells with periodicity:
#
#    * for a normal-incidence wave the exact solution has 𝖴y ≡ 0, and the wall
#      condition 𝖴y = 0 is EXACTLY consistent with it. It approximates nothing.
#      :periodic merely PERMITS 𝖴y ≡ 0 while also admitting a family of
#      y-periodic modes a true 1-D model does not have — an extra mode family
#      sitting in the same wavenumber band as the physics (at the old Ly=3.0 the
#      shortest such mode was 3.0 m against a 4.0 m carrier);
#    * ny=1 + :wall is 2.8x cheaper: 7215 free DOFs against 20202 for ny=3
#      periodic on the same 240-cell streamwise mesh. The wall pins the bottom
#      and top 𝖴y node layers (2886 constrained DOFs = 2 levels x 481 x-nodes
#      x Nσ), leaving only the middle layer free;
#    * a direct LU costs more than linearly in DOFs, so the wall-clock saving is
#      larger than 2.8x.
#
#  Use :periodic only when the case genuinely has oblique or short-crested
#  content, where a solid wall would reflect. All the run/local/run_1d_*.sh
#  cases are normal-incidence and therefore use the ny=1 + :wall default.

#  SIZING FOR A 12-RANK PARTITION (2026-08-06). The flume was 50 m at dx=0.5
#  (8 cells/λ) on 1 core; it is now 60 m at dx=0.25 (16 cells/λ) decomposed
#  12×1. The extra cores went into RESOLUTION rather than length on purpose:
#  at kd=5.5 the group velocity is only c_g=1.25 m/s, so every metre of extra
#  flume costs 0.8 s of simulated transit before the far field is usable —
#  lengthening the domain would have spent the new cores on waiting.

# ---- bathymetry: flat, or a y-invariant submerged bar --------------------
usebar = !flat_bed_flag(1)
hbar   = genv_f("BALFEM_HBAR", 1.0)
xbar   = genv_f("BALFEM_XBAR", 30.0)
wbar   = genv_f("BALFEM_WBAR", 5.0)
#  SHOULDER LENGTH. The bar is two back-to-back tanh shoulders at x = xbar ∓ wbar,
#  each with transition length `sramp`: the bed crosses ~90 % of hbar over 2.9·sramp.
#  The default wbar/3 is the smooth trapezoid every earlier case used; BALFEM_SBAR
#  shortens it into a SQUARE-CROSS-SECTION (box) bar — 0.5 m puts the shoulder in
#  ~1.5 m ≈ 6 cells at dx = 0.25.
#  ⚠ A vertical step is NOT admissible and must not be requested by driving sramp
#  to zero: the residual carries ∇h explicitly (rule 4), so a discontinuous bed is
#  not representable by the model — it would only be aliased by the mesh, and the
#  ∇h terms would then measure the mesh rather than the bathymetry.
sramp  = genv_f("BALFEM_SBAR", wbar / 3.0)
sramp > 0.0 || error("BALFEM_SBAR must be > 0 (a vertical step has no ∇h the model can carry)")
h_bathy = usebar ?
    (x -> d - 0.5*hbar*(tanh((x[1]-(xbar-wbar))/sramp) - tanh((x[1]-(xbar+wbar))/sramp))) :
    nothing
bedtag = usebar ? "bar" : "flat"

# ---- wave generation -----------------------------------------------------
# `wave_bc` decides the mechanism by TYPE (resolve_wave_gen): nothing ⇒ interior
# source; a WaveInput ⇒ Dirichlet boundary generation. The sea state is built
# from a seeded spectrum, so every MPI rank gets an identical component table.
vert0 = assemble_vertical_tensors(M, 1, cbdy_override() === nothing ?
                                  [0.0, 0.728, 1.0] : cbdy_override())
if wave_gen_kind == "inner"
    wave_bc  = nothing
    spL      = genv_f("BALFEM_SPONGE_L", 12.0)
    #  THE SOURCE MUST CLEAR THE SPONGE. The old default `0.2*Lx` happened to
    #  equal `sponge_wL` at both the previous (50 m / 10 m) and the re-sized
    #  (60 m / 12 m) geometry, putting the Gaussian line source exactly on the
    #  sponge edge: the run then reports max|η| pinned at x = x_wm with
    #  eta_max_damped/eta_max_int ≈ 0.95, i.e. the wave is being absorbed as
    #  fast as it is made and there is no clean plane wave to measure.
    #  Default is now the sponge edge + 1.5 wavelengths (λ = 4.0 m at kd=5.5).
    x_wm     = genv_f("BALFEM_XWM", spL + 6.0)
    use_relax = genv_b("BALFEM_RELAX", 0)
elseif wave_gen_kind == "bc"
    wave_bc  = WaveInput(vert0; A=Awave, T=Twave, d=d,
                         T_ramp=genv_f("BALFEM_TRAMP", 2*Twave), profile=bc_profile_sym())
    x_wm     = 0.0
    spL      = 0.0                       # the inflow boundary must not be sponged
    use_relax = genv_b("BALFEM_RELAX", 1)
else                                      # sea
    wave_bc  = WaveInput(vert0, build_airy_state(d); d=d,
                         T_ramp=tramp_val(), profile=bc_profile_sym())
    x_wm     = 0.0
    spL      = 0.0
    use_relax = genv_b("BALFEM_RELAX", 1)
end
spR = genv_f("BALFEM_SPONGE_R", 15.0)

# Refuse a geometry where the interior source sits inside (or within half a
# wavelength of) the left sponge — it silently destroys the case.
if wave_gen_kind == "inner" && x_wm < spL + 2.0
    error("BALFEM_XWM ($(x_wm) m) is inside or too close to the left sponge " *
          "(width $(spL) m). The source would be absorbed as fast as it radiates. " *
          "Move it to at least $(spL + 2.0) m, or shrink BALFEM_SPONGE_L.")
end

# Gauge rake (SEQUENTIAL ONLY — the distributed driver evaluates no points).
# Stations at 20/40/60/80 % of the flume plus a Goda-Suzuki pair at mid-length.
y_c    = Ly/2
gauges = use_mpi ? Tuple{Float64,Float64}[] :
         [(0.2Lx, y_c), (0.4Lx, y_c), (0.5Lx, y_c), (0.5Lx + 1.0, y_c),
          (0.6Lx, y_c), (0.8Lx, y_c)]

#  STANDARDISED OUTPUT NAME — markdown_files/OUTPUT_NAMING_PROPOSAL.md
#      <model>_<domain>_<wave>_<regime>_<nlp>_<bed>_<discr>_<amplitude>_<period>[_extra]
#      P1LFE-2_1d_bcplane_nl_full_flat_Q2Q1_A0.1_T1.6
#  One generator for every driver (output_dir_name in src/utilities.jl); the old
#  per-driver `flume_…`/`small2d_…`/`small_…` prefixes named the SCRIPT, not the case.
#  ⚠ The discretisation token is why this exists: nothing in the old names distinguished
#  Q2/Q1 from Q2/Q2, which is how equal-order runs were compared against the Taylor-Hood
#  MMS campaign for months (CLAUDE.md rule 12b).
_is_sea   = wave_gen_kind == "sea"
_wavekind = _is_sea ? "irr" : "plane"
_gen      = wave_gen_kind == "inner" ? :inner : :bc
_amp      = _is_sea ? genv_f("BALFEM_HS", 0.002) : Awave
_per      = _is_sea ? genv_f("BALFEM_TP", Twave) : Twave
#  Extras: only what a study DELIBERATELY varies, in the spec's order.
_extra = String[]
solver_sym() === :theta && push!(_extra, "theta")
(solver_sym() === :sdirk && tableau_sym() !== :SDIRK_2_2) &&
    push!(_extra, lowercase(replace(String(tableau_sym()), "_" => "")))
genv_b("BALFEM_USE_AD", 0) && push!(_extra, "ad")
genv_b("BALFEM_NLP_INLOOP", 0) && push!(_extra, "inloop")
genv_b("BALFEM_MIXED", 0) && push!(_extra, "mixed")
haskey(ENV, "BALFEM_P_AUX") && push!(_extra, "aux"*ENV["BALFEM_P_AUX"])
let m = lowercase(genv("BALFEM_C3_MASK","both")); m == "both" || push!(_extra, "c3"*m) end
genv_i("BALFEM_QUAD_EXTRA", 0) != 0 && push!(_extra, "q$(genv_i("BALFEM_QUAD_EXTRA",0))")

_name = output_dir_name(; M = M, p_vert = p_vert, ny = ny, y_wall_bc = ybc_sym,
                          wave_kind = _wavekind, wave_gen = _gen,
                          regime = regime_sym(), nl_pressure = nl_pressure_sym(),
                          bed = bedtag, p_u = feord, p_eta = p_eta,
                          amplitude = _amp, period = _per, irregular = _is_sea,
                          nx = nx, nx_in_name = haskey(ENV, "BALFEM_NX"),
                          extra = _extra)
#  unique_output_dir NEVER overwrites: it suffixes _v2, _v3 … A re-executed batch wrote
#  six finished runs over their own output on 2026-09-06 and destroyed them.
outdir = haskey(ENV, "BALFEM_OUTDIR") ? genv("BALFEM_OUTDIR", "") :
         unique_output_dir(joinpath(ROOT, "output", "local_1d"), _name)

if is_rank0()
    @printf("############################################################\n")
    @printf("# QUASI-1D FLUME | gen=%s | %s %s %s | A=%g T=%g\n",
            wave_gen_kind, regime_sym(), nl_pressure_sym(), bedtag, Awave, Twave)
    @printf("#   %s | M=%d | domain %.0f×%.0f m | mesh %d×%d (dx=%.3f) | %s\n",
            model_name, M, Lx, Ly, nx, ny, Lx/nx,
            use_mpi ? "MPI $(genv_i("BALFEM_PX",12))×1 ranks" : "sequential (+gauges)")
    @printf("#   dt=%g s | %g periods → T_final=%.1f s | sponge L/R = %.0f/%.0f, μ=%.0f\n",
            dt, periods, Tfinal, spL, spR, mumax)
    usebar && @printf("#   bar: h=%.2f m at x=%.1f m, half-width %.1f m, shoulder %.2f m → crest depth %.2f m\n",
                      hbar, xbar, wbar, sramp, d - hbar)
    @printf("#   out=%s\n", outdir)
    @printf("############################################################\n")
    flush(stdout)
end

common = (M=M, p_vertical=p_vert, c_bdy=cbdy_override(), p_u=feord, p_eta=p_eta, quad_extra=genv_i("BALFEM_QUAD_EXTRA", 0),
          h_val=d, h_bathy=h_bathy, T_wave=Twave, A_wave=Awave,
          x_wm=x_wm, y_wm=nothing,
          sponge_wL=spL, sponge_wR=spR, sponge_wB=0.0, sponge_wT=0.0, mu_max=mumax,
          T_final=Tfinal, dt=dt,
          regime=regime_sym(), nl_pressure=nl_pressure_sym(),
          flat_bed=flat_bed_flag(1),
          y_wall_bc=ybc_sym, x_wall_bc=false,
          wave_bc=wave_bc, bc_side=bc_side_sym(), bc_profile=bc_profile_sym(),
          relax_bc=use_relax, relax_width=relax_w_val(),
          #  BALFEM_USE_AD=1 swaps the hand Jacobians for AD of the SAME residual.
          #  jacobian_u omits every N derivative and the pressure eta-derivatives by
          #  design (problem.jl: "the N blocks add to the residual but not here"), so
          #  Newton is quasi-Newton. AD makes it exact and is therefore the direct
          #  test of whether a nonlinear failure is a SOLVE defect or an OPERATOR one.
          #  SEQUENTIAL ONLY -- there is no use_ad on the distributed path -- and
          #  substantially slower per assembly, so use it to diagnose, not to run.
          use_ad=genv_b("BALFEM_USE_AD", 0),
          #  BALFEM_NLP_INLOOP=1 evaluates the Class-III L2 projections from the CURRENT
          #  Newton iterate (static condensation) instead of freezing them one step behind.
          #  Removes the O(dt) lag error; only meaningful with nl_pressure=:full.
          nlp_inloop=genv_b("BALFEM_NLP_INLOOP", 0),
          #  BALFEM_MIXED=1 replaces the Class-III L2 projections with GENUINE UNKNOWNS
          #  (5 fields for grad-S only, 7 with grad-b). Projection-free and lag-free, but
          #  AD Jacobians and sequential only -- a DIAGNOSTIC, far slower per step.
          mixed=genv_b("BALFEM_MIXED", 0),
          #  BALFEM_P_AUX sets the FE order of the mixed auxiliary unknowns (𝖦, 𝖥).
          #  Unset → the velocity order (every mixed run before 2026-09-23).
          p_aux=(haskey(ENV, "BALFEM_P_AUX") ? genv_i("BALFEM_P_AUX", feord) : nothing),
          #  BALFEM_C3_MASK isolates WHICH Class-III object is assembled, to answer
          #  which of the two carries the :full instability:
          #    "both" (default) = ordinary :full
          #    "gs"   = the ∇𝖲 family only  (components {1,2,5}, collapsed + 𝓝²'s remainder)
          #    "gb"   = component 4 only    (the ∇𝖻 carrier)
          #    "none" = Class-III suppressed; NOT :native ({3,6,7,8} still assembled)
          c3_mask=(m -> m == "gs" ? (true,false) : m == "gb" ? (false,true) :
                        m == "none" ? (false,false) : (true,true))(lowercase(genv("BALFEM_C3_MASK","both"))),
          output_dir=outdir, save_every=save_ev,
          write_w=write_w_flag(), write_pressure=write_p_flag(), rho=rho_val(),
          solver_type=solver_sym(), tableau=tableau_sym(),
          nl_iter=nl_iter_val(), nl_tol=nl_tol_val(),
          print_every=genv_i("BALFEM_PRINT_EVERY", 10),
          check_every=genv_i("BALFEM_CHECK_EVERY", 0),
          diag_every=genv_i("BALFEM_DIAG_EVERY", 0),
          div_factor=genv_f("BALFEM_DIV_FACTOR", 20.0))

if use_mpi
    px = genv_i("BALFEM_PX", 12)
    nx % px == 0 || error("BALFEM_NX ($nx) must be divisible by BALFEM_PX ($px)")
    diags, vert, prob = setup_and_run_distributed(;
        cpu_grid=(px, 1), domain=(0.0, Lx, 0.0, Ly), partition=(nx, ny),
        ls_rtol=ls_rtol_val(), ls_maxiter=ls_maxiter_val(), krylov_m=krylov_m_val(), precond=precond_sym(),
        common...)
else
    diags, vert, prob = setup_and_run(;
        domain=((0.0, Lx), (0.0, Ly)), partition=(nx, ny),
        gauges=gauges, common...)
end

#  `tag` was a leftover from before the output_dir_name rename (rule 2c) and was defined
#  NOWHERE, so every SUCCESSFUL run threw UndefVarError here -- after writing every solve
#  result, VTK file and diagnostics row. A batch runner reading exit codes scored a completed
#  100-period run as a failure (rule 35). `_name` is the case identity this line wanted.
is_rank0() && @printf("flume_1d [%s] done: %d steps → %s\n", _name, length(diags), outdir)
