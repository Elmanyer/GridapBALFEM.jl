# ==============================================================
#  run_periodic_1d.jl — the x-PERIODIC closed flume: the INTERIOR discretisation alone
#
#  THE QUESTION (LaTeX chapter "Stability Analysis", §7.7). Every full-pressure failure on the
#  60 m flume pins at the inflow (x ≈ 0.12–0.5 m) or at the relaxation-zone edge. The flume
#  cannot tell apart
#     (a) an INTERIOR discretisation instability that the boundary merely seeds, from
#     (b) a boundary treatment (linear Dirichlet data + relaxation towards a LINEAR target,
#         forced onto a NONLINEAR model) that keeps injecting grid-scale content.
#  This driver removes every boundary mechanism: no inflow, no relaxation zone, no sponge,
#  no source. What remains is the interior semi-discrete operator, the time integrator and
#  the nonlinear state.
#
#  THE NULL IS KNOWN. The domain is exactly BALFEM_NLAMBDA model wavelengths long. With ONE
#  wavelength only the harmonics n·k0 exist, so Benjamin–Feir sidebands are excluded by
#  construction, and a Stokes wave at κa ≈ 0.16 is superharmonically stable: in the continuum
#  NOTHING should grow. Any sustained growth of the high-k band is either a discretisation
#  instability or ill-posedness of the continuum — separated by the dx ladder (growth rate
#  converging ⇒ physical/well-posed; growing like dx^-α ⇒ grid-scale).
#
#  THE STATE. The solver's own discrete linear eigenmode (WaveInput :model polarisation) at
#  amplitude A, hot-started at t = 0 with NO ramp, plus a deterministic broadband seed of
#  size BALFEM_NOISE on η and 𝖴x at every representable harmonic (golden-ratio phases —
#  no RNG dependency, identical on every run). The seed makes growth measurable from
#  ~1e-8 instead of round-off, so a Floquet rate shows within tens of periods.
#  ⚠ The linear eigenmode is not an exact nonlinear wave: it sheds BOUNDED free harmonics.
#  Read growth from the log-linear rise of the HIGH-k band energy
#  (postprocessing/examples/periodic_growth.jl), never from max η.
#
#  v2: BALFEM_NL_PRESSURE=0/1 — all eight 𝓝 components off/on (Class III by the broken
#  formulation; there is no other treatment and no component mask).
#
#  SKELETON STABILISATION (src/broken.jl):
#     BALFEM_STAB      jumpgrad | ghostvolume        BALFEM_CIP_ORDER  jumpgrad orders (≤ 2)
#     BALFEM_CIP_GU    γ_u (velocity)   BALFEM_CIP_GE  γ_η (surface)   BALFEM_CIP_HEXP  s (2)
#
#  KNOBS (all BALFEM_*, as run_flume_1d.jl): NL_PRESSURE, FE_ORDER, P_ETA, AWAVE, TWAVE, D,
#  DT, PERIODS, SOLVER/TABLEAU, USE_AD, plus
#     BALFEM_NLAMBDA   wavelengths in the box          (default 1)
#     BALFEM_NCELL     cells per wavelength            (default 16  → dx = λ/16 ≈ 0.25 m)
#     BALFEM_KA        steepness k₀A (model k₀); overrides BALFEM_AWAVE   (default unset)
#     BALFEM_NOISE     seed amplitude                  (default 1e-8; 0 = none)
#     BALFEM_BED_SLOPE max bed slope of a PERIODIC cosine bed (default 0 = flat; needs FLAT_BED=0)
#     BALFEM_BED_NW    bed undulations per box         (default 1)
#
#  TIME INTEGRATOR (src/timeloop.jl build_time_integrator; the output name carries the choice):
#     BALFEM_SOLVER    theta | sdirk (= rk: any Gridap tableau) | gen_alpha       (default sdirk)
#     BALFEM_TABLEAU   SDIRK_*/DIRK_* (implicit) or EXRK_* (explicit)            (default SDIRK_2_2)
#     BALFEM_THETA_T   θ of the θ-method, 0.5 = Crank–Nicolson                    (default 0.5)
#     BALFEM_RHO_INF   ρ∞ of generalised-α, 1 = no damping, 0 = strongest         (default 0.5)
#
#  VARIABLE BED (models 4/6). d(x) = D + δ·cos(2π·nw·x/Lx): smooth (C∞, so ∇h and the bed Hessian
#  are exact) and periodic with the box, so the seam is invisible. δ = slope/(2π·nw/Lx) is set
#  by the MAXIMUM SLOPE, the quantity that governs the ∇h terms. With nw = 1 and one wavelength
#  in the box the bed wavenumber is k0, NOT Bragg-resonant (class-I Bragg needs k_bed = 2k0, nw = 2):
#  the bed couples harmonics n ↔ n±1 without reflecting the carrier. The initial eigenmode is
#  the flat-bed one at the mean depth D; on a mild bed it sheds a small bounded free part.
#  RUN:  run/local/run_1dper_*.sh     (sequential: the case is a few hundred DOFs)
# ==============================================================

include(joinpath(@__DIR__, "..", "distributed", "_dist_common.jl"))
using Gridap.TensorValues: VectorValue     # the stacked 𝖴x initial condition

get!(ENV, "BALFEM_REGIME",      "nonlinear")
get!(ENV, "BALFEM_NL_PRESSURE", "0")
get!(ENV, "BALFEM_FLAT_BED",    "1")

M       = genv_i("BALFEM_M", 2)
p_vert  = genv_i("BALFEM_P_VERT", 1)
feord   = genv_i("BALFEM_FE_ORDER", 2)
p_eta   = genv_i("BALFEM_P_ETA", feord - 1)
d       = genv_f("BALFEM_D", 3.5)
Twave   = genv_f("BALFEM_TWAVE", 1.6)
Awave   = genv_f("BALFEM_AWAVE", 0.10)
dt      = genv_f("BALFEM_DT", 0.04)
nlam    = genv_i("BALFEM_NLAMBDA", 1)
ncell   = genv_i("BALFEM_NCELL", 16)
noise   = genv_f("BALFEM_NOISE", 1e-8)
periods = genv_f("BALFEM_PERIODS", 100.0)
Tfinal  = periods * Twave
bslope  = genv_f("BALFEM_BED_SLOPE", 0.0)
bnw     = genv_i("BALFEM_BED_NW", 1)
flatbed = flat_bed_flag(1)
flatbed == (bslope == 0) ||
    error("run_periodic_1d: BALFEM_FLAT_BED=$(Int(flatbed)) contradicts BALFEM_BED_SLOPE=$bslope " *
          "(a variable bed needs FLAT_BED=0 and a slope > 0; a flat bed needs FLAT_BED=1, slope 0)")

# ---- the base wave: the discrete linear eigenmode, no ramp ----------------------------
vert0 = assemble_vertical_tensors(M, p_vert, Vector{Float64}(resolve_cbdy(M, cbdy_override(), p_vert)))
Nσ    = vert0.N_dof
#  BALFEM_KA, the steepness k₀A with k₀ the MODEL's wavenumber at T, overrides BALFEM_AWAVE: a
#  campaign over T at fixed steepness needs no table of amplitudes (2026-11-11).
haskey(ENV, "BALFEM_KA") &&
    (Awave = genv_f("BALFEM_KA", 0.0) / GridapBALFEM.model_wavenumber(vert0, 2π / Twave, d, 9.81))
wi    = WaveInput(vert0; A = Awave, T = Twave, d = d, T_ramp = 0.0, profile = :model)
k0    = wi.ks[1]
λ0    = 2π / k0
Lx    = nlam * λ0                         # EXACTLY n model wavelengths ⇒ the wave is periodic
nx    = nlam * ncell
dx    = Lx / nx
Ly    = dx                                # one isotropic cell across, solid walls (rule 12)
inc   = incident_fields(wi)
kbed  = 2π * bnw / Lx
δbed  = bslope / kbed                     # bed amplitude from the maximum slope
hfun  = flatbed ? nothing : (x -> d + δbed * cos(kbed * x[1]))

# ---- deterministic broadband seed at every representable harmonic --------------------
#  harmonics n = 1 … p_u·nx of the BOX wavenumber 2π/Lx (up to the node Nyquist), equal
#  amplitude, golden-ratio phases; η and each 𝖴x component get independent phases.
nmax  = feord * nx
kL    = 2π / Lx
φg    = (n, j) -> 2π * mod(n * 0.6180339887498949 + j * 0.7548776662466927, 1.0)
anorm = noise / sqrt(nmax)
seedf = (x, j) -> anorm * sum(cos(n * kL * x + φg(n, j)) for n in 1:nmax)
eta0  = x -> inc.eta(x, 0.0) + (noise > 0 ? seedf(x[1], 0) : 0.0)
ux0   = x -> begin
    u = inc.ux(x, 0.0)
    noise > 0 || return u
    VectorValue(ntuple(j -> u[j] + seedf(x[1], j), Nσ)...)
end

# ---- output name (the standard grammar; wave token `icplane`, extras mark the box) ----
_extra = String["xper", "nlam$(nlam)", "n$(ncell)"]
noise > 0 && push!(_extra, @sprintf("seed%.0e", noise))
integrator_tag() === nothing || push!(_extra, integrator_tag())
genv_b("BALFEM_USE_AD", 0) && push!(_extra, "ad")
is_resume() && push!(_extra, "resume")
#  γ: unset ⇒ the driver default (DEFAULT_CIP_GAMMA with 𝓝, where the C⁰-CIP is mandatory; 0 otherwise)
_cgdef = nl_pressure_flag() ? DEFAULT_CIP_GAMMA : 0.0
_cgu = genv_f("BALFEM_CIP_GU", _cgdef); _cge = genv_f("BALFEM_CIP_GE", _cgdef)
_cgu > 0 && push!(_extra, @sprintf("cipu%g", _cgu))
_cge > 0 && push!(_extra, @sprintf("cipe%g", _cge))
haskey(ENV, "BALFEM_CIP_HEXP") && push!(_extra, "hexp" * ENV["BALFEM_CIP_HEXP"])
genv_i("BALFEM_CIP_ORDER", 1) > 1 && push!(_extra, "cipord" * ENV["BALFEM_CIP_ORDER"])
lowercase(genv("BALFEM_STAB", "jumpgrad")) == "ghostvolume" && push!(_extra, "ghost")
_name = output_dir_name(; M = M, p_vert = p_vert, ny = 1, y_wall_bc = :wall,
                          wave_kind = "plane", wave_gen = :ic,
                          regime = regime_sym(), nl_pressure = nl_pressure_flag(),
                          bed = flatbed ? "flat" : @sprintf("cosbed%gs%g", bnw, bslope), p_u = feord, p_eta = p_eta,
                          amplitude = Awave, period = Twave, irregular = false,
                          nx = nx, nx_in_name = false, extra = _extra)
outdir = haskey(ENV, "BALFEM_OUTDIR") ? genv("BALFEM_OUTDIR", "") :
         unique_output_dir(joinpath(ROOT, "output", "local_1d", "periodic"), _name)

@printf("############################################################\n")
@printf("# PERIODIC 1-D FLUME | %s %s %s | A=%g T=%g | %s\n",
        regime_sym(), nl_pressure_flag(), flatbed ? "flat" : "cosine bed", Awave, Twave, "P$(p_vert)LFE-$(M)")
@printf("#   box = %d × λ_model = %.4f m | k0 = %.4f (kd=%.2f, κa=%.3f) | %d cells/λ → dx=%.4f\n",
        nlam, Lx, k0, k0 * d, k0 * Awave, ncell, dx)
@printf("#   NO inflow, NO relaxation, NO sponge, NO source | seed %.1e on %d harmonics\n",
        noise, nmax)
flatbed || @printf("#   BED d(x) = %.3f + %.4f·cos(%.4f x): %d undulation(s)/box, max slope %.3f, d ∈ [%.3f, %.3f] m\n",
                   d, δbed, kbed, bnw, bslope, d - δbed, d + δbed)
@printf("#   dt=%g s | %g periods → T_final=%.1f s | out=%s\n", dt, periods, Tfinal, outdir)
@printf("############################################################\n")
flush(stdout)

diags, vert, prob = setup_and_run(;
    M = M, p_vertical = p_vert, c_bdy = cbdy_override(), p_u = feord, p_eta = p_eta,
    quad_extra = genv_i("BALFEM_QUAD_EXTRA", 0),
    domain = ((0.0, Lx), (0.0, Ly)), partition = (nx, 1),
    x_periodic = true, y_wall_bc = :wall, x_wall_bc = false,
    h_val = d, T_wave = Twave, A_wave = Awave,
    sponge_wL = 0.0, sponge_wR = 0.0, sponge_wB = 0.0, sponge_wT = 0.0,
    eta0_func = eta0, ux0_func = ux0,
    T_final = Tfinal, dt = dt,
    regime = regime_sym(), nl_pressure = nl_pressure_flag(), flat_bed = flatbed, h_bathy = hfun,
    use_ad = genv_b("BALFEM_USE_AD", 0),
    cip_gamma_u = _cgu, cip_gamma_eta = _cge,
    cip_hexp = genv_f("BALFEM_CIP_HEXP", 2.0),
    cip_order = genv_i("BALFEM_CIP_ORDER", 1),
    stabilization = Symbol(lowercase(genv("BALFEM_STAB", "jumpgrad"))),
    output_dir = outdir, save_every = genv_i("BALFEM_SAVE_EVERY", 5),
    #  sub-cell VTK sampling: 2·p_u samples per cell resolves every polynomial mode up to
    #  the NODE Nyquist (vertex-only output would alias all sub-element content away)
    vtk_nsubcells = genv_i("BALFEM_VTK_NSUB", 2 * feord),
    solver_type = solver_sym(), tableau = tableau_sym(), theta = theta_t_val(), rho_inf = rho_inf_val(),
    nl_iter = nl_iter_val(), nl_tol = nl_tol_val(),
    print_every = genv_i("BALFEM_PRINT_EVERY", 25),
    check_every = genv_i("BALFEM_CHECK_EVERY", 0),     # θ-scheme self-check; meaningless under SDIRK
    diag_every = genv_i("BALFEM_DIAG_EVERY", 5),
    eta_ref = Awave,
    div_factor = genv_f("BALFEM_DIV_FACTOR", 20.0),
    checkpoint_kw()...)

@printf("periodic_1d [%s] done: %d steps → %s\n", _name, length(diags), outdir)
