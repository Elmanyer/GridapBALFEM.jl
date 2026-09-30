# ==============================================================
#  run_periodic_1d.jl — the x-PERIODIC closed flume: the INTERIOR discretisation alone
#
#  THE QUESTION (LaTeX chapter "Stability Analysis", §7.7). Every :full failure on the
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
#  MIXED LAYOUT: the auxiliary 𝖦 is solved from its own constraint at t₀
#  (mixed_consistent_ic), so the hot start is consistent.
#
#  BROKEN / C⁰-IP (src/broken.jl, markdown_files/BROKEN_FORMULATION_PLAN.md):
#     BALFEM_BROKEN=1  Class-III 𝓚/𝓟 blocks by the distributional gradient (needs :full)
#     BALFEM_CIP_GU    C⁰-IP γ_u on ⟦∂ₙ𝖴⟧   BALFEM_CIP_GE  γ_η on ⟦∂ₙη⟧   BALFEM_CIP_HEXP  s (2)
#
#  KNOBS (all BALFEM_*, as run_flume_1d.jl): NL_PRESSURE, MIXED, P_AUX, C3_MASK, FE_ORDER,
#  P_ETA, AWAVE, TWAVE, D, DT, PERIODS, SOLVER/TABLEAU, USE_AD, NLP_INLOOP, plus
#     BALFEM_NLAMBDA   wavelengths in the box          (default 1)
#     BALFEM_NCELL     cells per wavelength            (default 16  → dx = λ/16 ≈ 0.25 m)
#     BALFEM_NOISE     seed amplitude                  (default 1e-8; 0 = none)
#  RUN:  run/local/run_1dper_*.sh     (sequential: the case is a few hundred DOFs)
# ==============================================================

include(joinpath(@__DIR__, "..", "distributed", "_dist_common.jl"))
using Gridap.TensorValues: VectorValue     # the stacked 𝖴x initial condition

get!(ENV, "BALFEM_REGIME",      "nonlinear")
get!(ENV, "BALFEM_NL_PRESSURE", "native")
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
flat_bed_flag(1) || error("run_periodic_1d: a periodic flume needs a FLAT bed (BALFEM_FLAT_BED=1)")

# ---- the base wave: the discrete linear eigenmode, no ramp ----------------------------
vert0 = assemble_vertical_tensors(M, p_vert, Vector{Float64}(resolve_cbdy(M, cbdy_override(), p_vert)))
Nσ    = vert0.N_dof
wi    = WaveInput(vert0; A = Awave, T = Twave, d = d, T_ramp = 0.0, profile = :model)
k0    = wi.ks[1]
λ0    = 2π / k0
Lx    = nlam * λ0                         # EXACTLY n model wavelengths ⇒ the wave is periodic
nx    = nlam * ncell
dx    = Lx / nx
Ly    = dx                                # one isotropic cell across, solid walls (rule 12)
inc   = incident_fields(wi)

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
solver_sym() === :theta && push!(_extra, "theta")
(solver_sym() === :sdirk && tableau_sym() !== :SDIRK_2_2) &&
    push!(_extra, lowercase(replace(String(tableau_sym()), "_" => "")))
genv_b("BALFEM_USE_AD", 0) && push!(_extra, "ad")
genv_b("BALFEM_NLP_INLOOP", 0) && push!(_extra, "inloop")
genv_b("BALFEM_MIXED", 0) && push!(_extra, "mixed")
genv_b("BALFEM_BROKEN", 0) && push!(_extra, "broken")
_cgu = genv_f("BALFEM_CIP_GU", 0.0); _cge = genv_f("BALFEM_CIP_GE", 0.0)
_cgu > 0 && push!(_extra, @sprintf("cipu%g", _cgu))
_cge > 0 && push!(_extra, @sprintf("cipe%g", _cge))
haskey(ENV, "BALFEM_CIP_HEXP") && push!(_extra, "hexp" * ENV["BALFEM_CIP_HEXP"])
haskey(ENV, "BALFEM_P_AUX") && push!(_extra, "aux" * ENV["BALFEM_P_AUX"])
let m = lowercase(genv("BALFEM_C3_MASK", "both")); m == "both" || push!(_extra, "c3" * m) end
_name = output_dir_name(; M = M, p_vert = p_vert, ny = 1, y_wall_bc = :wall,
                          wave_kind = "plane", wave_gen = :ic,
                          regime = regime_sym(), nl_pressure = nl_pressure_sym(),
                          bed = "flat", p_u = feord, p_eta = p_eta,
                          amplitude = Awave, period = Twave, irregular = false,
                          nx = nx, nx_in_name = false, extra = _extra)
outdir = haskey(ENV, "BALFEM_OUTDIR") ? genv("BALFEM_OUTDIR", "") :
         unique_output_dir(joinpath(ROOT, "output", "local_1d", "periodic"), _name)

@printf("############################################################\n")
@printf("# PERIODIC 1-D FLUME | %s %s flat | A=%g T=%g | %s\n",
        regime_sym(), nl_pressure_sym(), Awave, Twave, "P$(p_vert)LFE-$(M)")
@printf("#   box = %d × λ_model = %.4f m | k0 = %.4f (kd=%.2f, κa=%.3f) | %d cells/λ → dx=%.4f\n",
        nlam, Lx, k0, k0 * d, k0 * Awave, ncell, dx)
@printf("#   NO inflow, NO relaxation, NO sponge, NO source | seed %.1e on %d harmonics\n",
        noise, nmax)
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
    regime = regime_sym(), nl_pressure = nl_pressure_sym(), flat_bed = true,
    use_ad = genv_b("BALFEM_USE_AD", 0),
    nlp_inloop = genv_b("BALFEM_NLP_INLOOP", 0),
    mixed = genv_b("BALFEM_MIXED", 0),
    broken = genv_b("BALFEM_BROKEN", 0),
    cip_gamma_u = _cgu, cip_gamma_eta = _cge,
    cip_hexp = genv_f("BALFEM_CIP_HEXP", 2.0),
    p_aux = (haskey(ENV, "BALFEM_P_AUX") ? genv_i("BALFEM_P_AUX", feord) : nothing),
    c3_mask = (m -> m == "gs" ? (true, false) : m == "gb" ? (false, true) :
                    m == "none" ? (false, false) : (true, true))(lowercase(genv("BALFEM_C3_MASK", "both"))),
    output_dir = outdir, save_every = genv_i("BALFEM_SAVE_EVERY", 5),
    #  sub-cell VTK sampling: 2·p_u samples per cell resolves every polynomial mode up to
    #  the NODE Nyquist (vertex-only output would alias all sub-element content away)
    vtk_nsubcells = genv_i("BALFEM_VTK_NSUB", 2 * feord),
    solver_type = solver_sym(), tableau = tableau_sym(),
    nl_iter = nl_iter_val(), nl_tol = nl_tol_val(),
    print_every = genv_i("BALFEM_PRINT_EVERY", 25),
    check_every = genv_i("BALFEM_CHECK_EVERY", 0),     # θ-scheme self-check; meaningless under SDIRK
    diag_every = genv_i("BALFEM_DIAG_EVERY", 5),
    eta_ref = Awave,
    div_factor = genv_f("BALFEM_DIV_FACTOR", 20.0))

@printf("periodic_1d [%s] done: %d steps → %s\n", _name, length(diags), outdir)
