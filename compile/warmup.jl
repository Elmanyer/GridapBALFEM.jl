# ==============================================================
#  warmup.jl — trace-compile the COMPLETE BALFE-M solver into the sysimage
#
#  Run as PackageCompiler's precompile_execution_file (see compile.jl). Every case below is a TINY
#  run (4×4 cells, 2 time steps) whose only purpose is to make Julia compile — and PackageCompiler
#  record — every specialisation a production run will need, so that ranks LOAD the solver instead
#  of JIT-compiling it (~30–45 min per rank for models 1–4, ~3 h for the Q3/Q2 𝓝 path).
#
#  WHAT A CASE MUST DIFFER IN TO BE WORTH TRACING. Compiled code is keyed on TYPES, not values:
#    * the execution path — CellField (sequential) vs DistributedCellField + GMRES (distributed);
#    * the model — regime, nl_pressure (Class-III cellwise + facet terms + mandatory C⁰-CIP) and
#      flat_bed (the ∇h terms and the bed Hessian) select different residual/Jacobian expressions;
#    * the time integrator (RungeKutta / ThetaMethod / GeneralizedAlpha stage operators);
#    * wave generation — interior source vs Dirichlet boundary generation (a TransientTrialFESpace),
#      the WaveSpec sea-state conversion, the relaxation zone, the hot start;
#    * boundaries and ICs — y walls / open / periodic, x-periodic box, closed basin release;
#    * the skeleton stabiliser (:jumpgrad order 1–2, :ghostvolume) and, distributed, the GMRES
#      preconditioner (:jacobi / :schwarz / :gs);
#    * the vertical basis SIZE Nσ = M·p + 1 — a TYPE PARAMETER (VectorValue{Nσ}): every extra Nσ
#      compiles the whole solver again (BALFEM_WARMUP_VERT below).
#    Values (mesh size, FE order, dt, amplitudes) do not matter, so the cases are minimal. Every
#    case switches ON the per-step report, the field diagnostics (incl. the CIP power), the
#    governing-equation check, VTK output and the w/p reconstruction, so those are baked too.
#
#  ⚠ LIMIT OF ANY SYSIMAGE: the bathymetry is a user closure (`h_bathy`), and Gridap specialises
#  every term that contains H = d + η on that closure's TYPE. The flat bed is built inside the
#  package (stable type, fully covered); a variable bed defined in a launcher script has a type
#  of its own, so the H-dependent kernels of a variable-bed run still JIT. The variable-bed cases
#  below bake everything that does not depend on that type (the ∇h code paths, the drivers).
#
#  EVERY CASE IS GUARDED: a failure is reported (and listed in the summary at the end) but does
#  not throw away the build after the cold compile is paid. Read the summary in the build log.
#
#  KNOBS (environment):
#    BALFEM_WARMUP_VERT   vertical bases to trace, "M:p" list          default "2:1"  (P1LFE-2, Nσ=3,
#                         e.g. "2:1,3:1,2:2"                                            the production basis)
#    BALFEM_WARMUP_EXTRA  optional groups, comma list                  default ""     (none)
#                           mms  — the manufactured-solution drivers (verification), all 6 models
#                           ad   — the AD Jacobian path (use_ad=true; sequential, diagnosis only)
#    BALFEM_WARMUP_ONLY   restrict to case names containing this text  default ""     (all)
#    BALFEM_WARMUP_DRYRUN 1 = run nothing; check every case's keywords against its driver's
#                         signature (seconds, no compile) — do this before a multi-hour build
#
#  GridapBALFEM is loaded as a PACKAGE here, and compile.jl lists :GridapBALFEM among the baked
#  packages. That combination is what makes this trace stick: the specialisations compiled below
#  belong to the package, so they are retained in the image. Under the old `include()` loading
#  they belonged to a throw-away `Main.GridapBALFEM` and were discarded.
# ==============================================================

const ROOT = normpath(joinpath(@__DIR__, ".."))
using GridapBALFEM
using Gridap: VectorValue
using Printf

# ── knobs ──────────────────────────────────────────────────────────────────────
const VERT  = [Tuple(parse.(Int, split(strip(s), ':'))) for s in
               split(get(ENV, "BALFEM_WARMUP_VERT", "2:1"), ',') if !isempty(strip(s))]
const EXTRA = Set(strip.(split(lowercase(get(ENV, "BALFEM_WARMUP_EXTRA", "")), ',')))
const ONLY  = get(ENV, "BALFEM_WARMUP_ONLY", "")
const DRY   = get(ENV, "BALFEM_WARMUP_DRYRUN", "0") == "1"

# ── the six models (CLAUDE.md §4): regime, 𝓝, bed; 𝓝 needs p_eta ≥ 2 (Q3/Q2) ───────────────
const MODELS = [
    (id = 1, regime = :linear,    nlp = false, flat = true,  p_u = 2),
    (id = 2, regime = :linear,    nlp = false, flat = false, p_u = 2),
    (id = 3, regime = :nonlinear, nlp = false, flat = true,  p_u = 2),
    (id = 4, regime = :nonlinear, nlp = false, flat = false, p_u = 2),
    (id = 5, regime = :nonlinear, nlp = true,  flat = true,  p_u = 3),
    (id = 6, regime = :nonlinear, nlp = true,  flat = false, p_u = 3),
]

# ── geometry: a 4 m × 4 m tank, h = 3.5 m, a smooth bar for the variable-bed models ─────────
const D0 = 3.5
"Smooth submerged bar across the tank (a named function: one type for every variable-bed case)."
warmup_bar(x) = D0 - 0.5*(tanh((x[1] - 1.5)/0.4) - tanh((x[1] - 2.5)/0.4))
bed(m) = m.flat ? nothing : warmup_bar

# ── bookkeeping: every case is guarded and timed ───────────────────────────────────────────
const RESULTS = Tuple{String,Bool,Float64}[]
"Run `driver(; kw...)` as the case `name`, guarded and timed (DRY: only check the keywords)."
function trace(name::String, driver, kw::NamedTuple)
    isempty(ONLY) || occursin(ONLY, name) || return nothing
    if DRY
        bad = setdiff(keys(kw), Base.kwarg_decl(only(methods(driver))))
        isempty(bad) || @warn "warmup (dry run): case `$name` passes keywords `$driver` does not accept" bad
        push!(RESULTS, (name, isempty(bad), 0.0))
        return nothing
    end
    println("warmup ▸ ", name); flush(stdout)
    t0 = time()
    try
        driver(; kw...)
        push!(RESULTS, (name, true, time() - t0))
    catch e
        push!(RESULTS, (name, false, time() - t0))
        @warn "warmup: case `$name` FAILED and is NOT baked — runs of this kind will JIT it" exception = (e, catch_backtrace())
    end
    flush(stdout)
    return nothing
end
outdir(name) = joinpath(ROOT, "output", "warmup", replace(name, r"[^A-Za-z0-9_.-]" => "_"))

#  Options every case shares: a 2-step run with every per-step facility switched on.
common(M, p, m; name) = (
    M = M, p_vertical = p, c_bdy = nothing,
    p_u = m.p_u, p_eta = m.p_u - 1, h_val = D0, h_bathy = bed(m),
    regime = m.regime, nl_pressure = m.nlp, flat_bed = m.flat,
    T_wave = 1.6, A_wave = 0.001, T_final = 0.1, dt = 0.05,
    print_every = 1, diag_every = 1, check_every = 1, save_every = 1,
    write_w = true, write_pressure = true,
    output_dir = outdir(name))

#  Interior-source flume (the default generation): source at x = 1.5 m, x-sponges.
inner = (x_wm = 1.5, y_wm = nothing, sponge_wL = 0.5, sponge_wR = 0.5, mu_max = 5.0)
#  Dirichlet boundary generation on the left edge (no inflow sponge), with a relaxation zone.
bcgen = (wave_bc = :regular, bc_side = :left, sponge_wL = 0.0, sponge_wR = 0.5, mu_max = 5.0,
         relax_bc = true, relax_width = 1.0)

"Seeded JONSWAP AiryState (the production sea-state path; same construction as
examples/distributed/_dist_common.jl, with both seeds pinned)."
function warmup_sea_state(h)
    spec   = WaveSpec.ContinuousSpectrums.JONSWAP(0.02, 1.6, 3.3)
    dspec  = WaveSpec.SpectralSpreading.DiscreteSpectralSpreading(
                 spec, WaveSpec.SpectralSampling.UniformSampling(), 1/(2.5*1.6), 1/(0.75*1.6), 5;
                 domain = WaveSpec.SpectralSampling.Frequency, mess = false)
    spread = WaveSpec.AngularSpreading.change_seed!(
                 WaveSpec.AngularSpreading.DiscreteAngularSpreading(0.0; seed = 1), 1)
    state  = WaveSpec.AiryWaves.AiryState(dspec, spread, h)
    return WaveSpec.AiryWaves.change_seed!(state, 1)
end

zero_u(Nσ) = x -> VectorValue(ntuple(_ -> 0.0, Nσ)...)
eta_bump   = x -> 0.01 * cos(2π * x[1] / 4.0)

# ── SEQUENTIAL (CellField; LU + Newton) ────────────────────────────────────────────────────
function sequential_cases(M, p)
    Nσ  = M*p + 1
    tag = "P$(p)LFE-$(M)"
    seq(name; kw...) = trace(name, setup_and_run, values(kw))
    tank = (domain = ((0.0, 4.0), (0.0, 4.0)), partition = (4, 4), gauges = [(2.0, 2.0)])

    # every model, under both production integrators
    for m in MODELS, st in (:theta, :sdirk)
        n = "seq $tag model$(m.id) $st inner"
        seq(n; tank..., inner..., common(M, p, m; name = n)..., solver_type = st)
    end
    # the remaining integrators (time-loop specialisations), on a model without and with 𝓝
    for m in MODELS[[1, 5]], st in (:gen_alpha, :rk3)
        n = "seq $tag model$(m.id) $st inner"
        seq(n; tank..., inner..., common(M, p, m; name = n)..., solver_type = st)
    end
    # Dirichlet boundary generation + relaxation, y-periodic; regular wave and a WaveSpec sea
    for m in MODELS[[1, 3, 5]]
        n = "seq $tag model$(m.id) bcgen regular yperiodic"
        seq(n; tank..., bcgen..., common(M, p, m; name = n)..., y_wall_bc = :periodic, solver_type = :theta)
    end
    let m = MODELS[1], n = "seq $tag model1 bcgen sea"
        seq(n; tank..., bcgen..., common(M, p, m; name = n)..., wave_bc = warmup_sea_state(D0),
            solver_type = :theta)
    end
    let m = MODELS[1], n = "seq $tag model1 bcgen hotstart"
        seq(n; tank..., bcgen..., common(M, p, m; name = n)..., T_ramp = 0.0, ic_from_bc = true,
            solver_type = :theta)
    end
    # open lateral boundaries with y-sponges
    let m = MODELS[1], n = "seq $tag model1 yopen"
        seq(n; tank..., inner..., common(M, p, m; name = n)..., y_wall_bc = :open,
            sponge_wB = 0.5, sponge_wT = 0.5, solver_type = :theta)
    end
    # the closed x-periodic box (the stability study): models 3, 5, 6 under CN
    for m in MODELS[[3, 5, 6]]
        n = "seq $tag model$(m.id) xperiodic box"
        seq(n; domain = ((0.0, 4.0), (0.0, 1.0)), partition = (4, 1),
            common(M, p, m; name = n)..., A_wave = 0.0, x_periodic = true, x_wall_bc = false,
            sponge_wL = 0.0, sponge_wR = 0.0, eta0_func = eta_bump, ux0_func = zero_u(Nσ),
            solver_type = :theta)
    end
    # closed-basin release (x walls)
    for m in MODELS[[1, 5]]
        n = "seq $tag model$(m.id) release xwalls"
        seq(n; tank..., common(M, p, m; name = n)..., A_wave = 0.0, x_wall_bc = true,
            sponge_wL = 0.0, sponge_wR = 0.0, eta0_func = eta_bump, solver_type = :theta)
    end
    # skeleton stabilisers: optional CIP without 𝓝, order-2 jumps, and the ghost penalty
    let m = MODELS[3], n = "seq $tag model3 cip optional"
        seq(n; tank..., inner..., common(M, p, m; name = n)..., cip_gamma_u = 0.1, cip_gamma_eta = 0.1,
            solver_type = :theta)
    end
    for m in MODELS[[5, 6]]
        n = "seq $tag model$(m.id) jumpgrad order2"
        seq(n; tank..., inner..., common(M, p, m; name = n)..., cip_order = 2, solver_type = :theta)
        n = "seq $tag model$(m.id) ghostvolume"
        seq(n; tank..., inner..., common(M, p, m; name = n)..., stabilization = :ghostvolume,
            cip_gamma_u = 0.01, cip_gamma_eta = 0.01, solver_type = :theta)
    end
    # AD Jacobians (diagnosis only; opt-in)
    if "ad" in EXTRA
        for m in MODELS[[3, 5]]
            n = "seq $tag model$(m.id) use_ad"
            seq(n; tank..., inner..., common(M, p, m; name = n)..., use_ad = true, solver_type = :theta)
        end
    end
end

# ── DISTRIBUTED (DistributedCellField; GMRES + Newton), one rank ───────────────────────────
#  Types do not depend on the rank count, so cpu_grid = (1, 1) bakes what N ranks run.
function distributed_cases(M, p)
    tag = "P$(p)LFE-$(M)"
    dist(name; kw...) = trace(name, setup_and_run_distributed, values(kw))
    tank = (cpu_grid = (1, 1), domain = (0.0, 4.0, 0.0, 4.0), partition = (4, 4))
    #  setup_and_run_distributed has no gauges / ux0_func / use_ad / x_periodic.
    dcommon(m; name) = common(M, p, m; name = name)

    for m in MODELS, st in (:theta, :sdirk)
        n = "dist $tag model$(m.id) $st inner"
        dist(n; tank..., inner..., dcommon(m; name = n)..., solver_type = st)
    end
    # the GMRES preconditioners
    for pc in (:schwarz, :gs)
        m = MODELS[1]; n = "dist $tag model1 precond $pc"
        dist(n; tank..., inner..., dcommon(m; name = n)..., precond = pc, solver_type = :theta)
    end
    # boundary generation + relaxation, y-periodic; regular wave and a WaveSpec sea
    for m in MODELS[[1, 3, 5]]
        n = "dist $tag model$(m.id) bcgen regular yperiodic"
        dist(n; tank..., bcgen..., dcommon(m; name = n)..., y_wall_bc = :periodic, solver_type = :theta)
    end
    let m = MODELS[1], n = "dist $tag model1 bcgen sea"
        dist(n; tank..., bcgen..., dcommon(m; name = n)..., wave_bc = warmup_sea_state(D0),
             solver_type = :theta)
    end
    let m = MODELS[1], n = "dist $tag model1 yopen"
        dist(n; tank..., inner..., dcommon(m; name = n)..., y_wall_bc = :open,
             sponge_wB = 0.5, sponge_wT = 0.5, solver_type = :theta)
    end
    # closed-basin release (x walls)
    for m in MODELS[[1, 5]]
        n = "dist $tag model$(m.id) release xwalls"
        dist(n; tank..., dcommon(m; name = n)..., A_wave = 0.0, x_wall_bc = true,
             sponge_wL = 0.0, sponge_wR = 0.0, eta0_func = eta_bump, solver_type = :theta)
    end
    # order-2 jump penalty with 𝓝 (:ghostvolume is sequential only)
    let m = MODELS[5], n = "dist $tag model5 jumpgrad order2"
        dist(n; tank..., inner..., dcommon(m; name = n)..., cip_order = 2, solver_type = :theta)
    end
end

# ── MMS drivers (verification; opt-in) ──────────────────────────────────────────────────────
function mms_cases(M, p)
    tag = "P$(p)LFE-$(M)"
    hbed = (x, y) -> 1.0 + 0.1*sin(2π*x)*cos(2π*y)
    for m in MODELS
        kw = (nx = 4, ny = 4, dt = 0.05, T_final = 0.1, M = M, p_vert = p,
              p_u = m.p_u, p_eta = m.p_u - 1, regime = m.regime, nl_pressure = m.nlp,
              flat_bed = m.flat, hfun = m.flat ? nothing : hbed, solver_type = :theta)
        trace("mms seq $tag model$(m.id)", run_mms_case, kw)
        trace("mms dist $tag model$(m.id)", run_mms_case_distributed, (; kw..., cpu_grid = (1, 1), verbose = false))
    end
end

# ── run ────────────────────────────────────────────────────────────────────────────────────
println("warmup: vertical bases ", join(("P$(p)LFE-$(M)" for (M, p) in VERT), ", "),
        " | extra groups: ", isempty(filter(!isempty, EXTRA)) ? "none" : join(filter(!isempty, EXTRA), ","),
        isempty(ONLY) ? "" : " | only: \"$ONLY\"")
for (M, p) in VERT
    sequential_cases(M, p)
    distributed_cases(M, p)
    "mms" in EXTRA && mms_cases(M, p)
end

# ── summary: what is baked and what will JIT ───────────────────────────────────────────────
nok = count(r -> r[2], RESULTS); nfail = length(RESULTS) - nok
println("\n", "="^78)
@printf("  warmup summary: %d traced, %d FAILED (failed cases are not baked; their runs will JIT)\n", nok, nfail)
println("="^78)
for (name, ok, t) in RESULTS
    @printf("  %-6s %-55s %8.1f s\n", ok ? "ok" : "FAILED", name, t)
end
println("="^78)
all(r -> !r[2], filter(r -> startswith(r[1], "dist"), RESULTS)) && any(r -> startswith(r[1], "dist"), RESULTS) &&
    @warn "warmup: EVERY distributed case failed — most likely MPI is unavailable in the build " *
          "subprocess. The image carries the sequential solver only; distributed runs will JIT."
DRY && nfail > 0 && error("warmup (dry run): $nfail case(s) pass keywords their driver does not accept")
println(DRY ? "warmup dry run done — the keyword sets are valid." : "warmup done.")
