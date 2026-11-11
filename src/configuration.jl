# ==============================================================
#  configuration.jl — CONFIGURATION.txt: the human-readable identity card of a run (2026-11-11)
#
#  Every run of `setup_and_run` / `setup_and_run_distributed` writes `<output_dir>/CONFIGURATION.txt`
#  before the time loop starts. It is the companion of `run_manifest.toml`, written for a reader:
#    * the manifest is the complete machine record (every keyword, `repr`'d, plus the git state);
#    * CONFIGURATION.txt is the standardised summary a person scans to tell runs apart. The same
#      sections, keys and order in every run, a one-line ID at the top (grep it across runs), and
#      the DERIVED magnitudes the keywords do not show: the model's own wavelength, kd, kA, the
#      dispersion errors against Airy, the cells per wavelength, the steps per period, the
#      integrator's damping of the carrier (rule 15: a dissipative run states its cost).
#  Values are computed here from the configuration (`Base.@locals()` of the driver) and the objects
#  the driver has built (vertical tensors, FE spaces); nothing here feeds back into the run.
# ==============================================================

"Model number 1–6 of the (regime, nl_pressure, flat_bed) matrix (CLAUDE.md §4)."
model_number(regime::Symbol, nl_pressure::Bool, flat_bed::Bool) =
    regime === :linear ? (flat_bed ? 1 : 2) : nl_pressure ? (flat_bed ? 5 : 6) : (flat_bed ? 3 : 4)

"""
    carrier_amplification(solver_type, y; theta, tableau, rho_inf) → |R(iy)|

Amplitude factor per step of the integrator on a neutral oscillation `u̇ = iωu`, `y = ωΔt`
(1 = no damping). θ-method: `R = (1 + (1−θ)iy)/(1 − θiy)`; Runge–Kutta: the tableau's stability
function; generalised-α: the physical eigenvalue of its two-step recurrence (the one closest to
`exp(iy)`).
"""
function carrier_amplification(solver_type::Symbol, y::Float64; theta::Float64 = 0.5,
                               tableau::Symbol = :SDIRK_2_2, rho_inf::Float64 = 0.5)
    z = im * y
    if solver_type === :theta
        return abs((1 + (1 - theta) * z) / (1 - theta * z))
    elseif solver_type === :sdirk || solver_type === :rk
        return abs(tableau_R(Gridap.ODEs.ButcherTableau(getfield(Gridap.ODEs, tableau)()), z))
    else
        #  GeneralizedAlpha1 (Gridap): u₁ = u₀ + Δt((1−γ)v₀ + γv₁),  v_αm = λ u_αf  (λΔt = z)
        αf = 1 / (1 + rho_inf); αm = (3 - rho_inf) / (1 + rho_inf) / 2; γ = 1 / 2 + αm - αf
        L  = [1.0+0im -γ; -z*αf αm]                      # unknowns (u₁, Δt·v₁)
        Rr = [1.0+0im (1 - γ); z*(1 - αf) -(1 - αm)]     # knowns   (u₀, Δt·v₀)
        ev = eigvals(L \ Rr)
        return abs(ev[argmin(abs.(ev .- exp(z)))])
    end
end

"""
    integrator_lines(solver_type; theta, tableau, rho_inf) → Vector{Pair{String,String}}

The integrator for the configuration file, one property per line: the scheme, its parameter,
its kind, stages, order and stiff damping `|R(∞)|` (the same keys for every integrator).
"""
function integrator_lines(solver_type::Symbol; theta, tableau, rho_inf)
    if solver_type === :theta
        return ["time integrator" => "ThetaMethod (θ-method)",
                "integrator parameter" => @sprintf("θ = %.4g%s", theta, theta == 0.5 ? " (Crank–Nicolson)" : ""),
                "integrator kind" => "implicit",
                "integrator stages" => "1",
                "integrator order" => theta == 0.5 ? "2" : "1",
                "dissipative" => theta == 0.5 ? "no" : "yes",
                "|R(∞)|" => @sprintf("%.4g", abs((1 - theta) / theta))]
    elseif solver_type === :gen_alpha
        return ["time integrator" => "GeneralizedAlpha1 (generalised-α, first-order form)",
                "integrator parameter" => @sprintf("ρ∞ = %.4g", rho_inf),
                "integrator kind" => "implicit, carries (u, u̇)",
                "integrator stages" => "1",
                "integrator order" => "2",
                "dissipative" => rho_inf == 1 ? "no" : "yes",
                "|R(∞)|" => @sprintf("%.4g", rho_inf)]
    end
    tb   = Gridap.ODEs.ButcherTableau(getfield(Gridap.ODEs, tableau)())
    expl = Gridap.ODEs.TableauType(tb) == Gridap.ODEs.ExplicitTableau
    damp = 1 - abs(tableau_R(tb, 0.3im))
    return ["time integrator" => "RungeKutta" * (tableau in DEFECTIVE_TABLEAUS ? "  ⚠ defective Gridap tableau (rule 15)" : ""),
            "integrator parameter" => @sprintf("tableau %s", tableau),
            "integrator kind" => expl ? "explicit" : "diagonally implicit",
            "integrator stages" => string(length(Gridap.ODEs.get_weights(tb))),
            "integrator order" => @sprintf("%d (Gridap label)", Gridap.ODEs.get_order(tb)),
            "dissipative" => abs(damp) < 1e-12 ? "no" : damp > 0 ? "yes" : "no — AMPLIFIES oscillations",
            "|R(∞)|" => expl ? "— (explicit)" : @sprintf("%.4g", abs(tableau_R(tb, -1e12 + 0im)))]
end

"""
    try_write_run_configuration(output_dir, cfg; kw...)

`write_run_configuration` that never stops a run: a failure is reported as a warning (the
configuration file is a report, not part of the computation). Returns the path or `nothing`.
"""
function try_write_run_configuration(output_dir, cfg; kw...)
    try
        return write_run_configuration(output_dir, cfg; kw...)
    catch e
        @warn "CONFIGURATION.txt could not be written (the run continues)" exception = (e, catch_backtrace())
        return nothing
    end
end

"""
    write_run_configuration(output_dir, cfg; vert, c_bdy, n_free, git="", driver, nranks=1,
                            solver_desc=("", ""), t_start=0.0)

Write `<output_dir>/CONFIGURATION.txt` (see the file header). `cfg` is the driver's keyword set
(`Base.@locals()`, with the resolved CIP γ), `vert` the vertical tensors, `c_bdy` the σ-element
boundaries in use, `n_free` the free DOF count, `solver_desc` = (nonlinear, linear) solver
descriptions as the driver's banner prints them. Returns the path.
"""
function write_run_configuration(output_dir::AbstractString, cfg::AbstractDict; vert, c_bdy,
                                 n_free::Integer, git::AbstractString = "",
                                 driver::AbstractString = "setup_and_run", nranks::Integer = 1,
                                 solver_desc = ("", ""), t_start::Float64 = 0.0)
    c(k, default = nothing) = get(cfg, k, default)
    g = Float64(c(:g, 9.81)); d = Float64(c(:h_val, 3.5))
    regime = c(:regime); nlp = c(:nl_pressure, false); flat = c(:flat_bed, true)
    M = c(:M); p = c(:p_vertical, 1); Nσ = vert.N_dof
    p_u = c(:p_u); p_eta = c(:p_eta)
    (x0, x1), (y0, y1) = c(:domain); nx, ny = c(:partition)
    dx = (x1 - x0) / nx; dy = (y1 - y0) / ny
    T = Float64(c(:T_wave)); A = Float64(c(:A_wave)); ω = 2π / T
    dt = Float64(c(:dt)); T_final = Float64(c(:T_final))
    st = c(:solver_type); θ = Float64(c(:theta, 0.5)); tab = c(:tableau, :SDIRK_2_2)
    ρ = Float64(c(:rho_inf, 0.5))
    mid = model_number(regime, nlp, flat)

    # ---- derived wave magnitudes (carrier at T_wave, A_wave, reference depth h_val) ------------
    k_e = find_wavenumber(ω, d, g); λ_e = 2π / k_e
    k_m = model_wavenumber(vert, ω, d, g); λ_m = 2π / k_m
    kd_e = k_e * d; kd_m = k_m * d
    wp = wave_properties(vert, kd_e; g = g, d = d)
    err_C = wp.C / wp.Ce - 1; err_Cg = wp.Cg / wp.Cge - 1
    ar_C = applicable_range(vert; prop = :C, g = g); ar_Cg = applicable_range(vert; prop = :Cg, g = g)
    ω_inf = sqrt(g / d * dot(vert.Phi, (-vert.B) \ vert.Phi))
    H = 2A; miche = 0.142 * tanh(kd_m)
    dfn = c(:h_bathy) === nothing ? (x -> d) : c(:h_bathy)
    dsamp = [dfn(VectorValue(x, y)) for x in range(x0, x1; length = 41), y in range(y0, y1; length = 3)]

    # ---- time integration magnitudes --------------------------------------------------------------
    nsteps = round(Int, (T_final - t_start) / dt)
    y = ω * dt
    amp = carrier_amplification(st, y; theta = θ, tableau = tab, rho_inf = ρ)
    loss_run = 1 - amp^nsteps

    # ---- generation ----------------------------------------------------------------------------
    wave_bc = c(:wave_bc); x_per = c(:x_periodic, false)
    regular = wave_bc === nothing || wave_bc === :regular
    gen = if x_per
        "none (closed x-periodic domain; initial condition only)"
    elseif wave_bc !== nothing
        "Dirichlet boundary generation (" * (wave_bc isa Symbol ? string(wave_bc) : string(nameof(typeof(wave_bc)))) * ")"
    elseif c(:eta0_func) !== nothing && A == 0
        "none (initial-condition release)"
    elseif c(:y_wm) === nothing
        "interior line source (plane waves)"
    else
        "interior point source (ring waves)"
    end
    yb = c(:y_wall_bc, :wall) === :periodic ? "periodic" : c(:y_wall_bc, :wall) === :wall ? "solid walls" : "open (natural)"
    xb = x_per ? "periodic" : c(:x_wall_bc, false) ? "solid walls" : "open (natural)"
    γu = Float64(c(:cip_gamma_u, 0.0)); γη = Float64(c(:cip_gamma_eta, 0.0))
    itag = st === :theta ? @sprintf("theta(%.3g)", θ) : st === :gen_alpha ? @sprintf("gen_alpha(%.3g)", ρ) : string(tab)

    # ---- the one-line ID: what tells two runs apart at a glance ---------------------------------
    id = @sprintf("m%d | P%dLFE-%d | Q%d/Q%d | %d×%d cells | A=%.4g m T=%.4g s d=%.4g m | %s dt=%.4g s | t=%.4g–%.4g s | %s",
                  mid, p, M, p_u, p_eta, nx, ny, A, T, d, itag, dt, t_start, T_final,
                  nranks > 1 ? "$(nranks) ranks" : "sequential")

    #  one quantity per line; every run gets the same keys in the same order ("—" = not applicable)
    L = String[]
    sec(name) = (push!(L, ""); push!(L, "[" * name * "]"))
    kv(k, v) = push!(L, rpad(k, 34) * ": " * string(v))
    f(x; u = "") = @sprintf("%.5g", x) * (isempty(u) ? "" : " " * u)
    pct(x) = @sprintf("%+.3g %%", 100x)
    push!(L, "="^100)
    push!(L, "GridapBALFEM run configuration — CONFIGURATION.txt")
    push!(L, "="^100)
    kv("ID", id)
    kv("run", basename(normpath(output_dir)))
    kv("written", Libc.strftime("%Y-%m-%d %H:%M:%S", time()))
    kv("code", isempty(git) ? "(git state unknown)" : git)
    kv("driver", driver)

    sec("MODEL")
    kv("model", mid)
    kv("regime", regime)
    kv("nonlinear pressure 𝓝", nlp ? "on (all eight components)" : "off")
    kv("sea bed", flat ? "flat" : "variable")
    kv("vertical basis", @sprintf("P%dLFE-%d", p, M))
    kv("vertical elements M", M)
    kv("vertical order p", p)
    kv("vertical nodes Nσ", Nσ)
    kv("σ-element boundaries", "[" * join((@sprintf("%.4g", v) for v in c_bdy), ", ") * "]")
    kv("applicable kd, phase speed (2 %)", @sprintf("≤ %.4g", ar_C))
    kv("applicable kd, group speed (2 %)", @sprintf("≤ %.4g", ar_Cg))
    kv("Class III treatment", nlp ? "cellwise Hessians + kind-S facet terms" : "—")
    kv("skeleton stabilisation", (γu > 0 || γη > 0) ? string(c(:stabilization, :jumpgrad)) : "none")
    kv("CIP γ_u (velocity)", f(γu))
    kv("CIP γ_η (surface)", f(γη))
    kv("CIP order", (γu > 0 || γη > 0) ? c(:cip_order, 1) : "—")
    kv("CIP h exponent", (γu > 0 || γη > 0) ? f(Float64(c(:cip_hexp, 2.0))) : "—")

    sec("PHYSICS")
    kv("gravity g", f(g; u = "m/s²"))
    kv("still-water depth d", f(d; u = "m"))
    kv("bed depth range", flat ? "—" : @sprintf("%.5g – %.5g m", minimum(dsamp), maximum(dsamp)))
    kv("wave generation", gen)
    kv("wave input", regular ? "regular wave (A, T below)" : "sea state (A, T below are nominal)")
    kv("amplitude A", f(A; u = "m"))
    kv("wave height H = 2A", f(H; u = "m"))
    kv("period T", f(T; u = "s"))
    kv("angular frequency ω", f(ω; u = "rad/s"))
    kv("direction", f(rad2deg(Float64(c(:wave_dir, 0.0))); u = "°"))
    kv("wavelength λ (model)", f(λ_m; u = "m"))
    kv("wavelength λ (Airy)", f(λ_e; u = "m"))
    kv("wavenumber k (model)", f(k_m; u = "1/m"))
    kv("wavenumber k (Airy)", f(k_e; u = "1/m"))
    kv("kd (model)", f(kd_m))
    kv("kd (Airy)", f(kd_e) * (kd_e <= ar_Cg ? "" : kd_e <= ar_C ? "  ⚠ beyond the group-speed range" :
                                                                   "  ⚠ beyond both applicable ranges"))
    kv("steepness kA (model)", f(k_m * A))
    kv("steepness H/λ (model)", f(H / λ_m))
    kv("Miche limit H/λ = 0.142 tanh kd", f(miche))
    kv("H/λ relative to Miche", @sprintf("%.3g %%", 100 * (H / λ_m) / miche))
    kv("relative amplitude A/d", f(A / d))
    kv("relative height H/d", f(H / d))
    kv("Ursell number Hλ²/d³", f(H * λ_m^2 / d^3))
    kv("phase speed C (model)", f(wp.C; u = "m/s"))
    kv("phase-speed error vs Airy", pct(err_C))
    kv("group speed Cg (model)", f(wp.Cg; u = "m/s"))
    kv("group-speed error vs Airy", pct(err_Cg))
    kv("frequency ceiling ω∞ (model)", f(ω_inf; u = "rad/s"))
    kv("ω/ω∞", f(ω / ω_inf))
    kv("2ω/ω∞ (second harmonic)", f(2ω / ω_inf))

    sec("DOMAIN AND DISCRETISATION")
    kv("domain x", @sprintf("%.5g – %.5g m", x0, x1))
    kv("domain y", @sprintf("%.5g – %.5g m", y0, y1))
    kv("domain length / λ (model)", f((x1 - x0) / λ_m))
    kv("cells x × y", "$nx × $ny")
    kv("cell size dx", f(dx; u = "m"))
    kv("cell size dy", f(dy; u = "m"))
    kv("cells per wavelength λ/dx", f(λ_m / dx))
    kv("x boundaries", xb)
    kv("y boundaries", yb)
    kv("sponge width left", f(Float64(c(:sponge_wL, 0.0)); u = "m"))
    kv("sponge width right", f(Float64(c(:sponge_wR, 0.0)); u = "m"))
    kv("sponge width bottom", f(Float64(c(:sponge_wB, 0.0)); u = "m"))
    kv("sponge width top", f(Float64(c(:sponge_wT, 0.0)); u = "m"))
    kv("sponge strength μ_max", f(Float64(c(:mu_max, 0.0)); u = "1/s"))
    kv("source position", (wave_bc === nothing && !x_per && !(c(:eta0_func) !== nothing && A == 0)) ?
       (c(:y_wm) === nothing ? @sprintf("x = %.5g m", Float64(c(:x_wm, 0.0))) :
                               @sprintf("(x, y) = (%.5g, %.5g) m", Float64(c(:x_wm, 0.0)), Float64(c(:y_wm)))) : "—")
    kv("generation boundary", wave_bc === nothing ? "—" : string(c(:bc_side, :left)))
    kv("generation vertical profile", wave_bc === nothing ? "—" : string(c(:bc_profile, :model)))
    kv("relaxation zone width", (wave_bc !== nothing && c(:relax_bc, false)) ? f(Float64(c(:relax_width, 0.0)); u = "m") : "—")
    kv("FE pairing", "Taylor–Hood (Lagrange, C⁰)")
    kv("surface field η", "scalar, Q$(p_eta)")
    kv("velocity fields 𝖴x, 𝖴y", "VectorValue{$(Nσ)}, Q$(p_u)")
    kv("free DOFs", n_free)
    kv("quadrature degree", 2 * max(p_u, p_eta) + 2 + Int(c(:quad_extra, 0)))

    sec("TIME INTEGRATION")
    foreach(((k, v),) -> kv(k, v), integrator_lines(st; theta = θ, tableau = tab, rho_inf = ρ))
    kv("time step dt", f(dt; u = "s"))
    kv("start time", f(t_start; u = "s"))
    kv("end time", f(T_final; u = "s"))
    kv("steps", nsteps)
    kv("steps per period T/dt", f(T / dt))
    kv("ωΔt", f(y))
    kv("Courant number C·dt/dx", f(wp.C * dt / dx))
    kv("carrier damping per step", @sprintf("%.3e (1 − |R(iωΔt)|)", 1 - amp))
    kv("carrier amplitude lost over run", @sprintf("%.4g %%", 100 * loss_run))
    kv("nonlinear solver", isempty(solver_desc[1]) ? "Newton" : solver_desc[1])
    kv("linear solver", isempty(solver_desc[2]) ? "—" : solver_desc[2])
    kv("MPI ranks", nranks > 1 ? nranks : "1 (sequential)")
    kv("cpu grid", nranks > 1 ? string(c(:cpu_grid)) : "—")
    kv("restart from", c(:restart_from) === nothing ? "none (cold start)" : string(c(:restart_from)))

    sec("OUTPUT")
    kv("output directory", normpath(output_dir))
    kv("VTK snapshots", Int(c(:save_every, 0)) > 0 ? "every $(c(:save_every)) steps" : "off")
    kv("diagnostics", Int(c(:diag_every, 0)) > 0 ? "every $(c(:diag_every)) steps" : "off")
    kv("checkpoints", Int(c(:checkpoint_every, 0)) > 0 ? "every $(c(:checkpoint_every)) steps" : "off")
    kv("checkpoints kept", Int(c(:checkpoint_every, 0)) > 0 ? c(:checkpoint_keep, 2) : "—")
    kv("machine record", "run_manifest.toml (every keyword, git state)")
    push!(L, "="^100)

    mkpath(output_dir)
    path = joinpath(output_dir, "CONFIGURATION.txt")
    open(io -> foreach(l -> println(io, l), L), path, "w")
    return path
end
