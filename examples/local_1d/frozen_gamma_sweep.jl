# ==============================================================
#  frozen_gamma_sweep.jl — the frozen-state spectrum of the full model against the CIP penalty γ,
#  and what each time integrator does with it (pre-cluster check 2, 2026-11-12)
#
#  THE QUESTION. For model 5 (𝓝 + mandatory C⁰-CIP) in the closed box at n cells/λ: how does the
#  fastest frozen-state growth σ_max change with γ and with the carrier steepness kA, at which
#  wavenumber and frequency does it sit, and which (integrator, Δt) pairs damp it? The answer
#  sets the γ values of tier G of run/integrator_study (cases.py GAMMAS_G) around the predicted
#  boundary, before the cluster time is spent.
#
#  THE METHOD. stability_eig.jl's frozen-state linearisation, M* δu̇ + J* δu = 0 at the linear
#  eigenmode of amplitude A (Jacobians by coloured central FD of the assembled residual), with one
#  shortcut: the C⁰-CIP is LINEAR in u and in γ and has no u̇ term, so
#       J(γ) = J(γa) + (γ − γa)/(γb − γa) · (J(γb) − J(γa)),   M independent of γ.
#  Two FD Jacobians per kA give every γ; the script gates the shortcut at one direct γ.
#  Per (kA, γ): σ_max with its |ω| and dominant harmonic n, the growing modes, and the
#  scheme-and-Δt verdict  ρ_s(Δt) = max_λ ln|R_s(λΔt)| / Δt  (R_s the stability function of the
#  Gridap tableau, R = (1 + z/2)/(1 − z/2) for CN): > 0 means that integrator at that Δt does not
#  damp the frozen growth.
#
#  ⚠ WHAT IT CANNOT SEE. A frozen state is one instant of the wave. It measures instantaneous,
#  mostly high-k growth (the tool's design, stability_eig.jl header) and the damping budgets of the
#  CIP and of the integrators. It does NOT capture a resonant interaction pumped by the carrier's
#  time dependence (the m3 decay triad: pump at 2ω₀, a steady member and a 2ω₀ member). Those are
#  measured by the box runs (temporal spectrum), not here. Read σ_max as a screening of γ, never
#  as a stability verdict of the time-dependent run.
#
#  RUN:   julia --project=. examples/local_1d/frozen_gamma_sweep.jl
#  ENV:   FGS_NCELL (32), FGS_KA ("0.158,0.25,0.32"), FGS_GAMMAS ("0.1,0.2,0.3,0.5,1,2,3,5,10"),
#         FGS_MODEL (on | off), FGS_OUT (output/precluster/frozen_gamma_sweep)
#  OUT:   <FGS_OUT>/frozen_gamma_sweep_n<ncell>.csv and the log on stdout
# ==============================================================
include(joinpath(@__DIR__, "stability_eig.jl"))

const FGS_NCELL  = parse(Int, get(ENV, "FGS_NCELL", "32"))
const FGS_KA     = parse.(Float64, split(get(ENV, "FGS_KA", "0.158,0.25,0.32"), ","))
const FGS_GAMMAS = parse.(Float64, split(get(ENV, "FGS_GAMMAS", "0.1,0.2,0.3,0.5,1,2,3,5,10"), ","))
const FGS_MODEL  = Symbol(get(ENV, "FGS_MODEL", "on"))
const FGS_OUT    = get(ENV, "FGS_OUT", joinpath(@__DIR__, "..", "..", "output", "precluster", "frozen_gamma_sweep"))
const FGS_GA, FGS_GB, FGS_GGATE = 0.3, 3.0, 1.0          # the two FD γ and the gate γ
const FGS_SCHEMES = (:CN, :SDIRK_2_2, :SDIRK_4_3, :SDIRK_Crouzeix_3_4)
const FGS_SPP     = (80, 40, 20, 10)                       # Δt = T/spp

"""Frozen-state matrices of the box at amplitude A and penalty γ: the body of `stability_case`
up to J and M, returned with what the mode tagging needs. `want_M = false` skips M."""
function fgs_ops(; pu = 3, pe = pu - 1, model = FGS_MODEL, A, ncell, gamma, vert = se_vert(), want_M = true)
    Nσ = vert.N_dof
    ω  = 2π / SE_T
    k  = model_wavenumber(vert, ω, SE_D, SE_G)
    L  = 2π / k
    dx = L / ncell
    uamp = _component_uamp(vert, Float64[], :model, 1.0, ω, k, SE_D, SE_G)
    dm    = CartesianDiscreteModel((0.0, L, 0.0, dx), (ncell, 1); isperiodic = (true, false))
    trian = Triangulation(dm)
    qdeg  = 2 * max(pu, pe) + 2
    dΩh   = Measure(trian, qdeg)
    U, V  = build_fe_spaces(dm, pu, Nσ; y_wall_bc = :wall, x_wall_bc = false, p_eta = pe)
    prob  = build_problem(vert; h_bathy = (x -> SE_D), regime = :nonlinear,
                                nl_pressure = (model === :on), flat_bed = true,
                                model = dm, quad_degree = qdeg, p_eta = pe)
    gamma > 0 && attach_skeleton!(prob, dm; cip_gamma_u = gamma, cip_gamma_eta = gamma,
                                  stabilization = :jumpgrad, cip_order = 1,
                                  p_u = pu, p_eta = pe, degree = qdeg)
    zv   = VectorValue(ntuple(_ -> 0.0, Nσ)...)
    uh   = interpolate_everywhere([x -> A * cos(k * x[1]),
                                   x -> VectorValue(ntuple(j -> A * uamp[j] * cos(k * x[1]), Nσ)...),
                                   x -> zv], U)
    uth  = interpolate_everywhere([x -> A * ω * sin(k * x[1]),
                                   x -> VectorValue(ntuple(j -> A * uamp[j] * ω * sin(k * x[1]), Nσ)...),
                                   x -> zv], U)
    res  = (t, u, v) -> global_residual(t, u, v, prob, trian, dΩh)
    rvec(x, xt) = assemble_vector(v -> res(0.0, TransientCellField(FEFunction(U, x),
                                                                    (FEFunction(U, xt),)), v), V)
    x0  = copy(get_free_dof_values(uh)); xt0 = copy(get_free_dof_values(uth))
    pat = se_pattern(U, V, dΩh, Nσ, prob.skel[])
    col = se_colour(pat)
    J   = Matrix(se_fdjac(x -> rvec(x, xt0), x0, pat, col))
    M   = want_M ? Matrix(se_fdjac(w -> rvec(x0, w), xt0, pat, col)) : nothing
    nf   = [num_free_dofs(U.spaces[i]) for i in 1:length(U.spaces)]
    offs = cumsum([0; nf])
    return (J = J, M = M, U = U, trian = trian, pu = pu, ncell = ncell, L = L, k = k, ω = ω,
            Nσ = Nσ, rng2 = (offs[2] + 1):offs[3])
end

"Eigenpairs of −M⁻¹J with the dominant box harmonic n of each mode's 𝖴x field."
function fgs_spectrum(o, J, M, Fmap)
    E = eigen(-(M \ J))
    Bj, Fk, kgrid = Fmap
    Vux = E.vectors[o.rng2, :]
    Pw  = zeros(length(kgrid), length(E.values))
    for j in 1:o.Nσ
        Pw .+= abs2.(Fk * (Bj[j] * Vux))
    end
    ndom = [argmax(view(Pw, :, m)) - 1 for m in 1:length(E.values)]     # harmonic index n = k/k0
    return E.values, ndom
end

Rcn(z) = (1 + z / 2) / (1 - z / 2)
fgs_R(s) = s === :CN ? Rcn : (tb = Gridap.ODEs.ButcherTableau(s); z -> GridapBALFEM.tableau_R(tb, z))

"ρ_s(Δt) = max_λ ln|R_s(λΔt)|/Δt: the fastest growth that scheme s at Δt lets through (1/s)."
fgs_rate(R, λ, dt) = maximum(l -> log(abs(R(l * dt))), λ) / dt

function main()
    mkpath(FGS_OUT)
    csv = joinpath(FGS_OUT, "frozen_gamma_sweep_n$(FGS_NCELL)_$(FGS_MODEL).csv")
    Rs  = Dict(s => fgs_R(s) for s in FGS_SCHEMES)
    hdr = ["kA", "A", "gamma", "sigma_max", "omega_at_max", "n_at_max", "n_growing",
           "cip_carrier_damping"]
    for s in FGS_SCHEMES, spp in FGS_SPP; push!(hdr, "rate_$(s)_T$(spp)"); end
    open(csv, "w") do io; println(io, join(hdr, ",")); end
    @printf("frozen γ sweep: model %s, n%d, kA %s, γ %s\n", FGS_MODEL, FGS_NCELL,
            join(FGS_KA, " "), join(FGS_GAMMAS, " ")); flush(stdout)
    vert = se_vert()
    k0 = model_wavenumber(vert, 2π / SE_T, SE_D, SE_G)
    for (ika, kA) in enumerate(FGS_KA)
        A = kA / k0
        t0 = time()
        oa = fgs_ops(; A = A, ncell = FGS_NCELL, gamma = FGS_GA, vert = vert)
        ob = fgs_ops(; A = A, ncell = FGS_NCELL, gamma = FGS_GB, vert = vert, want_M = false)
        Fmap = se_fourier_map(oa.U.spaces[2], oa.trian, oa.pu, oa.ncell, oa.L, oa.Nσ)
        Jg(γ) = oa.J .+ ((γ - FGS_GA) / (FGS_GB - FGS_GA)) .* (ob.J .- oa.J)
        @printf("\n=== kA = %.3f (A = %.4f m): 2 FD Jacobians + M in %.0f s\n", kA, A, time() - t0)
        if ika == 1      # gate the linear-in-γ shortcut against a direct FD Jacobian
            oc = fgs_ops(; A = A, ncell = FGS_NCELL, gamma = FGS_GGATE, vert = vert, want_M = false)
            @printf("  gate: ‖J(γ=%g) direct − interpolated‖/‖J‖ = %.2e (FD noise ~1e-8)\n", FGS_GGATE,
                    norm(oc.J .- Jg(FGS_GGATE)) / norm(oc.J))
        end
        @printf("  %-6s %-11s %-8s %-4s %-5s %-10s | %s\n", "γ", "σmax /s", "|ω|", "n", "#grow",
                "CIP@carrier", "fastest growth let through, max ln|R|/Δt (1/s), Δt = T/" * join(FGS_SPP, ", T/"))
        for γ in FGS_GAMMAS
            λ, nd = fgs_spectrum(oa, Jg(γ), oa.M, Fmap)
            σ = real.(λ); im_ = argmax(σ)
            ngrow = count(>(1e-6), σ)
            # the CIP's damping of the carrier: the n = 1 modes nearest ω0 (frozen, one phase)
            c1 = findall(i -> nd[i] == 1, eachindex(λ))
            cipc = isempty(c1) ? NaN : -real(λ[c1[argmin(abs.(abs.(imag.(λ[c1])) .- oa.ω))]])
            rates = [fgs_rate(Rs[s], λ, SE_T / spp) for s in FGS_SCHEMES for spp in FGS_SPP]
            @printf("  %-6g %+.3e  %-8.3g %-4d %-5d %+.2e  |", γ, σ[im_], abs(imag(λ[im_])), nd[im_], ngrow, cipc)
            for (is, s) in enumerate(FGS_SCHEMES)
                @printf(" %s: %s", s, join([@sprintf("%+.2e", r) for r in rates[(is - 1) * length(FGS_SPP) .+ (1:length(FGS_SPP))]], " "))
            end
            println()
            # the growing modes, fastest first (at most 5)
            gi = sort(findall(>(1e-6), σ); by = i -> -σ[i])[1:min(5, ngrow)]
            isempty(gi) || println("         growing: ", join([@sprintf("σ %.3e |ω| %.3g n %d", σ[i], abs(imag(λ[i])), nd[i]) for i in gi], "; "))
            open(csv, "a") do io
                println(io, join(string.([kA, A, γ, σ[im_], abs(imag(λ[im_])), nd[im_], ngrow, cipc, rates...]), ","))
            end
            flush(stdout)
        end
    end
    println("\nwritten: ", csv)
end

main()
