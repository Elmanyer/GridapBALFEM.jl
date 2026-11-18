# ==============================================================
#  oss_cw_spectrum.jl — why does C_w make m5 mixed + OSS grow?  (branch mixed-vms-solver)
#
#  Box of oss_newton_trace.jl (m5 flat bed, 1 λ, 16 cells/λ, Q3/Q2, P1LFE-2, kA = 0.08, CN Δt = 0.04).
#  At the hot-start eigenmode scaled by a, it assembles A = ∂R/∂u (hand core + exact jet VMS part) and
#  B = ∂R/∂u̇ (Galerkin mass: OSS adds none), eliminates the algebraic fields (𝖦, 𝖥, Π, δ) by Schur
#  complement and reports the frozen spectrum of B u̇ = −A u: max Re λ, and for the most unstable mode
#  its character (vertical-structure share of the kinetic energy, η share, effective wavenumber
#  k_eff·h_e/π: 1 = the node-spacing Nyquist).
#
#  HYPOTHESIS UNDER TEST (H1). The stabilisation is ∫ L*(v)·τ P⊥(Rˢ). Its test operator L*(v) carries
#  only the transport parts (advection N_b∂_b W, PSPG Φ∂q), but Rˢ also carries 𝓝, whose strong form
#  has third derivatives of 𝗎 (∂ₐ(H² NP), NP ∋ 𝖦 ≈ ∇∇·(H𝗎)). The cross term τ_w ∇q·P⊥(e·𝓝) has no
#  sign and is ~ (U/c₀)(k h_e⁻¹ H)² larger than the positive part at the grid scale.
#  Discriminator: the variant `_noN` evaluates the SAME stabilisation with 𝓝 left out of Rˢ (a copy
#  of the problem with nl_pressure = false, used ONLY inside vms_contrib / oss_consistent_ic).
#  H1 predicts: wave1 unstable at a > 0 and stable at a = 0; wave1_noN stable; rate ∝ C_w and ∝ a.
#
#  ⚠ The core ∂R/∂u is the hand Jacobian (quasi-Newton in the 𝓝 blocks, rule 5): absolute rates are
#  indicative; the comparison between variants (same core) is the measurement. Frozen spectra are
#  not predictive of long runs (STABILITY.md §2).
#  Output: output/mixed_vms/oss_spectrum/report.txt
# ==============================================================
using GridapBALFEM, Gridap, Gridap.ODEs, LinearAlgebra, Printf

const ROOT_S = normpath(joinpath(@__DIR__, "..", ".."))
outdir = joinpath(ROOT_S, "output", "mixed_vms", "oss_spectrum"); mkpath(outdir)
M, p_vert, p_u, p_eta, d, Twave, kA, ncell, dt, θ = 2, 1, 3, 2, 3.5, 1.6, 0.08, 16, 0.04, 0.5
vert  = assemble_vertical_tensors(M, p_vert, Vector{Float64}(resolve_cbdy(M, nothing, p_vert)))
Nσ    = vert.N_dof
k0    = GridapBALFEM.model_wavenumber(vert, 2π / Twave, d, 9.81)
Awave = kA / k0
wi    = WaveInput(vert; A = Awave, T = Twave, d = d, T_ramp = 0.0, profile = :model)
Lx    = 2π / wi.ks[1]; dx = Lx / ncell
inc   = incident_fields(wi)
QDEG  = 2 * max(p_u, p_eta) + 2
he    = dx / p_u

model, trian = build_horizontal_model(((0.0, Lx), (0.0, dx)), (ncell, 1); x_periodic = true)
dΩ   = Measure(trian, QDEG)
prob = build_problem(vert; h_bathy = x -> d, regime = :nonlinear, nl_pressure = true, flat_bed = true,
                     model = model, quad_degree = QDEG, p_eta = p_eta, formulation = :mixed)
#  the same problem with 𝓝 off, SHARING the mixed/vms Refs: used only inside the stabilisation
fs    = Any[getfield(prob, f) for f in fieldnames(typeof(prob))]
fs[findfirst(==(:nl_pressure), fieldnames(typeof(prob)))] = false
probN = typeof(prob).name.wrapper(fs...)

U7, V7   = build_fe_spaces(model, p_u, Nσ; y_wall_bc = :wall, p_eta = p_eta, n_aux = 4, p_aux = p_u, n_oss = 0)
U10, V10 = build_fe_spaces(model, p_u, Nσ; y_wall_bc = :wall, p_eta = p_eta, n_aux = 4, p_aux = p_u, n_oss = 3)
nphys = sum(num_free_dofs(U7.spaces[i]) for i in 1:3)
nη    = num_free_dofs(U7.spaces[1])
@printf("DOFs: mixed %d, mixed+OSS %d, physics %d (η %d)   h_e = %.4f m, τ_w(C=1) = %.4e s\n",
        num_free_dofs(U7), num_free_dofs(U10), nphys, nη, he, he / (2 * sqrt(9.81 * d)))
flush(stdout)

"Hot-start state ×a on layout U, consistent auxiliary (and OSS) fields; exact u̇ of η, 𝖴x."
function box_state(a, U, probS; ε = 1e-5)
    pad = Any[GridapBALFEM.extra_field_zeros(prob)...]
    uh  = interpolate_everywhere(Any[x -> a * inc.eta(x, 0.0), x -> a * inc.ux(x, 0.0),
                                     x -> zero(inc.ux(x, 0.0)), pad...], U)
    ut  = interpolate_everywhere(Any[x -> a * (inc.eta(x, ε) - inc.eta(x, -ε)) / 2ε,
                                     x -> a * (inc.ux(x, ε) - inc.ux(x, -ε)) / 2ε,
                                     x -> zero(inc.ux(x, 0.0)), pad...], U)
    uh = mixed_consistent_ic(uh, prob, U, trian, dΩ)[1]
    has_oss(prob) && (uh = oss_consistent_ic(uh, probS, U, trian, dΩ)[1])
    return uh, ut
end

"A (hand core + exact VMS of probS) and B, dense."
function mats(uh, ut, U, V, probS)
    tu  = TransientCellField(uh, (ut,))
    rec = prob.vms[]
    prob.vms[] = nothing
    A = Matrix(assemble_matrix((du, v) -> jacobian_u(0.0, tu, du, v, prob, trian, dΩ), U, V))
    B = Matrix(assemble_matrix((dw, v) -> jacobian_u_t(0.0, tu, dw, v, prob, trian, dΩ), U, V))
    prob.vms[] = rec
    if rec !== nothing
        A .+= Matrix(assemble_matrix((du, v) -> vms_contrib(probS, 0.0, tu, v, trian, dΩ; du = du), U, V))
    end
    return A, B
end

"Mode character of the full vector y (physics + algebraic)."
function character(y, U)
    yh  = FEFunction(U, y)
    η, ux = yh[1], yh[2]
    Mv, Φ = prob.Mv, prob.Φ
    mq  = Operation(w -> w ⋅ (Mv ⋅ w))            # M-weighted square
    bq  = Operation(w -> (w ⋅ Φ)^2)               # depth-average square (M-energy of e ū)
    ek  = sum(∫(d * mq(ux)) * dΩ)
    ebt = sum(∫(d * bq(ux)) * dΩ)
    ep  = sum(∫(9.81 * η * η) * dΩ)
    dxu = alg_dx(ux)
    k2  = sum(∫(mq(dxu)) * dΩ) / max(sum(∫(mq(ux)) * dΩ), 1e-300)
    return (vs = (ek - ebt) / max(ek, 1e-300), eta = ep / max(ep + ek, 1e-300), kh = sqrt(max(k2, 0.0)) * he / π)
end

function spectrum(A, B, U)
    p = 1:nphys; q = (nphys + 1):size(A, 1)
    X  = A[q, q] \ A[q, p]
    Ar = A[p, p] - A[p, q] * X
    F  = eigen(-Ar, B[p, p])
    ok = findall(isfinite, F.values)
    λ  = F.values[ok]
    i  = ok[argmax(real.(λ))]
    yp = F.vectors[:, i]; yp ./= norm(yp)
    y  = vcat(yp, -X * yp)
    cr = character(real.(y), U); ci = character(imag.(y), U)
    c  = norm(real.(yp)) ≥ norm(imag.(yp)) ? cr : ci
    nunst = count(l -> real(l) > 1e-6, λ)
    amp = maximum(abs, (1 .+ (1 - θ) * dt .* λ) ./ (1 .- θ * dt .* λ))
    return (re = real(F.values[i]), im = imag(F.values[i]), nunst = nunst, maxabs = maximum(abs, λ),
            amp = amp, c...)
end

#  (tag, coefficients (C_a, C_w, C_η) or nothing, 𝓝 in the stabilised residual)
variants = [("none", nothing, true),
            ("wave1", (0.0, 1.0, 0.0), true), ("wave1_noN", (0.0, 1.0, 0.0), false),
            ("wave0.1", (0.0, 0.1, 0.0), true),
            ("adv_eta", (1.0, 0.0, 1.0), true), ("all1_noN", (1.0, 1.0, 1.0), false)]
AMPS = parse.(Float64, split(get(ENV, "BALFEM_SPEC_AMPS", "0,0.5,1,2"), ","))
rows = String[]
hdr  = @sprintf("%-10s %5s %11s %10s %6s %10s %8s | %6s %6s %6s %8s", "variant", "a", "max Re λ", "Im λ",
                "#unst", "max|λ|", "CN amp", "vs", "η", "k h/π", "t[s]")
println(hdr); push!(rows, hdr); flush(stdout)
for a in AMPS, (tag, c, withN) in variants
    t0 = time()
    if c === nothing
        prob.vms[] = nothing
        U, V, probS = U7, V7, prob
    else
        attach_vms!(prob, trian; p_u = p_u, c_adv = c[1], c_wave = c[2], c_eta = c[3], dt = dt, oss = true)
        U, V, probS = U10, V10, (withN ? prob : probN)
    end
    uh, ut = box_state(a, U, probS)
    A, B   = mats(uh, ut, U, V, probS)
    s      = spectrum(A, B, U)
    r = @sprintf("%-10s %5.2f %+11.4e %10.3e %6d %10.3e %8.5f | %6.3f %6.3f %6.3f %8.1f", tag, a, s.re, s.im,
                 s.nunst, s.maxabs, s.amp, s.vs, s.eta, s.kh, time() - t0)
    println(r); push!(rows, r); flush(stdout)
    open(joinpath(outdir, "report.txt"), "w") do io
        println(io, "oss_cw_spectrum.jl — m5 mixed flat, 1 λ, $(ncell) cells/λ, Q$(p_u)/Q$(p_eta), CN Δt=$(dt), τ_a Δt-capped; a = amplitude / run amplitude (kA=$(kA))")
        println(io, "max Re λ of the frozen B u̇ = −A u (hand core + exact VMS); mode of max Re λ: vs = vertical-structure share of KE, η = PE share, k h/π = effective wavenumber / node Nyquist")
        println(io, "_noN: the stabilisation's residual without 𝓝 (H1 discriminator)\n")
        foreach(x -> println(io, x), rows)
    end
end
println("\nreport → ", joinpath(outdir, "report.txt"))
