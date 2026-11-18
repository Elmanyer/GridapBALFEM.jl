# ==============================================================
#  vms_jacobian_gap.jl — why does the hand (quasi-Newton) Jacobian stop converging with VMS?
#  (branch mixed-vms-solver)
#
#  On the box state of compare_mixed_vms.jl (model 3, flat bed, hot-start eigenmode, 16 cells/λ,
#  Q3/Q2) it assembles the two matrices the integrator uses, by hand and by AD of the same residual:
#      A = ∂R/∂u,  B = ∂R/∂u̇,   CN Newton matrix  J = A + B/(θΔt)
#  and reports, per VMS variant:
#      relA, relB   ‖hand − AD‖/‖AD‖
#      ΔA_vms       ‖(A_hand − A_ad)_with VMS − (A_hand − A_ad)_without‖: the VMS part of the gap
#                   (0 ⇔ the VMS Jacobian is exact AND wired correctly; the rest is the core gap)
#      ρ            spectral radius of I − J_hand⁻¹ J_ad = the asymptotic contraction factor of the
#                   quasi-Newton iteration per iteration (ρ < 1 converges, ρ ≈ 1 stalls, ρ > 1 diverges)
#  and the frozen linearised spectrum of B u̇ = −A u (exact AD matrices; algebraic fields eliminated
#  by Schur complement): max Re λ, max |λ|, min |eig B|, the number of eigenvalues of B with a
#  negative real part (an indefinite effective mass), and the CN amplification.
#
#  BALFEM_GAP_OSS=1: the same with orthogonal subgrid scales (oss=true; 3 projection fields, the
#  state carries consistent projections). Output report_oss.txt instead of report.txt.
#  Output: output/mixed_vms/jacobian_gap/report[_oss].txt
# ==============================================================
using GridapBALFEM, Gridap, Gridap.ODEs, LinearAlgebra, Printf

const ROOT_G = normpath(joinpath(@__DIR__, "..", ".."))
outdir = joinpath(ROOT_G, "output", "mixed_vms", "jacobian_gap"); mkpath(outdir)
M, p_vert, p_u, p_eta, d, Twave, kA, ncell, dt, θ = 2, 1, 3, 2, 3.5, 1.6, 0.08, 16, 0.04, 0.5
vert  = assemble_vertical_tensors(M, p_vert, Vector{Float64}(resolve_cbdy(M, nothing, p_vert)))
Nσ    = vert.N_dof
k0    = GridapBALFEM.model_wavenumber(vert, 2π / Twave, d, 9.81)
Awave = kA / k0
wi    = WaveInput(vert; A = Awave, T = Twave, d = d, T_ramp = 0.0, profile = :model)
Lx    = 2π / wi.ks[1]; dx = Lx / ncell
inc   = incident_fields(wi)
QDEG  = 2p_u + 4
zv    = VectorValue(ntuple(_ -> 0.0, Nσ)...)

model, trian = build_horizontal_model(((0.0, Lx), (0.0, dx)), (ncell, 1); x_periodic = true)
dΩ   = Measure(trian, QDEG)
prob = build_problem(vert; h_bathy = x -> d, regime = :nonlinear, nl_pressure = false, flat_bed = true,
                     model = model, quad_degree = QDEG, p_eta = p_eta)
const OSS = get(ENV, "BALFEM_GAP_OSS", "0") == "1"
U, V = build_fe_spaces(model, p_u, Nσ; y_wall_bc = :wall, p_eta = p_eta, n_oss = OSS ? 3 : 0)
nphys = sum(num_free_dofs(U.spaces[i]) for i in 1:3)

"Box state at amplitude factor a (a = 1: the run's hot start), with its exact time derivative (FD in t)."
function box_state(a; t = 0.0, ε = 1e-5)
    pad = OSS ? Any[x -> zv, x -> zv, x -> 0.0] : Any[]
    uh = interpolate_everywhere(Any[x -> a * inc.eta(x, t), x -> a * inc.ux(x, t), x -> zv, pad...], U)
    ut = interpolate_everywhere(Any[x -> a * (inc.eta(x, t + ε) - inc.eta(x, t - ε)) / 2ε,
                                    x -> a * (inc.ux(x, t + ε) - inc.ux(x, t - ε)) / 2ε, x -> zv, pad...], U)
    has_oss(prob) && (uh = oss_consistent_ic(uh, prob, U, trian, dΩ)[1])
    return uh, ut
end

"Frozen spectrum of B u̇ = −A u on the physics block (algebraic fields eliminated)."
function spectrum(A, B)
    p = 1:nphys; q = (nphys + 1):size(A, 1)
    Ar = isempty(q) ? A[p, p] : A[p, p] - A[p, q] * (A[q, q] \ A[q, p])
    Bp = B[p, p]
    λ  = eigvals(-Ar, Bp)
    λ  = λ[isfinite.(λ)]
    eB = eigvals(Bp)
    amp = maximum(abs, (1 .+ (1 - θ) * dt .* λ) ./ (1 .- θ * dt .* λ))
    return maximum(real, λ), maximum(abs, λ), minimum(abs, eB), count(e -> real(e) < 0, eB), amp
end

function mats(uh, ut)
    tu = TransientCellField(uh, (ut,))
    Ah = assemble_matrix((du, v)  -> jacobian_u(0.0, tu, du, v, prob, trian, dΩ), U, V)
    Bh = assemble_matrix((dw, v)  -> jacobian_u_t(0.0, tu, dw, v, prob, trian, dΩ), U, V)
    Aa = jacobian(FEOperator((x, v) -> global_residual(0.0, TransientCellField(x, (ut,)), v, prob, trian, dΩ), U, V), uh)
    Ba = jacobian(FEOperator((w, v) -> global_residual(0.0, TransientCellField(uh, (w,)), v, prob, trian, dΩ), U, V), ut)
    return Matrix(Ah), Matrix(Bh), Matrix(Aa), Matrix(Ba)
end

variants = [("no VMS", nothing), ("all1", (1.0, 1.0, 1.0)), ("adv1", (1.0, 0.0, 0.0)),
            ("wave1", (0.0, 1.0, 0.0)), ("eta1", (0.0, 0.0, 1.0)), ("all0.1", (0.1, 0.1, 0.1))]
rows = String[]
push!(rows, @sprintf("%-8s %5s %10s %10s %10s %10s %8s  %-6s | %10s %10s %10s %5s %8s", "variant", "a", "relA", "relB",
                     "ΔA_vms", "ΔB_vms", "ρ", "gap", "max Re λ", "max |λ|", "min|eig B|", "B<0", "CN amp"))
nη = num_free_dofs(U.spaces[1]); nx_ = num_free_dofs(U.spaces[2])
blk(i) = i ≤ nη ? "η" : i ≤ nη + nx_ ? "Ux" : i ≤ nphys ? "Uy" : "Π/δ"
setv!(c) = OSS ? (c === nothing ? (prob.vms[] = nothing) :
                  attach_vms!(prob, trian; p_u = p_u, c_adv = c[1], c_wave = c[2], c_eta = c[3], oss = true)) :
                 attach_vms!(prob, trian; p_u = p_u, c_adv = c === nothing ? 0.0 : c[1],
                             c_wave = c === nothing ? 0.0 : c[2], c_eta = c === nothing ? 0.0 : c[3])
AMPS = parse.(Float64, split(get(ENV, "BALFEM_GAP_AMPS", "0.25,0.5,1,2"), ","))
for a in AMPS
    for (tag, c) in variants
        #  the state is rebuilt per variant: with OSS its consistent projections depend on the τ's
        setv!(c === nothing ? (OSS ? (1.0, 1.0, 1.0) : nothing) : c)
        uh, ut = box_state(a)
        setv!(nothing)
        Ah0, Bh0, Aa0, Ba0 = mats(uh, ut)
        if c === nothing
            Ah, Bh, Aa, Ba = Ah0, Bh0, Aa0, Ba0
        else
            setv!(c)
            Ah, Bh, Aa, Ba = mats(uh, ut)
            setv!(nothing)
        end
        Jh = Ah + Bh / (θ * dt); Ja = Aa + Ba / (θ * dt)
        #  no VMS on the OSS layout: the projection rows are absent, J is singular; ρ on the physics block
        ρ  = (OSS && c === nothing) ? maximum(abs, eigvals(I - Jh[1:nphys, 1:nphys] \ Ja[1:nphys, 1:nphys])) :
                                      maximum(abs, eigvals(I - Jh \ Ja))
        ΔA = norm((Ah - Aa) - (Ah0 - Aa0)) / norm(Aa); ΔB = norm((Bh - Ba) - (Bh0 - Ba0)) / norm(Ba)
        G  = abs.((Jh - Ja)); ij = argmax(G)
        sp = (OSS && c === nothing) ? spectrum(Aa[1:nphys, 1:nphys], Ba[1:nphys, 1:nphys]) : spectrum(Aa, Ba)
        push!(rows, @sprintf("%-8s %5.2f %10.3e %10.3e %10.3e %10.3e %8.4f  %-6s | %10.3e %10.3e %10.3e %5d %8.4f",
                             tag, a, norm(Ah - Aa) / norm(Aa), norm(Bh - Ba) / norm(Ba), ΔA, ΔB, ρ,
                             blk(ij[1]) * "×" * blk(ij[2]), sp...))
        println(rows[end]); flush(stdout)
    end
end
rname = get(ENV, "BALFEM_GAP_REPORT", OSS ? "report_oss.txt" : "report.txt")
open(joinpath(outdir, rname), "w") do io
    println(io, "vms_jacobian_gap.jl — model 3, flat, box 1 λ, $(ncell) cells/λ, Q$(p_u)/Q$(p_eta), CN Δt=$(dt); a = amplitude / run amplitude (kA=$(kA))")
    println(io, OSS ? "STABILISATION: orthogonal subgrid scales (oss=true)" : "STABILISATION: ASGS (full residual incl. u̇)")
    println(io, "ρ = spectral radius of I − J_hand⁻¹ J_AD (quasi-Newton contraction per iteration)\n")
    foreach(r -> println(io, r), rows)
end
println("\nreport → ", joinpath(outdir, rname))
