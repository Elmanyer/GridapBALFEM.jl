# ==============================================================
#  test_broken_formulation.jl — the broken (skeleton) Class-III path and the C⁰-IP penalty
#
#  Gates markdown_files/BROKEN_FORMULATION_PLAN.md §T6 (src/broken.jl; LaTeX §6.4).
#
#  G1  LAYER IDENTITY. cellwise ∂ₐ𝖲 (hand-expanded from ∇∇) MINUS the skeleton layer
#      ⟦𝖲⟧nₐ equals the mixed weak gradient −∫𝖲∂ₐψ + ∮𝖲ψnₐ — for 𝖲 AND 𝖻, in x AND y,
#      on a sloping bed, walled/open 2-D mesh, asymmetric state, smooth and rough. This
#      is an exact integration-by-parts identity: it fixes the SIGN of the layer, the
#      NORMAL convention and the HESSIAN slot order at once, to round-off. The cellwise
#      part alone is reported as the control (it must be O(1) wrong on a rough state).
#  G2  LAYER LIVE (rule 38d). The skeleton integral is a non-negligible part of the broken
#      Class-III contribution on a rough state.
#  G3  NO LEAKAGE. attach_skeleton! with everything off leaves prob.skel = nothing, and a
#      :native / :none residual is BITWISE unchanged; broken=true without :full is refused.
#  G4  CONSISTENCY. broken and projected (projection evaluated at the SAME state, no lag)
#      Class-III residuals converge to each other under refinement of a smooth state.
#  G5  C⁰-IP ALGEBRA. symmetric, positive semidefinite, zero on a global polynomial of
#      the trial degree (open mesh), and the η-penalty conserves mass (q ≡ 1).
#  G6  C⁰-IP SELECTIVITY. Rayleigh quotient on a λ = 2dx mode ≫ on a resolved mode.
#  G7  JACOBIAN. linear regime + C⁰-IP: the hand ∂R/∂u equals the FD oracle (the linear
#      Jacobian is exact, so the penalty block must be exact too); on the :full broken path
#      the hand↔FD gap is the deliberate Class-III omission and shrinks with amplitude.
#
#  RUN:  julia --project=. test/test_broken_formulation.jl
# ==============================================================
using GridapBALFEM, Gridap, Gridap.ODEs, LinearAlgebra, SparseArrays, Printf, Random

println("="^70); println("  test_broken_formulation.jl — broken Class-III + C⁰-IP"); println("="^70)
np = 0; nf = 0
chk(n, c) = (global np, nf; c ? (println("  PASS  $n"); np += 1) : (println("  FAIL  $n"); nf += 1))

const M, PV, PU, PE = 2, 1, 3, 2
const vert = assemble_vertical_tensors(M, PV, Vector{Float64}(resolve_cbdy(M, nothing, PV)))
const Nσ   = vert.N_dof
const Ex2  = VectorValue(1.0, 0.0)
const Ey2  = VectorValue(0.0, 1.0)
const QDEG = 2 * PU + 4
bed(x)  = 3.5 - 0.08 * x[1] + 0.05 * x[2]
flat(x) = 3.5

function setup(; nx = 6, ny = 3, Lx = 3.0, Ly = 1.5, periodic = false, bc = :open)
    model, trian = build_horizontal_model(((0.0, Lx), (0.0, Ly)), (nx, ny); x_periodic = periodic)
    dΩ   = Measure(trian, QDEG)
    U, V = build_fe_spaces(model, PU, Nσ; y_wall_bc = bc, p_eta = PE)
    return (model = model, trian = trian, dΩ = dΩ, U = U, V = V)
end

"Asymmetric smooth state at amplitude `a`, plus an optional rough (random) seed."
function state(S, a; rough = 0.0, seed = 1, kx = 1.0)
    uh = interpolate_everywhere(
        [x -> a * (0.05 * cos(1.1kx * x[1]) + 0.03 * sin(0.8x[2] + 0.3kx * x[1])),
         x -> VectorValue(ntuple(j -> a * 0.03j * sin(0.7kx * x[1] + 0.2j + 0.4x[2]), Nσ)...),
         x -> VectorValue(ntuple(j -> a * 0.02 * (Nσ + 1 - j) * cos(0.6x[2] + 0.5kx * x[1] - 0.1j), Nσ)...)],
        S.U)
    x = copy(get_free_dof_values(uh))
    rough > 0 && (x .+= rough .* randn(MersenneTwister(seed), length(x)))
    return FEFunction(S.U, x)
end

function mkprob(S; skel = :none, gu = 0.0, ge = 0.0, nlp = :full, regime = :nonlinear,
                bedf = bed, flat_bed = false)
    p = build_problem(vert; h_bathy = bedf, regime = regime, nl_pressure = nlp, flat_bed = flat_bed)
    skel == :broken && attach_skeleton!(p, S.model; broken = true, cip_gamma_u = gu,
                                        cip_gamma_eta = ge, degree = QDEG)
    skel == :cip    && attach_skeleton!(p, S.model; cip_gamma_u = gu, cip_gamma_eta = ge,
                                        degree = QDEG)
    return p
end

"The derived fields exactly as global_residual builds them."
function fields(prob, uh, trian)
    η, Ux, Uy = uh[1], uh[2], uh[3]
    d_cf = CellField(prob.h_bathy, trian)
    H    = d_cf + η
    dhx  = prob.flat_bed ? 0.0 * d_cf : alg_dx(d_cf)
    dhy  = prob.flat_bed ? 0.0 * d_cf : alg_dy(d_cf)
    dHx  = dhx + alg_dx(η); dHy = dhy + alg_dy(η)
    DU   = alg_dx(Ux) + alg_dy(Uy)
    bf   = dHx * Ux + dHy * Uy
    Sf   = H * DU + bf
    return (η = η, Ux = Ux, Uy = Uy, d_cf = d_cf, H = H, dHx = dHx, dHy = dHy, DU = DU, bf = bf, S = Sf)
end

"Residual VECTOR at FE function `uh`, with u̇ frozen (zero unless given)."
function resid(prob, S, uh; udot = nothing)
    ud = udot === nothing ? FEFunction(S.U, zeros(num_free_dofs(S.U))) : udot
    tu = TransientCellField(uh, (ud,))
    return assemble_vector(v -> global_residual(0.0, tu, v, prob, S.trian, S.dΩ), S.V)
end

_avgψ(ψ) = 0.5 * (ψ.plus + ψ.minus)

# ---------------------------------------------------------------------------------------
println("\n  G1 — layer identity: cellwise ∂ₐg − ⟦g⟧nₐ ≡ −∫g∂ₐψ + ∮gψnₐ  (g = 𝖲, 𝖻)")
function g1_errors(S, prob, uh)
    f = fields(prob, uh, S.trian)
    _, _, _, (dxS, dyS, dxb, dyb) =
        broken_class3_cell_fields(prob, f.d_cf, f.η, f.H, f.dHx, f.dHy, f.Ux, f.Uy, f.DU, f.S)
    sk  = prob.skel[]
    Vψ  = FESpace(S.model, ReferenceFE(lagrangian, VectorValue{Nσ,Float64}, PU); conformity = :H1)
    Γ   = BoundaryTriangulation(S.model); dΓ = Measure(Γ, QDEG); nΓ = get_normal_vector(Γ)
    nxΓ = Operation(n -> n ⋅ Ex2)(nΓ); nyΓ = Operation(n -> n ⋅ Ey2)(nΓ)
    errs = Float64[]; ctrl = Float64[]
    for (g, dxg, dyg, nm) in ((f.S, dxS, dyS, "𝖲"), (f.bf, dxb, dyb, "𝖻"))
        jg = g.plus - g.minus
        for (dg, na, nΓa, dψ, lbl) in ((dxg, sk.nx, nxΓ, alg_dx, "x"), (dyg, sk.ny, nyΓ, alg_dy, "y"))
            lay  = assemble_vector(ψ -> ∫(dg ⋅ ψ) * S.dΩ - ∫((jg ⋅ _avgψ(ψ)) * na) * sk.dΛ, Vψ)
            cell = assemble_vector(ψ -> ∫(dg ⋅ ψ) * S.dΩ, Vψ)
            weak = assemble_vector(ψ -> ∫(-(g ⋅ dψ(ψ))) * S.dΩ + ∫((g ⋅ ψ) * nΓa) * dΓ, Vψ)
            e = norm(lay - weak) / norm(weak); ec = norm(cell - weak) / norm(weak)
            @printf("    ∂%s%s  |cell+layer − weak|/|weak| = %.2e   (cellwise alone %.2e)\n", lbl, nm, e, ec)
            push!(errs, e); push!(ctrl, ec)
        end
    end
    return errs, ctrl
end
S1 = setup()
pb1 = mkprob(S1; skel = :broken)
for (lbl, uh) in (("smooth", state(S1, 1.0)), ("rough", state(S1, 1.0; rough = 0.02)))
    println("   $lbl state:")
    e, c = g1_errors(S1, pb1, uh)
    chk("G1 [$lbl] cellwise Hessians + skeleton layer = weak gradient (max $(round(maximum(e), sigdigits=2)))",
        maximum(e) < 1e-10)
    lbl == "rough" && chk("G1 [rough] the cellwise part ALONE is O(1) wrong — the layer is load-bearing (min $(round(minimum(c), sigdigits=2)))",
                          minimum(c) > 0.1)
end

# ---------------------------------------------------------------------------------------
println("\n  G2 — the skeleton layer is live in the assembled residual (rough state)")
let uh = state(S1, 1.0; rough = 0.01)
    r_none = resid(mkprob(S1), S1, uh)                  # projected, no projection state ⇒ no 𝓚/𝓟 Class III
    r_brk  = resid(pb1, S1, uh)
    f      = fields(pb1, uh, S1.trian)
    r_lay  = assemble_vector(v -> broken_class3_skeleton_contrib(pb1, pb1.skel[], f.H, f.dHx, f.dHy,
                                  f.Ux, f.Uy, f.S, f.bf, v[2], v[3], alg_dx(v[2]) + alg_dy(v[3])), S1.V)
    ratio  = norm(r_lay) / norm(r_brk - r_none)
    @printf("    |broken Class-III 𝓚+𝓟| = %.3e   |skeleton layer| = %.3e   ratio %.2f\n",
            norm(r_brk - r_none), norm(r_lay), ratio)
    chk("G2 the skeleton layer is a non-negligible part of the broken Class-III blocks (ratio $(round(ratio, sigdigits=2)))",
        0.05 < ratio < 20)
end

# ---------------------------------------------------------------------------------------
println("\n  G3 — no leakage")
let uh = state(S1, 1.0; rough = 0.01)
    pn = mkprob(S1; nlp = :native)
    r1 = resid(pn, S1, uh)
    attach_skeleton!(pn, S1.model; broken = false, cip_gamma_u = 0.0, cip_gamma_eta = 0.0, degree = QDEG)
    chk("G3 everything off ⇒ prob.skel stays nothing", pn.skel[] === nothing)
    chk("G3 :native residual bitwise unchanged", resid(pn, S1, uh) == r1)
    refused = try
        attach_skeleton!(mkprob(S1; nlp = :native), S1.model; broken = true, degree = QDEG); false
    catch
        true
    end
    chk("G3 broken=true without nl_pressure=:full is refused", refused)
end

# ---------------------------------------------------------------------------------------
println("\n  G4 — broken and projected (same state, no lag) converge to each other")
const VSTAR = [x -> 0.3 * sin(0.9x[1] + 0.2x[2]),
               x -> VectorValue(ntuple(j -> cos(0.5x[1] + 0.3j - 0.4x[2]), Nσ)...),
               x -> VectorValue(ntuple(j -> 0.5 * sin(0.7x[2] + 0.2x[1] + 0.1j), Nσ)...)]
function g4(nx)
    Sm = setup(nx = nx, ny = nx ÷ 2)
    u  = state(Sm, 1.0)
    vc = get_free_dof_values(interpolate_everywhere(VSTAR, Sm.U))
    r0 = resid(mkprob(Sm), Sm, u)
    rb = resid(mkprob(Sm; skel = :broken), Sm, u)
    pp = mkprob(Sm); update_nlp_state!(pp, build_nlp_ctx(Sm.model, PU, Nσ, Sm.trian, Sm.dΩ), u)
    rp = resid(pp, Sm, u)
    Jb = dot(rb - r0, vc); Jp = dot(rp - r0, vc)
    @printf("    nx=%3d  J_broken=% .7e  J_projected=% .7e  |diff|/|J| = %.3e\n", nx, Jb, Jp, abs(Jb - Jp) / abs(Jb))
    return abs(Jb - Jp) / abs(Jb)
end
let ds = [g4(n) for n in (6, 12, 24)]
    rates = log2.(ds[1:end-1] ./ ds[2:end])
    chk("G4 the two treatments converge to each other (pairwise orders $(round.(rates, digits=2)))",
        all(rates .> 2.0) && ds[end] < 1e-5)
end

# ---------------------------------------------------------------------------------------
println("\n  G5 — C⁰-IP algebra (open 2-D mesh, sloping bed)")
cipmat(Sm, p) = assemble_matrix((du, v) -> cip_contrib(p, p.skel[], du[1], du[2], du[3], v[1], v[2], v[3]), Sm.U, Sm.V)
let pc = mkprob(S1; skel = :cip, gu = 1.0, ge = 1.0, nlp = :none)
    A  = cipmat(S1, pc); Ad = Matrix(A)
    λ  = eigvals(Symmetric(0.5 * (Ad + Ad')))
    asym = norm(Ad - Ad') / norm(Ad)
    @printf("    asym %.2e  eigmin %.2e  eigmax %.2e\n", asym, minimum(λ), maximum(λ))
    chk("G5 symmetric", asym < 1e-12)
    chk("G5 positive semidefinite", minimum(λ) > -1e-10 * maximum(λ))
    poly = interpolate_everywhere([x -> 0.3 + 0.2x[1] - 0.1x[2] + 0.05x[1]^2 - 0.07x[1] * x[2] + 0.04x[2]^2,
        x -> VectorValue(ntuple(j -> 0.1j * x[1]^3 - 0.2x[1] * x[2]^2 + 0.3x[2]^3 + j * x[1] * x[2], Nσ)...),
        x -> VectorValue(ntuple(j -> 0.2x[2]^3 - 0.1j * x[1]^2 * x[2] + 0.05x[1]^3, Nσ)...)], S1.U)
    xp = get_free_dof_values(poly)
    e  = norm(A * xp) / (opnorm(Ad) * norm(xp))
    chk("G5 zero on a global polynomial of the trial degree ($(round(e, sigdigits=2)))", e < 1e-12)
    xr = randn(MersenneTwister(3), length(xp))
    z  = VectorValue(ntuple(_ -> 0.0, Nσ)...)
    q1 = get_free_dof_values(interpolate_everywhere([x -> 1.0, x -> z, x -> z], S1.U))
    m  = abs(dot(q1, A * xr)) / norm(A * xr)
    chk("G5 the η-penalty conserves mass: (q≡1)·J_h = $(round(m, sigdigits=2))", m < 1e-12)
end

# ---------------------------------------------------------------------------------------
println("\n  G6 — C⁰-IP is grid-scale selective (x-periodic box, Q3/Q2, 16 cells/λ)")
const SP = setup(nx = 16, ny = 1, Lx = 4.0, Ly = 0.25, periodic = true, bc = :wall)
let pcp = mkprob(SP; skel = :cip, gu = 1.0, nlp = :none, bedf = flat, flat_bed = true)
    Ap = cipmat(SP, pcp)
    Mm = assemble_matrix((du, v) -> ∫(du[2] ⋅ v[2]) * SP.dΩ, SP.U, SP.V)
    z  = VectorValue(ntuple(_ -> 0.0, Nσ)...)
    mode(k) = get_free_dof_values(interpolate_everywhere(
        [x -> 0.0, x -> VectorValue(ntuple(j -> cos(k * x[1]), Nσ)...), x -> z], SP.U))
    rq(x) = dot(x, Ap * x) / dot(x, Mm * x)
    rc = rq(mode(2π / 4)); rg = rq(mode(π / 0.25))
    @printf("    Rayleigh quotient: carrier %.3e   λ = 2dx %.3e   ratio %.2e\n", rc, rg, rg / rc)
    chk("G6 the λ = 2dx mode is penalised > 1e4 × the carrier", rg / rc > 1e4)
end

# ---------------------------------------------------------------------------------------
println("\n  G7 — Jacobian: hand ∂R/∂u vs central FD of the residual")
function jac_fd_gap(p, Sm, x0, d; h = 1e-7)
    ud  = FEFunction(Sm.U, zeros(length(x0)))
    res(x) = assemble_vector(v -> global_residual(0.0, TransientCellField(FEFunction(Sm.U, x), (ud,)), v,
                                                  p, Sm.trian, Sm.dΩ), Sm.V)
    J   = assemble_matrix((du, v) -> jacobian_u(0.0, TransientCellField(FEFunction(Sm.U, x0), (ud,)), du, v,
                                                p, Sm.trian, Sm.dΩ), Sm.U, Sm.V)
    fd  = (res(x0 .+ h .* d) .- res(x0 .- h .* d)) ./ (2h)
    return maximum(abs, J * d .- fd) / maximum(abs, fd)
end
let x0 = get_free_dof_values(state(S1, 1.0; rough = 0.01)), d = randn(MersenneTwister(7), length(x0))
    gl = jac_fd_gap(mkprob(S1; skel = :cip, gu = 0.3, ge = 0.3, nlp = :none, regime = :linear), S1, x0, d)
    @printf("    linear + C⁰-IP: rel gap %.2e\n", gl)
    chk("G7 linear regime + C⁰-IP: the Jacobian is exact (penalty block included)", gl < 1e-6)
    pbk = mkprob(S1; skel = :broken, gu = 0.3, ge = 0.3)
    gaps = Float64[]
    for a in (1.0, 0.5, 0.25)
        xa = get_free_dof_values(state(S1, a; rough = 0.01a))
        push!(gaps, jac_fd_gap(pbk, S1, xa, d))
        @printf("    :full broken + C⁰-IP, a = %.2f: rel gap %.3e\n", a, gaps[end])
    end
    chk("G7 :full broken: the quasi-Newton gap shrinks with amplitude (the deliberate omission)",
        all(gaps[1:end-1] ./ gaps[2:end] .> 1.2))
end

# ---------------------------------------------------------------------------------------
println("\n  G8 — hp-CIP, orders 1–2 (BROKEN_FORMULATION_PLAN.md §5.2)")
let S2 = setup()
    for ord in (1, 2)
        p = build_problem(vert; h_bathy = bed, regime = :nonlinear, nl_pressure = :none, flat_bed = false)
        attach_skeleton!(p, S2.model; cip_gamma_u = 1.0, cip_gamma_eta = 1.0, cip_order = ord,
                         p_u = PU, p_eta = PE, degree = QDEG)
        A  = cipmat(S2, p); Ad = Matrix(A)
        λ  = eigvals(Symmetric(0.5 * (Ad + Ad')))
        poly = interpolate_everywhere([x -> 0.3 + 0.2x[1] - 0.1x[2] + 0.05x[1]^2 - 0.07x[1] * x[2] + 0.04x[2]^2,
            x -> VectorValue(ntuple(j -> 0.1j * x[1]^3 - 0.2x[1] * x[2]^2 + 0.3x[2]^3 + j * x[1] * x[2], Nσ)...),
            x -> VectorValue(ntuple(j -> 0.2x[2]^3 - 0.1j * x[1]^2 * x[2] + 0.05x[1]^3, Nσ)...)], S2.U)
        xp = get_free_dof_values(poly)
        e  = norm(A * xp) / (opnorm(Ad) * norm(xp))
        rk = count(λ .> 1e-10 * maximum(λ))
        @printf("    order %d: orders (u,η) = (%d,%d)  asym %.1e  λmin/λmax %.1e  |A·poly| %.1e  rank %d/%d\n",
                ord, p.skel[].order_u, p.skel[].order_eta, norm(Ad - Ad') / norm(Ad),
                minimum(λ) / maximum(λ), e, rk, length(λ))
        chk("G8 [order $ord] symmetric PSD, zero on a global polynomial of the trial degree",
            norm(Ad - Ad') / norm(Ad) < 1e-12 && minimum(λ) > -1e-10 * maximum(λ) && e < 1e-12)
        ord == 2 && (global RANK2 = rk)
        ord == 1 && (global RANK1 = rk)
    end
    chk("G8 the order-2 terms are live: rank grows ($(RANK1) → $(RANK2))", RANK2 > RANK1)
    refused = try
        attach_skeleton!(mkprob(S2; nlp = :none), S2.model; cip_gamma_u = 1.0, cip_order = 3, degree = QDEG); false
    catch
        true
    end
    chk("G8 cip_order = 3 is refused (Gridap evaluates FE derivatives only up to order 2)", refused)
end

# =======================================================================================
#  G9 — :ghostvolume, the direct (volume) ghost penalty (GHOST_PENALTY_PLAN.md §3)
# =======================================================================================
println("\n  G9 — ghost-volume stabilisation")
"Problem with a given stabilisation attached (flat or sloping bed)."
function gprob(Sm; stab = :ghostvolume, gu = 1.0, ge = 1.0, regime = :nonlinear, nlp = :none,
               bedf = bed, flat_bed = false, pu = PU, pe = PE, order = 1)
    p = build_problem(vert; h_bathy = bedf, regime = regime, nl_pressure = nlp, flat_bed = flat_bed)
    attach_skeleton!(p, Sm.model; cip_gamma_u = gu, cip_gamma_eta = ge, stabilization = stab,
                     cip_order = order, p_u = pu, p_eta = pe, degree = 2 * pu + 4)
    return p
end
stabmat(Sm, p) = assemble_matrix((du, v) -> stab_contrib(p, p.skel[], du[1], du[2], du[3], v[1], v[2], v[3]), Sm.U, Sm.V)

let S2 = setup()                                       # open 2-D mesh, Q3/Q2, sloping bed
    p  = gprob(S2)
    A  = stabmat(S2, p); Ad = Matrix(A)
    λ  = eigvals(Symmetric(0.5 * (Ad + Ad')))
    asym = norm(Ad - Ad') / norm(Ad)
    poly = interpolate_everywhere([x -> 0.3 + 0.2x[1] - 0.1x[2] + 0.05x[1]^2 - 0.07x[1] * x[2] + 0.04x[2]^2,
        x -> VectorValue(ntuple(j -> 0.1j * x[1]^3 - 0.2x[1] * x[2]^2 + 0.3x[2]^3 + j * x[1] * x[2], Nσ)...),
        x -> VectorValue(ntuple(j -> 0.2x[2]^3 - 0.1j * x[1]^2 * x[2] + 0.05x[1]^3, Nσ)...)], S2.U)
    xp = get_free_dof_values(poly)
    e  = norm(A * xp) / (opnorm(Ad) * norm(xp))
    z  = VectorValue(ntuple(_ -> 0.0, Nσ)...)
    q1 = get_free_dof_values(interpolate_everywhere([x -> 1.0, x -> z, x -> z], S2.U))
    xr = randn(MersenneTwister(5), length(xp)); m = abs(dot(q1, A * xr)) / norm(A * xr)
    rk = count(λ .> 1e-10 * maximum(λ))
    @printf("    Q3/Q2 open mesh: asym %.1e  λmin/λmax %.1e  |A·poly| %.1e  mass %.1e  rank %d/%d\n",
            asym, minimum(λ) / maximum(λ), e, m, rk, length(λ))
    chk("G9a kernel: zero on a global polynomial of the trial degree ($(round(e, sigdigits=2)))", e < 1e-12)
    chk("G9b symmetric positive semidefinite", asym < 1e-12 && minimum(λ) > -1e-10 * maximum(λ))
    chk("G9d the η-penalty conserves mass ($(round(m, sigdigits=2)))", m < 1e-12)
    pj = gprob(S2; stab = :jumpgrad, order = 2); rj = count(eigvals(Symmetric(Matrix(stabmat(S2, pj)))) .> 1e-10 * maximum(λ))
    chk("G9b the ghost penalty sees MORE than the order-≤2 jump penalty (rank $rj → $rk): it reaches order 3", rk > rj)
end

#  G9c — the closed form (★): for p ≤ 2 every jump is computable, so the ghost matrix must equal
#  γ τ h^(s−3) Σ_{j,k} G_jk ∫_F ⟦∂ₙʲ·⟧⟦∂ₙᵏ·⟧,  G_jk = [1+(−1)^(j+k)] h^(j+k+1)/((j+k+1) j! k!)
function star_matrix(Sm, p, hx, hy)
    sk = p.skel[]; nx, ny = sk.nx, sk.ny
    hn = hx * (nx * nx) + hy * (ny * ny)
    Gjk(j, k) = Operation(h -> (1 + (-1)^(j + k)) * h^(j + k + 1) / ((j + k + 1) * factorial(j) * factorial(k)) *
                               h^(sk.hexp - 3))(hn)
    dn(f, j) = j == 1 ? (nx * (alg_dx(f).plus - alg_dx(f).minus) + ny * (alg_dy(f).plus - alg_dy(f).minus)) :
        ((nx * nx) * (alg_hess(f, 1, 1).plus - alg_hess(f, 1, 1).minus) +
         (2.0 * nx * ny) * (alg_hess(f, 1, 2).plus - alg_hess(f, 1, 2).minus) +
         (ny * ny) * (alg_hess(f, 2, 2).plus - alg_hess(f, 2, 2).minus))
    g = p.g; dval = 3.5
    τu = dval * sqrt(g * dval); τη = sqrt(g * dval)
    form(u, v) = ∫( sum(Gjk(j, k) * (τu * sk.gu) * ((alg_mul(p.Mv, dn(u[2], j)) ⋅ dn(v[2], k)) +
                                                     (alg_mul(p.Mv, dn(u[3], j)) ⋅ dn(v[3], k)))
                        for j in 1:2, k in 1:2) +
                    (Gjk(1, 1) * (τη * sk.ge)) * (dn(u[1], 1) * dn(v[1], 1)) ) * sk.dΛ
    return assemble_matrix(form, Sm.U, Sm.V)
end
for (lbl, Sm, hx, hy) in (("open 2-D mesh", nothing, 0.5, 0.5), ("x-periodic box (wrap facet)", nothing, 0.5, 0.5))
    Sq = lbl == "open 2-D mesh" ?
        (let (m, t) = build_horizontal_model(((0.0, 3.0), (0.0, 1.5)), (6, 3))
             U, V = build_fe_spaces(m, 2, Nσ; y_wall_bc = :open, p_eta = 1)
             (model = m, trian = t, dΩ = Measure(t, 8), U = U, V = V) end) :
        (let (m, t) = build_horizontal_model(((0.0, 4.0), (0.0, 0.5)), (8, 1); x_periodic = true)
             U, V = build_fe_spaces(m, 2, Nσ; y_wall_bc = :wall, p_eta = 1)
             (model = m, trian = t, dΩ = Measure(t, 8), U = U, V = V) end)
    p = gprob(Sq; bedf = flat, flat_bed = true, pu = 2, pe = 1, gu = 1.0, ge = 1.0)
    A = stabmat(Sq, p); B = star_matrix(Sq, p, hx, hy)
    e = norm(A - B) / norm(B)
    @printf("    Q2/Q1 %-28s |A_ghost − A_(★)|/|A_(★)| = %.2e\n", lbl, e)
    chk("G9c ghost ≡ closed form (★) on the $lbl ($(round(e, sigdigits=2)))", e < 1e-11)
end

let x0 = get_free_dof_values(state(S1, 1.0; rough = 0.01)), d = randn(MersenneTwister(9), length(x0))
    gl = jac_fd_gap(gprob(S1; regime = :linear, gu = 0.3, ge = 0.3), S1, x0, d)
    @printf("    linear + ghost: hand vs FD rel gap %.2e\n", gl)
    chk("G9e linear regime + ghost: the Jacobian is exact", gl < 1e-6)
end

let pcp = gprob(SP; gu = 1.0, ge = 0.0, bedf = flat, flat_bed = true)
    Ap = stabmat(SP, pcp)
    Mm = assemble_matrix((du, v) -> ∫(du[2] ⋅ v[2]) * SP.dΩ, SP.U, SP.V)
    z  = VectorValue(ntuple(_ -> 0.0, Nσ)...)
    mode(k) = get_free_dof_values(interpolate_everywhere(
        [x -> 0.0, x -> VectorValue(ntuple(j -> cos(k * x[1]), Nσ)...), x -> z], SP.U))
    rq(x) = dot(x, Ap * x) / dot(x, Mm * x)
    rc = rq(mode(2π / 4)); rg = rq(mode(π / 0.25))
    @printf("    ghost Rayleigh quotient: carrier %.3e   λ = 2dx %.3e   ratio %.2e\n", rc, rg, rg / rc)
    chk("G9f grid-scale selective (λ = 2dx / carrier > 1e4)", rg / rc > 1e4)
end

let uh = state(S1, 1.0; rough = 0.01)
    p_def = build_problem(vert; h_bathy = bed, regime = :nonlinear, nl_pressure = :native, flat_bed = false)
    attach_skeleton!(p_def, S1.model; cip_gamma_u = 0.3, cip_gamma_eta = 0.3, cip_order = 2, p_u = PU, p_eta = PE, degree = QDEG)
    p_jg = gprob(S1; stab = :jumpgrad, nlp = :native, gu = 0.3, ge = 0.3, order = 2)
    chk("G9g default stabilization is :jumpgrad, bitwise", resid(p_def, S1, uh) == resid(p_jg, S1, uh))
    r1 = try gprob(S1; order = 2); false catch; true end
    chk("G9g :ghostvolume with cip_order ≠ 1 is refused", r1)
    nm = CartesianDiscreteModel((0.0, 1.0, 0.0, 1.0), (4, 2); map = x -> VectorValue(x[1]^2, x[2]))
    r2 = try
        pp = build_problem(vert; h_bathy = flat, regime = :nonlinear, nl_pressure = :none, flat_bed = true)
        attach_skeleton!(pp, nm; cip_gamma_u = 1.0, stabilization = :ghostvolume, p_u = PU, p_eta = PE, degree = QDEG); false
    catch
        true
    end
    chk("G9g a non-uniform mesh is refused", r2)
end

# =======================================================================================
#  G10 — the EXACT Jacobian of the broken Class-III blocks (volume + skeleton layer)
#  Oracle: central FD of the Class-III residual ALONE, R_broken − R_(no Class III).
# =======================================================================================
println("\n  G10 — broken Class-III Jacobian vs FD")
function c3_gap(Sm, x0, d; h = 1e-7, bedf = bed, flat_bed = false)
    pb = build_problem(vert; h_bathy = bedf, regime = :nonlinear, nl_pressure = :full, flat_bed = flat_bed)
    attach_skeleton!(pb, Sm.model; broken = true, degree = QDEG)
    p0 = build_problem(vert; h_bathy = bedf, regime = :nonlinear, nl_pressure = :full, flat_bed = flat_bed)
    ud = FEFunction(Sm.U, zeros(length(x0)))
    R(p, x) = assemble_vector(v -> global_residual(0.0, TransientCellField(FEFunction(Sm.U, x), (ud,)), v,
                                                   p, Sm.trian, Sm.dΩ), Sm.V)
    Rc3(x) = R(pb, x) - R(p0, x)
    fd = (Rc3(x0 .+ h .* d) .- Rc3(x0 .- h .* d)) ./ (2h)
    uh = FEFunction(Sm.U, x0); d_cf = CellField(pb.h_bathy, Sm.trian)
    J  = assemble_matrix((du, v) -> broken_class3_jacobian(pb, pb.skel[], d_cf, uh[1], uh[2], uh[3],
                                                          du[1], du[2], du[3], v[2], v[3], Sm.dΩ), Sm.U, Sm.V)
    return maximum(abs, J * d .- fd) / maximum(abs, fd)
end
let x0 = get_free_dof_values(state(S1, 1.0; rough = 0.01)), d = randn(MersenneTwister(7), length(x0))
    g = c3_gap(S1, x0, d)
    @printf("    open 2-D mesh, sloping bed, Q3/Q2: rel gap %.2e\n", g)
    chk("G10 broken Class-III Jacobian = FD of the Class-III residual (sloping bed, 2-D)", g < 1e-6)
    gp = let xp = get_free_dof_values(state(SP, 1.0; rough = 0.01)), dp = randn(MersenneTwister(8), length(xp))
        c3_gap(SP, xp, dp; bedf = flat, flat_bed = true)
    end
    @printf("    x-periodic box, flat bed, Q3/Q2:   rel gap %.2e\n", gp)
    chk("G10 broken Class-III Jacobian = FD of the Class-III residual (periodic box, wrap facet)", gp < 1e-6)
end

println("\n" * "="^70); @printf("  %d passed, %d failed\n", np, nf); println("="^70)
nf == 0 || error("test_broken_formulation: $nf gate(s) failed")
