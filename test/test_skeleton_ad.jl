# ==============================================================
#  test_skeleton_ad.jl — automatic differentiation through SKELETON integrals
#
#  The v2 residual has skeleton (interior-facet) terms: the broken Class-III layer and the
#  skeleton stabilisers. Differentiating them by AD needs the Gridap fork's skeleton-AD patch
#  (Elmanyer/Gridap.jl @ fix-transient-multifield-ad, commit "Fix: automatic differentiation of
#  skeleton integrals on transient fields"): `DomainStyle(::Type{SkeletonCellFieldPair})` and
#  `TransientSingleFieldCellField` `.plus/.minus` delegating to the wrapped field.
#
#  S0  the patch is present (the type-level DomainStyle method exists).
#  S1  AD Jacobian of the broken Class-III residual (`broken_class3_residual`) = the hand
#      `broken_class3_jacobian`, to round-off — sloping 2-D mesh, and the x-periodic box (wrap
#      facet). Both are exact, so they must agree entry by entry.
#  S2  AD Jacobian of each skeleton stabiliser (:jumpgrad order 2, :ghostvolume) = its own form
#      on the increments (it is linear), to round-off.
#  S3  a transient run with use_ad=true (AD Jacobians of the WHOLE residual, skeleton included)
#      reaches the same answer as the hand Jacobians (rule 17b: the Jacobian sets the path, not
#      the root) — closed box, nl_pressure=true + ghost penalty, Crank–Nicolson.
#
#  RUN:  julia --project=. test/test_skeleton_ad.jl
# ==============================================================
using GridapBALFEM, Gridap, Gridap.ODEs, LinearAlgebra, SparseArrays, Printf, Random

println("="^70); println("  test_skeleton_ad.jl — AD through skeleton integrals"); println("="^70)
np = 0; nf = 0
chk(n, c) = (global np, nf; c ? (println("  PASS  $n"); np += 1) : (println("  FAIL  $n"); nf += 1); flush(stdout))

const PU, PE = 3, 2
const QDEG   = 2 * PU + 4
const vert   = assemble_vertical_tensors(2, 1, Vector{Float64}(resolve_cbdy(2, nothing, 1)))
const Nσ     = vert.N_dof
bed(x)  = 3.5 - 0.08 * x[1] + 0.05 * x[2] + 0.03 * sin(1.3 * x[1]) * cos(0.9 * x[2])
flat(x) = 3.5

function setup(; nx = 6, ny = 3, Lx = 3.0, Ly = 1.5, periodic = false, bc = :open)
    model, trian = build_horizontal_model(((0.0, Lx), (0.0, Ly)), (nx, ny); x_periodic = periodic)
    U, V = build_fe_spaces(model, PU, Nσ; y_wall_bc = bc, p_eta = PE)
    return (model = model, trian = trian, dΩ = Measure(trian, QDEG), U = U, V = V)
end

function state(S, a; seed)
    uh = interpolate_everywhere(
        [x -> a * (0.05 * cos(1.1x[1]) + 0.03 * sin(0.8x[2] + 0.3x[1])),
         x -> VectorValue(ntuple(j -> a * 0.03j * sin(0.7x[1] + 0.2j + 0.4x[2]), Nσ)...),
         x -> VectorValue(ntuple(j -> a * 0.02 * (Nσ + 1 - j) * cos(0.6x[2] + 0.5x[1] - 0.1j), Nσ)...)],
        S.U)
    x = copy(get_free_dof_values(uh)); x .+= 0.01a .* randn(MersenneTwister(seed), length(x))
    return FEFunction(S.U, x)
end

"The derived fields exactly as global_residual builds them (works on transient fields too)."
function fields(prob, u, trian)
    η, Ux, Uy = u[1], u[2], u[3]
    d_cf = CellField(prob.h_bathy, trian); H = d_cf + η
    dhx = prob.flat_bed ? 0.0 * d_cf : alg_dx(d_cf); dhy = prob.flat_bed ? 0.0 * d_cf : alg_dy(d_cf)
    dHx = dhx + alg_dx(η); dHy = dhy + alg_dy(η)
    DU = alg_dx(Ux) + alg_dy(Uy); bf = dHx * Ux + dHy * Uy
    return (η = η, Ux = Ux, Uy = Uy, d_cf = d_cf, H = H, dHx = dHx, dHy = dHy, DU = DU, bf = bf,
            S = H * DU + bf)
end

"AD Jacobian (w.r.t. u, u̇ frozen) of a residual form `res(tu, v)` wrapped in a TransientCellField —
the configuration the time integrator differentiates, and the one the patch repairs."
ad_jac(S, res, uh, ud) =
    Gridap.FESpaces.jacobian(Gridap.FESpaces.FEOperator((x, v) -> res(TransientCellField(x, (ud,)), v),
                                                        S.U, S.V), uh)
rel(A, B) = norm(Matrix(A) - Matrix(B)) / max(norm(Matrix(B)), 1e-300)

# ---------------------------------------------------------------------------------------
println("\n  S0 — the fork carries the skeleton-AD patch")
chk("S0 DomainStyle(::Type{SkeletonCellFieldPair}) is defined",
    any(m -> occursin("Type{<:Gridap.CellData.SkeletonCellFieldPair", string(m.sig)),
        methods(Gridap.CellData.DomainStyle)))

# ---------------------------------------------------------------------------------------
println("\n  S1 — AD of the broken Class-III residual = the hand broken_class3_jacobian")
for (lbl, S, bedf, fb) in (("sloping 2-D mesh", setup(), bed, false),
                           ("x-periodic box (wrap facet)", setup(nx = 8, ny = 1, Lx = 4.0, Ly = 0.5,
                                                                periodic = true, bc = :wall), flat, true))
    pb = build_problem(vert; h_bathy = bedf, regime = :nonlinear, nl_pressure = true, flat_bed = fb,
                       model = S.model, quad_degree = QDEG)
    uh = state(S, 1.0; seed = 1); ud = FEFunction(S.U, zeros(num_free_dofs(S.U)))
    res(tu, v) = (f = fields(pb, tu, S.trian);
                  broken_class3_residual(pb, pb.skel[], f.d_cf, f.η, f.H, f.dHx, f.dHy, f.Ux, f.Uy,
                                         f.DU, f.S, f.bf, v[2], v[3], alg_dx(v[2]) + alg_dy(v[3]), S.dΩ))
    t0 = time()
    Jad = ad_jac(S, res, uh, ud)
    d_cf = CellField(pb.h_bathy, S.trian)
    Jh  = assemble_matrix((du, v) -> broken_class3_jacobian(pb, pb.skel[], d_cf, uh[1], uh[2], uh[3],
                                                            du[1], du[2], du[3], v[2], v[3], S.dΩ), S.U, S.V)
    e = rel(Jad, Jh)
    @printf("    %-28s |J_AD − J_hand|/|J_hand| = %.2e   (|J| = %.3e, %.0f s)\n", lbl, e, norm(Matrix(Jh)), time() - t0)
    chk("S1 [$lbl] AD ≡ hand to round-off", e < 1e-10)
end

# ---------------------------------------------------------------------------------------
println("\n  S2 — AD of the skeleton stabilisers = their own (linear) form")
let S = setup(), uh = state(S, 1.0; seed = 2), ud = FEFunction(S.U, zeros(num_free_dofs(S.U)))
    for (stab, ord) in ((:jumpgrad, 2), (:ghostvolume, 1))
        p = build_problem(vert; h_bathy = bed, regime = :nonlinear, nl_pressure = false, flat_bed = false)
        attach_skeleton!(p, S.model; cip_gamma_u = 0.3, cip_gamma_eta = 0.3, stabilization = stab,
                         cip_order = ord, p_u = PU, p_eta = PE, degree = QDEG)
        res(tu, v) = stab_contrib(p, p.skel[], tu[1], tu[2], tu[3], v[1], v[2], v[3])
        Jad = ad_jac(S, res, uh, ud)
        Jh  = assemble_matrix((du, v) -> stab_contrib(p, p.skel[], du[1], du[2], du[3], v[1], v[2], v[3]),
                              S.U, S.V)
        e = rel(Jad, Jh)
        @printf("    %-12s |J_AD − J_form|/|J_form| = %.2e\n", stab, e)
        chk("S2 [$stab] AD ≡ the stabiliser form to round-off", e < 1e-10)
    end
end

# ---------------------------------------------------------------------------------------
println("\n  S3 — transient run: use_ad=true reaches the hand-Jacobian answer (rule 17b)")
function box_run(use_ad)
    λ = 4.0; n = 8
    vert0 = vert; k = 2π / λ
    eta0 = x -> 0.05 * cos(k * x[1])
    ux0  = x -> VectorValue(ntuple(j -> 0.02j * cos(k * x[1]), Nσ)...)
    out  = mktempdir()
    diags, _, _ = setup_and_run(; M = 2, p_vertical = 1, p_u = PU, p_eta = PE,
        domain = ((0.0, λ), (0.0, λ / n)), partition = (n, 1), x_periodic = true,
        y_wall_bc = :wall, x_wall_bc = false, h_val = 3.5, T_wave = 1.6, A_wave = 0.05,
        sponge_wL = 0.0, sponge_wR = 0.0, sponge_wB = 0.0, sponge_wT = 0.0,
        eta0_func = eta0, ux0_func = ux0, T_final = 0.12, dt = 0.04,
        regime = :nonlinear, nl_pressure = true, flat_bed = true,
        stabilization = :ghostvolume, cip_gamma_u = 0.01, cip_gamma_eta = 0.01,
        solver_type = :theta, nl_tol = 1e-10, use_ad = use_ad,
        output_dir = out, save_every = 0, print_every = 10^6, check_every = 0, diag_every = -1)
    return [d.eta_max for d in diags]
end
let h = box_run(false), a = box_run(true)
    d = maximum(abs.(h .- a)) / maximum(abs.(h))
    @printf("    max η: hand %s   AD %s   rel diff %.2e\n", string(round.(h, sigdigits = 8)),
            string(round.(a, sigdigits = 8)), d)
    chk("S3 use_ad=true runs and matches the hand-Jacobian trajectory (rel < 1e-8)",
        length(a) == length(h) && d < 1e-8)
end

println("\n" * "="^70); @printf("  %d passed, %d failed\n", np, nf); println("="^70)
nf == 0 || error("test_skeleton_ad: $nf gate(s) failed")
