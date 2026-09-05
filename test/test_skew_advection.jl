# ==============================================================
#  test_skew_advection.jl — the energy identity behind the skew-symmetric
#                           (energy-consistent) advection correction
#
#  WHAT IT PROVES. The advection block assembled by `global_residual` has an EXACT
#  energy production, derived in building_files/SKEW_SYMMETRIC_ADVECTION_PLAN.md §1:
#
#      ★★   n(U;U,U) = −½∫_Ω ∇·(Hū)(Σ_ij M_ij u_i·u_j) dΩ
#                      + ½∮_∂Ω H Σ𝓜_ikj (u_k·n)(u_j·u_i) dΓ        for ANY discrete U
#
#  so the operator is NOT energy-neutral: its production is exactly the continuity
#  defect, weighted by the kinetic-energy density. Discretely that defect cannot
#  cancel, because the cancellation needs the continuity equation tested against
#  q = ½Σ M_ij u_i·u_j — a QUARTIC in the velocity, which is not in the η space.
#  The correction adds that defect back as a consistent term.
#
#  ★★ rests on the σ-tensor identity ½(𝓖_ikj + 𝓖_jki) = ½𝓜_ikj − ½Φ_k M_ij, which is
#  gated separately (and on five vertical bases) in test_vertical.jl. THIS file gates
#  the horizontal half — the symmetrisation, the integration by parts, and the
#  assembled operator — on real discrete states.
#
#  ⚠ WHY THE BOUNDARY TERM IS IN THE ASSERTION AND NOT ASSUMED AWAY.
#  The first version of this test built its state with `interpolate_everywhere` and
#  asserted pi_res ≈ 0 in a walled basin, reasoning that u·n = 0 kills the flux. It
#  FAILED at ratio 2.885 — because `interpolate_everywhere` writes the DIRICHLET dofs
#  from the function too (that is what "everywhere" means; `interpolate` is the one
#  that takes them from the space). The state therefore never satisfied the walls and
#  the flux was large and real. Both forms are now tested, which is strictly stronger:
#  G5 pins the general identity INCLUDING the flux, G6 the closed-basin corollary.
#
#  RUN:  julia --project=. test/test_skew_advection.jl
# ==============================================================

using GridapBALFEM
using Gridap, Gridap.ODEs, Gridap.FESpaces, Gridap.Algebra
using LinearAlgebra, Printf

const A = GridapBALFEM

println("=" ^ 74)
println("  test_skew_advection.jl — the advection energy identity ★★")
println("=" ^ 74)

n_pass = 0; n_fail = 0
function check(name, cond, extra = "")
    global n_pass, n_fail
    cond ? (println("  PASS  $name $extra"); n_pass += 1) :
           (println("  FAIL  $name $extra"); n_fail += 1)
    flush(stdout)
end

vert     = assemble_vertical_tensors(2, 1, [0.0, 0.728, 1.0])
const Nσ = vert.N_dof
#  A SLOPING bed by default: ∇h ≠ 0 puts the bed-slope half of ∇·(H u_k) into the
#  identity, which a flat bed cannot test (rule: a flat-bed regression can never
#  test ∇h code). The flat bed is checked too, as the control.
h_sine(x) = 2.5 * (1 + 0.2*sin(1.3*x[1]))
h_flat(x) = 2.5

model, trian = build_horizontal_model(((0.0, 1.7), (0.0, 1.1)), (4, 4))
U, V = build_fe_spaces(model, 2, Nσ; y_wall_bc = :wall, x_wall_bc = true, p_eta = 1)
dΩh  = Measure(trian, 8)
Γ    = BoundaryTriangulation(model)
dΓ   = Measure(Γ, 8)
nΓ   = get_normal_vector(Γ)

zvv   = VectorValue(ntuple(_ -> 0.0, Nσ)...)
state = [ x -> 0.30*cos(1.1x[1])*cos(0.7x[2]),
          x -> VectorValue(ntuple(j -> (0.20+0.05j)*sin(0.9x[1]+0.3j)*cos(0.5x[2]), Nσ)...),
          x -> VectorValue(ntuple(j -> (0.15-0.03j)*cos(0.8x[1])*sin(0.6x[2]+0.2j), Nσ)...) ]

mkprob(hb; skew, regime = :nonlinear) =
    build_problem(vert; g = GridapBALFEM.g, h_bathy = hb, regime = regime,
                  nl_pressure = :none, flat_bed = false, skew_advection = skew,
                  mu_sponge = (x -> 0.0), wm_src = ((x, t) -> 0.0))

resvec(prob, uh, uth, t = 0.37) =
    assemble_vector(v -> A.global_residual(t, TransientCellField(uh, (uth,)), v,
                                           prob, trian, dΩh), V)

# ---------------------------------------------------------------------------
#  G1–G4 — the switch itself
# ---------------------------------------------------------------------------
uh  = interpolate_everywhere(state, U)
uth = interpolate_everywhere([x -> 0.55*state[1](x), x -> 0.55*state[2](x),
                              x -> 0.55*state[3](x)], U)

r_off  = resvec(mkprob(h_sine; skew = false), uh, uth)
r_on   = resvec(mkprob(h_sine; skew = true),  uh, uth)
r_off2 = resvec(mkprob(h_sine; skew = false), uh, uth)
check("G1  skew=false reproducible (bit-identical)", r_off == r_off2)

rel = norm(r_on - r_off) / norm(r_off)
check("G2  skew=true CHANGES the residual (the knob is live, not dead)",
      rel > 1e-8, @sprintf("(‖Δr‖/‖r‖ = %.3e)", rel))

#  The correction is CUBIC in the state, so it must vanish identically at rest —
#  and the discrete rest state is force-free by construction (gravity's baseline).
z = interpolate_everywhere([x -> 0.0, x -> zvv, x -> zvv], U)
check("G3  rest state untouched by the correction",
      resvec(mkprob(h_sine; skew = true), z, z) ==
      resvec(mkprob(h_sine; skew = false), z, z))

#  The correction lives INSIDE the advection block, which a linear model does not
#  assemble. It must therefore be inert — and build_problem warns, so a launcher
#  cannot set it on a linear run and read the null result as a treatment effect.
check("G4  linear regime inert (no advection block to correct)",
      resvec(mkprob(h_sine; skew = false, regime = :linear), uh, uth) ==
      resvec(mkprob(h_sine; skew = true,  regime = :linear), uh, uth))

# ---------------------------------------------------------------------------
#  G5/G6 — ★★ itself, measured on the assembled operator
# ---------------------------------------------------------------------------
"Return (pi_adv, pi_res, boundary_flux) for state `u` under bathymetry `hb`."
function production(hb, u)
    prob = mkprob(hb; skew = false)
    η, Ux, Uy = u[1], u[2], u[3]
    d_cf = CellField(hb, trian);  H = d_cf + η
    dHx  = A.alg_dx(d_cf) + A.alg_dx(η);   dHy = A.alg_dy(d_cf) + A.alg_dy(η)
    Sk   = H*(A.alg_dx(Ux) + A.alg_dy(Uy)) + (dHx*Ux + dHy*Uy)     # ∇·(H u_k)
    TMx  = A.alg_outer(Ux, A.alg_dx(Ux)) + A.alg_outer(Uy, A.alg_dy(Ux))
    TMy  = A.alg_outer(Ux, A.alg_dx(Uy)) + A.alg_outer(Uy, A.alg_dy(Uy))
    ke2  = (Ux ⋅ A.alg_mul(prob.Mv, Ux)) + (Uy ⋅ A.alg_mul(prob.Mv, Uy))
    dHu  = A.alg_dot(prob.Φ, Sk)                                    # ∇·(H ū)

    #  n(U;U,U): the residual's advection block, tested with W := U.
    pi_adv = sum(∫( H*((A.alg_dc3(prob.M3, TMx)) ⋅ Ux)
                  + H*((A.alg_dc3(prob.M3, TMy)) ⋅ Uy)
                  + ((A.alg_dc3(prob.G3, A.alg_outer(Sk, Ux))) ⋅ Ux)
                  + ((A.alg_dc3(prob.G3, A.alg_outer(Sk, Uy))) ⋅ Uy) ) * dΩh)
    pi_res = pi_adv + sum(∫( 0.5*dHu*ke2 ) * dΩh)

    #  ½∮ H Σ𝓜_ikj (u_k·n)(u_j·u_i) dΓ.  𝓜 is fully symmetric, so contracting
    #  P_ij = u_i·u_j over the trailing two indices leaves the k (flux) index.
    Hb  = CellField(hb, Γ) + η
    MP  = A.alg_dc3(prob.M3, A.alg_outer(Ux, Ux) + A.alg_outer(Uy, Uy))
    nx  = nΓ ⋅ VectorValue(1.0, 0.0);  ny = nΓ ⋅ VectorValue(0.0, 1.0)
    flux = sum(∫( Hb*( nx*(Ux ⋅ MP) + ny*(Uy ⋅ MP) ) ) * dΓ)
    return pi_adv, pi_res, flux
end

println()
for (nm, hb) in (("sloping bed", h_sine), ("flat bed", h_flat))
    pa, pr, fx = production(hb, uh)
    err = abs(pr - 0.5*fx) / max(abs(pa), 1e-300)
    @printf("    %-12s  pi_adv = %+.8e   pi_res = %+.8e   ½∮ = %+.8e\n", nm, pa, pr, 0.5*fx)
    check("G5  ★★ on a general state, $nm: pi_res = ½∮flux", err < 1e-10,
          @sprintf("(rel %.2e)", err))
end

#  Closed-basin corollary. `interpolate` (NOT `interpolate_everywhere`) takes the
#  Dirichlet dofs from the SPACE, so u·n = 0 on all four walls and the flux term
#  vanishes identically — leaving pi_adv equal to the continuity defect alone.
uw = interpolate(state, U)
pa, pr, fx = production(h_sine, uw)
@printf("\n    closed basin  pi_adv = %+.8e   pi_res = %+.8e   ½∮ = %+.8e\n", pa, pr, 0.5*fx)
check("G6  closed basin: the boundary flux vanishes", abs(fx) < 1e-10*max(abs(pa), 1.0),
      @sprintf("(∮ = %.2e)", fx))
check("G6  closed basin: pi_adv is EXACTLY the continuity defect (pi_res = 0)",
      abs(pr) < 1e-10*max(abs(pa), 1.0), @sprintf("(rel %.2e)", abs(pr)/abs(pa)))

println()
println("=" ^ 74)
@printf("  Results: %d PASS,  %d FAIL\n", n_pass, n_fail)
println("=" ^ 74)
n_fail > 0 ? error("test_skew_advection: $n_fail failed!") :
             println("  The advection energy identity ★★ holds on the assembled operator.")
