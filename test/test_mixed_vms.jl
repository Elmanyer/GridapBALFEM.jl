# ==============================================================
#  test_mixed_vms.jl — the mixed formulation and the residual-based (VMS) stabilisation
#
#  Gates of markdown_files/MIXED_VMS_PLAN.md §6 (src/mixed.jl, src/vms.jl, src/jets.jl).
#
#  G0  NO LEAKAGE. With formulation=:primal and no VMS attached, the residual and both Jacobians
#      are those of the v2 solver: attaching a VMS record with all coefficients 0 leaves them
#      BITWISE unchanged, and the refusals of the plan (§1) hold.
#  G1  AUXILIARY ROWS. 𝖦, 𝖥 from mixed_consistent_ic satisfy their rows to round-off.
#  G2  MIXED ↔ BROKEN. On a smooth state the mixed physics rows and the v2 broken rows (its CIP
#      subtracted) approximate the same operator: their difference falls under refinement.
#  G3  WEAK ↔ STRONG. ∫ v·R_strong equals the assembled weak residual of the physics rows up to
#      the inter-element jumps of a C⁰ field: the difference falls under refinement. This checks
#      every sign of the strong residuals independently. Models 3, 4, 5, 6.
#  G4  JACOBIANS. The jet tangents of mixed_contrib and vms_contrib (∂/∂u and ∂/∂u̇) against a
#      central finite difference of the same assembled terms.
#  G5  KNOBS LIVE (rule 38d). Each VMS coefficient changes the residual on its own.
#  G6  ORTHOGONAL SUBGRID SCALES (oss=true). (a) the projections from oss_consistent_ic satisfy
#      their rows; (b) the MASS IS GALERKIN: ∂R/∂u̇ with OSS attached is bitwise the matrix without
#      any stabilisation, and the OSS residual does not read u̇; (c) the jet ∂/∂u of the OSS terms
#      and rows against central FD; (d) the OSS terms fall under refinement on a smooth state.
#
#  RUN:  julia --project=. test/test_mixed_vms.jl          (BALFEM_TEST_FAST=1 skips model 6)
# ==============================================================
using GridapBALFEM, Gridap, Gridap.ODEs, LinearAlgebra, SparseArrays, Printf, Random

println("="^70); println("  test_mixed_vms.jl — mixed formulation + residual-based stabilisation"); println("="^70)
np = 0; nf = 0
chk(n, c) = (global np, nf; c ? (println("  PASS  $n"); np += 1) : (println("  FAIL  $n"); nf += 1); flush(stdout))

const M, PV, PU, PE = 2, 1, 3, 2
const vert = assemble_vertical_tensors(M, PV, Vector{Float64}(resolve_cbdy(M, nothing, PV)))
const Nσ   = vert.N_dof
const QDEG = 2 * PU + 4
const FAST = get(ENV, "BALFEM_TEST_FAST", "0") == "1"
const zvv  = VectorValue(ntuple(_ -> 0.0, Nσ)...)

"x-periodic strip with solid y-walls. `bed = nothing` ⇒ flat."
function setup(; nx = 8, ny = 1, Lx = 4.0, Ly = 0.5, form = :primal, nlp = false, bed = nothing, oss = false)
    model, trian = build_horizontal_model(((0.0, Lx), (0.0, Ly)), (nx, ny); x_periodic = true)
    dΩ   = Measure(trian, QDEG)
    hf   = bed === nothing ? (x -> 3.5) : bed
    prob = build_problem(vert; h_bathy = hf, regime = :nonlinear, nl_pressure = nlp,
                         flat_bed = bed === nothing, model = model, quad_degree = QDEG, p_eta = PE,
                         formulation = form)
    U, V = build_fe_spaces(model, PU, Nσ; y_wall_bc = :wall, p_eta = PE, n_aux = mixed_n_aux(prob),
                           n_oss = oss ? 3 : 0)
    oss && attach_vms!(prob, trian; p_u = PU, oss = true)
    return (; model, trian, dΩ, prob, U, V, Lx, Ly)
end

"Periodic cosine bed with a transverse slope (smooth, so every derivative is analytic)."
bedfun(Lx, Ly) = x -> 3.5 + 0.20 * cos(2π * x[1] / Lx) + 0.05 * cos(π * x[2] / Ly)

"Smooth state, periodic in x, 𝖴y = 0 on the walls. Auxiliary fields zero."
function state(S; a = 1.0)
    k = 2π / S.Lx; m = π / S.Ly
    fs = Any[x -> a * (0.05 * cos(k * x[1]) + 0.02 * sin(2k * x[1] + 0.3)) * (1 + 0.3 * cos(m * x[2])),
             x -> VectorValue(ntuple(j -> a * (0.03j * sin(k * x[1] + 0.2j) + 0.01 * cos(2k * x[1] - 0.1j)) *
                                          (1 + 0.2 * cos(m * x[2])), Nσ)...),
             x -> VectorValue(ntuple(j -> a * 0.02 * sin(m * x[2]) * cos(k * x[1] - 0.1j), Nσ)...)]
    append!(fs, extra_field_zeros(S.prob))
    return interpolate_everywhere(fs, S.U)
end

"Smooth u̇ of the same kind."
function dstate(S; a = 1.0)
    k = 2π / S.Lx; m = π / S.Ly
    fs = Any[x -> a * 0.04 * sin(k * x[1] + 0.1) * (1 + 0.2 * cos(m * x[2])),
             x -> VectorValue(ntuple(j -> a * 0.05j * cos(k * x[1] + 0.3j) * (1 + 0.1 * cos(m * x[2])), Nσ)...),
             x -> VectorValue(ntuple(j -> a * 0.03 * sin(m * x[2]) * sin(k * x[1] + 0.2j), Nσ)...)]
    append!(fs, extra_field_zeros(S.prob))
    return interpolate_everywhere(fs, S.U)
end

"State with consistent auxiliary fields (mixed 𝖦, 𝖥, then OSS Π, δ). Returns the larger relative row residual."
function cstate(S; a = 1.0)
    u = state(S; a = a); rel = 0.0
    if S.prob.formulation === :mixed
        u, r1 = mixed_consistent_ic(u, S.prob, S.U, S.trian, S.dΩ); rel = max(rel, r1)
    end
    if has_oss(S.prob)
        u, r2 = oss_consistent_ic(u, S.prob, S.U, S.trian, S.dΩ); rel = max(rel, r2)
    end
    return u, rel
end

tcf(uh, ud) = TransientCellField(uh, (ud,))
resid(S, uh, ud) = assemble_vector(v -> global_residual(0.0, tcf(uh, ud), v, S.prob, S.trian, S.dΩ), S.V)
n3(S) = sum(num_free_dofs(S.U.spaces[i]) for i in 1:3)        # DOFs of the three physics fields

# ---------------------------------------------------------------------------------------
println("\n  G0 — no leakage, refusals")
let S = setup(nx = 6), uh = state(S), ud = dstate(S)
    r0 = resid(S, uh, ud)
    J0 = assemble_matrix((du, v) -> jacobian_u(0.0, tcf(uh, ud), du, v, S.prob, S.trian, S.dΩ), S.U, S.V)
    attach_vms!(S.prob, S.trian; p_u = PU, c_adv = 0.0, c_wave = 0.0, c_eta = 0.0)
    chk("all three coefficients 0 ⇒ no VMS record", !has_vms(S.prob))
    chk("residual bitwise unchanged", resid(S, uh, ud) == r0)
    J1 = assemble_matrix((du, v) -> jacobian_u(0.0, tcf(uh, ud), du, v, S.prob, S.trian, S.dΩ), S.U, S.V)
    chk("jacobian_u bitwise unchanged", J1 == J0)
    chk("primal problem carries formulation=:primal and no mixed record",
        S.prob.formulation === :primal && S.prob.mixed[] === nothing)
end
refused(f) = try f(); false catch; true end
let model = build_horizontal_model(((0.0, 4.0), (0.0, 0.5)), (4, 1); x_periodic = true)[1]
    bp(; kw...) = build_problem(vert; h_bathy = x -> 3.5, flat_bed = true, model = model,
                                quad_degree = QDEG, p_eta = PE, kw...)
    chk("formulation=:mixed without 𝓝 is refused", refused(() -> bp(nl_pressure = false, formulation = :mixed)))
    chk("an unknown formulation is refused", refused(() -> bp(nl_pressure = true, formulation = :hybrid)))
    pm = bp(nl_pressure = true, formulation = :mixed)
    chk("mixed problem: no skeleton record, boundary record present",
        pm.skel[] === nothing && pm.mixed[] !== nothing)
    chk("mixed problem refuses a skeleton penalty",
        refused(() -> attach_skeleton!(pm, model; cip_gamma_u = 0.3, cip_gamma_eta = 0.3, p_u = PU,
                                       p_eta = PE, degree = QDEG)))
    pp = bp(nl_pressure = true)
    trian = Triangulation(model)
    chk("VMS on the primal 𝓝 solver is refused", refused(() -> attach_vms!(pp, trian; p_u = PU)))
    pl = build_problem(vert; h_bathy = x -> 3.5, flat_bed = true, regime = :linear)
    chk("VMS on a linear model is refused", refused(() -> attach_vms!(pl, trian; p_u = PU)))
    chk("a negative VMS coefficient is refused",
        refused(() -> attach_vms!(bp(nl_pressure = false), trian; p_u = PU, c_adv = -1.0)))
end

# ---------------------------------------------------------------------------------------
println("\n  G1 — auxiliary rows at the consistent state")
let S = setup(nx = 8, ny = 2, form = :mixed, nlp = true, bed = bedfun(4.0, 0.5))
    _, rel = cstate(S)
    @printf("    relative auxiliary residual = %.2e\n", rel)
    chk("𝖦, 𝖥 satisfy their rows (≤ 1e-10)", rel ≤ 1e-10)
end

# ---------------------------------------------------------------------------------------
println("\n  G2 — mixed ↔ broken: the physics rows agree under refinement")
function g2(nx)
    Sp = setup(nx = nx, form = :primal, nlp = true)
    Sm = setup(nx = nx, form = :mixed,  nlp = true)
    up = state(Sp); um, _ = cstate(Sm)
    udp = dstate(Sp); udm = dstate(Sm)
    rp  = resid(Sp, up, udp)
    sk  = Sp.prob.skel[]
    cip = assemble_vector(v -> stab_contrib(Sp.prob, sk, up[1], up[2], up[3], v[1], v[2], v[3]), Sp.V)
    rm  = resid(Sm, um, udm)[1:n3(Sm)]
    #  scale: the Class-III contribution itself (𝓝 on minus 𝓝 off would mix in {3,6,7,8})
    return norm(rm - (rp - cip)) / norm(rp - cip)
end
let e = [g2(nx) for nx in (8, 16)]
    @printf("    |mixed − broken|/|broken|:  nx=8 %.3e   nx=16 %.3e   (ratio %.1f)\n", e[1], e[2], e[1] / e[2])
    chk("difference falls under refinement (ratio > 2)", e[1] / e[2] > 2)
end

# ---------------------------------------------------------------------------------------
println("\n  G3 — weak ↔ strong: ∫ v·R_strong against the assembled weak residual")
function g3(nx; form, nlp, bed)
    #  refine BOTH directions: the state varies in y, and a fixed ny leaves a mesh-independent
    #  inter-element jump in ∂ᵧ that the comparison would read as a mismatch
    S = setup(nx = nx, ny = nx ÷ 3, form = form, nlp = nlp, bed = bed)
    uh, _ = cstate(S); ud = dstate(S)
    rw = resid(S, uh, ud)[1:n3(S)]
    vc = GridapBALFEM.vms_ctx(S.prob, S.trian; p_u = PU)
    rs = assemble_vector(v -> vms_strong_contrib(S.prob, 0.0, tcf(uh, ud), v, S.trian, S.dΩ; vc = vc), S.V)[1:n3(S)]
    return norm(rw - rs) / norm(rw)
end
g3cases = Any[("model 3 (flat, 𝓝 off, primal)",     (form = :primal, nlp = false, bed = nothing)),
              ("model 4 (variable bed, 𝓝 off)",     (form = :primal, nlp = false, bed = bedfun(4.0, 0.5))),
              ("model 5 (flat, 𝓝 on, mixed)",       (form = :mixed,  nlp = true,  bed = nothing))]
FAST || push!(g3cases, ("model 6 (variable bed, 𝓝 on, mixed)", (form = :mixed, nlp = true, bed = bedfun(4.0, 0.5))))
for (name, kw) in g3cases
    e = [g3(nx; kw...) for nx in (6, 12)]
    @printf("    %-40s nx=6 %.3e   nx=12 %.3e   (ratio %.1f)\n", name, e[1], e[2], e[1] / e[2])
    chk("$name: weak − strong falls under refinement (ratio > 2, fine < 5e-2)", e[1] / e[2] > 2 && e[2] < 5e-2)
end

# ---------------------------------------------------------------------------------------
println("\n  G4 — jet Jacobians against central finite differences")
"Relative error of the assembled Jacobian `J(u)` of the vector function `r(u)` in a random direction."
function fd_err(S, rfun, Jfun, x0; ε = 1e-6, seed = 7)
    δ  = randn(MersenneTwister(seed), length(x0)); δ ./= norm(δ)
    fe(x) = FEFunction(S.U, x)
    fd = (rfun(fe(x0 .+ ε .* δ)) .- rfun(fe(x0 .- ε .* δ))) ./ (2ε)
    Jδ = Jfun(fe(x0)) * δ
    return norm(Jδ - fd) / norm(fd)
end
function g4(; form, nlp, bed, label)
    S  = setup(nx = 4, ny = 2, form = form, nlp = nlp, bed = bed)
    uh, _ = cstate(S); ud = dstate(S)
    x0 = copy(get_free_dof_values(uh)); xd = copy(get_free_dof_values(ud))
    #  a rough perturbation: the Jacobian must be exact away from the smooth manifold too
    x0 .+= 1e-3 .* randn(MersenneTwister(3), length(x0))
    if form === :mixed
        e = fd_err(S, u -> assemble_vector(v -> mixed_contrib(S.prob, u, v, S.trian, S.dΩ), S.V),
                      u -> assemble_matrix((du, v) -> mixed_contrib(S.prob, u, v, S.trian, S.dΩ; du = du), S.U, S.V), x0)
        @printf("    %-22s mixed rows      ∂/∂u : %.2e\n", label, e)
        chk("$label: mixed_contrib Jacobian (≤ 1e-6)", e ≤ 1e-6)
    end
    attach_vms!(S.prob, S.trian; p_u = PU, c_adv = 1.0, c_wave = 1.0, c_eta = 1.0)
    eu = fd_err(S, u -> assemble_vector(v -> vms_contrib(S.prob, 0.0, tcf(u, ud), v, S.trian, S.dΩ), S.V),
                   u -> assemble_matrix((du, v) -> vms_contrib(S.prob, 0.0, tcf(u, ud), v, S.trian, S.dΩ; du = du), S.U, S.V), x0)
    u0 = FEFunction(S.U, x0)
    et = fd_err(S, w -> assemble_vector(v -> vms_contrib(S.prob, 0.0, tcf(u0, w), v, S.trian, S.dΩ), S.V),
                   w -> assemble_matrix((dw, v) -> vms_contrib(S.prob, 0.0, tcf(u0, w), v, S.trian, S.dΩ; dut = dw), S.U, S.V), xd)
    @printf("    %-22s VMS             ∂/∂u : %.2e    ∂/∂u̇ : %.2e\n", label, eu, et)
    chk("$label: vms_contrib ∂/∂u (≤ 1e-6)", eu ≤ 1e-6)
    chk("$label: vms_contrib ∂/∂u̇ (≤ 1e-6)", et ≤ 1e-6)
end
g4(form = :primal, nlp = false, bed = nothing,          label = "model 3")
g4(form = :mixed,  nlp = true,  bed = nothing,          label = "model 5 mixed")
FAST || g4(form = :mixed, nlp = true, bed = bedfun(4.0, 0.5), label = "model 6 mixed")

# ---------------------------------------------------------------------------------------
println("\n  G5 — each VMS coefficient is live")
let S = setup(nx = 4, ny = 2, form = :primal, nlp = false), uh = state(S), ud = dstate(S)
    vv(; kw...) = (attach_vms!(S.prob, S.trian; p_u = PU, c_adv = 0.0, c_wave = 0.0, c_eta = 0.0, kw...);
                   has_vms(S.prob) ? norm(assemble_vector(v -> vms_contrib(S.prob, 0.0, tcf(uh, ud), v, S.trian, S.dΩ), S.V)) : 0.0)
    na = vv(c_adv = 1.0); nw = vv(c_wave = 1.0); ne = vv(c_eta = 1.0)
    @printf("    |VMS term|:  C_a alone %.3e   C_w alone %.3e   C_η alone %.3e\n", na, nw, ne)
    chk("C_a alone contributes", na > 0)
    chk("C_w alone contributes", nw > 0)
    chk("C_η alone contributes", ne > 0)
    n2 = vv(c_wave = 2.0)
    chk("the C_w term is linear in C_w (state-independent τ_w)", isapprox(n2, 2nw; rtol = 1e-12))
end

# ---------------------------------------------------------------------------------------
println("\n  G6 — orthogonal subgrid scales")
function g6(; form, nlp, bed, label)
    S = setup(nx = 4, ny = 2, form = form, nlp = nlp, bed = bed, oss = true)
    uh, rel = cstate(S); ud = dstate(S)
    @printf("    %-22s relative projection-row residual = %.2e\n", label, rel)
    chk("$label: OSS projections satisfy their rows (≤ 1e-10)", rel ≤ 1e-10)
    #  (b) the mass: ∂R/∂u̇ with OSS == without any stabilisation, bitwise; no u̇ in the OSS residual
    Jt(S) = assemble_matrix((dw, v) -> jacobian_u_t(0.0, tcf(uh, ud), dw, v, S.prob, S.trian, S.dΩ), S.U, S.V)
    Jo  = Jt(S)
    rec = S.prob.vms[]; S.prob.vms[] = nothing
    Jn  = Jt(S)
    S.prob.vms[] = rec
    chk("$label: ∂R/∂u̇ with OSS is bitwise the unstabilised mass", Jo == Jn)
    ro(w) = assemble_vector(v -> vms_contrib(S.prob, 0.0, tcf(uh, w), v, S.trian, S.dΩ), S.V)
    chk("$label: the OSS residual does not read u̇", ro(ud) == ro(dstate(S; a = 3.7)))
    #  (c) exact ∂/∂u (terms + projection rows), away from the smooth manifold
    x0 = copy(get_free_dof_values(uh)); x0 .+= 1e-3 .* randn(MersenneTwister(3), length(x0))
    eu = fd_err(S, u -> assemble_vector(v -> vms_contrib(S.prob, 0.0, tcf(u, ud), v, S.trian, S.dΩ), S.V),
                   u -> assemble_matrix((du, v) -> vms_contrib(S.prob, 0.0, tcf(u, ud), v, S.trian, S.dΩ; du = du), S.U, S.V), x0)
    @printf("    %-22s OSS             ∂/∂u : %.2e\n", label, eu)
    chk("$label: OSS vms_contrib ∂/∂u (≤ 1e-6)", eu ≤ 1e-6)
end
g6(form = :primal, nlp = false, bed = nothing,          label = "model 3")
g6(form = :mixed,  nlp = true,  bed = nothing,          label = "model 5 mixed")
FAST || g6(form = :mixed, nlp = true, bed = bedfun(4.0, 0.5), label = "model 6 mixed")
#  (d) weak consistency: the OSS terms (physics rows) relative to the weak residual fall under
#  refinement; the ASGS terms are printed for comparison (they fall like τ ~ h only)
function g6d(nx; oss)
    S = setup(nx = nx, ny = nx ÷ 3, form = :primal, nlp = false, oss = oss)
    oss || attach_vms!(S.prob, S.trian; p_u = PU)
    uh, _ = cstate(S); ud = dstate(S)
    rs = assemble_vector(v -> vms_contrib(S.prob, 0.0, tcf(uh, ud), v, S.trian, S.dΩ), S.V)[1:n3(S)]
    rw = resid(S, uh, ud)[1:n3(S)] - rs
    return norm(rs) / norm(rw)
end
let eo = [g6d(nx; oss = true) for nx in (6, 12)], ea = [g6d(nx; oss = false) for nx in (6, 12)]
    @printf("    model 3, |stab|/|weak|:  OSS nx=6 %.3e  nx=12 %.3e (ratio %.1f)   ASGS %.3e → %.3e (ratio %.1f)\n",
            eo[1], eo[2], eo[1] / eo[2], ea[1], ea[2], ea[1] / ea[2])
    chk("model 3: the OSS terms fall under refinement (ratio > 2)", eo[1] / eo[2] > 2)
end

println("\n", "="^70)
@printf("  test_mixed_vms.jl: %d passed, %d failed\n", np, nf)
println("="^70)
nf == 0 || error("test_mixed_vms.jl: $nf gate(s) failed")
