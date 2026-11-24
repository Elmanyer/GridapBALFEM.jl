# ==============================================================
#  test_nlp_distributed_parity.jl — the SEQUENTIAL↔DISTRIBUTED gate of the full nonlinear model
#                                    (V2_SOLVER_PLAN.md step 10)
#
#  Models 5–6 (nl_pressure=true) carry two skeleton terms: the Class-III facet terms of blocks C
#  and L and the mandatory first-order C⁰-CIP penalty (LaTeX §6.4). Both live on the interior
#  facets, and on a DistributedDiscreteModel each facet is held once, by the rank that owns it,
#  with the ghost layer supplying the other cell. A facet on a partition boundary that is lost,
#  doubled, or evaluated with the wrong neighbour changes the assembled residual. This gate
#  compares the two execution paths on the SAME state.
#
#  WHAT IS COMPARED. DOF numberings differ between the paths, so vectors are not compared entry by
#  entry. Instead every operator is paired with a fixed smooth test field v* (and, for the
#  Jacobians, a fixed direction δu*): r(v*) = R(u; v*) and J(v*, δu*) are scalars, identical in
#  exact arithmetic, and every facet, the partition-boundary ones included, contributes to them.
#  They are compared to a relative 1e-10:
#    P1  facet count                     (each interior facet exactly once)
#    P2  the C⁰-CIP term alone           J_CIP(u; v*)
#    P3  the Class-III facet term alone  F_III(u; v*)   (depth via depth_cf: ghost-inclusive)
#    P4  the full residual               R(u, u̇; v*)
#    P5  ∂R/∂u paired                    v*·(∂R/∂u)·δu*
#    P6  ∂R/∂u̇ paired                    v*·(∂R/∂u̇)·δu*
#  G0 (negative control): P2 and P3 must be materially nonzero, or the parity is vacuous.
#
#  ⚠ THE REFERENCE IS COMPUTED IN-PROCESS (every rank redundantly), never pinned as a constant.
#
#  RUN (default backend: PartitionedArrays' in-process debug backend, 4 "ranks" in one process —
#  the distributed code path without MPI):
#    julia --project=. test/test_nlp_distributed_parity.jl
#  and under real MPI:
#    PARITY_BACKEND=mpi ~/.julia/bin/mpiexecjl --project=. -n 4 julia --project=. \
#        test/test_nlp_distributed_parity.jl
# ==============================================================

using GridapBALFEM
using Gridap
using Gridap.ODEs
using GridapDistributed
using PartitionedArrays
using Printf

const BACKEND  = get(ENV, "PARITY_BACKEND", "debug")
const CPU_GRID = (2, 2)
const NX, NY   = 8, 6
const LX, LY   = 1.7, 1.1
const PU, PE   = 3, 2
const QDEG     = 2*max(PU, PE) + 2
const T0       = 0.0
const D0, AB   = 1.0, 0.15
const TOL      = 1e-10

vert = assemble_vertical_tensors(2, 1, [0.0, 0.728, 1.0])
const Nσ = vert.N_dof

state_funs(a) = [
    x -> a * 0.30 * cos(1.1x[1]) * cos(0.7x[2]),
    x -> VectorValue(ntuple(j -> a * (0.20 + 0.05j) * sin(0.9x[1] + 0.3j) * cos(0.5x[2]), Nσ)...),
    x -> VectorValue(ntuple(j -> a * (0.15 - 0.03j) * cos(0.8x[1]) * sin(0.6x[2] + 0.2j), Nσ)...),
]
test_funs = [
    x -> 0.7 * sin(1.3x[1] + 0.2) * cos(0.9x[2]),
    x -> VectorValue(ntuple(j -> (0.4 + 0.1j) * cos(1.2x[1] - 0.1j) * sin(1.7x[2] + 0.3), Nσ)...),
    x -> VectorValue(ntuple(j -> (0.3 - 0.05j) * sin(0.7x[1] + 0.4) * sin(1.9x[2] + 0.1j), Nσ)...),
]
dir_funs = [
    x -> 0.2 * cos(2.1x[1]) * sin(1.1x[2] + 0.5),
    x -> VectorValue(ntuple(j -> (0.1 + 0.02j) * sin(1.5x[1]) * cos(1.3x[2] - 0.2j), Nσ)...),
    x -> VectorValue(ntuple(j -> (0.12 - 0.01j) * cos(1.1x[1] + 0.3j) * cos(0.8x[2]), Nσ)...),
]

"All pairings of one model on one (sequential or distributed) mesh."
function pairings(model, trian, flat)
    bed  = flat ? (x -> D0) : (x -> D0 * (1 + AB * sin(1.3x[1]) * cos(0.9x[2])))
    U, V = build_fe_spaces(model, PU, Nσ; y_wall_bc = :wall, x_wall_bc = true, p_eta = PE)
    dΩh  = Measure(trian, QDEG)
    prob = build_problem(vert; h_bathy = bed, regime = :nonlinear, nl_pressure = true,
                         flat_bed = flat, model = model, quad_degree = QDEG, p_eta = PE)
    sk   = prob.skel[]
    uh   = interpolate_everywhere(state_funs(1.0), U)
    uth  = interpolate_everywhere(state_funs(0.6), U)
    du   = interpolate_everywhere(dir_funs, U)
    vs   = interpolate_everywhere(test_funs, V)
    tu   = TransientCellField(uh, (uth,))
    η, Ux, Uy = uh[1], uh[2], uh[3]
    q, Wx, Wy = vs[1], vs[2], vs[3]
    d_cf = GridapBALFEM.depth_cf(prob, trian)      # ghost-inclusive when distributed
    H    = d_cf + η
    dhx  = flat ? 0.0*GridapBALFEM.alg_dx(d_cf) : GridapBALFEM.alg_dx(d_cf)
    dhy  = flat ? 0.0*GridapBALFEM.alg_dy(d_cf) : GridapBALFEM.alg_dy(d_cf)
    dHx  = dhx + GridapBALFEM.alg_dx(η);  dHy = dhy + GridapBALFEM.alg_dy(η)
    DU   = GridapBALFEM.alg_dx(Ux) + GridapBALFEM.alg_dy(Uy)
    bf   = dHx*Ux + dHy*Uy
    S    = H*DU + bf
    DW   = GridapBALFEM.alg_dx(Wx) + GridapBALFEM.alg_dy(Wy)
    return (nf   = sk.nfacets,
            cip  = sum(stab_contrib(prob, sk, η, Ux, Uy, q, Wx, Wy)),
            f3   = sum(broken_class3_skeleton_contrib(prob, sk, H, dHx, dHy, Ux, Uy, S, bf, Wx, Wy, DW)),
            res  = sum(global_residual(T0, tu, vs, prob, trian, dΩh)),
            ju   = sum(jacobian_u(T0, tu, du, vs, prob, trian, dΩh)),
            jut  = sum(jacobian_u_t(T0, tu, du, vs, prob, trian, dΩh)))
end

results = Tuple{String,Bool}[]
gate(name, ok) = (push!(results, (name, ok)); @printf("  %s  %s\n", ok ? "PASS" : "FAIL", name))
rel(a, b) = abs(a - b) / max(abs(b), 1e-300)

backend = BACKEND == "mpi" ? with_mpi : with_debug
backend() do distribute
    ranks = distribute(LinearIndices((prod(CPU_GRID),)))
    main  = i_am_main(ranks)
    main && println("=== test_nlp_distributed_parity: backend=$BACKEND, cpu_grid=$CPU_GRID, mesh $(NX)x$(NY), Q$PU/Q$PE ===")
    for (tag, flat) in (("M5 flat bed", true), ("M6 variable bed", false))
        ms, ts = build_horizontal_model(((0.0, LX), (0.0, LY)), (NX, NY))
        s = pairings(ms, ts, flat)
        md, td = build_horizontal_model_distributed(ranks, CPU_GRID, (0.0, LX, 0.0, LY), (NX, NY))
        d = pairings(md, td, flat)
        if main
            @printf("  %s: facets %d/%d | CIP %.6e/%.6e | F_III %.6e/%.6e\n", tag, s.nf, d.nf, s.cip, d.cip, s.f3, d.f3)
            @printf("  %s: R %.12e/%.12e | ∂R/∂u %.6e/%.6e | ∂R/∂u̇ %.6e/%.6e\n", tag, s.res, d.res, s.ju, d.ju, s.jut, d.jut)
            gate("G0 $tag: CIP and facet pairings materially nonzero", abs(s.cip) > 1e-8 && abs(s.f3) > 1e-8)
            gate("P1 $tag: interior facet count",                     s.nf == d.nf)
            gate("P2 $tag: C⁰-CIP term       (rel $(@sprintf("%.1e", rel(d.cip, s.cip))))", rel(d.cip, s.cip) < TOL)
            gate("P3 $tag: Class-III facets  (rel $(@sprintf("%.1e", rel(d.f3, s.f3))))",   rel(d.f3, s.f3) < TOL)
            gate("P4 $tag: full residual     (rel $(@sprintf("%.1e", rel(d.res, s.res))))", rel(d.res, s.res) < TOL)
            gate("P5 $tag: ∂R/∂u paired      (rel $(@sprintf("%.1e", rel(d.ju, s.ju))))",   rel(d.ju, s.ju) < TOL)
            gate("P6 $tag: ∂R/∂u̇ paired      (rel $(@sprintf("%.1e", rel(d.jut, s.jut))))", rel(d.jut, s.jut) < TOL)
        end
    end
    if main
        np = count(last, results); nfail = length(results) - np
        println("="^70); @printf("  %d passed, %d failed\n", np, nfail); println("="^70)
    end
end
