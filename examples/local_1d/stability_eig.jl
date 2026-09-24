# ==============================================================
#  stability_eig.jl — OPTION C: frozen-state linear stability of the semi-discrete system
#
#  THE QUESTION (CLAUDE.md §5.2e, 2026-09-24). Is the `:full` instability a property of the
#  CONTINUOUS model at finite amplitude (Hadamard ill-posedness: growth unbounded in k), or of
#  the DISCRETISATION (growth confined to each scheme's grid scale)? Yang & Liu run our exact
#  case (Reg03: A=0.10, T=1.6, kd=5.5) for 190 T with 5-point FD + RK4 and no filter, so a
#  stable discretisation of the model exists at least for their resolution.
#
#  THE METHOD (Fuhrman, Bingham, Madsen & Thomsen 2004, IJNMF 45, matrix-based analysis).
#  Linearise the semi-discrete system M(u) u̇ + r(u) = 0 about a frozen finite-amplitude state
#  u*:   M* δu̇ + J* δu = 0,   λ = eig(−M*⁻¹ J*).
#  Both matrices are EXACT (AD of the assembled residual), so no quasi-Newton gap enters.
#  Mixed layout: the auxiliary rows carry no time derivative, so 𝖦 is eliminated exactly by
#  its Schur complement (J_aa is the auxiliary Gram matrix):
#      M_pp δṗ + (J_pp − J_pa J_aa⁻¹ J_ap) δp = 0.
#
#  THE STATE. A periodic domain exactly ONE model wavelength long, carrying the solver's own
#  discrete linear eigenmode (waveinput.jl `_component_uamp`, :model) at amplitude A. One
#  wavelength EXCLUDES Benjamin–Feir sidebands by construction, and a Stokes wave at κa≈0.16
#  is superharmonically stable, so in the continuum no bounded-k growth is expected here at
#  all. u̇* is the same wave's time derivative; 𝖦* is solved from its own constraint.
#
#  THE READING. Each eigenpair is tagged with the dominant wavenumber of its velocity field.
#    * A = 0: every eigenvalue must be imaginary (linear core non-dissipative, no sponge) —
#      the sanity gate. A real part here is a defect in the harness or the linear core.
#    * σ_max vs dx at fixed A: bounded ⇒ well-posed; growing like dx^{-α} ⇒ the growth lives
#      at the grid scale.
#    * σ(k) across pairings: if the curves AGREE at resolved k and keep rising with k, the
#      growth is a property of the continuous model; if each scheme's growth sits at its own
#      grid scale and the curves disagree, it is discrete.
#  ⚠ FROZEN-STATE CAVEAT: the base wave is time-dependent, so λ is an instantaneous rate. It is
#  the right tool for the high-k question (those modes evolve much faster than the carrier),
#  not for O(ω) physics. And the periodic domain has no inflow, relaxation zone or sponge — a
#  boundary-driven mode will NOT appear here.
#
#  USE (from julia-mcp, env = this project):
#      include("examples/local_1d/stability_eig.jl")
#      r = stability_case(pu=2, pe=1, paux=2, model=:mixed, A=0.10, ncell=16)
# ==============================================================
using GridapBALFEM, Gridap, Gridap.TensorValues, LinearAlgebra, Printf, SparseArrays
using Gridap.ODEs: TransientCellField
import GridapBALFEM: _component_uamp, model_wavenumber

const SE_D, SE_T, SE_G = 3.5, 1.6, 9.81

"Vertical tensors for P1LFE-2 on the published mesh (the production basis)."
se_vert() = assemble_vertical_tensors(2, 1, Vector{Float64}(resolve_cbdy(2, nothing, 1)))

"""
    stability_case(; pu, pe, paux, model, A, ncell, t_phase=0.0, keep_vecs=false)

`model` ∈ (:none, :native, :mixed) — :mixed is `nl_pressure=:full` on the projection-free
5-field layout with c3_mask=(true,false), exactly the configuration of the c3*_mixed runs.
`ncell` = cells per wavelength. Returns a NamedTuple with the eigenvalues and, per eigenpair,
the dominant wavenumber of its 𝖴x field.
"""
function stability_case(; pu::Int, pe::Int = pu - 1, paux::Int = pu, model::Symbol,
                          A::Float64, ncell::Int, vert = se_vert(), quad_extra::Int = 0)
    Nσ    = vert.N_dof
    ω     = 2π / SE_T
    k     = model_wavenumber(vert, ω, SE_D, SE_G)
    L     = 2π / k
    dx    = L / ncell
    Ly    = dx
    uamp  = _component_uamp(vert, Float64[], :model, 1.0, ω, k, SE_D, SE_G)   # per unit A

    dm    = CartesianDiscreteModel((0.0, L, 0.0, Ly), (ncell, 1); isperiodic = (true, false))
    trian = Triangulation(dm)
    dΩh   = Measure(trian, 2 * max(pu, pe) + 2 + quad_extra)
    mixed = model === :mixed
    naux  = mixed ? 2 : 0
    U, V  = build_fe_spaces(dm, pu, Nσ; y_wall_bc = :wall, x_wall_bc = false, p_eta = pe,
                            n_aux = naux, p_aux = paux)

    nlp  = model === :none ? :none : model === :native ? :native : :full
    prob = build_problem(vert; h_bathy = (x -> SE_D), regime = :nonlinear, nl_pressure = nlp,
                               flat_bed = true, c3_mask = (true, false))

    zv   = VectorValue(ntuple(_ -> 0.0, Nσ)...)
    eta  = x -> A * cos(k * x[1])
    ux   = x -> VectorValue(ntuple(j -> A * uamp[j] * cos(k * x[1]), Nσ)...)
    etat = x -> A * ω * sin(k * x[1])
    uxt  = x -> VectorValue(ntuple(j -> A * uamp[j] * ω * sin(k * x[1]), Nσ)...)
    zf   = x -> zv
    uh   = interpolate_everywhere([eta,  ux,  zf, fill(zf, naux)...], U)
    uth  = interpolate_everywhere([etat, uxt, zf, fill(zf, naux)...], U)

    if mixed
        Γ  = BoundaryTriangulation(dm); dΓ = Measure(Γ, 10); nΓ = get_normal_vector(Γ)
        res = (t, u, v) -> global_residual_mixed(t, u, v, prob, trian, dΩh, dΓ, nΓ)
    else
        res = (t, u, v) -> global_residual(t, u, v, prob, trian, dΩh)
    end

    # field DOF ranges (ConsecutiveMultiFieldStyle: fields are contiguous, in order)
    nf   = [num_free_dofs(U.spaces[i]) for i in 1:length(U.spaces)]
    offs = cumsum([0; nf])
    rng(i) = (offs[i] + 1):offs[i + 1]
    ip   = 1:offs[4]                                   # physics block [η, 𝖴x, 𝖴y]
    ia   = mixed ? ((offs[4] + 1):offs[end]) : (1:0)

    #  ⚠ JACOBIANS BY COLOURED CENTRAL FINITE DIFFERENCES OF THE ASSEMBLED RESIDUAL, NOT AD.
    #  AD-compiling this residual took >60 min even at :none and the session was killed
    #  (2026-09-24). The FD Jacobian differentiates the SAME assembled residual the solver
    #  integrates, to ~1e-10 relative, reusing the residual path the solver already compiles.
    #  Its accuracy is GATED below: FD ∂R/∂u̇ vs the hand `jacobian_u_t` (exact, rule 5).
    rvec(x, xt) = assemble_vector(v -> res(0.0, TransientCellField(FEFunction(U, x),
                                                                    (FEFunction(U, xt),)), v), V)
    x0  = copy(get_free_dof_values(uh)); xt0 = copy(get_free_dof_values(uth))
    pat = se_pattern(U, V, dΩh, Nσ, naux)
    col = se_colour(pat)

    # 𝖦* consistent with (η*, 𝗎*): the aux rows are linear in 𝖦 with Gram block J_aa.
    if mixed
        J0 = se_fdjac(x -> rvec(x, xt0), x0, pat, col)
        r0 = rvec(x0, xt0)
        x0[ia] .-= Matrix(J0[ia, ia]) \ r0[ia]
        r1 = rvec(x0, xt0)
        aux_res = norm(r1[ia]) / max(norm(r0[ia]), 1e-300)
        uh = FEFunction(U, x0)
    else
        aux_res = 0.0
    end

    J = se_fdjac(x -> rvec(x, xt0), x0, pat, col)
    M = se_fdjac(w -> rvec(x0, w), xt0, pat, col)
    #  harness gate: FD ∂R/∂u̇ against the hand Jacobian (exact by rule 5), physics block only
    tu  = TransientCellField(uh, (uth,))
    Mh  = assemble_matrix((dut, v) -> jacobian_u_t(0.0, tu, dut, v, prob, trian, dΩh), U, V)
    m_gate = norm(Matrix(M)[ip, ip] - Matrix(Mh)[ip, ip]) / norm(Matrix(Mh)[ip, ip])

    Jd = Matrix(J); Md = Matrix(M)
    Jr = mixed ? Jd[ip, ip] - Jd[ip, ia] * (Jd[ia, ia] \ Jd[ia, ip]) : Jd
    Mr = Md[ip, ip]
    mix_leak = mixed ? norm(Md[ip, ia]) + norm(Md[ia, :]) : 0.0     # must be 0: no aux in u̇
    E  = eigen(-(Mr \ Jr))
    λ  = E.values

    # Dominant wavenumber of each eigenvector's 𝖴x field. A precomputed linear map
    # (𝖴x free DOFs → values at fine quadrature points) followed by a quadrature-weighted
    # Fourier transform, energy summed over the Nσ vertical components — all eigenvectors at
    # once as matrix products, never a point search.
    Bj, Fk, kgrid = se_fourier_map(U.spaces[2], trian, pu, ncell, L, Nσ)
    Vux  = E.vectors[rng(2), :]
    Pw   = zeros(length(kgrid), length(λ))
    for j in 1:Nσ
        Pw .+= abs2.(Fk * (Bj[j] * Vux))
    end
    kdom = [kgrid[argmax(view(Pw, :, m))] for m in 1:length(λ)]
    # fraction of each eigenvector's weight on 𝖴y (transverse, spurious in a 1-D flume)
    nrm2(v) = sum(abs2, v)
    uyfr = [nrm2(E.vectors[rng(3), m]) / max(nrm2(E.vectors[:, m]), 1e-300) for m in 1:length(λ)]

    σ    = real.(λ)
    imax = argmax(σ)
    return (pu = pu, pe = pe, paux = paux, model = model, A = A, ncell = ncell, dx = dx,
            L = L, k0 = k, ndof = length(ip), naux = length(ia), λ = λ, kdom = kdom,
            uyfr = uyfr, σmax = σ[imax], ωmax = abs(imag(λ[imax])), kmax = kdom[imax],
            uymax = uyfr[imax], rho = maximum(abs, λ), aux_res = aux_res, mix_leak = mix_leak,
            m_gate = m_gate, ncolour = maximum(col))
end

"""
Structural sparsity of the residual Jacobian: every residual term is a cell (or wall-face)
integral, so the pattern is the FE cell connectivity across ALL field pairs. Assembled from a
full-coupling bilinear form whose integrand is nonzero for every basis pair sharing a cell.
"""
function se_pattern(U, V, dΩh, Nσ, naux)
    one_ = VectorValue(ntuple(_ -> 1.0, Nσ)...)
    s(u) = u[1] + (u[2] ⋅ one_) + (u[3] ⋅ one_) +
           (naux ≥ 1 ? (u[4] ⋅ one_) : 0.0 * u[1]) + (naux ≥ 2 ? (u[5] ⋅ one_) : 0.0 * u[1])
    A = assemble_matrix((du, v) -> ∫(s(du) * s(v))dΩh, U, V)
    return spones(A)
end
spones(A::SparseMatrixCSC) = SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, ones(length(A.nzval)))

"Greedy column colouring: two columns share a colour only if they share no nonzero row."
function se_colour(P::SparseMatrixCSC)
    n = size(P, 2); col = zeros(Int, n)
    Pt = sparse(P')                                   # rows → columns
    for j in 1:n
        used = Set{Int}()
        for r in P.rowval[P.colptr[j]:(P.colptr[j + 1] - 1)]
            for c in Pt.rowval[Pt.colptr[r]:(Pt.colptr[r + 1] - 1)]
                col[c] > 0 && push!(used, col[c])
            end
        end
        c = 1; while c in used; c += 1; end
        col[j] = c
    end
    return col
end

"Central-difference Jacobian of f at x0 on pattern P with column colouring `col`."
function se_fdjac(f, x0, P::SparseMatrixCSC, col; h = 1e-6)
    J = similar(P); fill!(J.nzval, 0.0)
    for c in 1:maximum(col)
        js = findall(==(c), col)
        e  = zeros(length(x0)); e[js] .= h
        d  = (f(x0 .+ e) .- f(x0 .- e)) ./ (2h)
        for j in js, k in P.colptr[j]:(P.colptr[j + 1] - 1)
            J.nzval[k] = d[P.rowval[k]]
        end
    end
    return J
end

"""
Linear map from a VectorValue{Nσ} FE space's free DOFs to its values at fine quadrature points
(one dense matrix per vertical component), and the quadrature-weighted Fourier matrix at the
periodic wavenumbers 2πn/L, n = 0 … pu·ncell (the node Nyquist).
"""
function se_fourier_map(Ux, trian, pu, ncell, L, Nσ)
    dΩs = Measure(trian, 2 * pu + 6)
    pts = get_cell_points(dΩs)
    xq  = [p[1] for cellp in pts.cell_phys_point for p in cellp]
    wq  = [w for cellw in dΩs.quad.cell_weight for w in cellw]      # uniform mesh: |J| constant
    n   = num_free_dofs(Ux)
    Bj  = [zeros(length(xq), n) for _ in 1:Nσ]
    e   = zeros(n)
    for i in 1:n
        e .= 0.0; e[i] = 1.0
        vals = [v for cellv in fh_eval(FEFunction(Ux, e), pts) for v in cellv]
        for j in 1:Nσ
            Bj[j][:, i] .= getindex.(vals, j)
        end
    end
    kgrid = [2π * m / L for m in 0:(pu * ncell)]
    Fk    = [wq[q] * cis(-kk * xq[q]) for kk in kgrid, q in eachindex(xq)]
    return Bj, Fk, kgrid
end
fh_eval(fh, pts) = collect(fh(pts))

"One-line summary of a case (growth in 1/s, wavenumbers as multiples of the carrier k0)."
function se_line(r)
    @sprintf("%-6s Q%d/Q%d aux%-2s A=%.2f  n/λ=%3d dx=%.4f  ndof=%5d | σmax=%+.4e /s  (|ω|=%.3g, k/k0=%.1f, λ/dx=%.2f, uy=%.2f) | ρ=%.3g  auxres=%.1e leak=%.1e",
             r.model, r.pu, r.pe, r.model === :mixed ? string(r.paux) : "-", r.A, r.ncell, r.dx,
             r.ndof, r.σmax, r.ωmax, r.kmax / r.k0, r.kmax > 0 ? (2π / r.kmax) / r.dx : Inf,
             r.uymax, r.rho, r.aux_res, r.mix_leak) * @sprintf("  Mgate=%.1e ncol=%d", r.m_gate, r.ncolour)
end

"σ_max(k): the largest growth rate among eigenpairs whose dominant wavenumber falls in each bin."
function se_sigma_of_k(r; nbins = 12)
    kmax = maximum(r.kdom); edges = range(0, kmax + 1e-9, length = nbins + 1)
    out = Tuple{Float64,Float64}[]
    for b in 1:nbins
        sel = findall(k -> edges[b] <= k < edges[b + 1], r.kdom)
        isempty(sel) && continue
        push!(out, ((edges[b] + edges[b + 1]) / 2, maximum(real.(r.λ[sel]))))
    end
    return out
end
