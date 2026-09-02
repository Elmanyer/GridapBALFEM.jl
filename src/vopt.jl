# ==============================================================================
#  vopt.jl — VERTICAL GRID OPTIMISATION by the Yang & Liu total-relative-error
#  functional (StokesWaveFourierAnalysis.tex §sec: vertical grid optimisation,
#  eq: total error definition).
#
#  E_total = Ē_c + Ē_cg + Ē_shoal + Ē_u + Ē_w
#
#  where each term integrates a relative error over kd ∈ [0, Ω] against the
#  weight W(kd) = exp[(2^{-kd} − 2^{-π}) log 5], and the overbar is
#  normalisation by that term's MEDIAN over the sampled design population.
#
#  WHY A SELF-CONTAINED ANALYTIC BASIS, and not the Gridap assembly.
#  Two of the five terms — E_u and E_w — compare VERTICAL PROFILES pointwise:
#      E_u ~ ∫∫ |u_m(σ) − u_e(σ)| dσ W d(kd)
#  The absolute value inside the σ-integral means these CANNOT be formed from
#  the precomputed Gram tensors (Φ, M, B) the way C, C_g and γ can. They need
#  φ_j(σ) evaluated pointwise, and they need it inside an optimiser's innermost
#  loop. Calling Gridap per candidate mesh would also make the DOF ORDERING a
#  load-bearing assumption — the coefficient vector m must be matched to σ
#  positions, and Gridap's global numbering is not part of its public contract.
#
#  So the basis is built here in closed form, in a known ordering (nodes
#  ascending in σ), and the resulting (Φ, M, B) are CHECKED against
#  `assemble_dispersion_tensors` through the ordering-INVARIANT functional
#  R(μ) = Φᵀ(M+μ|B|)⁻¹Φ — see `vopt_selfcheck`. A permutation cannot hide there.
# ==============================================================================

"""
    SigmaBasis(M, p, c_bdy; nq_elem=…) → basis

Piecewise-Lagrange nodal basis of order `p` on the σ-mesh with element
interfaces `c_bdy`, tabulated at Gauss-Legendre points. Fields:

* `sig`, `wq`   — global quadrature points and weights on [0,1]
* `PHI[q,j]`    — φ_j(σ_q)
* `PIN[q,j]`    — φ_j_int(σ_q) = ∫₀^{σ_q} φ_j
* `nodes[j]`    — σ position of node j (ascending)

`N = M·p + 1`. Gauss-Legendre of `p+1` points per element integrates the degree-
`2p` products `φ_iφ_j` exactly; `φ_int` products reach `2p+2`, so `nq_elem`
defaults to `p+3` points per element, exact for both.
"""
struct SigmaBasis
    M      :: Int
    p      :: Int
    c_bdy  :: Vector{Float64}
    nodes  :: Vector{Float64}
    sig    :: Vector{Float64}
    wq     :: Vector{Float64}
    PHI    :: Matrix{Float64}
    PIN    :: Matrix{Float64}
end

#  Gauss-Legendre nodes/weights on [-1,1] via the Golub-Welsch eigenproblem of
#  the Jacobi matrix — a few lines, exact, and no extra dependency.
function _gauss_legendre(n::Int)
    i  = 1:(n-1)
    b  = @. i / sqrt(4.0*i^2 - 1.0)
    J  = diagm(-1 => collect(b), 1 => collect(b))
    F  = eigen(Symmetric(J))
    x  = F.values
    w  = 2.0 .* (F.vectors[1, :] .^ 2)
    return x, w
end

#  Lagrange basis of order p on p+1 EQUISPACED local nodes t ∈ [0,1], and its
#  antiderivative from the element's left endpoint. Both evaluated at `t`.
function _lagrange_local(p::Int, t::Float64)
    tn = collect(range(0.0, 1.0; length = p+1))
    L  = ones(Float64, p+1)
    for a in 1:(p+1), b in 1:(p+1)
        a == b && continue
        L[a] *= (t - tn[b]) / (tn[a] - tn[b])
    end
    return L
end

function SigmaBasis(M::Int, p::Int, c_bdy::AbstractVector{<:Real}; nq_elem::Int = 0)
    @assert length(c_bdy) == M + 1
    cb = collect(Float64, c_bdy)
    nq = nq_elem > 0 ? nq_elem : p + 3
    xg, wg = _gauss_legendre(nq)
    N = M*p + 1

    nodes = Float64[]
    for k in 1:M
        a, b = cb[k], cb[k+1]
        for i in 0:p
            s = a + (i/p)*(b - a)
            (k > 1 && i == 0) && continue        # shared interface node
            push!(nodes, s)
        end
    end
    @assert length(nodes) == N

    sig = Float64[]; wq = Float64[]
    for k in 1:M
        a, b = cb[k], cb[k+1]; h = b - a
        for q in 1:nq
            push!(sig, a + 0.5*h*(xg[q] + 1.0))
            push!(wq,  0.5*h*wg[q])
        end
    end
    nQ  = length(sig)
    PHI = zeros(Float64, nQ, N)
    PIN = zeros(Float64, nQ, N)

    #  φ: local Lagrange, scattered to global nodes. Global index of local node i
    #  of element k is (k-1)*p + i + 1.
    for k in 1:M
        a, b = cb[k], cb[k+1]; h = b - a
        for q in 1:nq
            row = (k-1)*nq + q
            t   = (sig[row] - a) / h
            L   = _lagrange_local(p, t)
            for i in 0:p
                PHI[row, (k-1)*p + i + 1] = L[i+1]
            end
        end
    end

    #  φ_int(σ) = ∫₀^σ φ. Accumulate element by element: the completed integral
    #  of each preceding element plus a within-element Gauss-Legendre partial
    #  integral from the element's left endpoint to σ_q. The partial integral is
    #  itself computed by a sub-quadrature, exact because φ is a polynomial.
    xs, ws = _gauss_legendre(nq + 2)
    carry  = zeros(Float64, N)
    for k in 1:M
        a, b = cb[k], cb[k+1]; h = b - a
        for q in 1:nq
            row = (k-1)*nq + q
            σq  = sig[row]; hq = σq - a
            acc = copy(carry)
            for r in eachindex(xs)
                tt = a + 0.5*hq*(xs[r] + 1.0)
                ww = 0.5*hq*ws[r]
                L  = _lagrange_local(p, (tt - a)/h)
                for i in 0:p
                    acc[(k-1)*p + i + 1] += ww * L[i+1]
                end
            end
            PIN[row, :] = acc
        end
        #  advance the carry by the FULL integral over element k
        for r in eachindex(xs)
            tt = a + 0.5*h*(xs[r] + 1.0)
            ww = 0.5*h*ws[r]
            L  = _lagrange_local(p, (tt - a)/h)
            for i in 0:p
                carry[(k-1)*p + i + 1] += ww * L[i+1]
            end
        end
    end
    return SigmaBasis(M, p, cb, nodes, sig, wq, PHI, PIN)
end

"""
    vopt_tensors(basis) → (Phi, Mmat, B)

The three dispersion tensors in this file's own node ordering.
`B = −∫φ_int_i φ_int_j` — the SAME stored-negative convention as `vertical.jl`.
"""
function vopt_tensors(b::SigmaBasis)
    W    = Diagonal(b.wq)
    Phi  = vec(sum(W * b.PHI; dims = 1))
    Mmat = b.PHI' * W * b.PHI
    B    = -(b.PIN' * W * b.PIN)
    return (Phi = Phi, Mmat = Mmat, B = B)
end

# ------------------------------------------------------------------------------
#  The five relative-error integrands
# ------------------------------------------------------------------------------

"""
    vopt_weight(kd) → W

Yang & Liu's weighting `W = exp[(2^{-kd} − 2^{-π}) log 5]`.

⚠ It does NOT decay to zero: `W(0)=4.17`, `W(∞)→5^{-2^{-π}}=0.834`. So it
re-weights shallow water by a factor ≈5 but does not make the kd-integral
convergent on its own — the upper limit `Ω` is a genuine DESIGN BAND, not a
numerical cutoff. That is why `vopt_Omega` exists.
"""
vopt_weight(kd::Float64) = exp((2.0^(-kd) - 2.0^(-pi)) * log(5.0))

"""
    polarization(T, kd) → m

The first-order polarization `m = (M + μ|B|)⁻¹Φ`, `μ = (kd)²`, SCALED so that
its depth average `Φᵀm` equals the exact profile's, `1/kd`.

That scaling is the physical one: continuity gives `−ωa + kd·Φᵀm = 0` for the
model and `∫₀¹ u_e dσ = aω/kd` for Airy, so matching the depth-integrated flux
compares the two SHAPES at equal volume transport. It also makes `w_m(1) = 1`
identically, matching `w_e(1) = 1` — see `profile_errors`.
"""
function polarization(T, kd::Float64)
    mu   = kd^2
    absB = -T.B
    m    = (T.Mmat .+ mu .* absB) \ T.Phi
    mean = dot(T.Phi, m)
    (isfinite(mean) && abs(mean) > 1e-300) || return fill(NaN, length(m))
    return m .* ((1.0/kd) / mean)
end

"""
    profile_errors(b, T, kd) → (Eu, Ew)

`∫₀¹|u_m−u_e|dσ / u_e(1)` and `∫₀¹|w_m−w_e|dσ / w_e(1)` at one `kd`, with

    u_e(σ) = cosh(kd σ)/sinh(kd),   u_e(1) = coth(kd)
    w_e(σ) = sinh(kd σ)/sinh(kd),   w_e(1) = 1
    u_m(σ) = Σ_j m_j φ_j(σ),        w_m(σ) = kd Σ_j m_j φ_int_j(σ)

`w_m` follows from `u_m` by continuity, so the single scaling in
`polarization` fixes both. Evaluated on the basis' own Gauss-Legendre rule; the
integrand is |piecewise-polynomial − analytic| so the rule is not exact, but the
kink locations are element interfaces and the sub-element error is far below the
design differences being compared.
"""
function profile_errors(b::SigmaBasis, T, kd::Float64)
    m = polarization(T, kd)
    any(!isfinite, m) && return (Inf, Inf)
    um = b.PHI * m
    wm = kd .* (b.PIN * m)
    s  = b.sig
    shk = sinh(kd)
    (isfinite(shk) && shk > 0) || return (Inf, Inf)
    ue = cosh.(kd .* s) ./ shk
    we = sinh.(kd .* s) ./ shk
    Eu = dot(b.wq, abs.(um .- ue)) / (cosh(kd)/shk)
    Ew = dot(b.wq, abs.(wm .- we))
    return (Eu, Ew)
end

"""
    raw_errors(M, p, c_bdy; Omega, nkd=…) → NamedTuple

The five RAW (un-normalised) error integrals of `eq: total error definition`
over `kd ∈ (0, Ω]`. Composite Gauss-Legendre in `kd`; the integrands are smooth
where the model is usable and merely large where it is not, so a moderate rule
suffices and the cost is dominated by the `N×N` solve per node.

`kd → 0` is a removable singularity in every term (all relative errors vanish),
but `sinh(kd)` and the `1/kd` scaling are numerically poor there, so the rule
starts at `kd_min = 1e-3` rather than 0.
"""
function raw_errors(M::Int, p::Int, c_bdy::AbstractVector{<:Real};
                    Omega::Float64, nkd::Int = 48, kd_min::Float64 = 1e-3)
    b = SigmaBasis(M, p, c_bdy)
    T = vopt_tensors(b)
    xg, wg = _gauss_legendre(nkd)
    a, c   = kd_min, Omega
    kds    = @. a + 0.5*(c - a)*(xg + 1.0)
    wk     = @. 0.5*(c - a)*wg

    Ec = Ecg = Esh = Eu = Ew = 0.0
    for (kd, w) in zip(kds, wk)
        W  = vopt_weight(kd)
        mu = kd^2
        Rm = model_R(T, mu)
        Re = airy_R(mu)
        (isfinite(Rm[1]) && Rm[1] > 0) || return (Ec=Inf,Ecg=Inf,Esh=Inf,Eu=Inf,Ew=Inf)
        #  C² = g d R  ⇒  (C_m²−C_e²)/C_e² = (R_m−R_e)/R_e
        Ec  += abs((Rm[1] - Re[1]) / Re[1]) * W * w
        #  C_g² ∝ (R+μR')²/R
        gm = (Rm[1] + mu*Rm[2])^2 / Rm[1]
        ge = (Re[1] + mu*Re[2])^2 / Re[1]
        Ecg += abs((gm - ge) / ge) * W * w
        #  shoaling gradient, from the SAME routine for model and exact so the
        #  finite-difference truncation is common and cancels (see wave_properties)
        wp   = wave_properties(T, kd)
        Esh += abs((wp.gamma_e - wp.gamma) / kd) * W * w
        pu, pw = profile_errors(b, T, kd)
        Eu  += pu * W * w
        Ew  += pw * W * w
    end
    #  the shoaling term is exponentiated, per eq: total error definition
    return (Ec = Ec, Ecg = Ecg, Esh = exp(Esh) - 1.0, Eu = Eu, Ew = Ew)
end

# ------------------------------------------------------------------------------
#  Median normalisation, and the optimiser
# ------------------------------------------------------------------------------

#  Unconstrained parametrisation of an ORDERED interface set. z ∈ R^M ↦ element
#  widths by softmax ↦ interfaces by cumulative sum. Ordering and the c_0=0,
#  c_M=1 endpoints hold by construction, so the optimiser never has to reject a
#  candidate — which a box-constrained search on c directly would have to do.
function _z_to_cbdy(z::AbstractVector{<:Real})
    e = exp.(z .- maximum(z))
    h = e ./ sum(e)
    return vcat(0.0, cumsum(h))
end

_median(v) = (s = sort(v); n = length(s);
              isodd(n) ? s[(n+1)÷2] : 0.5*(s[n÷2] + s[n÷2+1]))

"""
    vopt_medians(M, p; Omega, nsample) → (med, samples)

Medians of the five raw errors over a sampled design population — the `overbar`
of `eq: total error definition`. The medians are per-term SCALES: they are what
makes a 2 % celerity error and a 0.05 profile error commensurable, so the
argmin genuinely depends on them and they must be computed over a population
representative of the search, not over one design.

Sampling is a deterministic low-discrepancy sweep of the width simplex, so the
result is reproducible without seeding an RNG.
"""
function vopt_medians(M::Int, p::Int; Omega::Float64, nsample::Int = 64,
                      nkd::Int = 32)
    zs = Vector{Vector{Float64}}()
    if M == 1
        push!(zs, [0.0])
    else
        #  additive-recurrence (Kronecker) sequence in M−1 dims: low-discrepancy,
        #  deterministic, no RNG state
        g = 1.0 / (2.0^(1.0/M))
        for i in 1:nsample
            z = Float64[0.0]
            for k in 1:(M-1)
                frac = mod(i * g^k, 1.0)
                push!(z, 6.0*(frac - 0.5))       # z ∈ [−3,3] ⇒ wide width ratios
            end
            push!(zs, z)
        end
    end
    rows = NamedTuple[]
    for z in zs
        r = raw_errors(M, p, _z_to_cbdy(z); Omega = Omega, nkd = nkd)
        all(isfinite, (r.Ec, r.Ecg, r.Esh, r.Eu, r.Ew)) && push!(rows, r)
    end
    isempty(rows) && error("vopt_medians: no finite sample at M=$M p=$p Ω=$Omega")
    med = (Ec  = _median([r.Ec  for r in rows]),
           Ecg = _median([r.Ecg for r in rows]),
           Esh = _median([r.Esh for r in rows]),
           Eu  = _median([r.Eu  for r in rows]),
           Ew  = _median([r.Ew  for r in rows]))
    return (med = med, nsample = length(rows))
end

"""
    total_error(M, p, c_bdy; Omega, med, nkd) → E_total

`E_total = Ē_c + Ē_cg + Ē_shoal + Ē_u + Ē_w`, each term divided by its median.
"""
function total_error(M::Int, p::Int, c_bdy::AbstractVector{<:Real};
                     Omega::Float64, med, nkd::Int = 32)
    r = raw_errors(M, p, c_bdy; Omega = Omega, nkd = nkd)
    any(!isfinite, (r.Ec, r.Ecg, r.Esh, r.Eu, r.Ew)) && return Inf
    return r.Ec/med.Ec + r.Ecg/med.Ecg + r.Esh/med.Esh + r.Eu/med.Eu + r.Ew/med.Ew
end

#  Compass (pattern) search on z. Low-dimensional (M−1 ≤ 3 here) and the
#  objective is continuous but not smooth (absolute values inside the
#  integrals), so a derivative-free method is the right class; a pattern search
#  is chosen over Nelder-Mead because it cannot collapse its simplex onto a
#  kink and stall there.
function _pattern_search(f, z0::Vector{Float64}; step0 = 0.8, tol = 1e-6,
                         maxit = 4000)
    z = copy(z0); fz = f(z); step = step0; it = 0
    n = length(z)
    while step > tol && it < maxit
        improved = false
        for k in 2:n, s in (+1.0, -1.0)     # z[1] pinned: softmax is shift-invariant
            zt = copy(z); zt[k] += s*step
            ft = f(zt); it += 1
            if ft < fz
                z, fz = zt, ft; improved = true; break
            end
        end
        improved || (step *= 0.5)
    end
    return z, fz
end

"""
    optimise_cbdy(M, p; Omega, nsample, nkd, nrestart) → NamedTuple

Minimise `E_total` over the element interfaces at fixed `(M,p)` and design band
`Omega`. Returns `(c_bdy, E, med, raw)`.

Two stages, and the first is not optional: the sampling pass supplies the
MEDIANS, and the same samples then seed the local search. Multi-start from the
best `nrestart` samples guards against the objective's local minima, which are
real — the profile terms have kinks at element interfaces.
"""
function optimise_cbdy(M::Int, p::Int; Omega::Float64, nsample::Int = 64,
                       nkd::Int = 32, nrestart::Int = 4, verbose::Bool = false)
    M == 1 && return (c_bdy = [0.0, 1.0], E = NaN, med = nothing, raw = nothing)
    mm  = vopt_medians(M, p; Omega = Omega, nsample = nsample, nkd = nkd)
    med = mm.med
    obj(z) = total_error(M, p, _z_to_cbdy(z); Omega = Omega, med = med, nkd = nkd)

    g = 1.0 / (2.0^(1.0/M))
    cands = Vector{Tuple{Float64,Vector{Float64}}}()
    for i in 1:nsample
        z = Float64[0.0]
        for k in 1:(M-1); push!(z, 6.0*(mod(i*g^k, 1.0) - 0.5)); end
        push!(cands, (obj(z), z))
    end
    sort!(cands; by = first)

    best = (Inf, cands[1][2])
    for (_, z0) in cands[1:min(nrestart, length(cands))]
        z, fz = _pattern_search(obj, z0)
        fz < best[1] && (best = (fz, z))
    end
    cb = _z_to_cbdy(best[2])
    verbose && @info "optimise_cbdy M=$M p=$p Ω=$Omega → $(round.(cb; digits=4))  E=$(best[1])"
    return (c_bdy = cb, E = best[1], med = med,
            raw = raw_errors(M, p, cb; Omega = Omega, nkd = nkd))
end

"""
    vopt_selfcheck(M, p, c_bdy; mus) → max relative difference in R(μ)

Validates this file's analytic basis against `assemble_dispersion_tensors` (the
Gridap assembly) through `R(μ) = Φᵀ(M+μ|B|)⁻¹Φ`, which is INVARIANT under any
permutation of the DOFs. A node-ordering mismatch — the one failure mode a
hand-built basis is prone to — cannot hide in it.
"""
function vopt_selfcheck(M::Int, p::Int, c_bdy::AbstractVector{<:Real};
                        mus = [0.1, 1.0, 9.0, 100.0, 400.0])
    T  = vopt_tensors(SigmaBasis(M, p, c_bdy))
    Tg = assemble_dispersion_tensors(M, p, collect(Float64, c_bdy))
    return maximum(abs(model_R(T, mu)[1]/model_R(Tg, mu)[1] - 1.0) for mu in mus)
end

# ------------------------------------------------------------------------------
#  The design band Ω — a CALIBRATED FIXED POINT, not a free knob
# ------------------------------------------------------------------------------

"""
    VOPT_KAPPA

Surface-resolution constant `κ = Ω · h_top`, where `h_top = Δσ_top/p` is the
NODE spacing in the top element of the optimised mesh.

**Why this is the criterion.** `W` does not decay to zero (`vopt_weight`), so
`Ω` is a genuine design band and the optimum migrates with it. A single `Ω`
cannot serve every basis: calibrating on Yang & Liu Table 1 gives `Ω* = 7.60,
29.0, 100.0` for `M = 2,3,4` — the richer basis earns a wider band. What IS
invariant across those three is the product with the surface node spacing:

    M=2: 7.60 × 0.2716 = 2.06     M=3: 29.0 × 0.0771 = 2.24
    M=4: 100.0 × 0.0220 = 2.20                        mean 2.17 (±4 %)

That is the variational statement of `subsec: design rule` made into a design
rule: by (P2) the space must approximate `u_⋆ = cosh(kdσ)/cosh(kd)`, an
exponential boundary layer of thickness `1/(kd)` pinned to `σ=1`. Requiring the
surface node spacing to resolve that layer AT THE DESIGN BAND closes the loop —
`Ω` is then whatever band the mesh it produces can actually resolve.

Dividing by `p` rather than using `Δσ_top` directly is the basis-order-aware
form: a `p=2` top element carries an interior node, so it resolves a layer twice
as thick as a `p=1` element of the same size. ⚠ This generalisation is
REASONED, not calibrated — there is no published `p≥2` optimum to check it
against. It is stated here so a later measurement can refute it.
"""
const VOPT_KAPPA = 2.17

"""
    vopt_Omega(M, p; kappa=VOPT_KAPPA, …) → Ω

The design band as the fixed point of `Ω · Δσ_top(Ω)/p = κ`, where `Δσ_top(Ω)`
is the top-element thickness of the mesh `optimise_cbdy` returns at band `Ω`.

Solved by bisection on `log Ω`: `Δσ_top` decreases with `Ω` (a wider band drives
the interface towards the surface), so the residual `Ω·Δσ_top/p − κ` is
increasing and the bracket is reliable. `M=1` has no free interface and returns
its single-element band directly.
"""
function vopt_Omega(M::Int, p::Int; kappa::Float64 = VOPT_KAPPA,
                    lo::Float64 = 1.5, hi::Float64 = 400.0, tol::Float64 = 0.02,
                    nsample::Int = 48, nkd::Int = 24, maxit::Int = 24)
    M == 1 && return kappa * p / 1.0          # Δσ_top ≡ 1
    resid(Om) = (r = optimise_cbdy(M, p; Omega = Om, nsample = nsample,
                                   nkd = nkd, nrestart = 2);
                 Om * (1.0 - r.c_bdy[end-1]) / p - kappa)
    a, b = lo, hi
    fa, fb = resid(a), resid(b)
    fa > 0 && return a
    fb < 0 && return b
    for _ in 1:maxit
        m  = sqrt(a*b)                        # bisect in log Ω: the band spans 2 decades
        fm = resid(m)
        #  NOTE: written as an if/else, NOT `fm < 0 ? (a, fa = m, fm) : …`.
        #  Julia parses `(a, fa = m, fm)` as a TUPLE with a named field, not as
        #  the tuple assignment it looks like, so the bracket never moved and the
        #  bisection silently returned its first midpoint for every input.
        if fm < 0
            a = m; fa = fm
        else
            b = m; fb = fm
        end
        (b/a - 1.0) < tol && break
    end
    return sqrt(a*b)
end

"""
    optimised_cbdy(M, p; kappa=VOPT_KAPPA, …) → NamedTuple

The full design pipeline: resolve the band `Ω` by `vopt_Omega`, then minimise
`E_total` at that band. Returns `(c_bdy, Omega, E, raw)`.
"""
function optimised_cbdy(M::Int, p::Int; kappa::Float64 = VOPT_KAPPA,
                        nsample::Int = 80, nkd::Int = 32, nrestart::Int = 4)
    M == 1 && return (c_bdy = [0.0, 1.0], Omega = NaN, E = NaN, raw = nothing)
    Om = vopt_Omega(M, p; kappa = kappa)
    r  = optimise_cbdy(M, p; Omega = Om, nsample = nsample, nkd = nkd,
                       nrestart = nrestart)
    return (c_bdy = r.c_bdy, Omega = Om, E = r.E, raw = r.raw)
end
