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
#  ⚠ THE FIVE TERMS ARE NOT OF A COMMON FORM — Yang & Liu (2024) eq. (3.10):
#      E_c     = ∫ |[C_m² − C_e²] W / C_e²| d(kd)        SQUARED celerity, |·|
#      E_cg    = ∫ |[C_gm − C_ge] W / C_ge| d(kd)        FIRST power,     |·|
#      E_shoal = exp[ ∫ (γ_e − γ_m) W / kd  d(kd) ] − 1  SIGNED, no |·|
#      E_u     = ∫ [∫₀¹ |u_m−u_e| dσ / u_e(1)] W d(kd)   |·| inside dσ
#      E_w     = ∫ [∫₀¹ |w_m−w_e| dσ / w_e(1)] W d(kd)   |·| inside dσ
#  Squaring E_cg, or wrapping E_shoal in abs, are BOTH departures from the
#  reference. They were present here until 2026-09-08 and are corrected; the
#  measured effect on the optimum is small (≤0.003 in c at M=2) but the forms
#  are the specification, not an approximation to it.
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
numerical cutoff. Ω must therefore be CHOSEN and stated, not derived here.
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

⚠ `Esh` is returned as `exp(I) − 1` with `I` the SIGNED integral, so it is
negative in practice — unlike the other four, which are absolute-valued per
(3.10). That is intended and needs no repair: `exp(I) − 1` is SINGLE-SIGNED
across the design population, so the median `vopt_medians` divides by carries
the same sign and cancels it, making the normalised term positive and exactly
equal to the `|exp(I) − 1|` form. The cancellation requires the median to be
taken over the `exp(I) − 1` population — see the warning in `vopt_medians`.
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
        #  C² = g d R  ⇒  (C_m²−C_e²)/C_e² = (R_m−R_e)/R_e.  E_c is the one term
        #  (3.10) states in the SQUARED celerity, so this form is the paper's.
        Ec  += abs((Rm[1] - Re[1]) / Re[1]) * W * w
        #  C_g, and γ below, come from the SAME routine for model and exact, so
        #  the finite-difference truncation is common and cancels.
        wp   = wave_properties(T, kd)
        #  ⚠ (3.10) takes C_g to the FIRST power — NOT squared like E_c. The
        #  asymmetry is the paper's own; a C_g² form is ≈2× this to leading
        #  order and does not survive the median normalisation unchanged.
        Ecg += abs((wp.Cg - wp.Cge) / wp.Cge) * W * w
        #  ⚠ (3.10) carries NO |·| on the shoaling term — verified at 400 dpi
        #  against the PDF (the next line shows |u_m−u_e| with unmistakable
        #  bars, so the omission is deliberate, not a rendering artefact). It is
        #  SIGNED, and that is correct AND SELF-CONSISTENT, for a reason worth
        #  stating because it is easy to "fix" wrongly:
        #
        #    exp(I)−1 is SINGLE-SIGNED across the design population (γ_m > γ_e
        #    throughout), so its median is negative too — measured on the LFE-2
        #    Ω=8 scan: median(exp(I)−1) = −0.19158, median(|exp(I)−1|) = +0.19158,
        #    equal in magnitude. The normalised term E_sh/med is therefore
        #    POSITIVE, and identical to the |·|-outside form: both give c₂=0.7310.
        #    THE SIGN CANCELS AGAINST THE SAME-SIGNED MEDIAN.
        #
        #  ⚠ This holds ONLY if the median is taken over the exp(I)−1 population.
        #  Normalising by median(exp(I)) ≈ +0.81 breaks the cancellation, inverts
        #  the term's role, and drives the optimum to the scan edge (c₂→0). That
        #  was a live defect here on 2026-09-08; see `vopt_medians`.
        Esh += ((wp.gamma_e - wp.gamma) / kd) * W * w
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

#  DEFAULT SCAN INCREMENT. Yang & Liu scan "c2 values ranging from 0 to 1 with
#  an increment of 0.001" — exact for M=2 (999 designs). The scan is a product
#  over M−1 ordered interfaces, so a fixed increment costs C(1/dc−1, M−1) and
#  0.001 is unaffordable past M=2; these keep every M near ~5·10³ designs.
_default_dc(M::Int) = M <= 2 ? 0.001 : (M == 3 ? 0.01 : 0.03)

"""
    _scan_designs(M, dc) → Vector{Vector{Float64}}

Every design on the uniform SCAN of the ordered interface simplex
`0 < c_1 < … < c_{M-1} < 1` at increment `dc`, each returned as a full
`c_bdy = [0, c_1, …, c_{M-1}, 1]`. This is the population `vopt_medians`
normalises over and the set `optimise_cbdy` searches — the same set, because
that is what the reference does.
"""
function _scan_designs(M::Int, dc::Float64)
    M == 1 && return [[0.0, 1.0]]
    n = round(Int, 1.0/dc) - 1
    n >= M-1 || error("_scan_designs: dc=$dc too coarse for M=$M")
    g = [i*dc for i in 1:n]
    out = Vector{Vector{Float64}}()
    idx = collect(1:(M-1))                      # strictly increasing indices
    while true
        push!(out, vcat(0.0, g[idx], 1.0))
        k = M-1
        while k >= 1 && idx[k] == n - (M-1-k); k -= 1; end
        k == 0 && break
        idx[k] += 1
        for j in (k+1):(M-1); idx[j] = idx[j-1] + 1; end
    end
    return out
end

#  Inverse of `_z_to_cbdy`, so a scan point can seed the local refinement.
#  The softmax is shift-invariant, so z is pinned by z[1] = 0.
function _cbdy_to_z(c::AbstractVector{<:Real})
    h = diff(collect(Float64, c))
    all(>(0.0), h) || error("_cbdy_to_z: non-increasing c_bdy $c")
    z = log.(h); return z .- z[1]
end

"""
    vopt_medians(M, p; Omega, dc, nkd) → (med, designs, raws, nsample)

Medians of the five raw errors over the DESIGN SCAN — the `overbar` of
`eq: total error definition`. The medians are per-term SCALES: they are what
makes a 2 % celerity error and a 0.05 profile error commensurable, so the
argmin genuinely depends on them and they cannot be computed from one design.

⚠ THE POPULATION IS THE SCAN ITSELF, not an auxiliary sample. Yang & Liu take
the median over the sweep they optimise on ("c2 values ranging from 0 to 1 with
an increment of 0.001"), so the normalisation and the search see the same set.
Until 2026-09-08 this used a Kronecker low-discrepancy sweep of the WIDTH
simplex over z ∈ [−3,3] — a different population, and since the medians are the
trade-off weights, a different argmin. Measured against the published table 1
(LFE-2, Ω=8): the scan gives c₂ = 0.731 against 0.728, the Kronecker sweep
0.738. The gap grows with M, because the two populations diverge in shape as
the dimension rises.

The scan is deterministic — no RNG, and reproducible from `(M, dc)` alone.
"""
function vopt_medians(M::Int, p::Int; Omega::Float64, dc::Float64 = 0.0,
                      nkd::Int = 32)
    designs = _scan_designs(M, dc > 0 ? dc : _default_dc(M))
    rows = NamedTuple[]; keep = Vector{Vector{Float64}}()
    for c in designs
        r = raw_errors(M, p, c; Omega = Omega, nkd = nkd)
        if all(isfinite, (r.Ec, r.Ecg, r.Esh, r.Eu, r.Ew))
            push!(rows, r); push!(keep, c)
        end
    end
    isempty(rows) && error("vopt_medians: no finite design on the scan at M=$M p=$p Ω=$Omega")
    #  ⚠ THE SHOALING MEDIAN IS TAKEN OVER THE ERROR POPULATION exp(I)−1, like
    #  every other term — NOT over exp(I). The distinction is not cosmetic: I is
    #  small, so median(exp(I)) is POSITIVE (0.80842 on the LFE-2 Ω=8 scan) while
    #  the numerator exp(I)−1 is NEGATIVE (median −0.19158). Dividing by it flips
    #  the term's sign and inverts its role in the minimisation — the optimum
    #  then runs to the scan edge, c₂→0. Taking the median over the SAME
    #  population as the numerator makes the shared sign cancel, leaving a
    #  positive term of O(1) on the median design: commensurable with the other
    #  four, which is the whole purpose of the overbar. Measured ratio 4.2×.
    med = (Ec  = _median([r.Ec  for r in rows]),
           Ecg = _median([r.Ecg for r in rows]),
           Esh = _median([r.Esh for r in rows]),
           Eu  = _median([r.Eu  for r in rows]),
           Ew  = _median([r.Ew  for r in rows]))
    return (med = med, designs = keep, raws = rows, nsample = length(rows))
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
    optimise_cbdy(M, p; Omega, dc, nkd, nrestart, refine) → NamedTuple

Minimise `E_total` over the element interfaces at fixed `(M,p)` and design band
`Omega`. Returns `(c_bdy, E, med, raw)`.

ONE scan does both jobs, which is the reference's own procedure: it supplies the
MEDIANS and it is the set the argmin is taken over. The `raw_errors` evaluations
are therefore shared — following the reference here costs nothing over the
previous auxiliary-sample scheme.

`refine` then runs a local pattern search from the best `nrestart` scan points,
which is NOT in the reference and is needed only because `_default_dc` coarsens
with M (0.03 at M=4 against the published 3-decimal node sets). Multi-start
guards against the objective's local minima, which are real — the profile terms
have kinks at element interfaces. Set `refine=false` for the literal procedure.
"""
function optimise_cbdy(M::Int, p::Int; Omega::Float64, dc::Float64 = 0.0,
                       nkd::Int = 32, nrestart::Int = 4, refine::Bool = true,
                       verbose::Bool = false)
    M == 1 && return (c_bdy = [0.0, 1.0], E = NaN, med = nothing, raw = nothing)
    mm  = vopt_medians(M, p; Omega = Omega, dc = dc, nkd = nkd)
    med = mm.med
    Et  = [r.Ec/med.Ec + r.Ecg/med.Ecg + r.Esh/med.Esh + r.Eu/med.Eu + r.Ew/med.Ew
           for r in mm.raws]
    ord  = sortperm(Et)
    best = (Et[ord[1]], _cbdy_to_z(mm.designs[ord[1]]))
    if refine
        obj(z) = total_error(M, p, _z_to_cbdy(z); Omega = Omega, med = med, nkd = nkd)
        for i in ord[1:min(nrestart, length(ord))]
            z, fz = _pattern_search(obj, _cbdy_to_z(mm.designs[i]); step0 = 0.25)
            fz < best[1] && (best = (fz, z))
        end
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
#  THE DESIGN BAND Ω IS AN INPUT, NOT SOMETHING THIS FILE DERIVES
# ------------------------------------------------------------------------------
#  ⚠ REMOVED 2026-09-08: `VOPT_KAPPA`, `vopt_Omega` and `optimised_cbdy`.
#
#  They closed Ω for p ≥ 2 (where no published band exists) through the fixed
#  point Ω·Δσ_top/p = κ ≈ 2.17. That constant was a FIT, not a result: Ω was
#  first tuned to reproduce Yang & Liu table 1 (Ω* = 7.60, 29.0, 100.0), and the
#  near-constancy of the product was noticed afterwards and promoted to a design
#  rule. Three things make it unusable:
#
#    1. Those Ω* were fitted against the PRE-2026-09-08 objective, so κ is
#       downstream of the very code that has since been corrected.
#    2. The same boundary-layer argument, carried through analytically in
#       `StokesWaveFourierAnalysis.tex`, gives Δσ_top ≈ 2.94/kd_max — a 35 %
#       disagreement with the fitted 2.11-2.17. A physical constant should not
#       differ from its own derivation by a third.
#    3. The p-aware form (dividing by p) was, in its own words, "REASONED, not
#       calibrated", with no published p ≥ 2 optimum to test it against.
#
#  ⚠ DO NOT REINTRODUCE IT. To design a mesh, call `optimise_cbdy(M, p; Omega)`
#  with a band you can justify: it scans the parametric space, normalises each
#  relative error by the median over that same sweep, and returns the minimum.
#  It is κ-free. For p = 1 the bands are published (Ω = 8, 24, 80 for M = 2,3,4);
#  for p ≥ 2 choosing Ω is an OPEN QUESTION and must not be closed by a fit.
# ------------------------------------------------------------------------------
