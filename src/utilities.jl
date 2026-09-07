# ==============================================================
#  utilities.jl — physical setup helpers + the sequential driver
#
#  Everything needed to turn physical inputs into a running simulation on one
#  process: the dispersion-relation solver (`find_wavenumber`), the quadratic
#  sponge-layer profile (`make_sponge`), the internal Gaussian wavemakers
#  (`make_wavemaker_line/point`), and the top-level driver `setup_and_run`,
#  which wires the two stages together and marches the time loop.
# ==============================================================

"""
    find_wavenumber(omega, d, g)

Newton iteration for the linear (Airy) dispersion relation ω² = g k tanh(kd).
"""
function find_wavenumber(omega::Float64, d::Float64, g::Float64)
    k = omega^2 / g
    for _ in 1:100
        f  = g * k * tanh(k*d) - omega^2
        df = g * (tanh(k*d) + k*d*(1.0 - tanh(k*d)^2))
        dk = f / df
        k -= dk
        abs(dk) < 1e-13 * abs(k) && break
    end
    return k
end

"""
    dispersion_ratio(vert, g, d, kd_vals)

- vert:    the vertical tensor bundle from `assemble_vertical_tensors`.
- g:       gravitational acceleration [m/s²].
- d:       still-water depth [m].
- kd_vals: vector of kd values to evaluate.

Model-to-exact phase speed ratio Cm/Ce per kd:
  Ce = √(g tanh(kd)/k),   Cm² = g·d · Φᵀ (Mmat − B·kd²)⁻¹ Φ
(`B ≤ 0`, so Mmat − B·kd² = Mmat + |B|·kd² > 0.)
"""
function dispersion_ratio(vert, g::Float64, d::Float64,
                              kd_vals::AbstractVector{Float64})
    A, B, Phi = vert.Mmat, vert.B, vert.Phi
    ratios = zeros(length(kd_vals))
    for (idx, kd) in enumerate(kd_vals)
        k     = kd / d
        Ce    = sqrt(g * tanh(kd) / k)
        M_eff = A .- B .* kd^2
        try
            Cm_sq = g * d * dot(Phi, M_eff \ Phi)
            ratios[idx] = Cm_sq > 0 ? sqrt(Cm_sq) / Ce : NaN
        catch
            ratios[idx] = NaN
        end
    end
    return ratios
end

"""
    applicable_kd(vert, g, d; err=0.02, kd_max=200.0, n=2000)

Maximum kd with |Cm/Ce − 1| ≤ err.
"""
function applicable_kd(vert, g::Float64, d::Float64;
                           err::Float64=0.02, kd_max::Float64=200.0, n::Int=2000)
    kd_vals = LinRange(0.01, kd_max, n)
    ratios  = dispersion_ratio(vert, g, d, kd_vals)
    idx     = findlast(abs.(ratios .- 1.0) .<= err)
    return isnothing(idx) ? NaN : Float64(kd_vals[idx])
end

# ===========================================================================
#  LINEAR WAVE PROPERTIES — C, C_g, γ, for the model and for Airy
#
#  StokesWaveFourierAnalysis.tex §subsec: stokes numbers. These are what the
#  APPLICABLE RANGE is actually defined on, and the σ-mesh design problem is
#  posed over all three:
#
#      R(μ)  = Φᵀ(M + μ|B|)⁻¹Φ ,           μ = (kd)² ,  |B| = −B_stored
#      C     = √(g d R)
#      C_g   = C (1 + μ R'/R)              [eq: LFEM group velocity]
#      γ     = ∂lnA/∂ln d |_ω = −½ ∂ln C_g/∂ln d |_ω   [eq: shoaling gradient]
#
#  ⚠ **OPTIMISING C ALONE IS A TRAP, AND THE DERIVATION SAYS SO.** Driving the
#  interface towards the free surface raises kd_app^(C) monotonically while
#  kd_app^(Cg) COLLAPSES — for the two-element mesh, C rises to 15.5 at c₁=0.81
#  while C_g falls off a cliff from 10.5 to 3.3 between c₁=0.81 and 0.82. A design
#  excellent in celerity can be useless for energy propagation.
#
#  ⚠ **γ CHANGES SIGN near kd ≈ 1.2** (where C_g peaks), so a RELATIVE error
#  measure on it is singular and must not be used. The tolerance on γ is ABSOLUTE
#  (0.02); on C and C_g it is relative (2 %). This is a genuine methodological
#  subtlety and a known source of discrepancy between published tables.
# ===========================================================================

"""
    model_R(vert, mu) → (R, dR, d2R)

The dispersion functional `R(μ) = Φᵀ(M+μ|B|)⁻¹Φ` and its first two μ-derivatives,
`R' = −ΦᵀS|B|SΦ` and `R'' = 2ΦᵀS|B|S|B|SΦ` with `S = (M+μ|B|)⁻¹`. One factorisation
serves all three.
"""
function model_R(vert, mu::Float64)
    absB = -vert.B                       # B_stored = −B̃ ⇒ |B| = −B_stored
    S    = (vert.Mmat .+ mu .* absB) \ Matrix{Float64}(I, size(absB)...)
    SPhi = S * vert.Phi
    BS   = absB * SPhi
    R    = dot(vert.Phi, SPhi)
    dR   = -dot(SPhi, BS)
    d2R  = 2.0 * dot(BS, S * BS)
    return (R, dR, d2R)
end

#  Airy: R_e(μ) = tanh(x)/x with x = √μ. Its (1 + μR'/R) reduces EXACTLY to
#  ½(1 + 2x/sinh 2x), i.e. the textbook group-velocity factor — asserted in
#  test_dispersion_curve rather than assumed here.
function airy_R(mu::Float64)
    x = sqrt(mu)
    t = tanh(x); c2 = sech(x)^2
    R  = t/x
    dR = c2/(2x^2) - t/(2x^3)
    return (R, dR)
end

"""
    wave_properties(vert, kd; g=g, d=1.0) → (C, Cg, gamma, Ce, Cge, gamma_e)

The three linear wave properties for the model and for Airy at one `kd`.

`γ` is formed by a CENTRED DIFFERENCE in `ln d` at fixed `ω`, using the SAME
routine for both model and exact — so the finite-difference truncation is common
to the two and cancels in the difference `|γ_m − γ_e|` that the tolerance is
applied to. Doing it analytically would need `R''` chained through an implicit
`x(d)`; that derivation is exactly the kind that goes wrong silently, and the
quantity is validated against a published table either way.

`γ` depends only on `kd` (the `d`-dependence cancels), so `d` is a scale here.
"""
function wave_properties(vert, kd::Float64; g::Float64=g, d::Float64=1.0)
    #  ω² = g x² R(x²)/d, so at fixed ω a change in d moves x. Invert F(x)=x²R(x²)
    #  by NEWTON, using the derivative we already have:
    #      F'(x) = 2x(R + x²R') = 2xR·G,      G = 1 + μR'/R
    #  The only call site perturbs d by exp(±1e-4), so the guess x = kd is within
    #  1e-4 of the root and Newton lands in one or two steps. (Bisection to full
    #  precision costs ~200 evaluations here, and this function sits inside a
    #  minimax optimiser's innermost loop — the difference is hours.)
    #  F is increasing where the model is usable and SATURATES beyond its
    #  bandwidth (F → Φᵀ|B|⁻¹Φ), so a non-positive derivative means "no root in
    #  this direction" and the iteration stops rather than diverging.
    FdF_m(x) = (R = model_R(vert, x^2); (x^2*R[1], 2x*(R[1] + x^2*R[2])))
    FdF_e(x) = (R = airy_R(x^2);        (x^2*R[1], 2x*(R[1] + x^2*R[2])))
    function x_at(FdF, target, xguess)
        local x, F, dF, r, xn
        x = xguess
        for _ in 1:60
            F, dF = FdF(x)
            r = F - target
            (abs(r) <= 1e-14*max(abs(target), 1.0) || dF <= 0) && break
            xn = x - r/dF
            xn <= 0 && (xn = 0.5x)
            abs(xn - x) <= 1e-13*x && (x = xn; break)
            x = xn
        end
        return x
    end
    #  ⚠ `local` IS LOAD-BEARING — the third instance of this hazard in this
    #  codebase. In Julia a nested function assigning a name that is already local
    #  to the ENCLOSING function assigns the ENCLOSING variable. Both closures
    #  below naturally want to call their phase speed `C`, and `C` is also the
    #  value `wave_properties` returns. Without `local`, the last γ evaluation —
    #  `cg_airy`, at the perturbed depth — silently overwrote the returned `C`
    #  with the AIRY celerity, so `|C/Ce − 1|` collapsed to ~1e-12 for every mesh
    #  at every kd and the applicable range came out as the search cap instead of
    #  10.84. The failure was invisible in C_g and γ, which reproduced their
    #  published values exactly throughout.
    cg_model(dd, ω) = begin
        local x, R, dR, Cl
        x = x_at(FdF_m, ω^2*dd/g, kd)
        R, dR, _ = model_R(vert, x^2)
        Cl = sqrt(g*dd*R)
        (Cl, Cl*(1 + x^2*dR/R))
    end
    cg_airy(dd, ω) = begin
        local x, R, dR, Cl
        x = x_at(FdF_e, ω^2*dd/g, kd)
        R, dR = airy_R(x^2)
        Cl = sqrt(g*dd*R)
        (Cl, Cl*(1 + x^2*dR/R))
    end
    #  C and C_g at the REQUESTED kd are direct evaluations — evaluate them
    #  directly. Routing them through `x_at` (as this did) makes them depend on a
    #  root solve that DEGENERATES exactly where the model does: F(x) = x²R(x²)
    #  SATURATES at Φᵀ|B|⁻¹Φ for large kd, so bisecting F(x)=target on a flat
    #  function returns an arbitrary point of the flat region. The symptom was
    #  precise and misleading — |C/Ce − 1| decayed to 1e-16 as kd grew, i.e. the
    #  model looked PERFECT exactly where it is worst, and the applicable range
    #  came out as the search cap instead of 10.84.
    #  `x_at` is still needed for the γ difference, where d is perturbed at fixed ω
    #  and x genuinely moves — but there it stays in the neighbourhood of kd.
    Rm, dRm, _ = model_R(vert, kd^2)
    Re, dRe    = airy_R(kd^2)
    C   = sqrt(g*d*Rm);  Cg  = C  * (1 + kd^2*dRm/Rm)
    Ce  = sqrt(g*d*Re);  Cge = Ce * (1 + kd^2*dRe/Re)
    ωm = kd * C  / d          # ω = k C, k = kd/d
    ωe = kd * Ce / d
    h = 1e-4
    γm = -0.5*(log(cg_model(d*exp(h), ωm)[2]) - log(cg_model(d*exp(-h), ωm)[2]))/(2h)
    γe = -0.5*(log(cg_airy( d*exp(h), ωe)[2]) - log(cg_airy( d*exp(-h), ωe)[2]))/(2h)
    return (C=C, Cg=Cg, gamma=γm, Ce=Ce, Cge=Cge, gamma_e=γe)
end

"""
    property_errors(vert, kd; g=g) → (eC, eCg, egamma)

The three error measures the applicable range is defined on, **each already in the
units of its own tolerance**: relative for `C` and `C_g`, **ABSOLUTE for `γ`**
(which changes sign near `kd ≈ 1.2`, where a relative measure is singular).

Returns `(Inf, Inf, Inf)` on a degenerate mesh, so an optimiser is repelled rather
than crashed.
"""
function property_errors(vert, kd::Float64; g::Float64=g)
    w = try
        wave_properties(vert, kd; g=g)
    catch
        return (Inf, Inf, Inf)
    end
    (isfinite(w.C) && isfinite(w.Cg) && w.C > 0 && w.Cg != 0) || return (Inf, Inf, Inf)
    return (abs(w.C/w.Ce - 1.0), abs(w.Cg/w.Cge - 1.0), abs(w.gamma - w.gamma_e))
end

"""
    applicable_range(vert; prop=:C, tol=0.02, kd_max=400.0, n=4000, refine=1e-4)

`kd_app^(X) = max{K : |error_X(kd)| ≤ tol ∀ kd ≤ K}` — the definition of
`StokesWaveFourierAnalysis.tex` eq: applicable range definition.

> ⚠ **NOTE THE QUANTIFIER: the tolerance must hold THROUGHOUT `[0,K]`, not merely
> AT `K`.** The error curves are NOT monotone — for surface-clustered meshes `|e|`
> descends to an interior extremum, recovers, and descends again — so a `findlast`
> over a grid (which is what [`applicable_kd`](@ref) does) can step straight over a
> breach and report a range two to three times too large. That is not hypothetical:
> the two-element mesh `c₁ = 0.8696` has an interior dip that just breaches 2 % near
> `kd ≈ 8`, and its true range is 7.7 while a `findlast` reports > 20.
>
> **A coarse scan reproduces the same error**, because the dip can be narrower than
> the grid step. `n` here is deliberately large and the crossing is bisected.

`prop ∈ (:C, :Cg, :gamma)`; the `γ` tolerance is absolute.
"""
function applicable_range(vert; prop::Symbol=:C, tol::Float64=0.02, g::Float64=g,
                          kd_max::Float64=400.0, n::Int=4000, refine::Float64=1e-4)
    idx = prop === :C ? 1 : prop === :Cg ? 2 : prop === :gamma ? 3 :
          error("applicable_range: prop must be :C, :Cg or :gamma (got :$prop)")
    e(kd) = property_errors(vert, kd; g=g)[idx]
    kd0 = 1e-3
    e(kd0) > tol && return 0.0
    #  Geometric grid: the interesting structure is at small kd for gamma and at
    #  large kd for C, and a uniform grid resolves neither well at fixed cost.
    grid = exp.(range(log(kd0), log(kd_max); length=n))
    lo = kd0
    for kd in grid
        if e(kd) > tol
            hi = kd
            while hi - lo > refine
                mid = 0.5*(lo+hi)
                e(mid) > tol ? (hi = mid) : (lo = mid)
            end
            return lo
        end
        lo = kd
    end
    return kd_max
end

"""
    dispersion_error(vert, g, d, kd) → |Cm/Ce − 1|

The dispersion error at ONE `kd`, as a scalar and CONTINUOUS in the vertical mesh —
the objective a σ-mesh optimiser needs. `applicable_kd` below returns a grid index
and is therefore piecewise constant in `c_bdy`; this is not.

Returns `Inf` where `Cm²` is non-positive or the solve fails, so an optimiser
minimising it is repelled from degenerate meshes rather than crashing on them.
"""
function dispersion_error(vert, g::Float64, d::Float64, kd::Float64)
    k  = kd / d
    Ce = sqrt(g * tanh(kd) / k)
    M_eff = vert.Mmat .- vert.B .* kd^2
    Cm_sq = try
        g * d * dot(vert.Phi, M_eff \ vert.Phi)
    catch
        return Inf
    end
    (isfinite(Cm_sq) && Cm_sq > 0) || return Inf
    return abs(sqrt(Cm_sq) / Ce - 1.0)
end

"""
    applicable_kd_first(vert, g, d; err=0.02, kd_max=200.0, n=400, tol=1e-4) → Float64

The **FIRST-CROSSING** applicable `kd`: the smallest `kd` at which `|Cm/Ce − 1|`
exceeds `err`, located by a coarse scan and refined by bisection.

⚠ **This is NOT the same quantity as [`applicable_kd`](@ref)**, and the difference
matters whenever a mesh's dispersion error leaves the ±err band and re-enters.
`applicable_kd` takes `findlast(|ratio−1| ≤ err)` over its grid — the LAST grid
point inside the band, which happily jumps *past* an excursion — and, being a grid
index, it is piecewise constant in `c_bdy`.

Use this one to OPTIMISE a σ-mesh (continuous objective, and it is the honest
statement of the range over which the model is valid); use `applicable_kd` to
reproduce the published tables, which is how they were defined. Quote both when
they disagree — the disagreement IS the finding, namely that the mesh has an
excursion.

Returns `kd_max` if the error never exceeds `err` on `(0, kd_max]`, and `NaN` if it
is already exceeded at the smallest sampled `kd`.
"""
function applicable_kd_first(vert, g::Float64, d::Float64;
                             err::Float64=0.02, kd_max::Float64=200.0,
                             n::Int=400, tol::Float64=1e-4)
    kd0 = 0.01
    dispersion_error(vert, g, d, kd0) > err && return NaN   # bad already at the bottom
    grid = LinRange(kd0, kd_max, n)
    lo   = kd0
    for kd in grid
        e = dispersion_error(vert, g, d, Float64(kd))
        if e > err
            #  bracketed: [lo, kd] straddles the crossing. Bisect.
            hi = Float64(kd)
            while hi - lo > tol
                mid = 0.5*(lo + hi)
                dispersion_error(vert, g, d, mid) > err ? (hi = mid) : (lo = mid)
            end
            return lo
        end
        lo = Float64(kd)
    end
    return kd_max        # never left the band on this range
end

"""
    make_sponge(domain, wL, wR, wB, wT, mu_max)

- domain: ((x0,x1),(y0,y1)) or (x0,x1,y0,y1) rectangle.
- wL:     width of the LEFT sponge layer.
- wR:     width of the RIGHT sponge layer.
- wB:     width of the BOTTOM sponge layer.
- wT:     width of the TOP sponge layer.
- mu_max: maximum sponge strength.

Quadratic sponge μ(x) on up to four boundaries; corner regions clamp at
mu_max (max of the contributions, not the sum).
"""
function make_sponge(domain::Tuple, wL::Float64, wR::Float64,
                         wB::Float64, wT::Float64, mu_max::Float64)
    if domain isa Tuple{Tuple,Tuple}
        (x0,x1), (y0,y1) = domain
    else
        x0,x1,y0,y1 = domain
    end
    function mu_fn(x)
        xv = Float64(x[1]); yv = Float64(x[2])
        mu = 0.0
        wL > 0 && xv < x0 + wL && (mu = max(mu, mu_max * ((x0+wL-xv)/wL)^2))
        wR > 0 && xv > x1 - wR && (mu = max(mu, mu_max * ((xv-(x1-wR))/wR)^2))
        wB > 0 && yv < y0 + wB && (mu = max(mu, mu_max * ((y0+wB-yv)/wB)^2))
        wT > 0 && yv > y1 - wT && (mu = max(mu, mu_max * ((yv-(y1-wT))/wT)^2))
        return mu
    end
    return mu_fn
end

"""
    make_wavemaker_line(x_wm, A, T, k_wave; sigma_wm=1.5)

Gaussian line mass source at x = x_wm (long-crested plane waves):
  S(x,t) = 2Aω exp(−((x−x_wm)/σ)²) cos(ωt)
(factor 2: waves radiate in ±x). Enters continuity as −∫ q·S.
"""
function make_wavemaker_line(x_wm::Float64, A::Float64, T::Float64,
                                 k_wave::Float64; sigma_wm::Float64=1.5)
    omega = 2.0 * pi / T
    function wm_fn(x, t)
        xv = Float64(x[1])
        return 2.0 * A * omega * exp(-((xv - x_wm)/sigma_wm)^2) * cos(omega * t)
    end
    return wm_fn
end

"""
    make_wavemaker_point(x_wm, y_wm, A, T; sigma_wm=1.5)

Gaussian point mass source at (x_wm, y_wm) (ring waves).
"""
function make_wavemaker_point(x_wm::Float64, y_wm::Float64,
                                  A::Float64, T::Float64; sigma_wm::Float64=1.5)
    omega = 2.0 * pi / T
    function wm_fn(x, t)
        xv = Float64(x[1]); yv = Float64(x[2])
        r2 = ((xv - x_wm)^2 + (yv - y_wm)^2) / sigma_wm^2
        return 2.0 * A * omega * exp(-r2) * cos(omega * t)
    end
    return wm_fn
end

# Optimised vertical ELEMENT-BOUNDARY positions for the PIECEWISE-LINEAR family
# (Yang & Liu 2024, Table 1). These are element boundaries (M+1 of them), not
# σ-NODES: for p ≥ 2 each element carries p+1 nodes.
#
# ⚠ KEPT AS THE PUBLISHED VALUES, NOT REPLACED by our own optimiser output.
# `vopt.jl` reproduces them to ≤0.014 from the Yang & Liu total-relative-error
# functional (see VOPT_KAPPA), which is the validation of that implementation —
# but the published figures remain the reference for p = 1, so the table is not
# perturbed by our quadrature and search tolerances.
const DEFAULT_CBDY = Dict(
    1 => [0.0, 1.0],
    2 => [0.0, 0.728, 1.0],
    3 => [0.0, 0.726, 0.925, 1.0],
    4 => [0.0, 0.745, 0.923, 0.977, 1.0],
)

# Optimised element boundaries for p ≥ 2, keyed (M, p). DERIVED HERE by the same
# functional (`optimised_cbdy` in vopt.jl) — no published optimum exists for
# p ≥ 2. The design band Ω is the calibrated fixed point Ω·Δσ_top/p = VOPT_KAPPA.
# Regenerate with `optimised_cbdy(M, p)`; the values are recorded rather than
# recomputed because the search costs minutes and every MMS study needs them.
const DEFAULT_CBDY_P = Dict(
    (1, 2) => [0.0, 1.0],                    # single element: no free interface
    (2, 2) => [0.0, 0.8298, 1.0],
)

"""
    resolve_cbdy(M, c_bdy) → Vector{Float64}

The ONE place the σ-element boundaries are chosen: `c_bdy === nothing` picks the
optimised set for this `M` when one is tabulated and a uniform split of `[0,1]`
otherwise; anything else is passed through after a length/endpoint check.

Every entry point that builds vertical tensors (`setup_and_run`,
`setup_and_run_distributed`, and the MMS drivers) resolves through here, so a
default that is valid for one `M` cannot be hard-wired into a signature that
advertises `M` as a degree of freedom. That defect was live in
`run_conv_study`, where `c_bdy` was pinned to the **M=2** node set and any
`M ≠ 2` threw the `length(c_bdy) == M+1` assertion in `assemble_vertical_tensors`.
"""
function resolve_cbdy(M::Int, c_bdy::Union{Nothing,AbstractVector{<:Real}} = nothing,
                      p::Int = 1)
    M ≥ 1 || error("resolve_cbdy: M must be ≥ 1 (got $M)")
    #  p is a TRAILING POSITIONAL with a default, so every existing 2-argument
    #  call site keeps its meaning (p = 1, the published table). For p ≥ 2 the
    #  optimum is a different mesh — the p=1 table is a valid mesh at any p but
    #  is NOT optimal there, which is exactly what DEFAULT_CBDY_P fixes.
    cb = c_bdy === nothing ?
         (p == 1 ? get(DEFAULT_CBDY, M, collect(LinRange(0.0, 1.0, M + 1))) :
                   get(DEFAULT_CBDY_P, (M, p),
                       get(DEFAULT_CBDY, M, collect(LinRange(0.0, 1.0, M + 1))))) :
         collect(Float64, c_bdy)
    length(cb) == M + 1 || error(
        "resolve_cbdy: c_bdy has $(length(cb)) entries but M=$M needs M+1 = $(M+1) " *
        "element boundaries (they are boundaries, not σ-nodes — for p ≥ 2 the extra " *
        "nodes are interior to each element).")
    (cb[1] ≈ 0.0 && cb[end] ≈ 1.0) ||
        error("resolve_cbdy: c_bdy must start at 0 and end at 1 (got $(cb[1]) … $(cb[end])).")
    issorted(cb) || error("resolve_cbdy: c_bdy must be increasing (got $cb).")
    return cb
end

"""
    setup_and_run(; kwargs...) → (diags, vert, prob)

Top-level sequential driver. It executes the whole workflow in order —
Stage 1 (vertical tensors), the horizontal mesh and FE spaces, the forcing
(wavemaker/sponge or Dirichlet wave generation), the residual/operator, the
time integrator, the initial state, and finally the time loop — and returns the
per-step diagnostics `diags`, the vertical tensor bundle `vert`, and the problem
bundle `prob`. Every physical/numerical choice is a keyword argument (documented
inline on the signature below); the defaults describe a small plane-wave case.

Wave forcing comes from EITHER an internal wavemaker OR Dirichlet boundary
generation, chosen by the arguments:
  * internal wavemaker: `y_wm = nothing` → line source (long-crested plane
    wave); a number → point source at (x_wm, y_wm) (radial/ring waves);
  * Dirichlet boundary generation (waveinput.jl): pass `wave_bc` =
  * `:regular`             — monochromatic wave built from `A_wave`/`T_wave`;
  * a `WaveInput`          — any prebuilt component table;
  * a `WaveSpec.AiryWaves.AiryState` — stochastic sea state (auto-converted).
The interior wavemaker is then disabled and η/𝖴x (and 𝖴y for directional
seas) are prescribed on `bc_side` (`:left`/`:right`). Related kwargs:
`bc_profile` (`:model`/`:airy` vertical polarization), `T_ramp` (Hann ramp,
`nothing` → 2 peak periods), `ic_from_bc` (hot start from the incident field;
requires `T_ramp=0.0`), `relax_bc`+`relax_width` (generation/absorption
relaxation zone adjacent to the inflow, strength `mu_max`).
"""

"""
    check_flat_bed_consistency(h_bathy, domain, flat_bed; ns=7, rtol=1e-8) → Bool

Warn if the `flat_bed` switch disagrees with the prescribed bathymetry `h_bathy`
(`x → d(x)`) sampled on a small grid over `domain`: `flat_bed=true` over a varying
bed silently drops the ∇h (sloping-bed) physics, while `flat_bed=false` over a
constant bed assembles ∇h-terms that vanish anyway. Returns whether the bed is
constant. Sampling a pure function is cheap and rank-independent.
"""
function check_flat_bed_consistency(h_bathy, domain, flat_bed::Bool;
                                    ns::Int = 7, rtol::Float64 = 1e-8)
    (x0, x1), (y0, y1) = domain isa Tuple{Tuple,Tuple} ? domain :
                         ((domain[1], domain[2]), (domain[3], domain[4]))
    hs = [h_bathy(VectorValue(x, y)) for x in range(x0, x1; length=ns),
                                         y in range(y0, y1; length=ns)]
    hmin, hmax = extrema(hs); hmean = sum(hs) / length(hs)
    is_const = (hmax - hmin) <= rtol * max(abs(hmean), 1.0)
    if flat_bed && !is_const
        @warn "flat_bed=true but the bathymetry varies over the domain " *
              "(Δh=$(round(hmax - hmin; sigdigits=3)) m): the ∇h (sloping-bed) terms are being " *
              "dropped and the bed-slope physics is NOT represented — use flat_bed=false for variable bathymetry."
    elseif !flat_bed && is_const
        @info "flat_bed=false but h_bathy is constant (d=$(round(hmean; sigdigits=4)) m): " *
              "the ∇h terms vanish anyway; set flat_bed=true to skip assembling them."
    end
    return is_const
end

"""
    resolve_wave_gen(wave_gen, wave_bc) → Symbol

Map the user's `wave_gen` selector to one of the two mechanisms `:inner_res |
:bc_gen`, validating it against `wave_bc`:

  - `:inner_res` — interior Gaussian source (line ⇒ plane wave, point ⇒ ring wave);
                   requires `wave_bc === nothing`.
  - `:bc_gen`    — boundary Dirichlet generation. The boundary source — a
                   parametrised regular wave (from `A_wave`/`T_wave`/`wave_dir`),
                   a prebuilt `WaveInput`, or a WaveSpec `AiryState` — is dispatched
                   on the TYPE of `wave_bc` inside the driver; all feed the SAME
                   Dirichlet machinery (they differ only in how the `WaveInput`
                   component table is populated, not in the boundary mechanism).

`:auto` (the default) infers it: `wave_bc === nothing` → `:inner_res`, else `:bc_gen`.
"""
function resolve_wave_gen(wave_gen::Symbol, wave_bc)
    if wave_gen === :auto
        return wave_bc === nothing ? :inner_res : :bc_gen
    elseif wave_gen === :inner_res
        wave_bc === nothing ||
            error("wave_gen=:inner_res uses the interior source; leave wave_bc=nothing")
        return :inner_res
    elseif wave_gen === :bc_gen
        return :bc_gen
    else
        error("wave_gen must be :auto, :inner_res or :bc_gen (got :$wave_gen)")
    end
end

function setup_and_run(;
    # ---- Vertical discretisation ---------------------------------------------
    M            :: Int     = 2,          # number of vertical σ-elements (BALFE-M order: 2/3/4)
    p_vertical   :: Int     = 1,          # polynomial order of each σ-element (Nσ = M·p_vertical+1)
    c_bdy                   = nothing,    # σ-node boundary positions in [0,1]; nothing → optimised set
    # ---- Horizontal discretisation -------------------------------------------
    domain                  = ((0.0, 60.0), (0.0, 10.0)),  # ((x0,x1),(y0,y1)) extent [m]
    partition    :: Tuple   = (120, 20),  # (nx,ny) number of horizontal cells
    p_u :: Int     = 2,          # horizontal FE order (must be ≥2: Q1 zeroes the dispersion)
    quad_extra   :: Int     = 0,          # EXTRA quadrature degree on top of the default
                                          #   2·max(p_h,p_η)+2. The default integrates the LINEAR
                                          #   terms exactly but is one degree SHORT of the nonlinear
                                          #   advection integrand φᵢ·u_k·∇u_j·H (degree 3p+1 at
                                          #   equal order). Raise it to test aliasing hypotheses.
    p_eta        :: Int     = p_u - 1,    # surface FE order. MUST satisfy p_u = p_eta + 1
                                          #   (Taylor-Hood); check_taylor_hood ERRORS otherwise.
                                          #   Equal order is inf-sup deficient here — see rule 2b
                                          #   pairing: η enters momentum undifferentiated (via ∇·v
                                          #   after IBP), so it plays the pressure role of a Stokes
                                          #   system and equal-order continuous spaces are inf-sup
                                          #   deficient — the analytic MMS measures order p there
                                          #   rather than p+1, in BOTH fields.
                                          #   ⚠ A BETTER RATE IS NOT A BETTER ANSWER AT A GIVEN MESH:
                                          #   at nx=24 the equal-order Q3/Q3 was 40x MORE ACCURATE
                                          #   than Q3/Q2, because η sits in a richer space. Compare
                                          #   error-vs-DOF at your production resolution before
                                          #   switching. See MMS_CONVERGENCE_CAMPAIGN.md.
    # ---- Physical parameters -------------------------------------------------
    h_val        :: Float64 = 3.5,        # still-water depth [m] (flat bed unless h_bathy given)
    g            :: Float64 = g,          # gravitational acceleration [m/s²]
    T_wave       :: Float64 = 1.6,        # forcing wave period [s]
    A_wave       :: Float64 = 0.001,      # forcing wave amplitude [m] (keep small for stability)
    # ---- Internal wavemaker (used only when wave_bc is nothing) ---------------
    x_wm         :: Float64 = 12.0,       # wavemaker x-position [m]
    y_wm                    = nothing,    # nothing → line source (plane wave); number → point source
    # ---- Sponge layers (quadratic damping toward each boundary) ---------------
    sponge_wL    :: Float64 = 12.0,       # left-edge sponge width [m] (0 = none)
    sponge_wR    :: Float64 = 12.0,       # right-edge sponge width [m]
    sponge_wB    :: Float64 = 0.0,        # bottom-edge (y0) sponge width [m]
    sponge_wT    :: Float64 = 0.0,        # top-edge (y1) sponge width [m]
    mu_max       :: Float64 = 5.0,        # peak sponge strength (also the relaxation-zone strength)
    # ---- Time integration ----------------------------------------------------
    T_final      :: Float64 = 12.8,       # final simulated time [s]
    dt           :: Float64 = 0.02,       # time step [s]
    solver_type  :: Symbol  = :sdirk,     # integrator: :sdirk (default) | :theta | :gen_alpha | :rk3
    tableau      :: Symbol  = :SDIRK_2_2, # Runge–Kutta tableau when solver_type == :sdirk
    theta        :: Float64 = 0.5,        # θ for :theta (0.5 = Crank–Nicolson)
    rho_inf      :: Float64 = 0.5,        # high-frequency damping for :gen_alpha (ρ∞∈[0,1])
    # ---- Output --------------------------------------------------------------
    output_dir   :: String  = joinpath(@__DIR__, "..", "output", "seq_out"),  # VTK/pvd destination
    save_every   :: Int     = 0,          # write a VTK snapshot every N steps (0 = no VTK)
    gauges                  = [],         # list of (x,y) probe points; η is sampled there each step
    # ---- Boundary conditions -------------------------------------------------
    y_wall_bc    :: Symbol  = :wall,      # y-edge BC: :wall (𝖴y=0) | :open (natural) | :periodic (y-periodic)
    x_wall_bc    :: Bool    = false,      # solid walls on the x-edges (𝖴x=0); true for closed-basin IC
    # ---- Physics flags (switch individual residual terms on/off) --------------
    regime       :: Symbol  = :nonlinear, # :linear (linearised, no advection) | :nonlinear
    nl_pressure  :: Symbol  = :none,      # nonlinear pressure: :none | :native {3,6,7,8} | :full {+1,2,4,5}
    flat_bed     :: Bool    = false,      # sea-bed geometry: false = variable bathymetry (∇h≠0),
                                          #   true = flat bed (∇h≡0; every ∇h-term dropped, ∇η-terms kept)
    h_bathy                 = nothing,    # x → d(x): variable bathymetry (overrides h_val)
    eta0_func               = nothing,    # x → η₀(x): initial free surface (IC release: set x_wall_bc=true)
    # ---- Dirichlet boundary wave generation (waveinput.jl) --------------------
    wave_gen     :: Symbol  = :auto,      # :inner_res | :bc_gen  (:auto infers from wave_bc)
    wave_bc                 = nothing,    # nothing | :regular | WaveInput | WaveSpec AiryState
    wave_dir     :: Float64 = 0.0,        # propagation angle vs +x for a :bc_gen boundary wave [rad]
    bc_side      :: Symbol  = :left,      # generation boundary (:left/:right)
    bc_profile   :: Symbol  = :model,     # vertical polarization of the inflow (:model/:airy)
    T_ramp                  = nothing,    # Hann ramp-up time [s]; nothing → 2 peak periods
    ic_from_bc   :: Bool    = false,      # hot-start the field from the incident wave (needs T_ramp=0)
    relax_bc     :: Bool    = false,      # relaxation (generation/absorption) zone at the inflow
    relax_width  :: Float64 = 0.0,        # relaxation-zone width [m]; 0 → one peak wavelength
    # ---- Solver / diagnostics ------------------------------------------------
    use_ad       :: Bool    = false,      # build Jacobians by AD instead of the hand Jacobians
    show_trace   :: Bool    = false,      # print the Newton iteration trace
    nl_iter      :: Int     = 50,         # max Newton iterations per stage
    nl_tol       :: Float64 = 1e-5,       # Newton residual tolerance (‖r‖∞). Production default.
                                          #   MEASURED 2026-08-12, and it is NOT a null change:
                                          #   Newton takes an INTEGER number of iterations, so this
                                          #   tolerance acts as a step function. At 1e-6 Newton needs
                                          #   3 iterations/step; at 1e-5 (and at 1e-4 — bit-identical)
                                          #   it needs 2, and max|eta| shifts by 3.8e-5 relative.
                                          #   Accepted on the error budget, NOT on a null result: the
                                          #   dropped 3rd iteration polishes a residual (~3e-6) some
                                          #   600x smaller than the O(dt^2) time-discretisation error
                                          #   the answer already carries (‖R‖∞ ~ 1.8e-3, measured by
                                          #   the run's own residual check). Tests that need a sharper
                                          #   answer pin 1e-8 explicitly.
    print_every  :: Int     = 1,          # print a step report every N steps (1 = every step)
    check_every  :: Int     = 50,         # re-verify the governing equations every N steps (0 = off)
    check_tol    :: Float64 = 1e-8,       # tolerance for that verification (‖R‖∞)
    # ---- Field diagnostics (monitor.jl: max|η| location, invariants, RSS) -----
    diag_every   :: Int     = 0,          # sample every N steps (0 → = print_every; −1 = disabled)
    diag_csv     :: Bool    = true,       # write output_dir/diagnostics.csv (needs save_every≥0 dir)
    eta_ref                 = nothing,    # reference amplitude for the divergence guard (auto)
    div_factor   :: Float64 = 20.0,       # abort when max|η| > div_factor · eta_ref
    # ---- Reconstructed field output (at the Nσ vertical σ-nodes) --------------
    write_w        :: Bool    = false,    # also write vertical-velocity fields w_s<σ> to VTK
    write_pressure :: Bool    = false,    # also write total-pressure fields p_s<σ> to VTK
    rho            :: Float64 = rho,      # water density [kg/m³] (used for the pressure output)
)
    # Choose the σ-element boundaries: the paper's optimised set for this M when
    # available, otherwise a uniform split of [0,1]. One resolver, one definition.
    #  Pairing gate FIRST — before the σ-tensors, the mesh, the JIT and the hours.
    #  build_fe_spaces gates it too, but that is ~30 min of compilation away on a cold
    #  session, and a run should not get that far only to be told its element orders
    #  were wrong. CLAUDE.md rule 2b.
    check_taylor_hood(p_u, p_eta; where = "setup_and_run")
    c_bdy = resolve_cbdy(M, c_bdy)

    # --- STAGE 1: VERTICAL PRE-COMPUTATION (MESH INDEPENDENT, DONE ONCE) -------
    # Build the σ-basis and integrate every vertical tensor the residual needs.
    println("=== Vertical FE problem (algebraic BALFE-M) ===")
    vert = assemble_vertical_tensors(M, p_vertical, c_bdy)
    @printf("  Nσ=%d   ΣΦ=%.6f\n", vert.N_dof, sum(vert.Phi))   # ΣΦ=1 sanity check

    # Lateral (y) boundary condition: :wall / :open / :periodic.
    y_wall_bc in (:wall, :open, :periodic) ||
        error("setup_and_run: y_wall_bc must be :wall, :open or :periodic (got :$y_wall_bc)")
    y_periodic = y_wall_bc == :periodic  # set y_periodic = true   when  y_wall_bc == :periodic
    if y_periodic
        # Check for bottom/top sponge layers, incompatible with y-periodic BCs (they are redundant).
        (sponge_wB > 0 || sponge_wT > 0) &&
            @warn "y_wall_bc=:periodic — lateral sponges (sponge_wB/wT) are redundant on a " *
                  "periodic domain; set them to 0"
        # Check for a point source wavemaker, which is not y-periodic (the periodic image array is not a point source).
        !isnothing(y_wm) &&
            @warn "y_wall_bc=:periodic with a point source is not y-periodic (periodic image " *
                  "array); use a line source (y_wm=nothing)"
    end

    # --- STAGE 2 SETUP: HORIZONTAL MESH + INTEGRATION MEASURE ----------------- 
    # `y_periodic` glues the top/bottom edges when y_wall_bc == :periodic.
    model, trian = build_horizontal_model(domain, partition; y_periodic=y_periodic)
    # quadrature degree = 2·p_u+2 integrates the nonlinear (product) terms exactly enough.
    dΩh = Measure(trian, 2*max(p_u, p_eta) + 2 + quad_extra)

    # Forcing frequency and the matching wavenumber from the Airy relation
    # ω² = g k tanh(kd) (used to size the wavemaker and report kd).
    omega  = 2.0*pi/T_wave                      # compute wave frequency from the period
    k_wave = find_wavenumber(omega, h_val, g)   # compute corresponding wavenumber from the Airy dispersion relation

    # dfn: the bathymetry function: either a constant h_val or a user-supplied h_bathy(x).
    dfn    = isnothing(h_bathy) ? (x -> h_val) : h_bathy    # bathymetry: constant h_val or user d(x)

    # Unpack the domain corners (accept either nested or flat tuple form).
    if domain isa Tuple{Tuple,Tuple}
        (x0d, x1d), (y0d, y1d) = domain
    else
        x0d, x1d, y0d, y1d = domain
    end

    # ---- DIRICHLET BOUNDARY WAVE GENERATION (WAVEINPUT.JL) --------------------
    # If wave_bc is set, build the component table `wi` that drives the inflow
    # boundary; the interior wavemaker is disabled below. `wi` stays nothing for
    # the internal-wavemaker path.
    wg = resolve_wave_gen(wave_gen, wave_bc)
    @printf("  Wave generation: %s\n", string(wg))
    wi = nothing
    if wg === :bc_gen
        bc_side in (:left, :right) ||
            error("setup_and_run: bc_side must be :left or :right (got :$bc_side)")
        Tr = T_ramp === nothing ? nothing : Float64(T_ramp)

        # The boundary source is dispatched on the TYPE of wave_bc.
        if wave_bc isa WaveInput  # a prebuilt WaveInput passes through
            wi = wave_bc
        elseif wave_bc isa WaveSpec.AiryWaves.AiryState  # a WaveSpec AiryState is converted
            wi = WaveInput(vert, wave_bc; d=h_val, g=g, T_ramp=Tr, profile=bc_profile)
        else     # otherwise (nothing / :regular) a parametrised regular plane wave is built from A_wave/T_wave/wave_dir. 
            tr_val = Tr === nothing ? 2.0 * T_wave : Tr
            wi = WaveInput(vert; A=A_wave, T=T_wave, d=h_val, g=g, theta=wave_dir,
                        T_ramp=tr_val, profile=bc_profile)
        end
        # All three feed the same Dirichlet machinery, which needs to be feed the component table `wi`,
        # a WaveInput structure (the only difference is how the table is populated).

        println()
        waveinput_summary(wi)      # print Hs/Tp/components of the generated sea
        # Sanity: the generation boundary must sit over a constant depth equal to
        # the depth the boundary data was built for.
         
        # Identify the x-coordinate of the generation boundary (x0d or x1d) 
        xg = bc_side == :left ? x0d : x1d

        # Sample the bathymetry along the y-direction at that x coordinate (xg = generation BC position). The depth must be constant
        dsamp = [dfn(VectorValue(xg, y0d + s*(y1d - y0d))) for s in 0.0:0.25:1.0]

        # Check if all sampled depths are approximately equal to the WaveInput depth (wi.d) within a relative tolerance of 1e-8. If not, issue a warning.
        all(v -> isapprox(v, wi.d; rtol=1e-8), dsamp) ||
            @warn "wave_bc: depth along the generation boundary is not constant " *
                  "(or differs from the WaveInput depth $(wi.d) m)"

        # A directional (θ≠0) sea needs open/periodic lateral boundaries, not walls.
        wi.directional && y_wall_bc == :wall &&
            error("setup_and_run: a directional sea (θ≠0 components) requires " *
                  "y_wall_bc=:open (lateral sponges) or :periodic")
        wi.directional && y_periodic &&
            @warn "y_wall_bc=:periodic with a directional sea requires each component's " *
                  "transverse wavenumber k·sinθ to be a box harmonic 2π/Ly (not enforced)"

        # A hot start supplies the field at t=0, so it is incompatible with a ramp.
        ic_from_bc && wi.T_ramp > 0.0 &&
            error("setup_and_run: ic_from_bc=true requires T_ramp=0.0 " *
                  "(the hot start replaces the ramp)")

        # Warn if a sponge sits on the inflow (it would eat the incoming wave)…
        gen_w = bc_side == :left ? sponge_wL : sponge_wR
        gen_w > 0.0 && !relax_bc &&
            @warn "wave_bc: a plain sponge overlaps the generation boundary and " *
                  "damps the incident wave; set its width to 0 or use relax_bc=true"
        # …and if there is no sponge on the far side to absorb the outgoing wave.
        opp_w = bc_side == :left ? sponge_wR : sponge_wL
        opp_w > 0.0 ||
            @warn "wave_bc: no sponge opposite the inflow — expect reflections"
    end

    # --- STACKED FE SPACES [η,𝖴x,𝖴y] + BOUNDARY CONDITIONS wi --------------------
    println("\n=== 2D Horizontal FE problem (stacked [η,𝖴x,𝖴y]) ===")
    # For a generated sea, pass the time-varying Dirichlet data (η, 𝖴x, and 𝖴y
    # for directional seas) so build_fe_spaces makes the matching transient trials.

    # Build inflow BC data for build_fe_spaces. 
    # If wi is nothing, the interior wavemaker is used and no Dirichlet BCs are applied. 
    if wi === nothing
        inflow = nothing
    # Otherwise, the inflow BCs are built from the WaveInput structure wi.
    else
        if wi.directional
            uy_val = uy_bc(wi)
        else
            uy_val = nothing
        end
        inflow = (side = bc_side, eta = eta_bc(wi), ux = ux_bc(wi), uy = uy_val)
    end

    # Build the stacked FE spaces for the horizontal problem, applying the inflow BCs if provided.
    pe = p_eta
    U, V = build_fe_spaces(model, 
                                p_u,           # horizontal (velocity) FE order
                                vert.N_dof;             # number of vertical DOFs = number of stacked fields
                                y_wall_bc=y_wall_bc,    # lateral BC type
                                x_wall_bc=x_wall_bc,    # solid wall BC on x-edges
                                inflow=inflow,          # inflow BC data (η, 𝖴x, 𝖴y) if provided
                                p_eta=pe)               # surface FE order (see the kwarg note)

    @printf("  Fields: 3 (η + 2 stacked VectorValue{%d})   free DOFs: %d\n",
            vert.N_dof, num_free_dofs(U(0.0)))
    @printf("  Wave: λ=%.2f m, kd=%.2f\n", 2pi/k_wave, k_wave*h_val)

    # --- Forcing: sponge profile + internal wavemaker source ------------------
    # Sponge damping μ(x,y) grows quadratically toward the flagged boundaries.
    sponge = make_sponge(domain, sponge_wL, sponge_wR, sponge_wB, sponge_wT, mu_max)

    # Internal source S(x,t): none when generating at a boundary; a line source
    # (plane waves) when y_wm is nothing; a point source (ring waves) otherwise.
    if wi !== nothing
        wm = (x, t) -> 0.0  # Dirichlet generation disables the interior wavemaker
    elseif isnothing(y_wm)
        wm = make_wavemaker_line(x_wm, A_wave, T_wave, k_wave)  # line source (plane waves)
    else
        wm = make_wavemaker_point(x_wm, Float64(y_wm), A_wave, T_wave)  # point source (ring waves)
    end

    # --- Optional relaxation zone next to the inflow (generation + absorption) -
    # Blends the state toward the incident wave over a strip at the boundary.
    relax_mu_fn = x -> 0.0
    relax_tg    = nothing
    use_relax   = wi !== nothing && relax_bc
    if use_relax
        wrx = relax_width > 0.0 ? relax_width : 2.0*pi/wi.ks[argmax(wi.amps)]
        relax_mu_fn = bc_side == :left ?
            (x -> x[1] < x0d + wrx ? mu_max*((x0d + wrx - x[1])/wrx)^2 : 0.0) :
            (x -> x[1] > x1d - wrx ? mu_max*((x[1] - (x1d - wrx))/wrx)^2 : 0.0)
        relax_tg = incident_fields(wi)
        @printf("  Relaxation zone: width=%.2f m, μ_max=%.2f, %s boundary\n",
                wrx, mu_max, string(bc_side))
    end

    # Warn if the user has set flat_bed=true but the bathymetry varies, or vice versa.
    check_flat_bed_consistency(dfn, domain, flat_bed)   # warn on bed ↔ flat_bed mismatch

    # --- Assemble the problem bundle: vertical tensors → Gridap constants +
    #     the depth, forcing, and physics flags that define the residual. -------
    prob = build_problem(vert; g=g, 
                        h_bathy=dfn,                # bathymetry function (x → d(x))
                        regime=regime,              # linear/nonlinear physics
                        nl_pressure=nl_pressure,    # nonlinear pressure treatment
                        flat_bed=flat_bed,          # whether to drop ∇h terms (flat bed)
                        mu_sponge=sponge,           # sponge damping profile μ(x,y)
                        wm_src=wm,                  # internal wavemaker source S(x,t)
                        relax_bc=use_relax,         # whether to use a relaxation zone at the inflow
                        relax_mu=relax_mu_fn,       # relaxation-zone damping profile μ(x,y)
                        relax_tg=relax_tg)          # incident wave target for the relaxation zone

    # Build problem TransientFEOperator ->  Wrap the residual (+ Jacobians) into a Gridap operator
    # `use_ad` swaps the hand Jacobians for AD-generated ones (cross-checking only).
    op = use_ad ? build_ode_operator_ad(prob, U, V, trian, dΩh) :
                  build_ode_operator(prob, U, V, trian, dΩh)

    # The monitor transparently wraps the Newton solver to harvest per-step stats.
    monitor = SolverMonitor()

    # Build the time integrator (SDIRK by default) around a Newton+LU nonlinear solve.
    solver  = build_ode_solver(dt;  solver_type=solver_type,    # :sdirk | :theta | :gen_alpha | :rk3
                                    theta=theta,                # θ for :theta (0.5 = Crank–Nicolson)
                                    rho_inf=rho_inf,            # high-frequency damping for :gen_alpha (ρ∞∈[0,1])
                                    tableau=tableau,            # SDIRK tableau for :sdirk
                                    nl_iter=nl_iter,            # max Newton iterations per stage
                                    nl_tol=nl_tol,              # Newton residual tolerance (‖r‖∞)
                                    show_trace=show_trace,      # print the Newton iteration trace
                                    monitor=monitor)            # wrap the Newton solver to harvest per-step stats

    # Optional independent re-assembly of the governing equations for verification
    # -> check if the governing equations are satisfied at the current solution (‖R‖∞ < check_tol).
    # (its θ-scheme self-check is meaningful only for the :theta integrator).
    checker = check_every > 0 ?
              ResidualChecker(prob, U, V, trian, dΩh, dt, theta,
                                 solver_type == :theta) : nothing

    # Initial condition. Four cases: hot-start from the incident wave; rest state
    # for a generated sea; rest state (default); or a prescribed η₀(x) release.
    u0 = if wi !== nothing && ic_from_bc
        inc = incident_fields(wi)
        make_initial_conditions(U(0.0), vert.N_dof;
            eta0_func = x -> inc.eta(x, 0.0),
            ux0_func  = x -> inc.ux(x, 0.0),
            uy0_func  = x -> inc.uy(x, 0.0))
    elseif wi !== nothing
        make_initial_conditions(U(0.0), vert.N_dof; eta0_func=eta0_func)
    elseif isnothing(eta0_func)
        make_initial_conditions(U)
    else
        make_initial_conditions(U, vert.N_dof; eta0_func=eta0_func)
    end

    # For nl_pressure=:full, build the frozen-projection context (mass matrix
    # factorised once) used to evaluate the irreducible ∇H/𝓟 pressure halves.
    nlp = nl_pressure == :full ?
          (prob, build_nlp_ctx(model, p_u, vert.N_dof, trian, dΩh)) : nothing

    # Reconstruction context for optional w/p VTK output (nothing if both off).
    recon = build_field_recon(vert, dfn, g; rho=rho,
                                  write_w=write_w, write_pressure=write_pressure)
    if recon !== nothing
        @printf("  Field output: write_w=%s write_pressure=%s at σ-levels %s\n",
                string(write_w), string(write_pressure),
                string(round.(recon.levels; digits=3)))
    end

    # Field diagnostics: max|η| location + interior/damped split, |u|/|η|, mass
    # and energy invariants, process RSS, and the relative divergence guard.
    # `diag_every=0` follows `print_every` so every reported line is complete;
    # a negative value switches the whole block off.
    diag_n   = diag_every == 0 ? max(print_every, 1) : diag_every
    eta_ref_v = resolve_eta_ref(eta_ref, A_wave, wi, eta0_func, domain)
    rundiag  = diag_n > 0 ?
               build_run_diagnostics(prob, wi !== nothing ? U(0.0) : U, trian, dΩh;
                                     eta_ref=eta_ref_v, div_factor=div_factor,
                                     output_dir=output_dir, diag_csv=diag_csv,
                                     u0=u0) : nothing

    # Print the solver-configuration banner (integrator, tolerances, dt/steps).
    println()
    print_solver_banner(
        @sprintf("Newton (NLsolve, exact hand Jacobians) | max iters = %d | ftol (‖r‖∞) = %.1e",
                 nl_iter, nl_tol),
        "LU direct factorisation (sequential)";
        solver_type=solver_type, theta=theta, dt=dt, t0=0.0, T_final=T_final,
        print_every=print_every, check_every=check_every, check_tol=check_tol,
        monitor=monitor, diag_every=max(diag_n, 0), eta_ref=eta_ref_v,
        div_limit=rundiag === nothing ? NaN : rundiag.div_limit)

    # --- March the transient problem from 0 to T_final, collecting per-step
    #     diagnostics and writing VTK/gauge output as requested. ---------------
    println("\n=== Time loop (algebraic) ===")
    diags = run_time_loop(op, solver,                   # Gridap operator + time integrator, constructed at build_ode_solver
                            u0,                         # initial condition (Gridap FE function)
                            0.0,                        # initial time
                            T_final;                    # final time
                            output_dir=output_dir,      # VTK/pvd destination (output directory)
                            save_every=save_every,      # write a VTK snapshot every N steps (0 = no VTK)
                            trian=trian,                # horizontal mesh (for VTK output)
                            Nσ=vert.N_dof,              # number of vertical σ-levels = number of stacked horizontal fields (+ 1 for η) 
                            print_every=print_every,    # print a step report every N steps (1 = every step)
                            gauges=gauges,              # list of (x,y) probe points; η is sampled there each step
                            recon=recon,                # reconstruction context for optional w/p VTK output (nothing if both off)
                            trial_space=U,              # trial FE space (for VTK output)
                            dt=dt,                      # time step [s]
                            nlp=nlp,                    # frozen-projection context for nl_pressure=:full (nothing if not used)
                            monitor=monitor,            # wrap the Newton solver to harvest per-step stats
                            checker=checker,            # optional independent re-assembly of the governing equations for verification
                            check_every=check_every,    # re-verify the governing equations every N steps (0 = off)
                            check_tol=check_tol,        # tolerance for that verification (‖R‖∞)
                            rundiag=rundiag,            # field-diagnostics state (max|η| location, invariants, RSS)
                            diag_every=max(diag_n, 0))  # sample the field diagnostics every N steps (0 = off)

    return diags, vert, prob
end

# ==============================================================
#  Standardised output-directory naming
#  (building_files/OUTPUT_NAMING_PROPOSAL.md — the accepted spec)
#
#      <model>_<domain>_<wave>_<regime>_<nlp>_<bed>_<discr>_<amplitude>_<period>[_<extra>…]
#      P1LFE-2_1d_bcplane_nl_full_flat_Q2Q1_A0.10_T1.6
#
#  ONE generator for every driver, sequential and distributed. Three drivers building
#  names three different ways is the defect this replaces.
#
#  ⚠ FIELDS ARE NEVER OMITTED, and an absent field must never encode a default. The
#  discretisation field exists because nothing in the old names said Q2/Q1 vs Q2/Q2 —
#  which is exactly how the equal-order runs and the Taylor-Hood MMS campaign were
#  compared for months as though they were the same solver (CLAUDE.md rule 12b).
# ==============================================================

"Vertical basis token: `P{p_vert}LFE-{M}`. The legacy `M2` spelling is retired."
model_token(M::Int, p_vert::Int) = "P$(p_vert)LFE-$(M)"

"""
    domain_token(; ny, y_wall_bc) -> "1d" | "2d" | "2dper" | "2dopen"

`ny == 1` with solid walls is the narrow flume, i.e. a 1-D horizontal case
(CLAUDE.md rule 12). Otherwise the lateral boundary condition names the class.
"""
function domain_token(; ny::Int, y_wall_bc::Symbol)
    ny == 1 && y_wall_bc === :wall && return "1d"
    y_wall_bc === :periodic && return "2dper"
    y_wall_bc === :open     && return "2dopen"
    return "2d"
end

"""
    wave_token(kind, gen) -> e.g. "bcplane", "plane", "irr", "dir", "ichump"

The generation mechanism is a PREFIX on the wave type, not a separate field: `bcplane`
keeps the mechanism attached to the thing it generates and reads better than `plane_bc`.
`gen`: `:bc` | `:inner` | `:ic` | `:none`.
"""
function wave_token(kind::AbstractString, gen::Symbol)
    pre = gen === :bc ? "bc" : gen === :ic ? "ic" : ""
    return pre * kind
end

"`:linear`→`lin`, `:nonlinear`→`nl`."
regime_token(regime::Symbol) = regime === :linear ? "lin" : "nl"

"""
    discr_token(p_u, p_eta; nx=nothing, ny=nothing) -> "Q2Q1" | "Q2Q1-nx480" | …

⚠ `p_eta` is written even though Taylor-Hood makes it redundant in principle. A
discretisation that is not visible in the output is one nobody checks.
"""
function discr_token(p_u::Int, p_eta::Int; nx = nothing, ny = nothing)
    t = "Q$(p_u)Q$(p_eta)"
    nx === nothing && return t
    return ny === nothing || ny == 1 ? "$t-nx$(nx)" : "$t-nx$(nx)x$(ny)"
end

"Trim a float to a compact, unambiguous decimal (`0.1`→`0.10`, `1.6`→`1.6`, `0.001`→`0.001`)."
function _num(x::Real)
    s = @sprintf("%.10g", float(x))
    return s
end

"""
    output_dir_name(; M, p_vert, ny, y_wall_bc, wave_kind, wave_gen, regime,
                      nl_pressure, bed, p_u, p_eta, amplitude, period,
                      irregular=false, nx=nothing, nx_in_name=false, extra=String[]) -> String

Build the standardised directory name. **Refuses combinations that cannot exist**, rather
than labelling them: a name generator that can express an impossible case will eventually
be asked to label one.
"""
function output_dir_name(; M::Int, p_vert::Int,
                           ny::Int, y_wall_bc::Symbol,
                           wave_kind::AbstractString, wave_gen::Symbol,
                           regime::Symbol, nl_pressure::Symbol,
                           bed::AbstractString,
                           p_u::Int, p_eta::Int,
                           amplitude::Real, period::Real,
                           irregular::Bool = false,
                           nx = nothing, nx_in_name::Bool = false,
                           extra::AbstractVector{<:AbstractString} = String[])
    dom = domain_token(; ny = ny, y_wall_bc = y_wall_bc)
    #  ⚠ A 1-D domain carries one propagation direction; directional content has a
    #  transverse wavenumber a one-cell-wide flume cannot represent (rule 12). Refuse it.
    (dom == "1d" && occursin("dir", wave_kind)) && error(
        "output_dir_name: directional content ($wave_kind) is impossible on a 1-D domain " *
        "— a flume one cell across cannot represent k_y (CLAUDE.md rule 12).")
    regime === :linear && nl_pressure !== :none && error(
        "output_dir_name: regime=:linear with nl_pressure=:$nl_pressure — a linear model " *
        "carries no 𝓝 (mirrors resolve_physics).")
    check_taylor_hood(p_u, p_eta; where = "output_dir_name")

    amp = (irregular ? "Hs" : "A") * _num(amplitude)
    per = (irregular ? "Tp" : "T")  * _num(period)
    parts = [model_token(M, p_vert), dom,
             wave_token(wave_kind, wave_gen),
             regime_token(regime), String(nl_pressure), String(bed),
             discr_token(p_u, p_eta; nx = nx_in_name ? nx : nothing, ny = nothing),
             amp, per]
    append!(parts, extra)
    return join(parts, "_")
end

"""
    unique_output_dir(root, name) -> String

Join and, if the directory already exists and is non-empty, suffix `_v2`, `_v3`, …

⚠ **NEVER silently overwrite.** On 2026-09-06 a re-executed batch wrote six finished runs
over their own output and destroyed them. A path that cannot collide cannot do that.
"""
function unique_output_dir(root::AbstractString, name::AbstractString)
    path = joinpath(root, name)
    (!isdir(path) || isempty(readdir(path))) && return path
    v = 2
    while isdir("$(path)_v$(v)") && !isempty(readdir("$(path)_v$(v)"))
        v += 1
    end
    return "$(path)_v$(v)"
end
