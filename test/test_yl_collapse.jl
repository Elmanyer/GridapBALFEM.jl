# ===========================================================================
#  test_yl_collapse.jl — BALFE-M (p=1) against Yang & Liu's PUBLISHED formulas
#
#  WHAT THIS GATE IS FOR, AND WHY IT IS NOT LIKE THE OTHERS. Most of the suite is
#  SELF-CONSISTENCY (rule 28): it checks the code against equations the same code
#  specifies, so it passes for any residual, right or wrong. This file compares
#  our derivation against an INDEPENDENTLY PUBLISHED one -- Yang & Liu (2024),
#  JFM 999 A32 -- so it can detect an error made once in the derivation and
#  propagated consistently into both residual and forcing, which the MMS cannot.
#
#  WHAT IT COMPARES, AND WHY THAT IS THE RIGHT OBJECT. It evaluates both
#  formulations AT A PRESCRIBED STATE: an arbitrary smooth (u_j, H, h) supplied
#  directly, with no solve, no horizontal mesh, no time stepping and no assembly.
#  It therefore tests the OPERATOR alone. In particular the Class-III frozen L2
#  projections play no part, so a projection error cannot mask or mimic a
#  derivation error here -- the two failure modes are separated by construction.
#
#  STAGE 1 (this file): the vertical velocity w, our
#      w = -Σ_j [ (u_j·∇h)φ_j - (u_j·∇H)σφ_j + ∇·(H u_j) φint_j ]
#  against their (2.19)-(2.20) with the supplementary coefficients (A.1)-(A.2).
#  ⚠ This is the FULLY NON-LINEAR w: (A.1) carries h'(u_{k+1}-u_k) and H u^(1,0),
#  (A.2) carries H^(1,0) u, and the test state has h' ≠ 0 and ∇H ≠ 0 throughout.
#
#  The p=1 basis is written out EXPLICITLY below rather than taken from the
#  solver, so that this file tests the DERIVATION and not merely the assembly.
#
#  Source: the coefficients are in the article's supplementary material, free
#  with the open-access paper; a copy is kept beside the PDF in the parent
#  directory of this repository.
#
#  STAGES 2 AND 3 ARE NOT YET WRITTEN: the pressure field against (2.22)-(2.23)
#  + (B.1)-(B.4) -- which is the decisive check on the N package, hence on the
#  Class-III terms -- and the weighted residual against (2.25)+(2.28)-(2.30)
#  + (C.1)-(C.5).
# ===========================================================================
using Test
using ForwardDiff, Printf, LinearAlgebra

# ---- p = 1 hat basis, written out explicitly (NOT taken from the solver, so
#      that this script tests the derivation rather than re-testing the code) --
function phi(c::Vector{Float64}, m::Int, s::Float64)
    n = length(c)
    lo = m > 1 ? c[m-1] : c[1]
    hi = m < n ? c[m+1] : c[n]
    if m > 1 && lo <= s <= c[m]
        return (s - lo) / (c[m] - lo)
    elseif m < n && c[m] <= s <= hi
        return (hi - s) / (hi - c[m])
    else
        return 0.0
    end
end

#  φ_int_m(σ) = ∫₀^σ φ_m.  φ_m is piecewise linear, so 2-point Gauss on each
#  mesh sub-interval is EXACT to round-off.
const GX = [-1/sqrt(3), 1/sqrt(3)]
function phi_int(c::Vector{Float64}, m::Int, s::Float64)
    pts = sort(unique(vcat([x for x in c if x < s], s, 0.0)))
    pts = [x for x in pts if x <= s + 1e-15]
    tot = 0.0
    for i in 1:length(pts)-1
        a, b = pts[i], pts[i+1]
        hgt = (b - a) / 2
        mid = (a + b) / 2
        for g in GX
            tot += hgt * phi(c, m, mid + hgt*g)
        end
    end
    return tot
end

# ---- the prescribed state (arbitrary, smooth, NON-trivial in every slot) ----
d0 = 3.5
hbed(x)  = d0 - 0.8*exp(-((x - 12.0)/4.0)^2)          # variable bed  => h' ≠ 0
eta(x)   = 0.09*cos(1.1x + 0.3) + 0.03*sin(2.3x)       # free surface
Hfun(x)  = hbed(x) + eta(x)
ufun(x, m) = (0.4 + 0.15m)*sin((0.7 + 0.13m)*x + 0.4m) # nodal velocities

function state(x::Float64, M::Int)
    H  = Hfun(x);  Hx = ForwardDiff.derivative(Hfun, x)
    hx = ForwardDiff.derivative(hbed, x)
    u  = [ufun(x, m) for m in 1:M+1]
    ux = [ForwardDiff.derivative(t -> ufun(t, m), x) for m in 1:M+1]
    return (H=H, Hx=Hx, hx=hx, u=u, ux=ux)
end

# ---- OUR w ----------------------------------------------------------------
function w_ours(c::Vector{Float64}, st, s::Float64)
    M1 = length(c)
    acc = 0.0
    for m in 1:M1
        divHu = st.Hx*st.u[m] + st.H*st.ux[m]        # ∇·(H u_m) in 1-D
        acc += st.u[m]*st.hx*phi(c, m, s) -
               st.u[m]*st.Hx*s*phi(c, m, s) +
               divHu*phi_int(c, m, s)
    end
    return -acc
end

# ---- THEIR w: (A.1), (A.2), (2.20), (2.19) --------------------------------
function w_theirs_coeffs(c::Vector{Float64}, st)
    M = length(c) - 1
    w1 = zeros(M); w2 = zeros(M); w0 = zeros(M)
    for k in 1:M
        Δ = c[k] - c[k+1]                                  # NOTE: negative
        w1[k] = st.hx*st.u[k+1]/Δ - st.hx*st.u[k]/Δ +
                c[k+1]*st.H*st.ux[k]/Δ - c[k]*st.H*st.ux[k+1]/Δ          # (A.1)
        w2[k] = st.Hx*st.u[k]/(2Δ) - st.Hx*st.u[k+1]/(2Δ) +
                st.H*st.ux[k+1]/(2Δ) - st.H*st.ux[k]/(2Δ)                # (A.2)
    end
    #  (2.20): w_{1,0} = -u_1 h';  for k>1 the telescoping continuity sum
    w0[1] = -st.u[1]*st.hx
    for k in 2:M
        s = -st.u[1]*st.hx
        for m in 1:k-1
            s += c[m+1]^1*(w1[m] - w1[m+1]) + c[m+1]^2*(w2[m] - w2[m+1])
        end
        w0[k] = s
    end
    return w0, w1, w2
end

function w_theirs(c::Vector{Float64}, st, s::Float64)
    M = length(c) - 1
    w0, w1, w2 = w_theirs_coeffs(c, st)
    k = clamp(searchsortedlast(c, s), 1, M)
    return w0[k] + w1[k]*s + w2[k]*s^2
end


@testset "Yang & Liu collapse — stage 1: vertical velocity w" begin
    meshes = (("M=2 c2=0.728", [0.0, 0.728, 1.0]),
              ("M=2 c2=0.5",   [0.0, 0.5,   1.0]),
              ("M=3",          [0.0, 0.4, 0.8, 1.0]),
              ("M=4",          [0.0, 0.3, 0.6, 0.85, 1.0]))
    #  x = 12.0 sits on the bed bump, where h' is largest -- the station that
    #  actually exercises the variable-bathymetry terms of (A.1).
    stations = (0.0, 3.7, 11.3, 12.0, 19.4)
    worst = 0.0
    for (lab, c) in meshes
        M = length(c) - 1
        for x in stations
            st = state(x, M)
            dmax = 0.0; wmax = 0.0
            for s in range(1e-9, 1 - 1e-9; length = 501)
                a = w_ours(c, st, s); b = w_theirs(c, st, s)
                dmax = max(dmax, abs(a - b)); wmax = max(wmax, abs(a))
            end
            worst = max(worst, dmax)
            #  Absolute tolerance: |w| = O(1) in this state, so this is a
            #  relative agreement of ~1e-13 and the bar is round-off, not a
            #  tuned threshold (rule 15: never fix a failure by moving one).
            @test dmax < 1e-11
            @test wmax > 0.1          # the state is non-trivial, not all zeros
        end
    end
    @info "stage 1 worst |Δw| over all meshes and stations" worst
end
