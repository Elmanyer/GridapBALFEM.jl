# ==============================================================
#  cip_eigen_order_rule.jl — which jump orders must a penalty reach to leave no undamped
#  element-scale waves? The "weak-form order" rule against the "p − 1" rule.
#  (markdown_files/STABILISATION_ORDER_RULE.md)
#
#  THE QUESTION. A jump penalty of order ≤ m on Q_p fields leaves undamped exactly the
#  C^m piecewise polynomials (its kernel); oscillatory modes that fit there escape the
#  damping. Two candidate rules for the order that suffices:
#     (W) the order of the WEAK FORM: trial second derivatives ⇒ m = 2 for every field;
#     (P) the ELEMENT order: m = p − 1 per field, so the kernel is the maximally smooth
#         C^{p−1} splines (one wave dof per cell, no optical branches).
#  At Q3/Q2 the two coincide for the velocity (p − 1 = 2) — the case on record. Q4/Q3 separates
#  them: with m = 2 on both fields, (W) predicts every mid/high/sub-element mode damped; (P)
#  predicts undamped velocity modes left in the C² quartic splines (two wave dofs per cell).
#
#  METHOD. cip_eigen_analysis.jl: the LINEARISED flat-bed box operator −M_eff⁻¹(K + J_h),
#  one carrier wavelength, x-periodic, one cell across with walls, built from the solver's own
#  exact linear-regime Jacobians, diagonalised densely. Each oscillatory mode is labelled by its
#  dominant box harmonic (periodic_growth.jl bands); x-uniform (transverse) modes are excluded.
#  Damping = ENERGY rate 2·Re(−λ) (s⁻¹). ⚠ The linear operator has no instability: these are the
#  dissipation budgets available against the nonlinear growth, not stability verdicts.
#
#  CASES. Q4/Q3 and, as the same-batch control, Q3/Q2 — both at 16 cells per wavelength;
#  :jumpgrad orders 1 and 2, penalty on both fields and on the velocity only; γ ladder; plus the
#  ghost penalty (full order p) at Q4/Q3 as the complete-family reference.
#
#  OUTPUT (persistent; revisit without re-running):
#     output/local_1d/cip_eigen_order_rule/summary.csv        one row per case (band minima/maxima)
#     output/local_1d/cip_eigen_order_rule/modes_<case>.csv   every eigenvalue: re, im, n_dom, x_frac, band
#     output/local_1d/cip_eigen_order_rule/run_manifest.toml  commit, config, environment
#  RUN:  julia --project=. examples/local_1d/cip_eigen_order_rule.jl  (log next to the output)
# ==============================================================
include(joinpath(@__DIR__, "cip_eigen_analysis.jl"))      # box, bands, operator, dominant_n (no run)

#  EIG_OUT selects the output folder. The first pass (2026-10-06, folder root) classified modes on the
#  centre line only, which cannot see TRANSVERSE modes (η ∝ cos(πy/h) vanishes there); the second pass
#  (EIG_OUT=.../yfilter) adds the y-structure test below.
const OUT2 = get(ENV, "EIG_OUT", joinpath(@__DIR__, "..", "..", "output", "local_1d", "cip_eigen_order_rule"))
mkpath(OUT2)
LinearAlgebra.BLAS.set_num_threads(parse(Int, get(ENV, "EIG_BLAS", "2")))

const GAMMAS = (1e-4, 3e-4, 1e-3, 2e-3, 3e-3, 1e-2, 3e-2)
const NCELL  = 16

write_run_manifest(OUT2, Dict{String,Any}("script" => "cip_eigen_order_rule.jl", "pairings" => "Q4/Q3, Q3/Q2",
                          "ncell" => NCELL, "gammas" => collect(GAMMAS), "orders" => [1, 2],
                          "fields" => ["both", "u"], "ghost_reference" => "Q4/Q3",
                          "carrier_omega" => 3.913, "d" => D0, "Lx" => LX);
                   driver = "examples/local_1d/cip_eigen_order_rule.jl")

bandof(B, n) = (for (bn, r) in bands(B); n in r && return bn; end; "beyond")

"""y-structure of a mode on the one-cell-wide box: η and every 𝖴x component sampled on the lines
y = h/4 and 3h/4; `yodd` = energy share of the part ODD in y (top − bottom), `uy` = energy share of
𝖴y. A y-invariant run never excites a y-odd or 𝖴y-carried mode, and the x-facet penalties cannot
reach it (there are no interior y-facets), so such modes are classed `transverse`."""
sq(v) = sum(abs2, v.data)                          # |v|² of a Gridap VectorValue
function ystructure(B, vec)
    Ns = 2 * B.pu * B.ncell
    xs = [(m - 0.5) * LX / Ns for m in 1:Ns]
    odd = 0.0; tot = 0.0; ey = 0.0; ex = 0.0
    for part in (real(vec), imag(vec))
        maximum(abs, part) < 1e-14 && continue
        uh = FEFunction(B.U, part)
        for x in xs
            pb, pt = Point(x, B.h / 4), Point(x, 3B.h / 4)
            eb, et = uh[1](pb), uh[1](pt)
            ub, ut = uh[2](pb), uh[2](pt)
            vb, vt = uh[3](pb), uh[3](pt)
            odd += abs2(et - eb) + sq(ut - ub)
            tot += abs2(et - eb) + abs2(et + eb) + sq(ut - ub) + sq(ut + ub)
            ey  += sq(vb) + sq(vt)
            ex  += sq(ub) + sq(ut)
        end
    end
    return odd / max(tot, 1e-300), ey / max(ey + ex, 1e-300)
end

"Diagonalise one case, dump every mode, return the per-band summary."
function case!(summary, B, γu, γe, order, stab, fld)
    tag = @sprintf("Q%dQ%d_n%d_%s_%s_ord%d_g%g", B.pu, B.pe, B.ncell, stab, fld, order, γu)
    t0  = time()
    ev  = operator(B, γu, γe, order; stab = stab)
    keep = findall(imag.(ev.values) .> 1e-6)
    open(joinpath(OUT2, "modes_$tag.csv"), "w") do io
        println(io, "re,im,omega,energy_damping,n_dom,x_frac,y_odd,uy_frac,band")
        for i in keep
            n, xf = dominant_n(B, ev.vectors[:, i])
            yo, uy = ystructure(B, ev.vectors[:, i])
            λ = ev.values[i]
            band = (xf > 0.5 && yo < 0.5 && uy < 0.5) ? bandof(B, n) : "transverse"
            @printf(io, "%.10e,%.10e,%.8e,%.8e,%d,%.4f,%.4f,%.4f,%s\n", real(λ), imag(λ), abs(imag(λ)),
                    -2 * real(λ), n, xf, yo, uy, band)
        end
    end
    #  summary from the dump (same rules as analyse(): x-modes only, carrier = n=1 nearest 3.913)
    rows = [split(l, ",") for l in readlines(joinpath(OUT2, "modes_$tag.csv"))[2:end]]
    om   = [parse(Float64, r[3]) for r in rows]; dm = [parse(Float64, r[4]) for r in rows]
    bd   = [r[end] for r in rows]
    ic   = findall(==("carrier"), bd)
    car  = isempty(ic) ? NaN : dm[ic[argmin(abs.(om[ic] .- 3.913))]]
    nover = count((abs.(imag.(ev.values)) .< 1e-6) .& (real.(ev.values) .< -1e-8))
    stats(b) = (s = findall(==(b), bd); isempty(s) ? (0, NaN, NaN, NaN) :
                (length(s), minimum(dm[s]), maximum(dm[s]), om[s[argmin(dm[s])]]))
    mid, hi, sb = stats("mid"), stats("high"), stats("sub")
    @printf("%-40s | carrier %.2e | mid(%2d) min %.3e @ω=%.2f | high(%2d) min %.3e | sub(%2d) min %.3e | osc %3d over %3d [%.0fs]\n",
            tag, car, mid[1], mid[2], mid[4], hi[1], hi[2], sb[1], sb[2], length(keep), nover, time() - t0)
    println(summary, join((B.pu, B.pe, B.ncell, stab, fld, order, γu, γe, car,
                           mid..., hi..., sb..., length(keep), nover), ","))
    flush(summary); flush(stdout)
end

open(joinpath(OUT2, "summary.csv"), "w") do summary
    println(summary, "pu,pe,ncell,stab,fields,order,gamma_u,gamma_eta,carrier_damp," *
                     "mid_n,mid_min,mid_max,mid_omega_at_min,high_n,high_min,high_max,high_omega_at_min," *
                     "sub_n,sub_min,sub_max,sub_omega_at_min,n_osc,n_overdamped")
    for pu in (4, 3)                                  # Q4/Q3 first: it is the new question
        B = box(pu, NCELL)
        println("\n######## Q$(pu)/Q$(pu-1), $(NCELL) cells/λ — free DOFs $(num_free_dofs(B.U))"); flush(stdout)
        case!(summary, B, 0.0, 0.0, 1, :jumpgrad, "none")
        for fld in ("both", "u"), order in (1, 2), γ in GAMMAS
            case!(summary, B, γ, fld == "both" ? γ : 0.0, order, :jumpgrad, fld)
        end
        pu == 4 && for γ in GAMMAS
            case!(summary, B, γ, γ, 1, :ghostvolume, "both")
        end
    end
end
println("\ncip_eigen_order_rule: done → ", OUT2)
