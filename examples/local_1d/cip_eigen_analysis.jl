# ==============================================================
#  cip_eigen_analysis.jl — what does the hp C⁰-IP penalty do to the box spectrum?
#
#  markdown_files/BROKEN_FORMULATION_PLAN.md §5.7 step 2. The LINEARISED flat-bed operator of
#  the closed x-periodic box, −M_eff⁻¹ (K + J_h), assembled from the solver's own hand
#  Jacobians (`jacobian_u_t` → M_eff, `jacobian_u` → K + J_h; both exact in the linear regime)
#  and diagonalised densely. Each OSCILLATORY eigenmode is assigned to the wavenumber band of
#  `periodic_growth.jl` by the dominant harmonic of its sampled η and 𝖴x profiles, and its
#  ENERGY damping rate 2·Re(−λ) is reported per band, against γ and the penalty order.
#
#  THE WINDOW the screening needs (§5.2):
#     mid-band damping  ≳ 0.13 s⁻¹   (the late nonlinear-pressure growth rate it must beat)
#     carrier damping   ≪ 1/T_run ≈ 0.006 s⁻¹   (no locking)
#  ⚠ The linear operator has no instability: these are the dissipation budgets available
#  against the nonlinear growth, not stability verdicts.
#
#  RUN (julia-mcp or):  julia --project=. examples/local_1d/cip_eigen_analysis.jl
#  OUT: output/local_1d/cip_eigen/cip_eigen.csv
# ==============================================================
using GridapBALFEM, Gridap, Gridap.ODEs, LinearAlgebra, SparseArrays, Printf

const OUTD = joinpath(@__DIR__, "..", "..", "output", "local_1d", "cip_eigen")
mkpath(OUTD)
const VERT = assemble_vertical_tensors(2, 1, Vector{Float64}(resolve_cbdy(2, nothing, 1)))
const NS   = VERT.N_dof
const D0   = 3.5
const LX   = 4.0                      # one carrier wavelength (≈ the model λ at kd = 5.5)

function box(pu, ncell)
    pe = pu - 1
    h  = LX / ncell
    model, trian = build_horizontal_model(((0.0, LX), (0.0, h)), (ncell, 1); x_periodic = true)
    dΩ = Measure(trian, 2 * pu + 4)
    U, V = build_fe_spaces(model, pu, NS; y_wall_bc = :wall, p_eta = pe)
    return (pu = pu, pe = pe, ncell = ncell, h = h, model = model, trian = trian, dΩ = dΩ, U = U, V = V)
end

"Bands in box harmonics, as periodic_growth.jl: carrier, low, mid, high (node spectrum), sub-element."
function bands(B)
    nN = B.pu * B.ncell ÷ 2                      # node Nyquist
    return [("carrier", 1:1), ("low", 2:4), ("mid", 5:(nN ÷ 2)), ("high", (nN ÷ 2 + 1):nN),
            ("sub", (nN + 1):(B.pu * B.ncell))]
end

function operator(B, γu, γe, order; stab = :jumpgrad)
    p = build_problem(VERT; h_bathy = x -> D0, regime = :linear, nl_pressure = false, flat_bed = true)
    (γu > 0 || γe > 0) && attach_skeleton!(p, B.model; cip_gamma_u = γu, cip_gamma_eta = γe,
                                           cip_order = (stab === :ghostvolume ? 1 : order),
                                           stabilization = stab, p_u = B.pu, p_eta = B.pe,
                                           degree = 2 * B.pu + 4)
    z  = FEFunction(B.U, zeros(num_free_dofs(B.U)))
    tu = TransientCellField(z, (z,))
    K  = assemble_matrix((du, v) -> jacobian_u(0.0, tu, du, v, p, B.trian, B.dΩ), B.U, B.V)
    Me = assemble_matrix((du, v) -> jacobian_u_t(0.0, tu, du, v, p, B.trian, B.dΩ), B.U, B.V)
    return eigen(-(Matrix(Me) \ Matrix(K)))
end

"Dominant box harmonic of a mode: power spectra of η and all 𝖴x components sampled on a
uniform x-line (2p samples per cell), each normalised, summed over real and imaginary parts."
function dominant_n(B, vec)
    Ns = 2 * B.pu * B.ncell
    xs = [Point((m - 0.5) * LX / Ns, B.h / 2) for m in 1:Ns]
    P  = zeros(Ns ÷ 2 + 1)
    for part in (real(vec), imag(vec))
        maximum(abs, part) < 1e-14 && continue
        uh = FEFunction(B.U, part)
        ηs = [uh[1](x) for x in xs]
        us = [uh[2](x) for x in xs]
        for sig in (ηs, [u[j] for u in us, j in 1:NS])
            S = sig isa Vector ? reshape(sig, :, 1) : sig
            for c in 1:size(S, 2)
                F = abs2.([sum(S[m, c] * cis(-2π * n * (m - 1) / Ns) for m in 1:Ns) for n in 0:(Ns ÷ 2)])
                tot = sum(F)
                tot > 0 && (P .+= F ./ tot)
            end
        end
    end
    xfrac = 1 - P[1] / max(sum(P), 1e-300)        # share of the power that VARIES along x
    P[1] = 0.0                                    # the mean is not a wave
    return argmax(P) - 1, xfrac
end

function analyse(B, γu, γe, order; io = stdout, stab = :jumpgrad)
    ev  = operator(B, γu, γe, order; stab = stab)
    ω   = abs.(imag.(ev.values)); σ = real.(ev.values)
    #  each conjugate pair once (Im λ > 0); degenerate ±k travelling waves are KEPT as two modes
    keep = findall(imag.(ev.values) .> 1e-6)
    cls  = [dominant_n(B, ev.vectors[:, i]) for i in keep]
    #  ⚠ x-UNIFORM modes are skipped: on the one-cell-wide flume there is a TRANSVERSE standing
    #  mode (ω ≈ 7.55 rad/s, η varying in y only) that the x-facet penalty cannot see and a
    #  y-invariant problem never excites. Left in, it reads as an "undamped" mode with a random
    #  band label (BROKEN_FORMULATION_PLAN.md §5.8).
    xmode = [c[2] > 0.5 for c in cls]
    ns   = [c[1] for c in cls]
    damp = [-2 * σ[i] for i in keep]              # ENERGY damping rate (s⁻¹), comparable with σ_E
    nover = count((ω .< 1e-6) .& (σ .< -1e-8))
    #  the CARRIER is tracked as the x-wave nearest the carrier frequency, not by band label
    ic   = [k for k in eachindex(keep) if xmode[k] && ns[k] == 1]
    car  = isempty(ic) ? NaN : damp[ic[argmin(abs.(ω[keep[ic]] .- 3.913))]]
    row = Dict{String,Any}("carrier" => (length(ic), car, car))
    for (bn, r) in bands(B)
        bn == "carrier" && continue
        sel = [k for k in eachindex(keep) if xmode[k] && ns[k] in r]
        row[bn] = isempty(sel) ? (0, NaN, NaN) : (length(sel), minimum(damp[sel]), maximum(damp[sel]))
    end
    return row, length(keep), nover
end

function run_ladder(B; gammas = (0.0, 1e-4, 3e-4, 1e-3, 2e-3, 3e-3, 1e-2, 3e-2), orders = (1, 2),
                    fields = (:both, :u), csv = nothing)
    for fld in fields, order in orders, γ in gammas
        (γ == 0 && (order > 1 || fld != first(fields))) && continue
        γu = γ; γe = fld === :both ? γ : 0.0
        t0 = time()
        row, nk, nover = analyse(B, γu, γe, order)
        c = row["carrier"]; mid = row["mid"]; hi = row["high"]; sb = row["sub"]
        @printf("Q%d/Q%d n%-3d %-4s ord %d γ=%-6g | carrier %.2e | mid(%2d) min %.3e | high(%2d) min %.3e | sub(%2d) min %.3e | osc %3d overdamped %3d  [%.0fs]\n",
                B.pu, B.pe, B.ncell, fld, order, γ, c[2], mid[1], mid[2], hi[1], hi[2], sb[1], sb[2],
                nk, nover, time() - t0)
        flush(stdout)
        csv === nothing || println(csv, join((B.pu, B.ncell, fld, order, γ, c[2], mid[1], mid[2], mid[3],
                                             hi[1], hi[2], hi[3], sb[1], sb[2], sb[3], nk, nover), ","))
    end
end

if abspath(PROGRAM_FILE) == (@__FILE__) || get(ENV, "CIP_EIG_RUN", "0") == "1"
    open(joinpath(OUTD, "cip_eigen.csv"), "w") do csv
        println(csv, "pu,ncell,fields,order,gamma,carrier_damp,mid_n,mid_min,mid_max,high_n,high_min,high_max,sub_n,sub_min,sub_max,n_osc,n_overdamped")
        for (pu, nc) in ((3, 16), (2, 32))
            run_ladder(box(pu, nc); csv = csv)
            flush(csv)
        end
    end
    println("written ", joinpath(OUTD, "cip_eigen.csv"))
end
