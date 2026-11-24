# ==============================================================
#  reconstruct_3d_vtk.jl — the 3-D flow (u, v, w, p) of a 2-D run, for ParaView
#
#  A run stores only η and the horizontal velocity modes u{j}x, u{j}y (w and p are eliminated
#  from the model). This rebuilds, at every snapshot, the full field at the model's own nodes: the
#  VTK vertices × the Nσ vertical σ-nodes (P1LFE-2: σ = 0, 0.728, 1), on a curvilinear structured
#  grid with z = −d(x,y) + σ (d + η), and writes one .vts per
#  snapshot plus a recon3d.pvd collection. Point fields:
#    velocity  (u, v, w)        u at σ-node j is the mode uⱼ (and v); w from the model's kinematics
#    w                          w(σ) = −(φ(σ)·𝖺) + (σφ(σ)·𝖻) − (φ_int(σ)·𝖲)        [any bed]
#                                 𝖺ⱼ = uⱼ·∇d, 𝖻ⱼ = uⱼ·∇H, 𝖲ⱼ = H ∇·uⱼ + 𝖻ⱼ,  H = d + η
#    p, p_nh, p_dyn             p = ρg(1−σ)H + p_nh,  p_nh = −ρd² (π(σ)·∂ₜ∇·u),  p_dyn = ρgη + p_nh
#                                 (p_dyn = p minus the still-water hydrostatic ρg(−z))
#    eta, depth, sigma
#  The formulas are those of src/reconstruct.jl (the solver's own write_w / write_pressure) and
#  postprocessing/src/reconstruct.jl. DEFAULT: at the same Nσ σ-nodes only (no levels between them:
#  φⱼ(σₗ) = δⱼₗ, and φ_int, π enter only as the exact basis integrals of the model). With nσ ≥ 2 it
#  adds nσ uniform σ-levels, evaluated through the vertical shape functions (u piecewise linear,
#  w piecewise quadratic, p_nh piecewise cubic in σ).
#
#  ACCURACY (printed by the script):
#   * ∂x, ∂y: 4th-order central differences on the VTK vertex grid (Gridap writes Q2 fields at
#     the cell vertices only).
#   * ∂ₜ(∇·u) for p_nh: 8th-order central differences ACROSS SNAPSHOTS (the run saved every
#     save_every·dt), so p_nh is accurate only for frequencies well below the snapshot Nyquist;
#     the response table below says how much of each frequency survives. The solver's own p_s
#     uses consecutive time steps instead. w and the hydrostatic part need no time derivative.
#   * p_nh is the linear flat-bed expression (as in the solver): over a sloping bed it omits the
#     ∇d terms of the non-hydrostatic pressure.
#   * Check printed at the end: w at the surface against ∂ₜη (linear kinematic condition).
#
#  RUN (postprocessing environment):
#    julia --project=postprocessing postprocessing/examples/reconstruct_3d_vtk.jl <run dir> [nσ] [every]
#    nσ = 0 (default): the σ-nodes only, into <run dir>/recon3d/;
#    nσ ≥ 2: the σ-nodes plus nσ uniform levels, into <run dir>/recon3d_s<nσ>/
#  Bathymetry: a flat bed of depth BALFEM_D (3.5), or the linear-production bar when the run dir
#  name contains "_bar_" (d − ½ HBAR (tanh((x−XBAR+WBAR)/S) − tanh((x−XBAR−WBAR)/S)), S = WBAR/3;
#  HBAR, XBAR, WBAR read from run_manifest.toml). Assumes a zero initial state (prepended at t = 0).
# ==============================================================
include(joinpath(@__DIR__, "..", "src", "GridapBALFEMPost.jl"))
using .GridapBALFEMPost
using WriteVTK, Printf

run_dir = abspath(length(ARGS) ≥ 1 ? ARGS[1] : error("usage: reconstruct_3d_vtk.jl <run dir> [nσ] [every]"))
nσ_uni  = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 0            # 0: σ-nodes only; ≥ 2: + uniform levels
every   = length(ARGS) ≥ 3 ? parse(Int, ARGS[3]) : 1            # write every n-th snapshot
out_dir = joinpath(run_dir, nσ_uni ≥ 2 ? "recon3d_s$(nσ_uni)" : "recon3d"); mkpath(out_dir)   # σ-nodes: recon3d/
ρ, g    = 1025.0, 9.81
c_bdy   = [0.0, 0.728, 1.0]                                      # P1LFE-2 (resolve_cbdy, M = 2)

# ---- run parameters from the manifest ----------------------------------------
man = read(joinpath(run_dir, "run_manifest.toml"), String)
mval(key, default) = (m = match(Regex("\\b$key\\s*=\\s*\"([^\"]*)\""), man);
                      m === nothing ? default : parse(Float64, m.captures[1]))
d0 = mval("h_val", 3.5)
if occursin("_bar_", basename(run_dir))
    hbar, xbar, wbar = mval("BALFEM_HBAR", 1.5), mval("BALFEM_XBAR", 26.0), mval("BALFEM_WBAR", 6.0)
    S = wbar / 3
    depth  = x -> d0 - 0.5hbar * (tanh((x - (xbar - wbar)) / S) - tanh((x - (xbar + wbar)) / S))
    ddepth = x -> -0.5hbar / S * (sech((x - (xbar - wbar)) / S)^2 - sech((x - (xbar + wbar)) / S)^2)
    @printf("bed: bar, d = %.2f, HBAR %.2f, XBAR %.1f, WBAR %.1f (crest depth %.2f)\n",
            d0, hbar, xbar, wbar, depth(xbar))
else
    depth = x -> d0; ddepth = x -> 0.0
    @printf("bed: flat, d = %.2f\n", d0)
end

# ---- snapshots: .pvd → .pvtu → pieces ---------------------------------------
pvd = joinpath(run_dir, "solution.pvd")
txt = read(pvd, String)
sets = [(parse(Float64, m[1]), joinpath(run_dir, m[2]))
        for m in eachmatch(r"timestep=\"([^\"]+)\"[^>]*file=\"([^\"]+)\"", txt)]
sort!(sets; by = first)
times = first.(sets); nt = length(times)
Δt = times[2] - times[1]
maximum(abs, diff(times) .- Δt) < 1e-9 || error("snapshots are not equally spaced")

pieces(pvtu) = [joinpath(dirname(pvtu), m[1]) for m in eachmatch(r"<Piece Source=\"([^\"]+)\"", read(pvtu, String))]
function load_pvtu(pvtu)
    pts = Matrix{Float64}[]; flds = Dict{String,Vector{Float64}}[]
    for f in pieces(pvtu)
        p, fl = load_snapshot(f); push!(pts, p); push!(flds, fl)
    end
    P = reduce(vcat, pts)
    F = Dict(k => reduce(vcat, [fl[k] for fl in flds]) for k in keys(flds[1]))
    return P, F
end

P, F1 = load_pvtu(sets[1][2])
Nσ = count(k -> occursin(r"^u\d+x$", k), keys(F1))
Nσ == length(c_bdy) || error("run has $Nσ modes, the σ-mesh here has $(length(c_bdy)) nodes")
fieldlist = ["eta"; ["u$(j)$(c)" for j in 1:Nσ for c in ("x", "y")]]

# vertex grid: each (x, y) node ↦ one VTK point (duplicates carry the same CG value)
tol = 1e-6 * maximum(abs, P)
uniq(v) = (s = sort(v); o = [s[1]]; for x in s; x - o[end] > tol && push!(o, x); end; o)
xs, ys = uniq(P[:, 1]), uniq(P[:, 2]); Nx, Ny = length(xs), length(ys)
idx = zeros(Int, Nx, Ny)
for p in axes(P, 1)
    i = searchsortedlast(xs, P[p, 1] + tol); j = searchsortedlast(ys, P[p, 2] + tol)
    idx[i, j] = p
end
all(!=(0), idx) || error("the VTK points are not a full Cartesian grid")
hx, hy = xs[2] - xs[1], ys[2] - ys[1]
@printf("grid %d × %d (dx %.4f, dy %.4f), %d modes, %d snapshots every %.3f s (t = %.1f … %.1f)\n",
        Nx, Ny, hx, hy, Nσ, nt, Δt, times[1], times[end])

# all snapshots in memory: G[name][i, j, n], n = 1 is the zero initial state at t = 0
G = Dict(k => zeros(Nx, Ny, nt + 1) for k in fieldlist)
for (n, (_, f)) in enumerate(sets)
    Pn, Fn = n == 1 ? (P, F1) : load_pvtu(f)
    n == 2 && (Pn == P || error("snapshot 2 has a different point layout"))
    for k in fieldlist
        G[k][:, :, n + 1] .= Fn[k][idx]
    end
    n % 20 == 0 && println("  read $n / $nt")
end
T = [0.0; times]

# ---- derivatives ------------------------------------------------------------
"∂/∂(dim) on the uniform grid: 4th-order central inside, 2nd order next to and at the edges."
function dgrid(f::AbstractMatrix, h, dim)
    A = dim == 1 ? f : permutedims(f)
    D = similar(A); N = size(A, 1)
    @views begin
        D[3:N-2, :] .= (A[1:N-4, :] .- 8 .* A[2:N-3, :] .+ 8 .* A[4:N-1, :] .- A[5:N, :]) ./ (12h)
        D[2, :]   .= (A[3, :] .- A[1, :]) ./ (2h);     D[N-1, :] .= (A[N, :] .- A[N-2, :]) ./ (2h)
        D[1, :]   .= (-3 .* A[1, :] .+ 4 .* A[2, :] .- A[3, :]) ./ (2h)
        D[N, :]   .= (3 .* A[N, :] .- 4 .* A[N-1, :] .+ A[N-2, :]) ./ (2h)
    end
    return dim == 1 ? D : permutedims(D)
end

# central first-derivative weights of half-width m (orders 2m)
const CW = Dict(1 => [1/2], 2 => [2/3, -1/12], 3 => [3/4, -3/20, 1/60],
                4 => [4/5, -1/5, 4/105, -1/280])
"∂ₜ of a [Nx, Ny, nt+1] series at index n: the widest central stencil that fits, one-sided at the end."
function dtime(A, n, Δt)
    N = size(A, 3); m = min(4, n - 1, N - n)
    m == 0 && return (3 .* A[:, :, n] .- 4 .* A[:, :, n-1] .+ A[:, :, n-2]) ./ (2Δt)
    D = zeros(size(A, 1), size(A, 2))
    for (k, c) in enumerate(CW[m])
        D .+= c .* (A[:, :, n+k] .- A[:, :, n-k])
    end
    return D ./ Δt
end

# response of the 8th-order stencil: the fraction of ∂ₜ that survives at frequency f
resp(f) = (x = 2π * f * Δt; sum(2c * sin(k * x) for (k, c) in enumerate(CW[4])) / x)
println("∂ₜ across snapshots (Δt = $(Δt) s, Nyquist $(round(1/(2Δt); digits=3)) Hz): fraction kept")
for f in (0.25, 0.5, 0.75, 1.0)
    @printf("    f = %.2f Hz (T = %.2f s): %.3f\n", f, 1/f, resp(f))
end

# divergence of every mode at every snapshot (needed for ∂ₜ∇·u)
DIV = [zeros(Nx, Ny, nt + 1) for _ in 1:Nσ]
for n in 1:nt+1, j in 1:Nσ
    DIV[j][:, :, n] .= dgrid(G["u$(j)x"][:, :, n], hx, 1) .+ dgrid(G["u$(j)y"][:, :, n], hy, 2)
end

# ---- σ-levels and the vertical basis -----------------------------------------
basis = sigma_basis(c_bdy; p = 1)
σnodes = sort(unique(round.(reduce(vcat, basis.local_nodes); digits = 12)))
length(σnodes) == Nσ || error("$(length(σnodes)) σ-nodes for $Nσ modes")
σs  = nσ_uni ≥ 2 ? sort(unique(round.([collect(range(0, 1; length = nσ_uni)); σnodes]; digits = 12))) : σnodes
Nz  = length(σs)
PHI = [phi(basis, j, σ) for σ in σs, j in 1:Nσ]
PHI_INT = [phi_int(basis, j, σ) for σ in σs, j in 1:Nσ]
PI3 = [pi3(basis, j, σ) for σ in σs, j in 1:Nσ]
println("σ-levels: ", join(round.(σs; digits = 3), ", "))

dX = [depth(x) for x in xs, _ in ys]; dXx = [ddepth(x) for x in xs, _ in ys]   # ∂y d = 0

# ---- per snapshot: reconstruct and write -------------------------------------
X = Array{Float32}(undef, Nx, Ny, Nz); Y = similar(X); Z = similar(X)
for k in 1:Nz, j in 1:Ny
    X[:, j, k] .= xs; Y[:, j, k] .= ys[j]
end
U, V, W, Pt, Pnh, Pdyn, E = (similar(X) for _ in 1:7)
D3 = similar(X); S3 = similar(X)
for k in 1:Nz
    D3[:, :, k] .= dX; S3[:, :, k] .= σs[k]
end
check_num = 0.0; check_den = 0.0
xin = 0.15 * (xs[end] - xs[1]) .< xs .- xs[1] .< 0.6 * (xs[end] - xs[1])   # away from relaxation/sponge

paraview_collection(joinpath(out_dir, "recon3d")) do pvd
    for n in 2:nt+1
        (n - 2) % every == 0 || continue
        η = G["eta"][:, :, n]
        ηx, ηy = dgrid(η, hx, 1), dgrid(η, hy, 2)
        H = dX .+ η
        ux = [G["u$(j)x"][:, :, n] for j in 1:Nσ]; uy = [G["u$(j)y"][:, :, n] for j in 1:Nσ]
        𝖺 = [ux[j] .* dXx for j in 1:Nσ]                                   # u·∇d
        𝖻 = [ux[j] .* (dXx .+ ηx) .+ uy[j] .* ηy for j in 1:Nσ]             # u·∇H
        𝖲 = [H .* DIV[j][:, :, n] .+ 𝖻[j] for j in 1:Nσ]                   # ∇·(H u)
        dDIV = [dtime(DIV[j], n, Δt) for j in 1:Nσ]
        for k in 1:Nz
            σ = σs[k]
            u = sum(PHI[k, j] .* ux[j] for j in 1:Nσ); v = sum(PHI[k, j] .* uy[j] for j in 1:Nσ)
            w = sum(-PHI[k, j] .* 𝖺[j] .+ σ * PHI[k, j] .* 𝖻[j] .- PHI_INT[k, j] .* 𝖲[j] for j in 1:Nσ)
            pnh = -ρ .* dX .^ 2 .* sum(PI3[k, j] .* dDIV[j] for j in 1:Nσ)
            Z[:, :, k] .= -dX .+ σ .* H
            U[:, :, k] .= u; V[:, :, k] .= v; W[:, :, k] .= w
            Pnh[:, :, k] .= pnh; Pt[:, :, k] .= ρ * g * (1 - σ) .* H .+ pnh
            Pdyn[:, :, k] .= ρ * g .* η .+ pnh; E[:, :, k] .= η
            if k == Nz                                                     # w(1) vs ∂ₜη
                ηt = dtime(G["eta"], n, Δt)
                global check_num += sum(abs2, (w .- ηt)[xin, :]); global check_den += sum(abs2, ηt[xin, :])
            end
        end
        tag = replace(@sprintf("%.4f", T[n]), "." => "_")
        vtk_grid(joinpath(out_dir, "recon_t_$tag"), X, Y, Z; compress = true) do vtk
            vtk["velocity", VTKPointData()] = (U, V, W)
            vtk["w", VTKPointData()] = W
            vtk["p", VTKPointData()] = Pt
            vtk["p_nh", VTKPointData()] = Pnh
            vtk["p_dyn", VTKPointData()] = Pdyn
            vtk["eta", VTKPointData()] = E
            vtk["depth", VTKPointData()] = D3
            vtk["sigma", VTKPointData()] = S3
            pvd[T[n]] = vtk
        end
        (n - 1) % 20 == 0 && println("  wrote t = $(T[n])")
    end
end
@printf("check: ‖w(σ=1) − ∂ₜη‖ / ‖∂ₜη‖ = %.3e  (x in the middle 15–60 %% of the domain)\n",
        sqrt(check_num / check_den))
println("written: ", joinpath(out_dir, "recon3d.pvd"))
