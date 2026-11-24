# ==============================================================
#  periodic.jl — the closed x-periodic box (1-D stability runs)
#
#  A run of examples/local_1d/run_periodic_1d.jl is one carrier wavelength long, periodic in x,
#  one cell across in y. Every snapshot holds η(x) and the layer velocities u_j(x) (the horizontal
#  velocity at σ-node j; σ = 0 bed, σ = 1 surface). Inside the column u(x,σ) = Σ_j u_j(x) φ_j(σ).
#
#  Here: loading into [time × x] arrays, the split into box harmonics n (wavelength Lx/n), the
#  "parts" used by the plots (carrier, disturbance, steady, oscillating, any harmonic set), energies
#  per harmonic, the growth table, and frame selection for time marching.
# ==============================================================

using Serialization, TOML

"""A loaded periodic-box run.
`x` (Ns, periodic duplicate dropped, uniform), `t` (nt), `η[nt, Ns]`, `U[j, nt, Ns]` (layer j),
`basis` (the σ-basis), `D` (depth-average weights, ū = Σ D_j u_j), `A`, `B` (vertical mass matrices
∫φᵢφⱼ and ∫φᵢ_int φⱼ_int), `meta` (run parameters), `diag` (diagnostics.csv, or nothing)."""
struct PeriodicRun
    dir   :: String
    x     :: Vector{Float64}
    t     :: Vector{Float64}
    η     :: Matrix{Float64}
    U     :: Array{Float64,3}
    basis :: SigmaBasis
    D     :: Vector{Float64}
    A     :: Matrix{Float64}
    B     :: Matrix{Float64}
    meta  :: NamedTuple
    diag  :: Union{Nothing,CsvTable}
end
Base.show(io::IO, r::PeriodicRun) = @printf(io, "PeriodicRun(%s: T=%.2f s, d=%.2f m, %d cells, Nσ=%d, %d snapshots %.1f–%.1f s)",
                                             r.meta.name, r.meta.T, r.meta.d, r.meta.ncell, r.basis.Nσ, length(r.t), r.t[1], r.t[end])

const _PUBLISHED_CBDY = Dict(1 => [0.0, 1.0], 2 => [0.0, 0.728, 1.0], 3 => [0.0, 0.726, 0.925, 1.0],
                             4 => [0.0, 0.745, 0.923, 0.977, 1.0])

"Run parameters from run_manifest.toml (all values are stored as strings there)."
function _periodic_meta(dir)
    f = joinpath(dir, "run_manifest.toml")
    isfile(f) || error("load_periodic: no run_manifest.toml in $dir")
    cfg = TOML.parsefile(f)["config"]
    num(k, def = NaN) = haskey(cfg, k) ? something(tryparse(Float64, cfg[k]), def) : def
    tup(k) = Meta.parse(cfg[k])                     # "((0.0, 3.97), (0.0, 0.12))" → Expr
    dom = eval(tup("domain")); part = eval(tup("partition"))
    cb = cfg["c_bdy"] == "nothing" ? nothing : collect(Float64, eval(Meta.parse(cfg["c_bdy"])))
    (name = basename(normpath(dir)), T = num("T_wave"), d = num("h_val"), A = num("A_wave"), dt = num("dt"),
     g = num("g", 9.81), Lx = dom[1][2] - dom[1][1], Ly = dom[2][2] - dom[2][1], ncell = part[1],
     M = parse(Int, cfg["M"]), p_vert = parse(Int, get(cfg, "p_vertical", "1")), c_bdy = cb,
     p_u = parse(Int, get(cfg, "p_u", "0")), p_eta = parse(Int, get(cfg, "p_eta", "0")),
     regime = get(cfg, "regime", "?"), nl_pressure = get(cfg, "nl_pressure", "?"), flat_bed = get(cfg, "flat_bed", "?"),
     cip_gamma_u = num("cip_gamma_u", 0.0), cip_gamma_eta = num("cip_gamma_eta", 0.0), solver = get(cfg, "solver_type", "?"))
end

"""
    load_periodic(dir; c_bdy=nothing, cache=true) -> PeriodicRun

Load every `sol_t_*.vtu` of a periodic-box run (fields `eta`, `u<j>x`) into [time × x] arrays.
`c_bdy`: the σ element boundaries; default = the run's own (manifest), else the published p = 1 set
(at p ≥ 2 the solver default is a uniform split). The arrays are cached in `<dir>/_viz_fields.jls`
(delete it after the run has advanced).
"""
function load_periodic(dir::AbstractString; c_bdy = nothing, cache::Bool = true)
    meta = _periodic_meta(dir)
    cb = c_bdy !== nothing ? collect(float.(c_bdy)) : meta.c_bdy !== nothing ? meta.c_bdy :
         meta.p_vert == 1 ? _PUBLISHED_CBDY[meta.M] : collect(range(0.0, 1.0; length = meta.M + 1))
    b = sigma_basis(cb; p = meta.p_vert); N = b.Nσ
    cfile = joinpath(dir, "_viz_fields.jls")
    nfiles = count(f -> occursin(r"^sol_t_.*\.vtu$", f), readdir(dir))
    data = cache && isfile(cfile) ? deserialize(cfile) : nothing
    if data === nothing || length(data.t) != nfiles
        files = filter(f -> occursin(r"^sol_t_.*\.vtu$", f), readdir(dir))
        isempty(files) && error("load_periodic: no sol_t_*.vtu in $dir")
        ts = _time_from_name.(files); p = sortperm(ts); files = files[p]; ts = ts[p]
        P, F0 = load_snapshot(joinpath(dir, files[1]))
        haskey(F0, "u$(N)x") || error("load_periodic: no field u$(N)x: the basis (Nσ = $N) does not match the run")
        xr = round.(P[:, 1]; digits = 8); xs = sort(unique(xr))[1:end-1]       # drop x = Lx (≡ 0)
        idx = [findall(==(x), xr) for x in xs]                                 # the y-copies of each x
        η = zeros(length(ts), length(xs)); U = zeros(N, length(ts), length(xs))
        for (it, f) in enumerate(files)
            _, F = it == 1 ? (P, F0) : load_snapshot(joinpath(dir, f))
            η[it, :] = [mean(F["eta"][i]) for i in idx]
            for j in 1:N
                U[j, it, :] = [mean(F["u$(j)x"][i]) for i in idx]
            end
        end
        maximum(abs.(diff(diff(xs)))) < 1e-6 || @warn "load_periodic: x samples not uniform; the harmonic split assumes they are"
        data = (x = xs, t = ts, η = η, U = U)
        cache && serialize(cfile, data)
    end
    q(f) = _piecewise_int(b, f, 0.0, 1.0)
    D = [q(s -> phi(b, j, s)) for j in 1:N]
    A = [q(s -> phi(b, i, s) * phi(b, j, s)) for i in 1:N, j in 1:N]
    B = [q(s -> phi_int(b, i, s) * phi_int(b, j, s)) for i in 1:N, j in 1:N]
    dg = isfile(joinpath(dir, "diagnostics.csv")) ? load_csv(joinpath(dir, "diagnostics.csv")) : nothing
    return PeriodicRun(String(dir), data.x, data.t, data.η, data.U, b, D, A, B, meta, dg)
end

# ---- harmonics and parts ---------------------------------------------------------------------------

nlayers(r::PeriodicRun) = r.basis.Nσ
"Highest usable box harmonic (below the sample Nyquist)."
nharmonics(r::PeriodicRun) = length(r.x) ÷ 2 - 1
ω0(r::PeriodicRun) = 2π / r.meta.T
"Snapshots per wave period (for the period mean)."
period_samples(r::PeriodicRun) = max(1, round(Int, r.meta.T / (r.t[2] - r.t[1])))
"Indices of the snapshots with t ∈ [t0, t1]."
twindow(r::PeriodicRun, tw) = findall(t -> tw[1] ≤ t ≤ tw[2], r.t)
"Index of the snapshot closest to t."
tindex(r::PeriodicRun, t::Real) = argmin(abs.(r.t .- t))
"Depth-averaged velocity ū[nt, Ns] = Σ_j D_j u_j."
depth_average(r::PeriodicRun) = sum(r.D[j] .* r.U[j, :, :] for j in 1:nlayers(r))

"Keep only the box harmonics n ∈ ns of each row of `v` (rows = time, columns = x)."
function xband(v::AbstractMatrix, ns)
    Ns = size(v, 2); F = fft(v, 2); keep = falses(Ns)
    for n in ns
        0 ≤ n < Ns ÷ 2 || continue
        keep[n+1] = true; n > 0 && (keep[Ns-n+1] = true)
    end
    F[:, .!keep] .= 0
    return real.(ifft(F, 2))
end

"Complex amplitude of box harmonic n (coefficient of e^{+i k_n x}, k_n = 2πn/Lx) of every row."
xharmonic(v::AbstractMatrix, n::Int) = fft(v, 2)[:, n+1] ./ size(v, 2)

"Running mean over one wave period along time: what survives does not oscillate at ω0 or its multiples."
function steady_part(v::AbstractMatrix, m::Int)
    nt = size(v, 1); vs = similar(v)
    for i in 1:nt
        lo = clamp(i - m ÷ 2, 1, nt - m + 1); vs[i, :] = vec(mean(view(v, lo:lo+m-1, :); dims = 1))
    end
    return vs
end

"""
    part_rows(r, v, its, part; nmin=5, ncarrier=4) -> Matrix

The `part` of the field `v` [nt × Ns] at the snapshot indices `its`:
* `:total`; `:carrier` (n = 1…ncarrier); `:disturbance` (n ≥ nmin);
* `:steady` / `:oscillating`: the disturbance split by a running mean over one wave period;
* a vector or range of harmonics, e.g. `[16]` or `5:20`.
"""
function part_rows(r::PeriodicRun, v::AbstractMatrix, its, part; nmin::Int = 5, ncarrier::Int = 4)
    its = collect(its)
    part === :total && return v[its, :]
    part === :carrier && return xband(v[its, :], 1:ncarrier)
    part isa AbstractVector && return xband(v[its, :], part)
    part === :disturbance && return xband(v[its, :], nmin:nharmonics(r))
    if part in (:steady, :oscillating)
        m = period_samples(r); lo = max(1, minimum(its) - m); hi = min(length(r.t), maximum(its) + m)
        d = xband(v[lo:hi, :], nmin:nharmonics(r)); s = steady_part(d, m)
        k = its .- (lo - 1)
        return part === :steady ? s[k, :] : d[k, :] .- s[k, :]
    end
    error("part_rows: unknown part $part (use :total, :carrier, :disturbance, :steady, :oscillating or a harmonic list)")
end
part_label(part; nmin = 5, ncarrier = 4) =
    part === :total ? "total" : part === :carrier ? "carrier (n ≤ $ncarrier)" :
    part === :disturbance ? "disturbance (n ≥ $nmin)" : part === :steady ? "steady disturbance (period mean)" :
    part === :oscillating ? "oscillating disturbance" : "harmonics n = " * (length(part) == 1 ? string(first(part)) : "$(first(part))–$(last(part))")

"Layer fields of a part at one snapshot: (η[Ns], U[N, Ns], ū[Ns])."
function part_state(r::PeriodicRun, it::Int, part; kw...)
    N = nlayers(r)
    U = reduce(vcat, [part_rows(r, r.U[j, :, :], it:it, part; kw...) for j in 1:N])
    η = vec(part_rows(r, r.η, it:it, part; kw...))
    return η, U, vec(sum(r.D .* U; dims = 1))
end

"""u(x,σ) = Σ_j u_j(x) φ_j(σ) and the stream function ψ(x,σ) = d Σ_j u_j(x) φ_j_int(σ), on `nσ` levels.
Contours of ψ are streamlines of (u, w) with w = −∂ψ/∂x (linearised continuity, flat bed; for drawing)."""
function column_field(r::PeriodicRun, U::AbstractMatrix; nσ::Int = 100)
    σs = range(0, 1; length = nσ); N = nlayers(r)
    Φ = [phi(r.basis, j, σ) for σ in σs, j in 1:N]; Ψ = [phi_int(r.basis, j, σ) for σ in σs, j in 1:N]
    return σs, Φ * U, r.meta.d .* (Ψ * U)
end

# ---- frames ----------------------------------------------------------------------------------------

"""
    select_frames(r; trange=(t_first, t_last), dt=nothing, every=1, times=nothing, n=nothing,
                  phase=nothing) -> Vector{Int}

Snapshot indices for time marching. Precedence: `times` (explicit list, nearest snapshots) > `n`
(that many frames evenly over `trange`) > `dt` (one frame per `dt` seconds) > `every` (every k-th
snapshot). `phase ∈ [0,1)`: keep only frames at that phase of the wave period (stroboscopic view:
the carrier then looks frozen and only what does not move with it changes).
"""
function select_frames(r::PeriodicRun; trange = (r.t[1], r.t[end]), dt = nothing, every::Int = 1,
                       times = nothing, n = nothing, phase = nothing)
    times !== nothing && return unique([tindex(r, t) for t in times])
    its = twindow(r, trange)
    isempty(its) && error("select_frames: no snapshot in $trange")
    if phase !== nothing
        T = r.meta.T; t0 = r.t[its[1]]
        ts = (ceil((t0 / T) - phase) + phase) * T:T:r.t[its[end]]
        return unique([tindex(r, t) for t in ts])
    end
    n !== nothing && return unique(round.(Int, range(its[1], its[end]; length = n)))
    dt !== nothing && return unique([tindex(r, t) for t in r.t[its[1]]:dt:r.t[its[end]]])
    return its[1:every:end]
end

# ---- energy ----------------------------------------------------------------------------------------

"""
    harmonic_energy(r; its=all) -> NamedTuple of [nt × nh] matrices (energy per unit area, J/m²/ρ)

Linear energy of every box harmonic n = 1…nh, split into
* `pot`  : ½ g η²  (the surface);
* `mean` : ½ d ū²  (kinetic, depth-mean flow);
* `shear`: ½ d (uᵀA u − ū²) ≥ 0 (kinetic, vertical structure: what the depth average misses);
* `w`    : ½ d³ (∂ₓu)ᵀB(∂ₓu) (kinetic of the vertical velocity; NOT in the solver's monitor energy).
`pot + mean + shear` is the linearised (H → d) form of the monitor energy `½gη² + ½H uᵀM u`.
"""
function harmonic_energy(r::PeriodicRun; its = 1:length(r.t))
    its = collect(its); N = nlayers(r); nh = nharmonics(r); g = r.meta.g; d = r.meta.d
    Fη = fft(r.η[its, :], 2) ./ length(r.x)
    FU = [fft(r.U[j, its, :], 2) ./ length(r.x) for j in 1:N]
    pot = zeros(length(its), nh); mn = similar(pot); sh = similar(pot); ww = similar(pot)
    for n in 1:nh, i in eachindex(its)
        c = [FU[j][i, n+1] for j in 1:N]; ub = sum(r.D .* c); k = 2π * n / r.meta.Lx
        # a real field with coefficients c at ±n has mean square 2|c|²
        pot[i, n] = g * abs2(Fη[i, n+1])
        mn[i, n] = d * abs2(ub)
        sh[i, n] = d * (real(dot(c, r.A * c)) - abs2(ub))
        ww[i, n] = d^3 * k^2 * real(dot(c, r.B * c))
    end
    return (t = r.t[its], pot = pot, mean = mn, shear = sh, w = ww)
end

"""
    node_stride(r, p) -> Int

The snapshots hold `nsub = length(x)/ncell` uniform samples per cell (the VTK sub-cells), starting on a
cell edge; the Lagrange nodes of an order-`p` field sit at every `nsub ÷ p`-th sample. Errors if they
do not (nsub not a multiple of p).
"""
function node_stride(r::PeriodicRun, p::Integer)
    nsub, rm = divrem(length(r.x), r.meta.ncell)
    (rm == 0 && p ≥ 1 && nsub % p == 0 && abs(r.x[1]) < 1e-9) ||
        error("node_stride: $(length(r.x)) samples on $(r.meta.ncell) cells do not contain the order-$p nodes")
    return nsub ÷ p
end

"""
    mode_energy(r; its=all) -> (nη, pot, nu, mean, shear, w)

Linear energy of every Fourier mode of the FE solution, from its NODAL values, so every mode appears once
and no energy sits in images. (An FFT of the sub-cell samples, as in `harmonic_energy`, spreads a mode n
over n and its images ncell·q ± n (the cell-periodic shape inside each element); with p_u > 1 the images
of different modes share bins (mode 10 and the image of mode 42 on 32 cells), so they cannot be folded
back by n alone.)
* η on its ncell·p_eta nodes → modes `nη` = 0 … ncell·p_eta/2: `pot`;
* the velocities on their ncell·p_u nodes → modes `nu` = 0 … ncell·p_u/2: `mean`, `shear`, `w`.
Split as in `harmonic_energy`; rows are the snapshots `its`. n = 0 and the last mode (the node Nyquist)
carry half weight: a real sequence's mean square is |c₀|² + 2Σ|c_n|² + |c_N/2|².
"""
function mode_energy(r::PeriodicRun; its = 1:length(r.t))
    its = collect(its); N = nlayers(r); g = r.meta.g; d = r.meta.d
    sη = node_stride(r, r.meta.p_eta); su = node_stride(r, r.meta.p_u)
    Nη = length(r.x) ÷ sη; Nu = length(r.x) ÷ su
    Fη = fft(r.η[its, 1:sη:end], 2) ./ Nη
    FU = [fft(r.U[j, its, 1:su:end], 2) ./ Nu for j in 1:N]
    wt(n, Nn) = n == 0 || 2n == Nn ? 0.5 : 1.0
    nη = 0:Nη÷2; nu = 0:Nu÷2
    pot = [wt(n, Nη) * g * abs2(Fη[i, n+1]) for i in eachindex(its), n in nη]
    mn = zeros(length(its), length(nu)); sh = similar(mn); ww = similar(mn)
    for n in nu, i in eachindex(its)
        c = [FU[j][i, n+1] for j in 1:N]; ub = sum(r.D .* c); k = 2π * n / r.meta.Lx; f = wt(n, Nu)
        mn[i, n+1] = f * d * abs2(ub)
        sh[i, n+1] = f * d * (real(dot(c, r.A * c)) - abs2(ub))
        ww[i, n+1] = f * d^3 * k^2 * real(dot(c, r.B * c))
    end
    return (nη = nη, pot = pot, nu = nu, mean = mn, shear = sh, w = ww)
end

"""
    state_spectrum(r, it; avg=:period) -> (nη, pot, nu, mean, shear, w, its)

`mode_energy` at snapshot `it`: the surface energy per η mode (`nη`, `pot`) and the kinetic energy per
velocity mode (`nu`; `mean`, `shear`, `w`), n = 0 (the x-uniform flow) included. `avg = :period`
averages over one wave period centred on `it`: a free wave swaps its energy between η and u twice per
period, so the split of a single snapshot flickers while each mode's total does not. `avg = :instant`
uses the snapshot alone.
"""
function state_spectrum(r::PeriodicRun, it::Integer; avg::Symbol = :period)
    h = avg === :period ? period_samples(r) ÷ 2 : 0
    its = max(1, it - h):min(length(r.t), it + h)
    E = mode_energy(r; its); m(a) = vec(mean(a; dims = 1))
    return (nη = E.nη, pot = m(E.pot), nu = E.nu, mean = m(E.mean), shear = m(E.shear), w = m(E.w), its = its)
end

"""
    energy_history(r; modes=:groups, avg=:period, nmin=5, ncarrier=4) -> (t, labels, E, n, E0)

Total linear energy (surface + depth-mean + shear + w, from the nodal modes of `mode_energy`) over time,
`E[nt × nlines]`:
* `modes = :all`: every mode n = 0 … ncell·p_u/2, one column each (`n` lists them); the modes above
  ncell·p_eta/2 have no surface part (η has no such mode);
* `modes = :groups`: total, carrier (n = 1), bound harmonics (2 … ncarrier), disturbance (n ≥ nmin) and the
  x-uniform flow (n = 0);
* `modes = [14, 16, 32]` (any list): those modes.
`avg = :period` takes a running mean over one wave period (a free wave swaps its energy between η and u
twice per period); `:instant` keeps the snapshots. `E0` is the total at the first snapshot.
"""
function energy_history(r::PeriodicRun; modes = :groups, avg::Symbol = :period, nmin::Int = 5, ncarrier::Int = 4)
    E = mode_energy(r); nmax = last(E.nu)
    T = E.mean .+ E.shear .+ E.w
    T[:, 1:size(E.pot, 2)] .+= E.pot                        # η modes stop at ncell·p_eta/2
    E0 = sum(T[1, :])
    col(ns) = vec(sum(T[:, collect(ns) .+ 1]; dims = 2))
    if modes === :all
        lab = ["n = $n" for n in E.nu]; M = T; ns = collect(E.nu)
    elseif modes === :groups
        grp = ["total" => 0:nmax, "carrier (n = 1)" => 1:1, "bound harmonics (n = 2–$ncarrier)" => 2:ncarrier,
               "disturbance (n ≥ $nmin)" => nmin:nmax, "n = 0 (x-uniform flow)" => 0:0]
        lab = first.(grp); M = reduce(hcat, [col(last(g)) for g in grp]); ns = Int[]
    else
        ns = collect(Int, modes)
        all(n -> 0 ≤ n ≤ nmax, ns) || error("energy_history: modes must lie in 0 … $nmax (the velocity node limit)")
        lab = ["n = $n" for n in ns]; M = reduce(hcat, [col(n:n) for n in ns])
    end
    avg === :period && (M = steady_part(M, period_samples(r)))
    return (t = r.t, labels = lab, E = M, n = ns, E0 = E0)
end

"""Default bands of box harmonics, the same as `examples/periodic_growth.jl`: carrier (1), bound
harmonics (2–4), mid band (5 … nyq/2), HIGH band (nyq/2+1 … nyq) and sub-element (> nyq), with
nyq = ncell·p_u/2 the Nyquist harmonic of the velocity nodes."""
function default_bands(r::PeriodicRun; ncarrier::Int = 4)
    nh = nharmonics(r); nyq = r.meta.ncell * max(r.meta.p_u, 1) ÷ 2
    bands = [("carrier (n = 1)", 1:1), ("bound harmonics (2–$ncarrier)", 2:ncarrier)]
    rngs = [("mid band", ncarrier+1:nyq÷2), ("HIGH band (node spectrum)", nyq÷2+1:nyq), ("sub-element (> node Nyquist)", nyq+1:nh)]
    for (nm, rg) in rngs
        rg = first(rg):min(last(rg), nh)
        isempty(rg) || push!(bands, (@sprintf("%s (%d–%d)", nm, first(rg), last(rg)), rg))
    end
    return bands
end

# ---- growth table ----------------------------------------------------------------------------------

"""
    growing_harmonics(r; tw, nmin=2, top=12, layer=Nσ, io=stdout) -> Vector{NamedTuple}

Fit every harmonic of layer `layer` over the window `tw`: growth rate σ of |û| (1/s), oscillation
frequency ω/ω0 (phase rate; > 0 travels with the carrier; ⚠ folded into ±(π/Δt_snap)), the ratio
|u_layer|/|ū| (vertical structure: large = shear with no net flow) and η√(g/d)/|u_layer| (≈ 0 =
invisible at the surface). Prints the `top` first; returns all, sorted by σ (`sortby = :σ`) or by the final amplitude (`:amp`).
"""
function growing_harmonics(r::PeriodicRun; tw = (max(r.t[1], r.t[end] - 80), r.t[end]), nmin::Int = 2,
                           top::Int = 12, layer::Int = nlayers(r), io = stdout, sortby::Symbol = :σ)
    its = twindow(r, tw); t = r.t[its]; ub = depth_average(r)
    unwrap(p) = (q = copy(p); for i in 2:length(q); q[i] -= 2π * round((q[i] - q[i-1]) / 2π); end; q)
    rows = map(nmin:nharmonics(r)) do n
        c = xharmonic(r.U[layer, its, :], n); cb = xharmonic(ub[its, :], n); ce = xharmonic(r.η[its, :], n)
        (n = n, σ = cov(t, log.(abs.(c))) / var(t), ω = -cov(t, unwrap(angle.(c))) / var(t) / ω0(r),
         shear_ratio = mean(abs.(c)) / mean(abs.(cb)), surface_ratio = sqrt(r.meta.g / r.meta.d) * mean(abs.(ce)) / mean(abs.(c)),
         amp = abs(c[end]))
    end
    sort!(rows; by = x -> sortby === :amp ? -x.amp : -x.σ)
    if io !== nothing && top > 0
        @printf(io, "Fastest-growing harmonics of u%d, t ∈ [%.1f, %.1f] s:\n", layer, t[1], t[end])
        @printf(io, "  %4s  %10s  %8s  %12s  %12s  %10s\n", "n", "σ (1/s)", "ω/ω0", "|u|/|ū|", "η√(g/d)/|u|", "|û| at end")
        for x in rows[1:min(top, end)]
            @printf(io, "  %4d  %+10.4f  %+8.3f  %12.1f  %12.2f  %10.1e\n", x.n, x.σ, x.ω, x.shear_ratio, x.surface_ratio, x.amp)
        end
    end
    return rows
end

"""End of the clean (linear) growth: the last time at which the disturbance (n ≥ nmin) in the surface
layer is below `frac` × the carrier there. After it the disturbance changes the wave."""
function end_of_linear_growth(r::PeriodicRun; nmin::Int = 5, frac = 0.02)
    N = nlayers(r); v = r.U[N, :, :]
    dist = vec(sqrt.(mean(xband(v, nmin:nharmonics(r)) .^ 2; dims = 2)))
    car = vec(sqrt.(mean(xband(v, 1:1) .^ 2; dims = 2)))
    return r.t[something(findlast(dist .< frac .* car), length(r.t))]
end

"""The decay triad of the growing content: steady member = the LARGEST (at the end of the window)
growing harmonic n ≥ nmin with |ω| < 0.15 ω0 (the largest, not the fastest: a fast-growing harmonic
still at round-off level is not what the run is doing); pump = the carrier's 2nd harmonic, so the
fast member is n_s − 2 (travelling against the carrier at ≈ 2ω0). `nothing` if no steady harmonic
grows. ⚠ A guess from the table: confirm the pair (equal σ, ω ≈ 0 and ≈ 2ω0) before relying on it."""
function pick_triad(rows; nmin::Int = 5)
    c = filter(x -> x.n ≥ nmin && abs(x.ω) < 0.15 && x.σ > 0, rows)
    isempty(c) && return nothing
    x = c[argmax([x.amp for x in c])]
    return (x.n, x.n - 2)
end

"Linear dispersion of the layered model on a flat bed: ω(k) = √(g d k² Dᵀ(A + k²d²B)⁻¹D), and ω∞."
function model_dispersion(r::PeriodicRun)
    g, d = r.meta.g, r.meta.d
    ωk(k) = sqrt(g * d * k^2 * dot(r.D, (r.A + k^2 * d^2 * r.B) \ r.D))
    return ωk, sqrt(g / d * dot(r.D, r.B \ r.D))
end
