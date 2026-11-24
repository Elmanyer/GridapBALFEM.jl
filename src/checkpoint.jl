# ==============================================================
#  checkpoint.jl — write the solution state to disk and resume a run from it
#
#  A checkpoint is u_n (the free DOF values of the stacked [η, 𝖴x, 𝖴y] FEFunction), the time,
#  the step counter, the t = 0 mass/energy baselines of the diagnostics, and a FINGERPRINT of
#  everything that fixes the discrete problem. Both integrators in use (Crank–Nicolson,
#  SDIRK) are one-step methods, so u_n is the whole state.
#
#  ⚠ NOT BIT-EXACT (author's choice, 2026-10-31). Gridap starts each Newton solve from an
#  internal buffer that a checkpoint does not store: ThetaMethod's stage unknown starts from
#  the swap buffer holding u_{n−1} (ODESolutions.jl, ThetaMethod.jl:82); a DIRK stage from that
#  stage's slope at the previous step (RungeKuttaDIM.jl:82,128). Newton stops at nl_tol, so a
#  resumed run is a valid continuation that differs from the uninterrupted one at the Newton-
#  tolerance level (rule 24). test/test_restart.jl measures that difference.
#
#  FORMAT (no new dependency, readable by anything). In <output_dir>/checkpoints/:
#     ck_<step>.toml             header: `key = "value"` lines (valid TOML), [fingerprint]
#     ck_<step>.f64              sequential: the free values, raw little-endian Float64
#     ck_<step>.part<k>.f64      distributed: part k's OWN free values, in own order
#  The header is written LAST, so a checkpoint without its header is incomplete and ignored.
#  The header also stores Σv and Σv² of the values (%.17g), checked on load (relative 1e-12).
#
#  DISTRIBUTED. Gridap's distributed DOF numbering depends on the partition, so a distributed
#  checkpoint resumes only on the same cpu_grid (and mesh); anything else is refused.
# ==============================================================

#  The configuration keys that fix the discrete problem: a mismatch is refused.
const CHECKPOINT_HARD_KEYS = (:M, :p_vertical, :c_bdy, :domain, :partition, :cpu_grid, :p_u, :p_eta,
    :quad_extra, :h_val, :g, :regime, :nl_pressure, :flat_bed, :x_periodic, :y_wall_bc, :x_wall_bc,
    :cip_gamma_u, :cip_gamma_eta, :cip_hexp, :cip_order, :stabilization, :use_ad,
    :T_wave, :A_wave, :x_wm, :y_wm, :sponge_wL, :sponge_wR, :sponge_wB, :sponge_wT, :mu_max,
    :wave_gen, :wave_dir, :bc_side, :bc_profile, :relax_bc, :relax_width, :T_ramp)
#  Keys that change the time discretisation, not the problem: refused unless
#  restart_allow_change = true (the resumed run is then a different discretisation in time).
const CHECKPOINT_SOFT_KEYS = (:dt, :solver_type, :tableau, :theta, :rho_inf, :nl_tol, :nl_iter,
    :ls_rtol, :ls_maxiter, :krylov_m, :precond)

_ck_q(x) = repr(string(x))

#  A barrier only when MPI is live: the debug backend runs every part in one process.
_ck_barrier() = MPI.Initialized() ? MPI.Barrier(MPI.COMM_WORLD) : nothing

#  (count, Σv, Σv²) over the parts, each reduced separately (scalar reductions only).
function _ck_reduce_sums(parts)
    n  = reduce(+, map(p -> p[1], parts); init = 0)
    s1 = reduce(+, map(p -> p[2], parts); init = 0.0)
    s2 = reduce(+, map(p -> p[3], parts); init = 0.0)
    return n, s1, s2
end

"""
    checkpoint_fingerprint(config; n_free, bathy) → Vector{Pair{String,String}}

The fingerprint of a run: `repr` of every `CHECKPOINT_HARD_KEYS`/`CHECKPOINT_SOFT_KEYS` entry
present in `config` (the `Base.@locals()` keyword set of the driver), plus the free DOF count and
a sample of the bathymetry (`bathy`, a string). Functions and forcing objects (`eta0_func`,
`ux0_func`, `wave_bc`, `h_bathy`) are not compared by value: the bathymetry is compared through
its sample, and the initial-condition functions are irrelevant to a resumed run.
"""
function checkpoint_fingerprint(config::AbstractDict; n_free::Integer, bathy::AbstractString)
    fp = Pair{String,String}[]
    for k in (CHECKPOINT_HARD_KEYS..., CHECKPOINT_SOFT_KEYS...)
        #  ρ∞ (added 2026-11-11) belongs to the fingerprint of a generalised-α run only, so a
        #  checkpoint written before the key existed still resumes under the other integrators.
        k === :rho_inf && get(config, :solver_type, nothing) !== :gen_alpha && continue
        haskey(config, k) && push!(fp, string(k) => repr(config[k]))
    end
    push!(fp, "n_free" => string(n_free))
    push!(fp, "bathymetry" => bathy)
    return fp
end

"Sample `dfn` (x → depth) on a 7 × 3 grid over `domain` ((x0,x1),(y0,y1) or (x0,x1,y0,y1))."
function bathymetry_sample(dfn, domain)
    x0, x1, y0, y1 = length(domain) == 4 ? domain : (domain[1]..., domain[2]...)
    vals = [dfn(VectorValue(x, y)) for x in range(x0, x1; length = 7), y in range(y0, y1; length = 3)]
    return join((@sprintf("%.12g", v) for v in vals), ",")
end

_ck_dir(output_dir) = joinpath(output_dir, "checkpoints")
_ck_base(output_dir, step) = joinpath(_ck_dir(output_dir), @sprintf("ck_%09d", step))

function _ck_write_f64(path, v)
    open(path, "w") do io
        for x in v
            write(io, htol(Float64(x)))
        end
    end
end

function _ck_read_f64(path)
    n = filesize(path) ÷ 8
    v = Vector{Float64}(undef, n)
    open(path, "r") do io
        read!(io, v)
    end
    return ltoh.(v)
end

function _ck_write_header(base, step, t, nvals, s1, s2, nparts, fingerprint, mass0, energy0)
    commit = try
        strip(read(Cmd(["git", "-C", dirname(@__DIR__), "rev-parse", "HEAD"]), String))
    catch
        "unknown"
    end
    tmp = base * ".toml.tmp"
    open(tmp, "w") do io
        println(io, "# GridapBALFEM checkpoint — u_n only (not bit-exact on resume: src/checkpoint.jl)")
        println(io, "step     = ", step)
        @printf(io, "t        = %.17g\n", t)
        println(io, "n_values = ", nvals)
        println(io, "n_parts  = ", nparts)
        @printf(io, "sum      = %.17g\n", s1)
        @printf(io, "sum2     = %.17g\n", s2)
        @printf(io, "mass0    = %.17g\n", mass0)
        @printf(io, "energy0  = %.17g\n", energy0)
        println(io, "commit   = ", _ck_q(commit))
        println(io, "written  = ", _ck_q(Libc.strftime("%Y-%m-%dT%H:%M:%S", time())))
        println(io, "\n[fingerprint]")
        for (k, v) in fingerprint
            println(io, k, " = ", _ck_q(v))
        end
    end
    mv(tmp, base * ".toml"; force = true)       # the header lands last and atomically
end

"""
    write_checkpoint(output_dir, step, t, uh; fingerprint, rundiag=nothing, keep=2, ranks=nothing)

Write the checkpoint of `uh` (sequential `FEFunction` or distributed) at (`step`, `t`) into
`<output_dir>/checkpoints/`, then keep only the `keep` most recent checkpoints (`keep ≤ 0` keeps
all). The t = 0 diagnostics baselines are taken from `rundiag` so a resumed run keeps
reporting mass and energy relative to the ORIGINAL initial state.
"""
function write_checkpoint(output_dir::AbstractString, step::Integer, t::Real, uh;
                          fingerprint, rundiag = nothing, keep::Integer = 2, ranks = nothing)
    base  = _ck_base(output_dir, step)
    m0    = rundiag === nothing ? NaN : rundiag.mass0
    e0    = rundiag === nothing ? NaN : rundiag.energy0
    x     = get_free_dof_values(uh)
    if ranks === nothing
        mkpath(_ck_dir(output_dir))
        _ck_write_f64(base * ".f64", x)
        _ck_write_header(base, step, t, length(x), sum(x), sum(abs2, x), 1, fingerprint, m0, e0)
        _ck_prune(output_dir, keep)
    else
        i_am_main(ranks) && mkpath(_ck_dir(output_dir))
        _ck_barrier()
        parts = map(ranks, own_values(x)) do r, xo
            _ck_write_f64(@sprintf("%s.part%d.f64", base, r), xo)
            (length(xo), sum(xo), sum(abs2, xo))
        end
        tot = _ck_reduce_sums(parts)
        _ck_barrier()
        if i_am_main(ranks)
            _ck_write_header(base, step, t, tot[1], tot[2], tot[3], length(ranks), fingerprint, m0, e0)
            _ck_prune(output_dir, keep)
        end
        _ck_barrier()
    end
    return base * ".toml"
end

function _ck_prune(output_dir, keep)
    keep > 0 || return
    hdrs = sort(filter(f -> occursin(r"^ck_\d+\.toml$", f), readdir(_ck_dir(output_dir))))
    for h in hdrs[1:end-min(keep, length(hdrs))]
        stem = h[1:end-5]
        for f in readdir(_ck_dir(output_dir))
            startswith(f, stem * ".") && rm(joinpath(_ck_dir(output_dir), f))
        end
    end
end

"""
    read_checkpoint_header(path) → Dict{String,String}

`path` is a checkpoint header (`…/ck_<step>.toml`), a `checkpoints/` directory, or a run
directory; a directory resolves to its most recent complete checkpoint. The returned dict
has the top-level keys, the fingerprint entries as `"fp.<key>"`, and `"base"` (the path stem).
"""
function read_checkpoint_header(path::AbstractString)
    hdr = path
    if isdir(path)
        d = isdir(_ck_dir(path)) ? _ck_dir(path) : path
        hs = sort(filter(f -> occursin(r"^ck_\d+\.toml$", f), readdir(d)))
        isempty(hs) && error("read_checkpoint_header: no complete checkpoint in $d")
        hdr = joinpath(d, hs[end])
    end
    isfile(hdr) || error("read_checkpoint_header: not a checkpoint header: $hdr")
    out = Dict{String,String}("base" => hdr[1:end-5], "header" => hdr)
    sect = ""
    for line in eachline(hdr)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        if startswith(s, "[")
            sect = strip(s, ['[', ']']); continue
        end
        k, v = strip.(split(s, "="; limit = 2))
        v = startswith(v, "\"") ? unescape_string(v[2:end-1]) : v
        out[sect == "" ? k : sect == "fingerprint" ? "fp." * k : sect * "." * k] = v
    end
    return out
end

"""
    check_checkpoint_fingerprint(hdr, fingerprint; allow_change=false, where="restart")

Compare a checkpoint header with the current run's fingerprint. A HARD mismatch (the
discrete problem) always errors; a SOFT mismatch (the time discretisation) errors unless
`allow_change`, and is then listed. Returns the list of soft changes.
"""
function check_checkpoint_fingerprint(hdr::AbstractDict, fingerprint; allow_change::Bool = false,
                                      where::AbstractString = "restart")
    soft = String[]
    hard = String[]
    for (k, v) in fingerprint
        old = get(hdr, "fp." * k, nothing)
        old === nothing && (push!(hard, "$k: absent from the checkpoint, now $v"); continue)
        old == v && continue
        msg = "$k: checkpoint $old, now $v"
        Symbol(k) in CHECKPOINT_SOFT_KEYS ? push!(soft, msg) : push!(hard, msg)
    end
    isempty(hard) ||
        error("$where: the checkpoint belongs to a different discrete problem; refusing to resume.\n  " *
              join(hard, "\n  "))
    isempty(soft) || allow_change ||
        error("$where: the time discretisation differs from the checkpoint's (pass " *
              "restart_allow_change = true to resume anyway):\n  " * join(soft, "\n  "))
    return soft
end

function _ck_check_sums(hdr, n, s1, s2, where)
    n == parse(Int, hdr["n_values"]) ||
        error("$where: checkpoint holds $(hdr["n_values"]) values, read $n")
    e1, e2 = parse(Float64, hdr["sum"]), parse(Float64, hdr["sum2"])
    #  relative 1e-12: a distributed sum may be reduced in a different order on load
    ok(a, b) = abs(a - b) <= 1e-12 * max(abs(b), 1e-300) || a == b
    (ok(s1, e1) && ok(s2, e2)) ||
        error("$where: checkpoint values are corrupt (Σv $s1 ≠ $e1 or Σv² $s2 ≠ $e2)")
end

"""
    load_checkpoint_state(hdr, U; ranks=nothing, where="restart") → FEFunction

Rebuild u_n on the (already evaluated, if transient) trial space `U` from the checkpoint
`hdr` (`read_checkpoint_header`). Distributed: `ranks` given, same partition required.
"""
function load_checkpoint_state(hdr::AbstractDict, U; ranks = nothing, where::AbstractString = "restart")
    base = hdr["base"]
    if ranks === nothing
        parse(Int, hdr["n_parts"]) == 1 ||
            error("$where: a $(hdr["n_parts"])-part (distributed) checkpoint cannot resume a sequential run")
        v = _ck_read_f64(base * ".f64")
        length(v) == num_free_dofs(U) ||
            error("$where: checkpoint has $(length(v)) values, the space $(num_free_dofs(U)) free DOFs")
        _ck_check_sums(hdr, length(v), sum(v), sum(abs2, v), where)
        return FEFunction(U, v)
    end
    parse(Int, hdr["n_parts"]) == length(ranks) ||
        error("$where: checkpoint written on $(hdr["n_parts"]) parts, this run has $(length(ranks))")
    x = zero_free_values(U)
    parts = map(ranks, own_values(x)) do r, xo
        v = _ck_read_f64(@sprintf("%s.part%d.f64", base, r))
        length(v) == length(xo) ||
            error("$where: part $r holds $(length(v)) own values, this partition $(length(xo))")
        xo .= v
        (length(v), sum(v), sum(abs2, v))
    end
    _ck_check_sums(hdr, _ck_reduce_sums(parts)..., where)
    consistent!(x) |> wait
    return FEFunction(U, x)
end
