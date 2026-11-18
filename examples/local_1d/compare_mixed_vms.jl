# ==============================================================
#  compare_mixed_vms.jl — the new solver against the v2 solver, models 1–6, short box runs
#
#  Branch mixed-vms-solver (markdown_files/MIXED_VMS_PLAN.md). Same closed x-periodic box as
#  run_periodic_1d.jl (one model wavelength, ny = 1, solid y-walls, hot start from the discrete
#  linear eigenmode, no seed), Crank–Nicolson, Q3/Q2. For each model the REFERENCE is the v2
#  primal solver; the other arms are what this branch adds:
#
#      m1, m2   primal                      (linear: no mixed form, no VMS — the code path is v2's)
#      m3, m4   primal | primal + VMS
#      m5, m6   primal (Class III cellwise + facet terms + C⁰-CIP) | mixed | mixed + VMS
#
#  Compared: η at four gauges over the run (relative L2 difference of the series, against the
#  reference arm), the energy ratio E/E₀ at the end, max|η|, Newton iterations per step.
#
#  Output: output/mixed_vms/comparison/<case>/  (run_manifest.toml, CONFIGURATION.txt,
#  diagnostics.csv, gauges.csv) and output/mixed_vms/comparison/summary.txt.
#
#  KNOBS: BALFEM_CMP_MODELS (default "1,2,3,4,5,6"), BALFEM_CMP_PERIODS (5), BALFEM_CMP_NCELL (16),
#         BALFEM_CMP_KA (0.08), BALFEM_DT (0.04), BALFEM_CMP_SLOPE (0.05, models 2/4/6),
#         BALFEM_CMP_ARMS ("m5:mixed,…": only those arms, no summary), BALFEM_CMP_MODE (run | summary).
#  Parallel launcher: run/local/run_mixed_vms_compare.sh (one process per arm, memory-gated).
# ==============================================================
using GridapBALFEM, Printf
using Gridap.TensorValues: VectorValue

const ROOT_CMP = normpath(joinpath(@__DIR__, "..", ".."))
_e(k, d)   = get(ENV, k, string(d))
models   = parse.(Int, split(_e("BALFEM_CMP_MODELS", "1,2,3,4,5,6"), ","))
periods  = parse(Float64, _e("BALFEM_CMP_PERIODS", 5))
ncell    = parse(Int, _e("BALFEM_CMP_NCELL", 16))
kA       = parse(Float64, _e("BALFEM_CMP_KA", 0.08))
dt       = parse(Float64, _e("BALFEM_DT", 0.04))
slope    = parse(Float64, _e("BALFEM_CMP_SLOPE", 0.05))
outroot  = joinpath(ROOT_CMP, "output", "mixed_vms", "comparison")
M, p_vert, p_u, p_eta, d, Twave = 2, 1, 3, 2, 3.5, 1.6

vert0 = assemble_vertical_tensors(M, p_vert, Vector{Float64}(resolve_cbdy(M, nothing, p_vert)))
Nσ    = vert0.N_dof
k0    = GridapBALFEM.model_wavenumber(vert0, 2π / Twave, d, 9.81)
Awave = kA / k0
wi    = WaveInput(vert0; A = Awave, T = Twave, d = d, T_ramp = 0.0, profile = :model)
Lx    = 2π / wi.ks[1]
dx    = Lx / ncell
inc   = incident_fields(wi)
kbed  = 2π / Lx
bedf  = x -> d + (slope / kbed) * cos(kbed * x[1])
gauges = [(f * Lx, 0.5dx) for f in (0.125, 0.375, 0.625, 0.875)]

"The arms of one model: (tag, keyword set). The first is the reference."
function arms(m::Int)
    m in (1, 2) && return [("primal", (;))]
    #  `_ad`: the same operator with the exact AD Jacobian. With VMS the hand (quasi-Newton) Jacobian
    #  does not converge at step 2 (vms_newton_probe.jl, 2026-11-17: AD 3 iterations to 1e-14).
    #  `_oss`: orthogonal subgrid scales (2026-11-17, MIXED_VMS_PLAN.md §3.5), hand (exact) Jacobian.
    m in (3, 4) && return [("primal", (;)), ("primal_vms", (; vms = true)),
                           ("primal_vms_ad", (; vms = true, use_ad = true)),
                           ("primal_vms_oss", (; vms = true, vms_oss = true))]
    return [("primal_cip", (;)), ("mixed", (; formulation = :mixed)),
            ("mixed_vms", (; formulation = :mixed, vms = true)),
            ("mixed_vms_ad", (; formulation = :mixed, vms = true, use_ad = true)),
            ("mixed_vms_oss", (; formulation = :mixed, vms = true, vms_oss = true))]
end

function run_arm(m::Int, tag::String, kw)
    regime = m in (1, 2) ? :linear : :nonlinear
    nlp    = m in (5, 6)
    flat   = isodd(m)
    out    = joinpath(outroot, @sprintf("m%d_%s", m, tag))
    t0 = time()
    diags, _, _ = setup_and_run(;
        M = M, p_vertical = p_vert, p_u = p_u, p_eta = p_eta,
        domain = ((0.0, Lx), (0.0, dx)), partition = (ncell, 1),
        x_periodic = true, y_wall_bc = :wall,
        h_val = d, T_wave = Twave, A_wave = Awave,
        sponge_wL = 0.0, sponge_wR = 0.0,
        eta0_func = x -> inc.eta(x, 0.0), ux0_func = x -> inc.ux(x, 0.0),
        T_final = periods * Twave, dt = dt, solver_type = :theta, theta = 0.5,
        regime = regime, nl_pressure = nlp, flat_bed = flat, h_bathy = flat ? nothing : bedf,
        output_dir = out, save_every = 0, gauges = gauges,
        nl_tol = 1e-8, print_every = 25, check_every = 0, diag_every = 5, eta_ref = Awave,
        kw...)
    wall = time() - t0
    open(joinpath(out, "gauges.csv"), "w") do io
        println(io, "t,", join(("eta_g$i" for i in 1:length(gauges)), ","))
        for s in diags
            println(io, @sprintf("%.6f,", s.t), join((@sprintf("%.10e", v) for v in s.gauge_vals), ","))
        end
    end
    #  energy ratio at the last diagnostics sample
    lines = readlines(joinpath(out, "diagnostics.csv"))
    hdr   = split(lines[1], ","); last = split(lines[end], ",")
    eratio = parse(Float64, last[findfirst(==("energy_ratio"), hdr)])
    G = reduce(hcat, [s.gauge_vals for s in diags])           # gauges × steps
    return (; m, tag, G, eratio, emax = maximum(s.eta_max for s in diags),
              newton = sum(s.nl_iters for s in diags) / length(diags),
              nsteps = length(diags), expected = round(Int, periods * Twave / dt), wall)
end

"Write the arm's scalar results next to its run (read back by the summary mode)."
function write_arm_result(m, tag, r; status = "ok", msg = "")
    out = joinpath(outroot, @sprintf("m%d_%s", m, tag)); mkpath(out)
    open(joinpath(out, "arm_result.txt"), "w") do io
        println(io, "status=", status)
        r === nothing || for k in (:nsteps, :expected, :eratio, :emax, :newton, :wall)
            println(io, k, "=", getfield(r, k))
        end
        isempty(msg) || println(io, "message=", replace(msg, '\n' => ' '))
    end
end

"Read an arm back: (status, scalars, gauge matrix gauges × steps), or nothing if it never ran."
function read_arm(m, tag)
    out = joinpath(outroot, @sprintf("m%d_%s", m, tag))
    f = joinpath(out, "arm_result.txt"); isfile(f) || return nothing
    kv = Dict(split(l, "=", limit = 2)[1] => split(l, "=", limit = 2)[2] for l in readlines(f) if occursin('=', l))
    G = nothing
    if isfile(joinpath(out, "gauges.csv"))
        L = readlines(joinpath(out, "gauges.csv"))[2:end]
        G = reduce(hcat, [parse.(Float64, split(l, ",")[2:end]) for l in L])
    end
    return (; kv, G)
end

#  MODES. Default: run every arm of BALFEM_CMP_MODELS in this process, then summarise.
#  BALFEM_CMP_ARMS="m5:mixed,m6:mixed_vms" runs only those arms (one process per arm when launched by
#  run/local/run_mixed_vms_compare.sh). BALFEM_CMP_MODE=summary runs nothing and rebuilds summary.txt
#  from the arm_result.txt + gauges.csv of every arm on disk.
mode    = _e("BALFEM_CMP_MODE", "run")
armsel  = _e("BALFEM_CMP_ARMS", "")
mkpath(outroot)
if mode != "summary"
    for m in models, (tag, kw) in arms(m)
        isempty(armsel) || "m$m:$tag" in split(armsel, ",") || continue
        try
            r = run_arm(m, tag, kw)
            write_arm_result(m, tag, r)
            @printf("ARM m%d %-12s done: %d/%d steps, E/E0 %.6f, Newton %.2f, wall %.1f s\n",
                    m, tag, r.nsteps, r.expected, r.eratio, r.newton, r.wall)
        catch e
            write_arm_result(m, tag, nothing; status = "failed", msg = first(sprint(showerror, e), 200))
            @warn "compare_mixed_vms: m$m $tag failed" exception = (e, catch_backtrace())
            println("ARM m$m $tag FAILED")
        end
        flush(stdout)
    end
    isempty(armsel) || exit(0)          # a single-arm process does not write the summary
end

rows = String[]
push!(rows, @sprintf("%-4s %-12s %7s %12s %12s %10s %8s %9s", "m", "arm", "steps", "gauge_diff", "E/E0", "max|eta|/A", "Newton", "wall [s]"))
for m in models
    ref = nothing
    for (tag, _) in arms(m)
        a = read_arm(m, tag)
        if a === nothing
            push!(rows, @sprintf("%-4d %-12s not run", m, tag)); continue
        elseif a.kv["status"] != "ok"
            push!(rows, @sprintf("%-4d %-12s FAILED: %s", m, tag, get(a.kv, "message", "")[1:min(end, 150)])); continue
        end
        ref === nothing && (ref = a)
        n  = min(size(a.G, 2), size(ref.G, 2))
        gd = a === ref ? 0.0 : sqrt(sum(abs2, a.G[:, 1:n] .- ref.G[:, 1:n]) / sum(abs2, ref.G[:, 1:n]))
        p(k) = parse(Float64, a.kv[k])
        push!(rows, @sprintf("%-4d %-12s %3d/%-3d %12.3e %12.6f %10.4f %8.2f %9.1f", m, tag, Int(p("nsteps")), Int(p("expected")),
                             gd, p("eratio"), p("emax") / Awave, p("newton"), p("wall")))
    end
end
open(joinpath(outroot, "summary.txt"), "w") do io
    println(io, "compare_mixed_vms.jl — box: 1 λ, $(ncell) cells/λ, Q$(p_u)/Q$(p_eta), P$(p_vert)LFE-$(M), d=$(d) m, T=$(Twave) s, ",
                @sprintf("kA=%.3g (A=%.4g m), CN dt=%.3g s, %g periods; variable bed: cosine, max slope %.3g", kA, Awave, dt, periods, slope))
    println(io, "gauge_diff = relative L2 difference of η at 4 gauges over the run, against the first arm of each model")
    println(io)
    foreach(r -> println(io, r), rows)
end
println("\n", join(rows, "\n"))
println("\nsummary → ", joinpath(outroot, "summary.txt"))
