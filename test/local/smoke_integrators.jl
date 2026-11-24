# ==============================================================
#  smoke_integrators.jl — every time integrator through the closed-box driver, 2 steps each
#
#  WHY (2026-11-11). The time-integrator dissipation study (STABILITY.md §6b) will run on the cluster
#  from the system image, varying BALFEM_SOLVER / BALFEM_TABLEAU / BALFEM_THETA_T / BALFEM_RHO_INF.
#  Before any cluster hour is spent, every choice must compile and run on the real entry point,
#  examples/local_1d/run_periodic_1d.jl, for the models of the study.
#
#  WHAT. For each model, the driver is included with the environment of one case (Q3/Q2, 8 cells,
#  one wavelength, A = 0.10, T = 1.6 s, Δt = 0.04, 2 steps):
#    * θ-method θ = 0.5 (CN) and 0.6; generalised-α ρ∞ = 0, 0.5, 1;
#    * every diagonally implicit and every explicit Gridap tableau (Tableaus{DIM,EX}.jl, read from the
#      Gridap source, so a tableau added to the fork is tested too);
#    * refusals that must fire: an IMEX tableau, an unknown tableau name, ρ∞ = 1.5;
#    * checkpoints: a generalised-α run resumed from its own checkpoint (must pass), the same
#      checkpoint resumed with another ρ∞ without RESTART_ALLOW_CHANGE (must be refused), and a
#      checkpoint written BEFORE ρ∞ existed (campaign m3_q32_n16, t = 320 s) resumed under CN
#      (must pass: the fingerprint is backward compatible).
#  A case passes when the driver reports `done: 2 steps` and the last step is finite (or, for a
#  refusal, when it errors with the expected message). Gate lines start with `SMOKE `; the verdict is
#  the `SMOKE SUMMARY` line (rule 35), never the exit code.
#
#  KNOBS:  SMOKE_MODELS  comma list of m3,m5,m6           default "m3,m5,m6"
#          SMOKE_FULL    models that get EVERY tableau     default "m3,m5" (the others: one per kind)
#          SMOKE_ONLY    run only cases whose name contains this text
#  RUN (detached, ~1–3 h cold):
#    setsid nohup julia --project=. test/local/smoke_integrators.jl > output/smoke_integrators/_smoke.out 2>&1 < /dev/null &
# ==============================================================
using GridapBALFEM
using Gridap
using Logging
using Printf

const ROOT   = normpath(joinpath(@__DIR__, "..", ".."))
const DRIVER = joinpath(ROOT, "examples", "local_1d", "run_periodic_1d.jl")
const OUT    = joinpath(ROOT, "output", "smoke_integrators")
mkpath(OUT)

const MODELS = Dict(
    "m3" => Dict("BALFEM_NL_PRESSURE" => "0", "BALFEM_FLAT_BED" => "1", "BALFEM_BED_SLOPE" => "0"),
    "m5" => Dict("BALFEM_NL_PRESSURE" => "1", "BALFEM_FLAT_BED" => "1", "BALFEM_BED_SLOPE" => "0"),
    "m6" => Dict("BALFEM_NL_PRESSURE" => "1", "BALFEM_FLAT_BED" => "0", "BALFEM_BED_SLOPE" => "0.1"))
const BASE = Dict(
    "BALFEM_REGIME" => "nonlinear", "BALFEM_M" => "2", "BALFEM_P_VERT" => "1", "BALFEM_BED_NW" => "1",
    "BALFEM_FE_ORDER" => "3", "BALFEM_P_ETA" => "2", "BALFEM_NCELL" => "8", "BALFEM_NLAMBDA" => "1",
    "BALFEM_AWAVE" => "0.10", "BALFEM_TWAVE" => "1.6", "BALFEM_D" => "3.5", "BALFEM_DT" => "0.04",
    "BALFEM_PERIODS" => "0.05", "BALFEM_NOISE" => "1e-8", "BALFEM_USE_AD" => "0",
    "BALFEM_SAVE_EVERY" => "1", "BALFEM_DIAG_EVERY" => "1", "BALFEM_PRINT_EVERY" => "1",
    "BALFEM_CHECKPOINT_EVERY" => "0")
#  every key a case may set: cleared before each case so nothing leaks from the previous one
const CASE_KEYS = ["BALFEM_SOLVER", "BALFEM_TABLEAU", "BALFEM_THETA_T", "BALFEM_RHO_INF",
                   "BALFEM_RESTART_FROM", "BALFEM_RESTART_ALLOW_CHANGE", "BALFEM_CIP_GU", "BALFEM_CIP_GE"]

"Tableau names declared in one of Gridap's tableau files."
function gridap_tableaus(file)
    src = read(joinpath(pkgdir(Gridap), "src", "ODEs", "ODESolvers", file), String)
    return [m[1] for m in eachmatch(r"struct (\w+) <: TableauName", src)]
end

csv(s) = [strip(x) for x in split(s, ',') if !isempty(strip(x))]
const SEL_MODELS = csv(get(ENV, "SMOKE_MODELS", "m3,m5,m6"))
const FULL       = csv(get(ENV, "SMOKE_FULL", "m3,m5"))
const ONLY       = get(ENV, "SMOKE_ONLY", "")

#  (name, env, expect): expect = :run, or a substring the refusal message must contain
function cases_for(model)
    c = Tuple{String,Dict{String,String},Any}[]
    push!(c, ("theta0.50", Dict("BALFEM_SOLVER" => "theta", "BALFEM_THETA_T" => "0.5"), :run))
    push!(c, ("theta0.60", Dict("BALFEM_SOLVER" => "theta", "BALFEM_THETA_T" => "0.6"), :run))
    for r in ("0.0", "0.5", "1.0")
        push!(c, ("galpha$r", Dict("BALFEM_SOLVER" => "gen_alpha", "BALFEM_RHO_INF" => r), :run))
    end
    tabs = model in FULL ? vcat(gridap_tableaus("TableausDIM.jl"), gridap_tableaus("TableausEX.jl")) :
                           ["SDIRK_2_2", "DIRK_TRBDF_3_2", "EXRK_RungeKutta_4_4"]
    for tb in tabs
        push!(c, (tb, Dict("BALFEM_SOLVER" => "sdirk", "BALFEM_TABLEAU" => tb), :run))
    end
    push!(c, ("rk_alias", Dict("BALFEM_SOLVER" => "rk", "BALFEM_TABLEAU" => "SDIRK_4_3"), :run))
    push!(c, ("refuse_imex", Dict("BALFEM_SOLVER" => "sdirk", "BALFEM_TABLEAU" => "IMEXRK_2_2_2"), "IMEX"))
    push!(c, ("refuse_name", Dict("BALFEM_SOLVER" => "sdirk", "BALFEM_TABLEAU" => "SDIRK_9_9"),
              "not a Gridap Butcher tableau"))
    push!(c, ("refuse_rho", Dict("BALFEM_SOLVER" => "gen_alpha", "BALFEM_RHO_INF" => "1.5"), "rho_inf must lie"))
    return c
end

const RESULTS = Tuple{String,Bool,String}[]
function gate(name, ok, msg)
    push!(RESULTS, (name, ok, msg))
    println(ok ? "SMOKE PASS " : "SMOKE FAIL ", rpad(name, 40), msg); flush(stdout)
end

"Include the driver with `env` (on top of BASE and the model); return (log text, error or nothing)."
function run_driver(name, env; dir = joinpath(OUT, replace(name, r"[^A-Za-z0-9_.-]" => "_")))
    for k in CASE_KEYS; delete!(ENV, k); end
    for (k, v) in env; ENV[k] = v; end
    ENV["BALFEM_OUTDIR"] = dir
    mkpath(dir)
    logf = joinpath(dir, "smoke.log")
    err = nothing
    open(logf, "w") do io
        redirect_stdio(stdout = io, stderr = io) do
            with_logger(ConsoleLogger(io)) do
                try
                    Base.include(Main, DRIVER)
                catch e
                    err = e
                    showerror(io, e)
                end
            end
        end
    end
    return read(logf, String), err
end

"Last per-step line of a driver log: (finite?, seconds per step)."
function last_step(log)
    lines = filter(l -> startswith(l, "step "), split(log, '\n'))
    isempty(lines) && return (false, NaN)
    l = lines[end]
    m = match(r"eta_max=([-+0-9.eE]+|NaN|Inf)", l); s = match(r"\|\s+([0-9.]+) s \|", l)
    fin = m !== nothing && isfinite(parse(Float64, m[1]))
    return (fin, s === nothing ? NaN : parse(Float64, s[1]))
end

function judge(name, log, err, expect)
    if expect === :run
        err === nothing || return gate(name, false, "ERROR: " * first(sprint(showerror, err), 200))
        m = match(r"done: (\d+) steps", log)
        m === nothing && return gate(name, false, "no `done:` line")
        fin, sps = last_step(log)
        n = parse(Int, m[1])
        warn = occursin("amplifies neutral", log) ? "  [warned: amplifies oscillations]" :
               occursin("is defective", log) ? "  [warned: defective tableau]" : ""
        gate(name, n == 2 && fin, @sprintf("%d steps, finite=%s, %.2f s/step%s", n, fin, sps, warn))
    else
        msg = err === nothing ? "" : sprint(showerror, err)
        gate(name, err !== nothing && occursin(expect, msg),
             err === nothing ? "ran, but a refusal was expected" : "refused: " * first(msg, 120))
    end
end

for model in SEL_MODELS
    for (k, v) in merge(BASE, MODELS[model]); ENV[k] = v; end
    for (cname, env, expect) in cases_for(model)
        name = "$model $cname"
        isempty(ONLY) || occursin(ONLY, name) || continue
        t0 = time()
        log, err = run_driver(name, env)
        judge(name, log, err, expect)
        @printf("      (%.0f s wall)\n", time() - t0)
    end
end

# ---- checkpoints: generalised-α round trip, its refusal, and backward compatibility ----------
if "m3" in SEL_MODELS && (isempty(ONLY) || occursin("ck", ONLY))
    for (k, v) in merge(BASE, MODELS["m3"]); ENV[k] = v; end
    ga = Dict("BALFEM_SOLVER" => "gen_alpha", "BALFEM_RHO_INF" => "0.5")
    ENV["BALFEM_CHECKPOINT_EVERY"] = "1"
    log, err = run_driver("m3 ck_write galpha0.5", ga)
    judge("m3 ck_write galpha0.5", log, err, :run)
    ENV["BALFEM_CHECKPOINT_EVERY"] = "0"
    src = joinpath(OUT, "m3_ck_write_galpha0.5")
    ENV["BALFEM_PERIODS"] = "0.10"                                   # T_final is absolute: 2 more steps
    log, err = run_driver("m3 ck_resume galpha0.5", merge(ga, Dict("BALFEM_RESTART_FROM" => src)))
    judge("m3 ck_resume galpha0.5", log, err, :run)
    log, err = run_driver("m3 ck_resume galpha0.3 refused",
                          merge(ga, Dict("BALFEM_RHO_INF" => "0.3", "BALFEM_RESTART_FROM" => src)))
    judge("m3 ck_resume galpha0.3 refused", log, err, "time discretisation differs")
    #  a checkpoint written before ρ∞ existed: the campaign's m3 n16 box at t = 320 s, resumed under CN
    old = joinpath(ROOT, "output", "1D_stability_campaign", "m3_q32_n16")
    if isdir(joinpath(old, "checkpoints"))
        ENV["BALFEM_NCELL"] = "16"; ENV["BALFEM_PERIODS"] = string(200 + 2 * 0.04 / 1.6)
        log, err = run_driver("m3 ck_resume_pre_rho theta",
                              Dict("BALFEM_SOLVER" => "theta", "BALFEM_RESTART_FROM" => old))
        judge("m3 ck_resume_pre_rho theta", log, err, :run)
    else
        gate("m3 ck_resume_pre_rho theta", false, "no checkpoint at $old")
    end
end

npass = count(r -> r[2], RESULTS)
println("\nSMOKE SUMMARY: $npass / $(length(RESULTS)) passed")
for r in RESULTS
    r[2] || println("  FAILED: ", r[1], " — ", r[3])
end
