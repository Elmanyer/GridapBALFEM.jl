# ==============================================================
#  bench_depth_memo.jl — the depth-field memo (2026-11-11): bit-identical, and how much faster
#
#  WHAT. In ONE process, each box case runs twice through the real driver
#  (examples/local_1d/run_periodic_1d.jl): GridapBALFEM.DEPTH_MEMO[] = false (the depth fields rebuilt
#  at every call site, as before 2026-11-11) and = true (built once, `depth_fields`). A warm-up run
#  per case pays the compilation first, so both timed runs are warm. Each timed run writes a final
#  checkpoint, and the two checkpoints' raw Float64 state files must be byte-identical.
#  Cases: m5 (flat bed, the profiled case), m6 (cosine bed: ∇h and the bed Hessian live), m3.
#  Run it once in debug mode and once with Gridap in performance mode (LocalPreferences.toml
#  [Gridap] execution_mode = "performance"): BENCH_TAG labels the two.
#
#  Gate lines: `BENCH <case> identical=<true|false> off=<s/step> on=<s/step> speed-up=<x>`.
#  KNOBS: BENCH_STEPS (default 10), BENCH_MODELS (default "m5,m6,m3"), BENCH_TAG (default "debug").
#  RUN:   setsid nohup julia --project=. test/local/bench_depth_memo.jl > output/bench_depth_memo/<tag>.out 2>&1 < /dev/null &
# ==============================================================
using GridapBALFEM
using Logging
using Printf

const ROOT   = normpath(joinpath(@__DIR__, "..", ".."))
const DRIVER = joinpath(ROOT, "examples", "local_1d", "run_periodic_1d.jl")
const TAG    = get(ENV, "BENCH_TAG", "debug")
const OUT    = joinpath(ROOT, "output", "bench_depth_memo", TAG)
const STEPS  = parse(Int, get(ENV, "BENCH_STEPS", "10"))
mkpath(OUT)

const MODELS = Dict(
    "m3" => Dict("BALFEM_NL_PRESSURE" => "0", "BALFEM_FLAT_BED" => "1", "BALFEM_BED_SLOPE" => "0"),
    "m5" => Dict("BALFEM_NL_PRESSURE" => "1", "BALFEM_FLAT_BED" => "1", "BALFEM_BED_SLOPE" => "0"),
    "m6" => Dict("BALFEM_NL_PRESSURE" => "1", "BALFEM_FLAT_BED" => "0", "BALFEM_BED_SLOPE" => "0.1"))
const BASE = Dict(
    "BALFEM_REGIME" => "nonlinear", "BALFEM_M" => "2", "BALFEM_P_VERT" => "1", "BALFEM_BED_NW" => "1",
    "BALFEM_FE_ORDER" => "3", "BALFEM_P_ETA" => "2", "BALFEM_NCELL" => "8", "BALFEM_NLAMBDA" => "1",
    "BALFEM_AWAVE" => "0.10", "BALFEM_TWAVE" => "1.6", "BALFEM_D" => "3.5", "BALFEM_DT" => "0.04",
    "BALFEM_NOISE" => "1e-8", "BALFEM_USE_AD" => "0", "BALFEM_SOLVER" => "theta",
    "BALFEM_SAVE_EVERY" => "0", "BALFEM_DIAG_EVERY" => "5", "BALFEM_PRINT_EVERY" => "5",
    "BALFEM_CHECKPOINT_EVERY" => "100000", "BALFEM_CHECKPOINT_KEEP" => "1")

"Run the driver quietly into `dir` for `nsteps` steps; return the wall time of the whole run."
function run_case(dir, nsteps)
    ENV["BALFEM_PERIODS"] = string(nsteps * 0.04 / 1.6)
    ENV["BALFEM_OUTDIR"] = dir
    mkpath(dir)
    t0 = time()
    open(joinpath(dir, "run.log"), "w") do io
        redirect_stdio(stdout = io, stderr = io) do
            with_logger(ConsoleLogger(io)) do
                Base.include(Main, DRIVER)
            end
        end
    end
    return time() - t0
end

"The raw state file of the last checkpoint of a run."
final_state(dir) = (ck = joinpath(dir, "checkpoints");
                    f = sort(filter(endswith(".f64"), readdir(ck)))[end]; read(joinpath(ck, f)))

for model in split(get(ENV, "BENCH_MODELS", "m5,m6,m3"), ',')
    for (k, v) in merge(BASE, MODELS[model]); ENV[k] = v; end
    GridapBALFEM.DEPTH_MEMO[] = true
    run_case(joinpath(OUT, "$(model)_warmup"), 2)                     # compile everything first
    GridapBALFEM.DEPTH_MEMO[] = false
    t_off = run_case(joinpath(OUT, "$(model)_memo_off"), STEPS)
    GridapBALFEM.DEPTH_MEMO[] = true
    t_on  = run_case(joinpath(OUT, "$(model)_memo_on"), STEPS)
    same  = final_state(joinpath(OUT, "$(model)_memo_off")) == final_state(joinpath(OUT, "$(model)_memo_on"))
    @printf("BENCH %-3s %-11s identical=%s  off=%.2f s/step  on=%.2f s/step  speed-up=%.2fx\n",
            model, TAG, same, t_off / STEPS, t_on / STEPS, t_off / t_on)
    flush(stdout)
end
println("BENCH DONE ($TAG)")
