# ==============================================================
#  sponge_worker.jl — one worker of the sponge-reflection sweep (CH8_CAMPAIGN_TODO.md §7)
#
#  Started by run/local/run_1d_sponge_sweep.sh, never by hand. A worker runs up to MAX_CASES
#  cases of the queue in ONE Julia process, so the solver is compiled once per worker instead of
#  once per case (84 cases × ~15 min of JIT otherwise), and then exits: a bounded lifetime,
#  because long-lived Julia+Gridap processes degrade (rule 41). Workers pull cases from a shared
#  queue through claim directories (rule 44), so the cost-sorted queue is never sliced (rule 43).
#
#    julia --project=. run/local/sponge_worker.jl <cases.tsv> <out_dir> <log_dir> <max_cases>
#
#  Per case: every BALFEM_* variable is cleared (the flume driver `get!`s defaults into ENV, which
#  would otherwise leak into the next case), the common settings below and the case row are set,
#  and examples/local_1d/run_flume_1d.jl is included into a fresh module with stdout, stderr and
#  the logger sent to <log_dir>/<case>.log. Verdict from the log, never from an exit code
#  (rule 35): a finished case prints `flume_1d [...] done:`.
#  A claim is `<log_dir>/_claims/<case>/` holding `pid` (worker PID, host, start time). A case
#  whose claim exists is skipped; to re-run one, delete its claim directory, its log and its
#  output directory.
# ==============================================================
using Printf, Dates, Logging

length(ARGS) == 4 || error("usage: sponge_worker.jl <cases.tsv> <out_dir> <log_dir> <max_cases>")
const CASES_TSV, OUT_DIR, LOG_DIR = ARGS[1], ARGS[2], ARGS[3]
const MAX_CASES = parse(Int, ARGS[4])
const PROJ   = normpath(joinpath(@__DIR__, "..", ".."))
const DRIVER = joinpath(PROJ, "examples", "local_1d", "run_flume_1d.jl")
const CLAIMS = joinpath(LOG_DIR, "_claims")
mkpath(CLAIMS)

#  The Snellius linear suite (run/linear_production/_linear_common.sh, run_lin_bcplane_flat.sh),
#  sequential and 1-D: P1LFE-2, Q2/Q1, CN, dt = 0.02, d = 3.5, flat bed, Dirichlet generation of a
#  regular wave with the :model profile and the inflow relaxation zone, ny = 1 with walls,
#  dx = 0.25 (Ly = dx). VTK of η only, every 5 steps (0.1 s: a whole number of samples per period
#  for every T of the sweep).
const COMMON = Dict(
    "BALFEM_MPI" => "0", "BALFEM_WAVE_GEN" => "bc", "BALFEM_REGIME" => "linear",
    "BALFEM_NL_PRESSURE" => "0", "BALFEM_FLAT_BED" => "1",
    "BALFEM_M" => "2", "BALFEM_P_VERT" => "1", "BALFEM_FE_ORDER" => "2",
    "BALFEM_D" => "3.5", "BALFEM_DT" => "0.02", "BALFEM_SOLVER" => "theta",
    "BALFEM_AWAVE" => "0.1", "BALFEM_LY" => "0.25", "BALFEM_NY" => "1", "BALFEM_YBC" => "wall",
    "BALFEM_RELAX" => "1", "BALFEM_SAVE_EVERY" => "5", "BALFEM_PRINT_EVERY" => "250",
    "BALFEM_WRITE_W" => "0", "BALFEM_WRITE_PRESSURE" => "0",
)

function read_cases(path)
    lines = filter(!isempty, readlines(path))
    hdr = split(lines[1], '\t')
    return [Dict(zip(hdr, split(l, '\t'))) for l in lines[2:end]]
end

"Atomic claim: mkdir fails if another worker got there first."
function claim!(name)
    dir = joinpath(CLAIMS, name)
    try
        mkdir(dir)
    catch
        return false
    end
    write(joinpath(dir, "pid"), "$(getpid()) $(gethostname()) $(now())\n")
    return true
end

function run_case(c)
    name = c["name"]
    for k in collect(keys(ENV))
        startswith(k, "BALFEM_") && delete!(ENV, k)
    end
    merge!(ENV, COMMON)
    ENV["BALFEM_TWAVE"]    = c["T"]
    ENV["BALFEM_SPONGE_R"] = c["W"]
    ENV["BALFEM_MUMAX"]    = c["mu"]
    ENV["BALFEM_RELAX_W"]  = c["relax_w"]
    ENV["BALFEM_LX"]       = c["Lx"]
    ENV["BALFEM_NX"]       = c["nx"]
    ENV["BALFEM_TFINAL"]   = c["t_final"]
    ENV["BALFEM_OUTDIR"]   = joinpath(OUT_DIR, name)
    m = Module(Symbol("case_", replace(name, r"[^A-Za-z0-9_]" => "_")))
    Core.eval(m, :(include(p) = Base.include($m, p)))
    ok = false
    t0 = time()
    open(joinpath(LOG_DIR, "$name.log"), "w") do io
        println(io, "case $name | T=$(c["T"]) W=$(c["W"]) mu=$(c["mu"]) relax_w=$(c["relax_w"]) " *
                    "Lx=$(c["Lx"]) nx=$(c["nx"]) t_final=$(c["t_final"]) role=$(c["role"])")
        println(io, "started $(now()) | worker pid $(getpid())"); flush(io)
        redirect_stdio(stdout = io, stderr = io) do
            with_logger(ConsoleLogger(io)) do
                try
                    Base.include(m, DRIVER)
                    ok = true
                catch e
                    println(io, "ERROR in case $name:")
                    showerror(io, e, catch_backtrace())
                    println(io)
                end
            end
        end
        @printf(io, "--- case %s at %s (%.0f s) ---\n", ok ? "finished" : "FAILED", now(), time() - t0)
    end
    return ok
end

ndone = 0
for c in read_cases(CASES_TSV)
    ndone >= MAX_CASES && break
    isfile(joinpath(LOG_DIR, c["name"] * ".log")) && continue     # run before, by any worker
    claim!(c["name"]) || continue
    @printf("%s  worker %d: start %s\n", now(), getpid(), c["name"]); flush(stdout)
    ok = run_case(c)
    global ndone += 1
    @printf("%s  worker %d: %s %s\n", now(), getpid(), ok ? "done" : "FAILED", c["name"]); flush(stdout)
    GC.gc()
end
@printf("%s  worker %d: exiting after %d case(s)\n", now(), getpid(), ndone)
