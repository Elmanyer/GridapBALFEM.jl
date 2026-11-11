# Every figure and animation of a closed x-periodic box run (examples/local_1d/run_periodic_1d.jl),
# written to <run>/viz/ (or --out). The functions live in src/periodic.jl and src/periodic_plots.jl.
#
#   julia --project=postprocessing postprocessing/examples/periodic_viz.jl <run_dir> [options]
#
# options
#   --no-anim              figures only
#   --parts total,disturbance,steady,oscillating,16      what to animate (a number = that harmonic)
#   --dt 1.0               one animation frame per dt seconds of simulated time
#   --trange 150 246       animate only this time window
#   --phase 0.0            stroboscopic animation: one frame per wave period at this phase (0–1)
#   --fps 12               frames per second of the video
#   --out DIR              output folder (default <run_dir>/viz)
#
# From julia-mcp (postprocessing env):  empty!(ARGS); append!(ARGS, ["<run_dir>", "--no-anim"]);
#                                       include("postprocessing/examples/periodic_viz.jl")
(isdefined(Main, :GridapBALFEMPost) && isdefined(Main.GridapBALFEMPost, :load_periodic)) ||
    include(joinpath(@__DIR__, "..", "src", "GridapBALFEMPost.jl"))
using .GridapBALFEMPost      # calls below are qualified: they survive a re-include in a long-lived session

function periodic_viz_main(args)
    isempty(args) && error("usage: periodic_viz.jl <run_dir> [--no-anim] [--parts …] [--dt …] [--trange a b] [--phase φ] [--fps n] [--out DIR]")
    dir = args[1]; anim = true; parts = Any[:total, :disturbance]; dt = 1.0; trange = nothing; phase = nothing; fps = 12; out = nothing
    i = 2
    while i ≤ length(args)
        a = args[i]
        if a == "--no-anim"
            anim = false
        elseif a == "--parts"
            parts = Any[all(isdigit, s) ? [parse(Int, s)] : Symbol(s) for s in split(args[i+1], ",")]; i += 1
        elseif a == "--dt"
            dt = parse(Float64, args[i+1]); i += 1
        elseif a == "--trange"
            trange = (parse(Float64, args[i+1]), parse(Float64, args[i+2])); i += 2
        elseif a == "--phase"
            phase = parse(Float64, args[i+1]); i += 1
        elseif a == "--fps"
            fps = parse(Int, args[i+1]); i += 1
        elseif a == "--out"
            out = args[i+1]; i += 1
        else
            error("periodic_viz: unknown option $a")
        end
        i += 1
    end
    r = GridapBALFEMPost.load_periodic(dir); println(r)
    outdir = something(out, joinpath(dir, "viz"))
    GridapBALFEMPost.periodic_report(r; outdir, animate = false)
    if anim
        tr = something(trange, (r.t[1], r.t[end]))
        frames = phase === nothing ? GridapBALFEMPost.select_frames(r; trange = tr, dt) : GridapBALFEMPost.select_frames(r; trange = tr, phase)
        for part in parts
            tag = part isa Symbol ? String(part) : "n$(part[1])"
            f = GridapBALFEMPost.animate_run(r; frames, part, fps, out = joinpath(outdir, "state_$tag.mp4"))
            println("  ", basename(f), "  (", length(frames), " frames)")
        end
    end
    return outdir
end

periodic_viz_main(ARGS)
