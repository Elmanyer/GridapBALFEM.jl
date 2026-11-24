# abs_onset.jl — absolute-threshold burst onset (2026-11-23): the first period whose period-rms of the
# grid-band VS velocity exceeds 5 mm/s and 1 cm/s. Not usable at kA 0.32 or kd 1.7, where the bound
# grid-band content already exceeds 5 mm/s. Run after analyse.jl in the same session; writes
# $IS_OUT/abs_onset.tsv (pasted into events2.tsv by hand for pivot3.py).
open(joinpath(AN_DIR, "abs_onset.tsv"), "w") do io
    println(io, "name\tP_grid5mm\tP_grid1cm\tvs_grid_max\tper_end")
    for m in rows
        isfile(joinpath(AN_DIR, "_coeffs", m["name"] * ".jls")) || continue
        a = an_case(m["name"], m); np = floor(Int, a.per[end])
        np < 3 && (println(io, m["name"], "\tNaN\tNaN\tNaN\t", np); continue)
        s = [sqrt(mean(abs2, a.B["vs_grid"][(a.per .> p - 1) .& (a.per .≤ p)])) for p in 1:np]
        f(th) = (i = findfirst(>(th), s[3:end]); i === nothing ? "none" : string(i + 2))
        println(io, m["name"], '\t', f(5e-3), '\t', f(1e-2), '\t', @sprintf("%.3g", maximum(s[3:end])), '\t', np)
    end
end
