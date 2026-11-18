# ==============================================================
#  vms_newton_probe.jl — which VMS term stops Newton? (branch mixed-vms-solver)
#
#  Same box as compare_mixed_vms.jl, model 3 (BALFEM_PROBE_MODEL), a few steps, one VMS coefficient
#  at a time and scaled. Prints, per variant, the Newton iterations of each step or the failure.
#  Output: output/mixed_vms/newton_probe/<variant>/.
# ==============================================================
using GridapBALFEM, Printf

const ROOT_P = normpath(joinpath(@__DIR__, "..", ".."))
m       = parse(Int, get(ENV, "BALFEM_PROBE_MODEL", "3"))
nsteps  = parse(Int, get(ENV, "BALFEM_PROBE_STEPS", "3"))
ncell, kA, dt, slope = 16, 0.08, 0.04, 0.05
M, p_vert, p_u, p_eta, d, Twave = 2, 1, 3, 2, 3.5, 1.6
vert0 = assemble_vertical_tensors(M, p_vert, Vector{Float64}(resolve_cbdy(M, nothing, p_vert)))
k0    = GridapBALFEM.model_wavenumber(vert0, 2π / Twave, d, 9.81)
Awave = kA / k0
wi    = WaveInput(vert0; A = Awave, T = Twave, d = d, T_ramp = 0.0, profile = :model)
Lx    = 2π / wi.ks[1]; dx = Lx / ncell
inc   = incident_fields(wi)
bedf  = x -> d + (slope / (2π / Lx)) * cos(2π / Lx * x[1])

variants = [("all1",    (; vms_c_adv = 1.0, vms_c_wave = 1.0, vms_c_eta = 1.0)),
            ("all1_ad", (; vms_c_adv = 1.0, vms_c_wave = 1.0, vms_c_eta = 1.0, use_ad = true)),
            ("adv1", (; vms_c_adv = 1.0, vms_c_wave = 0.0, vms_c_eta = 0.0)),
            ("wave1", (; vms_c_adv = 0.0, vms_c_wave = 1.0, vms_c_eta = 0.0)),
            ("eta1",  (; vms_c_adv = 0.0, vms_c_wave = 0.0, vms_c_eta = 1.0)),
            ("wave0.1", (; vms_c_adv = 0.0, vms_c_wave = 0.1, vms_c_eta = 0.0)),
            ("eta0.1",  (; vms_c_adv = 0.0, vms_c_wave = 0.0, vms_c_eta = 0.1)),
            ("all0.1",  (; vms_c_adv = 0.1, vms_c_wave = 0.1, vms_c_eta = 0.1))]
rows = String[]
for (tag, kw) in variants
    out = joinpath(ROOT_P, "output", "mixed_vms", "newton_probe", @sprintf("m%d_%s", m, tag))
    line = try
        diags, _, _ = setup_and_run(;
            M = M, p_vertical = p_vert, p_u = p_u, p_eta = p_eta,
            domain = ((0.0, Lx), (0.0, dx)), partition = (ncell, 1), x_periodic = true, y_wall_bc = :wall,
            h_val = d, T_wave = Twave, A_wave = Awave, sponge_wL = 0.0, sponge_wR = 0.0,
            eta0_func = x -> inc.eta(x, 0.0), ux0_func = x -> inc.ux(x, 0.0),
            T_final = nsteps * dt, dt = dt, solver_type = :theta, theta = 0.5,
            regime = :nonlinear, nl_pressure = m in (5, 6), flat_bed = isodd(m),
            h_bathy = isodd(m) ? nothing : bedf, formulation = m in (5, 6) ? :mixed : :primal,
            output_dir = out, save_every = 0, nl_tol = 1e-8, print_every = 1, check_every = 0,
            diag_every = 1, eta_ref = Awave, vms = true, kw...)
        @sprintf("%-8s Newton per step: %s   res: %s", tag, join((s.nl_iters for s in diags), " "),
                 join((@sprintf("%.1e", s.res_nl) for s in diags), " "))
    catch e
        @sprintf("%-8s FAILED: %s", tag, first(sprint(showerror, e), 120))
    end
    push!(rows, line); println("PROBE ", line); flush(stdout)
end
println("\n", join(rows, "\n"))
