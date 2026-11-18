# ==============================================================
#  oss_newton_trace.jl — Newton trace of the first steps of mixed + OSS (branch mixed-vms-solver)
#
#  Same box as compare_mixed_vms.jl. m5/m6 mixed + OSS reach NaN at step 2 (2026-11-18, with the old
#  and with the completed jacobian_u). This prints the NLsolve trace of every iteration so the failure
#  can be read: Jacobian stall (slow contraction) or iterates running off (blow-up).
#  Knobs: BALFEM_PROBE_MODEL (5), BALFEM_PROBE_STEPS (2), BALFEM_PROBE_VARIANT (all1 | adv_eta | wave1 |
#  all0.1), BALFEM_PROBE_OSS (1), BALFEM_PROBE_NLTOL (1e-8), BALFEM_PROBE_TAG (output-name suffix),
#  BALFEM_VMS_DT (0; 1 adds (2/Δt)² to τ_a).
#  Output: output/mixed_vms/newton_probe/oss_m<m>_<variant>/.
# ==============================================================
using GridapBALFEM, Printf

const ROOT_P = normpath(joinpath(@__DIR__, "..", ".."))
m       = parse(Int, get(ENV, "BALFEM_PROBE_MODEL", "5"))
nsteps  = parse(Int, get(ENV, "BALFEM_PROBE_STEPS", "2"))
vname   = get(ENV, "BALFEM_PROBE_VARIANT", "all1")
oss     = get(ENV, "BALFEM_PROBE_OSS", "1") == "1"
nltol   = parse(Float64, get(ENV, "BALFEM_PROBE_NLTOL", "1e-8"))
suffix  = get(ENV, "BALFEM_PROBE_TAG", "")
dtterm  = get(ENV, "BALFEM_VMS_DT", "0") == "1"      # add (2/Δt)² to τ_a (bounds τ_a ≤ Δt/2)
ncell, kA, dt, slope = 16, 0.08, 0.04, 0.05
M, p_vert, p_u, p_eta, d, Twave = 2, 1, 3, 2, 3.5, 1.6
vert0 = assemble_vertical_tensors(M, p_vert, Vector{Float64}(resolve_cbdy(M, nothing, p_vert)))
k0    = GridapBALFEM.model_wavenumber(vert0, 2π / Twave, d, 9.81)
Awave = kA / k0
wi    = WaveInput(vert0; A = Awave, T = Twave, d = d, T_ramp = 0.0, profile = :model)
Lx    = 2π / wi.ks[1]; dx = Lx / ncell
inc   = incident_fields(wi)
bedf  = x -> d + (slope / (2π / Lx)) * cos(2π / Lx * x[1])

coeffs = Dict("all1"    => (1.0, 1.0, 1.0), "adv_eta" => (1.0, 0.0, 1.0), "wave1" => (0.0, 1.0, 0.0),
              "adv1"    => (1.0, 0.0, 0.0), "eta1" => (0.0, 0.0, 1.0), "all0.1" => (0.1, 0.1, 0.1))
ca, cw, ce = coeffs[vname]
out = joinpath(ROOT_P, "output", "mixed_vms", "newton_probe", @sprintf("%s_m%d_%s%s", oss ? "oss" : "asgs", m, vname, suffix))
println("PROBE m$m ", oss ? "OSS" : "ASGS", " $vname  C_a=$ca C_w=$cw C_η=$ce  Δt-term=$dtterm  → $out"); flush(stdout)
line = try
    diags, _, _ = setup_and_run(;
        M = M, p_vertical = p_vert, p_u = p_u, p_eta = p_eta,
        domain = ((0.0, Lx), (0.0, dx)), partition = (ncell, 1), x_periodic = true, y_wall_bc = :wall,
        h_val = d, T_wave = Twave, A_wave = Awave, sponge_wL = 0.0, sponge_wR = 0.0,
        eta0_func = x -> inc.eta(x, 0.0), ux0_func = x -> inc.ux(x, 0.0),
        T_final = nsteps * dt, dt = dt, solver_type = :theta, theta = 0.5,
        regime = :nonlinear, nl_pressure = m in (5, 6), flat_bed = isodd(m),
        h_bathy = isodd(m) ? nothing : bedf, formulation = m in (5, 6) ? :mixed : :primal,
        output_dir = out, save_every = 0, nl_tol = nltol, print_every = 1, check_every = 0,
        diag_every = 1, eta_ref = Awave, show_trace = true,
        vms = true, vms_oss = oss, vms_c_adv = ca, vms_c_wave = cw, vms_c_eta = ce, vms_dt_term = dtterm)
    @sprintf("Newton per step: %s   res: %s", join((s.nl_iters for s in diags), " "),
             join((@sprintf("%.1e", s.res_nl) for s in diags), " "))
catch e
    "FAILED: " * first(sprint(showerror, e), 160)
end
println("PROBE-RESULT m$m ", oss ? "OSS" : "ASGS", " $vname: ", line)
