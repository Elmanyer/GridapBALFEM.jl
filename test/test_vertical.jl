# ==============================================================
#  test_vertical.jl — vertical tensor set (FAST, no horizontal FEM)
#
#  Checks on assemble_vertical_tensors:
#    identities (Mmat sym/PSD, ΣΦ=1, φ_int BCs), the leading-pressure
#    identities P[:,:,3] = −B and Kcal = Pcal − ∫σΘφᵢ, and the dispersion
#    bridge vs Yang & Liu (2024) Table 1 applicable-kd values.
#
#  RUN:  julia --project=. GridapBALFEM.jl/test/test_vertical.jl
# ==============================================================

using GridapBALFEM
using Gridap.TensorValues   # VectorValue — used below to evaluate φ_int at σ=0,1.
                            # GridapBALFEM does not re-export Gridap's names, so a
                            # consumer that constructs Gridap values imports them.
using LinearAlgebra, Printf

println("=" ^ 60)
println("  test_vertical.jl — vertical tensor set")
println("=" ^ 60)

n_pass = 0; n_fail = 0
function check(name, cond)
    global n_pass, n_fail
    if cond; println("  PASS  $name"); n_pass += 1
    else;    println("  FAIL  $name"); n_fail += 1; end
end

g = 9.81; d = 3.5

# ---- P1LFE-2 -------------------------------------------------------------------
vert2 = assemble_vertical_tensors(2, 1, [0.0, 0.728, 1.0])
check("P1LFE-2: N_dof = 3", vert2.N_dof == 3)
check("P1LFE-2: ΣΦ = 1", abs(sum(vert2.Phi) - 1.0) < 1e-12)
check("P1LFE-2: Mmat symmetric", norm(vert2.Mmat - vert2.Mmat') < 1e-13)
check("P1LFE-2: Mmat positive definite", all(eigvals(Symmetric(vert2.Mmat)) .> 0))
check("P1LFE-2: B symmetric ≤ 0", norm(vert2.B - vert2.B') < 1e-13 &&
                                all(eigvals(Symmetric(vert2.B)) .< 1e-14))

# antiderivative boundary values: φⱼ_int(0)=0, φⱼ_int(1)=Φⱼ
p0 = VectorValue(0.0); p1 = VectorValue(1.0)
check("P1LFE-2: φⱼ_int(0) = 0",
      all(abs(vert2.phi_int_fns[j](p0)) < 1e-12 for j in 1:vert2.N_dof))
check("P1LFE-2: φⱼ_int(1) = Φⱼ",
      all(abs(vert2.phi_int_fns[j](p1) - vert2.Phi[j]) < 1e-12 for j in 1:vert2.N_dof))

# leading-pressure identities
check("P1LFE-2: P[:,:,3] = −B (dispersion carrier identity)",
      norm(vert2.P[:,:,3] + vert2.B) < 1e-12)
check("P1LFE-2: Mcal fully symmetric in (i,k,j)",
      all(abs(vert2.Mcal[i,k,j] - vert2.Mcal[k,i,j]) < 1e-13 &&
          abs(vert2.Mcal[i,k,j] - vert2.Mcal[i,j,k]) < 1e-13
          for i in 1:3, k in 1:3, j in 1:3))

# Fubini split: Kcal = Pcal − ∫σΘφᵢ  ⇒  Pcal − Kcal = ∫σΘφᵢ; for the strictly
# positive components c=1..5 (Θ ≥ 0, σφᵢ ≥ 0) this difference must be ≥ 0.
dPK = vert2.Pcal .- vert2.Kcal
check("P1LFE-2: Pcal − Kcal = ∫σΘφᵢ ≥ 0 for components 1–5",
      all(dPK[:, :, :, 1:5] .> -1e-12))
check("P1LFE-2: Pcal assembled (nonzero)", maximum(abs.(vert2.Pcal)) > 1e-6)

# dispersion bridge — Yang & Liu (2024) Table 1
kd2 = applicable_kd(vert2, g, d)
@printf("  P1LFE-2 applicable kd = %.1f  (Table 1: ~10.9)\n", kd2)
check("P1LFE-2: applicable kd ≈ 10.9 (±1.0)", abs(kd2 - 10.9) < 1.0)

# ---- P1LFE-3 -------------------------------------------------------------------
vert3 = assemble_vertical_tensors(3, 1, [0.0, 0.726, 0.925, 1.0])
kd3 = applicable_kd(vert3, g, d)
@printf("  P1LFE-3 applicable kd = %.1f  (Table 1: ~39.2)\n", kd3)
check("P1LFE-3: ΣΦ = 1", abs(sum(vert3.Phi) - 1.0) < 1e-12)
check("P1LFE-3: applicable kd ≈ 39.2 (±2.0)", abs(kd3 - 39.2) < 2.0)
check("P1LFE-3: P[:,:,3] = −B", norm(vert3.P[:,:,3] + vert3.B) < 1e-12)

# ---- ★ the skew-symmetry identity, on EVERY basis ------------------------------
#  ½(𝓖ᵢₖⱼ + 𝓖ⱼₖᵢ) = ½𝓜ᵢₖⱼ − ½Φₖ Mᵢⱼ
#
#  This is the algebraic core of the energy-consistent advection reformulation
#  (building_files/SKEW_SYMMETRIC_ADVECTION_PLAN.md §1.2): it is what makes the
#  advection block's exact energy production collapse to −½∫∇·(Hū)(ΣMᵢⱼuᵢ·uⱼ),
#  i.e. to a pure continuity defect, with the 𝓜 pieces cancelling identically.
#
#  It holds because ψₖ = σΦₖ − varphiₖ VANISHES AT BOTH ENDS of the water column
#  (ψₖ(0)=0 by construction, ψₖ(1)=Φₖ−varphiₖ(1)=0), so the σ-integration by parts
#  carries no boundary term — a statement about the vertical basis alone. It is
#  therefore BASIS-AGNOSTIC, and this test checks that claim rather than assuming
#  it: five bases spanning p=1,2 and Nσ=3,4,5. If it ever fails, the correction in
#  problem.jl is no longer energy-consistent on that basis and must not be used.
println()
println("  ★ skew-symmetry identity  ½(𝓖ᵢₖⱼ+𝓖ⱼₖᵢ) = ½𝓜ᵢₖⱼ − ½Φₖ Mᵢⱼ")
for (nm, M, p, cb) in (("P1LFE-2", 2, 1, [0.0, 0.728, 1.0]),
                       ("P1LFE-3", 3, 1, [0.0, 0.726, 0.925, 1.0]),
                       ("P1LFE-4", 4, 1, [0.0, 0.745, 0.923, 0.977, 1.0]),
                       ("P2LFE-1", 1, 2, [0.0, 1.0]),
                       ("P2LFE-2", 2, 2, [0.0, 0.8298, 1.0]))
    v  = assemble_vertical_tensors(M, p, cb)
    N  = v.N_dof
    Gc = v.Gcal; Mc = v.Mcal; Mm = v.Mmat; Ph = v.Phi
    star = maximum(abs(0.5*(Gc[i,k,j] + Gc[j,k,i]) - 0.5*Mc[i,k,j] + 0.5*Ph[k]*Mm[i,j])
                   for i in 1:N, k in 1:N, j in 1:N)
    #  𝓜 must be FULLY symmetric — the other half of the cancellation.
    asym = maximum(max(abs(Mc[i,k,j] - Mc[k,i,j]), abs(Mc[i,k,j] - Mc[i,j,k]))
                   for i in 1:N, k in 1:N, j in 1:N)
    @printf("    %-8s Nσ=%d   max|★| = %.2e   max|𝓜 asym| = %.2e\n", nm, N, star, asym)
    check("$nm: ★ identity holds (< 1e-12)", star < 1e-12)
    check("$nm: 𝓜 fully symmetric (< 1e-13)", asym < 1e-13)
end

println()
println("=" ^ 60)
@printf("  Results: %d PASS,  %d FAIL\n", n_pass, n_fail)
println("=" ^ 60)
n_fail > 0 ? error("test_vertical: $n_fail failed!") :
             println("  Vertical tensor set validated.")
