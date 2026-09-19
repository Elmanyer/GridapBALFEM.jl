# ==============================================================
#  test_class3_reduction.jl — the EXACT algebraic reduction of Class-III {1,2,5}
#
#  Gates the identity implemented by `alg_class3_weight` (src/tensors.jl) and used by
#  `nlp_gradH_reduced_contrib` / `nlp_P_reduced_contrib`:
#
#      Σ_{l∈{1,2,5}} T^(l) ⊙ N^(l)  =  W ⊙ (Σ_a G_a⊗U_a)  +  T^(2) ⊙ (S⊗DU)
#      W[i,k,j] = −T1[i,k,j] + T2[i,k,j] − T5[i,j,k]
#
#  Derivation and rationale: markdown_files/NEW_TREATMENT.md §A.2.
#
#  ⚠ G4 IS THE POINT OF THIS FILE. The reduction's one plausible failure mode is dropping
#  the (k,j) transpose on T5 — the result still type-checks, still runs, and is wrong.
#  G4 asserts the comparison can DETECT that (rule 34: a test validates a term only if it
#  can resolve that term's contribution). G2/G3 mean nothing without it.
#
#  Runs on the REAL vertical tensors from `assemble_vertical`, not random ones, over three
#  vertical bases — the identity is index algebra and must not depend on Nσ or on p.
# ==============================================================
module TestClass3Reduction

using Test
using Gridap
using Gridap.TensorValues
using GridapBALFEM
using Random

const NGATE = Ref(0)
gate(name, cond) = (NGATE[] += 1; @testset "$name" begin @test cond end)

"Direct form: the four-object assembly as `nlp_frozen_N` builds it."
function direct_sum(T1, T2, T5, Gx, Gy, Ux, Uy, S, DU)
    N1 = -1.0*(outer(Gx, Ux) + outer(Gy, Uy))
    N2 = -1.0*N1 + outer(S, DU)
    N5 = -1.0*(outer(Ux, Gx) + outer(Uy, Gy))
    return double_contraction(T1, N1) + double_contraction(T2, N2) +
           double_contraction(T5, N5)
end

"Reduced form: one W-contraction plus the admissible S⊗DU remainder."
function reduced_sum(W, T2, Gx, Gy, Ux, Uy, S, DU)
    GU = outer(Gx, Ux) + outer(Gy, Uy)
    return double_contraction(W, GU) + double_contraction(T2, outer(S, DU))
end

relerr(a, b, N) = maximum(abs.([a[i] - b[i] for i in 1:N])) /
                  max(maximum(abs.([a[i] for i in 1:N])), 1e-300)

@testset "Class-III algebraic reduction" begin

    # ---- G1: pin the contraction convention the derivation assumes -----------
    #  If Gridap ever contracted the LEADING two indices instead, every gate below
    #  would still pass against a consistently-wrong reduction. Pin it explicitly.
    let N = 3
        Random.seed!(11)
        T = ThirdOrderTensorValue{N,N,N}(rand(N,N,N)...)
        M = TensorValue{N,N}(rand(N,N)...)
        got  = double_contraction(T, M)
        want = [sum(T[i,k,j]*M[k,j] for k in 1:N, j in 1:N) for i in 1:N]
        gate("G1 double_contraction contracts the TRAILING two indices",
             maximum(abs.([got[i] - want[i] for i in 1:N])) < 1e-14)
    end

    # ---- G2/G3/G5: the real tensors, three vertical bases --------------------
    for (label, M, p) in (("P1LFE-2", 2, 1), ("P1LFE-3", 3, 1), ("P2LFE-2", 2, 2))
        vert = assemble_vertical_tensors(M, p, Vector{Float64}(resolve_cbdy(M, nothing, p)))
        N    = vert.N_dof
        Random.seed!(2026 + N)
        mkv() = VectorValue(rand(N)...)
        Gx, Gy, Ux, Uy, S, DU = mkv(), mkv(), mkv(), mkv(), mkv(), mkv()

        for (fam, arr) in (("K", vert.Kcal), ("P", vert.Pcal))
            T1 = alg_to_tensor3(arr[:, :, :, 1])
            T2 = alg_to_tensor3(arr[:, :, :, 2])
            T5 = alg_to_tensor3(arr[:, :, :, 5])
            W  = alg_to_tensor3(alg_class3_weight(arr[:, :, :, 1],
                                                  arr[:, :, :, 2],
                                                  arr[:, :, :, 5]))

            d = direct_sum(T1, T2, T5, Gx, Gy, Ux, Uy, S, DU)
            r = reduced_sum(W, T2, Gx, Gy, Ux, Uy, S, DU)
            e = relerr(d, r, N)
            gate("G2/G3 $label 𝓣=$fam reduced ≡ direct (rel $(round(e, sigdigits=3)))",
                 e < 1e-13)

            # ---- G4 THE CONTROL: without the (k,j) transpose it MUST fail ----
            Wbad = alg_to_tensor3([-arr[i,k,j,1] + arr[i,k,j,2] - arr[i,k,j,5]
                                   for i in 1:N, k in 1:N, j in 1:N])
            rbad = reduced_sum(Wbad, T2, Gx, Gy, Ux, Uy, S, DU)
            ebad = relerr(d, rbad, N)
            gate("G4 $label 𝓣=$fam CONTROL: no transpose ⇒ detectably wrong " *
                 "(rel $(round(ebad, sigdigits=3)))", ebad > 1e-3)

            # ---- G5: W is a constant, state-independent tensor ---------------
            Wb = alg_to_tensor3(alg_class3_weight(arr[:, :, :, 1],
                                                  arr[:, :, :, 2],
                                                  arr[:, :, :, 5]))
            gate("G5 $label 𝓣=$fam W is deterministic/constant",
                 maximum(abs.([W[i,k,j] - Wb[i,k,j]
                               for i in 1:N, k in 1:N, j in 1:N])) == 0.0)
        end
    end
end

println("test_class3_reduction: $(NGATE[]) gates")

end # module
