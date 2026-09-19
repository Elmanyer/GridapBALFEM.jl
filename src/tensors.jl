# ==============================================================
#  tensors.jl — constant-tensor constructors + pointwise algebra helpers
#
#  Constructors (build time only). Gridap TensorValue data is column-major
#  (first index fastest):
#    TensorValue{N,N}(data...)            → T[i,j]   = data[i + (j-1)N]
#    ThirdOrderTensorValue{N,N,N}(data..) → T[i,j,k] = data[i + (j-1)N + (k-1)N²]
#  (verified by test/test_primitives.jl)
#
#  Pointwise helpers wrap constants in Operation closures (closing over
#  constants is fine — only CellFields must not be closed over in loops).
#  Verified semantics (see markdown_files/ARCHITECTURE.md §1):
#    ∇ of a VectorValue{Nσ} field = TensorValue{2,Nσ} (spatial index FIRST)
#      ⇒ ∂x f = e_x ⋅ ∇f   (use e⋅∇f, never ∇f⋅e)
#    double_contraction(𝗧3, S) contracts the TRAILING two indices:
#      result_i = Σ_{k,j} 𝗧3[i,k,j] S[k,j]
# ==============================================================

"Vector{Float64} → VectorValue{N} constant."
alg_to_vec(v::AbstractVector{Float64}) = VectorValue{length(v),Float64}(v...)

"Matrix → TensorValue{N,N} constant (index-preserving)."
function alg_to_tensor2(M::AbstractMatrix{Float64})
    N = size(M, 1); @assert size(M, 2) == N
    data = ntuple(k -> M[(k-1) % N + 1, (k-1) ÷ N + 1], N * N)
    return TensorValue{N,N,Float64}(data...)
end

"Array{Float64,3} → ThirdOrderTensorValue{N,N,N} constant (index-preserving)."
function alg_to_tensor3(T::AbstractArray{Float64,3})
    N = size(T, 1); @assert size(T, 2) == N && size(T, 3) == N
    data = ntuple(N^3) do k
        k0 = k - 1
        i = k0 % N + 1
        j = (k0 ÷ N) % N + 1
        l = k0 ÷ (N * N) + 1
        T[i, j, l]
    end
    return ThirdOrderTensorValue{N,N,N,Float64}(data...)
end

"∂x of a scalar or VectorValue{Nσ} CellField: e_x ⋅ ∇f (spatial index first)."
alg_dx(f) = Operation(g -> Ex ⋅ g)(∇(f))

"∂y of a scalar or VectorValue{Nσ} CellField: e_y ⋅ ∇f."
alg_dy(f) = Operation(g -> Ey ⋅ g)(∇(f))

"Constant TensorValue{N,N} ⋅ VectorValue{N}-CellField → VectorValue{N}-CellField (matvec)."
alg_mul(A::TensorValue, u) = Operation(x -> A ⋅ x)(u)

"Constant VectorValue{N} ⋅ VectorValue{N}-CellField → scalar CellField."
alg_dot(a::VectorValue, u) = Operation(x -> a ⋅ x)(u)

"double_contraction(constant ThirdOrderTensorValue, TensorValue-CellField):
contracts the TRAILING two indices → VectorValue{N}-CellField."
alg_dc3(T::ThirdOrderTensorValue, S) =
    Operation(s -> Gridap.TensorValues.double_contraction(T, s))(S)

"Outer product of two VectorValue{N}-CellFields → TensorValue{N,N}-CellField,
(a⊗b)[k,j] = a[k]·b[j]."
alg_outer(a, b) = Operation(Gridap.TensorValues.outer)(a, b)

"Fuse two scalar CellFields into a VectorValue{2}-CellField."
alg_vec2(a, b) = Operation(VectorValue)(a, b)

"""
    alg_class3_weight(T1, T2, T5) -> Array{Float64,3}

Combined weight for the EXACT algebraic reduction of the Class-III components
`{1,2,5}` of `𝓝ₖⱼ` onto a single contraction (markdown_files/NEW_TREATMENT.md §A.2).

With `Gₐ = ∂ₐπ(𝖲)` and the code convention `(a⊗b)[k,j] = a[k]b[j]`,

    𝓝¹ = −Σₐ Gₐ⊗Uₐ        𝓝² = −𝓝¹ + 𝖲⊗DU        𝓝⁵ = −Σₐ Uₐ⊗Gₐ

and relabelling the dummy pair `k↔j` in the `𝓝⁵` term gives

    Σ_{ℓ∈{1,2,5}} T⁽ˡ⁾ ⊙ 𝓝⁽ˡ⁾  =  W ⊙ (Σₐ Gₐ⊗Uₐ)  +  T⁽²⁾ ⊙ (𝖲⊗DU)
    W[i,k,j] = −T1[i,k,j] + T2[i,k,j] − T5[i,j,k]

⚠ **`T5` enters TRANSPOSED IN ITS LAST TWO INDICES** — that transpose IS the content of
"𝓝⁵ is the (k,j)-transpose of 𝓝¹". The `𝓐/𝓚/𝓟` families are NOT symmetric in `(k,j)`,
so dropping it is a silent, plausible-looking error. `test_class3_reduction.jl` gate G4
is the control that fails if it is dropped.

Takes and returns plain arrays (called once, at problem construction).
"""
function alg_class3_weight(T1::AbstractArray{Float64,3},
                           T2::AbstractArray{Float64,3},
                           T5::AbstractArray{Float64,3})
    N = size(T1, 1)
    @assert size(T1) == size(T2) == size(T5) == (N, N, N)
    W = Array{Float64,3}(undef, N, N, N)
    for i in 1:N, k in 1:N, j in 1:N
        W[i, k, j] = -T1[i, k, j] + T2[i, k, j] - T5[i, j, k]
    end
    return W
end
