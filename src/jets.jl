# ==============================================================
#  jets.jl — value-and-tangent pairs: a residual integrand and its exact Jacobian from ONE formula
#
#  WHY (markdown_files/MIXED_VMS_PLAN.md §4). The mixed and VMS terms are long polynomials in the
#  fields and their derivatives. Writing each once for the residual and once more, by hand, for the
#  Jacobian is how the two drift apart. Here every such term is written once on `Jet`s:
#
#      Jet(v, d)      v = the quantity,   d = its directional derivative (the "tangent")
#
#  The value part is the residual integrand. Seeding the tangents of the unknowns with the trial
#  functions makes the tangent part the Jacobian integrand, by the product rule applied operation
#  by operation. Seeding the tangents of u̇ instead gives ∂R/∂u̇.
#
#  `JZero` is the exact zero. A term that is identically zero (the bed slope on a flat bed, a
#  switched-off coefficient, the tangent of a field that is not being differentiated) is carried
#  as `JZ` and PRUNED: nothing is assembled for it and no operation tree is built.
#
#  Rule 6 holds: nothing here differentiates in space. Spatial derivatives enter as separate jets
#  built from raw fields and raw bases (`alg_dx(u)` with tangent `alg_dx(du)`).
# ==============================================================

"The identically-zero term of the jet algebra (pruned, never assembled)."
struct JZero end
const JZ = JZero()

"""
    Jet(v, d)

A quantity `v` and its directional derivative `d`. Either part is a Gridap CellField expression,
a number, or `JZ` (identically zero). `Jet(v)` has a zero tangent.
"""
struct Jet{V,D}
    v :: V
    d :: D
end
Jet(v) = Jet(v, JZ)

"The zero jet."
jzero() = Jet(JZ, JZ)
jlift(a::Jet) = a
jlift(a)      = Jet(a, JZ)
jval(a::Jet)  = a.v
jtan(a::Jet)  = a.d

# ---- term algebra: a term is `JZ` or anything Gridap can operate on -------------------
tadd(::JZero, ::JZero) = JZ
tadd(a, ::JZero)       = a
tadd(::JZero, b)       = b
tadd(a, b)             = a + b
tsub(::JZero, ::JZero) = JZ
tsub(a, ::JZero)       = a
tsub(::JZero, b)       = (-1.0) * b
tsub(a, b)             = a - b
"Bilinear operation: zero if either factor is zero."
tbil(f, ::JZero, ::JZero) = JZ
tbil(f, ::JZero, b)       = JZ
tbil(f, a, ::JZero)       = JZ
tbil(f, a, b)             = f(a, b)
"Linear operation: zero on zero."
tlin(f, ::JZero) = JZ
tlin(f, a)       = f(a)

# ---- jet algebra ------------------------------------------------------------------------
"Product rule for a bilinear `f`: δf(a,b) = f(δa,b) + f(a,δb)."
jbil(f, a::Jet, b::Jet) = Jet(tbil(f, a.v, b.v), tadd(tbil(f, a.d, b.v), tbil(f, a.v, b.d)))
"A linear `f` acts on value and tangent alike."
jlin(f, a::Jet) = Jet(tlin(f, a.v), tlin(f, a.d))

Base.:+(a::Jet, b::Jet) = Jet(tadd(a.v, b.v), tadd(a.d, b.d))
Base.:+(a::Jet, b)      = a + jlift(b)
Base.:+(a, b::Jet)      = jlift(a) + b
Base.:-(a::Jet, b::Jet) = Jet(tsub(a.v, b.v), tsub(a.d, b.d))
Base.:-(a::Jet, b)      = a - jlift(b)
Base.:-(a, b::Jet)      = jlift(a) - b
Base.:-(a::Jet)         = Jet(tsub(JZ, a.v), tsub(JZ, a.d))
Base.:*(a::Jet, b::Jet) = jbil(*, a, b)
Base.:*(a::Jet, b)      = jbil(*, a, jlift(b))
Base.:*(a, b::Jet)      = jbil(*, jlift(a), b)
LinearAlgebra.dot(a::Jet, b::Jet) = jbil(⋅, a, b)
LinearAlgebra.dot(a::Jet, b)      = jbil(⋅, a, jlift(b))
LinearAlgebra.dot(a, b::Jet)      = jbil(⋅, jlift(a), b)

alg_outer(a::Jet, b::Jet) = jbil(alg_outer, a, b)
alg_outer(a::Jet, b)      = jbil(alg_outer, a, jlift(b))
alg_outer(a, b::Jet)      = jbil(alg_outer, jlift(a), b)
alg_mul(A::TensorValue, u::Jet)           = jlin(x -> alg_mul(A, x), u)
alg_dot(a::VectorValue, u::Jet)           = jlin(x -> alg_dot(a, x), u)
alg_dc3(T::ThirdOrderTensorValue, S::Jet) = jlin(x -> alg_dc3(T, x), S)

"Scalar CellField × constant vector → VectorValue-CellField (`s·a`). Linear in `s`."
alg_svec(a::VectorValue, s)      = Operation(x -> x * a)(s)
alg_svec(a::VectorValue, s::Jet) = jlin(x -> alg_svec(a, x), s)

"""
    jfun(f, df, a::Jet, params...)

A nonlinear scalar function of one jet: value `f(params..., a)`, tangent `df(params..., a)·δa`.
`params` are plain CellFields (no tangent).
"""
function jfun(f, df, a::Jet, params...)
    a.v isa JZero && error("jfun: the argument is identically zero; evaluate the function directly")
    return Jet(Operation(f)(params..., a.v), tbil(*, Operation(df)(params..., a.v), a.d))
end

"Integrate one part of a jet (`pick = jval` or `jtan`); `nothing` when that part is identically zero."
jint(pick, a::Jet, dΩ) = (x = pick(a); x isa JZero ? nothing : ∫(x) * dΩ)

"Accumulate DomainContributions, skipping the `nothing` of a pruned term."
cadd(r, ::Nothing)         = r
cadd(::Nothing, t)         = t
cadd(::Nothing, ::Nothing) = nothing
cadd(r, t)                 = r + t
