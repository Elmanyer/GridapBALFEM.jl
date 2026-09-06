# ==============================================================
#  horizontal.jl — Stage 2: 2D mesh and the STACKED FE spaces
#
#  Builds the horizontal discretisation the time loop solves on: the 2D mesh and
#  the MultiFieldFESpace with THREE fields
#      [η, 𝖴x, 𝖴y],   𝖴x,𝖴y ∈ VectorValue{Nσ}
#  Each velocity field carries all Nσ vertical modes in one vector-valued FE
#  function, so a solid-wall Dirichlet BC constrains the whole stacked field
#  (all Nσ components zero at once) with a single condition.
#
#  Boundary tags (CartesianDiscreteModel, 2D):
#    tag_1..tag_4 = corners, tag_5/6 = bottom/top edges, tag_7/8 = left/right.
#  Corner tags are MANDATORY in any wall BC (omitting them leaves corner DOFs
#  unconstrained → exponential instability; root CLAUDE.md rule 5).
#
#  All fields H1-conforming Lagrange, p_u ≥ 2 (linear elements zero the
#  dispersion term and disable all non-hydrostatic physics).
# ==============================================================

"""
    build_horizontal_model(domain, partition; y_periodic=false) → (model, trian)

2D Cartesian mesh on ((x0,x1),(y0,y1)) (or flat (x0,x1,y0,y1)) with
`partition = (nx, ny)` cells. Create the measure as `Measure(trian, 2*p_u+2)`.
`y_periodic=true` builds the mesh periodic in y (the matching top/bottom edge DOFs
are identified), for the `:periodic` lateral boundary condition.
"""
function build_horizontal_model(domain::Tuple, partition::Tuple; y_periodic::Bool=false)
    if domain isa Tuple{Tuple,Tuple}
        (x0,x1), (y0,y1) = domain
        dom_flat = (x0, x1, y0, y1)
    else
        dom_flat = domain
    end
    model = CartesianDiscreteModel(dom_flat, partition; isperiodic=(false, y_periodic))
    trian = Triangulation(model)
    return model, trian
end

"Boundary tags of one side of the CartesianDiscreteModel (corners INCLUDED)."
function side_tags(side::Symbol)
    side == :left   && return ["tag_1","tag_3","tag_7"]
    side == :right  && return ["tag_2","tag_4","tag_8"]
    side == :bottom && return ["tag_1","tag_2","tag_5"]
    side == :top    && return ["tag_3","tag_4","tag_6"]
    error("side_tags: side must be :left, :right, :bottom or :top (got :$side)")
end

"""
    check_taylor_hood(p_u, p_eta; where="")

**THE horizontal element-pairing gate. Every path that builds FE spaces goes through
it, and it RAISES rather than warns.**

The BALFE-M horizontal pairing is **Taylor-Hood by design**: `p_u = p_eta + 1`. This
is not a preference, it is a well-posedness requirement, and it is enforced here
because getting it wrong is silent, slow and expensive.

**Why.** `η` enters the momentum equation *undifferentiated* — it reaches the test
function only through `∇·v`, after the integration by parts that produces `R_P`. It
therefore plays exactly the role pressure plays in a Stokes system, and the pairing is
subject to the same **inf-sup (LBB)** condition. Equal-order continuous spaces
(`p_u = p_eta`) are inf-sup deficient: they admit a spurious checkerboard mode at
`λ ≈ 2·dx` which the discrete operator cannot see and refinement does not remove.

**What that cost, before this gate existed.** Every fully nonlinear run was launched
at equal order while the entire MMS verification campaign ran Taylor-Hood, so the
configuration that blew up was never the configuration that was verified. The result
looked like a physical instability of the model: an unbounded grid-scale mode that got
*worse* under refinement (blow-up at t≈24 s for `dx=0.25`, t≈10 s for `dx=0.125`). Ten
hypotheses were refuted against it — Jacobian, sponge, boundaries, domain length, CFL,
Benjamin–Feir, quadrature aliasing, advection energy conservation — over three days.
On Taylor-Hood the mode does not exist and the refinement signature inverts.
See CLAUDE.md rules 2b and 12b.

`p_eta ≥ 1` is required as well: `p_u = p_eta + 1 ≥ 2` then follows, which is rule 2
(`Q1` velocity zeroes `R_P` and disables all non-hydrostatic physics).
"""
function check_taylor_hood(p_u::Int, p_eta::Int; where::AbstractString = "")
    ctx = isempty(where) ? "" : "$where: "
    p_eta ≥ 1 || error(
        ctx * "p_eta = $p_eta is invalid — the surface order must be ≥ 1.\n" *
        "  A Q0 surface cannot represent ∇η at all.")
    p_u == p_eta + 1 || error(
        ctx * "NON-TAYLOR-HOOD horizontal pairing: p_u = $p_u, p_eta = $p_eta.\n" *
        "  BALFE-M requires p_u = p_eta + 1 (velocity exactly one order above the surface).\n" *
        (p_u == p_eta ?
          "  You asked for EQUAL ORDER, which is inf-sup deficient here: η plays the\n" *
          "  pressure role of a Stokes system, and equal-order continuous spaces admit a\n" *
          "  spurious checkerboard at λ ≈ 2·dx. This is the defect that produced the\n" *
          "  'nonlinear instability' of 2026-09 (CLAUDE.md rule 12b).\n" :
          "  A gap other than one is not a supported pairing; nothing in the verified\n" *
          "  scope was measured on it.\n") *
        "  Use p_u = $(p_eta + 1) with p_eta = $p_eta, or p_u = $p_u with p_eta = $(p_u - 1).")
    return nothing
end

"""
    build_fe_spaces(model, p_u, Nσ; y_wall_bc=:wall, x_wall_bc=false,
                        inflow=nothing) → (U, V)

Stacked 3-field MultiFieldFESpace `[η, 𝖴x, 𝖴y]`. The lateral (y) boundary
condition is selected by the symbol `y_wall_bc`:

- `:wall` (default): solid wall `u_j^y = 0 ∀j` on the y-boundaries — Dirichlet
  zero `VectorValue{Nσ}` on the whole 𝖴y field (tags 1–6).
- `:open`: natural (free / zero-flux) y-edges — no essential condition.
- `:periodic`: the y-edges are identified in the *mesh* (`build_horizontal_model`
  with `y_periodic=true`); at the FE-space level this behaves like `:open`
  (no y-Dirichlet), the DOF identification being carried by the periodic mesh.

Other arguments:
- `x_wall_bc=true`: closed basin — additionally `u_j^x = 0 ∀j` on x-boundaries
  (tags 1–4,7,8). REQUIRED for any initial-condition problem (soliton, sloshing,
  IC hump); flumes keep it false (open x-ends, sponge-absorbed).
- `inflow`: Dirichlet wave generation on one boundary —
  `(side=:left, eta=g_eta, ux=g_ux, uy=g_uy_or_nothing)` with `g(t) = x -> …`
  transient closures (from `eta_bc`/`ux_bc`/`uy_bc` of a `WaveInput`). η and
  𝖴x become `TransientTrialFESpace`s (the returned `U` is callable at `t`);
  `uy !== nothing` (directional sea) REQUIRES `y_wall_bc ≠ :wall` (`:open` or
  `:periodic`). With `x_wall_bc=true` the opposite x-side keeps its zero wall.
  Sides `:left` / `:right` are supported for generation.
"""

function build_fe_spaces(model, p_u::Int, Nσ::Int;
                             y_wall_bc::Symbol = :wall, x_wall_bc::Bool = false,
                             inflow = nothing,
                             p_eta::Int = p_u - 1)
    #  THE gate. Every FE space in this solver is built here, so validating the pairing
    #  at this one point makes a non-Taylor-Hood run impossible rather than merely
    #  discouraged. Default p_eta = p_u − 1 is Taylor-Hood by construction.
    check_taylor_hood(p_u, p_eta; where = "build_fe_spaces")
    y_wall_bc in (:wall, :open, :periodic) ||
        error("build_fe_spaces: y_wall_bc must be :wall, :open or :periodic (got :$y_wall_bc)")
    reffe_eta = ReferenceFE(lagrangian, Float64, p_eta)
    reffe_U   = ReferenceFE(lagrangian, VectorValue{Nσ,Float64}, p_u)
    zvv       = VectorValue(ntuple(_ -> 0.0, Nσ)...)
    # Dirichlet y-wall only for :wall; :open and :periodic impose no y-essential BC
    # (the periodic case identifies the y-edge DOFs at the mesh level instead).
    y_tags = y_wall_bc == :wall ? ["tag_1","tag_2","tag_3","tag_4","tag_5","tag_6"] : String[]
    x_tags = x_wall_bc ? ["tag_1","tag_2","tag_3","tag_4","tag_7","tag_8"] : String[]

    if inflow === nothing
        V_eta = FESpace(model, reffe_eta; conformity=:H1)
        U_eta = TrialFESpace(V_eta)
        V_Ux  = isempty(x_tags) ? FESpace(model, reffe_U; conformity=:H1) :
                                  FESpace(model, reffe_U; conformity=:H1, dirichlet_tags=x_tags)
        U_Ux  = isempty(x_tags) ? TrialFESpace(V_Ux) : TrialFESpace(V_Ux, zvv)
    else
        side = inflow.side
        side in (:left, :right) ||
            error("build_fe_spaces: wave-generation side must be :left or :right (got :$side)")
        gen_tags = side_tags(side)
        zfun     = t -> (x -> zvv)                       # static wall in transient format

        # η: prescribed on the generation side only
        V_eta = FESpace(model, reffe_eta; conformity=:H1, dirichlet_tags=gen_tags)
        U_eta = TransientTrialFESpace(V_eta, inflow.eta)

        # 𝖴x: generation side prescribed; opposite side stays a wall if x_wall_bc.
        # CONTRACT: TransientTrialFESpace(V, funs) pairs funs[i] with
        # dirichlet_tags[i] — the two vectors below must stay index-aligned,
        # and static walls must be wrapped in transient form (t -> x -> 0).
        # The corner tags of the two x-sides are disjoint ({1,3} vs {2,4}),
        # so no tag receives two different functions.
        if x_wall_bc
            opp_tags = side_tags(side == :left ? :right : :left)
            ux_tags  = vcat(gen_tags, opp_tags)
            ux_funs  = vcat(fill(inflow.ux, length(gen_tags)),
                            fill(zfun,      length(opp_tags)))
        else
            ux_tags = gen_tags
            ux_funs = fill(inflow.ux, length(gen_tags))
        end
        V_Ux = FESpace(model, reffe_U; conformity=:H1, dirichlet_tags=ux_tags)
        U_Ux = TransientTrialFESpace(V_Ux, ux_funs)
    end

    if inflow !== nothing && inflow.uy !== nothing
        # directional sea: v ≠ 0 on the generation side is incompatible with the
        # y-wall corner tags — lateral boundaries must be :open (sponge-absorbed)
        # or :periodic instead
        y_wall_bc == :wall &&
            error("build_fe_spaces: a directional inflow (uy prescribed) requires " *
                  "y_wall_bc=:open or :periodic (not :wall); the wall corner tags " *
                  "conflict with v ≠ 0 at the generation boundary")
        gen_tags = side_tags(inflow.side)
        V_Uy = FESpace(model, reffe_U; conformity=:H1, dirichlet_tags=gen_tags)
        U_Uy = TransientTrialFESpace(V_Uy, inflow.uy)
    else
        V_Uy  = isempty(y_tags) ? FESpace(model, reffe_U; conformity=:H1) :
                                  FESpace(model, reffe_U; conformity=:H1, dirichlet_tags=y_tags)
        U_Uy  = isempty(y_tags) ? TrialFESpace(V_Uy) : TrialFESpace(V_Uy, zvv)
    end

    # The element type of the space vectors must be an abstract FE-space supertype
    # so that a mixed-type vector (unconstrained + Dirichlet fields) is accepted by
    # the sequential Vector{<:SingleFieldFESpace} resp. distributed
    # Vector{<:DistributedSingleFieldFESpace} MultiFieldFESpace dispatch.
    # TransientTrialFESpace is NOT a SingleFieldFESpace — with an inflow the trial
    # vector is left untyped so the MultiFieldFESpace constructor can promote it,
    # and the resulting space is callable at t. Keep ConsecutiveMultiFieldStyle:
    # BlockMultiFieldStyle breaks the distributed Jacobi preconditioner's `diag`.
    if isa(V_eta, GridapDistributed.DistributedSingleFieldFESpace)
        # distributed transient trials ARE DistributedSingleFieldFESpaces
        # (GridapDistributed src/ODEs.jl) — the typed vector stays valid
        Vspaces = GridapDistributed.DistributedSingleFieldFESpace[V_eta, V_Ux, V_Uy]
        Uspaces = GridapDistributed.DistributedSingleFieldFESpace[U_eta, U_Ux, U_Uy]
    else
        Vspaces = Gridap.FESpaces.SingleFieldFESpace[V_eta, V_Ux, V_Uy]
        Uspaces = inflow === nothing ?
                  Gridap.FESpaces.SingleFieldFESpace[U_eta, U_Ux, U_Uy] :
                  [U_eta, U_Ux, U_Uy]
    end
    U = MultiFieldFESpace(Uspaces)
    V = MultiFieldFESpace(Vspaces)
    return U, V
end
