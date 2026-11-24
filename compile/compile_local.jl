# ==============================================================
#  compile_local.jl — build the GridapBALFEM system image for THIS WORKSTATION
#
#  The cluster build (compile.jl) refuses anything but the system OpenMPI. A workstation has
#  no LocalPreferences.toml, so MPI.jl uses the BUNDLED (JLL) MPI that ~/.julia/bin/mpiexecjl
#  launches (CLAUDE.md rule 21). This script bakes the same packages and runs the same
#  warmup.jl (the complete solver, sequential and distributed), but bound to that JLL MPI,
#  and writes a DIFFERENTLY NAMED image so the two can never be confused:
#      GridapBALFEM_sysimage_local.so   — this machine, mpiexecjl  (gitignored)
#      GridapBALFEM_sysimage.so         — Snellius, system OpenMPI (compile.jl)
#
#  PREREQUISITE: PackageCompiler in the shared environment (it stacks under --project):
#      julia --project=@v1.10 -e 'using Pkg; Pkg.add("PackageCompiler")'
#  BUILD (single rank, so the distributed warmup cases can initialise MPI):
#      ~/.julia/bin/mpiexecjl --project=. -n 1 julia --project=. compile/compile_local.jl
#  USE:
#      julia --project=. -J GridapBALFEM_sysimage_local.so script.jl
#      ~/.julia/bin/mpiexecjl --project=. -n 4 julia --project=. -J GridapBALFEM_sysimage_local.so script.jl
#  ⚠ Valid only for the src/ it was built from (rule 25): rebuild after any src/ edit.
# ==============================================================

using MPIPreferences
MPIPreferences.binary == "system" && error("""
    MPIPreferences.binary = "system": this project is pinned to a system MPI, which is the
    CLUSTER configuration. Use compile/compile.jl there; this script is for the JLL MPI.""")
using MPI
@info "MPI preflight: baking against the bundled MPI (mpiexecjl)" MPIPreferences.binary version = MPI.MPI_LIBRARY_VERSION_STRING

let id = Base.identify_package("GridapBALFEM")
    id === nothing && error("GridapBALFEM is not resolvable from $(Base.active_project()); " *
                            "build with --project=<GridapBALFEM.jl>")
    @info "package preflight OK — the solver will be baked into the image" uuid = id.uuid
end
using GridapBALFEM

using PackageCompiler
const IMAGE = joinpath(@__DIR__, "..", "GridapBALFEM_sysimage_local.so")
#  include_transitive_dependencies=false as in compile.jl: bake what the solver actually loads.
create_sysimage(
    [:GridapBALFEM, :Gridap, :GridapDistributed, :GridapSolvers, :PartitionedArrays,
     :MPI, :BlockArrays, :WaveSpec];
    sysimage_path                   = IMAGE,
    precompile_execution_file       = joinpath(@__DIR__, "warmup.jl"),
    include_transitive_dependencies = false,
)
@info "local sysimage written (compile_local.sh stamps it with the src/ hash)" IMAGE
