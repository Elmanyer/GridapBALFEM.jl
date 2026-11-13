#!/bin/bash
#SBATCH --job-name="compile_GridapBALFEM"
#SBATCH --partition=compute-p2
#SBATCH --account=innovation
#SBATCH --time=16:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=56G
#SBATCH --output=compile_GridapBALFEM.%j.out
#SBATCH --error=compile_GridapBALFEM.%j.err
#
# Build the GridapBALFEM system image on DelftBlue (the twin of compile_snellius.sh). Submit from
# the compile/ dir:
#     sbatch compile_blue.sh
# Produces  $HOME/GridapBALFEM.jl/GridapBALFEM_sysimage.so  (~1 GB).
# Build on the SAME partition you run on (compute-p2): the image is compiled for the build node's
# CPU (compile.jl sets no cpu_target), so an image built here does not load on compute-p1/memory.
# A different account: sbatch --account=<name> compile_blue.sh.
# Needs Julia 1.10 from juliaup (compile/load_modules_blue.sh says how to install it).

set -euo pipefail
PROJ=$HOME/GridapBALFEM.jl
export BALFEM_CLUSTER=blue
source "$PROJ/compile/load_modules_blue.sh"        # OpenMPI 5.0.10, juliaup Julia 1.10, MPI_LIB_DIR

# 0) The packages of Manifest.toml, PackageCompiler, and 1) MPI.jl pinned to the SYSTEM OpenMPI
#    and Gridap to performance mode (LocalPreferences.toml). These download (GitHub, registry).
#    If the compute nodes have no internet, run them once on the login node and submit with
#    BALFEM_PREFS_DONE=1:
#        source compile/load_modules_blue.sh
#        julia --project=. -e 'using Pkg; Pkg.instantiate(); haskey(Pkg.project().dependencies, "PackageCompiler") || Pkg.add("PackageCompiler")'
#        julia --project=. compile/set_preferences.jl
#        BALFEM_PREFS_DONE=1 sbatch compile_blue.sh
if [ "${BALFEM_PREFS_DONE:-0}" != "1" ]; then
    julia --project="$PROJ" -e 'using Pkg; Pkg.instantiate(); haskey(Pkg.project().dependencies, "PackageCompiler") || Pkg.add("PackageCompiler")'
    julia --project="$PROJ" "$PROJ/compile/set_preferences.jl"
fi

# 1b) Clean precompile with the new preference, then VERIFY that MPI binds the SYSTEM OpenMPI.
julia --project="$PROJ" -e 'using Pkg; Pkg.precompile()'
julia --project="$PROJ" -e '
    using MPIPreferences
    @assert MPIPreferences.binary == "system" "MPIPreferences.binary=$(MPIPreferences.binary) (expected system)"
    using MPI
    v = MPI.MPI_LIBRARY_VERSION_STRING
    println("MPI bound at build: ", v)
    @assert occursin("Open MPI", v) "MPI did not bind system Open MPI: $v"'

# 2) One-time: mpiexecjl in ~/.julia/bin.
command -v mpiexecjl >/dev/null || julia --project="$PROJ" -e 'using MPI; MPI.install_mpiexecjl()'

# 3) Build the sysimage. Single MPI rank so warmup.jl can trace the distributed path too.
mpiexecjl --project="$PROJ" -n 1 julia --project="$PROJ" "$PROJ/compile/compile.jl"

# 4) Stamp the image with a hash of src/*.jl (run/balfem_env.sh, balfem_check_sysimage_freshness).
source "$PROJ/run/balfem_env.sh"
balfem_write_sysimage_stamp "$PROJ/GridapBALFEM_sysimage.so" "$PROJ/src"

echo "Sysimage built: $PROJ/GridapBALFEM_sysimage.so"
