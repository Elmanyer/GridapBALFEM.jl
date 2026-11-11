# ==============================================================
#  set_preferences.jl — pin MPI.jl to the cluster's SYSTEM OpenMPI, Gridap to performance mode
#
#  WHY: without this, MPI.jl uses its bundled JLL binary (MPICH). A sysimage
#  built in that state bakes in MPICH, and launching it under the cluster's
#  OpenMPI mpiexec/PMIx then crashes (implementation mismatch — the "continuous
#  run kills"). This writes LocalPreferences.toml so MPI.jl binds to the system
#  OpenMPI at BOTH build and run time.
#
#  RUN IT ONCE, with the SAME --project used to build the sysimage AND to launch
#  runs (so the binding lands in the right LocalPreferences.toml):
#      source compile/load_modules_snellius.sh
#      julia --project=$HOME/GridapBALFEM.jl compile/set_preferences.jl
#
#  The libmpi path below MUST match the OpenMPI module in load_modules_snellius.sh
#  (Snellius: OpenMPI/5.0.3-GCC-13.3.0). If you change the module, change this.
# ==============================================================

using MPIPreferences, Pkg

### DelftBlue
#mpi_lib_dir = "/apps/arch/2024r1/software/linux-rhel8-cascadelake/gcc-11.3.0/openmpi-4.1.6-w6w5qi5ljesbctyoojlfialbynqt25jb/lib/"

# Snellius OpenMPI/5.0.3-GCC-13.3.0 (matches load_modules_snellius.sh)
mpi_lib_dir = "/sw/arch/RHEL9/EB_production/2024/software/OpenMPI/5.0.3-GCC-13.3.0/lib/"
libmpi = joinpath(mpi_lib_dir, "libmpi.so")

isfile(libmpi) || error(
    "libmpi.so not found in $mpi_lib_dir — load the OpenMPI module first " *
    "(`source compile/load_modules_snellius.sh`) and verify the path.")

# Use the OFFICIAL switch rather than a hand-rolled set_preferences! block. It
# dlopens the system libmpi, auto-detects the vendor/ABI, and writes ALL required
# keys (including the `_format` marker the hand-rolled block was missing — that
# omission makes MPIPreferences silently fall back to the JLL, which is exactly
# how a sysimage ends up baked against OpenMPI_jll). `force=true` overwrites any
# stale block.
MPIPreferences.use_system_binary(;
    library_names = [libmpi],
    mpiexec       = "mpiexec",
    force         = true,
)

# Add WaveSpec from repository, with the fix on change_seed function bug.
Pkg.add(url = "https://github.com/Elmanyer/WaveSpec.jl.git")

# Add Gridap from repository, with the fix for multifield transient calculation of Jacobians using AD.
Pkg.add(url = "https://github.com/Elmanyer/Gridap.jl.git", rev = "fix-transient-multifield-ad")

# Gridap in PERFORMANCE mode (author, 2026-11-11): Gridap's internal `@check` consistency checks
# are compiled away (≈ 20 % of a model-5 step in the 2026-10-31 profile). Real errors (`@assert`,
# `error`, `@notimplemented`) stay and the numerics are unchanged. The verification gates run in
# debug mode on the development machine; production runs in performance mode. Like the MPI
# binding, the mode is baked into the sysimage: set it BEFORE the build, rebuild after changing it.
# BALFEM_GRIDAP_MODE=debug keeps the checks (e.g. to chase an obscure crash on the cluster).
using Preferences
const GRIDAP_UUID = Base.UUID("56d4f2e9-7ea1-5844-9cf6-b9c51ca7ce8e")
gridap_mode = get(ENV, "BALFEM_GRIDAP_MODE", "performance")
gridap_mode in ("performance", "debug") || error("BALFEM_GRIDAP_MODE must be performance or debug")
set_preferences!(GRIDAP_UUID, "execution_mode" => gridap_mode; force = true)
@info "Gridap execution_mode = $gridap_mode (LocalPreferences.toml)"

# Rebuild MPI's precompile cache now so the system binding is in place BEFORE the
# sysimage build bakes it.
Pkg.precompile()
@info "MPIPreferences pinned to system OpenMPI via use_system_binary()" libmpi binary=MPIPreferences.binary abi=MPIPreferences.abi
