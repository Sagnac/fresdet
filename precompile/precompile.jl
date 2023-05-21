using Pkg
using PackageCompiler

Pkg.activate(".")

sysimage_path = "fresdet." * Libc.Libdl.dlext

precompile_execution_file = "precompile/precompile_script.jl"

create_sysimage(; sysimage_path, precompile_execution_file)

#=

When the plots appear interact with the figure so that all of the functions are
precompiled and then close the window to continue the precompilation process.

The sysimage can be loaded at start-up using the -J or --sysimage switch.

=#
