using PackageCompiler

# you can use create_sysimage(; kwargs...) if you've activated the project and
# instantiated the environment
create_sysimage([:Images, :FFTW, :StatsBase, :GLMakie];
                sysimage_path = "fresdet." * Libc.Libdl.dlext,
                precompile_execution_file = "precompile/precompile_script.jl")

#=

When the plots appear interact with the figure so that all of the functions are
precompiled and then close the window to continue the precompilation process.

The sysimage can be loaded at start-up using the -J or --sysimage switch.

=#
